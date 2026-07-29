#!/usr/bin/env python3
"""
######################################################################################################
#                                                                                                    #
# go.py - GO term id -> name / namespace mapping from the Gene Ontology (go-basic.obo).              #
#                                                                                                    #
# The InterProScan TSV this pipeline parses carries GO identifiers only (GO:0004672). go-basic.obo   #
# is fetched once into db/go/ (a copy already there is reused, so runs can be offline) and each      #
# [Term] stanza is parsed for the three fields that matter: id, name, namespace. compile then adds   #
# readable names and splits the terms by aspect (molecular_function / biological_process /           #
# cellular_component). If the ontology cannot be fetched, the run continues with ids only.           #
#                                                                                                    #
######################################################################################################
"""

from __future__ import annotations

import re

from . import external, net

# GO aspects, in the conventional reporting order.
ASPECTS = ("molecular_function", "biological_process", "cellular_component")

# InterProScan's go_terms field appends a source annotation to each id, e.g.
# "GO:0005515(InterPro)|GO:0004553(PANTHER)" -- strip it so the bare id both
# matches go-basic.obo and is what the pipeline reports.
_GO_ID = re.compile(r"GO:\d{7}")


def ensure_obo(cfg):
    """Return db/go/go-basic.obo, downloading it once if not already present."""
    dest = cfg.go_obo
    if dest.exists():
        external.log(f"[OK] GO ontology present: {dest.name}")
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    external.log(f"[..] Fetching GO ontology: {cfg.GO_OBO_URL}")
    # net.fetch sends a User-Agent (purl.obolibrary.org 403s requests without one)
    # and retries transient failures with backoff.
    dest.write_bytes(net.fetch(cfg.GO_OBO_URL, timeout=180))
    external.log(f"[OK] GO ontology cached: {dest}")
    return dest


def parse_obo(path) -> dict:
    """{GO id -> {'name': ..., 'namespace': ...}} from a go-basic.obo file.

    Only ``[Term]`` stanzas contribute: a stanza header opens a new record (and
    ``[Typedef]`` opens none, so ids like ``part_of`` never leak in); the id, name
    and namespace lines within it are collected in any order and stored when the
    stanza closes."""
    terms, cur = {}, None

    def _store(rec):
        if rec and rec.get("id") and rec.get("name"):
            terms[rec["id"]] = {"name": rec["name"], "namespace": rec.get("namespace", "")}

    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("["):
                _store(cur)
                cur = {} if line.strip() == "[Term]" else None
            elif cur is not None:
                if line.startswith("id: GO:"):
                    cur["id"] = line[4:].strip()
                elif line.startswith("name:"):
                    cur["name"] = line[5:].strip()
                elif line.startswith("namespace:"):
                    cur["namespace"] = line[10:].strip()
    _store(cur)
    return terms


def load_terms(cfg) -> dict:
    """GO id -> {'name', 'namespace'}, fetching the ontology if needed. Returns {}
    with a warning if it cannot be obtained, so a run never fails on GO names."""
    try:
        return parse_obo(ensure_obo(cfg))
    except Exception as e:  # noqa: BLE001 - names are optional enrichment
        external.log(f"[WARN] GO names unavailable ({e}); goTerms will carry ids only.")
        return {}


def clean_ids(go_field) -> list:
    """GO ids from a raw go_terms field, source annotations like '(InterPro)'
    stripped so each entry is a bare 'GO:0000000' id."""
    raw = [t.strip() for t in str(go_field).split("|")
           if t.strip() and t.strip() not in ("-", "nan")]
    ids = []
    for t in raw:
        m = _GO_ID.match(t)
        ids.append(m.group(0) if m else t)
    return ids


def names_for(go_field, terms) -> str:
    """'GO:0004672(InterPro)|GO:0005537(PANTHER)' -> 'protein kinase activity|mannose binding'.

    Order is preserved; an id missing from the ontology falls back to the bare id."""
    ids = clean_ids(go_field)
    if not ids or not terms:
        return "-"
    return "|".join(terms[i]["name"] if i in terms else i for i in ids)


def names_by_aspect(go_field, terms) -> dict:
    """Split a GO field into {aspect: 'name|name'} for the three GO aspects."""
    grouped = {aspect: [] for aspect in ASPECTS}
    for go_id in clean_ids(go_field):
        term = terms.get(go_id)
        if term and term.get("namespace") in grouped:
            grouped[term["namespace"]].append(term["name"])
    return {aspect: ("|".join(names) if names else "-") for aspect, names in grouped.items()}
