#!/usr/bin/env python3
"""
####################################################################################################
#                                                                                                  #
# test_interpro.py - InterProScan TSV header-assignment tests.                                     #
#                                                                                                  #
# The GO and pathways columns are optional (present only when --goterms / --pathways are on), so   #
# the header must be named by the actual column count.                                             #
#                                                                                                  #
####################################################################################################
"""

from gwiscan.features.interpro import (
    BASE_COLUMNS,
    _appl_list,
    _combine_gff,
    _header_for,
)


def test_appl_list_default_pfam():
    assert _appl_list("Pfam") == ["Pfam"]


def test_appl_list_cdd_optin():
    # CDD opt-in: comma-separated -> one appl per member db.
    assert _appl_list("Pfam,CDD") == ["Pfam", "CDD"]
    assert _appl_list("Pfam, CDD ") == ["Pfam", "CDD"]


def test_appl_list_accepts_list():
    assert _appl_list(["Pfam", "CDD", ""]) == ["Pfam", "CDD"]


# --- only Pfam-accession families are sent to InterProScan --------------------

from gwiscan.config import Config
from gwiscan.features.interpro import pfam_candidate_ids


def test_write_outputs_splits_by_member_db(tmp_path):
    from gwiscan.features.interpro import _write_outputs
    lines = [
        "P1\tmd5\t100\tPfam\tPF00139\td\t1\t50\t1e-9\tT\tdate\tIPR\td",
        "P1\tmd5\t100\tCDD\tcd00001\td\t1\t50\t1e-9\tT\tdate\tIPR\td",
        "P2\tmd5\t120\tPfam\tPF01453\td\t5\t80\t1e-8\tT\tdate\tIPR\td",
    ]
    _write_outputs(lines, [], tmp_path / "interproscan.tsv", tmp_path / "interproscan.gff3")

    assert (tmp_path / "interproscan.tsv").exists()               # combined
    pfam = tmp_path / "interproscan.pfam.tsv"
    cdd = tmp_path / "interproscan.cdd.tsv"
    assert sum(1 for _ in open(pfam)) == 3                        # header + 2 Pfam rows
    assert sum(1 for _ in open(cdd)) == 2                         # header + 1 CDD row


def test_pfam_candidate_ids_excludes_custom_and_blast_only(tmp_path):
    (tmp_path / "config").mkdir()
    (tmp_path / "config" / "family.tsv").write_text(
        "Family\tPfamModel\tBlastModel\n"
        "GNA\tPF01453\tgna.fasta\n"      # Pfam accession -> InterProScan
        "CRA\tCRA.hmm\tcra.fasta\n"      # custom HMM -> no InterProScan
        "EUL\t-\teul.fasta\n"            # BLAST-only -> no InterProScan
    )
    cfg = Config(root=tmp_path)
    # cfg.result() routes candidates_merged.tsv into intermediate/candidates/.
    cfg.result("candidates_merged.tsv").write_text(
        "proteinId\tfamily\taccession\tevalue\tbitscore\tstart\tend\tmethod\n"
        "P_gna\tGNA\tPF01453.1\t1e-40\t130\t10\t110\thmm\n"
        "P_cra\tCRA\t-\t1e-30\t100\t25\t300\thmm\n"
        "P_eul\tEUL\t-\t1e-20\t80\t5\t200\tblast\n"
    )
    assert pfam_candidate_ids(cfg) == {"P_gna"}


# --- local InterProScan command builder --------------------------------------

from gwiscan.config import Config
from gwiscan.features.interpro import _local_cmd


def test_local_cmd_offline_by_default(tmp_path):
    cfg = Config(root=tmp_path, INTERPROSCAN_BIN="interproscan.sh", THREADS=8)
    cmd = _local_cmd(cfg, tmp_path / "candidates.fasta", tmp_path / "out")
    assert cmd[0] == "interproscan.sh"
    assert "-f" in cmd and cmd[cmd.index("-f") + 1] == "TSV,GFF3"
    # No family table here, so the derivation falls back to Pfam.
    assert cmd[cmd.index("-appl") + 1] == "Pfam"
    assert cmd[cmd.index("-cpu") + 1] == "8"
    assert "-goterms" in cmd
    assert "-dp" in cmd                       # offline by default


def test_local_cmd_lookup_enabled_drops_dp(tmp_path):
    cfg = Config(root=tmp_path, INTERPRO_LOOKUP=True)
    cmd = _local_cmd(cfg, tmp_path / "candidates.fasta", tmp_path / "out")
    assert "-dp" not in cmd                    # online precalc lookup allowed

# A real InterProScan 6 row (pathways empty). 15 fields.
REAL_15 = ("Chr5g0402091\ta7e74c3153434db6cd2ce722d4275c30\t353\tPfam\tPF00069\t"
           "Protein kinase domain\t6\t265\t6.7E-52\tT\t20-06-2026\tIPR000719\t"
           "Protein kinase domain\t-\t-")


def test_base_only():
    assert _header_for(13) == BASE_COLUMNS


def test_with_goterms_only():
    # goterms on, pathways off -> 14 columns ending in go_terms
    h = _header_for(14)
    assert h[-1] == "go_terms"
    assert "pathways" not in h
    assert len(h) == 14


def test_with_goterms_and_pathways():
    h = _header_for(15)
    assert h[-2:] == ["go_terms", "pathways"]


def test_matches_real_row():
    row = REAL_15.split("\t")
    header = _header_for(len(row))
    assert len(header) == len(row)
    mapped = dict(zip(header, row))
    assert mapped["protein_id"] == "Chr5g0402091"
    assert mapped["analysis"] == "Pfam"
    assert mapped["ipr_acc"] == "IPR000719"
    assert mapped["ipr_desc"] == "Protein kinase domain"
    assert mapped["go_terms"] == "-"


def test_extra_columns_padded():
    # Never drop data if the format ever grows.
    h = _header_for(17)
    assert len(h) == 17
    assert h[15:] == ["col_16", "col_17"]


CHUNK_A = (
    "##gff-version 3\n"
    "Chr5g0402091\tPfam\tprotein_match\t6\t265\t6.7E-52\t+\t.\tName=PF00069\n"
    "##FASTA\n"
    ">Chr5g0402091\nMKT\n"
)
CHUNK_B = (
    "##gff-version 3\n"
    "Chr2g0319311\tPfam\tprotein_match\t220\t515\t4.4E-70\t+\t.\tName=PF01501\n"
    "##FASTA\n"
    ">Chr2g0319311\nAAA\n"
)


def test_combine_gff_merges_chunks():
    merged = _combine_gff([CHUNK_A, CHUNK_B])
    lines = merged.splitlines()
    # exactly one version pragma, at the top
    assert lines[0] == "##gff-version 3"
    assert lines.count("##gff-version 3") == 1
    # both features present
    assert any("PF00069" in ln for ln in lines)
    assert any("PF01501" in ln for ln in lines)
    # FASTA sections dropped
    assert "##FASTA" not in merged
    assert ">Chr5g0402091" not in merged


def test_combine_gff_empty():
    assert _combine_gff([]) == "##gff-version 3\n"


# --- network resilience: retrying transient EBI API failures ------------------

from gwiscan.features import interpro


def test_session_is_configured_with_retry_backoff():
    # A transient blip during a multi-hour poll shouldn't discard the run, so the
    # session's adapters carry a Retry with backoff covering POST (submit) too.
    interpro._session.cache_clear()
    session = interpro._session()
    retry = session.get_adapter("https://www.ebi.ac.uk").max_retries

    assert retry.total == interpro._RETRY_TOTAL
    assert retry.backoff_factor == interpro._RETRY_BACKOFF
    for status in (429, 500, 502, 503, 504):
        assert status in retry.status_forcelist
    # allowed_methods=None -> retry every method, including the POST submit
    assert retry.allowed_methods is None
    # let the caller's raise_for_status() produce the clean HTTP error message
    assert retry.raise_on_status is False


def test_session_is_cached_singleton():
    interpro._session.cache_clear()
    assert interpro._session() is interpro._session()


def test_submit_and_getters_go_through_the_retrying_session(monkeypatch):
    # Lock that the four network calls use _session() (which retries), not bare
    # requests.get/post (which don't).
    calls = {"post": 0, "get": 0}

    class _Resp:
        text = "JOBID"
        def raise_for_status(self):  # noqa: D401 - stub
            pass

    class _FakeSession:
        def post(self, *a, **k):
            calls["post"] += 1
            return _Resp()
        def get(self, *a, **k):
            calls["get"] += 1
            return _Resp()

    monkeypatch.setattr(interpro, "_session", lambda: _FakeSession())

    base = interpro.IPRSCAN_ENDPOINTS[5]
    assert interpro._submit(base, "e@x.z", ">s\nACD", ["Pfam"]) == "JOBID"
    interpro._result_types(base, "JOBID")
    interpro._fetch_result(base, "JOBID", "tsv")
    assert calls["post"] == 1
    assert calls["get"] == 2


# --- data integrity: an unfinished/incomplete chunk must not silently drop members ---


def _api_cfg(tmp_path):
    from gwiscan.config import Config
    return Config(root=tmp_path, EBI_EMAIL="e@x.z")


def _write_cand_fasta(tmp_path, n):
    fasta = tmp_path / "cand.fasta"
    fasta.write_text("".join(f">P{i}\nACDEFGHIK\n" for i in range(n)))
    return fasta


def test_run_api_raises_when_a_chunk_does_not_finish(tmp_path, monkeypatch):
    # A chunk that times out / fails leaves its sequences unannotated; confirm would
    # then drop those members silently. The stage must stop instead.
    monkeypatch.setattr(interpro, "_submit", lambda *a, **k: "JOB1")
    monkeypatch.setattr(interpro, "_poll", lambda *a, **k: "TIMEOUT")
    monkeypatch.setattr(interpro, "_result_types", lambda *a, **k: ["tsv", "gff"])
    monkeypatch.setattr(interpro, "_fetch_result", lambda *a, **k: "")

    import pytest
    with pytest.raises(RuntimeError, match="did not finish"):
        interpro._run_api(_api_cfg(tmp_path), _write_cand_fasta(tmp_path, 3))


def test_run_api_raises_when_finished_job_has_no_tsv(tmp_path, monkeypatch):
    # FINISHED but no TSV result type is the same silent-drop, so it is also fatal.
    monkeypatch.setattr(interpro, "_submit", lambda *a, **k: "JOB1")
    monkeypatch.setattr(interpro, "_poll", lambda *a, **k: "FINISHED")
    monkeypatch.setattr(interpro, "_result_types", lambda *a, **k: ["gff"])   # no tsv
    monkeypatch.setattr(interpro, "_fetch_result", lambda *a, **k: "")

    import pytest
    with pytest.raises(RuntimeError, match="no 'tsv' result type"):
        interpro._run_api(_api_cfg(tmp_path), _write_cand_fasta(tmp_path, 2))


def test_run_api_succeeds_when_all_chunks_finish(tmp_path, monkeypatch):
    # The happy path still returns the annotations; a missing GFF is only a warning.
    monkeypatch.setattr(interpro, "_submit", lambda *a, **k: "JOB1")
    monkeypatch.setattr(interpro, "_poll", lambda *a, **k: "FINISHED")
    monkeypatch.setattr(interpro, "_result_types", lambda *a, **k: ["tsv"])   # no gff
    row = "P0\tmd5\t9\tPfam\tPF00139\td\t1\t9\t1e-9\tT\tdate\tIPR\td"
    monkeypatch.setattr(interpro, "_fetch_result", lambda *a, **k: row)
    # Avoid a real XML fetch for the per-chunk version; keep it a single version.
    monkeypatch.setattr(interpro, "_job_version", lambda *a, **k: "5.78-109.0")

    lines, gff_texts, meta = interpro._run_api(_api_cfg(tmp_path), _write_cand_fasta(tmp_path, 1))
    assert lines == [row]
    assert gff_texts == []
    assert meta["interproscan_versions"] == ["5.78-109.0"]


def test_run_api_submits_in_waves_then_polls(tmp_path, monkeypatch):
    # EBI Job Dispatcher fair-use: jobs run in WAVES of at most MAX_CONCURRENT_JOBS.
    # Each wave is submitted together, then every job in it is polled/fetched before
    # the next wave is submitted -- this caps a species' concurrent jobs so several
    # species stay within EBI's limit while still overlapping jobs for speed.
    order = []
    ids = iter([f"JOB{i}" for i in range(10)])
    monkeypatch.setattr(interpro, "CHUNK_SIZE", 30)          # small chunks for the test
    monkeypatch.setattr(interpro, "MAX_CONCURRENT_JOBS", 2)  # small waves for the test
    monkeypatch.setattr(interpro.time, "sleep", lambda s: None)   # no real spacing
    monkeypatch.setattr(interpro, "_submit", lambda *a, **k: (order.append("submit"), next(ids))[1])
    monkeypatch.setattr(interpro, "_poll", lambda *a, **k: (order.append("poll"), "FINISHED")[1])
    monkeypatch.setattr(interpro, "_result_types", lambda *a, **k: ["tsv", "gff"])
    monkeypatch.setattr(interpro, "_fetch_result",
                        lambda *a, **k: "P\tmd5\t9\tPfam\tPF1\td\t1\t9\t1e-9\tT\td\tIPR\td")
    monkeypatch.setattr(interpro, "_job_version", lambda *a, **k: "5.78-109.0")

    # 65 sequences, CHUNK_SIZE 30 -> 3 chunks; MAX_CONCURRENT_JOBS 2 -> waves [2, 1]:
    # submit both of wave 1, poll both, then submit + poll wave 2.
    interpro._run_api(_api_cfg(tmp_path), _write_cand_fasta(tmp_path, 65))
    assert order == ["submit", "submit", "poll", "poll", "submit", "poll"]
