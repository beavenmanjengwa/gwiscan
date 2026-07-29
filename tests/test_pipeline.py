"""Stage-ordering guard for the linear pipeline.

Every stage's required input must be produced by an earlier stage. The key
regression this locks: domain-bed runs before extract-domains (domains.run
hard-requires results/domains.bed). This is what a per-stage unit test cannot
catch, and it is exactly what was missing before.
"""

import pytest

import gwiscan.pipeline as pipeline
from gwiscan.config import Config


def _patch_all(monkeypatch, calls, tools_available=True):
    # Ordering/resume tests are about which stages run in what order, not about
    # which external binaries happen to be installed on the test host -- force the
    # availability check so auto-skip doesn't drop optional-tool stages here.
    # (Auto-skip has its own dedicated tests that flip this off deliberately.)
    monkeypatch.setattr(pipeline.external, "available", lambda tool: tools_available)

    def recorder(name):
        def _fn(cfg, *args, **kwargs):
            calls.append(name)
        return _fn

    monkeypatch.setattr(pipeline.preflight, "run", recorder("preflight"))
    monkeypatch.setattr(pipeline.setupdb, "setup_shared", recorder("setup_shared"))
    monkeypatch.setattr(pipeline.setupdb, "setup_proteome", recorder("setup_proteome"))
    monkeypatch.setattr(pipeline.hmm, "run", recorder("hmm"))
    monkeypatch.setattr(pipeline.diamond, "run", recorder("diamond"))
    monkeypatch.setattr(pipeline.candidates, "run", recorder("merge"))
    monkeypatch.setattr(pipeline.score, "run", recorder("score"))
    monkeypatch.setattr(pipeline.interpro, "run", recorder("interpro"))
    monkeypatch.setattr(pipeline.confirm, "run", recorder("confirm"))
    monkeypatch.setattr(pipeline.protparam, "run", recorder("protparam"))
    monkeypatch.setattr(pipeline.targetp, "run", recorder("targetp"))
    monkeypatch.setattr(pipeline.deeptmhmm, "run", recorder("deeptmhmm"))
    monkeypatch.setattr(pipeline.deeploc, "run", recorder("deeploc"))
    monkeypatch.setattr(pipeline.compile_stage, "run", recorder("compile"))
    monkeypatch.setattr(pipeline.domain_bed, "run", recorder("domain_bed"))
    monkeypatch.setattr(pipeline.domains, "run", recorder("extract_domains"))
    monkeypatch.setattr(pipeline.mature, "run", recorder("extract_mature"))
    monkeypatch.setattr(pipeline.msa, "run", recorder("msa"))
    monkeypatch.setattr(pipeline.trimal, "run", recorder("trim"))
    monkeypatch.setattr(pipeline.logos, "run", recorder("weblogo"))
    monkeypatch.setattr(pipeline.meme, "run", recorder("meme"))
    monkeypatch.setattr(pipeline.iqtree, "run", recorder("iqtree"))
    monkeypatch.setattr(pipeline.provenance, "run", recorder("provenance"))


def test_domain_bed_precedes_extract_domains(tmp_path, monkeypatch):
    calls = []
    _patch_all(monkeypatch, calls)
    pipeline.run(Config(root=tmp_path))

    # domains.run requires results/domains.bed, which only domain_bed.run produces.
    assert "domain_bed" in calls, "domain-bed stage is missing from the pipeline"
    assert calls.index("domain_bed") < calls.index("extract_domains")


def test_stage_prerequisite_ordering(tmp_path, monkeypatch):
    calls = []
    _patch_all(monkeypatch, calls)
    pipeline.run(Config(root=tmp_path))

    order = {name: i for i, name in enumerate(calls)}
    # setup -> searches -> merge -> interpro -> confirm
    assert order["setup_proteome"] < order["hmm"] < order["merge"]
    assert order["setup_proteome"] < order["diamond"] < order["merge"]
    assert order["merge"] < order["interpro"] < order["confirm"]
    # detectability compares both searches, so both must have run
    assert order["hmm"] < order["score"] and order["diamond"] < order["score"]
    # extract-domains feeds msa/weblogo/meme/iqtree
    for downstream in ("msa", "weblogo", "meme", "iqtree"):
        assert order["extract_domains"] < order[downstream]
    # the mature track: targetp -> extract-mature -> msa (mature tree)
    assert order["targetp"] < order["extract_mature"] < order["msa"]
    # trimAl runs after the MSA and before IQ-TREE (it feeds the tree)
    assert order["msa"] < order["trim"] < order["iqtree"]


def test_from_stage_resumes_without_rerunning_earlier_stages(tmp_path, monkeypatch):
    # The exact scenario this was built for: a run died partway (e.g. an optional
    # tool missing) and everything through search-diamond already succeeded --
    # resuming from merge must not touch preflight/setup/hmm/diamond again.
    calls = []
    _patch_all(monkeypatch, calls)
    cfg = Config(root=tmp_path, FROM_STAGE="merge")
    pipeline.run(cfg)

    for earlier in ("preflight", "setup_shared", "setup_proteome", "hmm", "diamond"):
        assert earlier not in calls
    assert calls[0] == "merge"
    assert "provenance" in calls


def test_until_stops_after_named_stage(tmp_path, monkeypatch):
    calls = []
    _patch_all(monkeypatch, calls)
    cfg = Config(root=tmp_path, UNTIL_STAGE="compile")
    pipeline.run(cfg)

    assert calls[-1] == "compile"
    for later in ("domain_bed", "extract_domains", "msa", "weblogo", "meme", "iqtree", "provenance"):
        assert later not in calls


def test_skip_stages_are_not_invoked_but_run_continues(tmp_path, monkeypatch):
    # This is the failure mode that motivated the feature: MEME missing used to
    # abort the whole run. --skip meme must let iqtree/provenance still run.
    calls = []
    _patch_all(monkeypatch, calls)
    cfg = Config(root=tmp_path, SKIP_STAGES=["meme"])
    pipeline.run(cfg)

    assert "meme" not in calls
    assert "iqtree" in calls
    assert "provenance" in calls


def test_unknown_from_stage_raises_runtime_error(tmp_path, monkeypatch):
    calls = []
    _patch_all(monkeypatch, calls)
    cfg = Config(root=tmp_path, FROM_STAGE="not-a-stage")
    with pytest.raises(RuntimeError, match="not-a-stage"):
        pipeline.run(cfg)
    assert calls == []  # fails fast, before touching any stage


def test_missing_optional_tool_is_auto_skipped_not_fatal(tmp_path, monkeypatch):
    # The core of the fix: a missing optional tool (weblogo/meme/iqtree) must not
    # abort the run -- the stage is skipped, everything after it still runs.
    calls = []
    _patch_all(monkeypatch, calls)
    # meme's binary is unavailable; weblogo and iqtree are present.
    monkeypatch.setattr(pipeline.external, "available",
                        lambda tool: tool != "meme")
    pipeline.run(Config(root=tmp_path))  # MEME_BIN defaults to "meme"

    assert "meme" not in calls
    assert "weblogo" in calls
    assert "iqtree" in calls
    assert "provenance" in calls


def test_all_optional_tools_missing_still_completes(tmp_path, monkeypatch):
    calls = []
    _patch_all(monkeypatch, calls)
    # Every optional-tool binary is missing; a required-stage function is not.
    optional_bins = {"trimal", "weblogo", "meme", "iqtree"}
    monkeypatch.setattr(pipeline.external, "available",
                        lambda tool: tool not in optional_bins)
    pipeline.run(Config(root=tmp_path))

    for skipped in ("trim", "weblogo", "meme", "iqtree"):
        assert skipped not in calls
    # the core identification + annotation path still ran end to end
    for ran in ("hmm", "diamond", "merge", "interpro", "compile", "provenance"):
        assert ran in calls


def test_explicit_skip_message_wins_over_auto_skip(tmp_path, monkeypatch, capsys):
    calls = []
    _patch_all(monkeypatch, calls)
    # meme is BOTH explicitly --skip'd and unavailable: it should be reported as an
    # explicit [SKIP], not an [AUTO-SKIP] (no duplicate/contradictory message).
    monkeypatch.setattr(pipeline.external, "available",
                        lambda tool: tool != "meme")
    pipeline.run(Config(root=tmp_path, SKIP_STAGES=["meme"]))

    out = capsys.readouterr().out
    meme_label = next(lbl for key, lbl, *_ in pipeline.STAGES if key == "meme")
    assert f"[SKIP] {meme_label} (--skip meme)" in out
    assert f"[AUTO-SKIP] {meme_label}" not in out


def test_custom_bin_config_prevents_auto_skip(tmp_path, monkeypatch):
    # If the user points WEBLOGO_BIN at a real path, availability is checked against
    # THAT, not the bare name -- so the stage is not auto-skipped.
    calls = []
    _patch_all(monkeypatch, calls)
    real = "/opt/envs/x/bin/weblogo"
    monkeypatch.setattr(pipeline.external, "available", lambda tool: tool == real)
    pipeline.run(Config(root=tmp_path, WEBLOGO_BIN=real, MEME_BIN="meme",
                        IQTREE_BIN="iqtree"))

    assert "weblogo" in calls          # available at the configured path -> runs
    assert "meme" not in calls         # bare "meme" unavailable -> auto-skipped
    assert "iqtree" not in calls
