"""Tests for the run-only resume/stop/skip CLI flags and --list-stages."""

from gwiscan import cli, pipeline


def test_run_only_flags_not_on_other_subcommands():
    # --from-stage etc. are meaningless for a single-stage subcommand like
    # `gwiscan weblogo`; they must not be offered there.
    parser = cli.build_parser()
    args = parser.parse_args(["weblogo"])
    assert not hasattr(args, "FROM_STAGE")
    assert not hasattr(args, "UNTIL_STAGE")
    assert not hasattr(args, "SKIP_STAGES")


def test_run_parses_resume_stop_skip_flags():
    parser = cli.build_parser()
    args = parser.parse_args([
        "run", "--from-stage", "merge", "--until", "compile", "--skip", "meme,weblogo",
    ])
    assert args.FROM_STAGE == "merge"
    assert args.UNTIL_STAGE == "compile"
    assert args.SKIP_STAGES == ["meme", "weblogo"]


def test_config_from_args_tolerates_subcommands_without_resume_flags(tmp_path):
    # _config_from_args reads _OVERRIDE_KEYS via getattr(args, key, None); a
    # subcommand without --from-stage/--until/--skip must not raise AttributeError.
    parser = cli.build_parser()
    args = parser.parse_args(["preflight", "-C", str(tmp_path)])
    cfg = cli._config_from_args(args)
    assert cfg.FROM_STAGE == ""
    assert cfg.UNTIL_STAGE == ""
    assert cfg.SKIP_STAGES == []


def test_list_stages_prints_stage_keys_and_exits_zero(capsys):
    rc = cli.main(["run", "--list-stages"])
    assert rc == 0
    out = capsys.readouterr().out.splitlines()
    assert out == pipeline.STAGE_KEYS


def test_list_stages_does_not_require_a_project_dir(tmp_path, capsys, monkeypatch):
    # --list-stages must work from anywhere -- no config.yaml/project needed.
    monkeypatch.chdir(tmp_path)
    rc = cli.main(["run", "--list-stages"])
    assert rc == 0


def test_trim_subcommand_and_trimal_flags(tmp_path):
    parser = cli.build_parser()
    args = parser.parse_args(["trim", "--trimal-bin", "/opt/trimal",
                              "--trimal-method", "gappyout"])
    assert args.command == "trim"
    assert args.TRIMAL_BIN == "/opt/trimal"
    assert args.TRIMAL_METHOD == "gappyout"
    cfg = cli._config_from_args(args)
    assert cfg.TRIMAL_BIN == "/opt/trimal"
    assert cfg.TRIMAL_METHOD == "gappyout"


def test_trim_stage_in_list_between_msa_and_iqtree():
    keys = pipeline.STAGE_KEYS
    assert "trim" in keys
    assert keys.index("msa") < keys.index("trim") < keys.index("iqtree")


def test_threads_short_flag_p(tmp_path):
    parser = cli.build_parser()
    for argv in (["run", "-p", "6"], ["run", "--threads", "6"]):
        args = parser.parse_args(argv)
        assert args.THREADS == 6


def test_only_species_and_retry_failed_are_run_only(tmp_path):
    parser = cli.build_parser()
    args = parser.parse_args([
        "run", "--only-species", "Ath,Gma", "--retry-failed",
    ])
    assert args.SPECIES_ONLY == ["Ath", "Gma"]
    assert args.RETRY_FAILED is True
    # not offered on single-stage subcommands
    other = parser.parse_args(["weblogo"])
    assert not hasattr(other, "SPECIES_ONLY")
    assert not hasattr(other, "RETRY_FAILED")


def test_retry_failed_absent_is_none_so_config_default_wins(tmp_path):
    # store_true with default=None -> absent flag doesn't override config's False.
    parser = cli.build_parser()
    args = parser.parse_args(["run", "-C", str(tmp_path)])
    cfg = cli._config_from_args(args)
    assert cfg.RETRY_FAILED is False
    assert cfg.SPECIES_ONLY == []


def test_weblogo_and_meme_bin_flags_flow_into_config(tmp_path):
    parser = cli.build_parser()
    args = parser.parse_args([
        "weblogo", "-C", str(tmp_path),
        "--weblogo-bin", "/envs/x/bin/weblogo",
    ])
    cfg = cli._config_from_args(args)
    assert cfg.WEBLOGO_BIN == "/envs/x/bin/weblogo"

    args2 = parser.parse_args(["meme", "-C", str(tmp_path), "--meme-bin", "/opt/meme"])
    cfg2 = cli._config_from_args(args2)
    assert cfg2.MEME_BIN == "/opt/meme"


def test_annotation_flag_flows_into_config(tmp_path):
    parser = cli.build_parser()
    args = parser.parse_args(["coords", "-C", str(tmp_path), "--annotation", "/data/ann.gff3"])
    cfg = cli._config_from_args(args)
    assert cfg.ANNOTATION == "/data/ann.gff3"
    assert str(cfg.annotation) == "/data/ann.gff3"
