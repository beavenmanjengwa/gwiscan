"""Tests for `run`'s resume/stop/skip stage selection (pipeline.plan_stages)."""

import pytest

from gwiscan.pipeline import plan_stages

STAGES = [
    ("a", "Stage A", None, "a"),
    ("b", "Stage B", None, "b"),
    ("c", "Stage C", None, "c"),
    ("d", "Stage D", None, "d"),
]


def _keys(selected):
    return [k for k, *_ in selected]


def test_no_directives_runs_everything():
    assert _keys(plan_stages(STAGES)) == ["a", "b", "c", "d"]


def test_from_stage_skips_earlier_stages():
    assert _keys(plan_stages(STAGES, from_stage="c")) == ["c", "d"]


def test_from_stage_at_first_stage_is_a_noop():
    assert _keys(plan_stages(STAGES, from_stage="a")) == ["a", "b", "c", "d"]


def test_until_stops_after_named_stage_inclusive():
    assert _keys(plan_stages(STAGES, until_stage="b")) == ["a", "b"]


def test_from_stage_and_until_combine():
    assert _keys(plan_stages(STAGES, from_stage="b", until_stage="c")) == ["b", "c"]


def test_skip_stages_still_present_caller_decides_meaning():
    # plan_stages validates skip keys but doesn't drop them from the slice --
    # the caller (pipeline.run) is the one that turns membership into a [SKIP] log
    # line instead of executing the stage.
    selected = plan_stages(STAGES, skip={"b"})
    assert _keys(selected) == ["a", "b", "c", "d"]


def test_unknown_from_stage_raises_with_valid_list():
    with pytest.raises(ValueError, match=r"--from-stage 'nope'.*Valid stages: a, b, c, d"):
        plan_stages(STAGES, from_stage="nope")


def test_unknown_until_raises_with_valid_list():
    with pytest.raises(ValueError, match=r"--until 'nope'.*Valid stages: a, b, c, d"):
        plan_stages(STAGES, until_stage="nope")


def test_unknown_skip_raises_with_valid_list():
    with pytest.raises(ValueError, match=r"--skip has unknown stage\(s\) nope.*Valid stages"):
        plan_stages(STAGES, skip={"nope"})


def test_multiple_unknown_skips_all_named_sorted():
    with pytest.raises(ValueError, match=r"unknown stage\(s\) nope1, nope2"):
        plan_stages(STAGES, skip={"nope2", "nope1"})


def test_until_before_from_stage_raises():
    with pytest.raises(ValueError, match=r"occurs at or before --from-stage"):
        plan_stages(STAGES, from_stage="c", until_stage="a")


def test_until_equal_to_from_stage_is_valid_single_stage():
    assert _keys(plan_stages(STAGES, from_stage="b", until_stage="b")) == ["b"]
