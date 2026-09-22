from types import SimpleNamespace
from unittest.mock import MagicMock

import pytest
from pyroaring import BitMap64

import swanlab
from swanlab.sdk.cmd import temporary
from swanlab.sdk.cmd.init import init
from swanlab.sdk.cmd.run import finish


class _Metric:
    def __init__(self, min_step: int, steps: BitMap64):
        self.min_step = min_step
        self.steps = steps


def test_allow_backfill_only_relaxes_the_requested_metric_and_step():
    loss = _Metric(min_step=10, steps=BitMap64([3, 4]))
    accuracy = _Metric(min_step=10, steps=BitMap64([3]))
    metrics = SimpleNamespace(get={"train/loss": loss, "accuracy": accuracy}.get)
    run = SimpleNamespace(_ctx=SimpleNamespace(core=SimpleNamespace(_metrics=metrics)))

    temporary._allow_backfill(run, {"train": {"loss": 0.1}}, step=3)

    assert loss.min_step == 2
    assert 3 not in loss.steps
    assert 4 in loss.steps
    assert accuracy.min_step == 10
    assert 3 in accuracy.steps


def test_log_backfill_skips_existing_metrics_and_uploads_only_missing(monkeypatch):
    run = init(mode="disabled", resume=True, id="run-id")
    run._path = "/alice/project/run-id"
    run.log = MagicMock()
    api = MagicMock()
    api.run.return_value.metrics.return_value = {
        "list": [
            {"key": "loss", "metrics": [{"step": 3, "value": 0.1}]},
            {"key": "accuracy", "metrics": []},
        ]
    }
    monkeypatch.setattr(temporary, "Api", lambda: api)
    warning = MagicMock()
    monkeypatch.setattr(temporary.console, "warning", warning)
    try:
        with pytest.warns(FutureWarning, match="temporary compatibility API"):
            swanlab.log_backfill({"loss": 0.1, "accuracy": 0.9}, step=3)
        api.run.assert_called_once_with("alice/project/run-id")
        api.run.return_value.metrics.assert_called_once_with(
            keys=["loss", "accuracy"], range_query={"start": 3, "end": 3}
        )
        warning.assert_called_once_with("Skip backfill: metric 'loss' at step=3 already exists on SwanLab.")
        run.log.assert_called_once_with({"accuracy": 0.9}, step=3)
    finally:
        finish()
