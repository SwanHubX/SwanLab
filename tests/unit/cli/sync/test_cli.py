"""
@author: nexisato
@file: test_cli.py
@time: 2026/10/23
@description: swanlab sync CLI 范围参数测试
"""

import importlib

from click.testing import CliRunner

from swanlab.cli.sync import sync as sync_cli

sync_module = importlib.import_module("swanlab.cli.sync")


def test_sync_cli_forwards_range_options(tmp_path, monkeypatch):
    calls = []

    def fake_sync(path, settings, **kwargs):
        calls.append((path, settings, kwargs))

    monkeypatch.setattr(sync_module, "sync_cmd", fake_sync)

    result = CliRunner().invoke(
        sync_cli,
        [
            str(tmp_path),
            "--id",
            "run-1",
            "--step-start",
            "10",
            "--step-end",
            "20",
            "--start-time",
            "2026-10-20 10:00:00",
            "--end-time",
            "2026-10-20 12:00:00",
        ],
    )

    assert result.exit_code == 0, result.output
    assert len(calls) == 1
    path, settings, kwargs = calls[0]
    assert path == tmp_path
    assert settings.run.id == "run-1"
    assert kwargs == {
        "step_start": 10,
        "step_end": 20,
        "start_time": "2026-10-20 10:00:00",
        "end_time": "2026-10-20 12:00:00",
    }


def test_sync_cli_without_range_passes_none(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(sync_module, "sync_cmd", lambda path, settings, **kwargs: calls.append((path, kwargs)))

    result = CliRunner().invoke(sync_cli, [str(tmp_path)])

    assert result.exit_code == 0, result.output
    assert calls[0][1] == {
        "step_start": None,
        "step_end": None,
        "start_time": None,
        "end_time": None,
    }


def test_sync_cli_rejects_negative_step(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(sync_module, "sync_cmd", lambda *args, **kwargs: calls.append(args))

    result = CliRunner().invoke(sync_cli, [str(tmp_path), "--step-start", "-1"])

    assert result.exit_code != 0
    assert calls == []
