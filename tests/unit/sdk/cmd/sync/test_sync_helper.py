"""
@author: cunyue
@file: test_sync_helper.py.py
@time: 2026/5/17 19:30
@description: 测试sync的工具函数
"""

import pytest
from pydantic import ValidationError

from swanlab.sdk.cmd.sync import ensure_run_dir


def test_ensure_run_dir_ok(tmp_path):
    assert ensure_run_dir(tmp_path) == tmp_path.resolve()


def test_ensure_run_dir_requires_directory(tmp_path):
    file_path = tmp_path / "run.txt"
    file_path.write_text("", encoding="utf-8")

    with pytest.raises(ValidationError):
        ensure_run_dir(file_path)


def test_ensure_run_dir_requires_readable(tmp_path, monkeypatch):
    monkeypatch.setattr("swanlab.sdk.cmd.sync.os.access", lambda *_: False)

    with pytest.raises(PermissionError):
        ensure_run_dir(tmp_path)
