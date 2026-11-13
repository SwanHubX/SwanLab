"""
@author: Nexisato
@file: test_reap.py
@time: 2026/10/27
@description: 分级进程回收测试。
"""

import subprocess
import sys

from swanlab.sdk.internal.core_go import process as process_module
from swanlab.sdk.internal.core_go.process import reap_process


def test_reap_none_returns_none():
    assert reap_process(None) is None


def test_reap_returns_normal_exit_code():
    proc = subprocess.Popen([sys.executable, "-c", "raise SystemExit(0)"])

    assert reap_process(proc) == 0
    assert proc.poll() == 0


def test_reap_is_idempotent():
    proc = subprocess.Popen([sys.executable, "-c", "raise SystemExit(7)"])

    assert reap_process(proc) == 7
    assert reap_process(proc) == 7


def test_reap_escalates_to_kill_when_sigterm_is_ignored(monkeypatch):
    monkeypatch.setattr(process_module, "WAIT_TIMEOUT", 0.2)
    monkeypatch.setattr(process_module, "TERMINATE_TIMEOUT", 0.2)
    monkeypatch.setattr(process_module, "KILL_TIMEOUT", 2.0)
    proc = subprocess.Popen(
        [
            sys.executable,
            "-c",
            "import signal, time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(30)",
        ]
    )

    assert reap_process(proc) is not None
    assert proc.poll() is not None
