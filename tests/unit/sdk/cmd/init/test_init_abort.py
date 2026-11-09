"""init 后置步骤失败时的回滚。

覆盖 cmd/init.py 在 Run 构造完成后的失败路径：必须撤销 hook、清理全局 Run，
并让后续 init 仍可正常执行。
"""

import sys

import pytest

import swanlab
from swanlab.sdk.cmd.init import init
from swanlab.sdk.internal.core_python import CorePython
from swanlab.sdk.internal.pkg import helper
from swanlab.sdk.internal.run import has_run


def test_config_load_failure_cleans_run_and_hooks():
    original_excepthook = sys.excepthook

    with pytest.raises(FileNotFoundError):
        init(mode="disabled", config="/definitely/missing/swanlab-config.yaml")

    assert not has_run()
    assert swanlab.run is None
    assert sys.excepthook is original_excepthook


def test_run_constructor_failure_cleans_run_and_hooks(monkeypatch):
    original_excepthook = sys.excepthook

    def boom(*_args, **_kwargs):
        raise RuntimeError("welcome failed")

    monkeypatch.setattr("swanlab.sdk.internal.run.greeting.welcome", boom)

    with pytest.raises(RuntimeError, match="welcome failed"):
        init(mode="disabled")

    assert not has_run()
    assert swanlab.run is None
    assert sys.excepthook is original_excepthook


def test_init_can_retry_after_failure():
    with pytest.raises(FileNotFoundError):
        init(mode="disabled", config="/definitely/missing/swanlab-config.yaml")

    run = init(mode="disabled")
    assert has_run()
    run.finish()


def test_local_config_failure_reclaims_python_core(tmp_path, monkeypatch):
    """local init 已完成 deliver_run_start 后配置加载失败，CorePython 必须回收本地资源。"""
    monkeypatch.setattr(helper, "get_core_impl", lambda: "python")
    observed = []
    original_abort = CorePython.abort

    def record_abort(self, cause=None):
        original_abort(self, cause)
        observed.append((self._store, self._watcher, self._transport, self._heartbeat))

    monkeypatch.setattr(CorePython, "abort", record_abort)

    with pytest.raises(FileNotFoundError):
        init(
            mode="local",
            log_dir=str(tmp_path / "swanlog"),
            config=str(tmp_path / "missing-config.yaml"),
        )

    assert observed == [(None, None, None, None)]
    assert not has_run()
    assert swanlab.run is None
