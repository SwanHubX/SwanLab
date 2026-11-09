"""init 后置步骤失败时的回滚。

覆盖 cmd/init.py 在 Run 构造完成后的失败路径：必须撤销 hook、清理全局 Run，
并让后续 init 仍可正常执行。
"""

import sys

import pytest

import swanlab
from swanlab.sdk.cmd.init import init
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
