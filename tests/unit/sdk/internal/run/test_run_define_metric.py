"""Run.define_metric 输入校验分支的最小集成测试。

key 的清洗/截断/强转语义由 test_run_fmt.py 的 TestValidateKey 覆盖，
此处固化接线行为：清洗后的 key 进入事件、非法 key / 非法 glob 报错且不发出事件、
glob 分类结果（is_glob）随事件显式下发。
"""

import threading
from unittest.mock import MagicMock

import pytest

from swanlab.sdk.internal.bus.events import MetricDefineEvent
from swanlab.sdk.internal.pkg import fork
from swanlab.sdk.internal.run import Run, fmt


class _MockRun:
    """最小化的 Run 替身，供非绑定方法测试使用"""

    def __init__(self):
        self._api_lock = threading.RLock()
        self._init_pid = fork.current_pid()
        self.alive = True
        self._components = MagicMock()


@pytest.fixture(autouse=True)
def _reset_warned_keys():
    """fmt._WARNED_KEYS 是进程级缓存，用例前后重置避免跨用例污染。"""
    fmt._WARNED_KEYS.clear()
    yield
    fmt._WARNED_KEYS.clear()


def _emitted_key(mock_run: _MockRun) -> str:
    """断言事件已发出并返回其 key。"""
    mock_run._components.emitter.emit.assert_called_once()
    event = mock_run._components.emitter.emit.call_args[0][0]
    assert isinstance(event, MetricDefineEvent)
    return event.key


class TestDefineMetricKeyValidation:
    def test_sanitized_key_emitted(self, monkeypatch):
        """非法边缘字符被清洗后再进入事件，规则与 log 侧的规范 key 对齐。"""
        monkeypatch.setattr(fmt.console, "warning", MagicMock())
        mock_run = _MockRun()
        Run.define_metric(mock_run, "train/loss/")  # type: ignore
        assert _emitted_key(mock_run) == "train/loss"

    def test_invalid_key_rejected_without_emit(self, monkeypatch):
        """清洗后为空的 key 报错且不发出任何事件。"""
        error = MagicMock()
        monkeypatch.setattr(fmt.console, "error", error)
        mock_run = _MockRun()
        Run.define_metric(mock_run, "///")  # type: ignore
        error.assert_called_once()
        mock_run._components.emitter.emit.assert_not_called()


class TestDefineMetricGlobValidation:
    """glob 模式校验与 is_glob 分类（校验在 API 层同步完成，resolver 不再解析字符串）。"""

    @pytest.mark.parametrize("key", ["train/*", "*"], ids=lambda v: repr(v))
    def test_valid_glob_emitted_with_is_glob_true(self, key, monkeypatch):
        """合法 glob（末尾单个 '*'）原样进入事件，且 is_glob=True。"""
        monkeypatch.setattr(fmt.console, "error", MagicMock())
        mock_run = _MockRun()
        Run.define_metric(mock_run, key)  # type: ignore
        mock_run._components.emitter.emit.assert_called_once()
        event = mock_run._components.emitter.emit.call_args[0][0]
        assert isinstance(event, MetricDefineEvent)
        assert event.key == key
        assert event.is_glob is True

    def test_exact_key_emitted_with_is_glob_false(self, monkeypatch):
        """无 '*' 的 key 进入事件时 is_glob=False。"""
        monkeypatch.setattr(fmt.console, "error", MagicMock())
        mock_run = _MockRun()
        Run.define_metric(mock_run, "train/loss")  # type: ignore
        event = mock_run._components.emitter.emit.call_args[0][0]
        assert event.is_glob is False

    @pytest.mark.parametrize(
        "bad",
        [
            "*loss",  # '*' 在开头
            "train/*/loss",  # '*' 在中间
            "train/**",  # 末尾两个 '*'
            "**",  # 多个 '*'
            "a*b",  # '*' 在中间且非末尾
            "a***b",  # 多个 '*'
        ],
        ids=lambda v: repr(v),
    )
    def test_invalid_glob_rejected_without_emit(self, bad, monkeypatch):
        """非法 glob 报错且不发出任何事件。"""
        error = MagicMock()
        monkeypatch.setattr(fmt.console, "error", error)
        mock_run = _MockRun()
        Run.define_metric(mock_run, bad)  # type: ignore
        error.assert_called_once()
        mock_run._components.emitter.emit.assert_not_called()

    def test_invalid_glob_does_not_block_subsequent_valid_define(self, monkeypatch):
        """非法 define 被拒绝后，后续合法 define 仍正常发出事件。"""
        monkeypatch.setattr(fmt.console, "error", MagicMock())
        mock_run = _MockRun()
        Run.define_metric(mock_run, "*loss")  # type: ignore
        Run.define_metric(mock_run, "train/*")  # type: ignore
        assert mock_run._components.emitter.emit.call_count == 1
        event = mock_run._components.emitter.emit.call_args[0][0]
        assert event.key == "train/*"
        assert event.is_glob is True
