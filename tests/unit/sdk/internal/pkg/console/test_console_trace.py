"""
@author: cunyue
@file: test_console_trace.py
@time: 2026/4/11
@description: 测试 console.trace 异常栈打印
"""

from unittest.mock import MagicMock, patch

from swanlab.sdk.internal.pkg import console

# ==============================================================================
# 基础行为
# ==============================================================================


def test_trace_no_active_exception():
    """不在 except 块内调用 trace 时，应无操作（返回 None）"""
    mock_error = MagicMock()
    with patch.object(console, "error", mock_error):
        result = console.trace("should be ignored")
    assert result is None
    mock_error.assert_not_called()


def test_trace_inside_except():
    """在 except 块内调用 trace 时，应委托给 error 打印异常栈信息"""
    mock_error = MagicMock()
    with patch.object(console, "error", mock_error):
        try:
            _ = 1 / 0
        except ZeroDivisionError:
            console.trace("division failed")

    mock_error.assert_called_once()
    call_kwargs = mock_error.call_args
    # console_args 应包含异常类型和消息
    console_msg = call_kwargs[1]["console_args"][0]
    assert "ZeroDivisionError" in console_msg
    assert "division by zero" in console_msg
    # file_args 应包含完整错误栈
    file_msg = call_kwargs[1]["file_args"][0]
    assert "ZeroDivisionError" in file_msg
    assert "division by zero" in file_msg


def test_trace_includes_prefix_message():
    """trace 的前缀消息应出现在输出中"""
    mock_error = MagicMock()
    with patch.object(console, "error", mock_error):
        try:
            raise ValueError("bad value")
        except ValueError:
            console.trace("my prefix")

    console_msg = mock_error.call_args[1]["console_args"][0]
    assert console_msg.startswith("my prefix:")
    assert "ValueError" in console_msg
    assert "bad value" in console_msg


def test_trace_no_prefix():
    """不传前缀消息时，输出应直接以异常栈开头"""
    mock_error = MagicMock()
    with patch.object(console, "error", mock_error):
        try:
            raise RuntimeError("oops")
        except RuntimeError:
            console.trace()

    console_msg = mock_error.call_args[1]["console_args"][0]
    assert "RuntimeError" in console_msg
    assert "oops" in console_msg


# ==============================================================================
# max_frames 参数
# ==============================================================================


def _deep_error(n: int):
    """递归生成 n 层调用栈后抛出异常"""
    if n <= 0:
        raise ValueError("bottom")
    _deep_error(n - 1)


def test_trace_max_frames_default():
    """默认 max_frames=2，终端只显示最后 2 层栈帧"""
    mock_error = MagicMock()
    with patch.object(console, "error", mock_error):
        try:
            _deep_error(5)
        except ValueError:
            console.trace("frame test")

    console_msg = mock_error.call_args[1]["console_args"][0]
    # 短栈末尾应为: ... -> _deep_error:88 -> _deep_error:87 -> ValueError: bottom
    # 即 max_frames 个栈帧 + 1 个异常类型
    assert "_deep_error" in console_msg
    assert "ValueError" in console_msg
    # 栈帧间用 -> 连接，最后也用 -> 连接异常类型
    parts = console_msg.split(" -> ")
    # 去掉前缀部分后，栈帧 + 异常 = max_frames + 1
    # 前缀 "frame test: ..." 会被包含在第一个 part 中
    trace_parts = [p for p in parts if ":" in p and not p.startswith("frame test")]
    assert len(trace_parts) <= 3  # max_frames=2 个栈帧 + 1 个异常


def test_trace_max_frames_custom():
    """自定义 max_frames 时，终端只显示对应层数的栈帧"""
    mock_error = MagicMock()
    with patch.object(console, "error", mock_error):
        try:
            _deep_error(5)
        except ValueError:
            console.trace("frame test", max_frames=4)

    console_msg = mock_error.call_args[1]["console_args"][0]
    # max_frames=4 时，4 个栈帧 + 1 个异常
    assert "_deep_error" in console_msg
    assert "ValueError" in console_msg


# ==============================================================================
# level_name 参数
# ==============================================================================


def test_trace_level_error():
    """level_name="error" 时应委托给 error"""
    mock_error = MagicMock()
    mock_debug = MagicMock()
    with patch.object(console, "error", mock_error), patch.object(console, "debug", mock_debug):
        try:
            raise OSError("fail")
        except OSError:
            console.trace("err", level_name="error")

    mock_error.assert_called_once()
    mock_debug.assert_not_called()


def test_trace_level_debug():
    """level_name="debug" 时应委托给 debug"""
    mock_error = MagicMock()
    mock_debug = MagicMock()
    with patch.object(console, "error", mock_error), patch.object(console, "debug", mock_debug):
        try:
            raise OSError("fail")
        except OSError:
            console.trace("dbg", level_name="debug")

    mock_debug.assert_called_once()
    mock_error.assert_not_called()


# ==============================================================================
# write_to_file 参数
# ==============================================================================


def test_trace_write_to_file_default():
    """默认 write_to_file=True，应传递给 error/debug"""
    mock_error = MagicMock()
    with patch.object(console, "error", mock_error):
        try:
            raise ValueError("x")
        except ValueError:
            console.trace("msg")

    assert mock_error.call_args[1]["write_to_file"] is True


def test_trace_write_to_file_false():
    """write_to_file=False 时，应传递给 error/debug"""
    mock_error = MagicMock()
    with patch.object(console, "error", mock_error):
        try:
            raise ValueError("x")
        except ValueError:
            console.trace("msg", write_to_file=False)

    assert mock_error.call_args[1]["write_to_file"] is False


# ==============================================================================
# file_args 包含完整栈
# ==============================================================================


def test_trace_file_args_has_full_traceback():
    """file_args 应包含完整的多层错误栈，不受 max_frames 限制"""
    mock_error = MagicMock()
    with patch.object(console, "error", mock_error):
        try:
            _deep_error(5)
        except ValueError:
            console.trace("full stack", max_frames=1)

    file_msg = mock_error.call_args[1]["file_args"][0]
    # 完整栈应包含所有 _deep_error 调用（含 repeated 行中的引用）
    assert file_msg.count("_deep_error") >= 5
