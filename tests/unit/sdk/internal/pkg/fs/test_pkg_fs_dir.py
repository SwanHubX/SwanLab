"""
@author: cunyue
@file: test_fs_dir.py
@time: 2026/3/11 14:30
@description: SwanLab SDK 文件系统辅助函数测试
"""

import errno
import importlib
import sys
from pathlib import Path

import pytest

from swanlab.sdk.internal.pkg.fs import dir

sys.modules[".."] = type(sys)("mock_parent")
sys.modules[".."].console = type(sys)("mock_console")


def test_safe_mkdir_normal(tmp_path: Path):
    """测试常规单级目录创建"""
    target = tmp_path / "test_dir"
    result = dir.safe_mkdir(target)

    assert result == target
    assert target.exists()
    assert target.is_dir()


def test_safe_mkdir_nested(tmp_path: Path):
    """测试多级目录创建"""
    target = tmp_path / "a" / "b" / "c"
    dir.safe_mkdir(target)

    assert target.exists()
    assert target.is_dir()


def test_safe_mkdir_already_exists(tmp_path: Path):
    """测试目录已存在时的情况"""
    target = tmp_path / "exist_dir"
    target.mkdir()

    # 应该直接通过，不会报错
    result = dir.safe_mkdir(target)
    assert result == target


def test_safe_mkdir_nas_timeout(monkeypatch, tmp_path: Path):
    """测试极端 NAS 延迟/权限问题下，探针是否能正确抛出 TimeoutError"""
    target = tmp_path / "timeout_dir"

    time_calls = [0, 10, 20, 30]
    monkeypatch.setattr("swanlab.sdk.internal.pkg.fs.dir.time.time", lambda: time_calls.pop(0))

    def mock_mkstemp(*args, **kwargs):
        raise OSError("Simulated NAS IO Error")

    monkeypatch.setattr("swanlab.sdk.internal.pkg.fs.dir.tempfile.mkstemp", mock_mkstemp)

    with pytest.raises(TimeoutError, match="is not writable within") as exc_info:
        dir.safe_mkdir(target, timeout=5.0)

    # 超时错误应链上最后一次原始 OSError，便于诊断根因
    assert isinstance(exc_info.value.__cause__, OSError)


def test_safe_mkdir_permission_denied_raises_permission_error(monkeypatch, tmp_path: Path):
    """权限不足时应立即抛 PermissionError，而不是伪装成 NAS 超时。"""
    target = tmp_path / "readonly_dir"

    def mock_mkstemp(*args, **kwargs):
        raise PermissionError(errno.EACCES, "Permission denied")

    monkeypatch.setattr("swanlab.sdk.internal.pkg.fs.dir.tempfile.mkstemp", mock_mkstemp)

    with pytest.raises(PermissionError, match="Directory .* is not writable"):
        dir.safe_mkdir(target, timeout=5.0)


def test_safe_mkdir_does_not_depend_on_temporary_file(monkeypatch, tmp_path: Path):
    """依赖约束：可写性探测不得回退使用 tempfile.TemporaryFile。

    其匿名文件 / 「打开时 unlink」语义在部分 FUSE / NAS 上会返回 EIO（bpo-22326）；
    本测试只固化「不依赖它」这一约束，真实兼容性需在原 NAS 环境验证。
    """
    target = tmp_path / "fuse_dir"
    target.mkdir()

    def mock_temporary_file(*args, **kwargs):
        raise OSError(errno.EIO, "Input/output error")

    monkeypatch.setattr("swanlab.sdk.internal.pkg.fs.dir.tempfile.TemporaryFile", mock_temporary_file)

    result = dir.safe_mkdir(target, timeout=1.0)

    assert result == target
    # 探针不应留下任何垃圾文件
    assert not list(target.glob(dir.PROBE_PREFIX + "*"))


def test_probe_writable_no_leftover(tmp_path: Path):
    """探针正常路径：写入成功后应清理本次探针文件"""
    dir._probe_writable(tmp_path)

    assert not list(tmp_path.glob(dir.PROBE_PREFIX + "*"))


def test_probe_unlink_failure_does_not_fail_writability_probe(monkeypatch, tmp_path: Path):
    """unlink 失败只保留本次探针文件，不影响目录可写性判定"""
    target = tmp_path / "unlink_fail_dir"
    target.mkdir()

    real_unlink = dir.os.unlink
    unlink_calls = []
    trace_messages = []

    def mock_unlink(path, *args, **kwargs):
        if Path(path).name.startswith(dir.PROBE_PREFIX):
            unlink_calls.append(path)
            raise OSError(errno.EIO, "Input/output error")
        return real_unlink(path, *args, **kwargs)

    def mock_trace(message, *args, **kwargs):
        trace_messages.append((message, kwargs))

    monkeypatch.setattr("swanlab.sdk.internal.pkg.fs.dir.os.unlink", mock_unlink)
    monkeypatch.setattr("swanlab.sdk.internal.pkg.safe.console.trace", mock_trace)

    assert dir.safe_mkdir(target, timeout=5.0) == target
    assert len(unlink_calls) == 1
    assert len(list(target.glob(dir.PROBE_PREFIX + "*"))) == 1
    assert len(trace_messages) == 1
    assert "Failed to clean up writability probe file" in trace_messages[0][0]
    assert trace_messages[0][1]["write_to_tty"] is False


def test_probe_cleanup_failure_does_not_mask_original_error(monkeypatch, tmp_path: Path):
    """写入失败且清理（close/unlink）也失败时，原始写入异常不得被覆盖"""
    real_unlink = dir.os.unlink
    unlink_calls = []

    def mock_write(fd, data):
        raise OSError(errno.EIO, "simulated write failure")

    def mock_unlink(path, *args, **kwargs):
        if Path(path).name.startswith(dir.PROBE_PREFIX):
            unlink_calls.append(str(path))
            raise OSError(errno.EPERM, "simulated unlink failure")
        return real_unlink(path, *args, **kwargs)

    monkeypatch.setattr("swanlab.sdk.internal.pkg.fs.dir.os.write", mock_write)
    monkeypatch.setattr("swanlab.sdk.internal.pkg.fs.dir.os.unlink", mock_unlink)

    with pytest.raises(OSError, match="simulated write failure") as exc_info:
        dir._probe_writable(tmp_path)

    # 抛出的是原始写入错误而非清理错误，且清理确实被尝试
    assert exc_info.value.errno == errno.EIO
    assert len(unlink_calls) == 1


def test_timeout_env_invalid_string(monkeypatch):
    """测试环境变量为非法字符串时，是否能安全回退到默认值 5.0"""
    monkeypatch.setenv("SWANLAB_FS_TIMEOUT", "invalid_abc123")

    mock_warnings = []
    monkeypatch.setattr("swanlab.sdk.internal.pkg.fs.dir.console.warning", lambda msg: mock_warnings.append(msg))

    importlib.reload(dir)

    assert dir.TIMEOUT == 5.0
    assert len(mock_warnings) == 1
    assert "Invalid SWANLAB_FS_TIMEOUT" in mock_warnings[0]

    monkeypatch.delenv("SWANLAB_FS_TIMEOUT", raising=False)
    importlib.reload(dir)


def test_timeout_env_negative_number(monkeypatch):
    """测试环境变量为负数或0时，是否能安全回退到默认值 5.0"""
    monkeypatch.setenv("SWANLAB_FS_TIMEOUT", "-10.5")

    mock_warnings = []
    monkeypatch.setattr("swanlab.sdk.internal.pkg.fs.dir.console.warning", lambda msg: mock_warnings.append(msg))

    importlib.reload(dir)

    assert dir.TIMEOUT == 5.0
    assert len(mock_warnings) == 1
    assert "must be > 0" in mock_warnings[0]

    monkeypatch.delenv("SWANLAB_FS_TIMEOUT", raising=False)
    importlib.reload(dir)
