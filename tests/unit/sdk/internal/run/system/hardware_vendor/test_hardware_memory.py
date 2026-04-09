"""
@author: cunyue
@file: test_hardware_memory.py
@time: 2026/3/31
@description: Tests for generic Memory hardware vendor
"""

import sys
from unittest.mock import MagicMock, mock_open, patch

import pytest

from swanlab.sdk.internal.run.system.hardware_vendor.memory import Memory, _bytes_to_snapshot
from swanlab.sdk.typings.run.system import MemorySnapshot

_GB = 1024**3
_MB = 1024**2


class TestBytesToSnapshot:
    """Tests for _bytes_to_snapshot() — 纯函数，无需平台 mock"""

    def test_exact_1gb(self):
        """恰好 1 GB → total=1, unit=GB"""
        result = _bytes_to_snapshot(1 * _GB)
        assert result == MemorySnapshot(total=1, total_unit="GB")

    def test_16gb(self):
        """16 GB 整数转换"""
        result = _bytes_to_snapshot(16 * _GB)
        assert result == MemorySnapshot(total=16, total_unit="GB")

    def test_gb_truncates_remainder(self):
        """非整 GB 向下截断（1.9 GB → 1 GB）"""
        result = _bytes_to_snapshot(int(1.9 * _GB))
        assert result == MemorySnapshot(total=1, total_unit="GB")

    def test_less_than_1gb_uses_mb(self):
        """不足 1 GB 时降级到 MB"""
        result = _bytes_to_snapshot(512 * _MB)
        assert result == MemorySnapshot(total=512, total_unit="MB")

    def test_exactly_1023mb(self):
        """1023 MB → unit=MB"""
        result = _bytes_to_snapshot(1023 * _MB)
        assert result == MemorySnapshot(total=1023, total_unit="MB")

    def test_less_than_1mb_returns_none(self):
        """不足 1 MB（极端情况）→ None"""
        assert _bytes_to_snapshot(1024) is None

    def test_zero_returns_none(self):
        """0 字节 → None"""
        assert _bytes_to_snapshot(0) is None


class TestMemoryGet:
    """Tests for Memory.get()"""

    # ──────────────────────────────────────────────
    # Linux
    # ──────────────────────────────────────────────

    @pytest.mark.skipif(sys.platform != "linux", reason="Linux only")
    def test_get_linux_success(self):
        """Linux: 正常从 /proc/meminfo 采集总内存"""
        # 16 GB = 16777216 kB
        meminfo = "MemTotal:       16777216 kB\nMemFree:        8000000 kB\n"
        with patch("builtins.open", mock_open(read_data=meminfo)):
            result = Memory.get()

        assert isinstance(result, MemorySnapshot)
        assert result.total == 16
        assert result.total_unit == "GB"

    @pytest.mark.skipif(sys.platform != "linux", reason="Linux only")
    def test_get_linux_file_error(self):
        """Linux: /proc/meminfo 读取失败 → get() 返回 None"""
        with patch("builtins.open", side_effect=OSError("permission denied")):
            result = Memory.get()

        assert result is None

    # ──────────────────────────────────────────────
    # macOS
    # ──────────────────────────────────────────────

    @pytest.mark.skipif(sys.platform != "darwin", reason="macOS only")
    def test_get_darwin_success(self):
        """macOS: 正常从 sysctl hw.memsize 采集总内存"""
        # 16 GB in bytes
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout=str(16 * _GB) + "\n")
            result = Memory.get()

        assert isinstance(result, MemorySnapshot)
        assert result.total == 16
        assert result.total_unit == "GB"

    @pytest.mark.skipif(sys.platform != "darwin", reason="macOS only")
    def test_get_darwin_sysctl_fails(self):
        """macOS: sysctl 返回非零退出码 → get() 返回 None"""
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=1, stdout="")
            result = Memory.get()

        assert result is None

    # ──────────────────────────────────────────────
    # Windows
    # ──────────────────────────────────────────────

    @pytest.mark.skipif(sys.platform != "win32", reason="Windows only")
    def test_get_windows_success(self):
        """Windows: 正常从 wmic 采集总内存"""
        wmic_output = "TotalPhysicalMemory\n{}\n\n".format(16 * _GB)
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout=wmic_output)
            result = Memory.get()

        assert isinstance(result, MemorySnapshot)
        assert result.total == 16
        assert result.total_unit == "GB"

    @pytest.mark.skipif(sys.platform != "win32", reason="Windows only")
    def test_get_windows_wmic_fails(self):
        """Windows: wmic 返回非零退出码 → get() 返回 None"""
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=1, stdout="")
            result = Memory.get()

        assert result is None

    # ──────────────────────────────────────────────
    # 跨平台兜底逻辑
    # ──────────────────────────────────────────────

    def test_get_returns_none_when_total_bytes_is_none(self):
        """_get_total_bytes 返回 None 时，get() 应返回 None"""
        with patch.object(Memory, "_get_total_bytes", return_value=None):
            result = Memory.get()

        assert result is None

    def test_get_returns_none_when_total_bytes_is_zero(self):
        """_get_total_bytes 返回 0 时，get() 应返回 None"""
        with patch.object(Memory, "_get_total_bytes", return_value=0):
            result = Memory.get()

        assert result is None

    def test_get_returns_none_on_unexpected_exception(self):
        """_get_total_bytes 抛出异常时，get() 应返回 None 而不是抛出"""
        with patch.object(Memory, "_get_total_bytes", side_effect=RuntimeError("mocked failure")):
            result = Memory.get()

        assert result is None

    def test_get_returns_memory_snapshot_type(self):
        """正常路径下 get() 必须返回 MemorySnapshot 实例"""
        with patch.object(Memory, "_get_total_bytes", return_value=16 * _GB):
            result = Memory.get()

        assert isinstance(result, MemorySnapshot)


class TestMemoryGetTotalBytes:
    """Tests for Memory._get_total_bytes()"""

    # ──────────────────────────────────────────────
    # Linux
    # ──────────────────────────────────────────────

    @pytest.mark.skipif(sys.platform != "linux", reason="Linux only")
    def test_get_total_bytes_linux_success(self):
        """Linux: 正确解析 /proc/meminfo 的 MemTotal 行"""
        meminfo = "MemTotal:       16777216 kB\nMemFree:        8000000 kB\n"
        with patch("builtins.open", mock_open(read_data=meminfo)):
            result = Memory._get_total_bytes()

        assert result == 16777216 * 1024

    @pytest.mark.skipif(sys.platform != "linux", reason="Linux only")
    def test_get_total_bytes_linux_no_mem_total(self):
        """Linux: /proc/meminfo 中没有 MemTotal 行 → 返回 None"""
        meminfo = "MemFree:        8000000 kB\nBuffers:        512000 kB\n"
        with patch("builtins.open", mock_open(read_data=meminfo)):
            result = Memory._get_total_bytes()

        assert result is None

    @pytest.mark.skipif(sys.platform != "linux", reason="Linux only")
    def test_get_total_bytes_linux_file_error(self):
        """Linux: /proc/meminfo 打开失败 → 返回 None"""
        with patch("builtins.open", side_effect=OSError("no such file")):
            result = Memory._get_total_bytes()

        assert result is None

    # ──────────────────────────────────────────────
    # macOS
    # ──────────────────────────────────────────────

    @pytest.mark.skipif(sys.platform != "darwin", reason="macOS only")
    def test_get_total_bytes_darwin_success(self):
        """macOS: 正确解析 sysctl hw.memsize 输出（字节）"""
        bytes_val = 17179869184  # 16 GB
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout=str(bytes_val) + "\n")
            result = Memory._get_total_bytes()

        assert result == bytes_val

    @pytest.mark.skipif(sys.platform != "darwin", reason="macOS only")
    def test_get_total_bytes_darwin_sysctl_fails(self):
        """macOS: sysctl 返回非零 → 返回 None"""
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=1, stdout="")
            result = Memory._get_total_bytes()

        assert result is None

    # ──────────────────────────────────────────────
    # Windows
    # ──────────────────────────────────────────────

    @pytest.mark.skipif(sys.platform != "win32", reason="Windows only")
    def test_get_total_bytes_windows_success(self):
        """Windows: 正确解析 wmic TotalPhysicalMemory 输出（字节）"""
        bytes_val = 17179869184  # 16 GB
        wmic_output = "TotalPhysicalMemory\n{}\n\n".format(bytes_val)
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout=wmic_output)
            result = Memory._get_total_bytes()

        assert result == bytes_val

    @pytest.mark.skipif(sys.platform != "win32", reason="Windows only")
    def test_get_total_bytes_windows_wmic_fails(self):
        """Windows: wmic 返回非零 → 返回 None"""
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=1, stdout="")
            result = Memory._get_total_bytes()

        assert result is None

    @pytest.mark.skipif(sys.platform != "win32", reason="Windows only")
    def test_get_total_bytes_windows_empty_output(self):
        """Windows: wmic 仅有表头无数据行 → 返回 None"""
        wmic_output = "TotalPhysicalMemory\n\n"
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout=wmic_output)
            result = Memory._get_total_bytes()

        assert result is None

    # ──────────────────────────────────────────────
    # 跨平台
    # ──────────────────────────────────────────────

    def test_get_total_bytes_unknown_platform_returns_none(self):
        """未知平台 → 返回 None"""
        with patch("sys.platform", "freebsd12"):
            result = Memory._get_total_bytes()

        assert result is None
