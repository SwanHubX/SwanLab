"""
@author: cunyue
@file: test_cpu.py
@time: 2026/3/31
@description: Tests for generic CPU hardware vendor
"""

import sys
from unittest.mock import MagicMock, mock_open, patch

import pytest

from swanlab.sdk.internal.run.system.hardware_vendor.cpu import CPU
from swanlab.sdk.typings.run.system import CPUSnapshot


class TestCPUGet:
    """Tests for CPU.get()"""

    # ──────────────────────────────────────────────
    # Linux
    # ──────────────────────────────────────────────

    @pytest.mark.skipif(sys.platform != "linux", reason="Linux only")
    def test_get_linux_success(self):
        """Linux: 正常采集物理核心数（单路）"""
        lscpu_output = "# comment\n0,0\n1,0\n2,0\n3,0\n"
        with (
            patch("subprocess.run") as mock_run,
            patch("multiprocessing.cpu_count", return_value=8),
            patch.object(CPU, "_get_real_brand", return_value="Intel Core i7-9700K"),
        ):
            mock_run.return_value = MagicMock(returncode=0, stdout=lscpu_output)
            result = CPU().get()

        assert isinstance(result, CPUSnapshot)
        assert result.brand == "Intel Core i7-9700K"
        assert result.logical_count == 8
        assert result.physical_count == 4

    @pytest.mark.skipif(sys.platform != "linux", reason="Linux only")
    def test_get_linux_multi_socket(self):
        """Linux: 多路 CPU 场景，Core ID 在不同 Socket 间重复，需要 (Core,Socket) 对去重"""
        # 两路 CPU，每路 4 个核，Core ID 均为 0-3
        lscpu_output = "# comment\n0,0\n1,0\n2,0\n3,0\n0,1\n1,1\n2,1\n3,1\n"
        with (
            patch("subprocess.run") as mock_run,
            patch("multiprocessing.cpu_count", return_value=16),
            patch.object(CPU, "_get_real_brand", return_value="Intel Xeon"),
        ):
            mock_run.return_value = MagicMock(returncode=0, stdout=lscpu_output)
            result = CPU().get()

        assert result is not None
        assert result.physical_count == 8  # 4 cores × 2 sockets

    @pytest.mark.skipif(sys.platform != "linux", reason="Linux only")
    def test_get_linux_lscpu_fails(self):
        """Linux: lscpu 返回非零退出码，physical_count 应为 None"""
        with (
            patch("subprocess.run") as mock_run,
            patch("multiprocessing.cpu_count", return_value=4),
            patch.object(CPU, "_get_real_brand", return_value=None),
            patch("platform.processor", return_value="x86_64"),
        ):
            mock_run.return_value = MagicMock(returncode=1, stdout="")
            result = CPU().get()

        assert result is not None
        assert result.physical_count is None
        assert result.logical_count == 4

    # ──────────────────────────────────────────────
    # macOS
    # ──────────────────────────────────────────────

    @pytest.mark.skipif(sys.platform != "darwin", reason="macOS only")
    def test_get_darwin_success(self):
        """macOS: 正常采集物理核心数"""
        with (
            patch("subprocess.run") as mock_run,
            patch("multiprocessing.cpu_count", return_value=8),
            patch.object(CPU, "_get_real_brand", return_value="Apple M1 Pro"),
        ):
            mock_run.return_value = MagicMock(returncode=0, stdout="6\n")
            result = CPU().get()

        assert isinstance(result, CPUSnapshot)
        assert result.physical_count == 6
        assert result.logical_count == 8

    @pytest.mark.skipif(sys.platform != "darwin", reason="macOS only")
    def test_get_darwin_sysctl_fails(self):
        """macOS: sysctl 失败，physical_count 应为 None"""
        with (
            patch("subprocess.run") as mock_run,
            patch("multiprocessing.cpu_count", return_value=8),
            patch.object(CPU, "_get_real_brand", return_value="Some CPU"),
        ):
            mock_run.return_value = MagicMock(returncode=1, stdout="")
            result = CPU().get()

        assert result is not None
        assert result.physical_count is None

    # ──────────────────────────────────────────────
    # Windows
    # ──────────────────────────────────────────────

    @pytest.mark.skipif(sys.platform != "win32", reason="Windows only")
    def test_get_windows_success(self):
        """Windows: 正常采集物理核心数（单 CPU）"""
        wmic_output = "NumberOfCores\n8\n\n"
        with (
            patch("subprocess.run") as mock_run,
            patch("multiprocessing.cpu_count", return_value=16),
            patch.object(CPU, "_get_real_brand", return_value="Intel Core i9"),
        ):
            mock_run.return_value = MagicMock(returncode=0, stdout=wmic_output)
            result = CPU().get()

        assert isinstance(result, CPUSnapshot)
        assert result.physical_count == 8
        assert result.logical_count == 16

    @pytest.mark.skipif(sys.platform != "win32", reason="Windows only")
    def test_get_windows_multi_cpu(self):
        """Windows: 多 CPU 场景，累加所有物理核心数"""
        wmic_output = "NumberOfCores\n4\n6\n\n"
        with (
            patch("subprocess.run") as mock_run,
            patch("multiprocessing.cpu_count", return_value=20),
            patch.object(CPU, "_get_real_brand", return_value="Intel Xeon"),
        ):
            mock_run.return_value = MagicMock(returncode=0, stdout=wmic_output)
            result = CPU().get()

        assert result is not None
        assert result.physical_count == 10  # 4 + 6

    @pytest.mark.skipif(sys.platform != "win32", reason="Windows only")
    def test_get_windows_wmic_fails(self):
        """Windows: wmic 失败，physical_count 应为 None"""
        with (
            patch("subprocess.run") as mock_run,
            patch("multiprocessing.cpu_count", return_value=8),
            patch.object(CPU, "_get_real_brand", return_value="Some CPU"),
        ):
            mock_run.return_value = MagicMock(returncode=1, stdout="")
            result = CPU().get()

        assert result is not None
        assert result.physical_count is None

    # ──────────────────────────────────────────────
    # Brand fallback（跨平台）
    # ──────────────────────────────────────────────

    def test_get_brand_fallback_to_platform_processor(self):
        """_get_real_brand 返回 None 时，回退到 platform.processor()"""
        with (
            patch("sys.platform", "linux"),
            patch("subprocess.run") as mock_run,
            patch("multiprocessing.cpu_count", return_value=4),
            patch.object(CPU, "_get_real_brand", return_value=None),
            patch("platform.processor", return_value="x86_64"),
        ):
            mock_run.return_value = MagicMock(returncode=1, stdout="")
            result = CPU().get()

        assert result is not None
        assert result.brand == "x86_64"

    def test_get_brand_fallback_to_none_when_both_empty(self):
        """_get_real_brand 和 platform.processor() 均为空字符串时，brand 应为 None"""
        with (
            patch("sys.platform", "linux"),
            patch("subprocess.run") as mock_run,
            patch("multiprocessing.cpu_count", return_value=4),
            patch.object(CPU, "_get_real_brand", return_value=None),
            patch("platform.processor", return_value=""),
        ):
            mock_run.return_value = MagicMock(returncode=1, stdout="")
            result = CPU().get()

        assert result is not None
        assert result.brand is None

    def test_get_returns_cpu_snapshot_type(self):
        """get() 在正常情况下必须返回 CPUSnapshot 实例"""
        with (
            patch("sys.platform", "linux"),
            patch("subprocess.run") as mock_run,
            patch("multiprocessing.cpu_count", return_value=4),
            patch.object(CPU, "_get_real_brand", return_value="Test CPU"),
        ):
            mock_run.return_value = MagicMock(returncode=0, stdout="0,0\n1,0\n")
            result = CPU().get()

        assert isinstance(result, CPUSnapshot)

    def test_get_returns_none_on_unexpected_exception(self):
        """顶层异常（如 cpu_count 抛出）时，get() 应返回 None 而不是抛出"""
        with (
            patch("multiprocessing.cpu_count", side_effect=RuntimeError("mocked failure")),
        ):
            result = CPU().get()

        assert result is None


class TestCPUGetRealBrand:
    """Tests for CPU._get_real_brand()"""

    # ──────────────────────────────────────────────
    # Linux
    # ──────────────────────────────────────────────

    @pytest.mark.skipif(sys.platform != "linux", reason="Linux only")
    def test_get_real_brand_linux_success(self):
        """Linux: 从 /proc/cpuinfo 正确读取 model name"""
        cpuinfo = "processor\t: 0\nmodel name\t: Intel(R) Core(TM) i7-9700K CPU @ 3.60GHz\ncpu MHz\t\t: 3600.0\n"
        with patch("builtins.open", mock_open(read_data=cpuinfo)):
            result = CPU._get_real_brand()

        assert result == "Intel(R) Core(TM) i7-9700K CPU @ 3.60GHz"

    @pytest.mark.skipif(sys.platform != "linux", reason="Linux only")
    def test_get_real_brand_linux_no_model_name(self):
        """Linux: /proc/cpuinfo 中不含 model name 时，返回 None"""
        cpuinfo = "processor\t: 0\ncpu MHz\t\t: 3600.0\n"
        with patch("builtins.open", mock_open(read_data=cpuinfo)):
            result = CPU._get_real_brand()

        assert result is None

    @pytest.mark.skipif(sys.platform != "linux", reason="Linux only")
    def test_get_real_brand_linux_file_error(self):
        """Linux: /proc/cpuinfo 打开失败时，返回 None"""
        with patch("builtins.open", side_effect=OSError("permission denied")):
            result = CPU._get_real_brand()

        assert result is None

    # ──────────────────────────────────────────────
    # macOS
    # ──────────────────────────────────────────────

    @pytest.mark.skipif(sys.platform != "darwin", reason="macOS only")
    def test_get_real_brand_darwin_success(self):
        """macOS: 从 sysctl 正确读取 CPU 品牌"""
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stdout="Apple M2 Pro\n")
            result = CPU._get_real_brand()

        assert result == "Apple M2 Pro"

    @pytest.mark.skipif(sys.platform != "darwin", reason="macOS only")
    def test_get_real_brand_darwin_sysctl_fails(self):
        """macOS: sysctl 返回非零时，返回 None"""
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=1, stdout="")
            result = CPU._get_real_brand()

        assert result is None

    # ──────────────────────────────────────────────
    # Windows
    # ──────────────────────────────────────────────

    @pytest.mark.skipif(sys.platform != "win32", reason="Windows only")
    def test_get_real_brand_windows_success(self):
        """Windows: 从注册表正确读取 CPU 名称"""
        mock_key = MagicMock()
        with (
            patch("winreg.OpenKey", return_value=mock_key),
            patch("winreg.QueryValueEx", return_value=("Intel(R) Core(TM) i9-10900K CPU @ 3.70GHz  ", None)),
        ):
            result = CPU._get_real_brand()

        assert result == "Intel(R) Core(TM) i9-10900K CPU @ 3.70GHz"

    @pytest.mark.skipif(sys.platform != "win32", reason="Windows only")
    def test_get_real_brand_windows_registry_error(self):
        """Windows: 注册表访问失败时，返回 None"""
        with patch("winreg.OpenKey", side_effect=OSError("access denied")):
            result = CPU._get_real_brand()

        assert result is None

    # ──────────────────────────────────────────────
    # 跨平台：未知平台时返回 None
    # ──────────────────────────────────────────────

    def test_get_real_brand_unknown_platform_returns_none(self):
        """未知平台（如 freebsd）时，_get_real_brand 应返回 None"""
        with patch("sys.platform", "freebsd12"):
            result = CPU._get_real_brand()

        assert result is None
