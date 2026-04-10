"""
@author: cunyue
@file: test_apple.py
@time: 2026/3/31
@description: Tests for Apple Silicon hardware vendor
"""

import json
import sys
from unittest.mock import MagicMock, patch

import pytest

from swanlab.sdk.internal.run.system.hardware_vendor.apple import Apple
from swanlab.sdk.typings.run.system import AppleSiliconSnapshot, PlatformSlug, SystemShim


class TestAppleGet:
    """Test Apple.get() — 静态快照采集"""

    @pytest.mark.skipif(sys.platform != "darwin", reason="macOS only")
    def test_get_on_macos_with_apple_silicon(self):
        """Test successful detection on macOS with Apple Silicon"""
        mock_output = {"SPHardwareDataType": [{"chip_type": "Apple M3 Pro", "physical_memory": "18 GB"}]}
        with patch("subprocess.run") as mock_run, patch("multiprocessing.cpu_count", return_value=12):
            mock_run.return_value = MagicMock(stdout=json.dumps(mock_output))
            result = Apple.get()
            assert isinstance(result, AppleSiliconSnapshot)
            assert result.name == "Apple M3 Pro"
            assert result.memory == 18
            assert result.memory_unit == "GB"
            assert result.cpu_count == 12

    def test_get_on_non_macos(self):
        """Test returns None on non-macOS platforms"""
        with patch("sys.platform", "linux"):
            result = Apple.get()
            assert result is None

    @pytest.mark.skipif(sys.platform != "darwin", reason="macOS only")
    def test_get_with_intel_mac(self):
        """Test detection on Intel Mac (fallback to cpu_type)"""
        mock_output = {"SPHardwareDataType": [{"cpu_type": "Intel Core i7", "physical_memory": "16 GB"}]}
        with patch("subprocess.run") as mock_run, patch("multiprocessing.cpu_count", return_value=8):
            mock_run.return_value = MagicMock(stdout=json.dumps(mock_output))
            result = Apple.get()
            assert isinstance(result, AppleSiliconSnapshot)
            assert result.name == "Intel Core i7"
            assert result.memory == 16

    @pytest.mark.skipif(sys.platform != "darwin", reason="macOS only")
    def test_get_with_missing_chip_type(self):
        """Test returns None when chip_type is missing"""
        mock_output = {"SPHardwareDataType": [{"physical_memory": "16 GB"}]}
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=json.dumps(mock_output))
            result = Apple.get()
            assert result is None


class TestAppleNew:
    """Test Apple.new() — 监控采集器工厂"""

    @staticmethod
    def _make_shim(
        slug: PlatformSlug = "macos-arm",
    ) -> SystemShim:
        return SystemShim(slug=slug, enable_cpu=True, enable_memory=True)

    def test_returns_none_on_linux(self):
        """Linux 平台返回 None"""
        shim = self._make_shim(slug="linux")
        assert Apple.new(shim) is None

    def test_returns_none_on_macos_intel(self):
        """macOS Intel 平台返回 None"""
        shim = self._make_shim(slug="macos-intel")
        assert Apple.new(shim) is None

    def test_returns_none_on_windows(self):
        """Windows 平台返回 None"""
        shim = self._make_shim(slug="windows")
        assert Apple.new(shim) is None

    def test_returns_none_on_unknown(self):
        """未知平台返回 None"""
        shim = self._make_shim(slug="unknown")
        assert Apple.new(shim) is None

    def test_returns_tuple_on_macos_arm(self):
        """macOS-arm 平台返回 (Apple, SystemScalars) 元组"""
        shim = self._make_shim(slug="macos-arm")
        result = Apple.new(shim)
        assert result is not None
        collector, scalars = result
        assert isinstance(collector, Apple)
        assert isinstance(scalars, list)
        assert len(scalars) == 4

    def test_scalar_keys(self):
        """验证采集器注册了正确的 scalar key"""
        shim = self._make_shim(slug="macos-arm")
        result = Apple.new(shim)
        assert result is not None
        _, scalars = result
        keys = [s.key for s in scalars]
        assert keys == ["cpu.pct", "cpu.thds", "mem.pct", "mem.proc"]

    def test_scalar_chart_names(self):
        """验证 scalar 的 chart_name 分组正确"""
        shim = self._make_shim(slug="macos-arm")
        result = Apple.new(shim)
        assert result is not None
        _, scalars = result
        chart_names = [s.chart_name for s in scalars]
        assert chart_names == ["CPU Utilization", "Process CPU Threads", "System Memory", "Process Memory"]

    def test_collect_returns_four_metrics(self):
        """collect() 返回 4 个 (key, value) 采集结果"""
        mock_proc = MagicMock()
        mock_proc.num_threads.return_value = 8
        mock_proc.memory_info.return_value = MagicMock(rss=1024 * 1024 * 256)  # 256 MB
        with (
            patch("psutil.cpu_percent", return_value=42.5),
            patch("psutil.Process", return_value=mock_proc),
            patch("psutil.virtual_memory") as mock_vm,
        ):
            mock_vm.return_value = MagicMock(percent=67.3)
            shim = self._make_shim(slug="macos-arm")
            result = Apple.new(shim)

        assert result is not None
        collector, _ = result
        # collect 时 handler 闭包中已捕获 mock_proc，需再次 patch cpu_percent 和 virtual_memory
        with (
            patch("psutil.cpu_percent", return_value=42.5),
            patch("psutil.virtual_memory", return_value=MagicMock(percent=67.3)),
        ):
            metrics = collector.collect()

        assert len(metrics) == 4
        assert metrics[0] == ("cpu.pct", 42.5)
        assert metrics[1] == ("cpu.thds", 8)
        assert metrics[2] == ("mem.pct", 67.3)
        assert metrics[3] == ("mem.proc", 256.0)
