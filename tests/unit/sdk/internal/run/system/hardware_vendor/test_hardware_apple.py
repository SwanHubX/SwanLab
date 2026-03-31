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
from swanlab.sdk.typings.run.system import AppleSiliconSnapshot


class TestApple:
    """Test Apple Silicon hardware detection"""

    @pytest.mark.skipif(sys.platform != "darwin", reason="macOS only")
    def test_get_on_macos_with_apple_silicon(self):
        """Test successful detection on macOS with Apple Silicon"""
        mock_output = {"SPHardwareDataType": [{"chip_type": "Apple M3 Pro", "physical_memory": "18 GB"}]}
        with patch("subprocess.run") as mock_run, patch("multiprocessing.cpu_count", return_value=12):
            mock_run.return_value = MagicMock(stdout=json.dumps(mock_output))
            result = Apple().get()
            assert isinstance(result, AppleSiliconSnapshot)
            assert result.name == "Apple M3 Pro"
            assert result.memory == 18
            assert result.memory_unit == "GB"
            assert result.cpu_count == 12

    def test_get_on_non_macos(self):
        """Test returns None on non-macOS platforms"""
        with patch("sys.platform", "linux"):
            result = Apple().get()
            assert result is None

    @pytest.mark.skipif(sys.platform != "darwin", reason="macOS only")
    def test_get_with_intel_mac(self):
        """Test detection on Intel Mac (fallback to cpu_type)"""
        mock_output = {"SPHardwareDataType": [{"cpu_type": "Intel Core i7", "physical_memory": "16 GB"}]}
        with patch("subprocess.run") as mock_run, patch("multiprocessing.cpu_count", return_value=8):
            mock_run.return_value = MagicMock(stdout=json.dumps(mock_output))
            result = Apple().get()
            assert isinstance(result, AppleSiliconSnapshot)
            assert result.name == "Intel Core i7"
            assert result.memory == 16

    @pytest.mark.skipif(sys.platform != "darwin", reason="macOS only")
    def test_get_with_missing_chip_type(self):
        """Test returns None when chip_type is missing"""
        mock_output = {"SPHardwareDataType": [{"physical_memory": "16 GB"}]}
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(stdout=json.dumps(mock_output))
            result = Apple().get()
            assert result is None
