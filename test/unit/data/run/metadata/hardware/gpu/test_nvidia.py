"""
@author: cunyue
@file: test_nvidia.py
@time: 2024/12/5 13:27
@description: 测试NVIDIA GPU信息采集
"""

import pynvml
import pytest

from swanlab.data.run.metadata.hardware.gpu.nvidia import GpuCollector

try:
    pynvml.nvmlInit()
    count = pynvml.nvmlDeviceGetCount()
except Exception:  # noqa
    count = 0


@pytest.mark.skipif(count == 0, reason="No NVIDIA GPU found")
def test_before_impl():
    collector = GpuCollector(count)
    collector.before_collect_impl()
    assert len(collector.handles) == count


@pytest.mark.skipif(count == 0, reason="No NVIDIA GPU found")
def test_after_impl():
    collector = GpuCollector(count)
    collector.after_collect_impl()
    assert len(collector.handles) == 0


class TestGpuCollector:
    collector = GpuCollector(count)

    def setup_class(self):
        self.collector.before_collect_impl()

    def teardown_class(self):
        self.collector.after_collect_impl()

    @pytest.mark.skipif(count == 0, reason="No NVIDIA GPU found")
    def test_get_mem(self):
        # 获取handle
        idx = 0
        mem = self.collector.get_gpu_mem_pct(idx=idx)
        assert mem['name'] == "GPU 0 Memory Allocated (%)"
        assert mem['config'].y_range == (0, 100)
        assert mem['config'].metric_name == "GPU 0"
        assert 100 >= mem['value'] >= 0

    @pytest.mark.skipif(count == 0, reason="No NVIDIA GPU found")
    def test_get_utils(self):
        # 获取handle
        idx = 0
        utils = self.collector.get_gpu_util(idx=idx)
        assert utils['name'] == "GPU 0 Utilization (%)"
        assert utils['config'].y_range == (0, 100)
        assert utils['config'].metric_name == "GPU 0"
        assert 100 >= utils['value'] >= 0

    @pytest.mark.skipif(count == 0, reason="No NVIDIA GPU found")
    def test_get_temp(self):
        # 获取handle
        idx = 0
        temp = self.collector.get_gpu_temp(idx=idx)
        assert temp['name'] == "GPU 0 Temperature (℃)"
        assert temp['config'].y_range is None
        assert temp['config'].metric_name == "GPU 0"
        assert temp['value'] >= 0
        assert temp['config'].metric_name == "GPU 0"

    @pytest.mark.skipif(count == 0, reason="No NVIDIA GPU found")
    def test_get_power(self):
        # 获取handle
        idx = 0
        power = self.collector.get_gpu_power(idx=idx)
        assert power['name'] == "GPU 0 Power Usage (W)"
        assert power['config'].y_range is None
        assert power['config'].metric_name == "GPU 0"
        assert power['value'] >= 0
        assert power['config'].metric_name == "GPU 0"
