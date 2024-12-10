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
def test_get_mem():
    collector = GpuCollector()
    # 获取handle
    idx = 0
    handle = pynvml.nvmlDeviceGetHandleByIndex(idx)
    mem = collector.get_gpu_mem_pct(idx=idx, handle=handle)
    assert  mem['name'] == "GPU 0 Utilization (%)"
    assert  mem['config'].y_range == (0, 100)
    assert mem['config'].metric_name == "GPU 0"
    assert  100>=mem['value'] >= 0

@pytest.mark.skipif(count == 0, reason="No NVIDIA GPU found")
def test_get_temp():
    collector = GpuCollector()
    # 获取handle
    idx = 0
    handle = pynvml.nvmlDeviceGetHandleByIndex(idx)
    temp = collector.get_gpu_temp(idx=idx, handle=handle)
    assert temp['name'] == "GPU 0 Temperature (℃)"
    assert temp['config'].y_range is None
    assert temp['config'].metric_name == "GPU 0"
    assert  temp['value'] >= 0
    assert temp['config'].metric_name == "GPU 0"

@pytest.mark.skipif(count == 0, reason="No NVIDIA GPU found")
def test_get_power():
    collector = GpuCollector()
    # 获取handle
    idx = 0
    handle = pynvml.nvmlDeviceGetHandleByIndex(idx)
    power = collector.get_gpu_power(idx=idx, handle=handle)
    assert power['name'] == "GPU 0 Power Usage (W)"
    assert power['config'].y_range is None
    assert power['config'].metric_name == "GPU 0"
    assert  power['value'] >= 0
    assert power['config'].metric_name == "GPU 0"