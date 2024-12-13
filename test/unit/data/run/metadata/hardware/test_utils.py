"""
@author: cunyue
@file: test_utils.py
@time: 2024/12/9 20:55
@description: 测试硬件信息采集工具
"""

from swanlab.data.run.metadata.hardware.utils import random_index, CpuCollector, generate_key


def test_random_index():
    s = random_index()
    assert len(s) == 8
    assert s.isalnum()
    s = random_index(10)
    assert len(s) == 10
    assert s.isalnum()


def test_generate_key():
    s = generate_key("test")
    assert s == "__swanlab__.test"


def test_cpu_usage():
    c = CpuCollector()
    usage = c.get_cpu_usage()
    assert usage is not None
    assert 0 <= usage["value"] <= 100
    assert usage["key"].endswith("cpu.pct")
    assert usage["name"] == "CPU Utilization (%)"
    assert usage["config"].y_range == (0, 100)
    assert usage["config"].chart_name == "CPU Utilization (%)"
    assert usage["config"].chart_index is None


def test_per_cpu_usage():
    c = CpuCollector()
    usage = c.get_per_cpu_usage()
    assert usage is not None
    for idx, u in enumerate(usage):
        assert 0 <= u["value"] <= 100
        assert u["key"].endswith(f"cpu.{idx}.pct")
        assert u["name"] == f"CPU {idx} Utilization (%)"
        assert u["config"].y_range == (0, 100)
        assert u["config"].chart_name == f"CPU Utilization (per core) (%)"
        # 每个核心的index应该相同,因为必须要放在同一个图表中
        assert u["config"].metric_name == f"CPU {idx}"
