import pytest

from swanlab.data.run.metadata.hardware.dcu.hygon import DCUCollector, map_hygon_dcu

try:
    driver_version, dcu_map = map_hygon_dcu()
    max_mem_value = 0
    for dcu_id in dcu_map:
        mem_value = int(dcu_map[dcu_id]["memory"][:-2])
        max_mem_value = max(max_mem_value, mem_value)
    max_mem_value = max_mem_value * 1024
except Exception as e:
    driver_version = None
    dcu_map = None
    max_mem_value = None


@pytest.mark.skipif(driver_version is None, reason="Hygon DCU not found")
def test_hygon_dcu_info():
    collector = DCUCollector(dcu_map, max_mem_value)
    collector()
    assert collector.collect_num == 1
