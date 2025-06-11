import pytest

from swanlab.data.run.metadata.hardware.mlu.cambricon import CambriconCollector, map_mlu

try:
    driver, mlu_map = map_mlu()
    max_mem_value = 0
    for mlu_id in mlu_map:
        mlu_mem = int(mlu_map[mlu_id].get("memory", 0))
        max_mem_value = max(max_mem_value, mlu_mem)
    max_mem_value = max_mem_value * 1024
    collector = CambriconCollector(mlu_map, max_mem_value)
except Exception as e:
    driver = None
    mlu_map = None
    max_mem_value = None


@pytest.mark.skipif(driver is None, reason="Cambricon MLU not found")
def test_cambricon_mlu_info():
    collector = CambriconCollector(mlu_map, max_mem_value)
    collector()
    assert collector.collect_num == 1
