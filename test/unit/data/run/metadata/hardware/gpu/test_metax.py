import pytest

from swanlab.data.run.metadata.hardware.gpu.metax import MetaxCollector, map_metax_gpu

try:
    driver, maca_version, gpu_map = map_metax_gpu()
    max_mem_value = 0
    for gpu_id in gpu_map:
        mem_value = gpu_map[gpu_id]["memory"]
        if mem_value > max_mem_value:
            max_mem_value = mem_value
    max_mem_value = max_mem_value * 1024
except Exception as e:
    driver = None
    maca_version = None
    gpu_map = None
    max_mem_value = None


@pytest.mark.skipif(maca_version is None, reason="MetaX GPU not found")
def test_metax_gpu_info():
    collector = MetaxCollector(gpu_map, max_mem_value)
    collector()
    assert collector.collect_num == 1
