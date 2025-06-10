import pytest

from swanlab.data.run.metadata.hardware.gpu.moorethreads import MTTCollector, map_moorethreads_gpu

try:
    driver, gpu_map = map_moorethreads_gpu()
    max_mem_value = 0
    for gpu_id in gpu_map:
        gpu_mem = int(gpu_map[gpu_id].get("memory", 0))
        max_mem_value = max(max_mem_value, gpu_mem)
    max_mem_value = max_mem_value * 1024
except Exception as e:
    driver = None
    gpu_map = None
    max_mem_value = None


@pytest.mark.skipif(driver is None, reason="MooreThreads GPU not found")
def test_moorethreads_gpu_info():
    collector = MTTCollector(gpu_map, max_mem_value)
    collector()
    assert collector.collect_num == 1
