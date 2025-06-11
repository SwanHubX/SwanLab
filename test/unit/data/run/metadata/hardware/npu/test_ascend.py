import pytest

from swanlab.data.run.metadata.hardware.npu.ascend import (
    AscendCollector,
    get_cann_version,
    get_chip_usage,
    get_version,
    map_npu,
)

try:
    npu_version = get_version()
    cann_version = get_cann_version()
    npu_map = map_npu()
    hbm_value: int = 0
    for npu_id in npu_map:
        for chip_id in npu_map[npu_id]:
            chip_info = npu_map[npu_id][chip_id]
            usage = get_chip_usage(npu_id, chip_id)
            hbm_value = int(usage.get("hbm", 0)) if usage else 0
    hbm_value = hbm_value * 1024
except Exception:
    npu_version = None
    cann_version = None
    npu_map = {}
    hbm_value = 0


@pytest.mark.skipif(npu_version is None, reason="No Ascend NPU found")
def test_ascend_collector():
    collector = AscendCollector(npu_map, hbm_value)
    collector()
    assert collector.collect_num == 1
