import pytest
from swanlab.data.run.metadata.hardware.disk import (
    get_disk_info,
    DiskCollector,   
)

def test_get_disk_info():
    """
    获取磁盘信息
    """
    assert get_disk_info() is not None