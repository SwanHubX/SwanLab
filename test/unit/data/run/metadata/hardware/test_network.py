import pytest
from swanlab.data.run.metadata.hardware.network import (
    get_network_info,
    NetworkCollector,   
)

def test_get_network_info():
    """
    获取网络信息
    """
    assert get_network_info() is not None