"""
@author: nexisato
@file: test_cpu.py
@time: 2025-02-28 21:22:50
@description: 测试CPU信息采集
"""

import pytest
import platform
from swanlab.data.run.metadata.hardware.cpu import (
    get_cpu_brand_windows,
    get_cpu_brand_linux,
    get_cpu_info,
    CpuCollector,
)


@pytest.mark.skipif(platform.system() != "Linux", reason="Linux only")
def test_get_cpu_brand_linux():
    """
    获取不同linux系统 + 不同语言环境下的cpu品牌
    - Debian/Ubuntu/EulerOS
    """
    assert get_cpu_brand_linux() is not None


@pytest.mark.skipif(platform.system() != "Windows", reason="Windows only")
def test_get_cpu_brand_windows():
    """
    获取不同windows系统 + 不同语言环境下的cpu品牌
    - Windows 10/11
    """
    assert get_cpu_brand_windows() is not None


@pytest.mark.skipif(platform.system() != "Linux" or platform.system() != "Windows", reason="Linux or Windows only")
def test_get_cpu_info():
    """
    获取cpu信息
    """
    info, _ = get_cpu_info()
    assert info is not None and info["cores"] > 0
    assert info["brand"] is not None
