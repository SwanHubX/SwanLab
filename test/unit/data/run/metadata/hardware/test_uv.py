"""
@author: KashiwaByte
@file: test_uv.py
@time: 2025-04-14 13:28:37
@description: 测试uv
"""

import subprocess

import pytest
import time

try:
    has_uv = subprocess.run(["uv", "--version"], capture_output=True).returncode == 0
except FileNotFoundError:
    has_uv = False


@pytest.mark.skipif(not has_uv, reason="uv is not installed")
def test_uv():
    """
    获取uv信息，如果不存在则返回None
    Returns: str
    """
    start_time = time.time()
    result = subprocess.run(["uv", "pip", "list", "--format=freeze"], capture_output=True, text=True, timeout=15)
    end_time = time.time()
    elapsed_time = end_time - start_time

    # 检查命令是否成功运行
    if result.returncode == 0:
        # 使用 assert 来验证输出是否符合预期
        assert elapsed_time < 0.5
        assert isinstance(result.stdout, str)
    else:
        assert result.stdout is None
