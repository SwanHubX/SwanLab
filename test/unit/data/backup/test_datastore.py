"""
@author: cunyue
@file: test_datastore.py
@time: 2025/6/5 18:33
@description: 测试 DataStore 类的功能，包括数据写入、读取、校验和计算等
"""

import os.path

from nanoid import generate

from swanlab.data.porter.datastore import DataStore
from tutils import TEMP_PATH

logs = [generate(size=l) for l in range(1, 100001, 1000)]


def test_write(filename=os.path.join(TEMP_PATH, "backup.swanlab")):
    """
    测试文件写入功能
    """
    ds = DataStore()
    ds.open_for_write(filename)
    for log in logs:
        ds.write(log)


def test_scan(filename=os.path.join(TEMP_PATH, "backup.swanlab")):
    """
    测试文件扫描功能
    """
    # 首先写入
    test_write(filename)
    # 然后扫描
    ds = DataStore()
    ds.open_for_scan(filename)
    for i in range(len(logs)):
        log = ds.scan()
        assert log == logs[i], "Error: Scanned log does not match written log"
