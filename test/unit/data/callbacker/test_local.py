"""
@author: cunyue
@file: test_local.py
@time: 2025/5/18 17:24
@description: 测试本地回调器
"""

import os.path

from freezegun import freeze_time
from nanoid import generate
from swankit.core import SwanLabSharedSettings

import tutils as T
from swanlab.data.callbacker.local import LocalRunCallback
from swanlab.log.type import LogData


def test_local_write_handler(monkeypatch):
    """
    测试本地回调器的写入处理函数，当日期改变时，文件名也会改变
    """
    callback = LocalRunCallback()
    callback.settings = SwanLabSharedSettings(
        logdir=T.TEMP_PATH,
        run_id="test-local-write-handler",
        version="test-local-write-handler",
        should_save=True,
    )
    with freeze_time('2020-10-06'):
        a = generate()
        mockdata = LogData(
            type='stdout',
            contents=[{'message': a, "create_time": '12344', "epoch": 1}],
        )
        callback._write_handler(mockdata)
        filename = os.path.join(callback.settings.console_dir, "2020-10-06.log")
        assert os.path.exists(filename)
        with open(filename, "r") as f:
            content = f.readlines()
            assert content[-1] == a + '\n'
    with freeze_time('2020-10-07'):
        b = generate()
        mockdata = LogData(
            type='stdout',
            contents=[{'message': b, "create_time": '12344', "epoch": 1}],
        )
        callback._write_handler(mockdata)
        filename = os.path.join(callback.settings.console_dir, "2020-10-07.log")
        assert os.path.exists(filename)
        with open(filename, "r") as f:
            content = f.readlines()
            assert content[-1] == b + '\n'
