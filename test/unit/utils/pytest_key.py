#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
@DATE: 2024/4/3 01:21
@File: pytest_key.py
@IDE: pycharm
@Description:
    token测试
"""
import os.path
import pytest
from swanlab.error import KeyFileError
from swanlab.utils.key import save_key, get_key
from utils.test.config import TEMP_PATH
from typing import Tuple
import netrc
import nanoid

netrc_pah = os.path.join(TEMP_PATH, ".netrc")


# ---------------------------------- 工具函数 ----------------------------------

def create_key() -> Tuple[str, str, str]:
    """
    生成保存key所必须的一些参数，host，key、password
    """
    return nanoid.generate(size=8), nanoid.generate(size=11), nanoid.generate()


def _test_save_key(info, path=netrc_pah):
    save_key(path, *info)


# ---------------------------------- 测试 ----------------------------------

def test_save_key_success():
    """
    测试保存key成功
    """
    info = create_key()
    _test_save_key(info)
    # netrc读取
    nrc = netrc.netrc(netrc_pah).authenticators(info[0])
    assert nrc[0] == info[1]
    assert nrc[1] == ''
    assert nrc[2] == info[2]


def test_save_key_path_error():
    """
    测试保存key时，输入的path上级文件夹不存在的情况
    """
    _netrc_pah = os.path.join(os.path.dirname(netrc_pah), nanoid.generate(), os.path.basename(netrc_pah))
    info = create_key()
    with pytest.raises(KeyFileError):
        _test_save_key(info, _netrc_pah)
