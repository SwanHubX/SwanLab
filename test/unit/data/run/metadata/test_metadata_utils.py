"""
@author: cunyue
@file: test_metadata_utils.py
@time: 2025/9/9 12:39
@description: 测试工具函数
"""

import os

from swanlab.data.run.metadata.utils import check_env


@check_env('TEST_ENV', default_return='default_value')
def foo():
    return "bar"


def test_check_env():
    """
    测试check_env装饰器功能
    """
    value = foo()
    assert value == "bar"
    os.environ['TEST_ENV'] = 'true'
    value = foo()
    assert value == 'default_value'
    os.environ['TEST_ENV'] = 'TRUE'
    value = foo()
    assert value == 'default_value'
    os.environ['TEST_ENV'] = '0'
    value = foo()
    assert value == "bar"
    del os.environ['TEST_ENV']
    value = foo()
    assert value == "bar"
