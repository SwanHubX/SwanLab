"""
@author: cunyue
@file: test_helper.py
@time: 2024/12/6 13:52
@description: 测试工具函数
"""

from swanlab.data.run.helper import check_log_level


def test_check_log_level():
    assert check_log_level(None) == "info"
    assert check_log_level("debug") == "debug"
    assert check_log_level("info") == "info"
    assert check_log_level("warning") == "warning"
    assert check_log_level("error") == "error"
    assert check_log_level("critical") == "critical"
    assert check_log_level("not_exist") == "info"
    assert check_log_level("INFO") == "info"
