"""
@author: cunyue
@file: utils.py
@time: 2025/2/13 15:31
@description: 存放一些公用工具
"""

import traceback

from swanlab.log import swanlog


def error_print(tp):
    """
    错误打印
    """
    # 如果是KeyboardInterrupt异常
    if tp == KeyboardInterrupt:
        swanlog.info("KeyboardInterrupt by user")
    else:
        swanlog.info("Error happened while training")


def traceback_error(tb, val):
    """
    获取traceback信息
    """
    trace_list = traceback.format_tb(tb)
    html = ""
    for line in trace_list:
        html += line
    html += str(val)
    return html
