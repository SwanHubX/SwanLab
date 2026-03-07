"""
@author: cunyue
@file: helper.py
@time: 2026/3/7 17:49
@description: SwanLab 运行时客户端辅助函数
"""

import json
from typing import Dict, List, Union

import requests


def decode_response(resp: requests.Response) -> Union[Dict, List, str]:
    """
    解码响应，合并异常捕获
    """
    try:
        return resp.json()
    except (json.decoder.JSONDecodeError, requests.JSONDecodeError):
        return resp.text
