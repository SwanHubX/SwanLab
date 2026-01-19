"""
@author: Zhou QiYang
@file: utils.py
@time: 2026/1/19 15:23
@description: 用户相关接口的工具函数
"""

import requests

STATUS_OK = 200
STATUS_CREATED = 201


def check_deleted(resp: requests.Response):
    return resp.status_code == STATUS_OK


def check_created(resp: requests.Response):
    return resp.status_code == STATUS_CREATED
