"""
@author: cunyue
@file: webhook.py
@time: 2024/11/18 12:35
@description: 当swanlab初始化完毕时，会调用这个webhook函数，用于通知用户消息
"""

import os

import requests

from swanlab.data.run.metadata import get_cooperation_info
from swanlab.env import SwanLabEnv
from swanlab.log import swanlog


def try_send_webhook():
    """
    尝试发送消息，如果发送失败也不会抛出异常，但是会打印错误信息
    """
    webhook = os.getenv(SwanLabEnv.WEBHOOK.value)
    if not webhook:
        return
    data = get_cooperation_info()
    try:
        requests.post(webhook, json=data)
    except Exception as e:  # noqa
        swanlog.error(f"Failed to send message to webhook: {e}")
