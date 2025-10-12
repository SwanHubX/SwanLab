"""
@author: cunyue
@file: notification_bark.py
@time: 2025/10/12 13:37
@description: 使用bark通知实验运行情况
"""

import argparse
import random
import time

import swanlab
from swanlab.plugin import BarkCallback

parser = argparse.ArgumentParser(description="测试Bark通知功能")
parser.add_argument('--api-key', type=str, default=None, help="SwanLab的API Key")
parser.add_argument('--host', type=str, default=None, help="SwanLab的服务器地址")
parser.add_argument("--key", type=str, required=True, help="Bark的Key")
parser.add_argument("--url", type=str, default="https://api.day.app/", help="Bark的服务器地址")
args = parser.parse_args()

if args.api_key:
    swanlab.login(api_key=args.api_key, host=args.host)
# 集成BarkCallback
swanlab.register_callbacks([BarkCallback(key=args.key, url=args.url)])

epochs = 50
lr = 0.01
offset = random.random() / 5
# 初始化
swanlab.init(description="测试Bark通知功能", mode="cloud")
swanlab.config.epochs = epochs
swanlab.config.learning_rate = lr
# 模拟训练
for epoch in range(2, swanlab.config.epochs):
    acc = 1 - 2**-epoch - random.random() / epoch - offset
    loss = 2**-epoch + random.random() / epoch + offset
    print(f"epoch={epoch}, accuracy={acc}, loss={loss}")
    swanlab.log({"t/accuracy": acc, "loss": loss})
    time.sleep(0.5)
