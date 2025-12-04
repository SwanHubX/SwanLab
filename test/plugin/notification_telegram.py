"""
@author: yugangcao
@file: notification_telegram.py
@time: 2025/12/05
@description: 使用Telegram通知实验运行情况
"""

import argparse
import random
import time

import swanlab
from swanlab.plugin import TelegramCallback

parser = argparse.ArgumentParser(description="测试Telegram通知功能")
parser.add_argument('--api-key', type=str, default=None, help="SwanLab的API Key")
parser.add_argument('--host', type=str, default=None, help="SwanLab的服务器地址")
parser.add_argument("--bot-token", type=str, required=True, help="Telegram Bot Token (从@BotFather获取)")
parser.add_argument("--chat-id", type=str, required=True, help="Telegram Chat ID (从@userinfobot获取)")
parser.add_argument("--language", type=str, default="zh", choices=["zh", "en"], help="通知语言")
parser.add_argument("--notify-on-start", action="store_true", help="实验开始时也发送通知")
args = parser.parse_args()

if args.api_key:
    swanlab.login(api_key=args.api_key, host=args.host)
# 集成TelegramCallback
swanlab.register_callbacks([TelegramCallback(
    bot_token=args.bot_token,
    chat_id=args.chat_id,
    language=args.language,
    notify_on_start=args.notify_on_start,
)])

epochs = 50
lr = 0.01
offset = random.random() / 5
# 初始化
swanlab.init(description="测试Telegram通知功能", mode="cloud")
swanlab.config.epochs = epochs
swanlab.config.learning_rate = lr
# 模拟训练
for epoch in range(2, swanlab.config.epochs):
    acc = 1 - 2**-epoch - random.random() / epoch - offset
    loss = 2**-epoch + random.random() / epoch + offset
    print(f"epoch={epoch}, accuracy={acc}, loss={loss}")
    swanlab.log({"t/accuracy": acc, "loss": loss})
    time.sleep(0.5)
