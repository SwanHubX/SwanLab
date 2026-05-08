"""
SwanLab notification callbacks.

Each callback runs IO (HTTP/SMTP) in a shared background thread pool so the
training main process is never blocked.  Exceptions are isolated by the
CallbackManager's ``safe.block`` dispatcher — individual callback failures
do not affect other callbacks or the training loop.
"""

from swanlab.plugin.notification.bark import BarkCallback
from swanlab.plugin.notification.base import NotificationCallback
from swanlab.plugin.notification.dingtalk import DingTalkCallback
from swanlab.plugin.notification.discord import DiscordCallback
from swanlab.plugin.notification.email import EmailCallback
from swanlab.plugin.notification.lark import LarkCallback
from swanlab.plugin.notification.slack import SlackCallback
from swanlab.plugin.notification.telegram import TelegramCallback
from swanlab.plugin.notification.wecom import WeComCallback, WXWorkCallback

__all__ = [
    "LarkCallback",
    "DingTalkCallback",
    "WeComCallback",
    "WXWorkCallback",
    "DiscordCallback",
    "SlackCallback",
    "TelegramCallback",
    "EmailCallback",
    "BarkCallback",
    "NotificationCallback",
]
