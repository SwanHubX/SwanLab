from typing import Any, Dict, Optional

import requests

from swanlab.plugin.notification.base import NotificationCallback
from swanlab.sdk.internal.pkg import console, safe


class WeComCallback(NotificationCallback):
    """WeCom (企业微信) notification callback.

    Usage::

        from swanlab.plugin import WeComCallback

        swanlab.init(
            callbacks=[WeComCallback(webhook_url="https://...")]
        )
    """

    def __init__(self, webhook_url: str, language: str = "zh"):
        super().__init__(language=language)
        self._webhook_url = webhook_url

    def _send_notification(self, state: str, error: Optional[str]) -> None:
        with safe.block(requests.RequestException, message="❌ WeComBot sending failed"):
            content = self._build_content(state, error)
            payload: Dict[str, Any] = {"msgtype": "text", "text": {"content": content}}
            resp = requests.post(self._webhook_url, json=payload, timeout=10)
            resp.raise_for_status()
            result: Dict[str, Any] = resp.json()
            if result.get("errcode") and result["errcode"] != 0:
                console.warning(f"❌ WeComBot sending failed: {result.get('errmsg')}")
                return
            console.info("✅ WeComBot notification sent successfully")


# 前向兼容 WXWork 命名
WXWorkCallback = WeComCallback
