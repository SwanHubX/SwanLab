import base64
import hashlib
import hmac
from datetime import datetime
from typing import Any, Dict, Optional

import requests

from swanlab.plugin.notification.base import NotificationCallback
from swanlab.sdk.internal.pkg import console, safe


class _LarkBot:
    """Generate the signature required by Lark webhook verification."""

    def __init__(self, webhook_url: str, secret: Optional[str] = None):
        self.webhook_url = webhook_url
        self.secret = secret

    def gen_sign(self, timestamp: int) -> str:
        if not self.secret:
            raise ValueError("secret is required for signing")
        string_to_sign = f"{timestamp}\n{self.secret}"
        hmac_code = hmac.new(string_to_sign.encode("utf-8"), digestmod=hashlib.sha256).digest()
        return base64.b64encode(hmac_code).decode("utf-8")


class LarkCallback(NotificationCallback):
    """Lark (Feishu) notification callback.

    Usage::

        from swanlab.plugin import LarkCallback

        swanlab.init(
            callbacks=[LarkCallback(webhook_url="https://...", secret="...")]
        )
    """

    def __init__(self, webhook_url: str, secret: Optional[str] = None, language: str = "zh"):
        super().__init__(language=language)
        self._signer = _LarkBot(webhook_url, secret)

    def _send_notification(self, state: str, error: Optional[str]) -> None:
        with safe.block(requests.RequestException, message="❌ LarkBot sending failed"):
            content = self._build_content(state, error)
            timestamp = int(datetime.now().timestamp())
            payload: Dict[str, Any] = {
                "timestamp": timestamp,
                "msg_type": "text",
                "content": {"text": content},
            }
            if self._signer.secret:
                payload["sign"] = self._signer.gen_sign(timestamp)
            resp = requests.post(self._signer.webhook_url, json=payload, timeout=10)
            resp.raise_for_status()
            result: Dict[str, Any] = resp.json()
            if result.get("code") and result["code"] != 0:
                console.warning(f"❌ LarkBot sending failed: {result.get('msg')}")
                return
            console.info("✅ LarkBot notification sent successfully")
