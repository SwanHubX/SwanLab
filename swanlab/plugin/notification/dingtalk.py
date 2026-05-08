import base64
import hashlib
import hmac
import urllib.parse
from datetime import datetime
from typing import Any, Dict, Optional

import requests

from swanlab.plugin.notification.base import NotificationCallback
from swanlab.sdk.internal.pkg import console, safe


class _DingTalkBot:
    """DingTalk webhook bot with optional HMAC-SHA256 signing.

    DingTalk signs are embedded in the URL query params and expire after
    1 hour.  ``_check_sign()`` is called before each request to refresh
    the signature when necessary.
    """

    def __init__(self, webhook_url: str, secret: Optional[str] = None):
        self.webhook_url = webhook_url
        self.secret = secret
        self._sign_ts: float = 0.0
        if self.secret and self.secret.startswith("SEC"):
            self._refresh_sign()

    def _check_sign(self) -> None:
        if self.secret and self.secret.startswith("SEC"):
            if datetime.now().timestamp() - self._sign_ts >= 3600:
                self._refresh_sign()

    def _refresh_sign(self) -> None:
        if not self.secret:
            return
        self._sign_ts = datetime.now().timestamp()
        ts_ms = round(self._sign_ts * 1000)
        string_to_sign = f"{ts_ms}\n{self.secret}"
        hmac_code = hmac.new(self.secret.encode(), string_to_sign.encode(), digestmod=hashlib.sha256).digest()
        sign = base64.b64encode(hmac_code).decode("utf-8")
        parsed = urllib.parse.urlparse(self.webhook_url)
        params = dict(urllib.parse.parse_qsl(parsed.query))
        params.update({"timestamp": str(ts_ms), "sign": sign})
        self.webhook_url = urllib.parse.urlunparse(parsed._replace(query=urllib.parse.urlencode(params)))


class DingTalkCallback(NotificationCallback):
    """DingTalk (钉钉) notification callback.

    Usage::

        from swanlab.plugin import DingTalkCallback

        swanlab.init(
            callbacks=[DingTalkCallback(webhook_url="https://...", secret="SEC...")]
        )
    """

    def __init__(self, webhook_url: str, secret: Optional[str] = None, language: str = "zh"):
        super().__init__(language=language)
        self._bot = _DingTalkBot(webhook_url, secret)

    def _send_notification(self, state: str, error: Optional[str]) -> None:
        with safe.block(requests.RequestException, message="❌ DingTalkBot sending failed"):
            self._bot._check_sign()
            content = self._build_content(state, error)
            payload: Dict[str, Any] = {"msgtype": "text", "text": {"content": content}}
            resp = requests.post(self._bot.webhook_url, json=payload, timeout=10)
            resp.raise_for_status()
            result: Dict[str, Any] = resp.json()
            if result.get("errcode") and result["errcode"] != 0:
                console.warning(f"❌ DingTalkBot sending failed: {result.get('errmsg')}")
                return
            console.info("✅ DingTalkBot notification sent successfully")
