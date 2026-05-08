from typing import Any, Dict, Optional

import requests

from swanlab.plugin.notification.base import NotificationCallback
from swanlab.sdk.internal.pkg import console, safe


class TelegramCallback(NotificationCallback):
    """Telegram notification callback.

    Usage::

        from swanlab.plugin import TelegramCallback

        swanlab.init(
            callbacks=[TelegramCallback(bot_token="123456:ABC...", chat_id="123456789")]
        )
    """

    def __init__(self, bot_token: str, chat_id: str, language: str = "zh"):
        super().__init__(language=language)
        self._bot_token = bot_token
        self._chat_id = chat_id
        self._api_url = f"https://api.telegram.org/bot{bot_token}/sendMessage"

    def _send_notification(self, state: str, error: Optional[str]) -> None:
        with safe.block(requests.RequestException, message="❌ TelegramBot sending failed"):
            content = self._build_content(state, error)
            payload: Dict[str, Any] = {"chat_id": self._chat_id, "text": content}
            resp = requests.post(self._api_url, json=payload, timeout=10)
            resp.raise_for_status()
            result: Dict[str, Any] = resp.json()
            if not result.get("ok"):
                console.warning(f"❌ TelegramBot sending failed: {result.get('description')}")
                return
            console.info("✅ TelegramBot notification sent successfully")
