from typing import Any, Dict, Optional

import requests

from swanlab.plugin.notification.base import NotificationCallback
from swanlab.sdk.internal.pkg import console, safe


class DiscordCallback(NotificationCallback):
    """Discord notification callback.

    Usage::

        from swanlab.plugin import DiscordCallback

        swanlab.init(
            callbacks=[DiscordCallback(webhook_url="https://discord.com/api/webhooks/...")]
        )
    """

    def __init__(self, webhook_url: str, language: str = "zh"):
        super().__init__(language=language)
        self._webhook_url = webhook_url

    def _send_notification(self, state: str, error: Optional[str]) -> None:
        with safe.block(requests.RequestException, message="❌ DiscordBot sending failed"):
            content = self._build_content(state, error)
            payload: Dict[str, Any] = {"content": content}
            resp = requests.post(self._webhook_url, json=payload, timeout=10)
            resp.raise_for_status()
            if resp.status_code not in (200, 204):
                console.warning(f"❌ DiscordBot sending failed: {resp.text}")
                return
            console.info("✅ DiscordBot notification sent successfully")
