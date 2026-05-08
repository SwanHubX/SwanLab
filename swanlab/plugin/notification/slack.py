from typing import Any, Dict, Optional

import requests

from swanlab.plugin.notification.base import NotificationCallback
from swanlab.sdk.internal.pkg import console, safe


class SlackCallback(NotificationCallback):
    """Slack notification callback.

    Usage::

        from swanlab.plugin import SlackCallback

        swanlab.init(
            callbacks=[SlackCallback(webhook_url="https://hooks.slack.com/services/...")]
        )
    """

    def __init__(self, webhook_url: str, language: str = "zh"):
        super().__init__(language=language)
        self._webhook_url = webhook_url

    def _send_notification(self, state: str, error: Optional[str]) -> None:
        with safe.block(requests.RequestException, message="❌ SlackBot sending failed"):
            content = self._build_content(state, error)
            payload: Dict[str, Any] = {"text": content}
            resp = requests.post(self._webhook_url, json=payload, timeout=10)
            resp.raise_for_status()
            if resp.status_code not in (200, 204):
                console.warning(f"❌ SlackBot sending failed: {resp.text}")
                return
            console.info("✅ SlackBot notification sent successfully")
