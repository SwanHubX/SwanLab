import json
from typing import Any, Dict, Literal, Optional

import requests

from swanlab.plugin.notification.base import NotificationCallback
from swanlab.sdk.internal.pkg import console, safe


class BarkCallback(NotificationCallback):
    """Bark push notification callback for iOS devices.

    Usage::

        from swanlab.plugin import BarkCallback

        swanlab.init(
            callbacks=[BarkCallback(url="https://api.day.app")]
        )
    """

    def __init__(
        self,
        url: str = "https://api.day.app",
        title: str = "SwanLab",
        bark_level: Literal["critical", "active", "timeSensitive", "passive"] = "active",
        icon: str = "https://swanlab.cn/icon.png",
        group: Literal["exp_name", "project", "workspace", None] = None,
        click_jump: bool = True,
        language: str = "zh",
        key: Optional[str] = None,
    ):
        super().__init__(language=language)
        self._url = url.rstrip("/")
        self._title = title
        self._bark_level = bark_level
        self._icon = icon
        self._group = group
        self._click_jump = click_jump
        self._key = key if key else ""

    def _build_payload(self, error: Optional[str]) -> Dict[str, Any]:
        tpl = {
            "en": {
                "subtitle_success": "✅ Your experiment completed successfully",
                "subtitle_error": "❌ Your experiment encountered an error: {error}",
                "link_text": "Project: {project}\nWorkspace: {workspace}\nName: {exp_name}\nDescription: {description}",
                "offline_text": "Offline use of the project.",
            },
            "zh": {
                "subtitle_success": "✅ 您的实验已成功完成",
                "subtitle_error": "❌ 您的实验遇到错误: {error}",
                "link_text": "项目: {project}\n工作区: {workspace}\n实验名: {exp_name}\n描述: {description}",
                "offline_text": "项目离线使用。",
            },
        }
        t = tpl.get(self.language, tpl["zh"])

        if error:
            subtitle = t["subtitle_error"].format(error=error)
        else:
            subtitle = t["subtitle_success"]

        exp_url = self._get_url()
        if exp_url and self._settings:
            body = t["link_text"].format(
                project=self._settings.project.name or "",
                workspace=self._settings.project.workspace or "",
                exp_name=self._settings.experiment.name or "",
                description=self._settings.experiment.description or "",
            )
            if self._group == "exp_name":
                group = self._settings.experiment.name
            elif self._group == "project":
                group = self._settings.project.name
            elif self._group == "workspace":
                group = self._settings.project.workspace
            else:
                group = None
        else:
            body = t["offline_text"]
            group = None

        data: Dict[str, Any] = {
            "title": self._title,
            "subtitle": subtitle,
            "body": body,
            "level": self._bark_level,
            "icon": self._icon,
            "group": group,
            "device_key": self._key,
        }
        if exp_url and self._click_jump:
            data["url"] = exp_url
        return {k: v for k, v in data.items() if v is not None}

    def _send_notification(self, state: str, error: Optional[str]) -> None:  # noqa: ARG002
        with safe.block(requests.RequestException, message="❌ BarkBot sending failed"):
            payload = self._build_payload(error)
            resp = requests.post(
                url=f"{self._url}/push",
                headers={"Content-Type": "application/json; charset=utf-8"},
                data=json.dumps(payload, ensure_ascii=False).encode("utf-8"),
                timeout=10,
            )
            resp.raise_for_status()
            result = resp.json()
            if result.get("errcode") and result["errcode"] != 0:
                console.warning(f"❌ BarkBot sending failed: {result.get('errmsg')}")
                return
            console.info("✅ BarkBot notification sent successfully")
