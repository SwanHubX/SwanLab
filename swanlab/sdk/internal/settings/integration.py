"""
@author: cunyue
@file: integration.py
@time: 2026/3/5 20:25
@description: SwanLab 集成配置，配置 SwanLab 与外部系统的集成
"""

from typing import ClassVar, Type

from pydantic import BaseModel, Field, model_validator
from pydantic.config import ConfigDict

from .compat import (
    dashboard_host_factory,
    dashboard_port_factory,
    webhook_timeout_factory,
    webhook_url_factory,
    webhook_value_factory,
)


class WebhookSettings(BaseModel):
    url: str = Field(default_factory=webhook_url_factory)
    """
    Webhook URL for SwanLab notifications.
    """

    value: str = Field(default_factory=webhook_value_factory)
    """
    Webhook value for SwanLab notifications.
    """

    timeout: int = Field(default_factory=webhook_timeout_factory)
    """
    Webhook timeout for SwanLab notifications.
    """
    model_config = ConfigDict(frozen=True)


class DashBoardSettings(BaseModel):
    host: str = Field(default_factory=dashboard_host_factory)
    """
    Dashboard server host.
    """

    port: int = Field(default_factory=dashboard_port_factory, ge=1, le=65535, validate_default=True)
    """
    Dashboard server port.
    """
    model_config = ConfigDict(frozen=True)


class IntegrationSettings(BaseModel):
    """
    Configuration for SwanLab integrations.
    """

    Webhook: ClassVar[Type[WebhookSettings]] = WebhookSettings
    Dashboard: ClassVar[Type[DashBoardSettings]] = DashBoardSettings

    webhook: WebhookSettings = Field(default_factory=WebhookSettings)
    dashboard: DashBoardSettings = Field(default_factory=DashBoardSettings)

    @model_validator(mode="before")
    @classmethod
    def assemble_nested_env(cls, data: dict) -> dict:
        """
        拦截 Pydantic 因 max_split=1 截断生成的平铺环境变量，
        将其重新组装为嵌套字典，以适配内部结构。
        """
        if isinstance(data, dict):
            # 处理 webhook_xxx -> webhook: {xxx: ...}
            webhook_data = data.get("webhook", {})
            if isinstance(webhook_data, dict):
                has_update = False
                for key in list(data.keys()):
                    if key.startswith("webhook_"):
                        webhook_data[key[8:]] = data.pop(key)
                        has_update = True
                if has_update:
                    data["webhook"] = webhook_data

            # 处理 dashboard_xxx -> dashboard: {xxx: ...}
            dashboard_data = data.get("dashboard", {})
            if isinstance(dashboard_data, dict):
                has_update = False
                for key in list(data.keys()):
                    if key.startswith("dashboard_"):
                        dashboard_data[key[10:]] = data.pop(key)
                        has_update = True
                if has_update:
                    data["dashboard"] = dashboard_data

        return data

    model_config = ConfigDict(frozen=True)
