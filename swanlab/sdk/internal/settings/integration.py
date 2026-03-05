"""
@author: cunyue
@file: integration.py
@time: 2026/3/5 20:25
@description: SwanLab 集成配置，配置 SwanLab 与外部系统的集成
"""

import os
from typing import ClassVar, Type

from pydantic import BaseModel, Field


def webhook_url_factory() -> str:
    # 使用额外的 SWANLAB_WEBHOOK 环境变量，一方面是为了向下兼容（老版本是 SWANLAB_WEBHOOK），另一方面是自动生成的环境变量太长了
    return os.environ.get("SWANLAB_WEBHOOK", "")


def webhook_value_factory() -> str:
    # 使用额外的 SWANLAB_WEBHOOK_VALUE 环境变量，一方面是为了向下兼容（老版本是 SWANLAB_WEBHOOK_VALUE），另一方面是自动生成的环境变量太长了
    return os.environ.get("SWANLAB_WEBHOOK_VALUE", "")


def dashboard_host_factory() -> str:
    # 使用额外的 SWANLAB_DASHBOARD_HOST 环境变量，因为自动生成的环境变量太长了
    return os.environ.get("SWANLAB_DASHBOARD_HOST", "127.0.0.1")


def dashboard_port_factory() -> int:
    # 使用额外的 SWANLAB_DASHBOARD_PORT 环境变量，因为自动生成的环境变量太长了
    return int(os.environ.get("SWANLAB_DASHBOARD_PORT", "9090"))


class WebhookSettings(BaseModel):
    url: str = Field(default_factory=webhook_url_factory)
    """
    Webhook URL for SwanLab notifications.
    """

    value: str = Field(default_factory=webhook_value_factory)
    """
    Webhook value for SwanLab notifications.
    """


class DashBoardSettings(BaseModel):
    host: str = Field(default_factory=dashboard_host_factory)
    """
    Dashboard server host.
    """

    port: int = Field(default_factory=dashboard_port_factory)
    """
    Dashboard server port.
    """


class IntegrationSettings(BaseModel):
    """
    Configuration for SwanLab integrations.
    """

    Webhook: ClassVar[Type[WebhookSettings]] = WebhookSettings
    Dashboard: ClassVar[Type[DashBoardSettings]] = DashBoardSettings

    webhook: WebhookSettings = Field(default_factory=WebhookSettings)
    dashboard: DashBoardSettings = Field(default_factory=DashBoardSettings)
