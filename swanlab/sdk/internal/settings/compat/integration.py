"""集成域遗留环境变量兼容层。"""

from .env import getenv

__all__ = [
    "dashboard_host_factory",
    "dashboard_port_factory",
    "webhook_timeout_factory",
    "webhook_url_factory",
    "webhook_value_factory",
]


def webhook_url_factory() -> str:
    # 使用额外的 SWANLAB_WEBHOOK 环境变量，一方面是为了向下兼容（老版本是 SWANLAB_WEBHOOK），另一方面是自动生成的环境变量太长了
    return getenv("SWANLAB_WEBHOOK", "")


def webhook_value_factory() -> str:
    # 使用额外的 SWANLAB_WEBHOOK_VALUE 环境变量，一方面是为了向下兼容（老版本是 SWANLAB_WEBHOOK_VALUE），另一方面是自动生成的环境变量太长了
    return getenv("SWANLAB_WEBHOOK_VALUE", "")


def webhook_timeout_factory() -> int:
    # 使用额外的 SWANLAB_WEBHOOK_TIMEOUT 环境变量，因为自动生成的环境变量太长了
    return int(getenv("SWANLAB_WEBHOOK_TIMEOUT", "5"))


def dashboard_host_factory() -> str:
    # 使用额外的 SWANLAB_DASHBOARD_HOST 环境变量，因为自动生成的环境变量太长了
    return getenv("SWANLAB_DASHBOARD_HOST", "127.0.0.1")


def dashboard_port_factory() -> int:
    # 使用额外的 SWANLAB_DASHBOARD_PORT 环境变量，因为自动生成的环境变量太长了
    return int(getenv("SWANLAB_DASHBOARD_PORT", "5092"))
