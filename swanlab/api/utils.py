"""
@author: Zhou QiYang
@file: utils.py
@time: 2026/1/11 23:44
@description: OpenApi 中的基础对象与通用工具
"""

from dataclasses import dataclass
from functools import wraps

from swanlab.core_python.api.type import IdentityType
from swanlab.core_python.api.user import get_self_hosted_init
from swanlab.core_python.client import Client
from swanlab.error import ApiError
from swanlab.log import swanlog


@dataclass
class Label:
    """
    Project label object
    you can get the label name by str(label)
    """

    name: str

    def __str__(self) -> str:
        return self.name


def self_hosted(identity: IdentityType = "user"):
    """
    用于需要在私有化环境下使用的接口的装饰器。
    :param identity: 用户身份，默认为 "user"，如果为 "root"，则会额外验证是否为根用户。
    """

    def decorator(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            client = getattr(self, '_client', None)
            if not isinstance(client, Client):
                swanlog.warning("There is no SwanLab client instance.")
                return None

            # 1. 尝试获取私有化服务信息
            try:
                self_hosted_info = get_self_hosted_init(client)
            except ApiError:
                raise ValueError("You haven't launched a swanlab self-hosted instance. This usages are not available.")

            if not self_hosted_info.get("enabled", False):
                swanlog.warning("SwanLab self-hosted instance hasn't been ready yet.")
                return None

            if self_hosted_info.get("expired", True):
                swanlog.warning("SwanLab self-hosted instance has expired.")
                return None

            # 2. 检测用户权限（商业版root用户功能）
            if identity == "root" and not self_hosted_info.get("root", False):
                swanlog.warning("You don't have permission to perform this action. Please login as a root user")
                return None

            return func(self, *args, **kwargs)

        return wrapper

    return decorator


__all__ = ['Label', 'self_hosted']
