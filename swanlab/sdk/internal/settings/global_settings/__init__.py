"""
@author: cunyue
@description: 以全局单例的方式存储用户通过`merge_settings`方式保存的全局配置
"""

from typing import TYPE_CHECKING, Optional

if TYPE_CHECKING:
    from swanlab.sdk.internal.settings import Settings


global_instance: Optional["Settings"] = None


def get_global_settings() -> Optional["Settings"]:
    """Get the global settings instance.

    Returns:
        Optional[Settings]: The global settings instance if it exists, otherwise None.
    """
    return global_instance


def set_global_settings(settings: "Settings") -> None:
    """Set the global settings instance.

    Args:
        settings (Settings): The Settings object to set as the global instance.
    """
    global global_instance
    global_instance = settings
