"""
@author: cunyue
@file: merge_settings.py
@time: 2026/3/6 21:48
@description: swanlab.merge_settings 方法，合并自定义配置项
"""

from typing import Union

from swanlab.sdk.cmd.guard import with_cmd_lock, without_run
from swanlab.sdk.internal.settings import Settings
from swanlab.sdk.internal.settings import settings as swanlab_settings


@with_cmd_lock
@without_run("merge_settings")
def merge_settings(settings: Union[Settings, dict]) -> None:
    """Merge custom settings into the global SwanLab configuration.

    This function allows you to customize SwanLab's behavior before initializing a run.
    It must be called before `swanlab.init()`.

    :param settings: Custom settings to merge. Can be either a Settings object or a dict.

    Examples:

        Merge settings from a dictionary:

        >>> import swanlab
        >>> swanlab.merge_settings({"mode": "local", "logdir": "./my_logs"})
        >>> swanlab.init()

        Disable hardware monitoring:

        >>> import swanlab
        >>> swanlab.merge_settings({"hardware": {"monitor": False}})
        >>> swanlab.init()

        Use a Settings object:

        >>> import swanlab
        >>> from swanlab import Settings
        >>> custom_settings = Settings(mode="offline", logdir="./experiments")
        >>> swanlab.merge_settings(custom_settings)
        >>> swanlab.init()
    """
    swanlab_settings.merge_settings(settings)
