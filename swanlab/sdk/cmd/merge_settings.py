"""
@author: cunyue
@file: merge_settings.py
@time: 2026/3/6 21:48
@description: swanlab.merge_settings 方法，合并自定义配置项
"""

from typing import Union

from swanlab.sdk.cmd.helper import with_cmd_lock
from swanlab.sdk.internal.settings import Settings
from swanlab.sdk.internal.settings import settings as swanlab_settings


@with_cmd_lock
def merge_settings(settings: Union[Settings, dict]) -> None:
    """
    Merge custom settings into the current SwanLab settings.
    Note that this function must be called before calling `swanlab.init`.

    Examples:
    ---------
    >>> import swanlab
    >>> # 1. Merge from dict
    >>> swanlab.merge_settings({"hardware": {"monitor": False}})
    >>> # 2. Merge from SwanLabSettings instance
    >>> swanlab.merge_settings(Settings(hardware=Settings.Hardware(monitor=False)))
    >>>
    >>> # Now you can call swanlab.init(), which will use the merged settings
    >>> swanlab.init()

    :param settings: The custom settings to merge.
    """
    swanlab_settings.merge_settings(settings)
