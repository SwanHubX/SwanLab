"""
@author: cunyue
@file: __init__.py.py
@time: 2026/3/5 14:35
@description: SwanLab SDK，负责SwanLab库的核心指标上传功能
内部所有方法的基础是对 context 的正确初始化和维护
"""

from typing import Union

from .internal.run import SwanLabRun
from .internal.settings import Settings
from .internal.settings import settings as swanlab_settings

__all__ = ["Settings", "merge_settings", "init", "finish", "SwanLabRun"]


def merge_settings(settings: Union[Settings, dict]) -> None:
    """
    Merge custom settings into the current SwanLab settings.
    Note that this function must be called before calling `swanlab.init`.
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


def init(
    # project: str = None,
    # workspace: str = None,
    # experiment_name: str = None,
    # description: str = None,
    # job_type: str = None,
    # group: str = None,
    # tags: List[str] = None,
    # config: Union[dict, str] = None,
    # logdir: str = None,
    # mode: MODES = None,
    # load: str = None,
    # public: bool = None,
    # callbacks: List[SwanKitCallback] = None,
    # settings: Settings = None,
    # id: str = None,
    # resume: Union[ResumeType, bool] = None,
    # reinit: bool = None,
    **kwargs,
):
    """
    Start a new run to track and log. Once you have called this function, you can use 'swanlab.log' to log data to
    the current run. Meanwhile, you can use 'swanlab.finish' to finish the current run and close the current
    experiment. After calling this function, SwanLab will begin to record the console output of the current process,
    and register a callback function to the exit function.
    :param settings: The settings for the current experiment.
    :return: The SwanLabRun object.
    """
    # settings_with_env = Settings()


def finish():
    """
    Finish the current run and close the current experiment.
    """
    pass
