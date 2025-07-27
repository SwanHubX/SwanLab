# 导出初始化函数和log函数，以及一些数据处理模块
from .data import *
from .data.modules import (
    Audio,
    Image,
    Object3D,
    Molecule,
    Text,
    Video,
    echarts,
    roc_curve,
    pr_curve,
    confusion_matrix,
)
from .env import SwanLabEnv
from .package import get_package_version
from .swanlab_settings import Settings
from .sync import sync_wandb, sync_tensorboardX, sync_tensorboard_torch, sync_mlflow, sync

# 设置默认环境变量
SwanLabEnv.set_default()
# 检查当前需要检查的环境变量
SwanLabEnv.check()

# 导出 OpenApi 接口，必须要等待上述的 import 语句执行完毕以后才能导出，否则会触发循环引用
from .api import OpenApi

__version__ = get_package_version()

__all__ = [
    "login",
    "Settings",
    "merge_settings",
    "init",
    "log",
    "register_callbacks",
    "finish",
    "Audio",
    "Image",
    "echarts",
    "Object3D",
    "Molecule",
    "Text",
    "Video",
    "Run",
    "State",
    "get_run",
    "get_config",
    "config",
    "OpenApi",
    "sync_wandb",
    "sync_mlflow",
    "sync_tensorboardX",
    "sync_tensorboard_torch",
    "sync",
    "__version__",
]
