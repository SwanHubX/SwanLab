# 导出初始化函数和log函数，以及一些数据处理模块
from .data import *
from .env import SwanLabEnv
from .package import get_package_version
from .swanlab_settings import Settings
from .sync import sync_wandb, sync_tensorboardX, sync_tensorboard_torch, sync_mlflow

# 设置默认环境变量
SwanLabEnv.set_default()
# 检查当前需要检查的环境变量
SwanLabEnv.check()

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
    "Object3D",
    "Molecule",
    "Text",
    "Run",
    "State",
    "get_run",
    "get_config",
    "config",
    "__version__",
]
