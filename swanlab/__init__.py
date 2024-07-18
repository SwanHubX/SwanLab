# 导出初始化函数和log函数
from .data import (
    login,
    init,
    log,
    finish,
    Audio,
    Image,
    Text,
    Run,
    State,
    get_run,
    get_config,
)

from .data.run.main import config
from .package import get_package_version
from .env import SwanLabEnv

# 设置默认环境变量
SwanLabEnv.set_default()
# 检查当前需要检查的环境变量
SwanLabEnv.check()

__version__ = get_package_version()

__all__ = [
    "login",
    "init",
    "log",
    "finish",
    "Audio",
    "Image",
    "Text",
    "Run",
    "State",
    "get_run",
    "get_config",
    "config",
    "__version__",
]
