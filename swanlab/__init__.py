# 导出初始化函数和log函数
from .data import (
    login,
    init,
    log,
    finish,
    config,
    Audio,
    Image,
    Text,
    Run,
    State,
    get_run,
)


from .package import get_package_version


__version__ = get_package_version()
