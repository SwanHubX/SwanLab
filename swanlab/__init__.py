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


__version__ = get_package_version()
