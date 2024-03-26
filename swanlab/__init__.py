# 导出初始化函数和log函数
from .data import (
    init,
    log,
    finish,
    config,
    Audio,
    Image,
    Text,
    Run,
)

from .utils import get_package_version


__version__ = get_package_version()
