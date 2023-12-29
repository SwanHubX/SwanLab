# 导出初始化函数和log函数
from .database import (
    init,
    log,
    BaseType,
)
from .utils import get_package_version


__version__ = get_package_version()
