# 导出初始化函数和log函数
import swanlab.jupyter
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


def load_ipython_extension(ipython):
    ipython.register_magics(swanlab.jupyter.SwanLabMagics)


__version__ = get_package_version()
