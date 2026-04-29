"""
@author: cunyue
@file: env.py
@time: 2026/3/8 17:57
@description: SwanLab SDK 环境相关工具
"""

import io
import os
import sys

__all__ = ["is_jupyter", "is_interactive", "is_jupyter_supports_stdin", "DEBUG"]

# 设计上，DEBUG 变量独立于其他模块，包括 settings。但是我们又希望有一个环境变量去控制是否打印 debug 信息，所以这里额外绑定一个环境变量
DEBUG = os.getenv("SWANLAB_DEBUG", "false").lower() in ["true", "1", "yes", "on"]
"""
是否为调试模式
"""


def is_jupyter() -> bool:
    """判断当前是否在 Jupyter Notebook 环境中"""
    try:
        from IPython.core.getipython import get_ipython

        # get_ipython() 如果在普通 python 终端会返回 TerminalInteractiveShell
        # 如果在 jupyter 会返回 ZMQInteractiveShell
        shell = get_ipython().__class__.__name__
        if shell == "ZMQInteractiveShell":
            return True
        elif shell == "TerminalInteractiveShell":
            return False
        else:
            return False
    except NameError:
        return False
    except ImportError:
        return False


def is_jupyter_supports_stdin() -> bool:
    """检查 Jupyter kernel 是否支持 stdin 输入（如 nbconvert --execute 不支持）"""
    # noinspection PyBroadException
    try:
        from IPython.core.getipython import get_ipython

        ip = get_ipython()
        kernel = getattr(ip, "kernel", None)
        if kernel is not None and hasattr(kernel, "_allow_stdin"):
            return getattr(kernel, "_allow_stdin")
    except Exception:
        pass
    # 无法确定时保守返回 True（交互式 Jupyter 更常见）
    return True


def is_interactive() -> bool:
    """
    是否为可交互式环境（输入连接tty设备）
    特殊的环境：jupyter notebook（但排除 nbconvert 等不支持 stdin 的场景）
    """
    try:
        fd = sys.stdin.fileno()
        if os.isatty(fd):
            return True
        if is_jupyter():
            return is_jupyter_supports_stdin()
        return False
    # 当使用 capsys、capfd 或 monkeypatch 等 fixture 来捕获或修改标准输入输出时
    # 会抛出 io.UnsupportedOperation，多为测试情况，视为可交互
    except io.UnsupportedOperation:
        return True
    except AttributeError:
        # 有些特殊的运行环境连 sys.stdin 都没有
        return False
