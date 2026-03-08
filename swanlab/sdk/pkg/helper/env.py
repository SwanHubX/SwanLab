"""
@author: cunyue
@file: env.py
@time: 2026/3/8 17:57
@description: SwanLab SDK 环境相关工具
"""

import io
import os
import sys

__all__ = ["is_jupyter", "is_interactive"]


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


def is_interactive() -> bool:
    """
    是否为可交互式环境（输入连接tty设备）
    特殊的环境：jupyter notebook
    """
    try:
        fd = sys.stdin.fileno()
        return os.isatty(fd) or is_jupyter()
    # 当使用 capsys、capfd 或 monkeypatch 等 fixture 来捕获或修改标准输入输出时
    # 会抛出 io.UnsupportedOperation，多为测试情况，视为可交互
    except io.UnsupportedOperation:
        return True
    except AttributeError:
        # 有些特殊的运行环境连 sys.stdin 都没有
        return False
