"""
@author: cunyue
@file: requirements.py
@time: 2024/11/18 15:48
@description: 依赖包信息采集
"""

import subprocess


def get_requirements():
    """获取当前环境依赖"""
    try:
        # 运行pixi命令获取当前环境下的环境目录
        result = subprocess.run(["pixi", "list"], capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            return result.stdout

        # 运行uv命令获取当前环境下的环境目录
        result = subprocess.run(["uv", "pip", "list", "--format=freeze"], capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            return result.stdout
    except Exception:  # noqa: 捕获所有异常，避免因环境问题导致程序崩溃
        pass

    try:
        result = subprocess.run(["pip", "list", "--format=freeze"], capture_output=True, text=True, timeout=15)
        if result.returncode == 0:
            return result.stdout
    except Exception:  # noqa: 捕获所有异常，避免因环境问题导致程序崩溃
        pass

    return None
