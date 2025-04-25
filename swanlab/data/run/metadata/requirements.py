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
        # 运行uv命令获取当前环境下的环境目录
        result = subprocess.run(["uv", "pip", "list", "--format=freeze"], capture_output=True, text=True, timeout=15)
        if result.returncode != 0:
            raise FileNotFoundError("uv command run failed")
        return result.stdout
    except FileNotFoundError:
        # 运行pip命令获取当前环境下的环境目录
        result = subprocess.run(["pip", "list", "--format=freeze"], stdout=subprocess.PIPE, text=True)
        # 检查命令是否成功运行
        if result.returncode == 0:
            return result.stdout
        else:
            return None
