"""
@author: cunyue
@file: filename.py
@time: 2026/3/15 01:20
@description: 约定的文件名
"""

metadata = "swanlab-metadata.json"
"""
存放运行时产生的元数据文件，如指标、配置、状态等
"""

config = "config.yaml"
"""
存放运行时产生的配置文件，如超参数、模型等
"""

requirements = "requirements.txt"
"""
存放Python依赖文件，如requirements.txt
"""

conda = "conda.yaml"
"""
存放Conda环境文件，如conda.yaml
"""


def run(run_id: str) -> str:
    """
    根据运行ID生成运行文件名。
    """
    return f"run-{run_id}.swanlab"
