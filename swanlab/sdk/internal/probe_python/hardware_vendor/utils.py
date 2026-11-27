"""
@author: Nexisato
@file: utils.py
@time: 2026/11/26
@description: 厂商命令行调用的统一超时约定

probe 采集属于尽力而为：外部命令可能因驱动异常、GPU 挂起或容器限制长期不返回，
因此所有调用都要显式带上超时。超时抛出的 subprocess.TimeoutExpired 与“命令不存在”
一样由调用方的 safe.block / safe.decorator 吞掉，采集结果按缺失处理。
"""

__all__ = ["DEFAULT_TIMEOUT"]

# 单次外部命令的等待上限（秒）。取值需覆盖大节点上 nvidia-smi / npu-smi 的正常耗时。
DEFAULT_TIMEOUT = 10.0
