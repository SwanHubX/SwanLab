"""
@author: cunyue
@file: __init__.py.py
@time: 2026/3/31 00:23
@description: 系统环境信息采集模块
"""

# from swanlab.sdk.internal.settings import Settings
# from swanlab.sdk.typings.run.system import SystemEnvironment
#
# from . import conda, git, requirements, runtime
#
#
# def generate(settings: Settings) -> SystemEnvironment:
#     # 1. 获取一些常规环境信息
#     git_snapshot = git.get() if settings.environment.git else None
#     runtime_snapshot = runtime.get() if settings.environment.runtime else None
#     conda_snapshot = conda.get() if settings.environment.conda else None
#     requirements_snapshot = requirements.get() if settings.environment.requirements else None
#     # 2. 获取硬件信息
#     # 只有当选择不采集硬件信息也不监控硬件信息时，才不采集硬件信息
#
#     return SystemEnvironment()
