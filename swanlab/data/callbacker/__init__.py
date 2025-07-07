"""
@author: cunyue
@file: __init__.py
@time: 2025/6/21 17:10
@description: 回调器模块，四大回调器对应 swanlab 的四种运行模式。
1. local： 本地模式回调器，此回调不自动导出
2. cloud： 云端模式回调器
3. offline： 离线模式回调器
4. disabled： 禁用回调器
"""

from .cloud import CloudPyCallback
from .disabled import DisabledCallback
from .offline import OfflineCallback
