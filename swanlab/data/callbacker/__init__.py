"""
@author: cunyue
@file: __init__.py
@time: 2025/6/21 17:10
@description: 回调器模块，分为三大回调器：
1. local： 本地模式回调器
2. cloud： 云端模式回调器
3. offline： 离线模式回调器
"""

from .cloud import CloudPyCallback
from .local import LocalRunCallback
from .offline import OfflineCallback
