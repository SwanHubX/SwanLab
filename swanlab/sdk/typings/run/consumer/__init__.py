"""
@author: caddiesnew
@file: __init__.py
@time: 2026-11-03 17:26:33
@description: 后台消费者组件相关类型标注
"""

from typing import List, Optional, Union

# 单次 RPC 发送的 record 默认批次
DEFAULT_CONSUMER_RECORD_BATCH: int = 100
