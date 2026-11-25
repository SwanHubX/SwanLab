"""
@author: mx
@file: eval.py
@time: 2026/11/25
@description: SwanLab 运行时评估API类型
"""

from typing import Any, Dict, List, Optional, TypedDict


class InitEvalType(TypedDict):
    # 评估cuid
    cuid: str
    # 评估名称
    name: str
    # 评估标签
    tags: List[str]
    # 外部链接
    link: Optional[str]
    # 评估配置，格式为 {key: {value, desc, sort}}
    config: Optional[Dict[str, Any]]
    # ckpt 产生的训练步号，服务端序列化为十进制字符串
    step: Optional[str]
    # 关联实验的 cuid
    experimentCuid: Optional[str]
    # 创建时间
    createdAt: str
    # 更新时间
    updatedAt: str
