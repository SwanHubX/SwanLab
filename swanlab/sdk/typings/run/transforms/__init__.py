"""
@author: cunyue
@file: __init__.py
@time: 2026/3/17 19:35
@description: 数据处理模块类型标注，主要负责标注每个模块的init类型
"""

from typing import List, Optional, Union

CaptionType = Optional[str]

CaptionsType = Optional[Union[str, List[str]]]
