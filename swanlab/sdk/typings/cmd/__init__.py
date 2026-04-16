"""
@author: cunyue
@file: __init__.py
@time: 2026/4/15 23:35
@description: SwanLab SDK cmd 模块的类型定义和导出
"""

import os
from typing import Any, Dict, Literal, Union

LoginType = Union[bool, Literal["local", "root"]]
"""
登录类型
- `False`: 不保存登录信息
- `True` 或 `"root"`: 保存到用户目录下的 `.swanlab` 目录中
- `"local"`: 保存到当前目录下的 `.swanlab` 目录中
"""

ConfigLike = Union[Dict[str, Any], str, os.PathLike]
"""
swanlab.init 允许传入的配置类型
"""
