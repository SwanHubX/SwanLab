"""
@author: caddiesnew
@file: __init__.py
@time: 2026/10/07
@description: probe_python 工具函数
"""

from pathlib import Path
from typing import Optional

from swanlab.proto.swanlab.save.v1.save_pb2 import SaveRecord, SaveType
from swanlab.sdk.internal.pkg import fs


def make_save_record(name: str, save_type: SaveType, content: str, target: Optional[Path]) -> SaveRecord:
    """target 为 None 时构建内存记录，否则写入文件并引用路径。"""
    if target is None:
        return SaveRecord(name=name, type=save_type, payload=content.encode("utf-8"))
    fs.safe_write(target, content)
    return SaveRecord(name=name, type=save_type, source_path=target.absolute().as_posix())
