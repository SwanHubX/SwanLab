"""
@author: cunyue
@file: object3d.py
@time: 2026/4/30
@description: 3D 对象处理模块类型标注
"""

from pathlib import Path
from typing import TYPE_CHECKING, List, Union

if TYPE_CHECKING:
    from swanlab import vendor
    from swanlab.sdk.internal.run.transforms.object3d import Object3D

Object3DDataType = Union[
    "Object3D",
    str,
    Path,
    "vendor.np.ndarray",
    dict,
]

Object3DDatasType = Union[
    "Object3D",
    List["Object3D"],
    str,
    List[str],
    Path,
    List[Path],
    "vendor.np.ndarray",
    List["vendor.np.ndarray"],
    dict,
    List[dict],
]
