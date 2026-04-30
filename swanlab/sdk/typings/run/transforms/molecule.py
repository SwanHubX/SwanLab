"""
@author: cunyue
@file: molecule.py
@time: 2026/4/30
@description: 分子结构处理模块类型标注
"""

from pathlib import Path
from typing import TYPE_CHECKING, List, Union

if TYPE_CHECKING:
    from swanlab.sdk.internal.run.transforms.molecule import Molecule

MoleculeDataType = Union[
    "Molecule",
    str,
    Path,
]

MoleculeDatasType = Union[
    "Molecule",
    List["Molecule"],
    str,
    List[str],
    Path,
    List[Path],
]
