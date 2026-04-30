"""
@author: cunyue
@file: __init__.py
@time: 2026/4/30
@description: 分子结构处理模块，基于 RDKit 支持多种分子格式
"""

import functools
import hashlib
from io import BytesIO
from pathlib import Path
from typing import Optional

from swanlab import vendor
from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaItem
from swanlab.sdk.internal.context import TransformMedia
from swanlab.sdk.internal.pkg import fs
from swanlab.sdk.typings.run.transforms import CaptionType
from swanlab.sdk.typings.run.transforms.molecule import MoleculeDataType

# ---------- 文件扩展名 -> 加载方式映射 ----------

_FILE_EXT_MAP = {
    ".pdb": "pdb",
    ".sdf": "sdf",
    ".sd": "sdf",
    ".mol": "mol",
}


def _require_rdkit(func):
    """装饰器：检查 RDKit 是否可用，不可用时抛出 ImportError。"""

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        if vendor.rdkit is None:
            raise ImportError(
                "RDKit is required for Molecule. "
                "To enable it, please install the 'media' extra by running:\n"
                '    pip install "swanlab[media]"'
            )
        return func(*args, **kwargs)

    return wrapper


def _ensure_3d(mol: "vendor.rdkit.Chem.Mol"):
    if mol.GetNumConformers() == 0:
        from rdkit.Chem import AllChem

        AllChem.EmbedMolecule(mol)  # type: ignore[attr-defined]


class Molecule(TransformMedia):
    """分子结构转换类，依赖 RDKit。

    支持的输入类型:
    - str: SMILES 字符串 或 文件路径（.pdb/.sdf/.mol）
    - Path: 文件路径（.pdb/.sdf/.mol）
    - rdkit.Chem.Mol: RDKit 分子对象
    - Molecule: 套娃加载

    内部统一转换为 PDB 格式存储。
    """

    def __init__(self, data: MoleculeDataType, caption: CaptionType = None):
        super().__init__()

        # 套娃加载
        attrs = self._unwrap(data)
        if attrs:
            self.buffer: BytesIO = attrs["buffer"]
            self.caption: Optional[str] = caption if caption is not None else attrs.get("caption")
            return

        self.caption = caption
        pdb_data = self._resolve_input(data)
        self.buffer = BytesIO(pdb_data.encode())

    def _resolve_input(self, data) -> str:
        """根据输入类型解析为 PDB 字符串。"""
        Chem = vendor.rdkit.Chem if vendor.rdkit is not None else None
        Mol = Chem.Mol if Chem is not None else None

        # 1. RDKit Mol 对象
        if Mol is not None and isinstance(data, Mol):
            return self._mol_to_pdb(data)

        # 2. str/Path -> 先判断是否为文件路径，再尝试 SMILES
        if isinstance(data, (str, Path)):
            path = Path(data)
            if path.exists() and path.is_file():
                return self._handle_file(path)
            # 不是文件路径，尝试 SMILES（仅限 str）
            if isinstance(data, str):
                return self._handle_smiles(data)

        raise TypeError(
            f"Unsupported input type: {type(data)}. Expected rdkit.Chem.Mol, str (SMILES or file path), or Path."
        )

    @staticmethod
    @_require_rdkit
    def _mol_to_pdb(mol) -> str:
        """将 RDKit Mol 对象转换为 PDB 字符串。"""
        _ensure_3d(mol)
        return vendor.rdkit.Chem.MolToPDBBlock(mol)

    @staticmethod
    @_require_rdkit
    def _handle_smiles(smiles: str) -> str:
        """从 SMILES 字符串创建 Molecule。"""
        Chem = vendor.rdkit.Chem
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Could not parse SMILES string: {smiles}")
        _ensure_3d(mol)
        return Chem.MolToPDBBlock(mol)

    @staticmethod
    @_require_rdkit
    def _handle_file(path: Path) -> str:
        """根据文件扩展名加载分子文件。"""
        suffixes = path.suffixes
        if not suffixes:
            raise ValueError(f"File has no extension: {path}")

        all_tried = ["".join(suffixes[i:]).lower() for i in range(len(suffixes))]

        for suffix in all_tried:
            load_type = _FILE_EXT_MAP.get(suffix)
            if load_type is not None:
                if load_type == "pdb":
                    with open(path) as f:
                        return f.read()
                elif load_type == "sdf":
                    Chem = vendor.rdkit.Chem
                    suppl = Chem.SDMolSupplier(str(path))
                    mol = next(suppl)
                    if mol is None:
                        raise ValueError(f"Could not read molecule from SDF file: {path}")
                    _ensure_3d(mol)
                    return Chem.MolToPDBBlock(mol)
                elif load_type == "mol":
                    Chem = vendor.rdkit.Chem
                    mol = Chem.MolFromMolFile(str(path))
                    if mol is None:
                        raise ValueError(f"Could not read molecule from MOL file: {path}")
                    _ensure_3d(mol)
                    return Chem.MolToPDBBlock(mol)

        raise ValueError(
            f"Unsupported file type: {path.name}. Supported extensions: {', '.join(sorted(_FILE_EXT_MAP.keys()))}"
        )

    # ---------- 工厂方法 ----------

    @classmethod
    def from_mol(cls, mol, caption: CaptionType = None) -> "Molecule":
        """从 RDKit Mol 对象创建 Molecule。

        Args:
            mol: RDKit Chem.Mol 对象
            caption: 可选的说明文字
        """
        return cls(mol, caption=caption)

    @classmethod
    def from_smiles(cls, smiles: str, caption: CaptionType = None) -> "Molecule":
        """从 SMILES 字符串创建 Molecule。

        Args:
            smiles: SMILES 字符串，如 "CCO"
            caption: 可选的说明文字
        """
        return cls(smiles, caption=caption)

    @classmethod
    def from_pdb(cls, path, caption: CaptionType = None) -> "Molecule":
        """从 PDB 文件创建 Molecule。

        Args:
            path: PDB 文件路径
            caption: 可选的说明文字
        """
        return cls(path, caption=caption)

    @classmethod
    def from_sdf(cls, path, caption: CaptionType = None) -> "Molecule":
        """从 SDF 文件创建 Molecule。

        Args:
            path: SDF 文件路径
            caption: 可选的说明文字
        """
        return cls(path, caption=caption)

    @classmethod
    def from_mol_file(cls, path, caption: CaptionType = None) -> "Molecule":
        """从 MOL 文件创建 Molecule。

        Args:
            path: MOL 文件路径
            caption: 可选的说明文字
        """
        return cls(path, caption=caption)

    # ---------- 抽象方法实现 ----------

    @classmethod
    def column_type(cls) -> ColumnType:
        return ColumnType.COLUMN_TYPE_MOLECULE

    def transform(self, *, step: int, path: Path) -> MediaItem:
        content = self.buffer.getvalue()
        sha256 = hashlib.sha256(content).hexdigest()
        filename = f"{step:03d}-{sha256[:8]}.pdb"
        fs.safe_write(path / filename, content, mode="wb")
        return MediaItem(filename=filename, sha256=sha256, size=len(content), caption=self.caption or "")
