from dataclasses import dataclass
from typing import Dict, Optional, Tuple

from swankit.core.data import DataSuite as D
from swankit.core.data import MediaBuffer, MediaType

# 尝试导入 RDKit，如果 RDKit 不可用，则将其设置为 None
try:
    from rdkit import Chem
    from rdkit.Chem import Mol
    _has_rdkit = True
except ImportError:
    Chem = None
    Mol = None
    _has_rdkit = False


@dataclass()
class Molecule(MediaType):
    pdb_data: str
    caption: Optional[str] = None

    def __post_init__(self):
        """初始化后验证输入数据。"""
        if not isinstance(self.pdb_data, str):
            raise TypeError("pdb_data must be a string use RDKit.Chem.MolToPDBBlock.")

    @classmethod
    def from_mol(cls, mol: Mol, *, caption: Optional[str] = None, **kwargs) -> "Molecule":
        """从 RDKit Mol 对象创建 Molecule 实例。

        Args:
            mol: RDKit Mol 对象。
            caption: 可选的描述文本。

        Returns:
            Molecule: 一个新的 Molecule 实例。

        Raises:
            ValueError: 如果 RDKit 不可用。

        Examples:
            >>> from rdkit import Chem
            >>> mol = Chem.MolFromSmiles("CCO")
            >>> molecule = Molecule.from_mol(mol, caption="Ethanol")
        """
        if Chem is None or Mol is None:
            raise ValueError("RDKit is not available.")
        pdb_data = Chem.MolToPDBBlock(mol)
        return cls(pdb_data, caption=caption, **kwargs)

    # ---------------------------------- override ----------------------------------

    def parse(self) -> Tuple[str, MediaBuffer]:
        """Convert Molecule PDB to buffer for transmission.

        Returns:
            Tuple containing:
            - File name with format: molecule-step{step}-{hash}.pdb
            - MediaBuffer containing the molecule pdb data
        """

        data = self.pdb_data.encode()

        buffer = MediaBuffer()
        buffer.write(data)

        hash_name = D.get_hash_by_bytes(data)[:16]
        save_name = f"molecule-step{self.step}-{hash_name}.pdb"

        return save_name, buffer

    def get_chart(self) -> MediaType.Chart:
        """Return chart type for visualization"""
        return MediaType.Chart.OBJECT3D

    def get_section(self) -> str:
        """Return section name for organization"""
        return "Molecule"

    def get_more(self) -> Optional[Dict[str, str]]:
        """Return additional information (caption)"""
        return {"caption": self.caption} if self.caption else None
