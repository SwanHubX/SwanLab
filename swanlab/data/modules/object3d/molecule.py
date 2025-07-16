from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

from swanlab.toolkit import DataSuite as D, MediaBuffer, MediaType

# Attempt to import RDKit; if unavailable, set Chem and Mol to None
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Mol

    _has_rdkit = True
except ImportError:
    Chem = None
    Mol = None
    AllChem = None
    _has_rdkit = False


@dataclass()
class Molecule(MediaType):
    pdb_data: str
    caption: Optional[str] = None

    def __post_init__(self):
        """Validates input data after initialization."""
        if not isinstance(self.pdb_data, str):
            raise TypeError("pdb_data must be a string, use RDKit.Chem.MolToPDBBlock to convert.")

    @staticmethod
    def check_is_available():
        """Check if RDKit is available."""
        if not _has_rdkit:
            raise ImportError("RDKit is not available. You can install it by running 'pip install rdkit'.")

    @classmethod
    def from_mol(cls, mol: Mol, *, caption: Optional[str] = None, **kwargs) -> "Molecule":
        """Creates a Molecule instance from an RDKit Mol object.

        Args:
            mol: The RDKit Mol object.
            caption: Optional descriptive text.

        Returns:
            Molecule: A new Molecule instance.

        Raises:
            ImportError: If RDKit is not available.

        Examples:
            >>> from rdkit import Chem
            >>> mol = Chem.MolFromSmiles("CCO")
            >>> molecule = Molecule.from_mol(mol, caption="Ethanol")
        """
        cls.check_is_available()
        pdb_block = cls._convert_to_pdb_block(mol)
        return cls(pdb_block, caption=caption, **kwargs)

    @staticmethod
    def _convert_to_pdb_block(mol: Mol) -> str:
        """将分子转换为 PDB 字符串"""
        Molecule.check_is_available()
        if mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol)  # 生成 3D 坐标
        return Chem.MolToPDBBlock(mol)

    @classmethod
    def from_pdb_file(cls, pdb_file: Union[Path, str], *, caption: Optional[str] = None, **kwargs) -> "Molecule":
        """Creates a Molecule instance from a PDB file by reading the file content directly.

        Args:
            pdb_file: Path to the PDB file.
            caption: Optional descriptive text.

        Returns:
            Molecule: A new Molecule instance.

        Raises:
            ValueError: If RDKit is not available or the file cannot be read.
        """
        cls.check_is_available()
        try:
            with open(pdb_file) as f:
                pdb_data = f.read()
        except FileNotFoundError as e:
            raise FileNotFoundError(f"PDB file not found: {pdb_file}") from e
        except Exception as e:
            raise ValueError(f"Could not read PDB file: {pdb_file}. Error: {e}") from e

        # Directly create the Molecule instance with the pdb_data, skipping Mol object
        return cls(pdb_data, caption=caption, **kwargs)

    @classmethod
    def from_sdf_file(cls, sdf_file: Path, *, caption: Optional[str] = None, **kwargs) -> "Molecule":
        """Creates a Molecule instance from an SDF file.

        Args:
            sdf_file: Path to the SDF file.
            caption: Optional descriptive text.

        Returns:
            Molecule: A new Molecule instance.

        Raises:
            ImportError: If RDKit is not available.
        """
        cls.check_is_available()
        suppl = Chem.SDMolSupplier(str(sdf_file))
        mol = next(suppl)  # Assuming only one molecule in the SDF file, you can iterate if needed.
        if mol is None:
            raise ValueError(f"Could not read molecule from SDF file: {sdf_file}")
        return cls.from_mol(mol, caption=caption, **kwargs)

    @classmethod
    def from_smiles(cls, smiles: str, *, caption: Optional[str] = None, **kwargs) -> "Molecule":
        """Creates a Molecule instance from a SMILES string.

        Args:
            smiles: The SMILES string.
            caption: Optional descriptive text.

        Returns:
            Molecule: A new Molecule instance.

        Raises:
            ValueError: If RDKit is not available.
        """
        cls.check_is_available()
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Could not read molecule from SMILES string: {smiles}")
        return cls.from_mol(mol, caption=caption, **kwargs)

    @classmethod
    def from_mol_file(cls, mol_file: Path, *, caption: Optional[str] = None, **kwargs) -> "Molecule":
        """Creates a Molecule instance from a Mol file.

        Args:
            mol_file: Path to the Mol file.
            caption: Optional descriptive text.

        Returns:
            Molecule: A new Molecule instance.

        Raises:
            ValueError: If RDKit is not available or the file cannot be read.
        """
        cls.check_is_available()
        mol = Chem.MolFromMolFile(str(mol_file))
        if mol is None:
            raise ValueError(f"Could not read molecule from Mol file: {mol_file}")
        return cls.from_mol(mol, caption=caption, **kwargs)

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
        return MediaType.Chart.MOLECULE

    def get_section(self) -> str:
        """Return section name for organization"""
        return "Molecule"

    def get_more(self) -> Optional[Dict[str, str]]:
        """Return additional information (caption)"""
        return {"caption": self.caption} if self.caption else None
