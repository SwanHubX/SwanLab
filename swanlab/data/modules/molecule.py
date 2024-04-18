# -*- coding: utf-8 -*-
from .base import BaseType
from typing import Union, List, TYPE_CHECKING
import os
import shutil
import secrets
import string
import io

if TYPE_CHECKING:
    from typing import TextIO
    import rdkit.Chem

    RDKitDataType = Union[str, "rdkit.Chem.rdchem.Mol"]


def generate_id(length: int = 16) -> str:
    """Generate a random base-36 string of `length` digits."""
    # There are ~2.8T base-36 8-digit strings. If we generate 210k ids,
    # we'll have a ~1% chance of collision.
    alphabet = string.ascii_lowercase + string.digits
    return "".join(secrets.choice(alphabet) for _ in range(length))


class Molecule(BaseType):
    """SwanLab class for 3D Molecular data.

    Arguments:
        data_or_path: (string, io)
            Molecule can be initialized from a file name or an io object.
            Molecules 可以被一个文件路径或一个IO对象初始化, 如swanlab.Molecule("path/to/file")
        caption: (string)
            Caption associated with the molecule for display.
            与Molecule相关的标题, 用于在GUI上显示。swanlab.Molecule("path/to/file", caption="Ethanol")
        file_type: (string)
            Type of the file. If not provided, the file type will be inferred from the file extension.
            文件类型。如果未提供，文件类型将从文件扩展名中推断出来。
    """

    # IO流支持的文件类型
    SUPPORTED_TYPES = {
        "pdb",
        "pqr",
        "mmcif",
        "mcif",
        "cif",
        "sdf",
        "sd",
        "gro",
        "mol2",
        "mmtf",
    }

    # RDKIT支持的文件类型
    SUPPORTED_RDKIT_TYPES = {"mol", "sdf"}

    def __init__(
        self,
        data_or_path: Union[str, "TextIO", List["Molecule"]],
        caption: str = None,
        file_type: str = None,
    ) -> None:

        super().__init__(data_or_path)
        self.molecule_data = None
        self.caption = self.__convert_caption(caption)
        self.file_type = file_type

    def get_data(self):
        # 如果传入的是Molecule类列表
        if isinstance(self.value, list):
            return self.get_data_list()
        # 数据预处理
        ext = self.__preprocess(self.value)
        random_id = generate_id()

        # 生成保存路径
        save_dir = os.path.join(self.settings.static_dir, self.tag)
        save_name = f"Molecule-step{self.step}-{random_id}.{ext}"
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)
        save_path = os.path.join(save_dir, save_name)

        # 保存分子数据到指定目录
        self.__save(save_path)
        return save_name

    def __preprocess(self, data_or_path):
        if hasattr(data_or_path, "name"):
            data_or_path = data_or_path.name

        # 如果传入的是IO对象
        if hasattr(data_or_path, "read"):
            if hasattr(data_or_path, "seek"):
                data_or_path.seek(0)

            ext = self.file_type

            # 如果没有传入file_type参数
            if ext is None:
                raise TypeError("When using the io object, the file_type keyword argument must be passed.")

            if ext not in Molecule.SUPPORTED_TYPES:
                raise TypeError("Molecule 3D only supports files of the type: " + ", ".join(Molecule.SUPPORTED_TYPES))

            # 分子数据为IO对象读取的内容
            self.molecule_data = data_or_path

        # 如果传入的是文件路径
        elif isinstance(data_or_path, str):
            ext = os.path.splitext(data_or_path)[1][1:]
            if ext not in Molecule.SUPPORTED_TYPES:
                raise TypeError("Molecule 3D only supports files of the type: " + ", ".join(Molecule.SUPPORTED_TYPES))
            # 分子数据为文件路径
            self.molecule_data = data_or_path
        else:
            raise TypeError("The data passed to Melocule must be a file name or file object.")

        return ext

    def __save(self, save_path):

        try:
            # 如果传入的是IO对象
            if hasattr(self.molecule_data, "read"):
                # 将分子数据写入临时文件
                with open(save_path, "w") as f:
                    f.write(self.molecule_data.read())
            # 如果传入的是文件路径, 复制文件到指定目录
            else:
                shutil.copyfile(self.molecule_data, save_path)

        except Exception as e:
            raise TypeError(f"Could not save the Molecule to the path: {save_path}") from e

    @classmethod
    def from_rdkit(
        cls,
        data_or_path: "RDKitDataType",
        caption: str = None,
        convert_to_3d_and_optimize: bool = True,
        mmff_optimize_molecule_max_iterations: int = 200,
    ) -> "Molecule":
        """Convert RDKit-supported file/object types to swanlab.Molecule.
        将RDKit支持的文件/对象类型转换为swanlab.Molecule。

        Arguments:
            data_or_path: (string, rdkit.Chem.rdchem.Mol)
                Molecule can be initialized from a file name or an rdkit.Chem.rdchem.Mol object.
                Molecule可以从文件名或rdkit.Chem.rdchem.Mol对象初始化。
                如swanlab.Molecule.from_rdkit("path/to/file")

            caption: (string)
                Caption associated with the molecule for display.
                与Molecule相关的标题, 用于在GUI上显示。如swanlab.Molecule.from_rdkit("path/to/file", caption="Ethanol")

            convert_to_3d_and_optimize: (bool)
                Convert to rdkit.Chem.rdchem.Mol with 3D coordinates.
                This is an expensive operation that may take a long time for complicated molecules.
                将具有3D坐标的rdkit.Chem.rdchem.Mol转换。
                这是一个比较耗时的操作，对于复杂分子可能需要很长时间。

            mmff_optimize_molecule_max_iterations: (int)
                Number of iterations to use in rdkit.Chem.AllChem.MMFFOptimizeMolecule
                在rdkit.Chem.AllChem.MMFFOptimizeMolecule中要使用的迭代次数
        """

        try:
            import rdkit
        except ImportError as e:
            raise TypeError("swanlab.Molecule requires the rdkit pypi package. Install with 'pip install rdkit'.")

        from rdkit import Chem
        from rdkit.Chem import AllChem
        import pathlib

        # 如果传入的是文件路径
        if isinstance(data_or_path, str):
            path = pathlib.Path(data_or_path)
            ext = path.suffix.split(".")[-1]
            if ext not in Molecule.SUPPORTED_RDKIT_TYPES:
                raise TypeError(
                    "swanlab.Molecule.from_rdkit only supports files of the type: "
                    + ", ".join(Molecule.SUPPORTED_RDKIT_TYPES)
                )
            # 如果后缀是sdf, 则使用SDMolSupplier进行读取
            if ext == "sdf":
                with Chem.SDMolSupplier(data_or_path) as supplier:
                    molecule = next(supplier)  # 只获取第一个分子
            else:
                molecule = getattr(Chem, f"MolFrom{ext.capitalize()}File")(data_or_path)
        elif isinstance(data_or_path, Chem.rdchem.Mol):
            molecule = data_or_path
        else:
            raise TypeError("Data must be file name or an rdkit.Chem.rdchem.Mol object")

        if convert_to_3d_and_optimize:
            molecule = Chem.AddHs(molecule)
            AllChem.EmbedMolecule(molecule)
            AllChem.MMFFOptimizeMolecule(
                molecule,
                maxIters=mmff_optimize_molecule_max_iterations,
            )

        # 转换为Molecule支持的pdb格式
        pdb_block = Chem.rdmolfiles.MolToPDBBlock(molecule)

        return cls(io.StringIO(pdb_block), caption=caption, file_type="pdb")

    @classmethod
    def from_smiles(
        cls,
        data: str,
        caption: str = None,
        sanitize: bool = True,
        convert_to_3d_and_optimize: bool = True,
        mmff_optimize_molecule_max_iterations: int = 200,
    ) -> "Molecule":
        """Convert SMILES string to swanlab.Molecule.
        将SMILES字符串转换为swanlab.Molecule。

        Arguments:
            data: (string)
                SMILES string. 如swanlab.Molecule.from_smiles("CCO")

            caption: (string)
                Caption associated with the molecule for display
                与Molecule相关的标题, 用于在GUI上显示。

            sanitize: (bool)
                Check if the molecule is chemically reasonable by the RDKit's definition.
                通过RDKit的定义检查分子是否在化学上合理。

            convert_to_3d_and_optimize: (bool)
                Convert to rdkit.Chem.rdchem.Mol with 3D coordinates.
                This is an expensive operation that may take a long time for complicated molecules.
                将具有3D坐标的rdkit.Chem.rdchem.Mol转换。
                这是一个比较耗时的操作，对于复杂分子可能需要很长时间。

            mmff_optimize_molecule_max_iterations: (int)
                Number of iterations to use in rdkit.Chem.AllChem.MMFFOptimizeMolecule
                在rdkit.Chem.AllChem.MMFFOptimizeMolecule中要使用的迭代次数
        """

        try:
            import rdkit
        except ImportError as e:
            raise TypeError("swanlab.Molecule requires the rdkit pypi package. Install with 'pip install rdkit'.")

        from rdkit import Chem

        molecule = Chem.MolFromSmiles(data, sanitize=sanitize)
        if molecule is None:
            raise TypeError("Unable to parse the SMILES string.")

        return cls.from_rdkit(
            data_or_path=molecule,
            caption=caption,
            convert_to_3d_and_optimize=convert_to_3d_and_optimize,
            mmff_optimize_molecule_max_iterations=mmff_optimize_molecule_max_iterations,
        )

    def get_more(self, *args, **kwargs) -> dict:
        """返回config数据"""
        # 如果传入的是Molecule类列表
        if isinstance(self.value, list):
            return self.get_more_list()
        else:
            return (
                {
                    "caption": self.caption,
                }
                if self.caption is not None
                else None
            )

    def expect_types(self, *args, **kwargs) -> list:
        return ["str", "TextIO"]

    def __convert_caption(self, caption):
        """将caption转换为字符串"""
        # 如果类型是字符串，则不做转换
        if isinstance(caption, str):
            caption = caption
        # 如果类型是数字，则转换为字符串
        elif isinstance(caption, (int, float)):
            caption = str(caption)
        # 如果类型是None，则转换为默认字符串
        elif caption is None:
            caption = None
        else:
            raise TypeError("caption must be a string, int or float.")
        return caption.strip() if caption else None

    def get_namespace(self, *args, **kwargs) -> str:
        """设定分组名"""
        return "Molecule"

    def get_chart_type(self) -> str:
        """设定图表类型"""
        return self.chart.melocule
