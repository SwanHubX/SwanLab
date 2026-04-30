"""
@author: cunyue
@file: test_molecule.py
@time: 2026/4/30
@description: Molecule TransformMedia 单元测试
"""

import hashlib

import pytest
from google.protobuf.timestamp_pb2 import Timestamp

from swanlab.proto.swanlab.metric.column.v1.column_pb2 import ColumnType
from swanlab.proto.swanlab.metric.data.v1.data_pb2 import MediaItem
from swanlab.sdk.internal.run.transforms.molecule import Molecule

# RDKit 在 CI 环境中可能不可用，需要条件跳过
try:
    from rdkit import Chem

    _HAS_RDKIT = True
except ImportError:
    _HAS_RDKIT = False

requires_rdkit = pytest.mark.skipif(not _HAS_RDKIT, reason="RDKit not installed")


# ---------------------------------- Fixtures ----------------------------------


@pytest.fixture
def ethanol_smiles():
    return "CCO"


@pytest.fixture
def benzene_smiles():
    return "c1ccccc1"


@pytest.fixture
def rdkit_mol():
    """RDKit Mol 对象（苯）"""
    mol = Chem.MolFromSmiles("c1ccccc1")  # type: ignore[union-attr]
    return mol


@pytest.fixture
def pdb_file(tmp_path, ethanol_smiles):
    """生成一个临时 PDB 文件"""
    mol = Chem.MolFromSmiles(ethanol_smiles)  # type: ignore[union-attr]
    Chem.AllChem.EmbedMolecule(mol)  # type: ignore[union-attr]
    pdb_content = Chem.MolToPDBBlock(mol)  # type: ignore[union-attr]
    path = tmp_path / "test.pdb"
    path.write_text(pdb_content)
    return path


@pytest.fixture
def sdf_file(tmp_path, rdkit_mol):
    """生成一个临时 SDF 文件"""
    path = tmp_path / "test.sdf"
    writer = Chem.SDWriter(str(path))  # type: ignore[union-attr]
    Chem.AllChem.EmbedMolecule(rdkit_mol)  # type: ignore[union-attr]
    writer.write(rdkit_mol)
    writer.close()
    return path


@pytest.fixture
def mol_file(tmp_path, rdkit_mol):
    """生成一个临时 MOL 文件"""
    path = tmp_path / "test.mol"
    Chem.AllChem.EmbedMolecule(rdkit_mol)  # type: ignore[union-attr]
    Chem.MolToMolFile(rdkit_mol, str(path))  # type: ignore[union-attr]
    return path


# ---------------------------------- 构造测试 ----------------------------------


class TestMoleculeInit:
    @requires_rdkit
    def test_from_smiles(self, ethanol_smiles):
        """从 SMILES 字符串构造"""
        mol = Molecule(ethanol_smiles, caption="Ethanol")
        assert mol.caption == "Ethanol"
        assert len(mol.buffer.getvalue()) > 0

    @requires_rdkit
    def test_from_rdkit_mol(self, rdkit_mol):
        """从 RDKit Mol 对象构造"""
        mol = Molecule(rdkit_mol, caption="Benzene")
        assert mol.caption == "Benzene"
        assert len(mol.buffer.getvalue()) > 0

    @requires_rdkit
    def test_from_pdb_file(self, pdb_file):
        """从 PDB 文件路径构造"""
        mol = Molecule(str(pdb_file))
        assert len(mol.buffer.getvalue()) > 0

    @requires_rdkit
    def test_from_pdb_path_object(self, pdb_file):
        """从 Path 对象构造 PDB"""
        mol = Molecule(pdb_file)
        assert len(mol.buffer.getvalue()) > 0

    @requires_rdkit
    def test_from_sdf_file(self, sdf_file):
        """从 SDF 文件构造"""
        mol = Molecule(str(sdf_file))
        assert len(mol.buffer.getvalue()) > 0

    @requires_rdkit
    def test_from_mol_file(self, mol_file):
        """从 MOL 文件构造"""
        mol = Molecule(str(mol_file))
        assert len(mol.buffer.getvalue()) > 0

    @requires_rdkit
    def test_caption_stored(self, ethanol_smiles):
        """caption 被正确保存"""
        mol = Molecule(ethanol_smiles, caption="test caption")
        assert mol.caption == "test caption"

    @requires_rdkit
    def test_caption_none_by_default(self, ethanol_smiles):
        """caption 默认为 None"""
        mol = Molecule(ethanol_smiles)
        assert mol.caption is None


# ---------------------------------- 构造错误测试 ----------------------------------


class TestMoleculeInitErrors:
    def test_unsupported_type_raises(self):
        """传入不支持的类型应抛出 TypeError"""
        with pytest.raises(TypeError, match="Unsupported input type"):
            Molecule(12345)  # type: ignore

    @requires_rdkit
    def test_invalid_smiles_raises(self):
        """无效的 SMILES 字符串应抛出 ValueError"""
        with pytest.raises(ValueError, match="Could not parse SMILES"):
            Molecule("INVALID_SMILES_STRING!!!")

    def test_nonexistent_file_raises(self):
        """不存在的文件应抛出 TypeError（因为 str 既不是文件也不是 SMILES，走 SMILES 解析失败）"""
        with pytest.raises((ValueError, TypeError)):
            Molecule("/nonexistent/path/mol.pdb")

    @requires_rdkit
    def test_unsupported_file_type_raises(self, tmp_path):
        """不支持的文件扩展名应抛出 ValueError"""
        unsupported = tmp_path / "test.xyz"
        unsupported.write_text("data")
        with pytest.raises(ValueError, match="Unsupported file type"):
            Molecule(str(unsupported))


# ---------------------------------- 套娃加载测试 ----------------------------------


class TestMoleculeNesting:
    @requires_rdkit
    def test_wrap_copies_buffer(self, ethanol_smiles):
        """套娃加载复用内层 buffer"""
        inner = Molecule(ethanol_smiles, caption="inner")
        outer = Molecule(inner)
        assert outer.buffer is inner.buffer
        assert outer.caption == "inner"

    @requires_rdkit
    def test_outer_caption_overrides_inner(self, ethanol_smiles):
        """外层 caption 优先级高于内层"""
        inner = Molecule(ethanol_smiles, caption="inner")
        outer = Molecule(inner, caption="outer")
        assert outer.caption == "outer"

    @requires_rdkit
    def test_inner_caption_used_when_outer_none(self, ethanol_smiles):
        """外层 caption 为 None 时使用内层 caption"""
        inner = Molecule(ethanol_smiles, caption="inner")
        outer = Molecule(inner, caption=None)
        assert outer.caption == "inner"


# ---------------------------------- column_type 测试 ----------------------------------


class TestMoleculeColumnType:
    def test_column_type(self):
        assert Molecule.column_type() == ColumnType.COLUMN_TYPE_MOLECULE


# ---------------------------------- build_data_record 测试 ----------------------------------


class TestMoleculeBuildDataRecord:
    @requires_rdkit
    def test_build_data_record_structure(self, ethanol_smiles, tmp_path):
        """build_data_record 返回正确结构的 MediaRecord"""
        mol = Molecule(ethanol_smiles)
        item = mol.transform(step=1, path=tmp_path)
        ts = Timestamp()
        record = Molecule.build_data_record(key="mol", step=1, timestamp=ts, data=[item])
        assert record.key == "mol"
        assert record.step == 1
        assert record.type == ColumnType.COLUMN_TYPE_MOLECULE
        assert len(record.value.items) == 1


# ---------------------------------- transform 测试 ----------------------------------


class TestMoleculeTransform:
    @requires_rdkit
    def test_transform_returns_media_item(self, ethanol_smiles, tmp_path):
        """transform 返回 MediaItem"""
        item = Molecule(ethanol_smiles).transform(step=1, path=tmp_path)
        assert isinstance(item, MediaItem)

    @requires_rdkit
    def test_transform_pdb_filename(self, ethanol_smiles, tmp_path):
        """transform 输出 .pdb 文件"""
        item = Molecule(ethanol_smiles).transform(step=1, path=tmp_path)
        assert item.filename.endswith(".pdb")
        assert (tmp_path / item.filename).exists()

    @requires_rdkit
    def test_transform_sha256_correct(self, ethanol_smiles, tmp_path):
        """sha256 与落盘文件内容一致"""
        item = Molecule(ethanol_smiles).transform(step=1, path=tmp_path)
        content = (tmp_path / item.filename).read_bytes()
        assert item.sha256 == hashlib.sha256(content).hexdigest()

    @requires_rdkit
    def test_transform_size_correct(self, ethanol_smiles, tmp_path):
        """size 与落盘文件字节数一致"""
        item = Molecule(ethanol_smiles).transform(step=1, path=tmp_path)
        assert item.size == len((tmp_path / item.filename).read_bytes())

    @requires_rdkit
    def test_transform_caption_empty_when_none(self, ethanol_smiles, tmp_path):
        """caption 为 None 时返回空字符串"""
        item = Molecule(ethanol_smiles).transform(step=1, path=tmp_path)
        assert item.caption == ""

    @requires_rdkit
    def test_transform_caption_preserved(self, ethanol_smiles, tmp_path):
        """caption 有值时正确传递"""
        item = Molecule(ethanol_smiles, caption="hello").transform(step=1, path=tmp_path)
        assert item.caption == "hello"


# ---------------------------------- 工厂方法测试 ----------------------------------


class TestMoleculeFactoryMethods:
    @requires_rdkit
    def test_from_mol(self, rdkit_mol):
        """from_mol 从 RDKit Mol 对象创建"""
        mol = Molecule.from_mol(rdkit_mol, caption="Benzene")
        assert mol.caption == "Benzene"
        assert len(mol.buffer.getvalue()) > 0

    @requires_rdkit
    def test_from_smiles(self, ethanol_smiles):
        """from_smiles 从 SMILES 字符串创建"""
        mol = Molecule.from_smiles(ethanol_smiles, caption="Ethanol")
        assert mol.caption == "Ethanol"

    @requires_rdkit
    def test_from_pdb(self, pdb_file):
        """from_pdb 从 PDB 文件创建"""
        mol = Molecule.from_pdb(str(pdb_file), caption="From PDB")
        assert mol.caption == "From PDB"

    @requires_rdkit
    def test_from_sdf(self, sdf_file):
        """from_sdf 从 SDF 文件创建"""
        mol = Molecule.from_sdf(str(sdf_file), caption="From SDF")
        assert mol.caption == "From SDF"

    @requires_rdkit
    def test_from_mol_file(self, mol_file):
        """from_mol_file 从 MOL 文件创建"""
        mol = Molecule.from_mol_file(str(mol_file), caption="From MOL")
        assert mol.caption == "From MOL"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
