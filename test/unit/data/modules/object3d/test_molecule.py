import pytest
from rdkit.Chem import MolFromSmiles

from swanlab.data.modules.object3d import Molecule, Object3D


class TestMolecule:
    @pytest.fixture
    def mol(self):
        """Creates a simple RDKit Mol object."""
        return MolFromSmiles("CCO")

    @pytest.fixture
    def pdb_file(self, tmp_path):
        """Creates a dummy PDB file."""
        pdb_content = """
        ATOM      1  N     ALA A   1      -12.748   0.639   1.454  1.00  0.00           N
        ATOM      2  CA    ALA A   1      -11.372   0.337   1.839  1.00  0.00           C
        ATOM      3  C     ALA A   1      -10.488   1.542   2.310  1.00  0.00           C
        ATOM      4  O     ALA A   1      -10.920   2.633   2.779  1.00  0.00           O
        ATOM      5  CB    ALA A   1      -11.217  -0.843   2.873  1.00  0.00           C
        END
        """
        pdb_path = tmp_path / "test.pdb"
        pdb_path.write_text(pdb_content)
        return pdb_path

    @pytest.fixture
    def sdf_file(self, tmp_path):
        """Creates a dummy SDF file."""
        sdf_content = """
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  END
        """
        sdf_path = tmp_path / "test.sdf"
        sdf_path.write_text(sdf_content)
        return sdf_path

    @pytest.fixture
    def mol_file(self, tmp_path):
        """Creates a dummy Mol file."""
        mol_content = """
     RDKit          2D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.7500    1.2990    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  1  3  1  0
M  END
        """
        mol_path = tmp_path / "test.mol"
        mol_path.write_text(mol_content)
        return mol_path

    def test_from_mol(self, mol):
        """Tests creating a Molecule from an RDKit Mol object."""
        if mol:
            molecule = Molecule.from_mol(mol, caption="Ethanol")
            assert molecule.caption == "Ethanol"
            assert isinstance(molecule.pdb_data, str)

    def test_from_pdb_file(self, pdb_file):
        """Tests creating a Molecule from a PDB file."""
        molecule = Molecule.from_pdb_file(pdb_file, caption="Test PDB")
        assert molecule.caption == "Test PDB"
        assert isinstance(molecule.pdb_data, str)
        with open(pdb_file) as f:
            assert molecule.pdb_data == f.read()

    def test_from_sdf_file(self, sdf_file):
        """Tests creating a Molecule from an SDF file."""
        molecule = Molecule.from_sdf_file(sdf_file, caption="Test SDF")
        assert molecule.caption == "Test SDF"
        assert isinstance(molecule.pdb_data, str)

    def test_from_mol_file(self, mol_file):
        """Tests creating a Molecule from a Mol file."""
        molecule = Molecule.from_mol_file(mol_file, caption="Test Mol")
        assert molecule.caption == "Test Mol"
        assert isinstance(molecule.pdb_data, str)

    def test_from_smiles(self):
        """Tests creating a Molecule from a SMILES string."""
        molecule = Molecule.from_smiles("CCO", caption="Test SMILES")
        assert molecule.caption == "Test SMILES"
        assert isinstance(molecule.pdb_data, str)

    def test_parse(self, mol):
        """Tests the parse method."""
        if mol:
            molecule = Molecule.from_mol(mol)
            molecule.step = 1
            filename, buffer = molecule.parse()
            assert filename.startswith("molecule-step1-")
            assert filename.endswith(".pdb")
            assert buffer

    def test_get_chart(self, mol):
        """Tests the get_chart method."""
        if mol:
            molecule = Molecule.from_mol(mol)
            assert molecule.get_chart() == Molecule.Chart.MOLECULE

    def test_get_section(self, mol):
        """Tests the get_section method."""
        if mol:
            molecule = Molecule.from_mol(mol)
            assert molecule.get_section() == "Molecule"

    def test_get_more(self, mol):
        """Tests the get_more method."""
        if mol:
            molecule = Molecule.from_mol(mol, caption="Test")
            assert molecule.get_more() == {"caption": "Test"}
            molecule = Molecule.from_mol(mol)
            assert molecule.get_more() is None


class TestObject3DWithMolecule:
    @pytest.fixture
    def mol(self):
        """Creates a simple RDKit Mol object."""
        return MolFromSmiles("CCO")  # Ethanol

    def test_object3d_from_mol(self, mol):
        """Tests Object3D dispatch with a Molecule object."""
        if mol:
            obj = Object3D(mol, caption="Ethanol Object")
            assert isinstance(obj, Molecule)
            assert obj.caption == "Ethanol Object"
