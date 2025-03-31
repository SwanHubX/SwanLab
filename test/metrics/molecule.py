from rdkit import Chem

import swanlab

swanlab.init(project="molecule", public=True)

# from rdkit.Chem.Mol
chem = Chem.MolFromSmiles("CCO")
molecule = swanlab.Object3D(chem, caption="cco")
swanlab.log({"example": molecule})

# from file path
mol2 = swanlab.Object3D("./assets/molecule.example.pdb")
swanlab.log({"file": mol2})

# this file is from https://github.com/SwanHubX/SwanLab/pull/477
# You should download https://github.com/SwanHubX/SwanLab/files/15022006/5p21.pdb.zip
# and unzip it and rename to big_mol.example.pdb
#
# mol3 = swanlab.Object3D("./big_mol.example.pdb")
# swanlab.log({"big_mol": mol3})
