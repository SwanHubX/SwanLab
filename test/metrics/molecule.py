"""
@author: xj63
@file: molecule.py
@time: 2025/3/13 13:30
@description: 测试上传分子对象

NOTE 你需要下载下面的文件放到当前文件目录的assets文件夹下，才能运行这个测试
- 3D分子pdb文件: 此文件存放于 https://github.com/SwanHubX/SwanLab/pull/477 ，可通过 https://github.com/user-attachments/files/19605387/molecule.example.pdb.zip 下载
"""

# noinspection PyPackageRequirements
from rdkit import Chem

import swanlab

swanlab.init(project="molecule", public=True)

# from rdkit.Chem.Mol
chem = Chem.MolFromSmiles("CCO")
cco = swanlab.Molecule.from_mol(chem, caption="cco")

# from file path
file = swanlab.Molecule.from_pdb_file("./assets/molecule.example.pdb", caption="file")

# from pdb data
with open("./assets/molecule.example.pdb") as f:
    data = swanlab.Molecule(f.read(), caption="data")

# upload
swanlab.log({"file": file, "data": data, "example": cco})
