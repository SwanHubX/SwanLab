from rdkit import Chem

import swanlab

swanlab.init(project="molecule", public=True)

chem = Chem.MolFromSmiles("CCO")

molecule = swanlab.Object3D(chem, caption="cco")
swanlab.log({"example": molecule})

mol2 = swanlab.Object3D("./molecule.example.pdb")
swanlab.log({"file": mol2})
