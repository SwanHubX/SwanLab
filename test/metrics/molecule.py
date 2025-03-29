from rdkit import Chem

import swanlab

swanlab.init(project="molecule", public=True)

chem = Chem.MolFromSmiles("CCO")

molecule = swanlab.data.modules.object3d.Molecule.from_mol(chem, caption="cco")
swanlab.log({"example": molecule})
