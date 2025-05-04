from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, Lipinski
import json

smiles = input()
mol = Chem.MolFromSmiles(smiles)

output_json = '{"smiles" : "", "wight":"", "rings":"", "H donor":"", "H acceptors":"", "Lipinski":"", "Icon":""}'

output_json["wight"] = Descriptors.MolWt(mol)
output_json["rings"] = Lipinski.RingCount(mol)
output_json["H donor"] = Lipinski.NumHDonors(mol)
output_json["H Acceptors"] = Lipinski.NumHAcceptors(mol)

lipinski_rule = {
        'MolWt < 500': output_json["wight"] < 500,
        'NumHDonors <= 5': output_json["H donor"] <= 5,
        'NumHAcceptors <= 10': output_json["H Acceptors"] <= 10,
        'LogP <= 5': Lipinski.MolLogP(mol) <= 5
    }
output_json["Lipinski"] = all(lipinski_rule.values())

output_json["Icon"] = Draw.MolToImage(mol)