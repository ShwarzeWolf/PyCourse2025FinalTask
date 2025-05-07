import json
import base64
from io import BytesIO

from rdkit import Chem
from rdkit.Chem import (
    Descriptors,
    Draw, 
    Lipinski,
)

const MAX_molecular_wight = 500
const MAX_h_donors = 5
const MAX_h_acceptors = 10

def get_molecule_properties(smiles: str) -> dict:
    
    molecule_properties = {
        "smiles": smiles,
        "molecular_weight": round(Descriptors.MolWt(mol), 2),
        "ring_count": Lipinski.RingCount(mol),
        "h_donors": Lipinski.NumHDonors(mol),
        "h_acceptors": Lipinski.NumHAcceptors(mol),
    }

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")

    lipinski_rules = {
        "MolWt < 500": molecule_properties["molecular_weight"] < 500,
        "NumHDonors <= 5": molecule_properties["h_donors"] <= 5,
        "NumHAcceptors <= 10": molecule_properties["h_acceptors"] <= 10,
        "LogP <= 5": Lipinski.MolLogP(mol) <= 5
    }
    molecule_properties["Lipinski_rules"] = lipinski_rules
    molecule_properties["lipinski_compliant"] = all(lipinski_rules.values())

    img = Draw.MolToImage(mol)
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    molecule_properties["icon"] = base64.b64encode(buffered.getvalue()).decode("utf-8")

    return molecule_properties


if __name__ == "__main__":
    smiles = input("Enter SMILES string: ")
    result = molecule_properties(smiles)
    print(json.dumps(result, indent=2))
