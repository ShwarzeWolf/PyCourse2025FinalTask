from io import BytesIO
import json
import base64
from typing import Optional

from rdkit import Chem
from rdkit.Chem import (
    Descriptors,
    Draw,
    Lipinski,
)


MAX_MOLECULAR_WEIGHT = 500
PRECISION = 2


def get_molecule_properties(smiles: str) -> Optional[dict]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Invalid SMILES string")
        return

    molecular_properties = {
        "smiles": smiles,
        "molecular_weight": round(Descriptors.MolWt(mol), PRECISION),
        "ring_count": Lipinski.RingCount(mol),
        "h_donors": Lipinski.NumHDonors(mol),
        "h_acceptors": Lipinski.NumHAcceptors(mol),
    }

    # Evaluate Lipinski's Rule of Five
    lipinski_rules = {
        "MolWt < 500": molecular_properties["molecular_weight"] < MAX_MOLECULAR_WEIGHT,
        "NumHDonors <= 5": molecular_properties["h_donors"] <= 5,
        "NumHAcceptors <= 10": molecular_properties["h_acceptors"] <= 10,
        "LogP <= 5": Lipinski.MolLogP(mol) <= 5,
    }
    molecular_properties["lipinski_rules"] = lipinski_rules
    molecular_properties["lipinski_compliant"] = all(lipinski_rules.values())

    img = Draw.MolToImage(mol)
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    molecular_properties["icon"] = base64.b64encode(buffered.getvalue()).decode("utf-8")

    return molecular_properties


if __name__ == "__main__":
    smiles = input("Enter SMILES string: ")
    result = get_molecule_properties(smiles)
    print(json.dumps(result, indent=2))
