from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, Lipinski
import json
import base64
from io import BytesIO

def molecule_properties(smiles: str) -> dict:
    # Initialize output dictionary
    output_json = {
        "smiles": smiles,
        "molecular_weight": 0.0,
        "ring_count": 0,
        "h_donors": 0,
        "h_acceptors": 0,
        "lipinski_compliant": False,
        "lipinski_rules": {},
        "icon": ""
    }

    try:
        # Parse SMILES and create molecule
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")

        # Calculate molecular properties
        output_json["molecular_weight"] = round(Descriptors.MolWt(mol), 2)
        output_json["ring_count"] = Lipinski.RingCount(mol)
        output_json["h_donors"] = Lipinski.NumHDonors(mol)
        output_json["h_acceptors"] = Lipinski.NumHAcceptors(mol)

        # Evaluate Lipinski's Rule of Five
        lipinski_rules = {
            "MolWt < 500": output_json["molecular_weight"] < 500,
            "NumHDonors <= 5": output_json["h_donors"] <= 5,
            "NumHAcceptors <= 10": output_json["h_acceptors"] <= 10,
            "LogP <= 5": Lipinski.MolLogP(mol) <= 5
        }
        output_json["lipinski_rules"] = lipinski_rules
        output_json["lipinski_compliant"] = all(lipinski_rules.values())

        # Generate molecule image as base64 string
        img = Draw.MolToImage(mol)
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        output_json["icon"] = base64.b64encode(buffered.getvalue()).decode("utf-8")

        return output_json

    except Exception as e:
        return {"error": f"Failed to process SMILES: {str(e)}"}

# Example usage
if __name__ == "__main__":
    smiles = input("Enter SMILES string: ")
    result = molecule_properties(smiles)
    print(json.dumps(result, indent=2))
