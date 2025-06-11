import base64
from io import BytesIO
import csv

from rdkit import Chem
from rdkit.Chem import (
    Descriptors,
    Draw, 
    Lipinski,
    MolFromSmiles,
)

MAX_MOLECULAR_WEIGHT = 500
MAX_H_DONORS = 5
MAX_H_ACCEPTORS = 10


def get_molecule_properties(smiles: str) -> dict:

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")
    
    molecule_properties = {
        "smiles": smiles,
        "molecular_weight": round(Descriptors.MolWt(mol), 2),
        "ring_count": Lipinski.RingCount(mol),
        "h_donors": Lipinski.NumHDonors(mol),
        "h_acceptors": Lipinski.NumHAcceptors(mol),
    }

    lipinski_rules = {
        "MolWt < 500": molecule_properties["molecular_weight"] < MAX_MOLECULAR_WEIGHT,
        "NumHDonors <= 5": molecule_properties["h_donors"] <= MAX_H_DONORS,
        "NumHAcceptors <= 10": molecule_properties["h_acceptors"] <= MAX_H_ACCEPTORS,
        "LogP <= 5": Lipinski.MolLogP(mol) <= 5
    }
    molecule_properties["Lipinski_rules"] = lipinski_rules
    molecule_properties["lipinski_compliant"] = all(lipinski_rules.values())

    img = Draw.MolToImage(mol)
    buffered = BytesIO()
    img.save(buffered, format="PNG")
    molecule_properties["icon"] = base64.b64encode(buffered.getvalue()).decode("utf-8")

    return molecule_properties


class Dir:
    def get_directory(self, input_file="input_smiles.csv", output_file="library.csv"):
        header = ["SMILES", "Molecular weight", "Number of rings",
                  "Number of hydrogen bond donors", "Number of hydrogen bond acceptors",
                  "Lip_pass"]

        with open(input_file, newline="") as csvfile:
            reader = csv.DictReader(csvfile)
            molecules = [row["SMILES"] for row in reader]

        with open(output_file, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(header)

            for smi in molecules:
                mol = MolFromSmiles(smi)
                if not mol:
                    continue

                mw = Descriptors.MolWt(mol)
                rings = mol.GetRingInfo().NumRings()
                h_donors = Descriptors.NumHDonors(mol)
                h_acceptors = Descriptors.NumHAcceptors(mol)

                lipinski = (
                    mw < 500 and
                    h_donors <= 5 and
                    h_acceptors <= 10
                )
                writer.writerow([smi, mw, rings, h_donors, h_acceptors, "Yes" if lipinski else "No"])
