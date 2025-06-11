import csv

from rdkit.Chem import MolFromSmiles, AllChem
from rdkit.DataStructs import TanimotoSimilarity


def find_similar_mols(user_smiles, file_path, threshold=0.7):
    user_mol = MolFromSmiles(user_smiles)
    if not user_mol:
        return []

    user_fp = AllChem.GetMorganFingerprintAsBitVect(user_mol, 2, nBits=2048)

    similar_molecules = []
    with open(file_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            db_smiles = row["SMILES"]
            db_mol = MolFromSmiles(db_smiles)
            if not db_mol:
                continue

            db_fp = AllChem.GetMorganFingerprintAsBitVect(db_mol, 2, nBits=2048)
            sim = TanimotoSimilarity(user_fp, db_fp)
            if sim >= threshold:
                row["Similarity"] = sim
                similar_molecules.append(row)

    return sorted(similar_molecules, key=lambda x: x["Similarity"], reverse=True)
