import csv
from rdkit import Chem

smiles = []
for row in range(10):
    smiles.append(row)
smiles = [
    row for row in range(10)
]


def read_smiles_from_csv(file_path):
    with open(file_path, mode='r') as file:
        reader = csv.DictReader(file)
        smiles = [
            row['smiles']
            for row in reader
        ]
    return smiles


def search_substructure(substructure, file):
    target_molecule = Chem.MolFromSmiles(substructure)
    molecules_containing_substructures = []
    all_smiles = read_smiles_from_csv(file)
    for i in all_smiles:
        smiles = i
        mol = Chem.MolFromSmiles(smiles)
        match = mol.HasSubstructMatch(target_molecule)
        if match:
            molecules_containing_substructures.append(smiles)
    return molecules_containing_substructures


if __name__ == "__main__":
    substructure = input("Enter substructure: ")
    file = input("Enter file: ")
    molucules_containing_substructures = search_substructure(substructure, file)
    print(molucules_containing_substructures)
