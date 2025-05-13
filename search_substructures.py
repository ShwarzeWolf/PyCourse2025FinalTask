import rdkit
import csv
from keras.src.ops import append
from rdkit import Chem

def read_smiles_from_csv(file_path):
    smiles_list = []
    with open(file_path, mode='r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            smiles_list.append(row['smiles'])
    return smiles_list
def substructure_search(substructure, file):
    benzene = Chem.MolFromSmiles(substructure)
    Tmol = []
    all_smiles = read_smiles_from_csv(file)
    for i in range(len(all_smiles)):
        smiles = all_smiles[i]
        mol = Chem.MolFromSmiles(smiles)
        match = mol.HasSubstructMatch(benzene)
        if (match == True):
            Tmol.append(smiles)
    return Tmol

if __name__ == "__main__":
    substructure = input("Enter substructure: ")
    file = input("Enter file: ")
    molecules = substructure_search(substructure, file)
    print(molecules)