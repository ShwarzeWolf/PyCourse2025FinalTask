import csv
import settings
from rdkit.Chem import MolFromSmiles, Descriptors, AllChem
from rdkit.DataStructs import TanimotoSimilarity



class Dir:
    def __init__(self):
        pass

    def get_directory(self, input_file='input_smiles.csv', output_file='library.csv'):
        header = ["SMILES", "Molecular weight", "Number of rings",
                  "Number of hydrogen bond donors", "Number of hydrogen bond acceptors",
                  "Соответствует ли критерию Липински"]

        with open(input_file, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            molecules = [row['SMILES'] for row in reader]

        with open(output_file, 'w', newline='') as f:
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
                writer.writerow([smi, mw, rings, h_donors, h_acceptors, "Да" if lipinski else "Нет"])


if __name__ == '__main__':
    d = Dir()
    d.get_directory()