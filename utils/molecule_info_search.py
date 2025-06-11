import csv

from rdkit.Chem import MolFromSmiles, Descriptors


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
