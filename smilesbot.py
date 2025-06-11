import csv
from telebot import TeleBot
import settings
from rdkit.Chem import MolFromSmiles, Descriptors, AllChem
from rdkit.DataStructs import TanimotoSimilarity

class TGbot:
    def __init__(self):
        self.bot = TeleBot(settings.TOKEN)

        @self.bot.message_handler(commands=['start'])
        def greet(message):
            self.bot.send_message(message.chat.id, 'Hi, im Dmitriis bot')

        @self.bot.message_handler(commands=['get_similar_smile'])
        def get_smiles(message):
            self.bot.send_message(message.chat.id, 'Please, write a name of smile that you want me to find similar one to')
            self.bot.register_next_step_handler(message, self.compare_with_other_mols)

    def compare_with_other_mols(self, message):
        user_smiles = message.text
        similar = self.find_similar_mols(user_smiles, 'library.csv', threshold=0.7)

        if not similar:
            self.bot.send_message(message.chat.id, "No similar molecules were found.")
        else:
            reply = "Similar molecules found:\n\n"
            for row in similar:
                reply += f"SMILES: {row['SMILES']}, Similarity: {row['Similarity']:.2f}\n"
            self.bot.send_message(message.chat.id, reply)

    def find_similar_mols(self, user_smiles, file_path, threshold=0.7):
        user_mol = MolFromSmiles(user_smiles)
        if not user_mol:
            return []

        user_fp = AllChem.GetMorganFingerprintAsBitVect(user_mol, 2, nBits=2048)

        similar_molecules = []
        with open(file_path, newline='') as f:
            reader = csv.DictReader(f)
            for row in reader:
                db_smiles = row['SMILES']
                db_mol = MolFromSmiles(db_smiles)
                if not db_mol:
                    continue

                db_fp = AllChem.GetMorganFingerprintAsBitVect(db_mol, 2, nBits=2048)
                sim = TanimotoSimilarity(user_fp, db_fp)
                if sim >= threshold:
                    row['Similarity'] = sim
                    similar_molecules.append(row)

        return sorted(similar_molecules, key=lambda x: x['Similarity'], reverse=True)

    def run(self):
        self.bot.polling()


class Dir:
    def __init__(self):
        pass

    def get_directory(self, input_file='input_smiles.csv', output_file='library.csv'):
        header = ["SMILES", "Molecular weight", "Number of rings",
                  "Number of hydrogen bond donors", "Number of hydrogen bond acceptors",
                  "Lip_pass"]

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
                writer.writerow([smi, mw, rings, h_donors, h_acceptors, "Yes" if lipinski else "No"])


if __name__ == '__main__':
    d = Dir()
    d.get_directory()

    bot = TGbot()
    bot.run()


