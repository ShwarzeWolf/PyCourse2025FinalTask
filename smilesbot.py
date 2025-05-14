import csv

from telebot import TeleBot
from telegram import update
import pandas as pd
from rdkit.Chem import MolFromSmiles, AllChem
from rdkit.DataStructs import TanimotoSimilarity

import settings

bot = TeleBot(settings.TOKEN)


@bot.message_handler(commands=['start'])
def greet(message):
    print(message)
    bot.send_message(
        message.chat.id,
        'Hi, im Dmitriis bot'
    )


@bot.message_handler(commands=['help'])
def get_help(message):
    print(message)
    bot.send_message(
        message.chat.id,
        'Sorry, i cant help you yet'
    )
  
    
@bot.message_handler(commands=['get similar smile'])
def get_smiles_to_find_similar_one(message):
    bot.send_message(
        message.chat.id,
        'Please, write a name of smile that you want me to find similar one to'
    )
    users_smile = update.message.text
    
    directory = Library('data.csv')
    directory.write_file_header()
    directory.get_similar_smiles(users_smile)
    return users_smile


class Library:
    def __init__(self, filepath):
        self.filepath = filepath

    # def write_file_header(self):
    #     data_rows = [
    #         'SMILES', 'Molecular weight', 'Number of rings',
    #         'Number of hydrogen bond donors', 'Number of hydrogen bond acceptors',
    #         'Lipinski pass',
    #     ]
    #
    #     with open(self.filepath, 'w') as file:
    #         writer = csv.writer(file)
    #         writer.writerow(data_rows)

    def get_similar_smiles(self, target_smiles, threshold=0.5):
        df = pd.read_csv(self.filepath)

        # TODO add this functionality for calculations
        # convert SMILES into molecule object
        # targetMol = MolFromSmiles(target_SMILES)
        # mol = MolFromSmiles(smiles)
        #
        # # calculate ECFP4 fingerprints for 2 molecules
        # target_fps = AllChem.GetMorganFingerprintAsBitVect(targetMol, 2, nBits=2048)  # [001001000101001]
        # fps = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        #
        # # find a similarity rate using rdkit Tanimoto similarity function
        # similarity_score = TanimotoSimilarity(
        #     target_fps, fps
        # )  # [0, 1]
        #
        # print(similarity_score)

        # Понять, как правильно считать
        row = df[df['SMILES'] == target_smiles].iloc[0]
        v1 = (
            row['Molecular weight']
            + row['Number of rings']
            + row['Number of hydrogen bond donors']
            + row['Number of hydrogen bond acceptors']
        )

        similar_smiles = []
        for index, row in df.iterrows():
            row_sum = row['Molecular weight'] + row['Number of rings'] + row['Number of hydrogen bond donors'] + row['Number of hydrogen bond acceptors']
            diff = v1 - row_sum
            if abs(diff) <= 1:
                similar_smiles.append(row['name'])
        return similar_smiles

