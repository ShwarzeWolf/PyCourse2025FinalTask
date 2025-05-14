from symtable import Class
import csv
import requests
from telebot import TeleBot
import settings
from telegram import update
import pandas as pd
from rdkit.Chem import MolFromSmiles, Descriptors

class TGbot:


    def __init__(self):
        pass

    bot = TeleBot(settings.TOKEN)

    @bot.message_handler(commands=['start'])
    def greet(self, message):
        print(message)
        bot.send_message(message.chat.id,
        'Hi, im Dmitriis bot' )


    @bot.message_handler(commands=['help'])
    def help(self, message):
        print(message)
        bot.send_message(message.chat.id, 'Sorry, i cant help you yet')

    @bot.message_handler(commands=['get similar smile'])
    def get_smile_to_find_similar_one(self, message):
        bot.send_message(message.chat.id,
                         'Please, write a name of smile that you want me to find similar one to')
        users_smile = update.message.text
        return users_smile


class Dir:

    def __init__(self):
        pass


     def get_directory(self):
        data_rows = ["SMILES", "Molecular weight", "Number of rings",
                     "Number of hydrogen bond donors", "Number of hydrogen bond acceptors",
                     "Lipinski pass"]

        with open("data.csv", "w") as file:
            writer = csv.writer(file)
            writer.writerow(data_rows)


    def get_similar_smile(self):
        df = pd.read_csv(data.csv)
        row = df[df['name'] == TGbot.users_smile].iloc[0]
        v1 = row["Molecular weight"] + row["Number of rings"] + row["Number of hydrogen bond donors"] + row["Number of hydrogen bond acceptors"]
        similar_smiles = []
        for index, row in df.iterrows():
            row_sum = row["Molecular weight"] + row["Number of rings"] + row["Number of hydrogen bond donors"] + row["Number of hydrogen bond acceptors"]
            diff = v1 - row_sum
            if abs(diff) <= 1:
                similar_smiles.append(row['name'])
        return similar_smiles

