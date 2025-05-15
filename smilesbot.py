import os

from dotenv import load_dotenv
from telebot import TeleBot

import smiles_info

load_dotenv()
TG_BOT_TOKEN = os.getenv("TG_BOT_TOKEN")

bot = TeleBot(TG_BOT_TOKEN)


@bot.message_handler(commands=["start"])
def greet(message):
    print(message)
    bot.send_message(
        message.chat.id,
        "Hi, i\'m smiles bot",
    )

@bot.message_handler(commands=["work"])
def work(message):
    print(message)
    info = smiles_info.get_molecule_properties(message)
    bot.send_message(
        message.chat.id,
        f"Here some information about {message}:\nMolecular wight = {info["molecular_weight"]}\nCount of ring = {info["ring_count"]}\nCount of hydrogen bond donors = {info["h_donors"]}\nCount of hydrogen bond acceptors = {info["h_acceptors"]}\nWhether Lipinski rule = {info["lipinski_compliant"]}\n"
    )
    bot.send_photo(
        message.chat.id, 
        open(info['icon'], 'rb')
    )


@bot.message_handler(commands=["help"])
def get_help(message):
    print(message)
    bot.send_message(
        message.chat.id,
        "I can give you some information about organic compound:\n\t-Molecular wight\n\t-Count of ring\n\t-Count of hydrogen bond donors\n\t-Count of hydrogen bond acceptors\n\t-Whether Lipinski rule\n\t-picture of organic compound",
        # TODO Добавить более точное описание функционала
    )


bot.polling()
