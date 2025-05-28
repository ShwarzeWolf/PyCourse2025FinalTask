import os
import io

from dotenv import load_dotenv
from telebot import TeleBot
import pandas as pd

import smiles_info

load_dotenv()
TG_BOT_TOKEN = os.getenv("TG_BOT_TOKEN")

bot = TeleBot(TG_BOT_TOKEN)

with open("text_titles.txt", "r") as file:
    data = file.readlines()

@bot.message_handler(commands=["start"])
def greet(message):
    print(message)
    bot.send_message(
        message.chat.id,
        data[0],
    )

@bot.message_handler(commands=["work"])
def work(message):
    print(message)
    info = smiles_info.get_molecule_properties(message)
    bot.send_message(
        message.chat.id,
        data[1]
    )
    bot.send_photo(
        message.chat.id, 
        open(info['icon'], 'rb')
    )

@bot.message_handler(commands=["upload"])
def handle_csv(message):
    try:
        file_info = bot.get_file(message.document.file_id)
        downloaded_file = bot.download_file(file_info.file_path)
        csv_data = io.StringIO(downloaded_file.decode('utf-8'))
        df = pd.read_csv(csv_data)
    except:
        bot.send_message(
            message.chat.id,
            data[9]
        )
        raise ValueError("incorrect file type")

@bot.message_handler(commands=["help"])
def get_help(message):
    print(message)
    markup = TeleBot.types.ReplyKeyboardMarkup(row_width=2)
    btn1 = TeleBot.types.KeyboardButton(data[6], callback_data="button1")
    btn2 = TeleBot.types.KeyboardButton(data[7], callback_data="button2")
    markup.add(btn1, btn2)
    bot.send_message(
        message.chat.id,
        data[2],
        reply_markup=markup
    )
# TODO Добавить более точное описание функционала

@bot.callback_query_handler(func=lambda call: True)
def callback_handler(call):
    if call.data == "button1":
        bot.answer_callback_query(call.id, data[3])
    elif call.data == "button2":
        bot.answer_callback_query(call.id, data[4])

bot.polling()
