import os
import io
import json

from dotenv import load_dotenv
from telebot import TeleBot
import pandas as pd

import smiles_info


load_dotenv()
TG_BOT_TOKEN = os.getenv("TG_BOT_TOKEN")

bot = TeleBot(TG_BOT_TOKEN)


with open('text_titles.json', 'r', encoding='utf-8') as f:
    Text_for_message = json.load(f)


@bot.message_handler(commands=["start"])
def greet(message):
    print(message)
    bot.send_message(
        message.chat.id,
        Text_for_message['welcome'],
    )


@bot.message_handler(commands=["work"])
def send_smiles_information(message):
    print(message)
    info = smiles_info.get_molecule_properties(message)
    bot.send_message(
        message.chat.id,
        Text_for_message['send_info_smiles']
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
            Text_for_message["type_error"]
        )
        raise ValueError("incorrect file type")


@bot.message_handler(commands=["help"])
def get_help(message):
    print(message)
    markup = TeleBot.types.ReplyKeyboardMarkup(row_width=2)
    btn1 = TeleBot.types.KeyboardButton(Text_for_message["help_button_info"], callback_data="help_button_info")
    btn2 = TeleBot.types.KeyboardButton(Text_for_message["help_button_similar"], callback_data="help_button_similar")
    markup.add(btn1, btn2)
    bot.send_message(
        message.chat.id,
        Text_for_message["HELP_QUESTION"],
        reply_markup=markup
    )
# TODO Добавить более точное описание функционала


@bot.callback_query_handler(func=lambda call: True)
def callback_handler(call):
    if call.data == "help_button_info":
        bot.answer_callback_query(call.id, Text_for_message["info_about_info_search"])
    elif call.data == "help_button_similar":
        bot.answer_callback_query(call.id, Text_for_message["info_about_similar_search"])

bot.polling()
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
