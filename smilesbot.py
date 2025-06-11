import os
import json

from dotenv import load_dotenv
from telebot import TeleBot

import exploration
from smiles_info import get_molecule_properties



load_dotenv()
TG_BOT_TOKEN = os.getenv("TG_BOT_TOKEN")


bot = TeleBot(TG_BOT_TOKEN)


with open("text_titles.json", "r", encoding="utf-8") as f:
    TEXT_FOR_MESSAGE = json.load(f)


@bot.message_handler(commands=["start"])
def greet(message):
    print(message)
    bot.send_message(
        message.chat.id,
        TEXT_FOR_MESSAGE["welcome"],
    )


@bot.message_handler(commands=["smile info"])
def send_smiles_information(message):
    print(message)
    info = get_molecule_properties(message)
    bot.send_message(
        message.chat.id,
        TEXT_FOR_MESSAGE["send_info_smiles"]
    )
    bot.send_photo(
        message.chat.id, 
        open(info["icon"], "rb")
    )


@bot.message_handler(commands=["upload"])
def handle_csv(message):
    try:
        file_info = bot.get_file(message.document.file_id)
        downloaded_file = bot.download_file(file_info.file_path)
        try:
            exploration.check_file_corections(downloaded_file)
            bot.send_message(
                message.chat.id,
                TEXT_FOR_MESSAGE["file_correct"]
            )
            not_valid_molecules_df = exploration.cheak_content_corections(downloaded_file)
            if not_valid_molecules_df.empty:
                bot.send_message(
                    message.chat.id,
                    TEXT_FOR_MESSAGE["file_correct_content"]
                )
            else:
                bot.send_message(
                    message.chat.id,
                    TEXT_FOR_MESSAGE["no_parsing_error"]
                )
        except:
            bot.send_message(
                message.chat.id,
                TEXT_FOR_MESSAGE["content_error"]
            )
            raise ValueError("incorrect content in file")
    except:
        bot.send_message(
            message.chat.id,
            TEXT_FOR_MESSAGE["type_error"]
        )
        raise ValueError("incorrect file type")


@bot.message_handler(commands=["help"])
def get_help(message):
    print(message)
    markup = TeleBot.types.ReplyKeyboardMarkup(row_width=2)
    btn1 = TeleBot.types.KeyboardButton(
        TEXT_FOR_MESSAGE["help_button_info"],
        callback_data="help_button_info"
    )
    btn2 = TeleBot.types.KeyboardButton(
        TEXT_FOR_MESSAGE["help_button_similar"],
        callback_data="help_button_similar"
    )
    markup.add(btn1, btn2)
    bot.send_message(
        message.chat.id,
        TEXT_FOR_MESSAGE["HELP_QUESTION"],
        reply_markup=markup
    )
# TODO Добавить более точное описание функционала


@bot.callback_query_handler(func=lambda call: True)
def callback_handler(call):
    if call.data == "help_button_info":
        bot.answer_callback_query(call.id, TEXT_FOR_MESSAGE["info_about_info_search"])
    elif call.data == "help_button_similar":
        bot.answer_callback_query(call.id, TEXT_FOR_MESSAGE["info_about_similar_search"])


bot.polling()
