import os
import json

from dotenv import load_dotenv
from telebot import TeleBot

import exploration
from smiles_info import get_molecule_properties

from utils.search_substructures import search_substructure
from utils.similarity_search import find_similar_mols


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


@bot.message_handler(commands=["get_similar_smiles"])
def get_similar_molecules(message):
    bot.send_message(
        message.chat.id,
        "Please, write a name of smile that you want me to find similar one to"
    )
    bot.register_next_step_handler(message, compare_with_other_mols)


def compare_with_other_mols(message):
    user_smiles = message.text
    # TODO add the choice of library
    similar = find_similar_mols(user_smiles, "data/library.csv", threshold=0.7)

    if not similar:
        bot.send_message(message.chat.id, "No similar molecules were found.")
    else:
        reply = "Similar molecules found:\n\n"
        for row in similar:
            reply += f"SMILES: {row['SMILES']}, Similarity: {row['Similarity']:.2f}\n"
        bot.send_message(message.chat.id, reply)


@bot.message_handler(commands=["get_substructures"])
def get_substructures(message):
    bot.send_message(
        message.chat.id,
        "Please, write a name of smiles that you want me to find substructures."
    )
    bot.register_next_step_handler(message, find_substructures)


def find_substructures(message):
    user_smiles = message.text
    # TODO add the choice of library
    molecules_with_found_substructures = search_substructure(user_smiles, "data/library.csv")

    if not molecules_with_found_substructures:
        bot.send_message(message.chat.id, "No molecules were found.")
    else:
        reply = "Molecules with substructures found:\n\n"
        for row in molecules_with_found_substructures:
            reply += f"{row}\n"
        bot.send_message(message.chat.id, reply)


bot.polling()
