import requests
from telebot import TeleBot
import settings

bot = TeleBot(settings.TOKEN)


@bot.message_handler(commands=['start'])
def greet(message):
    print(message)
    bot.send_message(message.chat.id,
'Hi, im Dmitriis bot' )


@bot.message_handler(commands=['help'])
def help(message):
    print(message)
    bot.send_message(message.chat.id, 'Sorry, i cant help you yet')
