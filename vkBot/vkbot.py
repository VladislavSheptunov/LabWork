import time
import vk_api
import random

vk = vk_api.VkApi(token='3315a90175450ab77c6139a652f1c92f52faf96ccef6b7b2bd5219d9120dba5715f1a6b8445ccf0de2233')
#vk.auth(token_only=True)
vk._auth_token()

values = {'out': 0,'count': 100,'time_offset': 60}

cmdlist = {'cmdlist':['даймузыку', 'дайвидео']}
def gen_music(strmusic):
    strmusic = strmusic + str(random.randint(19, 22))
    return strmusic

def write_msg(user_id, s):
    vk.method('messages.send', {'user_id':user_id,'message':s})

def send_music(user_id):
    vk.method('messages.send', {
        'user_id':user_id,
        'attachment':gen_music('audio-153111771_4562390')
                                }
              )
while True:
    response = vk.method('messages.get', values)
    print(response)
    if response['items']:
        values['last_message_id'] = response['items'][0]['id']
    for item in response['items']:
        if "Привет" in response['items'][0]['body']:
            write_msg(response['items'][0]['user_id'],'Привет!')
        if "Как тебя зовут?" in response['items'][0]['body']:
            write_msg(response['items'][0]['user_id'], 'СемЁн!')
        if "Как дела?" in response['items'][0]['body']:
            write_msg(response['items'][0]['user_id'], 'Прекрасно :)')
        if "даймузыку" in response['items'][0]['body']:
            send_music(response['items'][0]['user_id'])
        else:
            write_msg(response['items'][0]['user_id'], 'Список команд:\n' + cmdlist['cmdlist'][0] + '\n' + cmdlist['cmdlist'][1])
    time.sleep(1)