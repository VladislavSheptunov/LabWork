import requests

class Weathet():
    #s_city = "Moscow"
    s_city = ""
    city_id = 0
    #appid = "cabdd5f5af57fadba82abb9d5f632d7c"
    appid = ""

    def __init__(self, s_city, city_id, appid):
        self.s_city = s_city
        self.appid = appid
        self.city_id = city_id

    def getWeather(self):
        try:
            res = requests.get("http://api.openweathermap.org/data/2.5/find",
                               params={'q': self.s_city, 'type': 'like', 'units': 'metric', 'APPID': self.appid})
            data = res.json()
            cities = ["{} ({})".format(d['name'], d['sys']['country'])
                      for d in data['list']]
            #print("city:", cities)
            self.city_id = data['list'][0]['id']
            #print('city_id=', city_id)
        except Exception as e:
            print("Exception (find):", e)
            pass
        try:
            res = requests.get("http://api.openweathermap.org/data/2.5/weather",
                               params={'id': self.city_id, 'units': 'metric', 'lang': 'ru', 'APPID': self.appid})
            data = res.json()
            #print("conditions:", data['weather'][0]['description'])
            #print("temp:", data['main']['temp'])
            #print("temp_min:", data['main']['temp_min'])
            #print("temp_max:", data['main']['temp_max'])
        except Exception as e:
            print("Exception (weather):", e)
            pass
        return data
