import requests, copy, json, numpy

base_url = "https://api.iextrading.com/"
version = "1.0/"
symbol = "INTU"
scope = "3m"
resp_str = requests.get(base_url + version + "/stock/" + symbol + "/chart/" + scope).json()

dailyLogReturn = []
for i in range(1, len(resp_str)):
	dailyLogReturn.append(numpy.log(resp_str[i] / resp_str[i-1]))

print(dailyLogReturn)

#dailyVol = numpy.std(dailyLogReturn)

#print(dailyVol)
#print(type(dailyVol))