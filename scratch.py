import optionValuation as ov

# July 21 2019, 220 Call = 2.50 (Put ~ 10.85)
# ~ aaplOption = ov.Option(211.95, 1, 220, (48/365), 0.024, 0.015, "c")

aaplOption = ov.Option(208.45, 220, (46/365), 0.024, 0.015, "p", "aapl", "3m")
# ~ (self, symbol, stockPrice, strikePrice, time, riskFree, dividendYield, optionType, volScope):


print(aaplOption.impliedVol(13.47))

