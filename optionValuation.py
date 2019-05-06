import json, requests, numpy, math
from scipy.stats import norm


# Basic structure forked from https://github.com/Matt-Brigida/fin_332_python_option_class


class Option:

	def __init__(self, stockPrice, strikePrice, time, riskFree, dividendYield, optionType, symbol, volScope):

		# Set core attributes
		self.stockPrice= stockPrice
		self.strikePrice= strikePrice
		self.time = time
		self.riskFree = riskFree
		self.dividendYield = dividendYield
		
		# Interpret if it's a call or a put option
		if any(x in optionType for x in ["C", "c"]):
			self.optionType = "CALL"
		elif any(x in optionType for x in ["P", "p"]):
			self.optionType = "PUT"
		else:
			raise Exception("'{}' is not a valid optionType. Valid values start with a P or a C, case insensitive".format(optionType))

		# Find the volatility/stdDev
		resp_str = requests.get("https://api.iextrading.com/1.0/stock/" + symbol + "/chart/" + volScope).json()
		dailyLogReturn = []
		for i in range(1, len(resp_str)):
			dailyLogReturn.append(numpy.log(resp_str[i]["vwap"] / resp_str[i-1]["vwap"]))
		vol = numpy.std(dailyLogReturn)*numpy.sqrt(252)
		self.vol = vol
		self.stdDev = vol

	################################################################
	### INTRINSIC VALUE ############################################
	################################################################
	### RETURNS THE INTRINSIC VALUE FOR THE OPTION #################
	################################################################

	def intrinsic(self):
		
		if self.optionType == "CALL":
			# Call Option
			return max(self.stockPrice - self.strikePrice, 0)
		elif self.optionType == "PUT":
			# Put Option
			return max(self.strikePrice - self.stockPrice, 0)

	################################################################
	### PUT-CALL PARITY ############################################
	################################################################
	### GIVEN THE PRICE OF AN OPTION AND IF THAT PRICE IS FOR A ####
	### CALL OR PUT, IT WILL RETURN THE PRICE FOR THE OPPOSITE  ####
	################################################################

	def parity(self, optionPrice, optionType):
		
		# C + PV(x) = P + S
		
		pvx = ((self.strikePrice) / ((1 + self.riskFree) ** self.time))
		
		if optionType == "CALL":
			# Given a call, find the put
			# P = C - S + PV(x)
			return (optionPrice - self.stockPrice + pvx)
			
		elif optionType == "PUT":
			# Given a put, find the call
			# C = P + S - PV(x)
			return (optionPrice + self.stockPrice - pvx)

	################################################################
	### BLACK SCHOLES VALUATION ####################################
	################################################################
	### APPLIES THE BLACK SCHOLES FORMULA FOR A SET OF GIVENS ######
	### THE VOLATILITY CAN BE DECLARED, BUT DEFAULTS TO WHAT  ######
	### IS FOUND IN __INIT__ THROUGH IEXTRADING.ORG           ######
	################################################################

	def blackScholes(self, stdDev=None):
		
		s = self.stockPrice
		x = self.strikePrice
		r = self.riskFree
		t = self.time
		
		if stdDev == None:
			stdDev = self.stdDev
		
		d1_num = numpy.log(s/x) + (r + stdDev**2 / 2) * t
		d1_den = stdDev * numpy.sqrt(t)
		d1 = d1_num / d1_den
		
		d2 = d1 - d1_den
		
		nd1 = norm.cdf(d1)
		nd2 = norm.cdf(d2)
		
		call = s * nd1 - x * numpy.exp(-1 * r * t) * nd2
		
		if self.optionType == "CALL":
			return call
		else:
			return self.parity(call, "CALL")


	################################################################
	### MONTE CARLO VALUATION ######################################
	################################################################
	### EACH "ITERATION" IS 250K LOOPS, AVG'D TOGETHER. ############
	### THEN, EACH ITERATION AVERAGE IS AVERAGED. ##################
	################################################################

	def monteCarlo(self, iterations=32):
		s_0 = self.stockPrice
		x = self.strikePrice
		r = self.riskFree
		t = self.time
		vol = self.stdDev
		
		iterList = []
		for j in range(iterations):
			
			payoffList = []
			for i in range(250000):
				payoffList.append(max(s_0*math.exp((r-0.5*vol**2)*t+vol*math.sqrt(t)*norm.rvs())-x, 0))
			
			iterList.append(numpy.mean(payoffList) / ((1+r)**t))
			# ~ print("Iteration {} of {} ::: Mean = {}".format(j+1,iterations,iterList[-1]))
			# ~ print("{} million possibilities have been simulated so far.".format((j+1)*0.25))
		
		call = numpy.mean(iterList)
			
		if self.optionType == "CALL":
			return call
			
		else:
			return self.parity(call, "CALL")
			
			
			
	################################################################
	### GREEKS #####################################################
	################################################################
	### MAKE SURE TO CALL Option.greeks() FIRST THEN ###############
	### CALL THE GREEK AS AN ATTRIBUTE OF THE CLASS  ###############
	################################################################
	
	def greeks(self, optionType = None):
		
		if optionType == None:
			optionType = self.optionType
		
		s = self.stockPrice
		x = self.strikePrice
		t = self.time
		r = self.riskFree
		stdDev = self.stdDev
		
		d1 = (numpy.log(s/x) + (r + stdDev**2 / 2) * t) / (stdDev * numpy.sqrt(t))
		d2 = d1 - (stdDev * numpy.sqrt(t))
		n_d1 = norm.cdf(d1)
		nPrime_d1 = norm.pdf(d1)
		sigma = self.stdDev
		
		
		# DELTA
		if optionType == "CALL":
			self.delta = norm.cdf(d1)
		else:
			self.delta = norm.cdf(d1) - 1
		
		
		# GAMMA
		self.gamma = nPrime_d1 / (s*sigma*(t**(1/2)))


		# THETA
		if optionType == "CALL":
			self.theta = ((-s*nPrime_d1*sigma)/(2*(t**0.5))) - (r*x*numpy.exp(-1 * r * t)*norm.cdf(d2))
		else:
			self.theta = ((-s*nPrime_d1*sigma)/(2*(t**.05))) - (r*x*numpy.exp(-1 * r * t)*norm.cdf(-1 * d2))


		# RHO
		if optionType == "CALL":
			self.rho = x * t * numpy.exp(-1 * r * t) * norm.cdf(d2)
		else:
			self.rho = x * t * numpy.exp(-1 * r * t) * norm.cdf(-1 * d2)
		
		
		# VEGA
		self.vega = s * (t**0.5) * nPrime_d1


	################################################################
	### IMPLIED VOLATILITY FINDER ##################################
	################################################################
	### PLUG IN THE ACTUAL PRICE, AND IT WILL FIND #################
	### IMPLIED VOLATILITY THROUGH GUESS AND CHECK #################
	################################################################
	
	def impliedVol(self, actualPrice):
		
		loopContinue = True
		i = 0
		
		highVol = 1
		midVol = 0.5
		lowVol = 0
		
		while loopContinue:
			highPrice = self.blackScholes(highVol)
			midPrice = self.blackScholes(midVol)
			lowPrice = self.blackScholes(lowVol)
			
			if (actualPrice > midPrice) and (actualPrice < highPrice):
				# ~ print("midPrice too low.")
				lowVol = midVol
				midVol = (highVol + midVol) / 2
			
			elif (actualPrice < midPrice) and (actualPrice > lowPrice):
				# ~ print("midPrice too high.")
				highVol = midVol
				midVol = (lowVol + midVol) / 2
			
			elif actualPrice == midPrice:
				# ~ print("yer a wizzard.")
				return midVol
			
			else:
				# ~ print("Something didn't go right.")
				pass
				
			# ~ print("===========")
			# ~ print("highVol: {}".format(highVol))
			# ~ print("midVol: {}".format(midVol))
			# ~ print("lowVol: {}".format(lowVol))
			# ~ print("highPrice: {}".format(highPrice))
			# ~ print("midPrice: {}".format(midPrice))
			# ~ print("lowPrice: {}".format(lowPrice))
			# ~ print("actualprice: {}".format(actualPrice))
			# ~ print("Iteration: i={}".format(i))
			# ~ print("===========")
			
			if abs(actualPrice - midPrice) < 0.00000000000001:
				# ~ print("Stopping loop because difference is less than one 1-trillionth of a penny.")
				loopContinue = False
				
			if i > 10000:
				# ~ print("Stopping loop because it's looped 10,000 times")
				loopContinue = False
				
			i = i + 1
			
		return midVol
