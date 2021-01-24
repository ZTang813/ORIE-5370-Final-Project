import datetime
import numpy as np
import pandas as pd

# specify timeperiod for computing RSI here
timeperiod = 30

# Load data
DF = pd.read_csv("output.csv")

# Stock names
stocks = list(DF.columns)[2:]
days = len(DF["Date"])
RSI = pd.DataFrame()
RSI["Date"] = DF["Date"]

for stock in stocks:
    prices = DF[stock]
    firstday = 0
    if np.isnan(prices[firstday]):
        print("NaN Detected for", stock)
        RSI[stock] = [0]*days
        continue

    diff = list(prices[firstday:].diff()[1:])
    diff.insert(0, 0)
    A = 0
    B = 0
    rsi=[0]*(timeperiod+1+firstday)
    aA=[0]*(timeperiod+1+firstday)
    aB=[0]*(timeperiod+1+firstday)
    for i in range (1,timeperiod+firstday):
        if diff[i]>=0:
            A=A+diff[i]
        else:
            B=B+diff[i]
    rsi[-1]=A/(A-B)*100  #计算得到第一个RSI
    aA[-1]=A/timeperiod
    aB[-1]=B/timeperiod


    for i in range(timeperiod+1+firstday, days):
        if diff[i]>=0:
            aA.append((aA[-1]*(timeperiod-1)+ diff[i])/timeperiod)
            aB.append(aB[-1]*(timeperiod-1)/timeperiod)
        else:
            aA.append(aA[-1]*(timeperiod-1)/timeperiod)
            aB.append((aB[-1]*(timeperiod-1)+ diff[i])/timeperiod)
        rsi.append(aA[-1]/(aA[-1]-aB[-1])*100)


    RSI[stock] = rsi

RSI.to_csv("stock_RSI.csv")
