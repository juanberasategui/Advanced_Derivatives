BSM = function(S,K,r,sigma,T) {
    d = (log(S/K) + (r + sigma^2/2)*T)/sqrt(T)/sigma
    return(S*pnorm(d) - K*exp(-r*T)*pnorm(d-sqrt(T)*sigma))
}

require('BatchGetSymbols') # install package to get time series data

first.date = '2000-01-01'
last.date = '2020-12-31'
freq.data = 'daily'; tic = c('MSFT')
PriceData = BatchGetSymbols(tickers = tic, first.date = first.date, last.date = last.date, freq = freq.data)

PriceData = PriceData$df.tickers
R = PriceData$ret.adjusted.prices[-1]

# Parameters
S0 = PriceData[length(PriceData[,1]),]$price.close; v =sd(R[(length(R)-251):length(R)])*sqrt(252)
r = 0.02; T = 0.5; K = 220

c0 = BSM(S0,K,r,v,T)

#Make bootstrap sample
n = 1000; set.seed(1); Rsample = sample(R,n,replace = TRUE)

#Scenarios for call price tomorrow and VaR in 10 days
cScenarios =BSM(S0*(1+Rsample),K,r,v,T-1/252)
VaR = quantile(c0-cScenarios,probs = 0.01)
VaR10 = VaR*sqrt(10); VaR10
