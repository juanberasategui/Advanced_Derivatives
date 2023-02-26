#state variables
r = 0.05; σ = 0.25; S0 = 100; T = 5; K = 165; N = 50000

#set seed
set.seed(1)

#calculate A
A = exp(-r*T) #bank account

#calculate BSN price
BSM = function(S0, r, σ, T, K){
  d1 = (log(S0/K)+(r+0.5*σ^2)*T)/(σ*sqrt(T))
  d2 = d1-σ*sqrt(T)
  price = S0*pnorm(d1)-K*exp(-r*T)*pnorm(d2)
  return(price)
}

#calculate BSM value
BSMval = BSM(S0, r, σ, T, K)


#make a loop to check that we get a 95% CI
CI_check = function(number){
  indicator_count = 0
  for (i in 1:number){
    X = rnorm(N)
    ST = S0*exp((r-0.5*σ^2)*T+σ*sqrt(T)*X)
    CT = as.vector(pmax(ST-K,0))
    BSMmc = A*mean(CT)
    sdBSMmc = sd(CT*A )/sqrt(N)
    uci = BSMmc + 1.96*sdBSMmc
    lci = BSMmc - 1.96*sdBSMmc
    indicator = ifelse(BSMval >= lci & BSMval <= uci, 1, 0)
    indicator_count = indicator_count + indicator
  }
  return(indicator_count)
}

#make a loop to check that we get a 95% CI, this time with antithetic variates
aCI_check = function(number){
  indicator_count = 0
  for (i in 1:number){
    X = rnorm(N)
    aX= as.vector(-X)
    aST = S0*exp((r-0.5*σ^2)*T+σ*sqrt(T)*aX)
    aCT = as.vector(pmax(aST-K,0))
    aBSMmc = A*mean(aCT)
    asdBSMmc = sd(aCT*A )/sqrt(N)
    auci = aBSMmc + 1.96*asdBSMmc
    alci = aBSMmc - 1.96*asdBSMmc
    indicator = ifelse(BSMval >= alci & BSMval <= auci, 1, 0)
    indicator_count = indicator_count + indicator
  }
  return(indicator_count)
}

#The same as before, however, we calculate the SE with the average of 2N simulations
avCI_check = function(number){
  indicator_count = 0
  for (i in 1:number){
    X = rnorm(N)
    aX= as.vector(-X)
    ST = S0*exp((r-0.5*σ^2)*T+σ*sqrt(T)*X)
    CT = as.vector(pmax(ST-K,0))
    BSMmc = A*mean(CT)
    aST = S0*exp((r-0.5*σ^2)*T+σ*sqrt(T)*aX)
    aCT = as.vector(pmax(aST-K,0))
    aBSMmc = A*mean(aCT)

    BSMav = (BSMmc+aBSMmc)/2
    #Calculate the antithetic variates payoff estimate
    CTav = (CT+aCT)/2
    #Calculate the antithetic variates standard error
    sdBSMav = sd(CTav*A)/sqrt(N); sdBSMav
    avuci = BSMav + 1.96*sdBSMav
    avlci = BSMav - 1.96*sdBSMav

    indicator = ifelse(BSMval >= avlci & BSMval <= avuci, 1, 0)
    indicator_count = indicator_count + indicator
  }
    return(indicator_count)
}

#check indicators
ctrBSMmc = CI_check(1000); ctrBSMmc
ctraBSMmc = aCI_check(1000); ctraBSMmc
ctravBSMmc = avCI_check(1000); ctravBSMmc

