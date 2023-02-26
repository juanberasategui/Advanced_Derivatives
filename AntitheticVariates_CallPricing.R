# ECO423 Project 3 solution in R

# Step 1
r=0.05; S0=100; sigma=0.25; K=165; T=5; N=50000

# Step 2
d1 = (log(S0/K)+(r+sigma^2/2)*T) / sigma/sqrt(T)
BSM = S0*pnorm(d1) - exp(-r*T)*K*pnorm(d1-sigma*sqrt(T)); BSM

# Step 3
set.seed(1)
X = rnorm(N)

# Step 4
ST = S0*exp((r-sigma^2/2)*T + sigma*sqrt(T)*X)

# Step 5. Note the use of 'pmax' rather than 'max'
CT = pmax(ST-K,0)

# Step 6
BSMmc = exp(-r*T)*mean(CT); BSMmc

# Step 7
sdBSMmc = sd(CT*exp(-r*T))/sqrt(N); sdBSMmc

# Step 8
uci = BSMmc + qnorm(1-0.05/2)*sdBSMmc
lci = BSMmc + qnorm(0.05/2)*sdBSMmc

# Step 9
if(BSM < uci & BSM > lci) 1 else 0

# Step 10
aX = -X

# Step 11
aST = S0*exp((r-sigma^2/2)*T + sigma*sqrt(T)*aX)
aCT = pmax(aST-K,0)
aBSMmc = exp(-r*T)*mean(aCT); aBSMmc
asdBSMmc = sd(aCT*exp(-r*T))/sqrt(N); asdBSMmc
auci = aBSMmc + qnorm(1-0.05/2)*asdBSMmc
alci = aBSMmc + qnorm(0.05/2)*asdBSMmc
if(BSM < auci & BSM > alci) 1 else 0

# Step 12
BSMav = (BSMmc+aBSMmc)/2

# Step 13
CTav = (CT+aCT)/2

# Step 14
sdBSMav = sd(CTav*exp(-r*T))/sqrt(N); sdBSMav
uciAv = BSMav + qnorm(1-0.05/2)*sdBSMav; uciAv
lciAv = BSMav + qnorm(0.05/2)*sdBSMav; lciAv
if(BSM < uciAv & BSM > lciAv) 1 else 0

# The % improvement in SEM from antithetic variates technique:
(sdBSMmc/sqrt(2)-sdBSMav) / (sdBSMmc/sqrt(2))

# Step 15
ctrBSMmc = 0; ctraBSMmc = 0; ctrBSMav = 0

for(i in 1:100) {
	X = rnorm(N)
	ST = S0*exp((r-sigma^2/2)*T + sigma*sqrt(T)*X)
	CT = pmax(ST-K,0)
	BSMmc = exp(-r*T)*mean(CT)
	sdBSMmc = sd(CT*exp(-r*T))/sqrt(N)
	uci = BSMmc + qnorm(1-0.05/2)*sdBSMmc
	lci = BSMmc + qnorm(0.05/2)*sdBSMmc
	if(BSM < uci & BSM > lci) {
		ctrBSMmc = ctrBSMmc + 1
	}
	
	aX = -X
	aST = S0*exp((r-sigma^2/2)*T + sigma*sqrt(T)*aX)
	aCT = pmax(aST-K,0)
	aBSMmc = exp(-r*T)*mean(aCT)
	asdBSMmc = sd(aCT*exp(-r*T))/sqrt(N)
	auci = aBSMmc + qnorm(1-0.05/2)*asdBSMmc
	alci = aBSMmc + qnorm(0.05/2)*asdBSMmc
	if(BSM < auci & BSM > alci) {
		ctraBSMmc = ctraBSMmc + 1
	}

	CTav = (CT+aCT)/2
	BSMav = exp(-r*T)*mean(CTav) # = (BSMmc+aBSMmc)/2
	sdBSMav = sd(CTav*exp(-r*T))/sqrt(N)
	uciAv = BSMav + qnorm(1-0.05/2)*sdBSMav
	lciAv = BSMav + qnorm(0.05/2)*sdBSMav
	if(BSM < uciAv & BSM > lciAv) {
		ctrBSMav = ctrBSMav + 1
	}
}

ctrBSMmc; ctraBSMmc; ctrBSMav
