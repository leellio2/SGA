install.packages('survival')
install.packages('Hmisc')
install.packages('coxed')

library(survival)
library(Hmisc)
library(coxed)

#################################################################
# Generate data
# 	#Generic survival data
#	#Treatement A has no causal effect on outcome (truth = 0)
#	#Treatment A is associated with characteristics X 
#################################################################

n=10000
simdata <- sim.survdata(N=n, T=5000, xvars=16, censor=.8, num.data.frames = 1)

coef = matrix(c(2,2,1.5,1.5,1,1,.5,.5,.25,.25,0,0,0,0,0,0))
linear = as.matrix(simdata$xdata) %*% coef 
propensity = exp(-3+linear)/(1+exp(-3+linear))
A = rbinom(n,1, propensity)

##################### Basic PS run and weights  	##########################

	ps_mod = glm(A ~ as.matrix(simdata$xdata),family="binomial")
	ps = ps_mod$fitted.values

	#1. Null 
	null.weight = rep(1/n,n)
	null.weight

	#2. IPTW stabilized then normalized (without normalization for bubble plots)
		# not stabilized - not used: IPTW = A/ps + (1-A)/(1-ps)
	IPTW = A*sum(A)/(ps*sum(A/ps))+ (1-A)*sum((1-A))/((1-ps)*sum((1-A)/(1-ps)))
	IPTW.n = IPTW/sum(IPTW)

	#3. IPTW with trimming, stabilized  then normalized (without normalization for bubble plots)
	#version with same length as overall sample, but 0 weight for the excluded
	ind3 = (ps>0.10 & ps<0.90)
	sum(ind3)
	zt=A*(ind3==1)
	pst=ps*(ind3==1) 
	Weight = zt/pst + (1-zt)/(1-pst)
	Weight[is.na(Weight)] <- 0
	IPTWT = zt*Weight*sum(zt)/sum(zt*Weight)  +
	 (1-zt)*Weight*sum(1-zt[ind3==1])/sum((1-zt)*Weight)
	IPTWT.n = IPTWT/sum(IPTWT)

	#4. Overlap weights stabilized then normalized (without normalization for bubble plots)
		# not stabilized - not used:OW = A*(1-ps) + (1-A)*ps
	OW = A*(1-ps)*(sum(A)/sum(A*(1-ps)))+
		   (1-A)*ps* (sum(1-A)/sum((1-A)*ps))
	OW.n = OW/sum(OW)		

	#Suppose we had a categorical variable of interest
	#creating for demo
	xcat = as_tibble(cbind(B1=(simdata$xdata$X1<0),B2=(simdata$xdata$X2<0)  ))
	xcat_ignore = as_tibble(cbind(B11=(simdata$xdata$X1>0),B12=(simdata$xdata$X2>0)  ))
	xcont = simdata$xdata[,-c(1,2)] 

	#making a test dataframe
	test_df <- as_tibble(cbind(xcont, xcat, xcat_ignore))
	head(test_df)
	
	#vectors of continuous and dichotomous variables
	catnames <- names(xcat)
	contnames <- names(xcont)

