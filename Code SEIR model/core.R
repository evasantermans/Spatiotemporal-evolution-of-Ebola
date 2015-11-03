########################################################
# Modelling the Ebola outbreak in West Africa 
# SEIR model & MCMC estimation
########################################################
# Author: Eva Santermans
# Last update: July 17, 2015
########################################################


####################
## CORE   
####################

###### ODE with parameters: 
### - Overdispersion parameters phi1, phi2 (0-1) [1-2] 
### - E(time=0) (2) [3]
### - proportion fatal cases phi (3) [4]
### - underreporting rate cases and deaths rho (4) [5]
### - Reproduction 1st time period R0 (5) [6]
### - Change in reproduction number ri (6-...) [7-...]

###### Fixed parameters: 
### - preinfperiod = 9.4 days
### - infperiod non-fatal cases = 16.4 days
### - infperiod fatal cases = 7.5 days


########################################################################################
## LIBRARIES

library(Rcpp)
library(inline)
library(LaplacesDemon)


## Piece-wise constant transmission beta
##################################################

Rpwc = function(pars,ndates){
  n.breaks = length(ndates)
  Rvec = pars[6:length(pars)]
  R = NULL
  R = c(R,rep(Rvec[1],ndates[1]))
  for (i in 2:length(Rvec)){
    R = c(R,rep(sum(Rvec[1:i]),ndates[i]))
  }
  beta = R/((1-pars[4])*16.4+pars[4]*7.5)   
  return(beta)
}


## ODE solver (Euler's method)
##################################################

odefunc <- '
Environment myEnv = Environment::global_env();
Function Rpwc = myEnv["Rpwc"]; 

// 1. Define Variables 
// 1.1. Input Variables from R;
Rcpp::NumericVector pars(x);
Rcpp::NumericVector co(a);
Rcpp::NumericVector time(t);
Rcpp::NumericVector ndates(nsel);

// 1.2. Parameters;
double popsize = co(0);
double precision = co(1);
double preinfperiod = 9.4;
double prop = 1-pars(3);
double infperiodr = 16.4;
double infperiodd = 7.5;
double sigmaparm = 1/preinfperiod;
double gammaparm = 1/infperiodr;
double alphaparm = 1/infperiodd;

Rcpp::NumericVector beta(time.size());
beta = Rpwc(pars,ndates);


// 1.3. State variables;
Rcpp::NumericVector N(time.size());
Rcpp::NumericVector S(time.size());
Rcpp::NumericVector E(time.size());
Rcpp::NumericVector Ir(time.size());
Rcpp::NumericVector Id(time.size());
Rcpp::NumericVector Inew(time.size());
Rcpp::NumericVector M(time.size());
Rcpp::NumericVector R(time.size());
N(0) = popsize;
S(0) = popsize-pars(2);
E(0) = pars(2);
Ir(0) = 0;
Id(0) = 0;
Inew(0) = 0;
M(0) = 0;
R(0) = 0;


// 2. System of ODEs
for (int i = 0; i < (time.size()-1); i++){
S(i+1) = S(i) - beta(i)*S(i)*(Ir(i)+Id(i))/N(i)*precision;
E(i+1) = E(i) + beta(i)*S(i)*(Ir(i)+Id(i))/N(i)*precision - sigmaparm*E(i)*precision;
Ir(i+1) = Ir(i) + prop*sigmaparm*E(i)*precision - gammaparm*Ir(i)*precision;
Id(i+1) = Id(i) + (1-prop)*sigmaparm*E(i)*precision - alphaparm*Id(i)*precision;
Inew(i+1) = Inew(i) + sigmaparm*E(i)*precision;
M(i+1) = M(i) + alphaparm*Id(i)*precision;
R(i+1) = R(i) + gammaparm*Ir(i)*precision;
N(i+1) = S(i+1) + E(i+1) + Ir(i+1) + Id(i) + R(i+1); 
}

return Rcpp::List::create(Rcpp::Named("time") = time,
Rcpp::Named("S") = S,
Rcpp::Named("E") = E,
Rcpp::Named("Ir") = Ir,
Rcpp::Named("Id") = Id,
Rcpp::Named("Inew") = Inew,
Rcpp::Named("M") = M,
Rcpp::Named("R") = R,
Rcpp::Named("beta") = beta);
'


Rcpp.odefunc.pwc = cxxfunction(signature(x="numeric",a="numeric",t="numeric",nsel="numeric"),
                               plugin="Rcpp",body=odefunc)


ODEfuncpwc = function(pars,n.cutoffs){
  a = c(popsize,precision)
  t = time
  res = Rcpp.odefunc.pwc(pars,a,t,n.cutoffs)
  return(res)
}


## Loglikelihood calculation
##################################################

Rlikelihood = function(pars,Inew,M){
  res1=diff(c(0,Inew[index1.cum==1]))
  res2=diff(c(0,M[index2.cum==1]))
  ll1=sum(log(1e-50+dnbinom(obs1,mu=as.vector(res1)*pars[5],size=pars[1])))
  ll2=sum(log(1e-50+dnbinom(obs2,mu=as.vector(res2)*pars[5],size=pars[2])))
  return(list(ll=ll1+ll2,res1=res1*pars[5],res2=res2*pars[5]))
}


## MCMC with LaplacesDemon
##################################################


Model = function(pars,Data){
  
  ## parameters
  overdisp = exp(pars[1:2])
  E0 = exp(pars[3])
  prop = exp(pars[4])/(exp(pars[4])+1)
  underreport = exp(pars[5])/(exp(pars[5])+1)
  R0 = exp(pars[6])
  r = pars[7:(5+n.breaks)]
  
  beta = c(overdisp,E0,prop,underreport,R0,r)
  
  ## log(prior densities)
  overdisp.prior = dhalfcauchy(overdisp,25,log=T)
  E0.prior = dunif(E0,0,1,log=T)
  underreport.prior = log(dtrunc(underreport,"norm",a=0,b=1,mean=1/3,sd=0.1))
  prop.prior = dbeta(prop,10,10,log=T)
  R0.prior = dunif(R0,0,10,log=T)
  r.prior = dunif(r,-2,2,log=T)
  
  beta.prior = sum(overdisp.prior) + E0.prior + prop.prior + underreport.prior + R0.prior + sum(r.prior)
  
  ## loglikelihood
  ode = ODEfuncpwc(beta,n.cutoffs)
  LL = Rlikelihood(beta,ode$Inew,ode$M)
  
  ## logposterior
  LP = LL$ll + sum(beta.prior)
  
  ## output
  modelout = list(LP=LP, Dev=-2*LL$ll, Monitor=c(LP,overdisp,E0,prop,underreport,R0,r),yhat=c(LL$res1,LL$res2),parm=pars)
  return(modelout)
}

