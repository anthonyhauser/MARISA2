#Package
library(deSolve)
library(zoo)
library("Rcpp")
library("BH")
library("tidyverse")
library("survival")
library(icenReg)
library(bbmle)
library(boot)

##############################################################################################################################
#Calibration 1: alpha1, alpha2, alpha4
#NNRTI, susc vs res: estimate impact of NNRTI resistance on treatment
#3 parameters: alpha1, alpha2 (decrease of efficacy due to NNRTI resistance compared with IeDEA estimates)
#and alpha4 (increase of efficacy due to absence of NNRTI resistance, compared with IeDEA estimates)
source("C:/Users/ahauser/Documents/Step2/Step2_revised/models/SIR_NRTI_res_function.R")
#try hazard ratio function
hazard_rct2(alpha1=2.07,alpha2=3.18,alpha4=1.61,t_limit=48,treat_pos=1)
hazard_rct2(alpha1=1.99,alpha2=3.25,alpha4=1.61,t_limit=48,treat_pos=1)
hazard_rct2(alpha1=1.97,alpha2=3.24,alpha4=1.62,t_limit=48,treat_pos=1)

#Likelihood function of hazard ratio
hazard_ml=function(alpha1,alpha2,alpha4,logsd){
  res=hazard_rct2(alpha1,alpha2,alpha4,48,1)
  res=-sum(dnorm(x=c(logit(0.88),log(3.13),log(3.13)),mean=c(logit(res[1]),log(res[2]),log(res[3])),sd=exp(logsd),log=TRUE)) 
  return(res)
}
mle=mle2(hazard_ml, start = list(alpha1=2,alpha2=2,alpha4=1,logsd=-5),  method = "Nelder-Mead")
print(mle) #2.07, 3.18, 1.61
hazard_rct2(coef(mle)[1],coef(mle)[2],coef(mle)[3],48,1)
mle_ci=confint(mle,parm=1:3)

#Likelihood function of hazard ratio
hazard_ml2=function(alpha){
  alpha1=alpha[1]
  alpha2=alpha[2]
  alpha4=alpha[3]
  res=hazard_rct2(alpha1,alpha2,alpha4,48,1)
  print(c(alpha1,alpha2,alpha4))
  res = sum((c(logit(0.88),log(3.13),log(3.13))-c(logit(res[1]),log(res[2]),log(res[3])))^2)
  return(res)
}

mle=optim(par= c(alpha1=1.97,alpha2=3.24,alpha4=1.62),fn=hazard_ml2,  method = "L-BFGS-B",lower = c(alpha1=1,alpha2=1,alpha4=1),upper=c(alpha1=5,alpha2=5,alpha4=5),control=list(reltol=1e-5))
print(mle) #1.97, 3.24, 1.62

#Sensitvity: OR between 1 and 5
hazard_rct2(1.97,3.24,1.62,48,1)
gamma=1.575
hazard_rct2(1.97*gamma,3.24*gamma,1.62,48,1) #OR of 5, after 12 months
hazard_rct2(1,1,1.62,48,1)
##############################################################################################################################
#Proportion of viral suppression among 1) susceptible, 2) 10% resistant, 3) resistant
#use SIR_ART2
source("C:/Users/ahauser/Documents/Step2/Step2_revised/models/SIR_NRTI_res_function.R")
params_init= c( art2_eff=1, #change in efficacy of second regimen
                res=5.0, #Time to acquire NNRTI resistance
                alpha1=1.97,
                alpha2=3.24,
                alpha4=1.62)
                
SIR_solve=function(t,y,parms){
  theta=parms[[1]]
  return(SIR_ART2(t,y,theta))
}
times=0:36

#susceptible
x_start=array(0,dim=c(6,4,2))
x_start[1,1:4,1]=rep(10,4)
x_start=as.vector(x_start)

params = params_init
data_art= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(6*4*2)]
data_art_t=array(data_art[length(times),],dim=c(6,4,2))
sum(data_art_t[2,1:4,1:2])/sum(data_art_t[1:3,1:4,1:2])

#10% resistant
x_start=array(0,dim=c(6,4,2))
x_start[1,1:4,1]=rep(10,4)*0.9
x_start[1,1:4,2]=rep(10,4)*0.1
x_start=as.vector(x_start)

params=params_init
data_art= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(6*4*2)]
data_art_t=array(data_art[length(times),],dim=c(6,4,2))
sum(data_art_t[2,1:4,1:2])/sum(data_art_t[1:3,1:4,1:2])

#resistant
x_start=array(0,dim=c(6,4,2))
x_start[1,1:4,2]=rep(10,4)
x_start=as.vector(x_start)

params=params_init
data_art= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(6*4*2)]
data_art_t=array(data_art[length(times),],dim=c(6,4,2))
sum(data_art_t[2,1:4,1:2])/sum(data_art_t[1:3,1:4,1:2])


#resistant
params_init= c( art2_eff=0.85, #change in efficacy of second regimen
                res=5.0, #Time to acquire NNRTI resistance
                alpha1=1.97,
                alpha2=3.24,
                alpha4=1.62)

x_start=array(0,dim=c(6,4,2))
x_start[1,1:4,1]=rep(10,4)
x_start[4,1:4,1]=rep(10,4)
x_start=as.vector(x_start)

params=params_init
data_art= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(6*4*2)]
data_art_t=array(data_art[length(times),],dim=c(6,4,2))
sum(data_art_t[2,1:4,1:2])/sum(data_art_t[1:3,1:4,1:2])
sum(data_art_t[5,1:4,1:2])/sum(data_art_t[4:6,1:4,1:2])
##############################################################################################################################
#Proportion of viral suppression among 1) susceptible, 2) 10% resistant, 3) resistant
#use SIR_ART3
source("C:/Users/ahauser/Documents/Step2/Step2_revised/models/SIR_NRTI_res_function.R")
params= c(res2=Inf, #Time to acquire NRTI resistance
          res1=5, #Time to acquire NNRTI resistance
          alpha1=1.97,
          alpha2=3.24,
          alpha4=1.62)

SIR_solve=function(t,y,parms){
  theta=parms[[1]]
  return(SIR_ART3(t,y,theta))
}
times=0:36

#susceptible
x_start=array(0,dim=c(3,4,2,2))
x_start[1,1:4,1,1]=c(10,10,10,10)
x_start=as.vector(x_start)

data_art= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(3*4*2*2)]
data_art_t=array(data_art[length(times),],dim=c(3,4,2,2))
sum(data_art_t[2,1:4,1:2,1:2])/sum(data_art_t[1:3,1:4,1:2,1:2])

#10% resistant
x_start=array(0,dim=c(3,4,2,2))
x_start[1,1:4,1,1]=c(10,10,10,10)*0.9
x_start[1,1:4,2,1]=c(10,10,10,10)*0.1
x_start=as.vector(x_start)

data_art= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(3*4*2*2)]
data_art_t=array(data_art[length(times),],dim=c(3,4,2,2))
sum(data_art_t[2,1:4,1:2,1:2])/sum(data_art_t[1:3,1:4,1:2,1:2])

#resistant
x_start=array(0,dim=c(3,4,2,2))
x_start[1,1:4,2,1]=c(10,10,10,10)
x_start=as.vector(x_start)

data_art= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(3*4*2*2)]
data_art_t=array(data_art[length(times),],dim=c(3,4,2,2))
sum(data_art_t[2,1:4,1:2,1:2])/sum(data_art_t[1:3,1:4,1:2,1:2])

##############################################################################################################################
#Calibration 2: alpha5 (called dtg_eff in the model)
#use results from NAMSAL study comparing  efficacy of DTG vs NNRTI
source("C:/Users/ahauser/Documents/Step2/Step2_revised/models/SIR_NRTI_res_function.R")

#namsal, table s8
p1 = (231+45+9)/293
p2 = (211+49+8)/279
OR_namsal = (1-p2)/p2 /((1-p1)/p1)

OR_dtg=function(nnrti_cd4,dtg_cd4,nnrti_res,alpha5){
  x_start=array(0,dim=c(6,4,2))
  x_start[1,1:4,1:2]=c((1-nnrti_res)*nnrti_cd4, nnrti_res*nnrti_cd4)
  x_start[4,1:4,1]=dtg_cd4
  x_start=as.vector(x_start)
  
  times=0:12
  params= c(art2_eff=alpha5, #change in efficacy of second regimen
            res=5.0, #Time to acquire NNRTI resistance
            alpha1=1.97,
            alpha2=3.24,
            alpha4=1.62)
  
  SIR_solve=function(t,y,parms){
    theta=parms[[1]]
    return(SIR_ART2(t,y,theta))
  }
  data_art= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(6*4*2)]
  tmax=dim(data_art)[1]
  data_tmax=array(data_art[tmax,],dim=c(6,4,2))
  res = (sum(data_tmax[5,1:4,1:2])/sum(data_tmax[6,1:4,1:2])) /  (sum(data_tmax[2,1:4,1:2])/sum(data_tmax[3,1:4,1:2])) 
  return(res)
}
namsal_eff_ml=function(alpha5,logsd,OR){
  nnrti_cd4=c(61,63,89,97) #namsal
  dtg_cd4=c(107,88,56,52) #namsal
  nnrti_res=0.03 #0.06
  
  res = OR_dtg(nnrti_cd4,dtg_cd4,nnrti_res,alpha5)
  print("-------------")
  print(res)
  res=-sum(dnorm(x=log(OR),mean=log(res),sd=exp(logsd),log=TRUE)) 
  print(alpha5)
  return(res)
}
namsal_eff_ml2=function(alpha5,OR){
  nnrti_cd4=c(61,63,89,97) #namsal
  dtg_cd4=c(107,88,56,52) #namsal
  nnrti_res=0.03 #0.06
  res = OR_dtg(nnrti_cd4,dtg_cd4,nnrti_res,alpha5)
  res=(log(OR)-log(res))^2 
  return(res)
}
dtg_eff_ml2=function(alpha5,OR){
  nnrti_cd4=rep(10,4)
  dtg_cd4=rep(10,4)
  nnrti_res=0
  res = OR_dtg(nnrti_cd4,dtg_cd4,nnrti_res,alpha5)
  res=(log(OR)-log(res))^2 
  return(res)
}

mle=mle2(namsal_eff_ml, start = list(alpha5=0.9, logsd=-3),  method = "L-BFGS-B",fixed=list(OR=1),lower = list(alpha5=0.1,logsd=-Inf),upper=list(alpha5=5,logsd=Inf))
print(mle)

mle=optim(par= c(alpha5=1),fn=namsal_eff_ml2,  method = "L-BFGS-B",OR=1.46,lower = c(alpha5=0.5),upper=c(alpha5=5))
print(mle) #0.85
OR_dtg(rep(10,4),rep(10,4),0,0.8484) #alpha5=0.85, OR=1.02
mle=optim(par= c(alpha5=1),fn=dtg_eff_ml2,  method = "L-BFGS-B",OR=2,lower = c(alpha3=0.5),upper=c(alpha3=5))
print(mle) #alpha5=1.25, OR=2
mle=optim(par= c(alpha5=1),fn=dtg_eff_ml2,  method = "L-BFGS-B",OR=5,lower = c(alpha3=0.5),upper=c(alpha3=5))
print(mle) #alpha5=2.02, OR=5

#Sensitivity analysis
mle=optim(par= c(alpha5=1),fn=dtg_eff_ml2,  method = "L-BFGS-B",OR=1,lower = c(alpha3=0.5),upper=c(alpha3=5))
print(mle) #alpha5=2.02, OR=5
######################################################################################################################################
#Calibration 3: nrti_res
#Time to acquire NRTI resistance
source("C:/Users/ahauser/Documents/Step2/Step2_revised/models/SIR_NRTI_res_function.R")
SIR_solve=function(t,y,parms){
  theta=parms[[1]]
  return(SIR_ART3(t,y,theta))
}
outcome=function(data){
  tmax=dim(data)[1]
  data_tmax=array(data[tmax,],dim=c(3,4,2,2))
  outcome = (sum(data_tmax[3,1:4,1:2,2])/sum(data_tmax[3,1:4,1:2,1:2])) 
  return(outcome)
}

#Initial values
x_start=array(0,dim=c(3,4,2,2))
x_start[1,1:4,1,1]=rep(10,4)*0.9
x_start[1,1:4,2,1]=rep(10,4)*0.1
x_start=as.vector(x_start)

times=0:36
params= c(res2=20, #Time to acquire NRTI resistance
          res1=5.0, #Time to acquire NNRTI resistance
          alpha1=1.97,
          alpha2=3.24,
          alpha4=1.62)
#Simulation
data_art= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(3*4*2*2)]
outcome(data_art)

#Likelihood function
nrti_res_ml=function(res2,logsd){
  params["res2"]=res2
  data_art= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(3*4*2*2)]
  res=outcome(data_art)
  print(res)
  res=-dnorm(x=logit(0.75*0.73),mean=logit(res),sd=exp(logsd),log=TRUE)
  return(res)
}
nrti_res_ml2=function(res2){
  params["res2"]=res2
  data_art= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(3*4*2*2)]
  res=outcome(data_art)
  res=(logit(0.75*0.73)-logit(res))^2
  return(res)
}
mle=mle2(nrti_res_ml, start = list(res2=20,logsd=-5),  method = "Nelder-Mead")
mle=optim(par= c(res2=40),fn=nrti_res_ml2,  method = "L-BFGS-B",lower = c(res2=1),upper=c(res2=100))
print(mle) #40.26

##############################################################################################################################
#Calibration 4: alpha3
#Impact of NRTI resistance on DTG-based regimen, alpha3. We assume that the OR of failing DTG when having/not having NRTI resistance is 2
source("C:/Users/ahauser/Documents/Step2/Step2_revised/models/SIR_NRTI_res_function.R")
SIR_solve=function(t,y,parms){
  theta=parms[[1]]
  return(SIR_ART3(t,y,theta))
}
outcome=function(data){
  tmax=dim(data)[1]
  data_tmax=array(data[tmax,],dim=c(3,4,2,2))
  outcome = (sum(data_tmax[3,1:4,1:2,1:2])/sum(data_tmax[1:3,1:4,1:2,1:2])) 
  return(outcome)
}

#Initial values
times=0:36
params= c(res2=Inf,
          res1=40, #Time to acquire NRTI resistance
          alpha1=1,
          alpha2=1,
          alpha4=1.62)

#susceptible
x_start=array(0,dim=c(3,4,2,2))
x_start[1,1:4,1,1]=rep(10,4)
x_start=as.vector(x_start)
data_susc= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(3*4*2*2)]

#resistant
x_start=array(0,dim=c(3,4,2,2))
x_start[1,1:4,2,1]=rep(10,4)
x_start=as.vector(x_start)
data_res= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(3*4*2*2)]

outcome(data_susc)
outcome(data_res)
(outcome(data_res)/(1-outcome(data_res))) / (outcome(data_susc)/(1-outcome(data_susc)))

#Likelihood function
nrti_res_impact_ml=function(alpha3,logsd,dtg_eff,OR){
  params[c("alpha1","alpha2")]=alpha3
  params["alpha4"]=1.62*dtg_eff
  
  x_start=array(0,dim=c(3,4,2,2))
  x_start[1,1:4,1,1]=rep(10,4)
  x_start=as.vector(x_start)
  data_susc= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(3*4*2*2)]
  data_art= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(3*4*2*2)]
  fail_susc=outcome(data_susc)
  
  x_start=array(0,dim=c(3,4,2,2))
  x_start[1,1:4,2,1]=rep(10,4)
  x_start=as.vector(x_start)
  data_res= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(3*4*2*2)]
  fail_res=outcome(data_res)
  
  res = (fail_res/(1-fail_res)) / ((fail_susc)/(1-fail_susc))
  print(res)
  res=-dnorm(x=log(OR),mean=log(res),sd=exp(logsd),log=TRUE)
  return(res)
}
nrti_res_impact_ml2=function(alpha3,dtg_eff,OR){
  params[c("alpha1","alpha2")]=alpha3
  params["alpha4"]=1.62*dtg_eff
  
  x_start=array(0,dim=c(3,4,2,2))
  x_start[1,1:4,1,1]=rep(10,4)
  x_start=as.vector(x_start)
  data_susc= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(3*4*2*2)]
  data_art= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(3*4*2*2)]
  fail_susc=outcome(data_susc)
  
  x_start=array(0,dim=c(3,4,2,2))
  x_start[1,1:4,2,1]=rep(10,4)
  x_start=as.vector(x_start)
  data_res= (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(3*4*2*2)]
  fail_res=outcome(data_res)
  
  res = (fail_res/(1-fail_res)) / ((fail_susc)/(1-fail_susc))
  #print(res)
  res=(log(OR)-log(res))^2
  return(res)
}
#old way, takes too much time
nrti_res_impact_ml(2,-5,0.85,1)
mle=mle2(nrti_res_impact_ml, start = list(alpha3=2, logsd=-3),  method = "L-BFGS-B",fixed=list(dtg_eff=0.85,OR=1),lower = list(alpha3=0.9,logsd=-Inf),upper=list(alpha3=10,logsd=Inf))
print(mle) #1

#Minimizing least square error
#dtg_eff=0.85, OR(DTG)=1.02 (NAMSAL)
mle=optim(par= c(alpha3=2),fn=nrti_res_impact_ml2,  method = "L-BFGS-B",dtg_eff=0.85,OR=1,lower = c(alpha3=0.9),upper=c(alpha3=10))
print(mle) #1
mle=optim(par= c(alpha3=2),fn=nrti_res_impact_ml2,  method = "L-BFGS-B",dtg_eff=0.85,OR=2,lower = c(alpha3=0.9),upper=c(alpha3=10))
print(mle) #1.46
mle=optim(par= c(alpha3=2),fn=nrti_res_impact_ml2,  method = "L-BFGS-B",dtg_eff=0.85,OR=5,lower = c(alpha3=0.9),upper=c(alpha3=10))
print(mle) #2.52

#dtg_eff=1.25, OR(DTG)=2
mle=optim(par= c(alpha3=2),fn=nrti_res_impact_ml2,  method = "L-BFGS-B",dtg_eff=1.25,OR=1,lower = c(alpha3=0.9),upper=c(alpha3=10))
print(mle) #1
mle=optim(par= c(alpha3=2),fn=nrti_res_impact_ml2,  method = "L-BFGS-B",dtg_eff=1.25,OR=2,lower = c(alpha3=0.9),upper=c(alpha3=10))
print(mle) #1.41
mle=optim(par= c(alpha3=2),fn=nrti_res_impact_ml2,  method = "L-BFGS-B",dtg_eff=1.25,OR=5,lower = c(alpha3=0.9),upper=c(alpha3=10))
print(mle) #2.35

#dtg_eff=2.02, OR(DTG)=5
mle=optim(par= c(alpha3=2),fn=nrti_res_impact_ml2,  method = "L-BFGS-B",dtg_eff=2.02,OR=1,lower = c(alpha3=0.9),upper=c(alpha3=10))
print(mle) #1
mle=optim(par= c(alpha3=2),fn=nrti_res_impact_ml2,  method = "L-BFGS-B",dtg_eff=2.02,OR=2,lower = c(alpha3=0.9),upper=c(alpha3=10))
print(mle) #1.38
mle=optim(par= c(alpha3=2),fn=nrti_res_impact_ml2,  method = "L-BFGS-B",dtg_eff=2.02,OR=5,lower = c(alpha3=0.9),upper=c(alpha3=10))
print(mle) #2.17



##############################################################################################################################
#Calibration 5: decrease of the switching rate to PI
#estimate so that 4% of patients are on PI on mid 2016, https://researchonline.lshtm.ac.uk/id/eprint/4646974/1/Cutting%20the%20cost%20of%20South_GOLD%20VoR.pdf

source("C:/Users/ahauser/Documents/Step2/Step2_revised/models/SIR_NRTI_res_function.R")
params_init= c(dtg_1st=0, #use of DTG as first-line, either 0 or 1
               dtg_switch=0, #use of DTG as switch, either 0 or 1
               p_dtg=0, #proportion of women using DTG
               dtg_eff=0.85, #DTG efficacy vs susceptible NNRTI
               nnrti_res=5, #Time (months) to acquire NNRTI resistance
               rev=125, #Time (months) for NNRTI-resistance reversion
               alpha1=1.97, #Impact of NNRTI resistance on NNRTI-based regimen
               alpha2=3.24, #Impact of NNRTI resistance on NNRTI-based regimen
               p_msm=0.05, #Proportion of MSM
               risk=8/3, #Increase risk of HIV transmission in MSM (e.g. due to anal intercourse)
               hiv_msm=1, #Increase in HIV prevalence in MSM
               nrti_res=40, #Time (months) of acquiring NRTI resistance
               alpha3=1, #Impact of NRTI resistance on DTG-based regimen
               var_switch_DTG_S=1, #variation in switching rate to DTG in suppressed individuals, only used in heatmap
               var_switch_DTG_F=1, #variation in switching rate to DTG in failing individuals, only used in heatmap
               alpha4=1.62, #impact of no resistance on NNRTI-based regimen
               DTG_start_year=2020,
               switch_decrease=NA)

x_start=x_start2005()
times=0:180

SIR_solve=function(t,y,parms){
  theta=parms[[1]]
  x_r=0
  x_i=0
  return(SIR_NRTI2(t,y,theta,x_r,x_ir))
}
data_2020 = (ode(c(x_start,rep(0,8)), times=times,func=SIR_solve,parms=list(params_init)))[,1+1:(24*4*2*2*2+2)]

tmax=dim(data)[1]
outcome=rowSums(data_2020[seq(1,11*12+6,1),posa(c(22:24),1:4,1:2,1:2,1:2)])/rowSums(data_2020[seq(1,11*12+6,1),posa(c(4:24),1:4,1:2,1:2,1:2)])

switch_ml = function(switch_decrease){
  params=params_init
  params["switch_decrease"] = switch_decrease
  x_start=x_start2005()
  times=0:180
  data_2020 = (ode(c(x_start,rep(0,8)), times=times,func=SIR_solve,parms=list(params)))[,1+1:(24*4*2*2*2+2)]
  res = sum(data_2020[11*12+6,posa(c(22:24),1:4,1:2,1:2,1:2)])/sum(data_2020[11*12+6,posa(c(4:24),1:4,1:2,1:2,1:2)])
  print(switch_decrease)
  print(res)
  print("-----------")
  res = (logit(res)-logit(0.04))^2
  return(res)
}

mle=optim(par= c(switch_decrease=4),fn=switch_ml,  method = "L-BFGS-B",lower = c(switch_decrease=1),upper=c(switch_decrease=20))
print(mle) #3.83