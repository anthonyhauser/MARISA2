#Package
library(deSolve)
library(zoo)
library("Rcpp")
library("BH")
library("tidyverse")
#source("C:/Users/ahauser/Documents/Step2/Step2_revised/models/SIR_NRTI_res_function.R")
source("SIR_NRTI_res_function.R")

#Load params.set and select row
load("params.set.RData")
#load("C:/Users/ahauser/Documents/Step2/Step2_revised/models/params.set.RData")
args=(commandArgs(TRUE))
print(args)
args=as.numeric(unlist(args))
print(args)
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
               switch_decrease=3.83)
if(args>1){#if args=1, params=params_init (baseline values), otherwise change values accroding to params.set
  params_init[c("nnrti_res","rev","alpha1","alpha2","p_msm","risk","hiv_msm","dtg_eff")] = params.set[args-1,]
}

#Run model, from 2005 to 2020
x_start=x_start2005()
times=0:180

SIR_solve=function(t,y,parms){
  theta=parms[[1]]
  x_r=0
  x_i=0
  return(SIR_NRTI2(t,y,theta,x_r,x_ir))
}
data_2020 = (ode(c(x_start,rep(0,8)), times=times,func=SIR_solve,parms=list(params_init)))[,1+1:(24*4*2*2*2+2)]
print(data_2020[181,])

#From 2020 to 2040
scen=c("0) all on NNRTI","1) all men on DTG","2) all men and women not of childbearing age on DTG",
       "3) all men and women not at risk of pregancy on DTG","4) all men and women on DTG","4) all men and 99% of women on DTG")
scenario = data.frame(name1=c("No DTG",rep("DTG as a first-line",5),rep("DTG as a first-line and switch",5),rep("DTG as a first-line and switch",5)),
                      name2=c(scen[1:6],scen[2:6],scen[2:6]),
                      dtg_1st = c(0,rep(1,15)),
                      dtg_switch = c(0,rep(0,5),rep(1,10)),
                      p_dtg = c(0,c(0,0.175,0.63,1,0.99),c(0,0.175,0.63,1,0.99),c(0,0.175,0.63,1,0.99)),
                      alpha3 = c(rep(1,11),rep(1.46,5))) 

times=180+0:240
data_scen=list()
data_scen_dummy=list()
for(i in 1:(dim(scenario)[1])){
  params=params_init
  params[c("dtg_1st","dtg_switch","p_dtg","alpha3")]=as.numeric(scenario[i,c("dtg_1st","dtg_switch","p_dtg","alpha3")])
  x_start=x_start2020(data_2020[181,],p_TDF = 1.0, p_DTG = scenario$p_dtg[i])
  data = (ode(c(x_start,rep(0,8)), times=times,func=SIR_solve,parms=list(params)))[,2:779] #remove time
  data_scen[[i]] = data[,1:(24*4*2*2*2+2)]
  
  data2=apply(data[,(2+24*4*2*2*2) + 1:8],2,function(x) diff(x))
  data2=cbind(rowSums(data[,posa(c(3),1,1:2,1,1:2)]),rowSums(data[,posa(c(3),2,1:2,1,1:2)]),rowSums(data[,posa(c(3),3,1:2,1,1:2)]),rowSums(data[,posa(c(3),4,1:2,1,1:2)]),
              rowSums(data[,posa(c(3),1,1:2,2,1:2)]),rowSums(data[,posa(c(3),2,1:2,2,1:2)]),rowSums(data[,posa(c(3),3,1:2,2,1:2)]),rowSums(data[,posa(c(3),4,1:2,2,1:2)]))
  data_scen_dummy[[i]] = data2
  print(i)
}

#################################################################################################################################################################
#Proportion of failure of NNRTI-based regimen after 12 and 24 months for people starting in 2030, 2035, 2040
time=0:24
params= c(res2=Inf, #Time to acquire NRTI resistance
          res1=5, #Time to acquire NNRTI resistance
          alpha1=1.97,
          alpha2=3.24,
          alpha4=1.62)

SIR_solve=function(t,y,parms){
  theta=parms[[1]]
  return(SIR_ART3(t,y,theta))
}
scen_fail=NULL
for(i in c(1,2,3,4,6,7,8,9,11,12,13,14,16)){
  #for(i in 1:9){
  data=data_scen_dummy[[i]]
  data2=data_scen[[i]]
  fail=NULL
  cd4_200=NULL
  for(year in c(2025,2030,2035,2040)){
    #xstart
    x_start=array(0,dim=c(3,4,2,2))
    res_cd4=data[(year-2020)*12+1,5:8]/(data[(year-2020)*12+1,1:4]+data[(year-2020)*12+1,5:8])
    x_start[1,1:4,1,1] = c(10,10,10,10) * (1 - res_cd4)
    x_start[1,1:4,2,1] = c(10,10,10,10) * res_cd4
    cd4_200=c(cd4_200,rep(sum(x_start[1,4,1:2,1:2])/sum(x_start[1,1:4,1:2,1:2]),2))
    x_start=as.vector(x_start)
    #simulation
    data_fail = (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(3*4*2*2)]
    data_12 = array(data_fail[12+1,],dim=c(3,4,2,2)) #12 months
    data_24 = array(data_fail[24+1,],dim=c(3,4,2,2)) #24 months
    fail=c(fail,c(sum(data_12[3,1:4,1:2,1:2])/sum(data_12[1:3,1:4,1:2,1:2]),
                  sum(data_24[3,1:4,1:2,1:2])/sum(data_24[1:3,1:4,1:2,1:2])))
  }
  scen_fail = rbind(scen_fail,data.frame(name1=scenario[i,"name1"],
                                         name2=scenario[i,"name2"],
                                         alpha3=scenario[i,"alpha3"],
                                         year=rep(c(2025,2030,2035,2040),each=2),
                                         time=rep(c(12,24),4),
                                         prop_fail=fail,
                                         cd4_200=cd4_200))
  print(scenario[i,"name2"])
  print(i)
}
#Save
save(data_2020,data_scen,scen_fail,file=paste("sens_",args,".RData",sep=""))