#heatmap for Ubelix array
#Heatmap, with x changing switching rate to dtg for supp or just start people

library(deSolve)
library(zoo)
library("Rcpp")
library("BH")

#Define the variation in DTG switch rate (x), and the proportion of DTG-eligible women (y)
args=(commandArgs(TRUE))
print(args)
args=as.numeric(unlist(args))
i=args[1]
j=args[2]

i_len=39
j_len=51
x=1/(0.5+(i-1)/(i_len-1)*9.5)
y=(j-1)/(j_len-1)

print(c(x,y))

#Load results in 2020 and model
#source("C:/Users/ahauser/Documents/Step2/Step2_revised/models/SIR_NRTI_res_function.R")
source("SIR_NRTI_res_function.R")
params_init= c(dtg_1st=1, #use of DTG as first-line, either 0 or 1
          dtg_switch=1, #use of DTG as switch, either 0 or 1
          p_dtg=y, #proportion of women using DTG
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
          alpha4=1.62,
          DTG_start_year=2020,
          switch_decrease=3.83)

x_start=x_start2005()
times=0:180

SIR_solve=function(t,y,parms){
  theta=parms[[1]]
  x_r=0
  x_i=0
  return(SIR_NRTI2(t,y,theta,x_r,x_ir))
}
#data_2020 = (ode(c(x_start,rep(0,8)), times=times,func=SIR_solve,parms=list(params_init)))[,1+1:(24*4*2*2*2+2)]
#data_2020 = readRDS("data_2020.rds")
load("data_2020.RData")

x_start=x_start2020(data_2020[181,],p_TDF = 1.0, p_DTG = y)

#2020
times=180+0:240
#########################################################################################################################
#Suppressed
params= c(dtg_1st=1, #use of DTG as first-line, either 0 or 1
               dtg_switch=1, #use of DTG as switch, either 0 or 1
               p_dtg=y, #proportion of women using DTG
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
               var_switch_DTG_S=x, #variation in switching rate to DTG in suppressed individuals, only used in heatmap
               var_switch_DTG_F=1, #variation in switching rate to DTG in failing individuals, only used in heatmap
               alpha4=1.62, #impact of no resistance on NNRTI-based regimen
               DTG_start_year=2020,
               switch_decrease=3.83)
#run_model
data_suppressed = (ode(c(x_start,rep(0,8)), times=times,func=SIR_solve,parms=list(params)))[,1+1:(24*4*2*2*2+2)]
res_suppressed = sum(data_suppressed[(2035-2020)*12+1,posa(c(2,3),1:4,1:2,2,1:2)]) / sum(data_suppressed[(2035-2020)*12+1,posa(c(2,3),1:4,1:2,1:2,1:2)])

#########################################################################################################################
#Failed
params= c(dtg_1st=1, #use of DTG as first-line, either 0 or 1
               dtg_switch=1, #use of DTG as switch, either 0 or 1
               p_dtg=y, #proportion of women using DTG
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
               var_switch_DTG_F=x, #variation in switching rate to DTG in failing individuals, only used in heatmap
               alpha4=1.62, #impact of no resistance on NNRTI-based regimen
               DTG_start_year=2020,
               switch_decrease=3.83)
#run_model
data_failed = (ode(c(x_start,rep(0,8)), times=times,func=SIR_solve,parms=list(params)))[,1+1:(24*4*2*2*2+2)]
res_failed = sum(data_failed[(2035-2020)*12+1,posa(c(2,3),1:4,1:2,2,1:2)]) / sum(data_failed[(2035-2020)*12+1,posa(c(2,3),1:4,1:2,1:2,1:2)])

#save data
print(c(res_suppressed,res_failed))
save(res_suppressed,res_failed,file=paste("heatmap_",i,"_",j,".RData",sep=""))