#Package
library(deSolve)
library(zoo)
library("Rcpp")
library("BH")
library("tidyverse")
library(cowplot)

#SIR_new_function.R: first model, without NRTI-resistance dimension, not used anymore
#From 2005 to 2020
#comp_24_4_2_2_time_2005_2020.cpp
source("C:/Users/ahauser/Documents/Step2/Step2_revised/models/SIR_new_function.R")
x_start=x_start2005()
times=0:(15*12)
params_init= c(dtg_1st=0, #use of DTG as first-line, either 0 or 1
               dtg_switch=0, #use of DTG as switch, either 0 or 1
               p_dtg=0, #proportion of women using DTG
               dtg_eff=1, #DTG efficacy vs susceptible NNRTI
               nnrti_res=5, #Time (months) to acquire NNRTI resistance
               rev=125, #Time (months) for NNRTI-resistance reversion
               alpha1=1.97, #Impact of NNRTI resistance on NNRTI-based regimen
               alpha2=3.24, #Impact of NNRTI resistance on NNRTI-based regimen
               p_msm=0.05, #Proportion of MSM
               risk=8/3, #Increase risk of HIV transmission in MSM (e.g. due to anal intercourse)
               hiv_msm=1, #Increase in HIV prevalence in MSM
               nrti_res=40, #Time of (months) of acquiring NRTI resistance
               alpha3=1, #Impact of NRTI resistance on DTG-based regimen
               var_switch_DTG_S=1, #variation in switching rate to DTG in suppressed individuals, only used in heatmap
               var_switch_DTG_F=1, #variation in switching rate to DTG in failing individuals, only used in heatmap
               alpha4=1.62, #impact of no resistance on NNRTI-based regimen
               DTG_start_year=2020,
               switch_decrease=3.83)

SIR_solve3=function(t,y,parms){
  theta=parms[[1]]
  x_r=0
  x_i=0
  return(SIR_new3(t,y,theta,x_r,x_ir))
}
data_2020 = (ode(x_start, times=times,func=SIR_solve3,parms=list(params_init)))[,2:387]

#From 2020 to 2040
scen=c("0) all on NNRTI","1) all men on DTG","2) all men and women not of childbearing age on DTG",
       "3) all men and women not at risk of pregancy on DTG","4) all men and women on DTG")
scenario = data.frame(name1=c("No DTG",rep("DTG as a first-line",4),rep("DTG as a first-line and switch",4),rep("DTG as a first-line and switch",4)),
                      name2=c(scen,scen[2:5],scen[2:5]),
                      dtg_1st = c(0,rep(1,12)),
                      dtg_switch = c(0,rep(0,4),rep(1,8)),
                      p_dtg = c(0,rep(c(0,0.175,0.63,1),3)),
                      alpha3 = c(rep(1,9),rep(1.53,4))) #Impact of NRTI resistance on DTG-based regimen
times=180+0:240
data_scen=list()
for(i in 1:(dim(scenario)[1])){
  params=params_init
  params[c("dtg_1st","dtg_switch","p_dtg")]=as.numeric(scenario[i,c("dtg_1st","dtg_switch","p_dtg")])
  params["alpha3"] = scenario[i,"alpha3"]
  x_start=x_start2020(data_2020[181,],p_TDF = 1.0, p_DTG = scenario$p_dtg[i])
  data_scen[[i]] = (ode(x_start, times=times,func=SIR_solve3,parms=list(params)))[,2:387]
  print(i)
}

##############################################################################################################################
#Output
outcome=function(data){
  tmax=dim(data)[1]
  outcome=rowSums(data[seq(1,tmax,1),posa(c(2,3),1:4,1:2,2)])/rowSums(data[seq(1,tmax,1),posa(c(2,3),1:4,1:2,1:2)])
  return(outcome)
}

#output: no DTG, from 2005 to 2020
outcome_2020=data.frame(year=seq(2005,2020,1/12))
outcome_2020$name1=scenario[1,"name1"]
outcome_2020$name2=scenario[1,"name2"]
outcome_2020$alpha3=scenario[1,"alpha3"]
outcome_2020$outcome1=outcome(data_2020)

#output: scenarios from 2020 to 2040
for(i in 1:(dim(scenario)[1])){
  #With simulations from ub2.job
  list_data=list()
  for(j in 1:1){
    #load(paste("C:/Users/ahauser/Documents/Step2/R_results/dtg_sens/dtg_sens17/data1_",index_list[i],"_",j,".RData",sep=""))
    list_data[[j]]=data_scen[[i]]
  }
  outcome_scen_i=data.frame(year=seq(2020,2040,1/12))
  outcome_scen_i$name1=scenario[i,"name1"]
  outcome_scen_i$name2=scenario[i,"name2"]
  outcome_scen_i$alpha3=scenario[i,"alpha3"]
  outcome_scen_i$outcome1=outcome(data_scen[[i]])
  #outcome_scen[,c("out_inf","out_sup")]=t(apply(sapply(list_data,function(x){outcome(x)}),1,function(y) quantile(y,probs=c(0.025,0.975))))
  if(i==1){ outcome_scen = outcome_scen_i }else{ outcome_scen = rbind(outcome_scen,outcome_scen_i) }
}

##############################################################################################################################
##############################################################################################################################
#SIR_NRTI_res_function.R: new model with NRTI-resistance dimension
#From 2005 to 2020
#comp_24_4_2_2_time_2005_2020.cpp
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
               switch_decrease=3.83)

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
save(data_2020,file="data_2020.RData")
save(data_2020,file="../R_results/data_2020.RData")

#From 2020 to 2040
scen=c("0) all on NNRTI","1) all men on DTG","2) all men and women not of childbearing age on DTG",
       "3) all men and women not at risk of pregancy on DTG","4) all men and women on DTG","4) all men and 99% of women on DTG")
scenario = data.frame(name1=c("No DTG",rep("DTG as a first-line",4),rep("DTG as a first-line and switch",5),rep("DTG as a first-line and switch",5)),
                      name2=c(scen[1:5],scen[2:6],scen[2:6]),
                      dtg_1st = c(0,rep(1,14)),
                      dtg_switch = c(0,rep(0,4),rep(1,10)),
                      p_dtg = c(0,c(0,0.175,0.63,1),c(0,0.175,0.63,1,0.99),c(0,0.175,0.63,1,0.99)),
                      alpha3 = c(rep(1,10),rep(1.46,5))) #Impact of NRTI resistance on DTG-based regimen

times=180+0:240
data_scen=list()
data_scen_dummy=list()
for(i in 1:(dim(scenario)[1])){
  params=params_init
  params[c("dtg_1st","dtg_switch","p_dtg","alpha3")]=as.numeric(scenario[i,c("dtg_1st","dtg_switch","p_dtg","alpha3")])
  x_start=x_start2020(data_2020[181,],p_TDF = 1.0, p_DTG = scenario$p_dtg[i])
  #array(x_start[3:770],dim=c(24,4,2,2,2))
  data = (ode(c(x_start,rep(0,8)), times=times,func=SIR_solve,parms=list(params)))[,2:779] #remove time
  #array(data[241,3:770],dim=c(24,4,2,2,2))
  data_scen[[i]] = data[,1:(24*4*2*2*2+2)]
  print(outcome(data_scen[[i]]))
  data2=apply(data[,(2+24*4*2*2*2) + 1:8],2,function(x) diff(x))
  data2=cbind(rowSums(data[,posa(c(3),1,1:2,1,1:2)]),rowSums(data[,posa(c(3),2,1:2,1,1:2)]),rowSums(data[,posa(c(3),3,1:2,1,1:2)]),rowSums(data[,posa(c(3),4,1:2,1,1:2)]),
              rowSums(data[,posa(c(3),1,1:2,2,1:2)]),rowSums(data[,posa(c(3),2,1:2,2,1:2)]),rowSums(data[,posa(c(3),3,1:2,2,1:2)]),rowSums(data[,posa(c(3),4,1:2,2,1:2)]))

  data_scen_dummy[[i]] = data2
  print(i)
  #print(sum(data[240,posa(c(3),1:4,1:2,2,1:2)])/sum(data[240,posa(c(3),1:4,1:2,1:2,1:2)]))
}
save(data_2020,data_scen,data_scen_dummy,file="../R_results/scen_main.RData")

#Proportion of failure of NNRTI-based regimen after 12 and 24 months for people starting in 2030, 2035, 2040
load("../R_results/scen_main.RData")
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
for(i in c(1,2,3,4,6,7,8)){
  #for(i in 1:9){
  data=data_scen_dummy[[i]]
  data2=data_scen[[i]]
  fail=NULL
  cd4_200=NULL
  print("-------------------------------")
  print(paste("Scenario", i))
  for(year in c(2030,2035,2040)){
    x_start=array(0,dim=c(3,4,2,2))
    #x_start[1,1:4,1:2,1]=data[(year-2020)*12+1,]
    res_cd4=data[(year-2020)*12+1,5:8]/(data[(year-2020)*12+1,1:4]+data[(year-2020)*12+1,5:8])
    x_start[1,1:4,1,1] = c(10,10,10,10) * (1 - res_cd4)
    x_start[1,1:4,2,1] = c(10,10,10,10) * res_cd4
    # x_start[1,1:4,1:2,1]=data[(year-2020)*12+1,]*1000/rep(c(sum(data[(year-2020)*12+1,c(1,5)]),
    #                                                     sum(data[(year-2020)*12+1,c(2,6)]),
    #                                                     sum(data[(year-2020)*12+1,c(3,7)]),
    #                                                     sum(data[(year-2020)*12+1,c(4,8)])),2)
    cd4_200=c(cd4_200,rep(sum(x_start[1,4,1:2,1:2])/sum(x_start[1,1:4,1:2,1:2]),2))
    print(res_cd4)
    print(sum(x_start[1,1:4,2,1])/sum(x_start[1,1:4,1:2,1]))
    
    x_start=as.vector(x_start)
    data_fail = (ode(c(x_start,0,0), times=times,func=SIR_solve,parms=list(params)))[,1+1:(3*4*2*2)]
    data_12 = array(data_fail[12+1,],dim=c(3,4,2,2)) #12 months
    data_24 = array(data_fail[24+1,],dim=c(3,4,2,2)) #24 months
    fail=c(fail,c(sum(data_12[3,1:4,1:2,1:2])/sum(data_12[1:3,1:4,1:2,1:2]),
                  sum(data_24[3,1:4,1:2,1:2])/sum(data_24[1:3,1:4,1:2,1:2])))
    print(c(sum(data_12[3,1:4,1:2,1:2])/sum(data_12[1:3,1:4,1:2,1:2]),
            sum(data_24[3,1:4,1:2,1:2])/sum(data_24[1:3,1:4,1:2,1:2])))
  }
  scen_fail = rbind(scen_fail,data.frame(name1=scenario[i,"name1"],name2=scenario[i,"name2"],alpha3=scenario[i,"alpha3"],
                                         year=rep(c(2030,2035,2040),each=2),time=rep(c(12,24),3),prop_fail=fail,cd4_200=cd4_200))
  print(i)
}


##############################################################################################################################
#Output
outcome=function(data){
  tmax=dim(data)[1]
  outcome=rowSums(data[seq(1,tmax,1),posa(c(2,3),1:4,1:2,2,1:2)])/rowSums(data[seq(1,tmax,1),posa(c(2,3),1:4,1:2,1:2,1:2)])
  return(outcome)
}

#output: no DTG, from 2005 to 2020
outcome_2020=data.frame(year=seq(2005,2020,1/12))
outcome_2020$name1=scenario[1,"name1"]
outcome_2020$name2=scenario[1,"name2"]
outcome_2020$alpha3=scenario[1,"alpha3"]
outcome_2020$outcome1=outcome(data_2020)

#output: scenarios from 2020 to 2040
for(i in 1:(dim(scenario)[1])){
#for(i in 1:9){
  #With simulations from ub2.job
  list_data=list()
  list_data[[1]]=data_scen[[i]]
  
  outcome_scen_i=data.frame(year=seq(2020,2040,1/12))
  outcome_scen_i$name1=scenario[i,"name1"]
  outcome_scen_i$name2=scenario[i,"name2"]
  outcome_scen_i$alpha3=scenario[i,"alpha3"]
  outcome_scen_i$outcome1=outcome(data_scen[[i]])
  #outcome_scen[,c("out_inf","out_sup")]=t(apply(sapply(list_data,function(x){outcome(x)}),1,function(y) quantile(y,probs=c(0.025,0.975))))
  if(i==1){ outcome_scen = outcome_scen_i }else{ outcome_scen = rbind(outcome_scen,outcome_scen_i) }
}

##############################################################################################################################
#Plot
color=c("black","blue","violet","orange","red")
#DTG as first-line
dtg_plot_s1<-filter(outcome_scen,name1=="No DTG" | name1=="DTG as a first-line") %>%
  filter(alpha3 == 1) %>%
  ggplot( mapping=aes(x=year,y=outcome1,
                      color=factor(name2,labels=c("All on NNRTI  ", "Only men on DTG  ",
                                                  "All men and women not of childbearing age on DTG  ",
                                                  "All men and women not at risk of pregnancy on DTG  ",
                                                  "All men and women on DTG  ")))) +
  geom_line(size=1.1)+
  theme_bw()+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2044))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.7))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("NNRTI TDR (in %)")+
  geom_vline(xintercept = 2020, linetype="dashed", 
             color = "black", size=1)+
  annotate("text",x=2021,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  scale_color_manual(values=rep(color,2))+
  geom_line(data=outcome_2020,color="black",size=1.1)+
  geom_line(data=outcome_scen[outcome_scen$name1=="No DTG",],color="black",size=1.1)
dtg_plot_s1

#DTG as switch
dtg_plot_s2<-filter(outcome_scen,name1=="No DTG" | (name1=="DTG as a first-line and switch" & alpha3 == 1 & name2!="4) all men and 99% of women on DTG")) %>%
  ggplot( mapping=aes(x=year,y=outcome1,
                      color=factor(name2,labels=c("All on NNRTI  ", "Only men on DTG  ",
                                                  "All men and women not of childbearing age on DTG  ",
                                                  "All men and women not at risk of pregnancy on DTG  ",
                                                  "All men and women on DTG  ")))) +
  geom_line(size=1.1)+
  theme_bw()+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2044))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.7))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("NNRTI TDR (in %)")+
  geom_vline(xintercept = 2020, linetype="dashed", 
             color = "black", size=1)+
  annotate("text",x=2021,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  scale_color_manual(values=rep(color,2))+
  geom_line(data=outcome_2020,color="black",size=1.1)+
  geom_line(data=outcome_scen[outcome_scen$name1=="No DTG",],color="black",size=1.1)
dtg_plot_s2

#DTG as switch
dtg_plot_s3<-filter(outcome_scen,name1=="No DTG" | (name1=="DTG as a first-line and switch" & alpha3 == 1.46 & name2!="4) all men and 99% of women on DTG")) %>%
  ggplot( mapping=aes(x=year,y=outcome1,
                      color=factor(name2,labels=c("All on NNRTI  ", "Only men on DTG  ",
                                                  "All men and women not of childbearing age on DTG  ",
                                                  "All men and women not at risk of pregnancy on DTG  ",
                                                  "All men and women on DTG  ")))) +
  geom_line(size=1.1)+
  theme_bw()+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2044))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.7))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("NNRTI TDR (in %)")+
  geom_vline(xintercept = 2020, linetype="dashed", 
             color = "black", size=1)+
  annotate("text",x=2021,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  scale_color_manual(values=rep(color,2))+
  geom_line(data=outcome_2020,color="black",size=1.1)+
  geom_line(data=outcome_scen[outcome_scen$name1=="No DTG",],color="black",size=1.1)
dtg_plot_s3

legend<-get_legend(
  dtg_plot_s1 + theme(legend.box.margin = margin(0, 0, 0, 0))+
    theme(legend.position = "bottom",)+
    #legend.spacing.x = unit(1.0, 'cm')) +
    theme(legend.direction = "horizontal")+
    theme(legend.title=element_blank())+
    guides(colour=guide_legend(nrow=2,byrow=FALSE,
                               keywidth=1.5,))+
    theme(legend.text=element_text(size=10))
)
dtg_plot<- plot_grid(dtg_plot_s1+theme(legend.position="none"),
                     dtg_plot_s2+theme(legend.position="none"),
                     dtg_plot_s3+theme(legend.position="none"),
                     align = 'vh',
                     labels = c("A", "B","C"),
                     hjust = -1,
                     nrow = 3)
plot_grid(dtg_plot+
            theme(plot.margin = unit(c(0.5,0.5,0.5,0.5)/10, "cm")),
          legend, rel_heights = c(1, 0.08),nrow=2)#8.6x7 inches

##############################################################################################################################
#Plot
color=c("black","blue","violet","orange","red")

pd = position_dodge(0.7)
mutate(scen_fail,type_alpha3=ifelse(alpha3==1,"No impact of NRTI resistance","Impact of NRTI resistance")) %>%
ggplot(aes(x =year,y =prop_fail,color = name2)) +
  geom_point(shape = 15, size = 4, position = pd) +
  facet_grid(. ~ time)+
  #geom_errorbar(aes(ymin = sup_inf,ymax = sup_sup),width = 0.2, size = 0.7, position = pd) +
  theme_bw() +
  theme(axis.title = element_text(size=15)) +
  ylab("Failure (%)")+
  xlab("Calendar year of ART start")+
  theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0))+
  theme(axis.text.x = element_text(size=12,angle = 0))+
  theme(axis.text.y = element_text(size=12,angle = 0))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.47))+
  scale_color_manual(values=rep(color[c(1,2,3,4,5)],2))+
  labs(color="DTG-prescription strategy")+
  theme(legend.box.margin = margin(0, 0, 0, 0))+
  theme(legend.position = "bottom",)+
  theme(legend.direction = "horizontal")+
  theme(legend.title=element_blank())+
  guides(colour=guide_legend(nrow=2,byrow=FALSE,keywidth=1.5,))+
  theme(legend.text=element_text(size=11))



##############################################################################################################################
##############################################################################################################################
#Check if SIR_NRTI same as SIR_new3 (when summing fifth dimension)
x_start=abs(rnorm(768,10,5))
x_start_array=array(x_start,dim=c(24,4,2,2,2))
x_start2=c(0,0,x_start)
x_start1=c(0,0, as.vector(x_start_array[1:24,1:4,1:2,1:2,1]+x_start_array[1:24,1:4,1:2,1:2,2]))
######################################################################################

source("C:/Users/ahauser/Documents/Step2/Step2_revised/models/SIR_new_function.R")
par_sens=c(1,5,125,2.64,4.9,0.05,8/3,1)
params=list(c(0,0,0,as.numeric(par_sens)))
SIR_solve3=function(t,y,parms){
  theta=parms[[1]]
  x_r=0
  x_i=0
  return(SIR_new3(t,y,theta,x_r,x_ir))
}
data_1= array(((SIR_solve3(190,x_start1,params))[[1]])[3:(24*4*2*2+2)],dim=c(24,4,2,2))

######################################################################################

source("C:/Users/ahauser/Documents/Step2/Step2_revised/models/SIR_NRTI_res_function.R")
par_sens=c(1,5,125,2.64,4.9,0.05,8/3,1,0,1)
params=list(c(0,0,0,as.numeric(par_sens)))

SIR_solve=function(t,y,parms){
  theta=parms[[1]]
  x_r=0
  x_i=0
  return(SIR_NRTI(t,y,theta,x_r,x_ir))
}
data_2 = array(((SIR_solve(190,x_start2,params))[[1]])[3:(24*4*2*2*2+2)],dim=c(24,4,2,2,2))
data_2 = data_2[1:24,1:4,1:2,1:2,1]+data_2[1:24,1:4,1:2,1:2,2]

(data_1-data_2)>0.000001