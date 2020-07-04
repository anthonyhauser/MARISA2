library(deSolve)
library("tidyverse")
library(cowplot)
source("C:/Users/ahauser/Documents/Step2/Step2_revised/models/SIR_NRTI_res_function.R")

#Main scenarios
#load
load("../R_results/scen_main.RData")
scen=c("0) all on NNRTI","1) all men on DTG","2) all men and women not of childbearing age on DTG",
       "3) all men and women not at risk of pregancy on DTG","4) all men and women on DTG","4) all men and 99% of women on DTG")
scenario = data.frame(name1=c("No DTG",rep("DTG as a first-line",4),rep("DTG as a first-line and switch",5),rep("DTG as a first-line and switch",5)),
                      name2=c(scen[1:5],scen[2:6],scen[2:6]),
                      dtg_1st = c(0,rep(1,14)),
                      dtg_switch = c(0,rep(0,4),rep(1,10)),
                      p_dtg = c(0,c(0,0.175,0.63,1),c(0,0.175,0.63,1,0.99),c(0,0.175,0.63,1,0.99)),
                      alpha3 = c(rep(1,10),rep(1.46,5))) #Impact of NRTI resistance on DTG-based regimen
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
dtg_plot_main1<-filter(outcome_scen,name1=="No DTG" | name1=="DTG as a first-line") %>%
  filter(alpha3 == 1) %>%
  ggplot( mapping=aes(x=year,y=outcome1,
                      color=factor(name2,labels=c("All on NNRTI  ", "Only men on DTG  ",
                                                  "All men and women not of childbearing age on DTG  ",
                                                  "All men and women not at risk of pregancy on DTG  ",
                                                  "All men and women on DTG  ")))) +
  geom_line(size=1.1)+
  theme_bw()+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2040))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.6))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("NNRTI TDR (in %)")+
  geom_vline(xintercept = 2020, linetype="dashed", 
             color = "black", size=0.8)+
  annotate("text",x=2021,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  scale_color_manual(values=rep(color,2))+
  geom_line(data=outcome_2020,color="black",size=1.1)+
  geom_line(data=outcome_scen[outcome_scen$name1=="No DTG",],color="black",size=1.1)

#DTG as switch
dtg_plot_main2<-filter(outcome_scen,name1=="No DTG" | (name1=="DTG as a first-line and switch" & alpha3 == 1 & name2!="4) all men and 99% of women on DTG")) %>%
  ggplot( mapping=aes(x=year,y=outcome1,
                      color=factor(name2,labels=c("All on NNRTI  ", "Only men on DTG  ",
                                                  "All men and women not of childbearing age on DTG  ",
                                                  "All men and women not at risk of pregancy on DTG  ",
                                                  "All men and women on DTG  ")))) +
  geom_line(size=1.1)+
  theme_bw()+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2040))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.6))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("NNRTI TDR (in %)")+
  geom_vline(xintercept = 2020, linetype="dashed", 
             color = "black", size=0.8)+
  annotate("text",x=2021,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  scale_color_manual(values=rep(color,2))+
  geom_line(data=outcome_2020,color="black",size=1.1)+
  geom_line(data=outcome_scen[outcome_scen$name1=="No DTG",],color="black",size=1.1)

#DTG as switch
dtg_plot_main3<-filter(outcome_scen,name1=="No DTG" | (name1=="DTG as a first-line and switch" & alpha3 == 1.46 & name2!="4) all men and 99% of women on DTG")) %>%
  ggplot( mapping=aes(x=year,y=outcome1,
                      color=factor(name2,labels=c("All on NNRTI  ", "Only men on DTG  ",
                                                  "All men and women not of childbearing age on DTG  ",
                                                  "All men and women not at risk of pregnancy on DTG  ",
                                                  "All men and women on DTG  ")))) +
  geom_line(size=1.1)+
  theme_bw()+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2040))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.6))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("NNRTI TDR (in %)")+
  geom_vline(xintercept = 2020, linetype="dashed", 
             color = "black", size=0.8)+
  annotate("text",x=2021,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  scale_color_manual(values=rep(color,2))+
  geom_line(data=outcome_2020,color="black",size=1.1)+
  geom_line(data=outcome_scen[outcome_scen$name1=="No DTG",],color="black",size=1.1)

legend_main<-get_legend(
  dtg_plot_main1 + theme(legend.box.margin = margin(0, 0, 0, 0))+
    theme(legend.position = "bottom",)+
    #legend.spacing.x = unit(1.0, 'cm')) +
    theme(legend.direction = "horizontal")+
    theme(legend.title=element_blank())+
    guides(colour=guide_legend(nrow=2,byrow=FALSE,
                               keywidth=1.5,))+
    theme(legend.text=element_text(size=10))
)
dtg_plot_main<- plot_grid(dtg_plot_main1+theme(legend.position="none")+theme(plot.margin = unit(c(0,0,0,10)/10, "cm")),
                     dtg_plot_main2+theme(legend.position="none"),
                     dtg_plot_main3+theme(legend.position="none"),
                     align = 'vh',
                     labels = c("A", "B","C"),
                     hjust = -1,
                     nrow = 3)
##############################################################################################################################
#No Treat-All policy
#From 2005 to 2020
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

SIR_solve=function(t,y,parms){
  theta=parms[[1]]
  x_r=0
  x_i=0
  return(SIR_NRTI2(t,y,theta,x_r,x_ir))
}


source("C:/Users/ahauser/Documents/Step2/Step2_revised/models/SIR_NRTI_res_function.R")

SIR_solve=function(t,y,parms){
  theta=parms[[1]]
  treat_all=0
  treat_interruption=0
  return(SIR_NRTI2_scenario(t,y,theta,treat_all,treat_interruption))
}

x_start=x_start2005()
times=0:(15*12)
a <- SIR_solve(1,c(x_start,rep(0,8)),list(params_init))

data_2020_notreatall = (ode(c(x_start,rep(0,8)), times=times,func=SIR_solve,parms=list(params_init)))[,1+1:(24*4*2*2*2+2)]
outcome(data_2020_notreatall)

#From 2020 to 2040
scen=c("0) all on NNRTI","1) all men on DTG","2) all men and women not of childbearing age on DTG",
       "3) all men and women not at risk of pregancy on DTG","4) all men and women on DTG","4) all men and 99% of women on DTG")
scenario = data.frame(name1=c("No DTG",rep("DTG as a first-line",4),rep("DTG as a first-line and switch",5),rep("DTG as a first-line and switch",5)),
                      name2=c(scen[1:5],scen[2:6],scen[2:6]),
                      dtg_1st = c(0,rep(1,14)),
                      dtg_switch = c(0,rep(0,4),rep(1,10)),
                      p_dtg = c(0,c(0,0.175,0.63,1),c(0,0.175,0.63,1,0.99),c(0,0.175,0.63,1,0.99)),
                      alpha3 = c(rep(1,10),rep(1.46,5))) #Impact of NRTI resistance on DTG-based regimen

times=181+0:240
data_scen_notreatall=list()
for(i in 1:(dim(scenario)[1])){
  params=params_init
  params[c("dtg_1st","dtg_switch","p_dtg","alpha3")]=as.numeric(scenario[i,c("dtg_1st","dtg_switch","p_dtg","alpha3")])
  x_start=x_start2020(data_2020_notreatall[181,],p_TDF = 1.0, p_DTG = scenario$p_dtg[i])
  data = (ode(c(x_start,rep(0,8)), times=times,func=SIR_solve,parms=list(params)))[,2:779] #remove time
  data_scen_notreatall[[i]] = data[,1:(24*4*2*2*2+2)]
  print(i)
  
  print(outcome(data_scen_notreatall[[i]]))
}

save(data_2020_notreatall,data_scen_notreatall,file="../R_results/scen_notreatall.RData")
##############################################################################################################################
#Treatment interruption
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


SIR_solve=function(t,y,parms){
  theta=parms[[1]]
  treat_all=1
  treat_interruption=1
  return(SIR_NRTI2_scenario(t,y,theta,treat_all,treat_interruption))
}
x_start=x_start2005()
times=0:(15*12)
data_2020_treatinterruption = (ode(c(x_start,rep(0,8)), times=times,func=SIR_solve,parms=list(params_init)))[,1+1:(24*4*2*2*2+2)]

#From 2020 to 2040
scen=c("0) all on NNRTI","1) all men on DTG","2) all men and women not of childbearing age on DTG",
       "3) all men and women not at risk of pregancy on DTG","4) all men and women on DTG","4) all men and 99% of women on DTG")
scenario = data.frame(name1=c("No DTG",rep("DTG as a first-line",4),rep("DTG as a first-line and switch",5),rep("DTG as a first-line and switch",5)),
                      name2=c(scen[1:5],scen[2:6],scen[2:6]),
                      dtg_1st = c(0,rep(1,14)),
                      dtg_switch = c(0,rep(0,4),rep(1,10)),
                      p_dtg = c(0,c(0,0.175,0.63,1),c(0,0.175,0.63,1,0.99),c(0,0.175,0.63,1,0.99)),
                      alpha3 = c(rep(1,10),rep(1.46,5))) #Impact of NRTI resistance on DTG-based regimen

times=181+0:240
data_scen_treatinterruption=list()
for(i in 1:(dim(scenario)[1])){
  params=params_init
  params[c("dtg_1st","dtg_switch","p_dtg","alpha3")]=as.numeric(scenario[i,c("dtg_1st","dtg_switch","p_dtg","alpha3")])
  x_start=x_start2020(data_2020_treatinterruption[181,],p_TDF = 1.0, p_DTG = scenario$p_dtg[i])
  data = (ode(c(x_start,rep(0,8)), times=times,func=SIR_solve,parms=list(params)))[,2:779] #remove time
  data_scen_treatinterruption[[i]] = data[,1:(24*4*2*2*2+2)]
  print(i)
}
save(data_2020_treatinterruption,data_scen_treatinterruption,file="../R_results/scen_treatinterruption.RData")
##############################################################################################################################
#Output for no Treat-All and treatment interruption
#data
load("../R_results/scen_notreatall.RData")
data_2020=data_2020_notreatall
data_scen=data_scen_notreatall
title_scen="No Treat-All policy"
filename="../figures/sens_notreatall.pdf"

load("../R_results/scen_treatinterruption.RData")
data_2020=data_2020_treatinterruption
data_scen=data_scen_treatinterruption
title_scen="Treatment interruption"
filename="../figures/sens_treatinterruption.pdf"

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
  outcome_scen_i=data.frame(year=seq(2020,2040,1/12))
  outcome_scen_i$name1=scenario[i,"name1"]
  outcome_scen_i$name2=scenario[i,"name2"]
  outcome_scen_i$alpha3=scenario[i,"alpha3"]
  outcome_scen_i$outcome1=outcome(data_scen[[i]])
  if(i==1){ outcome_scen = outcome_scen_i }else{ outcome_scen = rbind(outcome_scen,outcome_scen_i) }
}


##############################################################################################################################
#Plot for no Treat-All and treatment interruption
color=c("black","blue","violet","orange","red")
#DTG as first-line
dtg_plot_scen1<-filter(outcome_scen,name1=="No DTG" | name1=="DTG as a first-line") %>%
  filter(alpha3 == 1) %>%
  ggplot( mapping=aes(x=year,y=outcome1,
                      color=factor(name2,labels=c("All on NNRTI  ", "Only men on DTG  ",
                                                  "All men and women not of childbearing age on DTG  ",
                                                  "All men and women not of childbaring potential on DTG  ",
                                                  "All men and women on DTG  ")))) +
  geom_line(size=1.1)+
  theme_bw()+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2040))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.6))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("NNRTI TDR (in %)")+
  geom_vline(xintercept = 2020, linetype="dashed", 
             color = "black", size=0.8)+
  annotate("text",x=2021,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  scale_color_manual(values=rep(color,2))+
  geom_line(data=outcome_2020,color="black",size=1.1)+
  geom_line(data=outcome_scen[outcome_scen$name1=="No DTG",],color="black",size=1.1)

#DTG as switch
dtg_plot_scen2<-filter(outcome_scen,name1=="No DTG" | (name1=="DTG as a first-line and switch" & alpha3 == 1 & name2!="4) all men and 99% of women on DTG")) %>%
  ggplot( mapping=aes(x=year,y=outcome1,
                      color=factor(name2,labels=c("All on NNRTI  ", "Only men on DTG  ",
                                                  "All men and women not of childbearing age on DTG  ",
                                                  "All men and women not at risk of pregnancy on DTG  ",
                                                  "All men and women on DTG  ")))) +
  geom_line(size=1.1)+
  theme_bw()+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2040))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.6))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("NNRTI TDR (in %)")+
  geom_vline(xintercept = 2020, linetype="dashed", 
             color = "black", size=0.8)+
  annotate("text",x=2021,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  scale_color_manual(values=rep(color,2))+
  geom_line(data=outcome_2020,color="black",size=1.1)+
  geom_line(data=outcome_scen[outcome_scen$name1=="No DTG",],color="black",size=1.1)

#DTG as switch
dtg_plot_scen3<-filter(outcome_scen,name1=="No DTG" | (name1=="DTG as a first-line and switch" & alpha3 == 1.46 & name2!="4) all men and 99% of women on DTG")) %>%
  ggplot( mapping=aes(x=year,y=outcome1,
                      color=factor(name2,labels=c("All on NNRTI  ", "Only men on DTG  ",
                                                  "All men and women not of childbearing age on DTG  ",
                                                  "All men and women not at risk of pregnancy on DTG  ",
                                                  "All men and women on DTG  ")))) +
  geom_line(size=1.1)+
  theme_bw()+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2040))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.6))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("NNRTI TDR (in %)")+
  geom_vline(xintercept = 2020, linetype="dashed", 
             color = "black", size=0.8)+
  annotate("text",x=2021,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  scale_color_manual(values=rep(color,2))+
  geom_line(data=outcome_2020,color="black",size=1.1)+
  geom_line(data=outcome_scen[outcome_scen$name1=="No DTG",],color="black",size=1.1)

dtg_plot_scen<- plot_grid(dtg_plot_scen1+theme(legend.position="none")+theme(plot.margin = unit(c(0,0,0,10)/10, "cm")),
                          dtg_plot_scen2+theme(legend.position="none"),
                          dtg_plot_scen3+theme(legend.position="none"),
                          align = 'vh',
                          labels = c("A", "B","C"),
                          hjust = -1,
                          nrow = 3)
dtg_plot <- plot_grid(dtg_plot_main+theme(plot.margin = unit(c(10,0,0,0)/10, "cm")),
                      dtg_plot_scen+theme(plot.margin = unit(c(10,0,0,0)/10, "cm")),
                      labels=c("Baseline model",title_scen),
                      ncol=2)
dtg_plot <- plot_grid(dtg_plot,
              legend_main, rel_heights = c(1,0.1),nrow=2)#8.6x7 inches
dtg_plot
ggsave(filename=filename,width=22,height=20,unit="cm")
##############################################################################################################################
#Impact of NRTI-resistance and DTG-efficacy
#From 2005 to 2020
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


SIR_solve=function(t,y,parms){
  theta=parms[[1]]
  x_r=0
  x_i=0
  return(SIR_NRTI2(t,y,theta,x_r,x_ir))
}
#data_2020 = (ode(c(x_start,rep(0,8)), times=times,func=SIR_solve,parms=list(params_init)))[,1+1:(24*4*2*2*2+2)]
load("../R_results/data_2020.RData")

#From 2020 to 2040
scen=c("0) all on NNRTI","1) all men on DTG","2) all men and women not of childbearing age on DTG",
       "3) all men and women not at risk of pregancy on DTG","4) all men and women on DTG","4) all men and 99% of women on DTG")
scenario = data.frame(name1=c("No DTG",rep("DTG as a first-line and switch",4*9)),
                      name2=c(scen[1],rep(scen[2:5],9)),
                      dtg_1st = c(0,rep(1,4*9)),
                      dtg_switch = c(0,rep(1,4*9)),
                      alpha3 = c(1,rep(c(1,1.55,2.88, 1,1.49,2.71, 1,1.43,2.44),each=4)), #Impact of NRTI resistance on DTG-based regimen
                      alpha5 = c(1,rep(rep(c(0.85,1.25,2.02),each=3),each=4)), #Scaling factor for DTG
                      p_dtg = c(0,rep(c(0,0.175,0.63,1),9)))

times=181+0:240
data_scen_nrti_res=list()
for(i in 1:(dim(scenario)[1])){
  #for(i in 1:9){
  params=params_init
  params[c("dtg_1st","dtg_switch","p_dtg","alpha3","alpha5")]=as.numeric(scenario[i,c("dtg_1st","dtg_switch","p_dtg","alpha3","alpha5")])
  x_start=x_start2020(data_2020[181,],p_TDF = 1.0, p_DTG = scenario$p_dtg[i])
  data = (ode(c(x_start,rep(0,8)), times=times,func=SIR_solve,parms=list(params)))[,2:779] #remove time
  data_scen_nrti_res[[i]] = data[,1:(24*4*2*2*2+2)]
  print(i)
}
save(data_2020,data_scen_nrti_res,file="../R_results/scen_nrti_res.RData")
##################################################################################################################################
#Output
load("../R_results/scen_nrti_res.RData")

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
outcome_2020$alpha5=scenario[1,"alpha5"]
outcome_2020$outcome1=outcome(data_2020)

#output: scenarios from 2020 to 2040
for(i in 1:(dim(scenario)[1])){
  outcome_scen_i=data.frame(year=seq(2020,2040,1/12))
  outcome_scen_i$name1=scenario[i,"name1"]
  outcome_scen_i$name2=scenario[i,"name2"]
  outcome_scen_i$alpha3=scenario[i,"alpha3"]
  outcome_scen_i$alpha5=scenario[i,"alpha5"]
  outcome_scen_i$outcome1=outcome(data_scen_nrti_res[[i]])
  if(i==1){ outcome_scen = outcome_scen_i }else{ outcome_scen = rbind(outcome_scen,outcome_scen_i) }
}

#Transform scenario "No DTG" so that it will plotted for every combination of ORs
d_noDTG<-filter(outcome_scen,name1=="No DTG")
length_time=dim(d_noDTG)[1]
d_noDTG<-do.call("rbind", replicate(9, d_noDTG, simplify = FALSE))
d_noDTG<-mutate(d_noDTG,
          alpha3 = c(rep(c(1,1.55,2.88, 1,1.49,2.71, 1,1.43,2.44),each=length_time)), #Impact of NRTI resistance on DTG-based regimen
          alpha5 = c(rep(rep(c(0.85,1.25,2.02),each=3),each=length_time)))
outcome_scen2=rbind(d_noDTG,
  filter(outcome_scen,name1=="DTG as a first-line and switch"))
  
##################################################################################################################################
#Plot
color=c("black","blue","violet","orange","red")
plot1<-mutate(outcome_scen2,
              alpha3_type = factor(alpha3,levels=sort(c(unique(alpha3))),labels=c("OR(DTG failure|NRTI res) = 1",
                                                                                  rep("OR(DTG failure|NRTI res) = 2",3),
                                                                                  rep("OR(DTG failure|NRTI res) = 5",3))),
              alpha5_type = factor(alpha5,levels=sort(c(unique(alpha5))),labels=c("OR(failure| NNRTI vs DTG) = 1.02",
                                                                                  "OR(failure| NNRTI vs DTG) = 2",
                                                                                  "OR(failure| NNRTI vs DTG) = 5"))) %>%
    ggplot( mapping=aes(x=year,y=outcome1,
                        color=factor(name2,labels=c("No DTG", "Only men on DTG  ",
                                                    "All men and women not of childbearing age on DTG  ",
                                                    "All men and women not at risk of pregnancy on DTG  ",
                                                    "All men and women on DTG  ")))) +
    geom_line(size=1.1)+
    theme_bw()+
    facet_grid(alpha3_type ~ alpha5_type)+
    theme(legend.position="bottom")+
    scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2010,2040))+
    scale_y_continuous(labels = scales::percent,limits=c(0,0.6))+
    labs(color="Scenarios")+
    xlab("Calendar year")+
    ylab("NNRTI TDR (in %)")+
    geom_vline(xintercept = 2020, linetype="dashed", 
               color = "black", size=0.8)+
    annotate("text",x=2021,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
    scale_color_manual(values=rep(color,9))+
    geom_line(data=outcome_scen[outcome_scen$name1=="No DTG",],color="black",size=1.1)+
    geom_line(data=outcome_2020,color="black",size=1.1)

legend<-get_legend(
  plot1 + theme(legend.box.margin = margin(0, 0, 0, 0))+
    theme(legend.position = "bottom",)+
    #legend.spacing.x = unit(1.0, 'cm')) +
    theme(legend.direction = "horizontal")+
    theme(legend.title=element_blank())+
    guides(colour=guide_legend(nrow=2,byrow=FALSE,
                               keywidth=1.5,))+
    theme(legend.text=element_text(size=10))
)

plot1 <- plot_grid(plot1+
            theme(legend.position="none",plot.margin = unit(c(0.5,0.5,0.5,0.5)/10, "cm")),
          legend, rel_heights = c(1, 0.08),nrow=2)#8.6x7 inches
plot1
ggsave(filename="../figures/sens_nrti_res.pdf",width=22,height=20,unit="cm")
