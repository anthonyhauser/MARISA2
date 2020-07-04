library(plotly)
library("tidyverse")
library("cowplot")

##############################################################################################################################
##############################################################################################################################
#Scenarios
source("C:/Users/ahauser/Documents/Step2/Step2_revised/models/SIR_NRTI_res_function.R")
scen=c("0) all on NNRTI","1) all men on DTG","2) all men and women not of childbearing age on DTG",
       "3) all men and women not at risk of pregancy on DTG","4) all men and women on DTG","4) all men and 99% of women on DTG")
scenario = data.frame(name1=c("No DTG",rep("DTG as a first-line",5),rep("DTG as a first-line and switch",5),rep("DTG as a first-line and switch",5)),
                      name2=c(scen[1:6],scen[2:6],scen[2:6]),
                      dtg_1st = c(0,rep(1,15)),
                      dtg_switch = c(0,rep(0,5),rep(1,10)),
                      p_dtg = c(0,c(0,0.175,0.63,1,0.99),c(0,0.175,0.63,1,0.99),c(0,0.175,0.63,1,0.99)),
                      alpha3 = c(rep(1,11),rep(1.46,5))) #Impact of NRTI resistance on DTG-based regimen
##############################################################################################################################
##############################################################################################################################
#Outcome: ART use
#Load data from ubelix: data_2020, data_scen (list of 15 data frame), data_fail (list of 15 data frame)
load(paste("../R_results/sensitivity3/sens_",1,".RData",sep=""))
#Output from ubelix results
outcome1=function(data){
  tmax=dim(data)[1]
  nnrti=rowSums(data[seq(1,tmax,1),posa(c(4:15),1:4,1:2,1:2,1:2)])/rowSums(data[seq(1,tmax,1),posa(c(4:24),1:4,1:2,1:2,1:2)])
  return(nnrti)
}
outcome2=function(data){
  tmax=dim(data)[1]
  dtg=rowSums(data[seq(1,tmax,1),posa(c(16:21),1:4,1:2,1:2,1:2)])/rowSums(data[seq(1,tmax,1),posa(c(4:24),1:4,1:2,1:2,1:2)])
  return(dtg)
}

outcome_2020=data.frame(year=data.frame(year=rep(seq(2005,2020,1/12),2)))
outcome_2020$name1=scenario[1,"name1"]
outcome_2020$name2=scenario[1,"name2"]
outcome_2020$alpha3=scenario[1,"alpha3"]
outcome_2020$outcome=c(outcome1(data_2020),outcome2(data_2020))
outcome_2020$art=rep(c("NNRTI","DTG"),each=15*12+1)
outcome_2020<-rbind(mutate(outcome_2020,name1="DTG as a first-line"),
                    mutate(outcome_2020,name1="DTG as a first-line and switch"))

for(i in 1:(dim(scenario)[1])){
  data=data_scen[[i]]
  
  outcome_scen_i=data.frame(year=rep(seq(2020,2040,1/12),2))
  outcome_scen_i$name1=scenario[i,"name1"]
  outcome_scen_i$name2=scenario[i,"name2"]
  outcome_scen_i$alpha3=scenario[i,"alpha3"]
  outcome_scen_i$outcome=c(outcome1(data),outcome2(data))
  outcome_scen_i$art=rep(c("NNRTI","DTG"),each=20*12+1)
  if(i==1){ outcome_scen = outcome_scen_i }else{ outcome_scen = rbind(outcome_scen,outcome_scen_i) }
}

outcome_scen<-rbind( outcome_scen,
                 filter(outcome_scen,name1=="No DTG") %>%mutate(name1="DTG as a first-line"),
                 filter(outcome_scen,name1=="No DTG") %>%mutate(name1="DTG as a first-line and switch"),
                 filter(outcome_scen,name1=="No DTG") %>%mutate(name1="DTG as a first-line and switch, NRTI res.", alpha3=1.46))
##############################################################################################################################
#Plot
color=c("black","blue","violet","orange","red")
#DTG as first-line
data1 = filter(outcome_scen,name1 =="DTG as a first-line") %>%
  filter(name2 != "4) all men and 99% of women on DTG")
dtg_plot_s1<- ggplot(data=data1, mapping=aes(x=year,y=outcome,
                                            color=factor(name2,labels=c("All on NNRTI  ", "Only men on DTG  ",
                                                                        "All men and women not of childbearing age on DTG  ",
                                                                        "All men and women not at risk of pregnancy on DTG  ",
                                                                        "All men and women on DTG  ")))) +
  facet_grid( ~ art)+
  geom_line(size=1.1)+
  theme_bw()+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2040))+
  scale_y_continuous(labels = scales::percent,limits=c(0,1))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("Relative NNRTI/DTG use (in %)")+
  geom_vline(xintercept = 2020, linetype="dashed", 
             color = "black", size=1)+
  annotate("text",x=2011,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  scale_color_manual(values=rep(color,2))+
  geom_line(data=outcome_2020,aes(x=year,y=outcome),color="black",size=1.1)+
  geom_line(data=data1[data1$name2=="0) all on NNRTI",],aes(x=year,y=outcome),color="black",size=1.1)

data2 = filter(outcome_scen,name1 =="DTG as a first-line and switch") %>%
  filter(name2 != "4) all men and 99% of women on DTG")
dtg_plot_s2<- ggplot(data=data2, mapping=aes(x=year,y=outcome,
                                             color=factor(name2,labels=c("All on NNRTI  ", "Only men on DTG  ",
                                                                         "All men and women not of childbearing age on DTG  ",
                                                                         "All men and women not at risk of pregnancy on DTG  ",
                                                                         "All men and women on DTG  ")))) +
  facet_grid( ~ art)+
  geom_line(size=1.1)+
  theme_bw()+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2040))+
  scale_y_continuous(labels = scales::percent,limits=c(0,1))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("Relative NNRTI/DTG use (in %)")+
  geom_vline(xintercept = 2020, linetype="dashed", 
             color = "black", size=1)+
  annotate("text",x=2011,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  scale_color_manual(values=rep(color,2))+
  geom_line(data=outcome_2020,aes(x=year,y=outcome),color="black",size=1.1)+
  geom_line(data=data2[data2$name2=="0) all on NNRTI",],aes(x=year,y=outcome),color="black",size=1.1)


#Legend and aggregated plot
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
                     align = 'vh',
                     labels = c("A", "B"),
                     hjust = -1,
                     nrow = 2)

dtg_plot <- plot_grid(dtg_plot+
            theme(plot.margin = unit(c(0.5,0.5,0.5,0.5)/10, "cm")),
          legend, rel_heights = c(1, 0.1),nrow=2)#8.6x7 inches
dtg_plot
ggsave(filename="../figures/Fig2.pdf",width=20,height=20,unit="cm")
##############################################################################################################################
##############################################################################################################################
#Outcome: NNRTI resistance
#Load data from ubelix: data_2020, data_scen (list of 15 data frame), data_fail (list of 15 data frame)
#Output from ubelix results
list_data_2020=list()
list_data=list()
for(j in 1:201){
  load(paste("../R_results/sensitivity3/sens_",j,".RData",sep=""))
  list_data[[j]]=data_scen
  list_data_2020[[j]]=data_2020
  print(j)
}

#output: no DTG, from 2005 to 2020
outcome=function(data){
  tmax=dim(data)[1]
  outcome=rowSums(data[seq(1,tmax,1),posa(c(2,3),1:4,1:2,2,1:2)])/rowSums(data[seq(1,tmax,1),posa(c(2,3),1:4,1:2,1:2,1:2)])
  return(outcome)
}
outcome_2020=data.frame(year=seq(2005,2020,1/12))
outcome_2020$name1=scenario[1,"name1"]
outcome_2020$name2=scenario[1,"name2"]
outcome_2020$alpha3=scenario[1,"alpha3"]
outcome_2020$outcome1=outcome(data_2020)
outcome_2020[,c("out_inf","out_sup")]=t(apply(sapply(list_data_2020,function(x){outcome(x)}),1,function(y) quantile(y,probs=c(0.025,0.975))))

#output: scenarios from 2020 to 2040
for(i in 1:(dim(scenario)[1])){
 #With simulations from ub2.job
 outcome_scen_i=data.frame(year=seq(2020,2040,1/12))
 outcome_scen_i$name1=scenario[i,"name1"]
 outcome_scen_i$name2=scenario[i,"name2"]
 outcome_scen_i$alpha3=scenario[i,"alpha3"]
 outcome_scen_i$outcome1=outcome((list_data[[1]])[[i]])
 outcome_scen_i[,c("out_inf","out_sup")]=t(apply(sapply(list_data,function(x){outcome(x[[i]])}),1,function(y) quantile(y,probs=c(0.025,0.975))))
 if(i==1){ outcome_scen = outcome_scen_i }else{ outcome_scen = rbind(outcome_scen,outcome_scen_i) }
 print(i)
}

save(outcome_scen,outcome_2020,file="../R_results/nnrti_res_scen.RData")
##############################################################################################################################
#Write estimate into a .tex file
scen=c("0) all on NNRTI","1) all men on DTG","2) all men and women not of childbearing age on DTG",
       "3) all men and women not at risk of pregancy on DTG","4) all men and women on DTG","4) all men and 99% of women on DTG")
scenario = data.frame(name1=c("No DTG",rep("DTG as a first-line",5),rep("DTG as a first-line and switch",5),rep("DTG as a first-line and switch",5)),
                      name2=c(scen[1:6],scen[2:6],scen[2:6]),
                      dtg_1st = c(0,rep(1,15)),
                      dtg_switch = c(0,rep(0,5),rep(1,10)),
                      p_dtg = c(0,c(0,0.175,0.63,1,0.99),c(0,0.175,0.63,1,0.99),c(0,0.175,0.63,1,0.99)),
                      alpha3 = c(rep(1,11),rep(1.46,5))) #Impact of NRTI resistance on DTG-based regimen

load("../R_results/nnrti_res_scen.RData")
data<-filter(outcome_scen,year==2040 | year==2030)

get("{")
name1=c("a","b","c","da","ea","fa","fc","ga","gc","ha","hc","ia","ic")
value1=as.numeric(c(filter(data,year==2040,name1=="No DTG",name2==scen[1],alpha3==1) %>% select(outcome1),
         filter(data,year==2040,name1=="DTG as a first-line and switch",name2==scen[5],alpha3==1) %>% select(outcome1),
         filter(data,year==2040,name1=="DTG as a first-line and switch",name2==scen[4],alpha3==1) %>% select(outcome1),
         
         filter(data,year==2030,name1=="No DTG",name2==scen[1],alpha3==1) %>% select(outcome1),
         filter(data,year==2040,name1=="No DTG",name2==scen[1],alpha3==1) %>% select(outcome1),
         
         filter(data,year==2030,name1=="DTG as a first-line and switch",name2==scen[5],alpha3==1) %>% select(outcome1),
         filter(data,year==2030,name1=="DTG as a first-line and switch",name2==scen[5],alpha3==1.46) %>% select(outcome1),
         filter(data,year==2040,name1=="DTG as a first-line and switch",name2==scen[5],alpha3==1) %>% select(outcome1),
         filter(data,year==2040,name1=="DTG as a first-line and switch",name2==scen[5],alpha3==1.46) %>% select(outcome1),
         
         filter(data,year==2030,name1=="DTG as a first-line and switch",name2==scen[4],alpha3==1) %>% select(outcome1),
         filter(data,year==2030,name1=="DTG as a first-line and switch",name2==scen[4],alpha3==1.46) %>% select(outcome1),
         filter(data,year==2040,name1=="DTG as a first-line and switch",name2==scen[4],alpha3==1) %>% select(outcome1),
         filter(data,year==2040,name1=="DTG as a first-line and switch",name2==scen[4],alpha3==1.46) %>% select(outcome1)))
l=length(name1)
output=c()
for(i in 1:l){
  output[i] = paste0("\\newcommand{","\\number",name1[i],"}{",round(value1[i]*1000)/10,"\\% ","}")
}
cat(paste0(output,collapse="\n"),file="../R_results/paper_estimates1.tex")

name=c("db","eb","fb","fd","gb","gd","hb","hd","ib","id")
value1=as.numeric(c(filter(data,year==2030,name1=="No DTG",name2==scen[1],alpha3==1) %>% select(out_inf),
         filter(data,year==2040,name1=="No DTG",name2==scen[1],alpha3==1) %>% select(out_inf),
         
         filter(data,year==2030,name1=="DTG as a first-line and switch",name2==scen[5],alpha3==1) %>% select(out_inf),
         filter(data,year==2030,name1=="DTG as a first-line and switch",name2==scen[5],alpha3==1.46) %>% select(out_inf),
         filter(data,year==2040,name1=="DTG as a first-line and switch",name2==scen[5],alpha3==1) %>% select(out_inf),
         filter(data,year==2040,name1=="DTG as a first-line and switch",name2==scen[5],alpha3==1.46) %>% select(out_inf),
         
         filter(data,year==2030,name1=="DTG as a first-line and switch",name2==scen[4],alpha3==1) %>% select(out_inf),
         filter(data,year==2030,name1=="DTG as a first-line and switch",name2==scen[4],alpha3==1.46) %>% select(out_inf),
         filter(data,year==2040,name1=="DTG as a first-line and switch",name2==scen[4],alpha3==1) %>% select(out_inf),
         filter(data,year==2040,name1=="DTG as a first-line and switch",name2==scen[4],alpha3==1.46) %>% select(out_inf)))

value2=as.numeric(c(filter(data,year==2030,name1=="No DTG",name2==scen[1],alpha3==1) %>% select(out_sup),
         filter(data,year==2040,name1=="No DTG",name2==scen[1],alpha3==1) %>% select(out_sup),
         
         filter(data,year==2030,name1=="DTG as a first-line and switch",name2==scen[5],alpha3==1) %>% select(out_sup),
         filter(data,year==2030,name1=="DTG as a first-line and switch",name2==scen[5],alpha3==1.46) %>% select(out_sup),
         filter(data,year==2040,name1=="DTG as a first-line and switch",name2==scen[5],alpha3==1) %>% select(out_sup),
         filter(data,year==2040,name1=="DTG as a first-line and switch",name2==scen[5],alpha3==1.46) %>% select(out_sup),
         
         filter(data,year==2030,name1=="DTG as a first-line and switch",name2==scen[4],alpha3==1) %>% select(out_sup),
         filter(data,year==2030,name1=="DTG as a first-line and switch",name2==scen[4],alpha3==1.46) %>% select(out_sup),
         filter(data,year==2040,name1=="DTG as a first-line and switch",name2==scen[4],alpha3==1) %>% select(out_sup),
         filter(data,year==2040,name1=="DTG as a first-line and switch",name2==scen[4],alpha3==1.46) %>% select(out_sup)))

l=length(name)
output=c()
for(i in 1:l){
  output[i] = paste0("\\newcommand{","\\number",name[i],"}{",round(value1[i]*1000)/10,"\\%","-",round(value2[i]*1000)/10,"\\%","}")
}
cat(paste0(output,collapse="\n"),file="../R_results/paper_estimates2.tex")

##############################################################################################################################
#Plot
color=c("black","blue","violet","orange","red")
#DTG as first-line
data = filter(outcome_scen,name1=="No DTG" | (name1=="DTG as a first-line" & name2!="4) all men and 99% of women on DTG")) %>%
       filter(alpha3 == 1)
dtg_plot_s1<- ggplot(data=data, mapping=aes(x=year,y=outcome1,
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
  geom_rect(xmin=2040.5,xmax=2044.5,ymin=0,ymax=0.7,colour="grey50",fill="white",size=.3)+
  geom_segment(x=2040,xend=2040.5,y=-.03,yend=0,size=.3,colour="grey50")+
  geom_pointrange(data=data[data$year==2040,],aes(x=rep(seq(2041,2044,length.out = 5),1),y=outcome1,ymin=out_inf,ymax=out_sup),alpha=1,fill="white",size=.8,stroke=.4,shape=1)+
  geom_line(data=outcome_2020,color="black",size=1.1)+
  geom_line(data=outcome_scen[outcome_scen$name1=="No DTG",],color="black",size=1.1)

#DTG as first-line
data = filter(outcome_scen,name1=="No DTG" | (name1=="DTG as a first-line and switch" & alpha3 == 1 & name2!="4) all men and 99% of women on DTG"))
dtg_plot_s2<- ggplot(data=data, mapping=aes(x=year,y=outcome1,
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
  geom_rect(xmin=2040.5,xmax=2044.5,ymin=0,ymax=0.7,colour="grey50",fill="white",size=.3)+
  geom_segment(x=2040,xend=2040.5,y=-.03,yend=0,size=.3,colour="grey50")+
  geom_pointrange(data=data[data$year==2040,],aes(x=rep(seq(2041,2044,length.out = 5),1),y=outcome1,ymin=out_inf,ymax=out_sup),alpha=1,fill="white",size=.8,stroke=.4,shape=1)+
  geom_line(data=outcome_2020,color="black",size=1.1)+
  geom_line(data=outcome_scen[outcome_scen$name1=="No DTG",],color="black",size=1.1)


#DTG as first-line, NRTI resistance
data = filter(outcome_scen,name1=="No DTG" | (name1=="DTG as a first-line and switch" & alpha3 == 1.46 & name2!="4) all men and 99% of women on DTG"))
dtg_plot_s3<- ggplot(data=data, mapping=aes(x=year,y=outcome1,
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
  geom_rect(xmin=2040.5,xmax=2044.5,ymin=0,ymax=0.7,colour="grey50",fill="white",size=.3)+
  geom_segment(x=2040,xend=2040.5,y=-.03,yend=0,size=.3,colour="grey50")+
  geom_pointrange(data=data[data$year==2040,],aes(x=rep(seq(2041,2044,length.out = 5),1),y=outcome1,ymin=out_inf,ymax=out_sup),alpha=1,fill="white",size=.8,stroke=.4,shape=1)+
  geom_line(data=outcome_2020,color="black",size=1.1)+
  geom_line(data=outcome_scen[outcome_scen$name1=="No DTG",],color="black",size=1.1)

#Legend and aggregated plot
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
dtg_plot<-plot_grid(dtg_plot+
            theme(plot.margin = unit(c(0.5,0.5,0.5,0.5)/10, "cm")),
          legend, rel_heights = c(1, 0.08),nrow=2)#8.6x7 inches
dtg_plot
ggsave(filename="../figures/Fig3.pdf",width=20,height=22,unit="cm")
##############################################################################################################################
##############################################################################################################################
#Outcome: NNRTI failure

list_scen_fail=list()
for(j in 1:201){
  load(paste("../R_results/sensitivity3/sens_",j,".RData",sep=""))
  scen_fail$year=rep(rep(c(2025,2030,2035),each=2),dim(scen_fail)[1]/6)
  list_scen_fail[[j]]=scen_fail
  print(j)
}

scen_fail=list_scen_fail[[1]]
scen_fail[,c("out_inf","out_sup")]=t(apply(sapply(list_scen_fail,function(x){x[,"prop_fail"]}),1,function(y) quantile(y,probs=c(0.025,0.975))))
scen_fail <- filter(scen_fail,year %in% c(2035))
scen_fail$name1[scen_fail$alpha3==1.46]="DTG as a first-line and switch, NRTI res."
scen_fail$name1=factor(scen_fail$name1,levels=c("No DTG",
                                                "DTG as a first-line",
                                                "DTG as a first-line and switch",
                                                "DTG as a first-line and switch, NRTI res."),
                                       labels=c("No DTG",
                                                "DTG as a \nfirst-line",
                                               "DTG as a first-\nline and switch",
                                               "DTG as a first-\nline and switch\nwith NRTI-res."))
scen_fail$time=factor(scen_fail$time,levels=c(12,24),labels=c("After 12 months",
                                                "After 24 months"))
scen_fail$name2=factor(scen_fail$name2,labels=c("All on NNRTI  ", "Only men on DTG  ",
                                      "All men and women not of childbearing age on DTG  ",
                                      "All men and women not at risk of pregnancy on DTG  ",
                                      "All men and 99% of women on DTG  "))
##############################################################################################################################
get("{")
name1=c("k","l","lb","n","nb","m","mb")
value1=as.numeric(c(filter(data,name1=="No DTG",time=="After 24 months") %>% select(prop_fail),
                    filter(data,name1=="DTG as a first-\nline and switch",name2=="Only men on DTG  ",time=="After 24 months") %>% select(prop_fail),
                    filter(data,name1=="DTG as a first-\nline and switch\nwith NRTI-res.",name2=="Only men on DTG  ",time=="After 24 months") %>% select(prop_fail),
                    filter(data,name1=="DTG as a first-\nline and switch",name2=="All men and women not at risk of pregnancy on DTG  ",time=="After 24 months") %>% select(prop_fail),
                    filter(data,name1=="DTG as a first-\nline and switch\nwith NRTI-res.",name2=="All men and women not at risk of pregnancy on DTG  ",time=="After 24 months") %>% select(prop_fail),
                    filter(data,name1=="DTG as a first-\nline and switch",name2=="All men and 99% of women on DTG  ",time=="After 24 months") %>% select(prop_fail),
                    filter(data,name1=="DTG as a first-\nline and switch\nwith NRTI-res.",name2=="All men and 99% of women on DTG  ",time=="After 24 months") %>% select(prop_fail)))
l=length(name1)
output=c()
for(i in 1:l){
  output[i] = paste0("\\newcommand{","\\number",name1[i],"}{",round(value1[i]*1000)/10,"\\% ","}")
}
cat(paste0(output,collapse="\n"),file="../R_results/paper_estimates3.tex")

##############################################################################################################################
#Plot
color=c("black","blue","violet","orange","red")

pd = position_dodge(0.7)
data=filter(scen_fail,year==2035)
dtg_plot_fail = ggplot(data,
       aes(x =name1, 
           y =prop_fail,
           color = name2)) +
  facet_grid(. ~ time)+
  
  geom_pointrange(aes(ymin = out_inf,
                    ymax = out_sup),
                size = 0.7, 
                position = pd) +
  theme_bw() +
  theme(axis.title = element_text(size=15)) +
  ylab("NNRTI-failure (%)")+
  xlab("Strategy of DTG-introduction")+
  theme(strip.text.x = element_text(size = 12, colour = "black", angle = 0))+
  theme(axis.text.x = element_text(size=12,angle = 0))+
  theme(axis.text.y = element_text(size=12,angle = 0))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.21))+
  scale_color_manual(values=rep(color[c(1,2,3,4,5)],2))+
  labs(color="DTG-prescription strategy")+
  theme(legend.box.margin = margin(0, 0, 0, 0))+
  theme(legend.position = "bottom",)+
  theme(legend.direction = "horizontal")+
  theme(legend.title=element_blank())+
  guides(colour=guide_legend(nrow=2,byrow=FALSE,
                             keywidth=1.5,))+
  theme(legend.text=element_text(size=11))+
  theme(axis.title.x=element_blank())
dtg_plot_fail

ggsave(filename="../figures/Fig5.pdf",width=30,height=16,unit="cm")


##################################################################################################################
##################################################################################################################
#Heatmap
# with inverted x-axis, time is equally spaced

#Load data
i_len=39
j_len=51
heatmap_supp=matrix(NA,ncol=i_len,nrow=j_len)
heatmap_fail=matrix(NA,ncol=i_len,nrow=j_len)
for(i in 1:i_len){
  for(j in 1:j_len){
    load(paste("../R_results/heatmap2/heatmap_",i,"_",j,".RData",sep=""))
    heatmap_supp[j,i]=res_suppressed
    heatmap_fail[j,i]=res_failed
  }
  print(i)
}

min1=min(heatmap_supp)
min2=min(heatmap_fail)
max1=max(heatmap_supp)
max2=max(heatmap_fail)
min=min(min1,min2)
max=max(max1,max2)

##################################################################################################################
get("{")
name1=c("ja","jb")
value1=as.numeric(c(min,max))
l=length(name1)
output=c()
for(i in 1:l){
  output[i] = paste0("\\newcommand{","\\number",name1[i],"}{",round(value1[i]*1000)/10,"\\% ","}")
}
cat(paste0(output,collapse="\n"),file="../R_results/paper_estimates4.tex")
##################################################################################################################
#Plot
vals1 <- unique(scales::rescale(c(heatmap_supp),to=c((min1-min)/(max-min),(max1-min)/(max-min))))
vals1_01 <- unique(scales::rescale(c(heatmap_supp),to=c(0,1)))
o <- order(vals1_01, decreasing = FALSE)
cols <- scales::col_numeric("Blues", domain = c(0,1))(vals1)
colz1 <- setNames(data.frame((vals1_01)[o], cols[o]), NULL)

#which_x<-c(seq(1,i_len,by=5),i_len)
which_x<-which(seq(0.5,10,length.out = 39) %in% c(0.5,2,4,6,8,10))

plot_dtg<-plot_ly(z=heatmap_supp,colorscale = colz1,showscale=FALSE,colorbar = list(thickness=30,len=0.1,title = "",titlefont=list(size=2),tickfont=list(size=2)), type = "contour") %>%
  plotly::layout(title="    A. Suppressed ind.",
                 margin=c(1,1,1,1)/5,
                 xaxis=list(title="Time to switch to DTG (year)",list(thickness=1,len=0.1,title = "NNRTI PDR \n level in 2035",titlefont=list(size=1),tickfont=list(size=1)),titlefont=list(family="Arial",size=24),ticktext=seq(0.5,10,length.out = 39)[which_x],tickvals=which_x-1,tickmode="array",tickfont=list(size=18)),
                 yaxis=list(title="% women eligible for DTG",titlefont=list(family="Arial",size=24),ticktext=seq(0,100,by=10),tickvals=seq(0,50,by=5),tickmode="array",tickfont=list(size=18))
  )%>%
  colorbar(limits = c(min,max))

#Plot 2
vals2 <- unique(scales::rescale(c(heatmap_fail),to=c((min2-min)/(max-min),(max2-min)/(max-min))))
vals2_01 <- unique(scales::rescale(c(heatmap_fail),to=c(0,1)))
o <- order(vals2_01, decreasing = FALSE)
cols <- scales::col_numeric("Blues", domain = c(0,1))(vals2)
colz2 <- setNames(data.frame((vals2_01)[o], cols[o]), NULL)

plot_dtg_fail<-plot_ly(z=heatmap_fail,colorscale = colz1,colorbar=list(thickness=60,len=0.7,title = "NNRTI PDR \n level in 2035",titlefont=list(size=16),tickfont=list(size=18)), type = "contour") %>%
  plotly::layout(title="                                                    ",
                 margin=c(1,1,1,1)/5,
                 xaxis=list(title="Time to switch to DTG (year)",list(thickness=1,len=0.1,title = "NNRTI PDR \n level in 2035",titlefont=list(size=1),tickfont=list(size=1)),titlefont=list(family="Arial",size=24),ticktext=seq(0.5,10,length.out = 39)[which_x],tickvals=which_x-1,tickmode="array",tickfont=list(size=18)),
                 yaxis=list(title="% women eligible for DTG",titlefont=list(family="Arial",size=24),ticktext=seq(0,100,by=10),tickvals=seq(0,50,by=5),tickmode="array",tickfont=list(size=18))
  )%>%
  colorbar(limits = c(min,max))
#2 Plots
subplot(plot_dtg,plot_dtg_fail,margin=c(0.08,0.08,0.08,0.08),widths=c(0.5,0.5),titleX=TRUE,shareX = FALSE,titleY=TRUE,shareY = FALSE)
tiff("test.tiff", units="in", width=5, height=5, res=300)
#save as fig4.png
##############################################################################################################################
##############################################################################################################################
#Outcome: NNRTI resistance with uncertainty over time
#Output from ubelix results
load("../R_results/nnrti_res_scen.RData")

#No DTG
data = filter(outcome_scen,name1=="No DTG")
data = rbind(data,
             outcome_2020)
dtg_plot_1<- ggplot(data=data, mapping=aes(x=year,y=outcome1)) +
  geom_line(size=1.1,color="black")+
  theme_bw()+
  theme(legend.position="none")+
  geom_ribbon(aes(ymin=out_inf,ymax=out_sup),color=NA,fill="black",alpha=0.3)+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2040))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.7))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("NNRTI TDR (in %)")+
  geom_vline(xintercept = 2020, linetype="dashed", 
             color = "black", size=1)+
  annotate("text",x=2021,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)

legend<-get_legend(
  dtg_plot_s1 + theme(legend.box.margin = margin(0, 0, 0, 0))+
    theme(legend.position = "bottom",)+
    #legend.spacing.x = unit(1.0, 'cm')) +
    theme(legend.direction = "horizontal")+
    theme(legend.title=element_blank())+
    guides(colour=guide_legend(nrow=5,byrow=FALSE,
                               keywidth=1.5,))+
    theme(legend.text=element_text(size=10))
)


dtg_plot_p1<- plot_grid(dtg_plot_1+theme(plot.margin = unit(c(0,0,0,10)/10, "cm")),
                     legend,
                     rel_widths = c(1,1.8),
                     labels=c("A",""),
                     nrow = 1,
                     ncol=2)

#DTG as first-line
color=c("blue","violet","orange","red")
data = filter(outcome_scen,name1=="DTG as a first-line" & name2!="4) all men and 99% of women on DTG")
data = rbind(data,
             mutate(outcome_2020,name1="DTG as a first-line",name2=scen[2]),
             mutate(outcome_2020,name1="DTG as a first-line",name2=scen[3]),
             mutate(outcome_2020,name1="DTG as a first-line",name2=scen[4]),
             mutate(outcome_2020,name1="DTG as a first-line",name2=scen[5]))
dtg_plot_2<- ggplot(data=data, mapping=aes(x=year,y=outcome1,
                                            color=factor(name2,labels=c("Only men on DTG  ",
                                                                        "All men and women not of childbearing age on DTG  ",
                                                                        "All men and women not at risk of pregnancy on DTG  ",
                                                                        "All men and women on DTG  ")),
                                            fill=factor(name2,labels=c("Only men on DTG  ",
                                                                       "All men and women not of childbearing age on DTG  ",
                                                                       "All men and women not at risk of pregnancy on DTG  ",
                                                                       "All men and women on DTG  ")))) +
  geom_line(size=1.1)+
  facet_grid(name2 ~ .)+
  theme_bw()+
  theme(legend.position="none")+
  geom_ribbon(aes(ymin=out_inf,ymax=out_sup),color=NA,alpha=0.3)+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2040))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.7))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("NNRTI TDR (in %)")+
  geom_vline(xintercept = 2020, linetype="dashed", 
             color = "black", size=1)+
  annotate("text",x=2021,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  scale_color_manual(values=color)+
  scale_fill_manual(values=color)+
  theme(strip.background = element_blank(), strip.text.y = element_blank())

#DTG as first-line and switch
data = filter(outcome_scen,name1=="DTG as a first-line and switch" & name2!="4) all men and 99% of women on DTG" & alpha3==1)
data = rbind(data,
             mutate(outcome_2020,name1="DTG as a first-line and switch",name2=scen[2]),
             mutate(outcome_2020,name1="DTG as a first-line and switch",name2=scen[3]),
             mutate(outcome_2020,name1="DTG as a first-line and switch",name2=scen[4]),
             mutate(outcome_2020,name1="DTG as a first-line and switch",name2=scen[5]))
dtg_plot_3<- ggplot(data=data, mapping=aes(x=year,y=outcome1,
                                           color=factor(name2,labels=c("Only men on DTG  ",
                                                                       "All men and women not of childbearing age on DTG  ",
                                                                       "All men and women not at risk of pregnancy on DTG  ",
                                                                       "All men and women on DTG  ")),
                                           fill=factor(name2,labels=c("Only men on DTG  ",
                                                                      "All men and women not of childbearing age on DTG  ",
                                                                      "All men and women not at risk of pregnancy on DTG  ",
                                                                      "All men and women on DTG  ")))) +
  geom_line(size=1.1)+
  facet_grid(name2 ~ .)+
  theme_bw()+
  theme(legend.position="none")+
  geom_ribbon(aes(ymin=out_inf,ymax=out_sup),color=NA,alpha=0.3)+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2040))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.7))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("NNRTI TDR (in %)")+
  geom_vline(xintercept = 2020, linetype="dashed", 
             color = "black", size=1)+
  annotate("text",x=2021,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  scale_color_manual(values=color)+
  scale_fill_manual(values=color)+
  theme(strip.background = element_blank(), strip.text.y = element_blank())

#DTG as first-line
data = filter(outcome_scen,name1=="DTG as a first-line and switch" & name2!="4) all men and 99% of women on DTG" & alpha3==1.46)
data = rbind(data,
             mutate(outcome_2020,name1="DTG as a first-line and switch",name2=scen[2],alpha3=1.46),
             mutate(outcome_2020,name1="DTG as a first-line and switch",name2=scen[3],alpha3=1.46),
             mutate(outcome_2020,name1="DTG as a first-line and switch",name2=scen[4],alpha3=1.46),
             mutate(outcome_2020,name1="DTG as a first-line and switch",name2=scen[5],alpha3=1.46))
dtg_plot_4<- ggplot(data=data, mapping=aes(x=year,y=outcome1,
                                           color=factor(name2,labels=c("Only men on DTG  ",
                                                                       "All men and women not of childbearing age on DTG  ",
                                                                       "All men and women not at risk of pregnancy on DTG  ",
                                                                       "All men and women on DTG  ")),
                                           fill=factor(name2,labels=c("Only men on DTG  ",
                                                                      "All men and women not of childbearing age on DTG  ",
                                                                      "All men and women not at risk of pregnancy on DTG  ",
                                                                      "All men and women on DTG  ")))) +
  geom_line(size=1.1)+
  facet_grid(name2 ~ .)+
  theme_bw()+
  theme(legend.position="none")+
  geom_ribbon(aes(ymin=out_inf,ymax=out_sup),color=NA,alpha=0.3)+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2040))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.7))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("NNRTI TDR (in %)")+
  geom_vline(xintercept = 2020, linetype="dashed", 
             color = "black", size=1)+
  annotate("text",x=2021,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  scale_color_manual(values=color)+
  scale_fill_manual(values=color)+
  theme(strip.background = element_blank(), strip.text.y = element_blank())

dtg_plot_p2<- plot_grid(dtg_plot_2,
                        dtg_plot_3,
                        dtg_plot_4,
                        rel_widths = c(1,1,1),
                        labels=c("B","C","D"),
                        nrow = 1,
                        ncol=3)

plot_dtg <- plot_grid(dtg_plot_p1 + theme(plot.margin = unit(c(2,0,5,0)/10, "cm")),
          dtg_plot_p2,
          rel_heights = c(1,3.2),
          nrow = 2,
          ncol=1)
plot_dtg
ggsave(filename="../figures/sens_uncertainty.pdf",width=20,height=20,unit="cm")

##############################################################################################################################
##############################################################################################################################
#Outcome: NNRTI resistance difference
#Load data from ubelix: data_2020, data_scen (list of 15 data frame), data_fail (list of 15 data frame)
#Output from ubelix results
list_data_2020=list()
list_data=list()
for(j in 1:201){
  load(paste("../R_results/sensitivity3/sens_",j,".RData",sep=""))
  list_data[[j]]=data_scen
  list_data_2020[[j]]=data_2020
  print(j)
}

#output: no DTG, from 2005 to 2020
outcome=function(data){
  tmax=dim(data)[1]
  outcome=rowSums(data[seq(1,tmax,1),posa(c(2,3),1:4,1:2,2,1:2)])/rowSums(data[seq(1,tmax,1),posa(c(2,3),1:4,1:2,1:2,1:2)])
  return(outcome)
}
outcome_2020=data.frame(year=seq(2005,2020,1/12))
outcome_2020$name1=scenario[1,"name1"]
outcome_2020$name2=scenario[1,"name2"]
outcome_2020$alpha3=scenario[1,"alpha3"]
outcome_2020$outcome1=outcome(data_2020)
outcome_2020[,c("out_inf","out_sup")]=t(apply(sapply(list_data_2020,function(x){outcome(x)}),1,function(y) quantile(y,probs=c(0.025,0.975))))

#output: scenarios from 2020 to 2040
for(i in 2:(dim(scenario)[1])){
  #With simulations from ub2.job
  outcome_scen_i=data.frame(year=seq(2020,2040,1/12))
  outcome_scen_i$name1=scenario[i,"name1"]
  outcome_scen_i$name2=scenario[i,"name2"]
  outcome_scen_i$alpha3=scenario[i,"alpha3"]
  outcome_scen_i$outcome1=outcome((list_data[[1]])[[i]])-outcome((list_data[[1]])[[1]])
  outcome_scen_i[,c("out_inf","out_sup")]=t(apply(sapply(list_data,function(x){outcome(x[[i]])-outcome(x[[1]])}),1,function(y) quantile(y,probs=c(0.025,0.975))))
  if(i==2){ outcome_scen = outcome_scen_i }else{ outcome_scen = rbind(outcome_scen,outcome_scen_i) }
  print(i)
}

save(outcome_scen,file="../R_results/diff_nnrti_res_scen.RData")

##############################################################################################################################
#Plot of difference in NNRTI
load("../R_results/diff_nnrti_res_scen.RData")

#DTG as first-line
color=c("blue","violet","orange","red")
data = filter(outcome_scen,name1=="DTG as a first-line" & name2!="4) all men and 99% of women on DTG")
dtg_plot_1<- ggplot(data=data, mapping=aes(x=year,y=outcome1,
                                           color=factor(name2,labels=c("Only men on DTG  ",
                                                                       "All men and women not of childbearing age on DTG  ",
                                                                       "All men and women not at risk of pregnancy on DTG  ",
                                                                       "All men and women on DTG  ")),
                                           fill=factor(name2,labels=c("Only men on DTG  ",
                                                                      "All men and women not of childbearing age on DTG  ",
                                                                      "All men and women not at risk of pregnancy on DTG  ",
                                                                      "All men and women on DTG  ")))) +
  geom_line(size=1.1)+
  facet_grid(name2 ~ .)+
  theme_bw()+
  theme(legend.position="none")+
  geom_ribbon(aes(ymin=out_inf,ymax=out_sup),color=NA,alpha=0.3)+
  scale_x_continuous(breaks=c(2020,2030,2040),labels=c("2020","2030","2040"),limits=c(2020,2040))+
  scale_y_continuous(labels = scales::percent,limits=c(-0.52,0.1))+
  labs(color="Scenarios",fill="Scenarios")+
  xlab("Calendar year")+
  ylab("Difference in NNRTI TDR (in %)")+
  scale_color_manual(values=color)+
  scale_fill_manual(values=color)+
  theme(strip.background = element_blank(), strip.text.y = element_blank())

legend<-get_legend(
  dtg_plot_1 + theme(legend.box.margin = margin(0, 0, 0, 0))+
    theme(legend.position = "bottom",)+
    #legend.spacing.x = unit(1.0, 'cm')) +
    theme(legend.direction = "horizontal")+
    theme(legend.title=element_blank())+
    guides(colour=guide_legend(nrow=2,byrow=FALSE,
                               keywidth=1.5,))+
    theme(legend.text=element_text(size=10))
)

#DTG as first-line and switch
data = filter(outcome_scen,name1=="DTG as a first-line and switch" & name2!="4) all men and 99% of women on DTG" & alpha3==1)
dtg_plot_2<- ggplot(data=data, mapping=aes(x=year,y=outcome1,
                                           color=factor(name2,labels=c("Only men on DTG  ",
                                                                       "All men and women not of childbearing age on DTG  ",
                                                                       "All men and women not at risk of pregnancy on DTG  ",
                                                                       "All men and women on DTG  ")),
                                           fill=factor(name2,labels=c("Only men on DTG  ",
                                                                      "All men and women not of childbearing age on DTG  ",
                                                                      "All men and women not at risk of pregnancy on DTG  ",
                                                                      "All men and women on DTG  ")))) +
  geom_line(size=1.1)+
  facet_grid(name2 ~ .)+
  theme_bw()+
  theme(legend.position="none")+
  geom_ribbon(aes(ymin=out_inf,ymax=out_sup),color=NA,alpha=0.3)+
  scale_x_continuous(breaks=c(2020,2030,2040),labels=c("2020","2030","2040"),limits=c(2020,2040))+
  scale_y_continuous(labels = scales::percent,limits=c(-0.52,0.1))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("Difference in NNRTI TDR (in %)")+
  scale_color_manual(values=color)+
  scale_fill_manual(values=color)+
  theme(strip.background = element_blank(), strip.text.y = element_blank())

#DTG as first-line and switch, NRTI-res
data = filter(outcome_scen,name1=="DTG as a first-line and switch" & name2!="4) all men and 99% of women on DTG" & alpha3==1.46)
dtg_plot_3<- ggplot(data=data, mapping=aes(x=year,y=outcome1,
                                           color=factor(name2,labels=c("Only men on DTG  ",
                                                                       "All men and women not of childbearing age on DTG  ",
                                                                       "All men and women not at risk of pregnancy on DTG  ",
                                                                       "All men and women on DTG  ")),
                                           fill=factor(name2,labels=c("Only men on DTG  ",
                                                                      "All men and women not of childbearing age on DTG  ",
                                                                      "All men and women not at risk of pregnancy on DTG  ",
                                                                      "All men and women on DTG  ")))) +
  geom_line(size=1.1)+
  facet_grid(name2 ~ .)+
  theme_bw()+
  theme(legend.position="none")+
  geom_ribbon(aes(ymin=out_inf,ymax=out_sup),color=NA,alpha=0.3)+
  scale_x_continuous(breaks=c(2020,2030,2040),labels=c("2020","2030","2040"),limits=c(2020,2040))+
  scale_y_continuous(labels = scales::percent,limits=c(-0.52,0.1))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("Difference in NNRTI TDR (in %)")+
  scale_color_manual(values=color)+
  scale_fill_manual(values=color)+
  theme(strip.background = element_blank(), strip.text.y = element_blank())

dtg_plot_p1<- plot_grid(dtg_plot_1,
                        dtg_plot_2,
                        dtg_plot_3,
                        rel_widths = c(1,1,1),
                        labels=c("A","B","C"),
                        nrow = 1,
                        ncol=3)

plot_dtg <- plot_grid(dtg_plot_p1 + theme(plot.margin = unit(c(2,0,5,0)/10, "cm")),
                      legend,
                      rel_heights = c(1,0.1),
                      nrow = 2,
                      ncol=1)
plot_dtg
ggsave(filename="../figures/sens_diff_uncertainty.pdf",width=22,height=20,unit="cm")
