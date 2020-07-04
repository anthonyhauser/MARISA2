library(deSolve)
library(zoo)
library("Rcpp")
library("BH")
library("ggplot2")
library("dplyr")
library("cowplot")
##############################################################################################################################
#adapt params
#Assumptions
#load params from value_gender4_v3.R and p1, p2 from graphs.R
#source("C:/Users/ahauser/Documents/Step1/Rfiles/used scripts/marisa/value_gender4_v3.R")
#source("C:/Users/ahauser/Documents/Step2/dtg/value_gender4_v4.R")
source("C:/Users/ahauser/Documents/Step2/dtg/value_gender4_v5.R")
source("C:/Users/ahauser/Documents/Step2/dtg/parmin_dtg.R")
p1=c(rate1_inf=3.318999e+00, rate2_inf=5.440886e-01, rate1_diag=2.736363e+02, rate2_diag=7.710578e+00, rate_treat=1.879695e-03,rate_death=1.632000e-01)
p2=c(rate1_inf=NA, rate2_inf=NA, rate1_diag=NA, rate2_diag=NA, rate_treat=NA,rate_death=NA,q=0.05,rate_ratio=0.5,k1=1,k2=2,k3=2,alpha=2,rate_res=5,rate_susc=125)
#no treatment interruption
params<-within(params,{
  RateStopTreatFirst=4.35/c(Inf,Inf,Inf,Inf)
  RateStopSuppFirst=4.35/c(Inf,Inf,Inf,Inf)
  RateStopFailFirst=4.35/c(Inf,Inf,Inf,Inf)
  RateStopTreatSecond=4.35/c(Inf,Inf,Inf,Inf) #Assumption : T2 same as T1
  RateStopSuppSecond=4.35/c(Inf,Inf,Inf,Inf)
  RateStopFailSecond=4.35/c(Inf,Inf,Inf,Inf)
})
#no starting with second-line PI
params<-within(params,{
  #RateDirectTreatSecond=4.35/c(5000,12000,13000,1850)
  RateTreatSecond=4.35/c(531,427,294,189)*1/5
  RateDirectTreatSecond=c(0,0,0,0)
})
#Resistance parameters
# p2["alpha"]=2
# p2["alpha2"]=2

p2["alpha"]=2.64
p2["alpha2"]=4.9

##############################################################################################################################
#fix treat_dtg
dtg_women1=0.175 #above 50
dtg_women2=0.63 #above 50 or using contraception
treat_dtgb=c(p_w_start=0,p_w_switch=0,rate_first=0,rate_switch=0,p_treat=0,p_supp=0,p_fail=0,delay_first=0,delay_switch=0)
treat_dtg1=c(p_w_start=0,p_w_switch=0,rate_first=1,rate_switch=0,p_treat=0,p_supp=0,p_fail=0,delay_first=0,delay_switch=0)
treat_dtg2=c(p_w_start=dtg_women1,p_w_switch=0,rate_first=1,rate_switch=0,p_treat=0,p_supp=0,p_fail=0,delay_first=0,delay_switch=0)
treat_dtg3=c(p_w_start=dtg_women2,p_w_switch=0,rate_first=1,rate_switch=0,p_treat=0,p_supp=0,p_fail=0,delay_first=0,delay_switch=0)
treat_dtg4=c(p_w_start=1,p_w_switch=0,rate_first=1,rate_switch=0,p_treat=0,p_supp=0,p_fail=0,delay_first=0,delay_switch=0)
#switch also treat and suppressed
treat_dtg5=c(p_w_start=0,p_w_switch=0,rate_first=1,rate_switch=1,p_treat=1,p_supp=1,p_fail=1,delay_first=0,delay_switch=0)
treat_dtg6=c(p_w_start=dtg_women1,p_w_switch=dtg_women1,rate_first=1,rate_switch=1,p_treat=1,p_supp=1,p_fail=1,delay_first=0,delay_switch=0)
treat_dtg7=c(p_w_start=dtg_women2,p_w_switch=dtg_women2,rate_first=1,rate_switch=1,p_treat=1,p_supp=1,p_fail=1,delay_first=0,delay_switch=0)
treat_dtg8=c(p_w_start=1,p_w_switch=1,rate_first=1,rate_switch=1,p_treat=1,p_supp=1,p_fail=1,delay_first=0,delay_switch=0)
treat_dtg9=c(p_w_start=0.99,p_w_switch=0.99,rate_first=1,rate_switch=1,p_treat=1,p_supp=1,p_fail=1,delay_first=0,delay_switch=0)
treat_dtg_mat=rbind(treat_dtgb,treat_dtg1,treat_dtg2,treat_dtg3,treat_dtg4,treat_dtg5,treat_dtg6,treat_dtg7,treat_dtg8,treat_dtg9)
treat_dtg_mat[,"p_treat"]=0
##############################################################################################################################
#load function, run from 2005 to 2018
#Source and figures
figures_2_dtg_update=function(data,cex_t=1.5,cex_l=1.5,width=4,i=2,new_inf=FALSE){
  #Time axis
  tmax=dim(data)[1]
  time_start_abs=2005
  time_start_diff=2005+11/12
  time_stop_abs=time_start_abs+(tmax-1)/12
  time_stop_diff=time_start_diff+(tmax-1)/12-1
  
  #Data
  if(dim(data)[2]==250){
    new_infected=rollapply(diff(data[1:tmax,241]),12,sum)
    death_people=rollapply(diff(data[1:tmax,242]),12,sum)
    undiag_people=rowSums(data[seq(1,tmax,1),select_15(1,1:4,1:2,1:2)])
    treat_people=rowSums(data[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,1:2)])/rowSums(data[seq(1,tmax,1),select_15(1:15,1:4,1:2,1:2)])
    resis_people=rowSums(data[seq(1,tmax,1),select_15(c(5,15),1:4,2,1:2)])/rowSums(data[seq(1,tmax,1),select_15(c(5,15),1:4,1:2,1:2)])
    
    trandr_index=i
    if(i==2){trandr_index=c(i,12)}
    trandr=rowSums(data[seq(1,tmax,1),select_15(trandr_index,1:4,2,1:2)])/rowSums(data[seq(1,tmax,1),select_15(trandr_index,1:4,1:2,1:2)])
    if(new_inf){
      trandr=c(rep(NA,12),rollapply(diff(data[1:tmax,246]),12,sum)/rollapply(diff(data[1:tmax,241]),12,sum))
    }
  }else{
    new_infected=rollapply(diff(data[1:tmax,289]),12,sum)
    death_people=rollapply(diff(data[1:tmax,290]),12,sum)
    undiag_people=rowSums(data[seq(1,tmax,1),select_18(1,1:4,1:2,1:2)])
    treat_people=rowSums(data[seq(1,tmax,1),select_18(c(3:11,13:18),1:4,1:2,1:2)])/rowSums(data[seq(1,tmax,1),select_18(1:18,1:4,1:2,1:2)])
    resis_people=rowSums(data[seq(1,tmax,1),select_18(c(5,15,18),1:4,2,1:2)])/rowSums(data[seq(1,tmax,1),select_18(c(5,15,18),1:4,1:2,1:2)])
    
    trandr_index=i
    if(i==2){trandr_index=c(i,12)}
    trandr=rowSums(data[seq(1,tmax,1),select_18(trandr_index,1:4,2,1:2)])/rowSums(data[seq(1,tmax,1),select_18(trandr_index,1:4,1:2,1:2)])
    if(new_inf){
      trandr=c(rep(NA,12),rollapply(diff(data[1:tmax,246]),12,sum)/rollapply(diff(data[1:tmax,241]),12,sum))
    }
  }
  
  #Thembisa estimates
  par(mfrow = c(1,1),lwd=2,oma=c(0,0.1,0,0))
  m <- matrix(c(1,2,3,4,5,6,7,7),nrow = 4,ncol = 2,byrow = TRUE)
  graphics::layout(mat = m,heights = c(0.3,0.3,0.35,0.05),widths=c(0.5,0.5))
  
  #Plot1
  par(mai=c(0.2,0.8,0.5,0.4))
  plot(seq(2005+11/12,2015+11/12,length.out=11),newinf_value,xlim=c(time_start_abs,time_stop_diff),ylim=c(0,max(c(newinf_value,new_infected))),xaxt='n',xlab="",ylab="",main="A. New infected",pch=1,col="red",cex.main=2.5,cex.axis=2.2,cex.lab=2.2,cex=2,lwd=4)
  lines(seq(time_start_diff,time_stop_diff,1/12),new_infected,type="l",lwd=width)
  mtext(text="1000 people",side=2,line=4,cex=cex_t)
  
  #Plot2
  par(mai=c(0.2,0.7,0.5,0.4))
  plot(2005.5:2015.5,undiag_value,xlim=c(time_start_abs,time_stop_diff),ylim=c(0,max(undiag_value,undiag_people)),xaxt='n',xlab="",ylab="",main="B. Undiagnosed people",pch=1,col="red",cex.main=2.5,cex.axis=2.2,cex.lab=2.2,cex=2,lwd=4)
  lines(seq(time_start_abs,time_stop_abs,1/12),undiag_people,type="l",lwd=width)
  mtext(text="1000 people",side=2,line=4,cex=cex_t)
  
  #Plot3
  par(mai=c(0.2,0.8,0.5,0.4))
  plot(seq(2005+11/12,2015+11/12,length.out=11),death_value,xlim=c(time_start_abs,time_stop_diff),ylim=c(0,max(death_value,death_people)),xaxt='n',xlab="",ylab="",main="C. AIDS-related deaths",pch=1,col="red",cex.main=2.5,cex.axis=2.2,cex.lab=2.2,cex=2,lwd=4)
  lines(seq(time_start_diff,time_stop_diff,1/12),death_people,type="l",col="black",lwd=width)
  mtext(text="1000 people",side=2,line=4,cex=cex_t)
  
  #Plot4
  par(mai=c(0.2,0.7,0.5,0.4))
  plot(2005.5:2015.5,c(treat_percent*100),xlim=c(time_start_abs,time_stop_diff),ylim=c(0,max(treat_percent,treat_people)*100),xaxt='n',xlab="",ylab="",main="D. ART coverage",pch=1,col="red",cex.main=2.5,cex.axis=2.2,cex.lab=2.2,cex=2,lwd=4)
  lines(seq(time_start_abs,time_stop_abs,1/12),treat_people*100,type="l",col="black",lwd=width)
  mtext(text="percent",side=2,line=4,cex=cex_t)
  
  #Plot5
  par(mai=c(0.8,0.8,0.5,0.4))
  plot(2005.5:2016.5,c(rep(NA,5),80,NA,NA,NA,95,NA,NA),xlim=c(time_start_abs,time_stop_diff),ylim=c(min(resis_people[-1])*100,100),xaxt='n',xlab="",ylab="",main="E. ADR (Failing 1st-line)",pch=1,col="blue",cex.main=2.5,cex.axis=2.5,cex.lab=2.2,cex=2,lwd=4)
  lines(seq(time_start_abs,time_stop_abs,1/12),c(rep(NA,12),resis_people[13:length(resis_people)]*100),type="l",col="black",lwd=width)
  mtext(text="Calendar year",side=1,line=4,cex=cex_t)
  mtext(text="percent",side=2,line=4,cex=cex_t)
  axis(side=1,tick=TRUE,line=1,cex.axis=2.5,lwd=2,col="white")
  axis(side=1,labels=FALSE,tick=TRUE,line=0.5,cex.axis=3,lwd=2)
  
  #Plot6
  par(mai=c(0.8,0.7,0.5,0.4))
  #plot(2005.5:2016.5,c(1,1.25,3,2.5,3.5,4,4.75,6,8,10,9.5,8.5),ylim=c(min(trandr,trandr2,trandr3,trandr4)*100,max(trandr,trandr2,trandr3,trandr4)*100),xlab="Calendar year",ylab="percent",main="F. TDR Prevalence (diagnosed people)",pch=1,col="blue",cex.main=2,cex.axis=1.5,cex.lab=1.6,cex=2,lwd=4)
  plot(c(2005,2009,2011,2014),c(1.7,2.1,5.2,9),xlim=c(time_start_abs,time_stop_diff),ylim=c(0,max(0.1,trandr,na.rm=TRUE)*100),xaxt='n',xlab="",ylab="",main="F. TDR Prevalence (diagnosed people)",pch=1,col="blue",cex.main=2.5,cex.axis=2.5,cex.lab=2.5,cex=2,lwd=4)
  lines(seq(time_start_abs,time_stop_abs,1/12),trandr*100,type="l",col="black",lwd=width)
  mtext(text="Calendar year",side=1,line=4,cex=cex_t)
  mtext(text="percent",side=2,line=4,cex=cex_t)
  axis(side=1,tick=TRUE,line=1,cex.axis=2.5,lwd=2,tck=-0.1,col="white")
  axis(side=1,labels=FALSE,tick=TRUE,line=0.5,cex.axis=3,lwd=2,tck=-0.05)
  
  #Legend
  par(lwd=1,mai=c(0,0,0,0))
  plot.new()
  
  legend(x="left", inset =0,
         c("Simulated by the model","Estimated data (WHO)","Data from previous studies"),
         lty=c(1,NA,NA),pch=c(NA,"o","o"), lwd=c(2,3,3), col=c("black","red","blue"), box.col=NA,horiz=TRUE,cex=cex_l,text.width = c(0.15,0.15,0.15))
}
figures_dtg_4_update=function(list_data,scen=c("Scen. 1", "Scen. 2", "Scen. 3","Scen. 4"),cex_t=1.5,cex_l=2,width=4,label="Simulations",i=2,men_women=c(1:2),new_inf=FALSE,text_y="Level of NNRTI PDR among diagnosed individuals (%)    ",title_1="1. DTG as 1st-line",title_2="2. DTG as 1st-line and switch"){
  #Datasets 1,2,3,4
  data=list_data[[1]]
  data1=list_data[[2]]
  data2=list_data[[3]]
  data3=list_data[[4]]
  data4=list_data[[5]]
  data5=list_data[[6]]
  data6=list_data[[7]]
  data7=list_data[[8]]
  data8=list_data[[9]]
  
  ###################################################################################################################
  #Thembisa estimates
  par(mfrow = c(1,1),lwd=2,oma=c(0,0.1,0,0))
  
  m <- matrix(c(1,2,3,4),nrow = 4,ncol = 1,byrow = FALSE)
  graphics::layout(mat = m,heights = c(0.4,0.5,0.05,0.05),widths=c(0.5,0.5))
  
  #Time axis
  tmax=dim(data)[1]
  time_start_abs=2005
  time_start_diff=2005+11/12
  time_stop_abs=time_start_abs+(tmax-1)/12
  time_stop_diff=time_start_diff+(tmax-1)/12-1
  #Data
  trandr=c(rep(NA,11),rollapply(diff(data[1:tmax,244]),12,sum)/rollapply(diff(data[1:tmax,241]),12,sum))
  
  trandr_index=1
  if(i==2){trandr_index=c(i,12)}
  trandr=rowSums(data[seq(1,tmax,1),select_15(trandr_index,1:4,2,men_women)])/rowSums(data[seq(1,tmax,1),select_15(trandr_index,1:4,1:2,men_women)])
  trandr1=rowSums(data1[seq(1,tmax,1),select_15(trandr_index,1:4,2,men_women)])/rowSums(data1[seq(1,tmax,1),select_15(trandr_index,1:4,1:2,men_women)])
  trandr2=rowSums(data2[seq(1,tmax,1),select_15(trandr_index,1:4,2,men_women)])/rowSums(data2[seq(1,tmax,1),select_15(trandr_index,1:4,1:2,men_women)])
  trandr3=rowSums(data3[seq(1,tmax,1),select_15(trandr_index,1:4,2,men_women)])/rowSums(data3[seq(1,tmax,1),select_15(trandr_index,1:4,1:2,men_women)])
  trandr4=rowSums(data4[seq(1,tmax,1),select_15(trandr_index,1:4,2,men_women)])/rowSums(data4[seq(1,tmax,1),select_15(trandr_index,1:4,1:2,men_women)])
  trandr5=rowSums(data5[seq(1,tmax,1),select_15(trandr_index,1:4,2,men_women)])/rowSums(data5[seq(1,tmax,1),select_15(trandr_index,1:4,1:2,men_women)])
  trandr6=rowSums(data6[seq(1,tmax,1),select_15(trandr_index,1:4,2,men_women)])/rowSums(data6[seq(1,tmax,1),select_15(trandr_index,1:4,1:2,men_women)])
  trandr7=rowSums(data7[seq(1,tmax,1),select_15(trandr_index,1:4,2,men_women)])/rowSums(data7[seq(1,tmax,1),select_15(trandr_index,1:4,1:2,men_women)])
  trandr8=rowSums(data8[seq(1,tmax,1),select_15(trandr_index,1:4,2,men_women)])/rowSums(data8[seq(1,tmax,1),select_15(trandr_index,1:4,1:2,men_women)])
  
  if(new_inf){
    trandr=c(rep(NA,12),rollapply(diff(data[1:tmax,246]),12,sum)/rollapply(diff(data[1:tmax,241]),12,sum))
    trandr1=c(rep(NA,12),rollapply(diff(data1[1:tmax,246]),12,sum)/rollapply(diff(data1[1:tmax,241]),12,sum))
    trandr2=c(rep(NA,12),rollapply(diff(data2[1:tmax,246]),12,sum)/rollapply(diff(data2[1:tmax,241]),12,sum))
    trandr3=c(rep(NA,12),rollapply(diff(data3[1:tmax,246]),12,sum)/rollapply(diff(data3[1:tmax,241]),12,sum))
    trandr4=c(rep(NA,12),rollapply(diff(data4[1:tmax,246]),12,sum)/rollapply(diff(data4[1:tmax,241]),12,sum))
    trandr5=c(rep(NA,12),rollapply(diff(data5[1:tmax,246]),12,sum)/rollapply(diff(data5[1:tmax,241]),12,sum))
    trandr6=c(rep(NA,12),rollapply(diff(data6[1:tmax,246]),12,sum)/rollapply(diff(data6[1:tmax,241]),12,sum))
    trandr7=c(rep(NA,12),rollapply(diff(data7[1:tmax,246]),12,sum)/rollapply(diff(data7[1:tmax,241]),12,sum))
    trandr8=c(rep(NA,12),rollapply(diff(data8[1:tmax,246]),12,sum)/rollapply(diff(data8[1:tmax,241]),12,sum))
  }

  if(new_inf){text_y="Level of NNRTI PDR among newly infected individuals (%)"}
  #Plot1
  par(mai=c(0.2,0.8,0.5,0.4))
  #plot(2005.5:2016.5,c(1,1.25,3,2.5,3.5,4,4.75,6,8,10,9.5,8.5),ylim=c(min(trandr,trandr2,trandr3,trandr4)*100,max(trandr,trandr2,trandr3,trandr4)*100),xlab="Calendar year",ylab="percent",main="F. TDR Prevalence (diagnosed people)",pch=1,col="blue",cex.main=2,cex.axis=1.5,cex.lab=1.6,cex=2,lwd=4)
  #plot(c(2005,2009,2011,2014),c(1.7,2.1,5.2,9),xlim=c(time_start_abs,time_stop_diff),ylim=c(min(trandr,trandr2,trandr3,trandr4)*100,max(trandr,trandr2,trandr3,trandr4)*100),xaxt='n',xlab="",ylab="",main="1. PDR level (diagnosed people)",pch=1,col="blue",cex.main=2,cex.axis=2.2,cex.lab=2.5,cex=2,lwd=4)
  plot(seq(time_start_abs,time_stop_abs,1/12),trandr*100,type="l",col="black",lwd=width,xlim=c(time_start_abs,time_stop_diff),ylim=c(min(trandr,trandr1,trandr2,trandr3,trandr4,na.rm=TRUE)*100,max(trandr,trandr1,trandr2,trandr3,trandr4,na.rm=TRUE)*100),xaxt='n',xlab="",ylab="",main=ifelse(length(men_women)>1,title_1,paste(title_1,paste(c("men","women")[men_women],collapse=" and "),sep=", ")),pch=1,cex.main=2,cex.axis=2.2,cex.lab=2.5,cex=2)
  lines(seq(time_start_abs,time_stop_abs,1/12),trandr1*100,type="l",col="violet",lwd=width)
  lines(seq(time_start_abs,time_stop_abs,1/12),trandr2*100,type="l",col="blue",lwd=width)
  lines(seq(time_start_abs,time_stop_abs,1/12),trandr3*100,type="l",col="red",lwd=width)
  lines(seq(time_start_abs,time_stop_abs,1/12),trandr4*100,type="l",col="orange",lwd=width)
  abline(v=2018,lty=2,col="red")
  text(x=2015,y=max(trandr,trandr5,trandr6,trandr7,trandr8,na.rm=TRUE)*100*3/4,labels="DTG \n introduction",cex=2,col="red")
  mtext(text=text_y,side=2,line=4,cex=cex_t,adj=1)
  
  
  #Plot2
  par(mai=c(0.8,0.8,0.5,0.4))
  #plot(2005.5:2016.5,c(1,1.25,3,2.5,3.5,4,4.75,6,8,10,9.5,8.5),ylim=c(min(trandr,trandr2,trandr3,trandr4)*100,max(trandr,trandr2,trandr3,trandr4)*100),xlab="Calendar year",ylab="percent",main="F. TDR Prevalence (diagnosed people)",pch=1,col="blue",cex.main=2,cex.axis=1.5,cex.lab=1.6,cex=2,lwd=4)
  #plot(c(2005,2009,2011,2014),c(1.7,2.1,5.2,9),xlim=c(time_start_abs,time_stop_diff),ylim=c(min(trandr,trandr2,trandr3,trandr4)*100,max(trandr,trandr2,trandr3,trandr4)*100),xaxt='n',xlab="",ylab="",main="2. PDR level (diagnosed people)",pch=1,col="blue",cex.main=2,cex.axis=2.2,cex.lab=2.5,cex=2,lwd=4)
  plot(seq(time_start_abs,time_stop_abs,1/12),trandr*100,type="l",col="black",lwd=width,xlim=c(time_start_abs,time_stop_diff),ylim=c(min(trandr,trandr5,trandr6,trandr7,trandr8,na.rm=TRUE)*100,max(trandr,trandr5,trandr6,trandr7,trandr8,na.rm=TRUE)*100),xaxt='n',xlab="",ylab="",main=ifelse(length(men_women)>1,title_2,paste(title_2,paste(c("men","women")[men_women],collapse=" and "),sep=", ")),pch=1,cex.main=2,cex.axis=2.2,cex.lab=2.5,cex=2)
  lines(seq(time_start_abs,time_stop_abs,1/12),trandr5*100,type="l",col="violet",lwd=width)
  lines(seq(time_start_abs,time_stop_abs,1/12),trandr6*100,type="l",col="blue",lwd=width)
  lines(seq(time_start_abs,time_stop_abs,1/12),trandr7*100,type="l",col="red",lwd=width)
  lines(seq(time_start_abs,time_stop_abs,1/12),trandr8*100,type="l",col="orange",lwd=width)
  abline(v=2018,lty=2,col="red")
  text(x=2015,y=max(trandr,trandr5,trandr6,trandr7,trandr8,na.rm=TRUE)*100*3/4,labels="DTG \n introduction",cex=2,col="red")
  mtext(text="Calendar year",side=1,line=4,cex=cex_t)
  axis(side=1,tick=TRUE,line=1,cex.axis=2.5,lwd=2,col="white")
  axis(side=1,labels=FALSE,tick=TRUE,line=0.5,cex.axis=3,lwd=2)
  
  
  
  
  #Legend
  par(lwd=1,mai=c(0,0,0,0))
  plot.new()
  legend(x="left", inset =0,
         c("Baseline model",""),
         lty=c(1,NA),pch=c(NA,NA), lwd=c(2,3), col=c("black","blue"), box.col=NA,horiz=TRUE,cex=cex_l,text.width = c(0.15,0.15),seg.len=c(1.2,1.2))
  
  par(lwd=1,mai=c(0,0,0,0))
  plot.new()
  legend(x="left", inset =0,
         c(scen),
         lty=c(1,1,1,1),pch=c(NA,NA,NA,NA), lwd=c(2,2,2,2), col=c("violet","blue","red","orange"), box.col=NA,horiz=TRUE,cex=cex_l,text.width = c(0.22,0.10,0.15,0.17),seg.len=c(1,1,1,1))
  
}
figures_dtg_n_res=function(list_data,scen=c("Scen. 1", "Scen. 2", "Scen. 3","Scen. 4"),cex_t=1.5,cex_l=2,width=4,label="Simulations",i=2,men_women=c(1:2),new_inf=FALSE,text_y="Cumulative number of new NNRTI resistant cases (1000 people)    ",title_1="1. DTG as 1st-line",title_2="2. DTG as 1st-line and switch"){
  #Datasets 1,2,3,4
  data=list_data[[1]]
  data1=list_data[[2]]
  data2=list_data[[3]]
  data3=list_data[[4]]
  data4=list_data[[5]]
  data5=list_data[[6]]
  data6=list_data[[7]]
  data7=list_data[[8]]
  data8=list_data[[9]]
  
  ###################################################################################################################
  #Thembisa estimates
  par(mfrow = c(1,1),lwd=2,oma=c(0,0.1,0,0))
  
  m <- matrix(c(1,2,3,4),nrow = 4,ncol = 1,byrow = FALSE)
  graphics::layout(mat = m,heights = c(0.4,0.5,0.05,0.05),widths=c(0.5,0.5))
  
  #Time axis
  tmax=dim(data)[1]
  time_start_abs=2005
  time_start_diff=2005+11/12
  time_stop_abs=time_start_abs+(tmax-1)/12
  time_stop_diff=time_start_diff+(tmax-1)/12-1
  #Data
  trandr=c(rep(NA,11),rollapply(diff(data[1:tmax,244]),12,sum)/rollapply(diff(data[1:tmax,241]),12,sum))
  
  trandr_index=1
  if(i==2){trandr_index=c(i,12)}
  trandr=data[,243]/1000
  trandr1=data1[,243]/1000
  trandr2=data2[,243]/1000
  trandr3=data3[,243]/1000
  trandr4=data4[,243]/1000
  trandr5=data5[,243]/1000
  trandr6=data6[,243]/1000
  trandr7=data7[,243]/1000
  trandr8=data8[,243]/1000
  
  #Plot1
  par(mai=c(0.2,0.8,0.5,0.4))
  #plot(2005.5:2016.5,c(1,1.25,3,2.5,3.5,4,4.75,6,8,10,9.5,8.5),ylim=c(min(trandr,trandr2,trandr3,trandr4)*100,max(trandr,trandr2,trandr3,trandr4)*100),xlab="Calendar year",ylab="percent",main="F. TDR Prevalence (diagnosed people)",pch=1,col="blue",cex.main=2,cex.axis=1.5,cex.lab=1.6,cex=2,lwd=4)
  #plot(c(2005,2009,2011,2014),c(1.7,2.1,5.2,9),xlim=c(time_start_abs,time_stop_diff),ylim=c(min(trandr,trandr2,trandr3,trandr4)*100,max(trandr,trandr2,trandr3,trandr4)*100),xaxt='n',xlab="",ylab="",main="1. PDR level (diagnosed people)",pch=1,col="blue",cex.main=2,cex.axis=2.2,cex.lab=2.5,cex=2,lwd=4)
  plot(seq(time_start_abs,time_stop_abs,1/12),trandr*100,type="l",col="black",lwd=width,xlim=c(time_start_abs,time_stop_diff),ylim=c(min(trandr,trandr1,trandr2,trandr3,trandr4,na.rm=TRUE)*100,max(trandr,trandr1,trandr2,trandr3,trandr4,na.rm=TRUE)*100),xaxt='n',xlab="",ylab="",main=ifelse(length(men_women)>1,title_1,paste(title_1,paste(c("men","women")[men_women],collapse=" and "),sep=", ")),pch=1,cex.main=2,cex.axis=2.2,cex.lab=2.5,cex=2)
  lines(seq(time_start_abs,time_stop_abs,1/12),trandr1*100,type="l",col="violet",lwd=width)
  lines(seq(time_start_abs,time_stop_abs,1/12),trandr2*100,type="l",col="blue",lwd=width)
  lines(seq(time_start_abs,time_stop_abs,1/12),trandr3*100,type="l",col="red",lwd=width)
  lines(seq(time_start_abs,time_stop_abs,1/12),trandr4*100,type="l",col="orange",lwd=width)
  abline(v=2018,lty=2,col="red")
  text(x=2015,y=max(trandr,trandr5,trandr6,trandr7,trandr8,na.rm=TRUE)*100*3/4,labels="DTG \n introduction",cex=2,col="red")
  mtext(text=text_y,side=2,line=4,cex=cex_t,adj=1)
  
  
  #Plot2
  par(mai=c(0.8,0.8,0.5,0.4))
  #plot(2005.5:2016.5,c(1,1.25,3,2.5,3.5,4,4.75,6,8,10,9.5,8.5),ylim=c(min(trandr,trandr2,trandr3,trandr4)*100,max(trandr,trandr2,trandr3,trandr4)*100),xlab="Calendar year",ylab="percent",main="F. TDR Prevalence (diagnosed people)",pch=1,col="blue",cex.main=2,cex.axis=1.5,cex.lab=1.6,cex=2,lwd=4)
  #plot(c(2005,2009,2011,2014),c(1.7,2.1,5.2,9),xlim=c(time_start_abs,time_stop_diff),ylim=c(min(trandr,trandr2,trandr3,trandr4)*100,max(trandr,trandr2,trandr3,trandr4)*100),xaxt='n',xlab="",ylab="",main="2. PDR level (diagnosed people)",pch=1,col="blue",cex.main=2,cex.axis=2.2,cex.lab=2.5,cex=2,lwd=4)
  plot(seq(time_start_abs,time_stop_abs,1/12),trandr*100,type="l",col="black",lwd=width,xlim=c(time_start_abs,time_stop_diff),ylim=c(min(trandr,trandr5,trandr6,trandr7,trandr8,na.rm=TRUE)*100,max(trandr,trandr5,trandr6,trandr7,trandr8,na.rm=TRUE)*100),xaxt='n',xlab="",ylab="",main=ifelse(length(men_women)>1,title_2,paste(title_2,paste(c("men","women")[men_women],collapse=" and "),sep=", ")),pch=1,cex.main=2,cex.axis=2.2,cex.lab=2.5,cex=2)
  lines(seq(time_start_abs,time_stop_abs,1/12),trandr5*100,type="l",col="violet",lwd=width)
  lines(seq(time_start_abs,time_stop_abs,1/12),trandr6*100,type="l",col="blue",lwd=width)
  lines(seq(time_start_abs,time_stop_abs,1/12),trandr7*100,type="l",col="red",lwd=width)
  lines(seq(time_start_abs,time_stop_abs,1/12),trandr8*100,type="l",col="orange",lwd=width)
  abline(v=2018,lty=2,col="red")
  text(x=2015,y=max(trandr,trandr5,trandr6,trandr7,trandr8,na.rm=TRUE)*100*3/4,labels="DTG \n introduction",cex=2,col="red")
  mtext(text="Calendar year",side=1,line=4,cex=cex_t)
  axis(side=1,tick=TRUE,line=1,cex.axis=2.5,lwd=2,col="white")
  axis(side=1,labels=FALSE,tick=TRUE,line=0.5,cex.axis=3,lwd=2)
  
  
  
  
  #Legend
  par(lwd=1,mai=c(0,0,0,0))
  plot.new()
  legend(x="left", inset =0,
         c("Baseline model",""),
         lty=c(1,NA),pch=c(NA,NA), lwd=c(2,3), col=c("black","blue"), box.col=NA,horiz=TRUE,cex=cex_l,text.width = c(0.15,0.15),seg.len=c(1.2,1.2))
  
  par(lwd=1,mai=c(0,0,0,0))
  plot.new()
  legend(x="left", inset =0,
         c(scen),
         lty=c(1,1,1,1),pch=c(NA,NA,NA,NA), lwd=c(2,2,2,2), col=c("violet","blue","red","orange"), box.col=NA,horiz=TRUE,cex=cex_l,text.width = c(0.22,0.10,0.15,0.17),seg.len=c(1,1,1,1))
  
}



source("C:/Users/ahauser/Documents/Step2/dtg/parmin_dtg.R")
sourceCpp("C:/Users/ahauser/Documents/Step2/dtg/solvdiff_cpp_dtg_2_2005.cpp")
xstart_dtg_2_2005=xstart_f_dtg_2(p2,treat_dtgb)
theta=p1
#15 compartments
data=my_fun10_solver2_dtg_2(xstart_dtg_2_2005,theta,treat_dtgb,params,p2)
#figures_2_dtg_update(data)
save(data,file="C:/Users/ahauser/Documents/Step2/R_results/data_dtg_2005_2019.RData")
#load function, run from 2005 to 2018
#Counterfactual scenarios from 2018 to 2038
sourceCpp("C:/Users/ahauser/Documents/Step2/dtg/solvdiff_cpp_dtg_2_2019.cpp")
xstart_dtg_2_2018=data[dim(data)[1],]
dtg_list=list()
for(i in 1:10){
  xstart=xstart_dtg_2_2018
  
  xstart[select_15(12:15,1:4,1:2,2)]=(1-treat_dtg_mat[i,1])*(xstart_dtg_2_2018[select_15(2:5,1:4,1:2,2)]+xstart_dtg_2_2018[select_15(12:15,1:4,1:2,2)])
  xstart[select_15(2:5,1:4,1:2,2)]=treat_dtg_mat[i,1]*(xstart_dtg_2_2018[select_15(2:5,1:4,1:2,2)]+xstart_dtg_2_2018[select_15(12:15,1:4,1:2,2)])
  xstart=c(xstart[1:240],xstart_dtg_2_2018[241:250])
  
  data1=my_fun10_solver2_dtg_2(xstart,theta,treat_dtg_mat[i,],params,p2)
  rownames(data1)=1:(dim(data1)[1])
  data1=rbind(data[-dim(data)[1],],data1)
  dtg_list[[i]]<-data1
  print(i)
}
dtg_ref<-dtg_list

data<-dtg_ref[[9]]
tmax=dim(data)[1]
rowSums(data[seq(1,tmax,1),select_15(c(2,12),1:4,2,1:2)])/rowSums(data[seq(1,tmax,1),select_15(c(2,12),1:4,1:2,1:2)])
data<-dtg_ref[[1]]
tmax=dim(data)[1]
rowSums(data[seq(1,tmax,1),select_15(c(2,12),1:4,2,1:2)])/rowSums(data[seq(1,tmax,1),select_15(c(2,12),1:4,1:2,1:2)])



figures_dtg_4_update(dtg_ref)
figures_dtg_n_res(dtg_list)

d<-dtg_list[[2]]
(d[,241])[seq(1,421,12)]
d<-dtg_list[[5]]
(d[,241])[seq(1,421,12)]
d<-dtg_list[[1]]
(d[,243])[seq(1,421,12)]
d<-dtg_list[[5]]
(d[,243])[seq(1,421,12)]

#Unsupp level
d<-dtg_list[[2]]
rowSums(d[seq(1,421,12),select_15(c(1,2,12),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(1,2,12,3,5,6,8,9,11,13,15),1:4,1:2,1:2)])
rowSums(d[seq(1,421,12),select_15(c(3,5,13,15),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(1,2,12,3,5,6,8,9,11,13,15),1:4,1:2,1:2)])
rowSums(d[seq(1,421,12),select_15(c(6,8),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(1,2,12,3,5,6,8,9,11,13,15),1:4,1:2,1:2)])
rowSums(d[seq(1,421,12),select_15(c(9,11),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(1,2,12,3,5,6,8,9,11,13,15),1:4,1:2,1:2)])
#Unsupp level
d<-dtg_list[[5]]
rowSums(d[seq(1,421,12),select_15(c(1,2,12),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(1,2,12,3,5,6,8,9,11,13,15),1:4,1:2,1:2)])
rowSums(d[seq(1,421,12),select_15(c(3,5,13,15),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(1,2,12,3,5,6,8,9,11,13,15),1:4,1:2,1:2)])
rowSums(d[seq(1,421,12),select_15(c(6,8),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(1,2,12,3,5,6,8,9,11,13,15),1:4,1:2,1:2)])
rowSums(d[seq(1,421,12),select_15(c(9,11),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(1,2,12,3,5,6,8,9,11,13,15),1:4,1:2,1:2)])



#People on ART
d<-dtg_list_v5_noint[[2]]
rowSums(d[seq(1,421,12),select_15(c(3:11,13:15),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(1:15),1:4,1:2,1:2)])
d<-dtg_list_v5_int2[[2]]
rowSums(d[seq(1,421,12),select_15(c(3:11,13:15),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(1:15),1:4,1:2,1:2)])


#Unsupp level
d<-dtg_list[[2]]
rowSums(d[seq(1,421,12),select_15(c(3,5,13,15),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(3,5,6,8,9,11,13,15),1:4,1:2,1:2)])
rowSums(d[seq(1,421,12),select_15(c(9,11),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(3,5,6,8,9,11,13,15),1:4,1:2,1:2)])

d<-dtg_list[[5]]
rowSums(d[seq(1,421,12),select_15(c(3,5,6,8,9,11,13,15),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(3:11,13:15),1:4,1:2,1:2)])



#Unsupp level
d<-dtg_list[[2]]
rowSums(d[seq(1,421,12),select_15(c(3,5,6,8,9,11,13,15),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(3:11,13:15),1:4,1:2,1:2)])
rowSums(d[seq(1,421,12),select_15(c(3,5,13,15),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(3:5,13:15),1:4,1:2,1:2)])
rowSums(d[seq(1,421,12),select_15(c(9,11),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(9:11),1:4,1:2,1:2)])

d<-dtg_list[[5]]
rowSums(d[seq(1,421,12),select_15(c(3,5,6,8,9,11,13,15),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(3:11,13:15),1:4,1:2,1:2)])

#Res among treated people
d<-dtg_list[[2]]
rowSums(d[seq(1,421,12),select_15(c(3:11,13:15),1:4,2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(3:11,13:15),1:4,1:2,1:2)])
d<-dtg_list[[5]]
rowSums(d[seq(1,421,12),select_15(c(3:11,13:15),1:4,2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(3:11,13:15),1:4,1:2,1:2)])

#Res among unsup people
d<-dtg_list[[2]]
rowSums(d[seq(1,421,12),select_15(c(3,5,6,8,9,11,13,15),1:4,2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(3,5,6,8,9,11,13,15),1:4,1:2,1:2)])
d<-dtg_list[[5]]
rowSums(d[seq(1,421,12),select_15(c(3,5,6,8,9,11,13,15),1:4,2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(3,5,6,8,9,11,13,15),1:4,1:2,1:2)])

#Res among NNRTI unsup people
d<-dtg_list[[2]]
rowSums(d[seq(1,421,12),select_15(c(3,5,13,15),1:4,2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(3,5,13,15),1:4,1:2,1:2)])
d<-dtg_list[[5]]
rowSums(d[seq(1,421,12),select_15(c(3,5,13,15),1:4,2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(3,5,13,15),1:4,1:2,1:2)])

#Res among DTG unsup people
d<-dtg_list[[2]]
rowSums(d[seq(1,421,12),select_15(c(9,11),1:4,2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(9,11),1:4,1:2,1:2)])
d<-dtg_list[[5]]
rowSums(d[seq(1,421,12),select_15(c(9,11),1:4,2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(9,11),1:4,1:2,1:2)])

#prop of people on DTG, all
d<-dtg_list[[2]]
rowSums(d[seq(1,421,12),select_15(c(9:11),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(3:11,13:15),1:4,1:2,1:2)])
d<-dtg_list[[5]]
rowSums(d[seq(1,421,12),select_15(c(9:11),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(3:11,13:15),1:4,1:2,1:2)])

#prop of people on NNRTI, all
d<-dtg_list[[2]]
rowSums(d[seq(1,421,12),select_15(c(3:5,13:15),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(3:11,13:15),1:4,1:2,1:2)])
d<-dtg_list[[5]]
rowSums(d[seq(1,421,12),select_15(c(3:5,13:15),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(3:11,13:15),1:4,1:2,1:2)])

#prop of people on PI, all
d<-dtg_list[[2]]
rowSums(d[seq(1,421,12),select_15(c(6:8),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(3:11,13:15),1:4,1:2,1:2)])
d<-dtg_list[[5]]
rowSums(d[seq(1,421,12),select_15(c(6:8),1:4,1:2,1:2)])/rowSums(d[seq(1,421,12),select_15(c(3:11,13:15),1:4,1:2,1:2)])


#prop of people on DTG, women
d<-dtg_list[[5]]
rowSums(d[seq(1,421,12),select_15(c(9:11),1:4,1:2,2)])/rowSums(d[seq(1,421,12),select_15(c(3:11,13:15),1:4,1:2,2)])


load(paste("C:/Users/ahauser/Documents/Step2/R_results/dtg_sens/dtg_sens13/data1_",1,"_","1",".RData",sep=""))

##############################################################################################################################
#Data frame
scen=c("0) all on NNRTI","1) all men on DTG","2) all men and women not of childbearing age on DTG",
       "3) all men and women not at risk of pregancy on DTG","4) all men and women on DTG")

#plot use of regimens according to scenarios, gender
art_use=data.frame(year=rep(seq(2005,2039,1/12),90),
                   ART_use=NA,
                   main_scen=rep(c("DTG as first-line", "DTG as first-line and switch"),each=409*45),
                   sub_scen=rep(rep(scen,2),each=409*9),
                   gender=rep(rep(c("men and women", "men", "women"),10),each=409*3),
                   art=rep(rep(c("NNRTI","DTG","PI"),30),each=409))

# art_use=expand.grid(year=seq(2005,2038,1/12),
#             sub_scen=scen,
#             main_scen=c("DTG as first-line", "DTG as first-line and switch"),
#             sub_scen=scen,
#             gender=c("men and women", "men", "women"),
#             art=c("NNRTI","DTG","PI"))


art_use$gender<-factor(art_use$gender,levels=c("men and women", "men", "women"))

#Results from R script (above)
dtg_list2=c(dtg_list[1:5],dtg_list[1],dtg_list[6:9])
for(i in 1:10){
  data_dtg<-dtg_list2[[i]]
  tmax=dim(data_dtg)[1]
  v_all1<-rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:5,13:15),1:4,1:2,1:2)])/rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,1:2)])
  v_men1<-rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:5,13:15),1:4,1:2,1)])/rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,1)])
  v_women1<-rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:5,13:15),1:4,1:2,2)])/rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,2)])
  v_all2<-rowSums(data_dtg[seq(1,tmax,1),select_15(c(9:11),1:4,1:2,1:2)])/rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,1:2)])
  v_men2<-rowSums(data_dtg[seq(1,tmax,1),select_15(c(9:11),1:4,1:2,1)])/rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,1)])
  v_women2<-rowSums(data_dtg[seq(1,tmax,1),select_15(c(9:11),1:4,1:2,2)])/rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,2)])
  v_all3<-rowSums(data_dtg[seq(1,tmax,1),select_15(c(6:8),1:4,1:2,1:2)])/rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,1:2)])
  v_men3<-rowSums(data_dtg[seq(1,tmax,1),select_15(c(6:8),1:4,1:2,1)])/rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,1)])
  v_women3<-rowSums(data_dtg[seq(1,tmax,1),select_15(6:8,1:4,1:2,2)])/rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,2)])
  
  
  v=c(v_all1,v_all2,v_all3,v_men1,v_men2,v_men3,v_women1,v_women2,v_women3)
  l1=1+(409*9)*(i-1)
  l2=(409*9)*i
  
  l=seq(l1,l2,1)
  art_use$ART_use[l]<-v
}


#Results from Ubelix, dtg_sens
art_use=data.frame(year=rep(seq(2005,2040,1/12),90),
                   ART_use=NA,
                   main_scen=rep(c("DTG as first-line", "DTG as first-line and switch"),each=421*45),
                   sub_scen=rep(rep(scen,2),each=421*9),
                   gender=rep(rep(c("men and women", "men", "women"),10),each=421*3),
                   art=rep(rep(c("NNRTI","DTG","PI"),30),each=421))

index_list=c(1:5,1,6:9)
for(i in 1:10){
  #With simulations from ub2.job
  load(paste("C:/Users/ahauser/Documents/Step2/R_results/dtg_sens/dtg_sens17/data1_",index_list[i],"_","1",".RData",sep=""))
  data_dtg=data1
  #With simulations from ub.job
  #load(paste("C:/Users/ahauser/Documents/Step2/R_results/dtg_sens/dtg_sens6/dtg_sens",index_list[i],".RData",sep=""))
  #data_dtg<-list_data[[1]]
  tmax=dim(data_dtg)[1]
  v_all1<-rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:5,13:15),1:4,1:2,1:2)])/rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,1:2)])
  v_men1<-rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:5,13:15),1:4,1:2,1)])/rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,1)])
  v_women1<-rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:5,13:15),1:4,1:2,2)])/rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,2)])
  v_all2<-rowSums(data_dtg[seq(1,tmax,1),select_15(c(9:11),1:4,1:2,1:2)])/rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,1:2)])
  v_men2<-rowSums(data_dtg[seq(1,tmax,1),select_15(c(9:11),1:4,1:2,1)])/rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,1)])
  v_women2<-rowSums(data_dtg[seq(1,tmax,1),select_15(c(9:11),1:4,1:2,2)])/rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,2)])
  v_all3<-rowSums(data_dtg[seq(1,tmax,1),select_15(c(6:8),1:4,1:2,1:2)])/rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,1:2)])
  v_men3<-rowSums(data_dtg[seq(1,tmax,1),select_15(c(6:8),1:4,1:2,1)])/rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,1)])
  v_women3<-rowSums(data_dtg[seq(1,tmax,1),select_15(6:8,1:4,1:2,2)])/rowSums(data_dtg[seq(1,tmax,1),select_15(c(3:11,13:15),1:4,1:2,2)])
  
  v=c(v_all1,v_all2,v_all3,v_men1,v_men2,v_men3,v_women1,v_women2,v_women3)
  # l1=1+(397*9)*(i-1)
  # l2=(397*9)*i
  l1=1+(421*9)*(i-1)
  l2=(421*9)*i
  
  l=seq(l1,l2,1)
  art_use$ART_use[l]<-v
  print(i)
}


##############################################################################################################################
#Figure
#Plot with split by gender
art_use2 <- art_use[art_use$gender!="men and women",]
art_use2$gender<-factor(art_use2$gender,levels=c("men", "women"))
ggplot(data=art_use2, mapping=aes(x=year,y=ART_use,color=sub_scen))+
  geom_line(size=1.1,aes(linetype=art))+
  theme_bw()+
  facet_grid(main_scen+gender~sub_scen,switch="y")+
  theme(legend.position="bottom")

#Plot no split by gender
color=c("black","blue","violet","orange","red")

art_use3 <- art_use[art_use$gender=="men and women",]
art_use3 <- art_use3[art_use3$art %in% c("NNRTI","DTG","PI"),]
art_use3$art<-factor(art_use3$art,levels=c("NNRTI", "DTG","PI"))
art_use3 <- art_use3[art_use3$art %in% c("NNRTI","DTG"),]
art_use3$art<-factor(art_use3$art,levels=c("NNRTI", "DTG"))
  dtg_plot_a1<-ggplot(data=art_use3[art_use3$main_scen=="DTG as first-line",], mapping=aes(x=year,y=ART_use,
            color=factor(sub_scen,labels=c("All on NNRTI", "Only men on DTG",
              "All men and women not of childbearing age on DTG",
              "All men and women not at risk of pregnancy on DTG",
              "All men and women on DTG"))))+
  geom_line(size=1.1)+
  theme_bw()+
  facet_grid(.~art,switch="y")+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c(as.character(c(2010,2020,2030,2040))))+
  scale_y_continuous(labels = scales::percent)+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("Relative NNRTI/DTG use (in %)")+
  geom_vline(xintercept = 2019, linetype="dashed", 
             color = "black", size=1)+
  annotate("text",x=2010,y=0.8,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  scale_color_manual(values=rep(color,4))
  
  dtg_plot_a2<-ggplot(data=art_use3[art_use3$main_scen=="DTG as first-line and switch",], mapping=aes(x=year,y=ART_use,
                                   color=factor(sub_scen,labels=c("All on NNRTI", "Only men on DTG",
                                                                  "All men and women not of childbearing age on DTG",
                                                                  "All men and women not at risk of pregnancy on DTG",
                                                                  "All men and women on DTG"))))+
    geom_line(size=1.1)+
    theme_bw()+
    facet_grid(.~art,switch="y")+
    theme(legend.position="right")+
    scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c(as.character(c(2010,2020,2030,2040))))+
    scale_y_continuous(labels = scales::percent)+
    labs(color="Scenarios")+
    xlab("Calendar year")+
    ylab("Relative NNRTI/DTG use (in %)")+
    geom_vline(xintercept = 2019, linetype="dashed", 
               color = "black", size=1)+
    annotate("text",x=2010,y=0.8,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
    scale_color_manual(values=rep(color,4))
  
  legend<-get_legend(
    dtg_plot_a1 + theme(legend.box.margin = margin(0, 0, 0, 0))+
      theme(legend.position = "bottom",)+
      #legend.spacing.x = unit(1.0, 'cm')) +
      theme(legend.direction = "horizontal")+
      theme(legend.title=element_blank())+
      guides(colour=guide_legend(nrow=2,byrow=FALSE,
                                 keywidth=1.5,))+
      theme(legend.text=element_text(size=10))
  )
  dtg_plot<- plot_grid(dtg_plot_a1+theme(legend.position="none"),
                       dtg_plot_a2+theme(legend.position="none"),
                       align = 'vh',
                       labels = c("A", "B"),
                       hjust = -1,
                       nrow = 2)
  plot_grid(dtg_plot+
              theme(plot.margin = unit(c(0.5,0.5,0.5,0.5)/10, "cm")),
            legend, rel_heights = c(1, 0.08),nrow=2)#8.6x7 inches
##############################################################################################################################
##############################################################################################################################
#NNRTI resistance over time
#Database
nnrti_res=expand.grid(year=seq(2005,2040,1/12),
                      sub_scen=scen,
                      main_scen=c("DTG as first-line", "DTG as first-line and switch"))
nnrti_res$out=NA
nnrti_res$out_inf=NA
nnrti_res$out_sup=NA
#Function
outcome=function(data){
  tmax=dim(data)[1]
  rowSums(data[seq(1,tmax,1),select_15(c(2,12),1:4,2,1:2)])/rowSums(data[seq(1,tmax,1),select_15(c(2,12),1:4,1:2,1:2)])
}
#Computation: nnrti resistance level
index_list=c(1:5,1,6:9)
for(i in 1:10){
  #With simulations from ub2.job
  list_data=list()
  for(j in 1:201){
    load(paste("C:/Users/ahauser/Documents/Step2/R_results/dtg_sens/dtg_sens17/data1_",index_list[i],"_",j,".RData",sep=""))
    list_data[[j]]=data1
  }
  #With simulations from ub.job
  #load(paste("C:/Users/ahauser/Documents/Step2/R_results/dtg_sens/dtg_sens6/dtg_sens",index_list[i],".RData",sep=""))
  l1=1+(421)*(i-1)
  l2=(421)*i
  l=seq(l1,l2,1)
  
  nnrti_res[l,c("out")]=outcome(list_data[[1]])
  nnrti_res[l,c("out_inf","out_sup")]=t(apply(sapply(list_data,function(x){outcome(x)}),1,function(y) quantile(y,probs=c(0.025,0.975))))
  print(i)
}

nnrti_res[nnrti_res$year=="2040" & nnrti_res$sub_scen==scen[5] & nnrti_res$main_scen=="DTG as first-line and switch",]
View(nnrti_res[nnrti_res$year=="2030" | nnrti_res$year=="2040",])
View(nnrti_res[nnrti_res$year=="2040",])
View(nnrti_res[nnrti_res$year=="2035",])
##############################################################################################################################
#Check the prop of people with cd4<200 according to scenario. prop cd4<200 increase with DTG, probably because DTG reduces
#number of infected and therefore diagnosed people are infected for long time and, thus, have low cd4 counts
j=1
i=1
load(paste("C:/Users/ahauser/Documents/Step2/R_results/dtg_sens/dtg_sens12/data1_",i,"_",j,".RData",sep=""))
data_nodtg=data1

j=1
i=10
load(paste("C:/Users/ahauser/Documents/Step2/R_results/dtg_sens/dtg_sens12/data1_",i,"_",j,".RData",sep=""))
data_dtg=data1

outcome=function(data){
  tmax=dim(data)[1]
  rowSums(data[seq(1,tmax,1),select_15(c(2,12),1:4,2,1:2)])/rowSums(data[seq(1,tmax,1),select_15(c(2,12),1:4,1:2,1:2)])
}
outcome2=function(data){
  return(sum(data[1,select_15(13,1:4,2,1:2)])/sum(data[1,select_15(13,1:4,1:2,1:2)]))
}
outcome3=function(data){
  tmax=dim(data)[1]
  return(rowSums(data[1:tmax,select_15(13,1,1:2,1:2)])/rowSums(data[1:tmax,select_15(13,1:4,1:2,1:2)]))
}

outcome(data_nodtg)
outcome(data_dtg)

outcome2(data_nodtg)
outcome2(data_dtg)

outcome3(data_nodtg)
outcome3(data_dtg)

##############################################################################################################################
#Figure

color=c("black","blue","violet","orange","red")

dtg_plot2<-ggplot(data=nnrti_res, mapping=aes(x=year,y=out,
                                              color=factor(sub_scen,labels=c("All on NNRTI", "Only men on DTG",
                                                                             "All men and women not of\n childbearing age on DTG",
                                                                             "All men and women not at\n risk of pregnancy on DTG",
                                                                             "All men and women on DTG"))))+
geom_line(size=1.1)+
  theme_bw()+
  facet_grid(main_scen~.,switch="y")+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2044))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.7))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("Percentage of NNRTI resistance among diagnosed individuals")+
  geom_vline(xintercept = 2019, linetype="dashed", 
             color = "black", size=1)+
  annotate("text",x=2020,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  geom_rect(data=nnrti_res,xmin=2040.5,xmax=2044.5,ymin=0,ymax=0.7,colour="grey50",fill="white",size=.3)+
  geom_segment(data=nnrti_res,x=2040,xend=2040.5,y=-.03,yend=0,size=.3,colour="grey50")+
  geom_pointrange(data=nnrti_res[nnrti_res$year==2040,],aes(x=rep(seq(2041,2044,length.out = 5),2),y=out,ymin=out_inf,ymax=out_sup),alpha=1,fill="white",size=.8,stroke=.4,shape=1)+
  scale_color_manual(values=rep(color,2))


nnrti_res1<-nnrti_res[nnrti_res$main_scen=="DTG as first-line",]
dtg_plot_s1<-ggplot(data=nnrti_res1, mapping=aes(x=year,y=out,
                                              color=factor(sub_scen,labels=c("All on NNRTI          ", "Only men on DTG  ",
                                                                             "All men and women not of childbearing age on DTG  ",
                                                                             "All men and women not at risk of pregnancy on DTG  ",
                                                                             "All men and women on DTG  "))))+
  geom_line(size=1.1)+
  theme_bw()+
  #facet_grid(main_scen~.,switch="y")+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2044))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.7))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("NNRTI TDR (in %)")+
  geom_vline(xintercept = 2019, linetype="dashed", 
             color = "black", size=1)+
  annotate("text",x=2020,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  geom_rect(data=nnrti_res1,xmin=2040.5,xmax=2044.5,ymin=0,ymax=0.7,colour="grey50",fill="white",size=.3)+
  geom_segment(data=nnrti_res1,x=2040,xend=2040.5,y=-.03,yend=0,size=.3,colour="grey50")+
  geom_pointrange(data=nnrti_res1[nnrti_res1$year==2040,],aes(x=rep(seq(2041,2044,length.out = 5),1),y=out,ymin=out_inf,ymax=out_sup),alpha=1,fill="white",size=.8,stroke=.4,shape=1)+
  scale_color_manual(values=rep(color,2))


nnrti_res2<-nnrti_res[nnrti_res$main_scen=="DTG as first-line and switch",]
dtg_plot_s2<-ggplot(data=nnrti_res2, mapping=aes(x=year,y=out,
                                                       color=factor(sub_scen,labels=c("All on NNRTI          ", "Only men on DTG  ",
                                                                                                                          "All men and women not of\n childbearing age on DTG  ",
                                                                                                                          "All men and women not at\n risk of pregnancy on DTG  ",
                                                                                                                          "All men and women on DTG  "))))+
  geom_line(size=1.1)+
  theme_bw()+
  #facet_grid(main_scen~.,switch="y")+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2044))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.7))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("NNRTI TDR (in %)")+
  geom_vline(xintercept = 2019, linetype="dashed", 
             color = "black", size=1)+
  annotate("text",x=2020,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  geom_rect(data=nnrti_res2,xmin=2040.5,xmax=2044.5,ymin=0,ymax=0.7,colour="grey50",fill="white",size=.3)+
  geom_segment(data=nnrti_res2,x=2040,xend=2040.5,y=-.03,yend=0,size=.3,colour="grey50")+
  geom_pointrange(data=nnrti_res2[nnrti_res2$year==2040,],aes(x=rep(seq(2041,2044,length.out = 5),1),y=out,ymin=out_inf,ymax=out_sup),alpha=1,fill="white",size=.8,stroke=.4,shape=1)+
  scale_color_manual(values=rep(color,2))

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
dtg_plot<- plot_grid(dtg_plot_s1+theme(legend.position="none")+
                       theme(plot.margin=unit(c(1,0,0,0),"cm")),
                     dtg_plot_s2+theme(legend.position="none")+
                       theme(plot.margin=unit(c(1,0,0,0),"cm")),
                     align = 'vh',
                     labels = c("A. DTG as first-line",
                                "B. DTG as first-line and switch"),
                     hjust = 0,
                     label_x = 0.05,
                     nrow = 2)
plot_grid(dtg_plot+
          theme(plot.margin = unit(c(0.5,0.5,0.5,0.5)/10, "cm")),
          legend, rel_heights = c(1, 0.1),nrow=2)
##############################################################################################################################
##############################################################################################################################
#Sensitivity analysis
##############################################################################################################################
##############################################################################################################################
#Sensitivity: interruption rates
#Assumptions
#load params from value_gender4_v3.R and p1, p2 from graphs.R
#source("C:/Users/ahauser/Documents/Step1/Rfiles/used scripts/marisa/value_gender4_v3.R")
source("C:/Users/ahauser/Documents/Step2/dtg/value_gender4_v5.R")
source("C:/Users/ahauser/Documents/Step2/dtg/parmin_dtg.R")
p1=c(rate1_inf=3.318999e+00, rate2_inf=5.440886e-01, rate1_diag=2.736363e+02, rate2_diag=7.710578e+00, rate_treat=1.879695e-03,rate_death=1.632000e-01)
p2=c(rate1_inf=NA, rate2_inf=NA, rate1_diag=NA, rate2_diag=NA, rate_treat=NA,rate_death=NA,q=0.05,rate_ratio=0.5,k1=1,k2=2,k3=2,alpha=2,rate_res=5,rate_susc=125)
#no starting with second-line PI
params<-within(params,{
  #RateDirectTreatSecond=4.35/c(5000,12000,13000,1850)
  RateTreatSecond=4.35/c(531,427,294,189)*1/5
  RateDirectTreatSecond=c(0,0,0,0)
})
#Resistance parameters
p2["alpha"]=2.64
p2["alpha2"]=4.9


##################################################################################### 
#Simulations
source("C:/Users/ahauser/Documents/Step2/dtg/parmin_dtg.R")
sourceCpp("C:/Users/ahauser/Documents/Step2/dtg/solvdiff_cpp_dtg_2_2005.cpp")
xstart_dtg_2_2005=xstart_f_dtg_2(p2,treat_dtgb)
theta=p1
#15 compartments
data=my_fun10_solver2_dtg_2(xstart_dtg_2_2005,theta,treat_dtgb,params,p2)
figures_2_dtg_update(data)
save(data,file="C:/Users/ahauser/Documents/Step2/R_results/data_dtg_2005_2019.RData")
#load function, run from 2005 to 2018
#Counterfactual scenarios from 2018 to 2038
sourceCpp("C:/Users/ahauser/Documents/Step2/dtg/solvdiff_cpp_dtg_2_2019.cpp")
xstart_dtg_2_2018=data[dim(data)[1],]
dtg_list=list()
for(i in 1:10){
  xstart=xstart_dtg_2_2018
  
  xstart[select_15(12:15,1:4,1:2,2)]=(1-treat_dtg_mat[i,1])*(xstart_dtg_2_2018[select_15(2:5,1:4,1:2,2)]+xstart_dtg_2_2018[select_15(12:15,1:4,1:2,2)])
  xstart[select_15(2:5,1:4,1:2,2)]=treat_dtg_mat[i,1]*(xstart_dtg_2_2018[select_15(2:5,1:4,1:2,2)]+xstart_dtg_2_2018[select_15(12:15,1:4,1:2,2)])
  xstart=c(xstart[1:240],xstart_dtg_2_2018[241:250])
  
  data1=my_fun10_solver2_dtg_2(xstart,theta,treat_dtg_mat[i,],params,p2)
  rownames(data1)=1:(dim(data1)[1])
  data1=rbind(data[-dim(data)[1],],data1)
  dtg_list[[i]]<-data1
  print(i)
}
figures_dtg_4_update(dtg_list)
dtg_noint=dtg_list

##################################################################################### 
#NNRTI resistance over time
#Database
nnrti_res=expand.grid(year=seq(2005,2040,1/12),
                      sub_scen=scen,
                      main_scen=c("DTG as first-line", "DTG as first-line and switch"),
                      treat_int=c("No treatment interruption","Treatment interruption"))
nnrti_res$out=NA
#Function
outcome=function(data){
  tmax=dim(data)[1]
  rowSums(data[seq(1,tmax,1),select_15(c(2,12),1:4,2,1:2)])/rowSums(data[seq(1,tmax,1),select_15(c(2,12),1:4,1:2,1:2)])
}
#Computation: nnrti resistance level
index_list=c(1:5,1,6:9)
for(i in 1:10){
  l1=1+(421)*(i-1)
  l2=(421)*i
  l=seq(l1,l2,1)
  
  nnrti_res[l,c("out")]=outcome(dtg_noint[[index_list[i]]])
}
for(i in 1:10){
  l1=421*10+1+(421)*(i-1)
  l2=421*10+(421)*i
  l=seq(l1,l2,1)
  
  nnrti_res[l,c("out")]=outcome(dtg_ref[[index_list[i]]])
}


#####################################################################################  
dtg_plot2<-ggplot(data=nnrti_res, mapping=aes(x=year,y=out,
                                              color=factor(sub_scen,labels=c("All on NNRTI", "Only men on DTG",
                                                                             "All men and women not of childbearing age on DTG",
                                                                             "All men and women not at risk of pregnancy on DTG",
                                                                             "All men and women on DTG"))))+
  geom_line(aes(size=2))+
  theme_bw()+
  facet_grid(main_scen~treat_int,switch="y")+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2041))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.65))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("Percentage of NNRTI resistance among diagnosed individuals")+
  geom_vline(xintercept = 2019, linetype="dashed", 
             color = "black", size=1)+
  annotate("text",x=2020,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  scale_color_manual(values=rep(color,2))+ 
  scale_size(range=c(0.8, 1), guide=FALSE)+
  scale_alpha(range=c(0.3, 1), guide=FALSE)+
  theme(legend.box.margin = margin(0, 0, 0, 0))+
  theme(legend.position = "bottom",)+
  #legend.spacing.x = unit(1.0, 'cm')) +
  theme(legend.direction = "horizontal")+
  theme(legend.title=element_blank())+
  guides(colour=guide_legend(nrow=2,byrow=FALSE,
                             keywidth=1.5,))+
  theme(legend.text=element_text(size=10))
dtg_plot2

#####################################################################################  
dtg_plot2<-ggplot(data=nnrti_res, mapping=aes(x=year,y=out,group=interaction(sub_scen,treat_int),alpha=treat_int,size=treat_int,
                                              color=factor(sub_scen,labels=c("All on NNRTI", "Only men on DTG",
                                                                             "All men and women not of\n childbearing age on DTG",
                                                                             "All men and women not at\n risk of pregnancy on DTG",
                                                                             "All men and women on DTG"))))+
  geom_line()+
  theme_bw()+
  facet_grid(main_scen~.,switch="y")+
  theme(legend.position="right")+
  scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2041))+
  scale_y_continuous(labels = scales::percent,limits=c(0,0.65))+
  labs(color="Scenarios")+
  xlab("Calendar year")+
  ylab("Percentage of NNRTI resistance among diagnosed individuals")+
  geom_vline(xintercept = 2019, linetype="dashed", 
             color = "black", size=1)+
  annotate("text",x=2020,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
  scale_color_manual(values=rep(color,2))+ 
  scale_size(range=c(0.8, 1), guide=FALSE)+
  scale_alpha(range=c(0.3, 1), guide=FALSE)
dtg_plot2

##############################################################################################################################
##############################################################################################################################
#Sensitivity: no Treat-All policy
#Assumptions
#load params from value_gender4_v3.R and p1, p2 from graphs.R
#source("C:/Users/ahauser/Documents/Step1/Rfiles/used scripts/marisa/value_gender4_v3.R")
source("C:/Users/ahauser/Documents/Step2/dtg/value_gender4_v5.R")
source("C:/Users/ahauser/Documents/Step2/dtg/parmin_dtg.R")
p1=c(rate1_inf=3.318999e+00, rate2_inf=5.440886e-01, rate1_diag=2.736363e+02, rate2_diag=7.710578e+00, rate_treat=1.879695e-03,rate_death=1.632000e-01)
p2=c(rate1_inf=NA, rate2_inf=NA, rate1_diag=NA, rate2_diag=NA, rate_treat=NA,rate_death=NA,q=0.05,rate_ratio=0.5,k1=1,k2=2,k3=2,alpha=2,rate_res=5,rate_susc=125)
#no treatment interruption
params<-within(params,{
  RateStopTreatFirst=4.35/c(Inf,Inf,Inf,Inf)
  RateStopSuppFirst=4.35/c(Inf,Inf,Inf,Inf)
  RateStopFailFirst=4.35/c(Inf,Inf,Inf,Inf)
  RateStopTreatSecond=4.35/c(Inf,Inf,Inf,Inf) #Assumption : T2 same as T1
  RateStopSuppSecond=4.35/c(Inf,Inf,Inf,Inf)
  RateStopFailSecond=4.35/c(Inf,Inf,Inf,Inf)
})
#no starting with second-line PI
params<-within(params,{
  #RateDirectTreatSecond=4.35/c(5000,12000,13000,1850)
  RateTreatSecond=4.35/c(531,427,294,189)*1/5
  RateDirectTreatSecond=c(0,0,0,0)
})
#Resistance parameters
p2["alpha"]=2.64
p2["alpha2"]=4.9
#no treat all policy
params<-within(params,{
  TreatAll_frame$d=rep(1,4)
  })
  ##################################################################################### 
  #Simulations
  source("C:/Users/ahauser/Documents/Step2/dtg/parmin_dtg.R")
  sourceCpp("C:/Users/ahauser/Documents/Step2/dtg/solvdiff_cpp_dtg_2_2005.cpp")
  xstart_dtg_2_2005=xstart_f_dtg_2(p2,treat_dtgb)
  theta=p1
  #15 compartments
  data=my_fun10_solver2_dtg_2(xstart_dtg_2_2005,theta,treat_dtgb,params,p2)
  figures_2_dtg_update(data)
  #load function, run from 2005 to 2018
  #Counterfactual scenarios from 2018 to 2038
  sourceCpp("C:/Users/ahauser/Documents/Step2/dtg/solvdiff_cpp_dtg_2_2019.cpp")
  xstart_dtg_2_2018=data[dim(data)[1],]
  dtg_list=list()
  for(i in 1:10){
    xstart=xstart_dtg_2_2018
    
    xstart[select_15(12:15,1:4,1:2,2)]=(1-treat_dtg_mat[i,1])*(xstart_dtg_2_2018[select_15(2:5,1:4,1:2,2)]+xstart_dtg_2_2018[select_15(12:15,1:4,1:2,2)])
    xstart[select_15(2:5,1:4,1:2,2)]=treat_dtg_mat[i,1]*(xstart_dtg_2_2018[select_15(2:5,1:4,1:2,2)]+xstart_dtg_2_2018[select_15(12:15,1:4,1:2,2)])
    xstart=c(xstart[1:240],xstart_dtg_2_2018[241:250])
    
    data1=my_fun10_solver2_dtg_2(xstart,theta,treat_dtg_mat[i,],params,p2)
    rownames(data1)=1:(dim(data1)[1])
    data1=rbind(data[-dim(data)[1],],data1)
    dtg_list[[i]]<-data1
    print(i)
  }
  figures_dtg_4_update(dtg_list)
  dtg_treatall=dtg_list
  
  ##################################################################################### 
  #NNRTI resistance over time
  #Database
  nnrti_res=expand.grid(year=seq(2005,2040,1/12),
                        sub_scen=scen,
                        main_scen=c("DTG as first-line", "DTG as first-line and switch"),
                        treat_all=c("No Treat-All", "Treat-All"))
  nnrti_res$out=NA
  nnrti_res$out_inf=NA
  nnrti_res$out_sup=NA
  #Function
  outcome=function(data){
    tmax=dim(data)[1]
    rowSums(data[seq(1,tmax,1),select_15(c(2,12),1:4,2,1:2)])/rowSums(data[seq(1,tmax,1),select_15(c(2,12),1:4,1:2,1:2)])
  }
  #Computation: nnrti resistance level
  index_list=c(1:5,1,6:9)
  for(i in 1:10){
    l1=1+(421)*(i-1)
    l2=(421)*i
    l=seq(l1,l2,1)
    
    nnrti_res[l,c("out")]=outcome(dtg_treatall[[index_list[i]]])
  }
  for(i in 1:10){
    l1=421*10+1+(421)*(i-1)
    l2=421*10+(421)*i
    l=seq(l1,l2,1)
    
    nnrti_res[l,c("out")]=outcome(dtg_ref[[index_list[i]]])
  }
  
  #####################################################################################  
  dtg_plot2<-ggplot(data=nnrti_res, mapping=aes(x=year,y=out,
                                                color=factor(sub_scen,labels=c("All on NNRTI", "Only men on DTG",
                                                                               "All men and women not of childbearing age on DTG",
                                                                               "All men and women not at risk of pregnancy on DTG",
                                                                               "All men and women on DTG"))))+
    geom_line(aes(size=2))+
    theme_bw()+
    facet_grid(main_scen~treat_all,switch="y")+
    theme(legend.position="right")+
    scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2041))+
    scale_y_continuous(labels = scales::percent,limits=c(0,0.65))+
    labs(color="Scenarios")+
    xlab("Calendar year")+
    ylab("Percentage of NNRTI resistance among diagnosed individuals")+
    geom_vline(xintercept = 2019, linetype="dashed", 
               color = "black", size=1)+
    annotate("text",x=2020,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
    scale_color_manual(values=rep(color,2))+ 
    scale_size(range=c(0.8, 1), guide=FALSE)+
    scale_alpha(range=c(0.3, 1), guide=FALSE)+
    theme(legend.box.margin = margin(0, 0, 0, 0))+
    theme(legend.position = "bottom",)+
    #legend.spacing.x = unit(1.0, 'cm')) +
    theme(legend.direction = "horizontal")+
    theme(legend.title=element_blank())+
    guides(colour=guide_legend(nrow=2,byrow=FALSE,
                               keywidth=1.5,))+
    theme(legend.text=element_text(size=10))
  dtg_plot2 
  
  #####################################################################################  
  dtg_plot2<-ggplot(data=nnrti_res, mapping=aes(x=year,y=out,group=interaction(sub_scen,treat_all),alpha=treat_all,size=treat_all,
                                                color=factor(sub_scen,labels=c("All on NNRTI", "Only men on DTG",
                                                                               "All men and women not of\n childbearing age on DTG",
                                                                               "All men and women not at\n risk of pregnancy on DTG",
                                                                               "All men and women on DTG"))))+
    geom_line(aes(size=2))+
    theme_bw()+
    facet_grid(main_scen~.,switch="y")+
    theme(legend.position="right")+
    scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2041))+
    scale_y_continuous(labels = scales::percent,limits=c(0,0.65))+
    labs(color="Scenarios")+
    xlab("Calendar year")+
    ylab("Percentage of NNRTI resistance among diagnosed individuals")+
    geom_vline(xintercept = 2019, linetype="dashed", 
               color = "black", size=1)+
    annotate("text",x=2020,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
    scale_color_manual(values=rep(color,2))+ 
    scale_size(range=c(0.8, 1), guide=FALSE)+
    scale_alpha(range=c(0.3, 1), guide=FALSE)+
    theme(legend.box.margin = margin(0, 0, 0, 0))+
    theme(legend.position = "bottom",)+
    #legend.spacing.x = unit(1.0, 'cm')) +
    theme(legend.direction = "horizontal")+
    theme(legend.title=element_blank())+
    guides(colour=guide_legend(nrow=2,byrow=FALSE,
                               keywidth=1.5,))+
    theme(legend.text=element_text(size=10))
  dtg_plot2
##############################################################################################################################
##############################################################################################################################
  #Sensitivity: DTG efficacy
  #Assumptions
  #load params from value_gender4_v3.R and p1, p2 from graphs.R
  #source("C:/Users/ahauser/Documents/Step1/Rfiles/used scripts/marisa/value_gender4_v3.R")
  source("C:/Users/ahauser/Documents/Step2/dtg/value_gender4_v4.R")
  source("C:/Users/ahauser/Documents/Step2/dtg/parmin_dtg.R")
  p1=c(rate1_inf=3.318999e+00, rate2_inf=5.440886e-01, rate1_diag=2.736363e+02, rate2_diag=7.710578e+00, rate_treat=1.879695e-03,rate_death=1.632000e-01)
  p2=c(rate1_inf=NA, rate2_inf=NA, rate1_diag=NA, rate2_diag=NA, rate_treat=NA,rate_death=NA,q=0.05,rate_ratio=0.5,k1=1,k2=2,k3=2,alpha=2,rate_res=5,rate_susc=125)
  #no treatment interruption
  params<-within(params,{
    RateStopTreatFirst=4.35/c(Inf,Inf,Inf,Inf)
    RateStopSuppFirst=4.35/c(Inf,Inf,Inf,Inf)
    RateStopFailFirst=4.35/c(Inf,Inf,Inf,Inf)
    RateStopTreatSecond=4.35/c(Inf,Inf,Inf,Inf) #Assumption : T2 same as T1
    RateStopSuppSecond=4.35/c(Inf,Inf,Inf,Inf)
    RateStopFailSecond=4.35/c(Inf,Inf,Inf,Inf)
  })
  #no starting with second-line PI
  params<-within(params,{
    #RateDirectTreatSecond=4.35/c(5000,12000,13000,1850)
    RateTreatSecond=4.35/c(531,427,294,189)*1/5
    RateDirectTreatSecond=c(0,0,0,0)
  })
  #Resistance parameters
  # p2["alpha"]=2
  # p2["alpha2"]=2
  
  p2["alpha"]=2.64
  p2["alpha2"]=4.9
  
  #####################################################################################  
  #Simulations
  source("C:/Users/ahauser/Documents/Step2/dtg/parmin_dtg.R")
  sourceCpp("C:/Users/ahauser/Documents/Step2/dtg/solvdiff_cpp_dtg_2_2005.cpp")
  xstart_dtg_2_2005=xstart_f_dtg_2(p2,treat_dtgb)
  theta=p1
  #15 compartments
  data=my_fun10_solver2_dtg_2(xstart_dtg_2_2005,theta,treat_dtgb,params,p2)
  figures_2_dtg_update(data)
  save(data,file="C:/Users/ahauser/Documents/Step2/R_results/data_dtg_2005_2019.RData")
  #load function, run from 2005 to 2018
  #Counterfactual scenarios from 2018 to 2038
  sourceCpp("C:/Users/ahauser/Documents/Step2/dtg/solvdiff_cpp_dtg_2_2019.cpp")
  xstart_dtg_2_2018=data[dim(data)[1],]
  dtg_list=list()
  for(i in 1:10){
    xstart=xstart_dtg_2_2018
    
    xstart[select_15(12:15,1:4,1:2,2)]=(1-treat_dtg_mat[i,1])*(xstart_dtg_2_2018[select_15(2:5,1:4,1:2,2)]+xstart_dtg_2_2018[select_15(12:15,1:4,1:2,2)])
    xstart[select_15(2:5,1:4,1:2,2)]=treat_dtg_mat[i,1]*(xstart_dtg_2_2018[select_15(2:5,1:4,1:2,2)]+xstart_dtg_2_2018[select_15(12:15,1:4,1:2,2)])
    xstart=c(xstart[1:240],xstart_dtg_2_2018[241:250])
    
    list_eff=list()
    dtg_eff<-seq(1,2,length.out = 5)
    for(j in 1:5){
      params<-within(params,{
        RateTreatToFailFirstDTG=4.35/c(203,198,164,112)*1/dtg_eff[j]
        RateFailFirstDTG=4.35/c(767,582,270,96)*1/dtg_eff[j]
        RateSuppFirstDTG=4.35/c(30,30,31,34)*dtg_eff[j]
        RateFailToSuppTreatFirstDTG=4.35/c(10,56,24,51)*dtg_eff[j]
      })
      data1=my_fun10_solver2_dtg_2(xstart,theta,treat_dtg_mat[i,],params,p2)
      rownames(data1)=1:(dim(data1)[1])
      data1=rbind(data[-dim(data)[1],],data1)
      list_eff[[j]]<-data1
      
    }
    
    dtg_list[[i]]<-list_eff
    print(i)
  }
  dtg_eff=dtg_list
  
  #####################################################################################  
  #NNRTI resistance over time
  #Database
  nnrti_res=expand.grid(year=seq(2005,2040,1/12),
                        eff=seq(1,2,length.out=5),
                        sub_scen=scen,
                        main_scen=c("DTG as first-line", "DTG as first-line and switch"))
  nnrti_res$out=NA
  
  #Function
  outcome=function(data){
    tmax=dim(data)[1]
    rowSums(data[seq(1,tmax,1),select_15(c(2,12),1:4,2,1:2)])/rowSums(data[seq(1,tmax,1),select_15(c(2,12),1:4,1:2,1:2)])
  }
  #Computation: nnrti resistance level
  index_list=c(1:5,1,6:9)
  for(i in 1:10){
    list<-dtg_eff[[index_list[i]]]
    for(j in 1:5){
      n=(i-1)*5+j
      l1=1+(421)*(n-1)
      l2=(421)*n
      l=seq(l1,l2,1)
      nnrti_res[l,c("out")]=outcome(list[[j]])
      print(l1)
      print(l2)
    }
  }
  nnrti_res$gg_alpha<-1
  nnrti_res$gg_alpha[nnrti_res$eff>1]<-0
  #####################################################################################  
  dtg_plot2<-ggplot(data=nnrti_res, mapping=aes(x=year,y=out,group=interaction(sub_scen,eff),alpha=gg_alpha,size=gg_alpha,
                                                color=factor(sub_scen,labels=c("All on NNRTI", "Only men on DTG",
                                                                               "All men and women not of childbearing age on DTG",
                                                                               "All men and women not at risk of pregnancy on DTG",
                                                                               "All men and women on DTG"))))+
    geom_line(aes(size=2))+
    theme_bw()+
    facet_grid(main_scen~.,switch="y")+
    theme(legend.position="right")+
    scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2041))+
    scale_y_continuous(labels = scales::percent,limits=c(0,0.65))+
    labs(color="Scenarios")+
    xlab("Calendar year")+
    ylab("Percentage of NNRTI resistance among diagnosed individuals")+
    geom_vline(xintercept = 2019, linetype="dashed", 
               color = "black", size=1)+
    annotate("text",x=2020,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
    scale_color_manual(values=rep(color,2))+ 
    scale_size(range=c(0.8, 1), guide=FALSE)+
    scale_alpha(range=c(0.3, 1), guide=FALSE)+
    theme(legend.box.margin = margin(0, 0, 0, 0))+
    theme(legend.position = "bottom",)+
    #legend.spacing.x = unit(1.0, 'cm')) +
    theme(legend.direction = "horizontal")+
    theme(legend.title=element_blank())+
    guides(colour=guide_legend(nrow=2,byrow=FALSE,
                               keywidth=1.5,))+
    theme(legend.text=element_text(size=10))
  dtg_plot2
  #####################################################################################  
  #not important
  nnrti_res=expand.grid(year=seq(2005,2040,1/12),
                        eff=seq(1,2,length.out=11),
                        sub_scen=scen,
                        main_scen=c("DTG as first-line", "DTG as first-line and switch"))
  nnrti_res$out=NA
  
  #Function
  outcome=function(data){
    tmax=dim(data)[1]
    rowSums(data[seq(1,tmax,1),select_15(c(2,12),1:4,2,1:2)])
  }
  #Computation: nnrti resistance level
  index_list=c(1:5,1,6:9)
  for(i in 1:10){
    list<-dtg_eff[[index_list[i]]]
    for(j in 1:11){
      n=(i-1)*11+j
      l1=1+(421)*(n-1)
      l2=(421)*n
      l=seq(l1,l2,1)
      nnrti_res[l,c("out")]=outcome(list[[j]])
      print(l1)
      print(l2)
    }
  }
  nnrti_res$gg_alpha<-1
  nnrti_res$gg_alpha[nnrti_res$eff>1]<-0
  #####################################################################################  
  dtg_plot2<-ggplot(data=nnrti_res, mapping=aes(x=year,y=out,group=interaction(sub_scen,eff),alpha=gg_alpha,size=gg_alpha,
                                                color=factor(sub_scen,labels=c("All on NNRTI", "Only men on DTG",
                                                                               "All men and women not of\n childbearing age on DTG",
                                                                               "All men and women not at\n risk of pregnancy on DTG",
                                                                               "All men and women on DTG"))))+
    geom_line()+
    theme_bw()+
    facet_grid(main_scen~.,switch="y")+
    theme(legend.position="right")+
    scale_x_continuous(breaks=c(2010,2020,2030,2040),labels=c("2010","2020","2030","2040"),limits=c(2005,2041))+
    scale_y_continuous(limits=c(0,500))+
    labs(color="Scenarios")+
    xlab("Calendar year")+
    ylab("NNRTI resistance cases among diagnosed")+
    geom_vline(xintercept = 2019, linetype="dashed", 
               color = "black", size=1)+
    annotate("text",x=2020,y=0.5,hjust=0,label="DTG\nintroduction",color="black",size=3.5)+
    scale_color_manual(values=rep(color,2))+ 
    scale_size(range=c(0.8, 1), guide=FALSE)+
    scale_alpha(range=c(0.3, 1), guide=FALSE)
  dtg_plot2
##############################################################################################################################
##############################################################################################################################
  #Stage
  #Database
  scen=c("1) all men on DTG","2) all men and women not of childbearing age on DTG",
         "3) all men and women not at risk of pregancy on DTG","4) all men and women on DTG")
  res_stage=expand.grid(stage=c("Unsuppressed, on ART","Infected"),
                        gender=c("men","women"),
                        sub_scen=scen,
                        main_scen=c("DTG as first-line", "DTG as first-line and switch"),
                        res=NA,
                        res_inf=NA,
                        res_sup=NA)
  stage_ind=list(c(5,8,11,15),c(1))
  gender_ind=list(c(1),c(2))
  
  #Function
  outcome=function(data,stage,gender){
    stage_ind=as.numeric((list(c(5,8,11,15),c(1)))[[which(c("Unsuppressed, on ART","Infected")==as.character(stage))]])
    gender_ind=as.vector((list(c(1),c(2)))[[which(c("men","women")==as.character(gender))]])
    return(sum(data[20*12+1,select_15(stage_ind,1:4,2,gender_ind)])/sum(data[20*12+1,select_15(stage_ind,1:4,1:2,gender_ind)]))
  }
  #Computation: nnrti resistance level
  index_list=c(2:9)
  for(i in 1:8){
    #With simulations from ub2.job
    list_data=list()
    for(j in 1:201){
      load(paste("C:/Users/ahauser/Documents/Step2/R_results/dtg_sens/dtg_sens9/data1_",index_list[i],"_",j,".RData",sep=""))
      list_data[[j]]=data1
    }
    #With simulations from ub.job
    #load(paste("C:/Users/ahauser/Documents/Step2/R_results/dtg_sens/dtg_sens6/dtg_sens",i+1,".RData",sep=""))
    res_stage$res[(4*(i-1)+1):(4*i)]=apply(res_stage[(4*(i-1)+1):(4*i),],1,function(x) outcome(list_data[[1]],as.character(x["stage"]),as.character(x["gender"])))
    res_stage[(4*(i-1)+1):(4*i),c("res_inf","res_sup")]=t(apply(res_stage[(4*(i-1)+1):(4*i),],1,function(x){
                quantile(sapply(list_data,function(y) outcome(y,x["stage"],x["gender"])),probs=c(0.025,0.975))
              }))
    print(i)
  }
  
  res_stage=res_stage[res_stage$sub_scen=="2) all men and women not of childbearing age on DTG" |
                        res_stage$sub_scen== "4) all men and women on DTG",]
  res_stage=res_stage[res_stage$main_scen=="DTG as first-line and switch",]
  
##############################################################################################################################
#Figure
  pd = position_dodge(.4)
  ggplot(res_stage,
    aes(x =stage, 
        y =res,
        color = gender)) +
    geom_point(shape = 15, 
            size = 4, 
         position = pd) +
    facet_grid(. ~ sub_scen)+
  
  geom_errorbar(aes(ymin = res_inf,
                    ymax = res_sup),
                width = 0.2, 
                size = 0.7, 
                position = pd) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold")) +
  ylab("Level of NNRTI resistance")+
    theme(axis.text.x = element_text(angle = 0))+
    scale_y_continuous(labels = scales::percent,limits=c(0,1))

##############################################################################################################################
##############################################################################################################################
  #Suppression
  #Database
  scen=c("All on NNRTI","All men on DTG"," All men and women not of childbearing age on DTG",
         "All men and women not at risk of pregancy on DTG","All men and 99% of women on DTG")
  sup_level=expand.grid(year_treated=c(1,2,5),
                        year_start=c(2025,2030,2035),
                        sub_scen=scen,
                        main_scen=c("DTG as first-line", "DTG as first-line and switch"),
                        sup=NA,
                        sup_inf=NA,
                        sup_sup=NA,
                        nnrti_res=NA,
                        nnrti_res_inf=NA,
                        nnrti_res_sup=NA,
                        prop200=NA)
  #Function
  outcome=function(data,years,mort=FALSE){
    d1=0
    d2=0
    if(mort){
      d1=data[years,241]
      d2=data[years,242]
    }
    if(length(years)==1){
      return((sum(data[years,select_15(15,1:4,1:2,1:2)])+d2)/
               (sum(data[years,select_15(13:15,1:4,1:2,1:2)])+d1+d2))
    }else{
      return((rowSums(data[years,select_15(15,1:4,1:2,1:2)])+d2)/
               (rowSums(data[years,select_15(13:15,1:4,1:2,1:2)])+d1+d2))
    }
  }
  outcome2=function(data){
    return(sum(data[1,select_15(13,1:4,2,1:2)])/sum(data[1,select_15(13,1:4,1:2,1:2)]))
  }
  outcome3=function(data){
    return(sum(data[1,select_15(13,4,1:2,1:2)])/sum(data[1,select_15(13,1:4,1:2,1:2)]))
  }
  #Computation: suppression level
  #index_list=c(1:5,1,6:9)
  index_list=c(1:5,1,6:8,10) #choose 10 instead of 9 as 10 is 99% and therefore allow to calculate suppression
  #5 is taken but doesn't work as it has 100% of DTG-eligible people, so there are no people in 13
  for(i in 1:10){
    
    #Load job 2 arrays
    data_sup_list2025=as.list(rep(NA,201))
    data_sup_list2030=as.list(rep(NA,201))
    data_sup_list2035=as.list(rep(NA,201))
    for(j in 1:201){
      load(paste("C:/Users/ahauser/Documents/Step2/R_results/dtg_supp/dtg_supp8_2/data_sup_list2025","_",index_list[i],"_",j,".RData",sep=""))
      data_2025=matrix(data_2025,nrow=169)
      data_sup_list2025[[j]]=data_2025
      
      load(paste("C:/Users/ahauser/Documents/Step2/R_results/dtg_supp/dtg_supp8_2/data_sup_list2030","_",index_list[i],"_",j,".RData",sep=""))
      data_2030=matrix(data_2030,nrow=169)
      data_sup_list2030[[j]]=data_2030
      
      load(paste("C:/Users/ahauser/Documents/Step2/R_results/dtg_supp/dtg_supp8_2/data_sup_list2035","_",index_list[i],"_",j,".RData",sep=""))
      data_2035=matrix(data_2035,nrow=169)
      data_sup_list2035[[j]]=data_2035
    }
    #Load job 1 array
    #load(paste("C:/Users/ahauser/Documents/Step2/R_results/dtg_supp/dtg_supp1/data_sup_list2025","_",index_list[i],".RData",sep=""))
    
    
    supp2025=outcome(data_sup_list2025[[1]],c(13,25,61),mort=TRUE)
    supp2025_sens=t(sapply(as.list(c(13,25,61)),function(x){
      quantile(sapply(data_sup_list2025,function(y) outcome(y,x,mort=TRUE)),probs=c(0.025,0.975),na.rm=TRUE)
      }))
    
    supp2030=outcome(data_sup_list2030[[1]],c(13,25,61),mort=TRUE)
    supp2030_sens=t(sapply(as.list(c(13,25,61)),function(x){
      quantile(sapply(data_sup_list2030,function(y) outcome(y,x,mort=TRUE)),probs=c(0.025,0.975),na.rm=TRUE)
    }))
    
    supp2035=outcome(data_sup_list2035[[1]],c(13,25,61),mort=TRUE)
    supp2035_sens=t(sapply(as.list(c(13,25,61)),function(x){
      quantile(sapply(data_sup_list2035,function(y) outcome(y,x,mort=TRUE)),probs=c(0.025,0.975),na.rm=TRUE)
    }))
    
    res2025_sens=quantile(sapply(data_sup_list2025,function(y) outcome2(y)),probs=c(0.025,0.975),na.rm=TRUE)
    res2030_sens=quantile(sapply(data_sup_list2030,function(y) outcome2(y)),probs=c(0.025,0.975),na.rm=TRUE)
    res2035_sens=quantile(sapply(data_sup_list2035,function(y) outcome2(y)),probs=c(0.025,0.975),na.rm=TRUE)
    
    sup_level[((i-1)*9+1):(i*9),"sup"]=c(supp2025,supp2030,supp2035)
    sup_level[((i-1)*9+1):(i*9),c("sup_inf","sup_sup")]=rbind(supp2025_sens,supp2030_sens,supp2035_sens)
    sup_level[((i-1)*9+1):(i*9),"nnrti_res"]=rep(c(outcome2(data_sup_list2025[[1]]),
                                               outcome2(data_sup_list2030[[1]]),
                                               outcome2(data_sup_list2035[[1]])),each=3)
    sup_level[((i-1)*9+1):(i*9),c("nnrti_res_inf","nnrti_res_sup")]=matrix(c(rep(res2025_sens,3),
                                                                  rep(res2030_sens,3),
                                                                  rep(res2035_sens,3)),ncol=2,byrow = TRUE)
    sup_level[((i-1)*9+1):(i*9),"prop200"]=rep(c(outcome3(data_sup_list2025[[1]]),
                                                   outcome3(data_sup_list2030[[1]]),
                                                   outcome3(data_sup_list2035[[1]])),each=3)
    
    print(i)   
  }
  
  #Remove what is not needed
  # sup_level=sup_level[sup_level$sub_scen=="0) all on NNRTI" |
  #                     sup_level$sub_scen=="1) all men on DTG" |
  #                     sup_level$sub_scen=="3) all men and women not at risk of pregancy on DTG" |
  #                     sup_level$sub_scen== "4) all men and women on DTG",]
  sup_level=sup_level[sup_level$main_scen=="DTG as first-line and switch",]
  View(sup_level)
  View(sup_level[sup_level$year_start=="2035" & sup_level$year_treated==2,])
  View(sup_level[sup_level$year_start=="2035" & sup_level$year_treated=="After 2 years of ART",])
##############################################################################################################################
#Figure
  sup_level2=sup_level
  sup_level=sup_level[sup_level$year_treated==1 | sup_level$year_treated==2,]
  sup_level=sup_level[sup_level$year_start==2025 | sup_level$year_start==2035,]
  sup_level$year_start=as.factor(as.character(sup_level$year_start))
  sup_level$year_treated=as.factor(as.character(sup_level$year_treated))
  levels(sup_level$year_treated)=c("After 1 year of ART", "After 2 years of ART")
  
  pd = position_dodge(0.7)
  ggplot(sup_level,
         aes(x =year_start, 
             y =sup,
             color = sub_scen)) +
    geom_point(shape = 15, 
               size = 4, 
               position = pd) +
    facet_grid(. ~ year_treated)+
    
    geom_errorbar(aes(ymin = sup_inf,
                      ymax = sup_sup),
                  width = 0.2, 
                  size = 0.7, 
                  position = pd) +
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
    guides(colour=guide_legend(nrow=2,byrow=FALSE,
                               keywidth=1.5,))+
    theme(legend.text=element_text(size=11))
##############################################################################################################################
##############################################################################################################################
#heatmap  
  library(plotly)
  # with inverted x-axis, time is equally spaced
  
  #Load data
  i_len=39
  j_len=51
  data_heatmap=list()
  data_heatmap_fail=list()
  for(i in 1:i_len){
    for(j in 1:j_len){
      #load(paste("C:/Users/ahauser/Documents/Step2/65471464/heatmap",j_len*(i-1)+j,".RData",sep=""))
      #data_heatmap[[j_len*(i-1)+j]]=data_h
      #heatmap_4
      load(paste("C:/Users/ahauser/Documents/Step2/R_results/dtg_heatmap/heatmap_8/heatmap_",i,"_",j,".RData",sep=""))
      data_heatmap[[j_len*(i-1)+j]]=data_h
      #heatmap_fail_5
      load(paste("C:/Users/ahauser/Documents/Step2/R_results/dtg_heatmap/heatmap_fail_8/heatmap_",i,"_",j,".RData",sep=""))
      data_heatmap_fail[[j_len*(i-1)+j]]=data_h
    }
    print(i)
  }

  #Calculation of tdr among diagnosed people
  heatmap=matrix(sapply(data_heatmap,function(x) sum(x[35*12+1,select_15(c(2,12),1:4,2,1:2)])/
                               sum(x[35*12+1,select_15(c(2,12),1:4,1:2,1:2)])),c(j_len,i_len))
  heatmap_fail=matrix(sapply(data_heatmap_fail,function(x) sum(x[35*12+1,select_15(c(2,12),1:4,2,1:2)])/
                          sum(x[35*12+1,select_15(c(2,12),1:4,1:2,1:2)])),c(j_len,i_len))
  
##################################################################################################################
  
  min1=min(heatmap)
  min2=min(heatmap_fail)
  max1=max(heatmap)
  max2=max(heatmap_fail)
  min=min(min1,min2)
  max=max(max1,max2)
  
##################################################################################################################
  #Plot 1
  vals1 <- unique(scales::rescale(c(heatmap),to=c((min1-min)/(max-min),(max1-min)/(max-min))))
  vals1_01 <- unique(scales::rescale(c(heatmap),to=c(0,1)))
  o <- order(vals1_01, decreasing = FALSE)
  cols <- scales::col_numeric("Blues", domain = c(0,1))(vals1)
  colz1 <- setNames(data.frame((vals1_01)[o], cols[o]), NULL)
  
  #which_x<-c(seq(1,i_len,by=5),i_len)
  which_x<-which(seq(0.5,10,length.out = 39) %in% c(0.5,2,4,6,8,10))
  
  plot_dtg<-plot_ly(z=heatmap,colorscale = colz1,colorbar = list(thickness=30,len=0.1,title = "",titlefont=list(size=2),tickfont=list(size=2)), type = "contour") %>%
    plotly::layout(title="    A. Suppressed ind.",
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
                   xaxis=list(title="Time to switch to DTG (year)",list(thickness=1,len=0.1,title = "NNRTI PDR \n level in 2035",titlefont=list(size=1),tickfont=list(size=1)),titlefont=list(family="Arial",size=24),ticktext=seq(0.5,10,length.out = 39)[which_x],tickvals=which_x-1,tickmode="array",tickfont=list(size=18)),
                   yaxis=list(title="% women eligible for DTG",titlefont=list(family="Arial",size=24),ticktext=seq(0,100,by=10),tickvals=seq(0,50,by=5),tickmode="array",tickfont=list(size=18))
    )%>%
    colorbar(limits = c(min,max))
  #2 Plots
  subplot(plot_dtg,plot_dtg_fail,margin=0.08,titleX=TRUE,shareX = TRUE,titleY=TRUE,shareY = FALSE)
##################################################################################################################
  
  
  plot_dtg_fail<-plot_ly(z=heatmap_fail,colorscale = colz1,colorbar=list(thickness=60,len=0.7,title = "NNRTI PDR \n level in 2035",titlefont=list(size=16),tickfont=list(size=18)), type = "contour") %>%
    plotly::layout(title="                                                    ",
                   xaxis=list(title="Time to switch to DTG (year)",list(thickness=1,len=0.1,title = "NNRTI PDR \n level in 2035",titlefont=list(size=1),tickfont=list(size=1)),titlefont=list(family="Arial",size=24),ticktext=seq(0.5,10,length.out = 39)[c(seq(1,i_len,by=5),i_len)],tickvals=c(seq(1,i_len,by=5),i_len)-1,tickmode="array",tickfont=list(size=18)),
                   yaxis=list(title="% women eligible for DTG",titlefont=list(family="Arial",size=24),ticktext=seq(0,100,by=10),tickvals=seq(0,50,by=5),tickmode="array",tickfont=list(size=18))
    )%>%
    colorbar(limits = c(min,max))
  
  #Plot 1
  vals1 <- unique(scales::rescale(c(heatmap),to=c((min1-min)/(max-min),(max1-min)/(max-min))))
  vals1_01 <- unique(scales::rescale(c(heatmap),to=c(0,1)))
  o <- order(vals1_01, decreasing = FALSE)
  cols <- scales::col_numeric("Blues", domain = c(0,1))(vals1)
  colz1 <- setNames(data.frame((vals1_01)[o], cols[o]), NULL)
  
  plot_dtg<-plot_ly(z=heatmap,colorscale = colz1,colorbar = list(thickness=30,len=0.1,title = "",titlefont=list(size=2),tickfont=list(size=2)), type = "contour") %>%
    plotly::layout(title="    A. Suppressed ind.",
                   xaxis=list(title="Time to switch to DTG (year)",list(thickness=1,len=0.1,title = "NNRTI PDR \n level in 2035",titlefont=list(size=1),tickfont=list(size=1)),titlefont=list(family="Arial",size=24),ticktext=seq(0.5,10,length.out = 39)[c(seq(1,i_len,by=5),i_len)],tickvals=c(seq(1,i_len,by=5),i_len)-1,tickmode="array",tickfont=list(size=18)),
                   yaxis=list(title="% women eligible for DTG",titlefont=list(family="Arial",size=24),ticktext=seq(0,100,by=10),tickvals=seq(0,50,by=5),tickmode="array",tickfont=list(size=18))
    )%>%
    colorbar(limits = c(min,max))
  
  #Plot 2
  vals2 <- unique(scales::rescale(c(heatmap_fail),to=c((min2-min)/(max-min),(max2-min)/(max-min))))
  vals2_01 <- unique(scales::rescale(c(heatmap_fail),to=c(0,1)))
  o <- order(vals2_01, decreasing = FALSE)
  cols <- scales::col_numeric("Blues", domain = c(0,1))(vals2)
  colz2 <- setNames(data.frame((vals2_01)[o], cols[o]), NULL)
  
  
  plot_dtg_fail<-plot_ly(z=heatmap_fail,colorscale = colz2,colorbar=list(thickness=60,len=0.7,title = "NNRTI PDR \n level in 2035",titlefont=list(size=16),tickfont=list(size=18)), type = "contour") %>%
    plotly::layout(title="                                                    ",
                   xaxis=list(title="Time to switch to DTG (year)",titlefont=list(family="Arial",size=24),ticktext=c(1,5,10,15,18),tickvals=(c(1,5,10,15,18)-1)/17*38,tickmode="array",tickfont=list(size=18)),
                   yaxis=list(title="% women eligible for DTG",titlefont=list(family="Arial",size=24),ticktext=seq(0,100,by=10),tickvals=seq(0,50,by=5),tickmode="array",tickfont=list(size=18))
    )%>%
    colorbar(limits = c(min,max))

  
  #2 Plots
  subplot(plot_dtg,plot_dtg_fail,margin=0.08,titleX=TRUE,shareX = TRUE,titleY=TRUE,shareY = FALSE)
##################################################################################################################
  
  
  vals1 <- unique(scales::rescale(c(heatmap),to=c(0,1)))
  vals2 <- unique(scales::rescale(c(heatmap),to=c(0,0.1)))
  o <- order(vals1, decreasing = FALSE)
  cols <- scales::col_numeric("Blues", domain = c(0,1))(vals2)
  colz1 <- setNames(data.frame((vals1)[o], cols[o]), NULL)
  
  
  
  o <- order(vals, decreasing = FALSE)
  cols <- scales::col_numeric("Blues", domain = c(0,1))(vals1)
  colz1 <- setNames(data.frame((vals1_01)[o], cols[o]), NULL)
  
  
  
  vals1 <- unique(scales::rescale(c(heatmap),to=c(0,1)))
  vals2 <- unique(scales::rescale(c(heatmap),to=c(0,0.1)))
  o <- order(vals1, decreasing = FALSE)
  cols <- scales::col_numeric("Blues", domain = c(0,1))(vals2)
  colz1 <- setNames(data.frame((vals1)[o], cols[o]), NULL)
  

  #Figure
  vals <- unique(scales::rescale(c(heatmap)))
  o <- order(vals, decreasing = FALSE)
  cols <- scales::col_numeric("Blues", domain = NULL)(vals)
  colz <- setNames(data.frame(vals[o], cols[o]), NULL)
  plot_dtg<-plot_ly(z=heatmap,colorscale = colz1, type = "contour") %>%
    plotly::layout(title="    A. Suppressed ind.",
                   xaxis=list(title="Time to switch to DTG (year)",list(thickness=1,len=0.1,title = "NNRTI PDR \n level in 2035",titlefont=list(size=1),tickfont=list(size=1)),titlefont=list(family="Arial",size=24),ticktext=seq(0.5,10,length.out = 39)[c(seq(1,i_len,by=5),i_len)],tickvals=c(seq(1,i_len,by=5),i_len)-1,tickmode="array",tickfont=list(size=18)),
                   yaxis=list(title="% women eligible for DTG",titlefont=list(family="Arial",size=24),ticktext=seq(0,100,by=10),tickvals=seq(0,50,by=5),tickmode="array",tickfont=list(size=18))
    )
  plot_dtg

  
  vals <- unique(scales::rescale(c(heatmap_fail)))
  o <- order(vals, decreasing = FALSE)
  cols <- scales::col_numeric("Blues", domain = NULL)(vals)
  colz <- setNames(data.frame(vals[o], cols[o]), NULL)
  plot_dtg_fail<-plot_ly(z=heatmap_fail,colorscale = colz,colorbar = list(thickness=60,len=0.7,title = "NNRTI PDR \n level in 2035",titlefont=list(size=20),tickfont=list(size=18)), type = "contour") %>%
    plotly::layout(title="                                                    ",
                   xaxis=list(title="Time to switch to DTG (year)",titlefont=list(family="Arial",size=24),ticktext=c(1,5,10,15,18),tickvals=(c(1,5,10,15,18)-1)/17*38,tickmode="array",tickfont=list(size=18)),
                   yaxis=list(title="% women eligible for DTG",titlefont=list(family="Arial",size=24),ticktext=seq(0,100,by=10),tickvals=seq(0,50,by=5),tickmode="array",tickfont=list(size=18))
                   )
  plot_dtg_fail
  
  plot_grid(plot_dtg,plot_dtg_fail, labels = "AUTO")
  
  subplot(plot_dtg,plot_dtg_fail,margin=0.08)
  
  
  
  
  
  
  #Figure
  vals <- unique(scales::rescale(c(heatmap)))
  o <- order(vals, decreasing = FALSE)
  cols <- scales::col_numeric("Blues", domain = NULL)(vals)
  colz <- setNames(data.frame(vals[o], cols[o]), NULL)
  plot_dtg<-plot_ly(z=heatmap,colorscale = colz, type = "contour") %>%
    plotly::layout(xaxis=list(title="Time to switch to DTG (year)",list(thickness=1,len=0.1,title = "",titlefont=list(size=1),tickfont=list(size=1)),titlefont=list(family="Arial",size=24),ticktext=seq(0.5,10,length.out = 39)[c(seq(1,i_len,by=5),i_len)],tickvals=c(seq(1,i_len,by=5),i_len)-1,tickmode="array",tickfont=list(size=18)),
                   yaxis=list(title="% women eligible for DTG",titlefont=list(family="Arial",size=24),ticktext=seq(0,100,by=10),tickvals=seq(0,50,by=5),tickmode="array",tickfont=list(size=18))
    )%>%
    colorbar(limits = c(min,max))
  plot_dtg
  
  library(fields)
  contour(heatmap, add = TRUE)
  
  vals <- unique(scales::rescale(c(heatmap)))
  o <- order(vals, decreasing = FALSE)
  cols <- scales::col_numeric("Blues", domain = NULL)(vals)
  colz <- setNames(data.frame(vals[o], cols[o]), NULL)
  plot_dtg_fail<-plot_ly(z=heatmap_fail,colorscale = colz,colorbar = list(thickness=60,len=0.7,title = "NNRTI PDR \n level in 2035",titlefont=list(size=20),tickfont=list(size=18)), type = "contour") %>%
    plotly::layout(xaxis=list(title="Time to switch to DTG (year)",titlefont=list(family="Arial",size=24),ticktext=c(1,5,10,15,18),tickvals=(c(1,5,10,15,18)-1)/17*38,tickmode="array",tickfont=list(size=18)),
                   yaxis=list(title="% women eligible for DTG",titlefont=list(family="Arial",size=24),ticktext=seq(0,100,by=10),tickvals=seq(0,50,by=5),tickmode="array",tickfont=list(size=18))
    )%>%
    colorbar(limits = c(min,max))
  plot_dtg_fail
  
  subplot(plot_dtg,plot_dtg_fail)
  