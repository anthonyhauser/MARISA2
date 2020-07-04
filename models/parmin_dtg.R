###########################################################################################
#DTG, not more used

#11 care stages, not more used
select_dtg=function(r1,r2,r3,r4){
  l1=11
  l2=4
  l3=2
  l4=2
  ar=array(1:(l1*l2*l3*l4),dim=c(l1,l2,l3,l4))
  return(as.vector(ar[r1,r2,r3,r4]))
}
mod_cpp_dtg=function(t,x,p1,treat_dtg){
  #mod =function(t,x, params){
  with(as.list(params), {
    x=x[1:(88*2)]
    ##################################################################################################################
    #Optimization
    rate1_inf <- p1[1] #number of unprotected sexual acts/month
    rate2_inf <- p1[2]
    rate1_diag <- p1[3] #diagnosis rate
    rate2_diag <- p1[4] #diagnosis rate
    rate_treat <- p1[5] #scale parameter for overall treatment rate
    rate_death <- p1[6] #scale parameter for overall death
    q <- p2[7]
    #Range
    rate_ratio <- p2[8] #varying parameter to estimate death for treated (but neither suppressed nor failed)
    
    alpha <- p2[12] #ratio between outcome when resistant/sensitive
    rate_res <- p2[13] #resistant rate
    rate_susc <- p2[14] #reversion rate
    
    #1. Infection and Diagnosis
    RateInf1=array(rep(RateTransm,each=4),dim=c(4,2,2))*rate1_inf
    RateInf2=RateInf1*rate2_inf
    
    Diag_frame[,4]=Diag_frame[,4]*rate2_diag
    RateDiag=1/rate1_diag*apply(Diag_frame,1,function(x) smooth_st(2005+t/12,x["a"],x["b"],x["c"],x["d"]))
    
    RateDiag=array(rep(RateDiag,2)*rep(c(1,diag_women),each=4),dim=c(4,2))+
      array(rep(oi_inc,2)*rep(smooth_st(2005+t/12,oi_test["a"],oi_test["b"],oi_test["c"],oi_test["d"]),8),dim=c(4,2))+
      array(c(rep(0,4),preg_inc)*rep(smooth_st(2005+t/12,preg_test["a"],preg_test["b"],preg_test["c"],preg_test["d"]),8),dim=c(4,2))
    
    #2. Treatment theoretical (if no counterfactual scenario)
    RateTreatFirst_th=apply(Treat_frame,1,function(x) smooth_st(2005+t/12,x["a"],x["b"],x["c"],x["d"]))*smooth_st(2005+t/12,2001,2012,1,200/12)*rate_treat
    
    RateTreatFirst_f=function(y){
      if(y<2005){
        return(apply(Treat_frame_original,1,function(x) smooth_st(y,x["a"],x["b"],x["c"],x["d"]))*smooth_st(y,2001,2012,1,200/12)*rate_treat)
      }else{
        return(apply(
          data.frame(t1=apply(Treat_frame,1,function(x) lin_st(y,x["a"],x["b"],x["c"],x["d"])),
                     t2=apply(Treat_frame2,1,function(x) lin_st(y,x["a"],x["b"],x["c"],x["d"]))*as.numeric(y>Treat_frame2$a[4])),
          1,max)*smooth_st(y,2001,2012,1,200/12)*rate_treat)
      }
    }
    RateTreatFirst_th=RateTreatFirst_f(2005+t/12)
    
    RateDirectTreatSecond_th=RateDirectTreatSecond
    
    #now
    RateTreatFirst=matrix(c(RateTreatFirst_th,RateTreatFirst_th),nrow=2,ncol=4,byrow=T)
    RateDirectTreatSecond=matrix(c(RateDirectTreatSecond_th,RateDirectTreatSecond_th),nrow=2,ncol=4,byrow=T)
    
    #Resistance testing at baseline
    #RateTreatFirst=matrix(c(RateTreatFirst_th,c(0,0,0,0)),nrow=2,ncol=4,byrow=T)
    #RateDirectTreatSecond=matrix(c(RateDirectTreatSecond_th,RateTreatFirst_th+RateDirectTreatSecond_th),nrow=2,ncol=4,byrow=T)
    #RateTreatFirst=rep(rate_treat,4)
    
    #2. Treatment DTG
    #RateTreatFirstDTG[res,cd4,gender]
    RateTreatFirstDTG=array(rep(RateTreatFirst*treat_dtg[3],2)*rep(c(1,treat_dtg[1]),each=8),dim=c(2,4,2))*smooth_st(2005+t/12,2018,2018+treat_dtg["delay_first"],0,1)
    #RateSwitchDTG[care,cd4,res,gender]
    RateSwitchDTG=array(rep(rep(treat_dtg[4]*RateTreatSecond,each=3)*rep(c(treat_dtg[5],treat_dtg[6],treat_dtg[7]),4),4)*
                          #rep(c(1,treat_dtg[2]),each=24),dim=c(3,4,2,2))*
                          rep(c(1,1),each=24),dim=c(3,4,2,2))*
      smooth_st(2005+t/12,2018,2018+treat_dtg["delay_switch"],0,1)
    RateSwitchDTG[1:2,1:4,1:2,1:2]=rep(treat_dtg[4],32)*rep(rep(1/12,2)*c(treat_dtg[5],treat_dtg[6]),16)*
      #rep(c(1,treat_dtg[2]),each=16)*
      rep(c(1,1),each=16)*
      smooth_st(2005+t/12,2018,2018+treat_dtg["delay_switch"],0,1)
    #RateTreatFirst [res,cd4,gender]
    RateTreatFirst=array(rep(RateTreatFirst,2),dim=c(2,4,2))*apply(array(rep(1-(c(1,treat_dtg[1])*smooth_st(2005+t/12,2018,2018+treat_dtg["delay_first"],0,1)*treat_dtg[3]),each=8),dim=c(2,4,2)),c(1,2,3),function(x) max(x,0))
    #RateTreatSecond[CD4,gender]
    RateTreatSecond=rep(RateTreatSecond,2)*apply(array(rep(1-(c(1,treat_dtg[2])*smooth_st(2005+t/12,2018,2018+treat_dtg["delay_switch"],0,1)*treat_dtg[4]),each=4),dim=c(4,2)),c(1,2),function(x) max(x,0))
    # RateTreatFirst=matrix(rep(0,8),nrow=2,ncol=4,byrow=T)
    # RateDirectTreatSecond=matrix(rep(0,8),nrow=2,ncol=4,byrow=T)
    #3. Death
    #p1=p1_f(t,q)
    #p1=min(1,q*mean(sapply(2005+seq(t-36,t,1)/12,function(x) apply(Treat_frame[4,],1,function(y) smooth_st(x,y["a"],y["b"],y["c"],y["d"]))*smooth_st(x,2001,2012,1,200/12))))
    p1=min(1,q*mean(sapply(2005+seq(t-36,t,1)/12,function(x) RateTreatFirst_f(x)[4]/rate_treat)))
    #p1=p1/2+0.2
    mu[1:2,4,1]=p1*40.9+(1-p1)*134.4
    mu[4,4,1]=p1*8.3+(1-p1)*41.7
    mu[5,4,1]=p1*11.8+(1-p1)*59.7
    mu[3,1:4,1]=rate_ratio*mu[4,1:4,1]+(1-rate_ratio)*mu[5,1:4,1]
    mu[6:8,1:4,1]=mu[3:5,1:4,1]
    mu[1:8,1:4,2]=mu[1:8,1:4,1]
    mu=mu*rate_death/1000
    
    #Back and forth mutation and impact on treatment outcome
    RateResistant=1/rate_res
    RateSusceptible=1/rate_susc
    #From T to S
    RateSuppFirstSusc=RateSuppFirst/1
    RateSuppFirstResis=1/alpha*RateSuppFirst/1
    #From S to F
    RateFailFirstSusc=RateFailFirst*1
    RateFailFirstResis=alpha*RateFailFirst*1
    #From T to F
    RateTreatToFailFirstSusc=RateTreatToFailFirst*1
    RateTreatToFailFirstResis=alpha*RateTreatToFailFirst*1
    #From F to S
    RateFailToSuppTreatFirstSusc=RateFailToSuppTreatFirst/1
    RateFailToSuppTreatFirstResis=1/alpha*RateFailToSuppTreatFirst/1
    
    #For second-line regimen, no effect of NNRTI resistance on treatment as it is PI
    ##################################################################################################################
    dx=array(0,dim=c(11,4,2,2))
    x=array(x,dim=c(11,4,2,2))
    I=apply(x,4,sum)
    
    #First approach to model S: S=N-I
    N=N_value_f(t+1)
    S=N-I
    #Percentage getting NNRTI-based 1st-line[gender]
    p_NNRTI=sapply(1-(c(1,treat_dtg[1])*smooth_st(2005+t/12,2018,2018+treat_dtg["delay_first"],0,1)*treat_dtg[3]),function(x) max(x,0))
    #infected and treated children
    I_m15=rep(infected_m15_f(t+1)/2,2) #number of infected children reaching 15 years old at a given month
    t_m15=rep(treat_m15_f(t+1)/2,2)
    
    dx[2,1:4,1,1:2]=rep((1-res_m15)*I_m15/4,each=4)
    dx[2,1:4,2,1:2]=rep(res_m15*I_m15/4,each=4)
    dx[3,1:4,1,1:2]=rep((1-res_m15)*t_m15*p_NNRTI/4,each=4)
    dx[3,1:4,2,1:2]=rep(res_m15*t_m15*p_NNRTI/4,each=4)
    dx[9,1:4,1,1:2]=rep((1-res_m15)*t_m15*(1-p_NNRTI)/4,each=4)
    dx[9,1:4,2,1:2]=rep(res_m15*t_m15*(1-p_NNRTI)/4,each=4)
    
    for(j in 1:2){
      #dx[x1,x2,x3,x4] x1:care stages x2:disease progression x3:resistance x4:gender
      #dx with only interaction with the neighbours in the compartmental model
      dx[1,1,1,j]=S[j]/N[j]*(sum(RateInf1[1:4,1:2,j]*apply(x[c(1),1:4,1,1:2],c(1,2),sum))+sum(RateInf2[1:4,1:2,j]*apply(x[c(2,3,5,6,8,9,11),1:4,1,1:2],c(2,3),sum)))-
        (RateStageInf[1]+RateDiag[1,j]+mu[1,1,1])*x[1,1,1,j]
      dx[2,1,1,j]=dx[2,1,1,j]+RateDiag[1,j]*x[1,1,1,j]-(RateTreatFirst[1,1,j]+RateStageDiag[1]+mu[2,1,1])*x[2,1,1,j]
      dx[3,1,1,j]=dx[3,1,1,j]+RateTreatFirst[1,1,j]*x[2,1,1,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,1])*x[3,1,1,j]+RateStageTreatFirstL[1]*x[3,2,1,j]
      dx[4,1,1,j]=RateSuppFirstSusc[1]*x[3,1,1,j]-(RateFailFirstSusc[1]+mu[4,1,1])*x[4,1,1,j]+RateStageSuppFirst[1]*x[4,2,1,j]
      dx[5,1,1,j]=RateFailFirstSusc[1]*x[4,1,1,j]-(RateTreatSecond[1,j]+RateStageFailFirst[1]+mu[5,1,1])*x[5,1,1,j]
      dx[6,1,1,j]=RateTreatSecond[1,j]*x[5,1,1,j]-(RateSuppSecond[1]+RateStageTreatSecondR[1]+mu[6,1,1])*x[6,1,1,j]+RateStageTreatSecondL[1]*x[6,2,1,j]
      dx[7,1,1,j]=RateSuppSecond[1]*x[6,1,1,j]-(RateFailSecond[1]+mu[7,1,1])*x[7,1,1,j]+RateStageSuppSecond[1]*x[7,2,1,j]
      dx[8,1,1,j]=RateFailSecond[1]*x[7,1,1,j]-(RateStageFailSecond[1]+mu[8,1,1])*x[8,1,1,j]
      
      dx[1,2,1,j]=-(RateStageInf[2]+RateDiag[2,j]+mu[1,2,1])*x[1,2,1,j]+RateStageInf[1]*x[1,1,1,j]
      dx[2,2,1,j]=dx[2,2,1,j]+RateDiag[2,j]*x[1,2,1,j]-(RateTreatFirst[1,2,j]+RateStageDiag[2]+mu[2,2,1])*x[2,2,1,j]+RateStageDiag[1]*x[2,1,1,j]
      dx[3,2,1,j]=dx[3,2,1,j]+RateTreatFirst[1,2,j]*x[2,2,1,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,1])*x[3,2,1,j]+RateStageTreatFirstR[1]*x[3,1,1,j]+RateStageTreatFirstL[2]*x[3,3,1,j]
      dx[4,2,1,j]=RateSuppFirstSusc[2]*x[3,2,1,j]-(RateFailFirstSusc[2]+RateStageSuppFirst[1]+mu[4,2,1])*x[4,2,1,j]+RateStageSuppFirst[2]*x[4,3,1,j]
      dx[5,2,1,j]=RateFailFirstSusc[2]*x[4,2,1,j]-(RateTreatSecond[2,j]+RateStageFailFirst[2]+mu[5,2,1])*x[5,2,1,j]+RateStageFailFirst[1]*x[5,1,1,j]
      dx[6,2,1,j]=RateTreatSecond[2,j]*x[5,2,1,j]-(RateSuppSecond[2]+RateStageTreatSecondR[2]+RateStageTreatSecondL[1]+mu[6,2,1 ])*x[6,2,1,j]+RateStageTreatSecondR[1]*x[6,1,1,j]+RateStageTreatSecondL[2]*x[6,3,1,j]
      dx[7,2,1,j]=RateSuppSecond[2]*x[6,2,1,j]-(RateFailSecond[2]+RateStageSuppSecond[1]+mu[7,2,1])*x[7,2,1,j]+RateStageSuppSecond[2]*x[7,3,1,j]
      dx[8,2,1,j]=RateFailSecond[2]*x[7,2,1,j]-(RateStageFailSecond[2]+mu[8,2,1])*x[8,2,1,j]+RateStageFailSecond[1]*x[8,1,1,j]
      
      dx[1,3,1,j]=-(RateStageInf[3]+RateDiag[3,j]+mu[1,3,1])*x[1,3,1,j]+RateStageInf[2]*x[1,2,1,j]
      dx[2,3,1,j]=dx[2,3,1,j]+RateDiag[3,j]*x[1,3,1,j]-(RateTreatFirst[1,3,j]+RateStageDiag[3]+mu[2,3,1])*x[2,3,1,j]+RateStageDiag[2]*x[2,2,1,j]
      dx[3,3,1,j]=dx[3,3,1,j]+RateTreatFirst[1,3,j]*x[2,3,1,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,1])*x[3,3,1,j]+RateStageTreatFirstR[2]*x[3,2,1,j]+RateStageTreatFirstL[3]*x[3,4,1,j]
      dx[4,3,1,j]=RateSuppFirstSusc[3]*x[3,3,1,j]-(RateFailFirstSusc[3]+RateStageSuppFirst[2]+mu[4,3,1])*x[4,3,1,j]+RateStageSuppFirst[3]*x[4,4,1,j]
      dx[5,3,1,j]=RateFailFirstSusc[3]*x[4,3,1,j]-(RateTreatSecond[3,j]+RateStageFailFirst[3]+mu[5,3,1])*x[5,3,1,j]+RateStageFailFirst[2]*x[5,2,1,j]
      dx[6,3,1,j]=RateTreatSecond[3,j]*x[5,3,1,j]-(RateSuppSecond[3]+RateStageTreatSecondR[3]+RateStageTreatSecondL[2]+mu[6,3,1])*x[6,3,1,j]+RateStageTreatSecondR[2]*x[6,2,1,j]+RateStageTreatSecondL[3]*x[6,4,1,j]
      dx[7,3,1,j]=RateSuppSecond[3]*x[6,3,1,j]-(RateFailSecond[3]+RateStageSuppSecond[2]+mu[7,3,1])*x[7,3,1,j]+RateStageSuppSecond[3]*x[7,4,1,j]
      dx[8,3,1,j]=RateFailSecond[3]*x[7,3,1,j]-(RateStageFailSecond[3]+mu[8,3,1])*x[8,3,1,j]+RateStageFailSecond[2]*x[8,2,1,j]
      
      dx[1,4,1,j]=-(RateDiag[4,j]+mu[1,4,1])*x[1,4,1,j]+RateStageInf[3]*x[1,3,1,j]
      dx[2,4,1,j]=dx[2,4,1,j]+RateDiag[4,j]*x[1,4,1,j]-(RateTreatFirst[1,4,j]+mu[2,4,1])*x[2,4,1,j]+RateStageDiag[3]*x[2,3,1,j]
      dx[3,4,1,j]=dx[3,4,1,j]+RateTreatFirst[1,4,j]*x[2,4,1,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,1])*x[3,4,1,j]+RateStageTreatFirstR[3]*x[3,3,1,j]
      dx[4,4,1,j]=RateSuppFirstSusc[4]*x[3,4,1,j]-(RateFailFirstSusc[4]+RateStageSuppFirst[3]+mu[4,4,1])*x[4,4,1,j]
      dx[5,4,1,j]=RateFailFirstSusc[4]*x[4,4,1,j]-(RateTreatSecond[4,j]+mu[5,4,1])*x[5,4,1,j]+RateStageFailFirst[3]*x[5,3,1,j]
      dx[6,4,1,j]=RateTreatSecond[4,j]*x[5,4,1,j]-(RateSuppSecond[4]+RateStageTreatSecondL[3]+mu[6,4,1])*x[6,4,1,j]+RateStageTreatSecondR[3]*x[6,3,1,j]
      dx[7,4,1,j]=RateSuppSecond[4]*x[6,4,1,j]-(RateFailSecond[4]+RateStageSuppSecond[3]+mu[7,4,1])*x[7,4,1,j]
      dx[8,4,1,j]=RateFailSecond[4]*x[7,4,1,j]-mu[8,4,1 ]*x[8,4,1,j]+RateStageFailSecond[3]*x[8,3,1,j]
      
      ##############################################################################################################################################################
      
      dx[1,1,1,j]=dx[1,1,1,j]
      dx[2,1,1,j]=dx[2,1,1,j]-RateDirectTreatSecond[1,1]*x[2,1,1,j]+RateStopTreatFirst[1]*x[3,1,1,j]+RateStopSuppFirst[1]*x[4,1,1,j]+RateStopFailFirst[1]*x[5,1,1,j]+
        RateStopTreatSecond[1]*x[6,1,1,j]+RateStopSuppSecond[1]*x[7,1,1,j]+RateStopFailSecond[1]*x[8,1,1,j]
      dx[3,1,1,j]=dx[3,1,1,j]-RateStopTreatFirst[1]*x[3,1,1,j]-RateTreatToFailFirstSusc[1]*x[3,1,1,j]
      dx[4,1,1,j]=dx[4,1,1,j]-RateStopSuppFirst[1]*x[4,1,1,j]+RateFailToSuppTreatFirstSusc[1]*x[5,1,1,j]-RateSuppFirstToSecond[1]*x[4,1,1,j]
      dx[5,1,1,j]=dx[5,1,1,j]-RateStopFailFirst[1]*x[5,1,1,j]-RateFailToSuppTreatFirstSusc[1]*x[5,1,1,j]+RateTreatToFailFirstSusc[1]*x[3,1,1,j]
      dx[6,1,1,j]=dx[6,1,1,j]-RateStopTreatSecond[1]*x[6,1,1,j]+RateDirectTreatSecond[1,1]*x[2,1,1,j]-RateTreatToFailSecond[1]*x[6,1,1,j]
      dx[7,1,1,j]=dx[7,1,1,j]-RateStopSuppSecond[1]*x[7,1,1,j]+RateFailToSuppTreatSecond[1]*x[8,1,1,j]+RateSuppFirstToSecond[1]*x[4,1,1,j]
      dx[8,1,1,j]=dx[8,1,1,j]-RateStopFailSecond[1]*x[8,1,1,j]-RateFailToSuppTreatSecond[1]*x[8,1,1,j]+RateTreatToFailSecond[1]*x[6,1,1,j]
      
      dx[1,2,1,j]=dx[1,2,1,j]
      dx[2,2,1,j]=dx[2,2,1,j]-RateDirectTreatSecond[1,2]*x[2,2,1,j]+RateStopTreatFirst[2]*x[3,2,1,j]+RateStopSuppFirst[2]*x[4,2,1,j]+RateStopFailFirst[2]*x[5,2,1,j]+
        RateStopTreatSecond[2]*x[6,2,1,j]+RateStopSuppSecond[2]*x[7,2,1,j]+RateStopFailSecond[2]*x[8,2,1,j]
      dx[3,2,1,j]=dx[3,2,1,j]-RateStopTreatFirst[2]*x[3,2,1,j]-RateTreatToFailFirstSusc[2]*x[3,2,1,j]
      dx[4,2,1,j]=dx[4,2,1,j]-RateStopSuppFirst[2]*x[4,2,1,j]+RateFailToSuppTreatFirstSusc[2]*x[5,2,1,j]-RateSuppFirstToSecond[2]*x[4,2,1,j]
      dx[5,2,1,j]=dx[5,2,1,j]-RateStopFailFirst[2]*x[5,2,1,j]-RateFailToSuppTreatFirstSusc[2]*x[5,2,1,j]+RateTreatToFailFirstSusc[2]*x[3,2,1,j]
      dx[6,2,1,j]=dx[6,2,1,j]-RateStopTreatSecond[2]*x[6,2,1,j]+RateDirectTreatSecond[1,2]*x[2,2,1,j]-RateTreatToFailSecond[2]*x[6,2,1,j]
      dx[7,2,1,j]=dx[7,2,1,j]-RateStopSuppSecond[2]*x[7,2,1,j]+RateFailToSuppTreatSecond[2]*x[8,2,1,j]+RateSuppFirstToSecond[2]*x[4,2,1,j]
      dx[8,2,1,j]=dx[8,2,1,j]-RateStopFailSecond[2]*x[8,2,1,j]-RateFailToSuppTreatSecond[2]*x[8,2,1,j]+RateTreatToFailSecond[2]*x[6,2,1,j]
      
      dx[1,3,1,j]=dx[1,3,1,j]
      dx[2,3,1,j]=dx[2,3,1,j]-RateDirectTreatSecond[1,3]*x[2,3,1,j]+RateStopTreatFirst[3]*x[3,3,1,j]+RateStopSuppFirst[3]*x[4,3,1,j]+RateStopFailFirst[3]*x[5,3,1,j]+
        RateStopTreatSecond[3]*x[6,3,1,j]+RateStopSuppSecond[3]*x[7,3,1,j]+RateStopFailSecond[3]*x[8,3,1,j]
      dx[3,3,1,j]=dx[3,3,1,j]-RateStopTreatFirst[3]*x[3,3,1,j]-RateTreatToFailFirstSusc[3]*x[3,3,1,j]
      dx[4,3,1,j]=dx[4,3,1,j]-RateStopSuppFirst[3]*x[4,3,1,j]+RateFailToSuppTreatFirstSusc[3]*x[5,3,1,j]-RateSuppFirstToSecond[3]*x[4,3,1,j]
      dx[5,3,1,j]=dx[5,3,1,j]-RateStopFailFirst[3]*x[5,3,1,j]-RateFailToSuppTreatFirstSusc[3]*x[5,3,1,j]+RateTreatToFailFirstSusc[3]*x[3,3,1,j]
      dx[6,3,1,j]=dx[6,3,1,j]-RateStopTreatSecond[3]*x[6,3,1,j]+RateDirectTreatSecond[1,3]*x[2,3,1,j]-RateTreatToFailSecond[3]*x[6,3,1,j]
      dx[7,3,1,j]=dx[7,3,1,j]-RateStopSuppSecond[3]*x[7,3,1,j]+RateFailToSuppTreatSecond[3]*x[8,3,1,j]+RateSuppFirstToSecond[3]*x[4,3,1,j]
      dx[8,3,1,j]=dx[8,3,1,j]-RateStopFailSecond[3]*x[8,3,1,j]-RateFailToSuppTreatSecond[3]*x[8,3,1,j]+RateTreatToFailSecond[3]*x[6,3,1,j]
      
      dx[1,4,1,j]=dx[1,4,1,j]
      dx[2,4,1,j]=dx[2,4,1,j]-RateDirectTreatSecond[1,4]*x[2,4,1,j]+RateStopTreatFirst[4]*x[3,4,1,j]+RateStopSuppFirst[4]*x[4,4,1,j]+RateStopFailFirst[4]*x[5,4,1,j]+
        RateStopTreatSecond[4]*x[6,4,1,j]+RateStopSuppSecond[4]*x[7,4,1,j]+RateStopFailSecond[4]*x[8,4,1,j]
      dx[3,4,1,j]=dx[3,4,1,j]-RateStopTreatFirst[4]*x[3,4,1,j]-RateTreatToFailFirstSusc[4]*x[3,4,1,j]
      dx[4,4,1,j]=dx[4,4,1,j]-RateStopSuppFirst[4]*x[4,4,1,j]+RateFailToSuppTreatFirstSusc[4]*x[5,4,1,j]-RateSuppFirstToSecond[4]*x[4,4,1,j]
      dx[5,4,1,j]=dx[5,4,1,j]-RateStopFailFirst[4]*x[5,4,1,j]-RateFailToSuppTreatFirstSusc[4]*x[5,4,1,j]+RateTreatToFailFirstSusc[4]*x[3,4,1,j]
      dx[6,4,1,j]=dx[6,4,1,j]-RateStopTreatSecond[4]*x[6,4,1,j]+RateDirectTreatSecond[1,4]*x[2,4,1,j]-RateTreatToFailSecond[4]*x[6,4,1,j]
      dx[7,4,1,j]=dx[7,4,1,j]-RateStopSuppSecond[4]*x[7,4,1,j]+RateFailToSuppTreatSecond[4]*x[8,4,1,j]+RateSuppFirstToSecond[4]*x[4,4,1,j]
      dx[8,4,1,j]=dx[8,4,1,j]-RateStopFailSecond[4]*x[8,4,1,j]-RateFailToSuppTreatSecond[4]*x[8,4,1,j]+RateTreatToFailSecond[4]*x[6,4,1,j]
      #############################################################################################################################################################################
      #############################################################################################################################################################################
      dx[1,1,2,j]=S[j]/N[j]*(sum(RateInf1[1:4,1:2,j]*apply(x[c(1),1:4,2,1:2],c(1,2),sum))+sum(RateInf2[1:4,1:2,j]*apply(x[c(2,3,5,6,8,9,11),1:4,2,1:2],c(2,3),sum)))-
        (RateStageInf[1]+RateDiag[1,j]+mu[1,1,2 ])*x[1,1,2,j]
      dx[2,1,2,j]=dx[2,1,2,j]+RateDiag[1,j]*x[1,1,2,j]-(RateTreatFirst[2,1,j]+RateStageDiag[1]+mu[2,1,2 ])*x[2,1,2,j]
      dx[3,1,2,j]=dx[3,1,2,j]+RateTreatFirst[2,1,j]*x[2,1,2,j]-(RateSuppFirstResis[1]+RateStageTreatFirstR[1]+mu[3,1,2])*x[3,1,2,j]+RateStageTreatFirstL[1]*x[3,2,2,j]
      dx[4,1,2,j]=RateSuppFirstResis[1]*x[3,1,2,j]-(RateFailFirstResis[1]+mu[4,1,2])*x[4,1,2,j]+RateStageSuppFirst[1]*x[4,2,2,j]
      dx[5,1,2,j]=RateFailFirstResis[1]*x[4,1,2,j]-(RateTreatSecond[1,j]+RateStageFailFirst[1]+mu[5,1,2 ])*x[5,1,2,j]
      dx[6,1,2,j]=RateTreatSecond[1,j]*x[5,1,2,j]-(RateSuppSecond[1]+RateStageTreatSecondR[1]+mu[6,1,2 ])*x[6,1,2,j]+RateStageTreatSecondL[1]*x[6,2,2,j]
      dx[7,1,2,j]=RateSuppSecond[1]*x[6,1,2,j]-(RateFailSecond[1]+mu[7,1,2])*x[7,1,2,j]+RateStageSuppSecond[1]*x[7,2,2,j]
      dx[8,1,2,j]=RateFailSecond[1]*x[7,1,2,j]-(RateStageFailSecond[1]+mu[8,1,2])*x[8,1,2,j]
      
      dx[1,2,2,j]=-(RateStageInf[2]+RateDiag[2,j]+mu[1,2,2])*x[1,2,2,j]+RateStageInf[1]*x[1,1,2,j]
      dx[2,2,2,j]=dx[2,2,2,j]+RateDiag[2,j]*x[1,2,2,j]-(RateTreatFirst[2,2,j]+RateStageDiag[2]+mu[2,2,2])*x[2,2,2,j]+RateStageDiag[1]*x[2,1,2,j]
      dx[3,2,2,j]=dx[3,2,2,j]+RateTreatFirst[2,2,j]*x[2,2,2,j]-(RateSuppFirstResis[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,2 ])*x[3,2,2,j]+RateStageTreatFirstR[1]*x[3,1,2,j]+RateStageTreatFirstL[2]*x[3,3,2,j]
      dx[4,2,2,j]=RateSuppFirstResis[2]*x[3,2,2,j]-(RateFailFirstResis[2]+RateStageSuppFirst[2]+mu[4,2,2])*x[4,2,2,j]+RateStageSuppFirst[2]*x[4,3,2,j]
      dx[5,2,2,j]=RateFailFirstResis[2]*x[4,2,2,j]-(RateTreatSecond[2,j]+RateStageFailFirst[2]+mu[5,2,2])*x[5,2,2,j]+RateStageFailFirst[1]*x[5,1,2,j]
      dx[6,2,2,j]=RateTreatSecond[2,j]*x[5,2,2,j]-(RateSuppSecond[2]+RateStageTreatSecondR[2]+RateStageTreatSecondL[1]+mu[6,2,2])*x[6,2,2,j]+RateStageTreatSecondR[1]*x[6,1,2,j]+RateStageTreatSecondL[2]*x[6,3,2,j]
      dx[7,2,2,j]=RateSuppSecond[2]*x[6,2,2,j]-(RateFailSecond[2]+RateStageSuppSecond[1]+mu[7,2,2])*x[7,2,2,j]+RateStageSuppSecond[2]*x[7,3,2,j]
      dx[8,2,2,j]=RateFailSecond[2]*x[7,2,2,j]-(RateStageFailSecond[2]+mu[8,2,2])*x[8,2,2,j]+RateStageFailSecond[1]*x[8,1,2,j]
      
      dx[1,3,2,j]=-(RateStageInf[3]+RateDiag[3,j]+mu[1,3,2])*x[1,3,2,j]+RateStageInf[2]*x[1,2,2,j]
      dx[2,3,2,j]=dx[2,3,2,j]+RateDiag[3,j]*x[1,3,2,j]-(RateTreatFirst[2,3,j]+RateStageDiag[3]+mu[2,3,2])*x[2,3,2,j]+RateStageDiag[2]*x[2,2,2,j]
      dx[3,3,2,j]=dx[3,3,2,j]+RateTreatFirst[2,3,j]*x[2,3,2,j]-(RateSuppFirstResis[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,2])*x[3,3,2,j]+RateStageTreatFirstR[2]*x[3,2,2,j]+RateStageTreatFirstL[3]*x[3,4,2,j]
      dx[4,3,2,j]=RateSuppFirstResis[3]*x[3,3,2,j]-(RateFailFirstResis[3]+RateStageSuppFirst[2]+mu[4,3,2])*x[4,3,2,j]+RateStageSuppFirst[3]*x[4,4,2,j]
      dx[5,3,2,j]=RateFailFirstResis[3]*x[4,3,2,j]-(RateTreatSecond[3,j]+RateStageFailFirst[3]+mu[5,3,2])*x[5,3,2,j]+RateStageFailFirst[2]*x[5,2,2,j]
      dx[6,3,2,j]=RateTreatSecond[3,j]*x[5,3,2,j]-(RateSuppSecond[3]+RateStageTreatSecondR[3]+RateStageTreatSecondL[2]+mu[6,3,2])*x[6,3,2,j]+RateStageTreatSecondR[2]*x[6,2,2,j]+RateStageTreatSecondL[3]*x[6,4,2,j]
      dx[7,3,2,j]=RateSuppSecond[3]*x[6,3,2,j]-(RateFailSecond[3]+RateStageSuppSecond[2]+mu[7,3,2])*x[7,3,2,j]+RateStageSuppSecond[3]*x[7,4,2,j]
      dx[8,3,2,j]=RateFailSecond[3]*x[7,3,2,j]-(RateStageFailSecond[3]+mu[8,3,2])*x[8,3,2,j]+RateStageFailSecond[2]*x[8,2,2,j]
      
      dx[1,4,2,j]=-(RateDiag[4,j]+mu[1,4,2])*x[1,4,2,j]+RateStageInf[3]*x[1,3,2,j]
      dx[2,4,2,j]=dx[2,4,2,j]+RateDiag[4,j]*x[1,4,2,j]-(RateTreatFirst[2,4,j]+mu[2,4,2])*x[2,4,2,j]+RateStageDiag[3]*x[2,3,2,j]
      dx[3,4,2,j]=dx[3,4,2,j]+RateTreatFirst[2,4,j]*x[2,4,2,j]-(RateSuppFirstResis[4]+RateStageTreatFirstL[3]+mu[3,4,2])*x[3,4,2,j]+RateStageTreatFirstR[3]*x[3,3,2,j]
      dx[4,4,2,j]=RateSuppFirstResis[4]*x[3,4,2,j]-(RateFailFirstResis[4]+RateStageSuppFirst[3]+mu[4,4,2])*x[4,4,2,j]
      dx[5,4,2,j]=RateFailFirstResis[4]*x[4,4,2,j]-(RateTreatSecond[4,j]+mu[5,4,2])*x[5,4,2,j]+RateStageFailFirst[3]*x[5,3,2,j]
      dx[6,4,2,j]=RateTreatSecond[4,j]*x[5,4,2,j]-(RateSuppSecond[4]+RateStageTreatSecondL[3]+mu[6,4,2])*x[6,4,2,j]+RateStageTreatSecondR[3]*x[6,3,2,j]
      dx[7,4,2,j]=RateSuppSecond[4]*x[6,4,2,j]-(RateFailSecond[4]+RateStageSuppSecond[3]+mu[7,4,2])*x[7,4,2,j]
      dx[8,4,2,j]=RateFailSecond[4]*x[7,4,2,j]-mu[8,4,2]*x[8,4,2,j]+RateStageFailSecond[3]*x[8,3,2,j]
      
      ##############################################################################################################################################################
      
      dx[1,1,2,j]=dx[1,1,2,j]
      dx[2,1,2,j]=dx[2,1,2,j]-RateDirectTreatSecond[2,1]*x[2,1,2,j]+RateStopTreatFirst[1]*x[3,1,2,j]+RateStopSuppFirst[1]*x[4,1,2,j]+RateStopFailFirst[1]*x[5,1,2,j]+
        RateStopTreatSecond[1]*x[6,1,2,j]+RateStopSuppSecond[1]*x[7,1,2,j]+RateStopFailSecond[1]*x[8,1,2,j]
      dx[3,1,2,j]=dx[3,1,2,j]-RateStopTreatFirst[1]*x[3,1,2,j]-RateTreatToFailFirstResis[1]*x[3,1,2,j]
      dx[4,1,2,j]=dx[4,1,2,j]-RateStopSuppFirst[1]*x[4,1,2,j]+RateFailToSuppTreatFirstResis[1]*x[5,1,2,j]-RateSuppFirstToSecond[1]*x[4,1,2,j]
      dx[5,1,2,j]=dx[5,1,2,j]-RateStopFailFirst[1]*x[5,1,2,j]-RateFailToSuppTreatFirstResis[1]*x[5,1,2,j]+RateTreatToFailFirstResis[1]*x[3,1,2,j]
      dx[6,1,2,j]=dx[6,1,2,j]-RateStopTreatSecond[1]*x[6,1,2,j]+RateDirectTreatSecond[2,1]*x[2,1,2,j]-RateTreatToFailSecond[1]*x[6,1,2,j]
      dx[7,1,2,j]=dx[7,1,2,j]-RateStopSuppSecond[1]*x[7,1,2,j]+RateFailToSuppTreatSecond[1]*x[8,1,2,j]+RateSuppFirstToSecond[1]*x[4,1,2,j]
      dx[8,1,2,j]=dx[8,1,2,j]-RateStopFailSecond[1]*x[8,1,2,j]-RateFailToSuppTreatSecond[1]*x[8,1,2,j]+RateTreatToFailSecond[1]*x[6,1,2,j]
      
      dx[1,2,2,j]=dx[1,2,2,j]
      dx[2,2,2,j]=dx[2,2,2,j]-RateDirectTreatSecond[2,2]*x[2,2,2,j]+RateStopTreatFirst[2]*x[3,2,2,j]+RateStopSuppFirst[2]*x[4,2,2,j]+RateStopFailFirst[2]*x[5,2,2,j]+
        RateStopTreatSecond[2]*x[6,2,2,j]+RateStopSuppSecond[2]*x[7,2,2,j]+RateStopFailSecond[2]*x[8,2,2,j]
      dx[3,2,2,j]=dx[3,2,2,j]-RateStopTreatFirst[2]*x[3,2,2,j]-RateTreatToFailFirstResis[2]*x[3,2,2,j]
      dx[4,2,2,j]=dx[4,2,2,j]-RateStopSuppFirst[2]*x[4,2,2,j]+RateFailToSuppTreatFirstResis[2]*x[5,2,2,j]-RateSuppFirstToSecond[2]*x[4,2,2,j]
      dx[5,2,2,j]=dx[5,2,2,j]-RateStopFailFirst[2]*x[5,2,2,j]-RateFailToSuppTreatFirstResis[2]*x[5,2,2,j]+RateTreatToFailFirstResis[2]*x[3,2,2,j]
      dx[6,2,2,j]=dx[6,2,2,j]-RateStopTreatSecond[2]*x[6,2,2,j]+RateDirectTreatSecond[2,2]*x[2,2,2,j]-RateTreatToFailSecond[2]*x[6,2,2,j]
      dx[7,2,2,j]=dx[7,2,2,j]-RateStopSuppSecond[2]*x[7,2,2,j]+RateFailToSuppTreatSecond[2]*x[8,2,2,j]+RateSuppFirstToSecond[2]*x[4,2,2,j]
      dx[8,2,2,j]=dx[8,2,2,j]-RateStopFailSecond[2]*x[8,2,2,j]-RateFailToSuppTreatSecond[2]*x[8,2,2,j]+RateTreatToFailSecond[2]*x[6,2,2,j]
      
      dx[1,3,2,j]=dx[1,3,2,j]
      dx[2,3,2,j]=dx[2,3,2,j]-RateDirectTreatSecond[2,3]*x[2,3,2,j]+RateStopTreatFirst[3]*x[3,3,2,j]+RateStopSuppFirst[3]*x[4,3,2,j]+RateStopFailFirst[3]*x[5,3,2,j]+
        RateStopTreatSecond[3]*x[6,3,2,j]+RateStopSuppSecond[3]*x[7,3,2,j]+RateStopFailSecond[3]*x[8,3,2,j]
      dx[3,3,2,j]=dx[3,3,2,j]-RateStopTreatFirst[3]*x[3,3,2,j]-RateTreatToFailFirstResis[3]*x[3,3,2,j]
      dx[4,3,2,j]=dx[4,3,2,j]-RateStopSuppFirst[3]*x[4,3,2,j]+RateFailToSuppTreatFirstResis[3]*x[5,3,2,j]-RateSuppFirstToSecond[3]*x[4,3,2,j]
      dx[5,3,2,j]=dx[5,3,2,j]-RateStopFailFirst[3]*x[5,3,2,j]-RateFailToSuppTreatFirstResis[3]*x[5,3,2,j]+RateTreatToFailFirstResis[3]*x[3,3,2,j]
      dx[6,3,2,j]=dx[6,3,2,j]-RateStopTreatSecond[3]*x[6,3,2,j]+RateDirectTreatSecond[2,3]*x[2,3,2,j]-RateTreatToFailSecond[3]*x[6,3,2,j]
      dx[7,3,2,j]=dx[7,3,2,j]-RateStopSuppSecond[3]*x[7,3,2,j]+RateFailToSuppTreatSecond[3]*x[8,3,2,j]+RateSuppFirstToSecond[3]*x[4,3,2,j]
      dx[8,3,2,j]=dx[8,3,2,j]-RateStopFailSecond[3]*x[8,3,2,j]-RateFailToSuppTreatSecond[3]*x[8,3,2,j]+RateTreatToFailSecond[3]*x[6,3,2,j]
      
      dx[1,4,2,j]=dx[1,4,2,j]
      dx[2,4,2,j]=dx[2,4,2,j]-RateDirectTreatSecond[2,4]*x[2,4,2,j]+RateStopTreatFirst[4]*x[3,4,2,j]+RateStopSuppFirst[4]*x[4,4,2,j]+RateStopFailFirst[4]*x[5,4,2,j]+
        RateStopTreatSecond[4]*x[6,4,2,j]+RateStopSuppSecond[4]*x[7,4,2,j]+RateStopFailSecond[4]*x[8,4,2,j]
      dx[3,4,2,j]=dx[3,4,2,j]-RateStopTreatFirst[4]*x[3,4,2,j]-RateTreatToFailFirstResis[4]*x[3,4,2,j]
      dx[4,4,2,j]=dx[4,4,2,j]-RateStopSuppFirst[4]*x[4,4,2,j]+RateFailToSuppTreatFirstResis[4]*x[5,4,2,j]-RateSuppFirstToSecond[4]*x[4,4,2,j]
      dx[5,4,2,j]=dx[5,4,2,j]-RateStopFailFirst[4]*x[5,4,2,j]-RateFailToSuppTreatFirstResis[4]*x[5,4,2,j]+RateTreatToFailFirstResis[4]*x[3,4,2,j]
      dx[6,4,2,j]=dx[6,4,2,j]-RateStopTreatSecond[4]*x[6,4,2,j]+RateDirectTreatSecond[2,4]*x[2,4,2,j]-RateTreatToFailSecond[4]*x[6,4,2,j]
      dx[7,4,2,j]=dx[7,4,2,j]-RateStopSuppSecond[4]*x[7,4,2,j]+RateFailToSuppTreatSecond[4]*x[8,4,2,j]+RateSuppFirstToSecond[4]*x[4,4,2,j]
      dx[8,4,2,j]=dx[8,4,2,j]-RateStopFailSecond[4]*x[8,4,2,j]-RateFailToSuppTreatSecond[4]*x[8,4,2,j]+RateTreatToFailSecond[4]*x[6,4,2,j]
      
      #############################################################################################################################################################################
      #############################################################################################################################################################################
      #dx with interaction between resistant and susceptible compartments
      dx[1,1,1,j]=dx[1,1,1,j]+RateSusceptible*x[1,1,2,j]
      dx[2,1,1,j]=dx[2,1,1,j]+RateSusceptible*x[2,1,2,j]
      dx[3,1,1,j]=dx[3,1,1,j]
      dx[4,1,1,j]=dx[4,1,1,j]
      dx[5,1,1,j]=dx[5,1,1,j]-RateResistant*x[5,1,1,j]
      dx[6,1,1,j]=dx[6,1,1,j]
      dx[7,1,1,j]=dx[7,1,1,j]
      dx[8,1,1,j]=dx[8,1,1,j]
      
      dx[1,2,1,j]=dx[1,2,1,j]+RateSusceptible*x[1,2,2,j]
      dx[2,2,1,j]=dx[2,2,1,j]+RateSusceptible*x[2,2,2,j]
      dx[3,2,1,j]=dx[3,2,1,j]
      dx[4,2,1,j]=dx[4,2,1,j]
      dx[5,2,1,j]=dx[5,2,1,j]-RateResistant*x[5,2,1,j]
      dx[6,2,1,j]=dx[6,2,1,j]
      dx[7,2,1,j]=dx[7,2,1,j]
      dx[8,2,1,j]=dx[8,2,1,j]
      
      dx[1,3,1,j]=dx[1,3,1,j]+RateSusceptible*x[1,3,2,j]
      dx[2,3,1,j]=dx[2,3,1,j]+RateSusceptible*x[2,3,2,j]
      dx[3,3,1,j]=dx[3,3,1,j]
      dx[4,3,1,j]=dx[4,3,1,j]
      dx[5,3,1,j]=dx[5,3,1,j]-RateResistant*x[5,3,1,j]
      dx[6,3,1,j]=dx[6,3,1,j]
      dx[7,3,1,j]=dx[7,3,1,j]
      dx[8,3,1,j]=dx[8,3,1,j]
      
      dx[1,4,1,j]=dx[1,4,1,j]+RateSusceptible*x[1,4,2,j]
      dx[2,4,1,j]=dx[2,4,1,j]+RateSusceptible*x[2,4,2,j]
      dx[3,4,1,j]=dx[3,4,1,j]
      dx[4,4,1,j]=dx[4,4,1,j]
      dx[5,4,1,j]=dx[5,4,1,j]-RateResistant*x[5,4,1,j]
      dx[6,4,1,j]=dx[6,4,1,j]
      dx[7,4,1,j]=dx[7,4,1,j]
      dx[8,4,1,j]=dx[8,4,1,j]
      
      #############################################################################################################################################################################
      dx[1,1,2,j]=dx[1,1,2,j]-RateSusceptible*x[1,1,2,j]
      dx[2,1,2,j]=dx[2,1,2,j]-RateSusceptible*x[2,1,2,j]
      dx[3,1,2,j]=dx[3,1,2,j]
      dx[4,1,2,j]=dx[4,1,2,j]
      dx[5,1,2,j]=dx[5,1,2,j]+RateResistant*x[5,1,1,j]
      dx[6,1,2,j]=dx[6,1,2,j]
      dx[7,1,2,j]=dx[7,1,2,j]
      dx[8,1,2,j]=dx[8,1,2,j]
      
      dx[1,2,2,j]=dx[1,2,2,j]-RateSusceptible*x[1,2,2,j]
      dx[2,2,2,j]=dx[2,2,2,j]-RateSusceptible*x[2,2,2,j]
      dx[3,2,2,j]=dx[3,2,2,j]
      dx[4,2,2,j]=dx[4,2,2,j]
      dx[5,2,2,j]=dx[5,2,2,j]+RateResistant*x[5,2,1,j]
      dx[6,2,2,j]=dx[6,2,2,j]
      dx[7,2,2,j]=dx[7,2,2,j]
      dx[8,2,2,j]=dx[8,2,2,j]
      
      dx[1,3,2,j]=dx[1,3,2,j]-RateSusceptible*x[1,3,2,j]
      dx[2,3,2,j]=dx[2,3,2,j]-RateSusceptible*x[2,3,2,j]
      dx[3,3,2,j]=dx[3,3,2,j]
      dx[4,3,2,j]=dx[4,3,2,j]
      dx[5,3,2,j]=dx[5,3,2,j]+RateResistant*x[5,3,1,j]
      dx[6,3,2,j]=dx[6,3,2,j]
      dx[7,3,2,j]=dx[7,3,2,j]
      dx[8,3,2,j]=dx[8,3,2,j]
      
      dx[1,4,2,j]=dx[1,4,2,j]-RateSusceptible*x[1,4,2,j]
      dx[2,4,2,j]=dx[2,4,2,j]-RateSusceptible*x[2,4,2,j]
      dx[3,4,2,j]=dx[3,4,2,j]
      dx[4,4,2,j]=dx[4,4,2,j]
      dx[5,4,2,j]=dx[5,4,2,j]+RateResistant*x[5,4,1,j]
      dx[6,4,2,j]=dx[6,4,2,j]
      dx[7,4,2,j]=dx[7,4,2,j]
      dx[8,4,2,j]=dx[8,4,2,j]
      
      
      #RateTreatFirstDTG[res,cd4,gender]------------------------------
      #NNRTI susceptible
      dx[2,1,1,j]=dx[2,1,1,j]-RateTreatFirstDTG[1,1,j]*x[2,1,1,j]
      dx[2,2,1,j]=dx[2,2,1,j]-RateTreatFirstDTG[1,2,j]*x[2,2,1,j]
      dx[2,3,1,j]=dx[2,3,1,j]-RateTreatFirstDTG[1,3,j]*x[2,3,1,j]
      dx[2,4,1,j]=dx[2,4,1,j]-RateTreatFirstDTG[1,4,j]*x[2,4,1,j]
      
      dx[2,1,2,j]=dx[2,1,2,j]-RateTreatFirstDTG[2,1,j]*x[2,1,2,j]
      dx[2,2,2,j]=dx[2,2,2,j]-RateTreatFirstDTG[2,2,j]*x[2,2,2,j]
      dx[2,3,2,j]=dx[2,3,2,j]-RateTreatFirstDTG[2,3,j]*x[2,3,2,j]
      dx[2,4,2,j]=dx[2,4,2,j]-RateTreatFirstDTG[2,4,j]*x[2,4,2,j]
      
      #RateSwitchDTG[care,cd4,res,gender]----------------------------
      #NNRTI susceptible
      dx[3,1,1,j]=dx[3,1,1,j]-RateSwitchDTG[1,1,1,j]*x[3,1,1,j]
      dx[3,2,1,j]=dx[3,2,1,j]-RateSwitchDTG[1,2,1,j]*x[3,2,1,j]
      dx[3,3,1,j]=dx[3,3,1,j]-RateSwitchDTG[1,3,1,j]*x[3,3,1,j]
      dx[3,4,1,j]=dx[3,4,1,j]-RateSwitchDTG[1,4,1,j]*x[3,4,1,j]
      
      dx[4,1,1,j]=dx[4,1,1,j]-RateSwitchDTG[2,1,1,j]*x[4,1,1,j]
      dx[4,2,1,j]=dx[4,2,1,j]-RateSwitchDTG[2,2,1,j]*x[4,2,1,j]
      dx[4,3,1,j]=dx[4,3,1,j]-RateSwitchDTG[2,3,1,j]*x[4,3,1,j]
      dx[4,4,1,j]=dx[4,4,1,j]-RateSwitchDTG[2,4,1,j]*x[4,4,1,j]
      
      dx[5,1,1,j]=dx[5,1,1,j]-RateSwitchDTG[3,1,1,j]*x[5,1,1,j]
      dx[5,2,1,j]=dx[5,2,1,j]-RateSwitchDTG[3,2,1,j]*x[5,2,1,j]
      dx[5,3,1,j]=dx[5,3,1,j]-RateSwitchDTG[3,3,1,j]*x[5,3,1,j]
      dx[5,4,1,j]=dx[5,4,1,j]-RateSwitchDTG[3,4,1,j]*x[5,4,1,j]
      
      #NNRTI resistant
      dx[3,1,2,j]=dx[3,1,2,j]-RateSwitchDTG[1,1,2,j]*x[3,1,2,j]
      dx[3,2,2,j]=dx[3,2,2,j]-RateSwitchDTG[1,2,2,j]*x[3,2,2,j]
      dx[3,3,2,j]=dx[3,3,2,j]-RateSwitchDTG[1,3,2,j]*x[3,3,2,j]
      dx[3,4,2,j]=dx[3,4,2,j]-RateSwitchDTG[1,4,2,j]*x[3,4,2,j]
      
      dx[4,1,2,j]=dx[4,1,2,j]-RateSwitchDTG[2,1,2,j]*x[4,1,2,j]
      dx[4,2,2,j]=dx[4,2,2,j]-RateSwitchDTG[2,2,2,j]*x[4,2,2,j]
      dx[4,3,2,j]=dx[4,3,2,j]-RateSwitchDTG[2,3,2,j]*x[4,3,2,j]
      dx[4,4,2,j]=dx[4,4,2,j]-RateSwitchDTG[2,4,2,j]*x[4,4,2,j]
      
      dx[5,1,2,j]=dx[5,1,2,j]-RateSwitchDTG[3,1,2,j]*x[5,1,2,j]
      dx[5,2,2,j]=dx[5,2,2,j]-RateSwitchDTG[3,2,2,j]*x[5,2,2,j]
      dx[5,3,2,j]=dx[5,3,2,j]-RateSwitchDTG[3,3,2,j]*x[5,3,2,j]
      dx[5,4,2,j]=dx[5,4,2,j]-RateSwitchDTG[3,4,2,j]*x[5,4,2,j]
      
      #DTG compartments-------------------------------------
      #NNRTI susceptible
      dx[9,1,1,j]=dx[9,1,1,j]+RateTreatFirstDTG[1,1,j]*x[2,1,1,j]+RateSwitchDTG[3,1,1,j]*x[5,1,1,j]+RateSwitchDTG[1,1,1,j]*x[3,1,1,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,1])*x[9,1,1,j]+RateStageTreatFirstL[1]*x[9,2,1,j]
      dx[10,1,1,j]=RateSuppFirstSusc[1]*x[9,1,1,j]-(RateFailFirstSusc[1]+mu[4,1,1])*x[10,1,1,j]+RateStageSuppFirst[1]*x[10,2,1,j]+RateSwitchDTG[2,1,1,j]*x[4,1,1,j]
      dx[11,1,1,j]=RateFailFirstSusc[1]*x[10,1,1,j]-(RateStageFailFirst[1]+mu[5,1,1])*x[11,1,1,j]
      
      dx[9,2,1,j]=dx[9,2,1,j]+RateTreatFirstDTG[1,2,j]*x[2,2,1,j]+RateSwitchDTG[3,2,1,j]*x[5,2,1,j]+RateSwitchDTG[1,2,1,j]*x[3,2,1,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,1])*x[9,2,1,j]+RateStageTreatFirstR[1]*x[9,1,1,j]+RateStageTreatFirstL[2]*x[9,3,1,j]
      dx[10,2,1,j]=RateSuppFirstSusc[2]*x[9,2,1,j]-(RateFailFirstSusc[2]+RateStageSuppFirst[1]+mu[4,2,1])*x[10,2,1,j]+RateStageSuppFirst[2]*x[10,3,1,j]+RateSwitchDTG[2,2,1,j]*x[4,2,1,j]
      dx[11,2,1,j]=RateFailFirstSusc[2]*x[10,2,1,j]-(RateStageFailFirst[2]+mu[5,2,1])*x[11,2,1,j]+RateStageFailFirst[1]*x[11,1,1,j]
      
      dx[9,3,1,j]=dx[9,3,1,j]+RateTreatFirstDTG[1,3,j]*x[2,3,1,j]+RateSwitchDTG[3,3,1,j]*x[5,3,1,j]+RateSwitchDTG[1,3,1,j]*x[3,3,1,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,1])*x[9,3,1,j]+RateStageTreatFirstR[2]*x[9,2,1,j]+RateStageTreatFirstL[3]*x[9,4,1,j]
      dx[10,3,1,j]=RateSuppFirstSusc[3]*x[9,3,1,j]-(RateFailFirstSusc[3]+RateStageSuppFirst[2]+mu[4,3,1])*x[10,3,1,j]+RateStageSuppFirst[3]*x[10,4,1,j]+RateSwitchDTG[2,3,1,j]*x[4,3,1,j]
      dx[11,3,1,j]=RateFailFirstSusc[3]*x[10,3,1,j]-(RateStageFailFirst[3]+mu[5,3,1])*x[11,3,1,j]+RateStageFailFirst[2]*x[11,2,1,j]
      
      dx[9,4,1,j]=dx[9,4,1,j]+RateTreatFirstDTG[1,4,j]*x[2,4,1,j]+RateSwitchDTG[3,4,1,j]*x[5,4,1,j]+RateSwitchDTG[1,4,1,j]*x[3,4,1,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,1])*x[9,4,1,j]+RateStageTreatFirstR[3]*x[9,3,1,j]
      dx[10,4,1,j]=RateSuppFirstSusc[4]*x[9,4,1,j]-(RateFailFirstSusc[4]+RateStageSuppFirst[3]+mu[4,4,1])*x[10,4,1,j]+RateSwitchDTG[2,4,1,j]*x[4,4,1,j]
      dx[11,4,1,j]=RateFailFirstSusc[4]*x[10,4,1,j]-(mu[5,4,1])*x[11,4,1,j]+RateStageFailFirst[3]*x[11,3,1,j]
      #NNRTI resistant
      dx[9,1,2,j]=dx[9,1,2,j]+RateTreatFirstDTG[2,1,j]*x[2,1,2,j]+RateSwitchDTG[3,1,2,j]*x[5,1,2,j]+RateSwitchDTG[1,1,2,j]*x[3,1,2,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,2])*x[9,1,2,j]+RateStageTreatFirstL[1]*x[9,2,2,j]
      dx[10,1,2,j]=RateSuppFirstSusc[1]*x[9,1,2,j]-(RateFailFirstSusc[1]+mu[4,1,2])*x[10,1,2,j]+RateStageSuppFirst[1]*x[10,2,2,j]+RateSwitchDTG[2,1,2,j]*x[4,1,2,j]
      dx[11,1,2,j]=RateFailFirstSusc[1]*x[10,1,2,j]-(RateStageFailFirst[1]+mu[5,1,2])*x[11,1,2,j]
      
      dx[9,2,2,j]=dx[9,2,2,j]+RateTreatFirstDTG[2,2,j]*x[2,2,2,j]+RateSwitchDTG[3,2,2,j]*x[5,2,2,j]+RateSwitchDTG[1,2,2,j]*x[3,2,2,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,2])*x[9,2,2,j]+RateStageTreatFirstR[1]*x[9,1,2,j]+RateStageTreatFirstL[2]*x[9,3,2,j]
      dx[10,2,2,j]=RateSuppFirstSusc[2]*x[9,2,2,j]-(RateFailFirstSusc[2]+RateStageSuppFirst[1]+mu[4,2,2])*x[10,2,2,j]+RateStageSuppFirst[2]*x[10,3,2,j]+RateSwitchDTG[2,2,2,j]*x[4,2,2,j]
      dx[11,2,2,j]=RateFailFirstSusc[2]*x[10,2,2,j]-(RateStageFailFirst[2]+mu[5,2,2])*x[11,2,2,j]+RateStageFailFirst[1]*x[11,1,2,j]
      
      dx[9,3,2,j]=dx[9,3,2,j]+RateTreatFirstDTG[2,3,j]*x[2,3,2,j]+RateSwitchDTG[3,3,2,j]*x[5,3,2,j]+RateSwitchDTG[1,3,2,j]*x[3,3,2,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,2])*x[9,3,2,j]+RateStageTreatFirstR[2]*x[9,2,2,j]+RateStageTreatFirstL[3]*x[9,4,2,j]
      dx[10,3,2,j]=RateSuppFirstSusc[3]*x[9,3,2,j]-(RateFailFirstSusc[3]+RateStageSuppFirst[2]+mu[4,3,2])*x[10,3,2,j]+RateStageSuppFirst[3]*x[10,4,2,j]+RateSwitchDTG[2,3,2,j]*x[4,3,2,j]
      dx[11,3,2,j]=RateFailFirstSusc[3]*x[10,3,2,j]-(RateStageFailFirst[3]+mu[5,3,2])*x[11,3,2,j]+RateStageFailFirst[2]*x[11,2,2,j]
      
      dx[9,4,2,j]=dx[9,4,2,j]+RateTreatFirstDTG[2,4,j]*x[2,4,2,j]+RateSwitchDTG[3,4,2,j]*x[5,4,2,j]+RateSwitchDTG[1,4,2,j]*x[3,4,2,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,2])*x[9,4,2,j]+RateStageTreatFirstR[3]*x[9,3,2,j]
      dx[10,4,2,j]=RateSuppFirstSusc[4]*x[9,4,2,j]-(RateFailFirstSusc[4]+RateStageSuppFirst[3]+mu[4,4,2])*x[10,4,2,j]+RateSwitchDTG[2,4,2,j]*x[4,4,2,j]
      dx[11,4,2,j]=RateFailFirstSusc[4]*x[10,4,2,j]-(mu[5,4,2])*x[11,4,2,j]+RateStageFailFirst[3]*x[11,3,2,j]
      
      
      
      
      # #DTG compartments-------------------------------------
      # #NNRTI susceptible
      # dx[9,1,1,j]=dx[9,1,1,j]+RateTreatFirstDTG[1,1,j]*x[2,1,1,j]+RateSwitchDTG[1,1,1,j]*x[3,1,1,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,1])*x[9,1,1,j]+RateStageTreatFirstL[1]*x[9,2,1,j]
      # dx[10,1,1,j]=RateSwitchDTG[2,1,1,j]*x[4,1,1,j]+RateSuppFirstSusc[1]*x[9,1,1,j]-(RateFailFirstSusc[1]+mu[4,1,1])*x[10,1,1,j]+RateStageSuppFirst[1]*x[10,2,1,j]
      # dx[11,1,1,j]=RateSwitchDTG[3,1,1,j]*x[5,1,1,j]+RateFailFirstSusc[1]*x[10,1,1,j]-(RateStageFailFirst[1]+mu[5,1,1])*x[11,1,1,j]
      # 
      # dx[9,2,1,j]=dx[9,2,1,j]+RateTreatFirstDTG[1,2,j]*x[2,2,1,j]+RateSwitchDTG[1,2,1,j]*x[3,2,1,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,1])*x[9,2,1,j]+RateStageTreatFirstR[1]*x[9,1,1,j]+RateStageTreatFirstL[2]*x[9,3,1,j]
      # dx[10,2,1,j]=RateSwitchDTG[2,2,1,j]*x[4,2,1,j]+RateSuppFirstSusc[2]*x[9,2,1,j]-(RateFailFirstSusc[2]+RateStageSuppFirst[1]+mu[4,2,1])*x[10,2,1,j]+RateStageSuppFirst[2]*x[10,3,1,j]
      # dx[11,2,1,j]=RateSwitchDTG[3,2,1,j]*x[5,2,1,j]+RateFailFirstSusc[2]*x[10,2,1,j]-(RateStageFailFirst[2]+mu[5,2,1])*x[11,2,1,j]+RateStageFailFirst[1]*x[11,1,1,j]
      # 
      # dx[9,3,1,j]=dx[9,3,1,j]+RateTreatFirstDTG[1,3,j]*x[2,3,1,j]+RateSwitchDTG[1,3,1,j]*x[3,3,1,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,1])*x[9,3,1,j]+RateStageTreatFirstR[2]*x[9,2,1,j]+RateStageTreatFirstL[3]*x[9,4,1,j]
      # dx[10,3,1,j]=RateSwitchDTG[2,3,1,j]*x[4,3,1,j]+RateSuppFirstSusc[3]*x[9,3,1,j]-(RateFailFirstSusc[3]+RateStageSuppFirst[2]+mu[4,3,1])*x[10,3,1,j]+RateStageSuppFirst[3]*x[10,4,1,j]
      # dx[11,3,1,j]=RateSwitchDTG[3,3,1,j]*x[5,3,1,j]+RateFailFirstSusc[3]*x[10,3,1,j]-(RateStageFailFirst[3]+mu[5,3,1])*x[11,3,1,j]+RateStageFailFirst[2]*x[11,2,1,j]
      # 
      # dx[9,4,1,j]=dx[9,4,1,j]+RateTreatFirstDTG[1,4,j]*x[2,4,1,j]+RateSwitchDTG[1,4,1,j]*x[3,4,1,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,1])*x[9,4,1,j]+RateStageTreatFirstR[3]*x[9,3,1,j]
      # dx[10,4,1,j]=RateSwitchDTG[2,4,1,j]*x[4,4,1,j]+RateSuppFirstSusc[4]*x[9,4,1,j]-(RateFailFirstSusc[4]+RateStageSuppFirst[3]+mu[4,4,1])*x[10,4,1,j]
      # dx[11,4,1,j]=RateSwitchDTG[3,4,1,j]*x[5,4,1,j]+RateFailFirstSusc[4]*x[10,4,1,j]-(mu[5,4,1])*x[11,4,1,j]+RateStageFailFirst[3]*x[11,3,1,j]
      # ##############################################################################################################################################################
      
      dx[9,1,1,j]=dx[9,1,1,j]-RateStopTreatFirst[1]*x[9,1,1,j]-RateTreatToFailFirstSusc[1]*x[9,1,1,j]
      dx[10,1,1,j]=dx[10,1,1,j]-RateStopSuppFirst[1]*x[10,1,1,j]+RateFailToSuppTreatFirstSusc[1]*x[11,1,1,j]
      dx[11,1,1,j]=dx[11,1,1,j]-RateStopFailFirst[1]*x[11,1,1,j]-RateFailToSuppTreatFirstSusc[1]*x[11,1,1,j]+RateTreatToFailFirstSusc[1]*x[9,1,1,j]
      
      dx[9,2,1,j]=dx[9,2,1,j]-RateStopTreatFirst[2]*x[9,2,1,j]-RateTreatToFailFirstSusc[2]*x[9,2,1,j]
      dx[10,2,1,j]=dx[10,2,1,j]-RateStopSuppFirst[2]*x[10,2,1,j]+RateFailToSuppTreatFirstSusc[2]*x[11,2,1,j]-RateSuppFirstToSecond[2]*x[10,2,1,j]
      dx[11,2,1,j]=dx[11,2,1,j]-RateStopFailFirst[2]*x[11,2,1,j]-RateFailToSuppTreatFirstSusc[2]*x[11,2,1,j]+RateTreatToFailFirstSusc[2]*x[9,2,1,j]
      
      dx[9,3,1,j]=dx[9,3,1,j]-RateStopTreatFirst[3]*x[9,3,1,j]-RateTreatToFailFirstSusc[3]*x[9,3,1,j]
      dx[10,3,1,j]=dx[10,3,1,j]-RateStopSuppFirst[3]*x[10,3,1,j]+RateFailToSuppTreatFirstSusc[3]*x[11,3,1,j]
      dx[11,3,1,j]=dx[11,3,1,j]-RateStopFailFirst[3]*x[11,3,1,j]-RateFailToSuppTreatFirstSusc[3]*x[11,3,1,j]+RateTreatToFailFirstSusc[3]*x[9,3,1,j]
      
      dx[9,4,1,j]=dx[9,4,1,j]-RateStopTreatFirst[4]*x[9,4,1,j]-RateTreatToFailFirstSusc[4]*x[9,4,1,j]
      dx[10,4,1,j]=dx[10,4,1,j]-RateStopSuppFirst[4]*x[10,4,1,j]+RateFailToSuppTreatFirstSusc[4]*x[11,4,1,j]
      dx[11,4,1,j]=dx[11,4,1,j]-RateStopFailFirst[4]*x[11,4,1,j]-RateFailToSuppTreatFirstSusc[4]*x[11,4,1,j]+RateTreatToFailFirstSusc[4]*x[9,4,1,j]
      #############################################################################################################################################################################
      #############################################################################################################################################################################
      #NNRTI resistant
      # dx[9,1,2,j]=dx[9,1,2,j]+RateTreatFirstDTG[2,1,j]*x[2,1,2,j]+RateSwitchDTG[1,1,2,j]*x[3,1,2,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,2])*x[9,1,2,j]+RateStageTreatFirstL[1]*x[9,2,2,j]
      # dx[10,1,2,j]=RateSwitchDTG[2,1,2,j]*x[4,1,2,j]+RateSuppFirstSusc[1]*x[9,1,2,j]-(RateFailFirstSusc[1]+mu[4,1,2])*x[10,1,2,j]+RateStageSuppFirst[1]*x[10,2,2,j]
      # dx[11,1,2,j]=RateSwitchDTG[3,1,2,j]*x[5,1,2,j]+RateFailFirstSusc[1]*x[10,1,2,j]-(RateStageFailFirst[1]+mu[5,1,2])*x[11,1,2,j]
      # 
      # dx[9,2,2,j]=dx[9,2,2,j]+RateTreatFirstDTG[2,2,j]*x[2,2,2,j]+RateSwitchDTG[1,2,2,j]*x[3,2,2,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,2])*x[9,2,2,j]+RateStageTreatFirstR[1]*x[9,1,2,j]+RateStageTreatFirstL[2]*x[9,3,2,j]
      # dx[10,2,2,j]=RateSwitchDTG[2,2,2,j]*x[4,2,2,j]+RateSuppFirstSusc[2]*x[9,2,2,j]-(RateFailFirstSusc[2]+RateStageSuppFirst[1]+mu[4,2,2])*x[10,2,2,j]+RateStageSuppFirst[2]*x[10,3,2,j]
      # dx[11,2,2,j]=RateSwitchDTG[3,2,2,j]*x[5,2,2,j]+RateFailFirstSusc[2]*x[10,2,2,j]-(RateStageFailFirst[2]+mu[5,2,2])*x[11,2,2,j]+RateStageFailFirst[1]*x[11,1,2,j]
      # 
      # dx[9,3,2,j]=dx[9,3,2,j]+RateTreatFirstDTG[2,3,j]*x[2,3,2,j]+RateSwitchDTG[1,3,2,j]*x[3,3,2,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,2])*x[9,3,2,j]+RateStageTreatFirstR[2]*x[9,2,2,j]+RateStageTreatFirstL[3]*x[9,4,2,j]
      # dx[10,3,2,j]=RateSwitchDTG[2,3,2,j]*x[4,3,2,j]+RateSuppFirstSusc[3]*x[9,3,2,j]-(RateFailFirstSusc[3]+RateStageSuppFirst[2]+mu[4,3,2])*x[10,3,2,j]+RateStageSuppFirst[3]*x[10,4,2,j]
      # dx[11,3,2,j]=RateSwitchDTG[3,3,2,j]*x[5,3,2,j]+RateFailFirstSusc[3]*x[10,3,2,j]-(RateStageFailFirst[3]+mu[5,3,2])*x[11,3,2,j]+RateStageFailFirst[2]*x[11,2,2,j]
      # 
      # dx[9,4,2,j]=dx[9,4,2,j]+RateTreatFirstDTG[2,4,j]*x[2,4,2,j]+RateSwitchDTG[1,4,2,j]*x[3,4,2,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,2])*x[9,4,2,j]+RateStageTreatFirstR[3]*x[9,3,2,j]
      # dx[10,4,2,j]=RateSwitchDTG[2,4,2,j]*x[4,4,2,j]+RateSuppFirstSusc[4]*x[9,4,2,j]-(RateFailFirstSusc[4]+RateStageSuppFirst[3]+mu[4,4,2])*x[10,4,2,j]
      # dx[11,4,2,j]=RateSwitchDTG[3,4,2,j]*x[5,4,2,j]+RateFailFirstSusc[4]*x[10,4,2,j]-(mu[5,4,2])*x[11,4,2,j]+RateStageFailFirst[3]*x[11,3,2,j]
      ##############################################################################################################################################################
      
      dx[9,1,2,j]=dx[9,1,2,j]-RateStopTreatFirst[1]*x[9,1,2,j]-RateTreatToFailFirstSusc[1]*x[9,1,2,j]
      dx[10,1,2,j]=dx[10,1,2,j]-RateStopSuppFirst[1]*x[10,1,2,j]+RateFailToSuppTreatFirstSusc[1]*x[11,1,2,j]
      dx[11,1,2,j]=dx[11,1,2,j]-RateStopFailFirst[1]*x[11,1,2,j]-RateFailToSuppTreatFirstSusc[1]*x[11,1,2,j]+RateTreatToFailFirstSusc[1]*x[9,1,2,j]
      
      dx[9,2,2,j]=dx[9,2,2,j]-RateStopTreatFirst[2]*x[9,2,2,j]-RateTreatToFailFirstSusc[2]*x[9,2,2,j]
      dx[10,2,2,j]=dx[10,2,2,j]-RateStopSuppFirst[2]*x[10,2,2,j]+RateFailToSuppTreatFirstSusc[2]*x[11,2,2,j]
      dx[11,2,2,j]=dx[11,2,2,j]-RateStopFailFirst[2]*x[11,2,2,j]-RateFailToSuppTreatFirstSusc[2]*x[11,2,2,j]+RateTreatToFailFirstSusc[2]*x[9,2,2,j]
      
      dx[9,3,2,j]=dx[9,3,2,j]-RateStopTreatFirst[3]*x[9,3,2,j]-RateTreatToFailFirstSusc[3]*x[9,3,2,j]
      dx[10,3,2,j]=dx[10,3,2,j]-RateStopSuppFirst[3]*x[10,3,2,j]+RateFailToSuppTreatFirstSusc[3]*x[11,3,2,j]
      dx[11,3,2,j]=dx[11,3,2,j]-RateStopFailFirst[3]*x[11,3,2,j]-RateFailToSuppTreatFirstSusc[3]*x[11,3,2,j]+RateTreatToFailFirstSusc[3]*x[9,3,2,j]
      
      dx[9,4,2,j]=dx[9,4,2,j]-RateStopTreatFirst[4]*x[9,4,2,j]-RateTreatToFailFirstSusc[4]*x[9,4,2,j]
      dx[10,4,2,j]=dx[10,4,2,j]-RateStopSuppFirst[4]*x[10,4,2,j]+RateFailToSuppTreatFirstSusc[4]*x[11,4,2,j]
      dx[11,4,2,j]=dx[11,4,2,j]-RateStopFailFirst[4]*x[11,4,2,j]-RateFailToSuppTreatFirstSusc[4]*x[11,4,2,j]+RateTreatToFailFirstSusc[4]*x[9,4,2,j]
    }
    dx[x+dx<0]=0
    #new_infections
    new_inf=S[1]/N[1]*(sum(RateInf1[1:4,1:2,1]*apply(x[c(1),1:4,1:2,1:2],c(1,3),sum))+
                         sum(RateInf2[1:4,1:2,1]*apply(x[c(2,3,5,6,8,9,11),1:4,1:2,1:2],c(2,4),sum)))+
      S[2]/N[2]*(sum(RateInf1[1:4,1:2,2]*apply(x[c(1),1:4,1:2,1:2],c(1,3),sum))+
                   sum(RateInf2[1:4,1:2,2]*apply(x[c(2,3,5,6,8,9,11),1:4,1:2,1:2],c(2,4),sum)))
    #death
    death=sum(apply(x[1:8,1:4,1:2,1:2],c(1,2,3),sum)*mu)+sum(apply(x[9:11,1:4,1:2,1:2],c(1,2,3),sum)*mu[3:5,1:4,1:2])
    
    #new_infections resistant
    new_inf_res=S[1]/N[1]*(sum(RateInf1[1:4,1:2,1]*apply(x[c(1),1:4,2,1:2],c(1,2),sum))+
                             sum(RateInf2[1:4,1:2,1]*apply(x[c(2,3,5,6,8,9,11),1:4,2,1:2],c(2,3),sum)))+
      S[2]/N[2]*(sum(RateInf1[1:4,1:2,2]*apply(x[c(1),1:4,2,1:2],c(1,2),sum))+
                   sum(RateInf2[1:4,1:2,2]*apply(x[c(2,3,5,6,8,9,11),1:4,2,1:2],c(2,3),sum)))
    
    n1=S[1]/N[1]*sum(RateInf1[1:4,1:2,1]*apply(x[c(1),1:4,2,1:2],c(1,2),sum))#susc inf
    n2=S[1]/N[1]*sum(RateInf2[1:4,1:2,1]*apply(x[c(2,3,5,6,8,9,11),1:4,2,1:2],c(2,3),sum))#susc diag
    n3=S[2]/N[2]*sum(RateInf1[1:4,1:2,2]*apply(x[c(1),1:4,2,1:2],c(1,2),sum))#res inf
    n4=S[2]/N[2]*sum(RateInf2[1:4,1:2,2]*apply(x[c(2,3,5,6,8,9,11),1:4,2,1:2],c(2,3),sum))#res diag
    #death resistant
    death_res=sum(apply(x[1:8,1:4,2,1:2],c(1,2),sum)*mu[1:8,1:4,2])
    
    #new treated
    new_treat=sum(RateTreatFirst[1:2,1:4,1:2]*aperm(x[2,1:4,1:2,1:2],c(2,1,3)))+
      sum(RateDirectTreatSecond[1:2,1:4]*aperm(apply(x[2,1:4,1:2,1:2],c(1,2),sum),c(2,1)))
    #new treated resistant
    new_treat_res=sum(RateTreatFirst[2,1:4,1:2]*x[2,1:4,2,1:2])+
      sum(RateDirectTreatSecond[2,1:4]*apply(x[2,1:4,2,1:2],c(1),sum))
    
    # new_treat_susc=sum(RateTreatFirstDTG[1,1:4,1:2]*x[2,1:4,1,1:2])
    # new_treat_res=sum(RateTreatFirstDTG[2,1:4,1:2]*x[2,1:4,2,1:2])
    # new_diag_susc=sum(RateDiag[1:4]*apply(x[1,1:4,1,1:2],1,sum))
    # new_diag_res=sum(RateDiag[1:4]*apply(x[1,1:4,2,1:2],1,sum))
    # #new_diag_susc=sum(x[2,1:4,1,1:2])
    # #new_diag_res=sum(x[2,1:4,2,1:2])
    # mort_susc=sum(mu[2,1:4,1]*apply(x[2,1:4,1,1:2],1,sum))
    # mort_res=sum(mu[2,1:4,2]*apply(x[2,1:4,2,1:2],1,sum))
    
    #dx[3:8,1:4,1,1:2]=dx[3:8,1:4,1,1:2]+infected_m15_month[t+1]*x[3:8,1:4,1,1:2]/sum(x[3:8,1:4,1,1:2])
    ####################################################################################################################
    res=c(dx,new_inf,death,new_treat,new_inf_res,death_res,new_treat_res)
    #res=c(dx,new_inf,new_inf_res,n1,n2,n3,n4)
    #res=c(dx,new_treat_susc,new_treat_res,new_diag_susc,new_diag_res,mort_susc,mort_res)
    return(c(res))
    #return(list(c(res,newinf,newsu,newres,death)))
    #})
  })
}
xstart_f_dtg=function(par){
  k1=par[9]
  k2=par[10]
  k3=par[11]
  start_treat=107
  start_value=array(0,dim=c(11,4,2,2))
  start_value[1,1:4,1,1:2]=c(apply(undiag_value,1,sum)[1]*0.42*repart_f(k1),apply(undiag_value,1,sum)[1]*0.58*repart_f(k1))
  start_value[2,1:4,1,1:2]=c((apply(diag_value,1,sum)[1]-start_treat)*0.42*repart_f(k2),(apply(diag_value,1,sum)[1]-start_treat)*0.58*repart_f(k2))
  start_value[3,1:4,1,1:2]=c(start_treat*0.42*repart_f(k3),start_treat*0.58*repart_f(k3))
  # start_value[4,1:4,1,1:2]=start_value[3,1:4,1,1:2]*0.45
  # start_value[5,1:4,1,1:2]=start_value[3,1:4,1,1:2]*0.15
  # start_value[3,1:4,1,1:2]=start_value[3,1:4,1,1:2]*0.4
  
  start_value[4,1:4,1,1:2]=start_value[3,1:4,1,1:2]*0
  start_value[5,1:4,1,1:2]=start_value[3,1:4,1,1:2]*0
  start_value[3,1:4,1,1:2]=start_value[3,1:4,1,1:2]*1
  
  start_value[1:2,1:4,2,1:2]=array(rep(rep(c(0.01,0.01,0.01,0.01),each=2),2),dim=c(2,4,2))*start_value[1:2,1:4,1,1:2]
  start_value[1:2,1:4,1,1:2]=array(rep(rep(1-c(0.01,0.01,0.01,0.01),each=2),2),dim=c(2,4,2))*start_value[1:2,1:4,1,1:2]
  start_value[3,1:4,2,1:2]=0.01*start_value[3,1:4,1,1:2]
  start_value[3,1:4,1,1:2]=0.99*start_value[3,1:4,1,1:2]
  start_value[4,1:4,2,1:2]=0.01*start_value[4,1:4,1,1:2]
  start_value[4,1:4,1,1:2]=0.99*start_value[4,1:4,1,1:2]
  start_value[5,1:4,2,1:2]=0.7*start_value[5,1:4,1,1:2]
  start_value[5,1:4,1,1:2]=0.3*start_value[5,1:4,1,1:2]
  xstart=c(start_value,rep(0,6))
  return(xstart)
}

#15 care stages, not more used, see below
mod_cpp_dtg_2=function(t,x,p1,treat_dtg){
  #mod =function(t,x, params){
  with(as.list(params), {
    x=x[1:(240)] # 15*4*2*2=240
    ##################################################################################################################
    #Optimization
    rate1_inf <- p1[1] #number of unprotected sexual acts/month
    rate2_inf <- p1[2]
    rate1_diag <- p1[3] #diagnosis rate
    rate2_diag <- p1[4] #diagnosis rate
    rate_treat <- p1[5] #scale parameter for overall treatment rate
    rate_death <- p1[6] #scale parameter for overall death
    q <- p2[7]
    #Range
    rate_ratio <- p2[8] #varying parameter to estimate death for treated (but neither suppressed nor failed)
    
    alpha <- p2[12] #ratio between outcome when resistant/sensitive
    rate_res <- p2[13] #resistant rate
    rate_susc <- p2[14] #reversion rate
    
    #1. Infection and Diagnosis
    RateInf1=array(rep(RateTransm,each=4),dim=c(4,2,2))*rate1_inf
    RateInf2=RateInf1*rate2_inf
    
    Diag_frame[,4]=Diag_frame[,4]*rate2_diag
    RateDiag=1/rate1_diag*apply(Diag_frame,1,function(x) smooth_st(2005+t/12,x["a"],x["b"],x["c"],x["d"]))
    
    RateDiag=array(rep(RateDiag,2)*rep(c(1,diag_women),each=4),dim=c(4,2))+
      array(rep(oi_inc,2)*rep(smooth_st(2005+t/12,oi_test["a"],oi_test["b"],oi_test["c"],oi_test["d"]),8),dim=c(4,2))+
      array(c(rep(0,4),preg_inc)*rep(smooth_st(2005+t/12,preg_test["a"],preg_test["b"],preg_test["c"],preg_test["d"]),8),dim=c(4,2))
    
    #2. Treatment theoretical (if no counterfactual scenario)
    RateTreatFirst_th=apply(Treat_frame,1,function(x) smooth_st(2005+t/12,x["a"],x["b"],x["c"],x["d"]))*smooth_st(2005+t/12,2001,2012,1,200/12)*rate_treat
    
    RateTreatFirst_f=function(y){
      if(y<2005){
        return(apply(Treat_frame_original,1,function(x) smooth_st(y,x["a"],x["b"],x["c"],x["d"]))*smooth_st(y,2001,2012,1,200/12)*rate_treat)
      }else{
        return(apply(
          data.frame(t1=apply(Treat_frame,1,function(x) lin_st(y,x["a"],x["b"],x["c"],x["d"])),
                     t2=apply(Treat_frame2,1,function(x) lin_st(y,x["a"],x["b"],x["c"],x["d"]))*as.numeric(y>Treat_frame2$a[4])),
          1,max)*smooth_st(y,2001,2012,1,200/12)*rate_treat)
      }
    }
    RateTreatFirst_th=RateTreatFirst_f(2005+t/12)
    
    RateDirectTreatSecond_th=RateDirectTreatSecond
    
    #now
    RateTreatFirst=matrix(c(RateTreatFirst_th,RateTreatFirst_th),nrow=2,ncol=4,byrow=T)
    RateDirectTreatSecond=matrix(c(RateDirectTreatSecond_th,RateDirectTreatSecond_th),nrow=2,ncol=4,byrow=T)
    
    #Resistance testing at baseline
    #RateTreatFirst=matrix(c(RateTreatFirst_th,c(0,0,0,0)),nrow=2,ncol=4,byrow=T)
    #RateDirectTreatSecond=matrix(c(RateDirectTreatSecond_th,RateTreatFirst_th+RateDirectTreatSecond_th),nrow=2,ncol=4,byrow=T)
    #RateTreatFirst=ep(rate_treat,4)
    
    #2. Treatment DTG
    #percentage of women eligible for dtg
    if(treat_dtg[1]!=treat_dtg[2] & treat_dtg[2]!=0){print("percentage of women eligible for dtg not equal for FirstTreat and Switch.")}
    p_women_dtg=treat_dtg[1]
    #RateTreatFirstDTG[res,cd4,gender]
    RateTreatFirstDTG=array(rep(RateTreatFirst*treat_dtg[3],2)*rep(c(1,1),each=8),dim=c(2,4,2))*
      smooth_st(2005+t/12,2018,2018+treat_dtg["delay_first"],0,1) #delay
    #RateSwitchDTG[care,cd4,res,gender]
    #All people switch at a rate equal to RateTreatSecond
    RateSwitchDTG=array(rep(rep(treat_dtg[4]*RateTreatSecond,each=3)*rep(c(treat_dtg[5],treat_dtg[6],treat_dtg[7]),4),4),dim=c(3,4,2,2))*
      smooth_st(2005+t/12,2018,2018+treat_dtg["delay_switch"],0,1) #delay
    # RateSwitchDTG=array(rep(rep(treat_dtg[4]*rep(1/12,4),each=3)*rep(c(treat_dtg[5],treat_dtg[6],treat_dtg[7]),4),4),dim=c(3,4,2,2))*
    #    smooth_st(2005+t/12,2018,2018+treat_dtg["delay_switch"],0,1) #delay
    #Treated and suppressed people switch at a rate equal to 1year^-1, failed people at a rate still equal to RateTreatSecond
    RateSwitchDTG[1:2,1:4,1:2,1:2]=rep(treat_dtg[4],32)*rep(rep(1/12,2)*c(treat_dtg[5],treat_dtg[6]),16)*
      smooth_st(2005+t/12,2018,2018+treat_dtg["delay_switch"],0,1) #delay
    #RateTreatFirst [res,cd4,gender]
    RateTreatFirst_noDTG=array(rep(RateTreatFirst,2),dim=c(2,4,2))
    RateTreatFirst=apply(RateTreatFirst_noDTG-RateTreatFirstDTG,c(1,2,3),function(x) max(x,0))
    
    #RateTreatSecond[CD4,gender]
    RateTreatSecond_noDTG=array(rep(RateTreatSecond,2),dim=c(4,2))
    # if(2005+t/12<2018){RateTreatSecond_noDTG=array(rep(RateTreatSecond,2),dim=c(4,2))
    # }else{RateTreatSecond_noDTG=array(rep(1/12,8),dim=c(4,2))}
    RateTreatSecond=apply(RateTreatSecond_noDTG-apply(RateSwitchDTG[3,1:4,1:2,1:2],c(1,3),mean),c(1,2),function(x) max(x,0))
    #we take the mean here as RateTreatSecond is the same for men and women.
    #It has not any impact as long as RateSwitchDTG is also the same for men and women.
    
    # RateTreatFirst=matrix(rep(0,8),nrow=2,ncol=4,byrow=T)
    # RateDirectTreatSecond=matrix(rep(0,8),nrow=2,ncol=4,byrow=T)
    #3. Death
    #p1=p1_f(t,q)
    #p1=min(1,q*mean(sapply(2005+seq(t-36,t,1)/12,function(x) apply(Treat_frame[4,],1,function(y) smooth_st(x,y["a"],y["b"],y["c"],y["d"]))*smooth_st(x,2001,2012,1,200/12))))
    p1=min(1,q*mean(sapply(2005+seq(t-36,t,1)/12,function(x) RateTreatFirst_f(x)[4]/rate_treat)))
    print(p1)
    #p1=p1/2+0.2
    mu[1:2,4,1]=p1*40.9+(1-p1)*134.4
    mu[4,4,1]=p1*8.3+(1-p1)*41.7
    mu[5,4,1]=p1*11.8+(1-p1)*59.7
    mu[3,1:4,1]=rate_ratio*mu[4,1:4,1]+(1-rate_ratio)*mu[5,1:4,1]
    mu[6:8,1:4,1]=mu[3:5,1:4,1]
    mu[1:8,1:4,2]=mu[1:8,1:4,1]
    mu=mu*rate_death/1000
    print(p1*40.9+(1-p1)*134.4)
    print(mu)
    #Back and forth mutation and impact on treatment outcome
    RateResistant=1/rate_res
    RateSusceptible=1/rate_susc
    #From T to S
    RateSuppFirstSusc=RateSuppFirst/1
    RateSuppFirstResis=1/alpha*RateSuppFirst/1
    #From S to F
    RateFailFirstSusc=RateFailFirst*1
    RateFailFirstResis=alpha*RateFailFirst*1
    #From T to F
    RateTreatToFailFirstSusc=RateTreatToFailFirst*1
    RateTreatToFailFirstResis=alpha*RateTreatToFailFirst*1
    #From F to S
    RateFailToSuppTreatFirstSusc=RateFailToSuppTreatFirst/1
    RateFailToSuppTreatFirstResis=1/alpha*RateFailToSuppTreatFirst/1
    
    #For second-line regimen, no effect of NNRTI resistance on treatment as it is PI
    ##################################################################################################################
    dx=array(0,dim=c(15,4,2,2))
    x=array(x,dim=c(15,4,2,2))
    I=apply(x,4,sum)
    
    #First approach to model S: S=N-I
    N=N_value_f(t+1)
    S=N-I
    #Percentage getting NNRTI-based 1st-line[gender]
    p_NNRTI=sum(RateTreatFirst)/sum(RateTreatFirst+RateTreatFirstDTG)*c(1,treat_dtg[1])+c(0,1-treat_dtg[1])
    #infected and treated children
    I_m15=rep(infected_m15_f(t)/2,2) #number of infected children reaching 15 years old at a given month
    t_m15=rep(treat_m15_f(t)/2,2)
    
    # dx[2,1:4,1,1:2]=rep((1-res_m15)*I_m15/4,each=4)
    # dx[2,1:4,2,1:2]=rep(res_m15*I_m15/4,each=4)
    # dx[3,1:4,1,1]=rep((1-res_m15)*t_m15[1]*p_NNRTI[1]/4,each=4)
    # dx[3,1:4,2,1]=rep(res_m15*t_m15[1]*p_NNRTI[1]/4,each=4)
    # dx[13,1:4,1,2]=rep((1-res_m15)*t_m15[2]*p_NNRTI[2]/4,each=4)
    # dx[13,1:4,2,2]=rep(res_m15*t_m15[2]*p_NNRTI[2]/4,each=4)
    # dx[9,1:4,1,1:2]=rep((1-res_m15)*t_m15*(1-p_NNRTI)/4,each=4)
    # dx[9,1:4,2,1:2]=rep(res_m15*t_m15*(1-p_NNRTI)/4,each=4)
    
    # if(t>1){
    #   print("RateTreatFirst")
    #   print(RateTreatFirst)
    #   print("RateTreatSecond")
    #   print(RateTreatSecond)
    #   print("RateDirectTreatSecond")
    #   print(RateDirectTreatSecond)
    #   print("RateTreatSecond_noDTG")
    #   print(RateTreatSecond_noDTG)
    #   print("RateSwitchDTG")
    #   print(RateSwitchDTG)
    #   print("RateTreatFirstDTG")
    #   print(RateTreatFirstDTG)
    #   print("RateDiag")
    #   print(RateDiag)
    # }
    
    for(j in 1:2){
      #dx[x1,x2,x3,x4] x1:care stages x2:disease progression x3:resistance x4:gender
      #dx with only interaction with the neighbours in the compartmental model
      dx[1,1,1,j]=S[j]/N[j]*(sum(RateInf1[1:4,1:2,j]*apply(x[c(1),1:4,1,1:2],c(1,2),sum))+sum(RateInf2[1:4,1:2,j]*apply(x[c(2,3,5,6,8,9,11,12,13,15),1:4,1,1:2],c(2,3),sum)))-
        (RateStageInf[1]+RateDiag[1,j]+mu[1,1,1])*x[1,1,1,j]
      dx[2,1,1,j]=dx[2,1,1,j]+RateDiag[1,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,1,1,j]-(RateTreatFirst[1,1,j]+RateStageDiag[1]+mu[2,1,1])*x[2,1,1,j]
      dx[3,1,1,j]=dx[3,1,1,j]+RateTreatFirst[1,1,j]*x[2,1,1,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,1])*x[3,1,1,j]+RateStageTreatFirstL[1]*x[3,2,1,j]
      dx[4,1,1,j]=RateSuppFirstSusc[1]*x[3,1,1,j]-(RateFailFirstSusc[1]+mu[4,1,1])*x[4,1,1,j]+RateStageSuppFirst[1]*x[4,2,1,j]
      dx[5,1,1,j]=RateFailFirstSusc[1]*x[4,1,1,j]-(RateTreatSecond[1,j]+RateStageFailFirst[1]+mu[5,1,1])*x[5,1,1,j]
      dx[6,1,1,j]=RateTreatSecond_noDTG[1,j]*x[11,1,1,j]+RateTreatSecond_noDTG[1,j]*x[15,1,1,j]+RateTreatSecond[1,j]*x[5,1,1,j]-(RateSuppSecond[1]+RateStageTreatSecondR[1]+mu[6,1,1])*x[6,1,1,j]+RateStageTreatSecondL[1]*x[6,2,1,j]
      dx[7,1,1,j]=RateSuppSecond[1]*x[6,1,1,j]-(RateFailSecond[1]+mu[7,1,1])*x[7,1,1,j]+RateStageSuppSecond[1]*x[7,2,1,j]
      dx[8,1,1,j]=RateFailSecond[1]*x[7,1,1,j]-(RateStageFailSecond[1]+mu[8,1,1])*x[8,1,1,j]
      
      dx[1,2,1,j]=-(RateStageInf[2]+RateDiag[2,j]+mu[1,2,1])*x[1,2,1,j]+RateStageInf[1]*x[1,1,1,j]
      dx[2,2,1,j]=dx[2,2,1,j]+RateDiag[2,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,2,1,j]-(RateTreatFirst[1,2,j]+RateStageDiag[2]+mu[2,2,1])*x[2,2,1,j]+RateStageDiag[1]*x[2,1,1,j]
      dx[3,2,1,j]=dx[3,2,1,j]+RateTreatFirst[1,2,j]*x[2,2,1,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,1])*x[3,2,1,j]+RateStageTreatFirstR[1]*x[3,1,1,j]+RateStageTreatFirstL[2]*x[3,3,1,j]
      dx[4,2,1,j]=RateSuppFirstSusc[2]*x[3,2,1,j]-(RateFailFirstSusc[2]+RateStageSuppFirst[1]+mu[4,2,1])*x[4,2,1,j]+RateStageSuppFirst[2]*x[4,3,1,j]
      dx[5,2,1,j]=RateFailFirstSusc[2]*x[4,2,1,j]-(RateTreatSecond[2,j]+RateStageFailFirst[2]+mu[5,2,1])*x[5,2,1,j]+RateStageFailFirst[1]*x[5,1,1,j]
      dx[6,2,1,j]=RateTreatSecond_noDTG[2,j]*x[11,2,1,j]+RateTreatSecond_noDTG[2,j]*x[15,2,1,j]+RateTreatSecond[2,j]*x[5,2,1,j]-(RateSuppSecond[2]+RateStageTreatSecondR[2]+RateStageTreatSecondL[1]+mu[6,2,1])*x[6,2,1,j]+RateStageTreatSecondR[1]*x[6,1,1,j]+RateStageTreatSecondL[2]*x[6,3,1,j]
      dx[7,2,1,j]=RateSuppSecond[2]*x[6,2,1,j]-(RateFailSecond[2]+RateStageSuppSecond[1]+mu[7,2,1])*x[7,2,1,j]+RateStageSuppSecond[2]*x[7,3,1,j]
      dx[8,2,1,j]=RateFailSecond[2]*x[7,2,1,j]-(RateStageFailSecond[2]+mu[8,2,1])*x[8,2,1,j]+RateStageFailSecond[1]*x[8,1,1,j]
      
      dx[1,3,1,j]=-(RateStageInf[3]+RateDiag[3,j]+mu[1,3,1])*x[1,3,1,j]+RateStageInf[2]*x[1,2,1,j]
      dx[2,3,1,j]=dx[2,3,1,j]+RateDiag[3,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,3,1,j]-(RateTreatFirst[1,3,j]+RateStageDiag[3]+mu[2,3,1])*x[2,3,1,j]+RateStageDiag[2]*x[2,2,1,j]
      dx[3,3,1,j]=dx[3,3,1,j]+RateTreatFirst[1,3,j]*x[2,3,1,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,1])*x[3,3,1,j]+RateStageTreatFirstR[2]*x[3,2,1,j]+RateStageTreatFirstL[3]*x[3,4,1,j]
      dx[4,3,1,j]=RateSuppFirstSusc[3]*x[3,3,1,j]-(RateFailFirstSusc[3]+RateStageSuppFirst[2]+mu[4,3,1])*x[4,3,1,j]+RateStageSuppFirst[3]*x[4,4,1,j]
      dx[5,3,1,j]=RateFailFirstSusc[3]*x[4,3,1,j]-(RateTreatSecond[3,j]+RateStageFailFirst[3]+mu[5,3,1])*x[5,3,1,j]+RateStageFailFirst[2]*x[5,2,1,j]
      dx[6,3,1,j]=RateTreatSecond_noDTG[3,j]*x[11,3,1,j]+RateTreatSecond_noDTG[3,j]*x[15,3,1,j]+RateTreatSecond[3,j]*x[5,3,1,j]-(RateSuppSecond[3]+RateStageTreatSecondR[3]+RateStageTreatSecondL[2]+mu[6,3,1])*x[6,3,1,j]+RateStageTreatSecondR[2]*x[6,2,1,j]+RateStageTreatSecondL[3]*x[6,4,1,j]
      dx[7,3,1,j]=RateSuppSecond[3]*x[6,3,1,j]-(RateFailSecond[3]+RateStageSuppSecond[2]+mu[7,3,1])*x[7,3,1,j]+RateStageSuppSecond[3]*x[7,4,1,j]
      dx[8,3,1,j]=RateFailSecond[3]*x[7,3,1,j]-(RateStageFailSecond[3]+mu[8,3,1])*x[8,3,1,j]+RateStageFailSecond[2]*x[8,2,1,j]
      
      dx[1,4,1,j]=-(RateDiag[4,j]+mu[1,4,1])*x[1,4,1,j]+RateStageInf[3]*x[1,3,1,j]
      dx[2,4,1,j]=dx[2,4,1,j]+RateDiag[4,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,4,1,j]-(RateTreatFirst[1,4,j]+mu[2,4,1])*x[2,4,1,j]+RateStageDiag[3]*x[2,3,1,j]
      dx[3,4,1,j]=dx[3,4,1,j]+RateTreatFirst[1,4,j]*x[2,4,1,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,1])*x[3,4,1,j]+RateStageTreatFirstR[3]*x[3,3,1,j]
      dx[4,4,1,j]=RateSuppFirstSusc[4]*x[3,4,1,j]-(RateFailFirstSusc[4]+RateStageSuppFirst[3]+mu[4,4,1])*x[4,4,1,j]
      dx[5,4,1,j]=RateFailFirstSusc[4]*x[4,4,1,j]-(RateTreatSecond[4,j]+mu[5,4,1])*x[5,4,1,j]+RateStageFailFirst[3]*x[5,3,1,j]
      dx[6,4,1,j]=RateTreatSecond_noDTG[4,j]*x[11,4,1,j]+RateTreatSecond_noDTG[4,j]*x[15,4,1,j]+RateTreatSecond[4,j]*x[5,4,1,j]-(RateSuppSecond[4]+RateStageTreatSecondL[3]+mu[6,4,1])*x[6,4,1,j]+RateStageTreatSecondR[3]*x[6,3,1,j]
      dx[7,4,1,j]=RateSuppSecond[4]*x[6,4,1,j]-(RateFailSecond[4]+RateStageSuppSecond[3]+mu[7,4,1])*x[7,4,1,j]
      dx[8,4,1,j]=RateFailSecond[4]*x[7,4,1,j]-mu[8,4,1 ]*x[8,4,1,j]+RateStageFailSecond[3]*x[8,3,1,j]
      
      ##############################################################################################################################################################
      
      dx[1,1,1,j]=dx[1,1,1,j]
      dx[2,1,1,j]=dx[2,1,1,j]-RateDirectTreatSecond[1,1]*x[2,1,1,j]+RateStopTreatFirst[1]*x[3,1,1,j]+RateStopSuppFirst[1]*x[4,1,1,j]+RateStopFailFirst[1]*x[5,1,1,j]+
        (RateStopTreatSecond[1]*x[6,1,1,j]+RateStopSuppSecond[1]*x[7,1,1,j]+RateStopFailSecond[1]*x[8,1,1,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,1,1,j]=dx[3,1,1,j]-RateStopTreatFirst[1]*x[3,1,1,j]-RateTreatToFailFirstSusc[1]*x[3,1,1,j]
      dx[4,1,1,j]=dx[4,1,1,j]-RateStopSuppFirst[1]*x[4,1,1,j]+RateFailToSuppTreatFirstSusc[1]*x[5,1,1,j]-RateSuppFirstToSecond[1]*x[4,1,1,j]
      dx[5,1,1,j]=dx[5,1,1,j]-RateStopFailFirst[1]*x[5,1,1,j]-RateFailToSuppTreatFirstSusc[1]*x[5,1,1,j]+RateTreatToFailFirstSusc[1]*x[3,1,1,j]
      dx[6,1,1,j]=dx[6,1,1,j]-RateStopTreatSecond[1]*x[6,1,1,j]+RateDirectTreatSecond[1,1]*x[2,1,1,j]+RateDirectTreatSecond[1,1]*x[12,1,1,j]-RateTreatToFailSecond[1]*x[6,1,1,j]
      dx[7,1,1,j]=dx[7,1,1,j]-RateStopSuppSecond[1]*x[7,1,1,j]+RateFailToSuppTreatSecond[1]*x[8,1,1,j]+RateSuppFirstToSecond[1]*x[4,1,1,j]
      dx[8,1,1,j]=dx[8,1,1,j]-RateStopFailSecond[1]*x[8,1,1,j]-RateFailToSuppTreatSecond[1]*x[8,1,1,j]+RateTreatToFailSecond[1]*x[6,1,1,j]
      
      dx[1,2,1,j]=dx[1,2,1,j]
      dx[2,2,1,j]=dx[2,2,1,j]-RateDirectTreatSecond[1,2]*x[2,2,1,j]+RateStopTreatFirst[2]*x[3,2,1,j]+RateStopSuppFirst[2]*x[4,2,1,j]+RateStopFailFirst[2]*x[5,2,1,j]+
        (RateStopTreatSecond[2]*x[6,2,1,j]+RateStopSuppSecond[2]*x[7,2,1,j]+RateStopFailSecond[2]*x[8,2,1,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,2,1,j]=dx[3,2,1,j]-RateStopTreatFirst[2]*x[3,2,1,j]-RateTreatToFailFirstSusc[2]*x[3,2,1,j]
      dx[4,2,1,j]=dx[4,2,1,j]-RateStopSuppFirst[2]*x[4,2,1,j]+RateFailToSuppTreatFirstSusc[2]*x[5,2,1,j]-RateSuppFirstToSecond[2]*x[4,2,1,j]
      dx[5,2,1,j]=dx[5,2,1,j]-RateStopFailFirst[2]*x[5,2,1,j]-RateFailToSuppTreatFirstSusc[2]*x[5,2,1,j]+RateTreatToFailFirstSusc[2]*x[3,2,1,j]
      dx[6,2,1,j]=dx[6,2,1,j]-RateStopTreatSecond[2]*x[6,2,1,j]+RateDirectTreatSecond[1,2]*x[2,2,1,j]+RateDirectTreatSecond[1,2]*x[12,2,1,j]-RateTreatToFailSecond[2]*x[6,2,1,j]
      dx[7,2,1,j]=dx[7,2,1,j]-RateStopSuppSecond[2]*x[7,2,1,j]+RateFailToSuppTreatSecond[2]*x[8,2,1,j]+RateSuppFirstToSecond[2]*x[4,2,1,j]
      dx[8,2,1,j]=dx[8,2,1,j]-RateStopFailSecond[2]*x[8,2,1,j]-RateFailToSuppTreatSecond[2]*x[8,2,1,j]+RateTreatToFailSecond[2]*x[6,2,1,j]
      
      dx[1,3,1,j]=dx[1,3,1,j]
      dx[2,3,1,j]=dx[2,3,1,j]-RateDirectTreatSecond[1,3]*x[2,3,1,j]+RateStopTreatFirst[3]*x[3,3,1,j]+RateStopSuppFirst[3]*x[4,3,1,j]+RateStopFailFirst[3]*x[5,3,1,j]+
        (RateStopTreatSecond[3]*x[6,3,1,j]+RateStopSuppSecond[3]*x[7,3,1,j]+RateStopFailSecond[3]*x[8,3,1,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,3,1,j]=dx[3,3,1,j]-RateStopTreatFirst[3]*x[3,3,1,j]-RateTreatToFailFirstSusc[3]*x[3,3,1,j]
      dx[4,3,1,j]=dx[4,3,1,j]-RateStopSuppFirst[3]*x[4,3,1,j]+RateFailToSuppTreatFirstSusc[3]*x[5,3,1,j]-RateSuppFirstToSecond[3]*x[4,3,1,j]
      dx[5,3,1,j]=dx[5,3,1,j]-RateStopFailFirst[3]*x[5,3,1,j]-RateFailToSuppTreatFirstSusc[3]*x[5,3,1,j]+RateTreatToFailFirstSusc[3]*x[3,3,1,j]
      dx[6,3,1,j]=dx[6,3,1,j]-RateStopTreatSecond[3]*x[6,3,1,j]+RateDirectTreatSecond[1,3]*x[2,3,1,j]+RateDirectTreatSecond[1,3]*x[12,3,1,j]-RateTreatToFailSecond[3]*x[6,3,1,j]
      dx[7,3,1,j]=dx[7,3,1,j]-RateStopSuppSecond[3]*x[7,3,1,j]+RateFailToSuppTreatSecond[3]*x[8,3,1,j]+RateSuppFirstToSecond[3]*x[4,3,1,j]
      dx[8,3,1,j]=dx[8,3,1,j]-RateStopFailSecond[3]*x[8,3,1,j]-RateFailToSuppTreatSecond[3]*x[8,3,1,j]+RateTreatToFailSecond[3]*x[6,3,1,j]
      
      dx[1,4,1,j]=dx[1,4,1,j]
      dx[2,4,1,j]=dx[2,4,1,j]-RateDirectTreatSecond[1,4]*x[2,4,1,j]+RateStopTreatFirst[4]*x[3,4,1,j]+RateStopSuppFirst[4]*x[4,4,1,j]+RateStopFailFirst[4]*x[5,4,1,j]+
        (RateStopTreatSecond[4]*x[6,4,1,j]+RateStopSuppSecond[4]*x[7,4,1,j]+RateStopFailSecond[4]*x[8,4,1,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,4,1,j]=dx[3,4,1,j]-RateStopTreatFirst[4]*x[3,4,1,j]-RateTreatToFailFirstSusc[4]*x[3,4,1,j]
      dx[4,4,1,j]=dx[4,4,1,j]-RateStopSuppFirst[4]*x[4,4,1,j]+RateFailToSuppTreatFirstSusc[4]*x[5,4,1,j]-RateSuppFirstToSecond[4]*x[4,4,1,j]
      dx[5,4,1,j]=dx[5,4,1,j]-RateStopFailFirst[4]*x[5,4,1,j]-RateFailToSuppTreatFirstSusc[4]*x[5,4,1,j]+RateTreatToFailFirstSusc[4]*x[3,4,1,j]
      dx[6,4,1,j]=dx[6,4,1,j]-RateStopTreatSecond[4]*x[6,4,1,j]+RateDirectTreatSecond[1,4]*x[2,4,1,j]+RateDirectTreatSecond[1,4]*x[12,4,1,j]-RateTreatToFailSecond[4]*x[6,4,1,j]
      dx[7,4,1,j]=dx[7,4,1,j]-RateStopSuppSecond[4]*x[7,4,1,j]+RateFailToSuppTreatSecond[4]*x[8,4,1,j]+RateSuppFirstToSecond[4]*x[4,4,1,j]
      dx[8,4,1,j]=dx[8,4,1,j]-RateStopFailSecond[4]*x[8,4,1,j]-RateFailToSuppTreatSecond[4]*x[8,4,1,j]+RateTreatToFailSecond[4]*x[6,4,1,j]
      #############################################################################################################################################################################
      #############################################################################################################################################################################
      dx[1,1,2,j]=S[j]/N[j]*(sum(RateInf1[1:4,1:2,j]*apply(x[c(1),1:4,2,1:2],c(1,2),sum))+sum(RateInf2[1:4,1:2,j]*apply(x[c(2,3,5,6,8,9,11,12,13,15),1:4,2,1:2],c(2,3),sum)))-
        (RateStageInf[1]+RateDiag[1,j]+mu[1,1,2 ])*x[1,1,2,j]
      dx[2,1,2,j]=dx[2,1,2,j]+RateDiag[1,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,1,2,j]-(RateTreatFirst[2,1,j]+RateStageDiag[1]+mu[2,1,2])*x[2,1,2,j]
      dx[3,1,2,j]=dx[3,1,2,j]+RateTreatFirst[2,1,j]*x[2,1,2,j]-(RateSuppFirstResis[1]+RateStageTreatFirstR[1]+mu[3,1,2])*x[3,1,2,j]+RateStageTreatFirstL[1]*x[3,2,2,j]
      dx[4,1,2,j]=RateSuppFirstResis[1]*x[3,1,2,j]-(RateFailFirstResis[1]+mu[4,1,2])*x[4,1,2,j]+RateStageSuppFirst[1]*x[4,2,2,j]
      dx[5,1,2,j]=RateFailFirstResis[1]*x[4,1,2,j]-(RateTreatSecond[1,j]+RateStageFailFirst[1]+mu[5,1,2 ])*x[5,1,2,j]
      dx[6,1,2,j]=RateTreatSecond_noDTG[1,j]*x[11,1,2,j]+RateTreatSecond_noDTG[1,j]*x[15,1,2,j]+RateTreatSecond[1,j]*x[5,1,2,j]-(RateSuppSecond[1]+RateStageTreatSecondR[1]+mu[6,1,2])*x[6,1,2,j]+RateStageTreatSecondL[1]*x[6,2,2,j]
      dx[7,1,2,j]=RateSuppSecond[1]*x[6,1,2,j]-(RateFailSecond[1]+mu[7,1,2])*x[7,1,2,j]+RateStageSuppSecond[1]*x[7,2,2,j]
      dx[8,1,2,j]=RateFailSecond[1]*x[7,1,2,j]-(RateStageFailSecond[1]+mu[8,1,2])*x[8,1,2,j]
      
      dx[1,2,2,j]=-(RateStageInf[2]+RateDiag[2,j]+mu[1,2,2])*x[1,2,2,j]+RateStageInf[1]*x[1,1,2,j]
      dx[2,2,2,j]=dx[2,2,2,j]+RateDiag[2,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,2,2,j]-(RateTreatFirst[2,2,j]+RateStageDiag[2]+mu[2,2,2])*x[2,2,2,j]+RateStageDiag[1]*x[2,1,2,j]
      dx[3,2,2,j]=dx[3,2,2,j]+RateTreatFirst[2,2,j]*x[2,2,2,j]-(RateSuppFirstResis[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,2])*x[3,2,2,j]+RateStageTreatFirstR[1]*x[3,1,2,j]+RateStageTreatFirstL[2]*x[3,3,2,j]
      dx[4,2,2,j]=RateSuppFirstResis[2]*x[3,2,2,j]-(RateFailFirstResis[2]+RateStageSuppFirst[1]+mu[4,2,2])*x[4,2,2,j]+RateStageSuppFirst[2]*x[4,3,2,j]
      dx[5,2,2,j]=RateFailFirstResis[2]*x[4,2,2,j]-(RateTreatSecond[2,j]+RateStageFailFirst[2]+mu[5,2,2])*x[5,2,2,j]+RateStageFailFirst[1]*x[5,1,2,j]
      dx[6,2,2,j]=RateTreatSecond_noDTG[2,j]*x[11,2,2,j]+RateTreatSecond_noDTG[2,j]*x[15,2,2,j]+RateTreatSecond[2,j]*x[5,2,2,j]-(RateSuppSecond[2]+RateStageTreatSecondR[2]+RateStageTreatSecondL[1]+mu[6,2,2])*x[6,2,2,j]+RateStageTreatSecondR[1]*x[6,1,2,j]+RateStageTreatSecondL[2]*x[6,3,2,j]
      dx[7,2,2,j]=RateSuppSecond[2]*x[6,2,2,j]-(RateFailSecond[2]+RateStageSuppSecond[1]+mu[7,2,2])*x[7,2,2,j]+RateStageSuppSecond[2]*x[7,3,2,j]
      dx[8,2,2,j]=RateFailSecond[2]*x[7,2,2,j]-(RateStageFailSecond[2]+mu[8,2,2])*x[8,2,2,j]+RateStageFailSecond[1]*x[8,1,2,j]
      
      dx[1,3,2,j]=-(RateStageInf[3]+RateDiag[3,j]+mu[1,3,2])*x[1,3,2,j]+RateStageInf[2]*x[1,2,2,j]
      dx[2,3,2,j]=dx[2,3,2,j]+RateDiag[3,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,3,2,j]-(RateTreatFirst[2,3,j]+RateStageDiag[3]+mu[2,3,2])*x[2,3,2,j]+RateStageDiag[2]*x[2,2,2,j]
      dx[3,3,2,j]=dx[3,3,2,j]+RateTreatFirst[2,3,j]*x[2,3,2,j]-(RateSuppFirstResis[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,2])*x[3,3,2,j]+RateStageTreatFirstR[2]*x[3,2,2,j]+RateStageTreatFirstL[3]*x[3,4,2,j]
      dx[4,3,2,j]=RateSuppFirstResis[3]*x[3,3,2,j]-(RateFailFirstResis[3]+RateStageSuppFirst[2]+mu[4,3,2])*x[4,3,2,j]+RateStageSuppFirst[3]*x[4,4,2,j]
      dx[5,3,2,j]=RateFailFirstResis[3]*x[4,3,2,j]-(RateTreatSecond[3,j]+RateStageFailFirst[3]+mu[5,3,2])*x[5,3,2,j]+RateStageFailFirst[2]*x[5,2,2,j]
      dx[6,3,2,j]=RateTreatSecond_noDTG[3,j]*x[11,3,2,j]+RateTreatSecond_noDTG[3,j]*x[15,3,2,j]+RateTreatSecond[3,j]*x[5,3,2,j]-(RateSuppSecond[3]+RateStageTreatSecondR[3]+RateStageTreatSecondL[2]+mu[6,3,2])*x[6,3,2,j]+RateStageTreatSecondR[2]*x[6,2,2,j]+RateStageTreatSecondL[3]*x[6,4,2,j]
      dx[7,3,2,j]=RateSuppSecond[3]*x[6,3,2,j]-(RateFailSecond[3]+RateStageSuppSecond[2]+mu[7,3,2])*x[7,3,2,j]+RateStageSuppSecond[3]*x[7,4,2,j]
      dx[8,3,2,j]=RateFailSecond[3]*x[7,3,2,j]-(RateStageFailSecond[3]+mu[8,3,2])*x[8,3,2,j]+RateStageFailSecond[2]*x[8,2,2,j]
      
      dx[1,4,2,j]=-(RateDiag[4,j]+mu[1,4,2])*x[1,4,2,j]+RateStageInf[3]*x[1,3,2,j]
      dx[2,4,2,j]=dx[2,4,2,j]+RateDiag[4,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,4,2,j]-(RateTreatFirst[2,4,j]+mu[2,4,2])*x[2,4,2,j]+RateStageDiag[3]*x[2,3,2,j]
      dx[3,4,2,j]=dx[3,4,2,j]+RateTreatFirst[2,4,j]*x[2,4,2,j]-(RateSuppFirstResis[4]+RateStageTreatFirstL[3]+mu[3,4,2])*x[3,4,2,j]+RateStageTreatFirstR[3]*x[3,3,2,j]
      dx[4,4,2,j]=RateSuppFirstResis[4]*x[3,4,2,j]-(RateFailFirstResis[4]+RateStageSuppFirst[3]+mu[4,4,2])*x[4,4,2,j]
      dx[5,4,2,j]=RateFailFirstResis[4]*x[4,4,2,j]-(RateTreatSecond[4,j]+mu[5,4,2])*x[5,4,2,j]+RateStageFailFirst[3]*x[5,3,2,j]
      dx[6,4,2,j]=RateTreatSecond_noDTG[4,j]*x[11,4,2,j]+RateTreatSecond_noDTG[4,j]*x[15,4,2,j]+RateTreatSecond[4,j]*x[5,4,2,j]-(RateSuppSecond[4]+RateStageTreatSecondL[3]+mu[6,4,2])*x[6,4,2,j]+RateStageTreatSecondR[3]*x[6,3,2,j]
      dx[7,4,2,j]=RateSuppSecond[4]*x[6,4,2,j]-(RateFailSecond[4]+RateStageSuppSecond[3]+mu[7,4,2])*x[7,4,2,j]
      dx[8,4,2,j]=RateFailSecond[4]*x[7,4,2,j]-mu[8,4,2]*x[8,4,2,j]+RateStageFailSecond[3]*x[8,3,2,j]
      
      ##############################################################################################################################################################
      
      dx[1,1,2,j]=dx[1,1,2,j]
      dx[2,1,2,j]=dx[2,1,2,j]-RateDirectTreatSecond[2,1]*x[2,1,2,j]+RateStopTreatFirst[1]*x[3,1,2,j]+RateStopSuppFirst[1]*x[4,1,2,j]+RateStopFailFirst[1]*x[5,1,2,j]+
        (RateStopTreatSecond[1]*x[6,1,2,j]+RateStopSuppSecond[1]*x[7,1,2,j]+RateStopFailSecond[1]*x[8,1,2,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,1,2,j]=dx[3,1,2,j]-RateStopTreatFirst[1]*x[3,1,2,j]-RateTreatToFailFirstResis[1]*x[3,1,2,j]
      dx[4,1,2,j]=dx[4,1,2,j]-RateStopSuppFirst[1]*x[4,1,2,j]+RateFailToSuppTreatFirstResis[1]*x[5,1,2,j]-RateSuppFirstToSecond[1]*x[4,1,2,j]
      dx[5,1,2,j]=dx[5,1,2,j]-RateStopFailFirst[1]*x[5,1,2,j]-RateFailToSuppTreatFirstResis[1]*x[5,1,2,j]+RateTreatToFailFirstResis[1]*x[3,1,2,j]
      dx[6,1,2,j]=dx[6,1,2,j]-RateStopTreatSecond[1]*x[6,1,2,j]+RateDirectTreatSecond[2,1]*x[2,1,2,j]+RateDirectTreatSecond[2,1]*x[12,1,2,j]-RateTreatToFailSecond[1]*x[6,1,2,j]
      dx[7,1,2,j]=dx[7,1,2,j]-RateStopSuppSecond[1]*x[7,1,2,j]+RateFailToSuppTreatSecond[1]*x[8,1,2,j]+RateSuppFirstToSecond[1]*x[4,1,2,j]
      dx[8,1,2,j]=dx[8,1,2,j]-RateStopFailSecond[1]*x[8,1,2,j]-RateFailToSuppTreatSecond[1]*x[8,1,2,j]+RateTreatToFailSecond[1]*x[6,1,2,j]
      
      dx[1,2,2,j]=dx[1,2,2,j]
      dx[2,2,2,j]=dx[2,2,2,j]-RateDirectTreatSecond[2,2]*x[2,2,2,j]+RateStopTreatFirst[2]*x[3,2,2,j]+RateStopSuppFirst[2]*x[4,2,2,j]+RateStopFailFirst[2]*x[5,2,2,j]+
        (RateStopTreatSecond[2]*x[6,2,2,j]+RateStopSuppSecond[2]*x[7,2,2,j]+RateStopFailSecond[2]*x[8,2,2,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,2,2,j]=dx[3,2,2,j]-RateStopTreatFirst[2]*x[3,2,2,j]-RateTreatToFailFirstResis[2]*x[3,2,2,j]
      dx[4,2,2,j]=dx[4,2,2,j]-RateStopSuppFirst[2]*x[4,2,2,j]+RateFailToSuppTreatFirstResis[2]*x[5,2,2,j]-RateSuppFirstToSecond[2]*x[4,2,2,j]
      dx[5,2,2,j]=dx[5,2,2,j]-RateStopFailFirst[2]*x[5,2,2,j]-RateFailToSuppTreatFirstResis[2]*x[5,2,2,j]+RateTreatToFailFirstResis[2]*x[3,2,2,j]
      dx[6,2,2,j]=dx[6,2,2,j]-RateStopTreatSecond[2]*x[6,2,2,j]+RateDirectTreatSecond[2,2]*x[2,2,2,j]+RateDirectTreatSecond[2,2]*x[12,2,2,j]-RateTreatToFailSecond[2]*x[6,2,2,j]
      dx[7,2,2,j]=dx[7,2,2,j]-RateStopSuppSecond[2]*x[7,2,2,j]+RateFailToSuppTreatSecond[2]*x[8,2,2,j]+RateSuppFirstToSecond[2]*x[4,2,2,j]
      dx[8,2,2,j]=dx[8,2,2,j]-RateStopFailSecond[2]*x[8,2,2,j]-RateFailToSuppTreatSecond[2]*x[8,2,2,j]+RateTreatToFailSecond[2]*x[6,2,2,j]
      
      dx[1,3,2,j]=dx[1,3,2,j]
      dx[2,3,2,j]=dx[2,3,2,j]-RateDirectTreatSecond[2,3]*x[2,3,2,j]+RateStopTreatFirst[3]*x[3,3,2,j]+RateStopSuppFirst[3]*x[4,3,2,j]+RateStopFailFirst[3]*x[5,3,2,j]+
        (RateStopTreatSecond[3]*x[6,3,2,j]+RateStopSuppSecond[3]*x[7,3,2,j]+RateStopFailSecond[3]*x[8,3,2,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,3,2,j]=dx[3,3,2,j]-RateStopTreatFirst[3]*x[3,3,2,j]-RateTreatToFailFirstResis[3]*x[3,3,2,j]
      dx[4,3,2,j]=dx[4,3,2,j]-RateStopSuppFirst[3]*x[4,3,2,j]+RateFailToSuppTreatFirstResis[3]*x[5,3,2,j]-RateSuppFirstToSecond[3]*x[4,3,2,j]
      dx[5,3,2,j]=dx[5,3,2,j]-RateStopFailFirst[3]*x[5,3,2,j]-RateFailToSuppTreatFirstResis[3]*x[5,3,2,j]+RateTreatToFailFirstResis[3]*x[3,3,2,j]
      dx[6,3,2,j]=dx[6,3,2,j]-RateStopTreatSecond[3]*x[6,3,2,j]+RateDirectTreatSecond[2,3]*x[2,3,2,j]+RateDirectTreatSecond[2,3]*x[12,3,2,j]-RateTreatToFailSecond[3]*x[6,3,2,j]
      dx[7,3,2,j]=dx[7,3,2,j]-RateStopSuppSecond[3]*x[7,3,2,j]+RateFailToSuppTreatSecond[3]*x[8,3,2,j]+RateSuppFirstToSecond[3]*x[4,3,2,j]
      dx[8,3,2,j]=dx[8,3,2,j]-RateStopFailSecond[3]*x[8,3,2,j]-RateFailToSuppTreatSecond[3]*x[8,3,2,j]+RateTreatToFailSecond[3]*x[6,3,2,j]
      
      dx[1,4,2,j]=dx[1,4,2,j]
      dx[2,4,2,j]=dx[2,4,2,j]-RateDirectTreatSecond[2,4]*x[2,4,2,j]+RateStopTreatFirst[4]*x[3,4,2,j]+RateStopSuppFirst[4]*x[4,4,2,j]+RateStopFailFirst[4]*x[5,4,2,j]+
        (RateStopTreatSecond[4]*x[6,4,2,j]+RateStopSuppSecond[4]*x[7,4,2,j]+RateStopFailSecond[4]*x[8,4,2,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,4,2,j]=dx[3,4,2,j]-RateStopTreatFirst[4]*x[3,4,2,j]-RateTreatToFailFirstResis[4]*x[3,4,2,j]
      dx[4,4,2,j]=dx[4,4,2,j]-RateStopSuppFirst[4]*x[4,4,2,j]+RateFailToSuppTreatFirstResis[4]*x[5,4,2,j]-RateSuppFirstToSecond[4]*x[4,4,2,j]
      dx[5,4,2,j]=dx[5,4,2,j]-RateStopFailFirst[4]*x[5,4,2,j]-RateFailToSuppTreatFirstResis[4]*x[5,4,2,j]+RateTreatToFailFirstResis[4]*x[3,4,2,j]
      dx[6,4,2,j]=dx[6,4,2,j]-RateStopTreatSecond[4]*x[6,4,2,j]+RateDirectTreatSecond[2,4]*x[2,4,2,j]+RateDirectTreatSecond[2,4]*x[12,4,2,j]-RateTreatToFailSecond[4]*x[6,4,2,j]
      dx[7,4,2,j]=dx[7,4,2,j]-RateStopSuppSecond[4]*x[7,4,2,j]+RateFailToSuppTreatSecond[4]*x[8,4,2,j]+RateSuppFirstToSecond[4]*x[4,4,2,j]
      dx[8,4,2,j]=dx[8,4,2,j]-RateStopFailSecond[4]*x[8,4,2,j]-RateFailToSuppTreatSecond[4]*x[8,4,2,j]+RateTreatToFailSecond[4]*x[6,4,2,j]
      
      #############################################################################################################################################################################
      #############################################################################################################################################################################
      #dx with interaction between resistant and susceptible compartments
      dx[1,1,1,j]=dx[1,1,1,j]+RateSusceptible*x[1,1,2,j]
      dx[2,1,1,j]=dx[2,1,1,j]+RateSusceptible*x[2,1,2,j]
      dx[5,1,1,j]=dx[5,1,1,j]-RateResistant*x[5,1,1,j]
      
      dx[1,2,1,j]=dx[1,2,1,j]+RateSusceptible*x[1,2,2,j]
      dx[2,2,1,j]=dx[2,2,1,j]+RateSusceptible*x[2,2,2,j]
      dx[5,2,1,j]=dx[5,2,1,j]-RateResistant*x[5,2,1,j]
      
      dx[1,3,1,j]=dx[1,3,1,j]+RateSusceptible*x[1,3,2,j]
      dx[2,3,1,j]=dx[2,3,1,j]+RateSusceptible*x[2,3,2,j]
      dx[5,3,1,j]=dx[5,3,1,j]-RateResistant*x[5,3,1,j]
      
      dx[1,4,1,j]=dx[1,4,1,j]+RateSusceptible*x[1,4,2,j]
      dx[2,4,1,j]=dx[2,4,1,j]+RateSusceptible*x[2,4,2,j]
      dx[5,4,1,j]=dx[5,4,1,j]-RateResistant*x[5,4,1,j]
      
      #############################################################################################################################################################################
      dx[1,1,2,j]=dx[1,1,2,j]-RateSusceptible*x[1,1,2,j]
      dx[2,1,2,j]=dx[2,1,2,j]-RateSusceptible*x[2,1,2,j]
      dx[5,1,2,j]=dx[5,1,2,j]+RateResistant*x[5,1,1,j]
      
      dx[1,2,2,j]=dx[1,2,2,j]-RateSusceptible*x[1,2,2,j]
      dx[2,2,2,j]=dx[2,2,2,j]-RateSusceptible*x[2,2,2,j]
      dx[5,2,2,j]=dx[5,2,2,j]+RateResistant*x[5,2,1,j]
      
      dx[1,3,2,j]=dx[1,3,2,j]-RateSusceptible*x[1,3,2,j]
      dx[2,3,2,j]=dx[2,3,2,j]-RateSusceptible*x[2,3,2,j]
      dx[5,3,2,j]=dx[5,3,2,j]+RateResistant*x[5,3,1,j]
      
      dx[1,4,2,j]=dx[1,4,2,j]-RateSusceptible*x[1,4,2,j]
      dx[2,4,2,j]=dx[2,4,2,j]-RateSusceptible*x[2,4,2,j]
      dx[5,4,2,j]=dx[5,4,2,j]+RateResistant*x[5,4,1,j]
      
      
      #RateTreatFirstDTG[res,cd4,gender]------------------------------
      #NNRTI susceptible
      dx[2,1,1,j]=dx[2,1,1,j]-RateTreatFirstDTG[1,1,j]*x[2,1,1,j]
      dx[2,2,1,j]=dx[2,2,1,j]-RateTreatFirstDTG[1,2,j]*x[2,2,1,j]
      dx[2,3,1,j]=dx[2,3,1,j]-RateTreatFirstDTG[1,3,j]*x[2,3,1,j]
      dx[2,4,1,j]=dx[2,4,1,j]-RateTreatFirstDTG[1,4,j]*x[2,4,1,j]
      #NNRTI resistant
      dx[2,1,2,j]=dx[2,1,2,j]-RateTreatFirstDTG[2,1,j]*x[2,1,2,j]
      dx[2,2,2,j]=dx[2,2,2,j]-RateTreatFirstDTG[2,2,j]*x[2,2,2,j]
      dx[2,3,2,j]=dx[2,3,2,j]-RateTreatFirstDTG[2,3,j]*x[2,3,2,j]
      dx[2,4,2,j]=dx[2,4,2,j]-RateTreatFirstDTG[2,4,j]*x[2,4,2,j]
      
      #RateSwitchDTG[care,cd4,res,gender]----------------------------
      #NNRTI susceptible
      dx[3,1,1,j]=dx[3,1,1,j]-RateSwitchDTG[1,1,1,j]*x[3,1,1,j]
      dx[3,2,1,j]=dx[3,2,1,j]-RateSwitchDTG[1,2,1,j]*x[3,2,1,j]
      dx[3,3,1,j]=dx[3,3,1,j]-RateSwitchDTG[1,3,1,j]*x[3,3,1,j]
      dx[3,4,1,j]=dx[3,4,1,j]-RateSwitchDTG[1,4,1,j]*x[3,4,1,j]
      
      dx[4,1,1,j]=dx[4,1,1,j]-RateSwitchDTG[2,1,1,j]*x[4,1,1,j]
      dx[4,2,1,j]=dx[4,2,1,j]-RateSwitchDTG[2,2,1,j]*x[4,2,1,j]
      dx[4,3,1,j]=dx[4,3,1,j]-RateSwitchDTG[2,3,1,j]*x[4,3,1,j]
      dx[4,4,1,j]=dx[4,4,1,j]-RateSwitchDTG[2,4,1,j]*x[4,4,1,j]
      
      dx[5,1,1,j]=dx[5,1,1,j]-RateSwitchDTG[3,1,1,j]*x[5,1,1,j]
      dx[5,2,1,j]=dx[5,2,1,j]-RateSwitchDTG[3,2,1,j]*x[5,2,1,j]
      dx[5,3,1,j]=dx[5,3,1,j]-RateSwitchDTG[3,3,1,j]*x[5,3,1,j]
      dx[5,4,1,j]=dx[5,4,1,j]-RateSwitchDTG[3,4,1,j]*x[5,4,1,j]
      
      #NNRTI resistant
      dx[3,1,2,j]=dx[3,1,2,j]-RateSwitchDTG[1,1,2,j]*x[3,1,2,j]
      dx[3,2,2,j]=dx[3,2,2,j]-RateSwitchDTG[1,2,2,j]*x[3,2,2,j]
      dx[3,3,2,j]=dx[3,3,2,j]-RateSwitchDTG[1,3,2,j]*x[3,3,2,j]
      dx[3,4,2,j]=dx[3,4,2,j]-RateSwitchDTG[1,4,2,j]*x[3,4,2,j]
      
      dx[4,1,2,j]=dx[4,1,2,j]-RateSwitchDTG[2,1,2,j]*x[4,1,2,j]
      dx[4,2,2,j]=dx[4,2,2,j]-RateSwitchDTG[2,2,2,j]*x[4,2,2,j]
      dx[4,3,2,j]=dx[4,3,2,j]-RateSwitchDTG[2,3,2,j]*x[4,3,2,j]
      dx[4,4,2,j]=dx[4,4,2,j]-RateSwitchDTG[2,4,2,j]*x[4,4,2,j]
      
      dx[5,1,2,j]=dx[5,1,2,j]-RateSwitchDTG[3,1,2,j]*x[5,1,2,j]
      dx[5,2,2,j]=dx[5,2,2,j]-RateSwitchDTG[3,2,2,j]*x[5,2,2,j]
      dx[5,3,2,j]=dx[5,3,2,j]-RateSwitchDTG[3,3,2,j]*x[5,3,2,j]
      dx[5,4,2,j]=dx[5,4,2,j]-RateSwitchDTG[3,4,2,j]*x[5,4,2,j]
      
      #DTG compartments-------------------------------------
      #NNRTI susceptible
      dx[9,1,1,j]=dx[9,1,1,j]+RateTreatFirstDTG[1,1,j]*x[2,1,1,j]+RateSwitchDTG[3,1,1,j]*x[5,1,1,j]+RateSwitchDTG[1,1,1,j]*x[3,1,1,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,1])*x[9,1,1,j]+RateStageTreatFirstL[1]*x[9,2,1,j]
      dx[10,1,1,j]=RateSuppFirstSusc[1]*x[9,1,1,j]-(RateFailFirstDTG[1]+mu[4,1,1])*x[10,1,1,j]+RateStageSuppFirst[1]*x[10,2,1,j]+RateSwitchDTG[2,1,1,j]*x[4,1,1,j]
      dx[11,1,1,j]=RateFailFirstDTG[1]*x[10,1,1,j]-(RateTreatSecond_noDTG[1,j]+RateStageFailFirst[1]+mu[5,1,1])*x[11,1,1,j]
      dx[12,1,1,j]=dx[12,1,1,j]+RateDiag[1,j]*((j==2)*(1-p_women_dtg))*x[1,1,1,j]-(RateTreatFirst_noDTG[1,1,j]+RateStageDiag[1]+mu[2,1,1])*x[12,1,1,j]
      dx[13,1,1,j]=dx[13,1,1,j]+RateTreatFirst_noDTG[1,1,j]*x[12,1,1,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,1])*x[13,1,1,j]+RateStageTreatFirstL[1]*x[13,2,1,j]
      dx[14,1,1,j]=dx[14,1,1,j]+RateSuppFirstSusc[1]*x[13,1,1,j]-(RateFailFirstSusc[1]+mu[4,1,1])*x[14,1,1,j]+RateStageSuppFirst[1]*x[14,2,1,j]
      dx[15,1,1,j]=dx[15,1,1,j]+RateFailFirstSusc[1]*x[14,1,1,j]-(RateTreatSecond_noDTG[1,j]+RateStageFailFirst[1]+mu[5,1,1])*x[15,1,1,j]
      
      dx[9,2,1,j]=dx[9,2,1,j]+RateTreatFirstDTG[1,2,j]*x[2,2,1,j]+RateSwitchDTG[3,2,1,j]*x[5,2,1,j]+RateSwitchDTG[1,2,1,j]*x[3,2,1,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,1])*x[9,2,1,j]+RateStageTreatFirstR[1]*x[9,1,1,j]+RateStageTreatFirstL[2]*x[9,3,1,j]
      dx[10,2,1,j]=RateSuppFirstSusc[2]*x[9,2,1,j]-(RateFailFirstDTG[2]+RateStageSuppFirst[1]+mu[4,2,1])*x[10,2,1,j]+RateStageSuppFirst[2]*x[10,3,1,j]+RateSwitchDTG[2,2,1,j]*x[4,2,1,j]
      dx[11,2,1,j]=RateFailFirstDTG[2]*x[10,2,1,j]-(RateTreatSecond_noDTG[2,j]+RateStageFailFirst[2]+mu[5,2,1])*x[11,2,1,j]+RateStageFailFirst[1]*x[11,1,1,j]
      dx[12,2,1,j]=dx[12,2,1,j]+RateDiag[2,j]*((j==2)*(1-p_women_dtg))*x[1,2,1,j]-(RateTreatFirst_noDTG[1,2,j]+RateStageDiag[2]+mu[2,2,1])*x[12,2,1,j]+RateStageDiag[1]*x[12,1,1,j]
      dx[13,2,1,j]=dx[13,2,1,j]+RateTreatFirst_noDTG[1,2,j]*x[12,2,1,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,1])*x[13,2,1,j]+RateStageTreatFirstR[1]*x[13,1,1,j]+RateStageTreatFirstL[2]*x[13,3,1,j]
      dx[14,2,1,j]=dx[14,2,1,j]+RateSuppFirstSusc[2]*x[13,2,1,j]-(RateFailFirstSusc[2]+RateStageSuppFirst[1]+mu[4,2,1])*x[14,2,1,j]+RateStageSuppFirst[2]*x[14,3,1,j]
      dx[15,2,1,j]=dx[15,2,1,j]+RateFailFirstSusc[2]*x[14,2,1,j]-(RateTreatSecond_noDTG[2,j]+RateStageFailFirst[2]+mu[5,2,1])*x[15,2,1,j]+RateStageFailFirst[1]*x[15,1,1,j]
      
      dx[9,3,1,j]=dx[9,3,1,j]+RateTreatFirstDTG[1,3,j]*x[2,3,1,j]+RateSwitchDTG[3,3,1,j]*x[5,3,1,j]+RateSwitchDTG[1,3,1,j]*x[3,3,1,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,1])*x[9,3,1,j]+RateStageTreatFirstR[2]*x[9,2,1,j]+RateStageTreatFirstL[3]*x[9,4,1,j]
      dx[10,3,1,j]=RateSuppFirstSusc[3]*x[9,3,1,j]-(RateFailFirstDTG[3]+RateStageSuppFirst[2]+mu[4,3,1])*x[10,3,1,j]+RateStageSuppFirst[3]*x[10,4,1,j]+RateSwitchDTG[2,3,1,j]*x[4,3,1,j]
      dx[11,3,1,j]=RateFailFirstDTG[3]*x[10,3,1,j]-(RateTreatSecond_noDTG[3,j]+RateStageFailFirst[3]+mu[5,3,1])*x[11,3,1,j]+RateStageFailFirst[2]*x[11,2,1,j]
      dx[12,3,1,j]=dx[12,3,1,j]+RateDiag[3,j]*((j==2)*(1-p_women_dtg))*x[1,3,1,j]-(RateTreatFirst_noDTG[1,3,j]+RateStageDiag[3]+mu[2,3,1])*x[12,3,1,j]+RateStageDiag[2]*x[12,2,1,j]
      dx[13,3,1,j]=dx[13,3,1,j]+RateTreatFirst_noDTG[1,3,j]*x[12,3,1,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,1])*x[13,3,1,j]+RateStageTreatFirstR[2]*x[13,2,1,j]+RateStageTreatFirstL[3]*x[13,4,1,j]
      dx[14,3,1,j]=dx[14,3,1,j]+RateSuppFirstSusc[3]*x[13,3,1,j]-(RateFailFirstSusc[3]+RateStageSuppFirst[2]+mu[4,3,1])*x[14,3,1,j]+RateStageSuppFirst[3]*x[14,4,1,j]
      dx[15,3,1,j]=dx[15,3,1,j]+RateFailFirstSusc[3]*x[14,3,1,j]-(RateTreatSecond_noDTG[3,j]+RateStageFailFirst[3]+mu[5,3,1])*x[15,3,1,j]+RateStageFailFirst[2]*x[15,2,1,j]
      
      dx[9,4,1,j]=dx[9,4,1,j]+RateTreatFirstDTG[1,4,j]*x[2,4,1,j]+RateSwitchDTG[3,4,1,j]*x[5,4,1,j]+RateSwitchDTG[1,4,1,j]*x[3,4,1,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,1])*x[9,4,1,j]+RateStageTreatFirstR[3]*x[9,3,1,j]
      dx[10,4,1,j]=RateSuppFirstSusc[4]*x[9,4,1,j]-(RateFailFirstDTG[4]+RateStageSuppFirst[3]+mu[4,4,1])*x[10,4,1,j]+RateSwitchDTG[2,4,1,j]*x[4,4,1,j]
      dx[11,4,1,j]=RateFailFirstDTG[4]*x[10,4,1,j]-(RateTreatSecond_noDTG[4,j]+mu[5,4,1])*x[11,4,1,j]+RateStageFailFirst[3]*x[11,3,1,j]
      dx[12,4,1,j]=dx[12,4,1,j]+RateDiag[4,j]*((j==2)*(1-p_women_dtg))*x[1,4,1,j]-(RateTreatFirst_noDTG[1,4,j]+mu[2,4,1])*x[12,4,1,j]+RateStageDiag[3]*x[12,3,1,j]
      dx[13,4,1,j]=dx[13,4,1,j]+RateTreatFirst_noDTG[1,4,j]*x[12,4,1,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,1])*x[13,4,1,j]+RateStageTreatFirstR[3]*x[13,3,1,j]
      dx[14,4,1,j]=dx[14,4,1,j]+RateSuppFirstSusc[4]*x[13,4,1,j]-(RateFailFirstSusc[4]+RateStageSuppFirst[3]+mu[4,4,1])*x[14,4,1,j]
      dx[15,4,1,j]=dx[15,4,1,j]+RateFailFirstSusc[4]*x[14,4,1,j]-(RateTreatSecond_noDTG[4,j]+mu[5,4,1])*x[15,4,1,j]+RateStageFailFirst[3]*x[15,3,1,j]
      
      #NNRTI resistant
      dx[9,1,2,j]=dx[9,1,2,j]+RateTreatFirstDTG[2,1,j]*x[2,1,2,j]+RateSwitchDTG[3,1,2,j]*x[5,1,2,j]+RateSwitchDTG[1,1,2,j]*x[3,1,2,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,2])*x[9,1,2,j]+RateStageTreatFirstL[1]*x[9,2,2,j]
      dx[10,1,2,j]=RateSuppFirstSusc[1]*x[9,1,2,j]-(RateFailFirstDTG[1]+mu[4,1,2])*x[10,1,2,j]+RateStageSuppFirst[1]*x[10,2,2,j]+RateSwitchDTG[2,1,2,j]*x[4,1,2,j]
      dx[11,1,2,j]=RateFailFirstDTG[1]*x[10,1,2,j]-(RateTreatSecond_noDTG[1,j]+RateStageFailFirst[1]+mu[5,1,2])*x[11,1,2,j]
      dx[12,1,2,j]=dx[12,1,2,j]+RateDiag[1,j]*((j==2)*(1-p_women_dtg))*x[1,1,2,j]-(RateTreatFirst_noDTG[2,1,j]+RateStageDiag[1]+mu[2,1,2])*x[12,1,2,j]
      dx[13,1,2,j]=dx[13,1,2,j]+RateTreatFirst_noDTG[2,1,j]*x[12,1,2,j]-(RateSuppFirstResis[1]+RateStageTreatFirstR[1]+mu[3,1,2])*x[13,1,2,j]+RateStageTreatFirstL[1]*x[13,2,2,j]
      dx[14,1,2,j]=dx[14,1,2,j]+RateSuppFirstResis[1]*x[13,1,2,j]-(RateFailFirstResis[1]+mu[4,1,2])*x[14,1,2,j]+RateStageSuppFirst[1]*x[14,2,2,j]
      dx[15,1,2,j]=dx[15,1,2,j]+RateFailFirstResis[1]*x[14,1,2,j]-(RateTreatSecond_noDTG[1,j]+RateStageFailFirst[1]+mu[5,1,2])*x[15,1,2,j]
      
      dx[9,2,2,j]=dx[9,2,2,j]+RateTreatFirstDTG[2,2,j]*x[2,2,2,j]+RateSwitchDTG[3,2,2,j]*x[5,2,2,j]+RateSwitchDTG[1,2,2,j]*x[3,2,2,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,2])*x[9,2,2,j]+RateStageTreatFirstR[1]*x[9,1,2,j]+RateStageTreatFirstL[2]*x[9,3,2,j]
      dx[10,2,2,j]=RateSuppFirstSusc[2]*x[9,2,2,j]-(RateFailFirstDTG[2]+RateStageSuppFirst[1]+mu[4,2,2])*x[10,2,2,j]+RateStageSuppFirst[2]*x[10,3,2,j]+RateSwitchDTG[2,2,2,j]*x[4,2,2,j]
      dx[11,2,2,j]=RateFailFirstDTG[2]*x[10,2,2,j]-(RateTreatSecond_noDTG[2,j]+RateStageFailFirst[2]+mu[5,2,2])*x[11,2,2,j]+RateStageFailFirst[1]*x[11,1,2,j]
      dx[12,2,2,j]=dx[12,2,2,j]+RateDiag[2,j]*((j==2)*(1-p_women_dtg))*x[1,2,2,j]-(RateTreatFirst_noDTG[2,2,j]+RateStageDiag[2]+mu[2,2,2])*x[12,2,2,j]+RateStageDiag[1]*x[12,1,2,j]
      dx[13,2,2,j]=dx[13,2,2,j]+RateTreatFirst_noDTG[2,2,j]*x[12,2,2,j]-(RateSuppFirstResis[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,2])*x[13,2,2,j]+RateStageTreatFirstR[1]*x[13,1,2,j]+RateStageTreatFirstL[2]*x[13,3,2,j]
      dx[14,2,2,j]=dx[14,2,2,j]+RateSuppFirstResis[2]*x[13,2,2,j]-(RateFailFirstResis[2]+RateStageSuppFirst[1]+mu[4,2,2])*x[14,2,2,j]+RateStageSuppFirst[2]*x[14,3,2,j]
      dx[15,2,2,j]=dx[15,2,2,j]+RateFailFirstResis[2]*x[14,2,2,j]-(RateTreatSecond_noDTG[2,j]+RateStageFailFirst[2]+mu[5,2,2])*x[15,2,2,j]+RateStageFailFirst[1]*x[15,1,2,j]
      
      dx[9,3,2,j]=dx[9,3,2,j]+RateTreatFirstDTG[2,3,j]*x[2,3,2,j]+RateSwitchDTG[3,3,2,j]*x[5,3,2,j]+RateSwitchDTG[1,3,2,j]*x[3,3,2,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,2])*x[9,3,2,j]+RateStageTreatFirstR[2]*x[9,2,2,j]+RateStageTreatFirstL[3]*x[9,4,2,j]
      dx[10,3,2,j]=RateSuppFirstSusc[3]*x[9,3,2,j]-(RateFailFirstDTG[3]+RateStageSuppFirst[2]+mu[4,3,2])*x[10,3,2,j]+RateStageSuppFirst[3]*x[10,4,2,j]+RateSwitchDTG[2,3,2,j]*x[4,3,2,j]
      dx[11,3,2,j]=RateFailFirstDTG[3]*x[10,3,2,j]-(RateTreatSecond_noDTG[3,j]+RateStageFailFirst[3]+mu[5,3,2])*x[11,3,2,j]+RateStageFailFirst[2]*x[11,2,2,j]
      dx[12,3,2,j]=dx[12,3,2,j]+RateDiag[3,j]*((j==2)*(1-p_women_dtg))*x[1,3,2,j]-(RateTreatFirst_noDTG[2,3,j]+RateStageDiag[3]+mu[2,3,2])*x[12,3,2,j]+RateStageDiag[2]*x[12,2,2,j]
      dx[13,3,2,j]=dx[13,3,2,j]+RateTreatFirst_noDTG[2,3,j]*x[12,3,2,j]-(RateSuppFirstResis[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,2])*x[13,3,2,j]+RateStageTreatFirstR[2]*x[13,2,2,j]+RateStageTreatFirstL[3]*x[13,4,2,j]
      dx[14,3,2,j]=dx[14,3,2,j]+RateSuppFirstResis[3]*x[13,3,2,j]-(RateFailFirstResis[3]+RateStageSuppFirst[2]+mu[4,3,2])*x[14,3,2,j]+RateStageSuppFirst[3]*x[14,4,2,j]
      dx[15,3,2,j]=dx[15,3,2,j]+RateFailFirstResis[3]*x[14,3,2,j]-(RateTreatSecond_noDTG[3,j]+RateStageFailFirst[3]+mu[5,3,2])*x[15,3,2,j]+RateStageFailFirst[2]*x[15,2,2,j]
      
      dx[9,4,2,j]=dx[9,4,2,j]+RateTreatFirstDTG[2,4,j]*x[2,4,2,j]+RateSwitchDTG[3,4,2,j]*x[5,4,2,j]+RateSwitchDTG[1,4,2,j]*x[3,4,2,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,2])*x[9,4,2,j]+RateStageTreatFirstR[3]*x[9,3,2,j]
      dx[10,4,2,j]=RateSuppFirstSusc[4]*x[9,4,2,j]-(RateFailFirstDTG[4]+RateStageSuppFirst[3]+mu[4,4,2])*x[10,4,2,j]+RateSwitchDTG[2,4,2,j]*x[4,4,2,j]
      dx[11,4,2,j]=RateFailFirstDTG[4]*x[10,4,2,j]-(RateTreatSecond_noDTG[4,j]+mu[5,4,2])*x[11,4,2,j]+RateStageFailFirst[3]*x[11,3,2,j]
      dx[12,4,2,j]=dx[12,4,2,j]+RateDiag[4,j]*((j==2)*(1-p_women_dtg))*x[1,4,2,j]-(RateTreatFirst_noDTG[2,4,j]+mu[2,4,2])*x[12,4,2,j]+RateStageDiag[3]*x[12,3,2,j]
      dx[13,4,2,j]=dx[13,4,2,j]+RateTreatFirst_noDTG[2,4,j]*x[12,4,2,j]-(RateSuppFirstResis[4]+RateStageTreatFirstL[3]+mu[3,4,2])*x[13,4,2,j]+RateStageTreatFirstR[3]*x[13,3,2,j]
      dx[14,4,2,j]=dx[14,4,2,j]+RateSuppFirstResis[4]*x[13,4,2,j]-(RateFailFirstResis[4]+RateStageSuppFirst[3]+mu[4,4,2])*x[14,4,2,j]
      dx[15,4,2,j]=dx[15,4,2,j]+RateFailFirstResis[4]*x[14,4,2,j]-(RateTreatSecond_noDTG[4,j]+mu[5,4,2])*x[15,4,2,j]+RateStageFailFirst[3]*x[15,3,2,j]
      
      
      
      # #DTG compartments-------------------------------------
      # #NNRTI susceptible
      # dx[9,1,1,j]=dx[9,1,1,j]+RateTreatFirstDTG[1,1,j]*x[2,1,1,j]+RateSwitchDTG[1,1,1,j]*x[3,1,1,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,1])*x[9,1,1,j]+RateStageTreatFirstL[1]*x[9,2,1,j]
      # dx[10,1,1,j]=RateSwitchDTG[2,1,1,j]*x[4,1,1,j]+RateSuppFirstSusc[1]*x[9,1,1,j]-(RateFailFirstSusc[1]+mu[4,1,1])*x[10,1,1,j]+RateStageSuppFirst[1]*x[10,2,1,j]
      # dx[11,1,1,j]=RateSwitchDTG[3,1,1,j]*x[5,1,1,j]+RateFailFirstSusc[1]*x[10,1,1,j]-(RateStageFailFirst[1]+mu[5,1,1])*x[11,1,1,j]
      # 
      # dx[9,2,1,j]=dx[9,2,1,j]+RateTreatFirstDTG[1,2,j]*x[2,2,1,j]+RateSwitchDTG[1,2,1,j]*x[3,2,1,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,1])*x[9,2,1,j]+RateStageTreatFirstR[1]*x[9,1,1,j]+RateStageTreatFirstL[2]*x[9,3,1,j]
      # dx[10,2,1,j]=RateSwitchDTG[2,2,1,j]*x[4,2,1,j]+RateSuppFirstSusc[2]*x[9,2,1,j]-(RateFailFirstSusc[2]+RateStageSuppFirst[1]+mu[4,2,1])*x[10,2,1,j]+RateStageSuppFirst[2]*x[10,3,1,j]
      # dx[11,2,1,j]=RateSwitchDTG[3,2,1,j]*x[5,2,1,j]+RateFailFirstSusc[2]*x[10,2,1,j]-(RateStageFailFirst[2]+mu[5,2,1])*x[11,2,1,j]+RateStageFailFirst[1]*x[11,1,1,j]
      # 
      # dx[9,3,1,j]=dx[9,3,1,j]+RateTreatFirstDTG[1,3,j]*x[2,3,1,j]+RateSwitchDTG[1,3,1,j]*x[3,3,1,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,1])*x[9,3,1,j]+RateStageTreatFirstR[2]*x[9,2,1,j]+RateStageTreatFirstL[3]*x[9,4,1,j]
      # dx[10,3,1,j]=RateSwitchDTG[2,3,1,j]*x[4,3,1,j]+RateSuppFirstSusc[3]*x[9,3,1,j]-(RateFailFirstSusc[3]+RateStageSuppFirst[2]+mu[4,3,1])*x[10,3,1,j]+RateStageSuppFirst[3]*x[10,4,1,j]
      # dx[11,3,1,j]=RateSwitchDTG[3,3,1,j]*x[5,3,1,j]+RateFailFirstSusc[3]*x[10,3,1,j]-(RateStageFailFirst[3]+mu[5,3,1])*x[11,3,1,j]+RateStageFailFirst[2]*x[11,2,1,j]
      # 
      # dx[9,4,1,j]=dx[9,4,1,j]+RateTreatFirstDTG[1,4,j]*x[2,4,1,j]+RateSwitchDTG[1,4,1,j]*x[3,4,1,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,1])*x[9,4,1,j]+RateStageTreatFirstR[3]*x[9,3,1,j]
      # dx[10,4,1,j]=RateSwitchDTG[2,4,1,j]*x[4,4,1,j]+RateSuppFirstSusc[4]*x[9,4,1,j]-(RateFailFirstSusc[4]+RateStageSuppFirst[3]+mu[4,4,1])*x[10,4,1,j]
      # dx[11,4,1,j]=RateSwitchDTG[3,4,1,j]*x[5,4,1,j]+RateFailFirstSusc[4]*x[10,4,1,j]-(mu[5,4,1])*x[11,4,1,j]+RateStageFailFirst[3]*x[11,3,1,j]
      # ##############################################################################################################################################################
      
      dx[9,1,1,j]=dx[9,1,1,j]-RateStopTreatFirst[1]*x[9,1,1,j]-RateTreatToFailFirstDTG[1]*x[9,1,1,j]
      dx[10,1,1,j]=dx[10,1,1,j]-RateStopSuppFirst[1]*x[10,1,1,j]+RateFailToSuppTreatFirstSusc[1]*x[11,1,1,j]
      dx[11,1,1,j]=dx[11,1,1,j]-RateStopFailFirst[1]*x[11,1,1,j]-RateFailToSuppTreatFirstSusc[1]*x[11,1,1,j]+RateTreatToFailFirstDTG[1]*x[9,1,1,j]
      dx[12,1,1,j]=dx[12,1,1,j]-RateDirectTreatSecond[1,1]*x[12,1,1,j]+RateStopTreatFirst[1]*x[13,1,1,j]+RateStopSuppFirst[1]*x[14,1,1,j]+RateStopFailFirst[1]*x[15,1,1,j]+
        (RateStopTreatSecond[1]*x[6,1,1,j]+RateStopSuppSecond[1]*x[7,1,1,j]+RateStopFailSecond[1]*x[8,1,1,j])*((j==2)*(1-p_women_dtg))
      dx[13,1,1,j]=dx[13,1,1,j]-RateStopTreatFirst[1]*x[13,1,1,j]-RateTreatToFailFirstSusc[1]*x[13,1,1,j]
      dx[14,1,1,j]=dx[14,1,1,j]-RateStopSuppFirst[1]*x[14,1,1,j]+RateFailToSuppTreatFirstSusc[1]*x[15,1,1,j]-RateSuppFirstToSecond[1]*x[14,1,1,j]
      dx[15,1,1,j]=dx[15,1,1,j]-RateStopFailFirst[1]*x[15,1,1,j]-RateFailToSuppTreatFirstSusc[1]*x[15,1,1,j]+RateTreatToFailFirstSusc[1]*x[13,1,1,j]
      
      dx[9,2,1,j]=dx[9,2,1,j]-RateStopTreatFirst[2]*x[9,2,1,j]-RateTreatToFailFirstDTG[2]*x[9,2,1,j]
      dx[10,2,1,j]=dx[10,2,1,j]-RateStopSuppFirst[2]*x[10,2,1,j]+RateFailToSuppTreatFirstSusc[2]*x[11,2,1,j]-RateSuppFirstToSecond[2]*x[10,2,1,j]
      dx[11,2,1,j]=dx[11,2,1,j]-RateStopFailFirst[2]*x[11,2,1,j]-RateFailToSuppTreatFirstSusc[2]*x[11,2,1,j]+RateTreatToFailFirstDTG[2]*x[9,2,1,j]
      dx[12,2,1,j]=dx[12,2,1,j]-RateDirectTreatSecond[1,2]*x[12,2,1,j]+RateStopTreatFirst[2]*x[13,2,1,j]+RateStopSuppFirst[2]*x[14,2,1,j]+RateStopFailFirst[2]*x[15,2,1,j]+
        (RateStopTreatSecond[2]*x[6,2,1,j]+RateStopSuppSecond[2]*x[7,2,1,j]+RateStopFailSecond[2]*x[8,2,1,j])*((j==2)*(1-p_women_dtg))
      dx[13,2,1,j]=dx[13,2,1,j]-RateStopTreatFirst[2]*x[13,2,1,j]-RateTreatToFailFirstSusc[2]*x[13,2,1,j]
      dx[14,2,1,j]=dx[14,2,1,j]-RateStopSuppFirst[2]*x[14,2,1,j]+RateFailToSuppTreatFirstSusc[2]*x[15,2,1,j]-RateSuppFirstToSecond[2]*x[14,2,1,j]
      dx[15,2,1,j]=dx[15,2,1,j]-RateStopFailFirst[2]*x[15,2,1,j]-RateFailToSuppTreatFirstSusc[2]*x[15,2,1,j]+RateTreatToFailFirstSusc[2]*x[13,2,1,j]
      
      dx[9,3,1,j]=dx[9,3,1,j]-RateStopTreatFirst[3]*x[9,3,1,j]-RateTreatToFailFirstDTG[3]*x[9,3,1,j]
      dx[10,3,1,j]=dx[10,3,1,j]-RateStopSuppFirst[3]*x[10,3,1,j]+RateFailToSuppTreatFirstSusc[3]*x[11,3,1,j]
      dx[11,3,1,j]=dx[11,3,1,j]-RateStopFailFirst[3]*x[11,3,1,j]-RateFailToSuppTreatFirstSusc[3]*x[11,3,1,j]+RateTreatToFailFirstDTG[3]*x[9,3,1,j]
      dx[12,3,1,j]=dx[12,3,1,j]-RateDirectTreatSecond[1,3]*x[12,3,1,j]+RateStopTreatFirst[3]*x[13,3,1,j]+RateStopSuppFirst[3]*x[14,3,1,j]+RateStopFailFirst[3]*x[15,3,1,j]+
        (RateStopTreatSecond[3]*x[6,3,1,j]+RateStopSuppSecond[3]*x[7,3,1,j]+RateStopFailSecond[3]*x[8,3,1,j])*((j==2)*(1-p_women_dtg))
      dx[13,3,1,j]=dx[13,3,1,j]-RateStopTreatFirst[3]*x[13,3,1,j]-RateTreatToFailFirstSusc[3]*x[13,3,1,j]
      dx[14,3,1,j]=dx[14,3,1,j]-RateStopSuppFirst[3]*x[14,3,1,j]+RateFailToSuppTreatFirstSusc[3]*x[15,3,1,j]-RateSuppFirstToSecond[3]*x[14,3,1,j]
      dx[15,3,1,j]=dx[15,3,1,j]-RateStopFailFirst[3]*x[15,3,1,j]-RateFailToSuppTreatFirstSusc[3]*x[15,3,1,j]+RateTreatToFailFirstSusc[3]*x[13,3,1,j]
      
      dx[9,4,1,j]=dx[9,4,1,j]-RateStopTreatFirst[4]*x[9,4,1,j]-RateTreatToFailFirstDTG[4]*x[9,4,1,j]
      dx[10,4,1,j]=dx[10,4,1,j]-RateStopSuppFirst[4]*x[10,4,1,j]+RateFailToSuppTreatFirstSusc[4]*x[11,4,1,j]
      dx[11,4,1,j]=dx[11,4,1,j]-RateStopFailFirst[4]*x[11,4,1,j]-RateFailToSuppTreatFirstSusc[4]*x[11,4,1,j]+RateTreatToFailFirstDTG[4]*x[9,4,1,j]
      dx[12,4,1,j]=dx[12,4,1,j]-RateDirectTreatSecond[1,4]*x[12,4,1,j]+RateStopTreatFirst[4]*x[13,4,1,j]+RateStopSuppFirst[4]*x[14,4,1,j]+RateStopFailFirst[4]*x[15,4,1,j]+
        (RateStopTreatSecond[4]*x[6,4,1,j]+RateStopSuppSecond[4]*x[7,4,1,j]+RateStopFailSecond[4]*x[8,4,1,j])*((j==2)*(1-p_women_dtg))
      dx[13,4,1,j]=dx[13,4,1,j]-RateStopTreatFirst[4]*x[13,4,1,j]-RateTreatToFailFirstSusc[4]*x[13,4,1,j]
      dx[14,4,1,j]=dx[14,4,1,j]-RateStopSuppFirst[4]*x[14,4,1,j]+RateFailToSuppTreatFirstSusc[4]*x[15,4,1,j]-RateSuppFirstToSecond[4]*x[14,4,1,j]
      dx[15,4,1,j]=dx[15,4,1,j]-RateStopFailFirst[4]*x[15,4,1,j]-RateFailToSuppTreatFirstSusc[4]*x[15,4,1,j]+RateTreatToFailFirstSusc[4]*x[13,4,1,j]
      
      #############################################################################################################################################################################
      #############################################################################################################################################################################
      #NNRTI resistant
      # dx[9,1,2,j]=dx[9,1,2,j]+RateTreatFirstDTG[2,1,j]*x[2,1,2,j]+RateSwitchDTG[1,1,2,j]*x[3,1,2,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,2])*x[9,1,2,j]+RateStageTreatFirstL[1]*x[9,2,2,j]
      # dx[10,1,2,j]=RateSwitchDTG[2,1,2,j]*x[4,1,2,j]+RateSuppFirstSusc[1]*x[9,1,2,j]-(RateFailFirstSusc[1]+mu[4,1,2])*x[10,1,2,j]+RateStageSuppFirst[1]*x[10,2,2,j]
      # dx[11,1,2,j]=RateSwitchDTG[3,1,2,j]*x[5,1,2,j]+RateFailFirstSusc[1]*x[10,1,2,j]-(RateStageFailFirst[1]+mu[5,1,2])*x[11,1,2,j]
      # 
      # dx[9,2,2,j]=dx[9,2,2,j]+RateTreatFirstDTG[2,2,j]*x[2,2,2,j]+RateSwitchDTG[1,2,2,j]*x[3,2,2,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,2])*x[9,2,2,j]+RateStageTreatFirstR[1]*x[9,1,2,j]+RateStageTreatFirstL[2]*x[9,3,2,j]
      # dx[10,2,2,j]=RateSwitchDTG[2,2,2,j]*x[4,2,2,j]+RateSuppFirstSusc[2]*x[9,2,2,j]-(RateFailFirstSusc[2]+RateStageSuppFirst[1]+mu[4,2,2])*x[10,2,2,j]+RateStageSuppFirst[2]*x[10,3,2,j]
      # dx[11,2,2,j]=RateSwitchDTG[3,2,2,j]*x[5,2,2,j]+RateFailFirstSusc[2]*x[10,2,2,j]-(RateStageFailFirst[2]+mu[5,2,2])*x[11,2,2,j]+RateStageFailFirst[1]*x[11,1,2,j]
      # 
      # dx[9,3,2,j]=dx[9,3,2,j]+RateTreatFirstDTG[2,3,j]*x[2,3,2,j]+RateSwitchDTG[1,3,2,j]*x[3,3,2,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,2])*x[9,3,2,j]+RateStageTreatFirstR[2]*x[9,2,2,j]+RateStageTreatFirstL[3]*x[9,4,2,j]
      # dx[10,3,2,j]=RateSwitchDTG[2,3,2,j]*x[4,3,2,j]+RateSuppFirstSusc[3]*x[9,3,2,j]-(RateFailFirstSusc[3]+RateStageSuppFirst[2]+mu[4,3,2])*x[10,3,2,j]+RateStageSuppFirst[3]*x[10,4,2,j]
      # dx[11,3,2,j]=RateSwitchDTG[3,3,2,j]*x[5,3,2,j]+RateFailFirstSusc[3]*x[10,3,2,j]-(RateStageFailFirst[3]+mu[5,3,2])*x[11,3,2,j]+RateStageFailFirst[2]*x[11,2,2,j]
      # 
      # dx[9,4,2,j]=dx[9,4,2,j]+RateTreatFirstDTG[2,4,j]*x[2,4,2,j]+RateSwitchDTG[1,4,2,j]*x[3,4,2,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,2])*x[9,4,2,j]+RateStageTreatFirstR[3]*x[9,3,2,j]
      # dx[10,4,2,j]=RateSwitchDTG[2,4,2,j]*x[4,4,2,j]+RateSuppFirstSusc[4]*x[9,4,2,j]-(RateFailFirstSusc[4]+RateStageSuppFirst[3]+mu[4,4,2])*x[10,4,2,j]
      # dx[11,4,2,j]=RateSwitchDTG[3,4,2,j]*x[5,4,2,j]+RateFailFirstSusc[4]*x[10,4,2,j]-(mu[5,4,2])*x[11,4,2,j]+RateStageFailFirst[3]*x[11,3,2,j]
      ##############################################################################################################################################################
      
      dx[9,1,2,j]=dx[9,1,2,j]-RateStopTreatFirst[1]*x[9,1,2,j]-RateTreatToFailFirstDTG[1]*x[9,1,2,j]
      dx[10,1,2,j]=dx[10,1,2,j]-RateStopSuppFirst[1]*x[10,1,2,j]+RateFailToSuppTreatFirstSusc[1]*x[11,1,2,j]
      dx[11,1,2,j]=dx[11,1,2,j]-RateStopFailFirst[1]*x[11,1,2,j]-RateFailToSuppTreatFirstSusc[1]*x[11,1,2,j]+RateTreatToFailFirstDTG[1]*x[9,1,2,j]
      dx[12,1,2,j]=dx[12,1,2,j]-RateDirectTreatSecond[2,1]*x[12,1,2,j]+RateStopTreatFirst[1]*x[13,1,2,j]+RateStopSuppFirst[1]*x[14,1,2,j]+RateStopFailFirst[1]*x[15,1,2,j]+
        (RateStopTreatSecond[1]*x[6,1,2,j]+RateStopSuppSecond[1]*x[7,1,2,j]+RateStopFailSecond[1]*x[8,1,2,j])*((j==2)*(1-p_women_dtg))
      dx[13,1,2,j]=dx[13,1,2,j]-RateStopTreatFirst[1]*x[13,1,2,j]-RateTreatToFailFirstResis[1]*x[13,1,2,j]
      dx[14,1,2,j]=dx[14,1,2,j]-RateStopSuppFirst[1]*x[14,1,2,j]+RateFailToSuppTreatFirstResis[1]*x[15,1,2,j]-RateSuppFirstToSecond[1]*x[14,1,2,j]
      dx[15,1,2,j]=dx[15,1,2,j]-RateStopFailFirst[1]*x[15,1,2,j]-RateFailToSuppTreatFirstResis[1]*x[15,1,2,j]+RateTreatToFailFirstResis[1]*x[13,1,2,j]
      
      dx[9,2,2,j]=dx[9,2,2,j]-RateStopTreatFirst[2]*x[9,2,2,j]-RateTreatToFailFirstDTG[2]*x[9,2,2,j]
      dx[10,2,2,j]=dx[10,2,2,j]-RateStopSuppFirst[2]*x[10,2,2,j]+RateFailToSuppTreatFirstSusc[2]*x[11,2,2,j]
      dx[11,2,2,j]=dx[11,2,2,j]-RateStopFailFirst[2]*x[11,2,2,j]-RateFailToSuppTreatFirstSusc[2]*x[11,2,2,j]+RateTreatToFailFirstDTG[2]*x[9,2,2,j]
      dx[12,2,2,j]=dx[12,2,2,j]-RateDirectTreatSecond[2,2]*x[12,2,2,j]+RateStopTreatFirst[2]*x[13,2,2,j]+RateStopSuppFirst[2]*x[14,2,2,j]+RateStopFailFirst[2]*x[15,2,2,j]+
        (RateStopTreatSecond[2]*x[6,2,2,j]+RateStopSuppSecond[2]*x[7,2,2,j]+RateStopFailSecond[2]*x[8,2,2,j])*((j==2)*(1-p_women_dtg))
      dx[13,2,2,j]=dx[13,2,2,j]-RateStopTreatFirst[2]*x[13,2,2,j]-RateTreatToFailFirstResis[2]*x[13,2,2,j]
      dx[14,2,2,j]=dx[14,2,2,j]-RateStopSuppFirst[2]*x[14,2,2,j]+RateFailToSuppTreatFirstResis[2]*x[15,2,2,j]-RateSuppFirstToSecond[2]*x[14,2,2,j]
      dx[15,2,2,j]=dx[15,2,2,j]-RateStopFailFirst[2]*x[15,2,2,j]-RateFailToSuppTreatFirstResis[2]*x[15,2,2,j]+RateTreatToFailFirstResis[2]*x[13,2,2,j]
      
      dx[9,3,2,j]=dx[9,3,2,j]-RateStopTreatFirst[3]*x[9,3,2,j]-RateTreatToFailFirstDTG[3]*x[9,3,2,j]
      dx[10,3,2,j]=dx[10,3,2,j]-RateStopSuppFirst[3]*x[10,3,2,j]+RateFailToSuppTreatFirstSusc[3]*x[11,3,2,j]
      dx[11,3,2,j]=dx[11,3,2,j]-RateStopFailFirst[3]*x[11,3,2,j]-RateFailToSuppTreatFirstSusc[3]*x[11,3,2,j]+RateTreatToFailFirstDTG[3]*x[9,3,2,j]
      dx[12,3,2,j]=dx[12,3,2,j]-RateDirectTreatSecond[2,3]*x[12,3,2,j]+RateStopTreatFirst[3]*x[13,3,2,j]+RateStopSuppFirst[3]*x[14,3,2,j]+RateStopFailFirst[3]*x[15,3,2,j]+
        (RateStopTreatSecond[3]*x[6,3,2,j]+RateStopSuppSecond[3]*x[7,3,2,j]+RateStopFailSecond[3]*x[8,3,2,j])*((j==2)*(1-p_women_dtg))
      dx[13,3,2,j]=dx[13,3,2,j]-RateStopTreatFirst[3]*x[13,3,2,j]-RateTreatToFailFirstResis[3]*x[13,3,2,j]
      dx[14,3,2,j]=dx[14,3,2,j]-RateStopSuppFirst[3]*x[14,3,2,j]+RateFailToSuppTreatFirstResis[3]*x[15,3,2,j]-RateSuppFirstToSecond[3]*x[14,3,2,j]
      dx[15,3,2,j]=dx[15,3,2,j]-RateStopFailFirst[3]*x[15,3,2,j]-RateFailToSuppTreatFirstResis[3]*x[15,3,2,j]+RateTreatToFailFirstResis[3]*x[13,3,2,j]
      
      dx[9,4,2,j]=dx[9,4,2,j]-RateStopTreatFirst[4]*x[9,4,2,j]-RateTreatToFailFirstDTG[4]*x[9,4,2,j]
      dx[10,4,2,j]=dx[10,4,2,j]-RateStopSuppFirst[4]*x[10,4,2,j]+RateFailToSuppTreatFirstSusc[4]*x[11,4,2,j]
      dx[11,4,2,j]=dx[11,4,2,j]-RateStopFailFirst[4]*x[11,4,2,j]-RateFailToSuppTreatFirstSusc[4]*x[11,4,2,j]+RateTreatToFailFirstDTG[4]*x[9,4,2,j]
      dx[12,4,2,j]=dx[12,4,2,j]-RateDirectTreatSecond[2,4]*x[12,4,2,j]+RateStopTreatFirst[4]*x[13,4,2,j]+RateStopSuppFirst[4]*x[14,4,2,j]+RateStopFailFirst[4]*x[15,4,2,j]+
        (RateStopTreatSecond[4]*x[6,4,2,j]+RateStopSuppSecond[4]*x[7,4,2,j]+RateStopFailSecond[4]*x[8,4,2,j])*((j==2)*(1-p_women_dtg))
      dx[13,4,2,j]=dx[13,4,2,j]-RateStopTreatFirst[4]*x[13,4,2,j]-RateTreatToFailFirstResis[4]*x[13,4,2,j]
      dx[14,4,2,j]=dx[14,4,2,j]-RateStopSuppFirst[4]*x[14,4,2,j]+RateFailToSuppTreatFirstResis[4]*x[15,4,2,j]-RateSuppFirstToSecond[4]*x[14,4,2,j]
      dx[15,4,2,j]=dx[15,4,2,j]-RateStopFailFirst[4]*x[15,4,2,j]-RateFailToSuppTreatFirstResis[4]*x[15,4,2,j]+RateTreatToFailFirstResis[4]*x[13,4,2,j]
      
      
      
      #Resistance acquisition and reversion
      dx[12,1,2,j]=dx[12,1,2,j]-RateSusceptible*x[12,1,2,j]
      dx[15,1,2,j]=dx[15,1,2,j]+RateResistant*x[15,1,1,j]
      dx[12,2,2,j]=dx[12,2,2,j]-RateSusceptible*x[12,2,2,j]
      dx[15,2,2,j]=dx[15,2,2,j]+RateResistant*x[15,2,1,j]
      dx[12,3,2,j]=dx[12,3,2,j]-RateSusceptible*x[12,3,2,j]
      dx[15,3,2,j]=dx[15,3,2,j]+RateResistant*x[15,3,1,j]
      dx[12,4,2,j]=dx[12,4,2,j]-RateSusceptible*x[12,4,2,j]
      dx[15,4,2,j]=dx[15,4,2,j]+RateResistant*x[15,4,1,j]
      
      dx[12,1,1,j]=dx[12,1,1,j]+RateSusceptible*x[12,1,2,j]
      dx[15,1,1,j]=dx[15,1,1,j]-RateResistant*x[15,1,1,j]
      dx[12,2,1,j]=dx[12,2,1,j]+RateSusceptible*x[12,2,2,j]
      dx[15,2,1,j]=dx[15,2,1,j]-RateResistant*x[15,2,1,j]
      dx[12,3,1,j]=dx[12,3,1,j]+RateSusceptible*x[12,3,2,j]
      dx[15,3,1,j]=dx[15,3,1,j]-RateResistant*x[15,3,1,j]
      dx[12,4,1,j]=dx[12,4,1,j]+RateSusceptible*x[12,4,2,j]
      dx[15,4,1,j]=dx[15,4,1,j]-RateResistant*x[15,4,1,j]
      
    }
    dx[x+dx<0]=0
    #dx[12:15,1:4,1:2,1]=0
    #new_infections
    new_inf=S[1]/N[1]*(sum(RateInf1[1:4,1:2,1]*apply(x[c(1),1:4,1:2,1:2],c(1,3),sum))+
                         sum(RateInf2[1:4,1:2,1]*apply(x[c(2,3,5,6,8,9,11,12,13,15),1:4,1:2,1:2],c(2,4),sum)))+
      S[2]/N[2]*(sum(RateInf1[1:4,1:2,2]*apply(x[c(1),1:4,1:2,1:2],c(1,3),sum))+
                   sum(RateInf2[1:4,1:2,2]*apply(x[c(2,3,5,6,8,9,11,12,13,15),1:4,1:2,1:2],c(2,4),sum)))
    #death
    death=sum(apply(x[1:8,1:4,1:2,1:2],c(1,2,3),sum)*mu)+sum(apply(x[12,1:4,1:2,1:2],c(1,2),sum)*mu[2,1:4,1:2])+
      sum(apply(x[9:11,1:4,1:2,1:2],c(1,2,3),sum)*mu[3:5,1:4,1:2])+sum(apply(x[13:15,1:4,1:2,1:2],c(1,2,3),sum)*mu[3:5,1:4,1:2])
    death_w=sum(x[1:8,1:4,1:2,2]*mu)+sum(x[12,1:4,1:2,2]*mu[2,1:4,1:2])+
      sum(x[9:11,1:4,1:2,2]*mu[3:5,1:4,1:2])+sum(x[13:15,1:4,1:2,2]*mu[3:5,1:4,1:2])
    death_w_nnrti=sum(x[3:5,1:4,1:2,2]*mu[3:5,1:4,1:2])+sum(x[13:15,1:4,1:2,2]*mu[3:5,1:4,1:2])
    death_w_inel=sum(x[13:15,1:4,1:2,2]*mu[3:5,1:4,1:2])
    death2=sum(apply(x[3:5,1:4,1:2,1:2],c(1,2,3),sum)*mu[3:5,1:4,1:2])
    death3=sum(apply(x[13:15,1:4,1:2,1:2],c(1,2,3),sum)*mu[3:5,1:4,1:2])
    death4=sum(apply(x[9:11,1:4,1:2,1:2],c(1,2,3),sum)*mu[3:5,1:4,1:2])
    death5=sum(apply(x[2,1:4,1:2,1:2],c(1,2),sum)*mu[2,1:4,1:2])
    #new_infections resistant
    new_inf_res=S[1]/N[1]*(sum(RateInf1[1:4,1:2,1]*apply(x[c(1),1:4,2,1:2],c(1,2),sum))+
                             sum(RateInf2[1:4,1:2,1]*apply(x[c(2,3,5,6,8,9,11,12,13,15),1:4,2,1:2],c(2,3),sum)))+
      S[2]/N[2]*(sum(RateInf1[1:4,1:2,2]*apply(x[c(1),1:4,2,1:2],c(1,2),sum))+
                   sum(RateInf2[1:4,1:2,2]*apply(x[c(2,3,5,6,8,9,11,12,13,15),1:4,2,1:2],c(2,3),sum)))
    
    n1=S[1]/N[1]*sum(RateInf1[1:4,1:2,1]*apply(x[c(1),1:4,2,1:2],c(1,2),sum))#susc inf
    n2=S[1]/N[1]*sum(RateInf2[1:4,1:2,1]*apply(x[c(2,3,5,6,8,9,11),1:4,2,1:2],c(2,3),sum))#susc diag
    n3=S[2]/N[2]*sum(RateInf1[1:4,1:2,2]*apply(x[c(1),1:4,2,1:2],c(1,2),sum))#res inf
    n4=S[2]/N[2]*sum(RateInf2[1:4,1:2,2]*apply(x[c(2,3,5,6,8,9,11),1:4,2,1:2],c(2,3),sum))#res diag
    #death resistant
    death_res=sum(apply(x[1:8,1:4,2,1:2],c(1,2),sum)*mu[1:8,1:4,2])+sum(apply(x[12,1:4,2,1:2],c(1),sum)*mu[2,1:4,2])+
      sum(apply(x[9:11,1:4,2,1:2],c(1,2),sum)*mu[3:5,1:4,2])+sum(apply(x[13:15,1:4,2,1:2],c(1,2),sum)*mu[3:5,1:4,2])
    
    
    death_treat=sum(apply(x[3:8,1:4,1:2,1:2],c(1,2,3),sum)*mu[3:8,1:4,1:2])+
      sum(apply(x[9:11,1:4,1:2,1:2],c(1,2,3),sum)*mu[3:5,1:4,1:2])+sum(apply(x[13:15,1:4,1:2,1:2],c(1,2,3),sum)*mu[3:5,1:4,1:2])
    death_treat_res=sum(apply(x[3:8,1:4,2,1:2],c(1,2),sum)*mu[3:8,1:4,2])+
      sum(apply(x[9:11,1:4,2,1:2],c(1,2),sum)*mu[3:5,1:4,2])+sum(apply(x[13:15,1:4,2,1:2],c(1,2),sum)*mu[3:5,1:4,2])
    
    
    new_treat=sum(RateTreatFirst[1:2,1:4,1:2]*aperm(x[2,1:4,1:2,1:2],perm=c(2,1,3)))+
      sum(RateTreatFirstDTG[1:2,1:4,1:2]*aperm(x[2,1:4,1:2,1:2],perm=c(2,1,3)))+
      sum(RateDirectTreatSecond[1:2,1:4]*aperm(apply(x[2,1:4,1:2,1:2],c(1,2),sum),c(2,1)))+
      sum(RateTreatFirst[1:2,1:4,1:2]*aperm(x[12,1:4,1:2,1:2],perm=c(2,1,3)))+
      sum(RateDirectTreatSecond[1:2,1:4]*aperm(apply(x[12,1:4,1:2,1:2],c(1,2),sum),c(2,1)))
    
    
    new_treat_res=sum(RateTreatFirst[2,1:4,1:2]*x[2,1:4,2,1:2])+
      sum(RateTreatFirstDTG[2,1:4,1:2]*x[2,1:4,2,1:2])+
      sum(RateDirectTreatSecond[2,1:4]*apply(x[2,1:4,2,1:2],c(1),sum))+
      sum(RateTreatFirst[2,1:4,1:2]*x[12,1:4,2,1:2])+
      sum(RateDirectTreatSecond[2,1:4]*apply(x[12,1:4,2,1:2],c(1),sum))
    
    death_diag=sum(apply(x[12,1:4,1:2,1:2],c(1,2),sum)*mu[2,1:4,1:2])+sum(apply(x[2,1:4,1:2,1:2],c(1,2),sum)*mu[2,1:4,1:2])
    death_diag_res=sum(apply(x[12,1:4,2,1:2],c(1),sum)*mu[2,1:4,2])+sum(apply(x[2,1:4,2,1:2],c(1),sum)*mu[2,1:4,2])
    
    new_diag=sum(RateDiag[1:4,1:2]*apply(x[1,1:4,1:2,1:2],c(1,3),sum))
    new_diag_res=sum(RateDiag[1:4,1:2]*apply(x[1,1:4,2,1:2],c(1,2),sum))
    
    
    #why last cd4 class go down quicker in res than in susc
    new_pi_susc=sum(apply(RateDirectTreatSecond[1:2,1:4],c(2),sum)*apply(x[c(2,12),1:4,1,1:2],c(2),sum))+
      sum(RateTreatSecond[1:4,1:2]*x[c(5),1:4,1,1:2])+
      sum(RateTreatSecond_noDTG[1:4,1:2]*apply(x[c(5,15),1:4,1,1:2],c(2,3),sum))
    new_pi_res=sum(apply(RateDirectTreatSecond[1:2,1:4],c(2),sum)*apply(x[c(2,12),1:4,2,1:2],c(2),sum))+
      sum(RateTreatSecond[1:4,1:2]*x[c(5),1:4,2,1:2])+
      sum(RateTreatSecond_noDTG[1:4,1:2]*apply(x[c(5,15),1:4,2,1:2],c(2,3),sum))
    new_pi_susc_4=sum(RateDirectTreatSecond[1:2,4]*apply(x[c(2,12),4,1,1:2],c(2),sum))+
      sum(RateTreatSecond[4,1:2]*x[c(5),4,1,1:2])+
      sum(RateTreatSecond_noDTG[4,1:2]*apply(x[c(5,15),4,1,1:2],c(2),sum))
    new_pi_res_4=sum(RateDirectTreatSecond[1:2,4]*apply(x[c(2,12),4,2,1:2],c(2),sum))+
      sum(RateTreatSecond[4,1:2]*x[c(5),4,2,1:2])+
      sum(RateTreatSecond_noDTG[4,1:2]*apply(x[c(5,15),4,2,1:2],c(2),sum))
    
    
    m_susc=sum(RateDiag[1:4,1]*x[1,1:4,1,1]+RateDiag[1:4,2]*p_women_dtg*x[1,1:4,1,2])
    m_res=sum(RateDiag[1:4,1]*x[1,1:4,2,1]+RateDiag[1:4,2]*p_women_dtg*x[1,1:4,2,2])
    
    m_susc=m_susc-sum(RateTreatFirst[1,1:4,1]*x[2,1:4,1,1]+RateTreatFirst[1,1:4,2]*x[2,1:4,1,2])
    m_res=m_res-sum(RateTreatFirst[2,1:4,1]*x[2,1:4,2,1]+RateTreatFirst[2,1:4,2]*x[2,1:4,2,2])
    
    m_susc=m_susc-sum(RateDirectTreatSecond[1,1:4]*x[2,1:4,1,1]+RateDirectTreatSecond[1,1:4]*x[2,1:4,1,2])
    m_res=m_res-sum(RateDirectTreatSecond[2,1:4]*x[2,1:4,2,1]+RateDirectTreatSecond[2,1:4]*x[2,1:4,2,2])
    
    m_susc=m_susc-sum(RateTreatFirstDTG[1,1:4,1]*x[2,1:4,1,1]+RateTreatFirstDTG[1,1:4,2]*x[2,1:4,1,2])
    m_res=m_res-sum(RateTreatFirstDTG[2,1:4,1]*x[2,1:4,2,1]+RateTreatFirstDTG[2,1:4,2]*x[2,1:4,2,2])
    
    m_susc=m_susc-sum(mu[2,1:4,1]*x[2,1:4,1,1]+mu[2,1:4,1]*x[2,1:4,1,2])
    m_res=m_res-sum(mu[2,1:4,2]*x[2,1:4,2,1]+mu[2,1:4,2]*x[2,1:4,2,2])
    
    # new_treat_susc=sum(RateTreatFirstDTG[1,1:4,1:2]*x[2,1:4,1,1:2])
    # new_treat_res=sum(RateTreatFirstDTG[2,1:4,1:2]*x[2,1:4,2,1:2])
    
    
    mort_susc=sum(mu[1,1:4,1]*apply(x[1,1:4,1,1:2],1,sum))
    mort_res=sum(mu[1,1:4,2]*apply(x[1,1:4,2,1:2],1,sum))
    
    # #new_diag_susc=sum(x[2,1:4,1,1:2])
    # #new_diag_res=sum(x[2,1:4,2,1:2])
    # mort_susc=sum(mu[2,1:4,1]*apply(x[2,1:4,1,1:2],1,sum))
    # mort_res=sum(mu[2,1:4,2]*apply(x[2,1:4,2,1:2],1,sum))
    
    #dx[3:8,1:4,1,1:2]=dx[3:8,1:4,1,1:2]+infected_m15_month[t+1]*x[3:8,1:4,1,1:2]/sum(x[3:8,1:4,1,1:2])
    ####################################################################################################################
    res=c(dx,new_inf,death,new_diag,new_treat,death_treat,new_inf_res,death_res,death_w,death_w_nnrti,death_w_inel)
    #res=c(dx,new_inf,death,new_diag,new_treat,death_treat,new_inf_res,death_res,new_diag_res,new_treat_res,death_treat_res)
    #res=c(dx,new_inf,death,new_diag,new_pi_susc,new_pi_res,new_inf_res,death_res,new_diag_res,new_pi_susc_4,new_pi_res_4)
    #res=c(dx,new_treat_susc,new_treat_res,new_diag_susc,new_diag_res,mort_susc,mort_res)
    #res=c(dx,new_inf,death,new_diag,new_treat,death_diag,new_inf_res,death_res,new_diag_res,new_treat_res,death_diag_res)
    return(c(res))
    #return(list(c(res,newinf,newsu,newres,death)))
    #})
  })
}

###########################################################################################
#DTG
#Rates
#rates are defined in a separate function, unlike before
mod_rate=function(t,x,p1,treat_dtg,parms,p2){
  
  e=as.list(parms)
  e<-within(as.list(e), {
    rate1_inf <- p1[1] #number of unprotected sexual acts/month
    rate2_inf <- p1[2]
    rate1_diag <- p1[3] #diagnosis rate
    rate2_diag <- p1[4] #diagnosis rate
    rate_treat <- p1[5] #scale parameter for overall treatment rate
    rate_death <- p1[6] #scale parameter for overall death
    q <- p2[7]
    #Range
    rate_ratio <- p2[8] #varying parameter to estimate death for treated (but neither suppressed nor failed)
    rate_ratio <- 0.9
    
    alpha1 <- p2[12] #ratio between outcome when resistant/sensitive
    rate_res <- p2[13] #resistant rate
    rate_susc <- p2[14] #reversion rate
    alpha2 <- p2[12]
    
    #1. Infection and Diagnosis
    RateInf1<-array(rep(RateTransm,each=4),dim=c(4,2,2))*rate1_inf
    RateInf2<-RateInf1*rate2_inf
    
    Diag_frame[,4]<-Diag_frame[,4]*rate2_diag
    RateDiag<-1/rate1_diag*apply(Diag_frame,1,function(x) smooth_st(2005+t/12,x["a"],x["b"],x["c"],x["d"]))
    
    RateDiag<-array(rep(RateDiag,2)*rep(c(1,diag_women),each=4),dim=c(4,2))+
      array(rep(oi_inc,2)*rep(smooth_st(2005+t/12,oi_test["a"],oi_test["b"],oi_test["c"],oi_test["d"]),8),dim=c(4,2))+
      array(c(rep(0,4),preg_inc)*rep(smooth_st(2005+t/12,preg_test["a"],preg_test["b"],preg_test["c"],preg_test["d"]),8),dim=c(4,2))
    
    #2. Treatment theoretical (if no counterfactual scenario)
    #2.1 Rate to first- and second-line regimen, NNRTI and PI
    
    #not more used
    #RateTreatFirst_th<-apply(Treat_frame,1,function(x) smooth_st(2005+t/12,x["a"],x["b"],x["c"],x["d"]))*smooth_st(2005+t/12,2001,2012,1,200/12)*rate_treat
    
    #Rate treat first function, second part is only used with scenario of earlier treat-all policy (take the max to have continuous increase in each cd4 stata)
    RateTreatFirst_f<-function(y){
      if(y<2005){
        return(apply(Treat_frame_original,1,function(x) smooth_st(y,x["a"],x["b"],x["c"],x["d"]))*smooth_st(y,2001,2012,1,200/12)*rate_treat)
      }else{
        return(apply(
          data.frame(t1=apply(Treat_frame,1,function(x) lin_st(y,x["a"],x["b"],x["c"],x["d"])),
                     t2=apply(Treat_frame2,1,function(x) lin_st(y,x["a"],x["b"],x["c"],x["d"]))*as.numeric(y>Treat_frame2$a[4])),
          1,max)*smooth_st(y,2001,2012,1,200/12)*rate_treat)
      }
    }
    #Theoretical treatment rates for first- and second-line regimen
    RateTreatFirst_th<-RateTreatFirst_f(2005+t/12)
    RateDirectTreatSecond_th<-RateDirectTreatSecond
    
    #RateTreatFirst[res,cd4]
    RateTreatFirst<-matrix(c(RateTreatFirst_th,RateTreatFirst_th),nrow=2,ncol=4,byrow=T)
    #RateTreatSecond[res,cd4]
    RateDirectTreatSecond<-matrix(c(RateDirectTreatSecond_th,RateDirectTreatSecond_th),nrow=2,ncol=4,byrow=T)
    
    #2.2 Treatment rates in case of DTG introduction
    #percentage of women eligible for dtg
    if(treat_dtg[1]!=treat_dtg[2] & treat_dtg[2]!=0){print("percentage of women eligible for dtg not equal for FirstTreat and Switch.")}
    p_women_dtg<-treat_dtg[1]
    #RateTreatFirstDTG[res,cd4,gender]
    RateTreatFirstDTG<-array(rep(RateTreatFirst*treat_dtg[3],2)*rep(c(1,1),each=8),dim=c(2,4,2))*
      smooth_st(2005+t/12,2018,2018+treat_dtg["delay_first"],0,1) #delay
    #-------------------------------------------------------------------------------------------------------------
    #2.2.1 RateSwitchDTG[care,cd4,res,gender]
    #1) All people switch at a rate equal to RateTreatSecond
    #RateSwitchDTG<-array(rep(rep(treat_dtg[4]*RateTreatSecond,each=3)*rep(c(treat_dtg[5],treat_dtg[6],treat_dtg[7]),4),4),dim=c(3,4,2,2))*
    #smooth_st(2005+t/12,2018,2018+treat_dtg["delay_switch"],0,1) #delay
    RateSwitchDTG<-array(rep(rep(RateTreatSecond,each=3)*rep(c(treat_dtg[5],treat_dtg[6],treat_dtg[7]),4),4),dim=c(3,4,2,2))*
      smooth_st(2005+t/12,2018,2018+treat_dtg["delay_switch"],0,1) #delay
    #2) All people switch at a rate equal to 1/12 months^(-1)
    # RateSwitchDTG=array(rep(rep(treat_dtg[4]*rep(1/12,4),each=3)*rep(c(treat_dtg[5],treat_dtg[6],treat_dtg[7]),4),4),dim=c(3,4,2,2))*
    #    smooth_st(2005+t/12,2018,2018+treat_dtg["delay_switch"],0,1) #delay
    #3) Treated and suppressed people switch at a rate equal to 1year^-1, failed people at a rate still equal to 1) RateTreatSecond or 2) 1/12 months^(-1)
    RateSwitchDTG[1:2,1:4,1:2,1:2]<-rep(treat_dtg[4],32)*rep(rep(1/12,2)*c(treat_dtg[5],treat_dtg[6]),16)*
      smooth_st(2005+t/12,2018,2018+treat_dtg["delay_switch"],0,1) #delay
    #-------------------------------------------------------------------------------------------------------------
    #2.2.2 Modified treatment rates
    #RateTreatFirst [res,cd4,gender]
    RateTreatFirst_noDTG<-array(rep(RateTreatFirst,2),dim=c(2,4,2)) #Treatment rate to NNRTI for DTG-ineligible people
    RateTreatFirst<-apply(RateTreatFirst_noDTG-RateTreatFirstDTG,c(1,2,3),function(x) max(x,0))  #Treatment rate to NNRTI for DTG-eligible people
    
    #RateTreatSecond_noDTG[CD4,gender]
    #1) Treat rates to second for dtg-ineligible people stay the same
    RateTreatSecond_noDTG<-array(rep(RateTreatSecond,2),dim=c(4,2)) #Treatment rate to PI for DTG-ineligible people
    #2) Treat rates to second for dtg-ineligible people increase in 2018 (only do that when switching rate to dtg of failed people is 1/12, in order to have same switching rates for both dtg-ineligible and -eligible)
    # if(2005+t/12<2018){RateTreatSecond_noDTG=array(rep(RateTreatSecond,2),dim=c(4,2))
    # }else{RateTreatSecond_noDTG=array(rep(1/12,8),dim=c(4,2))}
    #RateTreatSecond[CD4,gender]
    RateTreatSecond<-apply(RateTreatSecond_noDTG-apply(RateSwitchDTG[3,1:4,1:2,1:2],c(1,3),mean),c(1,2),function(x) max(x,0)) #Treatment rate to PI for DTG-eligible people
    #we take the mean here as RateTreatSecond is the same for men and women. It has not any impact as long as RateSwitchDTG is also the same for men and women.
    ############
    #print(RateTreatSecond_noDTG)
    #print(RateTreatSecond)
    ############
    #3. Death
    p1<-min(1,q*mean(sapply(2005+seq(t-36,t,1)/12,function(x) RateTreatFirst_f(x)[4]/rate_treat)))
    p1=1-0.27*exp(-q*(mean(sapply(2005+seq(t-36,t,1)/12,function(x) RateTreatFirst_f(x)[4])-
                             sapply(2005+seq(-36,0,1)/12,function(x) RateTreatFirst_f(x)[4])))/rate_treat)
    mu[1:2,4,1]<-p1*40.9+(1-p1)*134.4
    mu[4,4,1]<-p1*8.3+(1-p1)*41.7
    mu[5,4,1]<-p1*11.8+(1-p1)*59.7
    mu[3,1:4,1]<-rate_ratio*mu[4,1:4,1]+(1-rate_ratio)*mu[5,1:4,1]
    mu[6:8,1:4,1]<-mu[3:5,1:4,1]
    mu[1:8,1:4,2]<-mu[1:8,1:4,1]
    mu<-mu*rate_death/1000
    
    #Back and forth mutation and impact on treatment outcome
    RateResistant<-1/rate_res
    RateSusceptible<-1/rate_susc
    #From T to S
    RateSuppFirstSusc<-RateSuppFirst
    RateSuppFirstResis<-1/alpha1*RateSuppFirst
    #From S to F
    RateFailFirstSusc<-RateFailFirst
    RateFailFirstResis<-alpha2*RateFailFirst
    #From T to F
    RateTreatToFailFirstSusc<-RateTreatToFailFirst
    RateTreatToFailFirstResis<-alpha1*RateTreatToFailFirst
    #From F to S
    RateFailToSuppTreatFirstSusc<-RateFailToSuppTreatFirst
    RateFailToSuppTreatFirstResis<-1/alpha2*RateFailToSuppTreatFirst
    
    #For second-line regimen, no effect of NNRTI resistance on treatment as it is PI
    ##################################################################################################################
    if(length(x)==250){
      I<-apply(array(x,dim=c(15,4,2,2)),4,sum)
    }else if(length(x)==298){
      I<-apply(array(x,dim=c(18,4,2,2)),4,sum)
    }else{
      stop("length of x !=250 or 298 (my comment)")
    }
    
    #First approach to model S: S=N-I
    N<-N_value_f(t+1)
    S<-N-I
    #Percentage getting NNRTI-based 1st-line[gender]
    p_NNRTI<-sum(RateTreatFirst)/sum(RateTreatFirst+RateTreatFirstDTG)*c(1,treat_dtg[1])+c(0,1-treat_dtg[1])
  })
  return(e)
  #return(list(e[["RateSwitchDTG"]],e[["RateTreatSecond"]]))
}
mod_rate_2res=function(t,x,p1,treat_dtg,parms,p2){
  
  e=as.list(parms)
  e<-within(as.list(e), {
    rate1_inf <- p1[1] #number of unprotected sexual acts/month
    rate2_inf <- p1[2]
    rate1_diag <- p1[3] #diagnosis rate
    rate2_diag <- p1[4] #diagnosis rate
    rate_treat <- p1[5] #scale parameter for overall treatment rate
    rate_death <- p1[6] #scale parameter for overall death
    q <- p2[7]
    #Range
    rate_ratio <- p2[8] #varying parameter to estimate death for treated (but neither suppressed nor failed)
    rate_ratio <- 0.9
    
    alpha1 <- p2[12] #ratio between outcome when resistant/sensitive
    rate_res <- p2[13] #resistant rate
    rate_susc <- p2[14] #reversion rate
    alpha2 <- p2[15]
    
    #1. Infection and Diagnosis
    RateInf1<-array(rep(RateTransm,each=4),dim=c(4,2,2))*rate1_inf
    RateInf2<-RateInf1*rate2_inf
    
    Diag_frame[,4]<-Diag_frame[,4]*rate2_diag
    RateDiag<-1/rate1_diag*apply(Diag_frame,1,function(x) smooth_st(2005+t/12,x["a"],x["b"],x["c"],x["d"]))
    
    RateDiag<-array(rep(RateDiag,2)*rep(c(1,diag_women),each=4),dim=c(4,2))+
      array(rep(oi_inc,2)*rep(smooth_st(2005+t/12,oi_test["a"],oi_test["b"],oi_test["c"],oi_test["d"]),8),dim=c(4,2))+
      array(c(rep(0,4),preg_inc)*rep(smooth_st(2005+t/12,preg_test["a"],preg_test["b"],preg_test["c"],preg_test["d"]),8),dim=c(4,2))
    
    #2. Treatment theoretical (if no counterfactual scenario)
    #2.1 Rate to first- and second-line regimen, NNRTI and PI
    
    #not more used
    #RateTreatFirst_th<-apply(Treat_frame,1,function(x) smooth_st(2005+t/12,x["a"],x["b"],x["c"],x["d"]))*smooth_st(2005+t/12,2001,2012,1,200/12)*rate_treat
    
    #Rate treat first function, second part is only used with scenario of earlier treat-all policy (take the max to have continuous increase in each cd4 stata)
    RateTreatFirst_f<-function(y){
      #apply(TreatAll_frame,1,function(x) lin_st(y,x["a"],x["b"],x["c"],x["d"]*(1/12)/(rate_treat*200/12)))*
      apply(TreatAll_frame,1,function(x) lin_st(y,x["a"],x["b"],x["c"],x["d"]))* #Treat-All
        apply(Treat_frame,1,function(x) lin_st(y,x["a"],x["b"],x["c"],x["d"]))* #before treat-all, higher rate for lower cd4 counts
        smooth_st(y,2001,2012,1,200/12)*rate_treat #increase in time and baseline rate in 2001
    }
    #Theoretical treatment rates for first- and second-line regimen
    RateTreatFirst_th<-RateTreatFirst_f(2005+t/12)
    RateDirectTreatSecond_th<-RateDirectTreatSecond
    #RateTreatFirst[res,cd4]
    RateTreatFirst<-matrix(c(RateTreatFirst_th,RateTreatFirst_th),nrow=2,ncol=4,byrow=T)
    
    #RateTreatSecond[res,cd4]
    RateDirectTreatSecond<-matrix(c(RateDirectTreatSecond_th,RateDirectTreatSecond_th),nrow=2,ncol=4,byrow=T)
    
    #2.2 Treatment rates in case of DTG introduction
    #percentage of women eligible for dtg
    if(treat_dtg[1]!=treat_dtg[2] & treat_dtg[2]!=0){print("percentage of women eligible for dtg not equal for FirstTreat and Switch.")}
    p_women_dtg<-treat_dtg[1]
    #RateTreatFirstDTG[res,cd4,gender]
    RateTreatFirstDTG<-array(rep(RateTreatFirst*treat_dtg[3],2)*rep(c(1,1),each=8),dim=c(2,4,2))*
      smooth_st(2005+t/12,2019,2019+treat_dtg["delay_first"],0,1) #delay
    #-------------------------------------------------------------------------------------------------------------
    #2.2.1 RateSwitchDTG[care,cd4,res,gender]
    #1) All people switch at a rate equal to RateTreatSecond
    #RateSwitchDTG<-array(rep(rep(treat_dtg[4]*RateTreatSecond,each=3)*rep(c(treat_dtg[5],treat_dtg[6],treat_dtg[7]),4),4),dim=c(3,4,2,2))*
    #smooth_st(2005+t/12,2019,2019+treat_dtg["delay_switch"],0,1) #delay
    RateSwitchDTG<-array(rep(rep(RateTreatSecond,each=3)*rep(c(treat_dtg[5],treat_dtg[6],treat_dtg[7]),4),4),dim=c(3,4,2,2))*
      smooth_st(2005+t/12,2019,2019+treat_dtg["delay_switch"],0,1) #delay
    #2) All people switch at a rate equal to 1/12 months^(-1)
    # RateSwitchDTG=array(rep(rep(treat_dtg[4]*rep(1/12,4),each=3)*rep(c(treat_dtg[5],treat_dtg[6],treat_dtg[7]),4),4),dim=c(3,4,2,2))*
    #    smooth_st(2005+t/12,2019,2019+treat_dtg["delay_switch"],0,1) #delay
    #3) Treated and suppressed people switch at a rate equal to 1year^-1, failed people at a rate still equal to 1) RateTreatSecond or 2) 1/12 months^(-1)
    RateSwitchDTG[1:2,1:4,1:2,1:2]<-rep(treat_dtg[4],32)*rep(rep(1/12,2)*c(treat_dtg[5],treat_dtg[6]),16)*
      smooth_st(2005+t/12,2019,2019+treat_dtg["delay_switch"],0,1) #delay
    #4) Failed people switch at a rate equal to 1year^-1 (if treat_dtg[3]==1 & treat_dtg[4]==1, DTG 1st-line and switch)
    #or at RateTreatSecond (if treat_dtg[3]==1 & treat_dtg[4]==0, DTG 1st-line), or 0 (if treat_dtg[3]==0 & treat_dtg[4]==0, no DTG)
    RateSwitchDTG[3,1:4,1:2,1:2]<-pmax(array(rep(1/(12),16)*treat_dtg[4]*treat_dtg[7],dim=c(4,2,2)),array(rep(RateTreatSecond,4),dim=c(4,2,2)))*
      treat_dtg[3]*smooth_st(2005+t/12,2019,2019+treat_dtg["delay_switch"],0,1)
    print(RateSwitchDTG)
    print(RateTreatSecond)
    #-------------------------------------------------------------------------------------------------------------
    #2.2.2 Modified treatment rates
    #RateTreatFirst [res,cd4,gender]
    RateTreatFirst_noDTG<-array(rep(RateTreatFirst,2),dim=c(2,4,2)) #Treatment rate to NNRTI for DTG-ineligible people
    RateTreatFirst<-apply(RateTreatFirst_noDTG-RateTreatFirstDTG,c(1,2,3),function(x) max(x,0))  #Treatment rate to NNRTI for DTG-eligible people
    
    # print("Treat NNRTI DTG inel --- [res,cd4,gender]")
    # print(RateTreatFirst_noDTG)
    # print("Treat NNRTI DTG elig --- [res,cd4,gender]")
    # print(RateTreatFirst)
    # print("Treat DTG --- [res,cd4,gender]")
    # print(RateTreatFirstDTG)
    #RateTreatSecond_noDTG[CD4,gender]
    #1) Treat rates to second for dtg-ineligible people stay the same
    RateTreatSecond_noDTG<-array(rep(RateTreatSecond,2),dim=c(4,2)) #Treatment rate to PI for DTG-ineligible people
    #2) Treat rates to second for dtg-ineligible people increase in 2019 (only do that when switching rate to dtg of failed people is 1/12, in order to have same switching rates for both dtg-ineligible and -eligible)
    # if(2005+t/12<2019){RateTreatSecond_noDTG=array(rep(RateTreatSecond,2),dim=c(4,2))
    # }else{RateTreatSecond_noDTG=array(rep(1/12,8),dim=c(4,2))}
    
    #RateTreatSecond[CD4,gender]
    RateTreatSecond<-apply(RateTreatSecond_noDTG-apply(RateSwitchDTG[3,1:4,1:2,1:2],c(1,3),mean),c(1,2),function(x) max(x,0)) #Treatment rate to PI for DTG-eligible people
    #we take the mean here as RateTreatSecond is the same for men and women. It has not any impact as long as RateSwitchDTG is also the same for men and women.
    #assume that second-line for DTG-eligible is DTG (not more PI) from 2019, even in the scenario of no switch to DTG
    ############
    #print(RateTreatSecond_noDTG)
    #print(RateTreatSecond)
    ############
    #3. Death
    p1<-min(1,q*mean(sapply(2005+seq(t-36,t,1)/12,function(x) RateTreatFirst_f(x)[4]/rate_treat)))
    p1=1-0.27*exp(-q*(mean(sapply(2005+seq(t-36,t,1)/12,function(x) RateTreatFirst_f(x)[4])-
                             sapply(2005+seq(-36,0,1)/12,function(x) RateTreatFirst_f(x)[4])))/rate_treat)

    mu[1:2,4,1]<-p1*40.9+(1-p1)*134.4
    mu[4,4,1]<-p1*8.3+(1-p1)*41.7
    mu[5,4,1]<-p1*11.8+(1-p1)*59.7
    mu[3,1:4,1]<-rate_ratio*mu[4,1:4,1]+(1-rate_ratio)*mu[5,1:4,1]
    mu[6:8,1:4,1]<-mu[3:5,1:4,1]
    mu[1:8,1:4,2]<-mu[1:8,1:4,1]
    mu<-mu*rate_death/1000
    
    #Back and forth mutation and impact on treatment outcome
    RateResistant<-1/rate_res
    RateSusceptible<-1/rate_susc
    #From T to S
    RateSuppFirstSusc<-RateSuppFirst
    RateSuppFirstResis<-1/alpha1*RateSuppFirst
    #From S to F
    RateFailFirstSusc<-RateFailFirst
    RateFailFirstResis<-alpha2*RateFailFirst
    #From T to F
    RateTreatToFailFirstSusc<-RateTreatToFailFirst
    RateTreatToFailFirstResis<-alpha1*RateTreatToFailFirst
    #From F to S
    RateFailToSuppTreatFirstSusc<-RateFailToSuppTreatFirst
    RateFailToSuppTreatFirstResis<-1/alpha2*RateFailToSuppTreatFirst
    
    #For second-line regimen, no effect of NNRTI resistance on treatment as it is PI
    ##################################################################################################################
    if(length(x)==250){
      I<-apply(array(x,dim=c(15,4,2,2)),4,sum)
    }else if(length(x)==298){
      I<-apply(array(x,dim=c(18,4,2,2)),4,sum)
    }else{
      stop(paste("length of x !=250 or 298 (my comment), length=", length(x)))
    }
    
    #First approach to model S: S=N-I
    N<-N_value_f(t+1)
    #N=c(1000000,1000000)
    
    # print("RateTreatFirst_noDTG")
    # print(RateTreatFirst_noDTG)
    # print("RateTreatFirst")
    # print(RateTreatFirst)
    # print("RateTreatFirstDTG")
    # print(RateTreatFirstDTG)
    # print("RateTreatSecond")
    # print(RateTreatSecond)
    # print("RateTreatSecond_noDTG")
    # print(RateTreatSecond_noDTG)
    # print("RateSwitchDTG")
    # print(RateSwitchDTG)
    
    S<-N-I
    #S=rep(0,0)
    #Percentage getting NNRTI-based 1st-line[gender]
    p_NNRTI<-sum(RateTreatFirst)/sum(RateTreatFirst+RateTreatFirstDTG)*c(1,treat_dtg[1])+c(0,1-treat_dtg[1])
  })
  return(e)
  #return(list(e[["RateSwitchDTG"]],e[["RateTreatSecond"]]))
}

#15 care stages
select_15=function(r1,r2,r3,r4){
  l1=15
  l2=4
  l3=2
  l4=2
  ar=array(1:(l1*l2*l3*l4),dim=c(l1,l2,l3,l4))
  return(as.vector(ar[r1,r2,r3,r4]))
}
xstart_f_dtg_2=function(par,treat_dtg){
  k1=par[9]
  k2=par[10]
  k3=par[11]
  
  start_value=array(0,dim=c(15,4,2,2))
  start_value[1,1:4,1,1:2]=c(undiag_start[1]*repart_f(k1),undiag_start[2]*repart_f(k1))
  start_value[2,1:4,1,1:2]=c(diag_start[1]*repart_f(k2),diag_start[2]*repart_f(k2))
  start_value[3,1:4,1,1:2]=c(treat_start[1]*repart_f(k3),treat_start[2]*repart_f(k3))
  
  start_value[4,1:4,1,1:2]=start_value[3,1:4,1,1:2]*0
  start_value[5,1:4,1,1:2]=start_value[3,1:4,1,1:2]*0
  start_value[3,1:4,1,1:2]=start_value[3,1:4,1,1:2]*1
  
  start_value[1:2,1:4,2,1:2]=array(rep(rep(c(0.01,0.01,0.01,0.01),each=2),2),dim=c(2,4,2))*start_value[1:2,1:4,1,1:2]
  start_value[1:2,1:4,1,1:2]=array(rep(rep(1-c(0.01,0.01,0.01,0.01),each=2),2),dim=c(2,4,2))*start_value[1:2,1:4,1,1:2]
  start_value[3,1:4,2,1:2]=0.01*start_value[3,1:4,1,1:2]
  start_value[3,1:4,1,1:2]=0.99*start_value[3,1:4,1,1:2]
  start_value[4,1:4,2,1:2]=0.01*start_value[4,1:4,1,1:2]
  start_value[4,1:4,1,1:2]=0.99*start_value[4,1:4,1,1:2]
  start_value[5,1:4,2,1:2]=0.7*start_value[5,1:4,1,1:2]
  start_value[5,1:4,1,1:2]=0.3*start_value[5,1:4,1,1:2]
  
  start_value[12:15,1:4,1:2,2]=(1-treat_dtg[1])*start_value[2:5,1:4,1:2,2]
  start_value[2:5,1:4,1:2,2]=treat_dtg[1]*start_value[2:5,1:4,1:2,2]
  
  xstart=c(start_value,rep(0,10))
  return(xstart)
}
mod_dtg=function(t,x,p1,treat_dtg,parms,p2){
  
  #parms2=mod_rate(t,x,p1,treat_dtg,parms,p2)
  parms2=mod_rate_2res(t,x,p1,treat_dtg,parms,p2) #2 rates to model the impact of res on treatment
  dx=array(0,dim=c(15,4,2,2))
  x=array(x,dim=c(15,4,2,2))
  
  #People reaching 15 not considered, because uncertainty about proportion of resistance,
  #because you have to take into account both mtct of res. and acquisition of resistance due to pmtct (was not needed to be considered in marisa model 2005-2016)
  
  # dx[1,1:4,1,1:2]=rep(c(predict(inf_m,1+t/12)$y/4,predict(inf_w,1+t/12)$y/4)/12,each=4)
  # dx[1,1:4,2,1:2]=0
  # dx[2,1:4,1,1:2]=rep(c(predict(diag_m,1+t/12)$y/4,predict(diag_w,1+t/12)$y/4)/12,each=4)
  # dx[2,1:4,2,1:2]=0
  # dx[3,1:4,1,1:2]=(1-0.9*0.2)*rep(c(predict(treat_m,1+t/12)$y/4,predict(treat_w,1+t/12)$y/4)/12,each=4)
  # dx[3,1:4,2,1:2]=0.9*0.2*rep(c(predict(treat_m,1+t/12)$y/4,predict(treat_w,1+t/12)$y/4)/12,each=4)
  
  with(as.list(parms2), {
    for(j in 1:2){
      #dx[x1,x2,x3,x4] x1:care stages x2:disease progression x3:resistance x4:gender
      #dx with only interaction with the neighbours in the compartmental model
      dx[1,1,1,j]=S[j]/N[j]*(sum(RateInf1[1:4,1:2,j]*apply(x[c(1),1:4,1,1:2],c(1,2),sum))+sum(RateInf2[1:4,1:2,j]*apply(x[c(2,3,5,6,8,9,11,12,13,15),1:4,1,1:2],c(2,3),sum)))-
        (RateStageInf[1]+RateDiag[1,j]+mu[1,1,1])*x[1,1,1,j]
      print("-----------")
      print(RateDiag[1,j]*x[1,1,1,j])
      dx[2,1,1,j]=dx[2,1,1,j]+RateDiag[1,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,1,1,j]-(RateTreatFirst[1,1,j]+RateStageDiag[1]+mu[2,1,1])*x[2,1,1,j]
      dx[3,1,1,j]=dx[3,1,1,j]+RateTreatFirst[1,1,j]*x[2,1,1,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,1])*x[3,1,1,j]+RateStageTreatFirstL[1]*x[3,2,1,j]
      dx[4,1,1,j]=RateSuppFirstSusc[1]*x[3,1,1,j]-(RateFailFirstSusc[1]+mu[4,1,1])*x[4,1,1,j]+RateStageSuppFirst[1]*x[4,2,1,j]
      dx[5,1,1,j]=RateFailFirstSusc[1]*x[4,1,1,j]-(RateTreatSecond[1,j]+RateStageFailFirst[1]+mu[5,1,1])*x[5,1,1,j]
      dx[6,1,1,j]=RateTreatSecond_noDTG[1,j]*x[11,1,1,j]+RateTreatSecond_noDTG[1,j]*x[15,1,1,j]+RateTreatSecond[1,j]*x[5,1,1,j]-(RateSuppSecond[1]+RateStageTreatSecondR[1]+mu[6,1,1])*x[6,1,1,j]+RateStageTreatSecondL[1]*x[6,2,1,j]
      dx[7,1,1,j]=RateSuppSecond[1]*x[6,1,1,j]-(RateFailSecond[1]+mu[7,1,1])*x[7,1,1,j]+RateStageSuppSecond[1]*x[7,2,1,j]
      dx[8,1,1,j]=RateFailSecond[1]*x[7,1,1,j]-(RateStageFailSecond[1]+mu[8,1,1])*x[8,1,1,j]
      
      dx[1,2,1,j]=-(RateStageInf[2]+RateDiag[2,j]+mu[1,2,1])*x[1,2,1,j]+RateStageInf[1]*x[1,1,1,j]
      dx[2,2,1,j]=dx[2,2,1,j]+RateDiag[2,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,2,1,j]-(RateTreatFirst[1,2,j]+RateStageDiag[2]+mu[2,2,1])*x[2,2,1,j]+RateStageDiag[1]*x[2,1,1,j]
      dx[3,2,1,j]=dx[3,2,1,j]+RateTreatFirst[1,2,j]*x[2,2,1,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,1])*x[3,2,1,j]+RateStageTreatFirstR[1]*x[3,1,1,j]+RateStageTreatFirstL[2]*x[3,3,1,j]
      dx[4,2,1,j]=RateSuppFirstSusc[2]*x[3,2,1,j]-(RateFailFirstSusc[2]+RateStageSuppFirst[1]+mu[4,2,1])*x[4,2,1,j]+RateStageSuppFirst[2]*x[4,3,1,j]
      dx[5,2,1,j]=RateFailFirstSusc[2]*x[4,2,1,j]-(RateTreatSecond[2,j]+RateStageFailFirst[2]+mu[5,2,1])*x[5,2,1,j]+RateStageFailFirst[1]*x[5,1,1,j]
      dx[6,2,1,j]=RateTreatSecond_noDTG[2,j]*x[11,2,1,j]+RateTreatSecond_noDTG[2,j]*x[15,2,1,j]+RateTreatSecond[2,j]*x[5,2,1,j]-(RateSuppSecond[2]+RateStageTreatSecondR[2]+RateStageTreatSecondL[1]+mu[6,2,1])*x[6,2,1,j]+RateStageTreatSecondR[1]*x[6,1,1,j]+RateStageTreatSecondL[2]*x[6,3,1,j]
      dx[7,2,1,j]=RateSuppSecond[2]*x[6,2,1,j]-(RateFailSecond[2]+RateStageSuppSecond[1]+mu[7,2,1])*x[7,2,1,j]+RateStageSuppSecond[2]*x[7,3,1,j]
      dx[8,2,1,j]=RateFailSecond[2]*x[7,2,1,j]-(RateStageFailSecond[2]+mu[8,2,1])*x[8,2,1,j]+RateStageFailSecond[1]*x[8,1,1,j]
      
      dx[1,3,1,j]=-(RateStageInf[3]+RateDiag[3,j]+mu[1,3,1])*x[1,3,1,j]+RateStageInf[2]*x[1,2,1,j]
      dx[2,3,1,j]=dx[2,3,1,j]+RateDiag[3,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,3,1,j]-(RateTreatFirst[1,3,j]+RateStageDiag[3]+mu[2,3,1])*x[2,3,1,j]+RateStageDiag[2]*x[2,2,1,j]
      dx[3,3,1,j]=dx[3,3,1,j]+RateTreatFirst[1,3,j]*x[2,3,1,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,1])*x[3,3,1,j]+RateStageTreatFirstR[2]*x[3,2,1,j]+RateStageTreatFirstL[3]*x[3,4,1,j]
      dx[4,3,1,j]=RateSuppFirstSusc[3]*x[3,3,1,j]-(RateFailFirstSusc[3]+RateStageSuppFirst[2]+mu[4,3,1])*x[4,3,1,j]+RateStageSuppFirst[3]*x[4,4,1,j]
      dx[5,3,1,j]=RateFailFirstSusc[3]*x[4,3,1,j]-(RateTreatSecond[3,j]+RateStageFailFirst[3]+mu[5,3,1])*x[5,3,1,j]+RateStageFailFirst[2]*x[5,2,1,j]
      dx[6,3,1,j]=RateTreatSecond_noDTG[3,j]*x[11,3,1,j]+RateTreatSecond_noDTG[3,j]*x[15,3,1,j]+RateTreatSecond[3,j]*x[5,3,1,j]-(RateSuppSecond[3]+RateStageTreatSecondR[3]+RateStageTreatSecondL[2]+mu[6,3,1])*x[6,3,1,j]+RateStageTreatSecondR[2]*x[6,2,1,j]+RateStageTreatSecondL[3]*x[6,4,1,j]
      dx[7,3,1,j]=RateSuppSecond[3]*x[6,3,1,j]-(RateFailSecond[3]+RateStageSuppSecond[2]+mu[7,3,1])*x[7,3,1,j]+RateStageSuppSecond[3]*x[7,4,1,j]
      dx[8,3,1,j]=RateFailSecond[3]*x[7,3,1,j]-(RateStageFailSecond[3]+mu[8,3,1])*x[8,3,1,j]+RateStageFailSecond[2]*x[8,2,1,j]
      
      dx[1,4,1,j]=-(RateDiag[4,j]+mu[1,4,1])*x[1,4,1,j]+RateStageInf[3]*x[1,3,1,j]
      dx[2,4,1,j]=dx[2,4,1,j]+RateDiag[4,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,4,1,j]-(RateTreatFirst[1,4,j]+mu[2,4,1])*x[2,4,1,j]+RateStageDiag[3]*x[2,3,1,j]
      dx[3,4,1,j]=dx[3,4,1,j]+RateTreatFirst[1,4,j]*x[2,4,1,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,1])*x[3,4,1,j]+RateStageTreatFirstR[3]*x[3,3,1,j]
      dx[4,4,1,j]=RateSuppFirstSusc[4]*x[3,4,1,j]-(RateFailFirstSusc[4]+RateStageSuppFirst[3]+mu[4,4,1])*x[4,4,1,j]
      dx[5,4,1,j]=RateFailFirstSusc[4]*x[4,4,1,j]-(RateTreatSecond[4,j]+mu[5,4,1])*x[5,4,1,j]+RateStageFailFirst[3]*x[5,3,1,j]
      dx[6,4,1,j]=RateTreatSecond_noDTG[4,j]*x[11,4,1,j]+RateTreatSecond_noDTG[4,j]*x[15,4,1,j]+RateTreatSecond[4,j]*x[5,4,1,j]-(RateSuppSecond[4]+RateStageTreatSecondL[3]+mu[6,4,1])*x[6,4,1,j]+RateStageTreatSecondR[3]*x[6,3,1,j]
      dx[7,4,1,j]=RateSuppSecond[4]*x[6,4,1,j]-(RateFailSecond[4]+RateStageSuppSecond[3]+mu[7,4,1])*x[7,4,1,j]
      dx[8,4,1,j]=RateFailSecond[4]*x[7,4,1,j]-mu[8,4,1 ]*x[8,4,1,j]+RateStageFailSecond[3]*x[8,3,1,j]
      
      ##############################################################################################################################################################
      
      dx[1,1,1,j]=dx[1,1,1,j]
      dx[2,1,1,j]=dx[2,1,1,j]-RateDirectTreatSecond[1,1]*x[2,1,1,j]+RateStopTreatFirst[1]*x[3,1,1,j]+RateStopSuppFirst[1]*x[4,1,1,j]+RateStopFailFirst[1]*x[5,1,1,j]+
        (RateStopTreatSecond[1]*x[6,1,1,j]+RateStopSuppSecond[1]*x[7,1,1,j]+RateStopFailSecond[1]*x[8,1,1,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,1,1,j]=dx[3,1,1,j]-RateStopTreatFirst[1]*x[3,1,1,j]-RateTreatToFailFirstSusc[1]*x[3,1,1,j]
      dx[4,1,1,j]=dx[4,1,1,j]-RateStopSuppFirst[1]*x[4,1,1,j]+RateFailToSuppTreatFirstSusc[1]*x[5,1,1,j]-RateSuppFirstToSecond[1]*x[4,1,1,j]
      dx[5,1,1,j]=dx[5,1,1,j]-RateStopFailFirst[1]*x[5,1,1,j]-RateFailToSuppTreatFirstSusc[1]*x[5,1,1,j]+RateTreatToFailFirstSusc[1]*x[3,1,1,j]
      dx[6,1,1,j]=dx[6,1,1,j]-RateStopTreatSecond[1]*x[6,1,1,j]+RateDirectTreatSecond[1,1]*x[2,1,1,j]+RateDirectTreatSecond[1,1]*x[12,1,1,j]-RateTreatToFailSecond[1]*x[6,1,1,j]
      dx[7,1,1,j]=dx[7,1,1,j]-RateStopSuppSecond[1]*x[7,1,1,j]+RateFailToSuppTreatSecond[1]*x[8,1,1,j]+RateSuppFirstToSecond[1]*x[4,1,1,j]
      dx[8,1,1,j]=dx[8,1,1,j]-RateStopFailSecond[1]*x[8,1,1,j]-RateFailToSuppTreatSecond[1]*x[8,1,1,j]+RateTreatToFailSecond[1]*x[6,1,1,j]
      
      dx[1,2,1,j]=dx[1,2,1,j]
      dx[2,2,1,j]=dx[2,2,1,j]-RateDirectTreatSecond[1,2]*x[2,2,1,j]+RateStopTreatFirst[2]*x[3,2,1,j]+RateStopSuppFirst[2]*x[4,2,1,j]+RateStopFailFirst[2]*x[5,2,1,j]+
        (RateStopTreatSecond[2]*x[6,2,1,j]+RateStopSuppSecond[2]*x[7,2,1,j]+RateStopFailSecond[2]*x[8,2,1,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,2,1,j]=dx[3,2,1,j]-RateStopTreatFirst[2]*x[3,2,1,j]-RateTreatToFailFirstSusc[2]*x[3,2,1,j]
      dx[4,2,1,j]=dx[4,2,1,j]-RateStopSuppFirst[2]*x[4,2,1,j]+RateFailToSuppTreatFirstSusc[2]*x[5,2,1,j]-RateSuppFirstToSecond[2]*x[4,2,1,j]
      dx[5,2,1,j]=dx[5,2,1,j]-RateStopFailFirst[2]*x[5,2,1,j]-RateFailToSuppTreatFirstSusc[2]*x[5,2,1,j]+RateTreatToFailFirstSusc[2]*x[3,2,1,j]
      dx[6,2,1,j]=dx[6,2,1,j]-RateStopTreatSecond[2]*x[6,2,1,j]+RateDirectTreatSecond[1,2]*x[2,2,1,j]+RateDirectTreatSecond[1,2]*x[12,2,1,j]-RateTreatToFailSecond[2]*x[6,2,1,j]
      dx[7,2,1,j]=dx[7,2,1,j]-RateStopSuppSecond[2]*x[7,2,1,j]+RateFailToSuppTreatSecond[2]*x[8,2,1,j]+RateSuppFirstToSecond[2]*x[4,2,1,j]
      dx[8,2,1,j]=dx[8,2,1,j]-RateStopFailSecond[2]*x[8,2,1,j]-RateFailToSuppTreatSecond[2]*x[8,2,1,j]+RateTreatToFailSecond[2]*x[6,2,1,j]
      
      dx[1,3,1,j]=dx[1,3,1,j]
      dx[2,3,1,j]=dx[2,3,1,j]-RateDirectTreatSecond[1,3]*x[2,3,1,j]+RateStopTreatFirst[3]*x[3,3,1,j]+RateStopSuppFirst[3]*x[4,3,1,j]+RateStopFailFirst[3]*x[5,3,1,j]+
        (RateStopTreatSecond[3]*x[6,3,1,j]+RateStopSuppSecond[3]*x[7,3,1,j]+RateStopFailSecond[3]*x[8,3,1,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,3,1,j]=dx[3,3,1,j]-RateStopTreatFirst[3]*x[3,3,1,j]-RateTreatToFailFirstSusc[3]*x[3,3,1,j]
      dx[4,3,1,j]=dx[4,3,1,j]-RateStopSuppFirst[3]*x[4,3,1,j]+RateFailToSuppTreatFirstSusc[3]*x[5,3,1,j]-RateSuppFirstToSecond[3]*x[4,3,1,j]
      dx[5,3,1,j]=dx[5,3,1,j]-RateStopFailFirst[3]*x[5,3,1,j]-RateFailToSuppTreatFirstSusc[3]*x[5,3,1,j]+RateTreatToFailFirstSusc[3]*x[3,3,1,j]
      dx[6,3,1,j]=dx[6,3,1,j]-RateStopTreatSecond[3]*x[6,3,1,j]+RateDirectTreatSecond[1,3]*x[2,3,1,j]+RateDirectTreatSecond[1,3]*x[12,3,1,j]-RateTreatToFailSecond[3]*x[6,3,1,j]
      dx[7,3,1,j]=dx[7,3,1,j]-RateStopSuppSecond[3]*x[7,3,1,j]+RateFailToSuppTreatSecond[3]*x[8,3,1,j]+RateSuppFirstToSecond[3]*x[4,3,1,j]
      dx[8,3,1,j]=dx[8,3,1,j]-RateStopFailSecond[3]*x[8,3,1,j]-RateFailToSuppTreatSecond[3]*x[8,3,1,j]+RateTreatToFailSecond[3]*x[6,3,1,j]
      
      dx[1,4,1,j]=dx[1,4,1,j]
      dx[2,4,1,j]=dx[2,4,1,j]-RateDirectTreatSecond[1,4]*x[2,4,1,j]+RateStopTreatFirst[4]*x[3,4,1,j]+RateStopSuppFirst[4]*x[4,4,1,j]+RateStopFailFirst[4]*x[5,4,1,j]+
        (RateStopTreatSecond[4]*x[6,4,1,j]+RateStopSuppSecond[4]*x[7,4,1,j]+RateStopFailSecond[4]*x[8,4,1,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,4,1,j]=dx[3,4,1,j]-RateStopTreatFirst[4]*x[3,4,1,j]-RateTreatToFailFirstSusc[4]*x[3,4,1,j]
      dx[4,4,1,j]=dx[4,4,1,j]-RateStopSuppFirst[4]*x[4,4,1,j]+RateFailToSuppTreatFirstSusc[4]*x[5,4,1,j]-RateSuppFirstToSecond[4]*x[4,4,1,j]
      dx[5,4,1,j]=dx[5,4,1,j]-RateStopFailFirst[4]*x[5,4,1,j]-RateFailToSuppTreatFirstSusc[4]*x[5,4,1,j]+RateTreatToFailFirstSusc[4]*x[3,4,1,j]
      dx[6,4,1,j]=dx[6,4,1,j]-RateStopTreatSecond[4]*x[6,4,1,j]+RateDirectTreatSecond[1,4]*x[2,4,1,j]+RateDirectTreatSecond[1,4]*x[12,4,1,j]-RateTreatToFailSecond[4]*x[6,4,1,j]
      dx[7,4,1,j]=dx[7,4,1,j]-RateStopSuppSecond[4]*x[7,4,1,j]+RateFailToSuppTreatSecond[4]*x[8,4,1,j]+RateSuppFirstToSecond[4]*x[4,4,1,j]
      dx[8,4,1,j]=dx[8,4,1,j]-RateStopFailSecond[4]*x[8,4,1,j]-RateFailToSuppTreatSecond[4]*x[8,4,1,j]+RateTreatToFailSecond[4]*x[6,4,1,j]
      #############################################################################################################################################################################
      #############################################################################################################################################################################
      dx[1,1,2,j]=S[j]/N[j]*(sum(RateInf1[1:4,1:2,j]*apply(x[c(1),1:4,2,1:2],c(1,2),sum))+sum(RateInf2[1:4,1:2,j]*apply(x[c(2,3,5,6,8,9,11,12,13,15),1:4,2,1:2],c(2,3),sum)))-
        (RateStageInf[1]+RateDiag[1,j]+mu[1,1,2 ])*x[1,1,2,j]
      dx[2,1,2,j]=dx[2,1,2,j]+RateDiag[1,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,1,2,j]-(RateTreatFirst[2,1,j]+RateStageDiag[1]+mu[2,1,2])*x[2,1,2,j]
      dx[3,1,2,j]=dx[3,1,2,j]+RateTreatFirst[2,1,j]*x[2,1,2,j]-(RateSuppFirstResis[1]+RateStageTreatFirstR[1]+mu[3,1,2])*x[3,1,2,j]+RateStageTreatFirstL[1]*x[3,2,2,j]
      dx[4,1,2,j]=RateSuppFirstResis[1]*x[3,1,2,j]-(RateFailFirstResis[1]+mu[4,1,2])*x[4,1,2,j]+RateStageSuppFirst[1]*x[4,2,2,j]
      dx[5,1,2,j]=RateFailFirstResis[1]*x[4,1,2,j]-(RateTreatSecond[1,j]+RateStageFailFirst[1]+mu[5,1,2 ])*x[5,1,2,j]
      dx[6,1,2,j]=RateTreatSecond_noDTG[1,j]*x[11,1,2,j]+RateTreatSecond_noDTG[1,j]*x[15,1,2,j]+RateTreatSecond[1,j]*x[5,1,2,j]-(RateSuppSecond[1]+RateStageTreatSecondR[1]+mu[6,1,2])*x[6,1,2,j]+RateStageTreatSecondL[1]*x[6,2,2,j]
      dx[7,1,2,j]=RateSuppSecond[1]*x[6,1,2,j]-(RateFailSecond[1]+mu[7,1,2])*x[7,1,2,j]+RateStageSuppSecond[1]*x[7,2,2,j]
      dx[8,1,2,j]=RateFailSecond[1]*x[7,1,2,j]-(RateStageFailSecond[1]+mu[8,1,2])*x[8,1,2,j]
      
      dx[1,2,2,j]=-(RateStageInf[2]+RateDiag[2,j]+mu[1,2,2])*x[1,2,2,j]+RateStageInf[1]*x[1,1,2,j]
      dx[2,2,2,j]=dx[2,2,2,j]+RateDiag[2,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,2,2,j]-(RateTreatFirst[2,2,j]+RateStageDiag[2]+mu[2,2,2])*x[2,2,2,j]+RateStageDiag[1]*x[2,1,2,j]
      dx[3,2,2,j]=dx[3,2,2,j]+RateTreatFirst[2,2,j]*x[2,2,2,j]-(RateSuppFirstResis[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,2])*x[3,2,2,j]+RateStageTreatFirstR[1]*x[3,1,2,j]+RateStageTreatFirstL[2]*x[3,3,2,j]
      dx[4,2,2,j]=RateSuppFirstResis[2]*x[3,2,2,j]-(RateFailFirstResis[2]+RateStageSuppFirst[1]+mu[4,2,2])*x[4,2,2,j]+RateStageSuppFirst[2]*x[4,3,2,j]
      dx[5,2,2,j]=RateFailFirstResis[2]*x[4,2,2,j]-(RateTreatSecond[2,j]+RateStageFailFirst[2]+mu[5,2,2])*x[5,2,2,j]+RateStageFailFirst[1]*x[5,1,2,j]
      dx[6,2,2,j]=RateTreatSecond_noDTG[2,j]*x[11,2,2,j]+RateTreatSecond_noDTG[2,j]*x[15,2,2,j]+RateTreatSecond[2,j]*x[5,2,2,j]-(RateSuppSecond[2]+RateStageTreatSecondR[2]+RateStageTreatSecondL[1]+mu[6,2,2])*x[6,2,2,j]+RateStageTreatSecondR[1]*x[6,1,2,j]+RateStageTreatSecondL[2]*x[6,3,2,j]
      dx[7,2,2,j]=RateSuppSecond[2]*x[6,2,2,j]-(RateFailSecond[2]+RateStageSuppSecond[1]+mu[7,2,2])*x[7,2,2,j]+RateStageSuppSecond[2]*x[7,3,2,j]
      dx[8,2,2,j]=RateFailSecond[2]*x[7,2,2,j]-(RateStageFailSecond[2]+mu[8,2,2])*x[8,2,2,j]+RateStageFailSecond[1]*x[8,1,2,j]
      
      dx[1,3,2,j]=-(RateStageInf[3]+RateDiag[3,j]+mu[1,3,2])*x[1,3,2,j]+RateStageInf[2]*x[1,2,2,j]
      dx[2,3,2,j]=dx[2,3,2,j]+RateDiag[3,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,3,2,j]-(RateTreatFirst[2,3,j]+RateStageDiag[3]+mu[2,3,2])*x[2,3,2,j]+RateStageDiag[2]*x[2,2,2,j]
      dx[3,3,2,j]=dx[3,3,2,j]+RateTreatFirst[2,3,j]*x[2,3,2,j]-(RateSuppFirstResis[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,2])*x[3,3,2,j]+RateStageTreatFirstR[2]*x[3,2,2,j]+RateStageTreatFirstL[3]*x[3,4,2,j]
      dx[4,3,2,j]=RateSuppFirstResis[3]*x[3,3,2,j]-(RateFailFirstResis[3]+RateStageSuppFirst[2]+mu[4,3,2])*x[4,3,2,j]+RateStageSuppFirst[3]*x[4,4,2,j]
      dx[5,3,2,j]=RateFailFirstResis[3]*x[4,3,2,j]-(RateTreatSecond[3,j]+RateStageFailFirst[3]+mu[5,3,2])*x[5,3,2,j]+RateStageFailFirst[2]*x[5,2,2,j]
      dx[6,3,2,j]=RateTreatSecond_noDTG[3,j]*x[11,3,2,j]+RateTreatSecond_noDTG[3,j]*x[15,3,2,j]+RateTreatSecond[3,j]*x[5,3,2,j]-(RateSuppSecond[3]+RateStageTreatSecondR[3]+RateStageTreatSecondL[2]+mu[6,3,2])*x[6,3,2,j]+RateStageTreatSecondR[2]*x[6,2,2,j]+RateStageTreatSecondL[3]*x[6,4,2,j]
      dx[7,3,2,j]=RateSuppSecond[3]*x[6,3,2,j]-(RateFailSecond[3]+RateStageSuppSecond[2]+mu[7,3,2])*x[7,3,2,j]+RateStageSuppSecond[3]*x[7,4,2,j]
      dx[8,3,2,j]=RateFailSecond[3]*x[7,3,2,j]-(RateStageFailSecond[3]+mu[8,3,2])*x[8,3,2,j]+RateStageFailSecond[2]*x[8,2,2,j]
      
      dx[1,4,2,j]=-(RateDiag[4,j]+mu[1,4,2])*x[1,4,2,j]+RateStageInf[3]*x[1,3,2,j]
      dx[2,4,2,j]=dx[2,4,2,j]+RateDiag[4,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,4,2,j]-(RateTreatFirst[2,4,j]+mu[2,4,2])*x[2,4,2,j]+RateStageDiag[3]*x[2,3,2,j]
      dx[3,4,2,j]=dx[3,4,2,j]+RateTreatFirst[2,4,j]*x[2,4,2,j]-(RateSuppFirstResis[4]+RateStageTreatFirstL[3]+mu[3,4,2])*x[3,4,2,j]+RateStageTreatFirstR[3]*x[3,3,2,j]
      dx[4,4,2,j]=RateSuppFirstResis[4]*x[3,4,2,j]-(RateFailFirstResis[4]+RateStageSuppFirst[3]+mu[4,4,2])*x[4,4,2,j]
      dx[5,4,2,j]=RateFailFirstResis[4]*x[4,4,2,j]-(RateTreatSecond[4,j]+mu[5,4,2])*x[5,4,2,j]+RateStageFailFirst[3]*x[5,3,2,j]
      dx[6,4,2,j]=RateTreatSecond_noDTG[4,j]*x[11,4,2,j]+RateTreatSecond_noDTG[4,j]*x[15,4,2,j]+RateTreatSecond[4,j]*x[5,4,2,j]-(RateSuppSecond[4]+RateStageTreatSecondL[3]+mu[6,4,2])*x[6,4,2,j]+RateStageTreatSecondR[3]*x[6,3,2,j]
      dx[7,4,2,j]=RateSuppSecond[4]*x[6,4,2,j]-(RateFailSecond[4]+RateStageSuppSecond[3]+mu[7,4,2])*x[7,4,2,j]
      dx[8,4,2,j]=RateFailSecond[4]*x[7,4,2,j]-mu[8,4,2]*x[8,4,2,j]+RateStageFailSecond[3]*x[8,3,2,j]
      
      ##############################################################################################################################################################
      
      dx[1,1,2,j]=dx[1,1,2,j]
      dx[2,1,2,j]=dx[2,1,2,j]-RateDirectTreatSecond[2,1]*x[2,1,2,j]+RateStopTreatFirst[1]*x[3,1,2,j]+RateStopSuppFirst[1]*x[4,1,2,j]+RateStopFailFirst[1]*x[5,1,2,j]+
        (RateStopTreatSecond[1]*x[6,1,2,j]+RateStopSuppSecond[1]*x[7,1,2,j]+RateStopFailSecond[1]*x[8,1,2,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,1,2,j]=dx[3,1,2,j]-RateStopTreatFirst[1]*x[3,1,2,j]-RateTreatToFailFirstResis[1]*x[3,1,2,j]
      dx[4,1,2,j]=dx[4,1,2,j]-RateStopSuppFirst[1]*x[4,1,2,j]+RateFailToSuppTreatFirstResis[1]*x[5,1,2,j]-RateSuppFirstToSecond[1]*x[4,1,2,j]
      dx[5,1,2,j]=dx[5,1,2,j]-RateStopFailFirst[1]*x[5,1,2,j]-RateFailToSuppTreatFirstResis[1]*x[5,1,2,j]+RateTreatToFailFirstResis[1]*x[3,1,2,j]
      dx[6,1,2,j]=dx[6,1,2,j]-RateStopTreatSecond[1]*x[6,1,2,j]+RateDirectTreatSecond[2,1]*x[2,1,2,j]+RateDirectTreatSecond[2,1]*x[12,1,2,j]-RateTreatToFailSecond[1]*x[6,1,2,j]
      dx[7,1,2,j]=dx[7,1,2,j]-RateStopSuppSecond[1]*x[7,1,2,j]+RateFailToSuppTreatSecond[1]*x[8,1,2,j]+RateSuppFirstToSecond[1]*x[4,1,2,j]
      dx[8,1,2,j]=dx[8,1,2,j]-RateStopFailSecond[1]*x[8,1,2,j]-RateFailToSuppTreatSecond[1]*x[8,1,2,j]+RateTreatToFailSecond[1]*x[6,1,2,j]
      
      dx[1,2,2,j]=dx[1,2,2,j]
      dx[2,2,2,j]=dx[2,2,2,j]-RateDirectTreatSecond[2,2]*x[2,2,2,j]+RateStopTreatFirst[2]*x[3,2,2,j]+RateStopSuppFirst[2]*x[4,2,2,j]+RateStopFailFirst[2]*x[5,2,2,j]+
        (RateStopTreatSecond[2]*x[6,2,2,j]+RateStopSuppSecond[2]*x[7,2,2,j]+RateStopFailSecond[2]*x[8,2,2,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,2,2,j]=dx[3,2,2,j]-RateStopTreatFirst[2]*x[3,2,2,j]-RateTreatToFailFirstResis[2]*x[3,2,2,j]
      dx[4,2,2,j]=dx[4,2,2,j]-RateStopSuppFirst[2]*x[4,2,2,j]+RateFailToSuppTreatFirstResis[2]*x[5,2,2,j]-RateSuppFirstToSecond[2]*x[4,2,2,j]
      dx[5,2,2,j]=dx[5,2,2,j]-RateStopFailFirst[2]*x[5,2,2,j]-RateFailToSuppTreatFirstResis[2]*x[5,2,2,j]+RateTreatToFailFirstResis[2]*x[3,2,2,j]
      dx[6,2,2,j]=dx[6,2,2,j]-RateStopTreatSecond[2]*x[6,2,2,j]+RateDirectTreatSecond[2,2]*x[2,2,2,j]+RateDirectTreatSecond[2,2]*x[12,2,2,j]-RateTreatToFailSecond[2]*x[6,2,2,j]
      dx[7,2,2,j]=dx[7,2,2,j]-RateStopSuppSecond[2]*x[7,2,2,j]+RateFailToSuppTreatSecond[2]*x[8,2,2,j]+RateSuppFirstToSecond[2]*x[4,2,2,j]
      dx[8,2,2,j]=dx[8,2,2,j]-RateStopFailSecond[2]*x[8,2,2,j]-RateFailToSuppTreatSecond[2]*x[8,2,2,j]+RateTreatToFailSecond[2]*x[6,2,2,j]
      
      dx[1,3,2,j]=dx[1,3,2,j]
      dx[2,3,2,j]=dx[2,3,2,j]-RateDirectTreatSecond[2,3]*x[2,3,2,j]+RateStopTreatFirst[3]*x[3,3,2,j]+RateStopSuppFirst[3]*x[4,3,2,j]+RateStopFailFirst[3]*x[5,3,2,j]+
        (RateStopTreatSecond[3]*x[6,3,2,j]+RateStopSuppSecond[3]*x[7,3,2,j]+RateStopFailSecond[3]*x[8,3,2,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,3,2,j]=dx[3,3,2,j]-RateStopTreatFirst[3]*x[3,3,2,j]-RateTreatToFailFirstResis[3]*x[3,3,2,j]
      dx[4,3,2,j]=dx[4,3,2,j]-RateStopSuppFirst[3]*x[4,3,2,j]+RateFailToSuppTreatFirstResis[3]*x[5,3,2,j]-RateSuppFirstToSecond[3]*x[4,3,2,j]
      dx[5,3,2,j]=dx[5,3,2,j]-RateStopFailFirst[3]*x[5,3,2,j]-RateFailToSuppTreatFirstResis[3]*x[5,3,2,j]+RateTreatToFailFirstResis[3]*x[3,3,2,j]
      dx[6,3,2,j]=dx[6,3,2,j]-RateStopTreatSecond[3]*x[6,3,2,j]+RateDirectTreatSecond[2,3]*x[2,3,2,j]+RateDirectTreatSecond[2,3]*x[12,3,2,j]-RateTreatToFailSecond[3]*x[6,3,2,j]
      dx[7,3,2,j]=dx[7,3,2,j]-RateStopSuppSecond[3]*x[7,3,2,j]+RateFailToSuppTreatSecond[3]*x[8,3,2,j]+RateSuppFirstToSecond[3]*x[4,3,2,j]
      dx[8,3,2,j]=dx[8,3,2,j]-RateStopFailSecond[3]*x[8,3,2,j]-RateFailToSuppTreatSecond[3]*x[8,3,2,j]+RateTreatToFailSecond[3]*x[6,3,2,j]
      
      dx[1,4,2,j]=dx[1,4,2,j]
      dx[2,4,2,j]=dx[2,4,2,j]-RateDirectTreatSecond[2,4]*x[2,4,2,j]+RateStopTreatFirst[4]*x[3,4,2,j]+RateStopSuppFirst[4]*x[4,4,2,j]+RateStopFailFirst[4]*x[5,4,2,j]+
        (RateStopTreatSecond[4]*x[6,4,2,j]+RateStopSuppSecond[4]*x[7,4,2,j]+RateStopFailSecond[4]*x[8,4,2,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,4,2,j]=dx[3,4,2,j]-RateStopTreatFirst[4]*x[3,4,2,j]-RateTreatToFailFirstResis[4]*x[3,4,2,j]
      dx[4,4,2,j]=dx[4,4,2,j]-RateStopSuppFirst[4]*x[4,4,2,j]+RateFailToSuppTreatFirstResis[4]*x[5,4,2,j]-RateSuppFirstToSecond[4]*x[4,4,2,j]
      dx[5,4,2,j]=dx[5,4,2,j]-RateStopFailFirst[4]*x[5,4,2,j]-RateFailToSuppTreatFirstResis[4]*x[5,4,2,j]+RateTreatToFailFirstResis[4]*x[3,4,2,j]
      dx[6,4,2,j]=dx[6,4,2,j]-RateStopTreatSecond[4]*x[6,4,2,j]+RateDirectTreatSecond[2,4]*x[2,4,2,j]+RateDirectTreatSecond[2,4]*x[12,4,2,j]-RateTreatToFailSecond[4]*x[6,4,2,j]
      dx[7,4,2,j]=dx[7,4,2,j]-RateStopSuppSecond[4]*x[7,4,2,j]+RateFailToSuppTreatSecond[4]*x[8,4,2,j]+RateSuppFirstToSecond[4]*x[4,4,2,j]
      dx[8,4,2,j]=dx[8,4,2,j]-RateStopFailSecond[4]*x[8,4,2,j]-RateFailToSuppTreatSecond[4]*x[8,4,2,j]+RateTreatToFailSecond[4]*x[6,4,2,j]
      
      #############################################################################################################################################################################
      #############################################################################################################################################################################
      #dx with interaction between resistant and susceptible compartments
      dx[1,1,1,j]=dx[1,1,1,j]+RateSusceptible*x[1,1,2,j]
      dx[2,1,1,j]=dx[2,1,1,j]+RateSusceptible*x[2,1,2,j]
      dx[5,1,1,j]=dx[5,1,1,j]-RateResistant*x[5,1,1,j]
      
      dx[1,2,1,j]=dx[1,2,1,j]+RateSusceptible*x[1,2,2,j]
      dx[2,2,1,j]=dx[2,2,1,j]+RateSusceptible*x[2,2,2,j]
      dx[5,2,1,j]=dx[5,2,1,j]-RateResistant*x[5,2,1,j]
      
      dx[1,3,1,j]=dx[1,3,1,j]+RateSusceptible*x[1,3,2,j]
      dx[2,3,1,j]=dx[2,3,1,j]+RateSusceptible*x[2,3,2,j]
      dx[5,3,1,j]=dx[5,3,1,j]-RateResistant*x[5,3,1,j]
      
      dx[1,4,1,j]=dx[1,4,1,j]+RateSusceptible*x[1,4,2,j]
      dx[2,4,1,j]=dx[2,4,1,j]+RateSusceptible*x[2,4,2,j]
      dx[5,4,1,j]=dx[5,4,1,j]-RateResistant*x[5,4,1,j]
      
      #############################################################################################################################################################################
      dx[1,1,2,j]=dx[1,1,2,j]-RateSusceptible*x[1,1,2,j]
      dx[2,1,2,j]=dx[2,1,2,j]-RateSusceptible*x[2,1,2,j]
      dx[5,1,2,j]=dx[5,1,2,j]+RateResistant*x[5,1,1,j]
      
      dx[1,2,2,j]=dx[1,2,2,j]-RateSusceptible*x[1,2,2,j]
      dx[2,2,2,j]=dx[2,2,2,j]-RateSusceptible*x[2,2,2,j]
      dx[5,2,2,j]=dx[5,2,2,j]+RateResistant*x[5,2,1,j]
      
      dx[1,3,2,j]=dx[1,3,2,j]-RateSusceptible*x[1,3,2,j]
      dx[2,3,2,j]=dx[2,3,2,j]-RateSusceptible*x[2,3,2,j]
      dx[5,3,2,j]=dx[5,3,2,j]+RateResistant*x[5,3,1,j]
      
      dx[1,4,2,j]=dx[1,4,2,j]-RateSusceptible*x[1,4,2,j]
      dx[2,4,2,j]=dx[2,4,2,j]-RateSusceptible*x[2,4,2,j]
      dx[5,4,2,j]=dx[5,4,2,j]+RateResistant*x[5,4,1,j]
      
      
      #RateTreatFirstDTG[res,cd4,gender]------------------------------
      #NNRTI susceptible
      dx[2,1,1,j]=dx[2,1,1,j]-RateTreatFirstDTG[1,1,j]*x[2,1,1,j]
      dx[2,2,1,j]=dx[2,2,1,j]-RateTreatFirstDTG[1,2,j]*x[2,2,1,j]
      dx[2,3,1,j]=dx[2,3,1,j]-RateTreatFirstDTG[1,3,j]*x[2,3,1,j]
      dx[2,4,1,j]=dx[2,4,1,j]-RateTreatFirstDTG[1,4,j]*x[2,4,1,j]
      #NNRTI resistant
      dx[2,1,2,j]=dx[2,1,2,j]-RateTreatFirstDTG[2,1,j]*x[2,1,2,j]
      dx[2,2,2,j]=dx[2,2,2,j]-RateTreatFirstDTG[2,2,j]*x[2,2,2,j]
      dx[2,3,2,j]=dx[2,3,2,j]-RateTreatFirstDTG[2,3,j]*x[2,3,2,j]
      dx[2,4,2,j]=dx[2,4,2,j]-RateTreatFirstDTG[2,4,j]*x[2,4,2,j]
      
      #RateSwitchDTG[care,cd4,res,gender]----------------------------
      #NNRTI susceptible
      dx[3,1,1,j]=dx[3,1,1,j]-RateSwitchDTG[1,1,1,j]*x[3,1,1,j]
      dx[3,2,1,j]=dx[3,2,1,j]-RateSwitchDTG[1,2,1,j]*x[3,2,1,j]
      dx[3,3,1,j]=dx[3,3,1,j]-RateSwitchDTG[1,3,1,j]*x[3,3,1,j]
      dx[3,4,1,j]=dx[3,4,1,j]-RateSwitchDTG[1,4,1,j]*x[3,4,1,j]
      
      dx[4,1,1,j]=dx[4,1,1,j]-RateSwitchDTG[2,1,1,j]*x[4,1,1,j]
      dx[4,2,1,j]=dx[4,2,1,j]-RateSwitchDTG[2,2,1,j]*x[4,2,1,j]
      dx[4,3,1,j]=dx[4,3,1,j]-RateSwitchDTG[2,3,1,j]*x[4,3,1,j]
      dx[4,4,1,j]=dx[4,4,1,j]-RateSwitchDTG[2,4,1,j]*x[4,4,1,j]
      
      dx[5,1,1,j]=dx[5,1,1,j]-RateSwitchDTG[3,1,1,j]*x[5,1,1,j]
      dx[5,2,1,j]=dx[5,2,1,j]-RateSwitchDTG[3,2,1,j]*x[5,2,1,j]
      dx[5,3,1,j]=dx[5,3,1,j]-RateSwitchDTG[3,3,1,j]*x[5,3,1,j]
      dx[5,4,1,j]=dx[5,4,1,j]-RateSwitchDTG[3,4,1,j]*x[5,4,1,j]
      
      #NNRTI resistant
      dx[3,1,2,j]=dx[3,1,2,j]-RateSwitchDTG[1,1,2,j]*x[3,1,2,j]
      dx[3,2,2,j]=dx[3,2,2,j]-RateSwitchDTG[1,2,2,j]*x[3,2,2,j]
      dx[3,3,2,j]=dx[3,3,2,j]-RateSwitchDTG[1,3,2,j]*x[3,3,2,j]
      dx[3,4,2,j]=dx[3,4,2,j]-RateSwitchDTG[1,4,2,j]*x[3,4,2,j]
      
      dx[4,1,2,j]=dx[4,1,2,j]-RateSwitchDTG[2,1,2,j]*x[4,1,2,j]
      dx[4,2,2,j]=dx[4,2,2,j]-RateSwitchDTG[2,2,2,j]*x[4,2,2,j]
      dx[4,3,2,j]=dx[4,3,2,j]-RateSwitchDTG[2,3,2,j]*x[4,3,2,j]
      dx[4,4,2,j]=dx[4,4,2,j]-RateSwitchDTG[2,4,2,j]*x[4,4,2,j]
      
      dx[5,1,2,j]=dx[5,1,2,j]-RateSwitchDTG[3,1,2,j]*x[5,1,2,j]
      dx[5,2,2,j]=dx[5,2,2,j]-RateSwitchDTG[3,2,2,j]*x[5,2,2,j]
      dx[5,3,2,j]=dx[5,3,2,j]-RateSwitchDTG[3,3,2,j]*x[5,3,2,j]
      dx[5,4,2,j]=dx[5,4,2,j]-RateSwitchDTG[3,4,2,j]*x[5,4,2,j]
      
      #DTG compartments-------------------------------------
      #NNRTI susceptible
      dx[9,1,1,j]=dx[9,1,1,j]+RateTreatFirstDTG[1,1,j]*x[2,1,1,j]+RateSwitchDTG[3,1,1,j]*x[5,1,1,j]+RateSwitchDTG[1,1,1,j]*x[3,1,1,j]-(RateSuppFirstDTG[1]+RateStageTreatFirstR[1]+mu[3,1,1])*x[9,1,1,j]+RateStageTreatFirstL[1]*x[9,2,1,j]
      dx[10,1,1,j]=RateSuppFirstDTG[1]*x[9,1,1,j]-(RateFailFirstDTG[1]+mu[4,1,1])*x[10,1,1,j]+RateStageSuppFirst[1]*x[10,2,1,j]+RateSwitchDTG[2,1,1,j]*x[4,1,1,j]
      dx[11,1,1,j]=RateFailFirstDTG[1]*x[10,1,1,j]-(RateTreatSecond_noDTG[1,j]+RateStageFailFirst[1]+mu[5,1,1])*x[11,1,1,j]
      dx[12,1,1,j]=dx[12,1,1,j]+RateDiag[1,j]*((j==2)*(1-p_women_dtg))*x[1,1,1,j]-(RateTreatFirst_noDTG[1,1,j]+RateStageDiag[1]+mu[2,1,1])*x[12,1,1,j]
      dx[13,1,1,j]=dx[13,1,1,j]+RateTreatFirst_noDTG[1,1,j]*x[12,1,1,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,1])*x[13,1,1,j]+RateStageTreatFirstL[1]*x[13,2,1,j]
      dx[14,1,1,j]=dx[14,1,1,j]+RateSuppFirstSusc[1]*x[13,1,1,j]-(RateFailFirstSusc[1]+mu[4,1,1])*x[14,1,1,j]+RateStageSuppFirst[1]*x[14,2,1,j]
      dx[15,1,1,j]=dx[15,1,1,j]+RateFailFirstSusc[1]*x[14,1,1,j]-(RateTreatSecond_noDTG[1,j]+RateStageFailFirst[1]+mu[5,1,1])*x[15,1,1,j]
      
      dx[9,2,1,j]=dx[9,2,1,j]+RateTreatFirstDTG[1,2,j]*x[2,2,1,j]+RateSwitchDTG[3,2,1,j]*x[5,2,1,j]+RateSwitchDTG[1,2,1,j]*x[3,2,1,j]-(RateSuppFirstDTG[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,1])*x[9,2,1,j]+RateStageTreatFirstR[1]*x[9,1,1,j]+RateStageTreatFirstL[2]*x[9,3,1,j]
      dx[10,2,1,j]=RateSuppFirstDTG[2]*x[9,2,1,j]-(RateFailFirstDTG[2]+RateStageSuppFirst[1]+mu[4,2,1])*x[10,2,1,j]+RateStageSuppFirst[2]*x[10,3,1,j]+RateSwitchDTG[2,2,1,j]*x[4,2,1,j]
      dx[11,2,1,j]=RateFailFirstDTG[2]*x[10,2,1,j]-(RateTreatSecond_noDTG[2,j]+RateStageFailFirst[2]+mu[5,2,1])*x[11,2,1,j]+RateStageFailFirst[1]*x[11,1,1,j]
      dx[12,2,1,j]=dx[12,2,1,j]+RateDiag[2,j]*((j==2)*(1-p_women_dtg))*x[1,2,1,j]-(RateTreatFirst_noDTG[1,2,j]+RateStageDiag[2]+mu[2,2,1])*x[12,2,1,j]+RateStageDiag[1]*x[12,1,1,j]
      dx[13,2,1,j]=dx[13,2,1,j]+RateTreatFirst_noDTG[1,2,j]*x[12,2,1,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,1])*x[13,2,1,j]+RateStageTreatFirstR[1]*x[13,1,1,j]+RateStageTreatFirstL[2]*x[13,3,1,j]
      dx[14,2,1,j]=dx[14,2,1,j]+RateSuppFirstSusc[2]*x[13,2,1,j]-(RateFailFirstSusc[2]+RateStageSuppFirst[1]+mu[4,2,1])*x[14,2,1,j]+RateStageSuppFirst[2]*x[14,3,1,j]
      dx[15,2,1,j]=dx[15,2,1,j]+RateFailFirstSusc[2]*x[14,2,1,j]-(RateTreatSecond_noDTG[2,j]+RateStageFailFirst[2]+mu[5,2,1])*x[15,2,1,j]+RateStageFailFirst[1]*x[15,1,1,j]
      
      dx[9,3,1,j]=dx[9,3,1,j]+RateTreatFirstDTG[1,3,j]*x[2,3,1,j]+RateSwitchDTG[3,3,1,j]*x[5,3,1,j]+RateSwitchDTG[1,3,1,j]*x[3,3,1,j]-(RateSuppFirstDTG[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,1])*x[9,3,1,j]+RateStageTreatFirstR[2]*x[9,2,1,j]+RateStageTreatFirstL[3]*x[9,4,1,j]
      dx[10,3,1,j]=RateSuppFirstDTG[3]*x[9,3,1,j]-(RateFailFirstDTG[3]+RateStageSuppFirst[2]+mu[4,3,1])*x[10,3,1,j]+RateStageSuppFirst[3]*x[10,4,1,j]+RateSwitchDTG[2,3,1,j]*x[4,3,1,j]
      dx[11,3,1,j]=RateFailFirstDTG[3]*x[10,3,1,j]-(RateTreatSecond_noDTG[3,j]+RateStageFailFirst[3]+mu[5,3,1])*x[11,3,1,j]+RateStageFailFirst[2]*x[11,2,1,j]
      dx[12,3,1,j]=dx[12,3,1,j]+RateDiag[3,j]*((j==2)*(1-p_women_dtg))*x[1,3,1,j]-(RateTreatFirst_noDTG[1,3,j]+RateStageDiag[3]+mu[2,3,1])*x[12,3,1,j]+RateStageDiag[2]*x[12,2,1,j]
      dx[13,3,1,j]=dx[13,3,1,j]+RateTreatFirst_noDTG[1,3,j]*x[12,3,1,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,1])*x[13,3,1,j]+RateStageTreatFirstR[2]*x[13,2,1,j]+RateStageTreatFirstL[3]*x[13,4,1,j]
      dx[14,3,1,j]=dx[14,3,1,j]+RateSuppFirstSusc[3]*x[13,3,1,j]-(RateFailFirstSusc[3]+RateStageSuppFirst[2]+mu[4,3,1])*x[14,3,1,j]+RateStageSuppFirst[3]*x[14,4,1,j]
      dx[15,3,1,j]=dx[15,3,1,j]+RateFailFirstSusc[3]*x[14,3,1,j]-(RateTreatSecond_noDTG[3,j]+RateStageFailFirst[3]+mu[5,3,1])*x[15,3,1,j]+RateStageFailFirst[2]*x[15,2,1,j]
      
      dx[9,4,1,j]=dx[9,4,1,j]+RateTreatFirstDTG[1,4,j]*x[2,4,1,j]+RateSwitchDTG[3,4,1,j]*x[5,4,1,j]+RateSwitchDTG[1,4,1,j]*x[3,4,1,j]-(RateSuppFirstDTG[4]+RateStageTreatFirstL[3]+mu[3,4,1])*x[9,4,1,j]+RateStageTreatFirstR[3]*x[9,3,1,j]
      dx[10,4,1,j]=RateSuppFirstDTG[4]*x[9,4,1,j]-(RateFailFirstDTG[4]+RateStageSuppFirst[3]+mu[4,4,1])*x[10,4,1,j]+RateSwitchDTG[2,4,1,j]*x[4,4,1,j]
      dx[11,4,1,j]=RateFailFirstDTG[4]*x[10,4,1,j]-(RateTreatSecond_noDTG[4,j]+mu[5,4,1])*x[11,4,1,j]+RateStageFailFirst[3]*x[11,3,1,j]
      dx[12,4,1,j]=dx[12,4,1,j]+RateDiag[4,j]*((j==2)*(1-p_women_dtg))*x[1,4,1,j]-(RateTreatFirst_noDTG[1,4,j]+mu[2,4,1])*x[12,4,1,j]+RateStageDiag[3]*x[12,3,1,j]
      dx[13,4,1,j]=dx[13,4,1,j]+RateTreatFirst_noDTG[1,4,j]*x[12,4,1,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,1])*x[13,4,1,j]+RateStageTreatFirstR[3]*x[13,3,1,j]
      dx[14,4,1,j]=dx[14,4,1,j]+RateSuppFirstSusc[4]*x[13,4,1,j]-(RateFailFirstSusc[4]+RateStageSuppFirst[3]+mu[4,4,1])*x[14,4,1,j]
      dx[15,4,1,j]=dx[15,4,1,j]+RateFailFirstSusc[4]*x[14,4,1,j]-(RateTreatSecond_noDTG[4,j]+mu[5,4,1])*x[15,4,1,j]+RateStageFailFirst[3]*x[15,3,1,j]
      
      #NNRTI resistant
      dx[9,1,2,j]=dx[9,1,2,j]+RateTreatFirstDTG[2,1,j]*x[2,1,2,j]+RateSwitchDTG[3,1,2,j]*x[5,1,2,j]+RateSwitchDTG[1,1,2,j]*x[3,1,2,j]-(RateSuppFirstDTG[1]+RateStageTreatFirstR[1]+mu[3,1,2])*x[9,1,2,j]+RateStageTreatFirstL[1]*x[9,2,2,j]
      dx[10,1,2,j]=RateSuppFirstDTG[1]*x[9,1,2,j]-(RateFailFirstDTG[1]+mu[4,1,2])*x[10,1,2,j]+RateStageSuppFirst[1]*x[10,2,2,j]+RateSwitchDTG[2,1,2,j]*x[4,1,2,j]
      dx[11,1,2,j]=RateFailFirstDTG[1]*x[10,1,2,j]-(RateTreatSecond_noDTG[1,j]+RateStageFailFirst[1]+mu[5,1,2])*x[11,1,2,j]
      dx[12,1,2,j]=dx[12,1,2,j]+RateDiag[1,j]*((j==2)*(1-p_women_dtg))*x[1,1,2,j]-(RateTreatFirst_noDTG[2,1,j]+RateStageDiag[1]+mu[2,1,2])*x[12,1,2,j]
      dx[13,1,2,j]=dx[13,1,2,j]+RateTreatFirst_noDTG[2,1,j]*x[12,1,2,j]-(RateSuppFirstResis[1]+RateStageTreatFirstR[1]+mu[3,1,2])*x[13,1,2,j]+RateStageTreatFirstL[1]*x[13,2,2,j]
      dx[14,1,2,j]=dx[14,1,2,j]+RateSuppFirstResis[1]*x[13,1,2,j]-(RateFailFirstResis[1]+mu[4,1,2])*x[14,1,2,j]+RateStageSuppFirst[1]*x[14,2,2,j]
      dx[15,1,2,j]=dx[15,1,2,j]+RateFailFirstResis[1]*x[14,1,2,j]-(RateTreatSecond_noDTG[1,j]+RateStageFailFirst[1]+mu[5,1,2])*x[15,1,2,j]
      
      dx[9,2,2,j]=dx[9,2,2,j]+RateTreatFirstDTG[2,2,j]*x[2,2,2,j]+RateSwitchDTG[3,2,2,j]*x[5,2,2,j]+RateSwitchDTG[1,2,2,j]*x[3,2,2,j]-(RateSuppFirstDTG[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,2])*x[9,2,2,j]+RateStageTreatFirstR[1]*x[9,1,2,j]+RateStageTreatFirstL[2]*x[9,3,2,j]
      dx[10,2,2,j]=RateSuppFirstDTG[2]*x[9,2,2,j]-(RateFailFirstDTG[2]+RateStageSuppFirst[1]+mu[4,2,2])*x[10,2,2,j]+RateStageSuppFirst[2]*x[10,3,2,j]+RateSwitchDTG[2,2,2,j]*x[4,2,2,j]
      dx[11,2,2,j]=RateFailFirstDTG[2]*x[10,2,2,j]-(RateTreatSecond_noDTG[2,j]+RateStageFailFirst[2]+mu[5,2,2])*x[11,2,2,j]+RateStageFailFirst[1]*x[11,1,2,j]
      dx[12,2,2,j]=dx[12,2,2,j]+RateDiag[2,j]*((j==2)*(1-p_women_dtg))*x[1,2,2,j]-(RateTreatFirst_noDTG[2,2,j]+RateStageDiag[2]+mu[2,2,2])*x[12,2,2,j]+RateStageDiag[1]*x[12,1,2,j]
      dx[13,2,2,j]=dx[13,2,2,j]+RateTreatFirst_noDTG[2,2,j]*x[12,2,2,j]-(RateSuppFirstResis[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,2])*x[13,2,2,j]+RateStageTreatFirstR[1]*x[13,1,2,j]+RateStageTreatFirstL[2]*x[13,3,2,j]
      dx[14,2,2,j]=dx[14,2,2,j]+RateSuppFirstResis[2]*x[13,2,2,j]-(RateFailFirstResis[2]+RateStageSuppFirst[1]+mu[4,2,2])*x[14,2,2,j]+RateStageSuppFirst[2]*x[14,3,2,j]
      dx[15,2,2,j]=dx[15,2,2,j]+RateFailFirstResis[2]*x[14,2,2,j]-(RateTreatSecond_noDTG[2,j]+RateStageFailFirst[2]+mu[5,2,2])*x[15,2,2,j]+RateStageFailFirst[1]*x[15,1,2,j]
      
      dx[9,3,2,j]=dx[9,3,2,j]+RateTreatFirstDTG[2,3,j]*x[2,3,2,j]+RateSwitchDTG[3,3,2,j]*x[5,3,2,j]+RateSwitchDTG[1,3,2,j]*x[3,3,2,j]-(RateSuppFirstDTG[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,2])*x[9,3,2,j]+RateStageTreatFirstR[2]*x[9,2,2,j]+RateStageTreatFirstL[3]*x[9,4,2,j]
      dx[10,3,2,j]=RateSuppFirstDTG[3]*x[9,3,2,j]-(RateFailFirstDTG[3]+RateStageSuppFirst[2]+mu[4,3,2])*x[10,3,2,j]+RateStageSuppFirst[3]*x[10,4,2,j]+RateSwitchDTG[2,3,2,j]*x[4,3,2,j]
      dx[11,3,2,j]=RateFailFirstDTG[3]*x[10,3,2,j]-(RateTreatSecond_noDTG[3,j]+RateStageFailFirst[3]+mu[5,3,2])*x[11,3,2,j]+RateStageFailFirst[2]*x[11,2,2,j]
      dx[12,3,2,j]=dx[12,3,2,j]+RateDiag[3,j]*((j==2)*(1-p_women_dtg))*x[1,3,2,j]-(RateTreatFirst_noDTG[2,3,j]+RateStageDiag[3]+mu[2,3,2])*x[12,3,2,j]+RateStageDiag[2]*x[12,2,2,j]
      dx[13,3,2,j]=dx[13,3,2,j]+RateTreatFirst_noDTG[2,3,j]*x[12,3,2,j]-(RateSuppFirstResis[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,2])*x[13,3,2,j]+RateStageTreatFirstR[2]*x[13,2,2,j]+RateStageTreatFirstL[3]*x[13,4,2,j]
      dx[14,3,2,j]=dx[14,3,2,j]+RateSuppFirstResis[3]*x[13,3,2,j]-(RateFailFirstResis[3]+RateStageSuppFirst[2]+mu[4,3,2])*x[14,3,2,j]+RateStageSuppFirst[3]*x[14,4,2,j]
      dx[15,3,2,j]=dx[15,3,2,j]+RateFailFirstResis[3]*x[14,3,2,j]-(RateTreatSecond_noDTG[3,j]+RateStageFailFirst[3]+mu[5,3,2])*x[15,3,2,j]+RateStageFailFirst[2]*x[15,2,2,j]
      
      dx[9,4,2,j]=dx[9,4,2,j]+RateTreatFirstDTG[2,4,j]*x[2,4,2,j]+RateSwitchDTG[3,4,2,j]*x[5,4,2,j]+RateSwitchDTG[1,4,2,j]*x[3,4,2,j]-(RateSuppFirstDTG[4]+RateStageTreatFirstL[3]+mu[3,4,2])*x[9,4,2,j]+RateStageTreatFirstR[3]*x[9,3,2,j]
      dx[10,4,2,j]=RateSuppFirstDTG[4]*x[9,4,2,j]-(RateFailFirstDTG[4]+RateStageSuppFirst[3]+mu[4,4,2])*x[10,4,2,j]+RateSwitchDTG[2,4,2,j]*x[4,4,2,j]
      dx[11,4,2,j]=RateFailFirstDTG[4]*x[10,4,2,j]-(RateTreatSecond_noDTG[4,j]+mu[5,4,2])*x[11,4,2,j]+RateStageFailFirst[3]*x[11,3,2,j]
      dx[12,4,2,j]=dx[12,4,2,j]+RateDiag[4,j]*((j==2)*(1-p_women_dtg))*x[1,4,2,j]-(RateTreatFirst_noDTG[2,4,j]+mu[2,4,2])*x[12,4,2,j]+RateStageDiag[3]*x[12,3,2,j]
      dx[13,4,2,j]=dx[13,4,2,j]+RateTreatFirst_noDTG[2,4,j]*x[12,4,2,j]-(RateSuppFirstResis[4]+RateStageTreatFirstL[3]+mu[3,4,2])*x[13,4,2,j]+RateStageTreatFirstR[3]*x[13,3,2,j]
      dx[14,4,2,j]=dx[14,4,2,j]+RateSuppFirstResis[4]*x[13,4,2,j]-(RateFailFirstResis[4]+RateStageSuppFirst[3]+mu[4,4,2])*x[14,4,2,j]
      dx[15,4,2,j]=dx[15,4,2,j]+RateFailFirstResis[4]*x[14,4,2,j]-(RateTreatSecond_noDTG[4,j]+mu[5,4,2])*x[15,4,2,j]+RateStageFailFirst[3]*x[15,3,2,j]
      
      
      
      # #DTG compartments-------------------------------------
      # #NNRTI susceptible
      # dx[9,1,1,j]=dx[9,1,1,j]+RateTreatFirstDTG[1,1,j]*x[2,1,1,j]+RateSwitchDTG[1,1,1,j]*x[3,1,1,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,1])*x[9,1,1,j]+RateStageTreatFirstL[1]*x[9,2,1,j]
      # dx[10,1,1,j]=RateSwitchDTG[2,1,1,j]*x[4,1,1,j]+RateSuppFirstSusc[1]*x[9,1,1,j]-(RateFailFirstSusc[1]+mu[4,1,1])*x[10,1,1,j]+RateStageSuppFirst[1]*x[10,2,1,j]
      # dx[11,1,1,j]=RateSwitchDTG[3,1,1,j]*x[5,1,1,j]+RateFailFirstSusc[1]*x[10,1,1,j]-(RateStageFailFirst[1]+mu[5,1,1])*x[11,1,1,j]
      # 
      # dx[9,2,1,j]=dx[9,2,1,j]+RateTreatFirstDTG[1,2,j]*x[2,2,1,j]+RateSwitchDTG[1,2,1,j]*x[3,2,1,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,1])*x[9,2,1,j]+RateStageTreatFirstR[1]*x[9,1,1,j]+RateStageTreatFirstL[2]*x[9,3,1,j]
      # dx[10,2,1,j]=RateSwitchDTG[2,2,1,j]*x[4,2,1,j]+RateSuppFirstSusc[2]*x[9,2,1,j]-(RateFailFirstSusc[2]+RateStageSuppFirst[1]+mu[4,2,1])*x[10,2,1,j]+RateStageSuppFirst[2]*x[10,3,1,j]
      # dx[11,2,1,j]=RateSwitchDTG[3,2,1,j]*x[5,2,1,j]+RateFailFirstSusc[2]*x[10,2,1,j]-(RateStageFailFirst[2]+mu[5,2,1])*x[11,2,1,j]+RateStageFailFirst[1]*x[11,1,1,j]
      # 
      # dx[9,3,1,j]=dx[9,3,1,j]+RateTreatFirstDTG[1,3,j]*x[2,3,1,j]+RateSwitchDTG[1,3,1,j]*x[3,3,1,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,1])*x[9,3,1,j]+RateStageTreatFirstR[2]*x[9,2,1,j]+RateStageTreatFirstL[3]*x[9,4,1,j]
      # dx[10,3,1,j]=RateSwitchDTG[2,3,1,j]*x[4,3,1,j]+RateSuppFirstSusc[3]*x[9,3,1,j]-(RateFailFirstSusc[3]+RateStageSuppFirst[2]+mu[4,3,1])*x[10,3,1,j]+RateStageSuppFirst[3]*x[10,4,1,j]
      # dx[11,3,1,j]=RateSwitchDTG[3,3,1,j]*x[5,3,1,j]+RateFailFirstSusc[3]*x[10,3,1,j]-(RateStageFailFirst[3]+mu[5,3,1])*x[11,3,1,j]+RateStageFailFirst[2]*x[11,2,1,j]
      # 
      # dx[9,4,1,j]=dx[9,4,1,j]+RateTreatFirstDTG[1,4,j]*x[2,4,1,j]+RateSwitchDTG[1,4,1,j]*x[3,4,1,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,1])*x[9,4,1,j]+RateStageTreatFirstR[3]*x[9,3,1,j]
      # dx[10,4,1,j]=RateSwitchDTG[2,4,1,j]*x[4,4,1,j]+RateSuppFirstSusc[4]*x[9,4,1,j]-(RateFailFirstSusc[4]+RateStageSuppFirst[3]+mu[4,4,1])*x[10,4,1,j]
      # dx[11,4,1,j]=RateSwitchDTG[3,4,1,j]*x[5,4,1,j]+RateFailFirstSusc[4]*x[10,4,1,j]-(mu[5,4,1])*x[11,4,1,j]+RateStageFailFirst[3]*x[11,3,1,j]
      # ##############################################################################################################################################################
      
      dx[9,1,1,j]=dx[9,1,1,j]-RateStopTreatFirst[1]*x[9,1,1,j]-RateTreatToFailFirstDTG[1]*x[9,1,1,j]
      dx[10,1,1,j]=dx[10,1,1,j]-RateStopSuppFirst[1]*x[10,1,1,j]+RateFailToSuppTreatFirstDTG[1]*x[11,1,1,j]
      dx[11,1,1,j]=dx[11,1,1,j]-RateStopFailFirst[1]*x[11,1,1,j]-RateFailToSuppTreatFirstDTG[1]*x[11,1,1,j]+RateTreatToFailFirstDTG[1]*x[9,1,1,j]
      dx[12,1,1,j]=dx[12,1,1,j]-RateDirectTreatSecond[1,1]*x[12,1,1,j]+RateStopTreatFirst[1]*x[13,1,1,j]+RateStopSuppFirst[1]*x[14,1,1,j]+RateStopFailFirst[1]*x[15,1,1,j]+
        (RateStopTreatSecond[1]*x[6,1,1,j]+RateStopSuppSecond[1]*x[7,1,1,j]+RateStopFailSecond[1]*x[8,1,1,j])*((j==2)*(1-p_women_dtg))
      dx[13,1,1,j]=dx[13,1,1,j]-RateStopTreatFirst[1]*x[13,1,1,j]-RateTreatToFailFirstSusc[1]*x[13,1,1,j]
      dx[14,1,1,j]=dx[14,1,1,j]-RateStopSuppFirst[1]*x[14,1,1,j]+RateFailToSuppTreatFirstSusc[1]*x[15,1,1,j]-RateSuppFirstToSecond[1]*x[14,1,1,j]
      dx[15,1,1,j]=dx[15,1,1,j]-RateStopFailFirst[1]*x[15,1,1,j]-RateFailToSuppTreatFirstSusc[1]*x[15,1,1,j]+RateTreatToFailFirstSusc[1]*x[13,1,1,j]
      
      dx[9,2,1,j]=dx[9,2,1,j]-RateStopTreatFirst[2]*x[9,2,1,j]-RateTreatToFailFirstDTG[2]*x[9,2,1,j]
      dx[10,2,1,j]=dx[10,2,1,j]-RateStopSuppFirst[2]*x[10,2,1,j]+RateFailToSuppTreatFirstDTG[2]*x[11,2,1,j]-RateSuppFirstToSecond[2]*x[10,2,1,j]
      dx[11,2,1,j]=dx[11,2,1,j]-RateStopFailFirst[2]*x[11,2,1,j]-RateFailToSuppTreatFirstDTG[2]*x[11,2,1,j]+RateTreatToFailFirstDTG[2]*x[9,2,1,j]
      dx[12,2,1,j]=dx[12,2,1,j]-RateDirectTreatSecond[1,2]*x[12,2,1,j]+RateStopTreatFirst[2]*x[13,2,1,j]+RateStopSuppFirst[2]*x[14,2,1,j]+RateStopFailFirst[2]*x[15,2,1,j]+
        (RateStopTreatSecond[2]*x[6,2,1,j]+RateStopSuppSecond[2]*x[7,2,1,j]+RateStopFailSecond[2]*x[8,2,1,j])*((j==2)*(1-p_women_dtg))
      dx[13,2,1,j]=dx[13,2,1,j]-RateStopTreatFirst[2]*x[13,2,1,j]-RateTreatToFailFirstSusc[2]*x[13,2,1,j]
      dx[14,2,1,j]=dx[14,2,1,j]-RateStopSuppFirst[2]*x[14,2,1,j]+RateFailToSuppTreatFirstSusc[2]*x[15,2,1,j]-RateSuppFirstToSecond[2]*x[14,2,1,j]
      dx[15,2,1,j]=dx[15,2,1,j]-RateStopFailFirst[2]*x[15,2,1,j]-RateFailToSuppTreatFirstSusc[2]*x[15,2,1,j]+RateTreatToFailFirstSusc[2]*x[13,2,1,j]
      
      dx[9,3,1,j]=dx[9,3,1,j]-RateStopTreatFirst[3]*x[9,3,1,j]-RateTreatToFailFirstDTG[3]*x[9,3,1,j]
      dx[10,3,1,j]=dx[10,3,1,j]-RateStopSuppFirst[3]*x[10,3,1,j]+RateFailToSuppTreatFirstDTG[3]*x[11,3,1,j]
      dx[11,3,1,j]=dx[11,3,1,j]-RateStopFailFirst[3]*x[11,3,1,j]-RateFailToSuppTreatFirstDTG[3]*x[11,3,1,j]+RateTreatToFailFirstDTG[3]*x[9,3,1,j]
      dx[12,3,1,j]=dx[12,3,1,j]-RateDirectTreatSecond[1,3]*x[12,3,1,j]+RateStopTreatFirst[3]*x[13,3,1,j]+RateStopSuppFirst[3]*x[14,3,1,j]+RateStopFailFirst[3]*x[15,3,1,j]+
        (RateStopTreatSecond[3]*x[6,3,1,j]+RateStopSuppSecond[3]*x[7,3,1,j]+RateStopFailSecond[3]*x[8,3,1,j])*((j==2)*(1-p_women_dtg))
      dx[13,3,1,j]=dx[13,3,1,j]-RateStopTreatFirst[3]*x[13,3,1,j]-RateTreatToFailFirstSusc[3]*x[13,3,1,j]
      dx[14,3,1,j]=dx[14,3,1,j]-RateStopSuppFirst[3]*x[14,3,1,j]+RateFailToSuppTreatFirstSusc[3]*x[15,3,1,j]-RateSuppFirstToSecond[3]*x[14,3,1,j]
      dx[15,3,1,j]=dx[15,3,1,j]-RateStopFailFirst[3]*x[15,3,1,j]-RateFailToSuppTreatFirstSusc[3]*x[15,3,1,j]+RateTreatToFailFirstSusc[3]*x[13,3,1,j]
      
      dx[9,4,1,j]=dx[9,4,1,j]-RateStopTreatFirst[4]*x[9,4,1,j]-RateTreatToFailFirstDTG[4]*x[9,4,1,j]
      dx[10,4,1,j]=dx[10,4,1,j]-RateStopSuppFirst[4]*x[10,4,1,j]+RateFailToSuppTreatFirstDTG[4]*x[11,4,1,j]
      dx[11,4,1,j]=dx[11,4,1,j]-RateStopFailFirst[4]*x[11,4,1,j]-RateFailToSuppTreatFirstDTG[4]*x[11,4,1,j]+RateTreatToFailFirstDTG[4]*x[9,4,1,j]
      dx[12,4,1,j]=dx[12,4,1,j]-RateDirectTreatSecond[1,4]*x[12,4,1,j]+RateStopTreatFirst[4]*x[13,4,1,j]+RateStopSuppFirst[4]*x[14,4,1,j]+RateStopFailFirst[4]*x[15,4,1,j]+
        (RateStopTreatSecond[4]*x[6,4,1,j]+RateStopSuppSecond[4]*x[7,4,1,j]+RateStopFailSecond[4]*x[8,4,1,j])*((j==2)*(1-p_women_dtg))
      dx[13,4,1,j]=dx[13,4,1,j]-RateStopTreatFirst[4]*x[13,4,1,j]-RateTreatToFailFirstSusc[4]*x[13,4,1,j]
      dx[14,4,1,j]=dx[14,4,1,j]-RateStopSuppFirst[4]*x[14,4,1,j]+RateFailToSuppTreatFirstSusc[4]*x[15,4,1,j]-RateSuppFirstToSecond[4]*x[14,4,1,j]
      dx[15,4,1,j]=dx[15,4,1,j]-RateStopFailFirst[4]*x[15,4,1,j]-RateFailToSuppTreatFirstSusc[4]*x[15,4,1,j]+RateTreatToFailFirstSusc[4]*x[13,4,1,j]
      
      #############################################################################################################################################################################
      #############################################################################################################################################################################
      #NNRTI resistant
      # dx[9,1,2,j]=dx[9,1,2,j]+RateTreatFirstDTG[2,1,j]*x[2,1,2,j]+RateSwitchDTG[1,1,2,j]*x[3,1,2,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,2])*x[9,1,2,j]+RateStageTreatFirstL[1]*x[9,2,2,j]
      # dx[10,1,2,j]=RateSwitchDTG[2,1,2,j]*x[4,1,2,j]+RateSuppFirstSusc[1]*x[9,1,2,j]-(RateFailFirstSusc[1]+mu[4,1,2])*x[10,1,2,j]+RateStageSuppFirst[1]*x[10,2,2,j]
      # dx[11,1,2,j]=RateSwitchDTG[3,1,2,j]*x[5,1,2,j]+RateFailFirstSusc[1]*x[10,1,2,j]-(RateStageFailFirst[1]+mu[5,1,2])*x[11,1,2,j]
      # 
      # dx[9,2,2,j]=dx[9,2,2,j]+RateTreatFirstDTG[2,2,j]*x[2,2,2,j]+RateSwitchDTG[1,2,2,j]*x[3,2,2,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,2])*x[9,2,2,j]+RateStageTreatFirstR[1]*x[9,1,2,j]+RateStageTreatFirstL[2]*x[9,3,2,j]
      # dx[10,2,2,j]=RateSwitchDTG[2,2,2,j]*x[4,2,2,j]+RateSuppFirstSusc[2]*x[9,2,2,j]-(RateFailFirstSusc[2]+RateStageSuppFirst[1]+mu[4,2,2])*x[10,2,2,j]+RateStageSuppFirst[2]*x[10,3,2,j]
      # dx[11,2,2,j]=RateSwitchDTG[3,2,2,j]*x[5,2,2,j]+RateFailFirstSusc[2]*x[10,2,2,j]-(RateStageFailFirst[2]+mu[5,2,2])*x[11,2,2,j]+RateStageFailFirst[1]*x[11,1,2,j]
      # 
      # dx[9,3,2,j]=dx[9,3,2,j]+RateTreatFirstDTG[2,3,j]*x[2,3,2,j]+RateSwitchDTG[1,3,2,j]*x[3,3,2,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,2])*x[9,3,2,j]+RateStageTreatFirstR[2]*x[9,2,2,j]+RateStageTreatFirstL[3]*x[9,4,2,j]
      # dx[10,3,2,j]=RateSwitchDTG[2,3,2,j]*x[4,3,2,j]+RateSuppFirstSusc[3]*x[9,3,2,j]-(RateFailFirstSusc[3]+RateStageSuppFirst[2]+mu[4,3,2])*x[10,3,2,j]+RateStageSuppFirst[3]*x[10,4,2,j]
      # dx[11,3,2,j]=RateSwitchDTG[3,3,2,j]*x[5,3,2,j]+RateFailFirstSusc[3]*x[10,3,2,j]-(RateStageFailFirst[3]+mu[5,3,2])*x[11,3,2,j]+RateStageFailFirst[2]*x[11,2,2,j]
      # 
      # dx[9,4,2,j]=dx[9,4,2,j]+RateTreatFirstDTG[2,4,j]*x[2,4,2,j]+RateSwitchDTG[1,4,2,j]*x[3,4,2,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,2])*x[9,4,2,j]+RateStageTreatFirstR[3]*x[9,3,2,j]
      # dx[10,4,2,j]=RateSwitchDTG[2,4,2,j]*x[4,4,2,j]+RateSuppFirstSusc[4]*x[9,4,2,j]-(RateFailFirstSusc[4]+RateStageSuppFirst[3]+mu[4,4,2])*x[10,4,2,j]
      # dx[11,4,2,j]=RateSwitchDTG[3,4,2,j]*x[5,4,2,j]+RateFailFirstSusc[4]*x[10,4,2,j]-(mu[5,4,2])*x[11,4,2,j]+RateStageFailFirst[3]*x[11,3,2,j]
      ##############################################################################################################################################################
      
      dx[9,1,2,j]=dx[9,1,2,j]-RateStopTreatFirst[1]*x[9,1,2,j]-RateTreatToFailFirstDTG[1]*x[9,1,2,j]
      dx[10,1,2,j]=dx[10,1,2,j]-RateStopSuppFirst[1]*x[10,1,2,j]+RateFailToSuppTreatFirstDTG[1]*x[11,1,2,j]
      dx[11,1,2,j]=dx[11,1,2,j]-RateStopFailFirst[1]*x[11,1,2,j]-RateFailToSuppTreatFirstDTG[1]*x[11,1,2,j]+RateTreatToFailFirstDTG[1]*x[9,1,2,j]
      dx[12,1,2,j]=dx[12,1,2,j]-RateDirectTreatSecond[2,1]*x[12,1,2,j]+RateStopTreatFirst[1]*x[13,1,2,j]+RateStopSuppFirst[1]*x[14,1,2,j]+RateStopFailFirst[1]*x[15,1,2,j]+
        (RateStopTreatSecond[1]*x[6,1,2,j]+RateStopSuppSecond[1]*x[7,1,2,j]+RateStopFailSecond[1]*x[8,1,2,j])*((j==2)*(1-p_women_dtg))
      dx[13,1,2,j]=dx[13,1,2,j]-RateStopTreatFirst[1]*x[13,1,2,j]-RateTreatToFailFirstResis[1]*x[13,1,2,j]
      dx[14,1,2,j]=dx[14,1,2,j]-RateStopSuppFirst[1]*x[14,1,2,j]+RateFailToSuppTreatFirstResis[1]*x[15,1,2,j]-RateSuppFirstToSecond[1]*x[14,1,2,j]
      dx[15,1,2,j]=dx[15,1,2,j]-RateStopFailFirst[1]*x[15,1,2,j]-RateFailToSuppTreatFirstResis[1]*x[15,1,2,j]+RateTreatToFailFirstResis[1]*x[13,1,2,j]
      
      dx[9,2,2,j]=dx[9,2,2,j]-RateStopTreatFirst[2]*x[9,2,2,j]-RateTreatToFailFirstDTG[2]*x[9,2,2,j]
      dx[10,2,2,j]=dx[10,2,2,j]-RateStopSuppFirst[2]*x[10,2,2,j]+RateFailToSuppTreatFirstDTG[2]*x[11,2,2,j]
      dx[11,2,2,j]=dx[11,2,2,j]-RateStopFailFirst[2]*x[11,2,2,j]-RateFailToSuppTreatFirstDTG[2]*x[11,2,2,j]+RateTreatToFailFirstDTG[2]*x[9,2,2,j]
      dx[12,2,2,j]=dx[12,2,2,j]-RateDirectTreatSecond[2,2]*x[12,2,2,j]+RateStopTreatFirst[2]*x[13,2,2,j]+RateStopSuppFirst[2]*x[14,2,2,j]+RateStopFailFirst[2]*x[15,2,2,j]+
        (RateStopTreatSecond[2]*x[6,2,2,j]+RateStopSuppSecond[2]*x[7,2,2,j]+RateStopFailSecond[2]*x[8,2,2,j])*((j==2)*(1-p_women_dtg))
      dx[13,2,2,j]=dx[13,2,2,j]-RateStopTreatFirst[2]*x[13,2,2,j]-RateTreatToFailFirstResis[2]*x[13,2,2,j]
      dx[14,2,2,j]=dx[14,2,2,j]-RateStopSuppFirst[2]*x[14,2,2,j]+RateFailToSuppTreatFirstResis[2]*x[15,2,2,j]-RateSuppFirstToSecond[2]*x[14,2,2,j]
      dx[15,2,2,j]=dx[15,2,2,j]-RateStopFailFirst[2]*x[15,2,2,j]-RateFailToSuppTreatFirstResis[2]*x[15,2,2,j]+RateTreatToFailFirstResis[2]*x[13,2,2,j]
      
      dx[9,3,2,j]=dx[9,3,2,j]-RateStopTreatFirst[3]*x[9,3,2,j]-RateTreatToFailFirstDTG[3]*x[9,3,2,j]
      dx[10,3,2,j]=dx[10,3,2,j]-RateStopSuppFirst[3]*x[10,3,2,j]+RateFailToSuppTreatFirstDTG[3]*x[11,3,2,j]
      dx[11,3,2,j]=dx[11,3,2,j]-RateStopFailFirst[3]*x[11,3,2,j]-RateFailToSuppTreatFirstDTG[3]*x[11,3,2,j]+RateTreatToFailFirstDTG[3]*x[9,3,2,j]
      dx[12,3,2,j]=dx[12,3,2,j]-RateDirectTreatSecond[2,3]*x[12,3,2,j]+RateStopTreatFirst[3]*x[13,3,2,j]+RateStopSuppFirst[3]*x[14,3,2,j]+RateStopFailFirst[3]*x[15,3,2,j]+
        (RateStopTreatSecond[3]*x[6,3,2,j]+RateStopSuppSecond[3]*x[7,3,2,j]+RateStopFailSecond[3]*x[8,3,2,j])*((j==2)*(1-p_women_dtg))
      dx[13,3,2,j]=dx[13,3,2,j]-RateStopTreatFirst[3]*x[13,3,2,j]-RateTreatToFailFirstResis[3]*x[13,3,2,j]
      dx[14,3,2,j]=dx[14,3,2,j]-RateStopSuppFirst[3]*x[14,3,2,j]+RateFailToSuppTreatFirstResis[3]*x[15,3,2,j]-RateSuppFirstToSecond[3]*x[14,3,2,j]
      dx[15,3,2,j]=dx[15,3,2,j]-RateStopFailFirst[3]*x[15,3,2,j]-RateFailToSuppTreatFirstResis[3]*x[15,3,2,j]+RateTreatToFailFirstResis[3]*x[13,3,2,j]
      
      dx[9,4,2,j]=dx[9,4,2,j]-RateStopTreatFirst[4]*x[9,4,2,j]-RateTreatToFailFirstDTG[4]*x[9,4,2,j]
      dx[10,4,2,j]=dx[10,4,2,j]-RateStopSuppFirst[4]*x[10,4,2,j]+RateFailToSuppTreatFirstDTG[4]*x[11,4,2,j]
      dx[11,4,2,j]=dx[11,4,2,j]-RateStopFailFirst[4]*x[11,4,2,j]-RateFailToSuppTreatFirstDTG[4]*x[11,4,2,j]+RateTreatToFailFirstDTG[4]*x[9,4,2,j]
      dx[12,4,2,j]=dx[12,4,2,j]-RateDirectTreatSecond[2,4]*x[12,4,2,j]+RateStopTreatFirst[4]*x[13,4,2,j]+RateStopSuppFirst[4]*x[14,4,2,j]+RateStopFailFirst[4]*x[15,4,2,j]+
        (RateStopTreatSecond[4]*x[6,4,2,j]+RateStopSuppSecond[4]*x[7,4,2,j]+RateStopFailSecond[4]*x[8,4,2,j])*((j==2)*(1-p_women_dtg))
      dx[13,4,2,j]=dx[13,4,2,j]-RateStopTreatFirst[4]*x[13,4,2,j]-RateTreatToFailFirstResis[4]*x[13,4,2,j]
      dx[14,4,2,j]=dx[14,4,2,j]-RateStopSuppFirst[4]*x[14,4,2,j]+RateFailToSuppTreatFirstResis[4]*x[15,4,2,j]-RateSuppFirstToSecond[4]*x[14,4,2,j]
      dx[15,4,2,j]=dx[15,4,2,j]-RateStopFailFirst[4]*x[15,4,2,j]-RateFailToSuppTreatFirstResis[4]*x[15,4,2,j]+RateTreatToFailFirstResis[4]*x[13,4,2,j]
      
      
      
      #Resistance acquisition and reversion
      dx[12,1,2,j]=dx[12,1,2,j]-RateSusceptible*x[12,1,2,j]
      dx[15,1,2,j]=dx[15,1,2,j]+RateResistant*x[15,1,1,j]
      dx[12,2,2,j]=dx[12,2,2,j]-RateSusceptible*x[12,2,2,j]
      dx[15,2,2,j]=dx[15,2,2,j]+RateResistant*x[15,2,1,j]
      dx[12,3,2,j]=dx[12,3,2,j]-RateSusceptible*x[12,3,2,j]
      dx[15,3,2,j]=dx[15,3,2,j]+RateResistant*x[15,3,1,j]
      dx[12,4,2,j]=dx[12,4,2,j]-RateSusceptible*x[12,4,2,j]
      dx[15,4,2,j]=dx[15,4,2,j]+RateResistant*x[15,4,1,j]
      
      dx[12,1,1,j]=dx[12,1,1,j]+RateSusceptible*x[12,1,2,j]
      dx[15,1,1,j]=dx[15,1,1,j]-RateResistant*x[15,1,1,j]
      dx[12,2,1,j]=dx[12,2,1,j]+RateSusceptible*x[12,2,2,j]
      dx[15,2,1,j]=dx[15,2,1,j]-RateResistant*x[15,2,1,j]
      dx[12,3,1,j]=dx[12,3,1,j]+RateSusceptible*x[12,3,2,j]
      dx[15,3,1,j]=dx[15,3,1,j]-RateResistant*x[15,3,1,j]
      dx[12,4,1,j]=dx[12,4,1,j]+RateSusceptible*x[12,4,2,j]
      dx[15,4,1,j]=dx[15,4,1,j]-RateResistant*x[15,4,1,j]
      
    }
    dx[x+dx<0]=0
    
    #################################################################
    #dtg_impact.R: Resistance parameters and estimate of time to suppression
    d1=sum(mu[3:4,1:4,1]*apply(x[3:4,1:4,1:2,1:2],c(1,2),sum))
    d2=sum(mu[5,1:4,1]*apply(x[5,1:4,1:2,1:2],c(1),sum))
    new_treat2=c(RateTreatFirst_noDTG[1,1:4,1]*x[12,1:4,1,2],
                RateTreatFirst_noDTG[2,1:4,1]*x[12,1:4,2,2])
    #dx[3:8,1:4,1,1:2]=dx[3:8,1:4,1,1:2]+infected_m15_month[t+1]*x[3:8,1:4,1,1:2]/sum(x[3:8,1:4,1,1:2])
    ####################################################################################################################
    #res=c(dx,new_inf,death,new_diag,new_treat,death_treat,new_inf_res,death_res,death_w,death_w_inel,death_w_inel)
    
    # d1<-sum(RateTreatFirst[1,1:4,2]*apply(x[2,1:4,1:2,2],c(1),sum))+
    #   sum(RateTreatFirst_noDTG[1,1:4,2]*apply(x[12,1:4,1:2,2],c(1),sum))+
    #   sum(RateTreatFirstDTG[1,1:4,2]*apply(x[2,1:4,1:2,2],c(1),sum))
    # 
    # d2<-sum(RateTreatFirst_noDTG[1,1:4,2]*apply(x[12,1:4,1:2,2],c(1),sum))
    # 
    # d3<-sum(RateTreatFirstDTG[1,1:4,2]*apply(x[2,1:4,1:2,2],c(1),sum))
    # 
    # d4<-RateTreatFirst_noDTG[1,1:4,2]
    # d5<-RateTreatFirstDTG[1,1:4,2]
    
    new_inf=S[1]/N[1]*(sum(RateInf1[1:4,1:2,1]*apply(x[c(1),1:4,1:2,1:2],c(1,3),sum))+
                            sum(RateInf2[1:4,1:2,1]*apply(x[c(2,3,5,6,8,9,11,12,13,15),1:4,1:2,1:2],c(2,4),sum)))+
      S[2]/N[2]*(sum(RateInf1[1:4,1:2,2]*apply(x[c(1),1:4,1:2,1:2],c(1,3),sum))+
                   sum(RateInf2[1:4,1:2,2]*apply(x[c(2,3,5,6,8,9,11,12,13,15),1:4,1:2,1:2],c(2,4),sum)))
    new_inf_susc=S[1]/N[1]*(sum(RateInf1[1:4,1:2,1]*apply(x[c(1),1:4,1,1:2],c(1,2),sum))+
                              sum(RateInf2[1:4,1:2,1]*apply(x[c(2,3,5,6,8,9,11,12,13,15),1:4,1,1:2],c(2,3),sum)))+
      S[2]/N[2]*(sum(RateInf1[1:4,1:2,2]*apply(x[c(1),1:4,1,1:2],c(1,2),sum))+
                   sum(RateInf2[1:4,1:2,2]*apply(x[c(2,3,5,6,8,9,11,12,13,15),1:4,1,1:2],c(2,3),sum)))
    new_inf_res=S[1]/N[1]*(sum(RateInf1[1:4,1:2,1]*apply(x[c(1),1:4,2,1:2],c(1,2),sum))+
                             sum(RateInf2[1:4,1:2,1]*apply(x[c(2,3,5,6,8,9,11,12,13,15),1:4,2,1:2],c(2,3),sum)))+
      S[2]/N[2]*(sum(RateInf1[1:4,1:2,2]*apply(x[c(1),1:4,2,1:2],c(1,2),sum))+
                   sum(RateInf2[1:4,1:2,2]*apply(x[c(2,3,5,6,8,9,11,12,13,15),1:4,2,1:2],c(2,3),sum)))
    
    
    d1=sum(mu[3:4,1:4,1]*apply(x[13:14,1:4,1:2,1:2],c(1,2),sum))
    d2=sum(mu[5,1:4,1]*apply(x[15,1:4,1:2,1:2],c(1),sum))
    new_treat2=c(RateTreatFirst_noDTG[1,1:4,1]*x[12,1:4,1,2],
                 RateTreatFirst_noDTG[2,1:4,1]*x[12,1:4,2,2])
    ####################################################################################################################
    #res_dtg_impact=c(dx,new_inf,new_inf_susc,new_inf_res,rep(0,7))
    res_dtg_impact=c(dx,d1,d2,new_treat2)
    return(res_dtg_impact)
    #})
  })
}
mod_dtg_lsoda=function(t,x,parms){
  theta=parms[[1]]
  treat_dtg=parms[[2]]
  parms1=parms[[3]]
  p2=parms[[4]]
  x2=rep(0,250)
  
  x2[select_15(3:5,1:4,1:2,1:2)]=x[1:48]
  x2[241:242]=x[49:50]
  data=mod_dtg(t,x2,theta,treat_dtgb,params,p2)
  x2=data[select_15(3:5,1:4,1:2,1:2)]
  x2=c(x2,data[241:242])
  
  return(list(x2))
}

#18 care stages
select_18=function(r1,r2,r3,r4){
  l1=18
  l2=4
  l3=2
  l4=2
  ar=array(1:(l1*l2*l3*l4),dim=c(l1,l2,l3,l4))
  return(as.vector(ar[r1,r2,r3,r4]))
}
xstart_2005=function(par,treat_dtg){
  k1=par[9]
  k2=par[10]
  k3=par[11]
  
  start_value=array(0,dim=c(15,4,2,2))
  start_value[1,1:4,1,1:2]=c(undiag_start[1]*repart_f(k1),undiag_start[2]*repart_f(k1))
  start_value[2,1:4,1,1:2]=c(diag_start[1]*repart_f(k2),diag_start[2]*repart_f(k2))
  start_value[3,1:4,1,1:2]=c(treat_start[1]*repart_f(k3),treat_start[2]*repart_f(k3))
  
  start_value[4,1:4,1,1:2]=start_value[3,1:4,1,1:2]*0
  start_value[5,1:4,1,1:2]=start_value[3,1:4,1,1:2]*0
  start_value[3,1:4,1,1:2]=start_value[3,1:4,1,1:2]*1
  
  start_value[1:2,1:4,2,1:2]=array(rep(rep(c(0.01,0.01,0.01,0.01),each=2),2),dim=c(2,4,2))*start_value[1:2,1:4,1,1:2]
  start_value[1:2,1:4,1,1:2]=array(rep(rep(1-c(0.01,0.01,0.01,0.01),each=2),2),dim=c(2,4,2))*start_value[1:2,1:4,1,1:2]
  start_value[3,1:4,2,1:2]=0.01*start_value[3,1:4,1,1:2]
  start_value[3,1:4,1,1:2]=0.99*start_value[3,1:4,1,1:2]
  start_value[4,1:4,2,1:2]=0.01*start_value[4,1:4,1,1:2]
  start_value[4,1:4,1,1:2]=0.99*start_value[4,1:4,1,1:2]
  start_value[5,1:4,2,1:2]=0.7*start_value[5,1:4,1,1:2]
  start_value[5,1:4,1,1:2]=0.3*start_value[5,1:4,1,1:2]
  
  start_value[12:15,1:4,1:2,2]=(1-treat_dtg[1])*start_value[2:5,1:4,1:2,2]
  start_value[2:5,1:4,1:2,2]=treat_dtg[1]*start_value[2:5,1:4,1:2,2]
  
  xstart=c(start_value,rep(0,10))
  return(xstart)
}
mod_dtg_18=function(t,x,p1,treat_dtg,parms,p2,comp_time){
  
  parms2=mod_rate(t,x,p1,treat_dtg,parms,p2)
  dx=array(0,dim=c(18,4,2,2))
  x=array(x,dim=c(18,4,2,2))
  
  #Define comp_time globally, used here to fix time from which dtg-ineligible ind go to compartments 16
  
  with(as.list(parms2), {
    for(j in 1:2){
      #dx[x1,x2,x3,x4] x1:care stages x2:disease progression x3:resistance x4:gender
      #dx with only interaction with the neighbours in the compartmental model
      dx[1,1,1,j]=S[j]/N[j]*(sum(RateInf1[1:4,1:2,j]*apply(x[c(1),1:4,1,1:2],c(1,2),sum))+sum(RateInf2[1:4,1:2,j]*apply(x[c(2,3,5,6,8,9,11,12,13,15,16,18),1:4,1,1:2],c(2,3),sum)))-
        (RateStageInf[1]+RateDiag[1,j]+mu[1,1,1])*x[1,1,1,j]
      dx[2,1,1,j]=dx[2,1,1,j]+RateDiag[1,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,1,1,j]-(RateTreatFirst[1,1,j]+RateStageDiag[1]+mu[2,1,1])*x[2,1,1,j]
      dx[3,1,1,j]=dx[3,1,1,j]+RateTreatFirst[1,1,j]*x[2,1,1,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,1])*x[3,1,1,j]+RateStageTreatFirstL[1]*x[3,2,1,j]
      dx[4,1,1,j]=RateSuppFirstSusc[1]*x[3,1,1,j]-(RateFailFirstSusc[1]+mu[4,1,1])*x[4,1,1,j]+RateStageSuppFirst[1]*x[4,2,1,j]
      dx[5,1,1,j]=RateFailFirstSusc[1]*x[4,1,1,j]-(RateTreatSecond[1,j]+RateStageFailFirst[1]+mu[5,1,1])*x[5,1,1,j]
      dx[6,1,1,j]=RateTreatSecond_noDTG[1,j]*(x[11,1,1,j]+x[15,1,1,j]+x[18,1,1,j])+RateTreatSecond[1,j]*x[5,1,1,j]-(RateSuppSecond[1]+RateStageTreatSecondR[1]+mu[6,1,1])*x[6,1,1,j]+RateStageTreatSecondL[1]*x[6,2,1,j]
      dx[7,1,1,j]=RateSuppSecond[1]*x[6,1,1,j]-(RateFailSecond[1]+mu[7,1,1])*x[7,1,1,j]+RateStageSuppSecond[1]*x[7,2,1,j]
      dx[8,1,1,j]=RateFailSecond[1]*x[7,1,1,j]-(RateStageFailSecond[1]+mu[8,1,1])*x[8,1,1,j]
      
      dx[1,2,1,j]=-(RateStageInf[2]+RateDiag[2,j]+mu[1,2,1])*x[1,2,1,j]+RateStageInf[1]*x[1,1,1,j]
      dx[2,2,1,j]=dx[2,2,1,j]+RateDiag[2,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,2,1,j]-(RateTreatFirst[1,2,j]+RateStageDiag[2]+mu[2,2,1])*x[2,2,1,j]+RateStageDiag[1]*x[2,1,1,j]
      dx[3,2,1,j]=dx[3,2,1,j]+RateTreatFirst[1,2,j]*x[2,2,1,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,1])*x[3,2,1,j]+RateStageTreatFirstR[1]*x[3,1,1,j]+RateStageTreatFirstL[2]*x[3,3,1,j]
      dx[4,2,1,j]=RateSuppFirstSusc[2]*x[3,2,1,j]-(RateFailFirstSusc[2]+RateStageSuppFirst[1]+mu[4,2,1])*x[4,2,1,j]+RateStageSuppFirst[2]*x[4,3,1,j]
      dx[5,2,1,j]=RateFailFirstSusc[2]*x[4,2,1,j]-(RateTreatSecond[2,j]+RateStageFailFirst[2]+mu[5,2,1])*x[5,2,1,j]+RateStageFailFirst[1]*x[5,1,1,j]
      dx[6,2,1,j]=RateTreatSecond_noDTG[2,j]*(x[11,2,1,j]+x[15,2,1,j]+x[18,2,1,j])+RateTreatSecond[2,j]*x[5,2,1,j]-(RateSuppSecond[2]+RateStageTreatSecondR[2]+RateStageTreatSecondL[1]+mu[6,2,1])*x[6,2,1,j]+RateStageTreatSecondR[1]*x[6,1,1,j]+RateStageTreatSecondL[2]*x[6,3,1,j]
      dx[7,2,1,j]=RateSuppSecond[2]*x[6,2,1,j]-(RateFailSecond[2]+RateStageSuppSecond[1]+mu[7,2,1])*x[7,2,1,j]+RateStageSuppSecond[2]*x[7,3,1,j]
      dx[8,2,1,j]=RateFailSecond[2]*x[7,2,1,j]-(RateStageFailSecond[2]+mu[8,2,1])*x[8,2,1,j]+RateStageFailSecond[1]*x[8,1,1,j]
      
      dx[1,3,1,j]=-(RateStageInf[3]+RateDiag[3,j]+mu[1,3,1])*x[1,3,1,j]+RateStageInf[2]*x[1,2,1,j]
      dx[2,3,1,j]=dx[2,3,1,j]+RateDiag[3,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,3,1,j]-(RateTreatFirst[1,3,j]+RateStageDiag[3]+mu[2,3,1])*x[2,3,1,j]+RateStageDiag[2]*x[2,2,1,j]
      dx[3,3,1,j]=dx[3,3,1,j]+RateTreatFirst[1,3,j]*x[2,3,1,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,1])*x[3,3,1,j]+RateStageTreatFirstR[2]*x[3,2,1,j]+RateStageTreatFirstL[3]*x[3,4,1,j]
      dx[4,3,1,j]=RateSuppFirstSusc[3]*x[3,3,1,j]-(RateFailFirstSusc[3]+RateStageSuppFirst[2]+mu[4,3,1])*x[4,3,1,j]+RateStageSuppFirst[3]*x[4,4,1,j]
      dx[5,3,1,j]=RateFailFirstSusc[3]*x[4,3,1,j]-(RateTreatSecond[3,j]+RateStageFailFirst[3]+mu[5,3,1])*x[5,3,1,j]+RateStageFailFirst[2]*x[5,2,1,j]
      dx[6,3,1,j]=RateTreatSecond_noDTG[3,j]*(x[11,3,1,j]+x[15,3,1,j]+x[18,3,1,j])+RateTreatSecond[3,j]*x[5,3,1,j]-(RateSuppSecond[3]+RateStageTreatSecondR[3]+RateStageTreatSecondL[2]+mu[6,3,1])*x[6,3,1,j]+RateStageTreatSecondR[2]*x[6,2,1,j]+RateStageTreatSecondL[3]*x[6,4,1,j]
      dx[7,3,1,j]=RateSuppSecond[3]*x[6,3,1,j]-(RateFailSecond[3]+RateStageSuppSecond[2]+mu[7,3,1])*x[7,3,1,j]+RateStageSuppSecond[3]*x[7,4,1,j]
      dx[8,3,1,j]=RateFailSecond[3]*x[7,3,1,j]-(RateStageFailSecond[3]+mu[8,3,1])*x[8,3,1,j]+RateStageFailSecond[2]*x[8,2,1,j]
      
      dx[1,4,1,j]=-(RateDiag[4,j]+mu[1,4,1])*x[1,4,1,j]+RateStageInf[3]*x[1,3,1,j]
      dx[2,4,1,j]=dx[2,4,1,j]+RateDiag[4,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,4,1,j]-(RateTreatFirst[1,4,j]+mu[2,4,1])*x[2,4,1,j]+RateStageDiag[3]*x[2,3,1,j]
      dx[3,4,1,j]=dx[3,4,1,j]+RateTreatFirst[1,4,j]*x[2,4,1,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,1])*x[3,4,1,j]+RateStageTreatFirstR[3]*x[3,3,1,j]
      dx[4,4,1,j]=RateSuppFirstSusc[4]*x[3,4,1,j]-(RateFailFirstSusc[4]+RateStageSuppFirst[3]+mu[4,4,1])*x[4,4,1,j]
      dx[5,4,1,j]=RateFailFirstSusc[4]*x[4,4,1,j]-(RateTreatSecond[4,j]+mu[5,4,1])*x[5,4,1,j]+RateStageFailFirst[3]*x[5,3,1,j]
      dx[6,4,1,j]=RateTreatSecond_noDTG[4,j]*(x[11,4,1,j]+x[15,4,1,j]+x[18,4,1,j])+RateTreatSecond[4,j]*x[5,4,1,j]-(RateSuppSecond[4]+RateStageTreatSecondL[3]+mu[6,4,1])*x[6,4,1,j]+RateStageTreatSecondR[3]*x[6,3,1,j]
      dx[7,4,1,j]=RateSuppSecond[4]*x[6,4,1,j]-(RateFailSecond[4]+RateStageSuppSecond[3]+mu[7,4,1])*x[7,4,1,j]
      dx[8,4,1,j]=RateFailSecond[4]*x[7,4,1,j]-mu[8,4,1 ]*x[8,4,1,j]+RateStageFailSecond[3]*x[8,3,1,j]
      
      ##############################################################################################################################################################
      
      dx[1,1,1,j]=dx[1,1,1,j]
      dx[2,1,1,j]=dx[2,1,1,j]-RateDirectTreatSecond[1,1]*x[2,1,1,j]+RateStopTreatFirst[1]*x[3,1,1,j]+RateStopSuppFirst[1]*x[4,1,1,j]+RateStopFailFirst[1]*x[5,1,1,j]+
        (RateStopTreatSecond[1]*x[6,1,1,j]+RateStopSuppSecond[1]*x[7,1,1,j]+RateStopFailSecond[1]*x[8,1,1,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,1,1,j]=dx[3,1,1,j]-RateStopTreatFirst[1]*x[3,1,1,j]-RateTreatToFailFirstSusc[1]*x[3,1,1,j]
      dx[4,1,1,j]=dx[4,1,1,j]-RateStopSuppFirst[1]*x[4,1,1,j]+RateFailToSuppTreatFirstSusc[1]*x[5,1,1,j]-RateSuppFirstToSecond[1]*x[4,1,1,j]
      dx[5,1,1,j]=dx[5,1,1,j]-RateStopFailFirst[1]*x[5,1,1,j]-RateFailToSuppTreatFirstSusc[1]*x[5,1,1,j]+RateTreatToFailFirstSusc[1]*x[3,1,1,j]
      dx[6,1,1,j]=dx[6,1,1,j]-RateStopTreatSecond[1]*x[6,1,1,j]+RateDirectTreatSecond[1,1]*x[2,1,1,j]+RateDirectTreatSecond[1,1]*x[12,1,1,j]-RateTreatToFailSecond[1]*x[6,1,1,j]
      dx[7,1,1,j]=dx[7,1,1,j]-RateStopSuppSecond[1]*x[7,1,1,j]+RateFailToSuppTreatSecond[1]*x[8,1,1,j]+RateSuppFirstToSecond[1]*(x[4,1,1,j]+x[14,1,1,j]+x[17,1,1,j])
      dx[8,1,1,j]=dx[8,1,1,j]-RateStopFailSecond[1]*x[8,1,1,j]-RateFailToSuppTreatSecond[1]*x[8,1,1,j]+RateTreatToFailSecond[1]*x[6,1,1,j]
      
      dx[1,2,1,j]=dx[1,2,1,j]
      dx[2,2,1,j]=dx[2,2,1,j]-RateDirectTreatSecond[1,2]*x[2,2,1,j]+RateStopTreatFirst[2]*x[3,2,1,j]+RateStopSuppFirst[2]*x[4,2,1,j]+RateStopFailFirst[2]*x[5,2,1,j]+
        (RateStopTreatSecond[2]*x[6,2,1,j]+RateStopSuppSecond[2]*x[7,2,1,j]+RateStopFailSecond[2]*x[8,2,1,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,2,1,j]=dx[3,2,1,j]-RateStopTreatFirst[2]*x[3,2,1,j]-RateTreatToFailFirstSusc[2]*x[3,2,1,j]
      dx[4,2,1,j]=dx[4,2,1,j]-RateStopSuppFirst[2]*x[4,2,1,j]+RateFailToSuppTreatFirstSusc[2]*x[5,2,1,j]-RateSuppFirstToSecond[2]*x[4,2,1,j]
      dx[5,2,1,j]=dx[5,2,1,j]-RateStopFailFirst[2]*x[5,2,1,j]-RateFailToSuppTreatFirstSusc[2]*x[5,2,1,j]+RateTreatToFailFirstSusc[2]*x[3,2,1,j]
      dx[6,2,1,j]=dx[6,2,1,j]-RateStopTreatSecond[2]*x[6,2,1,j]+RateDirectTreatSecond[1,2]*x[2,2,1,j]+RateDirectTreatSecond[1,2]*x[12,2,1,j]-RateTreatToFailSecond[2]*x[6,2,1,j]
      dx[7,2,1,j]=dx[7,2,1,j]-RateStopSuppSecond[2]*x[7,2,1,j]+RateFailToSuppTreatSecond[2]*x[8,2,1,j]+RateSuppFirstToSecond[2]*(x[4,2,1,j]+x[14,2,1,j]+x[17,2,1,j])
      dx[8,2,1,j]=dx[8,2,1,j]-RateStopFailSecond[2]*x[8,2,1,j]-RateFailToSuppTreatSecond[2]*x[8,2,1,j]+RateTreatToFailSecond[2]*x[6,2,1,j]
      
      dx[1,3,1,j]=dx[1,3,1,j]
      dx[2,3,1,j]=dx[2,3,1,j]-RateDirectTreatSecond[1,3]*x[2,3,1,j]+RateStopTreatFirst[3]*x[3,3,1,j]+RateStopSuppFirst[3]*x[4,3,1,j]+RateStopFailFirst[3]*x[5,3,1,j]+
        (RateStopTreatSecond[3]*x[6,3,1,j]+RateStopSuppSecond[3]*x[7,3,1,j]+RateStopFailSecond[3]*x[8,3,1,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,3,1,j]=dx[3,3,1,j]-RateStopTreatFirst[3]*x[3,3,1,j]-RateTreatToFailFirstSusc[3]*x[3,3,1,j]
      dx[4,3,1,j]=dx[4,3,1,j]-RateStopSuppFirst[3]*x[4,3,1,j]+RateFailToSuppTreatFirstSusc[3]*x[5,3,1,j]-RateSuppFirstToSecond[3]*x[4,3,1,j]
      dx[5,3,1,j]=dx[5,3,1,j]-RateStopFailFirst[3]*x[5,3,1,j]-RateFailToSuppTreatFirstSusc[3]*x[5,3,1,j]+RateTreatToFailFirstSusc[3]*x[3,3,1,j]
      dx[6,3,1,j]=dx[6,3,1,j]-RateStopTreatSecond[3]*x[6,3,1,j]+RateDirectTreatSecond[1,3]*x[2,3,1,j]+RateDirectTreatSecond[1,3]*x[12,3,1,j]-RateTreatToFailSecond[3]*x[6,3,1,j]
      dx[7,3,1,j]=dx[7,3,1,j]-RateStopSuppSecond[3]*x[7,3,1,j]+RateFailToSuppTreatSecond[3]*x[8,3,1,j]+RateSuppFirstToSecond[3]*(x[4,3,1,j]+x[14,3,1,j]+x[17,3,1,j])
      dx[8,3,1,j]=dx[8,3,1,j]-RateStopFailSecond[3]*x[8,3,1,j]-RateFailToSuppTreatSecond[3]*x[8,3,1,j]+RateTreatToFailSecond[3]*x[6,3,1,j]
      
      dx[1,4,1,j]=dx[1,4,1,j]
      dx[2,4,1,j]=dx[2,4,1,j]-RateDirectTreatSecond[1,4]*x[2,4,1,j]+RateStopTreatFirst[4]*x[3,4,1,j]+RateStopSuppFirst[4]*x[4,4,1,j]+RateStopFailFirst[4]*x[5,4,1,j]+
        (RateStopTreatSecond[4]*x[6,4,1,j]+RateStopSuppSecond[4]*x[7,4,1,j]+RateStopFailSecond[4]*x[8,4,1,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,4,1,j]=dx[3,4,1,j]-RateStopTreatFirst[4]*x[3,4,1,j]-RateTreatToFailFirstSusc[4]*x[3,4,1,j]
      dx[4,4,1,j]=dx[4,4,1,j]-RateStopSuppFirst[4]*x[4,4,1,j]+RateFailToSuppTreatFirstSusc[4]*x[5,4,1,j]-RateSuppFirstToSecond[4]*x[4,4,1,j]
      dx[5,4,1,j]=dx[5,4,1,j]-RateStopFailFirst[4]*x[5,4,1,j]-RateFailToSuppTreatFirstSusc[4]*x[5,4,1,j]+RateTreatToFailFirstSusc[4]*x[3,4,1,j]
      dx[6,4,1,j]=dx[6,4,1,j]-RateStopTreatSecond[4]*x[6,4,1,j]+RateDirectTreatSecond[1,4]*x[2,4,1,j]+RateDirectTreatSecond[1,4]*x[12,4,1,j]-RateTreatToFailSecond[4]*x[6,4,1,j]
      dx[7,4,1,j]=dx[7,4,1,j]-RateStopSuppSecond[4]*x[7,4,1,j]+RateFailToSuppTreatSecond[4]*x[8,4,1,j]+RateSuppFirstToSecond[4]*(x[4,4,1,j]+x[14,4,1,j]+x[17,4,1,j])
      dx[8,4,1,j]=dx[8,4,1,j]-RateStopFailSecond[4]*x[8,4,1,j]-RateFailToSuppTreatSecond[4]*x[8,4,1,j]+RateTreatToFailSecond[4]*x[6,4,1,j]
      #############################################################################################################################################################################
      #############################################################################################################################################################################
      dx[1,1,2,j]=S[j]/N[j]*(sum(RateInf1[1:4,1:2,j]*apply(x[c(1),1:4,2,1:2],c(1,2),sum))+sum(RateInf2[1:4,1:2,j]*apply(x[c(2,3,5,6,8,9,11,12,13,15,16,18),1:4,2,1:2],c(2,3),sum)))-
        (RateStageInf[1]+RateDiag[1,j]+mu[1,1,2 ])*x[1,1,2,j]
      dx[2,1,2,j]=dx[2,1,2,j]+RateDiag[1,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,1,2,j]-(RateTreatFirst[2,1,j]+RateStageDiag[1]+mu[2,1,2])*x[2,1,2,j]
      dx[3,1,2,j]=dx[3,1,2,j]+RateTreatFirst[2,1,j]*x[2,1,2,j]-(RateSuppFirstResis[1]+RateStageTreatFirstR[1]+mu[3,1,2])*x[3,1,2,j]+RateStageTreatFirstL[1]*x[3,2,2,j]
      dx[4,1,2,j]=RateSuppFirstResis[1]*x[3,1,2,j]-(RateFailFirstResis[1]+mu[4,1,2])*x[4,1,2,j]+RateStageSuppFirst[1]*x[4,2,2,j]
      dx[5,1,2,j]=RateFailFirstResis[1]*x[4,1,2,j]-(RateTreatSecond[1,j]+RateStageFailFirst[1]+mu[5,1,2 ])*x[5,1,2,j]
      dx[6,1,2,j]=RateTreatSecond_noDTG[1,j]*(x[11,1,2,j]+x[15,1,2,j]+x[18,1,2,j])+RateTreatSecond[1,j]*x[5,1,2,j]-(RateSuppSecond[1]+RateStageTreatSecondR[1]+mu[6,1,2])*x[6,1,2,j]+RateStageTreatSecondL[1]*x[6,2,2,j]
      dx[7,1,2,j]=RateSuppSecond[1]*x[6,1,2,j]-(RateFailSecond[1]+mu[7,1,2])*x[7,1,2,j]+RateStageSuppSecond[1]*x[7,2,2,j]
      dx[8,1,2,j]=RateFailSecond[1]*x[7,1,2,j]-(RateStageFailSecond[1]+mu[8,1,2])*x[8,1,2,j]
      
      dx[1,2,2,j]=-(RateStageInf[2]+RateDiag[2,j]+mu[1,2,2])*x[1,2,2,j]+RateStageInf[1]*x[1,1,2,j]
      dx[2,2,2,j]=dx[2,2,2,j]+RateDiag[2,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,2,2,j]-(RateTreatFirst[2,2,j]+RateStageDiag[2]+mu[2,2,2])*x[2,2,2,j]+RateStageDiag[1]*x[2,1,2,j]
      dx[3,2,2,j]=dx[3,2,2,j]+RateTreatFirst[2,2,j]*x[2,2,2,j]-(RateSuppFirstResis[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,2])*x[3,2,2,j]+RateStageTreatFirstR[1]*x[3,1,2,j]+RateStageTreatFirstL[2]*x[3,3,2,j]
      dx[4,2,2,j]=RateSuppFirstResis[2]*x[3,2,2,j]-(RateFailFirstResis[2]+RateStageSuppFirst[1]+mu[4,2,2])*x[4,2,2,j]+RateStageSuppFirst[2]*x[4,3,2,j]
      dx[5,2,2,j]=RateFailFirstResis[2]*x[4,2,2,j]-(RateTreatSecond[2,j]+RateStageFailFirst[2]+mu[5,2,2])*x[5,2,2,j]+RateStageFailFirst[1]*x[5,1,2,j]
      dx[6,2,2,j]=RateTreatSecond_noDTG[2,j]*(x[11,2,2,j]+x[15,2,2,j]+x[18,2,2,j])+RateTreatSecond[2,j]*x[5,2,2,j]-(RateSuppSecond[2]+RateStageTreatSecondR[2]+RateStageTreatSecondL[1]+mu[6,2,2])*x[6,2,2,j]+RateStageTreatSecondR[1]*x[6,1,2,j]+RateStageTreatSecondL[2]*x[6,3,2,j]
      dx[7,2,2,j]=RateSuppSecond[2]*x[6,2,2,j]-(RateFailSecond[2]+RateStageSuppSecond[1]+mu[7,2,2])*x[7,2,2,j]+RateStageSuppSecond[2]*x[7,3,2,j]
      dx[8,2,2,j]=RateFailSecond[2]*x[7,2,2,j]-(RateStageFailSecond[2]+mu[8,2,2])*x[8,2,2,j]+RateStageFailSecond[1]*x[8,1,2,j]
      
      dx[1,3,2,j]=-(RateStageInf[3]+RateDiag[3,j]+mu[1,3,2])*x[1,3,2,j]+RateStageInf[2]*x[1,2,2,j]
      dx[2,3,2,j]=dx[2,3,2,j]+RateDiag[3,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,3,2,j]-(RateTreatFirst[2,3,j]+RateStageDiag[3]+mu[2,3,2])*x[2,3,2,j]+RateStageDiag[2]*x[2,2,2,j]
      dx[3,3,2,j]=dx[3,3,2,j]+RateTreatFirst[2,3,j]*x[2,3,2,j]-(RateSuppFirstResis[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,2])*x[3,3,2,j]+RateStageTreatFirstR[2]*x[3,2,2,j]+RateStageTreatFirstL[3]*x[3,4,2,j]
      dx[4,3,2,j]=RateSuppFirstResis[3]*x[3,3,2,j]-(RateFailFirstResis[3]+RateStageSuppFirst[2]+mu[4,3,2])*x[4,3,2,j]+RateStageSuppFirst[3]*x[4,4,2,j]
      dx[5,3,2,j]=RateFailFirstResis[3]*x[4,3,2,j]-(RateTreatSecond[3,j]+RateStageFailFirst[3]+mu[5,3,2])*x[5,3,2,j]+RateStageFailFirst[2]*x[5,2,2,j]
      dx[6,3,2,j]=RateTreatSecond_noDTG[3,j]*(x[11,3,2,j]+x[15,3,2,j]+x[18,3,2,j])+RateTreatSecond[3,j]*x[5,3,2,j]-(RateSuppSecond[3]+RateStageTreatSecondR[3]+RateStageTreatSecondL[2]+mu[6,3,2])*x[6,3,2,j]+RateStageTreatSecondR[2]*x[6,2,2,j]+RateStageTreatSecondL[3]*x[6,4,2,j]
      dx[7,3,2,j]=RateSuppSecond[3]*x[6,3,2,j]-(RateFailSecond[3]+RateStageSuppSecond[2]+mu[7,3,2])*x[7,3,2,j]+RateStageSuppSecond[3]*x[7,4,2,j]
      dx[8,3,2,j]=RateFailSecond[3]*x[7,3,2,j]-(RateStageFailSecond[3]+mu[8,3,2])*x[8,3,2,j]+RateStageFailSecond[2]*x[8,2,2,j]
      
      dx[1,4,2,j]=-(RateDiag[4,j]+mu[1,4,2])*x[1,4,2,j]+RateStageInf[3]*x[1,3,2,j]
      dx[2,4,2,j]=dx[2,4,2,j]+RateDiag[4,j]*((j==2)*(p_women_dtg)+(j==1))*x[1,4,2,j]-(RateTreatFirst[2,4,j]+mu[2,4,2])*x[2,4,2,j]+RateStageDiag[3]*x[2,3,2,j]
      dx[3,4,2,j]=dx[3,4,2,j]+RateTreatFirst[2,4,j]*x[2,4,2,j]-(RateSuppFirstResis[4]+RateStageTreatFirstL[3]+mu[3,4,2])*x[3,4,2,j]+RateStageTreatFirstR[3]*x[3,3,2,j]
      dx[4,4,2,j]=RateSuppFirstResis[4]*x[3,4,2,j]-(RateFailFirstResis[4]+RateStageSuppFirst[3]+mu[4,4,2])*x[4,4,2,j]
      dx[5,4,2,j]=RateFailFirstResis[4]*x[4,4,2,j]-(RateTreatSecond[4,j]+mu[5,4,2])*x[5,4,2,j]+RateStageFailFirst[3]*x[5,3,2,j]
      dx[6,4,2,j]=RateTreatSecond_noDTG[4,j]*(x[11,4,2,j]+x[15,4,2,j]+x[18,4,2,j])+RateTreatSecond[4,j]*x[5,4,2,j]-(RateSuppSecond[4]+RateStageTreatSecondL[3]+mu[6,4,2])*x[6,4,2,j]+RateStageTreatSecondR[3]*x[6,3,2,j]
      dx[7,4,2,j]=RateSuppSecond[4]*x[6,4,2,j]-(RateFailSecond[4]+RateStageSuppSecond[3]+mu[7,4,2])*x[7,4,2,j]
      dx[8,4,2,j]=RateFailSecond[4]*x[7,4,2,j]-mu[8,4,2]*x[8,4,2,j]+RateStageFailSecond[3]*x[8,3,2,j]
      
      ##############################################################################################################################################################
      
      dx[1,1,2,j]=dx[1,1,2,j]
      dx[2,1,2,j]=dx[2,1,2,j]-RateDirectTreatSecond[2,1]*x[2,1,2,j]+RateStopTreatFirst[1]*x[3,1,2,j]+RateStopSuppFirst[1]*x[4,1,2,j]+RateStopFailFirst[1]*x[5,1,2,j]+
        (RateStopTreatSecond[1]*x[6,1,2,j]+RateStopSuppSecond[1]*x[7,1,2,j]+RateStopFailSecond[1]*x[8,1,2,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,1,2,j]=dx[3,1,2,j]-RateStopTreatFirst[1]*x[3,1,2,j]-RateTreatToFailFirstResis[1]*x[3,1,2,j]
      dx[4,1,2,j]=dx[4,1,2,j]-RateStopSuppFirst[1]*x[4,1,2,j]+RateFailToSuppTreatFirstResis[1]*x[5,1,2,j]-RateSuppFirstToSecond[1]*x[4,1,2,j]
      dx[5,1,2,j]=dx[5,1,2,j]-RateStopFailFirst[1]*x[5,1,2,j]-RateFailToSuppTreatFirstResis[1]*x[5,1,2,j]+RateTreatToFailFirstResis[1]*x[3,1,2,j]
      dx[6,1,2,j]=dx[6,1,2,j]-RateStopTreatSecond[1]*x[6,1,2,j]+RateDirectTreatSecond[2,1]*x[2,1,2,j]+RateDirectTreatSecond[2,1]*x[12,1,2,j]-RateTreatToFailSecond[1]*x[6,1,2,j]
      dx[7,1,2,j]=dx[7,1,2,j]-RateStopSuppSecond[1]*x[7,1,2,j]+RateFailToSuppTreatSecond[1]*x[8,1,2,j]+RateSuppFirstToSecond[1]*(x[4,1,2,j]+x[14,1,2,j]+x[17,1,2,j])
      dx[8,1,2,j]=dx[8,1,2,j]-RateStopFailSecond[1]*x[8,1,2,j]-RateFailToSuppTreatSecond[1]*x[8,1,2,j]+RateTreatToFailSecond[1]*x[6,1,2,j]
      
      dx[1,2,2,j]=dx[1,2,2,j]
      dx[2,2,2,j]=dx[2,2,2,j]-RateDirectTreatSecond[2,2]*x[2,2,2,j]+RateStopTreatFirst[2]*x[3,2,2,j]+RateStopSuppFirst[2]*x[4,2,2,j]+RateStopFailFirst[2]*x[5,2,2,j]+
        (RateStopTreatSecond[2]*x[6,2,2,j]+RateStopSuppSecond[2]*x[7,2,2,j]+RateStopFailSecond[2]*x[8,2,2,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,2,2,j]=dx[3,2,2,j]-RateStopTreatFirst[2]*x[3,2,2,j]-RateTreatToFailFirstResis[2]*x[3,2,2,j]
      dx[4,2,2,j]=dx[4,2,2,j]-RateStopSuppFirst[2]*x[4,2,2,j]+RateFailToSuppTreatFirstResis[2]*x[5,2,2,j]-RateSuppFirstToSecond[2]*x[4,2,2,j]
      dx[5,2,2,j]=dx[5,2,2,j]-RateStopFailFirst[2]*x[5,2,2,j]-RateFailToSuppTreatFirstResis[2]*x[5,2,2,j]+RateTreatToFailFirstResis[2]*x[3,2,2,j]
      dx[6,2,2,j]=dx[6,2,2,j]-RateStopTreatSecond[2]*x[6,2,2,j]+RateDirectTreatSecond[2,2]*x[2,2,2,j]+RateDirectTreatSecond[2,2]*x[12,2,2,j]-RateTreatToFailSecond[2]*x[6,2,2,j]
      dx[7,2,2,j]=dx[7,2,2,j]-RateStopSuppSecond[2]*x[7,2,2,j]+RateFailToSuppTreatSecond[2]*x[8,2,2,j]+RateSuppFirstToSecond[2]*(x[4,2,2,j]+x[14,2,2,j]+x[17,2,2,j])
      dx[8,2,2,j]=dx[8,2,2,j]-RateStopFailSecond[2]*x[8,2,2,j]-RateFailToSuppTreatSecond[2]*x[8,2,2,j]+RateTreatToFailSecond[2]*x[6,2,2,j]
      
      dx[1,3,2,j]=dx[1,3,2,j]
      dx[2,3,2,j]=dx[2,3,2,j]-RateDirectTreatSecond[2,3]*x[2,3,2,j]+RateStopTreatFirst[3]*x[3,3,2,j]+RateStopSuppFirst[3]*x[4,3,2,j]+RateStopFailFirst[3]*x[5,3,2,j]+
        (RateStopTreatSecond[3]*x[6,3,2,j]+RateStopSuppSecond[3]*x[7,3,2,j]+RateStopFailSecond[3]*x[8,3,2,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,3,2,j]=dx[3,3,2,j]-RateStopTreatFirst[3]*x[3,3,2,j]-RateTreatToFailFirstResis[3]*x[3,3,2,j]
      dx[4,3,2,j]=dx[4,3,2,j]-RateStopSuppFirst[3]*x[4,3,2,j]+RateFailToSuppTreatFirstResis[3]*x[5,3,2,j]-RateSuppFirstToSecond[3]*x[4,3,2,j]
      dx[5,3,2,j]=dx[5,3,2,j]-RateStopFailFirst[3]*x[5,3,2,j]-RateFailToSuppTreatFirstResis[3]*x[5,3,2,j]+RateTreatToFailFirstResis[3]*x[3,3,2,j]
      dx[6,3,2,j]=dx[6,3,2,j]-RateStopTreatSecond[3]*x[6,3,2,j]+RateDirectTreatSecond[2,3]*x[2,3,2,j]+RateDirectTreatSecond[2,3]*x[12,3,2,j]-RateTreatToFailSecond[3]*x[6,3,2,j]
      dx[7,3,2,j]=dx[7,3,2,j]-RateStopSuppSecond[3]*x[7,3,2,j]+RateFailToSuppTreatSecond[3]*x[8,3,2,j]+RateSuppFirstToSecond[3]*(x[4,3,2,j]+x[14,3,2,j]+x[17,3,2,j])
      dx[8,3,2,j]=dx[8,3,2,j]-RateStopFailSecond[3]*x[8,3,2,j]-RateFailToSuppTreatSecond[3]*x[8,3,2,j]+RateTreatToFailSecond[3]*x[6,3,2,j]
      
      dx[1,4,2,j]=dx[1,4,2,j]
      dx[2,4,2,j]=dx[2,4,2,j]-RateDirectTreatSecond[2,4]*x[2,4,2,j]+RateStopTreatFirst[4]*x[3,4,2,j]+RateStopSuppFirst[4]*x[4,4,2,j]+RateStopFailFirst[4]*x[5,4,2,j]+
        (RateStopTreatSecond[4]*x[6,4,2,j]+RateStopSuppSecond[4]*x[7,4,2,j]+RateStopFailSecond[4]*x[8,4,2,j])*((j==2)*(p_women_dtg)+(j==1))
      dx[3,4,2,j]=dx[3,4,2,j]-RateStopTreatFirst[4]*x[3,4,2,j]-RateTreatToFailFirstResis[4]*x[3,4,2,j]
      dx[4,4,2,j]=dx[4,4,2,j]-RateStopSuppFirst[4]*x[4,4,2,j]+RateFailToSuppTreatFirstResis[4]*x[5,4,2,j]-RateSuppFirstToSecond[4]*x[4,4,2,j]
      dx[5,4,2,j]=dx[5,4,2,j]-RateStopFailFirst[4]*x[5,4,2,j]-RateFailToSuppTreatFirstResis[4]*x[5,4,2,j]+RateTreatToFailFirstResis[4]*x[3,4,2,j]
      dx[6,4,2,j]=dx[6,4,2,j]-RateStopTreatSecond[4]*x[6,4,2,j]+RateDirectTreatSecond[2,4]*x[2,4,2,j]+RateDirectTreatSecond[2,4]*x[12,4,2,j]-RateTreatToFailSecond[4]*x[6,4,2,j]
      dx[7,4,2,j]=dx[7,4,2,j]-RateStopSuppSecond[4]*x[7,4,2,j]+RateFailToSuppTreatSecond[4]*x[8,4,2,j]+RateSuppFirstToSecond[4]*(x[4,4,2,j]+x[14,4,2,j]+x[17,4,2,j])
      dx[8,4,2,j]=dx[8,4,2,j]-RateStopFailSecond[4]*x[8,4,2,j]-RateFailToSuppTreatSecond[4]*x[8,4,2,j]+RateTreatToFailSecond[4]*x[6,4,2,j]
      
      #############################################################################################################################################################################
      #############################################################################################################################################################################
      #dx with interaction between resistant and susceptible compartments
      dx[1,1,1,j]=dx[1,1,1,j]+RateSusceptible*x[1,1,2,j]
      dx[2,1,1,j]=dx[2,1,1,j]+RateSusceptible*x[2,1,2,j]
      dx[5,1,1,j]=dx[5,1,1,j]-RateResistant*x[5,1,1,j]
      
      dx[1,2,1,j]=dx[1,2,1,j]+RateSusceptible*x[1,2,2,j]
      dx[2,2,1,j]=dx[2,2,1,j]+RateSusceptible*x[2,2,2,j]
      dx[5,2,1,j]=dx[5,2,1,j]-RateResistant*x[5,2,1,j]
      
      dx[1,3,1,j]=dx[1,3,1,j]+RateSusceptible*x[1,3,2,j]
      dx[2,3,1,j]=dx[2,3,1,j]+RateSusceptible*x[2,3,2,j]
      dx[5,3,1,j]=dx[5,3,1,j]-RateResistant*x[5,3,1,j]
      
      dx[1,4,1,j]=dx[1,4,1,j]+RateSusceptible*x[1,4,2,j]
      dx[2,4,1,j]=dx[2,4,1,j]+RateSusceptible*x[2,4,2,j]
      dx[5,4,1,j]=dx[5,4,1,j]-RateResistant*x[5,4,1,j]
      
      #############################################################################################################################################################################
      dx[1,1,2,j]=dx[1,1,2,j]-RateSusceptible*x[1,1,2,j]
      dx[2,1,2,j]=dx[2,1,2,j]-RateSusceptible*x[2,1,2,j]
      dx[5,1,2,j]=dx[5,1,2,j]+RateResistant*x[5,1,1,j]
      
      dx[1,2,2,j]=dx[1,2,2,j]-RateSusceptible*x[1,2,2,j]
      dx[2,2,2,j]=dx[2,2,2,j]-RateSusceptible*x[2,2,2,j]
      dx[5,2,2,j]=dx[5,2,2,j]+RateResistant*x[5,2,1,j]
      
      dx[1,3,2,j]=dx[1,3,2,j]-RateSusceptible*x[1,3,2,j]
      dx[2,3,2,j]=dx[2,3,2,j]-RateSusceptible*x[2,3,2,j]
      dx[5,3,2,j]=dx[5,3,2,j]+RateResistant*x[5,3,1,j]
      
      dx[1,4,2,j]=dx[1,4,2,j]-RateSusceptible*x[1,4,2,j]
      dx[2,4,2,j]=dx[2,4,2,j]-RateSusceptible*x[2,4,2,j]
      dx[5,4,2,j]=dx[5,4,2,j]+RateResistant*x[5,4,1,j]
      
      
      #RateTreatFirstDTG[res,cd4,gender]------------------------------
      #NNRTI susceptible
      dx[2,1,1,j]=dx[2,1,1,j]-RateTreatFirstDTG[1,1,j]*x[2,1,1,j]
      dx[2,2,1,j]=dx[2,2,1,j]-RateTreatFirstDTG[1,2,j]*x[2,2,1,j]
      dx[2,3,1,j]=dx[2,3,1,j]-RateTreatFirstDTG[1,3,j]*x[2,3,1,j]
      dx[2,4,1,j]=dx[2,4,1,j]-RateTreatFirstDTG[1,4,j]*x[2,4,1,j]
      #NNRTI resistant
      dx[2,1,2,j]=dx[2,1,2,j]-RateTreatFirstDTG[2,1,j]*x[2,1,2,j]
      dx[2,2,2,j]=dx[2,2,2,j]-RateTreatFirstDTG[2,2,j]*x[2,2,2,j]
      dx[2,3,2,j]=dx[2,3,2,j]-RateTreatFirstDTG[2,3,j]*x[2,3,2,j]
      dx[2,4,2,j]=dx[2,4,2,j]-RateTreatFirstDTG[2,4,j]*x[2,4,2,j]
      
      #RateSwitchDTG[care,cd4,res,gender]----------------------------
      #NNRTI susceptible
      dx[3,1,1,j]=dx[3,1,1,j]-RateSwitchDTG[1,1,1,j]*x[3,1,1,j]
      dx[3,2,1,j]=dx[3,2,1,j]-RateSwitchDTG[1,2,1,j]*x[3,2,1,j]
      dx[3,3,1,j]=dx[3,3,1,j]-RateSwitchDTG[1,3,1,j]*x[3,3,1,j]
      dx[3,4,1,j]=dx[3,4,1,j]-RateSwitchDTG[1,4,1,j]*x[3,4,1,j]
      
      dx[4,1,1,j]=dx[4,1,1,j]-RateSwitchDTG[2,1,1,j]*x[4,1,1,j]
      dx[4,2,1,j]=dx[4,2,1,j]-RateSwitchDTG[2,2,1,j]*x[4,2,1,j]
      dx[4,3,1,j]=dx[4,3,1,j]-RateSwitchDTG[2,3,1,j]*x[4,3,1,j]
      dx[4,4,1,j]=dx[4,4,1,j]-RateSwitchDTG[2,4,1,j]*x[4,4,1,j]
      
      dx[5,1,1,j]=dx[5,1,1,j]-RateSwitchDTG[3,1,1,j]*x[5,1,1,j]
      dx[5,2,1,j]=dx[5,2,1,j]-RateSwitchDTG[3,2,1,j]*x[5,2,1,j]
      dx[5,3,1,j]=dx[5,3,1,j]-RateSwitchDTG[3,3,1,j]*x[5,3,1,j]
      dx[5,4,1,j]=dx[5,4,1,j]-RateSwitchDTG[3,4,1,j]*x[5,4,1,j]
      
      #NNRTI resistant
      dx[3,1,2,j]=dx[3,1,2,j]-RateSwitchDTG[1,1,2,j]*x[3,1,2,j]
      dx[3,2,2,j]=dx[3,2,2,j]-RateSwitchDTG[1,2,2,j]*x[3,2,2,j]
      dx[3,3,2,j]=dx[3,3,2,j]-RateSwitchDTG[1,3,2,j]*x[3,3,2,j]
      dx[3,4,2,j]=dx[3,4,2,j]-RateSwitchDTG[1,4,2,j]*x[3,4,2,j]
      
      dx[4,1,2,j]=dx[4,1,2,j]-RateSwitchDTG[2,1,2,j]*x[4,1,2,j]
      dx[4,2,2,j]=dx[4,2,2,j]-RateSwitchDTG[2,2,2,j]*x[4,2,2,j]
      dx[4,3,2,j]=dx[4,3,2,j]-RateSwitchDTG[2,3,2,j]*x[4,3,2,j]
      dx[4,4,2,j]=dx[4,4,2,j]-RateSwitchDTG[2,4,2,j]*x[4,4,2,j]
      
      dx[5,1,2,j]=dx[5,1,2,j]-RateSwitchDTG[3,1,2,j]*x[5,1,2,j]
      dx[5,2,2,j]=dx[5,2,2,j]-RateSwitchDTG[3,2,2,j]*x[5,2,2,j]
      dx[5,3,2,j]=dx[5,3,2,j]-RateSwitchDTG[3,3,2,j]*x[5,3,2,j]
      dx[5,4,2,j]=dx[5,4,2,j]-RateSwitchDTG[3,4,2,j]*x[5,4,2,j]
      
      #DTG compartments-------------------------------------
      #NNRTI susceptible
      dx[9,1,1,j]=dx[9,1,1,j]+RateTreatFirstDTG[1,1,j]*x[2,1,1,j]+RateSwitchDTG[3,1,1,j]*x[5,1,1,j]+RateSwitchDTG[1,1,1,j]*x[3,1,1,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,1])*x[9,1,1,j]+RateStageTreatFirstL[1]*x[9,2,1,j]
      dx[10,1,1,j]=RateSuppFirstSusc[1]*x[9,1,1,j]-(RateFailFirstDTG[1]+mu[4,1,1])*x[10,1,1,j]+RateStageSuppFirst[1]*x[10,2,1,j]+RateSwitchDTG[2,1,1,j]*x[4,1,1,j]
      dx[11,1,1,j]=RateFailFirstDTG[1]*x[10,1,1,j]-(RateTreatSecond_noDTG[1,j]+RateStageFailFirst[1]+mu[5,1,1])*x[11,1,1,j]
      dx[12,1,1,j]=dx[12,1,1,j]+RateDiag[1,j]*((j==2)*(1-p_women_dtg))*x[1,1,1,j]-(RateTreatFirst_noDTG[1,1,j]+RateStageDiag[1]+mu[2,1,1])*x[12,1,1,j]
      dx[13,1,1,j]=dx[13,1,1,j]+RateTreatFirst_noDTG[1,1,j]*(t<comp_time)*x[12,1,1,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,1])*x[13,1,1,j]+RateStageTreatFirstL[1]*x[13,2,1,j]
      dx[14,1,1,j]=dx[14,1,1,j]+RateSuppFirstSusc[1]*x[13,1,1,j]-(RateFailFirstSusc[1]+mu[4,1,1])*x[14,1,1,j]+RateStageSuppFirst[1]*x[14,2,1,j]
      dx[15,1,1,j]=dx[15,1,1,j]+RateFailFirstSusc[1]*x[14,1,1,j]-(RateTreatSecond_noDTG[1,j]+RateStageFailFirst[1]+mu[5,1,1])*x[15,1,1,j]
      dx[16,1,1,j]=dx[16,1,1,j]+RateTreatFirst_noDTG[1,1,j]*(t>=comp_time)*x[12,1,1,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,1])*x[16,1,1,j]+RateStageTreatFirstL[1]*x[16,2,1,j]
      dx[17,1,1,j]=dx[17,1,1,j]+RateSuppFirstSusc[1]*x[16,1,1,j]-(RateFailFirstSusc[1]+mu[4,1,1])*x[17,1,1,j]+RateStageSuppFirst[1]*x[17,2,1,j]
      dx[18,1,1,j]=dx[18,1,1,j]+RateFailFirstSusc[1]*x[17,1,1,j]-(RateTreatSecond_noDTG[1,j]+RateStageFailFirst[1]+mu[5,1,1])*x[18,1,1,j]
      
      dx[9,2,1,j]=dx[9,2,1,j]+RateTreatFirstDTG[1,2,j]*x[2,2,1,j]+RateSwitchDTG[3,2,1,j]*x[5,2,1,j]+RateSwitchDTG[1,2,1,j]*x[3,2,1,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,1])*x[9,2,1,j]+RateStageTreatFirstR[1]*x[9,1,1,j]+RateStageTreatFirstL[2]*x[9,3,1,j]
      dx[10,2,1,j]=RateSuppFirstSusc[2]*x[9,2,1,j]-(RateFailFirstDTG[2]+RateStageSuppFirst[1]+mu[4,2,1])*x[10,2,1,j]+RateStageSuppFirst[2]*x[10,3,1,j]+RateSwitchDTG[2,2,1,j]*x[4,2,1,j]
      dx[11,2,1,j]=RateFailFirstDTG[2]*x[10,2,1,j]-(RateTreatSecond_noDTG[2,j]+RateStageFailFirst[2]+mu[5,2,1])*x[11,2,1,j]+RateStageFailFirst[1]*x[11,1,1,j]
      dx[12,2,1,j]=dx[12,2,1,j]+RateDiag[2,j]*((j==2)*(1-p_women_dtg))*x[1,2,1,j]-(RateTreatFirst_noDTG[1,2,j]+RateStageDiag[2]+mu[2,2,1])*x[12,2,1,j]+RateStageDiag[1]*x[12,1,1,j]
      dx[13,2,1,j]=dx[13,2,1,j]+RateTreatFirst_noDTG[1,2,j]*(t<comp_time)*x[12,2,1,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,1])*x[13,2,1,j]+RateStageTreatFirstR[1]*x[13,1,1,j]+RateStageTreatFirstL[2]*x[13,3,1,j]
      dx[14,2,1,j]=dx[14,2,1,j]+RateSuppFirstSusc[2]*x[13,2,1,j]-(RateFailFirstSusc[2]+RateStageSuppFirst[1]+mu[4,2,1])*x[14,2,1,j]+RateStageSuppFirst[2]*x[14,3,1,j]
      dx[15,2,1,j]=dx[15,2,1,j]+RateFailFirstSusc[2]*x[14,2,1,j]-(RateTreatSecond_noDTG[2,j]+RateStageFailFirst[2]+mu[5,2,1])*x[15,2,1,j]+RateStageFailFirst[1]*x[15,1,1,j]
      dx[16,2,1,j]=dx[16,2,1,j]+RateTreatFirst_noDTG[1,2,j]*(t>=comp_time)*x[12,2,1,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,1])*x[16,2,1,j]+RateStageTreatFirstR[1]*x[16,1,1,j]+RateStageTreatFirstL[2]*x[16,3,1,j]
      dx[17,2,1,j]=dx[17,2,1,j]+RateSuppFirstSusc[2]*x[16,2,1,j]-(RateFailFirstSusc[2]+RateStageSuppFirst[1]+mu[4,2,1])*x[17,2,1,j]+RateStageSuppFirst[2]*x[17,3,1,j]
      dx[18,2,1,j]=dx[18,2,1,j]+RateFailFirstSusc[2]*x[17,2,1,j]-(RateTreatSecond_noDTG[2,j]+RateStageFailFirst[2]+mu[5,2,1])*x[18,2,1,j]+RateStageFailFirst[1]*x[18,1,1,j]
      
      dx[9,3,1,j]=dx[9,3,1,j]+RateTreatFirstDTG[1,3,j]*x[2,3,1,j]+RateSwitchDTG[3,3,1,j]*x[5,3,1,j]+RateSwitchDTG[1,3,1,j]*x[3,3,1,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,1])*x[9,3,1,j]+RateStageTreatFirstR[2]*x[9,2,1,j]+RateStageTreatFirstL[3]*x[9,4,1,j]
      dx[10,3,1,j]=RateSuppFirstSusc[3]*x[9,3,1,j]-(RateFailFirstDTG[3]+RateStageSuppFirst[2]+mu[4,3,1])*x[10,3,1,j]+RateStageSuppFirst[3]*x[10,4,1,j]+RateSwitchDTG[2,3,1,j]*x[4,3,1,j]
      dx[11,3,1,j]=RateFailFirstDTG[3]*x[10,3,1,j]-(RateTreatSecond_noDTG[3,j]+RateStageFailFirst[3]+mu[5,3,1])*x[11,3,1,j]+RateStageFailFirst[2]*x[11,2,1,j]
      dx[12,3,1,j]=dx[12,3,1,j]+RateDiag[3,j]*((j==2)*(1-p_women_dtg))*x[1,3,1,j]-(RateTreatFirst_noDTG[1,3,j]+RateStageDiag[3]+mu[2,3,1])*x[12,3,1,j]+RateStageDiag[2]*x[12,2,1,j]
      dx[13,3,1,j]=dx[13,3,1,j]+RateTreatFirst_noDTG[1,3,j]*(t<comp_time)*x[12,3,1,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,1])*x[13,3,1,j]+RateStageTreatFirstR[2]*x[13,2,1,j]+RateStageTreatFirstL[3]*x[13,4,1,j]
      dx[14,3,1,j]=dx[14,3,1,j]+RateSuppFirstSusc[3]*x[13,3,1,j]-(RateFailFirstSusc[3]+RateStageSuppFirst[2]+mu[4,3,1])*x[14,3,1,j]+RateStageSuppFirst[3]*x[14,4,1,j]
      dx[15,3,1,j]=dx[15,3,1,j]+RateFailFirstSusc[3]*x[14,3,1,j]-(RateTreatSecond_noDTG[3,j]+RateStageFailFirst[3]+mu[5,3,1])*x[15,3,1,j]+RateStageFailFirst[2]*x[15,2,1,j]
      dx[16,3,1,j]=dx[16,3,1,j]+RateTreatFirst_noDTG[1,3,j]*(t>=comp_time)*x[12,3,1,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,1])*x[16,3,1,j]+RateStageTreatFirstR[2]*x[16,2,1,j]+RateStageTreatFirstL[3]*x[16,4,1,j]
      dx[17,3,1,j]=dx[17,3,1,j]+RateSuppFirstSusc[3]*x[16,3,1,j]-(RateFailFirstSusc[3]+RateStageSuppFirst[2]+mu[4,3,1])*x[17,3,1,j]+RateStageSuppFirst[3]*x[17,4,1,j]
      dx[18,3,1,j]=dx[18,3,1,j]+RateFailFirstSusc[3]*x[17,3,1,j]-(RateTreatSecond_noDTG[3,j]+RateStageFailFirst[3]+mu[5,3,1])*x[18,3,1,j]+RateStageFailFirst[2]*x[18,2,1,j]
      
      dx[9,4,1,j]=dx[9,4,1,j]+RateTreatFirstDTG[1,4,j]*x[2,4,1,j]+RateSwitchDTG[3,4,1,j]*x[5,4,1,j]+RateSwitchDTG[1,4,1,j]*x[3,4,1,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,1])*x[9,4,1,j]+RateStageTreatFirstR[3]*x[9,3,1,j]
      dx[10,4,1,j]=RateSuppFirstSusc[4]*x[9,4,1,j]-(RateFailFirstDTG[4]+RateStageSuppFirst[3]+mu[4,4,1])*x[10,4,1,j]+RateSwitchDTG[2,4,1,j]*x[4,4,1,j]
      dx[11,4,1,j]=RateFailFirstDTG[4]*x[10,4,1,j]-(RateTreatSecond_noDTG[4,j]+mu[5,4,1])*x[11,4,1,j]+RateStageFailFirst[3]*x[11,3,1,j]
      dx[12,4,1,j]=dx[12,4,1,j]+RateDiag[4,j]*((j==2)*(1-p_women_dtg))*x[1,4,1,j]-(RateTreatFirst_noDTG[1,4,j]+mu[2,4,1])*x[12,4,1,j]+RateStageDiag[3]*x[12,3,1,j]
      dx[13,4,1,j]=dx[13,4,1,j]+RateTreatFirst_noDTG[1,4,j]*(t<comp_time)*x[12,4,1,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,1])*x[13,4,1,j]+RateStageTreatFirstR[3]*x[13,3,1,j]
      dx[14,4,1,j]=dx[14,4,1,j]+RateSuppFirstSusc[4]*x[13,4,1,j]-(RateFailFirstSusc[4]+RateStageSuppFirst[3]+mu[4,4,1])*x[14,4,1,j]
      dx[15,4,1,j]=dx[15,4,1,j]+RateFailFirstSusc[4]*x[14,4,1,j]-(RateTreatSecond_noDTG[4,j]+mu[5,4,1])*x[15,4,1,j]+RateStageFailFirst[3]*x[15,3,1,j]
      dx[16,4,1,j]=dx[16,4,1,j]+RateTreatFirst_noDTG[1,4,j]*(t>=comp_time)*x[12,4,1,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,1])*x[16,4,1,j]+RateStageTreatFirstR[3]*x[16,3,1,j]
      dx[17,4,1,j]=dx[17,4,1,j]+RateSuppFirstSusc[4]*x[16,4,1,j]-(RateFailFirstSusc[4]+RateStageSuppFirst[3]+mu[4,4,1])*x[17,4,1,j]
      dx[18,4,1,j]=dx[18,4,1,j]+RateFailFirstSusc[4]*x[17,4,1,j]-(RateTreatSecond_noDTG[4,j]+mu[5,4,1])*x[18,4,1,j]+RateStageFailFirst[3]*x[18,3,1,j]
      
      #NNRTI resistant
      dx[9,1,2,j]=dx[9,1,2,j]+RateTreatFirstDTG[2,1,j]*x[2,1,2,j]+RateSwitchDTG[3,1,2,j]*x[5,1,2,j]+RateSwitchDTG[1,1,2,j]*x[3,1,2,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,2])*x[9,1,2,j]+RateStageTreatFirstL[1]*x[9,2,2,j]
      dx[10,1,2,j]=RateSuppFirstSusc[1]*x[9,1,2,j]-(RateFailFirstDTG[1]+mu[4,1,2])*x[10,1,2,j]+RateStageSuppFirst[1]*x[10,2,2,j]+RateSwitchDTG[2,1,2,j]*x[4,1,2,j]
      dx[11,1,2,j]=RateFailFirstDTG[1]*x[10,1,2,j]-(RateTreatSecond_noDTG[1,j]+RateStageFailFirst[1]+mu[5,1,2])*x[11,1,2,j]
      dx[12,1,2,j]=dx[12,1,2,j]+RateDiag[1,j]*((j==2)*(1-p_women_dtg))*x[1,1,2,j]-(RateTreatFirst_noDTG[2,1,j]+RateStageDiag[1]+mu[2,1,2])*x[12,1,2,j]
      dx[13,1,2,j]=dx[13,1,2,j]+RateTreatFirst_noDTG[2,1,j]*(t<comp_time)*x[12,1,2,j]-(RateSuppFirstResis[1]+RateStageTreatFirstR[1]+mu[3,1,2])*x[13,1,2,j]+RateStageTreatFirstL[1]*x[13,2,2,j]
      dx[14,1,2,j]=dx[14,1,2,j]+RateSuppFirstResis[1]*x[13,1,2,j]-(RateFailFirstResis[1]+mu[4,1,2])*x[14,1,2,j]+RateStageSuppFirst[1]*x[14,2,2,j]
      dx[15,1,2,j]=dx[15,1,2,j]+RateFailFirstResis[1]*x[14,1,2,j]-(RateTreatSecond_noDTG[1,j]+RateStageFailFirst[1]+mu[5,1,2])*x[15,1,2,j]
      dx[16,1,2,j]=dx[16,1,2,j]+RateTreatFirst_noDTG[2,1,j]*(t>=comp_time)*x[12,1,2,j]-(RateSuppFirstResis[1]+RateStageTreatFirstR[1]+mu[3,1,2])*x[16,1,2,j]+RateStageTreatFirstL[1]*x[16,2,2,j]
      dx[17,1,2,j]=dx[17,1,2,j]+RateSuppFirstResis[1]*x[16,1,2,j]-(RateFailFirstResis[1]+mu[4,1,2])*x[17,1,2,j]+RateStageSuppFirst[1]*x[17,2,2,j]
      dx[18,1,2,j]=dx[18,1,2,j]+RateFailFirstResis[1]*x[17,1,2,j]-(RateTreatSecond_noDTG[1,j]+RateStageFailFirst[1]+mu[5,1,2])*x[18,1,2,j]
      
      
      dx[9,2,2,j]=dx[9,2,2,j]+RateTreatFirstDTG[2,2,j]*x[2,2,2,j]+RateSwitchDTG[3,2,2,j]*x[5,2,2,j]+RateSwitchDTG[1,2,2,j]*x[3,2,2,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,2])*x[9,2,2,j]+RateStageTreatFirstR[1]*x[9,1,2,j]+RateStageTreatFirstL[2]*x[9,3,2,j]
      dx[10,2,2,j]=RateSuppFirstSusc[2]*x[9,2,2,j]-(RateFailFirstDTG[2]+RateStageSuppFirst[1]+mu[4,2,2])*x[10,2,2,j]+RateStageSuppFirst[2]*x[10,3,2,j]+RateSwitchDTG[2,2,2,j]*x[4,2,2,j]
      dx[11,2,2,j]=RateFailFirstDTG[2]*x[10,2,2,j]-(RateTreatSecond_noDTG[2,j]+RateStageFailFirst[2]+mu[5,2,2])*x[11,2,2,j]+RateStageFailFirst[1]*x[11,1,2,j]
      dx[12,2,2,j]=dx[12,2,2,j]+RateDiag[2,j]*((j==2)*(1-p_women_dtg))*x[1,2,2,j]-(RateTreatFirst_noDTG[2,2,j]+RateStageDiag[2]+mu[2,2,2])*x[12,2,2,j]+RateStageDiag[1]*x[12,1,2,j]
      dx[13,2,2,j]=dx[13,2,2,j]+RateTreatFirst_noDTG[2,2,j]*(t<comp_time)*x[12,2,2,j]-(RateSuppFirstResis[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,2])*x[13,2,2,j]+RateStageTreatFirstR[1]*x[13,1,2,j]+RateStageTreatFirstL[2]*x[13,3,2,j]
      dx[14,2,2,j]=dx[14,2,2,j]+RateSuppFirstResis[2]*x[13,2,2,j]-(RateFailFirstResis[2]+RateStageSuppFirst[1]+mu[4,2,2])*x[14,2,2,j]+RateStageSuppFirst[2]*x[14,3,2,j]
      dx[15,2,2,j]=dx[15,2,2,j]+RateFailFirstResis[2]*x[14,2,2,j]-(RateTreatSecond_noDTG[2,j]+RateStageFailFirst[2]+mu[5,2,2])*x[15,2,2,j]+RateStageFailFirst[1]*x[15,1,2,j]
      dx[16,2,2,j]=dx[16,2,2,j]+RateTreatFirst_noDTG[2,2,j]*(t>=comp_time)*x[12,2,2,j]-(RateSuppFirstResis[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,2])*x[16,2,2,j]+RateStageTreatFirstR[1]*x[16,1,2,j]+RateStageTreatFirstL[2]*x[16,3,2,j]
      dx[17,2,2,j]=dx[17,2,2,j]+RateSuppFirstResis[2]*x[16,2,2,j]-(RateFailFirstResis[2]+RateStageSuppFirst[1]+mu[4,2,2])*x[17,2,2,j]+RateStageSuppFirst[2]*x[17,3,2,j]
      dx[18,2,2,j]=dx[18,2,2,j]+RateFailFirstResis[2]*x[17,2,2,j]-(RateTreatSecond_noDTG[2,j]+RateStageFailFirst[2]+mu[5,2,2])*x[18,2,2,j]+RateStageFailFirst[1]*x[18,1,2,j]
      
      
      dx[9,3,2,j]=dx[9,3,2,j]+RateTreatFirstDTG[2,3,j]*x[2,3,2,j]+RateSwitchDTG[3,3,2,j]*x[5,3,2,j]+RateSwitchDTG[1,3,2,j]*x[3,3,2,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,2])*x[9,3,2,j]+RateStageTreatFirstR[2]*x[9,2,2,j]+RateStageTreatFirstL[3]*x[9,4,2,j]
      dx[10,3,2,j]=RateSuppFirstSusc[3]*x[9,3,2,j]-(RateFailFirstDTG[3]+RateStageSuppFirst[2]+mu[4,3,2])*x[10,3,2,j]+RateStageSuppFirst[3]*x[10,4,2,j]+RateSwitchDTG[2,3,2,j]*x[4,3,2,j]
      dx[11,3,2,j]=RateFailFirstDTG[3]*x[10,3,2,j]-(RateTreatSecond_noDTG[3,j]+RateStageFailFirst[3]+mu[5,3,2])*x[11,3,2,j]+RateStageFailFirst[2]*x[11,2,2,j]
      dx[12,3,2,j]=dx[12,3,2,j]+RateDiag[3,j]*((j==2)*(1-p_women_dtg))*x[1,3,2,j]-(RateTreatFirst_noDTG[2,3,j]+RateStageDiag[3]+mu[2,3,2])*x[12,3,2,j]+RateStageDiag[2]*x[12,2,2,j]
      dx[13,3,2,j]=dx[13,3,2,j]+RateTreatFirst_noDTG[2,3,j]*(t<comp_time)*x[12,3,2,j]-(RateSuppFirstResis[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,2])*x[13,3,2,j]+RateStageTreatFirstR[2]*x[13,2,2,j]+RateStageTreatFirstL[3]*x[13,4,2,j]
      dx[14,3,2,j]=dx[14,3,2,j]+RateSuppFirstResis[3]*x[13,3,2,j]-(RateFailFirstResis[3]+RateStageSuppFirst[2]+mu[4,3,2])*x[14,3,2,j]+RateStageSuppFirst[3]*x[14,4,2,j]
      dx[15,3,2,j]=dx[15,3,2,j]+RateFailFirstResis[3]*x[14,3,2,j]-(RateTreatSecond_noDTG[3,j]+RateStageFailFirst[3]+mu[5,3,2])*x[15,3,2,j]+RateStageFailFirst[2]*x[15,2,2,j]
      dx[16,3,2,j]=dx[16,3,2,j]+RateTreatFirst_noDTG[2,3,j]*(t>=comp_time)*x[12,3,2,j]-(RateSuppFirstResis[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,2])*x[16,3,2,j]+RateStageTreatFirstR[2]*x[16,2,2,j]+RateStageTreatFirstL[3]*x[16,4,2,j]
      dx[17,3,2,j]=dx[17,3,2,j]+RateSuppFirstResis[3]*x[16,3,2,j]-(RateFailFirstResis[3]+RateStageSuppFirst[2]+mu[4,3,2])*x[17,3,2,j]+RateStageSuppFirst[3]*x[17,4,2,j]
      dx[18,3,2,j]=dx[18,3,2,j]+RateFailFirstResis[3]*x[17,3,2,j]-(RateTreatSecond_noDTG[3,j]+RateStageFailFirst[3]+mu[5,3,2])*x[18,3,2,j]+RateStageFailFirst[2]*x[18,2,2,j]
      
      
      dx[9,4,2,j]=dx[9,4,2,j]+RateTreatFirstDTG[2,4,j]*x[2,4,2,j]+RateSwitchDTG[3,4,2,j]*x[5,4,2,j]+RateSwitchDTG[1,4,2,j]*x[3,4,2,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,2])*x[9,4,2,j]+RateStageTreatFirstR[3]*x[9,3,2,j]
      dx[10,4,2,j]=RateSuppFirstSusc[4]*x[9,4,2,j]-(RateFailFirstDTG[4]+RateStageSuppFirst[3]+mu[4,4,2])*x[10,4,2,j]+RateSwitchDTG[2,4,2,j]*x[4,4,2,j]
      dx[11,4,2,j]=RateFailFirstDTG[4]*x[10,4,2,j]-(RateTreatSecond_noDTG[4,j]+mu[5,4,2])*x[11,4,2,j]+RateStageFailFirst[3]*x[11,3,2,j]
      dx[12,4,2,j]=dx[12,4,2,j]+RateDiag[4,j]*((j==2)*(1-p_women_dtg))*x[1,4,2,j]-(RateTreatFirst_noDTG[2,4,j]+mu[2,4,2])*x[12,4,2,j]+RateStageDiag[3]*x[12,3,2,j]
      dx[13,4,2,j]=dx[13,4,2,j]+RateTreatFirst_noDTG[2,4,j]*x[12,4,2,j]*(t<comp_time)-(RateSuppFirstResis[4]+RateStageTreatFirstL[3]+mu[3,4,2])*x[13,4,2,j]+RateStageTreatFirstR[3]*x[13,3,2,j]
      dx[14,4,2,j]=dx[14,4,2,j]+RateSuppFirstResis[4]*x[13,4,2,j]-(RateFailFirstResis[4]+RateStageSuppFirst[3]+mu[4,4,2])*x[14,4,2,j]
      dx[15,4,2,j]=dx[15,4,2,j]+RateFailFirstResis[4]*x[14,4,2,j]-(RateTreatSecond_noDTG[4,j]+mu[5,4,2])*x[15,4,2,j]+RateStageFailFirst[3]*x[15,3,2,j]
      dx[16,4,2,j]=dx[16,4,2,j]+RateTreatFirst_noDTG[2,4,j]*x[12,4,2,j]*(t>=comp_time)-(RateSuppFirstResis[4]+RateStageTreatFirstL[3]+mu[3,4,2])*x[16,4,2,j]+RateStageTreatFirstR[3]*x[16,3,2,j]
      dx[17,4,2,j]=dx[17,4,2,j]+RateSuppFirstResis[4]*x[16,4,2,j]-(RateFailFirstResis[4]+RateStageSuppFirst[3]+mu[4,4,2])*x[17,4,2,j]
      dx[18,4,2,j]=dx[18,4,2,j]+RateFailFirstResis[4]*x[17,4,2,j]-(RateTreatSecond_noDTG[4,j]+mu[5,4,2])*x[18,4,2,j]+RateStageFailFirst[3]*x[18,3,2,j]
      
      
      
      ##############################################################################################################################################################
      
      dx[9,1,1,j]=dx[9,1,1,j]-RateStopTreatFirst[1]*x[9,1,1,j]-RateTreatToFailFirstDTG[1]*x[9,1,1,j]
      dx[10,1,1,j]=dx[10,1,1,j]-RateStopSuppFirst[1]*x[10,1,1,j]+RateFailToSuppTreatFirstSusc[1]*x[11,1,1,j]
      dx[11,1,1,j]=dx[11,1,1,j]-RateStopFailFirst[1]*x[11,1,1,j]-RateFailToSuppTreatFirstSusc[1]*x[11,1,1,j]+RateTreatToFailFirstDTG[1]*x[9,1,1,j]
      dx[12,1,1,j]=dx[12,1,1,j]-RateDirectTreatSecond[1,1]*x[12,1,1,j]+RateStopTreatFirst[1]*x[13,1,1,j]+RateStopSuppFirst[1]*x[14,1,1,j]+RateStopFailFirst[1]*x[15,1,1,j]+
        (RateStopTreatSecond[1]*x[6,1,1,j]+RateStopSuppSecond[1]*x[7,1,1,j]+RateStopFailSecond[1]*x[8,1,1,j])*((j==2)*(1-p_women_dtg))
      dx[13,1,1,j]=dx[13,1,1,j]-RateStopTreatFirst[1]*x[13,1,1,j]-RateTreatToFailFirstSusc[1]*x[13,1,1,j]
      dx[14,1,1,j]=dx[14,1,1,j]-RateStopSuppFirst[1]*x[14,1,1,j]+RateFailToSuppTreatFirstSusc[1]*x[15,1,1,j]-RateSuppFirstToSecond[1]*x[14,1,1,j]
      dx[15,1,1,j]=dx[15,1,1,j]-RateStopFailFirst[1]*x[15,1,1,j]-RateFailToSuppTreatFirstSusc[1]*x[15,1,1,j]+RateTreatToFailFirstSusc[1]*x[13,1,1,j]
      dx[16,1,1,j]=dx[16,1,1,j]-RateStopTreatFirst[1]*x[16,1,1,j]-RateTreatToFailFirstSusc[1]*x[16,1,1,j]
      dx[17,1,1,j]=dx[17,1,1,j]-RateStopSuppFirst[1]*x[17,1,1,j]+RateFailToSuppTreatFirstSusc[1]*x[18,1,1,j]-RateSuppFirstToSecond[1]*x[17,1,1,j]
      dx[18,1,1,j]=dx[18,1,1,j]-RateStopFailFirst[1]*x[18,1,1,j]-RateFailToSuppTreatFirstSusc[1]*x[18,1,1,j]+RateTreatToFailFirstSusc[1]*x[16,1,1,j]
      
      dx[9,2,1,j]=dx[9,2,1,j]-RateStopTreatFirst[2]*x[9,2,1,j]-RateTreatToFailFirstDTG[2]*x[9,2,1,j]
      dx[10,2,1,j]=dx[10,2,1,j]-RateStopSuppFirst[2]*x[10,2,1,j]+RateFailToSuppTreatFirstSusc[2]*x[11,2,1,j]-RateSuppFirstToSecond[2]*x[10,2,1,j]
      dx[11,2,1,j]=dx[11,2,1,j]-RateStopFailFirst[2]*x[11,2,1,j]-RateFailToSuppTreatFirstSusc[2]*x[11,2,1,j]+RateTreatToFailFirstDTG[2]*x[9,2,1,j]
      dx[12,2,1,j]=dx[12,2,1,j]-RateDirectTreatSecond[1,2]*x[12,2,1,j]+RateStopTreatFirst[2]*x[13,2,1,j]+RateStopSuppFirst[2]*x[14,2,1,j]+RateStopFailFirst[2]*x[15,2,1,j]+
        (RateStopTreatSecond[2]*x[6,2,1,j]+RateStopSuppSecond[2]*x[7,2,1,j]+RateStopFailSecond[2]*x[8,2,1,j])*((j==2)*(1-p_women_dtg))
      dx[13,2,1,j]=dx[13,2,1,j]-RateStopTreatFirst[2]*x[13,2,1,j]-RateTreatToFailFirstSusc[2]*x[13,2,1,j]
      dx[14,2,1,j]=dx[14,2,1,j]-RateStopSuppFirst[2]*x[14,2,1,j]+RateFailToSuppTreatFirstSusc[2]*x[15,2,1,j]-RateSuppFirstToSecond[2]*x[14,2,1,j]
      dx[15,2,1,j]=dx[15,2,1,j]-RateStopFailFirst[2]*x[15,2,1,j]-RateFailToSuppTreatFirstSusc[2]*x[15,2,1,j]+RateTreatToFailFirstSusc[2]*x[13,2,1,j]
      dx[16,2,1,j]=dx[16,2,1,j]-RateStopTreatFirst[2]*x[16,2,1,j]-RateTreatToFailFirstSusc[2]*x[16,2,1,j]
      dx[17,2,1,j]=dx[17,2,1,j]-RateStopSuppFirst[2]*x[17,2,1,j]+RateFailToSuppTreatFirstSusc[2]*x[18,2,1,j]-RateSuppFirstToSecond[2]*x[17,2,1,j]
      dx[18,2,1,j]=dx[18,2,1,j]-RateStopFailFirst[2]*x[18,2,1,j]-RateFailToSuppTreatFirstSusc[2]*x[18,2,1,j]+RateTreatToFailFirstSusc[2]*x[16,2,1,j]
      
      dx[9,3,1,j]=dx[9,3,1,j]-RateStopTreatFirst[3]*x[9,3,1,j]-RateTreatToFailFirstDTG[3]*x[9,3,1,j]
      dx[10,3,1,j]=dx[10,3,1,j]-RateStopSuppFirst[3]*x[10,3,1,j]+RateFailToSuppTreatFirstSusc[3]*x[11,3,1,j]
      dx[11,3,1,j]=dx[11,3,1,j]-RateStopFailFirst[3]*x[11,3,1,j]-RateFailToSuppTreatFirstSusc[3]*x[11,3,1,j]+RateTreatToFailFirstDTG[3]*x[9,3,1,j]
      dx[12,3,1,j]=dx[12,3,1,j]-RateDirectTreatSecond[1,3]*x[12,3,1,j]+RateStopTreatFirst[3]*x[13,3,1,j]+RateStopSuppFirst[3]*x[14,3,1,j]+RateStopFailFirst[3]*x[15,3,1,j]+
        (RateStopTreatSecond[3]*x[6,3,1,j]+RateStopSuppSecond[3]*x[7,3,1,j]+RateStopFailSecond[3]*x[8,3,1,j])*((j==2)*(1-p_women_dtg))
      dx[13,3,1,j]=dx[13,3,1,j]-RateStopTreatFirst[3]*x[13,3,1,j]-RateTreatToFailFirstSusc[3]*x[13,3,1,j]
      dx[14,3,1,j]=dx[14,3,1,j]-RateStopSuppFirst[3]*x[14,3,1,j]+RateFailToSuppTreatFirstSusc[3]*x[15,3,1,j]-RateSuppFirstToSecond[3]*x[14,3,1,j]
      dx[15,3,1,j]=dx[15,3,1,j]-RateStopFailFirst[3]*x[15,3,1,j]-RateFailToSuppTreatFirstSusc[3]*x[15,3,1,j]+RateTreatToFailFirstSusc[3]*x[13,3,1,j]
      dx[16,3,1,j]=dx[16,3,1,j]-RateStopTreatFirst[3]*x[16,3,1,j]-RateTreatToFailFirstSusc[3]*x[16,3,1,j]
      dx[17,3,1,j]=dx[17,3,1,j]-RateStopSuppFirst[3]*x[17,3,1,j]+RateFailToSuppTreatFirstSusc[3]*x[18,3,1,j]-RateSuppFirstToSecond[3]*x[17,3,1,j]
      dx[18,3,1,j]=dx[18,3,1,j]-RateStopFailFirst[3]*x[18,3,1,j]-RateFailToSuppTreatFirstSusc[3]*x[18,3,1,j]+RateTreatToFailFirstSusc[3]*x[16,3,1,j]
      
      dx[9,4,1,j]=dx[9,4,1,j]-RateStopTreatFirst[4]*x[9,4,1,j]-RateTreatToFailFirstDTG[4]*x[9,4,1,j]
      dx[10,4,1,j]=dx[10,4,1,j]-RateStopSuppFirst[4]*x[10,4,1,j]+RateFailToSuppTreatFirstSusc[4]*x[11,4,1,j]
      dx[11,4,1,j]=dx[11,4,1,j]-RateStopFailFirst[4]*x[11,4,1,j]-RateFailToSuppTreatFirstSusc[4]*x[11,4,1,j]+RateTreatToFailFirstDTG[4]*x[9,4,1,j]
      dx[12,4,1,j]=dx[12,4,1,j]-RateDirectTreatSecond[1,4]*x[12,4,1,j]+RateStopTreatFirst[4]*x[13,4,1,j]+RateStopSuppFirst[4]*x[14,4,1,j]+RateStopFailFirst[4]*x[15,4,1,j]+
        (RateStopTreatSecond[4]*x[6,4,1,j]+RateStopSuppSecond[4]*x[7,4,1,j]+RateStopFailSecond[4]*x[8,4,1,j])*((j==2)*(1-p_women_dtg))
      dx[13,4,1,j]=dx[13,4,1,j]-RateStopTreatFirst[4]*x[13,4,1,j]-RateTreatToFailFirstSusc[4]*x[13,4,1,j]
      dx[14,4,1,j]=dx[14,4,1,j]-RateStopSuppFirst[4]*x[14,4,1,j]+RateFailToSuppTreatFirstSusc[4]*x[15,4,1,j]-RateSuppFirstToSecond[4]*x[14,4,1,j]
      dx[15,4,1,j]=dx[15,4,1,j]-RateStopFailFirst[4]*x[15,4,1,j]-RateFailToSuppTreatFirstSusc[4]*x[15,4,1,j]+RateTreatToFailFirstSusc[4]*x[13,4,1,j]
      dx[16,4,1,j]=dx[16,4,1,j]-RateStopTreatFirst[4]*x[16,4,1,j]-RateTreatToFailFirstSusc[4]*x[16,4,1,j]
      dx[17,4,1,j]=dx[17,4,1,j]-RateStopSuppFirst[4]*x[17,4,1,j]+RateFailToSuppTreatFirstSusc[4]*x[18,4,1,j]-RateSuppFirstToSecond[4]*x[17,4,1,j]
      dx[18,4,1,j]=dx[18,4,1,j]-RateStopFailFirst[4]*x[18,4,1,j]-RateFailToSuppTreatFirstSusc[4]*x[18,4,1,j]+RateTreatToFailFirstSusc[4]*x[16,4,1,j]
      
      #############################################################################################################################################################################
      #############################################################################################################################################################################
      
      dx[9,1,2,j]=dx[9,1,2,j]-RateStopTreatFirst[1]*x[9,1,2,j]-RateTreatToFailFirstDTG[1]*x[9,1,2,j]
      dx[10,1,2,j]=dx[10,1,2,j]-RateStopSuppFirst[1]*x[10,1,2,j]+RateFailToSuppTreatFirstSusc[1]*x[11,1,2,j]
      dx[11,1,2,j]=dx[11,1,2,j]-RateStopFailFirst[1]*x[11,1,2,j]-RateFailToSuppTreatFirstSusc[1]*x[11,1,2,j]+RateTreatToFailFirstDTG[1]*x[9,1,2,j]
      dx[12,1,2,j]=dx[12,1,2,j]-RateDirectTreatSecond[2,1]*x[12,1,2,j]+RateStopTreatFirst[1]*x[13,1,2,j]+RateStopSuppFirst[1]*x[14,1,2,j]+RateStopFailFirst[1]*x[15,1,2,j]+
        (RateStopTreatSecond[1]*x[6,1,2,j]+RateStopSuppSecond[1]*x[7,1,2,j]+RateStopFailSecond[1]*x[8,1,2,j])*((j==2)*(1-p_women_dtg))
      dx[13,1,2,j]=dx[13,1,2,j]-RateStopTreatFirst[1]*x[13,1,2,j]-RateTreatToFailFirstResis[1]*x[13,1,2,j]
      dx[14,1,2,j]=dx[14,1,2,j]-RateStopSuppFirst[1]*x[14,1,2,j]+RateFailToSuppTreatFirstResis[1]*x[15,1,2,j]-RateSuppFirstToSecond[1]*x[14,1,2,j]
      dx[15,1,2,j]=dx[15,1,2,j]-RateStopFailFirst[1]*x[15,1,2,j]-RateFailToSuppTreatFirstResis[1]*x[15,1,2,j]+RateTreatToFailFirstResis[1]*x[13,1,2,j]
      dx[16,1,2,j]=dx[16,1,2,j]-RateStopTreatFirst[1]*x[16,1,2,j]-RateTreatToFailFirstResis[1]*x[16,1,2,j]
      dx[17,1,2,j]=dx[17,1,2,j]-RateStopSuppFirst[1]*x[17,1,2,j]+RateFailToSuppTreatFirstResis[1]*x[18,1,2,j]-RateSuppFirstToSecond[1]*x[17,1,2,j]
      dx[18,1,2,j]=dx[18,1,2,j]-RateStopFailFirst[1]*x[18,1,2,j]-RateFailToSuppTreatFirstResis[1]*x[18,1,2,j]+RateTreatToFailFirstResis[1]*x[16,1,2,j]
      
      dx[9,2,2,j]=dx[9,2,2,j]-RateStopTreatFirst[2]*x[9,2,2,j]-RateTreatToFailFirstDTG[2]*x[9,2,2,j]
      dx[10,2,2,j]=dx[10,2,2,j]-RateStopSuppFirst[2]*x[10,2,2,j]+RateFailToSuppTreatFirstSusc[2]*x[11,2,2,j]
      dx[11,2,2,j]=dx[11,2,2,j]-RateStopFailFirst[2]*x[11,2,2,j]-RateFailToSuppTreatFirstSusc[2]*x[11,2,2,j]+RateTreatToFailFirstDTG[2]*x[9,2,2,j]
      dx[12,2,2,j]=dx[12,2,2,j]-RateDirectTreatSecond[2,2]*x[12,2,2,j]+RateStopTreatFirst[2]*x[13,2,2,j]+RateStopSuppFirst[2]*x[14,2,2,j]+RateStopFailFirst[2]*x[15,2,2,j]+
        (RateStopTreatSecond[2]*x[6,2,2,j]+RateStopSuppSecond[2]*x[7,2,2,j]+RateStopFailSecond[2]*x[8,2,2,j])*((j==2)*(1-p_women_dtg))
      dx[13,2,2,j]=dx[13,2,2,j]-RateStopTreatFirst[2]*x[13,2,2,j]-RateTreatToFailFirstResis[2]*x[13,2,2,j]
      dx[14,2,2,j]=dx[14,2,2,j]-RateStopSuppFirst[2]*x[14,2,2,j]+RateFailToSuppTreatFirstResis[2]*x[15,2,2,j]-RateSuppFirstToSecond[2]*x[14,2,2,j]
      dx[15,2,2,j]=dx[15,2,2,j]-RateStopFailFirst[2]*x[15,2,2,j]-RateFailToSuppTreatFirstResis[2]*x[15,2,2,j]+RateTreatToFailFirstResis[2]*x[13,2,2,j]
      dx[16,2,2,j]=dx[16,2,2,j]-RateStopTreatFirst[2]*x[16,2,2,j]-RateTreatToFailFirstResis[2]*x[16,2,2,j]
      dx[17,2,2,j]=dx[17,2,2,j]-RateStopSuppFirst[2]*x[17,2,2,j]+RateFailToSuppTreatFirstResis[2]*x[18,2,2,j]-RateSuppFirstToSecond[2]*x[17,2,2,j]
      dx[18,2,2,j]=dx[18,2,2,j]-RateStopFailFirst[2]*x[18,2,2,j]-RateFailToSuppTreatFirstResis[2]*x[18,2,2,j]+RateTreatToFailFirstResis[2]*x[16,2,2,j]
      
      dx[9,3,2,j]=dx[9,3,2,j]-RateStopTreatFirst[3]*x[9,3,2,j]-RateTreatToFailFirstDTG[3]*x[9,3,2,j]
      dx[10,3,2,j]=dx[10,3,2,j]-RateStopSuppFirst[3]*x[10,3,2,j]+RateFailToSuppTreatFirstSusc[3]*x[11,3,2,j]
      dx[11,3,2,j]=dx[11,3,2,j]-RateStopFailFirst[3]*x[11,3,2,j]-RateFailToSuppTreatFirstSusc[3]*x[11,3,2,j]+RateTreatToFailFirstDTG[3]*x[9,3,2,j]
      dx[12,3,2,j]=dx[12,3,2,j]-RateDirectTreatSecond[2,3]*x[12,3,2,j]+RateStopTreatFirst[3]*x[13,3,2,j]+RateStopSuppFirst[3]*x[14,3,2,j]+RateStopFailFirst[3]*x[15,3,2,j]+
        (RateStopTreatSecond[3]*x[6,3,2,j]+RateStopSuppSecond[3]*x[7,3,2,j]+RateStopFailSecond[3]*x[8,3,2,j])*((j==2)*(1-p_women_dtg))
      dx[13,3,2,j]=dx[13,3,2,j]-RateStopTreatFirst[3]*x[13,3,2,j]-RateTreatToFailFirstResis[3]*x[13,3,2,j]
      dx[14,3,2,j]=dx[14,3,2,j]-RateStopSuppFirst[3]*x[14,3,2,j]+RateFailToSuppTreatFirstResis[3]*x[15,3,2,j]-RateSuppFirstToSecond[3]*x[14,3,2,j]
      dx[15,3,2,j]=dx[15,3,2,j]-RateStopFailFirst[3]*x[15,3,2,j]-RateFailToSuppTreatFirstResis[3]*x[15,3,2,j]+RateTreatToFailFirstResis[3]*x[13,3,2,j]
      dx[16,3,2,j]=dx[16,3,2,j]-RateStopTreatFirst[3]*x[16,3,2,j]-RateTreatToFailFirstResis[3]*x[16,3,2,j]
      dx[17,3,2,j]=dx[17,3,2,j]-RateStopSuppFirst[3]*x[17,3,2,j]+RateFailToSuppTreatFirstResis[3]*x[18,3,2,j]-RateSuppFirstToSecond[3]*x[17,3,2,j]
      dx[18,3,2,j]=dx[18,3,2,j]-RateStopFailFirst[3]*x[18,3,2,j]-RateFailToSuppTreatFirstResis[3]*x[18,3,2,j]+RateTreatToFailFirstResis[3]*x[16,3,2,j]
      
      dx[9,4,2,j]=dx[9,4,2,j]-RateStopTreatFirst[4]*x[9,4,2,j]-RateTreatToFailFirstDTG[4]*x[9,4,2,j]
      dx[10,4,2,j]=dx[10,4,2,j]-RateStopSuppFirst[4]*x[10,4,2,j]+RateFailToSuppTreatFirstSusc[4]*x[11,4,2,j]
      dx[11,4,2,j]=dx[11,4,2,j]-RateStopFailFirst[4]*x[11,4,2,j]-RateFailToSuppTreatFirstSusc[4]*x[11,4,2,j]+RateTreatToFailFirstDTG[4]*x[9,4,2,j]
      dx[12,4,2,j]=dx[12,4,2,j]-RateDirectTreatSecond[2,4]*x[12,4,2,j]+RateStopTreatFirst[4]*x[13,4,2,j]+RateStopSuppFirst[4]*x[14,4,2,j]+RateStopFailFirst[4]*x[15,4,2,j]+
        (RateStopTreatSecond[4]*x[6,4,2,j]+RateStopSuppSecond[4]*x[7,4,2,j]+RateStopFailSecond[4]*x[8,4,2,j])*((j==2)*(1-p_women_dtg))
      dx[13,4,2,j]=dx[13,4,2,j]-RateStopTreatFirst[4]*x[13,4,2,j]-RateTreatToFailFirstResis[4]*x[13,4,2,j]
      dx[14,4,2,j]=dx[14,4,2,j]-RateStopSuppFirst[4]*x[14,4,2,j]+RateFailToSuppTreatFirstResis[4]*x[15,4,2,j]-RateSuppFirstToSecond[4]*x[14,4,2,j]
      dx[15,4,2,j]=dx[15,4,2,j]-RateStopFailFirst[4]*x[15,4,2,j]-RateFailToSuppTreatFirstResis[4]*x[15,4,2,j]+RateTreatToFailFirstResis[4]*x[13,4,2,j]
      dx[16,4,2,j]=dx[16,4,2,j]-RateStopTreatFirst[4]*x[16,4,2,j]-RateTreatToFailFirstResis[4]*x[16,4,2,j]
      dx[17,4,2,j]=dx[17,4,2,j]-RateStopSuppFirst[4]*x[17,4,2,j]+RateFailToSuppTreatFirstResis[4]*x[18,4,2,j]-RateSuppFirstToSecond[4]*x[17,4,2,j]
      dx[18,4,2,j]=dx[18,4,2,j]-RateStopFailFirst[4]*x[18,4,2,j]-RateFailToSuppTreatFirstResis[4]*x[18,4,2,j]+RateTreatToFailFirstResis[4]*x[16,4,2,j]
      
      
      #Resistance acquisition and reversion
      dx[12,1,2,j]=dx[12,1,2,j]-RateSusceptible*x[12,1,2,j]
      dx[15,1,2,j]=dx[15,1,2,j]+RateResistant*x[15,1,1,j]
      dx[18,1,2,j]=dx[18,1,2,j]+RateResistant*x[18,1,1,j]
      dx[12,2,2,j]=dx[12,2,2,j]-RateSusceptible*x[12,2,2,j]
      dx[15,2,2,j]=dx[15,2,2,j]+RateResistant*x[15,2,1,j]
      dx[18,2,2,j]=dx[18,2,2,j]+RateResistant*x[18,2,1,j]
      dx[12,3,2,j]=dx[12,3,2,j]-RateSusceptible*x[12,3,2,j]
      dx[15,3,2,j]=dx[15,3,2,j]+RateResistant*x[15,3,1,j]
      dx[18,3,2,j]=dx[18,3,2,j]+RateResistant*x[18,3,1,j]
      dx[12,4,2,j]=dx[12,4,2,j]-RateSusceptible*x[12,4,2,j]
      dx[15,4,2,j]=dx[15,4,2,j]+RateResistant*x[15,4,1,j]
      dx[18,4,2,j]=dx[18,4,2,j]+RateResistant*x[18,4,1,j]
      
      dx[12,1,1,j]=dx[12,1,1,j]+RateSusceptible*x[12,1,2,j]
      dx[15,1,1,j]=dx[15,1,1,j]-RateResistant*x[15,1,1,j]
      dx[18,1,1,j]=dx[18,1,1,j]-RateResistant*x[18,1,1,j]
      dx[12,2,1,j]=dx[12,2,1,j]+RateSusceptible*x[12,2,2,j]
      dx[15,2,1,j]=dx[15,2,1,j]-RateResistant*x[15,2,1,j]
      dx[18,2,1,j]=dx[18,2,1,j]-RateResistant*x[18,2,1,j]
      dx[12,3,1,j]=dx[12,3,1,j]+RateSusceptible*x[12,3,2,j]
      dx[15,3,1,j]=dx[15,3,1,j]-RateResistant*x[15,3,1,j]
      dx[18,3,1,j]=dx[18,3,1,j]-RateResistant*x[18,3,1,j]
      dx[12,4,1,j]=dx[12,4,1,j]+RateSusceptible*x[12,4,2,j]
      dx[15,4,1,j]=dx[15,4,1,j]-RateResistant*x[15,4,1,j]
      dx[18,4,1,j]=dx[18,4,1,j]-RateResistant*x[18,4,1,j]
      
    }
    dx[x+dx<0]=0
    #dx[12:15,1:4,1:2,1]=0
    #new_infections
    new_inf=S[1]/N[1]*(sum(RateInf1[1:4,1:2,1]*apply(x[c(1),1:4,1:2,1:2],c(1,3),sum))+
                         sum(RateInf2[1:4,1:2,1]*apply(x[c(2,3,5,6,8,9,11,12,13,15,16,18),1:4,1:2,1:2],c(2,4),sum)))+
      S[2]/N[2]*(sum(RateInf1[1:4,1:2,2]*apply(x[c(1),1:4,1:2,1:2],c(1,3),sum))+
                   sum(RateInf2[1:4,1:2,2]*apply(x[c(2,3,5,6,8,9,11,12,13,15,16,18),1:4,1:2,1:2],c(2,4),sum)))
    #death
    death=sum(apply(x[1:8,1:4,1:2,1:2],c(1,2,3),sum)*mu)+sum(apply(x[12,1:4,1:2,1:2],c(1,2),sum)*mu[2,1:4,1:2])+
      sum(apply(x[9:11,1:4,1:2,1:2],c(1,2,3),sum)*mu[3:5,1:4,1:2])+
      sum(apply(x[13:15,1:4,1:2,1:2],c(1,2,3),sum)*mu[3:5,1:4,1:2])+
      sum(apply(x[16:18,1:4,1:2,1:2],c(1,2,3),sum)*mu[3:5,1:4,1:2])
    
    
    death_w=sum(x[1:8,1:4,1:2,2]*mu)+sum(x[12,1:4,1:2,2]*mu[2,1:4,1:2])+
      sum(x[9:11,1:4,1:2,2]*mu[3:5,1:4,1:2])+sum(x[13:15,1:4,1:2,2]*mu[3:5,1:4,1:2])
    death_w_nnrti=sum(x[3:5,1:4,1:2,2]*mu[3:5,1:4,1:2])+sum(x[13:15,1:4,1:2,2]*mu[3:5,1:4,1:2])
    death_w_inel=sum(x[13:15,1:4,1:2,2]*mu[3:5,1:4,1:2])
    death_w_inel_2018=ifelse(sum(x[16:18,1:4,1:2,2])==0,
                             0,sum(x[16:18,1:4,1:2,2]*mu[3:5,1:4,1:2])/sum(x[16:18,1:4,1:2,2]))
    
    treat_16=t(RateTreatFirst_noDTG[1:2,1:4,2])*x[12,1:4,1:2,2]
    
    
    death2=sum(apply(x[3:5,1:4,1:2,1:2],c(1,2,3),sum)*mu[3:5,1:4,1:2])
    death3=sum(apply(x[13:15,1:4,1:2,1:2],c(1,2,3),sum)*mu[3:5,1:4,1:2])
    death4=sum(apply(x[9:11,1:4,1:2,1:2],c(1,2,3),sum)*mu[3:5,1:4,1:2])
    death5=sum(apply(x[2,1:4,1:2,1:2],c(1,2),sum)*mu[2,1:4,1:2])
    #new_infections resistant
    new_inf_res=S[1]/N[1]*(sum(RateInf1[1:4,1:2,1]*apply(x[c(1),1:4,2,1:2],c(1,2),sum))+
                             sum(RateInf2[1:4,1:2,1]*apply(x[c(2,3,5,6,8,9,11,12,13,15),1:4,2,1:2],c(2,3),sum)))+
      S[2]/N[2]*(sum(RateInf1[1:4,1:2,2]*apply(x[c(1),1:4,2,1:2],c(1,2),sum))+
                   sum(RateInf2[1:4,1:2,2]*apply(x[c(2,3,5,6,8,9,11,12,13,15),1:4,2,1:2],c(2,3),sum)))
    
    n1=S[1]/N[1]*sum(RateInf1[1:4,1:2,1]*apply(x[c(1),1:4,2,1:2],c(1,2),sum))#susc inf
    n2=S[1]/N[1]*sum(RateInf2[1:4,1:2,1]*apply(x[c(2,3,5,6,8,9,11),1:4,2,1:2],c(2,3),sum))#susc diag
    n3=S[2]/N[2]*sum(RateInf1[1:4,1:2,2]*apply(x[c(1),1:4,2,1:2],c(1,2),sum))#res inf
    n4=S[2]/N[2]*sum(RateInf2[1:4,1:2,2]*apply(x[c(2,3,5,6,8,9,11),1:4,2,1:2],c(2,3),sum))#res diag
    #death resistant
    death_res=sum(apply(x[1:8,1:4,2,1:2],c(1,2),sum)*mu[1:8,1:4,2])+sum(apply(x[12,1:4,2,1:2],c(1),sum)*mu[2,1:4,2])+
      sum(apply(x[9:11,1:4,2,1:2],c(1,2),sum)*mu[3:5,1:4,2])+sum(apply(x[13:15,1:4,2,1:2],c(1,2),sum)*mu[3:5,1:4,2])
    
    
    death_treat=sum(apply(x[3:8,1:4,1:2,1:2],c(1,2,3),sum)*mu[3:8,1:4,1:2])+
      sum(apply(x[9:11,1:4,1:2,1:2],c(1,2,3),sum)*mu[3:5,1:4,1:2])+sum(apply(x[13:15,1:4,1:2,1:2],c(1,2,3),sum)*mu[3:5,1:4,1:2])
    death_treat_res=sum(apply(x[3:8,1:4,2,1:2],c(1,2),sum)*mu[3:8,1:4,2])+
      sum(apply(x[9:11,1:4,2,1:2],c(1,2),sum)*mu[3:5,1:4,2])+sum(apply(x[13:15,1:4,2,1:2],c(1,2),sum)*mu[3:5,1:4,2])
    
    
    new_treat=sum(RateTreatFirst[1:2,1:4,1:2]*aperm(x[2,1:4,1:2,1:2],perm=c(2,1,3)))+
      sum(RateTreatFirstDTG[1:2,1:4,1:2]*aperm(x[2,1:4,1:2,1:2],perm=c(2,1,3)))+
      sum(RateDirectTreatSecond[1:2,1:4]*aperm(apply(x[2,1:4,1:2,1:2],c(1,2),sum),c(2,1)))+
      sum(RateTreatFirst[1:2,1:4,1:2]*aperm(x[12,1:4,1:2,1:2],perm=c(2,1,3)))+
      sum(RateDirectTreatSecond[1:2,1:4]*aperm(apply(x[12,1:4,1:2,1:2],c(1,2),sum),c(2,1)))
    
    
    new_treat_res=sum(RateTreatFirst[2,1:4,1:2]*x[2,1:4,2,1:2])+
      sum(RateTreatFirstDTG[2,1:4,1:2]*x[2,1:4,2,1:2])+
      sum(RateDirectTreatSecond[2,1:4]*apply(x[2,1:4,2,1:2],c(1),sum))+
      sum(RateTreatFirst[2,1:4,1:2]*x[12,1:4,2,1:2])+
      sum(RateDirectTreatSecond[2,1:4]*apply(x[12,1:4,2,1:2],c(1),sum))
    
    death_diag=sum(apply(x[12,1:4,1:2,1:2],c(1,2),sum)*mu[2,1:4,1:2])+sum(apply(x[2,1:4,1:2,1:2],c(1,2),sum)*mu[2,1:4,1:2])
    death_diag_res=sum(apply(x[12,1:4,2,1:2],c(1),sum)*mu[2,1:4,2])+sum(apply(x[2,1:4,2,1:2],c(1),sum)*mu[2,1:4,2])
    
    new_diag=sum(RateDiag[1:4,1:2]*apply(x[1,1:4,1:2,1:2],c(1,3),sum))
    new_diag_res=sum(RateDiag[1:4,1:2]*apply(x[1,1:4,2,1:2],c(1,2),sum))
    
    
    #why last cd4 class go down quicker in res than in susc
    new_pi_susc=sum(apply(RateDirectTreatSecond[1:2,1:4],c(2),sum)*apply(x[c(2,12),1:4,1,1:2],c(2),sum))+
      sum(RateTreatSecond[1:4,1:2]*x[c(5),1:4,1,1:2])+
      sum(RateTreatSecond_noDTG[1:4,1:2]*apply(x[c(5,15),1:4,1,1:2],c(2,3),sum))
    new_pi_res=sum(apply(RateDirectTreatSecond[1:2,1:4],c(2),sum)*apply(x[c(2,12),1:4,2,1:2],c(2),sum))+
      sum(RateTreatSecond[1:4,1:2]*x[c(5),1:4,2,1:2])+
      sum(RateTreatSecond_noDTG[1:4,1:2]*apply(x[c(5,15),1:4,2,1:2],c(2,3),sum))
    new_pi_susc_4=sum(RateDirectTreatSecond[1:2,4]*apply(x[c(2,12),4,1,1:2],c(2),sum))+
      sum(RateTreatSecond[4,1:2]*x[c(5),4,1,1:2])+
      sum(RateTreatSecond_noDTG[4,1:2]*apply(x[c(5,15),4,1,1:2],c(2),sum))
    new_pi_res_4=sum(RateDirectTreatSecond[1:2,4]*apply(x[c(2,12),4,2,1:2],c(2),sum))+
      sum(RateTreatSecond[4,1:2]*x[c(5),4,2,1:2])+
      sum(RateTreatSecond_noDTG[4,1:2]*apply(x[c(5,15),4,2,1:2],c(2),sum))
    
    
    m_susc=sum(RateDiag[1:4,1]*x[1,1:4,1,1]+RateDiag[1:4,2]*p_women_dtg*x[1,1:4,1,2])
    m_res=sum(RateDiag[1:4,1]*x[1,1:4,2,1]+RateDiag[1:4,2]*p_women_dtg*x[1,1:4,2,2])
    
    m_susc=m_susc-sum(RateTreatFirst[1,1:4,1]*x[2,1:4,1,1]+RateTreatFirst[1,1:4,2]*x[2,1:4,1,2])
    m_res=m_res-sum(RateTreatFirst[2,1:4,1]*x[2,1:4,2,1]+RateTreatFirst[2,1:4,2]*x[2,1:4,2,2])
    
    m_susc=m_susc-sum(RateDirectTreatSecond[1,1:4]*x[2,1:4,1,1]+RateDirectTreatSecond[1,1:4]*x[2,1:4,1,2])
    m_res=m_res-sum(RateDirectTreatSecond[2,1:4]*x[2,1:4,2,1]+RateDirectTreatSecond[2,1:4]*x[2,1:4,2,2])
    
    m_susc=m_susc-sum(RateTreatFirstDTG[1,1:4,1]*x[2,1:4,1,1]+RateTreatFirstDTG[1,1:4,2]*x[2,1:4,1,2])
    m_res=m_res-sum(RateTreatFirstDTG[2,1:4,1]*x[2,1:4,2,1]+RateTreatFirstDTG[2,1:4,2]*x[2,1:4,2,2])
    
    m_susc=m_susc-sum(mu[2,1:4,1]*x[2,1:4,1,1]+mu[2,1:4,1]*x[2,1:4,1,2])
    m_res=m_res-sum(mu[2,1:4,2]*x[2,1:4,2,1]+mu[2,1:4,2]*x[2,1:4,2,2])
    
    # new_treat_susc=sum(RateTreatFirstDTG[1,1:4,1:2]*x[2,1:4,1,1:2])
    # new_treat_res=sum(RateTreatFirstDTG[2,1:4,1:2]*x[2,1:4,2,1:2])
    
    
    mort_susc=sum(mu[1,1:4,1]*apply(x[1,1:4,1,1:2],1,sum))
    mort_res=sum(mu[1,1:4,2]*apply(x[1,1:4,2,1:2],1,sum))
    
    # #new_diag_susc=sum(x[2,1:4,1,1:2])
    # #new_diag_res=sum(x[2,1:4,2,1:2])
    # mort_susc=sum(mu[2,1:4,1]*apply(x[2,1:4,1,1:2],1,sum))
    # mort_res=sum(mu[2,1:4,2]*apply(x[2,1:4,2,1:2],1,sum))
    
    #dx[3:8,1:4,1,1:2]=dx[3:8,1:4,1,1:2]+infected_m15_month[t+1]*x[3:8,1:4,1,1:2]/sum(x[3:8,1:4,1,1:2])
    ####################################################################################################################
    #res=c(dx,treat_16,0,death_w_inel_2018)
    #res=c(dx,new_inf,death,new_diag,new_treat,death_treat,new_inf_res,death_res,death_w,death_w_nnrti,death_w_inel_2018)
    
    res=c(dx,new_inf,death,new_diag,new_treat,death_treat,new_inf_res,death_res,new_diag_res,new_treat_res,death_treat_res)
    #res=c(dx,new_inf,death,new_diag,new_pi_susc,new_pi_res,new_inf_res,death_res,new_diag_res,new_pi_susc_4,new_pi_res_4)
    #res=c(dx,new_treat_susc,new_treat_res,new_diag_susc,new_diag_res,mort_susc,mort_res)
    #res=c(dx,new_inf,death,new_diag,new_treat,death_diag,new_inf_res,death_res,new_diag_res,new_treat_res,death_diag_res)
    return(c(res))
    #return(list(c(res,newinf,newsu,newres,death)))
    #})
  })
}

#Mortality computation
#15
# mod_dtg=function(t,x,p1,treat_dtg,parms){
#   
#   parms2=mod_rate(t,x,p1,treat_dtg,parms)
#   dx=array(0,dim=c(15,4,2,2))
#   x=array(x,dim=c(15,4,2,2))
#   
#   #Define comp_time globally, used here to fix time from which dtg-ineligible ind go to compartments 16
#   
#   with(as.list(parms2), {
#     for(j in 1:2){
#       dx[13,1,1,j]=dx[13,1,1,j]+RateTreatFirst_noDTG[1,1,j]*x[12,1,1,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,1])*x[13,1,1,j]+RateStageTreatFirstL[1]*x[13,2,1,j]
#       dx[14,1,1,j]=dx[14,1,1,j]+RateSuppFirstSusc[1]*x[13,1,1,j]-(RateFailFirstSusc[1]+mu[4,1,1])*x[14,1,1,j]+RateStageSuppFirst[1]*x[14,2,1,j]
#       dx[15,1,1,j]=dx[15,1,1,j]+RateFailFirstSusc[1]*x[14,1,1,j]-(RateTreatSecond_noDTG[1,j]+RateStageFailFirst[1]+mu[5,1,1])*x[15,1,1,j]
#       
#       dx[13,2,1,j]=dx[13,2,1,j]+RateTreatFirst_noDTG[1,2,j]*x[12,2,1,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,1])*x[13,2,1,j]+RateStageTreatFirstR[1]*x[13,1,1,j]+RateStageTreatFirstL[2]*x[13,3,1,j]
#       dx[14,2,1,j]=dx[14,2,1,j]+RateSuppFirstSusc[2]*x[13,2,1,j]-(RateFailFirstSusc[2]+RateStageSuppFirst[1]+mu[4,2,1])*x[14,2,1,j]+RateStageSuppFirst[2]*x[14,3,1,j]
#       dx[15,2,1,j]=dx[15,2,1,j]+RateFailFirstSusc[2]*x[14,2,1,j]-(RateTreatSecond_noDTG[2,j]+RateStageFailFirst[2]+mu[5,2,1])*x[15,2,1,j]+RateStageFailFirst[1]*x[15,1,1,j]
#       
#       dx[13,3,1,j]=dx[13,3,1,j]+RateTreatFirst_noDTG[1,3,j]*x[12,3,1,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,1])*x[13,3,1,j]+RateStageTreatFirstR[2]*x[13,2,1,j]+RateStageTreatFirstL[3]*x[13,4,1,j]
#       dx[14,3,1,j]=dx[14,3,1,j]+RateSuppFirstSusc[3]*x[13,3,1,j]-(RateFailFirstSusc[3]+RateStageSuppFirst[2]+mu[4,3,1])*x[14,3,1,j]+RateStageSuppFirst[3]*x[14,4,1,j]
#       dx[15,3,1,j]=dx[15,3,1,j]+RateFailFirstSusc[3]*x[14,3,1,j]-(RateTreatSecond_noDTG[3,j]+RateStageFailFirst[3]+mu[5,3,1])*x[15,3,1,j]+RateStageFailFirst[2]*x[15,2,1,j]
#       
#       dx[13,4,1,j]=dx[13,4,1,j]+RateTreatFirst_noDTG[1,4,j]*x[12,4,1,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,1])*x[13,4,1,j]+RateStageTreatFirstR[3]*x[13,3,1,j]
#       dx[14,4,1,j]=dx[14,4,1,j]+RateSuppFirstSusc[4]*x[13,4,1,j]-(RateFailFirstSusc[4]+RateStageSuppFirst[3]+mu[4,4,1])*x[14,4,1,j]
#       dx[15,4,1,j]=dx[15,4,1,j]+RateFailFirstSusc[4]*x[14,4,1,j]-(RateTreatSecond_noDTG[4,j]+mu[5,4,1])*x[15,4,1,j]+RateStageFailFirst[3]*x[15,3,1,j]
#       
#       #NNRTI resistant
#       dx[13,1,2,j]=dx[13,1,2,j]+RateTreatFirst_noDTG[2,1,j]*x[12,1,2,j]-(RateSuppFirstResis[1]+RateStageTreatFirstR[1]+mu[3,1,2])*x[13,1,2,j]+RateStageTreatFirstL[1]*x[13,2,2,j]
#       dx[14,1,2,j]=dx[14,1,2,j]+RateSuppFirstResis[1]*x[13,1,2,j]-(RateFailFirstResis[1]+mu[4,1,2])*x[14,1,2,j]+RateStageSuppFirst[1]*x[14,2,2,j]
#       dx[15,1,2,j]=dx[15,1,2,j]+RateFailFirstResis[1]*x[14,1,2,j]-(RateTreatSecond_noDTG[1,j]+RateStageFailFirst[1]+mu[5,1,2])*x[15,1,2,j]
#       
#       dx[13,2,2,j]=dx[13,2,2,j]+RateTreatFirst_noDTG[2,2,j]*x[12,2,2,j]-(RateSuppFirstResis[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,2])*x[13,2,2,j]+RateStageTreatFirstR[1]*x[13,1,2,j]+RateStageTreatFirstL[2]*x[13,3,2,j]
#       dx[14,2,2,j]=dx[14,2,2,j]+RateSuppFirstResis[2]*x[13,2,2,j]-(RateFailFirstResis[2]+RateStageSuppFirst[1]+mu[4,2,2])*x[14,2,2,j]+RateStageSuppFirst[2]*x[14,3,2,j]
#       dx[15,2,2,j]=dx[15,2,2,j]+RateFailFirstResis[2]*x[14,2,2,j]-(RateTreatSecond_noDTG[2,j]+RateStageFailFirst[2]+mu[5,2,2])*x[15,2,2,j]+RateStageFailFirst[1]*x[15,1,2,j]
#       
#       dx[13,3,2,j]=dx[13,3,2,j]+RateTreatFirst_noDTG[2,3,j]*x[12,3,2,j]-(RateSuppFirstResis[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,2])*x[13,3,2,j]+RateStageTreatFirstR[2]*x[13,2,2,j]+RateStageTreatFirstL[3]*x[13,4,2,j]
#       dx[14,3,2,j]=dx[14,3,2,j]+RateSuppFirstResis[3]*x[13,3,2,j]-(RateFailFirstResis[3]+RateStageSuppFirst[2]+mu[4,3,2])*x[14,3,2,j]+RateStageSuppFirst[3]*x[14,4,2,j]
#       dx[15,3,2,j]=dx[15,3,2,j]+RateFailFirstResis[3]*x[14,3,2,j]-(RateTreatSecond_noDTG[3,j]+RateStageFailFirst[3]+mu[5,3,2])*x[15,3,2,j]+RateStageFailFirst[2]*x[15,2,2,j]
#       
#       dx[13,4,2,j]=dx[13,4,2,j]+RateTreatFirst_noDTG[2,4,j]*x[12,4,2,j]-(RateSuppFirstResis[4]+RateStageTreatFirstL[3]+mu[3,4,2])*x[13,4,2,j]+RateStageTreatFirstR[3]*x[13,3,2,j]
#       dx[14,4,2,j]=dx[14,4,2,j]+RateSuppFirstResis[4]*x[13,4,2,j]-(RateFailFirstResis[4]+RateStageSuppFirst[3]+mu[4,4,2])*x[14,4,2,j]
#       dx[15,4,2,j]=dx[15,4,2,j]+RateFailFirstResis[4]*x[14,4,2,j]-(RateTreatSecond_noDTG[4,j]+mu[5,4,2])*x[15,4,2,j]+RateStageFailFirst[3]*x[15,3,2,j]
#       
#       ##############################################################################################################################################################
#       
#       dx[13,1,1,j]=dx[13,1,1,j]-RateStopTreatFirst[1]*x[13,1,1,j]-RateTreatToFailFirstSusc[1]*x[13,1,1,j]
#       dx[14,1,1,j]=dx[14,1,1,j]-RateStopSuppFirst[1]*x[14,1,1,j]+RateFailToSuppTreatFirstSusc[1]*x[15,1,1,j]-RateSuppFirstToSecond[1]*x[14,1,1,j]
#       dx[15,1,1,j]=dx[15,1,1,j]-RateStopFailFirst[1]*x[15,1,1,j]-RateFailToSuppTreatFirstSusc[1]*x[15,1,1,j]+RateTreatToFailFirstSusc[1]*x[13,1,1,j]
#       
#       dx[13,2,1,j]=dx[13,2,1,j]-RateStopTreatFirst[2]*x[13,2,1,j]-RateTreatToFailFirstSusc[2]*x[13,2,1,j]
#       dx[14,2,1,j]=dx[14,2,1,j]-RateStopSuppFirst[2]*x[14,2,1,j]+RateFailToSuppTreatFirstSusc[2]*x[15,2,1,j]-RateSuppFirstToSecond[2]*x[14,2,1,j]
#       dx[15,2,1,j]=dx[15,2,1,j]-RateStopFailFirst[2]*x[15,2,1,j]-RateFailToSuppTreatFirstSusc[2]*x[15,2,1,j]+RateTreatToFailFirstSusc[2]*x[13,2,1,j]
#       
#       dx[13,3,1,j]=dx[13,3,1,j]-RateStopTreatFirst[3]*x[13,3,1,j]-RateTreatToFailFirstSusc[3]*x[13,3,1,j]
#       dx[14,3,1,j]=dx[14,3,1,j]-RateStopSuppFirst[3]*x[14,3,1,j]+RateFailToSuppTreatFirstSusc[3]*x[15,3,1,j]-RateSuppFirstToSecond[3]*x[14,3,1,j]
#       dx[15,3,1,j]=dx[15,3,1,j]-RateStopFailFirst[3]*x[15,3,1,j]-RateFailToSuppTreatFirstSusc[3]*x[15,3,1,j]+RateTreatToFailFirstSusc[3]*x[13,3,1,j]
#       
#       dx[13,4,1,j]=dx[13,4,1,j]-RateStopTreatFirst[4]*x[13,4,1,j]-RateTreatToFailFirstSusc[4]*x[13,4,1,j]
#       dx[14,4,1,j]=dx[14,4,1,j]-RateStopSuppFirst[4]*x[14,4,1,j]+RateFailToSuppTreatFirstSusc[4]*x[15,4,1,j]-RateSuppFirstToSecond[4]*x[14,4,1,j]
#       dx[15,4,1,j]=dx[15,4,1,j]-RateStopFailFirst[4]*x[15,4,1,j]-RateFailToSuppTreatFirstSusc[4]*x[15,4,1,j]+RateTreatToFailFirstSusc[4]*x[13,4,1,j]
#       #############################################################################################################################################################################
#       #############################################################################################################################################################################
#       
#       dx[13,1,2,j]=dx[13,1,2,j]-RateStopTreatFirst[1]*x[13,1,2,j]-RateTreatToFailFirstResis[1]*x[13,1,2,j]
#       dx[14,1,2,j]=dx[14,1,2,j]-RateStopSuppFirst[1]*x[14,1,2,j]+RateFailToSuppTreatFirstResis[1]*x[15,1,2,j]-RateSuppFirstToSecond[1]*x[14,1,2,j]
#       dx[15,1,2,j]=dx[15,1,2,j]-RateStopFailFirst[1]*x[15,1,2,j]-RateFailToSuppTreatFirstResis[1]*x[15,1,2,j]+RateTreatToFailFirstResis[1]*x[13,1,2,j]
#       
#       dx[13,2,2,j]=dx[13,2,2,j]-RateStopTreatFirst[2]*x[13,2,2,j]-RateTreatToFailFirstResis[2]*x[13,2,2,j]
#       dx[14,2,2,j]=dx[14,2,2,j]-RateStopSuppFirst[2]*x[14,2,2,j]+RateFailToSuppTreatFirstResis[2]*x[15,2,2,j]-RateSuppFirstToSecond[2]*x[14,2,2,j]
#       dx[15,2,2,j]=dx[15,2,2,j]-RateStopFailFirst[2]*x[15,2,2,j]-RateFailToSuppTreatFirstResis[2]*x[15,2,2,j]+RateTreatToFailFirstResis[2]*x[13,2,2,j]
#       
#       dx[13,3,2,j]=dx[13,3,2,j]-RateStopTreatFirst[3]*x[13,3,2,j]-RateTreatToFailFirstResis[3]*x[13,3,2,j]
#       dx[14,3,2,j]=dx[14,3,2,j]-RateStopSuppFirst[3]*x[14,3,2,j]+RateFailToSuppTreatFirstResis[3]*x[15,3,2,j]-RateSuppFirstToSecond[3]*x[14,3,2,j]
#       dx[15,3,2,j]=dx[15,3,2,j]-RateStopFailFirst[3]*x[15,3,2,j]-RateFailToSuppTreatFirstResis[3]*x[15,3,2,j]+RateTreatToFailFirstResis[3]*x[13,3,2,j]
#       
#       dx[13,4,2,j]=dx[13,4,2,j]-RateStopTreatFirst[4]*x[13,4,2,j]-RateTreatToFailFirstResis[4]*x[13,4,2,j]
#       dx[14,4,2,j]=dx[14,4,2,j]-RateStopSuppFirst[4]*x[14,4,2,j]+RateFailToSuppTreatFirstResis[4]*x[15,4,2,j]-RateSuppFirstToSecond[4]*x[14,4,2,j]
#       dx[15,4,2,j]=dx[15,4,2,j]-RateStopFailFirst[4]*x[15,4,2,j]-RateFailToSuppTreatFirstResis[4]*x[15,4,2,j]+RateTreatToFailFirstResis[4]*x[13,4,2,j]
#       
#       #Resistance acquisition and reversion
#       dx[15,1,2,j]=dx[15,1,2,j]+RateResistant*x[15,1,1,j]
#       dx[15,2,2,j]=dx[15,2,2,j]+RateResistant*x[15,2,1,j]
#       dx[15,3,2,j]=dx[15,3,2,j]+RateResistant*x[15,3,1,j]
#       dx[15,4,2,j]=dx[15,4,2,j]+RateResistant*x[15,4,1,j]
#       
#       dx[15,1,1,j]=dx[15,1,1,j]-RateResistant*x[15,1,1,j]
#       dx[15,2,1,j]=dx[15,2,1,j]-RateResistant*x[15,2,1,j]
#       dx[15,3,1,j]=dx[15,3,1,j]-RateResistant*x[15,3,1,j]
#       dx[15,4,1,j]=dx[15,4,1,j]-RateResistant*x[15,4,1,j]
#     }
#     dx[x+dx<0]=0
#     
#     death_w_inel_2018=ifelse(sum(x[13:15,1:4,1:2,2])==0,
#                              0,sum(x[13:15,1:4,1:2,2]*mu[3:5,1:4,1:2])/sum(x[13:15,1:4,1:2,2]))
#     res=c(dx,death_w_inel_2018,rep(0,9))
#     
#     return(c(res))
#     #return(list(c(res,newinf,newsu,newres,death)))
#     #})
#   })
# }
# #18
# mod_dtg=function(t,x,p1,treat_dtg,parms){
#   
#   parms2=mod_rate(t,x,p1,treat_dtg,parms)
#   dx=array(0,dim=c(18,4,2,2))
#   x=array(x,dim=c(18,4,2,2))
#   
#   #Define comp_time globally, used here to fix time from which dtg-ineligible ind go to compartments 16
#   
#   with(as.list(parms2), {
#     for(j in 1:2){
#       dx[16,1,1,j]=dx[16,1,1,j]+RateTreatFirst_noDTG[1,1,j]*x[12,1,1,j]-(RateSuppFirstSusc[1]+RateStageTreatFirstR[1]+mu[3,1,1])*x[16,1,1,j]+RateStageTreatFirstL[1]*x[16,2,1,j]
#       dx[17,1,1,j]=dx[17,1,1,j]+RateSuppFirstSusc[1]*x[16,1,1,j]-(RateFailFirstSusc[1]+mu[4,1,1])*x[17,1,1,j]+RateStageSuppFirst[1]*x[17,2,1,j]
#       dx[18,1,1,j]=dx[18,1,1,j]+RateFailFirstSusc[1]*x[17,1,1,j]-(RateTreatSecond_noDTG[1,j]+RateStageFailFirst[1]+mu[5,1,1])*x[18,1,1,j]
#       
#       dx[16,2,1,j]=dx[16,2,1,j]+RateTreatFirst_noDTG[1,2,j]*x[12,2,1,j]-(RateSuppFirstSusc[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,1])*x[16,2,1,j]+RateStageTreatFirstR[1]*x[16,1,1,j]+RateStageTreatFirstL[2]*x[16,3,1,j]
#       dx[17,2,1,j]=dx[17,2,1,j]+RateSuppFirstSusc[2]*x[16,2,1,j]-(RateFailFirstSusc[2]+RateStageSuppFirst[1]+mu[4,2,1])*x[17,2,1,j]+RateStageSuppFirst[2]*x[17,3,1,j]
#       dx[18,2,1,j]=dx[18,2,1,j]+RateFailFirstSusc[2]*x[17,2,1,j]-(RateTreatSecond_noDTG[2,j]+RateStageFailFirst[2]+mu[5,2,1])*x[18,2,1,j]+RateStageFailFirst[1]*x[18,1,1,j]
#       
#       dx[16,3,1,j]=dx[16,3,1,j]+RateTreatFirst_noDTG[1,3,j]*x[12,3,1,j]-(RateSuppFirstSusc[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,1])*x[16,3,1,j]+RateStageTreatFirstR[2]*x[16,2,1,j]+RateStageTreatFirstL[3]*x[16,4,1,j]
#       dx[17,3,1,j]=dx[17,3,1,j]+RateSuppFirstSusc[3]*x[16,3,1,j]-(RateFailFirstSusc[3]+RateStageSuppFirst[2]+mu[4,3,1])*x[17,3,1,j]+RateStageSuppFirst[3]*x[17,4,1,j]
#       dx[18,3,1,j]=dx[18,3,1,j]+RateFailFirstSusc[3]*x[17,3,1,j]-(RateTreatSecond_noDTG[3,j]+RateStageFailFirst[3]+mu[5,3,1])*x[18,3,1,j]+RateStageFailFirst[2]*x[18,2,1,j]
#       
#       dx[16,4,1,j]=dx[16,4,1,j]+RateTreatFirst_noDTG[1,4,j]*x[12,4,1,j]-(RateSuppFirstSusc[4]+RateStageTreatFirstL[3]+mu[3,4,1])*x[16,4,1,j]+RateStageTreatFirstR[3]*x[16,3,1,j]
#       dx[17,4,1,j]=dx[17,4,1,j]+RateSuppFirstSusc[4]*x[16,4,1,j]-(RateFailFirstSusc[4]+RateStageSuppFirst[3]+mu[4,4,1])*x[17,4,1,j]
#       dx[18,4,1,j]=dx[18,4,1,j]+RateFailFirstSusc[4]*x[17,4,1,j]-(RateTreatSecond_noDTG[4,j]+mu[5,4,1])*x[18,4,1,j]+RateStageFailFirst[3]*x[18,3,1,j]
#       
#       #NNRTI resistant
#       dx[16,1,2,j]=dx[16,1,2,j]+RateTreatFirst_noDTG[2,1,j]*x[12,1,2,j]-(RateSuppFirstResis[1]+RateStageTreatFirstR[1]+mu[3,1,2])*x[16,1,2,j]+RateStageTreatFirstL[1]*x[16,2,2,j]
#       dx[17,1,2,j]=dx[17,1,2,j]+RateSuppFirstResis[1]*x[16,1,2,j]-(RateFailFirstResis[1]+mu[4,1,2])*x[17,1,2,j]+RateStageSuppFirst[1]*x[17,2,2,j]
#       dx[18,1,2,j]=dx[18,1,2,j]+RateFailFirstResis[1]*x[17,1,2,j]-(RateTreatSecond_noDTG[1,j]+RateStageFailFirst[1]+mu[5,1,2])*x[18,1,2,j]
#       
#       dx[16,2,2,j]=dx[16,2,2,j]+RateTreatFirst_noDTG[2,2,j]*x[12,2,2,j]-(RateSuppFirstResis[2]+RateStageTreatFirstR[2]+RateStageTreatFirstL[1]+mu[3,2,2])*x[16,2,2,j]+RateStageTreatFirstR[1]*x[16,1,2,j]+RateStageTreatFirstL[2]*x[16,3,2,j]
#       dx[17,2,2,j]=dx[17,2,2,j]+RateSuppFirstResis[2]*x[16,2,2,j]-(RateFailFirstResis[2]+RateStageSuppFirst[1]+mu[4,2,2])*x[17,2,2,j]+RateStageSuppFirst[2]*x[17,3,2,j]
#       dx[18,2,2,j]=dx[18,2,2,j]+RateFailFirstResis[2]*x[17,2,2,j]-(RateTreatSecond_noDTG[2,j]+RateStageFailFirst[2]+mu[5,2,2])*x[18,2,2,j]+RateStageFailFirst[1]*x[18,1,2,j]
#       
#       dx[16,3,2,j]=dx[16,3,2,j]+RateTreatFirst_noDTG[2,3,j]*x[12,3,2,j]-(RateSuppFirstResis[3]+RateStageTreatFirstR[3]+RateStageTreatFirstL[2]+mu[3,3,2])*x[16,3,2,j]+RateStageTreatFirstR[2]*x[16,2,2,j]+RateStageTreatFirstL[3]*x[16,4,2,j]
#       dx[17,3,2,j]=dx[17,3,2,j]+RateSuppFirstResis[3]*x[16,3,2,j]-(RateFailFirstResis[3]+RateStageSuppFirst[2]+mu[4,3,2])*x[17,3,2,j]+RateStageSuppFirst[3]*x[17,4,2,j]
#       dx[18,3,2,j]=dx[18,3,2,j]+RateFailFirstResis[3]*x[17,3,2,j]-(RateTreatSecond_noDTG[3,j]+RateStageFailFirst[3]+mu[5,3,2])*x[18,3,2,j]+RateStageFailFirst[2]*x[18,2,2,j]
#       
#       dx[16,4,2,j]=dx[16,4,2,j]+RateTreatFirst_noDTG[2,4,j]*x[12,4,2,j]-(RateSuppFirstResis[4]+RateStageTreatFirstL[3]+mu[3,4,2])*x[16,4,2,j]+RateStageTreatFirstR[3]*x[16,3,2,j]
#       dx[17,4,2,j]=dx[17,4,2,j]+RateSuppFirstResis[4]*x[16,4,2,j]-(RateFailFirstResis[4]+RateStageSuppFirst[3]+mu[4,4,2])*x[17,4,2,j]
#       dx[18,4,2,j]=dx[18,4,2,j]+RateFailFirstResis[4]*x[17,4,2,j]-(RateTreatSecond_noDTG[4,j]+mu[5,4,2])*x[18,4,2,j]+RateStageFailFirst[3]*x[18,3,2,j]
#       
#       ##############################################################################################################################################################
#       
#       dx[16,1,1,j]=dx[16,1,1,j]-RateStopTreatFirst[1]*x[16,1,1,j]-RateTreatToFailFirstSusc[1]*x[16,1,1,j]
#       dx[17,1,1,j]=dx[17,1,1,j]-RateStopSuppFirst[1]*x[17,1,1,j]+RateFailToSuppTreatFirstSusc[1]*x[18,1,1,j]-RateSuppFirstToSecond[1]*x[17,1,1,j]
#       dx[18,1,1,j]=dx[18,1,1,j]-RateStopFailFirst[1]*x[18,1,1,j]-RateFailToSuppTreatFirstSusc[1]*x[18,1,1,j]+RateTreatToFailFirstSusc[1]*x[16,1,1,j]
#       
#       dx[16,2,1,j]=dx[16,2,1,j]-RateStopTreatFirst[2]*x[16,2,1,j]-RateTreatToFailFirstSusc[2]*x[16,2,1,j]
#       dx[17,2,1,j]=dx[17,2,1,j]-RateStopSuppFirst[2]*x[17,2,1,j]+RateFailToSuppTreatFirstSusc[2]*x[18,2,1,j]-RateSuppFirstToSecond[2]*x[17,2,1,j]
#       dx[18,2,1,j]=dx[18,2,1,j]-RateStopFailFirst[2]*x[18,2,1,j]-RateFailToSuppTreatFirstSusc[2]*x[18,2,1,j]+RateTreatToFailFirstSusc[2]*x[16,2,1,j]
#       
#       dx[16,3,1,j]=dx[16,3,1,j]-RateStopTreatFirst[3]*x[16,3,1,j]-RateTreatToFailFirstSusc[3]*x[16,3,1,j]
#       dx[17,3,1,j]=dx[17,3,1,j]-RateStopSuppFirst[3]*x[17,3,1,j]+RateFailToSuppTreatFirstSusc[3]*x[18,3,1,j]-RateSuppFirstToSecond[3]*x[17,3,1,j]
#       dx[18,3,1,j]=dx[18,3,1,j]-RateStopFailFirst[3]*x[18,3,1,j]-RateFailToSuppTreatFirstSusc[3]*x[18,3,1,j]+RateTreatToFailFirstSusc[3]*x[16,3,1,j]
#       
#       dx[16,4,1,j]=dx[16,4,1,j]-RateStopTreatFirst[4]*x[16,4,1,j]-RateTreatToFailFirstSusc[4]*x[16,4,1,j]
#       dx[17,4,1,j]=dx[17,4,1,j]-RateStopSuppFirst[4]*x[17,4,1,j]+RateFailToSuppTreatFirstSusc[4]*x[18,4,1,j]-RateSuppFirstToSecond[4]*x[17,4,1,j]
#       dx[18,4,1,j]=dx[18,4,1,j]-RateStopFailFirst[4]*x[18,4,1,j]-RateFailToSuppTreatFirstSusc[4]*x[18,4,1,j]+RateTreatToFailFirstSusc[4]*x[16,4,1,j]
#       
#       #############################################################################################################################################################################
#       #############################################################################################################################################################################
#       
#       dx[16,1,2,j]=dx[16,1,2,j]-RateStopTreatFirst[1]*x[16,1,2,j]-RateTreatToFailFirstResis[1]*x[16,1,2,j]
#       dx[17,1,2,j]=dx[17,1,2,j]-RateStopSuppFirst[1]*x[17,1,2,j]+RateFailToSuppTreatFirstResis[1]*x[18,1,2,j]-RateSuppFirstToSecond[1]*x[17,1,2,j]
#       dx[18,1,2,j]=dx[18,1,2,j]-RateStopFailFirst[1]*x[18,1,2,j]-RateFailToSuppTreatFirstResis[1]*x[18,1,2,j]+RateTreatToFailFirstResis[1]*x[16,1,2,j]
#       
#       dx[16,2,2,j]=dx[16,2,2,j]-RateStopTreatFirst[2]*x[16,2,2,j]-RateTreatToFailFirstResis[2]*x[16,2,2,j]
#       dx[17,2,2,j]=dx[17,2,2,j]-RateStopSuppFirst[2]*x[17,2,2,j]+RateFailToSuppTreatFirstResis[2]*x[18,2,2,j]-RateSuppFirstToSecond[2]*x[17,2,2,j]
#       dx[18,2,2,j]=dx[18,2,2,j]-RateStopFailFirst[2]*x[18,2,2,j]-RateFailToSuppTreatFirstResis[2]*x[18,2,2,j]+RateTreatToFailFirstResis[2]*x[16,2,2,j]
#       
#       dx[16,3,2,j]=dx[16,3,2,j]-RateStopTreatFirst[3]*x[16,3,2,j]-RateTreatToFailFirstResis[3]*x[16,3,2,j]
#       dx[17,3,2,j]=dx[17,3,2,j]-RateStopSuppFirst[3]*x[17,3,2,j]+RateFailToSuppTreatFirstResis[3]*x[18,3,2,j]-RateSuppFirstToSecond[3]*x[17,3,2,j]
#       dx[18,3,2,j]=dx[18,3,2,j]-RateStopFailFirst[3]*x[18,3,2,j]-RateFailToSuppTreatFirstResis[3]*x[18,3,2,j]+RateTreatToFailFirstResis[3]*x[16,3,2,j]
#       
#       dx[16,4,2,j]=dx[16,4,2,j]-RateStopTreatFirst[4]*x[16,4,2,j]-RateTreatToFailFirstResis[4]*x[16,4,2,j]
#       dx[17,4,2,j]=dx[17,4,2,j]-RateStopSuppFirst[4]*x[17,4,2,j]+RateFailToSuppTreatFirstResis[4]*x[18,4,2,j]-RateSuppFirstToSecond[4]*x[17,4,2,j]
#       dx[18,4,2,j]=dx[18,4,2,j]-RateStopFailFirst[4]*x[18,4,2,j]-RateFailToSuppTreatFirstResis[4]*x[18,4,2,j]+RateTreatToFailFirstResis[4]*x[16,4,2,j]
#       
#       
#       #Resistance acquisition and reversion
#       dx[18,1,2,j]=dx[18,1,2,j]+RateResistant*x[18,1,1,j]
#       dx[18,2,2,j]=dx[18,2,2,j]+RateResistant*x[18,2,1,j]
#       dx[18,3,2,j]=dx[18,3,2,j]+RateResistant*x[18,3,1,j]
#       dx[18,4,2,j]=dx[18,4,2,j]+RateResistant*x[18,4,1,j]
#       
#       dx[18,1,1,j]=dx[18,1,1,j]-RateResistant*x[18,1,1,j]
#       dx[18,2,1,j]=dx[18,2,1,j]-RateResistant*x[18,2,1,j]
#       dx[18,3,1,j]=dx[18,3,1,j]-RateResistant*x[18,3,1,j]
#       dx[18,4,1,j]=dx[18,4,1,j]-RateResistant*x[18,4,1,j]
#     }
#     dx[x+dx<0]=0
#     
#     death_w_inel_2018=ifelse(sum(x[16:18,1:4,1:2,2])==0,
#                              0,sum(x[16:18,1:4,1:2,2]*mu[3:5,1:4,1:2])/sum(x[16:18,1:4,1:2,2]))
#     res=c(dx,death_w_inel_2018,rep(0,9))
#     
#     return(c(res))
#     #return(list(c(res,newinf,newsu,newres,death)))
#     #})
#   })
# }

# 
# xstart1=c(runif(240,0,1),rep(0,10))
# ru=runif(1,0,1)
# 
# xstart2=c(rep(0,288),rep(0,10))
# xstart2[select(1:15,1:4,1:2,1:2)]=xstart1[select_dtg_2(1:15,1:4,1:2,1:2)]
# xstart2[select(16:18,1:4,1:2,1:2)]=ru*xstart2[select(13:15,1:4,1:2,1:2)]
# xstart2[select(13:15,1:4,1:2,1:2)]=(1-ru)*xstart2[select(13:15,1:4,1:2,1:2)]
# 
# d1=mod_dtg(50,xstart1,theta,treat_dtgb,params)
# new_inf1=d1[241]
# d1=d1[1:240]
# d2=mod_dtg(50,xstart2,theta,treat_dtgb,params)
# new_inf2=d2[289]
# d2=d2[1:288]
# a2=array(0,dim=c(15,4,2,2))
# a2[13:15,1:4,1:2,1:2]=d2[select(16:18,1:4,1:2,1:2)]
# 
# array(d1[select_dtg_2(1:15,1:4,1:2,1:2)],dim=c(15,4,2,2))-
#   (array(d2[select(1:15,1:4,1:2,1:2)],dim=c(15,4,2,2))+a2)
# new_inf1-new_inf2