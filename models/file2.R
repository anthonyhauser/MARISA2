SIR_NRTI2_scenario=function( t , # time
                             y , # system state { susceptible , infected , recovered }
                             theta , # parameters { transmission rate , recovery rate }
                             treat_all,
                             treat_interruption ) { 
  dtg_1st=theta[1];
  dtg_switch=theta[2];
  dtg_eff=theta[4];
  nnrti_res=1.0/theta[5];
  rev=1.0/theta[6];
  alpha1=theta[7];
  alpha2=theta[8];
  p_msm=theta[9];
  risk_r=theta[10];
  hiv_msm=theta[11];
  nrti_res = 1.0/theta[12];  #time of NRTI resistance acquiring
  alpha3 = theta[13]; #effect of NRTI resistance on DTG-regimen
  var_switch_DTG_S=theta[14]; #Increase or decrease of switching rate (used in the heatmap)
  var_switch_DTG_F=theta[15]; #Increase or decrease of switching rate (used in the heatmap)
  alpha4=theta[16];
  DTG_start_year=theta[17]; #Year of DTG introduction
  switch_decrease=theta[18]; #decrease in switching rate compared to the one from IeDEA-SA
  
  # inf1  = c(0.008*0.05, 0.003*0.95, 0.003*0.95, 0.0);
  corr_factor_incidence=(p_msm*risk_r*hiv_msm+2.0*(1-p_msm))/(0.05*0.008/0.003 + 2.0*0.95)
  inf1 = c( 0.003*risk_r*p_msm*hiv_msm / corr_factor_incidence,
            0.003*(1-p_msm) / corr_factor_incidence,
            0.003*(1-p_msm) / corr_factor_incidence,
            0.0);
  inf2=3.318999;
  inf3=0.5440886;
  all_pos =c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24); #position of all HIV stages
  inf_pos =c(2,3,4,6,7,9,10,12,13,15,16,18,19,21,22,24); #position of infectious HIV stages knowing their HIV-positive status
  #Diagnosis
  diag1 = 273.6363; # number of months before diagnosis (without opportunistic disease - oi) in 2005: 22 years
  diag2 = smooth_f(2005.0 + t/12.0,2005.0,2015.0,1.0,7.710578); # fold increase between 2005 and 2015: 3 years
  diag_women=1.25; #increase of diagnosis rate (without oi) for women
  oi_inc=c(0.05/12.0, 0.12/12.0, 0.27/12.0, 0.9/12.0); #oi incidence by cd4
  oi_test=smooth_f(2005 + t/12.0,2005.0,2015.0,0.2,0.8); #proportion of oi test
  preg_inc =c(1.0 * 23.0/(12.0*1000.0),0.96 * 23.0/(12.0*1000.0),0.87 * 23.0/(12.0*1000.0),0.74 * 23.0/(12.0*1000.0)); #incidence of pregnancy by cd4
  preg_test=smooth_f(2005 + t/12.0,2005.0,2010.0,0.5,0.98); #proportion of pregnancy test
  #Treatment
  p_dtg = c(1.0, theta[3]); #proportion of women opting for DTG
  p_tdf = 1.0; #proportion of people opting/being prescribed tdf (over azt)
  rate_treat = 0.001879695; #free treatment parameter rates fixed in project 1
  t_cd4_elig = c(0.4, 0.5, 0.7, 1.0); #relative initiation rates by cd4 in 2020 (increase to become identical from 2022, Treat-All policy)
  if(treat_all==1){t_cd4_treat_all = c(1.0/0.4, 1.0/0.5, 1.0/0.7, 1.0)}else{t_cd4_treat_all = c(1.0, 1.0, 1.0, 1.0)};
  t_cd4_elig_year1 = c(2015.5, 2013.5, 2009.5, 2001.0);
  t_cd4_elig_year2 = c(2016.5, 2015.5, 2012.5, 2004.0);
  t_1st =c();
  for(i in 1:4){
    t_1st[i] = rate_treat *smooth_f(2005.0 + t/12.0, 2001.0, 2012.0, 1.0, 200.0/12.0) *  #Increase in treatment rate between 2001 and 2012
      lin_f(2005.0 + t/12.0, t_cd4_elig_year1[i],  t_cd4_elig_year2[i], 0.0, t_cd4_elig[i]) *  #Change in eligiblity criteria
      lin_f(2005.0 + t/12.0, 2017.0, 2022.0, 1.0, t_cd4_treat_all[i]);  #Increase in treatment rate due to Treat-All policy
  }
  t_switch = c(4.35/531.0 /switch_decrease, 4.35/427.0 /switch_decrease, 4.35/294.0 /switch_decrease, 4.35/189.0 /switch_decrease);#Switching rate to PI, failing individuals (either DTG-ineligible on NNRTI or DTG-eligible on DTG)
  t_switch_elig =c(); #Switch rate to PI, DTG-eligible
  t_1st_NNRTI_inel=c();
  t_1st_NNRTI_elig=c();
  t_1st_DTG=c();
  t_switch_DTG_S=c();
  t_switch_DTG_F=c();
  for(i in 1:4){
    t_1st_DTG [i]= t_1st[i] * smooth_f(2005.0 + t/12.0, DTG_start_year, DTG_start_year, 0.0, 1.0) * dtg_1st; #Treatment initiation rate, DTG
    t_1st_NNRTI_inel[i] = t_1st[i];
    t_1st_NNRTI_elig[i] = t_1st[i] - t_1st_DTG[i];
    t_switch_elig[i] = t_switch[i] * smooth_f(2005.0 + t/12.0, DTG_start_year, DTG_start_year, 1.0, 1.0 - dtg_1st); #switch always 0 after 2019, except if dtg_1st=0 (no DTG introduction) in which case it remains at 1
    t_switch_DTG_S[i] = 1.0/12.0 * dtg_switch * var_switch_DTG_S * smooth_f(2005.0 + t/12.0, DTG_start_year, DTG_start_year, 0.0, 1.0);
    if(dtg_1st==0.0) t_switch_DTG_F[i] = 0.0;
    if(dtg_1st==1.0 & dtg_switch==0.0) t_switch_DTG_F[i] = t_switch[i];
    if(dtg_switch==1.0) t_switch_DTG_F[i] = 1.0/12.0 * var_switch_DTG_F;
    t_switch_DTG_F[i] = t_switch_DTG_F[i] * smooth_f(2005.0 + t/12.0, DTG_start_year, DTG_start_year, 0.0, 1.0);
  }
  
  #Mortality
  rate_death =0.1632000; #Mortality rate per month per 1000 people (reference, supp people with cd4>500)
  #p1 =0.865181; #proportion of people with cd4>50 among people with cd4<200
  t_prov=c()
  mort_approx=c()
  for(l in 1:37){
    t_prov=t-l+1.0
    mort_approx[l] = rate_treat * smooth_f(2005.0 + t_prov/12.0, 2001.0, 2012.0, 1.0, 200.0/12.0) *  #Increase in treatment rate between 2001 and 2012
      lin_f(2005.0 + t_prov/12.0, t_cd4_elig_year1[4],  t_cd4_elig_year2[4], 0.0, t_cd4_elig[4]) *  #Change in eligiblity criteria
      lin_f(2005.0 + t_prov/12.0, 2017.0, 2022.0, 1.0, t_cd4_treat_all[4])
  }
  p1=1.0-0.27*exp(-0.05*(mean(mort_approx)-0.005219697)/rate_treat)
  
  mu = matrix(c(1.57 *rate_death/1000.0, 2.0 *rate_death/1000.0, 4.57 *rate_death/1000.0, (p1*40.9+(1-p1)*134.4) *rate_death/1000.0, #relative mortality for untreated by cd4
                (0.9 * 1 + 0.1 * 3.92) * rate_death/1000.0, (0.9 * 1.26 + 0.1 * 3.92) * rate_death/1000.0, #for people starting treatment
                (0.9 * 1.94 + 0.1 * 4.28) * rate_death/1000.0, (0.9 * (p1*8.3+(1-p1)*41.7) + 0.1 * (p1*11.8+(1-p1)*59.7)) * rate_death/1000.0,
                1.0 * rate_death/1000.0, 1.26 * rate_death/1000.0, 1.94 * rate_death/1000.0, (p1*8.3+(1-p1)*41.7) * rate_death/1000.0,  #for suppressed people
                3.92 * rate_death/1000.0, 3.92 * rate_death/1000.0, 4.28 * rate_death/1000.0, (p1*11.8+(1-p1)*59.7) * rate_death/1000.0),nrow=4,byrow=T);#for people failing treatment
  #CD4 progression, for untreated, treated with NNRTI, treated with PI
  cd4_untreated = c(4.35 / 260.0, 4.35 / 156.0, 4.35 / 182.0); #untreated
  cd4_j =c()
  cd4  = matrix(c(4.35 / 206.0, 4.35 / 133.0, 4.35 / 263.0, #NNRTI, start treatment, right
                  4.35 / 69.0, 4.35 / 70.0, 4.35 / 79.0,#start treatment, left
                  4.35 / 73.0, 4.35 / 61.0, 4.35 / 41.0,#suppressed
                  4.35 / 77.0, 4.35 / 66.0, 4.35 / 96.0,#failing
                  4.35 / 140.0, 4.35 / 98.0, 4.35 / 145.0, #PI,start treatment, right
                  4.35 / 68.0, 4.35 / 81.0, 4.35 / 177.0,#start treatment, left
                  4.35 / 72.0, 4.35 / 59.0, 4.35 / 32.0,#suppressed
                  4.35 / 61.0, 4.35 / 65.0, 4.35 / 69.0),nrow=2,byrow=T);#failing
  #Treatment suppression and failure rates, for NNRTI and PI, by cd4
  treat_j =c()
  treat = matrix(c(2.0 * 4.35 / 30.0, 2.0 * 4.35 / 30.0, 2.0 * 4.35 / 31.0, 2.0 * 4.35 / 34.0, #NNRTI.0, T to S
                   2.0 * 4.35 / 203.0, 2.0 * 4.35 / 198.0, 2.0 * 4.35 / 164.0, 2.0 * 4.35 / 112.0, #T to F
                   4.35 / 767.0, 4.35 / 582.0, 4.35 / 270.0, 4.35 / 96.0, #S to F
                   4.35 / 28.0, 4.35 / 56.0, 4.35 / 62.0, 4.35 / 79.0, #F to S
                   2.0 * 4.35 / 33.0, 2.0 * 4.35 / 33.0, 2.0 * 4.35 / 35.0, 2.0 * 4.35 / 43.0, #PI
                   2.0 * 4.35 / 124.0, 2.0 * 4.35 / 122.0, 2.0 * 4.35 / 103.0, 2.0 * 4.35 / 66.0, #T to F
                   4.35 / 267.0, 4.35 / 178.0, 4.35 / 174.0, 4.35 / 83.0, #S to F
                   4.35 / 10.0, 4.35 / 56.0, 4.35 / 24.0, 4.35 / 51.0),nrow=2,byrow=T); #F to S
  if(treat_interruption==0){
    treat_interruption_rate=matrix(rep(0,4*3),nrow=3,byrow=T)
  }else{
    treat_interruption_rate=matrix(c(4.35 / 1800, 4.35 / 1400, 4.35 / 750, 4.35 / 680,
                                     4.35 / 9000, 4.35 / 5400 , 4.35 / 3300, 4.35 / 1600,
                                     4.35 / 2700, 4.35 / 2080, 4.35 / 1240, 4.35 / 560),nrow=3,byrow=T)
  }
  
  #CD4 progression and ART efficacy by ART, 1: cd4 progression/efficacy equivalent to NNRTI, 2: equivalent to PI
  art_eq  = c(1, 1, 1, 1, 1, 1, 2);
  art_j =c()
  #Increase/decrease of efficacy due to NNRTI resistance by ART
  art_res_nnrti1 = matrix(c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                            alpha1, alpha1, alpha1, alpha1, 1.0, 1.0, 1.0),nrow=2,byrow=T);
  art_res_nnrti2 = matrix(c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                            alpha2, alpha2, alpha2, alpha2, 1.0, 1.0, 1.0),nrow=2,byrow=T);
  art_res_nrti = matrix(c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,  1.0, 1.0, 1.0, 1.0, alpha3, alpha3, 1.0),nrow=2,byrow=T);
  #Increase/decrease of efficacy relative to NNRTI (for NNRTI and DTG) and to PI (for PI)
  art_eff = c(alpha4, alpha4, alpha4, alpha4, alpha4 * dtg_eff, alpha4 * dtg_eff, 1.0);
  
  
  #N values by gender, used to defined susceptible by gender y[1], y[2]
  N_value=c(39845943.06,40029729.5,40544491.9,41221762.41,42376196.21,43678397.16,
            44530813.78,44878159.38,46527810.82,47609855.86,48395856.04,49250847.33,
            50105838.63,50960829.93)/1000
  
  N_v<-smooth.spline(1:14,N_value,df=5,all.knots=FALSE)
  N_value_f=function(x){
    return(predict(N_v,x/12)$y*c(0.5,0.5))
  }
  N<-N_value_f(t+1)
  
  #susceptible by gender
  y[1]=N[1]-sum(y[posa(all_pos,1:4,1,1:2,1:2)]);
  y[2]=N[2]-sum(y[posa(all_pos,1:4,2,1:2,1:2)]);
  
  
  #Differential equations
  dy_dt=rep(0,24*4*2*2*2+2)
  dy_dt_dummy = rep(0,4*2)
  #Susceptible, dydt=0 as values over time fixed by the spline function 
  dy_dt[1]=0;
  dy_dt[2]=0;
  
  for(l in 1:2){
    dy_dt[pos(1,1,1,l,1)]=dy_dt[pos(1,1,1,l,1)] +
      inf1[1] * inf2 * (sum(y[posa(1,1:4,1,l,1:2)]) + inf3 * sum(y[posa(inf_pos,1:4,1,l,1:2)])) * y[1] / (y[1] + sum(y[posa(all_pos,1:4,1,1:2,1:2)])) +
      inf1[2] * inf2 * (sum(y[posa(1,1:4,2,l,1:2)]) + inf3 * sum(y[posa(inf_pos,1:4,2,l,1:2)])) * y[1] / (y[1] + sum(y[posa(all_pos,1:4,1,1:2,1:2)]));
    dy_dt[pos(1,1,2,l,1)]=dy_dt[pos(1,1,2,l,1)] +
      inf1[3] * inf2 * (sum(y[posa(1,1:4,1,l,1:2)]) + inf3 * sum(y[posa(inf_pos,1:4,1,l,1:2)])) * y[2] / (y[2] + sum(y[posa(all_pos,1:4,2,1:2,1:2)])) +
      inf1[4] * inf2 * (sum(y[posa(1,1:4,2,l,1:2)]) + inf3 * sum(y[posa(inf_pos,1:4,2,l,1:2)])) * y[2] / (y[2] + sum(y[posa(all_pos,1:4,2,1:2,1:2)]));
  }
  
  for(j in 1:4){
    for(l in 1:2){
      for(m in 1:2){
        #Diagnosis
        #men
        dy_dt[pos(1,j,1,l,m)] = dy_dt[pos(1,j,1,l,m)] - (diag2 / diag1 + oi_inc[j] * oi_test) * y[pos(1,j,1,l,m)];
        #women
        dy_dt[pos(1,j,2,l,m)] = dy_dt[pos(1,j,2,l,m)] - (diag2 / diag1 * diag_women + oi_inc[j] * oi_test + preg_inc[j] * preg_test) * y[pos(1,j,2,l,m)];
        
        #DTG-eligible, men
        dy_dt[pos(2,j,1,l,m)] = dy_dt[pos(2,j,1,l,m)] + p_dtg[1] * (diag2 / diag1 + oi_inc[j] * oi_test) * y[pos(1,j,1,l,m)];
        #DTG-ineligible, men
        dy_dt[pos(3,j,1,l,m)] = dy_dt[pos(3,j,1,l,m)] + (1 - p_dtg[1]) * (diag2 / diag1 + oi_inc[j] * oi_test) * y[pos(1,j,1,l,m)];
        #DTG-eligible, women
        dy_dt[pos(2,j,2,l,m)] = dy_dt[pos(2,j,2,l,m)] + p_dtg[2] * (diag2 / diag1 * diag_women + oi_inc[j] * oi_test + preg_inc[j] * preg_test) * y[pos(1,j,2,l,m)];
        #DTG-ineligible, women
        dy_dt[pos(3,j,2,l,m)] = dy_dt[pos(3,j,2,l,m)] + (1 - p_dtg[2]) * (diag2 / diag1 * diag_women + oi_inc[j] * oi_test + preg_inc[j] * preg_test) * y[pos(1,j,2,l,m)];
        
        
        for(k in 1:2){
          #Treatment initiation
          #NNRTI: 1) NNRTI+TDF for DTG-inel, 2) NNRTI+TDF for DTG-elig, 3) NNRTI+AZT for DTG-inel, 4) NNRTI+AZT for DTG-elig
          dy_dt[pos(4,j,k,l,m)] = dy_dt[pos(4,j,k,l,m)] + p_tdf * t_1st_NNRTI_inel[j] * y[pos(3,j,k,l,m)];
          dy_dt[pos(7,j,k,l,m)] = dy_dt[pos(7,j,k,l,m)] + p_tdf * t_1st_NNRTI_elig[j] * y[pos(2,j,k,l,m)];
          dy_dt[pos(10,j,k,l,m)] = dy_dt[pos(10,j,k,l,m)] + (1-p_tdf) * t_1st_NNRTI_inel[j] * y[pos(3,j,k,l,m)];
          dy_dt[pos(13,j,k,l,m)] = dy_dt[pos(13,j,k,l,m)] + (1-p_tdf) * t_1st_NNRTI_elig[j] * y[pos(2,j,k,l,m)];
          #DTG: DTG+TDF
          dy_dt[pos(16,j,k,l,m)] = dy_dt[pos(16,j,k,l,m)] + t_1st_DTG[j] * y[pos(2,j,k,l,m)];
          
          dy_dt[pos(3,j,k,l,m)] = dy_dt[pos(3,j,k,l,m)] - t_1st_NNRTI_inel[j] * y[pos(3,j,k,l,m)];
          dy_dt[pos(2,j,k,l,m)] = dy_dt[pos(2,j,k,l,m)] - (t_1st_NNRTI_elig[j] + t_1st_DTG[j]) * y[pos(2,j,k,l,m)];
          
          #Dummy compartments, summing the number of people initiating DTG, according to cd4 and NNRTI resistance
          dy_dt_dummy[j + (l-1)*4] = dy_dt_dummy[j + (l-1)*4] + t_1st_NNRTI_inel[j] * y[pos(3,j,k,l,m)];
          
          #Treatment switch to PI after failure
          #Switch from NNRTI to PI, DTG-ineligible
          dy_dt[pos(22,j,k,l,m)] = dy_dt[pos(22,j,k,l,m)] + t_switch[j] * (y[pos(6,j,k,l,m)] + y[pos(12,j,k,l,m)]);
          dy_dt[pos(6,j,k,l,m)] = dy_dt[pos(6,j,k,l,m)] - t_switch[j] * y[pos(6,j,k,l,m)];
          dy_dt[pos(12,j,k,l,m)] = dy_dt[pos(12,j,k,l,m)] - t_switch[j] * y[pos(12,j,k,l,m)];
          #Switch from DTG to PI
          dy_dt[pos(22,j,k,l,m)] = dy_dt[pos(22,j,k,l,m)] + t_switch[j] * (y[pos(18,j,k,l,m)] + y[pos(21,j,k,l,m)]);
          dy_dt[pos(18,j,k,l,m)] = dy_dt[pos(18,j,k,l,m)] - t_switch[j] * y[pos(18,j,k,l,m)];
          dy_dt[pos(21,j,k,l,m)] = dy_dt[pos(21,j,k,l,m)] - t_switch[j] * y[pos(21,j,k,l,m)];
          #Switch from NNRTI to PI, DTG-eligible
          dy_dt[pos(22,j,k,l,m)] = dy_dt[pos(22,j,k,l,m)] + t_switch_elig[j] * (y[pos(9,j,k,l,m)] + y[pos(15,j,k,l,m)]);
          dy_dt[pos(9,j,k,l,m)] = dy_dt[pos(9,j,k,l,m)] - t_switch_elig[j] * y[pos(9,j,k,l,m)];
          dy_dt[pos(15,j,k,l,m)] = dy_dt[pos(15,j,k,l,m)] - t_switch_elig[j] * y[pos(15,j,k,l,m)];
          
          #Treatment switch to DTG, suppressed
          dy_dt[pos(17,j,k,l,m)] = dy_dt[pos(17,j,k,l,m)] + t_switch_DTG_S[j] * (y[pos(8,j,k,l,m)] + y[pos(14,j,k,l,m)]);
          dy_dt[pos(8,j,k,l,m)] = dy_dt[pos(8,j,k,l,m)] - t_switch_DTG_S[j] * y[pos(8,j,k,l,m)];
          dy_dt[pos(14,j,k,l,m)] = dy_dt[pos(14,j,k,l,m)] - t_switch_DTG_S[j] * y[pos(14,j,k,l,m)];
          
          #Treatment switch to DTG, failed
          dy_dt[pos(16,j,k,l,m)] = dy_dt[pos(16,j,k,l,m)] + t_switch_DTG_F[j] * y[pos(15,j,k,l,m)];
          dy_dt[pos(19,j,k,l,m)] = dy_dt[pos(19,j,k,l,m)] + t_switch_DTG_F[j] * y[pos(9,j,k,l,m)];
          
          dy_dt[pos(15,j,k,l,m)] = dy_dt[pos(15,j,k,l,m)] - t_switch_DTG_F[j] * y[pos(15,j,k,l,m)];
          dy_dt[pos(9,j,k,l,m)] = dy_dt[pos(9,j,k,l,m)] - t_switch_DTG_F[j] * y[pos(9,j,k,l,m)];
        }
      }
    }
  }
  
  #Treatment stages and CD4, mortality
  #Infected and diagnosed
  for(k in 1:2){
    for(l in 1:2){
      for(i in 1:3){
        for(m in 1:2){
          dy_dt[pos(i,1,k,l,m)] = dy_dt[pos(i,1,k,l,m)] - cd4_untreated[1] * y[pos(i,1,k,l,m)];
          dy_dt[pos(i,2,k,l,m)] = dy_dt[pos(i,2,k,l,m)] + cd4_untreated[1] * y[pos(i,1,k,l,m)] - cd4_untreated[2] * y[pos(i,2,k,l,m)];
          dy_dt[pos(i,3,k,l,m)] = dy_dt[pos(i,3,k,l,m)] + cd4_untreated[2] * y[pos(i,2,k,l,m)] - cd4_untreated[3] * y[pos(i,3,k,l,m)];
          dy_dt[pos(i,4,k,l,m)] = dy_dt[pos(i,4,k,l,m)] + cd4_untreated[3] * y[pos(i,3,k,l,m)];
          
          #Mortality
          for(j in 1:4){
            dy_dt[pos(i,j,k,l,m)] = dy_dt[pos(i,j,k,l,m)] - mu[1,j] * y[pos(i,j,k,l,m)];
          }
        }
      }
    }
  }
  
  #Treatment interruption
  for(j in 1:4){
    for(k in 1:2){
      for(l in 1:2){
        for(m in 1:2){
          #DTG-ineligible
          for(art in c(1,3)){
            dy_dt[pos(3,j,k,l,m)] = treat_interruption_rate[1,j]* y[pos(3*art+1,j,k,l,m)] +
              treat_interruption_rate[2,j]* y[pos(3*art+2,j,k,l,m)] +
              treat_interruption_rate[3,j]* y[pos(3*art+3,j,k,l,m)]
          }
          #DTG-eligible
          for(art in c(2,4,5,6)){
            dy_dt[pos(2,j,k,l,m)] = treat_interruption_rate[1,j]* y[pos(3*art+1,j,k,l,m)] +
              treat_interruption_rate[2,j]* y[pos(3*art+2,j,k,l,m)] +
              treat_interruption_rate[3,j]* y[pos(3*art+3,j,k,l,m)]
          }
          #PI-based regimen
          dy_dt[pos(3,j,k,l,m)] = (1 - p_dtg[k]) * treat_interruption_rate[1,j]* y[pos(3*7+1,j,k,l,m)]
          dy_dt[pos(2,j,k,l,m)] = p_dtg[k] * treat_interruption_rate[1,j]* y[pos(3*7+1,j,k,l,m)]
        }
      }
    }
  }
  
  #On ART
  for(art in 1:7){
    art_j=art_eq[art];
    for(j in 1:12){
      cd4_j[j] = cd4[art_j,j];
    }
    for(l in 1:2){
      for(m in 1:2){
        treat_j[1:4]=rep(0,4);
        for(i in 5:8){#T to F
          treat_j[i] =  treat[art_j,i] / art_eff[art] * art_res_nnrti1[l,art] * art_res_nrti[m,art];
        }
        for(i in 1:4){#T to S
          treat_j[i] =  1.0/3.0 - treat_j[i+4];
          #treat_j[i] =  treat[art_j,i] * art_eff[art] / art_res_nnrti1[l,art] / art_res_nrti[m,art];
        }
        for(i in 9:12){#S to F
          treat_j[i] =  treat[art_j,i] / art_eff[art] * art_res_nnrti2[l,art] * art_res_nrti[m,art];
        }
        for(i in 13:16){#F to S
          treat_j[i] =  treat[art_j,i] * art_eff[art] / art_res_nnrti2[l,art] / art_res_nrti[m,art];
        }
        
        for(k in 1:2){
          #Start
          dy_dt[pos(3*art+1,1,k,l,m)] = dy_dt[pos(3*art+1,1,k,l,m)] + cd4_j[4] * y[pos(3*art+1,2,k,l,m)] - cd4_j[1] * y[pos(3*art+1,1,k,l,m)]+
            -(treat_j[1] + treat_j[5]) * y[pos(3*art+1,1,k,l,m)]+
            - treat_interruption_rate[1,1]* y[pos(3*art+1,1,k,l,m)];
          dy_dt[pos(3*art+1,2,k,l,m)] = dy_dt[pos(3*art+1,2,k,l,m)] + cd4_j[1] * y[pos(3*art+1,1,k,l,m)] + cd4_j[5] * y[pos(3*art+1,3,k,l,m)] - (cd4_j[4] + cd4_j[2]) * y[pos(3*art+1,2,k,l,m)]+
            -(treat_j[2] + treat_j[6]) * y[pos(3*art+1,2,k,l,m)]+
            - treat_interruption_rate[1,2]* y[pos(3*art+1,2,k,l,m)];
          dy_dt[pos(3*art+1,3,k,l,m)] = dy_dt[pos(3*art+1,3,k,l,m)] + cd4_j[2] * y[pos(3*art+1,2,k,l,m)] + cd4_j[6] * y[pos(3*art+1,4,k,l,m)] - (cd4_j[5] + cd4_j[3]) * y[pos(3*art+1,3,k,l,m)]+
            -(treat_j[3] + treat_j[7]) * y[pos(3*art+1,3,k,l,m)]+
            - treat_interruption_rate[1,3]* y[pos(3*art+1,3,k,l,m)];
          dy_dt[pos(3*art+1,4,k,l,m)] = dy_dt[pos(3*art+1,4,k,l,m)] + cd4_j[3] * y[pos(3*art+1,3,k,l,m)] - cd4_j[6] * y[pos(3*art+1,4,k,l,m)]+
            -(treat_j[4] + treat_j[8]) * y[pos(3*art+1,4,k,l,m)]+
            - treat_interruption_rate[1,4]* y[pos(3*art+1,4,k,l,m)];
          #Suppressed
          dy_dt[pos(3*art+2,1,k,l,m)] = dy_dt[pos(3*art+2,1,k,l,m)] + cd4_j[7] * y[pos(3*art+2,2,k,l,m)]+
            + treat_j[1] * y[pos(3*art+1,1,k,l,m)] + treat_j[13] * y[pos(3*art+3,1,k,l,m)] - treat_j[9] * y[pos(3*art+2,1,k,l,m)]+
            - treat_interruption_rate[2,1]* y[pos(3*art+2,1,k,l,m)];
          dy_dt[pos(3*art+2,2,k,l,m)] = dy_dt[pos(3*art+2,2,k,l,m)] + cd4_j[8] * y[pos(3*art+2,3,k,l,m)] - cd4_j[7] * y[pos(3*art+2,2,k,l,m)]+
            + treat_j[2] * y[pos(3*art+1,2,k,l,m)] + treat_j[14] * y[pos(3*art+3,2,k,l,m)] - treat_j[10] * y[pos(3*art+2,2,k,l,m)]
          - treat_interruption_rate[2,2]* y[pos(3*art+2,2,k,l,m)];
          dy_dt[pos(3*art+2,3,k,l,m)] = dy_dt[pos(3*art+2,3,k,l,m)] + cd4_j[9] * y[pos(3*art+2,4,k,l,m)] - cd4_j[8] * y[pos(3*art+2,3,k,l,m)]+
            + treat_j[3] * y[pos(3*art+1,3,k,l,m)] + treat_j[15] * y[pos(3*art+3,3,k,l,m)] - treat_j[11] * y[pos(3*art+2,3,k,l,m)]
          - treat_interruption_rate[2,2]* y[pos(3*art+2,3,k,l,m)];
          dy_dt[pos(3*art+2,4,k,l,m)] = dy_dt[pos(3*art+2,4,k,l,m)] - cd4_j[9] * y[pos(3*art+2,4,k,l,m)]+
            + treat_j[4] * y[pos(3*art+1,4,k,l,m)] + treat_j[16] * y[pos(3*art+3,4,k,l,m)] - treat_j[12] * y[pos(3*art+2,4,k,l,m)]
          - treat_interruption_rate[2,2]* y[pos(3*art+2,4,k,l,m)];
          #Failed
          dy_dt[pos(3*art+3,1,k,l,m)] = dy_dt[pos(3*art+3,1,k,l,m)] - cd4_j[10] * y[pos(3*art+3,1,k,l,m)]+
            treat_j[5] * y[pos(3*art+1,1,k,l,m)] + treat_j[9] * y[pos(3*art+2,1,k,l,m)] - treat_j[13] * y[pos(3*art+3,1,k,l,m)]+
            - treat_interruption_rate[3,1]* y[pos(3*art+3,1,k,l,m)];
          dy_dt[pos(3*art+3,2,k,l,m)] = dy_dt[pos(3*art+3,2,k,l,m)] + cd4_j[10] * y[pos(3*art+3,1,k,l,m)] - cd4_j[11] * y[pos(3*art+3,2,k,l,m)]+
            treat_j[6] * y[pos(3*art+1,2,k,l,m)] + treat_j[10] * y[pos(3*art+2,2,k,l,m)] - treat_j[14] * y[pos(3*art+3,2,k,l,m)]+
            - treat_interruption_rate[3,2]* y[pos(3*art+3,2,k,l,m)];
          dy_dt[pos(3*art+3,3,k,l,m)] = dy_dt[pos(3*art+3,3,k,l,m)] + cd4_j[11] * y[pos(3*art+3,2,k,l,m)] - cd4_j[12] * y[pos(3*art+3,3,k,l,m)]+
            treat_j[7] * y[pos(3*art+1,3,k,l,m)] + treat_j[11] * y[pos(3*art+2,3,k,l,m)] - treat_j[15] * y[pos(3*art+3,3,k,l,m)]+
            - treat_interruption_rate[3,3]* y[pos(3*art+3,3,k,l,m)];
          dy_dt[pos(3*art+3,4,k,l,m)] = dy_dt[pos(3*art+3,4,k,l,m)] + cd4_j[12] * y[pos(3*art+3,3,k,l,m)]+
            treat_j[8] * y[pos(3*art+1,4,k,l,m)] + treat_j[12] * y[pos(3*art+2,4,k,l,m)] - treat_j[16] * y[pos(3*art+3,4,k,l,m)]+
            - treat_interruption_rate[3,4]* y[pos(3*art+3,4,k,l,m)];
          
          #Mortality 
          for(j in 1:4){
            dy_dt[pos(3*art+1,j,k,l,m)] = dy_dt[pos(3*art+1,j,k,l,m)] - mu[2,j] * y[pos(3*art+1,j,k,l,m)];
            dy_dt[pos(3*art+2,j,k,l,m)] = dy_dt[pos(3*art+2,j,k,l,m)] - mu[3,j] * y[pos(3*art+2,j,k,l,m)];
            dy_dt[pos(3*art+3,j,k,l,m)] = dy_dt[pos(3*art+3,j,k,l,m)] - mu[4,j] * y[pos(3*art+3,j,k,l,m)];
          }
        }
      }
    }
  }
  
  for(j in 1:4){
    for(k in 1:2){
      for(m in 1:2){
        #nnrti resistance acquiring
        dy_dt[pos(6,j,k,2,m)] = dy_dt[pos(6,j,k,2,m)] + nnrti_res * y[pos(6,j,k,1,m)];
        dy_dt[pos(9,j,k,2,m)] = dy_dt[pos(9,j,k,2,m)] + nnrti_res * y[pos(9,j,k,1,m)];
        dy_dt[pos(12,j,k,2,m)] = dy_dt[pos(12,j,k,2,m)] + nnrti_res * y[pos(12,j,k,1,m)];
        dy_dt[pos(15,j,k,2,m)] = dy_dt[pos(15,j,k,2,m)] + nnrti_res * y[pos(15,j,k,1,m)];
        
        dy_dt[pos(6,j,k,1,m)] = dy_dt[pos(6,j,k,1,m)] - nnrti_res * y[pos(6,j,k,1,m)];
        dy_dt[pos(9,j,k,1,m)] = dy_dt[pos(9,j,k,1,m)] - nnrti_res * y[pos(9,j,k,1,m)];
        dy_dt[pos(12,j,k,1,m)] = dy_dt[pos(12,j,k,1,m)] - nnrti_res * y[pos(12,j,k,1,m)];
        dy_dt[pos(15,j,k,1,m)] = dy_dt[pos(15,j,k,1,m)] - nnrti_res * y[pos(15,j,k,1,m)];
        #nnrti resistance reversion
        dy_dt[pos(1,j,k,1,m)] = dy_dt[pos(1,j,k,1,m)] + rev * y[pos(1,j,k,2,m)];
        dy_dt[pos(2,j,k,1,m)] = dy_dt[pos(2,j,k,1,m)] + rev * y[pos(2,j,k,2,m)];
        dy_dt[pos(3,j,k,1,m)] = dy_dt[pos(3,j,k,1,m)] + rev * y[pos(3,j,k,2,m)];
        
        dy_dt[pos(1,j,k,2,m)] = dy_dt[pos(1,j,k,2,m)] - rev * y[pos(1,j,k,2,m)];
        dy_dt[pos(2,j,k,2,m)] = dy_dt[pos(2,j,k,2,m)] - rev * y[pos(2,j,k,2,m)];
        dy_dt[pos(3,j,k,2,m)] = dy_dt[pos(3,j,k,2,m)] - rev * y[pos(3,j,k,2,m)];
      }
      for(l in 1:2){
        #nrti resistance acquiring
        dy_dt[pos(6,j,k,l,2)] = dy_dt[pos(6,j,k,l,2)] + nrti_res * y[pos(6,j,k,l,1)];
        dy_dt[pos(9,j,k,l,2)] = dy_dt[pos(9,j,k,l,2)] + nrti_res * y[pos(9,j,k,l,1)];
        dy_dt[pos(12,j,k,l,2)] = dy_dt[pos(12,j,k,l,2)] + nrti_res * y[pos(12,j,k,l,1)];
        dy_dt[pos(15,j,k,l,2)] = dy_dt[pos(15,j,k,l,2)] + nrti_res * y[pos(15,j,k,l,1)];
        
        dy_dt[pos(6,j,k,l,1)] = dy_dt[pos(6,j,k,l,1)] - nrti_res * y[pos(6,j,k,l,1)];
        dy_dt[pos(9,j,k,l,1)] = dy_dt[pos(9,j,k,l,1)] - nrti_res * y[pos(9,j,k,l,1)];
        dy_dt[pos(12,j,k,l,1)] = dy_dt[pos(12,j,k,l,1)] - nrti_res * y[pos(12,j,k,l,1)];
        dy_dt[pos(15,j,k,l,1)] = dy_dt[pos(15,j,k,l,1)] - nrti_res * y[pos(15,j,k,l,1)];
      }
    }
  }
  return(list(c(dy_dt,dy_dt_dummy)))
}