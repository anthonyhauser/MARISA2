smooth_f=function(x, a, b, c, d){
  if(a==b){
    if(x<a){
      return(c) 
    }else{
      return(d) 
    }
  }else{
    if(x<a){
      return(c) 
    }else if(x>b){
      return(d) 
    }else{
      return(c+(d-c)*(3*((x-a)/(b-a))^2-2*((x-a)/(b-a))^3))
    }
  }
}

lin_f=function(x,a,b,c,d){
  if(a==b){
    if(x<a){
      return(c)
    }else{
      return(d)
    }
  }else{
    if(x<a){
      return(c)
    }else if(x>b){
      return(d)
    }else{
      return(c+(d-c)*(x-a)/(b-a))
    }
  }
}

pos=function(i, j, k, l){
  return(2 + i+(j-1)*24+(k-1)*24*4+(l-1)*24*4*2)
}

posa=function(v1, v2, v3, v4){
  r = rep(0,length(v1)*length(v2)*length(v3)*length(v4))
  for(i in 1:length(v1)){
    for(j in 1:length(v2)){
      for(k in 1:length(v3)){
        for(l in 1:length(v4)){
          r[i+(j-1)*length(v1)+(k-1)*length(v1)*length(v2)+(l-1)*length(v1)*length(v2)*length(v3)]= 2 + v1[i]+(v2[j]-1)*24+(v3[k]-1)*24*4+(v4[l]-1)*24*4*2 
        }
      }
    }
  }
  return(r)
}

repart_f=function(k){
  x1=1/(1+k+k^2+k^3)
  x2=k*x1
  x3=k*x2
  x4=k*x3
  return(c(x1,x2,x3,x4))
}

x_start2005=function(){
  #Previously estimated parameters
  rate_ratio=0.5
  k1=1
  k2=2
  k3=2
  
  #Starting point
  undiag_start=c(1308090,1875356)/1000
  diag_start=c(355630,795094)/1000
  treat_start=c(32928,55927)/1000
  
  start_value=array(0,dim=c(24,4,2,2))
  start_value[1,1:4,1:2,1:2] = array(rep(repart_f(k1),2*2),dim=c(4,2,2)) *
    array(rep(rep(undiag_start,each=4),2),dim=c(4,2,2)) *
    array(rep(c(0.99, 0.01),each=4*2),dim=c(4,2,2))
  start_value[3,1:4,1:2,1:2] = array(rep(repart_f(k2),2*2),dim=c(4,2,2)) *
    array(rep(rep(diag_start,each=4),2),dim=c(4,2,2)) *
    array(rep(c(0.99, 0.01),each=4*2),dim=c(4,2,2))
  start_value[4,1:4,1:2,1:2] = array(rep(repart_f(k3),2*2),dim=c(4,2,2)) *
    array(rep(rep(treat_start,each=4),2),dim=c(4,2,2)) *
    array(rep(c(0.99, 0.01),each=4*2),dim=c(4,2,2))
  
  start_value=c(c(0,0),as.vector(start_value))
  return(start_value)
}

x_start2019=function(x,p_TDF,p_DTG){
  x <- array(x[-(1:2)], dim=c(24,4,2,2))
  #Diagnostic
  diag_m <- array(x[2,1:4,1,1:2]) + array(x[3,1:4,1,1:2])
  diag_w <- array(x[2,1:4,2,1:2]) + array(x[3,1:4,2,1:2])
  #NNRTI according to cd4, sex, and resistance
  nnrti_m <- array(x[4:6,1:4,1,1:2]) + array(x[7:9,1:4,1,1:2]) + array(x[10:12,1:4,1,1:2]) + array(x[13:15,1:4,1,1:2])
  nnrti_w <- array(x[4:6,1:4,2,1:2]) + array(x[7:9,1:4,2,1:2]) + array(x[10:12,1:4,2,1:2]) + array(x[13:15,1:4,2,1:2])
  
  #diagnostic
  x[2,1:4,1,1:2] = diag_m
  x[3,1:4,1,1:2] = 0
  x[2,1:4,2,1:2] = diag_w * p_DTG
  x[3,1:4,2,1:2] = diag_w * (1-p_DTG)
  
  #men: all eligible for DTG
  x[4:6,1:4,1,1:2] = 0
  x[7:9,1:4,1,1:2] = nnrti_m * p_TDF
  x[10:12,1:4,1,1:2] = 0
  x[13:15,1:4,1,1:2] = nnrti_m * (1-p_TDF)
  x[16:18,1:4,1,1:2] = 0 #0, as no DTG in 2019
  x[19:21,1:4,1,1:2] = 0 #0, as no DTG in 2019
  
  x[4:6,1:4,2,1:2] = nnrti_w * p_TDF * (1-p_DTG)
  x[7:9,1:4,2,1:2] = nnrti_w * p_TDF * p_DTG
  x[10:12,1:4,2,1:2] = nnrti_w * (1-p_TDF) * (1-p_DTG)
  x[13:15,1:4,2,1:2] = nnrti_w * (1-p_TDF) * p_DTG
  x[16:18,1:4,2,1:2] = 0 #0, as no DTG in 2019
  x[19:21,1:4,2,1:2] = 0 #0, as no DTG in 2019
  
  x=c(c(0,0),as.numeric(x))
  return(x)
}

SIR_new2=function( t , # time
                  y , # system state { susceptible , infected , recovered }
                  theta , # parameters { transmission rate , recovery rate }
                  x_r , # real valued fixed data
                  x_i ) {
  dtg_1st=theta[1];
  dtg_switch=theta[2];
  
  inf1  = c( 0.008*0.05, 0.003*0.95, 0.003*0.95, 0.0);
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
  t_cd4_treat_all = c(1.0/0.4, 1.0/0.5, 1.0/0.7, 1.0);
  t_cd4_elig_year1 = c(2015.5, 2013.5, 2009.5, 2001.0);
  t_cd4_elig_year2 = c(2016.5, 2015.5, 2012.5, 2004.0);
  t_1st =c();
  for(i in 1:4){
    t_1st[i] = rate_treat *smooth_f(2005.0 + t/12.0, 2001.0, 2012.0, 1.0, 200.0/12.0) *  #Increase in treatment rate between 2001 and 2012
    lin_f(2005.0 + t/12.0, t_cd4_elig_year1[i],  t_cd4_elig_year2[i], 0.0, t_cd4_elig[i]) *  #Change in eligiblity criteria
    lin_f(2005.0 + t/12.0, 2017.0, 2022.0, 1.0, t_cd4_treat_all[i]);  #Increase in treatment rate due to Treat-All policy
  }
  t_switch = c(4.35/531.0 /5.0, 4.35/427.0 /5.0, 4.35/294.0 /5.0, 4.35/189.0 /5.0);#Switching rate to PI, failing individuals (either DTG-ineligible on NNRTI or DTG-eligible on DTG)
  t_switch_elig =c(); #Switch rate to PI, DTG-eligible
  t_1st_NNRTI_inel=c();
  t_1st_NNRTI_elig=c();
  t_1st_DTG=c();
  t_switch_DTG_S=c();
  t_switch_DTG_F=c();
  for(i in 1:4){
    t_1st_DTG [i]= t_1st[i] * smooth_f(2005.0 + t/12.0, 2019.0, 2019.0, 0.0, 1.0) * dtg_1st; #Treatment initiation rate, DTG
    t_1st_NNRTI_inel[i] = t_1st[i];
    t_1st_NNRTI_elig[i] = t_1st[i] - t_1st_DTG[i];
    t_switch_elig[i] = t_switch[i] * smooth_f(2005.0 + t/12.0, 2019.0, 2019.0, 1.0, 1.0 - dtg_1st); #switch always 0 after 2019, except if dtg_1st=0 (no DTG introduction) in which case it remains at 1
    t_switch_DTG_S[i] = 1.0/12.0 * dtg_switch * smooth_f(2005.0 + t/12.0, 2019.0, 2019.0, 0.0, 1.0);
    if(dtg_1st==0.0) t_switch_DTG_F[i] = 0.0;
    if(dtg_1st==1.0 & dtg_switch==0.0) t_switch_DTG_F[i] = t_switch[i];
    if(dtg_switch==1.0) t_switch_DTG_F[i] = 1.0/12.0;
      t_switch_DTG_F[i] = t_switch_DTG_F[i] * smooth_f(2005.0 + t/12.0, 2019.0, 2019.0, 0.0, 1.0);
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
  #CD4 progression and ART efficacy by ART, 1: cd4 progression/efficacy equivalent to NNRTI, 2: equivalent to PI
  art_eq  = c(1, 1, 1, 1, 1, 1, 2);
  art_j =c()
  #Increase/decrease of efficacy due to NNRTI resistance by ART
  art_res1 = matrix(c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.64, 2.64, 2.64, 2.64, 1.0, 1.0, 1.0),nrow=2,byrow=T);
  art_res2 = matrix(c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 4.9, 4.9, 4.9, 4.9, 1.0, 1.0, 1.0),nrow=2,byrow=T);
  #Increase/decrease of efficacy relative to NNRTI (for NNRTI and DTG) and to PI (for PI)
  art_eff = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
  art_eff[5] = theta[4];
  art_eff[6] = theta[4];
  #Resistance
  res = 1/5.0;
  rev = 1/125.0;
  
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
  y[1]=N[1]-sum(y[posa(all_pos,1:4,1,1:2)]);
  y[2]=N[2]-sum(y[posa(all_pos,1:4,2,1:2)]);
  
  
  #Differential equations
  dy_dt=rep(0,24*16+2)
  #Susceptible, dydt=0 as values over time fixed by the spline function 
  dy_dt[1]=0;
  dy_dt[2]=0;
  
  for(l in 1:2){
    dy_dt[pos(1,1,1,l)]=dy_dt[pos(1,1,1,l)] +
      inf1[1] * inf2 * (sum(y[posa(1,1:4,1,l)]) + inf3 * sum(y[posa(inf_pos,1:4,1,l)])) * y[1] / (y[1] + sum(y[posa(all_pos,1:4,1,1:2)])) +
      inf1[2] * inf2 * (sum(y[posa(1,1:4,2,l)]) + inf3 * sum(y[posa(inf_pos,1:4,2,l)])) * y[1] / (y[1] + sum(y[posa(all_pos,1:4,1,1:2)]));
    dy_dt[pos(1,1,2,l)]=dy_dt[pos(1,1,2,l)] +
      inf1[3] * inf2 * (sum(y[posa(1,1:4,1,l)]) + inf3 * sum(y[posa(inf_pos,1:4,1,l)])) * y[2] / (y[2] + sum(y[posa(all_pos,1:4,2,1:2)])) +
      inf1[4] * inf2 * (sum(y[posa(1,1:4,2,l)]) + inf3 * sum(y[posa(inf_pos,1:4,2,l)])) * y[2] / (y[2] + sum(y[posa(all_pos,1:4,2,1:2)]));
  }
  
  for(j in 1:4){
    for(l in 1:2){
      #Diagnosis
      #men
      dy_dt[pos(1,j,1,l)] = dy_dt[pos(1,j,1,l)] - (diag2 / diag1 + oi_inc[j] * oi_test) * y[pos(1,j,1,l)];
      #women
      dy_dt[pos(1,j,2,l)] = dy_dt[pos(1,j,2,l)] - (diag2 / diag1 * diag_women + oi_inc[j] * oi_test + preg_inc[j] * preg_test) * y[pos(1,j,2,l)];
      
      #DTG-eligible, men
      dy_dt[pos(2,j,1,l)] = dy_dt[pos(2,j,1,l)] + p_dtg[1] * (diag2 / diag1 + oi_inc[j] * oi_test) * y[pos(1,j,1,l)];
      #DTG-ineligible, men
      dy_dt[pos(3,j,1,l)] = dy_dt[pos(3,j,1,l)] + (1 - p_dtg[1]) * (diag2 / diag1 + oi_inc[j] * oi_test) * y[pos(1,j,1,l)];
      #DTG-eligible, women
      dy_dt[pos(2,j,2,l)] = dy_dt[pos(2,j,2,l)] + p_dtg[2] * (diag2 / diag1 * diag_women + oi_inc[j] * oi_test + preg_inc[j] * preg_test) * y[pos(1,j,2,l)];
      #DTG-ineligible, women
      dy_dt[pos(3,j,2,l)] = dy_dt[pos(3,j,2,l)] + (1 - p_dtg[2]) * (diag2 / diag1 * diag_women + oi_inc[j] * oi_test + preg_inc[j] * preg_test) * y[pos(1,j,2,l)];
      
      
      for(k in 1:2){
        #Treatment initiation
        #NNRTI: 1) NNRTI+TDF for DTG-inel, 2) NNRTI+TDF for DTG-elig, 3) NNRTI+AZT for DTG-inel, 4) NNRTI+AZT for DTG-elig
        dy_dt[pos(4,j,k,l)] = dy_dt[pos(4,j,k,l)] + p_tdf * t_1st_NNRTI_inel[j] * y[pos(3,j,k,l)];
        dy_dt[pos(7,j,k,l)] = dy_dt[pos(7,j,k,l)] + p_tdf * t_1st_NNRTI_elig[j] * y[pos(2,j,k,l)];
        dy_dt[pos(10,j,k,l)] = dy_dt[pos(10,j,k,l)] + (1-p_tdf) * t_1st_NNRTI_inel[j] * y[pos(3,j,k,l)];
        dy_dt[pos(13,j,k,l)] = dy_dt[pos(13,j,k,l)] + (1-p_tdf) * t_1st_NNRTI_elig[j] * y[pos(2,j,k,l)];
        #DTG: DTG+TDF
        dy_dt[pos(16,j,k,l)] = dy_dt[pos(16,j,k,l)] + t_1st_DTG[j] * y[pos(2,j,k,l)];
        
        dy_dt[pos(3,j,k,l)] = dy_dt[pos(3,j,k,l)] - t_1st_NNRTI_inel[j] * y[pos(3,j,k,l)];
        dy_dt[pos(2,j,k,l)] = dy_dt[pos(2,j,k,l)] - (t_1st_NNRTI_elig[j] + t_1st_DTG[j]) * y[pos(2,j,k,l)];
        
        #Treatment switch to PI after failure
        #Switch from NNRTI to PI, DTG-ineligible
        dy_dt[pos(22,j,k,l)] = dy_dt[pos(22,j,k,l)] + t_switch[j] * (y[pos(6,j,k,l)] + y[pos(12,j,k,l)]);
        dy_dt[pos(6,j,k,l)] = dy_dt[pos(6,j,k,l)] - t_switch[j] * y[pos(6,j,k,l)];
        dy_dt[pos(12,j,k,l)] = dy_dt[pos(12,j,k,l)] - t_switch[j] * y[pos(12,j,k,l)];
        #Switch from DTG to PI
        dy_dt[pos(22,j,k,l)] = dy_dt[pos(22,j,k,l)] + t_switch[j] * (y[pos(18,j,k,l)] + y[pos(21,j,k,l)]);
        dy_dt[pos(18,j,k,l)] = dy_dt[pos(18,j,k,l)] - t_switch[j] * y[pos(18,j,k,l)];
        dy_dt[pos(21,j,k,l)] = dy_dt[pos(21,j,k,l)] - t_switch[j] * y[pos(21,j,k,l)];
        #Switch from NNRTI to PI, DTG-eligible
        dy_dt[pos(22,j,k,l)] = dy_dt[pos(22,j,k,l)] + t_switch_elig[j] * (y[pos(9,j,k,l)] + y[pos(15,j,k,l)]);
        dy_dt[pos(9,j,k,l)] = dy_dt[pos(9,j,k,l)] - t_switch_elig[j] * y[pos(9,j,k,l)];
        dy_dt[pos(15,j,k,l)] = dy_dt[pos(15,j,k,l)] - t_switch_elig[j] * y[pos(15,j,k,l)];
        
        #Treatment switch to DTG, suppressed
        dy_dt[pos(17,j,k,l)] = dy_dt[pos(17,j,k,l)] + t_switch_DTG_S[j] * (y[pos(8,j,k,l)] + y[pos(14,j,k,l)]);
        dy_dt[pos(8,j,k,l)] = dy_dt[pos(8,j,k,l)] - t_switch_DTG_S[j] * y[pos(8,j,k,l)];
        dy_dt[pos(14,j,k,l)] = dy_dt[pos(14,j,k,l)] - t_switch_DTG_S[j] * y[pos(14,j,k,l)];
        
        #Treatment switch to DTG, failed
        dy_dt[pos(16,j,k,l)] = dy_dt[pos(16,j,k,l)] + t_switch_DTG_F[j] * y[pos(15,j,k,l)];
        dy_dt[pos(19,j,k,l)] = dy_dt[pos(19,j,k,l)] + t_switch_DTG_F[j] * y[pos(9,j,k,l)];
        
        dy_dt[pos(15,j,k,l)] = dy_dt[pos(15,j,k,l)] - t_switch_DTG_F[j] * y[pos(15,j,k,l)];
        dy_dt[pos(9,j,k,l)] = dy_dt[pos(9,j,k,l)] - t_switch_DTG_F[j] * y[pos(9,j,k,l)];
      }
    }
  }
  
  #Treatment stages and CD4, mortality
  #Infected and diagnosed
  for(k in 1:2){
    for(l in 1:2){
      for(i in 1:3){
        dy_dt[pos(i,1,k,l)] = dy_dt[pos(i,1,k,l)] - cd4_untreated[1] * y[pos(i,1,k,l)];
        dy_dt[pos(i,2,k,l)] = dy_dt[pos(i,2,k,l)] + cd4_untreated[1] * y[pos(i,1,k,l)] - cd4_untreated[2] * y[pos(i,2,k,l)];
        dy_dt[pos(i,3,k,l)] = dy_dt[pos(i,3,k,l)] + cd4_untreated[2] * y[pos(i,2,k,l)] - cd4_untreated[3] * y[pos(i,3,k,l)];
        dy_dt[pos(i,4,k,l)] = dy_dt[pos(i,4,k,l)] + cd4_untreated[3] * y[pos(i,3,k,l)];
        
        #Mortality
        for(j in 1:4){
          dy_dt[pos(i,j,k,l)] = dy_dt[pos(i,j,k,l)] - mu[1,j] * y[pos(i,j,k,l)];
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
      for(i in 1:4){
        treat_j[i] =  treat[art_j,i] * art_eff[art] / art_res1[l,art];
      }
      for(i in 5:8){
        treat_j[i] =  treat[art_j,i] / art_eff[art] * art_res1[l,art];
      }
      for(i in 9:12){
        treat_j[i] =  treat[art_j,i] * art_eff[art] * art_res2[l,art];
      }
      for(i in 13:16){
        treat_j[i] =  treat[art_j,i] / art_eff[art] / art_res2[l,art];
      }
      
      for(k in 1:2){
        #Start
        dy_dt[pos(3*art+1,1,k,l)] = dy_dt[pos(3*art+1,1,k,l)] + cd4_j[4] * y[pos(3*art+1,2,k,l)] - cd4_j[1] * y[pos(3*art+1,1,k,l)]+
          -(treat_j[1] + treat_j[5]) * y[pos(3*art+1,1,k,l)];
        dy_dt[pos(3*art+1,2,k,l)] = dy_dt[pos(3*art+1,2,k,l)] + cd4_j[1] * y[pos(3*art+1,1,k,l)] + cd4_j[5] * y[pos(3*art+1,3,k,l)] - (cd4_j[4] + cd4_j[2]) * y[pos(3*art+1,2,k,l)]+
          -(treat_j[2] + treat_j[6]) * y[pos(3*art+1,2,k,l)];
        dy_dt[pos(3*art+1,3,k,l)] = dy_dt[pos(3*art+1,3,k,l)] + cd4_j[2] * y[pos(3*art+1,2,k,l)] + cd4_j[6] * y[pos(3*art+1,4,k,l)] - (cd4_j[5] + cd4_j[3]) * y[pos(3*art+1,3,k,l)]+
          -(treat_j[3] + treat_j[7]) * y[pos(3*art+1,3,k,l)];
        dy_dt[pos(3*art+1,4,k,l)] = dy_dt[pos(3*art+1,4,k,l)] + cd4_j[3] * y[pos(3*art+1,3,k,l)] - cd4_j[6] * y[pos(3*art+1,4,k,l)]+
          -(treat_j[4] + treat_j[8]) * y[pos(3*art+1,4,k,l)];
        #Suppressed
        dy_dt[pos(3*art+2,1,k,l)] = dy_dt[pos(3*art+2,1,k,l)] + cd4_j[7] * y[pos(3*art+2,2,k,l)]+
          + treat_j[1] * y[pos(3*art+1,1,k,l)] + treat_j[13] * y[pos(3*art+3,1,k,l)] - treat_j[9] * y[pos(3*art+2,1,k,l)];
        dy_dt[pos(3*art+2,2,k,l)] = dy_dt[pos(3*art+2,2,k,l)] + cd4_j[8] * y[pos(3*art+2,3,k,l)] - cd4_j[7] * y[pos(3*art+2,2,k,l)]+
          + treat_j[2] * y[pos(3*art+1,2,k,l)] + treat_j[14] * y[pos(3*art+3,2,k,l)] - treat_j[10] * y[pos(3*art+2,2,k,l)];
        dy_dt[pos(3*art+2,3,k,l)] = dy_dt[pos(3*art+2,3,k,l)] + cd4_j[9] * y[pos(3*art+2,4,k,l)] - cd4_j[8] * y[pos(3*art+2,3,k,l)]+
          + treat_j[3] * y[pos(3*art+1,3,k,l)] + treat_j[15] * y[pos(3*art+3,3,k,l)] - treat_j[11] * y[pos(3*art+2,3,k,l)];
        dy_dt[pos(3*art+2,4,k,l)] = dy_dt[pos(3*art+2,4,k,l)] - cd4_j[9] * y[pos(3*art+2,4,k,l)]+
          + treat_j[4] * y[pos(3*art+1,4,k,l)] + treat_j[16] * y[pos(3*art+3,4,k,l)] - treat_j[12] * y[pos(3*art+2,4,k,l)];
        #Failed
        dy_dt[pos(3*art+3,1,k,l)] = dy_dt[pos(3*art+3,1,k,l)] - cd4_j[10] * y[pos(3*art+3,1,k,l)]+
          treat_j[5] * y[pos(3*art+1,1,k,l)] + treat_j[9] * y[pos(3*art+2,1,k,l)] - treat_j[13] * y[pos(3*art+3,1,k,l)];
        dy_dt[pos(3*art+3,2,k,l)] = dy_dt[pos(3*art+3,2,k,l)] + cd4_j[10] * y[pos(3*art+3,1,k,l)] - cd4_j[11] * y[pos(3*art+3,2,k,l)]+
          treat_j[6] * y[pos(3*art+1,2,k,l)] + treat_j[10] * y[pos(3*art+2,2,k,l)] - treat_j[14] * y[pos(3*art+3,2,k,l)];
        dy_dt[pos(3*art+3,3,k,l)] = dy_dt[pos(3*art+3,3,k,l)] + cd4_j[11] * y[pos(3*art+3,2,k,l)] - cd4_j[12] * y[pos(3*art+3,3,k,l)]+
          treat_j[7] * y[pos(3*art+1,3,k,l)] + treat_j[11] * y[pos(3*art+2,3,k,l)] - treat_j[15] * y[pos(3*art+3,3,k,l)];
        dy_dt[pos(3*art+3,4,k,l)] = dy_dt[pos(3*art+3,4,k,l)] + cd4_j[12] * y[pos(3*art+3,3,k,l)]+
          treat_j[8] * y[pos(3*art+1,4,k,l)] + treat_j[12] * y[pos(3*art+2,4,k,l)] - treat_j[16] * y[pos(3*art+3,4,k,l)];
        
        #Mortality 
        for(j in 1:4){
          dy_dt[pos(3*art+1,j,k,l)] = dy_dt[pos(3*art+1,j,k,l)] - mu[2,j] * y[pos(3*art+1,j,k,l)];
          dy_dt[pos(3*art+2,j,k,l)] = dy_dt[pos(3*art+2,j,k,l)] - mu[3,j] * y[pos(3*art+2,j,k,l)];
          dy_dt[pos(3*art+3,j,k,l)] = dy_dt[pos(3*art+3,j,k,l)] - mu[4,j] * y[pos(3*art+3,j,k,l)];
        }
      }
    }
  }
  
  for(j in 1:4){
    for(k in 1:2){
      dy_dt[pos(6,j,k,2)] = dy_dt[pos(6,j,k,2)] + res * y[pos(6,j,k,1)];
      dy_dt[pos(9,j,k,2)] = dy_dt[pos(9,j,k,2)] + res * y[pos(9,j,k,1)];
      dy_dt[pos(12,j,k,2)] = dy_dt[pos(12,j,k,2)] + res * y[pos(12,j,k,1)];
      dy_dt[pos(15,j,k,2)] = dy_dt[pos(15,j,k,2)] + res * y[pos(15,j,k,1)];
      
      dy_dt[pos(6,j,k,1)] = dy_dt[pos(6,j,k,1)] - res * y[pos(6,j,k,1)];
      dy_dt[pos(9,j,k,1)] = dy_dt[pos(9,j,k,1)] - res * y[pos(9,j,k,1)];
      dy_dt[pos(12,j,k,1)] = dy_dt[pos(12,j,k,1)] - res * y[pos(12,j,k,1)];
      dy_dt[pos(15,j,k,1)] = dy_dt[pos(15,j,k,1)] - res * y[pos(15,j,k,1)];
      
      #resistance reversion
      dy_dt[pos(1,j,k,1)] = dy_dt[pos(1,j,k,1)] + rev * y[pos(1,j,k,2)];
      dy_dt[pos(2,j,k,1)] = dy_dt[pos(2,j,k,1)] + rev * y[pos(2,j,k,2)];
      dy_dt[pos(3,j,k,1)] = dy_dt[pos(3,j,k,1)] + rev * y[pos(3,j,k,2)];
      
      dy_dt[pos(1,j,k,2)] = dy_dt[pos(1,j,k,2)] - rev * y[pos(1,j,k,2)];
      dy_dt[pos(2,j,k,2)] = dy_dt[pos(2,j,k,2)] - rev * y[pos(2,j,k,2)];
      dy_dt[pos(3,j,k,2)] = dy_dt[pos(3,j,k,2)] - rev * y[pos(3,j,k,2)];
    }
  }
  return(list(dy_dt))
}
#different theta
SIR_new3=function( t , # time
                   y , # system state { susceptible , infected , recovered }
                   theta , # parameters { transmission rate , recovery rate }
                   x_r , # real valued fixed data
                   x_i ) {
  dtg_1st=theta[1];
  dtg_switch=theta[2];
  dtg_eff=theta[4];
  res=1.0/theta[5];
  rev=1.0/theta[6];
  alpha1=theta[7];
  alpha2=theta[8];
  p_msm=theta[9];
  risk_r=theta[10];
  hiv_msm=theta[11];
  var_switch_DTG_S=theta[12]; #Increase or decrease of switching rate (used in the heatmap)
  var_switch_DTG_F=theta[13]; #Increase or decrease of switching rate (used in the heatmap)
  
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
  t_cd4_treat_all = c(1.0/0.4, 1.0/0.5, 1.0/0.7, 1.0);
  t_cd4_elig_year1 = c(2015.5, 2013.5, 2009.5, 2001.0);
  t_cd4_elig_year2 = c(2016.5, 2015.5, 2012.5, 2004.0);
  t_1st =c();
  for(i in 1:4){
    t_1st[i] = rate_treat *smooth_f(2005.0 + t/12.0, 2001.0, 2012.0, 1.0, 200.0/12.0) *  #Increase in treatment rate between 2001 and 2012
      lin_f(2005.0 + t/12.0, t_cd4_elig_year1[i],  t_cd4_elig_year2[i], 0.0, t_cd4_elig[i]) *  #Change in eligiblity criteria
      lin_f(2005.0 + t/12.0, 2017.0, 2022.0, 1.0, t_cd4_treat_all[i]);  #Increase in treatment rate due to Treat-All policy
  }
  t_switch = c(4.35/531.0 /5.0, 4.35/427.0 /5.0, 4.35/294.0 /5.0, 4.35/189.0 /5.0);#Switching rate to PI, failing individuals (either DTG-ineligible on NNRTI or DTG-eligible on DTG)
  t_switch_elig =c(); #Switch rate to PI, DTG-eligible
  t_1st_NNRTI_inel=c();
  t_1st_NNRTI_elig=c();
  t_1st_DTG=c();
  t_switch_DTG_S=c();
  t_switch_DTG_F=c();
  for(i in 1:4){
    t_1st_DTG [i]= t_1st[i] * smooth_f(2005.0 + t/12.0, 2019.0, 2019.0, 0.0, 1.0) * dtg_1st; #Treatment initiation rate, DTG
    t_1st_NNRTI_inel[i] = t_1st[i];
    t_1st_NNRTI_elig[i] = t_1st[i] - t_1st_DTG[i];
    t_switch_elig[i] = t_switch[i] * smooth_f(2005.0 + t/12.0, 2019.0, 2019.0, 1.0, 1.0 - dtg_1st); #switch always 0 after 2019, except if dtg_1st=0 (no DTG introduction) in which case it remains at 1
    t_switch_DTG_S[i] = 1.0/12.0 * dtg_switch * var_switch_DTG_S * smooth_f(2005.0 + t/12.0, 2019.0, 2019.0, 0.0, 1.0);
    if(dtg_1st==0.0) t_switch_DTG_F[i] = 0.0;
    if(dtg_1st==1.0 & dtg_switch==0.0) t_switch_DTG_F[i] = t_switch[i];
    if(dtg_switch==1.0) t_switch_DTG_F[i] = 1.0/12.0 * var_switch_DTG_F;
    t_switch_DTG_F[i] = t_switch_DTG_F[i] * smooth_f(2005.0 + t/12.0, 2019.0, 2019.0, 0.0, 1.0);
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
  #CD4 progression and ART efficacy by ART, 1: cd4 progression/efficacy equivalent to NNRTI, 2: equivalent to PI
  art_eq  = c(1, 1, 1, 1, 1, 1, 2);
  art_j =c()
  #Increase/decrease of efficacy due to NNRTI resistance by ART
  art_res1 = matrix(c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, alpha1, alpha1, alpha1, alpha1, 1.0, 1.0, 1.0),nrow=2,byrow=T);
  art_res2 = matrix(c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, alpha2, alpha2, alpha2, alpha2, 1.0, 1.0, 1.0),nrow=2,byrow=T);
  #Increase/decrease of efficacy relative to NNRTI (for NNRTI and DTG) and to PI (for PI)
  art_eff = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0);
  art_eff[5] = dtg_eff;
  art_eff[6] = dtg_eff;

  
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
  y[1]=N[1]-sum(y[posa(all_pos,1:4,1,1:2)]);
  y[2]=N[2]-sum(y[posa(all_pos,1:4,2,1:2)]);
  
  
  #Differential equations
  dy_dt=rep(0,24*16+2)
  #Susceptible, dydt=0 as values over time fixed by the spline function 
  dy_dt[1]=0;
  dy_dt[2]=0;
  
  for(l in 1:2){
    dy_dt[pos(1,1,1,l)]=dy_dt[pos(1,1,1,l)] +
      inf1[1] * inf2 * (sum(y[posa(1,1:4,1,l)]) + inf3 * sum(y[posa(inf_pos,1:4,1,l)])) * y[1] / (y[1] + sum(y[posa(all_pos,1:4,1,1:2)])) +
      inf1[2] * inf2 * (sum(y[posa(1,1:4,2,l)]) + inf3 * sum(y[posa(inf_pos,1:4,2,l)])) * y[1] / (y[1] + sum(y[posa(all_pos,1:4,1,1:2)]));
    dy_dt[pos(1,1,2,l)]=dy_dt[pos(1,1,2,l)] +
      inf1[3] * inf2 * (sum(y[posa(1,1:4,1,l)]) + inf3 * sum(y[posa(inf_pos,1:4,1,l)])) * y[2] / (y[2] + sum(y[posa(all_pos,1:4,2,1:2)])) +
      inf1[4] * inf2 * (sum(y[posa(1,1:4,2,l)]) + inf3 * sum(y[posa(inf_pos,1:4,2,l)])) * y[2] / (y[2] + sum(y[posa(all_pos,1:4,2,1:2)]));
    }
  
  
  for(j in 1:4){
    for(l in 1:2){
      #Diagnosis
      #men
      dy_dt[pos(1,j,1,l)] = dy_dt[pos(1,j,1,l)] - (diag2 / diag1 + oi_inc[j] * oi_test) * y[pos(1,j,1,l)];
      #women
      dy_dt[pos(1,j,2,l)] = dy_dt[pos(1,j,2,l)] - (diag2 / diag1 * diag_women + oi_inc[j] * oi_test + preg_inc[j] * preg_test) * y[pos(1,j,2,l)];
      
      #DTG-eligible, men
      dy_dt[pos(2,j,1,l)] = dy_dt[pos(2,j,1,l)] + p_dtg[1] * (diag2 / diag1 + oi_inc[j] * oi_test) * y[pos(1,j,1,l)];
      #DTG-ineligible, men
      dy_dt[pos(3,j,1,l)] = dy_dt[pos(3,j,1,l)] + (1 - p_dtg[1]) * (diag2 / diag1 + oi_inc[j] * oi_test) * y[pos(1,j,1,l)];
      #DTG-eligible, women
      dy_dt[pos(2,j,2,l)] = dy_dt[pos(2,j,2,l)] + p_dtg[2] * (diag2 / diag1 * diag_women + oi_inc[j] * oi_test + preg_inc[j] * preg_test) * y[pos(1,j,2,l)];
      #DTG-ineligible, women
      dy_dt[pos(3,j,2,l)] = dy_dt[pos(3,j,2,l)] + (1 - p_dtg[2]) * (diag2 / diag1 * diag_women + oi_inc[j] * oi_test + preg_inc[j] * preg_test) * y[pos(1,j,2,l)];
      
      
      for(k in 1:2){
        #Treatment initiation
        #NNRTI: 1) NNRTI+TDF for DTG-inel, 2) NNRTI+TDF for DTG-elig, 3) NNRTI+AZT for DTG-inel, 4) NNRTI+AZT for DTG-elig
        dy_dt[pos(4,j,k,l)] = dy_dt[pos(4,j,k,l)] + p_tdf * t_1st_NNRTI_inel[j] * y[pos(3,j,k,l)];
        dy_dt[pos(7,j,k,l)] = dy_dt[pos(7,j,k,l)] + p_tdf * t_1st_NNRTI_elig[j] * y[pos(2,j,k,l)];
        dy_dt[pos(10,j,k,l)] = dy_dt[pos(10,j,k,l)] + (1-p_tdf) * t_1st_NNRTI_inel[j] * y[pos(3,j,k,l)];
        dy_dt[pos(13,j,k,l)] = dy_dt[pos(13,j,k,l)] + (1-p_tdf) * t_1st_NNRTI_elig[j] * y[pos(2,j,k,l)];
        #DTG: DTG+TDF
        dy_dt[pos(16,j,k,l)] = dy_dt[pos(16,j,k,l)] + t_1st_DTG[j] * y[pos(2,j,k,l)];
        
        dy_dt[pos(3,j,k,l)] = dy_dt[pos(3,j,k,l)] - t_1st_NNRTI_inel[j] * y[pos(3,j,k,l)];
        dy_dt[pos(2,j,k,l)] = dy_dt[pos(2,j,k,l)] - (t_1st_NNRTI_elig[j] + t_1st_DTG[j]) * y[pos(2,j,k,l)];
        
        #Treatment switch to PI after failure
        #Switch from NNRTI to PI, DTG-ineligible
        dy_dt[pos(22,j,k,l)] = dy_dt[pos(22,j,k,l)] + t_switch[j] * (y[pos(6,j,k,l)] + y[pos(12,j,k,l)]);
        dy_dt[pos(6,j,k,l)] = dy_dt[pos(6,j,k,l)] - t_switch[j] * y[pos(6,j,k,l)];
        dy_dt[pos(12,j,k,l)] = dy_dt[pos(12,j,k,l)] - t_switch[j] * y[pos(12,j,k,l)];
        #Switch from DTG to PI
        dy_dt[pos(22,j,k,l)] = dy_dt[pos(22,j,k,l)] + t_switch[j] * (y[pos(18,j,k,l)] + y[pos(21,j,k,l)]);
        dy_dt[pos(18,j,k,l)] = dy_dt[pos(18,j,k,l)] - t_switch[j] * y[pos(18,j,k,l)];
        dy_dt[pos(21,j,k,l)] = dy_dt[pos(21,j,k,l)] - t_switch[j] * y[pos(21,j,k,l)];
        #Switch from NNRTI to PI, DTG-eligible
        dy_dt[pos(22,j,k,l)] = dy_dt[pos(22,j,k,l)] + t_switch_elig[j] * (y[pos(9,j,k,l)] + y[pos(15,j,k,l)]);
        dy_dt[pos(9,j,k,l)] = dy_dt[pos(9,j,k,l)] - t_switch_elig[j] * y[pos(9,j,k,l)];
        dy_dt[pos(15,j,k,l)] = dy_dt[pos(15,j,k,l)] - t_switch_elig[j] * y[pos(15,j,k,l)];
        
        #Treatment switch to DTG, suppressed
        dy_dt[pos(17,j,k,l)] = dy_dt[pos(17,j,k,l)] + t_switch_DTG_S[j] * (y[pos(8,j,k,l)] + y[pos(14,j,k,l)]);
        dy_dt[pos(8,j,k,l)] = dy_dt[pos(8,j,k,l)] - t_switch_DTG_S[j] * y[pos(8,j,k,l)];
        dy_dt[pos(14,j,k,l)] = dy_dt[pos(14,j,k,l)] - t_switch_DTG_S[j] * y[pos(14,j,k,l)];
        
        #Treatment switch to DTG, failed
        dy_dt[pos(16,j,k,l)] = dy_dt[pos(16,j,k,l)] + t_switch_DTG_F[j] * y[pos(15,j,k,l)];
        dy_dt[pos(19,j,k,l)] = dy_dt[pos(19,j,k,l)] + t_switch_DTG_F[j] * y[pos(9,j,k,l)];
        
        dy_dt[pos(15,j,k,l)] = dy_dt[pos(15,j,k,l)] - t_switch_DTG_F[j] * y[pos(15,j,k,l)];
        dy_dt[pos(9,j,k,l)] = dy_dt[pos(9,j,k,l)] - t_switch_DTG_F[j] * y[pos(9,j,k,l)];
      }
    }
  }
  
  #Treatment stages and CD4, mortality
  #Infected and diagnosed
  for(k in 1:2){
    for(l in 1:2){
      for(i in 1:3){
        dy_dt[pos(i,1,k,l)] = dy_dt[pos(i,1,k,l)] - cd4_untreated[1] * y[pos(i,1,k,l)];
        dy_dt[pos(i,2,k,l)] = dy_dt[pos(i,2,k,l)] + cd4_untreated[1] * y[pos(i,1,k,l)] - cd4_untreated[2] * y[pos(i,2,k,l)];
        dy_dt[pos(i,3,k,l)] = dy_dt[pos(i,3,k,l)] + cd4_untreated[2] * y[pos(i,2,k,l)] - cd4_untreated[3] * y[pos(i,3,k,l)];
        dy_dt[pos(i,4,k,l)] = dy_dt[pos(i,4,k,l)] + cd4_untreated[3] * y[pos(i,3,k,l)];
        
        #Mortality
        for(j in 1:4){
          dy_dt[pos(i,j,k,l)] = dy_dt[pos(i,j,k,l)] - mu[1,j] * y[pos(i,j,k,l)];
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
      for(i in 1:4){
        treat_j[i] =  treat[art_j,i] * art_eff[art] / art_res1[l,art];
      }
      for(i in 5:8){
        treat_j[i] =  treat[art_j,i] / art_eff[art] * art_res1[l,art];
      }
      for(i in 9:12){
        treat_j[i] =  treat[art_j,i] * art_eff[art] * art_res2[l,art];
      }
      for(i in 13:16){
        treat_j[i] =  treat[art_j,i] / art_eff[art] / art_res2[l,art];
      }
      
      for(k in 1:2){
        #Start
        dy_dt[pos(3*art+1,1,k,l)] = dy_dt[pos(3*art+1,1,k,l)] + cd4_j[4] * y[pos(3*art+1,2,k,l)] - cd4_j[1] * y[pos(3*art+1,1,k,l)]+
          -(treat_j[1] + treat_j[5]) * y[pos(3*art+1,1,k,l)];
        dy_dt[pos(3*art+1,2,k,l)] = dy_dt[pos(3*art+1,2,k,l)] + cd4_j[1] * y[pos(3*art+1,1,k,l)] + cd4_j[5] * y[pos(3*art+1,3,k,l)] - (cd4_j[4] + cd4_j[2]) * y[pos(3*art+1,2,k,l)]+
          -(treat_j[2] + treat_j[6]) * y[pos(3*art+1,2,k,l)];
        dy_dt[pos(3*art+1,3,k,l)] = dy_dt[pos(3*art+1,3,k,l)] + cd4_j[2] * y[pos(3*art+1,2,k,l)] + cd4_j[6] * y[pos(3*art+1,4,k,l)] - (cd4_j[5] + cd4_j[3]) * y[pos(3*art+1,3,k,l)]+
          -(treat_j[3] + treat_j[7]) * y[pos(3*art+1,3,k,l)];
        dy_dt[pos(3*art+1,4,k,l)] = dy_dt[pos(3*art+1,4,k,l)] + cd4_j[3] * y[pos(3*art+1,3,k,l)] - cd4_j[6] * y[pos(3*art+1,4,k,l)]+
          -(treat_j[4] + treat_j[8]) * y[pos(3*art+1,4,k,l)];
        #Suppressed
        dy_dt[pos(3*art+2,1,k,l)] = dy_dt[pos(3*art+2,1,k,l)] + cd4_j[7] * y[pos(3*art+2,2,k,l)]+
          + treat_j[1] * y[pos(3*art+1,1,k,l)] + treat_j[13] * y[pos(3*art+3,1,k,l)] - treat_j[9] * y[pos(3*art+2,1,k,l)];
        dy_dt[pos(3*art+2,2,k,l)] = dy_dt[pos(3*art+2,2,k,l)] + cd4_j[8] * y[pos(3*art+2,3,k,l)] - cd4_j[7] * y[pos(3*art+2,2,k,l)]+
          + treat_j[2] * y[pos(3*art+1,2,k,l)] + treat_j[14] * y[pos(3*art+3,2,k,l)] - treat_j[10] * y[pos(3*art+2,2,k,l)];
        dy_dt[pos(3*art+2,3,k,l)] = dy_dt[pos(3*art+2,3,k,l)] + cd4_j[9] * y[pos(3*art+2,4,k,l)] - cd4_j[8] * y[pos(3*art+2,3,k,l)]+
          + treat_j[3] * y[pos(3*art+1,3,k,l)] + treat_j[15] * y[pos(3*art+3,3,k,l)] - treat_j[11] * y[pos(3*art+2,3,k,l)];
        dy_dt[pos(3*art+2,4,k,l)] = dy_dt[pos(3*art+2,4,k,l)] - cd4_j[9] * y[pos(3*art+2,4,k,l)]+
          + treat_j[4] * y[pos(3*art+1,4,k,l)] + treat_j[16] * y[pos(3*art+3,4,k,l)] - treat_j[12] * y[pos(3*art+2,4,k,l)];
        #Failed
        dy_dt[pos(3*art+3,1,k,l)] = dy_dt[pos(3*art+3,1,k,l)] - cd4_j[10] * y[pos(3*art+3,1,k,l)]+
          treat_j[5] * y[pos(3*art+1,1,k,l)] + treat_j[9] * y[pos(3*art+2,1,k,l)] - treat_j[13] * y[pos(3*art+3,1,k,l)];
        dy_dt[pos(3*art+3,2,k,l)] = dy_dt[pos(3*art+3,2,k,l)] + cd4_j[10] * y[pos(3*art+3,1,k,l)] - cd4_j[11] * y[pos(3*art+3,2,k,l)]+
          treat_j[6] * y[pos(3*art+1,2,k,l)] + treat_j[10] * y[pos(3*art+2,2,k,l)] - treat_j[14] * y[pos(3*art+3,2,k,l)];
        dy_dt[pos(3*art+3,3,k,l)] = dy_dt[pos(3*art+3,3,k,l)] + cd4_j[11] * y[pos(3*art+3,2,k,l)] - cd4_j[12] * y[pos(3*art+3,3,k,l)]+
          treat_j[7] * y[pos(3*art+1,3,k,l)] + treat_j[11] * y[pos(3*art+2,3,k,l)] - treat_j[15] * y[pos(3*art+3,3,k,l)];
        dy_dt[pos(3*art+3,4,k,l)] = dy_dt[pos(3*art+3,4,k,l)] + cd4_j[12] * y[pos(3*art+3,3,k,l)]+
          treat_j[8] * y[pos(3*art+1,4,k,l)] + treat_j[12] * y[pos(3*art+2,4,k,l)] - treat_j[16] * y[pos(3*art+3,4,k,l)];
        
        #Mortality 
        for(j in 1:4){
          dy_dt[pos(3*art+1,j,k,l)] = dy_dt[pos(3*art+1,j,k,l)] - mu[2,j] * y[pos(3*art+1,j,k,l)];
          dy_dt[pos(3*art+2,j,k,l)] = dy_dt[pos(3*art+2,j,k,l)] - mu[3,j] * y[pos(3*art+2,j,k,l)];
          dy_dt[pos(3*art+3,j,k,l)] = dy_dt[pos(3*art+3,j,k,l)] - mu[4,j] * y[pos(3*art+3,j,k,l)];
        }
      }
    }
  }
  
  for(j in 1:4){
    for(k in 1:2){
      dy_dt[pos(6,j,k,2)] = dy_dt[pos(6,j,k,2)] + res * y[pos(6,j,k,1)];
      dy_dt[pos(9,j,k,2)] = dy_dt[pos(9,j,k,2)] + res * y[pos(9,j,k,1)];
      dy_dt[pos(12,j,k,2)] = dy_dt[pos(12,j,k,2)] + res * y[pos(12,j,k,1)];
      dy_dt[pos(15,j,k,2)] = dy_dt[pos(15,j,k,2)] + res * y[pos(15,j,k,1)];
      
      dy_dt[pos(6,j,k,1)] = dy_dt[pos(6,j,k,1)] - res * y[pos(6,j,k,1)];
      dy_dt[pos(9,j,k,1)] = dy_dt[pos(9,j,k,1)] - res * y[pos(9,j,k,1)];
      dy_dt[pos(12,j,k,1)] = dy_dt[pos(12,j,k,1)] - res * y[pos(12,j,k,1)];
      dy_dt[pos(15,j,k,1)] = dy_dt[pos(15,j,k,1)] - res * y[pos(15,j,k,1)];
      
      #resistance reversion
      dy_dt[pos(1,j,k,1)] = dy_dt[pos(1,j,k,1)] + rev * y[pos(1,j,k,2)];
      dy_dt[pos(2,j,k,1)] = dy_dt[pos(2,j,k,1)] + rev * y[pos(2,j,k,2)];
      dy_dt[pos(3,j,k,1)] = dy_dt[pos(3,j,k,1)] + rev * y[pos(3,j,k,2)];
      
      dy_dt[pos(1,j,k,2)] = dy_dt[pos(1,j,k,2)] - rev * y[pos(1,j,k,2)];
      dy_dt[pos(2,j,k,2)] = dy_dt[pos(2,j,k,2)] - rev * y[pos(2,j,k,2)];
      dy_dt[pos(3,j,k,2)] = dy_dt[pos(3,j,k,2)] - rev * y[pos(3,j,k,2)];
    }
  }
  return(list(dy_dt))
}
