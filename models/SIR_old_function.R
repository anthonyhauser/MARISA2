SIR_old=function( t , # time
                  y , # system state { susceptible , infected , recovered }
                  theta , # parameters { transmission rate , recovery rate }
                  x_r , # real valued fixed data
                  x_i ) { # integer valued fixed data
  t_1st=rep(0,4)  
  treat_j=rep(0,16)
  
  #Infection
  inf1=c(0.008*0.05,0.003*0.95,0.003*0.95,0)  #Probability of infection, 1 unprotected sexual intercourse, M->M, M->F, F->M, F->F
  inf2=3.318999  #number of unprotected intercourse per month
  inf3=0.5440886  #reduced number of unprotected intercourse for people knowing their HIV-positive status
  all_pos=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)  #position of all HIV stages
  inf_pos=c(2,3,4,6,7,9,10,12,13,15,16,18,19,21,22,24) #position of infectious HIV stages
  #Diagnosis
  diag1 = 273.6363  # number of months before diagnosis (without opportunistic disease - oi) in 2005: 22 years
  diag2 = 7.710578  # fold increase between 2005 and 2015: 3 years
  diag_women=1.25  #increase of diagnosis rate (without oi) for women
  oi_inc=c(0.05/12.0, 0.12/12.0, 0.27/12.0, 0.9/12.0)  #oi incidence by cd4
  oi_test=0.8  #proportion of oi test
  preg_inc=c(1.0 * 23.0/(12.0*1000.0),0.96 * 23.0/(12.0*1000.0),0.87 * 23.0/(12.0*1000.0),0.74 * 23.0/(12.0*1000.0))  #incidence of pregnancy by cd4
  preg_test=0.98  #proportion of pregnancy test
  #Treatment
  p_dtg = c(1,1)  #proportion of women opting for DTG
  p_tdf = 1  #proportion of people opting/being prescribed tdf (over azt)
  rate_treat = 0.001879695  #free treatment parameter rates fixed in project 1
  #t_cd4 = c(0.76, 0.8, 0.88, 1.0)  #relative initiation rates by cd4 in 2020 (increase to become identical from 2022, Treat-All policy)
  #for(v in 1:4) t_1st[v] = rate_treat * 200.0/12.0 * smooth(2020.0 + t/12.0, 2020.0, 2022.0, t_cd4[v], 1.0)  #Treatment initiation rate, no DTG
  t_cd4 = c(0.4,0.5,0.7,1)
  for(v in 1:4) t_1st[v] = rate_treat * 200.0/12.0 * lin_st(2020.0 + t/12.0, 2017.0, 2022.0, t_cd4[v], 1.0)
  t_1st_DTG = t_1st  #Treatment initiation rate, DTG
  t_switch_DTG_S = c(1.0/12.0, 1.0/12.0, 1.0/12.0, 1.0/12.0)  #Switching rate to DTG, suppressed individuals
  t_switch_DTG_F = c(1.0/12.0, 1.0/12.0, 1.0/12.0, 1.0/12.0)  #Switching rate to DTG, failing individuals
  t_switch = c(4.35/531.0 /5.0, 4.35/427.0 /5.0, 4.35/294.0 /5.0, 4.35/189.0 /5.0)  #Switching rate to PI, failing individuals (either DTG-ineligible on NNRTI or DTG-eligible on DTG)
  #Mortality
  rate_death=0.1632000  #Mortality rate per month per 1000 people (reference, supp people with cd4>500)
  p1=0.865181  #proportion of people with cd4>50 among people with cd4<200
  mu = matrix(c(1.57 *rate_death/1000.0, 2 *rate_death/1000.0, 4.57 *rate_death/1000.0, (p1*40.9+(1-p1)*134.4) *rate_death/1000.0, #relative mortality for untreated by cd4
                (0.9 * 1 + 0.1 * 3.92) * rate_death/1000.0, (0.9 * 1.26 + 0.1 * 3.92) * rate_death/1000.0, #for people starting treatment
                (0.9 * 1.94 + 0.1 * 4.28) * rate_death/1000.0, (0.9 * (p1*8.3+(1-p1)*41.7) + 0.1 * (p1*11.8+(1-p1)*59.7)) * rate_death/1000.0,
                1 * rate_death/1000.0, 1.26 * rate_death/1000.0, 1.94 * rate_death/1000.0, (p1*8.3+(1-p1)*41.7) * rate_death/1000.0,#for people starting treatment
                3.92 * rate_death/1000.0, 3.92 * rate_death/1000.0, 4.28 * rate_death/1000.0, (p1*11.8+(1-p1)*59.7) * rate_death/1000.0),nrow=4,byrow =TRUE) #for people failing treatment
  #CD4 progression, for untreated, treated with NNRTI, treated with PI
  cd4_untreated = c(4.35 / 260.0, 4.35 / 156.0, 4.35 / 182.0)  #untreated
  cd4 = matrix(c(4.35 / 206.0, 4.35 / 133.0, 4.35 / 263.0, #NNRTI, start treatment, right
                 4.35 / 69.0, 4.35 / 70.0, 4.35 / 79.0,#start treatment, left
                 4.35 / 73.0, 4.35 / 61.0, 4.35 / 41.0,#suppressed
                 4.35 / 77.0, 4.35 / 66.0, 4.35 / 96.0,#failing
                 4.35 / 140.0, 4.35 / 98.0, 4.35 / 145.0, #PI,start treatment, right
                 4.35 / 68.0, 4.35 / 81.0, 4.35 / 177.0,#start treatment, left
                 4.35 / 72.0, 4.35 / 59.0, 4.35 / 32.0,#suppressed
                 4.35 / 61.0, 4.35 / 65.0, 4.35 / 69.0),nrow=2,byrow =TRUE) #failing
  #Treatment suppression and failure rates, for NNRTI and PI, by cd4
  treat = matrix(c(2.0 * 4.35 / 30.0, 2.0 * 4.35 / 30.0, 2.0 * 4.35 / 31.0, 2.0 * 4.35 / 34.0, #NNRTI.0, T to S
                   2.0 * 4.35 / 203.0, 2.0 * 4.35 / 198.0, 2.0 * 4.35 / 164.0, 2.0 * 4.35 / 112.0, #T to F
                   4.35 / 767.0, 4.35 / 582.0, 4.35 / 270.0, 4.35 / 96.0, #S to F
                   4.35 / 28.0, 4.35 / 56.0, 4.35 / 62.0, 4.35 / 79.0, #F to S
                   2.0 * 4.35 / 33.0, 2.0 * 4.35 / 33.0, 2.0 * 4.35 / 35.0, 2.0 * 4.35 / 43.0, #PI
                   2.0 * 4.35 / 124.0, 2.0 * 4.35 / 122.0, 2.0 * 4.35 / 103.0, 2.0 * 4.35 / 66.0, #T to F
                   4.35 / 267.0, 4.35 / 178.0, 4.35 / 174.0, 4.35 / 83.0, #S to F
                   4.35 / 10.0, 4.35 / 56.0, 4.35 / 24.0, 4.35 / 51.0), nrow=2,byrow =TRUE)  #F to S
  #CD4 progression and ART efficacy by ART, 1: cd4 progression/efficacy equivalent to NNRTI, 2: equivalent to PI
  art_eq = c(1, 1, 1, 1, 1, 1, 2) 
  #Increase/decrease of efficacy due to NNRTI resistance by ART
  art_res1 = matrix(c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.64, 2.64, 2.64, 2.64, 1.0, 1.0, 1.0),nrow=2,byrow =TRUE)
  art_res2 = matrix(c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 4.9, 4.9, 4.9, 4.9, 1.0, 1.0, 1.0),nrow=2,byrow =TRUE)
  #Increase/decrease of efficacy relative to NNRTI by ART
  art_eff  = c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0) 
  #Resistance
  res = 1/5.0 
  rev = 1/125.0
  
  dy_dt=rep(0,24*16+2)
  #Susceptible 
  dy_dt[1]=0;
  dy_dt[2]=0;
  y[1]=1000000;
  y[2]=1000000;
  
  
  #New infections
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
      dy_dt[pos(3,j,1,l)] = dy_dt[pos(3,j,1,l)] + (1- p_dtg[1]) * (diag2 / diag1 + oi_inc[j] * oi_test) * y[pos(1,j,1,l)];
      #DTG-eligible, women
      dy_dt[pos(2,j,2,l)] = dy_dt[pos(2,j,2,l)] + p_dtg[2] * (diag2 / diag1 * diag_women + oi_inc[j] * oi_test + preg_inc[j] * preg_test) * y[pos(1,j,2,l)];
      #DTG-ineligible, women
      dy_dt[pos(3,j,2,l)] = dy_dt[pos(3,j,2,l)] + (1 - p_dtg[2]) * (diag2 / diag1 * diag_women + oi_inc[j] * oi_test + preg_inc[j] * preg_test) * y[pos(1,j,2,l)];
      
      for(k in 1:2){
        #Treatment initiation
        dy_dt[pos(4,j,k,l)] = dy_dt[pos(4,j,k,l)] + p_tdf * t_1st[j] * y[pos(3,j,k,l)];
        dy_dt[pos(10,j,k,l)] = dy_dt[pos(10,j,k,l)] + (1-p_tdf) * t_1st[j] * y[pos(3,j,k,l)];
        dy_dt[pos(16,j,k,l)] = dy_dt[pos(16,j,k,l)] + t_1st_DTG[j] * y[pos(2,j,k,l)];
        
        dy_dt[pos(3,j,k,l)] = dy_dt[pos(3,j,k,l)] - t_1st[j] * y[pos(3,j,k,l)];
        dy_dt[pos(2,j,k,l)] = dy_dt[pos(2,j,k,l)] - t_1st_DTG[j] * y[pos(2,j,k,l)];
        
        #Treatment switch after failure (no DTG)
        dy_dt[pos(22,j,k,l)] = dy_dt[pos(22,j,k,l)] + t_switch[j] * (y[pos(6,j,k,l)] + y[pos(12,j,k,l)] + y[pos(18,j,k,l)] + y[pos(21,j,k,l)]);
        
        dy_dt[pos(6,j,k,l)] = dy_dt[pos(6,j,k,l)] - t_switch[j] * y[pos(6,j,k,l)];
        dy_dt[pos(12,j,k,l)] = dy_dt[pos(12,j,k,l)] - t_switch[j] * y[pos(12,j,k,l)];
        dy_dt[pos(18,j,k,l)] = dy_dt[pos(18,j,k,l)] - t_switch[j] * y[pos(18,j,k,l)];
        dy_dt[pos(21,j,k,l)] = dy_dt[pos(21,j,k,l)] - t_switch[j] * y[pos(21,j,k,l)];
        
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
    cd4_j = cd4[art_eq[art],];
    for(l in 1:2){
      for(i in 1:4){
        treat_j[i] =  treat[art_eq[art],i] * art_eff[art] / art_res1[l,art];
      }
      for(i in 5:8){
        treat_j[i] =  treat[art_eq[art],i] / art_eff[art] * art_res1[l,art];
      }
      for(i in 9:12){
        treat_j[i] =  treat[art_eq[art],i] * art_eff[art] * art_res2[l,art];
      }
      for(i in 13:16){
        treat_j[i] =  treat[art_eq[art],i] / art_eff[art] / art_res2[l,art];
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
  
  #NNRTI resistance
  #resistance acquisition
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
  
  return (dy_dt)
}
