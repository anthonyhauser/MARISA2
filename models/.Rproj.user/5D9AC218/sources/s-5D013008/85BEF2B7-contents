###################################################################################################################################
#Sensitivity analysis
#Resistance, transmission, and DTG-efficacy
#https://daphnia.ecology.uga.edu/drakelab/wp-content/uploads/2015/07/sensitivity-ebola.pdf
require(lhs) #add the lhs library
h <- 200 #choose number of points
set.seed(6242015)

#####################################################################################################
#Sensitivity variation
#Time to acquiring NNRTI resistance (1/(reversion rate)), in months
res.min <- 3
res.max <- 9
#Time to reversion from NNRTI resistance to wild-type (1/(reversion rate)), in months
rev.min <- 50
rev.max <- 200
#Impact of NNRTI resistance on the efficacy of NNRTI-based regimen
alpha.min <- 1
alpha.max <- 3.1
alpha2.min <- 1
alpha2.max <- 5.1
#Proportion of MSM
pmsm.min <- 0.01
pmsm.max <- 0.1
#Increased risk of transmission in MSM, relative to the general population
riskr.min <- 1
riskr.max <- 5
#Increased HIV prevalence in MSM, relative to the general population
hiv_msm.min <- 1
hiv_msm.max <- 3
#Increase in treatment efficacy of DTG, relative to NNRTI
dtg_eff.min <- 0.84
dtg_eff.max <- 1.25

#####################################################################################################
#Sensitivity matrix
lhs<-maximinLHS(h,8) #simulate
params.set <- cbind(
  nnrti_res = lhs[,1]*(res.max-res.min)+res.min,
  rev = lhs[,2]*(rev.max-rev.min)+rev.min,
  alpha1 = lhs[,3]*(alpha.max-alpha.min)+alpha.min,
  alpha2 = lhs[,3]*(alpha2.max-alpha2.min)+alpha2.min,
  p_msm = lhs[,4]*(pmsm.max-pmsm.min)+pmsm.min, 
  riskr = lhs[,5]*(riskr.max-riskr.min)+riskr.min,
  hiv_msm= lhs[,6]*(hiv_msm.max-hiv_msm.min)+hiv_msm.min,
  dtg_eff= lhs[,7]*(dtg_eff.max-dtg_eff.min)+dtg_eff.min)#increase in hiv prevalence in msm vs het

save(params.set,file = "C:/Users/ahauser/Documents/Step2/Step2_revised/models/params.set.RData")

