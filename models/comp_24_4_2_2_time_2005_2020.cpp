//Run simulation from 2005 to 2018, 15 care stages

// tell R you need Boost and Cpp11
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]

// include Rcpp, it takes care of most other headers you need
#include <Rcpp.h>

// include Boost's odeint and functional
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <functional>
#include <math.h>

using namespace Rcpp;
using namespace std;
using namespace boost::numeric::odeint;


typedef boost::array< double ,386 > state_type;

double t_start=0.0;
double t_end=120.0;
double t_step=1.0;
//int t_size = (int) (t_end-t_start)/t_step+1;
//int t_size = 91;


// Function from state_type to NumericVector
void boost_array_to_nvec2(state_type const& s,Rcpp::NumericVector & tmp) {
  for (size_t i = 0; i < s.size(); ++i) {
    tmp[i] = s[i];
  }
}

// Function from NumericVector to state_type
void nvec_to_boost_array2(Rcpp::NumericVector const& s, state_type & tmp) {
  for (size_t i = 0; i < 386; ++i) {
    tmp[i] = s[i];
  }
}

int pos(int i, int j, int k, int l){
  return 1 + i+(j-1)*24+(k-1)*24*4+(l-1)*24*4*2;
}
double sum_pos(state_type w, vector<int> vect1, vector<int> vect2, vector<int> vect3, vector<int> vect4){
  double r=0;
  for(unsigned int i=0; i<vect1.size(); i++){
    for(unsigned int j=0; j<vect2.size(); j++){
      for(unsigned int k=0; k<vect3.size(); k++){
        for(unsigned int l=0; l<vect4.size(); l++){
          r+=w[pos(vect1[i],vect2[j],vect3[k],vect4[l])];
        }
      }
    }
  }
  return r;
}

double smooth_f(double x, double a, double b, double c, double d){
  if(a==b){
    if(x<a){
      return c;
    }else{
      return d;
    }
  }else{
    if(x<a){
      return c;
    }else if(x>b){
      return d;
    }else{
      return c+(d-c)*(3.0*pow((x-a)/(b-a),2.0)-2.0*pow((x-a)/(b-a),3.0));
    }
  }
}


double lin_f(double x, double a, double b, double c, double d){
  if(a==b){
    if(x<a){
      return c;
    }else{
      return d;
    }
  }else{
    if(x<a){
      return c;
    }else if(x>b){
      return d;
    }else{
      return c+(d-c)*(x-a)/(b-a);
    }
  }
}


void ode( const state_type &y , state_type &dy_dt , const double t , const Rcpp::NumericVector theta) {
  double dtg_1st=theta[1];
  double dtg_switch=theta[2];
  
  double inf1 [] = { 0.008*0.05, 0.003*0.95, 0.003*0.95, 0.0};
  double inf2=3.318999;
  double inf3=0.5440886;
  std::vector<int> all_pos ={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}; //position of all HIV stages
  std::vector<int> inf_pos ={2,3,4,6,7,9,10,12,13,15,16,18,19,21,22,24}; //position of infectious HIV stages knowing their HIV-positive status
  //Diagnosis
  double diag1 = 273.6363; // number of months before diagnosis (without opportunistic disease - oi) in 2005: 22 years
  double diag2 = smooth_f(2005.0 + t/12.0,2005.0,2015.0,1.0,7.710578); // fold increase between 2005 and 2015: 3 years
  double diag_women=1.25; //increase of diagnosis rate (without oi) for women
  double oi_inc []={0.05/12.0, 0.12/12.0, 0.27/12.0, 0.9/12.0}; //oi incidence by cd4
  double oi_test=smooth_f(2005.0 + t/12.0,2005.0,2015.0,0.2,0.8); //proportion of oi test
  double preg_inc []={1.0 * 23.0/(12.0*1000.0),0.96 * 23.0/(12.0*1000.0),0.87 * 23.0/(12.0*1000.0),0.74 * 23.0/(12.0*1000.0)}; //incidence of pregnancy by cd4
  double preg_test=smooth_f(2005.0 + t/12.0,2005.0,2015.0,0.5,0.98); //proportion of pregnancy test
  //Treatment
  double p_dtg [] = {1.0, theta[3]}; //proportion of women opting for DTG
  double p_tdf = 1; //proportion of people opting/being prescribed tdf (over azt)
  double rate_treat = 0.001879695; //free treatment parameter rates fixed in project 1
  double t_cd4_elig [] = {0.4, 0.5, 0.7, 1.0}; //relative initiation rates by cd4 in 2020 (increase to become identical from 2022, Treat-All policy)
  double t_cd4_treat_all [] = {1.0/0.4, 1.0/0.5, 1.0/0.7, 1.0};
  double t_cd4_elig_year1 [] = {2015.5, 2013.5, 2009.5, 2001.0};
  double t_cd4_elig_year2 [] = {2016.5, 2015.5, 2012.5, 2004.0};
  double t_1st [4];
  for(int i=0;i<4;++i){
	  t_1st[i] = rate_treat *smooth_f(2005.0 + t/12.0, 2001.0, 2012.0, 1.0, 200.0/12.0) * //Increase in treatment rate between 2001 and 2012
	  						 lin_f(2005.0 + t/12.0, t_cd4_elig_year1[i],  t_cd4_elig_year1[i], 0.0, t_cd4_elig[i]) * //Change in eligiblity criteria
							 lin_f(2005.0 + t/12.0, 2017.0, 2022.0, 1.0, t_cd4_treat_all); //Increase in treatment rate due to Treat-All policy
  }
  double t_switch [] = {4.35/531.0 /5.0, 4.35/427.0 /5.0, 4.35/294.0 /5.0, 4.35/189.0 /5.0}; //Switching rate to PI, failing individuals (either DTG-ineligible on NNRTI or DTG-eligible on DTG)
  double t_switch_elig [4]; //Switch rate to PI, DTG-eligible
  double t_1st_NNRTI_inel[4];
  double t_1st_NNRTI_elig[4];
  double t_1st_DTG[4];
  double t_switch_DTG_S[4];
  double t_switch_DTG_F[4];
  for(int i=0;i<4;++i){
  	t_1st_DTG [i]= t_1st[i] * smooth_f(2005.0 + t/12.0, 2019.0, 2019.0, 0.0, 1.0) * dtg_1st; //Treatment initiation rate, DTG
	t_1st_NNRTI_inel[i] = t_1st[i];
	t_1st_NNRTI_elig[i] = t_1st[i] - t_1st_DTG[i];
	t_switch_elig[i] = t_switch[i] * smooth_f(2005.0 + t/12.0, 2019.0, 2019.0, 1.0, 1.0 - dtg_1st);
	t_switch_DTG_S[i] = 1.0/12.0 * dtg_switch * smooth_f(2005.0 + t/12.0, 2019.0, 2019.0, 0.0, 1.0);
	 if(dtg_1st==0.0) t_switch_DTG_F[i] = 0.0;
    if(dtg_1st==1.0 & dtg_switch==0.0) t_switch_DTG_F[i] = t_switch[i];
    if(dtg_switch==1.0) t_switch_DTG_F[i] = 1.0/12.0;
      t_switch_DTG_F[i] = t_switch_DTG_F[i] * smooth_f(2005.0 + t/12.0, 2019.0, 2019.0, 0.0, 1.0);
  }
  //Mortality
  double rate_death =0.1632000; //Mortality rate per month per 1000 people (reference, supp people with cd4>500)
  double t_prov;
  double mort_approx [37];
  for(int l=0;l<37;++l){
    t_prov=t-l;
    mort_approx[l] = rate_treat * smooth_f(2005.0 + t_prov/12.0, 2001.0, 2012.0, 1.0, 200.0/12.0) *  //Increase in treatment rate between 2001 and 2012
      lin_f(2005.0 + t_prov/12.0, t_cd4_elig_year1[4],  t_cd4_elig_year2[4], 0.0, t_cd4_elig[4]) *  //Change in eligiblity criteria
      lin_f(2005.0 + t_prov/12.0, 2017.0, 2022.0, 1.0, t_cd4_treat_all[4]);
  }
  double p1=1.0-0.27*exp(-0.05*(mean(mort_approx)-0.005219697)/rate_treat);
  double mu[4][4] = {{1.57 *rate_death/1000.0, 2.0 *rate_death/1000.0, 4.57 *rate_death/1000.0, (p1*40.9+(1-p1)*134.4) *rate_death/1000.0}, //relative mortality for untreated by cd4
  {(0.9 * 1 + 0.1 * 3.92) * rate_death/1000.0, (0.9 * 1.26 + 0.1 * 3.92) * rate_death/1000.0, //for people starting treatment
   (0.9 * 1.94 + 0.1 * 4.28) * rate_death/1000.0, (0.9 * (p1*8.3+(1-p1)*41.7) + 0.1 * (p1*11.8+(1-p1)*59.7)) * rate_death/1000.0},
   {1.0 * rate_death/1000.0, 1.26 * rate_death/1000.0, 1.94 * rate_death/1000.0, (p1*8.3+(1-p1)*41.7) * rate_death/1000.0},  //for suppressed people
   {3.92 * rate_death/1000.0, 3.92 * rate_death/1000.0, 4.28 * rate_death/1000.0, (p1*11.8+(1-p1)*59.7) * rate_death/1000.0}};//for people failing treatment
  //CD4 progression, for untreated, treated with NNRTI, treated with PI
  double cd4_untreated[] = {4.35 / 260.0, 4.35 / 156.0, 4.35 / 182.0}; //untreated
  double cd4_j [12];
  double cd4 [2][12] = { {4.35 / 206.0, 4.35 / 133.0, 4.35 / 263.0, //NNRTI, start treatment, right
                          4.35 / 69.0, 4.35 / 70.0, 4.35 / 79.0,//start treatment, left
                          4.35 / 73.0, 4.35 / 61.0, 4.35 / 41.0,//suppressed
                          4.35 / 77.0, 4.35 / 66.0, 4.35 / 96.0},//failing
                          {4.35 / 140.0, 4.35 / 98.0, 4.35 / 145.0, //PI,start treatment, right
                           4.35 / 68.0, 4.35 / 81.0, 4.35 / 177.0,//start treatment, left
                           4.35 / 72.0, 4.35 / 59.0, 4.35 / 32.0,//suppressed
                           4.35 / 61.0, 4.35 / 65.0, 4.35 / 69.0} };//failing
  //Treatment suppression and failure rates, for NNRTI and PI, by cd4
  double treat_j [16];
  double treat [2][16] = { {2.0 * 4.35 / 30.0, 2.0 * 4.35 / 30.0, 2.0 * 4.35 / 31.0, 2.0 * 4.35 / 34.0, //NNRTI.0, T to S
                            2.0 * 4.35 / 203.0, 2.0 * 4.35 / 198.0, 2.0 * 4.35 / 164.0, 2.0 * 4.35 / 112.0, //T to F
                            4.35 / 767.0, 4.35 / 582.0, 4.35 / 270.0, 4.35 / 96.0, //S to F
                            4.35 / 28.0, 4.35 / 56.0, 4.35 / 62.0, 4.35 / 79.0}, //F to S
                            {2.0 * 4.35 / 33.0, 2.0 * 4.35 / 33.0, 2.0 * 4.35 / 35.0, 2.0 * 4.35 / 43.0, //PI
                             2.0 * 4.35 / 124.0, 2.0 * 4.35 / 122.0, 2.0 * 4.35 / 103.0, 2.0 * 4.35 / 66.0, //T to F
                             4.35 / 267.0, 4.35 / 178.0, 4.35 / 174.0, 4.35 / 83.0, //S to F
                             4.35 / 10.0, 4.35 / 56.0, 4.35 / 24.0, 4.35 / 51.0} }; //F to S
  //CD4 progression and ART efficacy by ART, 1: cd4 progression/efficacy equivalent to NNRTI, 2: equivalent to PI
  int art_eq [] = {1, 1, 1, 1, 1, 1, 2};
  int art_j;
  //Increase/decrease of efficacy due to NNRTI resistance by ART
  double art_res1 [2][7] = {{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, {2.64, 2.64, 2.64, 2.64, 1.0, 1.0, 1.0}};
  double art_res2 [2][7] = {{1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}, {4.9, 4.9, 4.9, 4.9, 1.0, 1.0, 1.0}};
  //Increase/decrease of efficacy relative to NNRTI by ART
  double art_eff [] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  art_eff[5] = theta[4];
  art_eff[6] = theta[4];
  //Resistance
  double res = 1/5.0;
  double rev = 1/125.0;
  
  /////////////////////////////////////////////////////////////////
  //Add smooth spline function
  //
  //
  //
  //
  //
  //
  //
  /////////////////////////////////////////////////////////////////
  
  for(int i=0;i<386;++i){dy_dt[i] = 0;}
  
  // dy_dt[1] = -0.5*y[1];
  // dy_dt[2] = -0.5*y[2];
  // dy_dt[3]=0.5*y[1] + 0.5*y[2];
  
  // for(int l=1;l<3;++l){
  //   dy_dt[pos(1,1,1,l)]= inf1[1] * inf2 * (sum_pos(y,{1},{1,2,3,4},{1},{l})) * y[0] / (y[0] + sum_pos(y,all_pos,{1,2,3,4},{1},{1,2}));
  //   dy_dt[pos(1,1,2,l)]= inf1[1] * inf2 * (sum_pos(y,{1},{1,2,3,4},{1},{l})) * y[1] / (y[1] + sum_pos(y,all_pos,{1,2,3,4},{2},{1,2}));
  // }
  
  for(int l=1;l<3;++l){
    dy_dt[pos(1,1,1,l)]=dy_dt[pos(1,1,1,l)] +
      inf1[0] * inf2 * (sum_pos(y,{1},{1,2,3,4},{1},{l}) + inf3 * sum_pos(y,inf_pos,{1,2,3,4},{1},{l})) * y[0] / (y[0] + sum_pos(y,all_pos,{1,2,3,4},{1},{1,2})) +
      inf1[1] * inf2 * (sum_pos(y,{1},{1,2,3,4},{2},{l}) + inf3 * sum_pos(y,inf_pos,{1,2,3,4},{2},{l})) * y[0] / (y[0] + sum_pos(y,all_pos,{1,2,3,4},{1},{1,2}));
    dy_dt[pos(1,1,2,l)]=dy_dt[pos(1,1,2,l)] +
      inf1[2] * inf2 * (sum_pos(y,{1},{1,2,3,4},{1},{l}) + inf3 * sum_pos(y,inf_pos,{1,2,3,4},{1},{l})) * y[1] / (y[1] + sum_pos(y,all_pos,{1,2,3,4},{2},{1,2})) +
      inf1[3] * inf2 * (sum_pos(y,{1},{1,2,3,4},{2},{l}) + inf3 * sum_pos(y,inf_pos,{1,2,3,4},{2},{l})) * y[1] / (y[1] + sum_pos(y,all_pos,{1,2,3,4},{2},{1,2}));
  }
  
  for(int j=1;j<5;++j){
    for(int l=1;l<3;++l){
      //Diagnosis
      //men
      dy_dt[pos(1,j,1,l)] = dy_dt[pos(1,j,1,l)] - (diag2 / diag1 + oi_inc[j-1] * oi_test) * y[pos(1,j,1,l)];
      //women
      dy_dt[pos(1,j,2,l)] = dy_dt[pos(1,j,2,l)] - (diag2 / diag1 * diag_women + oi_inc[j-1] * oi_test + preg_inc[j-1] * preg_test) * y[pos(1,j,2,l)];
      
      //DTG-eligible, men
      dy_dt[pos(2,j,1,l)] = dy_dt[pos(2,j,1,l)] + p_dtg[0] * (diag2 / diag1 + oi_inc[j-1] * oi_test) * y[pos(1,j,1,l)];
      //DTG-ineligible, men
      dy_dt[pos(3,j,1,l)] = dy_dt[pos(3,j,1,l)] + (1 - p_dtg[0]) * (diag2 / diag1 + oi_inc[j-1] * oi_test) * y[pos(1,j,1,l)];
      //DTG-eligible, women
      dy_dt[pos(2,j,2,l)] = dy_dt[pos(2,j,2,l)] + p_dtg[1] * (diag2 / diag1 * diag_women + oi_inc[j-1] * oi_test + preg_inc[j-1] * preg_test) * y[pos(1,j,2,l)];
      //DTG-ineligible, women
      dy_dt[pos(3,j,2,l)] = dy_dt[pos(3,j,2,l)] + (1 - p_dtg[1]) * (diag2 / diag1 * diag_women + oi_inc[j-1] * oi_test + preg_inc[j-1] * preg_test) * y[pos(1,j,2,l)];
      
      
      for(int k=1;k<3;++k){
        //Treatment initiation
		//NNRTI: 1) NNRTI+TDF for DTG-inel, 2) NNRTI+TDF for DTG-elig, 3) NNRTI+AZT for DTG-inel, 4) NNRTI+AZT for DTG-elig
		dy_dt[pos(4,j,k,l)] = dy_dt[pos(4,j,k,l)] + p_tdf * t_1st_NNRTI_inel[j-1] * y[pos(3,j,k,l)];
        dy_dt[pos(7,j,k,l)] = dy_dt[pos(7,j,k,l)] + p_tdf * t_1st_NNRTI_elig[j-1] * y[pos(2,j,k,l)];
		dy_dt[pos(10,j,k,l)] = dy_dt[pos(10,j,k,l)] + (1-p_tdf) * t_1st_NNRTI_inel[j-1] * y[pos(3,j,k,l)];
        dy_dt[pos(13,j,k,l)] = dy_dt[pos(13,j,k,l)] + (1-p_tdf) * t_1st_NNRTI_elig[j-1] * y[pos(2,j,k,l)];
		//DTG: DTG+TDF
        dy_dt[pos(16,j,k,l)] = dy_dt[pos(16,j,k,l)] + t_1st_DTG[j-1] * y[pos(2,j,k,l)];
        
        dy_dt[pos(3,j,k,l)] = dy_dt[pos(3,j,k,l)] - t_1st_NNRTI_inel[j-1] * y[pos(3,j,k,l)];
        dy_dt[pos(2,j,k,l)] = dy_dt[pos(2,j,k,l)] - (t_1st_NNRTI_elig[j-1] + t_1st_DTG[j-1]) * y[pos(2,j,k,l)];
        
        //Treatment switch to PI after failure
        //Switch from NNRTI to PI, DTG-ineligible
        dy_dt[pos(22,j,k,l)] = dy_dt[pos(22,j,k,l)] + t_switch[j-1] * (y[pos(6,j,k,l)] + y[pos(12,j,k,l)]);
        dy_dt[pos(6,j,k,l)] = dy_dt[pos(6,j,k,l)] - t_switch[j-1] * y[pos(6,j,k,l)];
        dy_dt[pos(12,j,k,l)] = dy_dt[pos(12,j,k,l)] - t_switch[j-1] * y[pos(12,j,k,l)];
        //Switch from DTG to PI
        dy_dt[pos(22,j,k,l)] = dy_dt[pos(22,j,k,l)] + t_switch[j-1] * (y[pos(18,j,k,l)] + y[pos(21,j,k,l)]);
        dy_dt[pos(18,j,k,l)] = dy_dt[pos(18,j,k,l)] - t_switch[j-1] * y[pos(18,j,k,l)];
        dy_dt[pos(21,j,k,l)] = dy_dt[pos(21,j,k,l)] - t_switch[j-1] * y[pos(21,j,k,l)];
        //Switch from NNRTI to PI, DTG-eligible
        dy_dt[pos(22,j,k,l)] = dy_dt[pos(22,j,k,l)] + t_switch_elig[j-1] * (y[pos(9,j,k,l)] + y[pos(15,j,k,l)]);
        dy_dt[pos(9,j,k,l)] = dy_dt[pos(9,j,k,l)] - t_switch_elig[j-1] * y[pos(9,j,k,l)];
        dy_dt[pos(15,j,k,l)] = dy_dt[pos(15,j,k,l)] - t_switch_elig[j-1] * y[pos(15,j,k,l)];
        
        //Treatment switch to DTG, suppressed
        dy_dt[pos(17,j,k,l)] = dy_dt[pos(17,j,k,l)] + t_switch_DTG_S[j-1] * (y[pos(8,j,k,l)] + y[pos(14,j,k,l)]);
        dy_dt[pos(8,j,k,l)] = dy_dt[pos(8,j,k,l)] - t_switch_DTG_S[j-1] * y[pos(8,j,k,l)];
        dy_dt[pos(14,j,k,l)] = dy_dt[pos(14,j,k,l)] - t_switch_DTG_S[j-1] * y[pos(14,j,k,l)];
        
        //Treatment switch to DTG, failed
        dy_dt[pos(16,j,k,l)] = dy_dt[pos(16,j,k,l)] + t_switch_DTG_F[j-1] * y[pos(15,j,k,l)];
        dy_dt[pos(19,j,k,l)] = dy_dt[pos(19,j,k,l)] + t_switch_DTG_F[j-1] * y[pos(9,j,k,l)];
        
        dy_dt[pos(15,j,k,l)] = dy_dt[pos(15,j,k,l)] - t_switch_DTG_F[j-1] * y[pos(15,j,k,l)];
        dy_dt[pos(9,j,k,l)] = dy_dt[pos(9,j,k,l)] - t_switch_DTG_F[j-1] * y[pos(9,j,k,l)];
      }
    }
  }
  
  //Treatment stages and CD4, mortality
  //Infected and diagnosed
  for(int k=1;k<3;++k){
    for(int l=1;l<3;++l){
      for(int i=1;i<4;++i){
        dy_dt[pos(i,1,k,l)] = dy_dt[pos(i,1,k,l)] - cd4_untreated[1-1] * y[pos(i,1,k,l)];
        dy_dt[pos(i,2,k,l)] = dy_dt[pos(i,2,k,l)] + cd4_untreated[1-1] * y[pos(i,1,k,l)] - cd4_untreated[2-1] * y[pos(i,2,k,l)];
        dy_dt[pos(i,3,k,l)] = dy_dt[pos(i,3,k,l)] + cd4_untreated[2-1] * y[pos(i,2,k,l)] - cd4_untreated[3-1] * y[pos(i,3,k,l)];
        dy_dt[pos(i,4,k,l)] = dy_dt[pos(i,4,k,l)] + cd4_untreated[3-1] * y[pos(i,3,k,l)];
        
        //Mortality
        for(int j=1;j<5;++j){
          dy_dt[pos(i,j,k,l)] = dy_dt[pos(i,j,k,l)] - mu[1-1][j-1] * y[pos(i,j,k,l)];
        }
      }
    }
  }
  
  //On ART
  for(int art=1;art<8;++art){
    art_j=art_eq[art-1];
    for(int j=1;j<13;++j){
      cd4_j[j-1] = cd4[art_j-1][j-1];
    }
    for(int l=1;l<3;++l){
      for(int i=1;i<5;++i){
        treat_j[i-1] =  treat[art_j-1][i-1] * art_eff[art-1] / art_res1[l-1][art-1];
      }
      for(int i=5;i<9;++i){
        treat_j[i-1] =  treat[art_j-1][i-1] / art_eff[art-1] * art_res1[l-1][art-1];
      }
      for(int i=9;i<13;++i){
        treat_j[i-1] =  treat[art_j-1][i-1] * art_eff[art-1] * art_res2[l-1][art-1];
      }
      for(int i=13;i<17;++i){
        treat_j[i-1] =  treat[art_j-1][i-1] / art_eff[art-1] / art_res2[l-1][art-1];
      }
      
      for(int k=1;k<3;++k){
        //Start
        dy_dt[pos(3*art+1,1,k,l)] = dy_dt[pos(3*art+1,1,k,l)] + cd4_j[4-1] * y[pos(3*art+1,2,k,l)] - cd4_j[1-1] * y[pos(3*art+1,1,k,l)]+
          -(treat_j[1-1] + treat_j[5-1]) * y[pos(3*art+1,1,k,l)];
          dy_dt[pos(3*art+1,2,k,l)] = dy_dt[pos(3*art+1,2,k,l)] + cd4_j[1-1] * y[pos(3*art+1,1,k,l)] + cd4_j[5-1] * y[pos(3*art+1,3,k,l)] - (cd4_j[4-1] + cd4_j[2-1]) * y[pos(3*art+1,2,k,l)]+
          -(treat_j[2-1] + treat_j[6-1]) * y[pos(3*art+1,2,k,l)];
          dy_dt[pos(3*art+1,3,k,l)] = dy_dt[pos(3*art+1,3,k,l)] + cd4_j[2-1] * y[pos(3*art+1,2,k,l)] + cd4_j[6-1] * y[pos(3*art+1,4,k,l)] - (cd4_j[5-1] + cd4_j[3-1]) * y[pos(3*art+1,3,k,l)]+
          -(treat_j[3-1] + treat_j[7-1]) * y[pos(3*art+1,3,k,l)];
          dy_dt[pos(3*art+1,4,k,l)] = dy_dt[pos(3*art+1,4,k,l)] + cd4_j[3-1] * y[pos(3*art+1,3,k,l)] - cd4_j[6-1] * y[pos(3*art+1,4,k,l)]+
          -(treat_j[4-1] + treat_j[8-1]) * y[pos(3*art+1,4,k,l)];
          //Suppressed
          dy_dt[pos(3*art+2,1,k,l)] = dy_dt[pos(3*art+2,1,k,l)] + cd4_j[7-1] * y[pos(3*art+2,2,k,l)]+
          + treat_j[1-1] * y[pos(3*art+1,1,k,l)] + treat_j[13-1] * y[pos(3*art+3,1,k,l)] - treat_j[9-1] * y[pos(3*art+2,1,k,l)];
          dy_dt[pos(3*art+2,2,k,l)] = dy_dt[pos(3*art+2,2,k,l)] + cd4_j[8-1] * y[pos(3*art+2,3,k,l)] - cd4_j[7-1] * y[pos(3*art+2,2,k,l)]+
          + treat_j[2-1] * y[pos(3*art+1,2,k,l)] + treat_j[14-1] * y[pos(3*art+3,2,k,l)] - treat_j[10-1] * y[pos(3*art+2,2,k,l)];
          dy_dt[pos(3*art+2,3,k,l)] = dy_dt[pos(3*art+2,3,k,l)] + cd4_j[9-1] * y[pos(3*art+2,4,k,l)] - cd4_j[8-1] * y[pos(3*art+2,3,k,l)]+
          + treat_j[3-1] * y[pos(3*art+1,3,k,l)] + treat_j[15-1] * y[pos(3*art+3,3,k,l)] - treat_j[11-1] * y[pos(3*art+2,3,k,l)];
          dy_dt[pos(3*art+2,4,k,l)] = dy_dt[pos(3*art+2,4,k,l)] - cd4_j[9-1] * y[pos(3*art+2,4,k,l)]+
          + treat_j[4-1] * y[pos(3*art+1,4,k,l)] + treat_j[16-1] * y[pos(3*art+3,4,k,l)] - treat_j[12-1] * y[pos(3*art+2,4,k,l)];
          //Failed
          dy_dt[pos(3*art+3,1,k,l)] = dy_dt[pos(3*art+3,1,k,l)] - cd4_j[10-1] * y[pos(3*art+3,1,k,l)]+
          treat_j[5-1] * y[pos(3*art+1,1,k,l)] + treat_j[9-1] * y[pos(3*art+2,1,k,l)] - treat_j[13-1] * y[pos(3*art+3,1,k,l)];
          dy_dt[pos(3*art+3,2,k,l)] = dy_dt[pos(3*art+3,2,k,l)] + cd4_j[10-1] * y[pos(3*art+3,1,k,l)] - cd4_j[11-1] * y[pos(3*art+3,2,k,l)]+
            treat_j[6-1] * y[pos(3*art+1,2,k,l)] + treat_j[10-1] * y[pos(3*art+2,2,k,l)] - treat_j[14-1] * y[pos(3*art+3,2,k,l)];
          dy_dt[pos(3*art+3,3,k,l)] = dy_dt[pos(3*art+3,3,k,l)] + cd4_j[11-1] * y[pos(3*art+3,2,k,l)] - cd4_j[12-1] * y[pos(3*art+3,3,k,l)]+
            treat_j[7-1] * y[pos(3*art+1,3,k,l)] + treat_j[11-1] * y[pos(3*art+2,3,k,l)] - treat_j[15-1] * y[pos(3*art+3,3,k,l)];
          dy_dt[pos(3*art+3,4,k,l)] = dy_dt[pos(3*art+3,4,k,l)] + cd4_j[12-1] * y[pos(3*art+3,3,k,l)]+
            treat_j[8-1] * y[pos(3*art+1,4,k,l)] + treat_j[12-1] * y[pos(3*art+2,4,k,l)] - treat_j[16-1] * y[pos(3*art+3,4,k,l)];
          
          //Mortality 
          for(int j=1;j<5;++j){
            dy_dt[pos(3*art+1,j,k,l)] = dy_dt[pos(3*art+1,j,k,l)] - mu[2-1][j-1] * y[pos(3*art+1,j,k,l)];
            dy_dt[pos(3*art+2,j,k,l)] = dy_dt[pos(3*art+2,j,k,l)] - mu[3-1][j-1] * y[pos(3*art+2,j,k,l)];
            dy_dt[pos(3*art+3,j,k,l)] = dy_dt[pos(3*art+3,j,k,l)] - mu[4-1][j-1] * y[pos(3*art+3,j,k,l)];
          }
      }
    }
  }
  
  for(int j=1;j<5;++j){
    for(int k=1;k<3;++k){
      dy_dt[pos(6,j,k,2)] = dy_dt[pos(6,j,k,2)] + res * y[pos(6,j,k,1)];
      dy_dt[pos(9,j,k,2)] = dy_dt[pos(9,j,k,2)] + res * y[pos(9,j,k,1)];
      dy_dt[pos(12,j,k,2)] = dy_dt[pos(12,j,k,2)] + res * y[pos(12,j,k,1)];
      dy_dt[pos(15,j,k,2)] = dy_dt[pos(15,j,k,2)] + res * y[pos(15,j,k,1)];
      
      dy_dt[pos(6,j,k,1)] = dy_dt[pos(6,j,k,1)] - res * y[pos(6,j,k,1)];
      dy_dt[pos(9,j,k,1)] = dy_dt[pos(9,j,k,1)] - res * y[pos(9,j,k,1)];
      dy_dt[pos(12,j,k,1)] = dy_dt[pos(12,j,k,1)] - res * y[pos(12,j,k,1)];
      dy_dt[pos(15,j,k,1)] = dy_dt[pos(15,j,k,1)] - res * y[pos(15,j,k,1)];
      
      //resistance reversion
      dy_dt[pos(1,j,k,1)] = dy_dt[pos(1,j,k,1)] + rev * y[pos(1,j,k,2)];
      dy_dt[pos(2,j,k,1)] = dy_dt[pos(2,j,k,1)] + rev * y[pos(2,j,k,2)];
      dy_dt[pos(3,j,k,1)] = dy_dt[pos(3,j,k,1)] + rev * y[pos(3,j,k,2)];
      
      dy_dt[pos(1,j,k,2)] = dy_dt[pos(1,j,k,2)] - rev * y[pos(1,j,k,2)];
      dy_dt[pos(2,j,k,2)] = dy_dt[pos(2,j,k,2)] - rev * y[pos(2,j,k,2)];
      dy_dt[pos(3,j,k,2)] = dy_dt[pos(3,j,k,2)] - rev * y[pos(3,j,k,2)];
    }
  }
  
}


// Function to write the results
Rcpp::NumericMatrix data(121, 386);
Rcpp::NumericVector x2(386);
void write_cout_2( const state_type &x , const double t ) {
  int t_ind = t;
  boost_array_to_nvec2(x,x2);
  data(t_ind,_)= x2;
}


typedef runge_kutta_dopri5< state_type > stepper_type;




// Function that depends on theta, return a function that depends on x dxdt and t
state_type y_temp;
Rcpp::NumericVector nvec(386);
std::function<void(const state_type&, state_type&, const double)> ode2(const Rcpp::NumericVector theta) {
  return [&theta](const state_type &y, state_type &dydt, const double t) {
    ode(y,dydt,t,theta);
    //for(int i=0;i<386;++i) dydt[i]=y[i];
  };
}

// [[Rcpp::export]]
Rcpp::NumericMatrix boostExample(const Rcpp::NumericVector vs, const Rcpp::NumericVector theta) {
  nvec_to_boost_array2(vs,y_temp); // initial conditions
  integrate_const(make_dense_output( 1E-9 , 1E-9 , stepper_type () ) ,
                  ode2(theta) , y_temp ,t_start , t_end , t_step , write_cout_2 );
  return data;
}