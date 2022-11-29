





DATA_SECTION

!!cout << "Starting Beaufort Assessment Model" << endl;
!!cout << endl;
!!cout << "                BAM!" << endl;
!!cout << endl;






init_int styr;
init_int endyr;


init_int styr_rec_dev;

init_int endyr_rec_dev;

init_int endyr_rec_phase1;
init_int endyr_rec_phase2;


init_int endyr_selex_phase1;
init_int endyr_selex_phase2;


init_number sizelim1;
init_number sizelim2;


number nyrs;
number nyrs_rec;


 LOCAL_CALCS
   nyrs=endyr-styr+1.;
   nyrs_rec=endyr_rec_dev-styr_rec_dev+1.;
 END_CALCS
 

init_int nages;

init_vector agebins(1,nages);
 

init_int nages_agec;


init_vector agebins_agec(1,nages_agec);


init_int nlenbins;          
init_number lenbins_width;  


init_vector lenbins(1,nlenbins);
 

init_number max_F_spr_msy;

init_int n_iter_spr;

int n_iter_msy;
 LOCAL_CALCS
		n_iter_msy=n_iter_spr; 
 END_CALCS


init_int styr_rec_spr;

init_int endyr_rec_spr;
number nyrs_rec_spr;
 LOCAL_CALCS
   nyrs_rec_spr=endyr_rec_spr-styr_rec_spr+1.;
 END_CALCS 
 

init_int selpar_n_yrs_wgted;

init_number set_BiasCor;







init_int styr_cpue_sCT;
init_int endyr_cpue_sCT;
init_vector obs_cpue_sCT(styr_cpue_sCT,endyr_cpue_sCT);   
init_vector obs_cv_cpue_sCT(styr_cpue_sCT,endyr_cpue_sCT);    


init_int nyr_lenc_sCT;
init_ivector yrs_lenc_sCT(1,nyr_lenc_sCT);
init_vector nsamp_lenc_sCT(1,nyr_lenc_sCT);
init_vector nfish_lenc_sCT(1,nyr_lenc_sCT);
init_matrix obs_lenc_sCT(1,nyr_lenc_sCT,1,nlenbins);


init_int nyr_agec_sCT;
init_ivector yrs_agec_sCT(1,nyr_agec_sCT);
init_vector nsamp_agec_sCT(1,nyr_agec_sCT);
init_vector nfish_agec_sCT(1,nyr_agec_sCT);
init_matrix obs_agec_sCT(1,nyr_agec_sCT,1,nages_agec);




init_int styr_cpue_cHL;                                             
init_int endyr_cpue_cHL;                                            
init_vector obs_cpue_cHL(styr_cpue_cHL,endyr_cpue_cHL); 
init_vector obs_cv_cpue_cHL(styr_cpue_cHL,endyr_cpue_cHL);  


init_int styr_L_cHL;
init_int endyr_L_cHL;
init_vector obs_L_cHL(styr_L_cHL,endyr_L_cHL);
init_vector obs_cv_L_cHL(styr_L_cHL,endyr_L_cHL);


init_int nyr_lenc_cHL;
init_ivector yrs_lenc_cHL(1,nyr_lenc_cHL);
init_vector nsamp_lenc_cHL(1,nyr_lenc_cHL);
init_vector nfish_lenc_cHL(1,nyr_lenc_cHL);
init_matrix obs_lenc_cHL(1,nyr_lenc_cHL,1,nlenbins);
 

init_int nyr_agec_cHL;
init_ivector yrs_agec_cHL(1,nyr_agec_cHL);
init_vector nsamp_agec_cHL(1,nyr_agec_cHL);
init_vector nfish_agec_cHL(1,nyr_agec_cHL);
init_matrix obs_agec_cHL(1,nyr_agec_cHL,1,nages_agec);


init_int styr_D_cHL;
init_int endyr_D_cHL;
init_vector obs_released_cHL(styr_D_cHL,endyr_D_cHL);
init_vector obs_cv_D_cHL(styr_D_cHL,endyr_D_cHL);




init_int styr_L_cOT;
init_int endyr_L_cOT;
init_vector obs_L_cOT(styr_L_cOT,endyr_L_cOT);
init_vector obs_cv_L_cOT(styr_L_cOT,endyr_L_cOT);


init_int nyr_lenc_cOT;
init_ivector yrs_lenc_cOT(1,nyr_lenc_cOT);
init_vector nsamp_lenc_cOT(1,nyr_lenc_cOT);
init_vector nfish_lenc_cOT(1,nyr_lenc_cOT);
init_matrix obs_lenc_cOT(1,nyr_lenc_cOT,1,nlenbins);




init_int styr_cpue_rHB;                                             
init_int endyr_cpue_rHB;                                            
init_vector obs_cpue_rHB(styr_cpue_rHB,endyr_cpue_rHB);
init_vector obs_cv_cpue_rHB(styr_cpue_rHB,endyr_cpue_rHB); 


init_int styr_L_rHB;
init_int endyr_L_rHB;
init_vector obs_L_rHB(styr_L_rHB,endyr_L_rHB);   
init_vector obs_cv_L_rHB(styr_L_rHB,endyr_L_rHB);    


init_int nyr_lenc_rHB;
init_ivector yrs_lenc_rHB(1,nyr_lenc_rHB);
init_vector nsamp_lenc_rHB(1,nyr_lenc_rHB);
init_vector nfish_lenc_rHB(1,nyr_lenc_rHB);
init_matrix obs_lenc_rHB(1,nyr_lenc_rHB,1,nlenbins);


init_int nyr_agec_rHB;
init_ivector yrs_agec_rHB(1,nyr_agec_rHB);
init_vector nsamp_agec_rHB(1,nyr_agec_rHB);
init_vector nfish_agec_rHB(1,nyr_agec_rHB);
init_matrix obs_agec_rHB(1,nyr_agec_rHB,1,nages_agec);


init_int styr_D_rHB;
init_int endyr_D_rHB;
init_vector obs_released_rHB(styr_D_rHB,endyr_D_rHB);
init_vector obs_cv_D_rHB(styr_D_rHB,endyr_D_rHB);


init_int nyr_lenc_rHB_D;
init_ivector yrs_lenc_rHB_D(1,nyr_lenc_rHB_D);
init_vector nsamp_lenc_rHB_D(1,nyr_lenc_rHB_D);
init_vector nfish_lenc_rHB_D(1,nyr_lenc_rHB_D);
init_matrix obs_lenc_rHB_D(1,nyr_lenc_rHB_D,1,nlenbins);











init_int styr_L_rGN;
init_int endyr_L_rGN;
init_vector obs_L_rGN(styr_L_rGN,endyr_L_rGN);   
init_vector obs_cv_L_rGN(styr_L_rGN,endyr_L_rGN);    


init_int styr_D_rGN;
init_int endyr_D_rGN;
init_vector obs_released_rGN(styr_D_rGN,endyr_D_rGN);
init_vector obs_cv_D_rGN(styr_D_rGN,endyr_D_rGN);


init_int nyr_lenc_rGN;
init_ivector yrs_lenc_rGN(1,nyr_lenc_rGN);
init_vector nsamp_lenc_rGN(1,nyr_lenc_rGN);
init_vector nfish_lenc_rGN(1,nyr_lenc_rGN);
init_matrix obs_lenc_rGN(1,nyr_lenc_rGN,1,nlenbins);


init_int nyr_agec_rGN;
init_ivector yrs_agec_rGN(1,nyr_agec_rGN);
init_vector nsamp_agec_rGN(1,nyr_agec_rGN);
init_vector nfish_agec_rGN(1,nyr_agec_rGN);
init_matrix obs_agec_rGN(1,nyr_agec_rGN,1,nages_agec);








init_vector set_Linf(1,7);
init_vector set_K(1,7);
init_vector set_t0(1,7);
init_vector set_len_cv(1,7);




init_vector set_M_constant(1,7);     

init_vector set_steep(1,7);         
init_vector set_log_R0(1,7);        
init_vector set_R_autocorr(1,7);    
init_vector set_rec_sigma(1,7);     

init_vector set_log_dm_lenc_cHL(1,7);    
init_vector set_log_dm_lenc_cOT(1,7);    
init_vector set_log_dm_lenc_rHB(1,7);    
init_vector set_log_dm_lenc_rGN(1,7);    
init_vector set_log_dm_lenc_rHB_D(1,7);  
init_vector set_log_dm_lenc_sCT(1,7);   
init_vector set_log_dm_agec_cHL(1,7);    
init_vector set_log_dm_agec_rHB(1,7);    
init_vector set_log_dm_agec_rGN(1,7);    
init_vector set_log_dm_agec_sCT(1,7);   





init_vector set_selpar_A50_cHL2(1,7);
init_vector set_selpar_slope_cHL2(1,7);
init_vector set_selpar_A50_cHL3(1,7); 
init_vector set_selpar_slope_cHL3(1,7);


init_vector set_selpar_A50_cOT2(1,7);
init_vector set_selpar_A50_cOT3(1,7);
init_vector set_selpar_slope_cOT2(1,7);
init_vector set_selpar_A502_cOT2(1,7); 
init_vector set_selpar_slope2_cOT2(1,7);

init_vector set_selpar_A50_rHB1(1,7);
init_vector set_selpar_slope_rHB1(1,7);
init_vector set_selpar_A50_rHB2(1,7);
init_vector set_selpar_slope_rHB2(1,7);
init_vector set_selpar_A50_rHB3(1,7); 
init_vector set_selpar_slope_rHB3(1,7);

init_vector set_selpar_A50_rGN3(1,7); 
init_vector set_selpar_slope_rGN3(1,7);

init_vector set_selpar_A50_sCT(1,7); 
init_vector set_selpar_slope_sCT(1,7);
init_vector set_selpar_A502_sCT(1,7); 
init_vector set_selpar_slope2_sCT(1,7);

init_vector set_selpar_age1logit_D(1,7);



init_vector set_log_q_cpue_cHL(1,7);      
init_vector set_log_q_cpue_rHB(1,7);      

init_vector set_log_q_cpue_sCT(1,7);     


init_vector set_log_avg_F_L_cHL(1,7);
init_vector set_log_avg_F_L_cOT(1,7);
init_vector set_log_avg_F_L_rHB(1,7);
init_vector set_log_avg_F_L_rGN(1,7);
init_vector set_log_avg_F_D_cHL(1,7);
init_vector set_log_avg_F_D_rHB(1,7);
init_vector set_log_avg_F_D_rGN(1,7);



init_vector set_log_dev_F_L_cHL(1,3); 
init_vector set_log_dev_F_L_cOT(1,3); 
init_vector set_log_dev_F_L_rHB(1,3);
init_vector set_log_dev_F_L_rGN(1,3);
init_vector set_log_dev_F_D_cHL(1,3); 
init_vector set_log_dev_F_D_rHB(1,3);
init_vector set_log_dev_F_D_rGN(1,3);
init_vector set_log_dev_RWq(1,3);
init_vector set_log_dev_rec(1,3);
init_vector set_log_dev_Nage(1,3);

init_vector set_log_dev_vals_F_L__cHL(styr_L_cHL,endyr_L_cHL);
init_vector set_log_dev_vals_F_L__cOT(styr_L_cOT,endyr_L_cOT);
init_vector set_log_dev_vals_F_L__rHB(styr_L_rHB,endyr_L_rHB);
init_vector set_log_dev_vals_F_L__rGN(styr_L_rGN,endyr_L_rGN);
init_vector set_log_dev_vals_F_D__cHvals(styr_D_cHL,endyr_D_cHL);
init_vector set_log_dev_vals_F_D__HBvals(styr_D_rHB,endyr_D_rHB);
init_vector set_log_dev_vals_F_D__GRvals(styr_D_rGN,endyr_D_rGN);
init_vector set_log_dev_vals_rec(styr_rec_dev,endyr_rec_dev);
init_vector set_log_dev_vals_Nage(2,nages);             





init_number set_w_L;
init_number set_w_D;            
init_number set_w_cpue_cHL;         
init_number set_w_cpue_rHB;         

init_number set_w_cpue_sCT;        
init_number set_w_lenc_cHL;        
init_number set_w_lenc_cOT;        
init_number set_w_lenc_rHB;        
init_number set_w_lenc_rGN;		
init_number set_w_lenc_rHB_D;		
init_number set_w_lenc_sCT;		
init_number set_w_agec_cHL;        
init_number set_w_agec_rHB;        
init_number set_w_agec_rGN;		
init_number set_w_agec_sCT;		
init_number set_w_Nage_init;    
init_number set_w_rec;          
init_number set_w_rec_early;    
init_number set_w_rec_end;      
init_number set_w_fullF;        
init_number set_w_Ftune;        





init_number wgtpar_a;
init_number wgtpar_b;



init_vector obs_maturity_f(1,nages);            
init_vector obs_maturity_m(1,nages);            
init_vector obs_prop_f(1,nages);
init_number spawn_time_frac; 


init_vector set_M(1,nages);     


init_number set_Dmort_cHL;
init_number set_Dmort_rHB;
init_number set_Dmort_rGN;


init_int SR_switch;


init_int set_q_rate_phase;  
init_number set_q_rate;

init_int set_q_DD_phase;      
init_number set_q_DD_beta;    
init_number set_q_DD_beta_se;
init_int set_q_DD_stage;      


init_number set_RWq_var;     


init_number set_Ftune;
init_int set_Ftune_yr;


init_number minSS_lenc_cHL;
init_number minSS_lenc_cOT;
init_number minSS_lenc_rHB;
init_number minSS_lenc_rHB_D;
init_number minSS_lenc_rGN;
init_number minSS_lenc_sCT;


init_number minSS_agec_cHL;
init_number minSS_agec_rHB;
init_number minSS_agec_rGN;
init_number minSS_agec_sCT;


init_int endyr_proj; 
init_int styr_regs;  
init_int Fproj_switch; 
init_number Fproj_mult; 
int styr_proj;
 LOCAL_CALCS
   styr_proj=endyr+1;
 END_CALCS
 

init_matrix age_error(1,nages,1,nages);


int iyear;
int iage;
int ilen;
int ff;

number sqrt2pi;
number g2mt;                    
number g2kg;                    
number g2klb;                   
number mt2klb;                  
number mt2lb;                   
number dzero;                   
number huge_number;             

init_number end_of_data_file;

 LOCAL_CALCS
   if(end_of_data_file!=999)
   {
       cout << "*** WARNING: Data File NOT READ CORRECTLY ****" << endl;
       exit(0);  
   }
   else
   {cout << "Data File read correctly" << endl;} 
 END_CALCS   


PARAMETER_SECTION

 LOCAL_CALCS
  const double Linf_LO=set_Linf(2); const double Linf_HI=set_Linf(3); const double Linf_PH=set_Linf(4);
  const double K_LO=set_K(2); const double K_HI=set_K(3); const double K_PH=set_K(4);
  const double t0_LO=set_t0(2); const double t0_HI=set_t0(3); const double t0_PH=set_t0(4);  
  const double len_cv_LO=set_len_cv(2); const double len_cv_HI=set_len_cv(3); const double len_cv_PH=set_len_cv(4); 
     
  const double M_constant_LO=set_M_constant(2); const double M_constant_HI=set_M_constant(3); const double M_constant_PH=set_M_constant(4);        
  const double steep_LO=set_steep(2); const double steep_HI=set_steep(3); const double steep_PH=set_steep(4);
  const double log_R0_LO=set_log_R0(2); const double log_R0_HI=set_log_R0(3); const double log_R0_PH=set_log_R0(4);
  const double R_autocorr_LO=set_R_autocorr(2); const double R_autocorr_HI=set_R_autocorr(3); const double R_autocorr_PH=set_R_autocorr(4);
  const double rec_sigma_LO=set_rec_sigma(2); const double rec_sigma_HI=set_rec_sigma(3); const double rec_sigma_PH=set_rec_sigma(4);
  
  const double log_dm_cHL_lc_LO=set_log_dm_lenc_cHL(2); const double log_dm_cHL_lc_HI=set_log_dm_lenc_cHL(3); const double log_dm_cHL_lc_PH=set_log_dm_lenc_cHL(4);
  const double log_dm_cOT_lc_LO=set_log_dm_lenc_cOT(2); const double log_dm_cOT_lc_HI=set_log_dm_lenc_cOT(3); const double log_dm_cOT_lc_PH=set_log_dm_lenc_cOT(4);
  const double log_dm_rHB_lc_LO=set_log_dm_lenc_rHB(2); const double log_dm_rHB_lc_HI=set_log_dm_lenc_rHB(3); const double log_dm_rHB_lc_PH=set_log_dm_lenc_rHB(4);
  const double log_dm_rGN_lc_LO=set_log_dm_lenc_rGN(2); const double log_dm_rGN_lc_HI=set_log_dm_lenc_rGN(3); const double log_dm_rGN_lc_PH=set_log_dm_lenc_rGN(4);
  const double log_dm_rHB_D_lc_LO=set_log_dm_lenc_rHB_D(2); const double log_dm_rHB_D_lc_HI=set_log_dm_lenc_rHB_D(3); const double log_dm_rHB_D_lc_PH=set_log_dm_lenc_rHB_D(4);
  const double log_dm_sCT_lc_LO=set_log_dm_lenc_sCT(2); const double log_dm_sCT_lc_HI=set_log_dm_lenc_sCT(3); const double log_dm_sCT_lc_PH=set_log_dm_lenc_sCT(4);
  const double log_dm_cHL_ac_LO=set_log_dm_agec_cHL(2); const double log_dm_cHL_ac_HI=set_log_dm_agec_cHL(3); const double log_dm_cHL_ac_PH=set_log_dm_agec_cHL(4);
  const double log_dm_rHB_ac_LO=set_log_dm_agec_rHB(2); const double log_dm_rHB_ac_HI=set_log_dm_agec_rHB(3); const double log_dm_rHB_ac_PH=set_log_dm_agec_rHB(4);
  const double log_dm_rGN_ac_LO=set_log_dm_agec_rGN(2); const double log_dm_rGN_ac_HI=set_log_dm_agec_rGN(3); const double log_dm_rGN_ac_PH=set_log_dm_agec_rGN(4);
  const double log_dm_sCT_ac_LO=set_log_dm_agec_sCT(2); const double log_dm_sCT_ac_HI=set_log_dm_agec_sCT(3); const double log_dm_sCT_ac_PH=set_log_dm_agec_sCT(4);
  
  
  
  const double selpar_A50_cHL2_LO=set_selpar_A50_cHL2(2); const double selpar_A50_cHL2_HI=set_selpar_A50_cHL2(3); const double selpar_A50_cHL2_PH=set_selpar_A50_cHL2(4);
  const double selpar_slope_cHL2_LO=set_selpar_slope_cHL2(2); const double selpar_slope_cHL2_HI=set_selpar_slope_cHL2(3); const double selpar_slope_cHL2_PH=set_selpar_slope_cHL2(4);
  const double selpar_A50_cHL3_LO=set_selpar_A50_cHL3(2); const double selpar_A50_cHL3_HI=set_selpar_A50_cHL3(3); const double selpar_A50_cHL3_PH=set_selpar_A50_cHL3(4);
  const double selpar_slope_cHL3_LO=set_selpar_slope_cHL3(2); const double selpar_slope_cHL3_HI=set_selpar_slope_cHL3(3); const double selpar_slope_cHL3_PH=set_selpar_slope_cHL3(4);

  const double selpar_A50_cOT2_LO=set_selpar_A50_cOT2(2); const double selpar_A50_cOT2_HI=set_selpar_A50_cOT2(3); const double selpar_A50_cOT2_PH=set_selpar_A50_cOT2(4);
  const double selpar_A50_cOT3_LO=set_selpar_A50_cOT3(2); const double selpar_A50_cOT3_HI=set_selpar_A50_cOT3(3); const double selpar_A50_cOT3_PH=set_selpar_A50_cOT3(4);
  const double selpar_slope_cOT2_LO=set_selpar_slope_cOT2(2); const double selpar_slope_cOT2_HI=set_selpar_slope_cOT2(3); const double selpar_slope_cOT2_PH=set_selpar_slope_cOT2(4);
  const double selpar_A502_cOT2_LO=set_selpar_A502_cOT2(2); const double selpar_A502_cOT2_HI=set_selpar_A502_cOT2(3); const double selpar_A502_cOT2_PH=set_selpar_A502_cOT2(4);
  const double selpar_slope2_cOT2_LO=set_selpar_slope2_cOT2(2); const double selpar_slope2_cOT2_HI=set_selpar_slope2_cOT2(3); const double selpar_slope2_cOT2_PH=set_selpar_slope2_cOT2(4);
  
  const double selpar_A50_rHB1_LO=set_selpar_A50_rHB1(2); const double selpar_A50_rHB1_HI=set_selpar_A50_rHB1(3); const double selpar_A50_rHB1_PH=set_selpar_A50_rHB1(4);
  const double selpar_slope_rHB1_LO=set_selpar_slope_rHB1(2); const double selpar_slope_rHB1_HI=set_selpar_slope_rHB1(3); const double selpar_slope_rHB1_PH=set_selpar_slope_rHB1(4);
  const double selpar_A50_rHB2_LO=set_selpar_A50_rHB2(2); const double selpar_A50_rHB2_HI=set_selpar_A50_rHB2(3); const double selpar_A50_rHB2_PH=set_selpar_A50_rHB2(4);
  const double selpar_slope_rHB2_LO=set_selpar_slope_rHB2(2); const double selpar_slope_rHB2_HI=set_selpar_slope_rHB2(3); const double selpar_slope_rHB2_PH=set_selpar_slope_rHB2(4);
  const double selpar_A50_rHB3_LO=set_selpar_A50_rHB3(2); const double selpar_A50_rHB3_HI=set_selpar_A50_rHB3(3); const double selpar_A50_rHB3_PH=set_selpar_A50_rHB3(4);
  const double selpar_slope_rHB3_LO=set_selpar_slope_rHB3(2); const double selpar_slope_rHB3_HI=set_selpar_slope_rHB3(3); const double selpar_slope_rHB3_PH=set_selpar_slope_rHB3(4);
  

  const double selpar_A50_rGN3_LO=set_selpar_A50_rGN3(2); const double selpar_A50_rGN3_HI=set_selpar_A50_rGN3(3); const double selpar_A50_rGN3_PH=set_selpar_A50_rGN3(4);
  const double selpar_slope_rGN3_LO=set_selpar_slope_rGN3(2); const double selpar_slope_rGN3_HI=set_selpar_slope_rGN3(3); const double selpar_slope_rGN3_PH=set_selpar_slope_rGN3(4);
  
  const double selpar_A50_sCT_LO=set_selpar_A50_sCT(2); const double selpar_A50_sCT_HI=set_selpar_A50_sCT(3); const double selpar_A50_sCT_PH=set_selpar_A50_sCT(4);
  const double selpar_slope_sCT_LO=set_selpar_slope_sCT(2); const double selpar_slope_sCT_HI=set_selpar_slope_sCT(3); const double selpar_slope_sCT_PH=set_selpar_slope_sCT(4);
  const double selpar_A502_sCT_LO=set_selpar_A502_sCT(2); const double selpar_A502_sCT_HI=set_selpar_A502_sCT(3); const double selpar_A502_sCT_PH=set_selpar_A502_sCT(4);
  const double selpar_slope2_sCT_LO=set_selpar_slope2_sCT(2); const double selpar_slope2_sCT_HI=set_selpar_slope2_sCT(3); const double selpar_slope2_sCT_PH=set_selpar_slope2_sCT(4);
  
  const double selpar_age1logit_D_LO=set_selpar_age1logit_D(2); const double selpar_age1logit_D_HI=set_selpar_age1logit_D(3); const double selpar_age1logit_D_PH=set_selpar_age1logit_D(4);
  
  const double log_q_cHL_LO=set_log_q_cpue_cHL(2); const double log_q_cHL_HI=set_log_q_cpue_cHL(3); const double log_q_cHL_PH=set_log_q_cpue_cHL(4);
  const double log_q_rHB_LO=set_log_q_cpue_rHB(2); const double log_q_rHB_HI=set_log_q_cpue_rHB(3); const double log_q_rHB_PH=set_log_q_cpue_rHB(4);
  
  const double log_q_sCT_LO=set_log_q_cpue_sCT(2); const double log_q_sCT_HI=set_log_q_cpue_sCT(3); const double log_q_sCT_PH=set_log_q_cpue_sCT(4);
  
  const double log_avg_F_cHL_LO=set_log_avg_F_L_cHL(2); const double log_avg_F_cHL_HI=set_log_avg_F_L_cHL(3); const double log_avg_F_cHL_PH=set_log_avg_F_L_cHL(4);
  const double log_avg_F_cOT_LO=set_log_avg_F_L_cOT(2); const double log_avg_F_cOT_HI=set_log_avg_F_L_cOT(3); const double log_avg_F_cOT_PH=set_log_avg_F_L_cOT(4);
  const double log_avg_F_rHB_LO=set_log_avg_F_L_rHB(2); const double log_avg_F_rHB_HI=set_log_avg_F_L_rHB(3); const double log_avg_F_rHB_PH=set_log_avg_F_L_rHB(4); 
  const double log_avg_F_rGN_LO=set_log_avg_F_L_rGN(2); const double log_avg_F_rGN_HI=set_log_avg_F_L_rGN(3); const double log_avg_F_rGN_PH=set_log_avg_F_L_rGN(4); 
  const double log_avg_F_cHL_D_LO=set_log_avg_F_D_cHL(2); const double log_avg_F_cHL_D_HI=set_log_avg_F_D_cHL(3); const double log_avg_F_cHL_D_PH=set_log_avg_F_D_cHL(4);
  const double log_avg_F_rHB_D_LO=set_log_avg_F_D_rHB(2); const double log_avg_F_rHB_D_HI=set_log_avg_F_D_rHB(3); const double log_avg_F_rHB_D_PH=set_log_avg_F_D_rHB(4); 
  const double log_avg_F_rGN_D_LO=set_log_avg_F_D_rGN(2); const double log_avg_F_rGN_D_HI=set_log_avg_F_D_rGN(3); const double log_avg_F_rGN_D_PH=set_log_avg_F_D_rGN(4); 
  
  
  const double log_F_dev_cHL_LO=set_log_dev_F_L_cHL(1); const double log_F_dev_cHL_HI=set_log_dev_F_L_cHL(2); const double log_F_dev_cHL_PH=set_log_dev_F_L_cHL(3);  
  const double log_F_dev_cOT_LO=set_log_dev_F_L_cOT(1); const double log_F_dev_cOT_HI=set_log_dev_F_L_cOT(2); const double log_F_dev_cOT_PH=set_log_dev_F_L_cOT(3);   
  const double log_F_dev_rHB_LO=set_log_dev_F_L_rHB(1); const double log_F_dev_rHB_HI=set_log_dev_F_L_rHB(2); const double log_F_dev_rHB_PH=set_log_dev_F_L_rHB(3);   
  const double log_F_dev_rGN_LO=set_log_dev_F_L_rGN(1); const double log_F_dev_rGN_HI=set_log_dev_F_L_rGN(2); const double log_F_dev_rGN_PH=set_log_dev_F_L_rGN(3);   
  
  const double log_F_dev_cHL_D_LO=set_log_dev_F_D_cHL(1); const double log_F_dev_cHL_D_HI=set_log_dev_F_D_cHL(2); const double log_F_dev_cHL_D_PH=set_log_dev_F_D_cHL(3);   
  const double log_F_dev_rHB_D_LO=set_log_dev_F_D_rHB(1); const double log_F_dev_rHB_D_HI=set_log_dev_F_D_rHB(2); const double log_F_dev_rHB_D_PH=set_log_dev_F_D_rHB(3);   
  const double log_F_dev_rGN_D_LO=set_log_dev_F_D_rGN(1); const double log_F_dev_rGN_D_HI=set_log_dev_F_D_rGN(2); const double log_F_dev_rGN_D_PH=set_log_dev_F_D_rGN(3);   
  
  const double log_RWq_LO=set_log_dev_RWq(1); const double log_RWq_HI=set_log_dev_RWq(2); const double log_RWq_PH=set_log_dev_RWq(3);  
  
  const double log_rec_dev_LO=set_log_dev_rec(1); const double log_rec_dev_HI=set_log_dev_rec(2); const double log_rec_dev_PH=set_log_dev_rec(3);          
  const double log_Nage_dev_LO=set_log_dev_Nage(1); const double log_Nage_dev_HI=set_log_dev_Nage(2); const double log_Nage_dev_PH=set_log_dev_Nage(3);          
  
 END_CALCS
 

  
  init_bounded_number Linf(Linf_LO,Linf_HI,Linf_PH);
  init_bounded_number K(K_LO,K_HI,K_PH);
  init_bounded_number t0(t0_LO,t0_HI,t0_PH);
  init_bounded_number len_cv_val(len_cv_LO,len_cv_HI,len_cv_PH);  
  vector Linf_out(1,8);
  vector K_out(1,8);
  vector t0_out(1,8);
  vector len_cv_val_out(1,8);

  vector meanlen_TL(1,nages);   
   
  vector wgt_g(1,nages);        
  vector wgt_kg(1,nages);       
  vector wgt_mt(1,nages);       
  vector wgt_klb(1,nages);      
  vector wgt_lb(1,nages);       
    
    
  matrix len_cHL_mm(styr,endyr,1,nages);          
  matrix wholewgt_cHL_klb(styr,endyr,1,nages);    
  matrix len_cOT_mm(styr,endyr,1,nages);          
  matrix wholewgt_cOT_klb(styr,endyr,1,nages);    
  matrix len_rHB_mm(styr,endyr,1,nages);          
  matrix wholewgt_rHB_klb(styr,endyr,1,nages);    
  matrix len_rGN_mm(styr,endyr,1,nages);          
  matrix wholewgt_rGN_klb(styr,endyr,1,nages);    

  matrix len_cHL_D_mm(styr,endyr,1,nages);          
  matrix wholewgt_cHL_D_klb(styr,endyr,1,nages);    
  matrix len_rHB_D_mm(styr,endyr,1,nages);          
  matrix wholewgt_rHB_D_klb(styr,endyr,1,nages);    
  matrix len_rGN_D_mm(styr,endyr,1,nages);          
  matrix wholewgt_rGN_D_klb(styr,endyr,1,nages);    
  
  matrix lenprob(1,nages,1,nlenbins);           
  number zscore_len;                            
  vector cprob_lenvec(1,nlenbins);              
  number zscore_lzero;                          
  number cprob_lzero;                           
    
  
  matrix lenprob_cHL(1,nages,1,nlenbins);     
  matrix lenprob_cOT(1,nages,1,nlenbins);     
  matrix lenprob_rHB(1,nages,1,nlenbins);     
  matrix lenprob_rHB_D(1,nages,1,nlenbins);   
  matrix lenprob_rGN(1,nages,1,nlenbins);     
  matrix lenprob_sCT(1,nages,1,nlenbins);    
  
  
  vector len_sd(1,nages);
  vector len_cv(1,nages); 


  matrix pred_cHL_lenc(1,nyr_lenc_cHL,1,nlenbins);
  matrix pred_cOT_lenc(1,nyr_lenc_cOT,1,nlenbins);
  matrix pred_rHB_lenc(1,nyr_lenc_rHB,1,nlenbins);  
  matrix pred_rHB_D_lenc(1,nyr_lenc_rHB_D,1,nlenbins);  
  matrix pred_rGN_lenc(1,nyr_lenc_rGN,1,nlenbins);  
  matrix pred_sCT_lenc(1,nyr_lenc_sCT,1,nlenbins);  
  
  matrix pred_cHL_agec(1,nyr_agec_cHL,1,nages_agec);
  matrix pred_cHL_agec_allages(1,nyr_agec_cHL,1,nages);
  matrix ErrorFree_cHL_agec(1,nyr_agec_cHL,1,nages);    
  matrix pred_rHB_agec(1,nyr_agec_rHB,1,nages_agec);
  matrix pred_rHB_agec_allages(1,nyr_agec_rHB,1,nages);  
  matrix ErrorFree_rHB_agec(1,nyr_agec_rHB,1,nages);
  matrix pred_rGN_agec(1,nyr_agec_rGN,1,nages_agec);  
  matrix pred_rGN_agec_allages(1,nyr_agec_rGN,1,nages);  
  matrix ErrorFree_rGN_agec(1,nyr_agec_rGN,1,nages);  
  matrix pred_sCT_agec(1,nyr_agec_sCT,1,nages_agec);  
  matrix pred_sCT_agec_allages(1,nyr_agec_sCT,1,nages);  
  matrix ErrorFree_sCT_agec(1,nyr_agec_sCT,1,nages);  
 

  vector nsamp_cHL_lenc_allyr(styr,endyr);
  vector nsamp_cOT_lenc_allyr(styr,endyr);
  vector nsamp_rHB_lenc_allyr(styr,endyr);
  vector nsamp_rHB_D_lenc_allyr(styr,endyr);  
  vector nsamp_rGN_lenc_allyr(styr,endyr);
  vector nsamp_sCT_lenc_allyr(styr,endyr);
  
  vector nsamp_cHL_agec_allyr(styr,endyr);
  vector nsamp_rHB_agec_allyr(styr,endyr);
  vector nsamp_rGN_agec_allyr(styr,endyr);  
  vector nsamp_sCT_agec_allyr(styr,endyr);  
  

  vector nfish_cHL_lenc_allyr(styr,endyr);
  vector nfish_cOT_lenc_allyr(styr,endyr);
  vector nfish_rHB_lenc_allyr(styr,endyr);
  vector nfish_rHB_D_lenc_allyr(styr,endyr);  
  vector nfish_rGN_lenc_allyr(styr,endyr);
  vector nfish_sCT_lenc_allyr(styr,endyr);
  
  vector nfish_cHL_agec_allyr(styr,endyr);
  vector nfish_rHB_agec_allyr(styr,endyr);
  vector nfish_rGN_agec_allyr(styr,endyr);  
  vector nfish_sCT_agec_allyr(styr,endyr);  


  vector neff_cHL_lenc_allyr(styr,endyr);
  vector neff_cOT_lenc_allyr(styr,endyr);
  vector neff_rHB_lenc_allyr(styr,endyr);
  vector neff_rHB_D_lenc_allyr(styr,endyr);  
  vector neff_rGN_lenc_allyr(styr,endyr);
  vector neff_sCT_lenc_allyr(styr,endyr);
  
  vector neff_cHL_agec_allyr(styr,endyr);
  vector neff_rHB_agec_allyr(styr,endyr);
  vector neff_rGN_agec_allyr(styr,endyr);  
  vector neff_sCT_agec_allyr(styr,endyr);  


  matrix N(styr,endyr+1,1,nages);           
  matrix N_mdyr(styr,endyr,1,nages);        
  matrix N_spawn(styr,endyr,1,nages);       
  init_bounded_vector log_dev_Nage(2,nages,log_Nage_dev_LO,log_Nage_dev_HI,log_Nage_dev_PH);
  vector log_Nage_dev_output(1,nages);      
  matrix B(styr,endyr+1,1,nages);           
  vector totB(styr,endyr+1);                
  vector totN(styr,endyr+1);                
  vector SSB(styr,endyr);                   
  vector rec(styr,endyr+1);                 
  vector prop_f(1,nages);
  vector prop_m(1,nages);
  vector maturity_f(1,nages);
  vector maturity_m(1,nages);
  vector reprod(1,nages);
 

  init_bounded_number log_R0(log_R0_LO,log_R0_HI,log_R0_PH);        
  vector log_R0_out(1,8);
  number R0;                                  
  init_bounded_number steep(steep_LO,steep_HI,steep_PH); 
  vector steep_out(1,8);
  init_bounded_number rec_sigma(rec_sigma_LO,rec_sigma_HI,rec_sigma_PH);  
  vector rec_sigma_out(1,8);
  init_bounded_number R_autocorr(R_autocorr_LO,R_autocorr_HI,R_autocorr_PH);  
  vector R_autocorr_out(1,8);

  number rec_sigma_sq;                        
  number rec_logL_add;                        
 
  init_bounded_dev_vector log_dev_rec(styr_rec_dev,endyr_rec_dev,log_rec_dev_LO,log_rec_dev_HI,log_rec_dev_PH);
  vector log_rec_dev_output(styr,endyr+1);             
  vector log_rec_dev_out(styr_rec_dev,endyr_rec_dev);  
  
  number var_rec_dev;                           
  number sigma_rec_dev;                         
  number BiasCor;                               
  number S0;                                    
  number B0;                                    
  number R1;                                    
  number R_virgin;                              
  vector SdS0(styr,endyr);                      

  
  init_bounded_number log_dm_lenc_cHL(log_dm_cHL_lc_LO,log_dm_cHL_lc_HI,log_dm_cHL_lc_PH);
  init_bounded_number log_dm_lenc_cOT(log_dm_cOT_lc_LO,log_dm_cOT_lc_HI,log_dm_cOT_lc_PH);
  init_bounded_number log_dm_lenc_rHB(log_dm_rHB_lc_LO,log_dm_rHB_lc_HI,log_dm_rHB_lc_PH);
  init_bounded_number log_dm_lenc_rGN(log_dm_rGN_lc_LO,log_dm_rGN_lc_HI,log_dm_rGN_lc_PH);
  init_bounded_number log_dm_lenc_rHB_D(log_dm_rHB_D_lc_LO,log_dm_rHB_D_lc_HI,log_dm_rHB_D_lc_PH);
  init_bounded_number log_dm_lenc_sCT(log_dm_sCT_lc_LO,log_dm_sCT_lc_HI,log_dm_sCT_lc_PH);
  init_bounded_number log_dm_agec_cHL(log_dm_cHL_ac_LO,log_dm_cHL_ac_HI,log_dm_cHL_ac_PH);
  init_bounded_number log_dm_agec_rHB(log_dm_rHB_ac_LO,log_dm_rHB_ac_HI,log_dm_rHB_ac_PH);
  init_bounded_number log_dm_agec_rGN(log_dm_rGN_ac_LO,log_dm_rGN_ac_HI,log_dm_rGN_ac_PH);
  init_bounded_number log_dm_agec_sCT(log_dm_sCT_ac_LO,log_dm_sCT_ac_HI,log_dm_sCT_ac_PH);
  vector log_dm_cHL_lc_out(1,8);
  vector log_dm_cOT_lc_out(1,8);
  vector log_dm_rHB_lc_out(1,8);
  vector log_dm_rGN_lc_out(1,8);
  vector log_dm_rHB_D_lc_out(1,8);
  vector log_dm_sCT_lc_out(1,8);
  vector log_dm_cHL_ac_out(1,8);
  vector log_dm_rHB_ac_out(1,8);
  vector log_dm_rGN_ac_out(1,8);
  vector log_dm_sCT_ac_out(1,8);





  matrix sel_cHL(styr,endyr,1,nages);
  vector sel_cHL_block1(1,nages);
  vector sel_cHL_block2(1,nages);
  vector sel_cHL_block3(1,nages);
  
  init_bounded_number selpar_A50_cHL2(selpar_A50_cHL2_LO,selpar_A50_cHL2_HI,selpar_A50_cHL2_PH);
  init_bounded_number selpar_slope_cHL2(selpar_slope_cHL2_LO,selpar_slope_cHL2_HI,selpar_slope_cHL2_PH);
  init_bounded_number selpar_A50_cHL3(selpar_A50_cHL3_LO,selpar_A50_cHL3_HI,selpar_A50_cHL3_PH);
  init_bounded_number selpar_slope_cHL3(selpar_slope_cHL3_LO,selpar_slope_cHL3_HI,selpar_slope_cHL3_PH);
  
  vector selpar_A50_cHL2_out(1,8);
  vector selpar_slope_cHL2_out(1,8);
  vector selpar_A50_cHL3_out(1,8);
  vector selpar_slope_cHL3_out(1,8);

  
  matrix sel_cOT(styr,endyr,1,nages);
  vector sel_cOT_block1(1,nages);
  vector sel_cOT_block2(1,nages);
  vector sel_cOT_block3(1,nages);

  init_bounded_number selpar_A50_cOT2(selpar_A50_cOT2_LO,selpar_A50_cOT2_HI,selpar_A50_cOT2_PH);
  init_bounded_number selpar_A50_cOT3(selpar_A50_cOT3_LO,selpar_A50_cOT3_HI,selpar_A50_cOT3_PH);
  init_bounded_number selpar_slope_cOT2(selpar_slope_cOT2_LO,selpar_slope_cOT2_HI,selpar_slope_cOT2_PH);
  init_bounded_number selpar_A502_cOT2(selpar_A502_cOT2_LO,selpar_A502_cOT2_HI,selpar_A502_cOT2_PH);
  init_bounded_number selpar_slope2_cOT2(selpar_slope2_cOT2_LO,selpar_slope2_cOT2_HI,selpar_slope2_cOT2_PH);
  
  vector selpar_A50_cOT2_out(1,8);
  vector selpar_A50_cOT3_out(1,8);
  vector selpar_slope_cOT2_out(1,8);
  vector selpar_A502_cOT2_out(1,8);
  vector selpar_slope2_cOT2_out(1,8);
  
  
  matrix sel_rHB(styr,endyr,1,nages);
  vector sel_rHB_block1(1,nages);
  vector sel_rHB_block2(1,nages);
  vector sel_rHB_block3(1,nages);
  
  init_bounded_number selpar_A50_rHB1(selpar_A50_rHB1_LO,selpar_A50_rHB1_HI,selpar_A50_rHB1_PH);
  init_bounded_number selpar_slope_rHB1(selpar_slope_rHB1_LO,selpar_slope_rHB1_HI,selpar_slope_rHB1_PH);
  init_bounded_number selpar_A50_rHB2(selpar_A50_rHB2_LO,selpar_A50_rHB2_HI,selpar_A50_rHB2_PH);
  init_bounded_number selpar_slope_rHB2(selpar_slope_rHB2_LO,selpar_slope_rHB2_HI,selpar_slope_rHB2_PH);    
  init_bounded_number selpar_A50_rHB3(selpar_A50_rHB3_LO,selpar_A50_rHB3_HI,selpar_A50_rHB3_PH);
  init_bounded_number selpar_slope_rHB3(selpar_slope_rHB3_LO,selpar_slope_rHB3_HI,selpar_slope_rHB3_PH);
  
  vector selpar_A50_rHB1_out(1,8);
  vector selpar_slope_rHB1_out(1,8);
  vector selpar_A50_rHB2_out(1,8);
  vector selpar_slope_rHB2_out(1,8);
  vector selpar_A50_rHB3_out(1,8);
  vector selpar_slope_rHB3_out(1,8);

  
  matrix sel_rGN(styr,endyr,1,nages);
  vector sel_rGN_block1(1,nages);
  vector sel_rGN_block2(1,nages);
  vector sel_rGN_block3(1,nages);
 
  init_bounded_number selpar_A50_rGN3(selpar_A50_rGN3_LO,selpar_A50_rGN3_HI,selpar_A50_rGN3_PH);
  init_bounded_number selpar_slope_rGN3(selpar_slope_rGN3_LO,selpar_slope_rGN3_HI,selpar_slope_rGN3_PH);

  vector selpar_A50_rGN3_out(1,8);
  vector selpar_slope_rGN3_out(1,8);

  

  matrix sel_D(styr,endyr,1,nages);
  vector sel_D_block1(1,nages);
  vector sel_D_block2(1,nages);
  vector sel_D_block3(1,nages);
  vector prob_belowsizelim_block2(1,nages);
  vector prob_belowsizelim_block3(1,nages);
  number zscore_lsizelim1;
  number zscore_lsizelim2;
  number cprob_lsizelim1;
  number cprob_lsizelim2;

	
  init_bounded_number selpar_age1logit_D(selpar_age1logit_D_LO,selpar_age1logit_D_HI,selpar_age1logit_D_PH);
  number selpar_age1_D; 
  vector selpar_age1logit_D_out(1,8);
     

  matrix sel_sCT(styr,endyr,1,nages);
  vector sel_sCT_vec(1,nages);
  
  init_bounded_number selpar_A50_sCT(selpar_A50_sCT_LO,selpar_A50_sCT_HI,selpar_A50_sCT_PH);
  init_bounded_number selpar_slope_sCT(selpar_slope_sCT_LO,selpar_slope_sCT_HI,selpar_slope_sCT_PH);
  init_bounded_number selpar_A502_sCT(selpar_A502_sCT_LO,selpar_A502_sCT_HI,selpar_A502_sCT_PH);
  init_bounded_number selpar_slope2_sCT(selpar_slope2_sCT_LO,selpar_slope2_sCT_HI,selpar_slope2_sCT_PH);
  
  vector selpar_A50_sCT_out(1,8);
  vector selpar_slope_sCT_out(1,8);
  vector selpar_A502_sCT_out(1,8);
  vector selpar_slope2_sCT_out(1,8);
       

  
  vector sel_wgted_L(1,nages);  
  vector sel_wgted_D(1,nages);  
  vector sel_wgted_tot(1,nages);



  vector pred_cHL_cpue(styr_cpue_cHL,endyr_cpue_cHL);                   
  matrix N_cHL(styr_cpue_cHL,endyr_cpue_cHL,1,nages);                   
  vector pred_rHB_cpue(styr_cpue_rHB,endyr_cpue_rHB);                   
  matrix N_rHB(styr_cpue_rHB,endyr_cpue_rHB,1,nages);                   
  
  
  vector pred_sCT_cpue(styr_cpue_sCT,endyr_cpue_sCT);                
  matrix N_sCT(styr_cpue_sCT,endyr_cpue_sCT,1,nages);                
  
  

  init_bounded_number log_q_cpue_cHL(log_q_cHL_LO,log_q_cHL_HI,log_q_cHL_PH);
  init_bounded_number log_q_cpue_rHB(log_q_rHB_LO,log_q_rHB_HI,log_q_rHB_PH);
  
  init_bounded_number log_q_cpue_sCT(log_q_sCT_LO,log_q_sCT_HI,log_q_sCT_PH);  
  vector log_q_cHL_out(1,8);
  vector log_q_rHB_out(1,8);
  
  vector log_q_sCT_out(1,8);  
  
  number q_rate;
  vector q_rate_fcn_cHL(styr_cpue_cHL,endyr_cpue_cHL);         
  vector q_rate_fcn_rHB(styr_cpue_rHB,endyr_cpue_rHB);         
  
  

  number q_DD_beta;
  vector q_DD_fcn(styr,endyr);    
  number B0_q_DD;                 
  vector B_q_DD(styr,endyr);      


 
 
 init_bounded_vector q_RW_log_dev_cHL(styr_cpue_cHL,endyr_cpue_cHL-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
 init_bounded_vector q_RW_log_dev_rHB(styr_cpue_rHB,endyr_cpue_rHB-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
 
 vector q_RW_log_dev_sCT(styr_cpue_sCT,endyr_cpue_sCT-1); 


 vector q_cHL(styr_cpue_cHL,endyr_cpue_cHL); 
 vector q_rHB(styr_cpue_rHB,endyr_cpue_rHB); 
 
 vector q_sCT(styr_cpue_sCT,endyr_cpue_sCT); 



  matrix L_cHL_num(styr,endyr,1,nages);   
  matrix L_cHL_klb(styr,endyr,1,nages);   
  vector pred_cHL_L_knum(styr,endyr);     
  vector pred_cHL_L_klb(styr,endyr);      

  matrix L_cOT_num(styr,endyr,1,nages);   
  matrix L_cOT_klb(styr,endyr,1,nages);   
  vector pred_cOT_L_knum(styr,endyr);     
  vector pred_cOT_L_klb(styr,endyr);      

  matrix L_rHB_num(styr,endyr,1,nages);   
  matrix L_rHB_klb(styr,endyr,1,nages);   
  vector pred_rHB_L_knum(styr,endyr);     
  vector pred_rHB_L_klb(styr,endyr);      

  matrix L_rGN_num(styr,endyr,1,nages);   
  matrix L_rGN_klb(styr,endyr,1,nages);   
  vector pred_rGN_L_knum(styr,endyr);     
  vector pred_rGN_L_klb(styr,endyr);      

  matrix L_total_num(styr,endyr,1,nages);
  matrix L_total_klb(styr,endyr,1,nages);
  vector L_total_knum_yr(styr,endyr);    
  vector L_total_klb_yr(styr,endyr);     
  

  matrix D_cHL_num(styr,endyr,1,nages);            
  matrix D_cHL_klb(styr,endyr,1,nages);             
  
  vector obs_cHL_D(styr_D_cHL,endyr_D_cHL);           
  
  vector pred_cHL_D_knum(styr,endyr);     
  vector pred_cHL_D_klb(styr,endyr);      
  
  matrix D_rHB_num(styr,endyr,1,nages);             
  matrix D_rHB_klb(styr,endyr,1,nages);             
  
  vector obs_rHB_D(styr_D_rHB,endyr_D_rHB);           
  
  vector pred_rHB_D_knum(styr,endyr);     
  vector pred_rHB_D_klb(styr,endyr);      
  
  matrix D_rGN_num(styr,endyr,1,nages);             
  matrix D_rGN_klb(styr,endyr,1,nages);             
  
  vector obs_rGN_D(styr_D_rGN,endyr_D_rGN);           
  
  vector pred_rGN_D_knum(styr,endyr);     
  vector pred_rGN_D_klb(styr,endyr);      
  
  matrix D_total_num(styr,endyr,1,nages);          
  matrix D_total_klb(styr,endyr,1,nages);          
  vector D_total_knum_yr(styr,endyr);              
  vector D_total_klb_yr(styr,endyr);               

  number Dmort_cHL;
  number Dmort_rHB;
  number Dmort_rGN;  


  number F_cHL_prop;       
  number F_cOT_prop;       
  number F_rHB_prop;       
  number F_rGN_prop;       
  number F_cHL_D_prop;     
  number F_rHB_D_prop;     
  number F_rGN_D_prop;     
    
  number F_temp_sum;      

  vector F_end(1,nages);
  vector F_end_L(1,nages);  
  vector F_end_D(1,nages);    
  number F_end_apex;
  
  number SSB_msy_out;           
  number F_msy_out;             
  number msy_klb_out;           
  number msy_knum_out;          
  number D_msy_klb_out;         
  number D_msy_knum_out;        
  number B_msy_out;             
  number R_msy_out;             
  number spr_msy_out;           

  number F20_dum;				
  number F30_dum;				
  number F40_dum;				
  number F20_out;              	
  number F30_out;              	
  number F40_out;              	
  number SSB_F30_out;  
  number B_F30_out;
  number R_F30_out;
  number L_F30_knum_out;
  number L_F30_klb_out;
  number D_F30_knum_out;
  number D_F30_klb_out;    
  number rec_mean;  			

  vector N_age_msy(1,nages);         
  vector N_age_msy_spawn(1,nages);   
  vector L_age_msy(1,nages);         
  vector D_age_msy(1,nages);         
  vector Z_age_msy(1,nages);         
  vector F_L_age_msy(1,nages);       
  vector F_D_age_msy(1,nages);       
  vector F_msy(1,n_iter_msy);        
  vector spr_msy(1,n_iter_msy);      
  vector R_eq(1,n_iter_msy);         
  vector L_eq_klb(1,n_iter_msy);     
  vector L_eq_knum(1,n_iter_msy);    
  vector D_eq_klb(1,n_iter_msy);     
  vector D_eq_knum(1,n_iter_msy);    
  vector SSB_eq(1,n_iter_msy);       
  vector B_eq(1,n_iter_msy);         
  
  vector FdF_msy(styr,endyr);
  vector FdF30(styr,endyr);
  vector SdSSB_msy(styr,endyr);	 
  number SdSSB_msy_end;
  number FdF_msy_end;
  number FdF_msy_end_mean;           
  vector SdSSB_F30(styr,endyr);	 
  vector Sdmsst_F30(styr,endyr);	 
  number SdSSB_F30_end;
  number Sdmsst_F30_end;
  number FdF30_end_mean;             
  number Fend_mean_temp;			 
  number Fend_mean;					 
  vector L_age_F30(1,nages);         
  vector D_age_F30(1,nages);         
  
  vector wgt_wgted_L_klb(1,nages);   
  vector wgt_wgted_D_klb(1,nages);   
  number wgt_wgted_L_denom;          
  number wgt_wgted_D_denom;          

  number iter_inc_msy;               
  


  vector M(1,nages);                         
  init_bounded_number M_constant(M_constant_LO,M_constant_HI,M_constant_PH);   
  vector M_constant_out(1,8);
  number smsy2msstM;                         
  number smsy2msst75;                        
  
  matrix F(styr,endyr,1,nages);
  vector Fsum(styr,endyr);                   
  vector Fapex(styr,endyr);                  
  matrix Z(styr,endyr,1,nages);

  init_bounded_number log_avg_F_L_cHL(log_avg_F_cHL_LO,log_avg_F_cHL_HI,log_avg_F_cHL_PH);
  vector log_avg_F_cHL_out(1,8);  
  init_bounded_dev_vector log_dev_F_L_cHL(styr_L_cHL,endyr_L_cHL,log_F_dev_cHL_LO,log_F_dev_cHL_HI,log_F_dev_cHL_PH);
  vector log_F_dev_cHL_out(styr_L_cHL,endyr_L_cHL);
  matrix F_cHL(styr,endyr,1,nages);
  vector F_cHL_out(styr,endyr);        
  number log_F_dev_init_cHL;
  number log_F_dev_end_cHL; 

  init_bounded_number log_avg_F_L_cOT(log_avg_F_cOT_LO,log_avg_F_cOT_HI,log_avg_F_cOT_PH);
  vector log_avg_F_cOT_out(1,8);  
  init_bounded_dev_vector log_dev_F_L_cOT(styr_L_cOT,endyr_L_cOT,log_F_dev_cOT_LO,log_F_dev_cOT_HI,log_F_dev_cOT_PH);
  vector log_F_dev_cOT_out(styr_L_cOT,endyr_L_cOT);
  matrix F_cOT(styr,endyr,1,nages);
  vector F_cOT_out(styr,endyr);        
  number log_F_dev_init_cOT;
  number log_F_dev_end_cOT; 

  init_bounded_number log_avg_F_L_rHB(log_avg_F_rHB_LO,log_avg_F_rHB_HI,log_avg_F_rHB_PH);
  vector log_avg_F_rHB_out(1,8); 
  init_bounded_dev_vector log_dev_F_L_rHB(styr_L_rHB,endyr_L_rHB,log_F_dev_rHB_LO,log_F_dev_rHB_HI,log_F_dev_rHB_PH);    
  vector log_F_dev_rHB_out(styr_L_rHB,endyr_L_rHB);
  matrix F_rHB(styr,endyr,1,nages);
  vector F_rHB_out(styr,endyr); 
  number log_F_dev_init_rHB;    
  number log_F_dev_end_rHB;  

  init_bounded_number log_avg_F_L_rGN(log_avg_F_rGN_LO,log_avg_F_rGN_HI,log_avg_F_rGN_PH);
  vector log_avg_F_rGN_out(1,8); 
  init_bounded_dev_vector log_dev_F_L_rGN(styr_L_rGN,endyr_L_rGN,log_F_dev_rGN_LO,log_F_dev_rGN_HI,log_F_dev_rGN_PH);    
  vector log_F_dev_rGN_out(styr_L_rGN,endyr_L_rGN);
  matrix F_rGN(styr,endyr,1,nages);
  vector F_rGN_out(styr,endyr); 
  number log_F_dev_init_rGN;    
  number log_F_dev_end_rGN;  

  init_bounded_number log_avg_F_D_cHL(log_avg_F_cHL_D_LO,log_avg_F_cHL_D_HI,log_avg_F_cHL_D_PH);
  vector log_avg_F_cHL_D_out(1,8);  
  init_bounded_dev_vector log_dev_F_D_cHL(styr_D_cHL,endyr_D_cHL,log_F_dev_cHL_D_LO,log_F_dev_cHL_D_HI,log_F_dev_cHL_D_PH);
  vector log_F_dev_cHL_D_out(styr_D_cHL,endyr_D_cHL);
  matrix F_cHL_D(styr,endyr,1,nages);
  vector F_cHL_D_out(styr,endyr);        
  number log_F_dev_end_cHL_D; 
  number log_F_avgdev_cHL_D;

  init_bounded_number log_avg_F_D_rHB(log_avg_F_rHB_D_LO,log_avg_F_rHB_D_HI,log_avg_F_rHB_D_PH);
  vector log_avg_F_rHB_D_out(1,8);  
  init_bounded_dev_vector log_dev_F_D_rHB(styr_D_rHB,endyr_D_rHB,log_F_dev_rHB_D_LO,log_F_dev_rHB_D_HI,log_F_dev_rHB_D_PH);
  vector log_F_dev_rHB_D_out(styr_D_rHB,endyr_D_rHB);
  matrix F_rHB_D(styr,endyr,1,nages);
  vector F_rHB_D_out(styr,endyr);        
  number log_F_dev_end_rHB_D; 
  number log_F_avgdev_rHB_D;
  
  init_bounded_number log_avg_F_D_rGN(log_avg_F_rGN_D_LO,log_avg_F_rGN_D_HI,log_avg_F_rGN_D_PH);
  vector log_avg_F_rGN_D_out(1,8);  
  init_bounded_dev_vector log_dev_F_D_rGN(styr_D_rGN,endyr_D_rGN,log_F_dev_rGN_D_LO,log_F_dev_rGN_D_HI,log_F_dev_rGN_D_PH);
  vector log_F_dev_rGN_D_out(styr_D_rGN,endyr_D_rGN);
  matrix F_rGN_D(styr,endyr,1,nages);
  vector F_rGN_D_out(styr,endyr);        
  number log_F_dev_end_rGN_D; 
  number log_F_avgdev_rGN_D;

  
  
  
  
  


  vector N_age_spr(1,nages);         
  vector N_age_spr_spawn(1,nages);   
  vector L_age_spr(1,nages);         
  vector Z_age_spr(1,nages);         
  vector spr_static(styr,endyr);     
  vector F_L_age_spr(1,nages);       
  vector F_D_age_spr(1,nages);       
  vector F_spr(1,n_iter_spr);        
  vector spr_spr(1,n_iter_spr);      
  vector spr_ratio(1,n_iter_spr);    
  vector L_spr(1,n_iter_spr);        

  vector N_spr_F0(1,nages);          
  vector N_bpr_F0(1,nages);          
  vector N_spr_initial(1,nages);     
  vector N_initial_eq(1,nages);      
  vector F_initial(1,nages);         
  vector Z_initial(1,nages);         
  number spr_initial;                
  number spr_F0;                     
  number bpr_F0;                     

  number iter_inc_spr;               


 
  number sdnr_lc_cHL;
  number sdnr_lc_cOT;
  number sdnr_lc_rHB;
  number sdnr_lc_rHB_D;   
  number sdnr_lc_rGN;  
  number sdnr_lc_sCT; 
 
  number sdnr_ac_cHL;
  number sdnr_ac_rHB;
  number sdnr_ac_rGN;
  number sdnr_ac_sCT;  
  
  number sdnr_I_cHL;
  number sdnr_I_rHB;
  
  number sdnr_I_sCT;  
    

  number w_L;
  number w_D;
  
  number w_cpue_cHL;
  number w_cpue_rHB;
  
  number w_cpue_sCT; 
   
  number w_lenc_cHL;
  number w_lenc_cOT; 
  number w_lenc_rHB;
  number w_lenc_rHB_D;  
  number w_lenc_rGN;
  number w_lenc_sCT;  
  
  number w_agec_cHL; 
  number w_agec_rHB;
  number w_agec_rGN;  
  number w_agec_sCT; 
  
  number w_Nage_init;  
  number w_rec;
  number w_rec_early;
  number w_rec_end;
  number w_fullF;  
  number w_Ftune;

  number f_cHL_L; 
  number f_cOT_L; 
  number f_rHB_L; 
  number f_rGN_L; 

  number f_cHL_D; 
  number f_rHB_D; 
  number f_rGN_D; 

  number f_cHL_cpue;
  number f_rHB_cpue;
  
  number f_sCT_cpue;  
 
  number f_rHB_RWq_cpue;
  number f_cHL_RWq_cpue;
  
  number f_cHL_lenc;
  number f_cOT_lenc;  
  number f_rHB_lenc;
  number f_rHB_D_lenc;  
  number f_rGN_lenc;  
  number f_sCT_lenc;  
  
  number f_cHL_agec;
  number f_rHB_agec;   
  number f_rGN_agec; 
  number f_sCT_agec; 
  

  number f_Nage_init;              
  number f_rec_dev;                
  number f_rec_dev_early;          
  number f_rec_dev_end;            
  number f_fullF_constraint;       
  number f_Ftune;                  
  number f_priors;                 
  
   
   objective_function_value fval;
   number fval_data;
   number grad_max;
  
 
   number denom;                   
   number numer;                   

  
  number F_reg_proj;						   
  vector F_proj(styr_proj,endyr_proj);         
  vector L_knum_proj(styr_proj,endyr_proj);    
  vector L_klb_proj(styr_proj,endyr_proj);     
  vector D_knum_proj(styr_proj,endyr_proj);    
  vector D_klb_proj(styr_proj,endyr_proj);     

  vector B_proj(styr_proj,endyr_proj);         
  vector SSB_proj(styr_proj,endyr_proj);       
  vector R_proj(styr_proj,endyr_proj);     	   
  vector FL_age_proj(1,nages);      		   
  vector FD_age_proj(1,nages);      		   
  
  matrix N_proj(styr_proj,endyr_proj,1,nages);           
  matrix N_spawn_proj(styr_proj,endyr_proj,1,nages);     
  matrix Z_proj(styr_proj,endyr_proj,1,nages);           
  matrix L_age_proj(styr_proj,endyr_proj,1,nages);       
  matrix D_age_proj(styr_proj,endyr_proj,1,nages);       
  





GLOBALS_SECTION
  #include "admodel.h"          
  #include "admb2r.cpp"         
  #include <time.h>
	time_t start,finish;
	long hour,minute,second;	
	double elapsed_time;
	

RUNTIME_SECTION
 maximum_function_evaluations 1000, 2000,3000, 5000, 10000, 10000, 10000;
 convergence_criteria 1e-2, 1e-2,1e-3, 1e-3, 1e-3, 1e-4, 1e-4;
 


PRELIMINARY_CALCS_SECTION


  Dmort_cHL=set_Dmort_cHL; Dmort_rHB=set_Dmort_rHB; Dmort_rGN=set_Dmort_rGN;

  for(iyear=styr_D_cHL; iyear<=endyr_D_cHL; iyear++)
	{obs_cHL_D(iyear)=Dmort_cHL*obs_released_cHL(iyear);
	}
	
  for(iyear=styr_D_rHB; iyear<=endyr_D_rHB; iyear++)
	{obs_rHB_D(iyear)=Dmort_rHB*obs_released_rHB(iyear);
	}
	
  for(iyear=styr_D_rGN; iyear<=endyr_D_rGN; iyear++)
	{obs_rGN_D(iyear)=Dmort_rGN*obs_released_rGN(iyear);
	}
 
 
  Linf=set_Linf(1);
  K=set_K(1);
  t0=set_t0(1);
  len_cv_val=set_len_cv(1);
  
  M=set_M; 
  M_constant=set_M_constant(1);
  smsy2msstM=1.0-M_constant;
  smsy2msst75=0.75;  
  
  log_R0=set_log_R0(1);
  steep=set_steep(1);
  R_autocorr=set_R_autocorr(1);
  rec_sigma=set_rec_sigma(1);
  
  log_dm_lenc_cHL=set_log_dm_lenc_cHL(1);
  log_dm_lenc_cOT=set_log_dm_lenc_cOT(1);
  log_dm_lenc_rHB=set_log_dm_lenc_rHB(1);
  log_dm_lenc_rGN=set_log_dm_lenc_rGN(1);
  log_dm_lenc_rHB_D=set_log_dm_lenc_rHB_D(1);
  log_dm_lenc_sCT=set_log_dm_lenc_sCT(1);
  log_dm_agec_cHL=set_log_dm_agec_cHL(1);
  log_dm_agec_rHB=set_log_dm_agec_rHB(1);
  log_dm_agec_rGN=set_log_dm_agec_rGN(1);
  log_dm_agec_sCT=set_log_dm_agec_sCT(1);
  
  log_q_cpue_cHL=set_log_q_cpue_cHL(1);
  log_q_cpue_rHB=set_log_q_cpue_rHB(1);
  
  log_q_cpue_sCT=set_log_q_cpue_sCT(1);
  
  q_rate=set_q_rate;
  q_rate_fcn_cHL=1.0;   
  q_rate_fcn_rHB=1.0;   
  
  q_DD_beta=set_q_DD_beta;
  q_DD_fcn=1.0;

  q_RW_log_dev_cHL.initialize(); 
  q_RW_log_dev_rHB.initialize(); 
  
  q_RW_log_dev_sCT.initialize();
  
   if (set_q_rate_phase<0 & q_rate!=0.0)
  {
    for (iyear=styr_cpue_cHL; iyear<=endyr_cpue_cHL; iyear++)
      {   if (iyear>styr_cpue_cHL & iyear <=2003) 
          {
             q_rate_fcn_cHL(iyear)=(1.0+(iyear-styr_cpue_cHL)*q_rate)*q_rate_fcn_cHL(styr_cpue_cHL);  
          }
          if (iyear>2003) {q_rate_fcn_cHL(iyear)=q_rate_fcn_cHL(iyear-1);} 
      }   
    for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++)
      {   if (iyear>styr_cpue_rHB & iyear <=2003) 
          {
             q_rate_fcn_rHB(iyear)=(1.0+(iyear-styr_cpue_rHB)*q_rate)*q_rate_fcn_rHB(styr_cpue_rHB);  
          }
          if (iyear>2003) {q_rate_fcn_rHB(iyear)=q_rate_fcn_rHB(iyear-1);} 
      }   
    
      
          
             
          
          
      
	  
  } 

  w_L=set_w_L;
  w_D=set_w_D;
  
  w_cpue_cHL=set_w_cpue_cHL;
  w_cpue_rHB=set_w_cpue_rHB;
  
  w_cpue_sCT=set_w_cpue_sCT;
  
  w_lenc_cHL=set_w_lenc_cHL;
  w_lenc_cOT=set_w_lenc_cOT;  
  w_lenc_rHB=set_w_lenc_rHB;
  w_lenc_rHB_D=set_w_lenc_rHB_D; 
  w_lenc_rGN=set_w_lenc_rGN;
  w_lenc_sCT=set_w_lenc_sCT;    
    
  w_agec_cHL=set_w_agec_cHL;
  w_agec_rHB=set_w_agec_rHB;
  w_agec_rGN=set_w_agec_rGN; 
  w_agec_sCT=set_w_agec_sCT;  
 
  
  w_Nage_init=set_w_Nage_init;
  w_rec=set_w_rec;
  w_rec_early=set_w_rec_early;
  w_rec_end=set_w_rec_end;
  w_fullF=set_w_fullF;
  w_Ftune=set_w_Ftune;

  log_avg_F_L_cHL=set_log_avg_F_L_cHL(1);
  log_avg_F_L_cOT=set_log_avg_F_L_cOT(1);  
  log_avg_F_L_rHB=set_log_avg_F_L_rHB(1); 
  log_avg_F_L_rGN=set_log_avg_F_L_rGN(1); 
  log_avg_F_D_cHL=set_log_avg_F_D_cHL(1);
  log_avg_F_D_rHB=set_log_avg_F_D_rHB(1); 
  log_avg_F_D_rGN=set_log_avg_F_D_rGN(1); 
    
  log_dev_F_L_cHL=set_log_dev_vals_F_L__cHL;
  log_dev_F_L_cOT=set_log_dev_vals_F_L__cOT;
  log_dev_F_L_rHB=set_log_dev_vals_F_L__rHB;
  log_dev_F_L_rGN=set_log_dev_vals_F_L__rGN;
  log_dev_F_D_cHL=set_log_dev_vals_F_D__cHvals;
  log_dev_F_D_rHB=set_log_dev_vals_F_D__HBvals;
  log_dev_F_D_rGN=set_log_dev_vals_F_D__GRvals;
 
  
  
  selpar_A50_cHL2=set_selpar_A50_cHL2(1);
  selpar_slope_cHL2=set_selpar_slope_cHL2(1);
  selpar_A50_cHL3=set_selpar_A50_cHL3(1);
  selpar_slope_cHL3=set_selpar_slope_cHL3(1);

  selpar_A50_cOT2=set_selpar_A50_cOT2(1);
  selpar_A50_cOT3=set_selpar_A50_cOT3(1);
  selpar_slope_cOT2=set_selpar_slope_cOT2(1);
  selpar_A502_cOT2=set_selpar_A502_cOT2(1);
  selpar_slope2_cOT2=set_selpar_slope2_cOT2(1);
  
  selpar_A50_rHB1=set_selpar_A50_rHB1(1);
  selpar_slope_rHB1=set_selpar_slope_rHB1(1);
  selpar_A50_rHB2=set_selpar_A50_rHB2(1);
  selpar_slope_rHB2=set_selpar_slope_rHB2(1);
  selpar_A50_rHB3=set_selpar_A50_rHB3(1);
  selpar_slope_rHB3=set_selpar_slope_rHB3(1);

  selpar_A50_rGN3=set_selpar_A50_rGN3(1);
  selpar_slope_rGN3=set_selpar_slope_rGN3(1);
  
  selpar_age1logit_D=set_selpar_age1logit_D(1);
    
  selpar_A50_sCT=set_selpar_A50_sCT(1);
  selpar_slope_sCT=set_selpar_slope_sCT(1);
  selpar_A502_sCT=set_selpar_A502_sCT(1);
  selpar_slope2_sCT=set_selpar_slope2_sCT(1);

 sqrt2pi=sqrt(2.*3.14159265);
 g2mt=0.000001;         
 g2kg=0.001;            
 mt2klb=2.20462;        
 mt2lb=mt2klb*1000.0;   
 g2klb=g2mt*mt2klb;     
 dzero=0.00001;         
 huge_number=1.0e+10;   
 
 SSB_msy_out=0.0;

 iter_inc_msy=max_F_spr_msy/(n_iter_msy-1);
 iter_inc_spr=max_F_spr_msy/(n_iter_spr-1); 

 maturity_f=obs_maturity_f;
 maturity_m=obs_maturity_m;
 
 prop_f=obs_prop_f;
 prop_m=1.0-obs_prop_f;
  



      nsamp_cHL_lenc_allyr=missing;
	  nsamp_cOT_lenc_allyr=missing;  
	  nsamp_rHB_lenc_allyr=missing;
      nsamp_rHB_D_lenc_allyr=missing;  
	  nsamp_rGN_lenc_allyr=missing;
	  nsamp_sCT_lenc_allyr=missing;	  
      nsamp_cHL_agec_allyr=missing;
      nsamp_rHB_agec_allyr=missing;
	  nsamp_rGN_agec_allyr=missing;  
	  nsamp_sCT_agec_allyr=missing;  	  
      
      nfish_cHL_lenc_allyr=missing;
	  nfish_cOT_lenc_allyr=missing;  
	  nfish_rHB_lenc_allyr=missing;
      nfish_rHB_D_lenc_allyr=missing;  
	  nfish_rGN_lenc_allyr=missing;
	  nfish_sCT_lenc_allyr=missing;	  
      nfish_cHL_agec_allyr=missing;
      nfish_rHB_agec_allyr=missing;
	  nfish_rGN_agec_allyr=missing;  
	  nfish_sCT_agec_allyr=missing;  	  
   
      for (iyear=1; iyear<=nyr_lenc_cHL; iyear++)
         {if (nsamp_lenc_cHL(iyear)>=minSS_lenc_cHL)
           {nsamp_cHL_lenc_allyr(yrs_lenc_cHL(iyear))=nsamp_lenc_cHL(iyear);
            nfish_cHL_lenc_allyr(yrs_lenc_cHL(iyear))=nfish_lenc_cHL(iyear);}}
      for (iyear=1; iyear<=nyr_lenc_cOT; iyear++)
         {if (nsamp_lenc_cOT(iyear)>=minSS_lenc_cOT)
           {nsamp_cOT_lenc_allyr(yrs_lenc_cOT(iyear))=nsamp_lenc_cOT(iyear);
            nfish_cOT_lenc_allyr(yrs_lenc_cOT(iyear))=nfish_lenc_cOT(iyear);}}			
	  for (iyear=1; iyear<=nyr_lenc_rHB; iyear++)                           
         {if (nsamp_lenc_rHB(iyear)>=minSS_lenc_rHB)
            {nsamp_rHB_lenc_allyr(yrs_lenc_rHB(iyear))=nsamp_lenc_rHB(iyear);
             nfish_rHB_lenc_allyr(yrs_lenc_rHB(iyear))=nfish_lenc_rHB(iyear);}}
      for (iyear=1; iyear<=nyr_lenc_rHB_D; iyear++)                           
         {if (nsamp_lenc_rHB_D(iyear)>=minSS_lenc_rHB_D)
            {nsamp_rHB_D_lenc_allyr(yrs_lenc_rHB_D(iyear))=nsamp_lenc_rHB_D(iyear);
             nfish_rHB_D_lenc_allyr(yrs_lenc_rHB_D(iyear))=nfish_lenc_rHB_D(iyear);}}
	  for (iyear=1; iyear<=nyr_lenc_rGN; iyear++)                           
         {if (nsamp_lenc_rGN(iyear)>=minSS_lenc_rGN)
            {nsamp_rGN_lenc_allyr(yrs_lenc_rGN(iyear))=nsamp_lenc_rGN(iyear);
             nfish_rGN_lenc_allyr(yrs_lenc_rGN(iyear))=nfish_lenc_rGN(iyear);}}
	  for (iyear=1; iyear<=nyr_lenc_sCT; iyear++)                           
         {if (nsamp_lenc_sCT(iyear)>=minSS_lenc_sCT)
            {nsamp_sCT_lenc_allyr(yrs_lenc_sCT(iyear))=nsamp_lenc_sCT(iyear);
             nfish_sCT_lenc_allyr(yrs_lenc_sCT(iyear))=nfish_lenc_sCT(iyear);}}

	  for (iyear=1; iyear<=nyr_agec_cHL; iyear++)
         {if (nsamp_agec_cHL(iyear)>=minSS_agec_cHL)
           {nsamp_cHL_agec_allyr(yrs_agec_cHL(iyear))=nsamp_agec_cHL(iyear);
            nfish_cHL_agec_allyr(yrs_agec_cHL(iyear))=nfish_agec_cHL(iyear);}}
      for (iyear=1; iyear<=nyr_agec_rHB; iyear++)
          {if (nsamp_agec_rHB(iyear)>=minSS_agec_rHB)
            {nsamp_rHB_agec_allyr(yrs_agec_rHB(iyear))=nsamp_agec_rHB(iyear);
             nfish_rHB_agec_allyr(yrs_agec_rHB(iyear))=nfish_agec_rHB(iyear);}} 
      for (iyear=1; iyear<=nyr_agec_rGN; iyear++)  
         {if (nsamp_agec_rGN(iyear)>=minSS_agec_rGN)
           {nsamp_rGN_agec_allyr(yrs_agec_rGN(iyear))=nsamp_agec_rGN(iyear);
             nfish_rGN_agec_allyr(yrs_agec_rGN(iyear))=nfish_agec_rGN(iyear);}}  
	  for (iyear=1; iyear<=nyr_agec_sCT; iyear++)  
          {if (nsamp_agec_sCT(iyear)>=minSS_agec_sCT)
            {nsamp_sCT_agec_allyr(yrs_agec_sCT(iyear))=nsamp_agec_sCT(iyear);
             nfish_sCT_agec_allyr(yrs_agec_sCT(iyear))=nfish_agec_sCT(iyear);}} 

             

  F_msy(1)=0.0;  
  for (ff=2;ff<=n_iter_msy;ff++) {F_msy(ff)=F_msy(ff-1)+iter_inc_msy;}
  F_spr(1)=0.0;  
  for (ff=2;ff<=n_iter_spr;ff++) {F_spr(ff)=F_spr(ff-1)+iter_inc_spr;}



  F_cHL.initialize(); L_cHL_num.initialize();
  F_cOT.initialize(); L_cOT_num.initialize();
  F_rHB.initialize(); L_rHB_num.initialize();
  F_rGN.initialize(); L_rGN_num.initialize();
  F_cHL_D.initialize(); D_cHL_num.initialize();
  F_rHB_D.initialize(); D_rHB_num.initialize();
  F_rGN_D.initialize(); D_rGN_num.initialize();

  F_cHL_out.initialize();
  F_cOT_out.initialize();
  F_rHB_out.initialize();
  F_rGN_out.initialize();
  F_cHL_D_out.initialize();
  F_rHB_D_out.initialize();
  F_rGN_D_out.initialize();

  sel_cHL.initialize();
  sel_cOT.initialize();  
  sel_rHB.initialize();
  sel_rGN.initialize();
  sel_D.initialize();
  sel_sCT.initialize();

  sel_cHL_block1.initialize(); 
  sel_cHL_block2.initialize();
  sel_cHL_block3.initialize();  
  sel_cOT_block1.initialize();
  sel_cOT_block2.initialize();
  sel_cOT_block3.initialize();    
  sel_rHB_block1.initialize();
  sel_rHB_block2.initialize();
  sel_rHB_block3.initialize();
  sel_rGN_block1.initialize();
  sel_rGN_block2.initialize();
  sel_rGN_block3.initialize();
  sel_D_block1.initialize();   
  sel_D_block2.initialize();  
  sel_D_block3.initialize();
  sel_sCT_vec.initialize();
  prob_belowsizelim_block2.initialize();
  prob_belowsizelim_block3.initialize();
  
  log_rec_dev_output.initialize();  
  log_dev_rec=set_log_dev_vals_rec;
  log_Nage_dev_output.initialize();
  log_dev_Nage=set_log_dev_vals_Nage;
 
 


TOP_OF_MAIN_SECTION
  time(&start);
  arrmblsize=20000000;
  gradient_structure::set_MAX_NVAR_OFFSET(1600);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(2000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(10000);



PROCEDURE_SECTION
  
 
 
 
 
  get_length_weight_at_age(); 
  
  get_reprod();
  
  get_length_at_age_dist(); 
  
  get_weight_at_age_landings();
  
  get_spr_F0();
  
  get_selectivity(); 
  
  get_mortality(); 
  
  get_bias_corr(); 
  
  get_numbers_at_age(); 
  
  get_landings_numbers();
  
  get_landings_wgt();
  
  get_dead_discards(); 
  
  get_catchability_fcns(); 
  
  get_indices();
  
  get_length_comps();
  
  get_age_comps();
  
   evaluate_objective_function();
   
 
 
FUNCTION get_length_weight_at_age
	
    
    meanlen_TL=Linf*(1.0-mfexp(-K*(agebins-t0+0.5)));     
    wgt_kg=wgtpar_a*pow(meanlen_TL,wgtpar_b);             
    wgt_g=wgt_kg/g2kg;                                    
    wgt_mt=wgt_g*g2mt;                                    
    wgt_klb=mt2klb*wgt_mt;                                
    wgt_lb=mt2lb*wgt_mt;                                  
    
FUNCTION get_reprod 
 
    reprod=elem_prod((elem_prod(prop_f,maturity_f)+elem_prod(prop_m,maturity_m)),wgt_mt);
 
FUNCTION get_length_at_age_dist
  
    
	for (iage=1;iage<=nages;iage++)
   {len_cv(iage)=len_cv_val;
    len_sd(iage)=meanlen_TL(iage)*len_cv(iage);
	zscore_lzero=(0.0-meanlen_TL(iage))/len_sd(iage); 
	cprob_lzero=cumd_norm(zscore_lzero);
	    
    
	
    zscore_len=((lenbins(1)+0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage);
    cprob_lenvec(1)=cumd_norm(zscore_len);          
    lenprob(iage,1)=cprob_lenvec(1)-cprob_lzero;    
    
	zscore_lsizelim1=(sizelim1-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim1=cumd_norm(zscore_lsizelim1);                 
	prob_belowsizelim_block2(iage)=	cprob_lsizelim1-cprob_lzero; 

	zscore_lsizelim2=(sizelim2-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim2=cumd_norm(zscore_lsizelim2);                 
	prob_belowsizelim_block3(iage)=	cprob_lsizelim2-cprob_lzero; 
	
    
    
    for (ilen=2;ilen<nlenbins;ilen++)
      {
        zscore_len=((lenbins(ilen)+0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage); 
		cprob_lenvec(ilen)=cumd_norm(zscore_len);
        lenprob(iage,ilen)=cprob_lenvec(ilen)-cprob_lenvec(ilen-1);
      }
	
    
	
    zscore_len=((lenbins(nlenbins)-0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage); 
	lenprob(iage,nlenbins)=1.0-cumd_norm(zscore_len);
      lenprob(iage)=lenprob(iage)/(1.0-cprob_lzero);  
   }
   
  
  lenprob_cHL=lenprob;
  lenprob_cOT=lenprob; 
  lenprob_rHB=lenprob;
  lenprob_rHB_D=lenprob; 
  lenprob_rGN=lenprob;  
  lenprob_sCT=lenprob; 
  
  
FUNCTION get_weight_at_age_landings  
  
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    len_cHL_mm(iyear)=meanlen_TL;  
    wholewgt_cHL_klb(iyear)=wgt_klb; 
    len_cOT_mm(iyear)=meanlen_TL;  
    wholewgt_cOT_klb(iyear)=wgt_klb; 	
    len_rHB_mm(iyear)=meanlen_TL;
    wholewgt_rHB_klb(iyear)=wgt_klb;
    len_rGN_mm(iyear)=meanlen_TL;
    wholewgt_rGN_klb(iyear)=wgt_klb;

    len_cHL_D_mm(iyear)=meanlen_TL;  
    wholewgt_cHL_D_klb(iyear)=wgt_klb;
    len_rHB_D_mm(iyear)=meanlen_TL;
    wholewgt_rHB_D_klb(iyear)=wgt_klb;
    len_rGN_D_mm(iyear)=meanlen_TL;
    wholewgt_rGN_D_klb(iyear)=wgt_klb;
  }  
 
 
FUNCTION get_spr_F0
  
  N_spr_F0(1)=1.0*mfexp(-1.0*M(1)*spawn_time_frac); 
  N_bpr_F0(1)=1.0;      
  for (iage=2; iage<=nages; iage++)
  { N_spr_F0(iage)=N_spr_F0(iage-1)*mfexp(-1.0*(M(iage-1)*(1.0-spawn_time_frac) + M(iage)*spawn_time_frac)); 
    N_bpr_F0(iage)=N_bpr_F0(iage-1)*mfexp(-1.0*(M(iage-1)));    
  }
  N_spr_F0(nages)=N_spr_F0(nages)/(1.0-mfexp(-1.0*M(nages))); 
  N_bpr_F0(nages)=N_bpr_F0(nages)/(1.0-mfexp(-1.0*M(nages)));
  
  spr_F0=sum(elem_prod(N_spr_F0,reprod)); 
  bpr_F0=sum(elem_prod(N_bpr_F0,wgt_mt));    



FUNCTION get_selectivity

  sel_cHL_block2=logistic(agebins, selpar_A50_cHL2, selpar_slope_cHL2);
  sel_cHL_block3=logistic(agebins, selpar_A50_cHL3, selpar_slope_cHL3);
  sel_cHL_block1=sel_cHL_block2;
  sel_cOT_block2=logistic_double(agebins, selpar_A50_cOT2, selpar_slope_cOT2, selpar_A502_cOT2,selpar_slope2_cOT2);
  sel_cOT_block3=logistic_double(agebins, selpar_A50_cOT3, selpar_slope_cOT2, selpar_A502_cOT2,selpar_slope2_cOT2);
  sel_cOT_block1=sel_cOT_block2;
  sel_rHB_block1=logistic(agebins, selpar_A50_rHB1, selpar_slope_rHB1);
  sel_rHB_block2=logistic(agebins, selpar_A50_rHB2, selpar_slope_rHB2);
  sel_rHB_block3=logistic(agebins, selpar_A50_rHB3, selpar_slope_rHB3);
  sel_rGN_block1=sel_rHB_block1;
  sel_rGN_block2=sel_rHB_block2;
  sel_rGN_block3=logistic(agebins, selpar_A50_rGN3, selpar_slope_rGN3);
  sel_sCT_vec=logistic_double(agebins, selpar_A50_sCT, selpar_slope_sCT, selpar_A502_sCT,selpar_slope2_sCT);
	
  selpar_age1_D=mfexp(selpar_age1logit_D)/(1.0+mfexp(selpar_age1logit_D));
  sel_D_block2(1)=selpar_age1_D; sel_D_block2(2)=1.0; sel_D_block2(3,nages)=prob_belowsizelim_block2(3,nages);
  sel_D_block3(1)=selpar_age1_D; sel_D_block3(2)=1.0; sel_D_block3(3,nages)=prob_belowsizelim_block3(3,nages);
  sel_D_block1=sel_D_block2;
    
  
  for (iyear=styr; iyear<=endyr_selex_phase1; iyear++)
   {     
    sel_cHL(iyear)=sel_cHL_block1;
	sel_cOT(iyear)=sel_cOT_block1;
    sel_rHB(iyear)=sel_rHB_block1;
    sel_rGN(iyear)=sel_rGN_block1;
	sel_sCT(iyear)=sel_sCT_vec;
 	sel_D(iyear)=sel_D_block1;
   }
   
  
  
  for (iyear=(endyr_selex_phase1+1); iyear<=endyr_selex_phase2; iyear++)
   {
    sel_cHL(iyear)=sel_cHL_block2;
	sel_cOT(iyear)=sel_cOT_block2;
    sel_rHB(iyear)=sel_rHB_block2;
    sel_rGN(iyear)=sel_rGN_block2;
	sel_sCT(iyear)=sel_sCT_vec;
 	sel_D(iyear)=sel_D_block2;
   }
  
  
   for (iyear=(endyr_selex_phase2+1); iyear<=endyr; iyear++)
   {   
    sel_cHL(iyear)=sel_cHL_block3;
	sel_cOT(iyear)=sel_cOT_block3;
    sel_rHB(iyear)=sel_rHB_block3;
    sel_rGN(iyear)=sel_rGN_block3;
	sel_sCT(iyear)=sel_sCT_vec;
 	sel_D(iyear)=sel_D_block3;
   }  

   
FUNCTION get_mortality
  Fsum.initialize();
  Fapex.initialize();
  F.initialize();
  
  log_F_dev_init_cHL=sum(log_dev_F_L_cHL(styr_L_cHL,(styr_L_cHL+2)))/3.0;  
  log_F_dev_init_cOT=sum(log_dev_F_L_cOT(styr_L_cOT,(styr_L_cOT+2)))/3.0;  
  log_F_dev_init_rHB=sum(log_dev_F_L_rHB(styr_L_rHB,(styr_L_rHB+2)))/3.0;         
  log_F_dev_init_rGN=sum(log_dev_F_L_rGN(styr_L_rGN,(styr_L_rGN+2)))/3.0;         
  
  for (iyear=styr; iyear<=endyr; iyear++) 
  {
    if(iyear>=styr_L_cHL & iyear<=endyr_L_cHL) 
		{F_cHL_out(iyear)=mfexp(log_avg_F_L_cHL+log_dev_F_L_cHL(iyear));}     
    F_cHL(iyear)=sel_cHL(iyear)*F_cHL_out(iyear);
    Fsum(iyear)+=F_cHL_out(iyear);
    

	if(iyear>=styr_L_cOT & iyear<=endyr_L_cOT) 
		{F_cOT_out(iyear)=mfexp(log_avg_F_L_cOT+log_dev_F_L_cOT(iyear));}     
    F_cOT(iyear)=sel_cOT(iyear)*F_cOT_out(iyear);
    Fsum(iyear)+=F_cOT_out(iyear);
    

    if(iyear>=styr_L_rHB & iyear<=endyr_L_rHB) 
		{F_rHB_out(iyear)=mfexp(log_avg_F_L_rHB+log_dev_F_L_rHB(iyear));}     
    F_rHB(iyear)=sel_rHB(iyear)*F_rHB_out(iyear);
    Fsum(iyear)+=F_rHB_out(iyear);
    
   
    if(iyear>=styr_L_rGN & iyear<=endyr_L_rGN) 
		{F_rGN_out(iyear)=mfexp(log_avg_F_L_rGN+log_dev_F_L_rGN(iyear));}    
    if (iyear<styr_L_rGN)
		{F_rGN_out(iyear)=mfexp(log_avg_F_L_rGN+log_F_dev_init_rGN);}
	F_rGN(iyear)=sel_rGN(iyear)*F_rGN_out(iyear); 
    Fsum(iyear)+=F_rGN_out(iyear);

	
    log_F_avgdev_cHL_D=sum(log_dev_F_D_cHL(styr_D_cHL,endyr_D_cHL))/(endyr_D_cHL-styr_D_cHL+1.0);
    log_F_avgdev_rHB_D=sum(log_dev_F_D_rHB(styr_D_rHB,endyr_D_rHB))/(endyr_D_rHB-styr_D_rHB+1.0);   
    log_F_avgdev_rGN_D=sum(log_dev_F_D_rGN(styr_D_rGN,endyr_D_rGN))/(endyr_D_rGN-styr_D_rGN+1.0);	

    if(iyear>=styr_D_cHL & iyear<=endyr_D_cHL)
		{F_cHL_D_out(iyear)=mfexp(log_avg_F_D_cHL+log_dev_F_D_cHL(iyear));}    
    if(iyear>endyr_selex_phase1 & iyear<styr_D_cHL)
		{F_cHL_D_out(iyear)=mfexp(log_avg_F_D_cHL+log_F_avgdev_cHL_D);} 	
    F_cHL_D(iyear)=sel_D(iyear)*F_cHL_D_out(iyear);
    Fsum(iyear)+=F_cHL_D_out(iyear);
    
    if(iyear>=styr_D_rHB & iyear<=endyr_D_rHB)
		{F_rHB_D_out(iyear)=mfexp(log_avg_F_D_rHB+log_dev_F_D_rHB(iyear));}    
    if(iyear>endyr_selex_phase1 & iyear<styr_D_rHB)
      {F_rHB_D_out(iyear)=mfexp(log_avg_F_D_rHB+log_F_avgdev_rHB_D);} 	
    F_rHB_D(iyear)=sel_D(iyear)*F_rHB_D_out(iyear);
    Fsum(iyear)+=F_rHB_D_out(iyear);
 
    if(iyear>=styr_D_rGN & iyear<=endyr_D_rGN)
		{F_rGN_D_out(iyear)=mfexp(log_avg_F_D_rGN+log_dev_F_D_rGN(iyear));}  
	F_rGN_D(iyear)=sel_D(iyear)*F_rGN_D_out(iyear); 
	Fsum(iyear)+=F_rGN_D_out(iyear);
    
 
    
    F(iyear)=F_cHL(iyear);  
    F(iyear)+=F_cOT(iyear);
    F(iyear)+=F_rHB(iyear);
    F(iyear)+=F_rGN(iyear);
    F(iyear)+=F_cHL_D(iyear);
    F(iyear)+=F_rHB_D(iyear);
    F(iyear)+=F_rGN_D(iyear);
    
    Fapex(iyear)=max(F(iyear));
    Z(iyear)=M+F(iyear);
    
   }  

FUNCTION get_bias_corr
  var_rec_dev=norm2(log_dev_rec(styr_rec_dev,endyr_rec_dev)-
              sum(log_dev_rec(styr_rec_dev,endyr_rec_dev))/nyrs_rec)
              /(nyrs_rec-1.0);                           
  
  rec_sigma_sq=square(rec_sigma);
  if (set_BiasCor <= 0.0) {BiasCor=mfexp(rec_sigma_sq/2.0);}   
  else {BiasCor=set_BiasCor;}

FUNCTION get_numbers_at_age

  R0=mfexp(log_R0);
  S0=spr_F0*R0;
  R_virgin=SR_eq_func(R0, steep, spr_F0, spr_F0, BiasCor, SR_switch);
 
  B0=bpr_F0*R_virgin;   
  B0_q_DD=R_virgin*sum(elem_prod(N_bpr_F0(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages))); 
  			
  F_initial=sel_cHL(styr)*mfexp(log_avg_F_L_cHL+log_F_dev_init_cHL)+
            sel_cOT(styr)*mfexp(log_avg_F_L_cOT+log_F_dev_init_cOT)+
            sel_rHB(styr)*mfexp(log_avg_F_L_rHB+log_F_dev_init_rHB)+
            sel_rGN(styr)*mfexp(log_avg_F_L_rGN+log_F_dev_init_rGN);			
    
  Z_initial=M+F_initial;


  N_spr_initial(1)=1.0*mfexp(-1.0*Z_initial(1)*spawn_time_frac); 
  for (iage=2; iage<=nages; iage++)
    {
      N_spr_initial(iage)=N_spr_initial(iage-1)*
                   mfexp(-1.0*(Z_initial(iage-1)*(1.0-spawn_time_frac) + Z_initial(iage)*spawn_time_frac)); 
    }
  N_spr_initial(nages)=N_spr_initial(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); 
  spr_initial=sum(elem_prod(N_spr_initial,reprod));
  if (styr==styr_rec_dev) {R1=SR_eq_func(R0, steep, spr_F0, spr_initial, 1.0, SR_switch);} 
  else {R1=SR_eq_func(R0, steep, spr_F0, spr_initial, BiasCor, SR_switch);} 
  if(R1<10.0) {R1=10.0;} 


  N_initial_eq(1)=R1;
  for (iage=2; iage<=nages; iage++)
  {
    N_initial_eq(iage)=N_initial_eq(iage-1)*
        mfexp(-1.0*(Z_initial(iage-1)));    
  }
  
  N_initial_eq(nages)=N_initial_eq(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); 
  

  N(styr)(2,nages)=elem_prod(N_initial_eq(2,nages),mfexp(log_dev_Nage));
   
  if (styr==styr_rec_dev) {N(styr,1)=N_initial_eq(1)*mfexp(log_dev_rec(styr_rec_dev));}
  else {N(styr,1)=N_initial_eq(1);}
  
  N_mdyr(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*0.5))); 
  N_spawn(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*spawn_time_frac))); 

  SSB(styr)=sum(elem_prod(N_spawn(styr),reprod));
  B_q_DD(styr)=sum(elem_prod(N(styr)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
    

  for (iyear=styr; iyear<endyr; iyear++)
  {
    if(iyear<(styr_rec_dev-1)||iyear>(endyr_rec_dev-1)) 
    {
        N(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch);
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); 
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); 
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));
		B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));       
    }
    else   
    {
        N(iyear+1,1)=SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch)*mfexp(log_dev_rec(iyear+1));
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); 
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); 
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));
        B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
    }
  }
  
  
  N(endyr+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(endyr),SR_switch);
  N(endyr+1)(2,nages)=++elem_prod(N(endyr)(1,nages-1),(mfexp(-1.*Z(endyr)(1,nages-1))));
  N(endyr+1,nages)+=N(endyr,nages)*mfexp(-1.*Z(endyr,nages)); 

  
FUNCTION get_landings_numbers 
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      L_cHL_num(iyear,iage)=N(iyear,iage)*F_cHL(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_cOT_num(iyear,iage)=N(iyear,iage)*F_cOT(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);		
      L_rHB_num(iyear,iage)=N(iyear,iage)*F_rHB(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_rGN_num(iyear,iage)=N(iyear,iage)*F_rGN(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);        
    }          
    pred_cHL_L_knum(iyear)=sum(L_cHL_num(iyear))/1000.0;
	pred_cOT_L_knum(iyear)=sum(L_cOT_num(iyear))/1000.0;
    pred_rHB_L_knum(iyear)=sum(L_rHB_num(iyear))/1000.0;
    pred_rGN_L_knum(iyear)=sum(L_rGN_num(iyear))/1000.0;
  }

 
FUNCTION get_landings_wgt
  for (iyear=styr; iyear<=endyr; iyear++)
  {    
    L_cHL_klb(iyear)=elem_prod(L_cHL_num(iyear),wholewgt_cHL_klb(iyear));     
	L_cOT_klb(iyear)=elem_prod(L_cOT_num(iyear),wholewgt_cOT_klb(iyear));     
    L_rHB_klb(iyear)=elem_prod(L_rHB_num(iyear),wholewgt_rHB_klb(iyear));     
    L_rGN_klb(iyear)=elem_prod(L_rGN_num(iyear),wholewgt_rGN_klb(iyear));     
    
    pred_cHL_L_klb(iyear)=sum(L_cHL_klb(iyear));
	pred_cOT_L_klb(iyear)=sum(L_cOT_klb(iyear));
    pred_rHB_L_klb(iyear)=sum(L_rHB_klb(iyear));
    pred_rGN_L_klb(iyear)=sum(L_rGN_klb(iyear));    
  }
 
FUNCTION get_dead_discards
  

  
  for (iyear=styr; iyear<=endyr; iyear++)	  
  {
    for (iage=1; iage<=nages; iage++)
    {
      D_cHL_num(iyear,iage)=N(iyear,iage)*F_cHL_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
    }
    pred_cHL_D_knum(iyear)=sum(D_cHL_num(iyear))/1000.0;            
    pred_cHL_D_klb(iyear)=sum(elem_prod(D_cHL_num(iyear),wholewgt_cHL_D_klb(iyear))); 
  }

  
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      D_rHB_num(iyear,iage)=N(iyear,iage)*F_rHB_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
    }
    pred_rHB_D_knum(iyear)=sum(D_rHB_num(iyear))/1000.0;            
    pred_rHB_D_klb(iyear)=sum(elem_prod(D_rHB_num(iyear),wholewgt_rHB_D_klb(iyear))); 
  }

  
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      D_rGN_num(iyear,iage)=N(iyear,iage)*F_rGN_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
    }
    pred_rGN_D_knum(iyear)=sum(D_rGN_num(iyear))/1000.0;            
    pred_rGN_D_klb(iyear)=sum(elem_prod(D_rGN_num(iyear),wholewgt_rGN_D_klb(iyear))); 
  }
   
FUNCTION get_catchability_fcns    
 
  if (set_q_rate_phase>0.0)
  {

      for (iyear=styr_cpue_cHL; iyear<=endyr_cpue_cHL; iyear++)
      {   if (iyear>styr_cpue_cHL & iyear <=2003) 
          {
             q_rate_fcn_cHL(iyear)=(1.0+(iyear-styr_cpue_cHL)*q_rate)*q_rate_fcn_cHL(styr_cpue_cHL);  
          }
          if (iyear>2003) {q_rate_fcn_cHL(iyear)=q_rate_fcn_cHL(iyear-1);} 
      }   
      for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++)
      {   if (iyear>styr_cpue_rHB & iyear <=2003) 
          {
             q_rate_fcn_rHB(iyear)=(1.0+(iyear-styr_cpue_rHB)*q_rate)*q_rate_fcn_rHB(styr_cpue_rHB);  
          }
          if (iyear>2003) {q_rate_fcn_rHB(iyear)=q_rate_fcn_rHB(iyear-1);} 
      }   
      
      
          
             
          
          
      
      
  } 

 
  if (q_DD_beta>0.0) 
  {
    B_q_DD+=dzero;
    for (iyear=styr;iyear<=endyr;iyear++)
        {q_DD_fcn(iyear)=pow(B0_q_DD,q_DD_beta)*pow(B_q_DD(iyear),-q_DD_beta);}
          
  }  

     
FUNCTION get_indices


 
  q_cHL(styr_cpue_cHL)=mfexp(log_q_cpue_cHL); 
  for (iyear=styr_cpue_cHL; iyear<=endyr_cpue_cHL; iyear++)
  {
      N_cHL(iyear)=elem_prod(elem_prod(N_mdyr(iyear),sel_cHL(iyear)),wholewgt_cHL_klb(iyear));   
      pred_cHL_cpue(iyear)=q_cHL(iyear)*q_rate_fcn_cHL(iyear)*q_DD_fcn(iyear)*sum(N_cHL(iyear));
      if (iyear<endyr_cpue_cHL){q_cHL(iyear+1)=q_cHL(iyear)*mfexp(q_RW_log_dev_cHL(iyear));}
  }

 
  q_rHB(styr_cpue_rHB)=mfexp(log_q_cpue_rHB); 
  for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++)
  {   
      N_rHB(iyear)=elem_prod(N_mdyr(iyear),sel_rHB(iyear)); 
      pred_rHB_cpue(iyear)=q_rHB(iyear)*q_rate_fcn_rHB(iyear)*q_DD_fcn(iyear)*sum(N_rHB(iyear));
      if (iyear<endyr_cpue_rHB){q_rHB(iyear+1)=q_rHB(iyear)*mfexp(q_RW_log_dev_rHB(iyear));}
  }

 
  
  
  
  
      
      
      
      
      
  
 
 
  q_sCT(styr_cpue_sCT)=mfexp(log_q_cpue_sCT); 
  for (iyear=styr_cpue_sCT; iyear<=endyr_cpue_sCT; iyear++)
  {   
    N_sCT(iyear)=elem_prod(N_mdyr(iyear),sel_sCT(iyear)); 
	pred_sCT_cpue(iyear)=q_sCT(iyear)*q_DD_fcn(iyear)*sum(N_sCT(iyear));
    if (iyear<endyr_cpue_sCT){q_sCT(iyear+1)=q_sCT(iyear)*mfexp(q_RW_log_dev_sCT(iyear));}
  } 
  
FUNCTION get_length_comps

  
  for (iyear=1;iyear<=nyr_lenc_cHL;iyear++)
  {pred_cHL_lenc(iyear)=(L_cHL_num(yrs_lenc_cHL(iyear))*lenprob_cHL)/sum(L_cHL_num(yrs_lenc_cHL(iyear)));}
 
  
  for (iyear=1;iyear<=nyr_lenc_cOT;iyear++)
  {pred_cOT_lenc(iyear)=(L_cOT_num(yrs_lenc_cOT(iyear))*lenprob_cOT)/sum(L_cOT_num(yrs_lenc_cOT(iyear)));}

  
  for (iyear=1;iyear<=nyr_lenc_rHB;iyear++) 
  {pred_rHB_lenc(iyear)=(L_rHB_num(yrs_lenc_rHB(iyear))*lenprob_rHB)/sum(L_rHB_num(yrs_lenc_rHB(iyear)));}
  
  
  for (iyear=1;iyear<=nyr_lenc_rHB_D;iyear++) 
  {pred_rHB_D_lenc(iyear)=(D_rHB_num(yrs_lenc_rHB_D(iyear))*lenprob_rHB_D)/sum(D_rHB_num(yrs_lenc_rHB_D(iyear)));}

 
  for (iyear=1;iyear<=nyr_lenc_rGN;iyear++) 
  {pred_rGN_lenc(iyear)=(L_rGN_num(yrs_lenc_rGN(iyear))*lenprob_rGN)/sum(L_rGN_num(yrs_lenc_rGN(iyear)));}
 
  
  for (iyear=1;iyear<=nyr_lenc_sCT;iyear++) 
  {pred_sCT_lenc(iyear)=(N_sCT(yrs_lenc_sCT(iyear))*lenprob_sCT)/sum(N_sCT(yrs_lenc_sCT(iyear)));}

  
FUNCTION get_age_comps
   
  
  for (iyear=1;iyear<=nyr_agec_cHL;iyear++) 
  {
    ErrorFree_cHL_agec(iyear)=L_cHL_num(yrs_agec_cHL(iyear))/sum(L_cHL_num(yrs_agec_cHL(iyear)));  
    pred_cHL_agec_allages(iyear)=age_error*(ErrorFree_cHL_agec(iyear)/sum(ErrorFree_cHL_agec(iyear)));   
    for (iage=1; iage<=nages_agec; iage++) {pred_cHL_agec(iyear,iage)=pred_cHL_agec_allages(iyear,iage);} 
    
  }
 
  
 for (iyear=1;iyear<=nyr_agec_rHB;iyear++)
  {
    ErrorFree_rHB_agec(iyear)=L_rHB_num(yrs_agec_rHB(iyear))/sum(L_rHB_num(yrs_agec_rHB(iyear)));
    pred_rHB_agec_allages(iyear)=age_error*ErrorFree_rHB_agec(iyear); 
    for (iage=1; iage<=nages_agec; iage++) {pred_rHB_agec(iyear,iage)=pred_rHB_agec_allages(iyear,iage);} 
    
  }
  
  
   
 for (iyear=1;iyear<=nyr_agec_rGN;iyear++)
  {
    ErrorFree_rGN_agec(iyear)=L_rGN_num(yrs_agec_rGN(iyear))/sum(L_rGN_num(yrs_agec_rGN(iyear)));
    pred_rGN_agec_allages(iyear)=age_error*ErrorFree_rGN_agec(iyear); 
    for (iage=1; iage<=nages_agec; iage++) {pred_rGN_agec(iyear,iage)=pred_rGN_agec_allages(iyear,iage);} 
    
  }
 
   
 for (iyear=1;iyear<=nyr_agec_sCT;iyear++)
  {
    ErrorFree_sCT_agec(iyear)=N_sCT(yrs_agec_sCT(iyear))/sum(N_sCT(yrs_agec_sCT(iyear)));
    pred_sCT_agec_allages(iyear)=age_error*ErrorFree_sCT_agec(iyear); 
    for (iage=1; iage<=nages_agec; iage++) {pred_sCT_agec(iyear,iage)=pred_sCT_agec_allages(iyear,iage);} 
    
  }
  

FUNCTION get_weighted_current 
  F_temp_sum=0.0;
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cHL+
        sum(log_dev_F_L_cHL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cOT+
        sum(log_dev_F_L_cOT((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rHB+
        sum(log_dev_F_L_rHB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rGN+
        sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_D_cHL+
        sum(log_dev_F_D_cHL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_D_rHB+
        sum(log_dev_F_D_rHB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_D_rGN+
        sum(log_dev_F_D_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);

  F_cHL_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cHL+
        sum(log_dev_F_L_cHL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_cOT_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cOT+
        sum(log_dev_F_L_cOT((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_rHB_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rHB+
        sum(log_dev_F_L_rHB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_rGN_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rGN+
        sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_cHL_D_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_D_cHL+
        sum(log_dev_F_D_cHL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_rHB_D_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_D_rHB+
        sum(log_dev_F_D_rHB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_rGN_D_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_D_rGN+
        sum(log_dev_F_D_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  
  log_F_dev_end_cHL=sum(log_dev_F_L_cHL((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_cOT=sum(log_dev_F_L_cOT((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_rHB=sum(log_dev_F_L_rHB((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  
  log_F_dev_end_rGN=sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  

  log_F_dev_end_cHL_D=sum(log_dev_F_D_cHL((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_rHB_D=sum(log_dev_F_D_rHB((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_rGN_D=sum(log_dev_F_D_rGN((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;

  F_end_L=sel_cHL(endyr)*mfexp(log_avg_F_L_cHL+log_F_dev_end_cHL)+
          sel_cOT(endyr)*mfexp(log_avg_F_L_cOT+log_F_dev_end_cOT)+
          sel_rHB(endyr)*mfexp(log_avg_F_L_rHB+log_F_dev_end_rHB)+
          sel_rGN(endyr)*mfexp(log_avg_F_L_rGN+log_F_dev_end_rGN);   

  F_end_D=sel_D(endyr)*mfexp(log_avg_F_D_cHL+log_F_dev_end_cHL_D)+
          sel_D(endyr)*mfexp(log_avg_F_D_rHB+log_F_dev_end_rHB_D)+
          sel_D(endyr)*mfexp(log_avg_F_D_rGN+log_F_dev_end_rGN_D);   
    
  F_end=F_end_L+F_end_D;
  F_end_apex=max(F_end);
  
  sel_wgted_tot=F_end/F_end_apex;
  sel_wgted_L=elem_prod(sel_wgted_tot, elem_div(F_end_L,F_end));
  sel_wgted_D=elem_prod(sel_wgted_tot, elem_div(F_end_D,F_end));
  
  wgt_wgted_L_denom=F_cHL_prop+F_cOT_prop+F_rHB_prop+F_rGN_prop;  
  wgt_wgted_L_klb=F_cHL_prop/wgt_wgted_L_denom*wholewgt_cHL_klb(endyr)+ 
				  F_cOT_prop/wgt_wgted_L_denom*wholewgt_cOT_klb(endyr)+ 
                  F_rHB_prop/wgt_wgted_L_denom*wholewgt_rHB_klb(endyr)+
                  F_rGN_prop/wgt_wgted_L_denom*wholewgt_rGN_klb(endyr);                          

  wgt_wgted_D_denom=F_cHL_D_prop+F_rHB_D_prop+F_rGN_D_prop;  
  wgt_wgted_D_klb=F_cHL_D_prop/wgt_wgted_D_denom*wholewgt_cHL_D_klb(endyr)+ 
                  F_rHB_D_prop/wgt_wgted_D_denom*wholewgt_rHB_D_klb(endyr)+
                  F_rGN_D_prop/wgt_wgted_D_denom*wholewgt_rGN_D_klb(endyr);                
  
FUNCTION get_msy
  
  
  for(ff=1; ff<=n_iter_msy; ff++)
  {
    
    Z_age_msy=0.0;
    F_L_age_msy=0.0;
    F_D_age_msy=0.0;
          
    F_L_age_msy=F_msy(ff)*sel_wgted_L;
    F_D_age_msy=F_msy(ff)*sel_wgted_D;    
    Z_age_msy=M+F_L_age_msy+F_D_age_msy;         
    
    N_age_msy(1)=1.0;
    for (iage=2; iage<=nages; iage++)
      {N_age_msy(iage)=N_age_msy(iage-1)*mfexp(-1.*Z_age_msy(iage-1));}
    N_age_msy(nages)=N_age_msy(nages)/(1.0-mfexp(-1.*Z_age_msy(nages)));
    N_age_msy_spawn(1,(nages-1))=elem_prod(N_age_msy(1,(nages-1)),
                                   mfexp((-1.*Z_age_msy(1,(nages-1)))*spawn_time_frac));                 
    N_age_msy_spawn(nages)=(N_age_msy_spawn(nages-1)*(mfexp(-1.*(Z_age_msy(nages-1)*(1.0-spawn_time_frac) + 
                            Z_age_msy(nages)*spawn_time_frac) )))/(1.0-mfexp(-1.*Z_age_msy(nages)));
                     
    spr_msy(ff)=sum(elem_prod(N_age_msy_spawn,reprod));
	        
    R_eq(ff)=SR_eq_func(R0, steep, spr_msy(1), spr_msy(ff), BiasCor, SR_switch);
    
    if (R_eq(ff)<dzero) {R_eq(ff)=dzero;}    
    N_age_msy*=R_eq(ff);
    N_age_msy_spawn*=R_eq(ff);
    
    for (iage=1; iage<=nages; iage++)
    {
      L_age_msy(iage)=N_age_msy(iage)*(F_L_age_msy(iage)/Z_age_msy(iage))*
                      (1.-mfexp(-1.*Z_age_msy(iage)));
      D_age_msy(iage)=N_age_msy(iage)*(F_D_age_msy(iage)/Z_age_msy(iage))*
                      (1.-mfexp(-1.0*Z_age_msy(iage)));                      
    }
    
    SSB_eq(ff)=sum(elem_prod(N_age_msy_spawn,reprod));
	B_eq(ff)=sum(elem_prod(N_age_msy,wgt_mt));
    L_eq_klb(ff)=sum(elem_prod(L_age_msy,wgt_wgted_L_klb)); 
    L_eq_knum(ff)=sum(L_age_msy)/1000.0;  
    D_eq_klb(ff)=sum(elem_prod(D_age_msy,wgt_wgted_D_klb)); 
    D_eq_knum(ff)=sum(D_age_msy)/1000.0;    
  }  
  
  msy_klb_out=max(L_eq_klb); 
  
  for(ff=1; ff<=n_iter_msy; ff++)
  {
   if(L_eq_klb(ff) == msy_klb_out) 
      {    
        SSB_msy_out=SSB_eq(ff);
        B_msy_out=B_eq(ff);
        R_msy_out=R_eq(ff);
        msy_knum_out=L_eq_knum(ff);
        D_msy_knum_out=D_eq_knum(ff);
        D_msy_klb_out=D_eq_klb(ff);        
        F_msy_out=F_msy(ff);  
        spr_msy_out=spr_msy(ff);      
      }
  }


FUNCTION get_per_recruit_stuff

  
 
  for(iyear=styr; iyear<=endyr; iyear++)
  {
    N_age_spr(1)=1.0;
    for(iage=2; iage<=nages; iage++)
     {N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z(iyear,iage-1));}
    N_age_spr(nages)=N_age_spr(nages)/(1.0-mfexp(-1.*Z(iyear,nages)));    
    N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                mfexp(-1.*Z(iyear)(1,(nages-1))*spawn_time_frac));
    N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z(iyear)(nages-1)*(1.0-spawn_time_frac) + Z(iyear)(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z(iyear)(nages)));           
    spr_static(iyear)=sum(elem_prod(N_age_spr_spawn,reprod))/spr_F0;
  }
 
  
  for(ff=1; ff<=n_iter_spr; ff++)
  {
    
    Z_age_spr=0.0;
    F_L_age_spr=0.0;

    F_L_age_spr=F_spr(ff)*sel_wgted_L;
    F_D_age_spr=F_spr(ff)*sel_wgted_D;
    Z_age_spr=M+F_L_age_spr+F_D_age_spr;

    N_age_spr(1)=1.0;
    for (iage=2; iage<=nages; iage++)
     {N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z_age_spr(iage-1));}
    N_age_spr(nages)=N_age_spr(nages)/(1-mfexp(-1.*Z_age_spr(nages)));
    N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                   mfexp((-1.*Z_age_spr(1,(nages-1)))*spawn_time_frac));                 
    N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_spr(nages)));
    spr_spr(ff)=sum(elem_prod(N_age_spr_spawn,reprod));
    L_spr(ff)=0.0;
    for (iage=1; iage<=nages; iage++)
    {
      L_age_spr(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));
      L_spr(ff)+=L_age_spr(iage)*wgt_wgted_L_klb(iage)*1000.0; 
    }   
  }
  spr_ratio=spr_spr/spr_F0;
  F20_dum=min(fabs(spr_ratio-0.2));
  F30_dum=min(fabs(spr_ratio-0.3));
  F40_dum=min(fabs(spr_ratio-0.4));
  for(ff=1; ff<=n_iter_spr; ff++)
  {   
      if (fabs(spr_ratio(ff)-0.2)==F20_dum) {F20_out=F_spr(ff);}	  
	  if (fabs(spr_ratio(ff)-0.3)==F30_dum) {F30_out=F_spr(ff);}	  
	  if (fabs(spr_ratio(ff)-0.4)==F40_dum) {F40_out=F_spr(ff);}
  }
  rec=column(N,1);
  rec_mean=sum(rec(styr_rec_spr, endyr_rec_spr))/nyrs_rec_spr;
  R_F30_out=rec_mean;
  F_L_age_spr=F30_out*sel_wgted_L;
  F_D_age_spr=F30_out*sel_wgted_D;
  Z_age_spr=M+F_L_age_spr+F_D_age_spr;

  N_age_spr(1)=R_F30_out;
  for (iage=2; iage<=nages; iage++)
     {N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z_age_spr(iage-1));}
  N_age_spr(nages)=N_age_spr(nages)/(1-mfexp(-1.*Z_age_spr(nages)));
  N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                   mfexp((-1.*Z_age_spr(1,(nages-1)))*spawn_time_frac));                 
  N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_spr(nages)));

  for (iage=1; iage<=nages; iage++)
    {
      L_age_F30(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));
      D_age_F30(iage)=N_age_spr(iage)*(F_D_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.0*Z_age_spr(iage)));                      
    }
  SSB_F30_out=sum(elem_prod(N_age_spr_spawn,reprod));
  B_F30_out=sum(elem_prod(N_age_spr,wgt_mt));
  L_F30_klb_out=sum(elem_prod(L_age_F30,wgt_wgted_L_klb)); 
  L_F30_knum_out=sum(L_age_F30)/1000.0;  
  D_F30_klb_out=sum(elem_prod(D_age_F30,wgt_wgted_D_klb)); 
  D_F30_knum_out=sum(D_age_F30)/1000.0;    

	

FUNCTION get_miscellaneous_stuff


  if(var_rec_dev>0.0)
   {sigma_rec_dev=sqrt(var_rec_dev);} 
   else{sigma_rec_dev=0.0;}

  len_cv=elem_div(len_sd,meanlen_TL);
  
  
  L_total_num.initialize();
  L_total_klb.initialize();
  L_total_knum_yr.initialize();
  L_total_klb_yr.initialize();  
  D_total_num.initialize();
  D_total_klb.initialize();
  D_total_knum_yr.initialize();
  D_total_klb_yr.initialize();
  D_cHL_klb.initialize();
  D_rHB_klb.initialize();
  D_rGN_klb.initialize();
  
  for(iyear=styr; iyear<=endyr; iyear++)
  {
        L_total_klb_yr(iyear)=pred_cHL_L_klb(iyear)+pred_cOT_L_klb(iyear)+pred_rHB_L_klb(iyear)+pred_rGN_L_klb(iyear);
        L_total_knum_yr(iyear)=pred_cHL_L_knum(iyear)+pred_cOT_L_knum(iyear)+pred_rHB_L_knum(iyear)+pred_rGN_L_knum(iyear);
                
        B(iyear)=elem_prod(N(iyear),wgt_mt);
        totN(iyear)=sum(N(iyear));
        totB(iyear)=sum(B(iyear));   
        
        
	    if (iyear>endyr_selex_phase1 && iyear<=endyr)			
        {
         D_total_knum_yr(iyear)+=pred_cHL_D_knum(iyear);
         D_total_klb_yr(iyear)+=pred_cHL_D_klb(iyear);
         D_cHL_klb(iyear)=elem_prod(D_cHL_num(iyear),wholewgt_cHL_D_klb(iyear));     
        }
                
        
		if (iyear>endyr_selex_phase1 && iyear<=endyr)
        {
         D_total_knum_yr(iyear)+=pred_rHB_D_knum(iyear);
         D_total_klb_yr(iyear)+=pred_rHB_D_klb(iyear);
         D_rHB_klb(iyear)=elem_prod(D_rHB_num(iyear),wholewgt_rHB_D_klb(iyear));     
        }    
    	
        
		if (iyear>=styr_D_rGN && iyear<=endyr)	
        {
         D_total_knum_yr(iyear)+=pred_rGN_D_knum(iyear);
         D_total_klb_yr(iyear)+=pred_rGN_D_klb(iyear);
         D_rGN_klb(iyear)=elem_prod(D_rGN_num(iyear),wholewgt_rGN_D_klb(iyear));     
        }                          
  }
  
  L_total_num=L_cHL_num+L_cOT_num+L_rHB_num+L_rGN_num;   
  L_total_klb=L_cHL_klb+L_cOT_klb+L_rHB_klb+L_rGN_klb;   
 
  D_total_num=(D_cHL_num+D_rHB_num+D_rGN_num);          
  D_total_klb=D_cHL_klb+D_rHB_klb+D_rGN_klb;            
  
  
  B(endyr+1)=elem_prod(N(endyr+1),wgt_mt);
  totN(endyr+1)=sum(N(endyr+1));
  totB(endyr+1)=sum(B(endyr+1));  
  SdS0=SSB/S0;
  
  Fend_mean_temp=1.0;
  for (iyear=1; iyear<=selpar_n_yrs_wgted; iyear++) {Fend_mean_temp*=Fapex(endyr-iyear+1);}
  Fend_mean=pow(Fend_mean_temp,(1.0/selpar_n_yrs_wgted));	  
  if(F_msy_out>0)
    {
      FdF_msy=Fapex/F_msy_out;
      FdF_msy_end=FdF_msy(endyr);
      FdF_msy_end_mean=Fend_mean/F_msy_out;
    }
  if(SSB_msy_out>0)
    {
      SdSSB_msy=SSB/SSB_msy_out;
      SdSSB_msy_end=SdSSB_msy(endyr);
    }  

	if(F30_out>0)
    {
	  FdF30=Fapex/F30_out;
	  FdF30_end_mean=Fend_mean/F30_out;
	}
  if(SSB_F30_out>0)
    {
      SdSSB_F30=SSB/SSB_F30_out;
	  Sdmsst_F30=SSB/(smsy2msst75*SSB_F30_out);
      SdSSB_F30_end=SdSSB_F30(endyr);
	  Sdmsst_F30_end=Sdmsst_F30(endyr);
    }  	
   
   for(iyear=styr_rec_dev; iyear<=endyr_rec_dev; iyear++)
     {log_rec_dev_output(iyear)=log_dev_rec(iyear);}
   
   for(iage=2; iage<=nages; iage++)
     {log_Nage_dev_output(iage)=log_dev_Nage(iage);}


FUNCTION get_projection
  
    switch(Fproj_switch){
       case 1: 
          F_reg_proj=Fend_mean;
          break;
       case 2: 
          F_reg_proj=F_msy_out;
          break;
       case 3: 
          F_reg_proj=F30_out;
          break;     
       case 4: 
          F_reg_proj=F40_out;
          break;          		  
       default: 
          cout << "Error in input: Projection switch Fproj_switch must be set to 1, 2, 3, or 4." << endl;
          cout << "Presently it is set to " << Fproj_switch <<"."<< endl;
          exit(0);          
   }

  N_proj(styr_proj)=N(endyr+1); 
 
  for (iyear=styr_proj; iyear<=endyr_proj; iyear++) 
  {     
        if (iyear<styr_regs) {F_proj(iyear)=Fend_mean;}
		else {F_proj(iyear)=Fproj_mult*F_reg_proj;}
		
		FL_age_proj=sel_wgted_L*F_proj(iyear);
		FD_age_proj=sel_wgted_D*F_proj(iyear);
		
        Z_proj(iyear)=M+FL_age_proj+FD_age_proj;
        N_spawn_proj(iyear)(1,nages)=elem_prod(N_proj(iyear)(1,nages),(mfexp(-1.*(Z_proj(iyear)(1,nages))*spawn_time_frac))); 
		SSB_proj(iyear)= sum(elem_prod(N_spawn_proj(iyear),reprod));
        B_proj(iyear)=sum(elem_prod(N_proj(iyear),wgt_mt)); 
	     
		for (iage=1; iage<=nages; iage++)
			{L_age_proj(iyear,iage)=N_proj(iyear,iage)*FL_age_proj(iage)*(1.-mfexp(-1.*Z_proj(iyear,iage)))/Z_proj(iyear,iage);
		     D_age_proj(iyear,iage)=N_proj(iyear,iage)*FD_age_proj(iage)*(1.-mfexp(-1.*Z_proj(iyear,iage)))/Z_proj(iyear,iage);
			}          
        L_knum_proj(iyear)=sum(L_age_proj(iyear))/1000.0;
	    D_knum_proj(iyear)=sum(D_age_proj(iyear))/1000.0;
	    L_klb_proj(iyear)=sum(elem_prod(L_age_proj(iyear),wgt_wgted_L_klb));     
        D_klb_proj(iyear)=sum(elem_prod(D_age_proj(iyear),wgt_wgted_D_klb));     
		
		if (iyear<endyr_proj) {
			N_proj(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB_proj(iyear),SR_switch);
			N_proj(iyear+1)(2,nages)=++elem_prod(N_proj(iyear)(1,nages-1),(mfexp(-1.*Z_proj(iyear)(1,nages-1))));
			N_proj(iyear+1,nages)+=N_proj(iyear,nages)*mfexp(-1.*Z_proj(iyear,nages)); 
		}
  }
   R_proj=column(N_proj,1);                          



FUNCTION evaluate_objective_function
  
  
  fval=0.0;
  fval_data=0.0;  




  f_cHL_cpue=0.0;
  f_cHL_cpue=lk_lognormal(pred_cHL_cpue, obs_cpue_cHL, obs_cv_cpue_cHL, w_cpue_cHL);
  fval+=f_cHL_cpue;
  fval_data+=f_cHL_cpue;  

  f_rHB_cpue=0.0;
  f_rHB_cpue=lk_lognormal(pred_rHB_cpue, obs_cpue_rHB, obs_cv_cpue_rHB, w_cpue_rHB);
  fval+=f_rHB_cpue;
  fval_data+=f_rHB_cpue;  

  
  
  
  
  
  f_sCT_cpue=0.0;
  f_sCT_cpue=lk_lognormal(pred_sCT_cpue, obs_cpue_sCT, obs_cv_cpue_sCT, w_cpue_sCT);
  fval+=f_sCT_cpue;
  fval_data+=f_sCT_cpue;  
  

  
  
  f_cHL_L=lk_lognormal(pred_cHL_L_klb(styr_L_cHL,endyr_L_cHL), obs_L_cHL(styr_L_cHL,endyr_L_cHL),
                      obs_cv_L_cHL(styr_L_cHL,endyr_L_cHL), w_L);
  fval+=f_cHL_L;
  fval_data+=f_cHL_L;

  
  f_cOT_L=lk_lognormal(pred_cOT_L_klb(styr_L_cOT,endyr_L_cOT), obs_L_cOT(styr_L_cOT,endyr_L_cOT),
                      obs_cv_L_cOT(styr_L_cOT,endyr_L_cOT), w_L);
  fval+=f_cOT_L;
  fval_data+=f_cOT_L;
  
  
  f_rHB_L=lk_lognormal(pred_rHB_L_knum(styr_L_rHB,endyr_L_rHB), obs_L_rHB(styr_L_rHB,endyr_L_rHB), 
                      obs_cv_L_rHB(styr_L_rHB,endyr_L_rHB), w_L);
  fval+=f_rHB_L;
  fval_data+=f_rHB_L;  

  
  f_rGN_L=lk_lognormal(pred_rGN_L_knum(styr_L_rGN,endyr_L_rGN), obs_L_rGN(styr_L_rGN,endyr_L_rGN), 
                      obs_cv_L_rGN(styr_L_rGN,endyr_L_rGN), w_L);
  fval+=f_rGN_L;
  fval_data+=f_rGN_L;  


  
  
  f_cHL_D=lk_lognormal(pred_cHL_D_knum(styr_D_cHL,endyr_D_cHL), obs_cHL_D(styr_D_cHL,endyr_D_cHL), 
                      obs_cv_D_cHL(styr_D_cHL,endyr_D_cHL), w_D);
  fval+=f_cHL_D;
  fval_data+=f_cHL_D;  

  
  f_rHB_D=lk_lognormal(pred_rHB_D_knum(styr_D_rHB,endyr_D_rHB), obs_rHB_D(styr_D_rHB,endyr_D_rHB), 
                      obs_cv_D_rHB(styr_D_rHB,endyr_D_rHB), w_D);
  fval+=f_rHB_D;
  fval_data+=f_rHB_D;  

  
  f_rGN_D=lk_lognormal(pred_rGN_D_knum(styr_D_rGN,endyr_D_rGN), obs_rGN_D(styr_D_rGN,endyr_D_rGN), 
                      obs_cv_D_rGN(styr_D_rGN,endyr_D_rGN), w_D);
  fval+=f_rGN_D;
  fval_data+=f_rGN_D;  
  


  
  
  
  f_cHL_lenc=lk_dirichlet_multinomial(nsamp_lenc_cHL, pred_cHL_lenc, obs_lenc_cHL, nyr_lenc_cHL, double(nlenbins), minSS_lenc_cHL, log_dm_lenc_cHL);
  fval+=f_cHL_lenc;
  fval_data+=f_cHL_lenc;

  
  
  
  f_cOT_lenc=lk_dirichlet_multinomial(nsamp_lenc_cOT, pred_cOT_lenc, obs_lenc_cOT, nyr_lenc_cOT, double(nlenbins), minSS_lenc_cOT, log_dm_lenc_cOT);
  fval+=f_cOT_lenc;
  fval_data+=f_cOT_lenc;

  
  
  
  f_rHB_lenc=lk_dirichlet_multinomial(nsamp_lenc_rHB, pred_rHB_lenc, obs_lenc_rHB, nyr_lenc_rHB, double(nlenbins), minSS_lenc_rHB, log_dm_lenc_rHB);
  fval+=f_rHB_lenc;
  fval_data+=f_rHB_lenc;
  
  
  
  
  f_rGN_lenc=lk_dirichlet_multinomial(nsamp_lenc_rGN, pred_rGN_lenc, obs_lenc_rGN, nyr_lenc_rGN, double(nlenbins), minSS_lenc_rGN, log_dm_lenc_rGN);
  fval+=f_rGN_lenc;
  fval_data+=f_rGN_lenc;
  
  
  
  
  f_rHB_D_lenc=lk_dirichlet_multinomial(nsamp_lenc_rHB_D, pred_rHB_D_lenc, obs_lenc_rHB_D, nyr_lenc_rHB_D, double(nlenbins), minSS_lenc_rHB_D, log_dm_lenc_rHB_D);
  fval+=f_rHB_D_lenc;
  fval_data+=f_rHB_D_lenc;
  
  
  
  
  f_sCT_lenc=lk_dirichlet_multinomial(nsamp_lenc_sCT, pred_sCT_lenc, obs_lenc_sCT, nyr_lenc_sCT, double(nlenbins), minSS_lenc_sCT, log_dm_lenc_sCT);
  fval+=f_sCT_lenc;
  fval_data+=f_sCT_lenc;
   


  
  
  
  f_cHL_agec=lk_dirichlet_multinomial(nsamp_agec_cHL, pred_cHL_agec, obs_agec_cHL, nyr_agec_cHL, double(nages_agec), minSS_agec_cHL, log_dm_agec_cHL);
  fval+=f_cHL_agec;
  fval_data+=f_cHL_agec;

  
  
  
  f_rHB_agec=lk_dirichlet_multinomial(nsamp_agec_rHB, pred_rHB_agec, obs_agec_rHB, nyr_agec_rHB, double(nages_agec), minSS_agec_rHB, log_dm_agec_rHB);
  fval+=f_rHB_agec;
  fval_data+=f_rHB_agec;

  
  
  
  f_rGN_agec=lk_dirichlet_multinomial(nsamp_agec_rGN, pred_rGN_agec, obs_agec_rGN, nyr_agec_rGN, double(nages_agec), minSS_agec_rGN, log_dm_agec_rGN);
  fval+=f_rGN_agec;
  fval_data+=f_rGN_agec;
  
  
  
  
  f_sCT_agec=lk_dirichlet_multinomial(nsamp_agec_sCT, pred_sCT_agec, obs_agec_sCT, nyr_agec_sCT, double(nages_agec), minSS_agec_sCT, log_dm_agec_sCT);
  fval+=f_sCT_agec;
  fval_data+=f_sCT_agec;
  



  
  
  f_Nage_init=norm2(log_dev_Nage);        
  fval+=w_Nage_init*f_Nage_init;
  
  f_rec_dev=0.0;
  
  rec_logL_add=nyrs_rec*log(rec_sigma);
  f_rec_dev=(square(log_dev_rec(styr_rec_dev) + rec_sigma_sq/2.0)/(2.0*rec_sigma_sq));
  for(iyear=(styr_rec_dev+1); iyear<=endyr_rec_dev; iyear++)
  {f_rec_dev+=(square(log_dev_rec(iyear)-R_autocorr*log_dev_rec(iyear-1) + rec_sigma_sq/2.0)/
               (2.0*rec_sigma_sq));}
  f_rec_dev+=rec_logL_add;            
  fval+=w_rec*f_rec_dev;
   
  f_rec_dev_early=0.0; 
  if (w_rec_early>0.0)
    { if (styr_rec_dev<endyr_rec_phase1)
        {  
          for(iyear=styr_rec_dev; iyear<=endyr_rec_phase1; iyear++)
          
          
          {f_rec_dev_early+=square(log_dev_rec(iyear));}
        }
  fval+=w_rec_early*f_rec_dev_early;
  }
  
  f_rec_dev_end=0.0; 
  if (w_rec_end>0.0)
  { if (endyr_rec_phase2<endyr_rec_dev)
        {  
          for(iyear=(endyr_rec_phase2+1); iyear<=endyr_rec_dev; iyear++)
          
          
          {f_rec_dev_end+=square(log_dev_rec(iyear));}
        }
      fval+=w_rec_end*f_rec_dev_end;
   }  

  
  f_Ftune=0.0; 
  if (w_Ftune>0.0)
  {if (set_Ftune>0.0 && !last_phase()) {f_Ftune=square(Fapex(set_Ftune_yr)-set_Ftune);}
   fval+=w_Ftune*f_Ftune;
  }

  
  f_fullF_constraint=0.0;
  if (w_fullF>0.0)
  {for (iyear=styr; iyear<=endyr; iyear++)
        {if(Fapex(iyear)>3.0) {f_fullF_constraint+=(mfexp(Fapex(iyear)-3.0)-1.0);}}
   fval+=w_fullF*f_fullF_constraint;
  }
  
 
 f_rHB_RWq_cpue=0.0;
 for (iyear=styr_cpue_rHB; iyear<endyr_cpue_rHB; iyear++)
     {f_rHB_RWq_cpue+=square(q_RW_log_dev_rHB(iyear))/(2.0*set_RWq_var);}
 fval+=f_rHB_RWq_cpue;   

 f_cHL_RWq_cpue=0.0;
 for (iyear=styr_cpue_cHL; iyear<endyr_cpue_cHL; iyear++)
     {f_cHL_RWq_cpue+=square(q_RW_log_dev_cHL(iyear))/(2.0*set_RWq_var);}
 fval+=f_cHL_RWq_cpue;   
  




  f_priors=0.0; 
  f_priors+=neg_log_prior(len_cv_val,set_len_cv(5),set_len_cv(6),set_len_cv(7));
   
  f_priors+=neg_log_prior(steep,set_steep(5),set_log_R0(6),set_log_R0(7)); 
  f_priors+=neg_log_prior(log_R0,set_log_R0(5),set_log_R0(6),set_log_R0(7)); 
  f_priors+=neg_log_prior(R_autocorr,set_R_autocorr(5),set_R_autocorr(6),set_R_autocorr(7));
  f_priors+=neg_log_prior(rec_sigma,set_rec_sigma(5),set_rec_sigma(6),set_rec_sigma(7));
 
  
  
  f_priors+=neg_log_prior(selpar_A50_cHL2,set_selpar_A50_cHL2(5), set_selpar_A50_cHL2(6), set_selpar_A50_cHL2(7));
  f_priors+=neg_log_prior(selpar_slope_cHL2,set_selpar_slope_cHL2(5), set_selpar_slope_cHL2(6), set_selpar_slope_cHL2(7));
  f_priors+=neg_log_prior(selpar_A50_cHL3,set_selpar_A50_cHL3(5), set_selpar_A50_cHL3(6), set_selpar_A50_cHL3(7));
  f_priors+=neg_log_prior(selpar_slope_cHL3,set_selpar_slope_cHL3(5), set_selpar_slope_cHL3(6), set_selpar_slope_cHL3(7));

  f_priors+=neg_log_prior(selpar_A50_cOT2,set_selpar_A50_cOT2(5), set_selpar_A50_cOT2(6), set_selpar_A50_cOT2(7));
  f_priors+=neg_log_prior(selpar_A50_cOT3,set_selpar_A50_cOT3(5), set_selpar_A50_cOT3(6), set_selpar_A50_cOT3(7));
  f_priors+=neg_log_prior(selpar_slope_cOT2,set_selpar_slope_cOT2(5), set_selpar_slope_cOT2(6), set_selpar_slope_cOT2(7));
  f_priors+=neg_log_prior(selpar_A502_cOT2,set_selpar_A502_cOT2(5), set_selpar_A502_cOT2(6), set_selpar_A502_cOT2(7));
  f_priors+=neg_log_prior(selpar_slope2_cOT2,set_selpar_slope2_cOT2(5), set_selpar_slope2_cOT2(6), set_selpar_slope2_cOT2(7));
  
  f_priors+=neg_log_prior(selpar_A50_rHB1,set_selpar_A50_rHB1(5), set_selpar_A50_rHB1(6), set_selpar_A50_rHB1(7));
  f_priors+=neg_log_prior(selpar_slope_rHB1,set_selpar_slope_rHB1(5), set_selpar_slope_rHB1(6), set_selpar_slope_rHB1(7));
  f_priors+=neg_log_prior(selpar_A50_rHB2,set_selpar_A50_rHB2(5), set_selpar_A50_rHB2(6), set_selpar_A50_rHB2(7));
  f_priors+=neg_log_prior(selpar_slope_rHB2,set_selpar_slope_rHB2(5), set_selpar_slope_rHB2(6), set_selpar_slope_rHB2(7));
  f_priors+=neg_log_prior(selpar_A50_rHB3,set_selpar_A50_rHB3(5), set_selpar_A50_rHB3(6), set_selpar_A50_rHB3(7));
  f_priors+=neg_log_prior(selpar_slope_rHB3,set_selpar_slope_rHB3(5), set_selpar_slope_rHB3(6), set_selpar_slope_rHB3(7));

  f_priors+=neg_log_prior(selpar_A50_rGN3,set_selpar_A50_rGN3(5), set_selpar_A50_rGN3(6), set_selpar_A50_rGN3(7));
  f_priors+=neg_log_prior(selpar_slope_rGN3,set_selpar_slope_rGN3(5), set_selpar_slope_rGN3(6), set_selpar_slope_rGN3(7));
  
  f_priors+=neg_log_prior(selpar_A50_sCT,set_selpar_A50_sCT(5), set_selpar_A50_sCT(6), set_selpar_A50_sCT(7));
  f_priors+=neg_log_prior(selpar_slope_sCT,set_selpar_slope_sCT(5), set_selpar_slope_sCT(6), set_selpar_slope_sCT(7));
  f_priors+=neg_log_prior(selpar_A502_sCT,set_selpar_A502_sCT(5), set_selpar_A502_sCT(6), set_selpar_A502_sCT(7));
  f_priors+=neg_log_prior(selpar_slope2_sCT,set_selpar_slope2_sCT(5), set_selpar_slope2_sCT(6), set_selpar_slope2_sCT(7));
  
  f_priors+=neg_log_prior(selpar_age1logit_D,set_selpar_age1logit_D(5), set_selpar_age1logit_D(6), set_selpar_age1logit_D(7));

  f_priors+=neg_log_prior(log_q_cpue_cHL,set_log_q_cpue_cHL(5),set_log_q_cpue_cHL(6),set_log_q_cpue_cHL(7));
  f_priors+=neg_log_prior(log_q_cpue_rHB,set_log_q_cpue_rHB(5),set_log_q_cpue_rHB(6),set_log_q_cpue_rHB(7));
  
  f_priors+=neg_log_prior(log_q_cpue_sCT,set_log_q_cpue_sCT(5),set_log_q_cpue_sCT(6),set_log_q_cpue_sCT(7));
      
  f_priors+=neg_log_prior(log_dm_lenc_cHL,set_log_dm_lenc_cHL(5),set_log_dm_lenc_cHL(6),set_log_dm_lenc_cHL(7));
  f_priors+=neg_log_prior(log_dm_lenc_cOT,set_log_dm_lenc_cOT(5),set_log_dm_lenc_cOT(6),set_log_dm_lenc_cOT(7));
  f_priors+=neg_log_prior(log_dm_lenc_rHB,set_log_dm_lenc_rHB(5),set_log_dm_lenc_rHB(6),set_log_dm_lenc_rHB(7));
  f_priors+=neg_log_prior(log_dm_lenc_rGN,set_log_dm_lenc_rGN(5),set_log_dm_lenc_rGN(6),set_log_dm_lenc_rGN(7));
  f_priors+=neg_log_prior(log_dm_lenc_rHB_D,set_log_dm_lenc_rHB_D(5),set_log_dm_lenc_rHB_D(6),set_log_dm_lenc_rHB_D(7));
  f_priors+=neg_log_prior(log_dm_lenc_sCT,set_log_dm_lenc_sCT(5),set_log_dm_lenc_sCT(6),set_log_dm_lenc_sCT(7));
  f_priors+=neg_log_prior(log_dm_agec_cHL,set_log_dm_agec_cHL(5),set_log_dm_agec_cHL(6),set_log_dm_agec_cHL(7));
  f_priors+=neg_log_prior(log_dm_agec_rHB,set_log_dm_agec_rHB(5),set_log_dm_agec_rHB(6),set_log_dm_agec_rHB(7));
  f_priors+=neg_log_prior(log_dm_agec_rGN,set_log_dm_agec_rGN(5),set_log_dm_agec_rGN(6),set_log_dm_agec_rGN(7));
  f_priors+=neg_log_prior(log_dm_agec_sCT,set_log_dm_agec_sCT(5),set_log_dm_agec_sCT(6),set_log_dm_agec_sCT(7));
  
  fval+=f_priors;



FUNCTION dvar_vector logistic(const dvar_vector& ages, const dvariable& A50, const dvariable& slope)
  
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1./(1.+mfexp(-1.*slope*(ages-A50))); 
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;



FUNCTION dvar_vector logistic_exponential(const dvar_vector& ages, const dvariable& A50, const dvariable& slope, const dvariable& sigma, const dvariable& joint)
  
  
  RETURN_ARRAYS_INCREMENT();                                                                                                                                    
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());                                                                                                         
  Sel_Tmp=1.0;                                                                                                                                                  
  for (iage=1; iage<=nages; iage++)                                                                                                                             
  {                                                                                                                                                             
   if (ages(iage)<joint) {Sel_Tmp(iage)=1./(1.+mfexp(-1.*slope*(ages(iage)-A50)));}                                                                             
   if (ages(iage)>joint){Sel_Tmp(iage)=mfexp(-1.*square((ages(iage)-joint)/sigma));}                                                                            
  }                                                                                                                                                             
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);                                                                                                                                 
  RETURN_ARRAYS_DECREMENT();                                                                                                                                    
  return Sel_Tmp;   



FUNCTION dvar_vector logistic_double(const dvar_vector& ages, const dvariable& A501, const dvariable& slope1, const dvariable& A502, const dvariable& slope2)
  
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=elem_prod( (1./(1.+mfexp(-1.*slope1*(ages-A501)))),(1.-(1./(1.+mfexp(-1.*slope2*(ages-(A501+A502)))))) );     
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;



FUNCTION dvar_vector logistic_joint(const dvar_vector& ages, const dvariable& A501, const dvariable& slope1, const dvariable& A502, const dvariable& slope2, const dvariable& satval, const dvariable& joint)
  
  
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1.0; 
  for (iage=1; iage<=nages; iage++)
  {
   if (double(iage)<joint) {Sel_Tmp(iage)=1./(1.+mfexp(-1.*slope1*(ages(iage)-A501)));}  
   if (double(iage)>joint){Sel_Tmp(iage)=1.0-(1.0-satval)/(1.+mfexp(-1.*slope2*(ages(iage)-A502)));}  
  }  
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;



FUNCTION dvar_vector gaussian_double(const dvar_vector& ages, const dvariable& peak, const dvariable& top, const dvariable& ascwid, const dvariable& deswid, const dvariable& init, const dvariable& final)
  
  
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step1(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step2(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step3(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step4(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step5(ages.indexmin(),ages.indexmax());
  dvar_vector sel_step6(ages.indexmin(),ages.indexmax());
  dvar_vector pars_tmp(1,6); dvar_vector sel_tmp_iq(1,2);
  
  pars_tmp(1)=peak;
  pars_tmp(2)=peak+1.0+(0.99*ages(nages)-peak-1.0)/(1.0+mfexp(-top));
  pars_tmp(3)=mfexp(ascwid);
  pars_tmp(4)=mfexp(deswid);
  pars_tmp(5)=1.0/(1.0+mfexp(-init));
  pars_tmp(6)=1.0/(1.0+mfexp(-final));
       
  sel_tmp_iq(1)=mfexp(-(square(ages(1)-pars_tmp(1))/pars_tmp(3)));
  sel_tmp_iq(2)=mfexp(-(square(ages(nages)-pars_tmp(2))/pars_tmp(4)));
  
  sel_step1=mfexp(-(square(ages-pars_tmp(1))/pars_tmp(3)));
  sel_step2=pars_tmp(5)+(1.0-pars_tmp(5))*(sel_step1-sel_tmp_iq(1))/(1.0-sel_tmp_iq(1));  
  sel_step3=mfexp(-(square(ages-pars_tmp(2))/pars_tmp(4)));
  sel_step4=1.0+(pars_tmp(6)-1.0)*(sel_step3-1.0)/(sel_tmp_iq(2)-1.0);
  sel_step5=1.0/ (1.0+mfexp(-(20.0* elem_div((ages-pars_tmp(1)), (1.0+sfabs(ages-pars_tmp(1)))) )));
  sel_step6=1.0/(1.0+mfexp(-(20.0*elem_div((ages-pars_tmp(2)),(1.0+sfabs(ages-pars_tmp(2)))) )));  

  Sel_Tmp=elem_prod(sel_step2,(1.0-sel_step5))+ 
          elem_prod(sel_step5,((1.0-sel_step6)+ elem_prod(sel_step4,sel_step6)) ); 
 
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
    


FUNCTION dvariable SR_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& SSB, int func)
  
  
  RETURN_ARRAYS_INCREMENT();
  dvariable Recruits_Tmp;
  switch(func) {
    case 1: 
      Recruits_Tmp=((0.8*R0*h*SSB)/(0.2*R0*spr_F0*(1.0-h)+(h-0.2)*SSB));       
    break;
    case 2: 
      Recruits_Tmp=((SSB/spr_F0)*mfexp(h*(1-SSB/(R0*spr_F0))));       
    break;
  }
  RETURN_ARRAYS_DECREMENT();
  return Recruits_Tmp;
  


FUNCTION dvariable SR_eq_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& spr_F, const dvariable& BC, int func)
  
  
  RETURN_ARRAYS_INCREMENT();
  dvariable Recruits_Tmp;
  switch(func) {
    case 1: 
      Recruits_Tmp=(R0/((5.0*h-1.0)*spr_F))*(BC*4.0*h*spr_F-spr_F0*(1.0-h));    
    break;
    case 2: 
      Recruits_Tmp=R0/(spr_F/spr_F0)*(1.0+log(BC*spr_F/spr_F0)/h);      
    break;
  }
  RETURN_ARRAYS_DECREMENT();
  return Recruits_Tmp;



FUNCTION dvariable multinom_eff_N(const dvar_vector& pred_comp, const dvar_vector& obs_comp)
  
  dvariable EffN_Tmp; dvariable numer; dvariable denom;
  RETURN_ARRAYS_INCREMENT();
  numer=sum( elem_prod(pred_comp,(1.0-pred_comp)) );
  denom=sum( square(obs_comp-pred_comp) );
  if (denom>0.0) {EffN_Tmp=numer/denom;}
  else {EffN_Tmp=-missing;}                            
  RETURN_ARRAYS_DECREMENT();
  return EffN_Tmp;



FUNCTION dvariable lk_lognormal(const dvar_vector& pred, const dvar_vector& obs, const dvar_vector& cv, const dvariable& wgt_dat)
  
  
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.0001;
  dvar_vector var(cv.indexmin(),cv.indexmax()); 
  var=log(1.0+square(cv/wgt_dat));   
  LkvalTmp=sum(0.5*elem_div(square(log(elem_div((pred+small_number),(obs+small_number)))),var) );
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;



FUNCTION dvariable lk_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const double& minSS, const dvariable& wgt_dat)
  
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.0001;
  LkvalTmp=0.0;
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {LkvalTmp-=wgt_dat*nsamp(ii)*sum(elem_prod((obs_comp(ii)+small_number),
               log(elem_div((pred_comp(ii)+small_number), (obs_comp(ii)+small_number)))));
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;



FUNCTION dvariable lk_robust_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS, const dvariable& wgt_dat)
  
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.0001;
  LkvalTmp=0.0;
  dvar_matrix Eprime=elem_prod((1.0-obs_comp), obs_comp)+0.1/mbin; 
  dvar_vector nsamp_wgt=nsamp*wgt_dat;
  
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {LkvalTmp+= sum(0.5*log(Eprime(ii))-log(small_number+mfexp(elem_div((-square(obs_comp(ii)-pred_comp(ii))) , (Eprime(ii)*2.0/nsamp_wgt(ii)) ))) );
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;
  


FUNCTION dvariable lk_dirichlet_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS, const dvariable& log_dir_par)
  
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  LkvalTmp=0.0; 
  dvar_vector nsamp_adjust=nsamp*mfexp(log_dir_par);
  
  for (int ii=1; ii<=ncomp; ii++)
  {
	if (nsamp(ii)>=minSS)
    {
		LkvalTmp-=gammln(nsamp_adjust(ii))-gammln(nsamp(ii)+nsamp_adjust(ii));
		LkvalTmp-=sum(gammln(nsamp(ii)*obs_comp(ii)+nsamp_adjust(ii)*pred_comp(ii)));
        LkvalTmp+=sum(gammln(nsamp_adjust(ii)*pred_comp(ii)));		
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;



FUNCTION dvariable lk_logistic_normal(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS)
  
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.0001;
  LkvalTmp=0.0;
  dvar_matrix nu=pred_comp+0.0;
  dvar_matrix pred_plus=pred_comp+small_number;
  dvar_matrix obs_plus=obs_comp+small_number;

  dvariable nu_mean;
  dvariable nu_sum_sq;
  dvariable tau_hat_sq;
  dvariable year_count; 
   
  LkvalTmp=0.0;
  nu_sum_sq=0.0;
  year_count=0.0;
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {
		year_count+=1.0;
		nu_mean=sum( log(obs_plus(ii))-log(pred_plus(ii))  )/mbin;	
		for (int jj=1; jj<=mbin;jj++)
		{
			nu(ii,jj) = log(obs_plus(ii,jj)) - log(pred_plus(ii,jj)) - nu_mean;
			nu_sum_sq += square(nu(ii,jj));
		}
    }
  }  
  if (year_count>0.0)
  {
	  tau_hat_sq = nu_sum_sq/((mbin-1.0)*year_count);
	  LkvalTmp = (mbin-1.0)*year_count*log(tau_hat_sq);
  }
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;
  



FUNCTION  dvariable neg_log_prior(dvariable pred, const double& prior, dvariable var, int pdf)
  
    dvariable LkvalTmp;
    dvariable alpha, beta, ab_iq;
    dvariable big_number=1e10;
    LkvalTmp=0.0;
    
    switch(pdf) {
        case 1: 
          LkvalTmp=0.0;
          break;
        case 2: 
          if(prior<=0.0) cout << "YIKES: Don't use a lognormal distn for a negative prior" << endl;
          else if(pred<=0) LkvalTmp=big_number=1e10;
          else {
            if(var<0.0) var=log(1.0+var*var) ;      
            LkvalTmp= 0.5*( square(log(pred/prior))/var + log(var) );
          }
        break;
        case 3: 
          if(var<0.0 && prior!=0.0) var=square(var*prior);       
          else if(var<0.0 && prior==0.0) var=-var;               
          LkvalTmp= 0.5*( square(pred-prior)/var + log(var) );
          break;
        case 4: 
          if(var<0.0) var=square(var*prior);          
          if(prior<=0.0 || prior>=1.0) cout << "YIKES: Don't use a beta distn for a prior outside (0,1)" << endl;
          ab_iq=prior*(1.0-prior)/var - 1.0; alpha=prior*ab_iq; beta=(1.0-prior)*ab_iq;
          if(pred>=0 && pred<=1) LkvalTmp= (1.0-alpha)*log(pred)+(1.0-beta)*log(1.0-pred)-gammln(alpha+beta)+gammln(alpha)+gammln(beta);
          else LkvalTmp=big_number;
          break;
        default: 
          cout << "The prior must be either 1(lognormal), 2(normal), or 3(beta)." << endl;
          cout << "Presently it is " << pdf << endl;
          exit(0);
    }
    return LkvalTmp;



FUNCTION dvariable sdnr_multinomial(const double& ncomp, const dvar_vector& ages, const dvar_vector& nsamp, 
                                    const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const dvariable& wgt_dat)
  
  
  RETURN_ARRAYS_INCREMENT();
  dvariable SdnrTmp;
  dvar_vector o(1,ncomp);  
  dvar_vector p(1,ncomp);  
  dvar_vector ose(1,ncomp);  
  dvar_vector res(1,ncomp);
  SdnrTmp=0.0;
  for (int ii=1; ii<=ncomp; ii++)
  {
    o(ii)=sum(elem_prod(ages,obs_comp(ii)));
    p(ii)=sum(elem_prod(ages,pred_comp(ii)));
    ose(ii)=sqrt((sum(elem_prod(square(ages),pred_comp(ii)))-square(p(ii)))/(nsamp(ii)*wgt_dat));
  }
  res=elem_div((o-p),ose); 
  SdnrTmp=sqrt(sum(square(res-(sum(res)/ncomp))/(ncomp-1.0))); 
  RETURN_ARRAYS_DECREMENT();
  return SdnrTmp;



FUNCTION dvariable sdnr_lognormal(const dvar_vector& pred, const dvar_vector& obs, const dvar_vector& cv, const dvariable& wgt_dat)
  
  RETURN_ARRAYS_INCREMENT();
  dvariable SdnrTmp;
  dvariable small_number=0.00001;
  dvariable n;
  dvar_vector res(cv.indexmin(),cv.indexmax());
  SdnrTmp=0.0;
  res=elem_div(log(elem_div(obs+small_number,pred+small_number)),sqrt(log(1+square(cv/wgt_dat))));
  n=cv.indexmax()-cv.indexmin()+1;
  SdnrTmp=sqrt(sum(square(res-(sum(res)/n))/(n-1.0))); 
  RETURN_ARRAYS_DECREMENT();
  return SdnrTmp; 


REPORT_SECTION

  if (last_phase())  
  {
      cout<<"start report"<<endl;
      
	  get_weighted_current();
      
      get_msy();
      
      get_per_recruit_stuff();
      
      get_miscellaneous_stuff();
      
	  get_projection();
	  
	  
      grad_max=objective_function_value::pobjfun->gmax;
      time(&finish);
	  elapsed_time=difftime(finish,start);
	  hour=long(elapsed_time)/3600;
	  minute=long(elapsed_time)%3600/60;
	  second=(long(elapsed_time)%3600)%60;
	  cout<<endl<<endl<<"*******************************************"<<endl;
	  cout<<"--Start time: "<<ctime(&start)<<endl;
	  cout<<"--Finish time: "<<ctime(&finish)<<endl;
	  cout<<"--Runtime: ";
	  cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds"<<endl;
	  cout << "--TotalLikelihood: " << fval << endl;
      cout<<"--Final gradient: "<<objective_function_value::pobjfun->gmax << endl;
	  cout<<"*******************************************"<<endl;
      
      cout <<endl;     
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;
      
      cout <<"F status="<<FdF_msy_end<<endl;
      cout <<"Pop status="<<SdSSB_msy_end<<endl;
      cout << "h="<<steep<<"   R0="<<R0<<endl;
      
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;  
      
      report << "TotalLikelihood " << fval << endl;
      report << "N" << endl;
      report << N<<endl;
      report << "F" << endl;
      report << F <<endl;
      report << "prob_belowsizelim_block3" << endl;
	  report<<prob_belowsizelim_block3<<endl;
	  
      sdnr_lc_cHL=sdnr_multinomial(nyr_lenc_cHL, lenbins, nsamp_lenc_cHL, pred_cHL_lenc, obs_lenc_cHL, w_lenc_cHL);  
      sdnr_lc_cOT=sdnr_multinomial(nyr_lenc_cOT, lenbins, nsamp_lenc_cOT, pred_cOT_lenc, obs_lenc_cOT, w_lenc_cOT);  	  
      sdnr_lc_rHB=sdnr_multinomial(nyr_lenc_rHB, lenbins, nsamp_lenc_rHB, pred_rHB_lenc, obs_lenc_rHB, w_lenc_rHB); 
	  sdnr_lc_rGN=sdnr_multinomial(nyr_lenc_rGN, lenbins, nsamp_lenc_rGN, pred_rGN_lenc, obs_lenc_rGN, w_lenc_rGN); 
      sdnr_lc_rHB_D=sdnr_multinomial(nyr_lenc_rHB_D, lenbins, nsamp_lenc_rHB_D, pred_rHB_D_lenc, obs_lenc_rHB_D, w_lenc_rHB_D); 
	  sdnr_lc_sCT=sdnr_multinomial(nyr_lenc_sCT, lenbins, nsamp_lenc_sCT, pred_sCT_lenc, obs_lenc_sCT, w_lenc_sCT); 
       
      sdnr_ac_cHL=sdnr_multinomial(nyr_agec_cHL, agebins_agec, nsamp_agec_cHL, pred_cHL_agec, obs_agec_cHL, w_agec_cHL);  
      sdnr_ac_rHB=sdnr_multinomial(nyr_agec_rHB, agebins_agec, nsamp_agec_rHB, pred_rHB_agec, obs_agec_rHB, w_agec_rHB);  
	  sdnr_ac_rGN=sdnr_multinomial(nyr_agec_rGN, agebins_agec, nsamp_agec_rGN, pred_rGN_agec, obs_agec_rGN, w_agec_rGN);  
      sdnr_ac_sCT=sdnr_multinomial(nyr_agec_sCT, agebins_agec, nsamp_agec_sCT, pred_sCT_agec, obs_agec_sCT, w_agec_sCT);  
     
      sdnr_I_cHL=sdnr_lognormal(pred_cHL_cpue, obs_cpue_cHL, obs_cv_cpue_cHL, w_cpue_cHL);
      sdnr_I_rHB=sdnr_lognormal(pred_rHB_cpue, obs_cpue_rHB, obs_cv_cpue_rHB, w_cpue_rHB);
      
      sdnr_I_sCT=sdnr_lognormal(pred_sCT_cpue, obs_cpue_sCT, obs_cv_cpue_sCT, w_cpue_sCT);  
       
      
      
      
      
       Linf_out(8)=Linf; Linf_out(1,7)=set_Linf; 
       K_out(8)=K; K_out(1,7)=set_K;
       t0_out(8)=t0; t0_out(1,7)=set_t0;
       len_cv_val_out(8)=len_cv_val; len_cv_val_out(1,7)=set_len_cv;
	   	   
       log_R0_out(8)=log_R0; log_R0_out(1,7)=set_log_R0;
       M_constant_out(8)=M_constant; M_constant_out(1,7)=set_M_constant;
       steep_out(8)=steep; steep_out(1,7)=set_steep;
       rec_sigma_out(8)=rec_sigma; rec_sigma_out(1,7)=set_rec_sigma;
       R_autocorr_out(8)=R_autocorr; R_autocorr_out(1,7)=set_R_autocorr;
	   
	   log_dm_cHL_lc_out(8)=log_dm_lenc_cHL; log_dm_cHL_lc_out(1,7)=set_log_dm_lenc_cHL;
	   log_dm_cOT_lc_out(8)=log_dm_lenc_cOT; log_dm_cOT_lc_out(1,7)=set_log_dm_lenc_cOT;
	   log_dm_rHB_lc_out(8)=log_dm_lenc_rHB; log_dm_rHB_lc_out(1,7)=set_log_dm_lenc_rHB;
	   log_dm_rGN_lc_out(8)=log_dm_lenc_rGN; log_dm_rGN_lc_out(1,7)=set_log_dm_lenc_rGN;
	   log_dm_rHB_D_lc_out(8)=log_dm_lenc_rHB_D; log_dm_rHB_D_lc_out(1,7)=set_log_dm_lenc_rHB_D;
	   log_dm_sCT_lc_out(8)=log_dm_lenc_sCT; log_dm_sCT_lc_out(1,7)=set_log_dm_lenc_sCT;
	   log_dm_cHL_ac_out(8)=log_dm_agec_cHL; log_dm_cHL_ac_out(1,7)=set_log_dm_agec_cHL;
	   log_dm_rHB_ac_out(8)=log_dm_agec_rHB; log_dm_rHB_ac_out(1,7)=set_log_dm_agec_rHB;
	   log_dm_rGN_ac_out(8)=log_dm_agec_rGN; log_dm_rGN_ac_out(1,7)=set_log_dm_agec_rGN;
	   log_dm_sCT_ac_out(8)=log_dm_agec_sCT; log_dm_sCT_ac_out(1,7)=set_log_dm_agec_sCT;
       
       
       
       selpar_A50_cHL2_out(8)=selpar_A50_cHL2; selpar_A50_cHL2_out(1,7)=set_selpar_A50_cHL2;
       selpar_slope_cHL2_out(8)=selpar_slope_cHL2; selpar_slope_cHL2_out(1,7)=set_selpar_slope_cHL2;
       selpar_A50_cHL3_out(8)=selpar_A50_cHL3; selpar_A50_cHL3_out(1,7)=set_selpar_A50_cHL3;
       selpar_slope_cHL3_out(8)=selpar_slope_cHL3; selpar_slope_cHL3_out(1,7)=set_selpar_slope_cHL3;

	   selpar_A50_cOT2_out(8)=selpar_A50_cOT2; selpar_A50_cOT2_out(1,7)=set_selpar_A50_cOT2;
       selpar_A50_cOT3_out(8)=selpar_A50_cOT3; selpar_A50_cOT3_out(1,7)=set_selpar_A50_cOT3;
       selpar_slope_cOT2_out(8)=selpar_slope_cOT2; selpar_slope_cOT2_out(1,7)=set_selpar_slope_cOT2;
       selpar_A502_cOT2_out(8)=selpar_A502_cOT2; selpar_A502_cOT2_out(1,7)=set_selpar_A502_cOT2;
       selpar_slope2_cOT2_out(8)=selpar_slope2_cOT2; selpar_slope2_cOT2_out(1,7)=set_selpar_slope2_cOT2;	  
	   
       selpar_A50_rHB1_out(8)=selpar_A50_rHB1; selpar_A50_rHB1_out(1,7)=set_selpar_A50_rHB1;
       selpar_slope_rHB1_out(8)=selpar_slope_rHB1; selpar_slope_rHB1_out(1,7)=set_selpar_slope_rHB1;
       selpar_A50_rHB2_out(8)=selpar_A50_rHB2; selpar_A50_rHB2_out(1,7)=set_selpar_A50_rHB2;
       selpar_slope_rHB2_out(8)=selpar_slope_rHB2; selpar_slope_rHB2_out(1,7)=set_selpar_slope_rHB2;
       selpar_A50_rHB3_out(8)=selpar_A50_rHB3; selpar_A50_rHB3_out(1,7)=set_selpar_A50_rHB3;
       selpar_slope_rHB3_out(8)=selpar_slope_rHB3; selpar_slope_rHB3_out(1,7)=set_selpar_slope_rHB3;
	   
	   selpar_A50_rGN3_out(8)=selpar_A50_rGN3; selpar_A50_rGN3_out(1,7)=set_selpar_A50_rGN3;
       selpar_slope_rGN3_out(8)=selpar_slope_rGN3; selpar_slope_rGN3_out(1,7)=set_selpar_slope_rGN3;
       
       selpar_A50_sCT_out(8)=selpar_A50_sCT; selpar_A50_sCT_out(1,7)=set_selpar_A50_sCT;
       selpar_slope_sCT_out(8)=selpar_slope_sCT; selpar_slope_sCT_out(1,7)=set_selpar_slope_sCT;
       selpar_A502_sCT_out(8)=selpar_A502_sCT; selpar_A502_sCT_out(1,7)=set_selpar_A502_sCT;
       selpar_slope2_sCT_out(8)=selpar_slope2_sCT; selpar_slope2_sCT_out(1,7)=set_selpar_slope2_sCT;
       
	   selpar_age1logit_D_out(8)=selpar_age1logit_D; selpar_age1logit_D_out(1,7)=set_selpar_age1logit_D;

       log_q_cHL_out(8)=log_q_cpue_cHL; log_q_cHL_out(1,7)=set_log_q_cpue_cHL;
       log_q_rHB_out(8)=log_q_cpue_rHB; log_q_rHB_out(1,7)=set_log_q_cpue_rHB;
       
       log_q_sCT_out(8)=log_q_cpue_sCT; log_q_sCT_out(1,7)=set_log_q_cpue_sCT;
                     
       log_avg_F_cHL_out(8)=log_avg_F_L_cHL; log_avg_F_cHL_out(1,7)=set_log_avg_F_L_cHL;
	   log_avg_F_cOT_out(8)=log_avg_F_L_cOT; log_avg_F_cOT_out(1,7)=set_log_avg_F_L_cOT;
       log_avg_F_rHB_out(8)=log_avg_F_L_rHB; log_avg_F_rHB_out(1,7)=set_log_avg_F_L_rHB;
       log_avg_F_rGN_out(8)=log_avg_F_L_rGN; log_avg_F_rGN_out(1,7)=set_log_avg_F_L_rGN;       
       log_avg_F_cHL_D_out(8)=log_avg_F_D_cHL; log_avg_F_cHL_D_out(1,7)=set_log_avg_F_D_cHL;
       log_avg_F_rHB_D_out(8)=log_avg_F_D_rHB; log_avg_F_rHB_D_out(1,7)=set_log_avg_F_D_rHB;
       log_avg_F_rGN_D_out(8)=log_avg_F_D_rGN; log_avg_F_rGN_D_out(1,7)=set_log_avg_F_D_rGN;
        
       log_rec_dev_out(styr_rec_dev, endyr_rec_dev)=log_dev_rec;
       log_F_dev_cHL_out(styr_L_cHL,endyr_L_cHL)=log_dev_F_L_cHL;
	   log_F_dev_cOT_out(styr_L_cOT,endyr_L_cOT)=log_dev_F_L_cOT;
       log_F_dev_rHB_out(styr_L_rHB,endyr_L_rHB)=log_dev_F_L_rHB;
       log_F_dev_rGN_out(styr_L_rGN,endyr_L_rGN)=log_dev_F_L_rGN;
       log_F_dev_cHL_D_out(styr_D_cHL,endyr_D_cHL)=log_dev_F_D_cHL;
       log_F_dev_rHB_D_out(styr_D_rHB,endyr_D_rHB)=log_dev_F_D_rHB;
       log_F_dev_rGN_D_out(styr_D_rGN,endyr_D_rGN)=log_dev_F_D_rGN;
           
   #include "RedGrouper.cxx"   

  } 
  
 
