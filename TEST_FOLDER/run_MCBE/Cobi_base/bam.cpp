#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
  #include "admodel.h"          // Include AD class definitions
  #include "admb2r.cpp"         // Include S-compatible output functions (needs preceding)
  #include <time.h>
	time_t start,finish;
	long hour,minute,second;	
	double elapsed_time;
#ifdef DEBUG
  #include <chrono>
#endif
#include <admodel.h>
#ifdef USE_ADMB_CONTRIBS
#include <contrib.h>

#endif
  extern "C"  {
    void ad_boundf(int i);
  }
#include <bam.htp>

model_data::model_data(int argc,char * argv[]) : ad_comm(argc,argv)
{
  adstring tmpstring;
  tmpstring=adprogram_name + adstring(".dat");
  if (argc > 1)
  {
    int on=0;
    if ( (on=option_match(argc,argv,"-ind"))>-1)
    {
      if (on>argc-2 || argv[on+1][0] == '-')
      {
        cerr << "Invalid input data command line option"
                " -- ignored" << endl;
      }
      else
      {
        tmpstring = adstring(argv[on+1]);
      }
    }
  }
  global_datafile = new cifstream(tmpstring);
  if (!global_datafile)
  {
    cerr << "Error: Unable to allocate global_datafile in model_data constructor.";
    ad_exit(1);
  }
  if (!(*global_datafile))
  {
    delete global_datafile;
    global_datafile=NULL;
  }
cout << "Starting Beaufort Assessment Model" << endl;
cout << endl;
cout << "                BAM!" << endl;
cout << endl;
  styr.allocate("styr");
  endyr.allocate("endyr");
  styr_rec_dev.allocate("styr_rec_dev");
  endyr_rec_dev.allocate("endyr_rec_dev");
  endyr_rec_phase1.allocate("endyr_rec_phase1");
  endyr_rec_phase2.allocate("endyr_rec_phase2");
  endyr_selex_phase1_rGN.allocate("endyr_selex_phase1_rGN");
   nyrs=endyr-styr+1.;
   nyrs_rec=endyr_rec_dev-styr_rec_dev+1.;
  nages.allocate("nages");
  agebins.allocate(1,nages,"agebins");
  nages_agec.allocate("nages_agec");
  agebins_agec.allocate(1,nages_agec,"agebins_agec");
  nlenbins.allocate("nlenbins");
  lenbins_width.allocate("lenbins_width");
  lenbins.allocate(1,nlenbins,"lenbins");
  max_F_spr_msy.allocate("max_F_spr_msy");
  n_iter_spr.allocate("n_iter_spr");
		n_iter_msy=n_iter_spr; 
  styr_rec_spr.allocate("styr_rec_spr");
  endyr_rec_spr.allocate("endyr_rec_spr");
   nyrs_rec_spr=endyr_rec_spr-styr_rec_spr+1.;
  selpar_n_yrs_wgted.allocate("selpar_n_yrs_wgted");
  set_BiasCor.allocate("set_BiasCor");
  styr_L_cGN.allocate("styr_L_cGN");
  endyr_L_cGN.allocate("endyr_L_cGN");
  obs_L_cGN.allocate(styr_L_cGN,endyr_L_cGN,"obs_L_cGN");
  obs_cv_L_cGN.allocate(styr_L_cGN,endyr_L_cGN,"obs_cv_L_cGN");
  nyr_lenc_pool_cGN.allocate("nyr_lenc_pool_cGN");
  yrs_lenc_pool_cGN.allocate(1,nyr_lenc_pool_cGN,"yrs_lenc_pool_cGN");
  nfish_lenc_pool_cGN.allocate(1,nyr_lenc_pool_cGN,"nfish_lenc_pool_cGN");
  nyr_lenc_cGN.allocate("nyr_lenc_cGN");
  yrs_lenc_cGN.allocate(1,nyr_lenc_cGN,"yrs_lenc_cGN");
  nsamp_lenc_cGN.allocate(1,nyr_lenc_cGN,"nsamp_lenc_cGN");
  nfish_lenc_cGN.allocate(1,nyr_lenc_cGN,"nfish_lenc_cGN");
  obs_lenc_cGN.allocate(1,nyr_lenc_cGN,1,nlenbins,"obs_lenc_cGN");
  styr_cpue_rHB.allocate("styr_cpue_rHB");
  endyr_cpue_rHB.allocate("endyr_cpue_rHB");
  obs_cpue_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"obs_cpue_rHB");
  obs_cv_cpue_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"obs_cv_cpue_rHB");
  styr_L_rGN.allocate("styr_L_rGN");
  endyr_L_rGN.allocate("endyr_L_rGN");
  obs_L_rGN.allocate(styr_L_rGN,endyr_L_rGN,"obs_L_rGN");
  obs_cv_L_rGN.allocate(styr_L_rGN,endyr_L_rGN,"obs_cv_L_rGN");
  nyr_agec_rGN.allocate("nyr_agec_rGN");
  yrs_agec_rGN.allocate(1,nyr_agec_rGN,"yrs_agec_rGN");
  nsamp_agec_rGN.allocate(1,nyr_agec_rGN,"nsamp_agec_rGN");
  nfish_agec_rGN.allocate(1,nyr_agec_rGN,"nfish_agec_rGN");
  obs_agec_rGN.allocate(1,nyr_agec_rGN,1,nages_agec,"obs_agec_rGN");
  set_Linf.allocate(1,7,"set_Linf");
  set_K.allocate(1,7,"set_K");
  set_t0.allocate(1,7,"set_t0");
  set_len_cv.allocate(1,7,"set_len_cv");
  set_Linf_L.allocate(1,7,"set_Linf_L");
  set_K_L.allocate(1,7,"set_K_L");
  set_t0_L.allocate(1,7,"set_t0_L");
  set_len_cv_L.allocate(1,7,"set_len_cv_L");
  set_Linf_F.allocate(1,7,"set_Linf_F");
  set_K_F.allocate(1,7,"set_K_F");
  set_t0_F.allocate(1,7,"set_t0_F");
  set_len_cv_F.allocate(1,7,"set_len_cv_F");
  set_M_constant.allocate(1,7,"set_M_constant");
  set_steep.allocate(1,7,"set_steep");
  set_log_R0.allocate(1,7,"set_log_R0");
  set_R_autocorr.allocate(1,7,"set_R_autocorr");
  set_rec_sigma.allocate(1,7,"set_rec_sigma");
  set_log_dm_lenc_cGN.allocate(1,7,"set_log_dm_lenc_cGN");
  set_log_dm_agec_rGN.allocate(1,7,"set_log_dm_agec_rGN");
  set_selpar_A50_cGN1.allocate(1,7,"set_selpar_A50_cGN1");
  set_selpar_slope_cGN1.allocate(1,7,"set_selpar_slope_cGN1");
  set_selpar_A50_rGN1.allocate(1,7,"set_selpar_A50_rGN1");
  set_selpar_slope_rGN1.allocate(1,7,"set_selpar_slope_rGN1");
  set_selpar_A50_rGN2.allocate(1,7,"set_selpar_A50_rGN2");
  set_selpar_slope_rGN2.allocate(1,7,"set_selpar_slope_rGN2");
  set_log_q_cpue_rHB.allocate(1,7,"set_log_q_cpue_rHB");
  set_F_init.allocate(1,7,"set_F_init");
  set_log_avg_F_L_cGN.allocate(1,7,"set_log_avg_F_L_cGN");
  set_log_avg_F_L_rGN.allocate(1,7,"set_log_avg_F_L_rGN");
  set_log_dev_F_L_cGN.allocate(1,3,"set_log_dev_F_L_cGN");
  set_log_dev_F_L_rGN.allocate(1,3,"set_log_dev_F_L_rGN");
  set_log_dev_RWq.allocate(1,3,"set_log_dev_RWq");
  set_log_dev_rec.allocate(1,3,"set_log_dev_rec");
  set_log_dev_Nage.allocate(1,3,"set_log_dev_Nage");
  set_log_dev_vals_F_L_cGN.allocate(styr_L_cGN,endyr_L_cGN,"set_log_dev_vals_F_L_cGN");
  set_log_dev_vals_F_L_rGN.allocate(styr_L_rGN,endyr_L_rGN,"set_log_dev_vals_F_L_rGN");
  set_log_dev_vals_rec.allocate(styr_rec_dev,endyr_rec_dev,"set_log_dev_vals_rec");
  set_log_dev_vals_Nage.allocate(2,nages,"set_log_dev_vals_Nage");
  set_w_L.allocate("set_w_L");
  set_w_cpue_rHB.allocate("set_w_cpue_rHB");
  set_w_lenc_cGN.allocate("set_w_lenc_cGN");
  set_w_agec_rGN.allocate("set_w_agec_rGN");
  set_w_Nage_init.allocate("set_w_Nage_init");
  set_w_rec.allocate("set_w_rec");
  set_w_rec_early.allocate("set_w_rec_early");
  set_w_rec_end.allocate("set_w_rec_end");
  set_w_fullF.allocate("set_w_fullF");
  set_w_Ftune.allocate("set_w_Ftune");
  wgtpar_a.allocate("wgtpar_a");
  wgtpar_b.allocate("wgtpar_b");
  obs_maturity_f.allocate(1,nages,"obs_maturity_f");
  obs_prop_f.allocate(1,nages,"obs_prop_f");
  spawn_time_frac.allocate("spawn_time_frac");
  set_M.allocate(1,nages,"set_M");
  SR_switch.allocate("SR_switch");
  set_q_rate_phase.allocate("set_q_rate_phase");
  set_q_rate.allocate("set_q_rate");
  set_q_DD_phase.allocate("set_q_DD_phase");
  set_q_DD_beta.allocate("set_q_DD_beta");
  set_q_DD_beta_se.allocate("set_q_DD_beta_se");
  set_q_DD_stage.allocate("set_q_DD_stage");
  set_RWq_var.allocate("set_RWq_var");
  set_Ftune.allocate("set_Ftune");
  set_Ftune_yr.allocate("set_Ftune_yr");
  minSS_lenc_cGN.allocate("minSS_lenc_cGN");
  minSS_agec_rGN.allocate("minSS_agec_rGN");
  endyr_proj.allocate("endyr_proj");
  styr_regs.allocate("styr_regs");
  Fproj_switch.allocate("Fproj_switch");
  Fproj_mult.allocate("Fproj_mult");
   styr_proj=endyr+1;
  age_error.allocate(1,nages,1,nages,"age_error");
  end_of_data_file.allocate("end_of_data_file");
   if(end_of_data_file!=999)
   {
       cout << "*** WARNING: Data File NOT READ CORRECTLY ****" << endl;
       exit(0);  
   }
   else
   {cout << "Data File read correctly" << endl;} 
  if (global_datafile)
  {
    delete global_datafile;
    global_datafile = NULL;
  }
}

model_parameters::model_parameters(int sz,int argc,char * argv[]) : 
 model_data(argc,argv) , function_minimizer(sz)
{
  initializationfunction();
  const double Linf_LO=set_Linf(2); const double Linf_HI=set_Linf(3); const double Linf_PH=set_Linf(4);
  const double K_LO=set_K(2); const double K_HI=set_K(3); const double K_PH=set_K(4);
  const double t0_LO=set_t0(2); const double t0_HI=set_t0(3); const double t0_PH=set_t0(4);  
  const double len_cv_LO=set_len_cv(2); const double len_cv_HI=set_len_cv(3); const double len_cv_PH=set_len_cv(4); 
   
  const double Linf_L_LO=set_Linf_L(2); const double Linf_L_HI=set_Linf_L(3); const double Linf_L_PH=set_Linf_L(4);
  const double K_L_LO=set_K_L(2); const double K_L_HI=set_K_L(3); const double K_L_PH=set_K_L(4);
  const double t0_L_LO=set_t0_L(2); const double t0_L_HI=set_t0_L(3); const double t0_L_PH=set_t0_L(4);  
  const double len_cv_L_LO=set_len_cv_L(2); const double len_cv_L_HI=set_len_cv_L(3); const double len_cv_L_PH=set_len_cv_L(4);
  
  const double Linf_F_LO=set_Linf_F(2); const double Linf_F_HI=set_Linf_F(3); const double Linf_F_PH=set_Linf_F(4);
  const double K_F_LO=set_K_F(2); const double K_F_HI=set_K_F(3); const double K_F_PH=set_K_F(4);
  const double t0_F_LO=set_t0_F(2); const double t0_F_HI=set_t0_F(3); const double t0_F_PH=set_t0_F(4);  
  const double len_cv_F_LO=set_len_cv_F(2); const double len_cv_F_HI=set_len_cv_F(3); const double len_cv_F_PH=set_len_cv_F(4);
  
  const double M_constant_LO=set_M_constant(2); const double M_constant_HI=set_M_constant(3); const double M_constant_PH=set_M_constant(4);        
  const double steep_LO=set_steep(2); const double steep_HI=set_steep(3); const double steep_PH=set_steep(4);
  const double log_R0_LO=set_log_R0(2); const double log_R0_HI=set_log_R0(3); const double log_R0_PH=set_log_R0(4);
  const double R_autocorr_LO=set_R_autocorr(2); const double R_autocorr_HI=set_R_autocorr(3); const double R_autocorr_PH=set_R_autocorr(4);
  const double rec_sigma_LO=set_rec_sigma(2); const double rec_sigma_HI=set_rec_sigma(3); const double rec_sigma_PH=set_rec_sigma(4);
  
  const double log_dm_cGN_lc_LO=set_log_dm_lenc_cGN(2); const double log_dm_cGN_lc_HI=set_log_dm_lenc_cGN(3); const double log_dm_cGN_lc_PH=set_log_dm_lenc_cGN(4);
  //const double log_dm_cLL_lc_LO=set_log_dm_cLL_lc(2); const double log_dm_cLL_lc_HI=set_log_dm_cLL_lc(3); const double log_dm_cLL_lc_PH=set_log_dm_cLL_lc(4);
  //const double log_dm_rGN_lc_LO=set_log_dm_rGN_lc(2); const double log_dm_rGN_lc_HI=set_log_dm_rGN_lc(3); const double log_dm_rGN_lc_PH=set_log_dm_rGN_lc(4);
  const double log_dm_rGN_ac_LO=set_log_dm_agec_rGN(2); const double log_dm_rGN_ac_HI=set_log_dm_agec_rGN(3); const double log_dm_rGN_ac_PH=set_log_dm_agec_rGN(4);
  
  const double selpar_A50_cGN1_LO=set_selpar_A50_cGN1(2); const double selpar_A50_cGN1_HI=set_selpar_A50_cGN1(3); const double selpar_A50_cGN1_PH=set_selpar_A50_cGN1(4);
  const double selpar_slope_cGN1_LO=set_selpar_slope_cGN1(2); const double selpar_slope_cGN1_HI=set_selpar_slope_cGN1(3); const double selpar_slope_cGN1_PH=set_selpar_slope_cGN1(4);
  //const double selpar_A50_cGN2_LO=set_selpar_A50_cGN2(2); const double selpar_A50_cGN2_HI=set_selpar_A50_cGN2(3); const double selpar_A50_cGN2_PH=set_selpar_A50_cGN2(4);
  //const double selpar_slope_cGN2_LO=set_selpar_slope_cGN2(2); const double selpar_slope_cGN2_HI=set_selpar_slope_cGN2(3); const double selpar_slope_cGN2_PH=set_selpar_slope_cGN2(4);
  // const double selpar_A502_cGN2_LO=set_selpar_A502_cGN2(2); const double selpar_A502_cGN2_HI=set_selpar_A502_cGN2(3); const double selpar_A502_cGN2_PH=set_selpar_A502_cGN2(4);
  // const double selpar_slope2_cGN2_LO=set_selpar_slope2_cGN2(2); const double selpar_slope2_cGN2_HI=set_selpar_slope2_cGN2(3); const double selpar_slope2_cGN2_PH=set_selpar_slope2_cGN2(4); 
  //const double selpar_A50_cGN3_LO=set_selpar_A50_cGN3(2); const double selpar_A50_cGN3_HI=set_selpar_A50_cGN3(3); const double selpar_A50_cGN3_PH=set_selpar_A50_cGN3(4);
  //const double selpar_slope_cGN3_LO=set_selpar_slope_cGN3(2); const double selpar_slope_cGN3_HI=set_selpar_slope_cGN3(3); const double selpar_slope_cGN3_PH=set_selpar_slope_cGN3(4);  
  
  const double selpar_A50_rGN1_LO=set_selpar_A50_rGN1(2); const double selpar_A50_rGN1_HI=set_selpar_A50_rGN1(3); const double selpar_A50_rGN1_PH=set_selpar_A50_rGN1(4);
  const double selpar_slope_rGN1_LO=set_selpar_slope_rGN1(2); const double selpar_slope_rGN1_HI=set_selpar_slope_rGN1(3); const double selpar_slope_rGN1_PH=set_selpar_slope_rGN1(4);
  const double selpar_A50_rGN2_LO=set_selpar_A50_rGN2(2); const double selpar_A50_rGN2_HI=set_selpar_A50_rGN2(3); const double selpar_A50_rGN2_PH=set_selpar_A50_rGN2(4);
  const double selpar_slope_rGN2_LO=set_selpar_slope_rGN2(2); const double selpar_slope_rGN2_HI=set_selpar_slope_rGN2(3); const double selpar_slope_rGN2_PH=set_selpar_slope_rGN2(4);
  // const double selpar_A502_rGN2_LO=set_selpar_A502_rGN2(2); const double selpar_A502_rGN2_HI=set_selpar_A502_rGN2(3); const double selpar_A502_rGN2_PH=set_selpar_A502_rGN2(4);
  // const double selpar_slope2_rGN2_LO=set_selpar_slope2_rGN2(2); const double selpar_slope2_rGN2_HI=set_selpar_slope2_rGN2(3); const double selpar_slope2_rGN2_PH=set_selpar_slope2_rGN2(4);  
  //const double selpar_A50_rGN3_LO=set_selpar_A50_rGN3(2); const double selpar_A50_rGN3_HI=set_selpar_A50_rGN3(3); const double selpar_A50_rGN3_PH=set_selpar_A50_rGN3(4);
  //const double selpar_slope_rGN3_LO=set_selpar_slope_rGN3(2); const double selpar_slope_rGN3_HI=set_selpar_slope_rGN3(3); const double selpar_slope_rGN3_PH=set_selpar_slope_rGN3(4);  
  
  
  //const double log_q_cGN_LO=set_log_q_cGN(2); const double log_q_cGN_HI=set_log_q_cGN(3); const double log_q_cGN_PH=set_log_q_cGN(4);
  //const double log_q_cLL_LO=set_log_q_cLL(2); const double log_q_cLL_HI=set_log_q_cLL(3); const double log_q_cLL_PH=set_log_q_cLL(4);  
  const double log_q_rHB_LO=set_log_q_cpue_rHB(2); const double log_q_rHB_HI=set_log_q_cpue_rHB(3); const double log_q_rHB_PH=set_log_q_cpue_rHB(4);
  
  const double F_init_LO=set_F_init(2); const double F_init_HI=set_F_init(3); const double F_init_PH=set_F_init(4);
  const double log_avg_F_cGN_LO=set_log_avg_F_L_cGN(2); const double log_avg_F_cGN_HI=set_log_avg_F_L_cGN(3); const double log_avg_F_cGN_PH=set_log_avg_F_L_cGN(4);
  //const double log_avg_F_cLL_LO=set_log_avg_F_cLL(2); const double log_avg_F_cLL_HI=set_log_avg_F_cLL(3); const double log_avg_F_cLL_PH=set_log_avg_F_cLL(4);
  const double log_avg_F_rGN_LO=set_log_avg_F_L_rGN(2); const double log_avg_F_rGN_HI=set_log_avg_F_L_rGN(3); const double log_avg_F_rGN_PH=set_log_avg_F_L_rGN(4); 
  
  //-dev vectors-----------------------------------------------------------------------------------------------------------  
  const double log_F_dev_cGN_LO=set_log_dev_F_L_cGN(1); const double log_F_dev_cGN_HI=set_log_dev_F_L_cGN(2); const double log_F_dev_cGN_PH=set_log_dev_F_L_cGN(3);  
  //const double log_F_dev_cLL_LO=set_log_F_dev_cLL(1); const double log_F_dev_cLL_HI=set_log_F_dev_cLL(2); const double log_F_dev_cLL_PH=set_log_F_dev_cLL(3);      
  const double log_F_dev_rGN_LO=set_log_dev_F_L_rGN(1); const double log_F_dev_rGN_HI=set_log_dev_F_L_rGN(2); const double log_F_dev_rGN_PH=set_log_dev_F_L_rGN(3);     
  
  const double log_RWq_LO=set_log_dev_RWq(1); const double log_RWq_HI=set_log_dev_RWq(2); const double log_RWq_PH=set_log_dev_RWq(3);  
  
  const double log_rec_dev_LO=set_log_dev_rec(1); const double log_rec_dev_HI=set_log_dev_rec(2); const double log_rec_dev_PH=set_log_dev_rec(3);          
  const double log_Nage_dev_LO=set_log_dev_Nage(1); const double log_Nage_dev_HI=set_log_dev_Nage(2); const double log_Nage_dev_PH=set_log_dev_Nage(3);          
  
  Linf.allocate(Linf_LO,Linf_HI,Linf_PH,"Linf");
  K.allocate(K_LO,K_HI,K_PH,"K");
  t0.allocate(t0_LO,t0_HI,t0_PH,"t0");
  len_cv_val.allocate(len_cv_LO,len_cv_HI,len_cv_PH,"len_cv_val");
  Linf_out.allocate(1,8,"Linf_out");
  #ifndef NO_AD_INITIALIZE
    Linf_out.initialize();
  #endif
  K_out.allocate(1,8,"K_out");
  #ifndef NO_AD_INITIALIZE
    K_out.initialize();
  #endif
  t0_out.allocate(1,8,"t0_out");
  #ifndef NO_AD_INITIALIZE
    t0_out.initialize();
  #endif
  len_cv_val_out.allocate(1,8,"len_cv_val_out");
  #ifndef NO_AD_INITIALIZE
    len_cv_val_out.initialize();
  #endif
  meanlen_TL.allocate(1,nages,"meanlen_TL");
  #ifndef NO_AD_INITIALIZE
    meanlen_TL.initialize();
  #endif
  wgt_g.allocate(1,nages,"wgt_g");
  #ifndef NO_AD_INITIALIZE
    wgt_g.initialize();
  #endif
  wgt_kg.allocate(1,nages,"wgt_kg");
  #ifndef NO_AD_INITIALIZE
    wgt_kg.initialize();
  #endif
  wgt_mt.allocate(1,nages,"wgt_mt");
  #ifndef NO_AD_INITIALIZE
    wgt_mt.initialize();
  #endif
  wgt_klb.allocate(1,nages,"wgt_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_klb.initialize();
  #endif
  wgt_lb.allocate(1,nages,"wgt_lb");
  #ifndef NO_AD_INITIALIZE
    wgt_lb.initialize();
  #endif
  Linf_L.allocate(Linf_L_LO,Linf_L_HI,Linf_L_PH,"Linf_L");
  K_L.allocate(K_L_LO,K_L_HI,K_L_PH,"K_L");
  t0_L.allocate(t0_L_LO,t0_L_HI,t0_L_PH,"t0_L");
  len_cv_val_L.allocate(len_cv_L_LO,len_cv_L_HI,len_cv_L_PH,"len_cv_val_L");
  Linf_L_out.allocate(1,8,"Linf_L_out");
  #ifndef NO_AD_INITIALIZE
    Linf_L_out.initialize();
  #endif
  K_L_out.allocate(1,8,"K_L_out");
  #ifndef NO_AD_INITIALIZE
    K_L_out.initialize();
  #endif
  t0_L_out.allocate(1,8,"t0_L_out");
  #ifndef NO_AD_INITIALIZE
    t0_L_out.initialize();
  #endif
  len_cv_val_L_out.allocate(1,8,"len_cv_val_L_out");
  #ifndef NO_AD_INITIALIZE
    len_cv_val_L_out.initialize();
  #endif
  meanlen_TL_L.allocate(1,nages,"meanlen_TL_L");
  #ifndef NO_AD_INITIALIZE
    meanlen_TL_L.initialize();
  #endif
  wgt_g_L.allocate(1,nages,"wgt_g_L");
  #ifndef NO_AD_INITIALIZE
    wgt_g_L.initialize();
  #endif
  wgt_kg_L.allocate(1,nages,"wgt_kg_L");
  #ifndef NO_AD_INITIALIZE
    wgt_kg_L.initialize();
  #endif
  wgt_mt_L.allocate(1,nages,"wgt_mt_L");
  #ifndef NO_AD_INITIALIZE
    wgt_mt_L.initialize();
  #endif
  wgt_klb_L.allocate(1,nages,"wgt_klb_L");
  #ifndef NO_AD_INITIALIZE
    wgt_klb_L.initialize();
  #endif
  wgt_lb_L.allocate(1,nages,"wgt_lb_L");
  #ifndef NO_AD_INITIALIZE
    wgt_lb_L.initialize();
  #endif
  wgt_klb_gut_L.allocate(1,nages,"wgt_klb_gut_L");
  #ifndef NO_AD_INITIALIZE
    wgt_klb_gut_L.initialize();
  #endif
  wgt_lb_gut_L.allocate(1,nages,"wgt_lb_gut_L");
  #ifndef NO_AD_INITIALIZE
    wgt_lb_gut_L.initialize();
  #endif
  Linf_F.allocate(Linf_F_LO,Linf_F_HI,Linf_F_PH,"Linf_F");
  K_F.allocate(K_F_LO,K_F_HI,K_F_PH,"K_F");
  t0_F.allocate(t0_F_LO,t0_F_HI,t0_F_PH,"t0_F");
  len_cv_val_F.allocate(len_cv_F_LO,len_cv_F_HI,len_cv_F_PH,"len_cv_val_F");
  Linf_F_out.allocate(1,8,"Linf_F_out");
  #ifndef NO_AD_INITIALIZE
    Linf_F_out.initialize();
  #endif
  K_F_out.allocate(1,8,"K_F_out");
  #ifndef NO_AD_INITIALIZE
    K_F_out.initialize();
  #endif
  t0_F_out.allocate(1,8,"t0_F_out");
  #ifndef NO_AD_INITIALIZE
    t0_F_out.initialize();
  #endif
  len_cv_val_F_out.allocate(1,8,"len_cv_val_F_out");
  #ifndef NO_AD_INITIALIZE
    len_cv_val_F_out.initialize();
  #endif
  meanlen_TL_F.allocate(1,nages,"meanlen_TL_F");
  #ifndef NO_AD_INITIALIZE
    meanlen_TL_F.initialize();
  #endif
  wgt_g_F.allocate(1,nages,"wgt_g_F");
  #ifndef NO_AD_INITIALIZE
    wgt_g_F.initialize();
  #endif
  wgt_kg_F.allocate(1,nages,"wgt_kg_F");
  #ifndef NO_AD_INITIALIZE
    wgt_kg_F.initialize();
  #endif
  wgt_mt_F.allocate(1,nages,"wgt_mt_F");
  #ifndef NO_AD_INITIALIZE
    wgt_mt_F.initialize();
  #endif
  wgt_klb_F.allocate(1,nages,"wgt_klb_F");
  #ifndef NO_AD_INITIALIZE
    wgt_klb_F.initialize();
  #endif
  wgt_lb_F.allocate(1,nages,"wgt_lb_F");
  #ifndef NO_AD_INITIALIZE
    wgt_lb_F.initialize();
  #endif
  len_cGN_mm.allocate(styr,endyr,1,nages,"len_cGN_mm");
  #ifndef NO_AD_INITIALIZE
    len_cGN_mm.initialize();
  #endif
  wholewgt_cGN_klb.allocate(styr,endyr,1,nages,"wholewgt_cGN_klb");
  #ifndef NO_AD_INITIALIZE
    wholewgt_cGN_klb.initialize();
  #endif
  len_rHB_mm.allocate(styr,endyr,1,nages,"len_rHB_mm");
  #ifndef NO_AD_INITIALIZE
    len_rHB_mm.initialize();
  #endif
  wholewgt_rHB_klb.allocate(styr,endyr,1,nages,"wholewgt_rHB_klb");
  #ifndef NO_AD_INITIALIZE
    wholewgt_rHB_klb.initialize();
  #endif
  len_rGN_mm.allocate(styr,endyr,1,nages,"len_rGN_mm");
  #ifndef NO_AD_INITIALIZE
    len_rGN_mm.initialize();
  #endif
  wholewgt_rGN_klb.allocate(styr,endyr,1,nages,"wholewgt_rGN_klb");
  #ifndef NO_AD_INITIALIZE
    wholewgt_rGN_klb.initialize();
  #endif
  lenprob.allocate(1,nages,1,nlenbins,"lenprob");
  #ifndef NO_AD_INITIALIZE
    lenprob.initialize();
  #endif
  zscore_len.allocate("zscore_len");
  #ifndef NO_AD_INITIALIZE
  zscore_len.initialize();
  #endif
  cprob_lenvec.allocate(1,nlenbins,"cprob_lenvec");
  #ifndef NO_AD_INITIALIZE
    cprob_lenvec.initialize();
  #endif
  zscore_lzero.allocate("zscore_lzero");
  #ifndef NO_AD_INITIALIZE
  zscore_lzero.initialize();
  #endif
  cprob_lzero.allocate("cprob_lzero");
  #ifndef NO_AD_INITIALIZE
  cprob_lzero.initialize();
  #endif
  lenprob_L.allocate(1,nages,1,nlenbins,"lenprob_L");
  #ifndef NO_AD_INITIALIZE
    lenprob_L.initialize();
  #endif
  zscore_len_L.allocate("zscore_len_L");
  #ifndef NO_AD_INITIALIZE
  zscore_len_L.initialize();
  #endif
  cprob_lenvec_L.allocate(1,nlenbins,"cprob_lenvec_L");
  #ifndef NO_AD_INITIALIZE
    cprob_lenvec_L.initialize();
  #endif
  zscore_lzero_L.allocate("zscore_lzero_L");
  #ifndef NO_AD_INITIALIZE
  zscore_lzero_L.initialize();
  #endif
  cprob_lzero_L.allocate("cprob_lzero_L");
  #ifndef NO_AD_INITIALIZE
  cprob_lzero_L.initialize();
  #endif
  lenprob_F.allocate(1,nages,1,nlenbins,"lenprob_F");
  #ifndef NO_AD_INITIALIZE
    lenprob_F.initialize();
  #endif
  zscore_len_F.allocate("zscore_len_F");
  #ifndef NO_AD_INITIALIZE
  zscore_len_F.initialize();
  #endif
  cprob_lenvec_F.allocate(1,nlenbins,"cprob_lenvec_F");
  #ifndef NO_AD_INITIALIZE
    cprob_lenvec_F.initialize();
  #endif
  zscore_lzero_F.allocate("zscore_lzero_F");
  #ifndef NO_AD_INITIALIZE
  zscore_lzero_F.initialize();
  #endif
  cprob_lzero_F.allocate("cprob_lzero_F");
  #ifndef NO_AD_INITIALIZE
  cprob_lzero_F.initialize();
  #endif
  lenprob_cGN.allocate(1,nages,1,nlenbins,"lenprob_cGN");
  #ifndef NO_AD_INITIALIZE
    lenprob_cGN.initialize();
  #endif
  lenprob_rHB.allocate(1,nages,1,nlenbins,"lenprob_rHB");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB.initialize();
  #endif
  lenprob_rGN.allocate(1,nages,1,nlenbins,"lenprob_rGN");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN.initialize();
  #endif
  len_sd.allocate(1,nages,"len_sd");
  #ifndef NO_AD_INITIALIZE
    len_sd.initialize();
  #endif
  len_cv.allocate(1,nages,"len_cv");
  #ifndef NO_AD_INITIALIZE
    len_cv.initialize();
  #endif
  len_sd_L.allocate(1,nages,"len_sd_L");
  #ifndef NO_AD_INITIALIZE
    len_sd_L.initialize();
  #endif
  len_cv_L.allocate(1,nages,"len_cv_L");
  #ifndef NO_AD_INITIALIZE
    len_cv_L.initialize();
  #endif
  len_sd_F.allocate(1,nages,"len_sd_F");
  #ifndef NO_AD_INITIALIZE
    len_sd_F.initialize();
  #endif
  len_cv_F.allocate(1,nages,"len_cv_F");
  #ifndef NO_AD_INITIALIZE
    len_cv_F.initialize();
  #endif
  pred_cGN_lenc.allocate(1,nyr_lenc_cGN,1,nlenbins,"pred_cGN_lenc");
  #ifndef NO_AD_INITIALIZE
    pred_cGN_lenc.initialize();
  #endif
  pred_cGN_lenc_yr.allocate(1,nyr_lenc_pool_cGN,1,nlenbins,"pred_cGN_lenc_yr");
  #ifndef NO_AD_INITIALIZE
    pred_cGN_lenc_yr.initialize();
  #endif
  pred_rGN_agec.allocate(1,nyr_agec_rGN,1,nages_agec,"pred_rGN_agec");
  #ifndef NO_AD_INITIALIZE
    pred_rGN_agec.initialize();
  #endif
  pred_rGN_agec_allages.allocate(1,nyr_agec_rGN,1,nages,"pred_rGN_agec_allages");
  #ifndef NO_AD_INITIALIZE
    pred_rGN_agec_allages.initialize();
  #endif
  ErrorFree_rGN_agec.allocate(1,nyr_agec_rGN,1,nages,"ErrorFree_rGN_agec");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_rGN_agec.initialize();
  #endif
  nsamp_cGN_lenc_allyr.allocate(styr,endyr,"nsamp_cGN_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_cGN_lenc_allyr.initialize();
  #endif
  nsamp_rGN_agec_allyr.allocate(styr,endyr,"nsamp_rGN_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_rGN_agec_allyr.initialize();
  #endif
  nfish_cGN_lenc_allyr.allocate(styr,endyr,"nfish_cGN_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_cGN_lenc_allyr.initialize();
  #endif
  nfish_rGN_agec_allyr.allocate(styr,endyr,"nfish_rGN_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_rGN_agec_allyr.initialize();
  #endif
  neff_cGN_lenc_allyr.allocate(styr,endyr,"neff_cGN_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    neff_cGN_lenc_allyr.initialize();
  #endif
  neff_rGN_agec_allyr.allocate(styr,endyr,"neff_rGN_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    neff_rGN_agec_allyr.initialize();
  #endif
  N.allocate(styr,endyr+1,1,nages,"N");
  #ifndef NO_AD_INITIALIZE
    N.initialize();
  #endif
  N_mdyr.allocate(styr,endyr,1,nages,"N_mdyr");
  #ifndef NO_AD_INITIALIZE
    N_mdyr.initialize();
  #endif
  N_spawn.allocate(styr,endyr,1,nages,"N_spawn");
  #ifndef NO_AD_INITIALIZE
    N_spawn.initialize();
  #endif
  log_dev_Nage.allocate(2,nages,log_Nage_dev_LO,log_Nage_dev_HI,log_Nage_dev_PH,"log_dev_Nage");
  log_Nage_dev_output.allocate(1,nages,"log_Nage_dev_output");
  #ifndef NO_AD_INITIALIZE
    log_Nage_dev_output.initialize();
  #endif
  B.allocate(styr,endyr+1,1,nages,"B");
  #ifndef NO_AD_INITIALIZE
    B.initialize();
  #endif
  totB.allocate(styr,endyr+1,"totB");
  #ifndef NO_AD_INITIALIZE
    totB.initialize();
  #endif
  totN.allocate(styr,endyr+1,"totN");
  #ifndef NO_AD_INITIALIZE
    totN.initialize();
  #endif
  SSB.allocate(styr,endyr,"SSB");
  #ifndef NO_AD_INITIALIZE
    SSB.initialize();
  #endif
  SSB_knum.allocate(styr,endyr,"SSB_knum");
  #ifndef NO_AD_INITIALIZE
    SSB_knum.initialize();
  #endif
  rec.allocate(styr,endyr+1,"rec");
  #ifndef NO_AD_INITIALIZE
    rec.initialize();
  #endif
  prop_f.allocate(1,nages,"prop_f");
  #ifndef NO_AD_INITIALIZE
    prop_f.initialize();
  #endif
  maturity_f.allocate(1,nages,"maturity_f");
  #ifndef NO_AD_INITIALIZE
    maturity_f.initialize();
  #endif
  reprod.allocate(1,nages,"reprod");
  #ifndef NO_AD_INITIALIZE
    reprod.initialize();
  #endif
  reprodknum.allocate(1,nages,"reprodknum");
  #ifndef NO_AD_INITIALIZE
    reprodknum.initialize();
  #endif
  log_R0.allocate(log_R0_LO,log_R0_HI,log_R0_PH,"log_R0");
  log_R0_out.allocate(1,8,"log_R0_out");
  #ifndef NO_AD_INITIALIZE
    log_R0_out.initialize();
  #endif
  R0.allocate("R0");
  #ifndef NO_AD_INITIALIZE
  R0.initialize();
  #endif
  steep.allocate(steep_LO,steep_HI,steep_PH,"steep");
  steep_out.allocate(1,8,"steep_out");
  #ifndef NO_AD_INITIALIZE
    steep_out.initialize();
  #endif
  rec_sigma.allocate(rec_sigma_LO,rec_sigma_HI,rec_sigma_PH,"rec_sigma");
  rec_sigma_out.allocate(1,8,"rec_sigma_out");
  #ifndef NO_AD_INITIALIZE
    rec_sigma_out.initialize();
  #endif
  R_autocorr.allocate(R_autocorr_LO,R_autocorr_HI,R_autocorr_PH,"R_autocorr");
  R_autocorr_out.allocate(1,8,"R_autocorr_out");
  #ifndef NO_AD_INITIALIZE
    R_autocorr_out.initialize();
  #endif
  rec_sigma_sq.allocate("rec_sigma_sq");
  #ifndef NO_AD_INITIALIZE
  rec_sigma_sq.initialize();
  #endif
  rec_logL_add.allocate("rec_logL_add");
  #ifndef NO_AD_INITIALIZE
  rec_logL_add.initialize();
  #endif
  log_dev_rec.allocate(styr_rec_dev,endyr_rec_dev,log_rec_dev_LO,log_rec_dev_HI,log_rec_dev_PH,"log_dev_rec");
  log_rec_dev_output.allocate(styr,endyr+1,"log_rec_dev_output");
  #ifndef NO_AD_INITIALIZE
    log_rec_dev_output.initialize();
  #endif
  log_rec_dev_out.allocate(styr_rec_dev,endyr_rec_dev,"log_rec_dev_out");
  #ifndef NO_AD_INITIALIZE
    log_rec_dev_out.initialize();
  #endif
  var_rec_dev.allocate("var_rec_dev");
  #ifndef NO_AD_INITIALIZE
  var_rec_dev.initialize();
  #endif
  sigma_rec_dev.allocate("sigma_rec_dev");
  #ifndef NO_AD_INITIALIZE
  sigma_rec_dev.initialize();
  #endif
  BiasCor.allocate("BiasCor");
  #ifndef NO_AD_INITIALIZE
  BiasCor.initialize();
  #endif
  S0.allocate("S0");
  #ifndef NO_AD_INITIALIZE
  S0.initialize();
  #endif
  B0.allocate("B0");
  #ifndef NO_AD_INITIALIZE
  B0.initialize();
  #endif
  R1.allocate("R1");
  #ifndef NO_AD_INITIALIZE
  R1.initialize();
  #endif
  R_virgin.allocate("R_virgin");
  #ifndef NO_AD_INITIALIZE
  R_virgin.initialize();
  #endif
  SdS0.allocate(styr,endyr,"SdS0");
  #ifndef NO_AD_INITIALIZE
    SdS0.initialize();
  #endif
  log_dm_lenc_cGN.allocate(log_dm_cGN_lc_LO,log_dm_cGN_lc_HI,log_dm_cGN_lc_PH,"log_dm_lenc_cGN");
  log_dm_agec_rGN.allocate(log_dm_rGN_ac_LO,log_dm_rGN_ac_HI,log_dm_rGN_ac_PH,"log_dm_agec_rGN");
  log_dm_cGN_lc_out.allocate(1,8,"log_dm_cGN_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_cGN_lc_out.initialize();
  #endif
  log_dm_rGN_ac_out.allocate(1,8,"log_dm_rGN_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_rGN_ac_out.initialize();
  #endif
  sel_cGN.allocate(styr,endyr,1,nages,"sel_cGN");
  #ifndef NO_AD_INITIALIZE
    sel_cGN.initialize();
  #endif
  sel_cGN_vec.allocate(1,nages,"sel_cGN_vec");
  #ifndef NO_AD_INITIALIZE
    sel_cGN_vec.initialize();
  #endif
  selpar_A50_cGN1.allocate(selpar_A50_cGN1_LO,selpar_A50_cGN1_HI,selpar_A50_cGN1_PH,"selpar_A50_cGN1");
  selpar_slope_cGN1.allocate(selpar_slope_cGN1_LO,selpar_slope_cGN1_HI,selpar_slope_cGN1_PH,"selpar_slope_cGN1");
  selpar_A50_cGN1_out.allocate(1,8,"selpar_A50_cGN1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cGN1_out.initialize();
  #endif
  selpar_slope_cGN1_out.allocate(1,8,"selpar_slope_cGN1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cGN1_out.initialize();
  #endif
  sel_rHB.allocate(styr,endyr,1,nages,"sel_rHB");
  #ifndef NO_AD_INITIALIZE
    sel_rHB.initialize();
  #endif
  sel_rHB_block1.allocate(1,nages,"sel_rHB_block1");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_block1.initialize();
  #endif
  sel_rHB_block2.allocate(1,nages,"sel_rHB_block2");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_block2.initialize();
  #endif
  sel_rGN.allocate(styr,endyr,1,nages,"sel_rGN");
  #ifndef NO_AD_INITIALIZE
    sel_rGN.initialize();
  #endif
  sel_rGN_block1.allocate(1,nages,"sel_rGN_block1");
  #ifndef NO_AD_INITIALIZE
    sel_rGN_block1.initialize();
  #endif
  sel_rGN_block2.allocate(1,nages,"sel_rGN_block2");
  #ifndef NO_AD_INITIALIZE
    sel_rGN_block2.initialize();
  #endif
  selpar_A50_rGN1.allocate(selpar_A50_rGN1_LO,selpar_A50_rGN1_HI,selpar_A50_rGN1_PH,"selpar_A50_rGN1");
  selpar_slope_rGN1.allocate(selpar_slope_rGN1_LO,selpar_slope_rGN1_HI,selpar_slope_rGN1_PH,"selpar_slope_rGN1");
  selpar_A50_rGN2.allocate(selpar_A50_rGN2_LO,selpar_A50_rGN2_HI,selpar_A50_rGN2_PH,"selpar_A50_rGN2");
  selpar_slope_rGN2.allocate(selpar_slope_rGN2_LO,selpar_slope_rGN2_HI,selpar_slope_rGN2_PH,"selpar_slope_rGN2");
  selpar_A50_rGN1_out.allocate(1,8,"selpar_A50_rGN1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rGN1_out.initialize();
  #endif
  selpar_slope_rGN1_out.allocate(1,8,"selpar_slope_rGN1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rGN1_out.initialize();
  #endif
  selpar_A50_rGN2_out.allocate(1,8,"selpar_A50_rGN2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rGN2_out.initialize();
  #endif
  selpar_slope_rGN2_out.allocate(1,8,"selpar_slope_rGN2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rGN2_out.initialize();
  #endif
  sel_wgted_L.allocate(1,nages,"sel_wgted_L");
  #ifndef NO_AD_INITIALIZE
    sel_wgted_L.initialize();
  #endif
  sel_wgted_tot.allocate(1,nages,"sel_wgted_tot");
  #ifndef NO_AD_INITIALIZE
    sel_wgted_tot.initialize();
  #endif
  pred_rHB_cpue.allocate(styr_cpue_rHB,endyr_cpue_rHB,"pred_rHB_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_rHB_cpue.initialize();
  #endif
  N_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,1,nages,"N_rHB");
  #ifndef NO_AD_INITIALIZE
    N_rHB.initialize();
  #endif
  log_q_cpue_rHB.allocate(log_q_rHB_LO,log_q_rHB_HI,log_q_rHB_PH,"log_q_cpue_rHB");
  log_q_rHB_out.allocate(1,8,"log_q_rHB_out");
  #ifndef NO_AD_INITIALIZE
    log_q_rHB_out.initialize();
  #endif
  q_rate.allocate("q_rate");
  #ifndef NO_AD_INITIALIZE
  q_rate.initialize();
  #endif
  q_rate_fcn_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"q_rate_fcn_rHB");
  #ifndef NO_AD_INITIALIZE
    q_rate_fcn_rHB.initialize();
  #endif
  q_DD_beta.allocate("q_DD_beta");
  #ifndef NO_AD_INITIALIZE
  q_DD_beta.initialize();
  #endif
  q_DD_fcn.allocate(styr,endyr,"q_DD_fcn");
  #ifndef NO_AD_INITIALIZE
    q_DD_fcn.initialize();
  #endif
  B0_q_DD.allocate("B0_q_DD");
  #ifndef NO_AD_INITIALIZE
  B0_q_DD.initialize();
  #endif
  B_q_DD.allocate(styr,endyr,"B_q_DD");
  #ifndef NO_AD_INITIALIZE
    B_q_DD.initialize();
  #endif
  q_RW_log_dev_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB-1,log_RWq_LO,log_RWq_HI,log_RWq_PH,"q_RW_log_dev_rHB");
  q_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"q_rHB");
  #ifndef NO_AD_INITIALIZE
    q_rHB.initialize();
  #endif
  L_cGN_num.allocate(styr,endyr,1,nages,"L_cGN_num");
  #ifndef NO_AD_INITIALIZE
    L_cGN_num.initialize();
  #endif
  L_cGN_klb.allocate(styr,endyr,1,nages,"L_cGN_klb");
  #ifndef NO_AD_INITIALIZE
    L_cGN_klb.initialize();
  #endif
  pred_cGN_L_knum.allocate(styr,endyr,"pred_cGN_L_knum");
  #ifndef NO_AD_INITIALIZE
    pred_cGN_L_knum.initialize();
  #endif
  pred_cGN_L_klb.allocate(styr,endyr,"pred_cGN_L_klb");
  #ifndef NO_AD_INITIALIZE
    pred_cGN_L_klb.initialize();
  #endif
  L_rGN_num.allocate(styr,endyr,1,nages,"L_rGN_num");
  #ifndef NO_AD_INITIALIZE
    L_rGN_num.initialize();
  #endif
  L_rGN_klb.allocate(styr,endyr,1,nages,"L_rGN_klb");
  #ifndef NO_AD_INITIALIZE
    L_rGN_klb.initialize();
  #endif
  pred_rGN_L_knum.allocate(styr,endyr,"pred_rGN_L_knum");
  #ifndef NO_AD_INITIALIZE
    pred_rGN_L_knum.initialize();
  #endif
  pred_rGN_L_klb.allocate(styr,endyr,"pred_rGN_L_klb");
  #ifndef NO_AD_INITIALIZE
    pred_rGN_L_klb.initialize();
  #endif
  L_total_num.allocate(styr,endyr,1,nages,"L_total_num");
  #ifndef NO_AD_INITIALIZE
    L_total_num.initialize();
  #endif
  L_total_klb.allocate(styr,endyr,1,nages,"L_total_klb");
  #ifndef NO_AD_INITIALIZE
    L_total_klb.initialize();
  #endif
  L_total_knum_yr.allocate(styr,endyr,"L_total_knum_yr");
  #ifndef NO_AD_INITIALIZE
    L_total_knum_yr.initialize();
  #endif
  L_total_klb_yr.allocate(styr,endyr,"L_total_klb_yr");
  #ifndef NO_AD_INITIALIZE
    L_total_klb_yr.initialize();
  #endif
  F_cGN_prop.allocate("F_cGN_prop");
  #ifndef NO_AD_INITIALIZE
  F_cGN_prop.initialize();
  #endif
  F_rGN_prop.allocate("F_rGN_prop");
  #ifndef NO_AD_INITIALIZE
  F_rGN_prop.initialize();
  #endif
  F_init_cGN_prop.allocate("F_init_cGN_prop");
  #ifndef NO_AD_INITIALIZE
  F_init_cGN_prop.initialize();
  #endif
  F_init_rGN_prop.allocate("F_init_rGN_prop");
  #ifndef NO_AD_INITIALIZE
  F_init_rGN_prop.initialize();
  #endif
  F_temp_sum.allocate("F_temp_sum");
  #ifndef NO_AD_INITIALIZE
  F_temp_sum.initialize();
  #endif
  F_end.allocate(1,nages,"F_end");
  #ifndef NO_AD_INITIALIZE
    F_end.initialize();
  #endif
  F_end_L.allocate(1,nages,"F_end_L");
  #ifndef NO_AD_INITIALIZE
    F_end_L.initialize();
  #endif
  F_end_apex.allocate("F_end_apex");
  #ifndef NO_AD_INITIALIZE
  F_end_apex.initialize();
  #endif
  SSB_msy_out.allocate("SSB_msy_out");
  #ifndef NO_AD_INITIALIZE
  SSB_msy_out.initialize();
  #endif
  F_msy_out.allocate("F_msy_out");
  #ifndef NO_AD_INITIALIZE
  F_msy_out.initialize();
  #endif
  msy_klb_out.allocate("msy_klb_out");
  #ifndef NO_AD_INITIALIZE
  msy_klb_out.initialize();
  #endif
  msy_knum_out.allocate("msy_knum_out");
  #ifndef NO_AD_INITIALIZE
  msy_knum_out.initialize();
  #endif
  B_msy_out.allocate("B_msy_out");
  #ifndef NO_AD_INITIALIZE
  B_msy_out.initialize();
  #endif
  R_msy_out.allocate("R_msy_out");
  #ifndef NO_AD_INITIALIZE
  R_msy_out.initialize();
  #endif
  spr_msy_out.allocate("spr_msy_out");
  #ifndef NO_AD_INITIALIZE
  spr_msy_out.initialize();
  #endif
  F20_dum.allocate("F20_dum");
  #ifndef NO_AD_INITIALIZE
  F20_dum.initialize();
  #endif
  F30_dum.allocate("F30_dum");
  #ifndef NO_AD_INITIALIZE
  F30_dum.initialize();
  #endif
  F40_dum.allocate("F40_dum");
  #ifndef NO_AD_INITIALIZE
  F40_dum.initialize();
  #endif
  F20_out.allocate("F20_out");
  #ifndef NO_AD_INITIALIZE
  F20_out.initialize();
  #endif
  F30_out.allocate("F30_out");
  #ifndef NO_AD_INITIALIZE
  F30_out.initialize();
  #endif
  F40_out.allocate("F40_out");
  #ifndef NO_AD_INITIALIZE
  F40_out.initialize();
  #endif
  SSB_F30_out.allocate("SSB_F30_out");
  #ifndef NO_AD_INITIALIZE
  SSB_F30_out.initialize();
  #endif
  SSB_F30_knum_out.allocate("SSB_F30_knum_out");
  #ifndef NO_AD_INITIALIZE
  SSB_F30_knum_out.initialize();
  #endif
  B_F30_out.allocate("B_F30_out");
  #ifndef NO_AD_INITIALIZE
  B_F30_out.initialize();
  #endif
  R_F30_out.allocate("R_F30_out");
  #ifndef NO_AD_INITIALIZE
  R_F30_out.initialize();
  #endif
  L_F30_knum_out.allocate("L_F30_knum_out");
  #ifndef NO_AD_INITIALIZE
  L_F30_knum_out.initialize();
  #endif
  L_F30_klb_out.allocate("L_F30_klb_out");
  #ifndef NO_AD_INITIALIZE
  L_F30_klb_out.initialize();
  #endif
  SSB_F40_out.allocate("SSB_F40_out");
  #ifndef NO_AD_INITIALIZE
  SSB_F40_out.initialize();
  #endif
  SSB_F40_knum_out.allocate("SSB_F40_knum_out");
  #ifndef NO_AD_INITIALIZE
  SSB_F40_knum_out.initialize();
  #endif
  B_F40_out.allocate("B_F40_out");
  #ifndef NO_AD_INITIALIZE
  B_F40_out.initialize();
  #endif
  R_F40_out.allocate("R_F40_out");
  #ifndef NO_AD_INITIALIZE
  R_F40_out.initialize();
  #endif
  L_F40_knum_out.allocate("L_F40_knum_out");
  #ifndef NO_AD_INITIALIZE
  L_F40_knum_out.initialize();
  #endif
  L_F40_klb_out.allocate("L_F40_klb_out");
  #ifndef NO_AD_INITIALIZE
  L_F40_klb_out.initialize();
  #endif
  rec_mean.allocate("rec_mean");
  #ifndef NO_AD_INITIALIZE
  rec_mean.initialize();
  #endif
  N_age_msy.allocate(1,nages,"N_age_msy");
  #ifndef NO_AD_INITIALIZE
    N_age_msy.initialize();
  #endif
  N_age_msy_spawn.allocate(1,nages,"N_age_msy_spawn");
  #ifndef NO_AD_INITIALIZE
    N_age_msy_spawn.initialize();
  #endif
  L_age_msy.allocate(1,nages,"L_age_msy");
  #ifndef NO_AD_INITIALIZE
    L_age_msy.initialize();
  #endif
  Z_age_msy.allocate(1,nages,"Z_age_msy");
  #ifndef NO_AD_INITIALIZE
    Z_age_msy.initialize();
  #endif
  F_L_age_msy.allocate(1,nages,"F_L_age_msy");
  #ifndef NO_AD_INITIALIZE
    F_L_age_msy.initialize();
  #endif
  F_msy.allocate(1,n_iter_msy,"F_msy");
  #ifndef NO_AD_INITIALIZE
    F_msy.initialize();
  #endif
  spr_msy.allocate(1,n_iter_msy,"spr_msy");
  #ifndef NO_AD_INITIALIZE
    spr_msy.initialize();
  #endif
  R_eq.allocate(1,n_iter_msy,"R_eq");
  #ifndef NO_AD_INITIALIZE
    R_eq.initialize();
  #endif
  L_eq_klb.allocate(1,n_iter_msy,"L_eq_klb");
  #ifndef NO_AD_INITIALIZE
    L_eq_klb.initialize();
  #endif
  L_eq_knum.allocate(1,n_iter_msy,"L_eq_knum");
  #ifndef NO_AD_INITIALIZE
    L_eq_knum.initialize();
  #endif
  SSB_eq.allocate(1,n_iter_msy,"SSB_eq");
  #ifndef NO_AD_INITIALIZE
    SSB_eq.initialize();
  #endif
  SSB_eq_knum.allocate(1,n_iter_msy,"SSB_eq_knum");
  #ifndef NO_AD_INITIALIZE
    SSB_eq_knum.initialize();
  #endif
  B_eq.allocate(1,n_iter_msy,"B_eq");
  #ifndef NO_AD_INITIALIZE
    B_eq.initialize();
  #endif
  FdF_msy.allocate(styr,endyr,"FdF_msy");
  #ifndef NO_AD_INITIALIZE
    FdF_msy.initialize();
  #endif
  FdF30.allocate(styr,endyr,"FdF30");
  #ifndef NO_AD_INITIALIZE
    FdF30.initialize();
  #endif
  FdF40.allocate(styr,endyr,"FdF40");
  #ifndef NO_AD_INITIALIZE
    FdF40.initialize();
  #endif
  SdSSB_msy.allocate(styr,endyr,"SdSSB_msy");
  #ifndef NO_AD_INITIALIZE
    SdSSB_msy.initialize();
  #endif
  SdSSB_msy_end.allocate("SdSSB_msy_end");
  #ifndef NO_AD_INITIALIZE
  SdSSB_msy_end.initialize();
  #endif
  FdF_msy_end.allocate("FdF_msy_end");
  #ifndef NO_AD_INITIALIZE
  FdF_msy_end.initialize();
  #endif
  FdF_msy_end_mean.allocate("FdF_msy_end_mean");
  #ifndef NO_AD_INITIALIZE
  FdF_msy_end_mean.initialize();
  #endif
  SdSSB_F30.allocate(styr,endyr,"SdSSB_F30");
  #ifndef NO_AD_INITIALIZE
    SdSSB_F30.initialize();
  #endif
  Sdmsst_F30.allocate(styr,endyr,"Sdmsst_F30");
  #ifndef NO_AD_INITIALIZE
    Sdmsst_F30.initialize();
  #endif
  SdSSB_F30_end.allocate("SdSSB_F30_end");
  #ifndef NO_AD_INITIALIZE
  SdSSB_F30_end.initialize();
  #endif
  Sdmsst_F30_end.allocate("Sdmsst_F30_end");
  #ifndef NO_AD_INITIALIZE
  Sdmsst_F30_end.initialize();
  #endif
  FdF30_end_mean.allocate("FdF30_end_mean");
  #ifndef NO_AD_INITIALIZE
  FdF30_end_mean.initialize();
  #endif
  L_age_F30.allocate(1,nages,"L_age_F30");
  #ifndef NO_AD_INITIALIZE
    L_age_F30.initialize();
  #endif
  SdSSB_F40.allocate(styr,endyr,"SdSSB_F40");
  #ifndef NO_AD_INITIALIZE
    SdSSB_F40.initialize();
  #endif
  Sdmsst_F40.allocate(styr,endyr,"Sdmsst_F40");
  #ifndef NO_AD_INITIALIZE
    Sdmsst_F40.initialize();
  #endif
  SdSSB_F40_end.allocate("SdSSB_F40_end");
  #ifndef NO_AD_INITIALIZE
  SdSSB_F40_end.initialize();
  #endif
  Sdmsst_F40_end.allocate("Sdmsst_F40_end");
  #ifndef NO_AD_INITIALIZE
  Sdmsst_F40_end.initialize();
  #endif
  FdF40_end_mean.allocate("FdF40_end_mean");
  #ifndef NO_AD_INITIALIZE
  FdF40_end_mean.initialize();
  #endif
  Fend_mean_temp.allocate("Fend_mean_temp");
  #ifndef NO_AD_INITIALIZE
  Fend_mean_temp.initialize();
  #endif
  Fend_mean.allocate("Fend_mean");
  #ifndef NO_AD_INITIALIZE
  Fend_mean.initialize();
  #endif
  L_age_F40.allocate(1,nages,"L_age_F40");
  #ifndef NO_AD_INITIALIZE
    L_age_F40.initialize();
  #endif
  wgt_wgted_L_klb.allocate(1,nages,"wgt_wgted_L_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_wgted_L_klb.initialize();
  #endif
  wgt_wgted_L_denom.allocate("wgt_wgted_L_denom");
  #ifndef NO_AD_INITIALIZE
  wgt_wgted_L_denom.initialize();
  #endif
  iter_inc_msy.allocate("iter_inc_msy");
  #ifndef NO_AD_INITIALIZE
  iter_inc_msy.initialize();
  #endif
  M.allocate(1,nages,"M");
  #ifndef NO_AD_INITIALIZE
    M.initialize();
  #endif
  M_constant.allocate(M_constant_LO,M_constant_HI,M_constant_PH,"M_constant");
  M_constant_out.allocate(1,8,"M_constant_out");
  #ifndef NO_AD_INITIALIZE
    M_constant_out.initialize();
  #endif
  smsy2msstM.allocate("smsy2msstM");
  #ifndef NO_AD_INITIALIZE
  smsy2msstM.initialize();
  #endif
  smsy2msst75.allocate("smsy2msst75");
  #ifndef NO_AD_INITIALIZE
  smsy2msst75.initialize();
  #endif
  F.allocate(styr,endyr,1,nages,"F");
  #ifndef NO_AD_INITIALIZE
    F.initialize();
  #endif
  Fsum.allocate(styr,endyr,"Fsum");
  #ifndef NO_AD_INITIALIZE
    Fsum.initialize();
  #endif
  Fapex.allocate(styr,endyr,"Fapex");
  #ifndef NO_AD_INITIALIZE
    Fapex.initialize();
  #endif
  Z.allocate(styr,endyr,1,nages,"Z");
  #ifndef NO_AD_INITIALIZE
    Z.initialize();
  #endif
  log_avg_F_L_cGN.allocate(log_avg_F_cGN_LO,log_avg_F_cGN_HI,log_avg_F_cGN_PH,"log_avg_F_L_cGN");
  log_avg_F_cGN_out.allocate(1,8,"log_avg_F_cGN_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_cGN_out.initialize();
  #endif
  log_dev_F_L_cGN.allocate(styr_L_cGN,endyr_L_cGN,log_F_dev_cGN_LO,log_F_dev_cGN_HI,log_F_dev_cGN_PH,"log_dev_F_L_cGN");
  log_F_dev_cGN_out.allocate(styr_L_cGN,endyr_L_cGN,"log_F_dev_cGN_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_cGN_out.initialize();
  #endif
  F_cGN.allocate(styr,endyr,1,nages,"F_cGN");
  #ifndef NO_AD_INITIALIZE
    F_cGN.initialize();
  #endif
  F_cGN_out.allocate(styr,endyr,"F_cGN_out");
  #ifndef NO_AD_INITIALIZE
    F_cGN_out.initialize();
  #endif
  log_F_dev_init_cGN.allocate("log_F_dev_init_cGN");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_init_cGN.initialize();
  #endif
  log_F_dev_end_cGN.allocate("log_F_dev_end_cGN");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_cGN.initialize();
  #endif
  log_avg_F_L_rGN.allocate(log_avg_F_rGN_LO,log_avg_F_rGN_HI,log_avg_F_rGN_PH,"log_avg_F_L_rGN");
  log_avg_F_rGN_out.allocate(1,8,"log_avg_F_rGN_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_rGN_out.initialize();
  #endif
  log_dev_F_L_rGN.allocate(styr_L_rGN,endyr_L_rGN,log_F_dev_rGN_LO,log_F_dev_rGN_HI,log_F_dev_rGN_PH,"log_dev_F_L_rGN");
  log_F_dev_rGN_out.allocate(styr_L_rGN,endyr_L_rGN,"log_F_dev_rGN_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_rGN_out.initialize();
  #endif
  F_rGN.allocate(styr,endyr,1,nages,"F_rGN");
  #ifndef NO_AD_INITIALIZE
    F_rGN.initialize();
  #endif
  F_rGN_out.allocate(styr,endyr,"F_rGN_out");
  #ifndef NO_AD_INITIALIZE
    F_rGN_out.initialize();
  #endif
  log_F_dev_init_rGN.allocate("log_F_dev_init_rGN");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_init_rGN.initialize();
  #endif
  log_F_dev_end_rGN.allocate("log_F_dev_end_rGN");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_rGN.initialize();
  #endif
  F_init.allocate(F_init_LO,F_init_HI,F_init_PH,"F_init");
  F_init_out.allocate(1,8,"F_init_out");
  #ifndef NO_AD_INITIALIZE
    F_init_out.initialize();
  #endif
  F_init_denom.allocate("F_init_denom");
  #ifndef NO_AD_INITIALIZE
  F_init_denom.initialize();
  #endif
  N_age_spr.allocate(1,nages,"N_age_spr");
  #ifndef NO_AD_INITIALIZE
    N_age_spr.initialize();
  #endif
  N_age_spr_spawn.allocate(1,nages,"N_age_spr_spawn");
  #ifndef NO_AD_INITIALIZE
    N_age_spr_spawn.initialize();
  #endif
  L_age_spr.allocate(1,nages,"L_age_spr");
  #ifndef NO_AD_INITIALIZE
    L_age_spr.initialize();
  #endif
  Z_age_spr.allocate(1,nages,"Z_age_spr");
  #ifndef NO_AD_INITIALIZE
    Z_age_spr.initialize();
  #endif
  spr_static.allocate(styr,endyr,"spr_static");
  #ifndef NO_AD_INITIALIZE
    spr_static.initialize();
  #endif
  F_L_age_spr.allocate(1,nages,"F_L_age_spr");
  #ifndef NO_AD_INITIALIZE
    F_L_age_spr.initialize();
  #endif
  F_spr.allocate(1,n_iter_spr,"F_spr");
  #ifndef NO_AD_INITIALIZE
    F_spr.initialize();
  #endif
  spr_spr.allocate(1,n_iter_spr,"spr_spr");
  #ifndef NO_AD_INITIALIZE
    spr_spr.initialize();
  #endif
  spr_ratio.allocate(1,n_iter_spr,"spr_ratio");
  #ifndef NO_AD_INITIALIZE
    spr_ratio.initialize();
  #endif
  L_spr.allocate(1,n_iter_spr,"L_spr");
  #ifndef NO_AD_INITIALIZE
    L_spr.initialize();
  #endif
  N_spr_F0.allocate(1,nages,"N_spr_F0");
  #ifndef NO_AD_INITIALIZE
    N_spr_F0.initialize();
  #endif
  N_bpr_F0.allocate(1,nages,"N_bpr_F0");
  #ifndef NO_AD_INITIALIZE
    N_bpr_F0.initialize();
  #endif
  N_spr_initial.allocate(1,nages,"N_spr_initial");
  #ifndef NO_AD_INITIALIZE
    N_spr_initial.initialize();
  #endif
  N_initial_eq.allocate(1,nages,"N_initial_eq");
  #ifndef NO_AD_INITIALIZE
    N_initial_eq.initialize();
  #endif
  F_initial.allocate(1,nages,"F_initial");
  #ifndef NO_AD_INITIALIZE
    F_initial.initialize();
  #endif
  Z_initial.allocate(1,nages,"Z_initial");
  #ifndef NO_AD_INITIALIZE
    Z_initial.initialize();
  #endif
  spr_initial.allocate("spr_initial");
  #ifndef NO_AD_INITIALIZE
  spr_initial.initialize();
  #endif
  spr_F0.allocate("spr_F0");
  #ifndef NO_AD_INITIALIZE
  spr_F0.initialize();
  #endif
  bpr_F0.allocate("bpr_F0");
  #ifndef NO_AD_INITIALIZE
  bpr_F0.initialize();
  #endif
  iter_inc_spr.allocate("iter_inc_spr");
  #ifndef NO_AD_INITIALIZE
  iter_inc_spr.initialize();
  #endif
  sdnr_lc_cGN.allocate("sdnr_lc_cGN");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_cGN.initialize();
  #endif
  sdnr_ac_rGN.allocate("sdnr_ac_rGN");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_rGN.initialize();
  #endif
  sdnr_I_rHB.allocate("sdnr_I_rHB");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_rHB.initialize();
  #endif
  w_L.allocate("w_L");
  #ifndef NO_AD_INITIALIZE
  w_L.initialize();
  #endif
  w_cpue_rHB.allocate("w_cpue_rHB");
  #ifndef NO_AD_INITIALIZE
  w_cpue_rHB.initialize();
  #endif
  w_lenc_cGN.allocate("w_lenc_cGN");
  #ifndef NO_AD_INITIALIZE
  w_lenc_cGN.initialize();
  #endif
  w_agec_rGN.allocate("w_agec_rGN");
  #ifndef NO_AD_INITIALIZE
  w_agec_rGN.initialize();
  #endif
  w_Nage_init.allocate("w_Nage_init");
  #ifndef NO_AD_INITIALIZE
  w_Nage_init.initialize();
  #endif
  w_rec.allocate("w_rec");
  #ifndef NO_AD_INITIALIZE
  w_rec.initialize();
  #endif
  w_rec_early.allocate("w_rec_early");
  #ifndef NO_AD_INITIALIZE
  w_rec_early.initialize();
  #endif
  w_rec_end.allocate("w_rec_end");
  #ifndef NO_AD_INITIALIZE
  w_rec_end.initialize();
  #endif
  w_fullF.allocate("w_fullF");
  #ifndef NO_AD_INITIALIZE
  w_fullF.initialize();
  #endif
  w_Ftune.allocate("w_Ftune");
  #ifndef NO_AD_INITIALIZE
  w_Ftune.initialize();
  #endif
  f_cGN_L.allocate("f_cGN_L");
  #ifndef NO_AD_INITIALIZE
  f_cGN_L.initialize();
  #endif
  f_rGN_L.allocate("f_rGN_L");
  #ifndef NO_AD_INITIALIZE
  f_rGN_L.initialize();
  #endif
  f_rHB_cpue.allocate("f_rHB_cpue");
  #ifndef NO_AD_INITIALIZE
  f_rHB_cpue.initialize();
  #endif
  f_rHB_RWq_cpue.allocate("f_rHB_RWq_cpue");
  #ifndef NO_AD_INITIALIZE
  f_rHB_RWq_cpue.initialize();
  #endif
  f_cGN_lenc.allocate("f_cGN_lenc");
  #ifndef NO_AD_INITIALIZE
  f_cGN_lenc.initialize();
  #endif
  f_rGN_lenc.allocate("f_rGN_lenc");
  #ifndef NO_AD_INITIALIZE
  f_rGN_lenc.initialize();
  #endif
  f_rGN_agec.allocate("f_rGN_agec");
  #ifndef NO_AD_INITIALIZE
  f_rGN_agec.initialize();
  #endif
  f_Nage_init.allocate("f_Nage_init");
  #ifndef NO_AD_INITIALIZE
  f_Nage_init.initialize();
  #endif
  f_rec_dev.allocate("f_rec_dev");
  #ifndef NO_AD_INITIALIZE
  f_rec_dev.initialize();
  #endif
  f_rec_dev_early.allocate("f_rec_dev_early");
  #ifndef NO_AD_INITIALIZE
  f_rec_dev_early.initialize();
  #endif
  f_rec_dev_end.allocate("f_rec_dev_end");
  #ifndef NO_AD_INITIALIZE
  f_rec_dev_end.initialize();
  #endif
  f_fullF_constraint.allocate("f_fullF_constraint");
  #ifndef NO_AD_INITIALIZE
  f_fullF_constraint.initialize();
  #endif
  f_Ftune.allocate("f_Ftune");
  #ifndef NO_AD_INITIALIZE
  f_Ftune.initialize();
  #endif
  f_priors.allocate("f_priors");
  #ifndef NO_AD_INITIALIZE
  f_priors.initialize();
  #endif
  fval.allocate("fval");
  prior_function_value.allocate("prior_function_value");
  likelihood_function_value.allocate("likelihood_function_value");
  fval_data.allocate("fval_data");
  #ifndef NO_AD_INITIALIZE
  fval_data.initialize();
  #endif
  grad_max.allocate("grad_max");
  #ifndef NO_AD_INITIALIZE
  grad_max.initialize();
  #endif
  denom.allocate("denom");
  #ifndef NO_AD_INITIALIZE
  denom.initialize();
  #endif
  numer.allocate("numer");
  #ifndef NO_AD_INITIALIZE
  numer.initialize();
  #endif
  F_reg_proj.allocate("F_reg_proj");
  #ifndef NO_AD_INITIALIZE
  F_reg_proj.initialize();
  #endif
  F_proj.allocate(styr_proj,endyr_proj,"F_proj");
  #ifndef NO_AD_INITIALIZE
    F_proj.initialize();
  #endif
  L_knum_proj.allocate(styr_proj,endyr_proj,"L_knum_proj");
  #ifndef NO_AD_INITIALIZE
    L_knum_proj.initialize();
  #endif
  L_klb_proj.allocate(styr_proj,endyr_proj,"L_klb_proj");
  #ifndef NO_AD_INITIALIZE
    L_klb_proj.initialize();
  #endif
  B_proj.allocate(styr_proj,endyr_proj,"B_proj");
  #ifndef NO_AD_INITIALIZE
    B_proj.initialize();
  #endif
  SSB_proj.allocate(styr_proj,endyr_proj,"SSB_proj");
  #ifndef NO_AD_INITIALIZE
    SSB_proj.initialize();
  #endif
  R_proj.allocate(styr_proj,endyr_proj,"R_proj");
  #ifndef NO_AD_INITIALIZE
    R_proj.initialize();
  #endif
  FL_age_proj.allocate(1,nages,"FL_age_proj");
  #ifndef NO_AD_INITIALIZE
    FL_age_proj.initialize();
  #endif
  N_proj.allocate(styr_proj,endyr_proj,1,nages,"N_proj");
  #ifndef NO_AD_INITIALIZE
    N_proj.initialize();
  #endif
  N_spawn_proj.allocate(styr_proj,endyr_proj,1,nages,"N_spawn_proj");
  #ifndef NO_AD_INITIALIZE
    N_spawn_proj.initialize();
  #endif
  Z_proj.allocate(styr_proj,endyr_proj,1,nages,"Z_proj");
  #ifndef NO_AD_INITIALIZE
    Z_proj.initialize();
  #endif
  L_age_proj.allocate(styr_proj,endyr_proj,1,nages,"L_age_proj");
  #ifndef NO_AD_INITIALIZE
    L_age_proj.initialize();
  #endif
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{1000, 2000,3000, 5000, 10000;//, 10000, 10000;}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-2, 1e-2,1e-3, 1e-3, 1e-4;//, 1e-4, 1e-4;}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
 
 //Population		
  Linf=set_Linf(1);
  K=set_K(1);
  t0=set_t0(1);
  len_cv_val=set_len_cv(1);
  
  //All fisheries
  Linf_L=set_Linf_L(1);  
  K_L=set_K_L(1);
  t0_L=set_t0_L(1);
  len_cv_val_L=set_len_cv_L(1);
  //Females
  Linf_F=set_Linf_F(1);  
  K_F=set_K_F(1);
  t0_F=set_t0_F(1);
  len_cv_val_F=set_len_cv_F(1);
  
  M=set_M; 
  M_constant=set_M_constant(1);
  smsy2msstM=1.0-M_constant;
  smsy2msst75=0.75;  
  
  log_R0=set_log_R0(1);
  steep=set_steep(1);
  R_autocorr=set_R_autocorr(1);
  rec_sigma=set_rec_sigma(1);
  
  log_dm_lenc_cGN=set_log_dm_lenc_cGN(1);
  //log_dm_cLL_lc=set_log_dm_cLL_lc(1);
  // log_dm_rHB_lc=set_log_dm_rHB_lc(1);
  //log_dm_rGN_lc=set_log_dm_rGN_lc(1);
  log_dm_agec_rGN=set_log_dm_agec_rGN(1);
  
 // log_q_cGN=set_log_q_cGN(1);
  //log_q_cLL=set_log_q_cLL(1);  
  log_q_cpue_rHB=set_log_q_cpue_rHB(1);
  
  q_rate=set_q_rate;
  //q_rate_fcn_cGN=1.0;
  //q_rate_fcn_cLL=1.0;    
  q_rate_fcn_rHB=1.0;   
  q_DD_beta=set_q_DD_beta;
  q_DD_fcn=1.0;
  //q_RW_log_dev_cGN.initialize();
 // q_RW_log_dev_cLL.initialize();   
  q_RW_log_dev_rHB.initialize(); 
  
   if (set_q_rate_phase<0 & q_rate!=0.0)
  {
      
    for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++)
      {   if (iyear>styr_cpue_rHB & iyear <=2003) 
          {//q_rate_fcn_rHB(iyear)=(1.0+q_rate)*q_rate_fcn_rHB(iyear-1); //compound
             q_rate_fcn_rHB(iyear)=(1.0+(iyear-styr_cpue_rHB)*q_rate)*q_rate_fcn_rHB(styr_cpue_rHB);  //linear
          }
          if (iyear>2003) {q_rate_fcn_rHB(iyear)=q_rate_fcn_rHB(iyear-1);} 
      }   
  } //end q_rate conditional      
  w_L=set_w_L;
  
 // w_I_cGN=set_w_I_cGN;
 // w_I_cLL=set_w_I_cLL;
  w_cpue_rHB=set_w_cpue_rHB;
  
  w_lenc_cGN=set_w_lenc_cGN;
  //w_lc_cLL=set_w_lc_cLL;  
  //w_lc_rHB=set_w_lc_rHB;
  //w_lc_rGN=set_w_lc_rGN;  
  w_agec_rGN=set_w_agec_rGN; 
  
  w_Nage_init=set_w_Nage_init;
  w_rec=set_w_rec;
  w_rec_early=set_w_rec_early;
  w_rec_end=set_w_rec_end;
  w_fullF=set_w_fullF;
  w_Ftune=set_w_Ftune;
  
  F_init=set_F_init(1);  
  log_avg_F_L_cGN=set_log_avg_F_L_cGN(1);
  //log_avg_F_cLL=set_log_avg_F_cLL(1);  
  // log_avg_F_rHB=set_log_avg_F_rHB(1); 
  log_avg_F_L_rGN=set_log_avg_F_L_rGN(1); 
    
  log_dev_F_L_cGN=set_log_dev_vals_F_L_cGN;
  //log_F_dev_cLL=set_log_F_dev_cLL_vals;
  // log_F_dev_rHB=set_log_F_dev_rHB_vals;
  log_dev_F_L_rGN=set_log_dev_vals_F_L_rGN;
 
  selpar_A50_cGN1=set_selpar_A50_cGN1(1);
  selpar_slope_cGN1=set_selpar_slope_cGN1(1);
  //selpar_A50_cGN2=set_selpar_A50_cGN2(1);
  //selpar_slope_cGN2=set_selpar_slope_cGN2(1);
  // selpar_A502_cGN2=set_selpar_A502_cGN2(1);
  // selpar_slope2_cGN2=set_selpar_slope2_cGN2(1);  
  //selpar_A50_cGN3=set_selpar_A50_cGN3(1);
  //selpar_slope_cGN3=set_selpar_slope_cGN3(1);  
  
  selpar_A50_rGN1=set_selpar_A50_rGN1(1);
  selpar_slope_rGN1=set_selpar_slope_rGN1(1);
  selpar_A50_rGN2=set_selpar_A50_rGN2(1);
  selpar_slope_rGN2=set_selpar_slope_rGN2(1);
  // selpar_A502_rGN2=set_selpar_A502_rGN2(1);
  // selpar_slope2_rGN2=set_selpar_slope2_rGN2(1);  
  //selpar_A50_rGN3=set_selpar_A50_rGN3(1);
  //selpar_slope_rGN3=set_selpar_slope_rGN3(1);  
 sqrt2pi=sqrt(2.*3.14159265);
 g2mt=0.000001;         //conversion of grams to metric tons
 g2kg=0.001;            //conversion of grams to kg 
 mt2klb=2.20462;        //conversion of metric tons to 1000 lb 
 mt2lb=mt2klb*1000.0;   //conversion of metric tons to lb
 g2klb=g2mt*mt2klb;     //conversion of grams to 1000 lb 
 dzero=0.00001;         
 huge_number=1.0e+10;   
 
 SSB_msy_out=0.0;
 iter_inc_msy=max_F_spr_msy/(n_iter_msy-1);
 iter_inc_spr=max_F_spr_msy/(n_iter_spr-1); 
 maturity_f=obs_maturity_f;
 
 prop_f=obs_prop_f;
      nsamp_cGN_lenc_allyr=missing;
	  //nsamp_cLL_lenc_allyr=missing;   
	  //nsamp_rGN_lenc_allyr=missing;	  
	  nsamp_rGN_agec_allyr=missing;  
	  
      nfish_cGN_lenc_allyr=missing;
	  //nfish_cLL_lenc_allyr=missing;  
	  //nfish_rGN_lenc_allyr=missing;
	  nfish_rGN_agec_allyr=missing;
   
      for (iyear=1; iyear<=nyr_lenc_cGN; iyear++)
         {if (nsamp_lenc_cGN(iyear)>=minSS_lenc_cGN)
           {nsamp_cGN_lenc_allyr(yrs_lenc_cGN(iyear))=nsamp_lenc_cGN(iyear);
            nfish_cGN_lenc_allyr(yrs_lenc_cGN(iyear))=nfish_lenc_cGN(iyear);}}
      //for (iyear=1; iyear<=nyr_cLL_lenc; iyear++)
      //   {if (nsamp_cLL_lenc(iyear)>=minSS_cLL_lenc)
      //     {nsamp_cLL_lenc_allyr(yrs_cLL_lenc(iyear))=nsamp_cLL_lenc(iyear);
      //      nfish_cLL_lenc_allyr(yrs_cLL_lenc(iyear))=nfish_cLL_lenc(iyear);}}			
	  //for (iyear=1; iyear<=nyr_rGN_lenc; iyear++)                           
      //   {if (nsamp_rGN_lenc(iyear)>=minSS_rGN_lenc)
      //      {nsamp_rGN_lenc_allyr(yrs_rGN_lenc(iyear))=nsamp_rGN_lenc(iyear);
      //       nfish_rGN_lenc_allyr(yrs_rGN_lenc(iyear))=nfish_rGN_lenc(iyear);}}
	  for (iyear=1; iyear<=nyr_agec_rGN; iyear++)  
         {if (nsamp_agec_rGN(iyear)>=minSS_agec_rGN)
           {nsamp_rGN_agec_allyr(yrs_agec_rGN(iyear))=nsamp_agec_rGN(iyear);
             nfish_rGN_agec_allyr(yrs_agec_rGN(iyear))=nfish_agec_rGN(iyear);}}  
	  
  F_msy(1)=0.0;  
  for (ff=2;ff<=n_iter_msy;ff++) {F_msy(ff)=F_msy(ff-1)+iter_inc_msy;}
  F_spr(1)=0.0;  
  for (ff=2;ff<=n_iter_spr;ff++) {F_spr(ff)=F_spr(ff-1)+iter_inc_spr;}
  F_cGN.initialize(); L_cGN_num.initialize();
  //F_cLL.initialize(); L_cLL_num.initialize();
  // F_rHB.initialize(); L_rHB_num.initialize();
  F_rGN.initialize(); L_rGN_num.initialize();
  F_cGN_out.initialize();
  //F_cLL_out.initialize();
  // F_rHB_out.initialize();
  F_rGN_out.initialize();
  sel_cGN.initialize();
  sel_cGN_vec.initialize();
  sel_rHB.initialize();
  sel_rGN.initialize();
  
  //sel_cGN_block1.initialize();  
  //sel_cGN_block2.initialize();  
  //sel_cLL_block1.initialize();    
  //sel_cLL_block2.initialize(); 
  sel_rHB_block1.initialize();
  sel_rHB_block2.initialize();
  sel_rGN_block1.initialize();
  sel_rGN_block2.initialize();
  
  
  log_rec_dev_output.initialize();  
  log_dev_rec=set_log_dev_vals_rec;
  log_Nage_dev_output.initialize();
  log_dev_Nage=set_log_dev_vals_Nage;
 
 
}

void model_parameters::userfunction(void)
{
  fval =0.0;
 //cout<<"start"<<endl;
 //get_M_at_age(); //Needed only if M is estimated
  get_length_weight_at_age(); 
  //cout << "got length, weight, fecundity transitions" <<endl;
  get_reprod();
  //cout << "got reprod" << endl;
  get_length_at_age_dist(); 
  //cout<< "got predicted length at age distribution"<<endl;
  get_weight_at_age_landings();
  //cout<< "got weight at age of landings"<<endl; 
  get_spr_F0();
  //cout << "got F0 spr" << endl;
  get_selectivity(); 
  //cout << "got selectivity" << endl;
  get_mortality(); 
  // cout << "got mortalities" << endl;
  get_bias_corr(); 
  //cout<< "got recruitment bias correction" << endl;
  get_numbers_at_age(); 
  //cout << "got numbers at age" << endl;
  get_landings_numbers();
  //cout << "got landings in numbers" << endl;
  get_landings_wgt();
  //cout << "got landings in wgt" << endl;
  // get_dead_discards(); 
  //cout << "got dead discards in num and wgt" << endl;
  get_catchability_fcns(); 
  //cout << "got catchability_fcns" << endl;
  get_indices();
  //cout << "got indices" << endl;
  get_length_comps();
  // cout<< "got length comps"<< endl;
   get_age_comps();
  //cout<< "got age comps"<< endl;
   evaluate_objective_function();
   //cout << "objective function calculations complete" << endl;
}

void model_parameters::get_length_weight_at_age(void)
{
	//population total length in mm
    //compute mean length (mm TL) and weight (whole) at age
    meanlen_TL=Linf*(1.0-mfexp(-K*(agebins-t0+0.5))); //Actually fork length     
    wgt_kg=wgtpar_a*pow(meanlen_TL,wgtpar_b);             //whole wgt in kg 
    wgt_g=wgt_kg/g2kg;                                    //convert wgt in kg to weight in g    
    wgt_mt=wgt_g*g2mt;                                    //convert weight in g to weight in mt
    wgt_klb=mt2klb*wgt_mt;                                //1000 lb of whole wgt
    wgt_lb=mt2lb*wgt_mt;                                  //lb of whole wgt
	//All fisheries
	meanlen_TL_L=Linf_L*(1.0-mfexp(-K_L*(agebins-t0_L+0.5)));    //Landings total length in mm
    wgt_kg_L=wgtpar_a*pow(meanlen_TL_L,wgtpar_b);             //whole wgt in kg 
    wgt_g_L=wgt_kg_L/g2kg;                                    //convert wgt in kg to weight in g    
    wgt_mt_L=wgt_g_L*g2mt;                                    //convert weight in g to weight in mt
    wgt_klb_L=mt2klb*wgt_mt_L;                               //1000 lb of whole wgt
    wgt_lb_L=mt2lb*wgt_mt_L;                                 //1000 lb of whole wgt
	//Females
	meanlen_TL_F=Linf_F*(1.0-mfexp(-K_F*(agebins-t0_F+0.5)));    //Landings total length in mm
    wgt_kg_F=wgtpar_a*pow(meanlen_TL_F,wgtpar_b);             //whole wgt in kg 
    wgt_g_F=wgt_kg_F/g2kg;                                    //convert wgt in kg to weight in g    
    wgt_mt_F=wgt_g_F*g2mt;                                    //convert weight in g to weight in mt
    wgt_klb_F=mt2klb*wgt_mt_F;                               //1000 lb of whole wgt
    wgt_lb_F=mt2lb*wgt_mt_F;                                 //1000 lb of whole wgt
	//batchfec = mfexp(batchfecpar_a + batchfecpar_b*meanlen_TL);  // batch fecundity at length [should be batchfec = exp(a+bL) based on Harris 2004]
	//fec = batchfec*nbatch/fecpar_scale;                          // annual fecundity at length scaled to fecpar_scale units
}

void model_parameters::get_reprod(void)
{
    //reprod=elem_prod(prop_f,elem_prod(maturity_f,fec));
	reprod=elem_prod(elem_prod(prop_f,maturity_f),wgt_mt_F);
	reprodknum=elem_prod(prop_f,maturity_f)/1000.0;
	//+elem_prod(prop_m,maturity_m)),wg_mt);
}

void model_parameters::get_length_at_age_dist(void)
{
  //compute matrix of length at age, based on the normal distribution
    //population
	for (iage=1;iage<=nages;iage++)
   {len_cv(iage)=len_cv_val;
    len_sd(iage)=meanlen_TL(iage)*len_cv(iage);
	zscore_lzero=(0.0-meanlen_TL(iage))/len_sd(iage); 
	cprob_lzero=cumd_norm(zscore_lzero);
    //All fishery dependent
	//len_cv_L(iage)=mfexp(log_len_cv_L+log_len_cv_dev_L(iage));
    len_cv_L(iage)=len_cv_val_L;
    len_sd_L(iage)=meanlen_TL_L(iage)*len_cv_L(iage);
	zscore_lzero_L=(0.0-meanlen_TL_L(iage))/len_sd_L(iage); 
	cprob_lzero_L=cumd_norm(zscore_lzero_L);    
	//Females
	//len_cv_L(iage)=mfexp(log_len_cv_L+log_len_cv_dev_L(iage));
    len_cv_F(iage)=len_cv_val_F;
    len_sd_F(iage)=meanlen_TL_F(iage)*len_cv_F(iage);
	zscore_lzero_F=(0.0-meanlen_TL_F(iage))/len_sd_F(iage); 
	cprob_lzero_F=cumd_norm(zscore_lzero_F);
    //first length bin
	//population
    zscore_len=((lenbins(1)+0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage);
    cprob_lenvec(1)=cumd_norm(zscore_len);          //includes any probability mass below zero
    lenprob(iage,1)=cprob_lenvec(1)-cprob_lzero;    //removes any probability mass below zero
	//All fishery dependent
    zscore_len_L=((lenbins(1)+0.5*lenbins_width)-meanlen_TL_L(iage)) / len_sd_L(iage);
    cprob_lenvec_L(1)=cumd_norm(zscore_len_L);          //includes any probability mass below zero
    lenprob_L(iage,1)=cprob_lenvec_L(1)-cprob_lzero_L;    //removes any probability mass below zero
	//Females
    zscore_len_F=((lenbins(1)+0.5*lenbins_width)-meanlen_TL_F(iage)) / len_sd_F(iage);
    cprob_lenvec_F(1)=cumd_norm(zscore_len_F);          //includes any probability mass below zero
    lenprob_F(iage,1)=cprob_lenvec_F(1)-cprob_lzero_F;    //removes any probability mass below zero
    //most other length bins  
    //population
    for (ilen=2;ilen<nlenbins;ilen++)
      {
        zscore_len=((lenbins(ilen)+0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage); 
		cprob_lenvec(ilen)=cumd_norm(zscore_len);
        lenprob(iage,ilen)=cprob_lenvec(ilen)-cprob_lenvec(ilen-1);
      }
	//All fishery dependent	
    for (ilen=2;ilen<nlenbins;ilen++)
      {
        zscore_len_L=((lenbins(ilen)+0.5*lenbins_width)-meanlen_TL_L(iage)) / len_sd_L(iage); 
		cprob_lenvec_L(ilen)=cumd_norm(zscore_len_L);
        lenprob_L(iage,ilen)=cprob_lenvec_L(ilen)-cprob_lenvec_L(ilen-1);
      }
	//Females	
    for (ilen=2;ilen<nlenbins;ilen++)
      {
        zscore_len_F=((lenbins(ilen)+0.5*lenbins_width)-meanlen_TL_F(iage)) / len_sd_F(iage); 
		cprob_lenvec_F(ilen)=cumd_norm(zscore_len_F);
        lenprob_F(iage,ilen)=cprob_lenvec_F(ilen)-cprob_lenvec_F(ilen-1);
      }
    //last length bin is a plus group
	//population
    zscore_len=((lenbins(nlenbins)-0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage); 
	lenprob(iage,nlenbins)=1.0-cumd_norm(zscore_len);
      lenprob(iage)=lenprob(iage)/(1.0-cprob_lzero);  //renormalize to account for any prob mass below size=0
	//All fishery dependent
    zscore_len_L=((lenbins(nlenbins)-0.5*lenbins_width)-meanlen_TL_L(iage)) / len_sd_L(iage); 
	lenprob_L(iage,nlenbins)=1.0-cumd_norm(zscore_len_L);
      lenprob_L(iage)=lenprob_L(iage)/(1.0-cprob_lzero_L);  //renormalize to account for any prob mass below size=0
	//Females
    zscore_len_F=((lenbins(nlenbins)-0.5*lenbins_width)-meanlen_TL_F(iage)) / len_sd_F(iage); 
	lenprob_F(iage,nlenbins)=1.0-cumd_norm(zscore_len_F);
      lenprob_F(iage)=lenprob_F(iage)/(1.0-cprob_lzero_F);  //renormalize to account for any prob mass below size=0
   }
  //fleet and survey specific length probs, all assumed here to equal the popn
  lenprob_cGN=lenprob_L;
  //lenprob_cLL=lenprob; 
  lenprob_rHB=lenprob;
  //lenprob_rGN=lenprob;  
}

void model_parameters::get_weight_at_age_landings(void)
{
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    len_cGN_mm(iyear)=meanlen_TL_L;  
    wholewgt_cGN_klb(iyear)=wgt_klb_L; 
    //len_cLL_mm(iyear)=meanlen_TL;  
    //wholewgt_cLL_klb(iyear)=wgt_klb; 	
    len_rHB_mm(iyear)=meanlen_TL_L;
    wholewgt_rHB_klb(iyear)=wgt_klb_L;
    len_rGN_mm(iyear)=meanlen_TL_L;
    wholewgt_rGN_klb(iyear)=wgt_klb_L;
  }  
}

void model_parameters::get_spr_F0(void)
{
  //at mdyr, apply half this yr's mortality, half next yr's
  N_spr_F0(1)=1.0*mfexp(-1.0*M(1)*spawn_time_frac); //at peak spawning time
  N_bpr_F0(1)=1.0;      //at start of year
  for (iage=2; iage<=nages; iage++)
  { N_spr_F0(iage)=N_spr_F0(iage-1)*mfexp(-1.0*(M(iage-1)*(1.0-spawn_time_frac) + M(iage)*spawn_time_frac)); 
    N_bpr_F0(iage)=N_bpr_F0(iage-1)*mfexp(-1.0*(M(iage-1)));    
  }
  N_spr_F0(nages)=N_spr_F0(nages)/(1.0-mfexp(-1.0*M(nages))); //plus group (sum of geometric series)
  N_bpr_F0(nages)=N_bpr_F0(nages)/(1.0-mfexp(-1.0*M(nages)));
  spr_F0=sum(elem_prod(N_spr_F0,reprod)); 
  bpr_F0=sum(elem_prod(N_bpr_F0,wgt_mt));    
}

void model_parameters::get_selectivity(void)
{
  sel_cGN_vec=logistic(agebins, selpar_A50_cGN1, selpar_slope_cGN1);
  sel_rGN_block1=logistic(agebins, selpar_A50_rGN1, selpar_slope_rGN1);
  sel_rGN_block2=logistic(agebins, selpar_A50_rGN2, selpar_slope_rGN2);
  sel_rHB_block1=sel_rGN_block1; // Use rGN selectivity for rHB
  sel_rHB_block2=sel_rGN_block1; // Use rGN selectivity for rHB
  //-------- cGN --------// 
  for (iyear=styr; iyear<=endyr; iyear++)
	{sel_cGN(iyear) = sel_cGN_vec;}
  //---- rGN and rHB ----//   
  //BLOCK 1 for selex  
  for (iyear=styr; iyear<=endyr_selex_phase1_rGN; iyear++)
   {     
    sel_rHB(iyear)=sel_rHB_block1;
    sel_rGN(iyear)=sel_rGN_block1;
   }
  //BLOCK 2 for selex  
  for (iyear=(endyr_selex_phase1_rGN+1); iyear<=endyr; iyear++)//iyear<=endyr_selex_phase2_rGN; iyear++)
   {    
    sel_rHB(iyear)=sel_rHB_block2;
    sel_rGN(iyear)=sel_rGN_block2;
   }
}

void model_parameters::get_mortality(void)
{
  Fsum.initialize();
  Fapex.initialize();
  F.initialize();
  //initialization F is avg from first 3 yrs of observed landings
  log_F_dev_init_cGN=sum(log_dev_F_L_cGN(styr_L_cGN,(styr_L_cGN+2)))/3.0;  
  //log_F_dev_init_cLL=sum(log_F_dev_cLL(styr_cLL_L,(styr_cLL_L+2)))/3.0;          
  log_F_dev_init_rGN=sum(log_dev_F_L_rGN(styr_L_rGN,(styr_L_rGN+2)))/3.0;         
  for (iyear=styr; iyear<=endyr; iyear++) 
  {
    if(iyear>=styr_L_cGN & iyear<=endyr_L_cGN) //spans full time series
		{F_cGN_out(iyear)=mfexp(log_avg_F_L_cGN+log_dev_F_L_cGN(iyear));}     
    F_cGN(iyear)=sel_cGN(iyear)*F_cGN_out(iyear);
    Fsum(iyear)+=F_cGN_out(iyear);
	if(iyear>=styr_L_rGN & iyear<=endyr_L_rGN) //starts in 1981
		{F_rGN_out(iyear)=mfexp(log_avg_F_L_rGN+log_dev_F_L_rGN(iyear));}    
    if (iyear<styr_L_rGN)
		{F_rGN_out(iyear)=mfexp(log_avg_F_L_rGN+log_F_dev_init_rGN);}
	F_rGN(iyear)=sel_rGN(iyear)*F_rGN_out(iyear); 
    Fsum(iyear)+=F_rGN_out(iyear);
    //Total F at age
    F(iyear)=F_cGN(iyear);  //first in additive series (NO +=)
    //F(iyear)+=F_cLL(iyear);
    // F(iyear)+=F_rHB(iyear);
    F(iyear)+=F_rGN(iyear);
    Fapex(iyear)=max(F(iyear));
    Z(iyear)=M+F(iyear);
   }  //end iyear 
}

void model_parameters::get_bias_corr(void)
{
  var_rec_dev=norm2(log_dev_rec(styr_rec_dev,endyr_rec_dev)-
              sum(log_dev_rec(styr_rec_dev,endyr_rec_dev))/nyrs_rec)
              /(nyrs_rec-1.0);                           
  //if (set_BiasCor <= 0.0) {BiasCor=mfexp(var_rec_dev/2.0);}   //bias correction based on empirical residuals
  rec_sigma_sq=square(rec_sigma);
  if (set_BiasCor <= 0.0) {BiasCor=mfexp(rec_sigma_sq/2.0);}   //bias correction based on Rsigma               
  else {BiasCor=set_BiasCor;}
}

void model_parameters::get_numbers_at_age(void)
{
  R0=mfexp(log_R0);
  S0=spr_F0*R0;
  //R_virgin=SR_eq_func(R0, steep, spr_F0, spr_F0, BiasCor, SR_switch);
  R_virgin=BiasCor*R0; //changed to move away from an SR relationship
  B0=bpr_F0*R_virgin;   
  B0_q_DD=R_virgin*sum(elem_prod(N_bpr_F0(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages))); 
  // Commented out code block from Erik's ASPM for red snapper  
  F_init_denom=mfexp(log_avg_F_L_cGN+log_F_dev_init_cGN)+mfexp(log_avg_F_L_rGN+log_F_dev_init_rGN); //+mfexp(log_avg_F_cLL+log_F_dev_init_cLL)
  F_init_cGN_prop= mfexp(log_avg_F_L_cGN+log_F_dev_init_cGN)/F_init_denom;
  //F_init_cLL_prop= mfexp(log_avg_F_cLL+log_F_dev_init_cLL)/F_init_denom;
  F_init_rGN_prop= mfexp(log_avg_F_L_rGN+log_F_dev_init_rGN)/F_init_denom;
  F_initial=sel_cGN(styr)*F_init*F_init_cGN_prop+
            //sel_cLL(styr)*F_init*F_init_cLL_prop+
            sel_rGN(styr)*F_init*F_init_rGN_prop;
	//F_initial=sel_initial*F_init;		
  Z_initial=M+F_initial;  
  N_spr_initial(1)=1.0*mfexp(-1.0*Z_initial(1)*spawn_time_frac); //at peak spawning time;
  for (iage=2; iage<=nages; iage++)
    {
      N_spr_initial(iage)=N_spr_initial(iage-1)*
                   mfexp(-1.0*(Z_initial(iage-1)*(1.0-spawn_time_frac) + Z_initial(iage)*spawn_time_frac)); 
    }
  N_spr_initial(nages)=N_spr_initial(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
  spr_initial=sum(elem_prod(N_spr_initial,reprod));
  //if (styr==styr_rec_dev) {R1=SR_eq_func(R0, steep, spr_F0, spr_initial, 1.0, SR_switch);} //without bias correction (deviation added later)
  //else {R1=SR_eq_func(R0, steep, spr_F0, spr_initial, BiasCor, SR_switch);} //with bias correction
  if (styr==styr_rec_dev) {R1=R0;} //without bias correction (deviation added later)
  else {R1=BiasCor*R0;} //with bias correction
  if(R1<10.0) {R1=10.0;} //Avoid unrealistically low popn sizes during search algorithm
  N_initial_eq(1)=R1;
  for (iage=2; iage<=nages; iage++)
  {
    N_initial_eq(iage)=N_initial_eq(iage-1)*
        mfexp(-1.0*(Z_initial(iage-1)));    
  }
  //plus group calculation
  N_initial_eq(nages)=N_initial_eq(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
  N(styr)(2,nages)=elem_prod(N_initial_eq(2,nages),mfexp(log_dev_Nage));
  if (styr==styr_rec_dev) {N(styr,1)=N_initial_eq(1)*mfexp(log_dev_rec(styr_rec_dev));}
  else {N(styr,1)=N_initial_eq(1);}
  N_mdyr(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*0.5))); //mid year 
  N_spawn(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*spawn_time_frac))); //peak spawning time 
  SSB(styr)=sum(elem_prod(N_spawn(styr),reprod));
  SSB_knum(styr)=sum(elem_prod(N_spawn(styr),reprodknum));
  B_q_DD(styr)=sum(elem_prod(N(styr)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
  for (iyear=styr; iyear<endyr; iyear++)
  {
    if(iyear<(styr_rec_dev-1)||iyear>(endyr_rec_dev-1)) //recruitment follows S-R curve (with bias correction) exactly
    {
		N(iyear+1,1)=BiasCor*R0; //Changed to use ave rec instead of SR relationship
		//N(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch);
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); //plus group
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));
		SSB_knum(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprodknum));
		B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));       
    }
    else   //recruitment follows S-R curve with lognormal deviation
    {
        N(iyear+1,1)=R0*mfexp(log_dev_rec(iyear+1)); //Changed to use ave rec instead of SR relationship
        //N(iyear+1,1)=SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch)*mfexp(log_dev_rec(iyear+1));
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); //plus group
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));
		SSB_knum(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprodknum));
        B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
    }
  }
  //last year (projection) has no recruitment variability
  N(endyr+1,1)=BiasCor*R0; //Changed to use ave rec instead of SR relationship
  //N(endyr+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(endyr),SR_switch);
  N(endyr+1)(2,nages)=++elem_prod(N(endyr)(1,nages-1),(mfexp(-1.*Z(endyr)(1,nages-1))));
  N(endyr+1,nages)+=N(endyr,nages)*mfexp(-1.*Z(endyr,nages)); //plus group
}

void model_parameters::get_landings_numbers(void)
{
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      L_cGN_num(iyear,iage)=N(iyear,iage)*F_cGN(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      //L_cLL_num(iyear,iage)=N(iyear,iage)*F_cLL(iyear,iage)*
        //(1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);		
      L_rGN_num(iyear,iage)=N(iyear,iage)*F_rGN(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);        
    }          
    pred_cGN_L_knum(iyear)=sum(L_cGN_num(iyear))/1000.0;
	//pred_cLL_L_knum(iyear)=sum(L_cLL_num(iyear))/1000.0;
    pred_rGN_L_knum(iyear)=sum(L_rGN_num(iyear))/1000.0;
  }
}

void model_parameters::get_landings_wgt(void)
{
  for (iyear=styr; iyear<=endyr; iyear++)
  {    
    L_cGN_klb(iyear)=elem_prod(L_cGN_num(iyear),wholewgt_cGN_klb(iyear));     //in 1000 lb whole weight
	//L_cLL_klb(iyear)=elem_prod(L_cLL_num(iyear),wholewgt_cLL_klb(iyear));     //in 1000 lb whole weight
    // L_rHB_klb(iyear)=elem_prod(L_rHB_num(iyear),wholewgt_rHB_klb(iyear));     //in 1000 lb whole weight
    L_rGN_klb(iyear)=elem_prod(L_rGN_num(iyear),wholewgt_rGN_klb(iyear));     //in 1000 lb whole weight
    pred_cGN_L_klb(iyear)=sum(L_cGN_klb(iyear));
	//pred_cLL_L_klb(iyear)=sum(L_cLL_klb(iyear));
    // pred_rHB_L_klb(iyear)=sum(L_rHB_klb(iyear));
    pred_rGN_L_klb(iyear)=sum(L_rGN_klb(iyear));    
  }
}

void model_parameters::get_catchability_fcns(void)
{
 //Get rate increase if estimated, otherwise fixed above
  if (set_q_rate_phase>0.0)
  {
      for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++)
      {   if (iyear>styr_cpue_rHB & iyear <=2003) 
          {//q_rate_fcn_rHB(iyear)=(1.0+q_rate)*q_rate_fcn_rHB(iyear-1); //compound
             q_rate_fcn_rHB(iyear)=(1.0+(iyear-styr_cpue_rHB)*q_rate)*q_rate_fcn_rHB(styr_cpue_rHB);  //linear
          }
          if (iyear>2003) {q_rate_fcn_rHB(iyear)=q_rate_fcn_rHB(iyear-1);} 
      }      
  } //end q_rate conditional      
 //Get density dependence scalar (=1.0 if density independent model is used)   
  if (q_DD_beta>0.0) 
  {
    B_q_DD+=dzero;
    for (iyear=styr;iyear<=endyr;iyear++)
        {q_DD_fcn(iyear)=pow(B0_q_DD,q_DD_beta)*pow(B_q_DD(iyear),-q_DD_beta);}
          //{q_DD_fcn(iyear)=1.0+4.0/(1.0+mfexp(0.75*(B_q_DD(iyear)-0.1*B0_q_DD))); }
  }  
}

void model_parameters::get_indices(void)
{
 //rHB  cpue
  q_rHB(styr_cpue_rHB)=mfexp(log_q_cpue_rHB); 
  for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++)
  {   
      N_rHB(iyear)=elem_prod(N_mdyr(iyear),sel_rHB(iyear)); 
      pred_rHB_cpue(iyear)=q_rHB(iyear)*q_rate_fcn_rHB(iyear)*q_DD_fcn(iyear)*sum(N_rHB(iyear));
      if (iyear<endyr_cpue_rHB){q_rHB(iyear+1)=q_rHB(iyear)*mfexp(q_RW_log_dev_rHB(iyear));}
  }
}

void model_parameters::get_length_comps(void)
{
  //cGN lines 
  for (iyear=1;iyear<=nyr_lenc_pool_cGN;iyear++)
   {pred_cGN_lenc_yr(iyear)=(L_cGN_num(yrs_lenc_pool_cGN(iyear))*lenprob_cGN)/sum(L_cGN_num(yrs_lenc_pool_cGN(iyear)));}
  pred_cGN_lenc.initialize();
  for (iyear=1;iyear<=nyr_lenc_pool_cGN;iyear++)
    {pred_cGN_lenc(1) += nfish_lenc_pool_cGN(iyear) * pred_cGN_lenc_yr(iyear);}
  pred_cGN_lenc(1)=pred_cGN_lenc(1)/sum(nfish_lenc_pool_cGN);
  ////cGN longline 
  //for (iyear=1;iyear<=nyr_cLL_lenc;iyear++)
  //{pred_cLL_lenc(iyear)=(L_cLL_num(yrs_cLL_lenc(iyear))*lenprob_cLL)/sum(L_cLL_num(yrs_cLL_lenc(iyear)));}
  //general rec   
  //for (iyear=1;iyear<=nyr_rGN_lenc;iyear++) 
  //{pred_rGN_lenc(iyear)=(L_rGN_num(yrs_rGN_lenc(iyear))*lenprob_rGN)/sum(L_rGN_num(yrs_rGN_lenc(iyear)));}
}

void model_parameters::get_age_comps(void)
{
  //Recreational 
  for (iyear=1;iyear<=nyr_agec_rGN;iyear++)
  {
    ErrorFree_rGN_agec(iyear)=L_rGN_num(yrs_agec_rGN(iyear))/sum(L_rGN_num(yrs_agec_rGN(iyear)));
    pred_rGN_agec_allages(iyear)=age_error*ErrorFree_rGN_agec(iyear); 
    for (iage=1; iage<=nages_agec; iage++) {pred_rGN_agec(iyear,iage)=pred_rGN_agec_allages(iyear,iage);} 
    //for (iage=(nages_agec+1); iage<=nages; iage++) {pred_rGN_agec(iyear,nages_agec)+=pred_rGN_agec_allages(iyear,iage);} //plus group                        
  }
}

void model_parameters::get_weighted_current(void)
{
  F_temp_sum=0.0;
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cGN+
        sum(log_dev_F_L_cGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rGN+
        sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_cGN_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cGN+
        sum(log_dev_F_L_cGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_rGN_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rGN+
        sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  log_F_dev_end_cGN=sum(log_dev_F_L_cGN((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  //log_F_dev_end_cLL=sum(log_F_dev_cLL((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted; 
  log_F_dev_end_rGN=sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  
  F_end_L=sel_cGN(endyr)*mfexp(log_avg_F_L_cGN+log_F_dev_end_cGN)+
          //sel_cLL(endyr)*mfexp(log_avg_F_cLL+log_F_dev_end_cLL)+
		  sel_rGN(endyr)*mfexp(log_avg_F_L_rGN+log_F_dev_end_rGN);   
  F_end=F_end_L;
  F_end_apex=max(F_end);
  sel_wgted_tot=F_end/F_end_apex;
  sel_wgted_L=elem_prod(sel_wgted_tot, elem_div(F_end_L,F_end));
  wgt_wgted_L_denom=F_cGN_prop+F_rGN_prop;  //+F_rHB_prop+F_cLL_prop
  wgt_wgted_L_klb=F_cGN_prop/wgt_wgted_L_denom*wholewgt_cGN_klb(endyr)+ 
				  //F_cLL_prop/wgt_wgted_L_denom*wholewgt_cLL_klb(endyr)+ 
				  F_rGN_prop/wgt_wgted_L_denom*wholewgt_rGN_klb(endyr);                                      
}

void model_parameters::get_msy(void)
{
  //compute values as functions of F
  for(ff=1; ff<=n_iter_msy; ff++)
  {
    //uses fishery-weighted F's
    Z_age_msy=0.0;
    F_L_age_msy=0.0;
    F_L_age_msy=F_msy(ff)*sel_wgted_L; 
    Z_age_msy=M+F_L_age_msy;        
    N_age_msy(1)=1.0;
    for (iage=2; iage<=nages; iage++)
      {N_age_msy(iage)=N_age_msy(iage-1)*mfexp(-1.*Z_age_msy(iage-1));}
    N_age_msy(nages)=N_age_msy(nages)/(1.0-mfexp(-1.*Z_age_msy(nages)));
    N_age_msy_spawn(1,(nages-1))=elem_prod(N_age_msy(1,(nages-1)),
                                   mfexp((-1.*Z_age_msy(1,(nages-1)))*spawn_time_frac));                 
    N_age_msy_spawn(nages)=(N_age_msy_spawn(nages-1)*(mfexp(-1.*(Z_age_msy(nages-1)*(1.0-spawn_time_frac) + 
                            Z_age_msy(nages)*spawn_time_frac) )))/(1.0-mfexp(-1.*Z_age_msy(nages)));
    spr_msy(ff)=sum(elem_prod(N_age_msy_spawn,reprod));
    //R_eq(ff)=SR_eq_func(R0, steep, spr_msy(1), spr_msy(ff), BiasCor, SR_switch);
    R_eq(ff)=BiasCor*R0;
    if (R_eq(ff)<dzero) {R_eq(ff)=dzero;}    
    N_age_msy*=R_eq(ff);
    N_age_msy_spawn*=R_eq(ff);
    for (iage=1; iage<=nages; iage++)
    {
      L_age_msy(iage)=N_age_msy(iage)*(F_L_age_msy(iage)/Z_age_msy(iage))*
                      (1.-mfexp(-1.*Z_age_msy(iage)));             
    }
    SSB_eq(ff)=sum(elem_prod(N_age_msy_spawn,reprod));
	SSB_eq_knum(ff)=sum(elem_prod(N_age_msy_spawn,reprodknum));
	B_eq(ff)=sum(elem_prod(N_age_msy,wgt_mt));
    L_eq_klb(ff)=sum(elem_prod(L_age_msy,wgt_wgted_L_klb)); //in whole weight
    L_eq_knum(ff)=sum(L_age_msy)/1000.0;   
  }  
  msy_klb_out=max(L_eq_klb); //msy in whole weight 
  for(ff=1; ff<=n_iter_msy; ff++)
  {
   if(L_eq_klb(ff) == msy_klb_out) 
      {    
        SSB_msy_out=SSB_eq(ff);
        B_msy_out=B_eq(ff);
        R_msy_out=R_eq(ff);
        msy_knum_out=L_eq_knum(ff);      
        F_msy_out=F_msy(ff);  
        spr_msy_out=spr_msy(ff);      
      }
  }
}

void model_parameters::get_per_recruit_stuff(void)
{
  //static per-recruit stuff
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
  //compute SSB/R and YPR as functions of F
  for(ff=1; ff<=n_iter_spr; ff++)
  {
    //uses fishery-weighted F's, same as in MSY calculations
    Z_age_spr=0.0;
    F_L_age_spr=0.0;
    F_L_age_spr=F_spr(ff)*sel_wgted_L;
    Z_age_spr=M+F_L_age_spr;
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
      L_spr(ff)+=L_age_spr(iage)*wgt_wgted_L_klb(iage)*1000.0; //in lb whole wgt
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
  Z_age_spr=M+F_L_age_spr;
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
    }
  SSB_F30_out=sum(elem_prod(N_age_spr_spawn,reprod));
  SSB_F30_knum_out=sum(elem_prod(N_age_spr_spawn,reprodknum));
  B_F30_out=sum(elem_prod(N_age_spr,wgt_mt));
  L_F30_klb_out=sum(elem_prod(L_age_F30,wgt_wgted_L_klb)); //in whole weight
  L_F30_knum_out=sum(L_age_F30)/1000.0; 
 //F40 calcs
  rec=column(N,1);
  rec_mean=sum(rec(styr_rec_spr, endyr_rec_spr))/nyrs_rec_spr;
  R_F40_out=rec_mean;
  F_L_age_spr=F40_out*sel_wgted_L;
  Z_age_spr=M+F_L_age_spr;
  N_age_spr(1)=R_F40_out;
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
      L_age_F40(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));                 
    }
  SSB_F40_out=sum(elem_prod(N_age_spr_spawn,reprod));
  SSB_F40_knum_out=sum(elem_prod(N_age_spr_spawn,reprodknum));
  B_F40_out=sum(elem_prod(N_age_spr,wgt_mt));
  L_F40_klb_out=sum(elem_prod(L_age_F40,wgt_wgted_L_klb)); //in whole weight
  L_F40_knum_out=sum(L_age_F40)/1000.0; 
}

void model_parameters::get_miscellaneous_stuff(void)
{
  if(var_rec_dev>0.0)
   {sigma_rec_dev=sqrt(var_rec_dev);} //sample SD of predicted residuals (may not equal rec_sigma)  
   else{sigma_rec_dev=0.0;}
  len_cv=elem_div(len_sd,meanlen_TL);
  len_cv_L=elem_div(len_sd_L,meanlen_TL_L);
  len_cv_F=elem_div(len_sd_F,meanlen_TL_F);
  //compute total landings- and discards-at-age in 1000 fish and klb whole weight
  L_total_num.initialize();
  L_total_klb.initialize();
  L_total_knum_yr.initialize();
  L_total_klb_yr.initialize();  
  for(iyear=styr; iyear<=endyr; iyear++)
  {
        L_total_klb_yr(iyear)=pred_cGN_L_klb(iyear)+pred_rGN_L_klb(iyear);//+pred_rHB_L_klb(iyear)+pred_cLL_L_klb(iyear)
        L_total_knum_yr(iyear)=pred_cGN_L_knum(iyear)+pred_rGN_L_knum(iyear);//+pred_rHB_L_knum(iyear)+pred_cLL_L_knum(iyear)
        B(iyear)=elem_prod(N(iyear),wgt_mt);
        totN(iyear)=sum(N(iyear));
        totB(iyear)=sum(B(iyear));                           
  }
  L_total_num=L_cGN_num+L_rGN_num;//+L_rHB_num+L_cLL_num   //landings at age in number fish
  L_total_klb=L_cGN_klb+L_rGN_klb;//+L_rHB_klb+L_cLL_klb   //landings at age in klb whole weight
  //Time series of interest
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
	if(F40_out>0)
    {
	  FdF40=Fapex/F40_out;
	  FdF40_end_mean=Fend_mean/F40_out;
	}
  if(SSB_F40_out>0)
    {
      SdSSB_F40=SSB/SSB_F40_out;
	  Sdmsst_F40=SSB/(smsy2msst75*SSB_F40_out);
      SdSSB_F40_end=SdSSB_F40(endyr);
	  Sdmsst_F40_end=Sdmsst_F40(endyr);
    }  	
   //fill in log recruitment deviations for yrs they are nonzero
   for(iyear=styr_rec_dev; iyear<=endyr_rec_dev; iyear++)
     {log_rec_dev_output(iyear)=log_dev_rec(iyear);}
   //fill in log Nage deviations for ages they are nonzero (ages2+)
   for(iage=2; iage<=nages; iage++)
     {log_Nage_dev_output(iage)=log_dev_Nage(iage);}
}

void model_parameters::get_projection(void)
{
    switch(Fproj_switch){
       case 1: //F=Fcurrent
          F_reg_proj=Fend_mean;
          break;
       case 2: //F=Fmsy
          F_reg_proj=F_msy_out;
          break;
       case 3: //F=F30
          F_reg_proj=F30_out;
          break;     
       case 4: //F=F40
          F_reg_proj=F40_out;
          break;          		  
       default: // no such switch available
          cout << "Error in input: Projection switch Fproj_switch must be set to 1, 2, 3, or 4." << endl;
          cout << "Presently it is set to " << Fproj_switch <<"."<< endl;
          exit(0);          
   }
  N_proj(styr_proj)=N(endyr+1); //initial conditions computed previously
  for (iyear=styr_proj; iyear<=endyr_proj; iyear++) //recruitment follows S-R curve (with bias correction) exactly
  {     
        if (iyear<styr_regs) {F_proj(iyear)=Fend_mean;}
		else {F_proj(iyear)=Fproj_mult*F_reg_proj;}
		FL_age_proj=sel_wgted_L*F_proj(iyear);
        Z_proj(iyear)=M+FL_age_proj;//+FD_age_proj;
        N_spawn_proj(iyear)(1,nages)=elem_prod(N_proj(iyear)(1,nages),(mfexp(-1.*(Z_proj(iyear)(1,nages))*spawn_time_frac))); //peak spawning time
		SSB_proj(iyear)= sum(elem_prod(N_spawn_proj(iyear),reprod));
        B_proj(iyear)=sum(elem_prod(N_proj(iyear),wgt_mt)); //uses spawning weight
		for (iage=1; iage<=nages; iage++)
			{L_age_proj(iyear,iage)=N_proj(iyear,iage)*FL_age_proj(iage)*(1.-mfexp(-1.*Z_proj(iyear,iage)))/Z_proj(iyear,iage);
			}          
        L_knum_proj(iyear)=sum(L_age_proj(iyear))/1000.0;
	    L_klb_proj(iyear)=sum(elem_prod(L_age_proj(iyear),wgt_wgted_L_klb));     //in 1000 lb
		if (iyear<endyr_proj) {
			N_proj(iyear+1,1)=BiasCor*R0; //Changed to move away from an SR relationship
			//N_proj(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB_proj(iyear),SR_switch);
			N_proj(iyear+1)(2,nages)=++elem_prod(N_proj(iyear)(1,nages-1),(mfexp(-1.*Z_proj(iyear)(1,nages-1))));
			N_proj(iyear+1,nages)+=N_proj(iyear,nages)*mfexp(-1.*Z_proj(iyear,nages)); //plus group		
		}
  }
   R_proj=column(N_proj,1);                          
}

void model_parameters::evaluate_objective_function(void)
{
  //fval=square(xdum-9.0);
  fval=0.0;
  fval_data=0.0;  
  f_rHB_cpue=0.0;
  f_rHB_cpue=lk_lognormal(pred_rHB_cpue, obs_cpue_rHB, obs_cv_cpue_rHB, w_cpue_rHB);
  fval+=f_rHB_cpue;
  fval_data+=f_rHB_cpue;   
  //f_cGN_L in 1000 lb whole wgt  
  f_cGN_L=lk_lognormal(pred_cGN_L_klb(styr_L_cGN,endyr_L_cGN), obs_L_cGN(styr_L_cGN,endyr_L_cGN),
                      obs_cv_L_cGN(styr_L_cGN,endyr_L_cGN), w_L);
  fval+=f_cGN_L;
  fval_data+=f_cGN_L;
  //f_rGN_L in 1000 fish
  f_rGN_L=lk_lognormal(pred_rGN_L_knum(styr_L_rGN,endyr_L_rGN), obs_L_rGN(styr_L_rGN,endyr_L_rGN), 
                      obs_cv_L_rGN(styr_L_rGN,endyr_L_rGN), w_L);
  fval+=f_rGN_L;
  fval_data+=f_rGN_L;  
  //f_cGN_lenc
  //f_cGN_lenc=lk_robust_multinomial(nsamp_lenc_cGN, pred_cGN_lenc, obs_lenc_cGN, nyr_lenc_cGN, double(nlenbins), minSS_lenc_cGN, w_lenc_cGN);
  //f_cGN_lenc=lk_logistic_normal(nsamp_lenc_cGN, pred_cGN_lenc, obs_lenc_cGN, nyr_lenc_cGN, double(nlenbins), minSS_lenc_cGN);
  f_cGN_lenc=lk_dirichlet_multinomial(nsamp_lenc_cGN, pred_cGN_lenc, obs_lenc_cGN, nyr_lenc_cGN, double(nlenbins), minSS_lenc_cGN, log_dm_lenc_cGN);
  fval+=f_cGN_lenc;
  fval_data+=f_cGN_lenc;
  //f_rGN_agec=lk_robust_multinomial(nsamp_agec_rGN, pred_rGN_agec, obs_agec_rGN, nyr_agec_rGN, double(nages_agec), minSS_agec_rGN, w_agec_rGN);
  //f_rGN_agec=lk_logistic_normal(nsamp_agec_rGN, pred_rGN_agec, obs_agec_rGN, nyr_agec_rGN, double(nages_agec), minSS_agec_rGN);
  f_rGN_agec=lk_dirichlet_multinomial(nsamp_agec_rGN, pred_rGN_agec, obs_agec_rGN, nyr_agec_rGN, double(nages_agec), minSS_agec_rGN, log_dm_agec_rGN);
  fval+=f_rGN_agec;
  fval_data+=f_rGN_agec; 
  //Light penalty applied to log_dev_Nage for deviation from zero. If not estimated, this penalty equals zero.
  f_Nage_init=norm2(log_dev_Nage);        
  fval+=w_Nage_init*f_Nage_init;
  f_rec_dev=0.0;
  //rec_sigma_sq=square(rec_sigma);
  rec_logL_add=nyrs_rec*log(rec_sigma);
  f_rec_dev=(square(log_dev_rec(styr_rec_dev) + rec_sigma_sq/2.0)/(2.0*rec_sigma_sq));
  for(iyear=(styr_rec_dev+1); iyear<=endyr_rec_dev; iyear++)
  {f_rec_dev+=(square(log_dev_rec(iyear)-R_autocorr*log_dev_rec(iyear-1) + rec_sigma_sq/2.0)/
               (2.0*rec_sigma_sq));}
  f_rec_dev+=rec_logL_add;            
  fval+=w_rec*f_rec_dev;
  f_rec_dev_early=0.0; //possible extra constraint on early rec deviations
  if (w_rec_early>0.0)
    { if (styr_rec_dev<endyr_rec_phase1)
        {  
          for(iyear=styr_rec_dev; iyear<=endyr_rec_phase1; iyear++)
          //{f_rec_dev_early+=(square(log_dev_rec(iyear)-R_autocorr*log_dev_rec(iyear-1) + rec_sigma_sq/2.0)/
          //                  (2.0*rec_sigma_sq)) + rec_logL_add;}
          {f_rec_dev_early+=square(log_dev_rec(iyear));}
        }
  fval+=w_rec_early*f_rec_dev_early;
  }
  f_rec_dev_end=0.0; //possible extra constraint on ending rec deviations
  if (w_rec_end>0.0)
  { if (endyr_rec_phase2<endyr_rec_dev)
        {  
          for(iyear=(endyr_rec_phase2+1); iyear<=endyr_rec_dev; iyear++)
          //{f_rec_dev_end+=(square(log_dev_rec(iyear)-R_autocorr*log_dev_rec(iyear-1) + rec_sigma_sq/2.0)/
          //                 (2.0*rec_sigma_sq)) + rec_logL_add;}
          {f_rec_dev_end+=square(log_dev_rec(iyear));}
        }
      fval+=w_rec_end*f_rec_dev_end;
   }  
  //Ftune penalty: does not apply in last phase
  f_Ftune=0.0; 
  if (w_Ftune>0.0)
  {if (set_Ftune>0.0 && !last_phase()) {f_Ftune=square(Fapex(set_Ftune_yr)-set_Ftune);}
   fval+=w_Ftune*f_Ftune;
  }
  //Penalty if apical F exceeds 3.0
  f_fullF_constraint=0.0;
  if (w_fullF>0.0)
  {for (iyear=styr; iyear<=endyr; iyear++)
        {if(Fapex(iyear)>3.0) {f_fullF_constraint+=(mfexp(Fapex(iyear)-3.0)-1.0);}}
   fval+=w_fullF*f_fullF_constraint;
  }
 //Random walk components of fishery dependent indices
 //f_cGN_RWq_cpue=0.0;
 //for (iyear=styr_cGN_cpue; iyear<endyr_cGN_cpue; iyear++)
 //    {f_cGN_RWq_cpue+=square(q_RW_log_dev_cGN(iyear))/(2.0*set_RWq_var);}
 //fval+=f_cGN_RWq_cpue;
 //
 //f_cLL_RWq_cpue=0.0;
 //for (iyear=styr_cLL_cpue; iyear<endyr_cLL_cpue; iyear++)
 //    {f_cLL_RWq_cpue+=square(q_RW_log_dev_cLL(iyear))/(2.0*set_RWq_var);}
 //fval+=f_cLL_RWq_cpue; 
 //
 f_rHB_RWq_cpue=0.0;
 for (iyear=styr_cpue_rHB; iyear<endyr_cpue_rHB; iyear++)
     {f_rHB_RWq_cpue+=square(q_RW_log_dev_rHB(iyear))/(2.0*set_RWq_var);}
 fval+=f_rHB_RWq_cpue;     
  f_priors=0.0; 
  f_priors+=neg_log_prior(len_cv_val,set_len_cv(5),set_len_cv(6),set_len_cv(7));
  f_priors+=neg_log_prior(steep,set_steep(5),set_steep(6),set_steep(7)); 
  f_priors+=neg_log_prior(log_R0,set_log_R0(5),set_log_R0(6),set_log_R0(7)); 
  f_priors+=neg_log_prior(R_autocorr,set_R_autocorr(5),set_R_autocorr(6),set_R_autocorr(7));
  f_priors+=neg_log_prior(rec_sigma,set_rec_sigma(5),set_rec_sigma(6),set_rec_sigma(7));
  f_priors+=neg_log_prior(selpar_A50_cGN1,set_selpar_A50_cGN1(5), set_selpar_A50_cGN1(6), set_selpar_A50_cGN1(7));
  f_priors+=neg_log_prior(selpar_slope_cGN1,set_selpar_slope_cGN1(5), set_selpar_slope_cGN1(6), set_selpar_slope_cGN1(7));
  //f_priors+=neg_log_prior(selpar_A50_cGN2,set_selpar_A50_cGN2(5), set_selpar_A50_cGN2(6), set_selpar_A50_cGN2(7));
  //f_priors+=neg_log_prior(selpar_slope_cGN2,set_selpar_slope_cGN2(5), set_selpar_slope_cGN2(6), set_selpar_slope_cGN2(7));
  // f_priors+=neg_log_prior(selpar_A502_cGN2,set_selpar_A502_cGN2(5), set_selpar_A502_cGN2(6), set_selpar_A502_cGN2(7));
  // f_priors+=neg_log_prior(selpar_slope2_cGN2,set_selpar_slope2_cGN2(5), set_selpar_slope2_cGN2(6), set_selpar_slope2_cGN2(7));  
  //f_priors+=neg_log_prior(selpar_A50_cGN3,set_selpar_A50_cGN3(5), set_selpar_A50_cGN3(6), set_selpar_A50_cGN3(7));
  //f_priors+=neg_log_prior(selpar_slope_cGN3,set_selpar_slope_cGN3(5), set_selpar_slope_cGN3(6), set_selpar_slope_cGN3(7));  
  f_priors+=neg_log_prior(selpar_A50_rGN1,set_selpar_A50_rGN1(5), set_selpar_A50_rGN1(6), set_selpar_A50_rGN1(7));
  f_priors+=neg_log_prior(selpar_slope_rGN1,set_selpar_slope_rGN1(5), set_selpar_slope_rGN1(6), set_selpar_slope_rGN1(7));
  f_priors+=neg_log_prior(selpar_A50_rGN2,set_selpar_A50_rGN2(5), set_selpar_A50_rGN2(6), set_selpar_A50_rGN2(7));
  f_priors+=neg_log_prior(selpar_slope_rGN2,set_selpar_slope_rGN2(5), set_selpar_slope_rGN2(6), set_selpar_slope_rGN2(7));
   f_priors+=neg_log_prior(log_q_cpue_rHB,set_log_q_cpue_rHB(5),set_log_q_cpue_rHB(6),set_log_q_cpue_rHB(7));
  f_priors+=neg_log_prior(log_dm_lenc_cGN,set_log_dm_lenc_cGN(5),set_log_dm_lenc_cGN(6),set_log_dm_lenc_cGN(7));
  //f_priors+=neg_log_prior(log_dm_cLL_lc,set_log_dm_cLL_lc(5),set_log_dm_cLL_lc(6),set_log_dm_cLL_lc(7));
  f_priors+=neg_log_prior(log_dm_agec_rGN,set_log_dm_agec_rGN(5),set_log_dm_agec_rGN(6),set_log_dm_agec_rGN(7));
  //f_priors+=neg_log_prior(log_dm_rGN_lc,set_log_dm_rGN_lc(5),set_log_dm_rGN_lc(6),set_log_dm_rGN_lc(7));
  f_priors+=neg_log_prior(F_init,set_F_init(5),set_F_init(6),set_F_init(7));  
  fval+=f_priors;
}

dvar_vector model_parameters::logistic(const dvar_vector& ages, const dvariable& A50, const dvariable& slope)
{
  //ages=vector of ages, A50=age at 50% selectivity, slope=rate of increase
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1./(1.+mfexp(-1.*slope*(ages-A50))); //logistic;  
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
}

dvar_vector model_parameters::logistic_exponential(const dvar_vector& ages, const dvariable& A50, const dvariable& slope, const dvariable& sigma, const dvariable& joint)
{
  //ages=vector of ages, A50=age at 50% sel (ascending limb), slope=rate of increase, sigma=controls rate of descent (descending)                               
  //joint=age to join curves                                                                                                                                    
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
}

dvar_vector model_parameters::logistic_double(const dvar_vector& ages, const dvariable& A501, const dvariable& slope1, const dvariable& A502, const dvariable& slope2)
{
  //ages=vector of ages, A50=age at 50% selectivity, slope=rate of increase, A502=age at 50% decrease additive to A501, slope2=slope of decrease
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=elem_prod( (1./(1.+mfexp(-1.*slope1*(ages-A501)))),(1.-(1./(1.+mfexp(-1.*slope2*(ages-(A501+A502)))))) );     
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
}

dvar_vector model_parameters::logistic_joint(const dvar_vector& ages, const dvariable& A501, const dvariable& slope1, const dvariable& A502, const dvariable& slope2, const dvariable& satval, const dvariable& joint)
{
  //ages=vector of ages, A501=age at 50% sel (ascending limb), slope1=rate of increase,A502=age at 50% sel (descending), slope1=rate of increase (ascending), 
  //satval=saturation value of descending limb, joint=location in age vector to join curves (may equal age or age + 1 if age-0 is included)
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
}

dvar_vector model_parameters::gaussian_double(const dvar_vector& ages, const dvariable& peak, const dvariable& top, const dvariable& ascwid, const dvariable& deswid, const dvariable& init, const dvariable& final)
{
  //ages=vector of ages, peak=ascending inflection location (as logistic), top=width of plateau, ascwid=ascent width (as log(width))
  //deswid=descent width (as log(width))
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
}

dvariable model_parameters::SR_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& SSB, int func)
{
  //R0=virgin recruitment, h=steepness, spr_F0=spawners per recruit @ F=0, SSB=spawning biomass
  //func=1 for Beverton-Holt, 2 for Ricker
  RETURN_ARRAYS_INCREMENT();
  dvariable Recruits_Tmp;
  switch(func) {
    case 1: //Beverton-Holt
      Recruits_Tmp=((0.8*R0*h*SSB)/(0.2*R0*spr_F0*(1.0-h)+(h-0.2)*SSB));       
    break;
    case 2: //Ricker
      Recruits_Tmp=((SSB/spr_F0)*mfexp(h*(1-SSB/(R0*spr_F0))));       
    break;
	case 3: //None
      Recruits_Tmp=R0;      
    break;
  }
  RETURN_ARRAYS_DECREMENT();
  return Recruits_Tmp;
}

dvariable model_parameters::SR_eq_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& spr_F, const dvariable& BC, int func)
{
  //R0=virgin recruitment, h=steepness, spr_F0=spawners per recruit @ F=0, spr_F=spawners per recruit @ F, BC=bias correction
  //func=1 for Beverton-Holt, 2 for Ricker
  RETURN_ARRAYS_INCREMENT();
  dvariable Recruits_Tmp;
  switch(func) {
    case 1: //Beverton-Holt
      Recruits_Tmp=(R0/((5.0*h-1.0)*spr_F))*(BC*4.0*h*spr_F-spr_F0*(1.0-h));    
    break;
    case 2: //Ricker
      Recruits_Tmp=R0/(spr_F/spr_F0)*(1.0+log(BC*spr_F/spr_F0)/h);      
    break;
	case 3: //None
      Recruits_Tmp=BC*R0;      
    break;
  }
  RETURN_ARRAYS_DECREMENT();
  return Recruits_Tmp;
}

dvariable model_parameters::multinom_eff_N(const dvar_vector& pred_comp, const dvar_vector& obs_comp)
{
  //pred_comp=vector of predicted comps, obscomp=vector of observed comps
  dvariable EffN_Tmp; dvariable numer; dvariable denom;
  RETURN_ARRAYS_INCREMENT();
  numer=sum( elem_prod(pred_comp,(1.0-pred_comp)) );
  denom=sum( square(obs_comp-pred_comp) );
  if (denom>0.0) {EffN_Tmp=numer/denom;}
  else {EffN_Tmp=-missing;}                            
  RETURN_ARRAYS_DECREMENT();
  return EffN_Tmp;
}

dvariable model_parameters::lk_lognormal(const dvar_vector& pred, const dvar_vector& obs, const dvar_vector& cv, const dvariable& wgt_dat)
{
  //pred=vector of predicted vals, obs=vector of observed vals, cv=vector of CVs in arithmetic space, wgt_dat=constant scaling of CVs
  //small_number is small value to avoid log(0) during search
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.0001;
  dvar_vector var(cv.indexmin(),cv.indexmax()); //variance in log space
  var=log(1.0+square(cv/wgt_dat));   // convert cv in arithmetic space to variance in log space
  LkvalTmp=sum(0.5*elem_div(square(log(elem_div((pred+small_number),(obs+small_number)))),var) );
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;
}

dvariable model_parameters::lk_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const double& minSS, const dvariable& wgt_dat)
{
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, minSS=min N threshold, wgt_dat=scaling of N's
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
}

dvariable model_parameters::lk_robust_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS, const dvariable& wgt_dat)
{
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, mbin=number of bins, minSS=min N threshold, wgt_dat=scaling of N's
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.0001;
  LkvalTmp=0.0;
  dvar_matrix Eprime=elem_prod((1.0-obs_comp), obs_comp)+0.1/mbin; //E' of Francis 2011, p.1131  
  dvar_vector nsamp_wgt=nsamp*wgt_dat;
  //cout<<nsamp_wgt<<endl;
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {LkvalTmp+= sum(0.5*log(Eprime(ii))-log(small_number+mfexp(elem_div((-square(obs_comp(ii)-pred_comp(ii))) , (Eprime(ii)*2.0/nsamp_wgt(ii)) ))) );
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;
}

dvariable model_parameters::lk_dirichlet_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS, const dvariable& log_dir_par)
{
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, mbin=number of bins, minSS=min N threshold, wgt_dat=scaling of N's
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvariable small_number=0.00001;
  LkvalTmp=0.0; 
  dvar_vector nsamp_adjust=nsamp*mfexp(log_dir_par);
   //dvar_vector nsamp_adjust=mfexp(log_dir_par);
  for (int ii=1; ii<=ncomp; ii++)
  {
	if (nsamp(ii)>=minSS)
    {
		LkvalTmp-=gammln(nsamp_adjust(ii))-gammln(nsamp(ii)+nsamp_adjust(ii));
		LkvalTmp-=sum(gammln(nsamp(ii)*obs_comp(ii)+nsamp_adjust(ii)*pred_comp(ii)+small_number));
		LkvalTmp+=sum(gammln(nsamp_adjust(ii)*pred_comp(ii)+small_number));
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;
  // //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, mbin=number of bins, minSS=min N threshold, wgt_dat=scaling of N's
  // RETURN_ARRAYS_INCREMENT();
  // dvariable LkvalTmp;
  // LkvalTmp=0.0; 
  // dvar_vector nsamp_adjust=nsamp*mfexp(log_dir_par);
  // //dvar_vector nsamp_adjust=mfexp(log_dir_par);
  // for (int ii=1; ii<=ncomp; ii++)
  // {
	// if (nsamp(ii)>=minSS)
    // {
		// LkvalTmp-=gammln(nsamp_adjust(ii))-gammln(nsamp(ii)+nsamp_adjust(ii));
		// LkvalTmp-=sum(gammln(nsamp(ii)*obs_comp(ii)+nsamp_adjust(ii)*pred_comp(ii)));
        // LkvalTmp+=sum(gammln(nsamp_adjust(ii)*pred_comp(ii)));		
    // }
  // }  
  // RETURN_ARRAYS_DECREMENT();
  // return LkvalTmp;
}

dvariable model_parameters::lk_logistic_normal(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS)
{
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, mbin=number of bins, minSS=min N threshold
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
  dvariable year_count; //keeps track of years included in likelihood (i.e., that meet the sample size requirement)
  LkvalTmp=0.0;
  nu_sum_sq=0.0;
  year_count=0.0;
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {
		year_count+=1.0;
		nu_mean=sum( log(obs_plus(ii))-log(pred_plus(ii))  )/mbin;	//year-specific mean log residual
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
}

dvariable model_parameters::neg_log_prior(dvariable pred, const double& prior, dvariable var, int pdf)
{
  //prior=prior point estimate, var=variance (if negative, treated as CV in arithmetic space), pred=predicted value, pdf=prior type (1=none, 2=lognormal, 3=normal, 4=beta)
    dvariable LkvalTmp;
    dvariable alpha, beta, ab_iq;
    dvariable big_number=1e10;
    LkvalTmp=0.0;
    // compute generic pdf's
    switch(pdf) {
        case 1: //option to turn off prior
          LkvalTmp=0.0;
          break;
        case 2: // lognormal
          if(prior<=0.0) cout << "YIKES: Don't use a lognormal distn for a negative prior" << endl;
          else if(pred<=0) LkvalTmp=big_number=1e10;
          else {
            if(var<0.0) var=log(1.0+var*var) ;      // convert cv to variance on log scale
            LkvalTmp= 0.5*( square(log(pred/prior))/var + log(var) );
          }
        break;
        case 3: // normal
          if(var<0.0 && prior!=0.0) var=square(var*prior);       // convert cv to variance on observation scale
          else if(var<0.0 && prior==0.0) var=-var;               // cv not really appropriate if prior value equals zero
          LkvalTmp= 0.5*( square(pred-prior)/var + log(var) );
          break;
        case 4: // beta
          if(var<0.0) var=square(var*prior);          // convert cv to variance on observation scale
          if(prior<=0.0 || prior>=1.0) cout << "YIKES: Don't use a beta distn for a prior outside (0,1)" << endl;
          ab_iq=prior*(1.0-prior)/var - 1.0; alpha=prior*ab_iq; beta=(1.0-prior)*ab_iq;
          if(pred>=0 && pred<=1) LkvalTmp= (1.0-alpha)*log(pred)+(1.0-beta)*log(1.0-pred)-gammln(alpha+beta)+gammln(alpha)+gammln(beta);
          else LkvalTmp=big_number;
          break;
        default: // no such prior pdf currently available
          cout << "The prior must be either 1(lognormal), 2(normal), or 3(beta)." << endl;
          cout << "Presently it is " << pdf << endl;
          exit(0);
    }
    return LkvalTmp;
}

dvariable model_parameters::sdnr_multinomial(const double& ncomp, const dvar_vector& ages, const dvar_vector& nsamp, 
                                    const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const dvariable& wgt_dat)
{
  //ncomp=number of years of data, ages=vector of ages, nsamp=vector of N's, 
  //pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, wgt_dat=likelihood weight for data source
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
}

dvariable model_parameters::sdnr_lognormal(const dvar_vector& pred, const dvar_vector& obs, const dvar_vector& cv, const dvariable& wgt_dat)
{
  //nyr=number of years of data, pred=vector of predicted data, obs=vector of observed data, cv=vector of cv's, wgt_dat=likelihood weight for data source
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
}

void model_parameters::report(const dvector& gradients)
{
 adstring ad_tmp=initial_params::get_reportfile_name();
  ofstream report((char*)(adprogram_name + ad_tmp));
  if (!report)
  {
    cerr << "error trying to open report file"  << adprogram_name << ".rep";
    return;
  }
  if (last_phase())  
  {
      cout<<"start report"<<endl;
       //cout<<"xdum = "<<xdum<<endl;
	  get_weighted_current();
       cout<<"got weighted"<<endl;
      get_msy();
       cout<<"got msy"<<endl;
      get_per_recruit_stuff();
       cout<<"got per recruit"<<endl;  
      get_miscellaneous_stuff();
       cout<<"got misc stuff"<<endl;
	  get_projection();
	   cout<<"got projection"<<endl;
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
      cout << "F40=" << F40_out<< "   SSB.F40=" << SSB_F40_out <<endl;
      cout <<"F status="<<FdF_msy_end<<endl;
      cout <<"Pop status="<<SdSSB_msy_end<<endl;
      cout << "R0="<<R0<<endl;
      //cout << "xdum " << xdum << endl;
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;  
      report << "TotalLikelihood " << fval << endl;
      report << "N" << endl;
      report << N<<endl;
      report << "F" << endl;
      report << F <<endl;
      // report << "prob_belowsizelim_block3" << endl;
	  // report<<prob_belowsizelim_block3<<endl;
      sdnr_I_rHB=sdnr_lognormal(pred_rHB_cpue, obs_cpue_rHB, obs_cv_cpue_rHB, w_cpue_rHB);      
      //#################################################################################################
      //##  Passing parameters to vector for bounds check plotting
      //################################################################################################# 
       Linf_out(8)=Linf; Linf_out(1,7)=set_Linf; 
       K_out(8)=K; K_out(1,7)=set_K;
       t0_out(8)=t0; t0_out(1,7)=set_t0;
       len_cv_val_out(8)=len_cv_val; len_cv_val_out(1,7)=set_len_cv;
	   Linf_L_out(8)=Linf_L; Linf_L_out(1,7)=set_Linf_L; 
       K_L_out(8)=K_L; K_L_out(1,7)=set_K_L;
       t0_L_out(8)=t0_L; t0_L_out(1,7)=set_t0_L;
       len_cv_val_L_out(8)=len_cv_val_L; len_cv_val_L_out(1,7)=set_len_cv_L;
	   Linf_F_out(8)=Linf_F; Linf_F_out(1,7)=set_Linf_F; 
       K_F_out(8)=K_F; K_F_out(1,7)=set_K_F;
       t0_F_out(8)=t0_F; t0_F_out(1,7)=set_t0_F;
       len_cv_val_F_out(8)=len_cv_val_F; len_cv_val_F_out(1,7)=set_len_cv_F;
       log_R0_out(8)=log_R0; log_R0_out(1,7)=set_log_R0;
       M_constant_out(8)=M_constant; M_constant_out(1,7)=set_M_constant;
       steep_out(8)=steep; steep_out(1,7)=set_steep;
       rec_sigma_out(8)=rec_sigma; rec_sigma_out(1,7)=set_rec_sigma;
       R_autocorr_out(8)=R_autocorr; R_autocorr_out(1,7)=set_R_autocorr;
	   log_dm_cGN_lc_out(8)=log_dm_lenc_cGN; log_dm_cGN_lc_out(1,7)=set_log_dm_lenc_cGN;
	   //log_dm_cLL_lc_out(8)=log_dm_cLL_lc; log_dm_cLL_lc_out(1,7)=set_log_dm_cLL_lc;
	   //log_dm_rGN_lc_out(8)=log_dm_rGN_lc; log_dm_rGN_lc_out(1,7)=set_log_dm_rGN_lc;
       log_dm_rGN_ac_out(8)=log_dm_agec_rGN; log_dm_rGN_ac_out(1,7)=set_log_dm_agec_rGN;
       selpar_A50_cGN1_out(8)=selpar_A50_cGN1; selpar_A50_cGN1_out(1,7)=set_selpar_A50_cGN1;
       selpar_slope_cGN1_out(8)=selpar_slope_cGN1; selpar_slope_cGN1_out(1,7)=set_selpar_slope_cGN1;
	   selpar_A50_rGN1_out(8)=selpar_A50_rGN1; selpar_A50_rGN1_out(1,7)=set_selpar_A50_rGN1;
       selpar_slope_rGN1_out(8)=selpar_slope_rGN1; selpar_slope_rGN1_out(1,7)=set_selpar_slope_rGN1;
	   selpar_A50_rGN2_out(8)=selpar_A50_rGN2; selpar_A50_rGN2_out(1,7)=set_selpar_A50_rGN2;
       selpar_slope_rGN2_out(8)=selpar_slope_rGN2; selpar_slope_rGN2_out(1,7)=set_selpar_slope_rGN2;
       log_q_rHB_out(8)=log_q_cpue_rHB; log_q_rHB_out(1,7)=set_log_q_cpue_rHB;
       log_avg_F_cGN_out(8)=log_avg_F_L_cGN; log_avg_F_cGN_out(1,7)=set_log_avg_F_L_cGN;
	   //log_avg_F_cLL_out(8)=log_avg_F_cLL; log_avg_F_cLL_out(1,7)=set_log_avg_F_cLL;
       log_avg_F_rGN_out(8)=log_avg_F_L_rGN; log_avg_F_rGN_out(1,7)=set_log_avg_F_L_rGN;
       F_init_out(8)=F_init; F_init_out(1,7)=set_F_init;       
       log_rec_dev_out(styr_rec_dev, endyr_rec_dev)=log_dev_rec;
       log_F_dev_cGN_out(styr_L_cGN,endyr_L_cGN)=log_dev_F_L_cGN;
	   //log_F_dev_cLL_out(styr_cLL_L,endyr_cLL_L)=log_F_dev_cLL;
       log_F_dev_rGN_out(styr_L_rGN,endyr_L_rGN)=log_dev_F_L_rGN;
   #include "bam.cxx"   // write the R-compatible report
  } //endl last phase loop     
}

model_data::~model_data()
{}

model_parameters::~model_parameters()
{}

void model_parameters::final_calcs(void){}

#ifdef _BORLANDC_
  extern unsigned _stklen=10000U;
#endif


#ifdef __ZTC__
  extern unsigned int _stack=10000U;
#endif

  long int arrmblsize=0;

int main(int argc,char * argv[])
{
    ad_set_new_handler();
  ad_exit=&ad_boundf;
  time(&start);
  arrmblsize=20000000;
  gradient_structure::set_MAX_NVAR_OFFSET(1600);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(2000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(10000);
    gradient_structure::set_NO_DERIVATIVES();
#ifdef DEBUG
  #ifndef __SUNPRO_C
std::feclearexcept(FE_ALL_EXCEPT);
  #endif
  auto start = std::chrono::high_resolution_clock::now();
#endif
    gradient_structure::set_YES_SAVE_VARIABLES_VALUES();
    if (!arrmblsize) arrmblsize=15000000;
    model_parameters mp(arrmblsize,argc,argv);
    mp.iprint=10;
    mp.preliminary_calculations();
    mp.computations(argc,argv);
#ifdef DEBUG
  std::cout << endl << argv[0] << " elapsed time is " << std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() << " microseconds." << endl;
  #ifndef __SUNPRO_C
bool failedtest = false;
if (std::fetestexcept(FE_DIVBYZERO))
  { cerr << "Error: Detected division by zero." << endl; failedtest = true; }
if (std::fetestexcept(FE_INVALID))
  { cerr << "Error: Detected invalid argument." << endl; failedtest = true; }
if (std::fetestexcept(FE_OVERFLOW))
  { cerr << "Error: Detected overflow." << endl; failedtest = true; }
if (std::fetestexcept(FE_UNDERFLOW))
  { cerr << "Error: Detected underflow." << endl; }
if (failedtest) { std::abort(); } 
  #endif
#endif
    return 0;
}

extern "C"  {
  void ad_boundf(int i)
  {
    /* so we can stop here */
    exit(i);
  }
}
