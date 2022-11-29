#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
  #include "admodel.h"     // Include AD class definitions
  #include "admb2r.cpp"    // Include S-compatible output functions (needs preceding)
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
  endyr_selex_phase1.allocate("endyr_selex_phase1");
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
  selpar_n_yrs_wgted.allocate("selpar_n_yrs_wgted");
  set_BiasCor.allocate("set_BiasCor");
  styr_cpue_cHL.allocate("styr_cpue_cHL");
  endyr_cpue_cHL.allocate("endyr_cpue_cHL");
  obs_cpue_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL,"obs_cpue_cHL");
  obs_cv_cpue_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL,"obs_cv_cpue_cHL");
  styr_L_cHL.allocate("styr_L_cHL");
  endyr_L_cHL.allocate("endyr_L_cHL");
  obs_L_cHL.allocate(styr_L_cHL,endyr_L_cHL,"obs_L_cHL");
  obs_cv_L_cHL.allocate(styr_L_cHL,endyr_L_cHL,"obs_cv_L_cHL");
  nyr_lenc_cHL.allocate("nyr_lenc_cHL");
  yrs_lenc_cHL.allocate(1,nyr_lenc_cHL,"yrs_lenc_cHL");
  nsamp_lenc_cHL.allocate(1,nyr_lenc_cHL,"nsamp_lenc_cHL");
  nfish_lenc_cHL.allocate(1,nyr_lenc_cHL,"nfish_lenc_cHL");
  obs_lenc_cHL.allocate(1,nyr_lenc_cHL,1,nlenbins,"obs_lenc_cHL");
  nyr_agec_cHL.allocate("nyr_agec_cHL");
  yrs_agec_cHL.allocate(1,nyr_agec_cHL,"yrs_agec_cHL");
  nsamp_agec_cHL.allocate(1,nyr_agec_cHL,"nsamp_agec_cHL");
  nfish_agec_cHL.allocate(1,nyr_agec_cHL,"nfish_agec_cHL");
  obs_agec_cHL.allocate(1,nyr_agec_cHL,1,nages_agec,"obs_agec_cHL");
  styr_cpue_sTV.allocate("styr_cpue_sTV");
  endyr_cpue_sTV.allocate("endyr_cpue_sTV");
  obs_cpue_sTV.allocate(styr_cpue_sTV,endyr_cpue_sTV,"obs_cpue_sTV");
  obs_cv_cpue_sTV.allocate(styr_cpue_sTV,endyr_cpue_sTV,"obs_cv_cpue_sTV");
  nyr_lenc_sTV.allocate("nyr_lenc_sTV");
  yrs_lenc_sTV.allocate(1,nyr_lenc_sTV,"yrs_lenc_sTV");
  nsamp_lenc_sTV.allocate(1,nyr_lenc_sTV,"nsamp_lenc_sTV");
  nfish_lenc_sTV.allocate(1,nyr_lenc_sTV,"nfish_lenc_sTV");
  obs_lenc_sTV.allocate(1,nyr_lenc_sTV,1,nlenbins,"obs_lenc_sTV");
  nyr_agec_sTV.allocate("nyr_agec_sTV");
  yrs_agec_sTV.allocate(1,nyr_agec_sTV,"yrs_agec_sTV");
  nsamp_agec_sTV.allocate(1,nyr_agec_sTV,"nsamp_agec_sTV");
  nfish_agec_sTV.allocate(1,nyr_agec_sTV,"nfish_agec_sTV");
  obs_agec_sTV.allocate(1,nyr_agec_sTV,1,nages_agec,"obs_agec_sTV");
  styr_cpue_rHB.allocate("styr_cpue_rHB");
  endyr_cpue_rHB.allocate("endyr_cpue_rHB");
  obs_cpue_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"obs_cpue_rHB");
  obs_cv_cpue_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"obs_cv_cpue_rHB");
  styr_L_rHB.allocate("styr_L_rHB");
  endyr_L_rHB.allocate("endyr_L_rHB");
  obs_L_rHB.allocate(styr_L_rHB,endyr_L_rHB,"obs_L_rHB");
  obs_cv_L_rHB.allocate(styr_L_rHB,endyr_L_rHB,"obs_cv_L_rHB");
  nyr_lenc_rHB.allocate("nyr_lenc_rHB");
  yrs_lenc_rHB.allocate(1,nyr_lenc_rHB,"yrs_lenc_rHB");
  nsamp_lenc_rHB.allocate(1,nyr_lenc_rHB,"nsamp_lenc_rHB");
  nfish_lenc_rHB.allocate(1,nyr_lenc_rHB,"nfish_lenc_rHB");
  obs_lenc_rHB.allocate(1,nyr_lenc_rHB,1,nlenbins,"obs_lenc_rHB");
  nyr_agec_rHB.allocate("nyr_agec_rHB");
  yrs_agec_rHB.allocate(1,nyr_agec_rHB,"yrs_agec_rHB");
  nsamp_agec_rHB.allocate(1,nyr_agec_rHB,"nsamp_agec_rHB");
  nfish_agec_rHB.allocate(1,nyr_agec_rHB,"nfish_agec_rHB");
  obs_agec_rHB.allocate(1,nyr_agec_rHB,1,nages_agec,"obs_agec_rHB");
  styr_D_rHB.allocate("styr_D_rHB");
  endyr_D_rHB.allocate("endyr_D_rHB");
  obs_released_rHB.allocate(styr_D_rHB,endyr_D_rHB,"obs_released_rHB");
  obs_cv_D_rHB.allocate(styr_D_rHB,endyr_D_rHB,"obs_cv_D_rHB");
  nyr_lenc_rHB_D.allocate("nyr_lenc_rHB_D");
  yrs_lenc_rHB_D.allocate(1,nyr_lenc_rHB_D,"yrs_lenc_rHB_D");
  nsamp_lenc_rHB_D.allocate(1,nyr_lenc_rHB_D,"nsamp_lenc_rHB_D");
  nfish_lenc_rHB_D.allocate(1,nyr_lenc_rHB_D,"nfish_lenc_rHB_D");
  obs_lenc_rHB_D.allocate(1,nyr_lenc_rHB_D,1,nlenbins,"obs_lenc_rHB_D");
  styr_cpue_rGN.allocate("styr_cpue_rGN");
  endyr_cpue_rGN.allocate("endyr_cpue_rGN");
  obs_cpue_rGN.allocate(styr_cpue_rGN,endyr_cpue_rGN,"obs_cpue_rGN");
  obs_cv_cpue_rGN.allocate(styr_cpue_rGN,endyr_cpue_rGN,"obs_cv_cpue_rGN");
  styr_L_rGN.allocate("styr_L_rGN");
  endyr_L_rGN.allocate("endyr_L_rGN");
  obs_L_rGN.allocate(styr_L_rGN,endyr_L_rGN,"obs_L_rGN");
  obs_cv_L_rGN.allocate(styr_L_rGN,endyr_L_rGN,"obs_cv_L_rGN");
  nyr_lenc_rGN.allocate("nyr_lenc_rGN");
  yrs_lenc_rGN.allocate(1,nyr_lenc_rGN,"yrs_lenc_rGN");
  nsamp_lenc_rGN.allocate(1,nyr_lenc_rGN,"nsamp_lenc_rGN");
  nfish_lenc_rGN.allocate(1,nyr_lenc_rGN,"nfish_lenc_rGN");
  obs_lenc_rGN.allocate(1,nyr_lenc_rGN,1,nlenbins,"obs_lenc_rGN");
  styr_D_rGN.allocate("styr_D_rGN");
  endyr_D_rGN.allocate("endyr_D_rGN");
  obs_released_rGN.allocate(styr_D_rGN,endyr_D_rGN,"obs_released_rGN");
  obs_cv_D_rGN.allocate(styr_D_rGN,endyr_D_rGN,"obs_cv_D_rGN");
  set_Linf.allocate(1,7,"set_Linf");
  set_K.allocate(1,7,"set_K");
  set_t0.allocate(1,7,"set_t0");
  set_len_cv.allocate(1,7,"set_len_cv");
  set_Linf_L.allocate(1,7,"set_Linf_L");
  set_K_L.allocate(1,7,"set_K_L");
  set_t0_L.allocate(1,7,"set_t0_L");
  set_len_cv_L.allocate(1,7,"set_len_cv_L");
  set_Linf_sTV.allocate(1,7,"set_Linf_sTV");
  set_K_sTV.allocate(1,7,"set_K_sTV");
  set_t0_sTV.allocate(1,7,"set_t0_sTV");
  set_len_cv_sTV.allocate(1,7,"set_len_cv_sTV");
  set_M_constant.allocate(1,7,"set_M_constant");
  set_steep.allocate(1,7,"set_steep");
  set_log_R0.allocate(1,7,"set_log_R0");
  set_R_autocorr.allocate(1,7,"set_R_autocorr");
  set_rec_sigma.allocate(1,7,"set_rec_sigma");
  set_selpar_L50_cHL1.allocate(1,7,"set_selpar_L50_cHL1");
  set_selpar_slope_cHL1.allocate(1,7,"set_selpar_slope_cHL1");
  set_selpar_L50_sTV.allocate(1,7,"set_selpar_L50_sTV");
  set_selpar_slope_sTV.allocate(1,7,"set_selpar_slope_sTV");
  set_selpar_L51_rHB1.allocate(1,7,"set_selpar_L51_rHB1");
  set_selpar_slope1_rHB1.allocate(1,7,"set_selpar_slope1_rHB1");
  set_selpar_L52_rHB1.allocate(1,7,"set_selpar_L52_rHB1");
  set_selpar_slope2_rHB1.allocate(1,7,"set_selpar_slope2_rHB1");
  set_selpar_L50_rHB2.allocate(1,7,"set_selpar_L50_rHB2");
  set_selpar_slope_rHB2.allocate(1,7,"set_selpar_slope_rHB2");
  set_selpar_afull_rHB2.allocate(1,7,"set_selpar_afull_rHB2");
  set_selpar_sigma_rHB2.allocate(1,7,"set_selpar_sigma_rHB2");
  set_selpar_L50_rHB_D.allocate(1,7,"set_selpar_L50_rHB_D");
  set_selpar_slope_rHB_D.allocate(1,7,"set_selpar_slope_rHB_D");
  set_selpar_afull_rHB_D.allocate(1,7,"set_selpar_afull_rHB_D");
  set_selpar_sigma_rHB_D.allocate(1,7,"set_selpar_sigma_rHB_D");
  set_selpar_L51_rGN.allocate(1,7,"set_selpar_L51_rGN");
  set_selpar_slope1_rGN.allocate(1,7,"set_selpar_slope1_rGN");
  set_selpar_L52_rGN.allocate(1,7,"set_selpar_L52_rGN");
  set_selpar_slope2_rGN.allocate(1,7,"set_selpar_slope2_rGN");
  set_log_q_cpue_cHL.allocate(1,7,"set_log_q_cpue_cHL");
  set_log_q_cpue_sTV.allocate(1,7,"set_log_q_cpue_sTV");
  set_log_q_cpue_rHB.allocate(1,7,"set_log_q_cpue_rHB");
  set_log_q_cpue_rGN.allocate(1,7,"set_log_q_cpue_rGN");
  set_F_init.allocate(1,7,"set_F_init");
  set_log_avg_F_L_cHL.allocate(1,7,"set_log_avg_F_L_cHL");
  set_log_avg_F_L_rHB.allocate(1,7,"set_log_avg_F_L_rHB");
  set_log_avg_F_L_rGN.allocate(1,7,"set_log_avg_F_L_rGN");
  set_log_avg_F_D_rHB.allocate(1,7,"set_log_avg_F_D_rHB");
  set_log_avg_F_D_rGN.allocate(1,7,"set_log_avg_F_D_rGN");
  set_log_dev_F_L_cHL.allocate(1,3,"set_log_dev_F_L_cHL");
  set_log_dev_F_L_rHB.allocate(1,3,"set_log_dev_F_L_rHB");
  set_log_dev_F_L_rGN.allocate(1,3,"set_log_dev_F_L_rGN");
  set_log_dev_F_D_rHB.allocate(1,3,"set_log_dev_F_D_rHB");
  set_log_dev_F_D_rGN.allocate(1,3,"set_log_dev_F_D_rGN");
  set_log_dev_rec.allocate(1,3,"set_log_dev_rec");
  set_log_dev_Nage.allocate(1,3,"set_log_dev_Nage");
  set_log_dev_vals_F_L_cHL.allocate(styr_L_cHL,endyr_L_cHL,"set_log_dev_vals_F_L_cHL");
  set_log_dev_vals_F_L_rHB.allocate(styr_L_rHB,endyr_L_rHB,"set_log_dev_vals_F_L_rHB");
  set_log_dev_vals_F_L_rGN.allocate(styr_L_rGN,endyr_L_rGN,"set_log_dev_vals_F_L_rGN");
  set_log_dev_vals_F_D_rHB.allocate(styr_D_rHB,endyr_D_rHB,"set_log_dev_vals_F_D_rHB");
  set_log_dev_vals_F_D_rGN.allocate(styr_D_rGN,endyr_D_rGN,"set_log_dev_vals_F_D_rGN");
  set_log_dev_vals_rec.allocate(styr_rec_dev,endyr_rec_dev,"set_log_dev_vals_rec");
  set_log_dev_vals_Nage.allocate(2,nages,"set_log_dev_vals_Nage");
  set_w_L.allocate("set_w_L");
  set_w_D.allocate("set_w_D");
  set_w_cpue_cHL.allocate("set_w_cpue_cHL");
  set_w_cpue_sTV.allocate("set_w_cpue_sTV");
  set_w_cpue_rHB.allocate("set_w_cpue_rHB");
  set_w_cpue_rGN.allocate("set_w_cpue_rGN");
  set_w_lenc_cHL.allocate("set_w_lenc_cHL");
  set_w_lenc_sTV.allocate("set_w_lenc_sTV");
  set_w_lenc_rHB.allocate("set_w_lenc_rHB");
  set_w_lenc_rHB_D.allocate("set_w_lenc_rHB_D");
  set_w_lenc_rGN.allocate("set_w_lenc_rGN");
  set_w_agec_cHL.allocate("set_w_agec_cHL");
  set_w_agec_sTV.allocate("set_w_agec_sTV");
  set_w_agec_rHB.allocate("set_w_agec_rHB");
  set_w_Nage_init.allocate("set_w_Nage_init");
  set_w_rec.allocate("set_w_rec");
  set_w_rec_early.allocate("set_w_rec_early");
  set_w_rec_end.allocate("set_w_rec_end");
  set_w_fullF.allocate("set_w_fullF");
  set_w_Ftune.allocate("set_w_Ftune");
  wgtpar_a.allocate("wgtpar_a");
  wgtpar_b.allocate("wgtpar_b");
  fecpar_a.allocate("fecpar_a");
  fecpar_b.allocate("fecpar_b");
  fecpar_batches.allocate(1,nages,"fecpar_batches");
  fecpar_scale.allocate("fecpar_scale");
  obs_maturity_f.allocate(1,nages,"obs_maturity_f");
  obs_prop_f.allocate(1,nages,"obs_prop_f");
  spawn_time_frac.allocate("spawn_time_frac");
  set_M.allocate(1,nages,"set_M");
  max_obs_age.allocate("max_obs_age");
  set_Dmort_rHB.allocate("set_Dmort_rHB");
  set_Dmort_rGN.allocate("set_Dmort_rGN");
  SR_switch.allocate("SR_switch");
  set_q_rate_phase.allocate("set_q_rate_phase");
  set_q_rate.allocate("set_q_rate");
  set_q_DD_phase.allocate("set_q_DD_phase");
  set_q_DD_beta.allocate("set_q_DD_beta");
  set_q_DD_beta_se.allocate("set_q_DD_beta_se");
  set_q_DD_stage.allocate("set_q_DD_stage");
  set_q_RW_phase.allocate("set_q_RW_phase");
  set_q_RW_rec_var.allocate("set_q_RW_rec_var");
  set_Ftune.allocate("set_Ftune");
  set_Ftune_yr.allocate("set_Ftune_yr");
  minSS_lenc_cHL.allocate("minSS_lenc_cHL");
  minSS_lenc_sTV.allocate("minSS_lenc_sTV");
  minSS_lenc_rHB.allocate("minSS_lenc_rHB");
  minSS_lenc_rHB_D.allocate("minSS_lenc_rHB_D");
  minSS_lenc_rGN.allocate("minSS_lenc_rGN");
  minSS_agec_cHL.allocate("minSS_agec_cHL");
  minSS_agec_sTV.allocate("minSS_agec_sTV");
  minSS_agec_rHB.allocate("minSS_agec_rHB");
  age_error.allocate(1,nages,1,nages,"age_error");
  use_landings_growth.allocate("use_landings_growth");
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
  //POPULATION growth
  const double Linf_LO=set_Linf(2); const double Linf_HI=set_Linf(3); const double Linf_PH=set_Linf(4);
  const double K_LO=set_K(2); const double K_HI=set_K(3); const double K_PH=set_K(4);
  const double t0_LO=set_t0(2); const double t0_HI=set_t0(3); const double t0_PH=set_t0(4);  
  const double len_cv_LO=set_len_cv(2); const double len_cv_HI=set_len_cv(3); const double len_cv_PH=set_len_cv(4); 
  //FISHERY LANDINGS growth
  const double Linf_L_LO=set_Linf_L(2); const double Linf_L_HI=set_Linf_L(3); const double Linf_L_PH=set_Linf_L(4);
  const double K_L_LO=set_K_L(2); const double K_L_HI=set_K_L(3); const double K_L_PH=set_K_L(4);
  const double t0_L_LO=set_t0_L(2); const double t0_L_HI=set_t0_L(3); const double t0_L_PH=set_t0_L(4);  
  const double len_cv_L_LO=set_len_cv_L(2); const double len_cv_L_HI=set_len_cv_L(3); const double len_cv_L_PH=set_len_cv_L(4);
  //FISHERY INDEPENDENT growth (SERFS chevron trap/video)
  const double Linf_sTV_LO=set_Linf_sTV(2); const double Linf_sTV_HI=set_Linf_sTV(3); const double Linf_sTV_PH=set_Linf_sTV(4);
  const double K_sTV_LO=set_K_sTV(2); const double K_sTV_HI=set_K_sTV(3); const double K_sTV_PH=set_K_sTV(4);
  const double t0_sTV_LO=set_t0_sTV(2); const double t0_sTV_HI=set_t0_sTV(3); const double t0_sTV_PH=set_t0_sTV(4);  
  const double len_cv_sTV_LO=set_len_cv_sTV(2); const double len_cv_sTV_HI=set_len_cv_sTV(3); const double len_cv_sTV_PH=set_len_cv_sTV(4);
  
  const double M_constant_LO=set_M_constant(2); const double M_constant_HI=set_M_constant(3); const double M_constant_PH=set_M_constant(4);        
  const double steep_LO=set_steep(2); const double steep_HI=set_steep(3); const double steep_PH=set_steep(4);
  const double log_R0_LO=set_log_R0(2); const double log_R0_HI=set_log_R0(3); const double log_R0_PH=set_log_R0(4);
  const double R_autocorr_LO=set_R_autocorr(2); const double R_autocorr_HI=set_R_autocorr(3); const double R_autocorr_PH=set_R_autocorr(4);
  const double rec_sigma_LO=set_rec_sigma(2); const double rec_sigma_HI=set_rec_sigma(3); const double rec_sigma_PH=set_rec_sigma(4);
  const double selpar_L50_cHL1_LO=set_selpar_L50_cHL1(2); const double selpar_L50_cHL1_HI=set_selpar_L50_cHL1(3); const double selpar_L50_cHL1_PH=set_selpar_L50_cHL1(4);
  const double selpar_slope_cHL1_LO=set_selpar_slope_cHL1(2); const double selpar_slope_cHL1_HI=set_selpar_slope_cHL1(3); const double selpar_slope_cHL1_PH=set_selpar_slope_cHL1(4);
  const double selpar_L50_sTV_LO=set_selpar_L50_sTV(2); const double selpar_L50_sTV_HI=set_selpar_L50_sTV(3); const double selpar_L50_sTV_PH=set_selpar_L50_sTV(4);
  const double selpar_slope_sTV_LO=set_selpar_slope_sTV(2); const double selpar_slope_sTV_HI=set_selpar_slope_sTV(3); const double selpar_slope_sTV_PH=set_selpar_slope_sTV(4); 
  const double selpar_L51_rHB1_LO=set_selpar_L51_rHB1(2); const double selpar_L51_rHB1_HI=set_selpar_L51_rHB1(3); const double selpar_L51_rHB1_PH=set_selpar_L51_rHB1(4);
  const double selpar_slope1_rHB1_LO=set_selpar_slope1_rHB1(2); const double selpar_slope1_rHB1_HI=set_selpar_slope1_rHB1(3); const double selpar_slope1_rHB1_PH=set_selpar_slope1_rHB1(4);
  const double selpar_L52_rHB1_LO=set_selpar_L52_rHB1(2); const double selpar_L52_rHB1_HI=set_selpar_L52_rHB1(3); const double selpar_L52_rHB1_PH=set_selpar_L52_rHB1(4);
  const double selpar_slope2_rHB1_LO=set_selpar_slope2_rHB1(2); const double selpar_slope2_rHB1_HI=set_selpar_slope2_rHB1(3); const double selpar_slope2_rHB1_PH=set_selpar_slope2_rHB1(4);
  const double selpar_L50_rHB2_LO=set_selpar_L50_rHB2(2); const double selpar_L50_rHB2_HI=set_selpar_L50_rHB2(3); const double selpar_L50_rHB2_PH=set_selpar_L50_rHB2(4);
  const double selpar_slope_rHB2_LO=set_selpar_slope_rHB2(2); const double selpar_slope_rHB2_HI=set_selpar_slope_rHB2(3); const double selpar_slope_rHB2_PH=set_selpar_slope_rHB2(4);
  const double selpar_afull_rHB2_LO=set_selpar_afull_rHB2(2); const double selpar_afull_rHB2_HI=set_selpar_afull_rHB2(3); const double selpar_afull_rHB2_PH=set_selpar_afull_rHB2(4);
  const double selpar_sigma_rHB2_LO=set_selpar_sigma_rHB2(2); const double selpar_sigma_rHB2_HI=set_selpar_sigma_rHB2(3); const double selpar_sigma_rHB2_PH=set_selpar_sigma_rHB2(4);
  const double selpar_L50_rHB_D_LO=set_selpar_L50_rHB_D(2); const double selpar_L50_rHB_D_HI=set_selpar_L50_rHB_D(3); const double selpar_L50_rHB_D_PH=set_selpar_L50_rHB_D(4);
  const double selpar_slope_rHB_D_LO=set_selpar_slope_rHB_D(2); const double selpar_slope_rHB_D_HI=set_selpar_slope_rHB_D(3); const double selpar_slope_rHB_D_PH=set_selpar_slope_rHB_D(4);
  const double selpar_afull_rHB_D_LO=set_selpar_afull_rHB_D(2); const double selpar_afull_rHB_D_HI=set_selpar_afull_rHB_D(3); const double selpar_afull_rHB_D_PH=set_selpar_afull_rHB_D(4);
  const double selpar_sigma_rHB_D_LO=set_selpar_sigma_rHB_D(2); const double selpar_sigma_rHB_D_HI=set_selpar_sigma_rHB_D(3); const double selpar_sigma_rHB_D_PH=set_selpar_sigma_rHB_D(4);
  const double selpar_L51_rGN_LO=set_selpar_L51_rGN(2); const double selpar_L51_rGN_HI=set_selpar_L51_rGN(3); const double selpar_L51_rGN_PH=set_selpar_L51_rGN(4); 
  const double selpar_slope1_rGN_LO=set_selpar_slope1_rGN(2); const double selpar_slope1_rGN_HI=set_selpar_slope1_rGN(3); const double selpar_slope1_rGN_PH=set_selpar_slope1_rGN(4);  
  const double selpar_L52_rGN_LO=set_selpar_L52_rGN(2); const double selpar_L52_rGN_HI=set_selpar_L52_rGN(3); const double selpar_L52_rGN_PH=set_selpar_L52_rGN(4);  
  const double selpar_slope2_rGN_LO=set_selpar_slope2_rGN(2); const double selpar_slope2_rGN_HI=set_selpar_slope2_rGN(3); const double selpar_slope2_rGN_PH=set_selpar_slope2_rGN(4); 
  const double log_q_cHL_LO=set_log_q_cpue_cHL(2); const double log_q_cHL_HI=set_log_q_cpue_cHL(3); const double log_q_cHL_PH=set_log_q_cpue_cHL(4);
  const double log_q_sTV_LO=set_log_q_cpue_sTV(2); const double log_q_sTV_HI=set_log_q_cpue_sTV(3); const double log_q_sTV_PH=set_log_q_cpue_sTV(4);  
  const double log_q_rHB_LO=set_log_q_cpue_rHB(2); const double log_q_rHB_HI=set_log_q_cpue_rHB(3); const double log_q_rHB_PH=set_log_q_cpue_rHB(4);
  const double log_q_rGN_LO=set_log_q_cpue_rGN(2); const double log_q_rGN_HI=set_log_q_cpue_rGN(3); const double log_q_rGN_PH=set_log_q_cpue_rGN(4);
  const double F_init_LO=set_F_init(2); const double F_init_HI=set_F_init(3); const double F_init_PH=set_F_init(4);
  const double log_avg_F_cHL_LO=set_log_avg_F_L_cHL(2); const double log_avg_F_cHL_HI=set_log_avg_F_L_cHL(3); const double log_avg_F_cHL_PH=set_log_avg_F_L_cHL(4);
  const double log_avg_F_rHB_LO=set_log_avg_F_L_rHB(2); const double log_avg_F_rHB_HI=set_log_avg_F_L_rHB(3); const double log_avg_F_rHB_PH=set_log_avg_F_L_rHB(4); 
  const double log_avg_F_rGN_LO=set_log_avg_F_L_rGN(2); const double log_avg_F_rGN_HI=set_log_avg_F_L_rGN(3); const double log_avg_F_rGN_PH=set_log_avg_F_L_rGN(4); 
  const double log_avg_F_rHB_D_LO=set_log_avg_F_D_rHB(2); const double log_avg_F_rHB_D_HI=set_log_avg_F_D_rHB(3); const double log_avg_F_rHB_D_PH=set_log_avg_F_D_rHB(4); 
  const double log_avg_F_rGN_D_LO=set_log_avg_F_D_rGN(2); const double log_avg_F_rGN_D_HI=set_log_avg_F_D_rGN(3); const double log_avg_F_rGN_D_PH=set_log_avg_F_D_rGN(4); 
  //-dev vectors-----------------------------------------------------------------------------------------------------------  
  const double log_F_dev_cHL_LO=set_log_dev_F_L_cHL(1); const double log_F_dev_cHL_HI=set_log_dev_F_L_cHL(2); const double log_F_dev_cHL_PH=set_log_dev_F_L_cHL(3);   
  const double log_F_dev_rHB_LO=set_log_dev_F_L_rHB(1); const double log_F_dev_rHB_HI=set_log_dev_F_L_rHB(2); const double log_F_dev_rHB_PH=set_log_dev_F_L_rHB(3);   
  const double log_F_dev_rGN_LO=set_log_dev_F_L_rGN(1); const double log_F_dev_rGN_HI=set_log_dev_F_L_rGN(2); const double log_F_dev_rGN_PH=set_log_dev_F_L_rGN(3);   
  const double log_F_dev_rHB_D_LO=set_log_dev_F_D_rHB(1); const double log_F_dev_rHB_D_HI=set_log_dev_F_D_rHB(2); const double log_F_dev_rHB_D_PH=set_log_dev_F_D_rHB(3);   
  const double log_F_dev_rGN_D_LO=set_log_dev_F_D_rGN(1); const double log_F_dev_rGN_D_HI=set_log_dev_F_D_rGN(2); const double log_F_dev_rGN_D_PH=set_log_dev_F_D_rGN(3);   
  const double log_rec_dev_LO=set_log_dev_rec(1); const double log_rec_dev_HI=set_log_dev_rec(2); const double log_rec_dev_PH=set_log_dev_rec(3);          
  const double log_Nage_dev_LO=set_log_dev_Nage(1); const double log_Nage_dev_HI=set_log_dev_Nage(2); const double log_Nage_dev_PH=set_log_dev_Nage(3);          
  Linf.allocate(Linf_LO,Linf_HI,Linf_PH,"Linf");
  K.allocate(K_LO,K_HI,K_PH,"K");
  t0.allocate(t0_LO,t0_HI,t0_PH,"t0");
  len_cv_val.allocate(len_cv_LO,len_cv_HI,len_cv_PH,"len_cv_val");
  Linf_L.allocate(Linf_L_LO,Linf_L_HI,Linf_L_PH,"Linf_L");
  K_L.allocate(K_L_LO,K_L_HI,K_L_PH,"K_L");
  t0_L.allocate(t0_L_LO,t0_L_HI,t0_L_PH,"t0_L");
  len_cv_val_L.allocate(len_cv_L_LO,len_cv_L_HI,len_cv_L_PH,"len_cv_val_L");
  Linf_sTV.allocate(Linf_sTV_LO,Linf_sTV_HI,Linf_sTV_PH,"Linf_sTV");
  K_sTV.allocate(K_sTV_LO,K_sTV_HI,K_sTV_PH,"K_sTV");
  t0_sTV.allocate(t0_sTV_LO,t0_sTV_HI,t0_sTV_PH,"t0_sTV");
  len_cv_val_sTV.allocate(len_cv_sTV_LO,len_cv_sTV_HI,len_cv_sTV_PH,"len_cv_val_sTV");
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
  Linf_sTV_out.allocate(1,8,"Linf_sTV_out");
  #ifndef NO_AD_INITIALIZE
    Linf_sTV_out.initialize();
  #endif
  K_sTV_out.allocate(1,8,"K_sTV_out");
  #ifndef NO_AD_INITIALIZE
    K_sTV_out.initialize();
  #endif
  t0_sTV_out.allocate(1,8,"t0_sTV_out");
  #ifndef NO_AD_INITIALIZE
    t0_sTV_out.initialize();
  #endif
  len_cv_val_sTV_out.allocate(1,8,"len_cv_val_sTV_out");
  #ifndef NO_AD_INITIALIZE
    len_cv_val_sTV_out.initialize();
  #endif
  meanlen_FL.allocate(1,nages,"meanlen_FL");
  #ifndef NO_AD_INITIALIZE
    meanlen_FL.initialize();
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
  fecundity.allocate(1,nages,"fecundity");
  #ifndef NO_AD_INITIALIZE
    fecundity.initialize();
  #endif
  meanlen_FL_L.allocate(1,nages,"meanlen_FL_L");
  #ifndef NO_AD_INITIALIZE
    meanlen_FL_L.initialize();
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
  meanlen_FL_sTV.allocate(1,nages,"meanlen_FL_sTV");
  #ifndef NO_AD_INITIALIZE
    meanlen_FL_sTV.initialize();
  #endif
  len_cHL_mm.allocate(styr,endyr,1,nages,"len_cHL_mm");
  #ifndef NO_AD_INITIALIZE
    len_cHL_mm.initialize();
  #endif
  wholewgt_cHL_klb.allocate(styr,endyr,1,nages,"wholewgt_cHL_klb");
  #ifndef NO_AD_INITIALIZE
    wholewgt_cHL_klb.initialize();
  #endif
  len_sTV_mm.allocate(styr,endyr,1,nages,"len_sTV_mm");
  #ifndef NO_AD_INITIALIZE
    len_sTV_mm.initialize();
  #endif
  wholewgt_sTV_klb.allocate(styr,endyr,1,nages,"wholewgt_sTV_klb");
  #ifndef NO_AD_INITIALIZE
    wholewgt_sTV_klb.initialize();
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
  len_rHB_D_mm.allocate(styr,endyr,1,nages,"len_rHB_D_mm");
  #ifndef NO_AD_INITIALIZE
    len_rHB_D_mm.initialize();
  #endif
  wholewgt_rHB_D_klb.allocate(styr,endyr,1,nages,"wholewgt_rHB_D_klb");
  #ifndef NO_AD_INITIALIZE
    wholewgt_rHB_D_klb.initialize();
  #endif
  len_rGN_D_mm.allocate(styr,endyr,1,nages,"len_rGN_D_mm");
  #ifndef NO_AD_INITIALIZE
    len_rGN_D_mm.initialize();
  #endif
  wholewgt_rGN_D_klb.allocate(styr,endyr,1,nages,"wholewgt_rGN_D_klb");
  #ifndef NO_AD_INITIALIZE
    wholewgt_rGN_D_klb.initialize();
  #endif
  lenprob.allocate(1,nages,1,nlenbins,"lenprob");
  #ifndef NO_AD_INITIALIZE
    lenprob.initialize();
  #endif
  lenprob2.allocate(1,nages,1,nlenbins,"lenprob2");
  #ifndef NO_AD_INITIALIZE
    lenprob2.initialize();
  #endif
  zscore_len.allocate("zscore_len");
  #ifndef NO_AD_INITIALIZE
  zscore_len.initialize();
  #endif
  zscore_len2.allocate("zscore_len2");
  #ifndef NO_AD_INITIALIZE
  zscore_len2.initialize();
  #endif
  cprob_lenvec.allocate(1,nlenbins,"cprob_lenvec");
  #ifndef NO_AD_INITIALIZE
    cprob_lenvec.initialize();
  #endif
  cprob_lenvec2.allocate(1,nlenbins,"cprob_lenvec2");
  #ifndef NO_AD_INITIALIZE
    cprob_lenvec2.initialize();
  #endif
  zscore_lzero.allocate("zscore_lzero");
  #ifndef NO_AD_INITIALIZE
  zscore_lzero.initialize();
  #endif
  cprob_lzero.allocate("cprob_lzero");
  #ifndef NO_AD_INITIALIZE
  cprob_lzero.initialize();
  #endif
  zscore_lzero2.allocate("zscore_lzero2");
  #ifndef NO_AD_INITIALIZE
  zscore_lzero2.initialize();
  #endif
  cprob_lzero2.allocate("cprob_lzero2");
  #ifndef NO_AD_INITIALIZE
  cprob_lzero2.initialize();
  #endif
  lenprob_cHL.allocate(1,nages,1,nlenbins,"lenprob_cHL");
  #ifndef NO_AD_INITIALIZE
    lenprob_cHL.initialize();
  #endif
  lenprob_sTV.allocate(1,nages,1,nlenbins,"lenprob_sTV");
  #ifndef NO_AD_INITIALIZE
    lenprob_sTV.initialize();
  #endif
  lenprob_rHB.allocate(1,nages,1,nlenbins,"lenprob_rHB");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB.initialize();
  #endif
  lenprob_rHB_D.allocate(1,nages,1,nlenbins,"lenprob_rHB_D");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB_D.initialize();
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
  len_sd_sTV.allocate(1,nages,"len_sd_sTV");
  #ifndef NO_AD_INITIALIZE
    len_sd_sTV.initialize();
  #endif
  len_cv_sTV.allocate(1,nages,"len_cv_sTV");
  #ifndef NO_AD_INITIALIZE
    len_cv_sTV.initialize();
  #endif
  pred_cHL_lenc.allocate(1,nyr_lenc_cHL,1,nlenbins,"pred_cHL_lenc");
  #ifndef NO_AD_INITIALIZE
    pred_cHL_lenc.initialize();
  #endif
  pred_sTV_lenc.allocate(1,nyr_lenc_sTV,1,nlenbins,"pred_sTV_lenc");
  #ifndef NO_AD_INITIALIZE
    pred_sTV_lenc.initialize();
  #endif
  pred_rHB_lenc.allocate(1,nyr_lenc_rHB,1,nlenbins,"pred_rHB_lenc");
  #ifndef NO_AD_INITIALIZE
    pred_rHB_lenc.initialize();
  #endif
  pred_rHB_D_lenc.allocate(1,nyr_lenc_rHB_D,1,nlenbins,"pred_rHB_D_lenc");
  #ifndef NO_AD_INITIALIZE
    pred_rHB_D_lenc.initialize();
  #endif
  pred_rGN_lenc.allocate(1,nyr_lenc_rGN,1,nlenbins,"pred_rGN_lenc");
  #ifndef NO_AD_INITIALIZE
    pred_rGN_lenc.initialize();
  #endif
  pred_cHL_agec.allocate(1,nyr_agec_cHL,1,nages_agec,"pred_cHL_agec");
  #ifndef NO_AD_INITIALIZE
    pred_cHL_agec.initialize();
  #endif
  pred_cHL_agec_allages.allocate(1,nyr_agec_cHL,1,nages,"pred_cHL_agec_allages");
  #ifndef NO_AD_INITIALIZE
    pred_cHL_agec_allages.initialize();
  #endif
  ErrorFree_cHL_agec.allocate(1,nyr_agec_cHL,1,nages,"ErrorFree_cHL_agec");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_cHL_agec.initialize();
  #endif
  pred_sTV_agec.allocate(1,nyr_agec_sTV,1,nages_agec,"pred_sTV_agec");
  #ifndef NO_AD_INITIALIZE
    pred_sTV_agec.initialize();
  #endif
  pred_sTV_agec_allages.allocate(1,nyr_agec_sTV,1,nages,"pred_sTV_agec_allages");
  #ifndef NO_AD_INITIALIZE
    pred_sTV_agec_allages.initialize();
  #endif
  ErrorFree_sTV_agec.allocate(1,nyr_agec_sTV,1,nages,"ErrorFree_sTV_agec");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_sTV_agec.initialize();
  #endif
  pred_rHB_agec.allocate(1,nyr_agec_rHB,1,nages_agec,"pred_rHB_agec");
  #ifndef NO_AD_INITIALIZE
    pred_rHB_agec.initialize();
  #endif
  pred_rHB_agec_allages.allocate(1,nyr_agec_rHB,1,nages,"pred_rHB_agec_allages");
  #ifndef NO_AD_INITIALIZE
    pred_rHB_agec_allages.initialize();
  #endif
  ErrorFree_rHB_agec.allocate(1,nyr_agec_rHB,1,nages,"ErrorFree_rHB_agec");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_rHB_agec.initialize();
  #endif
  nsamp_cHL_lenc_allyr.allocate(styr,endyr,"nsamp_cHL_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_cHL_lenc_allyr.initialize();
  #endif
  nsamp_sTV_lenc_allyr.allocate(styr,endyr,"nsamp_sTV_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_sTV_lenc_allyr.initialize();
  #endif
  nsamp_rHB_lenc_allyr.allocate(styr,endyr,"nsamp_rHB_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_rHB_lenc_allyr.initialize();
  #endif
  nsamp_rHB_D_lenc_allyr.allocate(styr,endyr,"nsamp_rHB_D_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_rHB_D_lenc_allyr.initialize();
  #endif
  nsamp_rGN_lenc_allyr.allocate(styr,endyr,"nsamp_rGN_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_rGN_lenc_allyr.initialize();
  #endif
  nsamp_cHL_agec_allyr.allocate(styr,endyr,"nsamp_cHL_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_cHL_agec_allyr.initialize();
  #endif
  nsamp_sTV_agec_allyr.allocate(styr,endyr,"nsamp_sTV_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_sTV_agec_allyr.initialize();
  #endif
  nsamp_rHB_agec_allyr.allocate(styr,endyr,"nsamp_rHB_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_rHB_agec_allyr.initialize();
  #endif
  nfish_cHL_lenc_allyr.allocate(styr,endyr,"nfish_cHL_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_cHL_lenc_allyr.initialize();
  #endif
  nfish_sTV_lenc_allyr.allocate(styr,endyr,"nfish_sTV_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_sTV_lenc_allyr.initialize();
  #endif
  nfish_rHB_lenc_allyr.allocate(styr,endyr,"nfish_rHB_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_rHB_lenc_allyr.initialize();
  #endif
  nfish_rHB_D_lenc_allyr.allocate(styr,endyr,"nfish_rHB_D_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_rHB_D_lenc_allyr.initialize();
  #endif
  nfish_rGN_lenc_allyr.allocate(styr,endyr,"nfish_rGN_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_rGN_lenc_allyr.initialize();
  #endif
  nfish_cHL_agec_allyr.allocate(styr,endyr,"nfish_cHL_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_cHL_agec_allyr.initialize();
  #endif
  nfish_sTV_agec_allyr.allocate(styr,endyr,"nfish_sTV_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_sTV_agec_allyr.initialize();
  #endif
  nfish_rHB_agec_allyr.allocate(styr,endyr,"nfish_rHB_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_rHB_agec_allyr.initialize();
  #endif
  neff_cHL_lenc_allyr_out.allocate(styr,endyr,"neff_cHL_lenc_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_cHL_lenc_allyr_out.initialize();
  #endif
  neff_sTV_lenc_allyr_out.allocate(styr,endyr,"neff_sTV_lenc_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_sTV_lenc_allyr_out.initialize();
  #endif
  neff_rHB_lenc_allyr_out.allocate(styr,endyr,"neff_rHB_lenc_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_rHB_lenc_allyr_out.initialize();
  #endif
  neff_rGN_lenc_allyr_out.allocate(styr,endyr,"neff_rGN_lenc_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_rGN_lenc_allyr_out.initialize();
  #endif
  neff_rHB_D_lenc_allyr_out.allocate(styr,endyr,"neff_rHB_D_lenc_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_rHB_D_lenc_allyr_out.initialize();
  #endif
  neff_cHL_agec_allyr_out.allocate(styr,endyr,"neff_cHL_agec_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_cHL_agec_allyr_out.initialize();
  #endif
  neff_sTV_agec_allyr_out.allocate(styr,endyr,"neff_sTV_agec_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_sTV_agec_allyr_out.initialize();
  #endif
  neff_rHB_agec_allyr_out.allocate(styr,endyr,"neff_rHB_agec_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_rHB_agec_allyr_out.initialize();
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
  MatFemB.allocate(styr,endyr,"MatFemB");
  #ifndef NO_AD_INITIALIZE
    MatFemB.initialize();
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
  reprod2.allocate(1,nages,"reprod2");
  #ifndef NO_AD_INITIALIZE
    reprod2.initialize();
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
  sel_cHL.allocate(styr,endyr,1,nages,"sel_cHL");
  #ifndef NO_AD_INITIALIZE
    sel_cHL.initialize();
  #endif
  selpar_L50_cHL1.allocate(selpar_L50_cHL1_LO,selpar_L50_cHL1_HI,selpar_L50_cHL1_PH,"selpar_L50_cHL1");
  selpar_slope_cHL1.allocate(selpar_slope_cHL1_LO,selpar_slope_cHL1_HI,selpar_slope_cHL1_PH,"selpar_slope_cHL1");
  selpar_L50_cHL1_out.allocate(1,8,"selpar_L50_cHL1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_L50_cHL1_out.initialize();
  #endif
  selpar_slope_cHL1_out.allocate(1,8,"selpar_slope_cHL1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cHL1_out.initialize();
  #endif
  sel_sTV.allocate(styr,endyr,1,nages,"sel_sTV");
  #ifndef NO_AD_INITIALIZE
    sel_sTV.initialize();
  #endif
  selpar_L50_sTV.allocate(selpar_L50_sTV_LO,selpar_L50_sTV_HI,selpar_L50_sTV_PH,"selpar_L50_sTV");
  selpar_slope_sTV.allocate(selpar_slope_sTV_LO,selpar_slope_sTV_HI,selpar_slope_sTV_PH,"selpar_slope_sTV");
  selpar_L50_sTV_out.allocate(1,8,"selpar_L50_sTV_out");
  #ifndef NO_AD_INITIALIZE
    selpar_L50_sTV_out.initialize();
  #endif
  selpar_slope_sTV_out.allocate(1,8,"selpar_slope_sTV_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_sTV_out.initialize();
  #endif
  sel_rHB.allocate(styr,endyr,1,nages,"sel_rHB");
  #ifndef NO_AD_INITIALIZE
    sel_rHB.initialize();
  #endif
  selpar_L51_rHB1.allocate(selpar_L51_rHB1_LO,selpar_L51_rHB1_HI,selpar_L51_rHB1_PH,"selpar_L51_rHB1");
  selpar_slope1_rHB1.allocate(selpar_slope1_rHB1_LO,selpar_slope1_rHB1_HI,selpar_slope1_rHB1_PH,"selpar_slope1_rHB1");
  selpar_L52_rHB1.allocate(selpar_L52_rHB1_LO,selpar_L52_rHB1_HI,selpar_L52_rHB1_PH,"selpar_L52_rHB1");
  selpar_slope2_rHB1.allocate(selpar_slope2_rHB1_LO,selpar_slope2_rHB1_HI,selpar_slope2_rHB1_PH,"selpar_slope2_rHB1");
  selpar_L50_rHB2.allocate(selpar_L50_rHB2_LO,selpar_L50_rHB2_HI,selpar_L50_rHB2_PH,"selpar_L50_rHB2");
  selpar_slope_rHB2.allocate(selpar_slope_rHB2_LO,selpar_slope_rHB2_HI,selpar_slope_rHB2_PH,"selpar_slope_rHB2");
  selpar_afull_rHB2.allocate(selpar_afull_rHB2_LO,selpar_afull_rHB2_HI,selpar_afull_rHB2_PH,"selpar_afull_rHB2");
  selpar_sigma_rHB2.allocate(selpar_sigma_rHB2_LO,selpar_sigma_rHB2_HI,selpar_sigma_rHB2_PH,"selpar_sigma_rHB2");
  selpar_L51_rHB1_out.allocate(1,8,"selpar_L51_rHB1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_L51_rHB1_out.initialize();
  #endif
  selpar_slope1_rHB1_out.allocate(1,8,"selpar_slope1_rHB1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope1_rHB1_out.initialize();
  #endif
  selpar_L52_rHB1_out.allocate(1,8,"selpar_L52_rHB1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_L52_rHB1_out.initialize();
  #endif
  selpar_slope2_rHB1_out.allocate(1,8,"selpar_slope2_rHB1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_rHB1_out.initialize();
  #endif
  selpar_L50_rHB2_out.allocate(1,8,"selpar_L50_rHB2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_L50_rHB2_out.initialize();
  #endif
  selpar_slope_rHB2_out.allocate(1,8,"selpar_slope_rHB2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rHB2_out.initialize();
  #endif
  selpar_afull_rHB2_out.allocate(1,8,"selpar_afull_rHB2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_afull_rHB2_out.initialize();
  #endif
  selpar_sigma_rHB2_out.allocate(1,8,"selpar_sigma_rHB2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_sigma_rHB2_out.initialize();
  #endif
  sel_rHB_D.allocate(styr,endyr,1,nages,"sel_rHB_D");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_D.initialize();
  #endif
  selpar_L50_rHB_D.allocate(selpar_L50_rHB_D_LO,selpar_L50_rHB_D_HI,selpar_L50_rHB_D_PH,"selpar_L50_rHB_D");
  selpar_slope_rHB_D.allocate(selpar_slope_rHB_D_LO,selpar_slope_rHB_D_HI,selpar_slope_rHB_D_PH,"selpar_slope_rHB_D");
  selpar_afull_rHB_D.allocate(selpar_afull_rHB_D_LO,selpar_afull_rHB_D_HI,selpar_afull_rHB_D_PH,"selpar_afull_rHB_D");
  selpar_sigma_rHB_D.allocate(selpar_sigma_rHB_D_LO,selpar_sigma_rHB_D_HI,selpar_sigma_rHB_D_PH,"selpar_sigma_rHB_D");
  selpar_L50_rHB_D_out.allocate(1,8,"selpar_L50_rHB_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_L50_rHB_D_out.initialize();
  #endif
  selpar_slope_rHB_D_out.allocate(1,8,"selpar_slope_rHB_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rHB_D_out.initialize();
  #endif
  selpar_afull_rHB_D_out.allocate(1,8,"selpar_afull_rHB_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_afull_rHB_D_out.initialize();
  #endif
  selpar_sigma_rHB_D_out.allocate(1,8,"selpar_sigma_rHB_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_sigma_rHB_D_out.initialize();
  #endif
  sel_rGN.allocate(styr,endyr,1,nages,"sel_rGN");
  #ifndef NO_AD_INITIALIZE
    sel_rGN.initialize();
  #endif
  selpar_L51_rGN.allocate(selpar_L51_rGN_LO,selpar_L51_rGN_HI,selpar_L51_rGN_PH,"selpar_L51_rGN");
  selpar_slope1_rGN.allocate(selpar_slope1_rGN_LO,selpar_slope1_rGN_HI,selpar_slope1_rGN_PH,"selpar_slope1_rGN");
  selpar_L52_rGN.allocate(selpar_L52_rGN_LO,selpar_L52_rGN_HI,selpar_L52_rGN_PH,"selpar_L52_rGN");
  selpar_slope2_rGN.allocate(selpar_slope2_rGN_LO,selpar_slope2_rGN_HI,selpar_slope2_rGN_PH,"selpar_slope2_rGN");
  selpar_L51_rGN_out.allocate(1,8,"selpar_L51_rGN_out");
  #ifndef NO_AD_INITIALIZE
    selpar_L51_rGN_out.initialize();
  #endif
  selpar_slope1_rGN_out.allocate(1,8,"selpar_slope1_rGN_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope1_rGN_out.initialize();
  #endif
  selpar_L52_rGN_out.allocate(1,8,"selpar_L52_rGN_out");
  #ifndef NO_AD_INITIALIZE
    selpar_L52_rGN_out.initialize();
  #endif
  selpar_slope2_rGN_out.allocate(1,8,"selpar_slope2_rGN_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_rGN_out.initialize();
  #endif
  sel_rGN_D.allocate(styr,endyr,1,nages,"sel_rGN_D");
  #ifndef NO_AD_INITIALIZE
    sel_rGN_D.initialize();
  #endif
  sel_wgted_L.allocate(1,nages,"sel_wgted_L");
  #ifndef NO_AD_INITIALIZE
    sel_wgted_L.initialize();
  #endif
  sel_wgted_D.allocate(1,nages,"sel_wgted_D");
  #ifndef NO_AD_INITIALIZE
    sel_wgted_D.initialize();
  #endif
  sel_wgted_tot.allocate(1,nages,"sel_wgted_tot");
  #ifndef NO_AD_INITIALIZE
    sel_wgted_tot.initialize();
  #endif
  pred_cHL_cpue.allocate(styr_cpue_cHL,endyr_cpue_cHL,"pred_cHL_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_cHL_cpue.initialize();
  #endif
  N_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL,1,nages,"N_cHL");
  #ifndef NO_AD_INITIALIZE
    N_cHL.initialize();
  #endif
  pred_sTV_cpue.allocate(styr_cpue_sTV,endyr_cpue_sTV,"pred_sTV_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_sTV_cpue.initialize();
  #endif
  N_sTV.allocate(styr_cpue_sTV,endyr_cpue_sTV,1,nages,"N_sTV");
  #ifndef NO_AD_INITIALIZE
    N_sTV.initialize();
  #endif
  pred_rHB_cpue.allocate(styr_cpue_rHB,endyr_cpue_rHB,"pred_rHB_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_rHB_cpue.initialize();
  #endif
  N_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,1,nages,"N_rHB");
  #ifndef NO_AD_INITIALIZE
    N_rHB.initialize();
  #endif
  pred_rGN_cpue.allocate(styr_cpue_rGN,endyr_cpue_rGN,"pred_rGN_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_rGN_cpue.initialize();
  #endif
  N_rGN.allocate(styr_cpue_rGN,endyr_cpue_rGN,1,nages,"N_rGN");
  #ifndef NO_AD_INITIALIZE
    N_rGN.initialize();
  #endif
  log_q_cpue_cHL.allocate(log_q_cHL_LO,log_q_cHL_HI,log_q_cHL_PH,"log_q_cpue_cHL");
  log_q_cpue_sTV.allocate(log_q_sTV_LO,log_q_sTV_HI,log_q_sTV_PH,"log_q_cpue_sTV");
  log_q_cpue_rHB.allocate(log_q_rHB_LO,log_q_rHB_HI,log_q_rHB_PH,"log_q_cpue_rHB");
  log_q_cpue_rGN.allocate(log_q_rGN_LO,log_q_rGN_HI,log_q_rGN_PH,"log_q_cpue_rGN");
  log_q_cHL_out.allocate(1,8,"log_q_cHL_out");
  #ifndef NO_AD_INITIALIZE
    log_q_cHL_out.initialize();
  #endif
  log_q_sTV_out.allocate(1,8,"log_q_sTV_out");
  #ifndef NO_AD_INITIALIZE
    log_q_sTV_out.initialize();
  #endif
  log_q_rHB_out.allocate(1,8,"log_q_rHB_out");
  #ifndef NO_AD_INITIALIZE
    log_q_rHB_out.initialize();
  #endif
  log_q_rGN_out.allocate(1,8,"log_q_rGN_out");
  #ifndef NO_AD_INITIALIZE
    log_q_rGN_out.initialize();
  #endif
  q_rate.allocate("q_rate");
  #ifndef NO_AD_INITIALIZE
  q_rate.initialize();
  #endif
  q_rate_fcn_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL,"q_rate_fcn_cHL");
  #ifndef NO_AD_INITIALIZE
    q_rate_fcn_cHL.initialize();
  #endif
  q_rate_fcn_sTV.allocate(styr_cpue_sTV,endyr_cpue_sTV,"q_rate_fcn_sTV");
  #ifndef NO_AD_INITIALIZE
    q_rate_fcn_sTV.initialize();
  #endif
  q_rate_fcn_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"q_rate_fcn_rHB");
  #ifndef NO_AD_INITIALIZE
    q_rate_fcn_rHB.initialize();
  #endif
  q_rate_fcn_rGN.allocate(styr_cpue_rGN,endyr_cpue_rGN,"q_rate_fcn_rGN");
  #ifndef NO_AD_INITIALIZE
    q_rate_fcn_rGN.initialize();
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
  q_RW_log_dev_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL-1,"q_RW_log_dev_cHL");
  #ifndef NO_AD_INITIALIZE
    q_RW_log_dev_cHL.initialize();
  #endif
  q_RW_log_dev_sTV.allocate(styr_cpue_sTV,endyr_cpue_sTV-1,"q_RW_log_dev_sTV");
  #ifndef NO_AD_INITIALIZE
    q_RW_log_dev_sTV.initialize();
  #endif
  q_RW_log_dev_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB-1,"q_RW_log_dev_rHB");
  #ifndef NO_AD_INITIALIZE
    q_RW_log_dev_rHB.initialize();
  #endif
  q_RW_log_dev_rGN.allocate(styr_cpue_rGN,endyr_cpue_rGN-1,"q_RW_log_dev_rGN");
  #ifndef NO_AD_INITIALIZE
    q_RW_log_dev_rGN.initialize();
  #endif
  q_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL,"q_cHL");
  #ifndef NO_AD_INITIALIZE
    q_cHL.initialize();
  #endif
  q_sTV.allocate(styr_cpue_sTV,endyr_cpue_sTV,"q_sTV");
  #ifndef NO_AD_INITIALIZE
    q_sTV.initialize();
  #endif
  q_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"q_rHB");
  #ifndef NO_AD_INITIALIZE
    q_rHB.initialize();
  #endif
  q_rGN.allocate(styr_cpue_rGN,endyr_cpue_rGN,"q_rGN");
  #ifndef NO_AD_INITIALIZE
    q_rGN.initialize();
  #endif
  L_cHL_num.allocate(styr,endyr,1,nages,"L_cHL_num");
  #ifndef NO_AD_INITIALIZE
    L_cHL_num.initialize();
  #endif
  L_cHL_klb.allocate(styr,endyr,1,nages,"L_cHL_klb");
  #ifndef NO_AD_INITIALIZE
    L_cHL_klb.initialize();
  #endif
  pred_cHL_L_knum.allocate(styr,endyr,"pred_cHL_L_knum");
  #ifndef NO_AD_INITIALIZE
    pred_cHL_L_knum.initialize();
  #endif
  pred_cHL_L_klb.allocate(styr,endyr,"pred_cHL_L_klb");
  #ifndef NO_AD_INITIALIZE
    pred_cHL_L_klb.initialize();
  #endif
  L_rHB_num.allocate(styr,endyr,1,nages,"L_rHB_num");
  #ifndef NO_AD_INITIALIZE
    L_rHB_num.initialize();
  #endif
  L_rHB_klb.allocate(styr,endyr,1,nages,"L_rHB_klb");
  #ifndef NO_AD_INITIALIZE
    L_rHB_klb.initialize();
  #endif
  pred_rHB_L_knum.allocate(styr,endyr,"pred_rHB_L_knum");
  #ifndef NO_AD_INITIALIZE
    pred_rHB_L_knum.initialize();
  #endif
  pred_rHB_L_klb.allocate(styr,endyr,"pred_rHB_L_klb");
  #ifndef NO_AD_INITIALIZE
    pred_rHB_L_klb.initialize();
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
  D_rHB_num.allocate(styr,endyr,1,nages,"D_rHB_num");
  #ifndef NO_AD_INITIALIZE
    D_rHB_num.initialize();
  #endif
  D_rHB_klb.allocate(styr,endyr,1,nages,"D_rHB_klb");
  #ifndef NO_AD_INITIALIZE
    D_rHB_klb.initialize();
  #endif
  pred_rHB_D_knum.allocate(styr_D_rHB,endyr_D_rHB,"pred_rHB_D_knum");
  #ifndef NO_AD_INITIALIZE
    pred_rHB_D_knum.initialize();
  #endif
  obs_rHB_D.allocate(styr_D_rHB,endyr_D_rHB,"obs_rHB_D");
  #ifndef NO_AD_INITIALIZE
    obs_rHB_D.initialize();
  #endif
  pred_rHB_D_klb.allocate(styr_D_rHB,endyr_D_rHB,"pred_rHB_D_klb");
  #ifndef NO_AD_INITIALIZE
    pred_rHB_D_klb.initialize();
  #endif
  D_rGN_num.allocate(styr,endyr,1,nages,"D_rGN_num");
  #ifndef NO_AD_INITIALIZE
    D_rGN_num.initialize();
  #endif
  D_rGN_klb.allocate(styr,endyr,1,nages,"D_rGN_klb");
  #ifndef NO_AD_INITIALIZE
    D_rGN_klb.initialize();
  #endif
  pred_rGN_D_knum.allocate(styr_D_rGN,endyr_D_rGN,"pred_rGN_D_knum");
  #ifndef NO_AD_INITIALIZE
    pred_rGN_D_knum.initialize();
  #endif
  obs_rGN_D.allocate(styr_D_rGN,endyr_D_rGN,"obs_rGN_D");
  #ifndef NO_AD_INITIALIZE
    obs_rGN_D.initialize();
  #endif
  pred_rGN_D_klb.allocate(styr_D_rGN,endyr_D_rGN,"pred_rGN_D_klb");
  #ifndef NO_AD_INITIALIZE
    pred_rGN_D_klb.initialize();
  #endif
  D_total_num.allocate(styr,endyr,1,nages,"D_total_num");
  #ifndef NO_AD_INITIALIZE
    D_total_num.initialize();
  #endif
  D_total_klb.allocate(styr,endyr,1,nages,"D_total_klb");
  #ifndef NO_AD_INITIALIZE
    D_total_klb.initialize();
  #endif
  D_total_knum_yr.allocate(styr,endyr,"D_total_knum_yr");
  #ifndef NO_AD_INITIALIZE
    D_total_knum_yr.initialize();
  #endif
  D_total_klb_yr.allocate(styr,endyr,"D_total_klb_yr");
  #ifndef NO_AD_INITIALIZE
    D_total_klb_yr.initialize();
  #endif
  Dmort_rHB.allocate("Dmort_rHB");
  #ifndef NO_AD_INITIALIZE
  Dmort_rHB.initialize();
  #endif
  Dmort_rGN.allocate("Dmort_rGN");
  #ifndef NO_AD_INITIALIZE
  Dmort_rGN.initialize();
  #endif
  F_cHL_prop.allocate("F_cHL_prop");
  #ifndef NO_AD_INITIALIZE
  F_cHL_prop.initialize();
  #endif
  F_rHB_prop.allocate("F_rHB_prop");
  #ifndef NO_AD_INITIALIZE
  F_rHB_prop.initialize();
  #endif
  F_rGN_prop.allocate("F_rGN_prop");
  #ifndef NO_AD_INITIALIZE
  F_rGN_prop.initialize();
  #endif
  F_rHB_D_prop.allocate("F_rHB_D_prop");
  #ifndef NO_AD_INITIALIZE
  F_rHB_D_prop.initialize();
  #endif
  F_rGN_D_prop.allocate("F_rGN_D_prop");
  #ifndef NO_AD_INITIALIZE
  F_rGN_D_prop.initialize();
  #endif
  F_init_cHL_prop.allocate("F_init_cHL_prop");
  #ifndef NO_AD_INITIALIZE
  F_init_cHL_prop.initialize();
  #endif
  F_init_rHB_prop.allocate("F_init_rHB_prop");
  #ifndef NO_AD_INITIALIZE
  F_init_rHB_prop.initialize();
  #endif
  F_init_rGN_prop.allocate("F_init_rGN_prop");
  #ifndef NO_AD_INITIALIZE
  F_init_rGN_prop.initialize();
  #endif
  F_init_rHB_D_prop.allocate("F_init_rHB_D_prop");
  #ifndef NO_AD_INITIALIZE
  F_init_rHB_D_prop.initialize();
  #endif
  F_init_rGN_D_prop.allocate("F_init_rGN_D_prop");
  #ifndef NO_AD_INITIALIZE
  F_init_rGN_D_prop.initialize();
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
  F_end_D.allocate(1,nages,"F_end_D");
  #ifndef NO_AD_INITIALIZE
    F_end_D.initialize();
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
  D_msy_klb_out.allocate("D_msy_klb_out");
  #ifndef NO_AD_INITIALIZE
  D_msy_klb_out.initialize();
  #endif
  D_msy_knum_out.allocate("D_msy_knum_out");
  #ifndef NO_AD_INITIALIZE
  D_msy_knum_out.initialize();
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
  D_F30_knum_out.allocate("D_F30_knum_out");
  #ifndef NO_AD_INITIALIZE
  D_F30_knum_out.initialize();
  #endif
  D_F30_klb_out.allocate("D_F30_klb_out");
  #ifndef NO_AD_INITIALIZE
  D_F30_klb_out.initialize();
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
  D_age_msy.allocate(1,nages,"D_age_msy");
  #ifndef NO_AD_INITIALIZE
    D_age_msy.initialize();
  #endif
  Z_age_msy.allocate(1,nages,"Z_age_msy");
  #ifndef NO_AD_INITIALIZE
    Z_age_msy.initialize();
  #endif
  F_L_age_msy.allocate(1,nages,"F_L_age_msy");
  #ifndef NO_AD_INITIALIZE
    F_L_age_msy.initialize();
  #endif
  F_D_age_msy.allocate(1,nages,"F_D_age_msy");
  #ifndef NO_AD_INITIALIZE
    F_D_age_msy.initialize();
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
  D_eq_klb.allocate(1,n_iter_msy,"D_eq_klb");
  #ifndef NO_AD_INITIALIZE
    D_eq_klb.initialize();
  #endif
  D_eq_knum.allocate(1,n_iter_msy,"D_eq_knum");
  #ifndef NO_AD_INITIALIZE
    D_eq_knum.initialize();
  #endif
  SSB_eq.allocate(1,n_iter_msy,"SSB_eq");
  #ifndef NO_AD_INITIALIZE
    SSB_eq.initialize();
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
  Fend_mean_temp.allocate("Fend_mean_temp");
  #ifndef NO_AD_INITIALIZE
  Fend_mean_temp.initialize();
  #endif
  Fend_mean.allocate("Fend_mean");
  #ifndef NO_AD_INITIALIZE
  Fend_mean.initialize();
  #endif
  wgt_wgted_L_klb.allocate(1,nages,"wgt_wgted_L_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_wgted_L_klb.initialize();
  #endif
  wgt_wgted_D_klb.allocate(1,nages,"wgt_wgted_D_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_wgted_D_klb.initialize();
  #endif
  wgt_wgted_L_denom.allocate("wgt_wgted_L_denom");
  #ifndef NO_AD_INITIALIZE
  wgt_wgted_L_denom.initialize();
  #endif
  wgt_wgted_D_denom.allocate("wgt_wgted_D_denom");
  #ifndef NO_AD_INITIALIZE
  wgt_wgted_D_denom.initialize();
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
  smsy2msst.allocate("smsy2msst");
  #ifndef NO_AD_INITIALIZE
  smsy2msst.initialize();
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
  log_avg_F_L_cHL.allocate(log_avg_F_cHL_LO,log_avg_F_cHL_HI,log_avg_F_cHL_PH,"log_avg_F_L_cHL");
  log_avg_F_cHL_out.allocate(1,8,"log_avg_F_cHL_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_cHL_out.initialize();
  #endif
  log_dev_F_L_cHL.allocate(styr_L_cHL,endyr_L_cHL,log_F_dev_cHL_LO,log_F_dev_cHL_HI,log_F_dev_cHL_PH,"log_dev_F_L_cHL");
  log_F_dev_cHL_out.allocate(styr_L_cHL,endyr_L_cHL,"log_F_dev_cHL_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_cHL_out.initialize();
  #endif
  F_cHL.allocate(styr,endyr,1,nages,"F_cHL");
  #ifndef NO_AD_INITIALIZE
    F_cHL.initialize();
  #endif
  F_cHL_out.allocate(styr,endyr,"F_cHL_out");
  #ifndef NO_AD_INITIALIZE
    F_cHL_out.initialize();
  #endif
  log_F_dev_init_cHL.allocate("log_F_dev_init_cHL");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_init_cHL.initialize();
  #endif
  log_F_dev_end_cHL.allocate("log_F_dev_end_cHL");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_cHL.initialize();
  #endif
  log_avg_F_L_rHB.allocate(log_avg_F_rHB_LO,log_avg_F_rHB_HI,log_avg_F_rHB_PH,"log_avg_F_L_rHB");
  log_avg_F_rHB_out.allocate(1,8,"log_avg_F_rHB_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_rHB_out.initialize();
  #endif
  log_dev_F_L_rHB.allocate(styr_L_rHB,endyr_L_rHB,log_F_dev_rHB_LO,log_F_dev_rHB_HI,log_F_dev_rHB_PH,"log_dev_F_L_rHB");
  log_F_dev_rHB_out.allocate(styr_L_rHB,endyr_L_rHB,"log_F_dev_rHB_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_rHB_out.initialize();
  #endif
  F_rHB.allocate(styr,endyr,1,nages,"F_rHB");
  #ifndef NO_AD_INITIALIZE
    F_rHB.initialize();
  #endif
  F_rHB_out.allocate(styr,endyr,"F_rHB_out");
  #ifndef NO_AD_INITIALIZE
    F_rHB_out.initialize();
  #endif
  log_F_dev_init_rHB.allocate("log_F_dev_init_rHB");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_init_rHB.initialize();
  #endif
  log_F_dev_end_rHB.allocate("log_F_dev_end_rHB");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_rHB.initialize();
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
  log_avg_F_D_rHB.allocate(log_avg_F_rHB_D_LO,log_avg_F_rHB_D_HI,log_avg_F_rHB_D_PH,"log_avg_F_D_rHB");
  log_avg_F_rHB_D_out.allocate(1,8,"log_avg_F_rHB_D_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_rHB_D_out.initialize();
  #endif
  log_dev_F_D_rHB.allocate(styr_D_rHB,endyr_D_rHB,log_F_dev_rHB_D_LO,log_F_dev_rHB_D_HI,log_F_dev_rHB_D_PH,"log_dev_F_D_rHB");
  log_F_dev_rHB_D_out.allocate(styr_D_rHB,endyr_D_rHB,"log_F_dev_rHB_D_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_rHB_D_out.initialize();
  #endif
  F_rHB_D.allocate(styr,endyr,1,nages,"F_rHB_D");
  #ifndef NO_AD_INITIALIZE
    F_rHB_D.initialize();
  #endif
  F_rHB_D_out.allocate(styr,endyr,"F_rHB_D_out");
  #ifndef NO_AD_INITIALIZE
    F_rHB_D_out.initialize();
  #endif
  log_F_dev_init_rHB_D.allocate("log_F_dev_init_rHB_D");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_init_rHB_D.initialize();
  #endif
  log_F_dev_end_rHB_D.allocate("log_F_dev_end_rHB_D");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_rHB_D.initialize();
  #endif
  log_avg_F_D_rGN.allocate(log_avg_F_rGN_D_LO,log_avg_F_rGN_D_HI,log_avg_F_rGN_D_PH,"log_avg_F_D_rGN");
  log_avg_F_rGN_D_out.allocate(1,8,"log_avg_F_rGN_D_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_rGN_D_out.initialize();
  #endif
  log_dev_F_D_rGN.allocate(styr_D_rGN,endyr_D_rGN,log_F_dev_rGN_D_LO,log_F_dev_rGN_D_HI,log_F_dev_rGN_D_PH,"log_dev_F_D_rGN");
  log_F_dev_rGN_D_out.allocate(styr_D_rGN,endyr_D_rGN,"log_F_dev_rGN_D_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_rGN_D_out.initialize();
  #endif
  F_rGN_D.allocate(styr,endyr,1,nages,"F_rGN_D");
  #ifndef NO_AD_INITIALIZE
    F_rGN_D.initialize();
  #endif
  F_rGN_D_out.allocate(styr,endyr,"F_rGN_D_out");
  #ifndef NO_AD_INITIALIZE
    F_rGN_D_out.initialize();
  #endif
  log_F_dev_init_rGN_D.allocate("log_F_dev_init_rGN_D");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_init_rGN_D.initialize();
  #endif
  log_F_dev_end_rGN_D.allocate("log_F_dev_end_rGN_D");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_rGN_D.initialize();
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
  sdnr_lc_cHL.allocate("sdnr_lc_cHL");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_cHL.initialize();
  #endif
  sdnr_lc_sTV.allocate("sdnr_lc_sTV");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_sTV.initialize();
  #endif
  sdnr_lc_rHB.allocate("sdnr_lc_rHB");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_rHB.initialize();
  #endif
  sdnr_lc_rHB_D.allocate("sdnr_lc_rHB_D");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_rHB_D.initialize();
  #endif
  sdnr_lc_rGN.allocate("sdnr_lc_rGN");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_rGN.initialize();
  #endif
  sdnr_ac_cHL.allocate("sdnr_ac_cHL");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_cHL.initialize();
  #endif
  sdnr_ac_sTV.allocate("sdnr_ac_sTV");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_sTV.initialize();
  #endif
  sdnr_ac_rHB.allocate("sdnr_ac_rHB");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_rHB.initialize();
  #endif
  sdnr_ac_rGN.allocate("sdnr_ac_rGN");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_rGN.initialize();
  #endif
  sdnr_I_cHL.allocate("sdnr_I_cHL");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_cHL.initialize();
  #endif
  sdnr_I_sTV.allocate("sdnr_I_sTV");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_sTV.initialize();
  #endif
  sdnr_I_rHB.allocate("sdnr_I_rHB");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_rHB.initialize();
  #endif
  sdnr_I_rGN.allocate("sdnr_I_rGN");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_rGN.initialize();
  #endif
  w_L.allocate("w_L");
  #ifndef NO_AD_INITIALIZE
  w_L.initialize();
  #endif
  w_D.allocate("w_D");
  #ifndef NO_AD_INITIALIZE
  w_D.initialize();
  #endif
  w_cpue_cHL.allocate("w_cpue_cHL");
  #ifndef NO_AD_INITIALIZE
  w_cpue_cHL.initialize();
  #endif
  w_cpue_sTV.allocate("w_cpue_sTV");
  #ifndef NO_AD_INITIALIZE
  w_cpue_sTV.initialize();
  #endif
  w_cpue_rHB.allocate("w_cpue_rHB");
  #ifndef NO_AD_INITIALIZE
  w_cpue_rHB.initialize();
  #endif
  w_cpue_rGN.allocate("w_cpue_rGN");
  #ifndef NO_AD_INITIALIZE
  w_cpue_rGN.initialize();
  #endif
  w_lenc_cHL.allocate("w_lenc_cHL");
  #ifndef NO_AD_INITIALIZE
  w_lenc_cHL.initialize();
  #endif
  w_lenc_sTV.allocate("w_lenc_sTV");
  #ifndef NO_AD_INITIALIZE
  w_lenc_sTV.initialize();
  #endif
  w_lenc_rHB.allocate("w_lenc_rHB");
  #ifndef NO_AD_INITIALIZE
  w_lenc_rHB.initialize();
  #endif
  w_lenc_rHB_D.allocate("w_lenc_rHB_D");
  #ifndef NO_AD_INITIALIZE
  w_lenc_rHB_D.initialize();
  #endif
  w_lenc_rGN.allocate("w_lenc_rGN");
  #ifndef NO_AD_INITIALIZE
  w_lenc_rGN.initialize();
  #endif
  w_agec_cHL.allocate("w_agec_cHL");
  #ifndef NO_AD_INITIALIZE
  w_agec_cHL.initialize();
  #endif
  w_agec_sTV.allocate("w_agec_sTV");
  #ifndef NO_AD_INITIALIZE
  w_agec_sTV.initialize();
  #endif
  w_agec_rHB.allocate("w_agec_rHB");
  #ifndef NO_AD_INITIALIZE
  w_agec_rHB.initialize();
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
  f_cHL_L.allocate("f_cHL_L");
  #ifndef NO_AD_INITIALIZE
  f_cHL_L.initialize();
  #endif
  f_rHB_L.allocate("f_rHB_L");
  #ifndef NO_AD_INITIALIZE
  f_rHB_L.initialize();
  #endif
  f_rGN_L.allocate("f_rGN_L");
  #ifndef NO_AD_INITIALIZE
  f_rGN_L.initialize();
  #endif
  f_rHB_D.allocate("f_rHB_D");
  #ifndef NO_AD_INITIALIZE
  f_rHB_D.initialize();
  #endif
  f_rGN_D.allocate("f_rGN_D");
  #ifndef NO_AD_INITIALIZE
  f_rGN_D.initialize();
  #endif
  f_cHL_cpue.allocate("f_cHL_cpue");
  #ifndef NO_AD_INITIALIZE
  f_cHL_cpue.initialize();
  #endif
  f_sTV_cpue.allocate("f_sTV_cpue");
  #ifndef NO_AD_INITIALIZE
  f_sTV_cpue.initialize();
  #endif
  f_rHB_cpue.allocate("f_rHB_cpue");
  #ifndef NO_AD_INITIALIZE
  f_rHB_cpue.initialize();
  #endif
  f_rGN_cpue.allocate("f_rGN_cpue");
  #ifndef NO_AD_INITIALIZE
  f_rGN_cpue.initialize();
  #endif
  f_cHL_lenc.allocate("f_cHL_lenc");
  #ifndef NO_AD_INITIALIZE
  f_cHL_lenc.initialize();
  #endif
  f_sTV_lenc.allocate("f_sTV_lenc");
  #ifndef NO_AD_INITIALIZE
  f_sTV_lenc.initialize();
  #endif
  f_rHB_lenc.allocate("f_rHB_lenc");
  #ifndef NO_AD_INITIALIZE
  f_rHB_lenc.initialize();
  #endif
  f_rHB_D_lenc.allocate("f_rHB_D_lenc");
  #ifndef NO_AD_INITIALIZE
  f_rHB_D_lenc.initialize();
  #endif
  f_rGN_lenc.allocate("f_rGN_lenc");
  #ifndef NO_AD_INITIALIZE
  f_rGN_lenc.initialize();
  #endif
  f_cHL_agec.allocate("f_cHL_agec");
  #ifndef NO_AD_INITIALIZE
  f_cHL_agec.initialize();
  #endif
  f_sTV_agec.allocate("f_sTV_agec");
  #ifndef NO_AD_INITIALIZE
  f_sTV_agec.initialize();
  #endif
  f_rHB_agec.allocate("f_rHB_agec");
  #ifndef NO_AD_INITIALIZE
  f_rHB_agec.initialize();
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
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{1000, 2000, 3000, 10000;}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-2, 1e-2, 1e-3, 1e-4;}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  Dmort_rHB=set_Dmort_rHB;
  obs_rHB_D=Dmort_rHB*obs_released_rHB;
  Dmort_rGN=set_Dmort_rGN;
  obs_rGN_D=Dmort_rGN*obs_released_rGN;
  Linf=set_Linf(1);
  K=set_K(1);
  t0=set_t0(1);
  len_cv_val=set_len_cv(1);
  // FISHERY LANDINGS Growth
  Linf_L=set_Linf_L(1);
  K_L=set_K_L(1);
  t0_L=set_t0_L(1);
  len_cv_val_L=set_len_cv_L(1);
  // FISHERY INDEPENDENT (SERFS chevron trap/video) Growth
  Linf_sTV=set_Linf_sTV(1);
  K_sTV=set_K_sTV(1);
  t0_sTV=set_t0_sTV(1);
  len_cv_val_sTV=set_len_cv_sTV(1);
  M=set_M; 
  M_constant=set_M_constant(1);
  smsy2msst=1.0-M_constant;
  smsy2msst75=0.75;  
  
  log_R0=set_log_R0(1);
  steep=set_steep(1);
  R_autocorr=set_R_autocorr(1);
  rec_sigma=set_rec_sigma(1);
  
  log_q_cpue_cHL=set_log_q_cpue_cHL(1);
  log_q_cpue_sTV=set_log_q_cpue_sTV(1);  
  log_q_cpue_rHB=set_log_q_cpue_rHB(1);
  log_q_cpue_rGN=set_log_q_cpue_rGN(1);
  
  q_rate=set_q_rate;
  q_rate_fcn_cHL=1.0;   
  q_rate_fcn_sTV=1.0;   
  q_rate_fcn_rHB=1.0;   
  q_rate_fcn_rGN=1.0; 
  q_DD_beta=set_q_DD_beta;
  q_DD_fcn=1.0;
  q_RW_log_dev_cHL.initialize(); 
  q_RW_log_dev_sTV.initialize();  
  q_RW_log_dev_rHB.initialize(); 
  q_RW_log_dev_rGN.initialize(); 
      
  if (set_q_rate_phase<0 & q_rate!=0.0)
  {
    for (iyear=styr_cpue_cHL; iyear<=endyr_cpue_cHL; iyear++)
      {   if (iyear>styr_cpue_cHL & iyear <=2003) 
          {//q_rate_fcn_cHL(iyear)=(1.0+q_rate)*q_rate_fcn_cHL(iyear-1); //compound
             q_rate_fcn_cHL(iyear)=(1.0+(iyear-styr_cpue_cHL)*q_rate)*q_rate_fcn_cHL(styr_cpue_cHL);  //linear
          }
          if (iyear>2003) {q_rate_fcn_cHL(iyear)=q_rate_fcn_cHL(iyear-1);} 
      }   
    for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++)
      {   if (iyear>styr_cpue_rHB & iyear <=2003) 
          {//q_rate_fcn_rHB(iyear)=(1.0+q_rate)*q_rate_fcn_rHB(iyear-1); //compound
             q_rate_fcn_rHB(iyear)=(1.0+(iyear-styr_cpue_rHB)*q_rate)*q_rate_fcn_rHB(styr_cpue_rHB);  //linear
          }
          if (iyear>2003) {q_rate_fcn_rHB(iyear)=q_rate_fcn_rHB(iyear-1);} 
      }   
    for (iyear=styr_cpue_rGN; iyear<=endyr_cpue_rGN; iyear++)
      {   if (iyear>styr_cpue_rGN & iyear <=2003) 
          {//q_rate_fcn_rGN(iyear)=(1.0+q_rate)*q_rate_fcn_rGN(iyear-1); //compound
             q_rate_fcn_rGN(iyear)=(1.0+(iyear-styr_cpue_rGN)*q_rate)*q_rate_fcn_rGN(styr_cpue_rGN);  //linear
          }
          if (iyear>2003) {q_rate_fcn_rGN(iyear)=q_rate_fcn_rGN(iyear-1);} 
      }   
  } //end q_rate conditional      
  w_L=set_w_L;
  w_D=set_w_D;
  
  w_cpue_cHL=set_w_cpue_cHL;
  w_cpue_sTV=set_w_cpue_sTV;
  w_cpue_rHB=set_w_cpue_rHB;
  w_cpue_rGN=set_w_cpue_rGN;
  w_lenc_cHL=set_w_lenc_cHL;
  w_lenc_sTV=set_w_lenc_sTV;
  w_lenc_rHB=set_w_lenc_rHB;
  w_lenc_rHB_D=set_w_lenc_rHB_D;  
  w_lenc_rGN=set_w_lenc_rGN;   
  
  w_agec_cHL=set_w_agec_cHL;
  w_agec_sTV=set_w_agec_sTV;
  w_agec_rHB=set_w_agec_rHB;
  //w_ac_rGN=set_w_ac_rGN;
  w_Nage_init=set_w_Nage_init;
  w_rec=set_w_rec;
  w_rec_early=set_w_rec_early;
  w_rec_end=set_w_rec_end;
  w_fullF=set_w_fullF;
  w_Ftune=set_w_Ftune;
  F_init=set_F_init(1);
  log_avg_F_L_cHL=set_log_avg_F_L_cHL(1);
  log_avg_F_L_rHB=set_log_avg_F_L_rHB(1); 
  log_avg_F_L_rGN=set_log_avg_F_L_rGN(1); 
  log_avg_F_D_rHB=set_log_avg_F_D_rHB(1); 
  log_avg_F_D_rGN=set_log_avg_F_D_rGN(1); 
    
  log_dev_F_L_cHL=set_log_dev_vals_F_L_cHL;
  log_dev_F_L_rHB=set_log_dev_vals_F_L_rHB;
  log_dev_F_L_rGN=set_log_dev_vals_F_L_rGN;
  log_dev_F_D_rHB=set_log_dev_vals_F_D_rHB;
  log_dev_F_D_rGN=set_log_dev_vals_F_D_rGN;
 
  selpar_L50_cHL1=set_selpar_L50_cHL1(1);
  selpar_slope_cHL1=set_selpar_slope_cHL1(1);
  selpar_L50_sTV=set_selpar_L50_sTV(1);
  selpar_slope_sTV=set_selpar_slope_sTV(1);
  selpar_L51_rHB1=set_selpar_L51_rHB1(1);
  selpar_slope1_rHB1=set_selpar_slope1_rHB1(1);
  selpar_L52_rHB1=set_selpar_L52_rHB1(1);
  selpar_slope2_rHB1=set_selpar_slope2_rHB1(1);
  selpar_L50_rHB2=set_selpar_L50_rHB2(1);
  selpar_slope_rHB2=set_selpar_slope_rHB2(1);
  selpar_afull_rHB2=set_selpar_afull_rHB2(1);
  selpar_sigma_rHB2=set_selpar_sigma_rHB2(1);
  
  selpar_L50_rHB_D=set_selpar_L50_rHB_D(1);
  selpar_slope_rHB_D=set_selpar_slope_rHB_D(1);
  selpar_afull_rHB_D=set_selpar_afull_rHB_D(1);
  selpar_sigma_rHB_D=set_selpar_sigma_rHB_D(1);
  selpar_L51_rGN=set_selpar_L51_rGN(1);    
  selpar_slope1_rGN=set_selpar_slope1_rGN(1);  
  selpar_L52_rGN=set_selpar_L52_rGN(1);  
  selpar_slope2_rGN=set_selpar_slope2_rGN(1);  
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
 //lbins=lenbins; //NOT NEEDED
 
      nsamp_cHL_lenc_allyr=missing;
      nsamp_sTV_lenc_allyr=missing;
      nsamp_rHB_lenc_allyr=missing;
      nsamp_rHB_D_lenc_allyr=missing;  
      nsamp_rGN_lenc_allyr=missing;  
      nsamp_cHL_agec_allyr=missing;
      nsamp_sTV_agec_allyr=missing;
      nsamp_rHB_agec_allyr=missing;
      //nsamp_rGN_agec_allyr=missing;
      nfish_cHL_lenc_allyr=missing;
      nfish_sTV_lenc_allyr=missing;
      nfish_rHB_lenc_allyr=missing;
      nfish_rHB_D_lenc_allyr=missing; 
      nfish_rGN_lenc_allyr=missing;   
      nfish_cHL_agec_allyr=missing;
      nfish_sTV_agec_allyr=missing;
      nfish_rHB_agec_allyr=missing;
      //nfish_rGN_agec_allyr=missing;
      for (iyear=1; iyear<=nyr_lenc_cHL; iyear++)
         {if (nsamp_lenc_cHL(iyear)>=minSS_lenc_cHL)
           {nsamp_cHL_lenc_allyr(yrs_lenc_cHL(iyear))=nsamp_lenc_cHL(iyear);
            nfish_cHL_lenc_allyr(yrs_lenc_cHL(iyear))=nfish_lenc_cHL(iyear);}}
      for (iyear=1; iyear<=nyr_lenc_sTV; iyear++)
         {if (nsamp_lenc_sTV(iyear)>=minSS_lenc_sTV)
           {nsamp_sTV_lenc_allyr(yrs_lenc_sTV(iyear))=nsamp_lenc_sTV(iyear);
            nfish_sTV_lenc_allyr(yrs_lenc_sTV(iyear))=nfish_lenc_sTV(iyear);}}
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
      for (iyear=1; iyear<=nyr_agec_cHL; iyear++)
         {if (nsamp_agec_cHL(iyear)>=minSS_agec_cHL)
           {nsamp_cHL_agec_allyr(yrs_agec_cHL(iyear))=nsamp_agec_cHL(iyear);
            nfish_cHL_agec_allyr(yrs_agec_cHL(iyear))=nfish_agec_cHL(iyear);}}
      for (iyear=1; iyear<=nyr_agec_sTV; iyear++)
         {if (nsamp_agec_sTV(iyear)>=minSS_agec_sTV)
           {nsamp_sTV_agec_allyr(yrs_agec_sTV(iyear))=nsamp_agec_sTV(iyear);
            nfish_sTV_agec_allyr(yrs_agec_sTV(iyear))=nfish_agec_sTV(iyear);}}
      for (iyear=1; iyear<=nyr_agec_rHB; iyear++)
         {if (nsamp_agec_rHB(iyear)>=minSS_agec_rHB)
           {nsamp_rHB_agec_allyr(yrs_agec_rHB(iyear))=nsamp_agec_rHB(iyear);
            nfish_rHB_agec_allyr(yrs_agec_rHB(iyear))=nfish_agec_rHB(iyear);}} 
             
  F_msy(1)=0.0;  
  for (ff=2;ff<=n_iter_msy;ff++) {F_msy(ff)=F_msy(ff-1)+iter_inc_msy;}
  F_spr(1)=0.0;  
  for (ff=2;ff<=n_iter_spr;ff++) {F_spr(ff)=F_spr(ff-1)+iter_inc_spr;}
  F_cHL.initialize(); L_cHL_num.initialize();
  F_rHB.initialize(); L_rHB_num.initialize();
  F_rGN.initialize(); L_rGN_num.initialize();
  F_rHB_D.initialize(); D_rHB_num.initialize();
  F_rGN_D.initialize(); D_rGN_num.initialize();
  F_cHL_out.initialize();
  F_rHB_out.initialize();
  F_rGN_out.initialize();
  F_rHB_D_out.initialize();
  F_rGN_D_out.initialize();
      
  sel_cHL.initialize();
  sel_sTV.initialize();
  sel_rHB.initialize();
  sel_rHB_D.initialize();
  sel_rHB.initialize();  
  
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
 //cout << "got length, weight, fecundity transitions" << endl;
 get_reprod();
 //cout << "got repro stuff" << endl;
 get_length_at_age_dist(); 
 //cout << "got predicted length at age distribution" << endl;
 get_weight_at_age_landings();
 //cout << "got weight at age of landings" << endl; 
 get_spr_F0();
 //cout << "got F0 spr" << endl;
 get_selectivity(); 
 //cout << "got selectivity" << endl;
 get_mortality(); 
 //cout << "got mortalities" << endl;
 get_bias_corr(); 
 //cout << "got recruitment bias correction" << endl;
 get_numbers_at_age(); 
 //cout << "got numbers at age" << endl;
 get_landings_numbers();
 //cout << "got landings in numbers" << endl;
 get_landings_wgt();
 //cout << "got landings in wgt" << endl;
 get_dead_discards(); 
 get_catchability_fcns(); 
 //cout << "got catchability_fcns" << endl;
 get_indices();
 //cout << "got indices" << endl;
 get_length_comps();
 //cout << "got length comps" << endl;
 get_age_comps();
 //cout << "got age comps" << endl;
 evaluate_objective_function();
 //cout << "objective function calculations complete" << endl;
}

void model_parameters::get_length_weight_at_age(void)
{
  //compute mean length (mm TL) and weight (whole) at age
    meanlen_FL=Linf*(1.0-mfexp(-K*(agebins-t0+0.5)));     //total length in mm
    wgt_kg=wgtpar_a*pow(meanlen_FL,wgtpar_b);             //whole wgt in kg 
    wgt_g=wgt_kg/g2kg;                                    //convert wgt in kg to weight in g    
    wgt_mt=wgt_g*g2mt;                                    //convert weight in g to weight in mt
    wgt_klb=mt2klb*wgt_mt;                                //1000 lb of whole wgt
    wgt_lb=mt2lb*wgt_mt;                                  //lb of whole wgt
    fecundity=(fecpar_a+fecpar_b*meanlen_FL)/fecpar_scale; //annual egg production of a mature female at age in units of fecpar_scale
  //THESE CALCULATIONS ASSUME THE FISHERY LANDINGS GROWTH CURVE
    meanlen_FL_L=Linf_L*(1.0-mfexp(-K_L*(agebins-t0_L+0.5))); //fork length in mm
    wgt_kg_L=wgtpar_a*pow(meanlen_FL_L,wgtpar_b);             //wgt in kg KC changed from wgt in metric tons (wgt_mt) because L-W relationship is mm FL and kg whole weight
    wgt_g_L=wgt_kg_L/g2kg;                                    //KC convert wgt in kg from L-W relationship to weight in g    
    wgt_mt_L=wgt_g_L*g2mt;                                    //KC convert weight in g to weight in mt
    wgt_klb_L=mt2klb*wgt_mt_L;                                //1000 lb of whole wgt
    wgt_lb_L=mt2lb*wgt_mt_L;                                  //1000 lb of whole wgt
  //THESE CALCULATIONS ASSUME THE FISHERY (SERFS chevron trap/video) GROWTH CURVE
    meanlen_FL_sTV=Linf_sTV*(1.0-mfexp(-K_sTV*(agebins-t0_sTV+0.5)));    //fork length in mm
}

void model_parameters::get_reprod(void)
{
   //reprod is product of stuff going into reproductive capacity calcs
   reprod=elem_prod(elem_prod(elem_prod(prop_f,maturity_f),fecundity),fecpar_batches); 
   reprod2=elem_prod(elem_prod(prop_f,maturity_f),wgt_mt);
}

void model_parameters::get_length_at_age_dist(void)
{
  //compute matrix of length at age, based on the normal distribution
 dvar_vector length(1,nages);
  dvariable cvlen;
  if (use_landings_growth==1.0)														//LC added to allow for use of Landings Growth Curve
  	{
  	length=meanlen_FL_L;
  	cvlen=len_cv_val_L;
        len_sd=len_sd_L;
  	}
   if (use_landings_growth==0.0) 	
  	{
  	length=meanlen_FL;
 	cvlen=len_cv_val;
 //     sdlen=len_sd_val;
        len_sd=len_sd;
	}
  for (iage=1;iage<=nages;iage++)
  {
    len_cv(iage)=cvlen; 
    len_sd(iage)=length(iage)*len_cv(iage); 
    zscore_lzero=(0.0-length(iage))/len_sd(iage);  
    cprob_lzero=cumd_norm(zscore_lzero); 
    //first length bin
    zscore_len=((lenbins(1)+0.5*lenbins_width)-length(iage)) / len_sd(iage);  
    cprob_lenvec(1)=cumd_norm(zscore_len);          //includes any probability mass below zero
    lenprob(iage,1)=cprob_lenvec(1)-cprob_lzero;    //removes any probability mass below zero
    //most other length bins     
    for (ilen=2;ilen<nlenbins;ilen++)
      {
        zscore_len=((lenbins(ilen)+0.5*lenbins_width)-length(iage)) / len_sd(iage);
        cprob_lenvec(ilen)=cumd_norm(zscore_len);
        lenprob(iage,ilen)=cprob_lenvec(ilen)-cprob_lenvec(ilen-1);
      }
    //last length bin is a plus group
    zscore_len=((lenbins(nlenbins)-0.5*lenbins_width)-length(iage)) / len_sd(iage);
    lenprob(iage,nlenbins)=1.0-cumd_norm(zscore_len);
    lenprob(iage)=lenprob(iage)/(1.0-cprob_lzero);  //renormalize to account for any prob mass below size=0
  }
  //fleet and survey specific length probs, all assumed here to equal the popn
  lenprob_cHL=lenprob;
  lenprob_rHB=lenprob;
  lenprob_rGN=lenprob;
  lenprob_rHB_D=lenprob; 
    ///////////// THIS SECTION FOR SERFS chevron trap/video (sTV)  ////////////////////////// 
  lenprob2.initialize();
  length=meanlen_FL_sTV;
  cvlen=len_cv_val_sTV;
  len_sd=len_sd_sTV;
  for (iage=1;iage<=nages;iage++)
    {
    len_cv(iage)=cvlen;															
    len_sd(iage)=length(iage)*len_cv(iage);											
      zscore_lzero2=(0.0-length(iage))/len_sd(iage);
      cprob_lzero2=cumd_norm(zscore_lzero2);
      zscore_len2=((lenbins(1)+0.5*lenbins_width)-length(iage)) / len_sd(iage);  					
      cprob_lenvec2(1)=cumd_norm(zscore_len2);         
      lenprob2(iage,1)=cprob_lenvec2(1)-cprob_lzero2;   
      //most other length bins     
      for (ilen=2;ilen<nlenbins;ilen++)  
        {
          zscore_len2=((lenbins(ilen)+0.5*lenbins_width)-length(iage)) / len_sd(iage);  				
          cprob_lenvec2(ilen)=cumd_norm(zscore_len2);
          lenprob2(iage,ilen)=cprob_lenvec2(ilen)-cprob_lenvec2(ilen-1);
        }
      //last length bin is a plus group
      zscore_len2=((lenbins(nlenbins)-0.5*lenbins_width)-length(iage)) / len_sd(iage); 
          lenprob2(iage,nlenbins)=1.0-cumd_norm(zscore_len2); 
      lenprob2(iage)=lenprob2(iage)/(1.0-cprob_lzero2);  //renormalize to account for any prob mass below size=0
    }
    //fleet and survey specific length probs
    lenprob_sTV=lenprob2;
}

void model_parameters::get_weight_at_age_landings(void)
{
  //KC5--modified based on Lew code
   if(use_landings_growth==1.0){ 
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    len_cHL_mm(iyear)=meanlen_FL_L;
    wholewgt_cHL_klb(iyear)=wgt_klb_L; //wholeweight used to match index 
    len_rHB_mm(iyear)=meanlen_FL_L;
    wholewgt_rHB_klb(iyear)=wgt_klb_L; //KC6  change _out to _Lc
    len_rGN_mm(iyear)=meanlen_FL_L;
    wholewgt_rGN_klb(iyear)=wgt_klb_L; //KC6   change _out to _L
    len_rHB_D_mm(iyear)=meanlen_FL_L;
    wholewgt_rHB_D_klb(iyear)=wgt_klb_L;  //KC6; change _out to _L
    len_rGN_D_mm(iyear)=meanlen_FL_L;
    wholewgt_rGN_D_klb(iyear)=wgt_klb_L;  //KC6  change _out to _L
  } 
  }
   if(use_landings_growth==0.0){ 
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    len_cHL_mm(iyear)=meanlen_FL; 
    wholewgt_cHL_klb(iyear)=wgt_klb; //wholeweight used to match index  
    len_rHB_mm(iyear)=meanlen_FL;
    wholewgt_rHB_klb(iyear)=wgt_klb;  //KC6  KC6-got rid _out
    len_rGN_mm(iyear)=meanlen_FL;
    wholewgt_rGN_klb(iyear)=wgt_klb;  //KC6  KC6-got rid _out
    len_rHB_D_mm(iyear)=meanlen_FL;
    wholewgt_rHB_D_klb(iyear)=wgt_klb;  //KC6  KC6-got rid _out
    len_rGN_D_mm(iyear)=meanlen_FL;
    wholewgt_rGN_D_klb(iyear)=wgt_klb;  //KC6  KC6-got rid _out
  } 
  } 
}

void model_parameters::get_spr_F0(void)
{
  //at mdyr, apply half this yr's mortality, half next yr's
  N_spr_F0(1)=1.0*mfexp(-1.0*M(1)*spawn_time_frac); //at peak spawning time
  N_bpr_F0(1)=1.0;      //at start of year
  for (iage=2; iage<=nages; iage++)
  {
    //N_spr_F0(iage)=N_spr_F0(iage-1)*mfexp(-1.0*(M(iage-1));
    N_spr_F0(iage)=N_spr_F0(iage-1)*mfexp(-1.0*(M(iage-1)*(1.0-spawn_time_frac) + M(iage)*spawn_time_frac)); 
    N_bpr_F0(iage)=N_bpr_F0(iage-1)*mfexp(-1.0*(M(iage-1)));    
  }
  N_spr_F0(nages)=N_spr_F0(nages)/(1.0-mfexp(-1.0*M(nages))); //plus group (sum of geometric series)
  N_bpr_F0(nages)=N_bpr_F0(nages)/(1.0-mfexp(-1.0*M(nages)));
  spr_F0=sum(elem_prod(N_spr_F0,reprod));  
  bpr_F0=sum(elem_prod(N_bpr_F0,wgt_mt));    
}

void model_parameters::get_selectivity(void)
{
  //BLOCK 1 for selex. 
  for (iyear=styr; iyear<=endyr_selex_phase1; iyear++)
  {     
     sel_cHL(iyear)=logistic(agebins, selpar_L50_cHL1, selpar_slope_cHL1);
     sel_sTV(iyear)=logistic(agebins, selpar_L50_sTV, selpar_slope_sTV);
     //sel_sTV(iyear)=logistic_exponential(agebins, selpar_L50_sTV, selpar_slope_sTV, selpar_sigma_sTV, selpar_afull_sTV);   
     //sel_vid(iyear)=logistic_exponential(agebins, selpar_L50_vid, selpar_slope_vid, selpar_sigma_vid, selpar_afull_vid);      
     //sel_vid(iyear)=logistic(agebins, selpar_L50_vid, selpar_slope_vid);
     //sel_vid(iyear)=sel_sTV(iyear);
     sel_rHB(iyear)=logistic_double(agebins, selpar_L51_rHB1, selpar_slope1_rHB1, selpar_L52_rHB1, selpar_slope2_rHB1);
    //sel_rHB(iyear)=logistic(agebins, selpar_L50_rHB1, selpar_slope_rHB1);
    //sel_rHB(iyear)=logistic_exponential(agebins, selpar_L50_rHB1, selpar_slope_rHB1, selpar_sigma_rHB1, selpar_afull_rHB1);
     sel_rHB_D(iyear)=logistic_exponential(agebins, selpar_L50_rHB_D, selpar_slope_rHB_D, selpar_sigma_rHB_D, selpar_afull_rHB_D); 
     sel_rGN(iyear)=logistic_double(agebins, selpar_L51_rGN, selpar_slope1_rGN, selpar_L52_rGN, selpar_slope2_rGN);
    //sel_rGN(iyear)=logistic(agebins, selpar_L50_rGN, selpar_slope_rGN);
     sel_rGN_D(iyear)=sel_rHB_D(iyear); 
  }
}

void model_parameters::get_mortality(void)
{
  Fsum.initialize();
  Fapex.initialize();
  F.initialize();
  //initialization F is avg from first 3 yrs of observed landings
  log_F_dev_init_cHL=sum(log_dev_F_L_cHL(styr_L_cHL,(styr_L_cHL+2)))/3.0;         
  log_F_dev_init_rHB=sum(log_dev_F_L_rHB(styr_L_rHB,(styr_L_rHB+2)))/3.0;         
  log_F_dev_init_rGN=sum(log_dev_F_L_rGN(styr_L_rGN,(styr_L_rGN+2)))/3.0;         
  for (iyear=styr; iyear<=endyr; iyear++) 
  {
    if(iyear>=styr_L_cHL & iyear<=endyr_L_cHL)
    {  F_cHL_out(iyear)=mfexp(log_avg_F_L_cHL+log_dev_F_L_cHL(iyear)); //}    
    // if (iyear<styr_L_cHL){F_cHL_out(iyear)=mfexp(log_avg_F_L_cHL+log_F_dev_init_cHL);}        
       F_cHL(iyear)=sel_cHL(iyear)*F_cHL_out(iyear);
       Fsum(iyear)+=F_cHL_out(iyear);
    }
    if(iyear>=styr_L_rHB & iyear<=endyr_L_rHB)
    {  F_rHB_out(iyear)=mfexp(log_avg_F_L_rHB+log_dev_F_L_rHB(iyear)); //}    
    // if (iyear<styr_L_rHB){F_rHB_out(iyear)=mfexp(log_avg_F_L_rHB+log_F_dev_init_rHB);}        
       F_rHB(iyear)=sel_rHB(iyear)*F_rHB_out(iyear);
       Fsum(iyear)+=F_rHB_out(iyear);
    }
    if(iyear>=styr_L_rGN & iyear<=endyr_L_rGN) 
    {  F_rGN_out(iyear)=mfexp(log_avg_F_L_rGN+log_dev_F_L_rGN(iyear)); //}    
    // if (iyear<styr_L_rGN){F_rGN_out(iyear)=mfexp(log_avg_F_L_rGN+log_F_dev_init_rGN);}        
       F_rGN(iyear)=sel_rGN(iyear)*F_rGN_out(iyear); //general rec shares headboat selex 
       Fsum(iyear)+=F_rGN_out(iyear);
    }
    if(iyear>=styr_D_rHB & iyear<=endyr_D_rHB)
    {  F_rHB_D_out(iyear)=mfexp(log_avg_F_D_rHB+log_dev_F_D_rHB(iyear)); //}    
    // if (iyear<styr_D_rHB){F_rHB_D_out(iyear)=mfexp(log_avg_F_D_rHB+log_F_dev_init_rHB_D);}        
       F_rHB_D(iyear)=sel_rHB_D(iyear)*F_rHB_D_out(iyear);
       Fsum(iyear)+=F_rHB_D_out(iyear);
    }
    if(iyear>=styr_D_rGN & iyear<=endyr_D_rGN)
    {  F_rGN_D_out(iyear)=mfexp(log_avg_F_D_rGN+log_dev_F_D_rGN(iyear)); //}    
    // if (iyear<styr_D_rGN){F_rGN_D_out(iyear)=mfexp(log_avg_F_D_rGN+log_F_dev_init_rGN_D);}        
       F_rGN_D(iyear)=sel_rHB_D(iyear)*F_rGN_D_out(iyear); //general rec discards shares headboat discards selex
       Fsum(iyear)+=F_rGN_D_out(iyear);
    }
    //Total F at age
    F(iyear)=F_cHL(iyear);  //first in additive series (NO +=)
    F(iyear)+=F_rHB(iyear);
    F(iyear)+=F_rGN(iyear);
    F(iyear)+=F_rHB_D(iyear);
    F(iyear)+=F_rGN_D(iyear);
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
  //R_virgin=(R0/((5.0*steep-1.0)*spr_F0))*
  //               (BiasCor*4.0*steep*spr_F0-spr_F0*(1.0-steep));
  //R_virgin=R0/(spr_F0/spr_F0)*BiasCor*(1.0+log(spr_F0/spr_F0)/steep); //Ricker
  R_virgin=SR_eq_func(R0, steep, spr_F0, spr_F0, BiasCor, SR_switch);
  B0=bpr_F0*R_virgin;   
  B0_q_DD=R_virgin*sum(elem_prod(N_bpr_F0(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages))); 
  F_init_denom=mfexp(log_avg_F_L_cHL+log_F_dev_init_cHL)+
               mfexp(log_avg_F_L_rHB+log_F_dev_init_rHB)+
               mfexp(log_avg_F_L_rGN+log_F_dev_init_rGN)+
               mfexp(log_avg_F_D_rHB+log_F_dev_init_rHB_D)+
               mfexp(log_avg_F_D_rGN+log_F_dev_init_rGN_D);
  F_init_cHL_prop=mfexp(log_avg_F_L_cHL+log_F_dev_init_cHL)/F_init_denom;
  F_init_rHB_prop=mfexp(log_avg_F_L_rHB+log_F_dev_init_rHB)/F_init_denom;
  F_init_rGN_prop=mfexp(log_avg_F_L_rGN+log_F_dev_init_rGN)/F_init_denom;
  F_init_rHB_D_prop=mfexp(log_avg_F_D_rHB+log_F_dev_init_rHB_D)/F_init_denom;
  F_init_rGN_D_prop=mfexp(log_avg_F_D_rGN+log_F_dev_init_rGN_D)/F_init_denom;
  F_initial=sel_cHL(styr)*F_init*F_init_cHL_prop+
            sel_rHB(styr)*F_init*F_init_rHB_prop+
            sel_rGN(styr)*F_init*F_init_rGN_prop+
            sel_rHB_D(styr)*F_init*F_init_rHB_prop+
            sel_rGN_D(styr)*F_init*F_init_rGN_prop;
  Z_initial=M+F_initial;
  N_spr_initial(1)=1.0*mfexp(-1.0*Z_initial(1)*spawn_time_frac); //at peak spawning time;
  for (iage=2; iage<=nages; iage++)
    {
      N_spr_initial(iage)=N_spr_initial(iage-1)*
                   mfexp(-1.0*(Z_initial(iage-1)*(1.0-spawn_time_frac) + Z_initial(iage)*spawn_time_frac)); 
    }
  N_spr_initial(nages)=N_spr_initial(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
  spr_initial=sum(elem_prod(N_spr_initial,reprod));  
  //if (styr==styr_rec_dev) {R1=(R0/((5.0*steep-1.0)*spr_initial))*
  //               (4.0*steep*spr_initial-spr_F0*(1.0-steep));} //without bias correction (deviation added later)
  //else {R1=(R0/((5.0*steep-1.0)*spr_initial))*
  //               (BiasCor*4.0*steep*spr_initial-spr_F0*(1.0-steep));} //with bias correction                 
  if (styr==styr_rec_dev) {R1=SR_eq_func(R0, steep, spr_F0, spr_initial, 1.0, SR_switch);} //without bias correction (deviation added later)
  else {R1=SR_eq_func(R0, steep, spr_F0, spr_initial, BiasCor, SR_switch);} //with bias correction
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
  SSB(styr)=sum(elem_prod(N_spawn(styr),reprod));  //KC4
  MatFemB(styr)=sum(elem_prod(N_spawn(styr),reprod2));
  B_q_DD(styr)=sum(elem_prod(N(styr)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
  for (iyear=styr; iyear<endyr; iyear++)
  {
    if(iyear<(styr_rec_dev-1)||iyear>(endyr_rec_dev-1)) //recruitment follows S-R curve (with bias correction) exactly
    {
        N(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch);
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); //plus group
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod)); //KC
        MatFemB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod2));  //KC  
        B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));       
    }
    else   //recruitment follows S-R curve with lognormal deviation
    {
        N(iyear+1,1)=SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch)*mfexp(log_dev_rec(iyear+1));
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); //plus group
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));  //KC4
        MatFemB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod2));
        B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
    }
  }
  //last year (projection) has no recruitment variability
  N(endyr+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(endyr),SR_switch);
  N(endyr+1)(2,nages)=++elem_prod(N(endyr)(1,nages-1),(mfexp(-1.*Z(endyr)(1,nages-1))));
  N(endyr+1,nages)+=N(endyr,nages)*mfexp(-1.*Z(endyr,nages)); //plus group
}

void model_parameters::get_landings_numbers(void)
{
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      L_cHL_num(iyear,iage)=N(iyear,iage)*F_cHL(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_rHB_num(iyear,iage)=N(iyear,iage)*F_rHB(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_rGN_num(iyear,iage)=N(iyear,iage)*F_rGN(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);        
    }          
    pred_cHL_L_knum(iyear)=sum(L_cHL_num(iyear))/1000.0;    
    pred_rHB_L_knum(iyear)=sum(L_rHB_num(iyear))/1000.0;
    pred_rGN_L_knum(iyear)=sum(L_rGN_num(iyear))/1000.0;
  }
}

void model_parameters::get_landings_wgt(void)
{
  for (iyear=styr; iyear<=endyr; iyear++)
  {    
    L_cHL_klb(iyear)=elem_prod(L_cHL_num(iyear),wholewgt_cHL_klb(iyear));     //in 1000 lb whole weight  
    L_rHB_klb(iyear)=elem_prod(L_rHB_num(iyear),wholewgt_rHB_klb(iyear));     //in 1000 lb whole weight 
    L_rGN_klb(iyear)=elem_prod(L_rGN_num(iyear),wholewgt_rGN_klb(iyear));     //in 1000 lb whole weight 
    pred_cHL_L_klb(iyear)=sum(L_cHL_klb(iyear));
    pred_rHB_L_klb(iyear)=sum(L_rHB_klb(iyear));
    pred_rGN_L_klb(iyear)=sum(L_rGN_klb(iyear));    
  }
}

void model_parameters::get_dead_discards(void)
{
  //dead discards at age (number fish) 
  for (iyear=styr_D_rHB; iyear<=endyr_D_rHB; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      D_rHB_num(iyear,iage)=N(iyear,iage)*F_rHB_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
    }
    pred_rHB_D_knum(iyear)=sum(D_rHB_num(iyear))/1000.0;            //pred annual dead discards in 1000s (for matching data)
      pred_rHB_D_klb(iyear)=sum(elem_prod(D_rHB_num(iyear),wholewgt_rHB_D_klb(iyear))); //annual dead discards in 1000 lb whole (for output only) 
  }
  for (iyear=styr_D_rGN; iyear<=endyr_D_rGN; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      D_rGN_num(iyear,iage)=N(iyear,iage)*F_rGN_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
    }
    pred_rGN_D_knum(iyear)=sum(D_rGN_num(iyear))/1000.0;            //pred annual dead discards in 1000s (for matching data)
    pred_rGN_D_klb(iyear)=sum(elem_prod(D_rGN_num(iyear),wholewgt_rGN_D_klb(iyear))); //annual dead discards in 1000 lb whole (for output only) 
  }
}

void model_parameters::get_catchability_fcns(void)
{
 //Get rate increase if estimated, otherwise fixed above
  if (set_q_rate_phase>0.0)
  {
      for (iyear=styr_cpue_cHL; iyear<=endyr_cpue_cHL; iyear++)
      {   if (iyear>styr_cpue_cHL & iyear <=2003) 
          {//q_rate_fcn_cHL(iyear)=(1.0+q_rate)*q_rate_fcn_cHL(iyear-1); //compound
             q_rate_fcn_cHL(iyear)=(1.0+(iyear-styr_cpue_cHL)*q_rate)*q_rate_fcn_cHL(styr_cpue_cHL);  //linear
          }
          if (iyear>2003) {q_rate_fcn_cHL(iyear)=q_rate_fcn_cHL(iyear-1);} 
      }   
      for (iyear=styr_cpue_sTV; iyear<=endyr_cpue_sTV; iyear++) //KC added this section
      {   if (iyear>styr_cpue_sTV & iyear <=2003) 
          {//q_rate_fcn_sTV(iyear)=(1.0+q_rate)*q_rate_fcn_sTV(iyear-1); //compound
             q_rate_fcn_sTV(iyear)=(1.0+(iyear-styr_cpue_sTV)*q_rate)*q_rate_fcn_sTV(styr_cpue_sTV);  //linear
          }
          if (iyear>2003) {q_rate_fcn_sTV(iyear)=q_rate_fcn_sTV(iyear-1);} 
      } 
      for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++)
      {   if (iyear>styr_cpue_rHB & iyear <=2003) 
          {//q_rate_fcn_rHB(iyear)=(1.0+q_rate)*q_rate_fcn_rHB(iyear-1); //compound
             q_rate_fcn_rHB(iyear)=(1.0+(iyear-styr_cpue_rHB)*q_rate)*q_rate_fcn_rHB(styr_cpue_rHB);  //linear
          }
          if (iyear>2003) {q_rate_fcn_rHB(iyear)=q_rate_fcn_rHB(iyear-1);} 
      }   
      for (iyear=styr_cpue_rGN; iyear<=endyr_cpue_rGN; iyear++)
      {   if (iyear>styr_cpue_rGN & iyear <=2003) 
          {//q_rate_fcn_rGN(iyear)=(1.0+q_rate)*q_rate_fcn_rGN(iyear-1); //compound
             q_rate_fcn_rGN(iyear)=(1.0+(iyear-styr_cpue_rGN)*q_rate)*q_rate_fcn_rGN(styr_cpue_rGN);  //linear
          }
          if (iyear>2003) {q_rate_fcn_rGN(iyear)=q_rate_fcn_rGN(iyear-1);} 
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
 //cHL  cpue
  q_cHL(styr_cpue_cHL)=mfexp(log_q_cpue_cHL); 
  for (iyear=styr_cpue_cHL; iyear<=endyr_cpue_cHL; iyear++)
  {//index in weight units. original index in lb and re-scaled. predicted in klb whole weight, but difference in lb and klb is absorbed by q
      N_cHL(iyear)=elem_prod(elem_prod(N_mdyr(iyear),sel_cHL(iyear)),wholewgt_cHL_klb(iyear));   
      //N_cHL(iyear)=elem_prod(N_mdyr(iyear),sel_cHL(iyear)); 
      pred_cHL_cpue(iyear)=q_cHL(iyear)*q_rate_fcn_cHL(iyear)*q_DD_fcn(iyear)*sum(N_cHL(iyear));
      if (iyear<endyr_cpue_cHL){q_cHL(iyear+1)=q_cHL(iyear)*mfexp(q_RW_log_dev_cHL(iyear));}
  }
 //sTV  cpue 
  q_sTV(styr_cpue_sTV)=mfexp(log_q_cpue_sTV); 
  for (iyear=styr_cpue_sTV; iyear<=endyr_cpue_sTV; iyear++)
  {//index in weight units. original index in lb and re-scaled. predicted in klb whole weight, but difference in lb and klb is absorbed by q
      //N_sTV(iyear)=elem_prod(elem_prod(N_mdyr(iyear),sel_sTV(iyear)),wholewgt_sTV_klb(iyear));   
      N_sTV(iyear)=elem_prod(N_mdyr(iyear),sel_sTV(iyear));  
      pred_sTV_cpue(iyear)=q_sTV(iyear)*q_rate_fcn_sTV(iyear)*q_DD_fcn(iyear)*sum(N_sTV(iyear));
      if (iyear<endyr_cpue_sTV){q_sTV(iyear+1)=q_sTV(iyear)*mfexp(q_RW_log_dev_sTV(iyear));}
  }
 //rHB  cpue
  q_rHB(styr_cpue_rHB)=mfexp(log_q_cpue_rHB); 
  for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++)
  {   
      N_rHB(iyear)=elem_prod(N_mdyr(iyear),sel_rHB(iyear)); 
      pred_rHB_cpue(iyear)=q_rHB(iyear)*q_rate_fcn_rHB(iyear)*q_DD_fcn(iyear)*sum(N_rHB(iyear));
      if (iyear<endyr_cpue_rHB){q_rHB(iyear+1)=q_rHB(iyear)*mfexp(q_RW_log_dev_rHB(iyear));}
  }
 //rGN  cpue
  q_rGN(styr_cpue_rGN)=mfexp(log_q_cpue_rGN); 
  for (iyear=styr_cpue_rGN; iyear<=endyr_cpue_rGN; iyear++)
  {   
      N_rGN(iyear)=elem_prod(N_mdyr(iyear),sel_rHB(iyear)); //rGN uses rHB selex
      pred_rGN_cpue(iyear)=q_rGN(iyear)*q_rate_fcn_rGN(iyear)*q_DD_fcn(iyear)*sum(N_rGN(iyear));
      if (iyear<endyr_cpue_rGN){q_rGN(iyear+1)=q_rGN(iyear)*mfexp(q_RW_log_dev_rGN(iyear));}
  } 
}

void model_parameters::get_length_comps(void)
{
 //cGN handline
  for (iyear=1;iyear<=nyr_lenc_cHL;iyear++) 
  {pred_cHL_lenc(iyear)=(L_cHL_num(yrs_lenc_cHL(iyear))*lenprob_cHL)/sum(L_cHL_num(yrs_lenc_cHL(iyear)));}
 //SERFS chevron trap/video
  for (iyear=1;iyear<=nyr_lenc_sTV;iyear++) 
  {pred_sTV_lenc(iyear)=(N_sTV(yrs_lenc_sTV(iyear))*lenprob_sTV)/sum(N_sTV(yrs_lenc_sTV(iyear)));}
 //headboat
  for (iyear=1;iyear<=nyr_lenc_rHB;iyear++) 
  {pred_rHB_lenc(iyear)=(L_rHB_num(yrs_lenc_rHB(iyear))*lenprob_rHB)/sum(L_rHB_num(yrs_lenc_rHB(iyear)));}
 //headboat discards
  for (iyear=1;iyear<=nyr_lenc_rHB_D;iyear++) 
  {pred_rHB_D_lenc(iyear)=(D_rHB_num(yrs_lenc_rHB_D(iyear))*lenprob_rHB_D)/sum(D_rHB_num(yrs_lenc_rHB_D(iyear)));} 
 //rGN (MRIP)
  for (iyear=1;iyear<=nyr_lenc_rGN;iyear++) 
  {pred_rGN_lenc(iyear)=(L_rGN_num(yrs_lenc_rGN(iyear))*lenprob_rGN)/sum(L_rGN_num(yrs_lenc_rGN(iyear)));}
       // cout << pred_rHB_D_lenc <<endl;
}

void model_parameters::get_age_comps(void)
{
  //Commercial handline
  for (iyear=1;iyear<=nyr_agec_cHL;iyear++) 
  {
    ErrorFree_cHL_agec(iyear)=L_cHL_num(yrs_agec_cHL(iyear))/sum(L_cHL_num(yrs_agec_cHL(iyear)));  
    //ErrorFree_cHL_agec(iyear)=elem_prod(N(yrs_agec_cHL(iyear)),sel_cHL(yrs_agec_cHL(iyear)));
    pred_cHL_agec_allages(iyear)=age_error*(ErrorFree_cHL_agec(iyear)/sum(ErrorFree_cHL_agec(iyear)));   
    for (iage=1; iage<=nages_agec; iage++) {pred_cHL_agec(iyear,iage)=pred_cHL_agec_allages(iyear,iage);} 
    for (iage=(nages_agec+1); iage<=nages; iage++) {pred_cHL_agec(iyear,nages_agec)+=pred_cHL_agec_allages(iyear,iage);} //plus group                             
  }
  //SERFS chevron trap/video
  for (iyear=1;iyear<=nyr_agec_sTV;iyear++) 
  {
    ErrorFree_sTV_agec(iyear)=N_sTV(yrs_agec_sTV(iyear))/sum(N_sTV(yrs_agec_sTV(iyear)));
    //ErrorFree_sTV_agec(iyear)=elem_prod(N(yrs_agec_sTV(iyear)),sel_sTV(yrs_agec_sTV(iyear)));
    pred_sTV_agec_allages(iyear)=age_error*(ErrorFree_sTV_agec(iyear)/sum(ErrorFree_sTV_agec(iyear)));     
    for (iage=1; iage<=nages_agec; iage++) {pred_sTV_agec(iyear,iage)=pred_sTV_agec_allages(iyear,iage);} 
    for (iage=(nages_agec+1); iage<=nages; iage++) {pred_sTV_agec(iyear,nages_agec)+=pred_sTV_agec_allages(iyear,iage);} //plus group                           
  }
 //Headboat
 for (iyear=1;iyear<=nyr_agec_rHB;iyear++)
  {
    ErrorFree_rHB_agec(iyear)=L_rHB_num(yrs_agec_rHB(iyear))/sum(L_rHB_num(yrs_agec_rHB(iyear)));
    pred_rHB_agec_allages(iyear)=age_error*ErrorFree_rHB_agec(iyear); 
    for (iage=1; iage<=nages_agec; iage++) {pred_rHB_agec(iyear,iage)=pred_rHB_agec_allages(iyear,iage);} 
    for (iage=(nages_agec+1); iage<=nages; iage++) {pred_rHB_agec(iyear,nages_agec)+=pred_rHB_agec_allages(iyear,iage);} //plus group                        
  }
 // //rGN (MRIP)
 // for (iyear=1;iyear<=nyr_rGN_agec;iyear++)
  // {
    // ErrorFree_rGN_agec(iyear)=L_rGN_num(yrs_rGN_agec(iyear))/sum(L_rGN_num(yrs_rGN_agec(iyear)));
    // pred_rGN_agec_allages(iyear)=age_error*ErrorFree_rGN_agec(iyear); 
    // for (iage=1; iage<=nages_agec; iage++) {pred_rGN_agec(iyear,iage)=pred_rGN_agec_allages(iyear,iage);} 
    // for (iage=(nages_agec+1); iage<=nages; iage++) {pred_rGN_agec(iyear,nages_agec)+=pred_rGN_agec_allages(iyear,iage);} //plus group                        
  // }  
}

void model_parameters::get_weighted_current(void)
{
  F_temp_sum=0.0;
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cHL+
        sum(log_dev_F_L_cHL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rHB+
        sum(log_dev_F_L_rHB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rGN+
        sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_D_rHB+
        sum(log_dev_F_D_rHB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_D_rGN+
        sum(log_dev_F_D_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_cHL_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cHL+
        sum(log_dev_F_L_cHL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_rHB_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rHB+
        sum(log_dev_F_L_rHB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_rGN_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rGN+
        sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_rHB_D_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_D_rHB+
        sum(log_dev_F_D_rHB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_rGN_D_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_D_rGN+
        sum(log_dev_F_D_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  log_F_dev_end_cHL=sum(log_dev_F_L_cHL((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_rHB=sum(log_dev_F_L_rHB((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  
  log_F_dev_end_rGN=sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  
  log_F_dev_end_rHB_D=sum(log_dev_F_D_rHB((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_rGN_D=sum(log_dev_F_D_rGN((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  F_end_L=sel_cHL(endyr)*mfexp(log_avg_F_L_cHL+log_F_dev_end_cHL)+
          sel_rHB(endyr)*mfexp(log_avg_F_L_rHB+log_F_dev_end_rHB)+
          sel_rGN(endyr)*mfexp(log_avg_F_L_rGN+log_F_dev_end_rGN);     
  F_end_D=sel_rHB_D(endyr)*mfexp(log_avg_F_D_rHB+log_F_dev_end_rHB_D)+
          sel_rGN_D(endyr)*mfexp(log_avg_F_D_rGN+log_F_dev_end_rGN_D);  
  F_end=F_end_L+F_end_D;
  F_end_apex=max(F_end);
  sel_wgted_tot=F_end/F_end_apex;
  sel_wgted_L=elem_prod(sel_wgted_tot, elem_div(F_end_L,F_end));
  sel_wgted_D=elem_prod(sel_wgted_tot, elem_div(F_end_D,F_end));
  wgt_wgted_L_denom=F_cHL_prop+F_rHB_prop+F_rGN_prop; 
  wgt_wgted_L_klb=F_cHL_prop/wgt_wgted_L_denom*wholewgt_cHL_klb(endyr)+      
                  F_rHB_prop/wgt_wgted_L_denom*wholewgt_rHB_klb(endyr)+  
                  F_rGN_prop/wgt_wgted_L_denom*wholewgt_rGN_klb(endyr);                       
  wgt_wgted_D_denom=F_rHB_D_prop+F_rGN_D_prop;
  wgt_wgted_D_klb=F_rHB_D_prop/wgt_wgted_D_denom*wholewgt_rHB_D_klb(endyr)+  
                 F_rGN_D_prop/wgt_wgted_D_denom*wholewgt_rGN_D_klb(endyr); 
}

void model_parameters::get_msy(void)
{
  //compute values as functions of F
  for(ff=1; ff<=n_iter_msy; ff++)
  {
    //uses fishery-weighted F's
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
    L_eq_klb(ff)=sum(elem_prod(L_age_msy,wgt_wgted_L_klb)); //in whole weight
    L_eq_knum(ff)=sum(L_age_msy)/1000.0;  
    D_eq_klb(ff)=sum(elem_prod(D_age_msy,wgt_wgted_D_klb)); //in whole weight   
    D_eq_knum(ff)=sum(D_age_msy)/1000.0;    
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
        D_msy_knum_out=D_eq_knum(ff);
        D_msy_klb_out=D_eq_klb(ff);        
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
    Z_age_spr=M+F_L_age_spr+F_spr(ff)*sel_wgted_D;
    N_age_spr(1)=1.0;
    for (iage=2; iage<=nages; iage++)
     {N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z_age_spr(iage-1));}
    N_age_spr(nages)=N_age_spr(nages)/(1-mfexp(-1.*Z_age_spr(nages)));
    N_age_spr_spawn(1,(nages-1))=elem_prod(N_age_spr(1,(nages-1)),
                                   mfexp((-1.*Z_age_spr(1,(nages-1)))*spawn_time_frac));                 
    N_age_spr_spawn(nages)=(N_age_spr_spawn(nages-1)*
                          (mfexp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_spr(nages)));
    spr_spr(ff)=sum(elem_prod(N_age_spr_spawn,reprod));  //KC4
    L_spr(ff)=0.0;
    for (iage=1; iage<=nages; iage++)
    {
      L_age_spr(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));
      L_spr(ff)+=L_age_spr(iage)*wgt_wgted_L_klb(iage)*1000.0; //in lb gutted wgt
    }   
  }
  spr_ratio=spr_spr/spr_F0;	//added 1-19-16 here to...
  F20_dum=min(fabs(spr_ratio-0.2));
  F30_dum=min(fabs(spr_ratio-0.3));
  F40_dum=min(fabs(spr_ratio-0.4));
  for(ff=1; ff<=n_iter_spr; ff++)
  {   
      if (fabs(spr_ratio(ff)-0.2)==F20_dum) {F20_out=F_spr(ff);
	  }
	  if (fabs(spr_ratio(ff)-0.3)==F30_dum) {
		F30_out=F_spr(ff);
        SSB_F30_out=SSB_eq(ff);  //NOTE, this works bc F grid for msy calcs is the same as for spr calcs
		B_F30_out=B_eq(ff);
        R_F30_out=R_eq(ff);
        L_F30_knum_out=L_eq_knum(ff);
        L_F30_klb_out=L_eq_klb(ff);
		D_F30_knum_out=D_eq_knum(ff);
        D_F30_klb_out=D_eq_klb(ff);          	  
	  }
	  if (fabs(spr_ratio(ff)-0.4)==F40_dum) {
		F40_out=F_spr(ff);
	  }
  }  //here.
}

void model_parameters::get_miscellaneous_stuff(void)
{
  if(var_rec_dev>0.0)
   {sigma_rec_dev=sqrt(var_rec_dev);} //pow(var_rec_dev,0.5);  //sample SD of predicted residuals (may not equal rec_sigma)  
   else{sigma_rec_dev=0.0;}
  len_sd=elem_prod(len_cv,meanlen_FL);
 for(iage=1;iage<=nages;iage++)		//Added this loop to allow plotting of both Population and Landings length at age
    	{
    	len_cv(iage)=len_cv_val;
    	len_sd(iage)=meanlen_FL(iage)*len_cv(iage);
    	len_cv_L(iage)=len_cv_val_L;
	len_sd_L(iage)=meanlen_FL_L(iage)*len_cv_L(iage);
        len_cv_sTV(iage)=len_cv_val_sTV;
        len_sd_sTV(iage)=meanlen_FL_sTV(iage)*len_cv_sTV(iage);
  	}
  //compute total landings- and discards-at-age in 1000 fish and klb whole weight
  L_total_num.initialize();
  L_total_klb.initialize();
  L_total_knum_yr.initialize();
  L_total_klb_yr.initialize();  
  D_total_num.initialize();
  D_total_klb.initialize();
  D_total_knum_yr.initialize();
  D_total_klb_yr.initialize();
  for(iyear=styr; iyear<=endyr; iyear++)
  {
             L_total_klb_yr(iyear)=pred_cHL_L_klb(iyear)+pred_rHB_L_klb(iyear)+pred_rGN_L_klb(iyear);  
             L_total_knum_yr(iyear)=pred_cHL_L_knum(iyear)+pred_rHB_L_knum(iyear)+pred_rGN_L_knum(iyear); 
        B(iyear)=elem_prod(N(iyear),wgt_mt);
        totN(iyear)=sum(N(iyear));
        totB(iyear)=sum(B(iyear));   
        if (iyear>=styr_D_rHB && iyear<=endyr_D_rHB)
        {
         D_total_knum_yr(iyear)+=pred_rHB_D_knum(iyear);
         D_total_klb_yr(iyear)+=pred_rHB_D_klb(iyear);
         D_rHB_klb(iyear)=elem_prod(D_rHB_num(iyear),wholewgt_rHB_D_klb(iyear));     
        }    
        if (iyear>=styr_D_rGN && iyear<=endyr_D_rGN)
        {
         D_total_knum_yr(iyear)+=pred_rGN_D_knum(iyear);
         D_total_klb_yr(iyear)+=pred_rGN_D_klb(iyear);
         D_rGN_klb(iyear)=elem_prod(D_rGN_num(iyear),wholewgt_rGN_D_klb(iyear));    
        }                          
       }
  L_total_num=L_cHL_num+L_rHB_num+L_rGN_num;   //added landings at age in number fish
  L_total_klb=L_cHL_klb+L_rHB_klb+L_rGN_klb;   //landings at age in klb whole weight
  D_total_num=(D_rHB_num+D_rGN_num);          //discards at age in number fish
  D_total_klb=D_rHB_klb+D_rGN_klb;            //discards at age in klb whole weight
  //Time series of interest
  B(endyr+1)=elem_prod(N(endyr+1),wgt_mt);
  totN(endyr+1)=sum(N(endyr+1));
  totB(endyr+1)=sum(B(endyr+1));  
  //N_spawn(endyr+1)=N(endyr+1);
  //SSB(endyr+1)=sum(elem_prod(N_spawn(endyr+1),reprod));  
  //MatFemB(endyr+1)=sum(elem_prod(N_spawn(endyr+1),reprod2)); 
  rec=column(N,1);
  SdS0=SSB/S0;
  Fend_mean_temp=1.0;   //added 1-19-16 here to...
  for (iyear=1; iyear<=selpar_n_yrs_wgted; iyear++) {Fend_mean_temp*=Fapex(endyr-iyear+1);}
  Fend_mean=pow(Fend_mean_temp,(1.0/selpar_n_yrs_wgted));	 //here. 
  if(F_msy_out>0)
    {
      FdF_msy=Fapex/F_msy_out;
      FdF_msy_end=FdF_msy(endyr);
	  FdF_msy_end_mean=Fend_mean/F_msy_out;  //added 1-19-16
      //FdF_msy_end_mean=pow((FdF_msy(endyr)*FdF_msy(endyr-1)*FdF_msy(endyr-2)),(1.0/3.0));
    }
  if(SSB_msy_out>0)
    {
      SdSSB_msy=SSB/SSB_msy_out;
      SdSSB_msy_end=SdSSB_msy(endyr);
    }  
  if(F30_out>0)                  //added 1-19-16  here to...
    {
	  FdF30=Fapex/F30_out;
	  FdF30_end_mean=Fend_mean/F30_out;
	}
  if(SSB_F30_out>0)
    {
      SdSSB_F30=SSB/SSB_F30_out;
	  Sdmsst_F30=SSB/(smsy2msst*SSB_F30_out);  //KC-changed from smsy2msst75
      SdSSB_F30_end=SdSSB_F30(endyr);
	  Sdmsst_F30_end=Sdmsst_F30(endyr);
    }  							//here
   //fill in log recruitment deviations for yrs they are nonzero
   for(iyear=styr_rec_dev; iyear<=endyr_rec_dev; iyear++)
     {log_rec_dev_output(iyear)=log_dev_rec(iyear);}
   //fill in log Nage deviations for ages they are nonzero (ages2+)
   for(iage=2; iage<=nages; iage++)
     {log_Nage_dev_output(iage)=log_dev_Nage(iage);}
}

void model_parameters::get_effective_sample_sizes(void)
{
      neff_cHL_lenc_allyr_out=missing;
      neff_sTV_lenc_allyr_out=missing;
      neff_rHB_lenc_allyr_out=missing;
      neff_rHB_D_lenc_allyr_out=missing;  
      neff_rGN_lenc_allyr_out=missing;  
      neff_cHL_agec_allyr_out=missing;
      neff_sTV_agec_allyr_out=missing;
      neff_rHB_agec_allyr_out=missing;
      //neff_rGN_agec_allyr_out=missing;
      for (iyear=1; iyear<=nyr_lenc_cHL; iyear++)
         {if (nsamp_lenc_cHL(iyear)>=minSS_lenc_cHL)
            {neff_cHL_lenc_allyr_out(yrs_lenc_cHL(iyear))=multinom_eff_N(pred_cHL_lenc(iyear),obs_lenc_cHL(iyear));} 
          else {neff_cHL_lenc_allyr_out(yrs_lenc_cHL(iyear))=-99;}
         }
      for (iyear=1; iyear<=nyr_lenc_sTV; iyear++)
         {if (nsamp_lenc_sTV(iyear)>=minSS_lenc_sTV)
            {neff_sTV_lenc_allyr_out(yrs_lenc_sTV(iyear))=multinom_eff_N(pred_sTV_lenc(iyear),obs_lenc_sTV(iyear));} 
          else {neff_sTV_lenc_allyr_out(yrs_lenc_sTV(iyear))=-99;}
         }
      for (iyear=1; iyear<=nyr_lenc_rHB; iyear++)
         {if (nsamp_lenc_rHB(iyear)>=minSS_lenc_rHB)
            {neff_rHB_lenc_allyr_out(yrs_lenc_rHB(iyear))=multinom_eff_N(pred_rHB_lenc(iyear),obs_lenc_rHB(iyear));} 
          else {neff_rHB_lenc_allyr_out(yrs_lenc_rHB(iyear))=-99;}
         }
      for (iyear=1; iyear<=nyr_lenc_rHB_D; iyear++)  
         {if (nsamp_lenc_rHB_D(iyear)>=minSS_lenc_rHB_D)  
            {neff_rHB_D_lenc_allyr_out(yrs_lenc_rHB_D(iyear))=multinom_eff_N(pred_rHB_D_lenc(iyear),obs_lenc_rHB_D(iyear));}  
          else {neff_rHB_D_lenc_allyr_out(yrs_lenc_rHB_D(iyear))=-99;}  
         }  
      for (iyear=1; iyear<=nyr_lenc_rGN; iyear++)
         {if (nsamp_lenc_rGN(iyear)>=minSS_lenc_rGN)
            {neff_rGN_lenc_allyr_out(yrs_lenc_rGN(iyear))=multinom_eff_N(pred_rGN_lenc(iyear),obs_lenc_rGN(iyear));} 
          else {neff_rGN_lenc_allyr_out(yrs_lenc_rGN(iyear))=-99;}
         }
      for (iyear=1; iyear<=nyr_agec_cHL; iyear++)
         {if (nsamp_agec_cHL(iyear)>=minSS_agec_cHL)
            {neff_cHL_agec_allyr_out(yrs_agec_cHL(iyear))=multinom_eff_N(pred_cHL_agec(iyear),obs_agec_cHL(iyear));}                            
          else {neff_cHL_agec_allyr_out(yrs_agec_cHL(iyear))=-99;}
         }    
      for (iyear=1; iyear<=nyr_lenc_sTV; iyear++)
         {if (nsamp_lenc_sTV(iyear)>=minSS_lenc_sTV)
            {neff_sTV_lenc_allyr_out(yrs_lenc_sTV(iyear))=multinom_eff_N(pred_sTV_lenc(iyear),obs_lenc_sTV(iyear));} 
          else {neff_sTV_lenc_allyr_out(yrs_lenc_sTV(iyear))=-99;}
         }
      for (iyear=1; iyear<=nyr_agec_rHB; iyear++)
         {if (nsamp_agec_rHB(iyear)>=minSS_agec_rHB)
            {neff_rHB_agec_allyr_out(yrs_agec_rHB(iyear))=multinom_eff_N(pred_rHB_agec(iyear),obs_agec_rHB(iyear));}                            
          else {neff_rHB_agec_allyr_out(yrs_agec_rHB(iyear))=-99;}
         }
      // for (iyear=1; iyear<=nyr_rGN_agec; iyear++)
         // {if (nsamp_rGN_agec(iyear)>=minSS_rGN_agec)
            // {neff_rGN_agec_allyr_out(yrs_rGN_agec(iyear))=multinom_eff_N(pred_rGN_agec(iyear),obs_rGN_agec(iyear));}                            
          // else {neff_rGN_agec_allyr_out(yrs_rGN_agec(iyear))=-99;}
         // }    
}

void model_parameters::evaluate_objective_function(void)
{
  //fval=square(xdum-9.0);
  fval=0.0;
  fval_data=0.0;  
  f_sTV_cpue=0.0;
  f_sTV_cpue=lk_lognormal(pred_sTV_cpue, obs_cpue_sTV, obs_cv_cpue_sTV, w_cpue_sTV);
  fval+=f_sTV_cpue;
  fval_data+=f_sTV_cpue;
  //f_cHL_L in 1000 lb whole wgt
  f_cHL_L=lk_lognormal(pred_cHL_L_klb(styr_L_cHL,endyr_L_cHL), obs_L_cHL(styr_L_cHL,endyr_L_cHL),
                      obs_cv_L_cHL(styr_L_cHL,endyr_L_cHL), w_L);
  fval+=f_cHL_L;
  fval_data+=f_cHL_L;
  //f_rHB_L in 1000 fish
  f_rHB_L=lk_lognormal(pred_rHB_L_knum(styr_L_rHB,endyr_L_rHB), obs_L_rHB(styr_L_rHB,endyr_L_rHB), 
                      obs_cv_L_rHB(styr_L_rHB,endyr_L_rHB), w_L);
  fval+=f_rHB_L;
  fval_data+=f_rHB_L;  
  //f_rGN_L in 1000 fish
  f_rGN_L=lk_lognormal(pred_rGN_L_knum(styr_L_rGN,endyr_L_rGN), obs_L_rGN(styr_L_rGN,endyr_L_rGN), 
                      obs_cv_L_rGN(styr_L_rGN,endyr_L_rGN), w_L);
  fval+=f_rGN_L;
  fval_data+=f_rGN_L;  
  f_rHB_D=lk_lognormal(pred_rHB_D_knum(styr_D_rHB,endyr_D_rHB), obs_rHB_D(styr_D_rHB,endyr_D_rHB), 
                      obs_cv_D_rHB(styr_D_rHB,endyr_D_rHB), w_D);
  fval+=f_rHB_D;
  fval_data+=f_rHB_D;  
  //f_rGN_D in 1000 fish
  f_rGN_D=lk_lognormal(pred_rGN_D_knum(styr_D_rGN,endyr_D_rGN), obs_rGN_D(styr_D_rGN,endyr_D_rGN), 
                      obs_cv_D_rGN(styr_D_rGN,endyr_D_rGN), w_D);
  fval+=f_rGN_D;
  fval_data+=f_rGN_D;  
  f_cHL_lenc=lk_robust_multinomial(nsamp_lenc_cHL, pred_cHL_lenc, obs_lenc_cHL, nyr_lenc_cHL, double(nlenbins), minSS_lenc_cHL, w_lenc_cHL);
  //f_cHL_lenc=lk_multinomial(nsamp_lenc_cHL, pred_cHL_lenc, obs_lenc_cHL, nyr_lenc_cHL, minSS_lenc_cHL, w_lenc_cHL);
  fval+=f_cHL_lenc;
  fval_data+=f_cHL_lenc;
  f_sTV_lenc=lk_robust_multinomial(nsamp_lenc_sTV, pred_sTV_lenc, obs_lenc_sTV, nyr_lenc_sTV, double(nlenbins), minSS_lenc_sTV, w_lenc_sTV);
  //f_sTV_lenc=lk_multinomial(nsamp_lenc_sTV, pred_sTV_lenc, obs_lenc_sTV, nyr_lenc_sTV, minSS_lenc_sTV, w_lenc_sTV);
  fval+=f_sTV_lenc;
  fval_data+=f_sTV_lenc;
  f_rHB_lenc=lk_robust_multinomial(nsamp_lenc_rHB, pred_rHB_lenc, obs_lenc_rHB, nyr_lenc_rHB, double(nlenbins), minSS_lenc_rHB, w_lenc_rHB);
  //f_rHB_lenc=lk_multinomial(nsamp_lenc_rHB, pred_rHB_lenc, obs_lenc_rHB, nyr_lenc_rHB, minSS_lenc_rHB, w_lenc_rHB);
  fval+=f_rHB_lenc;
  fval_data+=f_rHB_lenc;
  f_rHB_D_lenc=lk_robust_multinomial(nsamp_lenc_rHB_D, pred_rHB_D_lenc, obs_lenc_rHB_D, nyr_lenc_rHB_D, double(nlenbins), minSS_lenc_rHB_D, w_lenc_rHB_D);  
  //f_rHB_D_lenc=lk_multinomial(nsamp_lenc_rHB_D, pred_rHB_D_lenc, obs_lenc_rHB_D, nyr_lenc_rHB_D, minSS_lenc_rHB_D, w_lc_rHB_d);  
  fval+=f_rHB_D_lenc;  
  fval_data+=f_rHB_D_lenc;    
  f_rGN_lenc=lk_robust_multinomial(nsamp_lenc_rGN, pred_rGN_lenc, obs_lenc_rGN, nyr_lenc_rGN, double(nlenbins), minSS_lenc_rGN, w_lenc_rGN);
  //f_rGN_lenc=lk_multinomial(nsamp_lenc_rGN, pred_rGN_lenc, obs_lenc_rGN, nyr_lenc_rGN, minSS_lenc_rGN, w_lenc_rGN);
  fval+=f_rGN_lenc;
  fval_data+=f_rGN_lenc; 
  f_cHL_agec=lk_robust_multinomial(nsamp_agec_cHL, pred_cHL_agec, obs_agec_cHL, nyr_agec_cHL, double(nages_agec), minSS_agec_cHL, w_agec_cHL);
  //f_cHL_agec=lk_multinomial(nsamp_agec_cHL, pred_cHL_agec, obs_agec_cHL, nyr_agec_cHL, minSS_agec_cHL, w_agec_cHL);
  fval+=f_cHL_agec;
  fval_data+=f_cHL_agec;
  f_sTV_agec=lk_robust_multinomial(nsamp_agec_sTV, pred_sTV_agec, obs_agec_sTV, nyr_agec_sTV, double(nages_agec), minSS_agec_sTV, w_agec_sTV);
  //f_sTV_agec=lk_multinomial(nsamp_agec_sTV, pred_sTV_agec, obs_agec_sTV, nyr_agec_sTV, minSS_agec_sTV, w_agec_sTV);
  fval+=f_sTV_agec;
  fval_data+=f_sTV_agec;
  f_rHB_agec=lk_robust_multinomial(nsamp_agec_rHB, pred_rHB_agec, obs_agec_rHB, nyr_agec_rHB, double(nages_agec), minSS_agec_rHB, w_agec_rHB);
  //f_rHB_agec=lk_multinomial(nsamp_agec_rHB, pred_rHB_agec, obs_agec_rHB, nyr_agec_rHB, minSS_agec_rHB, w_agec_rHB);
  fval+=f_rHB_agec;
  fval_data+=f_rHB_agec;
  // f_rGN_agec=lk_robust_multinomial(nsamp_rGN_agec, pred_rGN_agec, obs_rGN_agec, nyr_rGN_agec, double(nages_agec), minSS_rGN_agec, w_ac_rGN);
  // //f_rGN_agec=lk_multinomial(nsamp_rGN_agec, pred_rGN_agec, obs_rGN_agec, nyr_rGN_agec, minSS_rGN_agec, w_ac_rGN);
  // fval+=f_rGN_agec;
  // fval_data+=f_rGN_agec;
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
  f_priors=0.0; 
  f_priors+=neg_log_prior(Linf,set_Linf(5),set_Linf(6),set_Linf(7));
  f_priors+=neg_log_prior(K,set_K(5),set_K(6),set_K(7));
  f_priors+=neg_log_prior(t0,set_t0(5),set_t0(6),set_t0(7));
  f_priors+=neg_log_prior(len_cv_val,set_len_cv(5),set_len_cv(6),set_len_cv(7));
  f_priors+=neg_log_prior(M_constant,set_M_constant(5),set_M_constant(6),set_M_constant(7));
  f_priors+=neg_log_prior(Linf_L,set_Linf_L(5),set_Linf_L(6),set_Linf_L(7));
  f_priors+=neg_log_prior(K_L,set_K_L(5),set_K_L(6),set_K_L(7));
  f_priors+=neg_log_prior(t0_L,set_t0_L(5),set_t0_L(6),set_t0_L(7));
  f_priors+=neg_log_prior(len_cv_val_L,set_len_cv_L(5),set_len_cv_L(6),set_len_cv_L(7));
  f_priors+=neg_log_prior(Linf_sTV,set_Linf_sTV(5),set_Linf_sTV(6),set_Linf_sTV(7));
  f_priors+=neg_log_prior(K_sTV,set_K_sTV(5),set_K_sTV(6),set_K_sTV(7));
  f_priors+=neg_log_prior(t0_sTV,set_t0_sTV(5),set_t0_sTV(6),set_t0_sTV(7));
  f_priors+=neg_log_prior(len_cv_val_sTV,set_len_cv_sTV(5),set_len_cv_sTV(6),set_len_cv_sTV(7));  
  f_priors+=neg_log_prior(M_constant,set_M_constant(5),set_M_constant(6),set_M_constant(7));
  f_priors+=neg_log_prior(steep,set_steep(5),set_log_R0(6),set_log_R0(7)); 
  f_priors+=neg_log_prior(log_R0,set_log_R0(5),set_log_R0(6),set_log_R0(7)); 
  f_priors+=neg_log_prior(R_autocorr,set_R_autocorr(5),set_R_autocorr(6),set_R_autocorr(7));
  f_priors+=neg_log_prior(rec_sigma,set_rec_sigma(5),set_rec_sigma(6),set_rec_sigma(7));
  f_priors+=neg_log_prior(selpar_L50_cHL1,set_selpar_L50_cHL1(5), set_selpar_L50_cHL1(6), set_selpar_L50_cHL1(7));
  f_priors+=neg_log_prior(selpar_slope_cHL1,set_selpar_slope_cHL1(5), set_selpar_slope_cHL1(6), set_selpar_slope_cHL1(7));
  f_priors+=neg_log_prior(selpar_L50_sTV,set_selpar_L50_sTV(5), set_selpar_L50_sTV(6), set_selpar_L50_sTV(7));
  f_priors+=neg_log_prior(selpar_slope_sTV,set_selpar_slope_sTV(5), set_selpar_slope_sTV(6), set_selpar_slope_sTV(7));
  f_priors+=neg_log_prior(selpar_L51_rHB1,set_selpar_L51_rHB1(5), set_selpar_L51_rHB1(6), set_selpar_L51_rHB1(7));
  f_priors+=neg_log_prior(selpar_slope1_rHB1,set_selpar_slope1_rHB1(5), set_selpar_slope1_rHB1(6), set_selpar_slope1_rHB1(7));
  f_priors+=neg_log_prior(selpar_L52_rHB1,set_selpar_L52_rHB1(5), set_selpar_L52_rHB1(6), set_selpar_L52_rHB1(7));
  f_priors+=neg_log_prior(selpar_slope2_rHB1,set_selpar_slope2_rHB1(5), set_selpar_slope2_rHB1(6), set_selpar_slope2_rHB1(7));
  f_priors+=neg_log_prior(selpar_L50_rHB2,set_selpar_L50_rHB2(5), set_selpar_L50_rHB2(6), set_selpar_L50_rHB2(7));
  f_priors+=neg_log_prior(selpar_slope_rHB2,set_selpar_slope_rHB2(5), set_selpar_slope_rHB2(6), set_selpar_slope_rHB2(7));
  f_priors+=neg_log_prior(selpar_afull_rHB2,set_selpar_afull_rHB2(5), set_selpar_afull_rHB2(6), set_selpar_afull_rHB2(7));
  f_priors+=neg_log_prior(selpar_sigma_rHB2,set_selpar_sigma_rHB2(5), set_selpar_sigma_rHB2(6), set_selpar_sigma_rHB2(7));
  f_priors+=neg_log_prior(selpar_L50_rHB_D,set_selpar_L50_rHB_D(5), set_selpar_L50_rHB_D(6), set_selpar_L50_rHB_D(7));
  f_priors+=neg_log_prior(selpar_slope_rHB_D,set_selpar_slope_rHB_D(5), set_selpar_slope_rHB_D(6), set_selpar_slope_rHB_D(7));
  f_priors+=neg_log_prior(selpar_afull_rHB_D,set_selpar_afull_rHB_D(5), set_selpar_afull_rHB_D(6), set_selpar_afull_rHB_D(7));
  f_priors+=neg_log_prior(selpar_sigma_rHB_D,set_selpar_sigma_rHB_D(5), set_selpar_sigma_rHB_D(6), set_selpar_sigma_rHB_D(7));
  f_priors+=neg_log_prior(selpar_L51_rGN,set_selpar_L51_rGN(5), set_selpar_L51_rGN(6), set_selpar_L51_rGN(7));
  f_priors+=neg_log_prior(selpar_slope1_rGN,set_selpar_slope1_rGN(5), set_selpar_slope1_rGN(6), set_selpar_slope1_rGN(7));
  f_priors+=neg_log_prior(selpar_L52_rGN,set_selpar_L52_rGN(5), set_selpar_L52_rGN(6), set_selpar_L52_rGN(7));
  f_priors+=neg_log_prior(selpar_slope2_rGN,set_selpar_slope2_rGN(5), set_selpar_slope2_rGN(6), set_selpar_slope2_rGN(7));
  f_priors+=neg_log_prior(log_q_cpue_cHL,set_log_q_cpue_cHL(5),set_log_q_cpue_cHL(6),set_log_q_cpue_cHL(7));
  f_priors+=neg_log_prior(log_q_cpue_sTV,set_log_q_cpue_sTV(5),set_log_q_cpue_sTV(6),set_log_q_cpue_sTV(7)); 
  f_priors+=neg_log_prior(log_q_cpue_rHB,set_log_q_cpue_rHB(5),set_log_q_cpue_rHB(6),set_log_q_cpue_rHB(7));
  f_priors+=neg_log_prior(log_q_cpue_rGN,set_log_q_cpue_rGN(5),set_log_q_cpue_rGN(6),set_log_q_cpue_rGN(7));
  f_priors+=neg_log_prior(F_init,set_F_init(5),set_F_init(6),set_F_init(7));
  fval+=f_priors;
  //cout << "fval = " << fval << "  fval_data = " << fval_data << endl;
  //cout << endl;
}

dvar_vector model_parameters::logistic(const dvar_vector& ages, const dvariable& L50, const dvariable& slope)
{
  //ages=vector of ages, L50=age at 50% selectivity, slope=rate of increase
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1./(1.+mfexp(-1.*slope*(ages-L50))); //logistic;  
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
}

dvar_vector model_parameters::logistic_exponential(const dvar_vector& ages, const dvariable& L50, const dvariable& slope, const dvariable& sigma, const dvariable& joint)
{
  //ages=vector of ages, L50=age at 50% sel (ascending limb), slope=rate of increase, sigma=controls rate of descent (descending)                               
  //joint=age to join curves                                                                                                                                    
  RETURN_ARRAYS_INCREMENT();                                                                                                                                    
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());                                                                                                         
  Sel_Tmp=1.0;                                                                                                                                                  
  for (iage=1; iage<=nages; iage++)                                                                                                                             
  {                                                                                                                                                             
   if (ages(iage)<joint) {Sel_Tmp(iage)=1./(1.+mfexp(-1.*slope*(ages(iage)-L50)));}                                                                             
   if (ages(iage)>joint){Sel_Tmp(iage)=mfexp(-1.*square((ages(iage)-joint)/sigma));}                                                                            
  }                                                                                                                                                             
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);                                                                                                                                 
  RETURN_ARRAYS_DECREMENT();                                                                                                                                    
  return Sel_Tmp;   
}

dvar_vector model_parameters::logistic_double(const dvar_vector& ages, const dvariable& L501, const dvariable& slope1, const dvariable& L502, const dvariable& slope2)
{
  //ages=vector of ages, L50=age at 50% selectivity, slope=rate of increase, L502=age at 50% decrease additive to L501, slope2=slope of decrease
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=elem_prod( (1./(1.+mfexp(-1.*slope1*(ages-L501)))),(1.-(1./(1.+mfexp(-1.*slope2*(ages-(L501+L502)))))) );     
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
}

dvar_vector model_parameters::logistic_joint(const dvar_vector& ages, const dvariable& L501, const dvariable& slope1, const dvariable& L502, const dvariable& slope2, const dvariable& satval, const dvariable& joint)
{
  //ages=vector of ages, L501=age at 50% sel (ascending limb), slope1=rate of increase,L502=age at 50% sel (descending), slope1=rate of increase (ascending), 
  //satval=saturation value of descending limb, joint=location in age vector to join curves (may equal age or age + 1 if age-0 is included)
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1.0; 
  for (iage=1; iage<=nages; iage++)
  {
   if (double(iage)<joint) {Sel_Tmp(iage)=1./(1.+mfexp(-1.*slope1*(ages(iage)-L501)));}  
   if (double(iage)>joint){Sel_Tmp(iage)=1.0-(1.0-satval)/(1.+mfexp(-1.*slope2*(ages(iage)-L502)));}  
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
  dvariable small_number=0.00001;
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
  dvariable small_number=0.00001;   //KC  
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
  dvariable small_number=0.00001;   //KC 
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
      get_weighted_current();
      cout<<"got weighted"<<endl;
      get_msy();
      cout<<"got msy"<<endl;
 get_per_recruit_stuff();
      get_miscellaneous_stuff();
      cout<<"got misc stuff"<<endl;
      //cout<<"got per recruit"<<endl;  
      get_effective_sample_sizes();
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
      //cout << "BC Fmsy=" << F_msy_out<< "   BC SSBmsy=" << SSB_msy_out <<endl;
      cout <<"F status="<<FdF_msy_end<<endl;
      cout <<"Pop status="<<SdSSB_msy_end<<endl;
      cout << "h="<<steep<<"   R0="<<R0<<endl;
      //cout << "len_cv = "<<len_cv_val<<endl;
      //cout << "xdum " << xdum << endl;
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;  
     // cout << F_initial << endl;
      report << "TotalLikelihood " << fval << endl;
      report << "N" << endl;
      report << N<<endl;
      report << "F" << endl;
      report << F <<endl;
      sdnr_lc_cHL=sdnr_multinomial(nyr_lenc_cHL, lenbins, nsamp_lenc_cHL, pred_cHL_lenc, obs_lenc_cHL, w_lenc_cHL);  
      sdnr_lc_sTV=sdnr_multinomial(nyr_lenc_sTV, lenbins, nsamp_lenc_sTV, pred_sTV_lenc, obs_lenc_sTV, w_lenc_sTV);  
      sdnr_lc_rHB=sdnr_multinomial(nyr_lenc_rHB, lenbins, nsamp_lenc_rHB, pred_rHB_lenc, obs_lenc_rHB, w_lenc_rHB); 
      sdnr_lc_rHB_D=sdnr_multinomial(nyr_lenc_rHB_D, lenbins, nsamp_lenc_rHB_D, pred_rHB_D_lenc, obs_lenc_rHB_D, w_lenc_rHB_D);   
      sdnr_lc_rGN=sdnr_multinomial(nyr_lenc_rGN, lenbins, nsamp_lenc_rGN, pred_rGN_lenc, obs_lenc_rGN, w_lenc_rGN);  
      sdnr_ac_cHL=sdnr_multinomial(nyr_agec_cHL, agebins_agec, nsamp_agec_cHL, pred_cHL_agec, obs_agec_cHL, w_agec_cHL);  
      sdnr_ac_sTV=sdnr_multinomial(nyr_agec_sTV, agebins_agec, nsamp_agec_sTV, pred_sTV_agec, obs_agec_sTV, w_agec_sTV);  
      sdnr_ac_rHB=sdnr_multinomial(nyr_agec_rHB, agebins_agec, nsamp_agec_rHB, pred_rHB_agec, obs_agec_rHB, w_agec_rHB);  
      sdnr_I_sTV=sdnr_lognormal(pred_sTV_cpue, obs_cpue_sTV, obs_cv_cpue_sTV, w_cpue_sTV);  
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
       Linf_sTV_out(8)=Linf_sTV; Linf_sTV_out(1,7)=set_Linf_sTV; 
       K_sTV_out(8)=K_sTV; K_sTV_out(1,7)=set_K_sTV;
       t0_sTV_out(8)=t0_sTV; t0_sTV_out(1,7)=set_t0_sTV;
       len_cv_val_sTV_out(8)=len_cv_val_sTV; len_cv_val_sTV_out(1,7)=set_len_cv_sTV;
       log_R0_out(8)=log_R0; log_R0_out(1,7)=set_log_R0;
       M_constant_out(8)=M_constant; M_constant_out(1,7)=set_M_constant;
       steep_out(8)=steep; steep_out(1,7)=set_steep;
       rec_sigma_out(8)=rec_sigma; rec_sigma_out(1,7)=set_rec_sigma;
       R_autocorr_out(8)=R_autocorr; R_autocorr_out(1,7)=set_R_autocorr;
       selpar_L50_cHL1_out(8)=selpar_L50_cHL1; selpar_L50_cHL1_out(1,7)=set_selpar_L50_cHL1;
       selpar_slope_cHL1_out(8)=selpar_slope_cHL1; selpar_slope_cHL1_out(1,7)=set_selpar_slope_cHL1;
       selpar_L50_sTV_out(8)=selpar_L50_sTV; selpar_L50_sTV_out(1,7)=set_selpar_L50_sTV;
       selpar_slope_sTV_out(8)=selpar_slope_sTV; selpar_slope_sTV_out(1,7)=set_selpar_slope_sTV;
       selpar_L51_rHB1_out(8)=selpar_L51_rHB1; selpar_L51_rHB1_out(1,7)=set_selpar_L51_rHB1;
       selpar_slope1_rHB1_out(8)=selpar_slope1_rHB1; selpar_slope1_rHB1_out(1,7)=set_selpar_slope1_rHB1;
       selpar_L52_rHB1_out(8)=selpar_L52_rHB1; selpar_L52_rHB1_out(1,7)=set_selpar_L52_rHB1;
       selpar_slope2_rHB1_out(8)=selpar_slope2_rHB1; selpar_slope2_rHB1_out(1,7)=set_selpar_slope2_rHB1;
       selpar_L50_rHB2_out(8)=selpar_L50_rHB2; selpar_L50_rHB2_out(1,7)=set_selpar_L50_rHB2;
       selpar_slope_rHB2_out(8)=selpar_slope_rHB2; selpar_slope_rHB2_out(1,7)=set_selpar_slope_rHB2;
       selpar_afull_rHB2_out(8)=selpar_afull_rHB2; selpar_afull_rHB2_out(1,7)=set_selpar_afull_rHB2;
       selpar_sigma_rHB2_out(8)=selpar_sigma_rHB2; selpar_sigma_rHB2_out(1,7)=set_selpar_sigma_rHB2;
       selpar_L50_rHB_D_out(8)=selpar_L50_rHB_D; selpar_L50_rHB_D_out(1,7)=set_selpar_L50_rHB_D;
       selpar_slope_rHB_D_out(8)=selpar_slope_rHB_D; selpar_slope_rHB_D_out(1,7)=set_selpar_slope_rHB_D;
       selpar_afull_rHB_D_out(8)=selpar_afull_rHB_D; selpar_afull_rHB_D_out(1,7)=set_selpar_afull_rHB_D;
       selpar_sigma_rHB_D_out(8)=selpar_sigma_rHB_D; selpar_sigma_rHB_D_out(1,7)=set_selpar_sigma_rHB_D;
       selpar_L51_rGN_out(8)=selpar_L51_rGN; selpar_L51_rGN_out(1,7)=set_selpar_L51_rGN;
       selpar_slope1_rGN_out(8)=selpar_slope1_rGN; selpar_slope1_rGN_out(1,7)=set_selpar_slope1_rGN;
       selpar_L52_rGN_out(8)=selpar_L52_rGN; selpar_L52_rGN_out(1,7)=set_selpar_L52_rGN;
       selpar_slope2_rGN_out(8)=selpar_slope2_rGN; selpar_slope2_rGN_out(1,7)=set_selpar_slope2_rGN;
       log_q_cHL_out(8)=log_q_cpue_cHL; log_q_cHL_out(1,7)=set_log_q_cpue_cHL;
       log_q_sTV_out(8)=log_q_cpue_sTV; log_q_sTV_out(1,7)=set_log_q_cpue_sTV;  
       log_q_rHB_out(8)=log_q_cpue_rHB; log_q_rHB_out(1,7)=set_log_q_cpue_rHB;
       log_q_rGN_out(8)=log_q_cpue_rGN; log_q_rGN_out(1,7)=set_log_q_cpue_rGN;
       log_avg_F_cHL_out(8)=log_avg_F_L_cHL; log_avg_F_cHL_out(1,7)=set_log_avg_F_L_cHL;
       log_avg_F_rHB_out(8)=log_avg_F_L_rHB; log_avg_F_rHB_out(1,7)=set_log_avg_F_L_rHB;
       log_avg_F_rGN_out(8)=log_avg_F_L_rGN; log_avg_F_rGN_out(1,7)=set_log_avg_F_L_rGN;       
       log_avg_F_rHB_D_out(8)=log_avg_F_D_rHB; log_avg_F_rHB_D_out(1,7)=set_log_avg_F_D_rHB;
       log_avg_F_rGN_D_out(8)=log_avg_F_D_rGN; log_avg_F_rGN_D_out(1,7)=set_log_avg_F_D_rGN;
       F_init_out(8)=F_init; F_init_out(1,7)=set_F_init;
       log_rec_dev_out(styr_rec_dev, endyr_rec_dev)=log_dev_rec;
       log_F_dev_cHL_out(styr_L_cHL,endyr_L_cHL)=log_dev_F_L_cHL;
       log_F_dev_rHB_out(styr_L_rHB,endyr_L_rHB)=log_dev_F_L_rHB;
       log_F_dev_rGN_out(styr_L_rGN,endyr_L_rGN)=log_dev_F_L_rGN;
       log_F_dev_rHB_D_out(styr_D_rHB,endyr_D_rHB)=log_dev_F_D_rHB;
       log_F_dev_rGN_D_out(styr_D_rGN,endyr_D_rGN)=log_dev_F_D_rGN;
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
