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
  endyr_selex_phase1.allocate("endyr_selex_phase1");
  endyr_selex_phase2.allocate("endyr_selex_phase2");
  sizelim1.allocate("sizelim1");
  sizelim2.allocate("sizelim2");
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
  styr_L_cHL.allocate("styr_L_cHL");
  endyr_L_cHL.allocate("endyr_L_cHL");
  obs_L_cHL.allocate(styr_L_cHL,endyr_L_cHL,"obs_L_cHL");
  obs_cv_L_cHL.allocate(styr_L_cHL,endyr_L_cHL,"obs_cv_L_cHL");
  styr_L_cTW.allocate("styr_L_cTW");
  endyr_L_cTW.allocate("endyr_L_cTW");
  obs_L_cTW.allocate(styr_L_cTW,endyr_L_cTW,"obs_L_cTW");
  obs_cv_L_cTW.allocate(styr_L_cTW,endyr_L_cTW,"obs_cv_L_cTW");
  styr_L_rHB.allocate("styr_L_rHB");
  endyr_L_rHB.allocate("endyr_L_rHB");
  obs_L_rHB.allocate(styr_L_rHB,endyr_L_rHB,"obs_L_rHB");
  obs_cv_L_rHB.allocate(styr_L_rHB,endyr_L_rHB,"obs_cv_L_rHB");
  styr_L_rGN.allocate("styr_L_rGN");
  endyr_L_rGN.allocate("endyr_L_rGN");
  obs_L_rGN.allocate(styr_L_rGN,endyr_L_rGN,"obs_L_rGN");
  obs_cv_L_rGN.allocate(styr_L_rGN,endyr_L_rGN,"obs_cv_L_rGN");
  styr_D_cHL.allocate("styr_D_cHL");
  endyr_D_cHL.allocate("endyr_D_cHL");
  obs_released_cHL.allocate(styr_D_cHL,endyr_D_cHL,"obs_released_cHL");
  obs_cv_D_cHL.allocate(styr_D_cHL,endyr_D_cHL,"obs_cv_D_cHL");
  styr_D_rHB.allocate("styr_D_rHB");
  endyr_D_rHB.allocate("endyr_D_rHB");
  obs_released_rHB.allocate(styr_D_rHB,endyr_D_rHB,"obs_released_rHB");
  obs_cv_D_rHB.allocate(styr_D_rHB,endyr_D_rHB,"obs_cv_D_rHB");
  styr_D_rGN.allocate("styr_D_rGN");
  endyr_D_rGN.allocate("endyr_D_rGN");
  obs_released_rGN.allocate(styr_D_rGN,endyr_D_rGN,"obs_released_rGN");
  obs_cv_D_rGN.allocate(styr_D_rGN,endyr_D_rGN,"obs_cv_D_rGN");
  styr_cpue_sCT.allocate("styr_cpue_sCT");
  endyr_cpue_sCT.allocate("endyr_cpue_sCT");
  obs_cpue_sCT.allocate(styr_cpue_sCT,endyr_cpue_sCT,"obs_cpue_sCT");
  obs_cv_cpue_sCT.allocate(styr_cpue_sCT,endyr_cpue_sCT,"obs_cv_cpue_sCT");
  styr_cpue_rHB.allocate("styr_cpue_rHB");
  endyr_cpue_rHB.allocate("endyr_cpue_rHB");
  obs_cpue_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"obs_cpue_rHB");
  obs_cv_cpue_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"obs_cv_cpue_rHB");
  nyr_agec_sCT.allocate("nyr_agec_sCT");
  yrs_agec_sCT.allocate(1,nyr_agec_sCT,"yrs_agec_sCT");
  nsamp_agec_sCT.allocate(1,nyr_agec_sCT,"nsamp_agec_sCT");
  nfish_agec_sCT.allocate(1,nyr_agec_sCT,"nfish_agec_sCT");
  obs_agec_sCT.allocate(1,nyr_agec_sCT,1,nages_agec,"obs_agec_sCT");
  nyr_agec_cHL.allocate("nyr_agec_cHL");
  yrs_agec_cHL.allocate(1,nyr_agec_cHL,"yrs_agec_cHL");
  nsamp_agec_cHL.allocate(1,nyr_agec_cHL,"nsamp_agec_cHL");
  nfish_agec_cHL.allocate(1,nyr_agec_cHL,"nfish_agec_cHL");
  obs_agec_cHL.allocate(1,nyr_agec_cHL,1,nages_agec,"obs_agec_cHL");
  nyr_agec_rHB.allocate("nyr_agec_rHB");
  yrs_agec_rHB.allocate(1,nyr_agec_rHB,"yrs_agec_rHB");
  nsamp_agec_rHB.allocate(1,nyr_agec_rHB,"nsamp_agec_rHB");
  nfish_agec_rHB.allocate(1,nyr_agec_rHB,"nfish_agec_rHB");
  obs_agec_rHB.allocate(1,nyr_agec_rHB,1,nages_agec,"obs_agec_rHB");
  nyr_lenc_pool_cTW.allocate("nyr_lenc_pool_cTW");
  yrs_lenc_pool_cTW.allocate(1,nyr_lenc_pool_cTW,"yrs_lenc_pool_cTW");
  nsamp_lenc_pool_cTW.allocate(1,nyr_lenc_pool_cTW,"nsamp_lenc_pool_cTW");
  nfish_lenc_pool_cTW.allocate(1,nyr_lenc_pool_cTW,"nfish_lenc_pool_cTW");
  nyr_lenc_cTW.allocate("nyr_lenc_cTW");
  yrs_lenc_cTW.allocate(1,nyr_lenc_cTW,"yrs_lenc_cTW");
  nsamp_lenc_cTW.allocate(1,nyr_lenc_cTW,"nsamp_lenc_cTW");
  nfish_lenc_cTW.allocate(1,nyr_lenc_cTW,"nfish_lenc_cTW");
  obs_lenc_cTW.allocate(1,nyr_lenc_cTW,1,nlenbins,"obs_lenc_cTW");
  set_Linf.allocate(1,7,"set_Linf");
  set_K.allocate(1,7,"set_K");
  set_t0.allocate(1,7,"set_t0");
  set_len_cv_val.allocate(1,7,"set_len_cv_val");
  set_M_constant.allocate(1,7,"set_M_constant");
  set_steep.allocate(1,7,"set_steep");
  set_log_R0.allocate(1,7,"set_log_R0");
  set_R_autocorr.allocate(1,7,"set_R_autocorr");
  set_rec_sigma.allocate(1,7,"set_rec_sigma");
  set_log_dm_lenc_cTW.allocate(1,7,"set_log_dm_lenc_cTW");
  set_log_dm_agec_cHL.allocate(1,7,"set_log_dm_agec_cHL");
  set_log_dm_agec_rHB.allocate(1,7,"set_log_dm_agec_rHB");
  set_log_dm_agec_sCT.allocate(1,7,"set_log_dm_agec_sCT");
  set_A50_sel_cHL2.allocate(1,7,"set_A50_sel_cHL2");
  set_slope_sel_cHL2.allocate(1,7,"set_slope_sel_cHL2");
  set_A50_sel_cHL3.allocate(1,7,"set_A50_sel_cHL3");
  set_slope_sel_cHL3.allocate(1,7,"set_slope_sel_cHL3");
  set_A50_sel_cTW.allocate(1,7,"set_A50_sel_cTW");
  set_slope_sel_cTW.allocate(1,7,"set_slope_sel_cTW");
  set_A50_sel_rHB1.allocate(1,7,"set_A50_sel_rHB1");
  set_slope_sel_rHB1.allocate(1,7,"set_slope_sel_rHB1");
  set_A50_sel_rHB2.allocate(1,7,"set_A50_sel_rHB2");
  set_slope_sel_rHB2.allocate(1,7,"set_slope_sel_rHB2");
  set_A50_sel_rHB3.allocate(1,7,"set_A50_sel_rHB3");
  set_slope_sel_rHB3.allocate(1,7,"set_slope_sel_rHB3");
  set_A50_sel_sCT.allocate(1,7,"set_A50_sel_sCT");
  set_slope_sel_sCT.allocate(1,7,"set_slope_sel_sCT");
  set_log_q_cpue_rHB.allocate(1,7,"set_log_q_cpue_rHB");
  set_log_q_cpue_sCT.allocate(1,7,"set_log_q_cpue_sCT");
  set_log_avg_F_L_cHL.allocate(1,7,"set_log_avg_F_L_cHL");
  set_log_avg_F_L_cTW.allocate(1,7,"set_log_avg_F_L_cTW");
  set_log_avg_F_L_rHB.allocate(1,7,"set_log_avg_F_L_rHB");
  set_log_avg_F_L_rGN.allocate(1,7,"set_log_avg_F_L_rGN");
  set_log_avg_F_D_cHL.allocate(1,7,"set_log_avg_F_D_cHL");
  set_log_avg_F_D_rHB.allocate(1,7,"set_log_avg_F_D_rHB");
  set_log_avg_F_D_rGN.allocate(1,7,"set_log_avg_F_D_rGN");
  set_log_dev_F_L_cHL.allocate(1,3,"set_log_dev_F_L_cHL");
  set_log_dev_F_L_cTW.allocate(1,3,"set_log_dev_F_L_cTW");
  set_log_dev_F_L_rHB.allocate(1,3,"set_log_dev_F_L_rHB");
  set_log_dev_F_L_rGN.allocate(1,3,"set_log_dev_F_L_rGN");
  set_log_dev_F_D_cHL.allocate(1,3,"set_log_dev_F_D_cHL");
  set_log_dev_F_D_rHB.allocate(1,3,"set_log_dev_F_D_rHB");
  set_log_dev_F_D_rGN.allocate(1,3,"set_log_dev_F_D_rGN");
  set_log_dev_RWq.allocate(1,3,"set_log_dev_RWq");
  set_log_dev_rec.allocate(1,3,"set_log_dev_rec");
  set_log_dev_Nage.allocate(1,3,"set_log_dev_Nage");
  set_log_dev_vals_F_L_cHL.allocate(styr_L_cHL,endyr_L_cHL,"set_log_dev_vals_F_L_cHL");
  set_log_dev_vals_F_L_cTW.allocate(styr_L_cTW,endyr_L_cTW,"set_log_dev_vals_F_L_cTW");
  set_log_dev_vals_F_L_rHB.allocate(styr_L_rHB,endyr_L_rHB,"set_log_dev_vals_F_L_rHB");
  set_log_dev_vals_F_L_rGN.allocate(styr_L_rGN,endyr_L_rGN,"set_log_dev_vals_F_L_rGN");
  set_log_dev_vals_F_D_cHL.allocate(styr_D_cHL,endyr_D_cHL,"set_log_dev_vals_F_D_cHL");
  set_log_dev_vals_F_D_rHB.allocate(styr_D_rHB,endyr_D_rHB,"set_log_dev_vals_F_D_rHB");
  set_log_dev_vals_F_D_rGN.allocate(styr_D_rGN,endyr_D_rGN,"set_log_dev_vals_F_D_rGN");
  set_log_dev_vals_rec.allocate(styr_rec_dev,endyr_rec_dev,"set_log_dev_vals_rec");
  set_log_dev_vals_Nage.allocate(2,nages,"set_log_dev_vals_Nage");
  set_w_L.allocate("set_w_L");
  set_w_D.allocate("set_w_D");
  set_w_cpue_rHB.allocate("set_w_cpue_rHB");
  set_w_cpue_sCT.allocate("set_w_cpue_sCT");
  set_w_lenc_cTW.allocate("set_w_lenc_cTW");
  set_w_agec_cHL.allocate("set_w_agec_cHL");
  set_w_agec_rHB.allocate("set_w_agec_rHB");
  set_w_agec_sCT.allocate("set_w_agec_sCT");
  set_w_Nage_init.allocate("set_w_Nage_init");
  set_w_rec.allocate("set_w_rec");
  set_w_rec_early.allocate("set_w_rec_early");
  set_w_rec_end.allocate("set_w_rec_end");
  set_w_fullF.allocate("set_w_fullF");
  set_w_Ftune.allocate("set_w_Ftune");
  wgtpar_a.allocate("wgtpar_a");
  wgtpar_b.allocate("wgtpar_b");
  obs_maturity_f.allocate(styr,endyr,1,nages,"obs_maturity_f");
  obs_maturity_m.allocate(1,nages,"obs_maturity_m");
  obs_prop_f.allocate(1,nages,"obs_prop_f");
  spawn_time_frac.allocate("spawn_time_frac");
  set_M.allocate(1,nages,"set_M");
  set_Dmort_cHL.allocate("set_Dmort_cHL");
  set_Dmort_rHB.allocate("set_Dmort_rHB");
  set_Dmort_rGN.allocate("set_Dmort_rGN");
  SPR_rec_switch.allocate("SPR_rec_switch");
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
  minSS_lenc_cTW.allocate("minSS_lenc_cTW");
  minSS_agec_cHL.allocate("minSS_agec_cHL");
  minSS_agec_rHB.allocate("minSS_agec_rHB");
  minSS_agec_sCT.allocate("minSS_agec_sCT");
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
  const double len_cv_LO=set_len_cv_val(2); const double len_cv_HI=set_len_cv_val(3); const double len_cv_PH=set_len_cv_val(4); 
     
  const double M_constant_LO=set_M_constant(2); const double M_constant_HI=set_M_constant(3); const double M_constant_PH=set_M_constant(4);        
  const double steep_LO=set_steep(2); const double steep_HI=set_steep(3); const double steep_PH=set_steep(4);
  const double log_R0_LO=set_log_R0(2); const double log_R0_HI=set_log_R0(3); const double log_R0_PH=set_log_R0(4);
  const double R_autocorr_LO=set_R_autocorr(2); const double R_autocorr_HI=set_R_autocorr(3); const double R_autocorr_PH=set_R_autocorr(4);
  const double rec_sigma_LO=set_rec_sigma(2); const double rec_sigma_HI=set_rec_sigma(3); const double rec_sigma_PH=set_rec_sigma(4);
  
  const double log_dm_lenc_cTW_LO=set_log_dm_lenc_cTW(2); const double log_dm_lenc_cTW_HI=set_log_dm_lenc_cTW(3); const double log_dm_lenc_cTW_PH=set_log_dm_lenc_cTW(4);
  const double log_dm_agec_cHL_LO=set_log_dm_agec_cHL(2); const double log_dm_agec_cHL_HI=set_log_dm_agec_cHL(3); const double log_dm_agec_cHL_PH=set_log_dm_agec_cHL(4);
  const double log_dm_agec_rHB_LO=set_log_dm_agec_rHB(2); const double log_dm_agec_rHB_HI=set_log_dm_agec_rHB(3); const double log_dm_agec_rHB_PH=set_log_dm_agec_rHB(4);
  const double log_dm_agec_sCT_LO=set_log_dm_agec_sCT(2); const double log_dm_agec_sCT_HI=set_log_dm_agec_sCT(3); const double log_dm_agec_sCT_PH=set_log_dm_agec_sCT(4);
  
  const double A50_sel_cHL2_LO=set_A50_sel_cHL2(2); const double A50_sel_cHL2_HI=set_A50_sel_cHL2(3); const double A50_sel_cHL2_PH=set_A50_sel_cHL2(4);
  const double slope_sel_cHL2_LO=set_slope_sel_cHL2(2); const double slope_sel_cHL2_HI=set_slope_sel_cHL2(3); const double slope_sel_cHL2_PH=set_slope_sel_cHL2(4);
  const double A50_sel_cHL3_LO=set_A50_sel_cHL3(2); const double A50_sel_cHL3_HI=set_A50_sel_cHL3(3); const double A50_sel_cHL3_PH=set_A50_sel_cHL3(4);
  const double slope_sel_cHL3_LO=set_slope_sel_cHL3(2); const double slope_sel_cHL3_HI=set_slope_sel_cHL3(3); const double slope_sel_cHL3_PH=set_slope_sel_cHL3(4);
  const double A50_sel_cTW_LO=set_A50_sel_cTW(2); const double A50_sel_cTW_HI=set_A50_sel_cTW(3); const double A50_sel_cTW_PH=set_A50_sel_cTW(4);
  const double slope_sel_cTW_LO=set_slope_sel_cTW(2); const double slope_sel_cTW_HI=set_slope_sel_cTW(3); const double slope_sel_cTW_PH=set_slope_sel_cTW(4);
  
  const double A50_sel_rHB1_LO=set_A50_sel_rHB1(2); const double A50_sel_rHB1_HI=set_A50_sel_rHB1(3); const double A50_sel_rHB1_PH=set_A50_sel_rHB1(4);
  const double slope_sel_rHB1_LO=set_slope_sel_rHB1(2); const double slope_sel_rHB1_HI=set_slope_sel_rHB1(3); const double slope_sel_rHB1_PH=set_slope_sel_rHB1(4);
  const double A50_sel_rHB2_LO=set_A50_sel_rHB2(2); const double A50_sel_rHB2_HI=set_A50_sel_rHB2(3); const double A50_sel_rHB2_PH=set_A50_sel_rHB2(4);
  const double slope_sel_rHB2_LO=set_slope_sel_rHB2(2); const double slope_sel_rHB2_HI=set_slope_sel_rHB2(3); const double slope_sel_rHB2_PH=set_slope_sel_rHB2(4);
  const double A50_sel_rHB3_LO=set_A50_sel_rHB3(2); const double A50_sel_rHB3_HI=set_A50_sel_rHB3(3); const double A50_sel_rHB3_PH=set_A50_sel_rHB3(4);
  const double slope_sel_rHB3_LO=set_slope_sel_rHB3(2); const double slope_sel_rHB3_HI=set_slope_sel_rHB3(3); const double slope_sel_rHB3_PH=set_slope_sel_rHB3(4);
  
  
  const double A50_sel_sCT_LO=set_A50_sel_sCT(2); const double A50_sel_sCT_HI=set_A50_sel_sCT(3); const double A50_sel_sCT_PH=set_A50_sel_sCT(4);
  const double slope_sel_sCT_LO=set_slope_sel_sCT(2); const double slope_sel_sCT_HI=set_slope_sel_sCT(3); const double slope_sel_sCT_PH=set_slope_sel_sCT(4);
  
  
  const double log_q_L_rHB_LO=set_log_q_cpue_rHB(2); const double log_q_cpue_rHB_HI=set_log_q_cpue_rHB(3); const double log_q_cpue_rHB_PH=set_log_q_cpue_rHB(4);
  const double log_q_cpue_sCT_LO=set_log_q_cpue_sCT(2); const double log_q_cpue_sCT_HI=set_log_q_cpue_sCT(3); const double log_q_cpue_sCT_PH=set_log_q_cpue_sCT(4);
  
  const double log_avg_F_L_cHL_LO=set_log_avg_F_L_cHL(2); const double log_avg_F_L_cHL_HI=set_log_avg_F_L_cHL(3); const double log_avg_F_L_cHL_PH=set_log_avg_F_L_cHL(4);
  const double log_avg_F_L_cTW_LO=set_log_avg_F_L_cTW(2); const double log_avg_F_L_cTW_HI=set_log_avg_F_L_cTW(3); const double log_avg_F_L_cTW_PH=set_log_avg_F_L_cTW(4);
  const double log_avg_F_L_rHB_LO=set_log_avg_F_L_rHB(2); const double log_avg_F_L_rHB_HI=set_log_avg_F_L_rHB(3); const double log_avg_F_L_rHB_PH=set_log_avg_F_L_rHB(4); 
  const double log_avg_F_L_rGN_LO=set_log_avg_F_L_rGN(2); const double log_avg_F_L_rGN_HI=set_log_avg_F_L_rGN(3); const double log_avg_F_L_rGN_PH=set_log_avg_F_L_rGN(4); 
  const double log_avg_F_D_cHL_LO=set_log_avg_F_D_cHL(2); const double log_avg_F_D_cHL_HI=set_log_avg_F_D_cHL(3); const double log_avg_F_D_cHL_PH=set_log_avg_F_D_cHL(4);
  const double log_avg_F_D_rHB_LO=set_log_avg_F_D_rHB(2); const double log_avg_F_D_rHB_HI=set_log_avg_F_D_rHB(3); const double log_avg_F_D_rHB_PH=set_log_avg_F_D_rHB(4); 
  const double log_avg_F_D_rGN_LO=set_log_avg_F_D_rGN(2); const double log_avg_F_D_rGN_HI=set_log_avg_F_D_rGN(3); const double log_avg_F_D_rGN_PH=set_log_avg_F_D_rGN(4); 
  
  //-dev vectors-----------------------------------------------------------------------------------------------------------  
  const double log_F_dev_L_cHL_LO=set_log_dev_F_L_cHL(1); const double log_F_dev_L_cHL_HI=set_log_dev_F_L_cHL(2); const double log_F_dev_L_cHL_PH=set_log_dev_F_L_cHL(3);  
  const double log_F_dev_L_cTW_LO=set_log_dev_F_L_cTW(1); const double log_F_dev_L_cTW_HI=set_log_dev_F_L_cTW(2); const double log_F_dev_L_cTW_PH=set_log_dev_F_L_cTW(3);   
  const double log_F_dev_L_rHB_LO=set_log_dev_F_L_rHB(1); const double log_F_dev_L_rHB_HI=set_log_dev_F_L_rHB(2); const double log_F_dev_L_rHB_PH=set_log_dev_F_L_rHB(3);   
  const double log_F_dev_L_rGN_LO=set_log_dev_F_L_rGN(1); const double log_F_dev_L_rGN_HI=set_log_dev_F_L_rGN(2); const double log_F_dev_L_rGN_PH=set_log_dev_F_L_rGN(3);   
  
  const double log_F_dev_D_cHL_LO=set_log_dev_F_D_cHL(1); const double log_F_dev_D_cHL_HI=set_log_dev_F_D_cHL(2); const double log_F_dev_D_cHL_PH=set_log_dev_F_D_cHL(3);   
  const double log_F_dev_D_rHB_LO=set_log_dev_F_D_rHB(1); const double log_F_dev_D_rHB_HI=set_log_dev_F_D_rHB(2); const double log_F_dev_D_rHB_PH=set_log_dev_F_D_rHB(3);   
  const double log_F_dev_D_rGN_LO=set_log_dev_F_D_rGN(1); const double log_F_dev_D_rGN_HI=set_log_dev_F_D_rGN(2); const double log_F_dev_D_rGN_PH=set_log_dev_F_D_rGN(3);   
  
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
  len_cHL_mm.allocate(styr,endyr,1,nages,"len_cHL_mm");
  #ifndef NO_AD_INITIALIZE
    len_cHL_mm.initialize();
  #endif
  wholewgt_cHL_klb.allocate(styr,endyr,1,nages,"wholewgt_cHL_klb");
  #ifndef NO_AD_INITIALIZE
    wholewgt_cHL_klb.initialize();
  #endif
  len_cTW_mm.allocate(styr,endyr,1,nages,"len_cTW_mm");
  #ifndef NO_AD_INITIALIZE
    len_cTW_mm.initialize();
  #endif
  wholewgt_cTW_klb.allocate(styr,endyr,1,nages,"wholewgt_cTW_klb");
  #ifndef NO_AD_INITIALIZE
    wholewgt_cTW_klb.initialize();
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
  len_D_cHL_mm.allocate(styr,endyr,1,nages,"len_D_cHL_mm");
  #ifndef NO_AD_INITIALIZE
    len_D_cHL_mm.initialize();
  #endif
  wholewgt_D_cHL_klb.allocate(styr,endyr,1,nages,"wholewgt_D_cHL_klb");
  #ifndef NO_AD_INITIALIZE
    wholewgt_D_cHL_klb.initialize();
  #endif
  len_D_rHB_mm.allocate(styr,endyr,1,nages,"len_D_rHB_mm");
  #ifndef NO_AD_INITIALIZE
    len_D_rHB_mm.initialize();
  #endif
  wholewgt_D_rHB_klb.allocate(styr,endyr,1,nages,"wholewgt_D_rHB_klb");
  #ifndef NO_AD_INITIALIZE
    wholewgt_D_rHB_klb.initialize();
  #endif
  len_D_rGN_mm.allocate(styr,endyr,1,nages,"len_D_rGN_mm");
  #ifndef NO_AD_INITIALIZE
    len_D_rGN_mm.initialize();
  #endif
  wholewgt_D_rGN_klb.allocate(styr,endyr,1,nages,"wholewgt_D_rGN_klb");
  #ifndef NO_AD_INITIALIZE
    wholewgt_D_rGN_klb.initialize();
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
  lenprob_lenc_cTW.allocate(1,nages,1,nlenbins,"lenprob_lenc_cTW");
  #ifndef NO_AD_INITIALIZE
    lenprob_lenc_cTW.initialize();
  #endif
  len_sd.allocate(1,nages,"len_sd");
  #ifndef NO_AD_INITIALIZE
    len_sd.initialize();
  #endif
  len_cv.allocate(1,nages,"len_cv");
  #ifndef NO_AD_INITIALIZE
    len_cv.initialize();
  #endif
  pred_lenc_cTW.allocate(1,nyr_lenc_cTW,1,nlenbins,"pred_lenc_cTW");
  #ifndef NO_AD_INITIALIZE
    pred_lenc_cTW.initialize();
  #endif
  pred_lenc_cTW_yr.allocate(1,nyr_lenc_pool_cTW,1,nlenbins,"pred_lenc_cTW_yr");
  #ifndef NO_AD_INITIALIZE
    pred_lenc_cTW_yr.initialize();
  #endif
  pred_agec_cHL.allocate(1,nyr_agec_cHL,1,nages_agec,"pred_agec_cHL");
  #ifndef NO_AD_INITIALIZE
    pred_agec_cHL.initialize();
  #endif
  pred_agec_cHL_allages.allocate(1,nyr_agec_cHL,1,nages,"pred_agec_cHL_allages");
  #ifndef NO_AD_INITIALIZE
    pred_agec_cHL_allages.initialize();
  #endif
  ErrorFree_agec_cHL.allocate(1,nyr_agec_cHL,1,nages,"ErrorFree_agec_cHL");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_agec_cHL.initialize();
  #endif
  pred_agec_rHB.allocate(1,nyr_agec_rHB,1,nages_agec,"pred_agec_rHB");
  #ifndef NO_AD_INITIALIZE
    pred_agec_rHB.initialize();
  #endif
  pred_agec_rHB_allages.allocate(1,nyr_agec_rHB,1,nages,"pred_agec_rHB_allages");
  #ifndef NO_AD_INITIALIZE
    pred_agec_rHB_allages.initialize();
  #endif
  ErrorFree_agec_rHB.allocate(1,nyr_agec_rHB,1,nages,"ErrorFree_agec_rHB");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_agec_rHB.initialize();
  #endif
  pred_agec_sCT.allocate(1,nyr_agec_sCT,1,nages_agec,"pred_agec_sCT");
  #ifndef NO_AD_INITIALIZE
    pred_agec_sCT.initialize();
  #endif
  pred_agec_sCT_allages.allocate(1,nyr_agec_sCT,1,nages,"pred_agec_sCT_allages");
  #ifndef NO_AD_INITIALIZE
    pred_agec_sCT_allages.initialize();
  #endif
  ErrorFree_agec_sCT.allocate(1,nyr_agec_sCT,1,nages,"ErrorFree_agec_sCT");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_agec_sCT.initialize();
  #endif
  nsamp_lenc_cTW_allyr.allocate(styr,endyr,"nsamp_lenc_cTW_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_lenc_cTW_allyr.initialize();
  #endif
  nsamp_agec_cHL_allyr.allocate(styr,endyr,"nsamp_agec_cHL_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_agec_cHL_allyr.initialize();
  #endif
  nsamp_agec_rHB_allyr.allocate(styr,endyr,"nsamp_agec_rHB_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_agec_rHB_allyr.initialize();
  #endif
  nsamp_agec_sCT_allyr.allocate(styr,endyr,"nsamp_agec_sCT_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_agec_sCT_allyr.initialize();
  #endif
  nfish_lenc_cTW_allyr.allocate(styr,endyr,"nfish_lenc_cTW_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_lenc_cTW_allyr.initialize();
  #endif
  nfish_agec_cHL_allyr.allocate(styr,endyr,"nfish_agec_cHL_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_agec_cHL_allyr.initialize();
  #endif
  nfish_agec_rHB_allyr.allocate(styr,endyr,"nfish_agec_rHB_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_agec_rHB_allyr.initialize();
  #endif
  nfish_agec_sCT_allyr.allocate(styr,endyr,"nfish_agec_sCT_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_agec_sCT_allyr.initialize();
  #endif
  neff_lenc_cTW_allyr.allocate(styr,endyr,"neff_lenc_cTW_allyr");
  #ifndef NO_AD_INITIALIZE
    neff_lenc_cTW_allyr.initialize();
  #endif
  neff_agec_cHL_allyr.allocate(styr,endyr,"neff_agec_cHL_allyr");
  #ifndef NO_AD_INITIALIZE
    neff_agec_cHL_allyr.initialize();
  #endif
  neff_agec_rHB_allyr.allocate(styr,endyr,"neff_agec_rHB_allyr");
  #ifndef NO_AD_INITIALIZE
    neff_agec_rHB_allyr.initialize();
  #endif
  neff_agec_sCT_allyr.allocate(styr,endyr,"neff_agec_sCT_allyr");
  #ifndef NO_AD_INITIALIZE
    neff_agec_sCT_allyr.initialize();
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
  rec.allocate(styr,endyr+1,"rec");
  #ifndef NO_AD_INITIALIZE
    rec.initialize();
  #endif
  prop_f.allocate(1,nages,"prop_f");
  #ifndef NO_AD_INITIALIZE
    prop_f.initialize();
  #endif
  prop_m.allocate(1,nages,"prop_m");
  #ifndef NO_AD_INITIALIZE
    prop_m.initialize();
  #endif
  maturity_f.allocate(styr,endyr,1,nages,"maturity_f");
  #ifndef NO_AD_INITIALIZE
    maturity_f.initialize();
  #endif
  maturity_m.allocate(1,nages,"maturity_m");
  #ifndef NO_AD_INITIALIZE
    maturity_m.initialize();
  #endif
  reprod.allocate(styr,endyr,1,nages,"reprod");
  #ifndef NO_AD_INITIALIZE
    reprod.initialize();
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
  log_dm_lenc_cTW.allocate(log_dm_lenc_cTW_LO,log_dm_lenc_cTW_HI,log_dm_lenc_cTW_PH,"log_dm_lenc_cTW");
  log_dm_agec_cHL.allocate(log_dm_agec_cHL_LO,log_dm_agec_cHL_HI,log_dm_agec_cHL_PH,"log_dm_agec_cHL");
  log_dm_agec_rHB.allocate(log_dm_agec_rHB_LO,log_dm_agec_rHB_HI,log_dm_agec_rHB_PH,"log_dm_agec_rHB");
  log_dm_agec_sCT.allocate(log_dm_agec_sCT_LO,log_dm_agec_sCT_HI,log_dm_agec_sCT_PH,"log_dm_agec_sCT");
  log_dm_lenc_cTW_out.allocate(1,8,"log_dm_lenc_cTW_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_lenc_cTW_out.initialize();
  #endif
  log_dm_agec_cHL_out.allocate(1,8,"log_dm_agec_cHL_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_agec_cHL_out.initialize();
  #endif
  log_dm_agec_rHB_out.allocate(1,8,"log_dm_agec_rHB_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_agec_rHB_out.initialize();
  #endif
  log_dm_agec_sCT_out.allocate(1,8,"log_dm_agec_sCT_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_agec_sCT_out.initialize();
  #endif
  sel_cHL.allocate(styr,endyr,1,nages,"sel_cHL");
  #ifndef NO_AD_INITIALIZE
    sel_cHL.initialize();
  #endif
  sel_cHL_block1.allocate(1,nages,"sel_cHL_block1");
  #ifndef NO_AD_INITIALIZE
    sel_cHL_block1.initialize();
  #endif
  sel_cHL_block2.allocate(1,nages,"sel_cHL_block2");
  #ifndef NO_AD_INITIALIZE
    sel_cHL_block2.initialize();
  #endif
  sel_cHL_block3.allocate(1,nages,"sel_cHL_block3");
  #ifndef NO_AD_INITIALIZE
    sel_cHL_block3.initialize();
  #endif
  A50_sel_cHL2.allocate(A50_sel_cHL2_LO,A50_sel_cHL2_HI,A50_sel_cHL2_PH,"A50_sel_cHL2");
  slope_sel_cHL2.allocate(slope_sel_cHL2_LO,slope_sel_cHL2_HI,slope_sel_cHL2_PH,"slope_sel_cHL2");
  A50_sel_cHL2_out.allocate(1,8,"A50_sel_cHL2_out");
  #ifndef NO_AD_INITIALIZE
    A50_sel_cHL2_out.initialize();
  #endif
  slope_sel_cHL2_out.allocate(1,8,"slope_sel_cHL2_out");
  #ifndef NO_AD_INITIALIZE
    slope_sel_cHL2_out.initialize();
  #endif
  A50_sel_cHL3.allocate(A50_sel_cHL3_LO,A50_sel_cHL3_HI,A50_sel_cHL3_PH,"A50_sel_cHL3");
  slope_sel_cHL3.allocate(slope_sel_cHL3_LO,slope_sel_cHL3_HI,slope_sel_cHL3_PH,"slope_sel_cHL3");
  A50_sel_cHL3_out.allocate(1,8,"A50_sel_cHL3_out");
  #ifndef NO_AD_INITIALIZE
    A50_sel_cHL3_out.initialize();
  #endif
  slope_sel_cHL3_out.allocate(1,8,"slope_sel_cHL3_out");
  #ifndef NO_AD_INITIALIZE
    slope_sel_cHL3_out.initialize();
  #endif
  sel_cTW.allocate(styr,endyr,1,nages,"sel_cTW");
  #ifndef NO_AD_INITIALIZE
    sel_cTW.initialize();
  #endif
  sel_cTW_vec.allocate(1,nages,"sel_cTW_vec");
  #ifndef NO_AD_INITIALIZE
    sel_cTW_vec.initialize();
  #endif
  A50_sel_cTW.allocate(A50_sel_cTW_LO,A50_sel_cTW_HI,A50_sel_cTW_PH,"A50_sel_cTW");
  slope_sel_cTW.allocate(slope_sel_cTW_LO,slope_sel_cTW_HI,slope_sel_cTW_PH,"slope_sel_cTW");
  A50_sel_cTW_out.allocate(1,8,"A50_sel_cTW_out");
  #ifndef NO_AD_INITIALIZE
    A50_sel_cTW_out.initialize();
  #endif
  slope_sel_cTW_out.allocate(1,8,"slope_sel_cTW_out");
  #ifndef NO_AD_INITIALIZE
    slope_sel_cTW_out.initialize();
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
  sel_rHB_block3.allocate(1,nages,"sel_rHB_block3");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_block3.initialize();
  #endif
  A50_sel_rHB1.allocate(A50_sel_rHB1_LO,A50_sel_rHB1_HI,A50_sel_rHB1_PH,"A50_sel_rHB1");
  slope_sel_rHB1.allocate(slope_sel_rHB1_LO,slope_sel_rHB1_HI,slope_sel_rHB1_PH,"slope_sel_rHB1");
  A50_sel_rHB2.allocate(A50_sel_rHB2_LO,A50_sel_rHB2_HI,A50_sel_rHB2_PH,"A50_sel_rHB2");
  slope_sel_rHB2.allocate(slope_sel_rHB2_LO,slope_sel_rHB2_HI,slope_sel_rHB2_PH,"slope_sel_rHB2");
  A50_sel_rHB3.allocate(A50_sel_rHB3_LO,A50_sel_rHB3_HI,A50_sel_rHB3_PH,"A50_sel_rHB3");
  slope_sel_rHB3.allocate(slope_sel_rHB3_LO,slope_sel_rHB3_HI,slope_sel_rHB3_PH,"slope_sel_rHB3");
  A50_sel_rHB1_out.allocate(1,8,"A50_sel_rHB1_out");
  #ifndef NO_AD_INITIALIZE
    A50_sel_rHB1_out.initialize();
  #endif
  slope_sel_rHB1_out.allocate(1,8,"slope_sel_rHB1_out");
  #ifndef NO_AD_INITIALIZE
    slope_sel_rHB1_out.initialize();
  #endif
  A50_sel_rHB2_out.allocate(1,8,"A50_sel_rHB2_out");
  #ifndef NO_AD_INITIALIZE
    A50_sel_rHB2_out.initialize();
  #endif
  slope_sel_rHB2_out.allocate(1,8,"slope_sel_rHB2_out");
  #ifndef NO_AD_INITIALIZE
    slope_sel_rHB2_out.initialize();
  #endif
  A50_sel_rHB3_out.allocate(1,8,"A50_sel_rHB3_out");
  #ifndef NO_AD_INITIALIZE
    A50_sel_rHB3_out.initialize();
  #endif
  slope_sel_rHB3_out.allocate(1,8,"slope_sel_rHB3_out");
  #ifndef NO_AD_INITIALIZE
    slope_sel_rHB3_out.initialize();
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
  sel_rGN_block3.allocate(1,nages,"sel_rGN_block3");
  #ifndef NO_AD_INITIALIZE
    sel_rGN_block3.initialize();
  #endif
  sel_D_cHL.allocate(styr,endyr,1,nages,"sel_D_cHL");
  #ifndef NO_AD_INITIALIZE
    sel_D_cHL.initialize();
  #endif
  sel_D_rHB.allocate(styr,endyr,1,nages,"sel_D_rHB");
  #ifndef NO_AD_INITIALIZE
    sel_D_rHB.initialize();
  #endif
  sel_D_rGN.allocate(styr,endyr,1,nages,"sel_D_rGN");
  #ifndef NO_AD_INITIALIZE
    sel_D_rGN.initialize();
  #endif
  prob_belowsizelim_block2.allocate(1,nages,"prob_belowsizelim_block2");
  #ifndef NO_AD_INITIALIZE
    prob_belowsizelim_block2.initialize();
  #endif
  prob_belowsizelim_block3.allocate(1,nages,"prob_belowsizelim_block3");
  #ifndef NO_AD_INITIALIZE
    prob_belowsizelim_block3.initialize();
  #endif
  zscore_lsizelim1.allocate("zscore_lsizelim1");
  #ifndef NO_AD_INITIALIZE
  zscore_lsizelim1.initialize();
  #endif
  zscore_lsizelim2.allocate("zscore_lsizelim2");
  #ifndef NO_AD_INITIALIZE
  zscore_lsizelim2.initialize();
  #endif
  cprob_lsizelim1.allocate("cprob_lsizelim1");
  #ifndef NO_AD_INITIALIZE
  cprob_lsizelim1.initialize();
  #endif
  cprob_lsizelim2.allocate("cprob_lsizelim2");
  #ifndef NO_AD_INITIALIZE
  cprob_lsizelim2.initialize();
  #endif
  sel_sCT.allocate(styr,endyr,1,nages,"sel_sCT");
  #ifndef NO_AD_INITIALIZE
    sel_sCT.initialize();
  #endif
  sel_sCT_vec.allocate(1,nages,"sel_sCT_vec");
  #ifndef NO_AD_INITIALIZE
    sel_sCT_vec.initialize();
  #endif
  A50_sel_sCT.allocate(A50_sel_sCT_LO,A50_sel_sCT_HI,A50_sel_sCT_PH,"A50_sel_sCT");
  slope_sel_sCT.allocate(slope_sel_sCT_LO,slope_sel_sCT_HI,slope_sel_sCT_PH,"slope_sel_sCT");
  A50_sel_sCT_out.allocate(1,8,"A50_sel_sCT_out");
  #ifndef NO_AD_INITIALIZE
    A50_sel_sCT_out.initialize();
  #endif
  slope_sel_sCT_out.allocate(1,8,"slope_sel_sCT_out");
  #ifndef NO_AD_INITIALIZE
    slope_sel_sCT_out.initialize();
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
  pred_cpue_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"pred_cpue_rHB");
  #ifndef NO_AD_INITIALIZE
    pred_cpue_rHB.initialize();
  #endif
  N_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,1,nages,"N_rHB");
  #ifndef NO_AD_INITIALIZE
    N_rHB.initialize();
  #endif
  pred_cpue_sCT.allocate(styr_cpue_sCT,endyr_cpue_sCT,"pred_cpue_sCT");
  #ifndef NO_AD_INITIALIZE
    pred_cpue_sCT.initialize();
  #endif
  N_sCT.allocate(styr_cpue_sCT,endyr_cpue_sCT,1,nages,"N_sCT");
  #ifndef NO_AD_INITIALIZE
    N_sCT.initialize();
  #endif
  log_q_cpue_rHB.allocate(log_q_L_rHB_LO,log_q_cpue_rHB_HI,log_q_cpue_rHB_PH,"log_q_cpue_rHB");
  log_q_cpue_sCT.allocate(log_q_cpue_sCT_LO,log_q_cpue_sCT_HI,log_q_cpue_sCT_PH,"log_q_cpue_sCT");
  log_q_cpue_rHB_out.allocate(1,8,"log_q_cpue_rHB_out");
  #ifndef NO_AD_INITIALIZE
    log_q_cpue_rHB_out.initialize();
  #endif
  log_q_cpue_sCT_out.allocate(1,8,"log_q_cpue_sCT_out");
  #ifndef NO_AD_INITIALIZE
    log_q_cpue_sCT_out.initialize();
  #endif
  q_rate.allocate("q_rate");
  #ifndef NO_AD_INITIALIZE
  q_rate.initialize();
  #endif
  q_rate_fcn_cpue_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"q_rate_fcn_cpue_rHB");
  #ifndef NO_AD_INITIALIZE
    q_rate_fcn_cpue_rHB.initialize();
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
  q_RW_log_dev_cpue_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB-1,log_RWq_LO,log_RWq_HI,log_RWq_PH,"q_RW_log_dev_cpue_rHB");
  q_RW_log_dev_sCT.allocate(styr_cpue_sCT,endyr_cpue_sCT-1,"q_RW_log_dev_sCT");
  #ifndef NO_AD_INITIALIZE
    q_RW_log_dev_sCT.initialize();
  #endif
  q_cpue_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"q_cpue_rHB");
  #ifndef NO_AD_INITIALIZE
    q_cpue_rHB.initialize();
  #endif
  q_cpue_sCT.allocate(styr_cpue_sCT,endyr_cpue_sCT,"q_cpue_sCT");
  #ifndef NO_AD_INITIALIZE
    q_cpue_sCT.initialize();
  #endif
  L_cHL_num.allocate(styr,endyr,1,nages,"L_cHL_num");
  #ifndef NO_AD_INITIALIZE
    L_cHL_num.initialize();
  #endif
  L_cHL_klb.allocate(styr,endyr,1,nages,"L_cHL_klb");
  #ifndef NO_AD_INITIALIZE
    L_cHL_klb.initialize();
  #endif
  pred_L_cHL_knum.allocate(styr,endyr,"pred_L_cHL_knum");
  #ifndef NO_AD_INITIALIZE
    pred_L_cHL_knum.initialize();
  #endif
  pred_L_cHL_klb.allocate(styr,endyr,"pred_L_cHL_klb");
  #ifndef NO_AD_INITIALIZE
    pred_L_cHL_klb.initialize();
  #endif
  L_cTW_num.allocate(styr,endyr,1,nages,"L_cTW_num");
  #ifndef NO_AD_INITIALIZE
    L_cTW_num.initialize();
  #endif
  L_cTW_klb.allocate(styr,endyr,1,nages,"L_cTW_klb");
  #ifndef NO_AD_INITIALIZE
    L_cTW_klb.initialize();
  #endif
  pred_L_cTW_knum.allocate(styr,endyr,"pred_L_cTW_knum");
  #ifndef NO_AD_INITIALIZE
    pred_L_cTW_knum.initialize();
  #endif
  pred_L_cTW_klb.allocate(styr,endyr,"pred_L_cTW_klb");
  #ifndef NO_AD_INITIALIZE
    pred_L_cTW_klb.initialize();
  #endif
  L_rHB_num.allocate(styr,endyr,1,nages,"L_rHB_num");
  #ifndef NO_AD_INITIALIZE
    L_rHB_num.initialize();
  #endif
  L_rHB_klb.allocate(styr,endyr,1,nages,"L_rHB_klb");
  #ifndef NO_AD_INITIALIZE
    L_rHB_klb.initialize();
  #endif
  pred_L_rHB_knum.allocate(styr,endyr,"pred_L_rHB_knum");
  #ifndef NO_AD_INITIALIZE
    pred_L_rHB_knum.initialize();
  #endif
  pred_L_rHB_klb.allocate(styr,endyr,"pred_L_rHB_klb");
  #ifndef NO_AD_INITIALIZE
    pred_L_rHB_klb.initialize();
  #endif
  L_rGN_num.allocate(styr,endyr,1,nages,"L_rGN_num");
  #ifndef NO_AD_INITIALIZE
    L_rGN_num.initialize();
  #endif
  L_rGN_klb.allocate(styr,endyr,1,nages,"L_rGN_klb");
  #ifndef NO_AD_INITIALIZE
    L_rGN_klb.initialize();
  #endif
  pred_L_rGN_knum.allocate(styr,endyr,"pred_L_rGN_knum");
  #ifndef NO_AD_INITIALIZE
    pred_L_rGN_knum.initialize();
  #endif
  pred_L_rGN_klb.allocate(styr,endyr,"pred_L_rGN_klb");
  #ifndef NO_AD_INITIALIZE
    pred_L_rGN_klb.initialize();
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
  D_cHL_num.allocate(styr,endyr,1,nages,"D_cHL_num");
  #ifndef NO_AD_INITIALIZE
    D_cHL_num.initialize();
  #endif
  D_cHL_klb.allocate(styr,endyr,1,nages,"D_cHL_klb");
  #ifndef NO_AD_INITIALIZE
    D_cHL_klb.initialize();
  #endif
  obs_D_cHL.allocate(styr_D_cHL,endyr_D_cHL,"obs_D_cHL");
  #ifndef NO_AD_INITIALIZE
    obs_D_cHL.initialize();
  #endif
  pred_D_cHL_knum.allocate(styr,endyr,"pred_D_cHL_knum");
  #ifndef NO_AD_INITIALIZE
    pred_D_cHL_knum.initialize();
  #endif
  pred_D_cHL_klb.allocate(styr,endyr,"pred_D_cHL_klb");
  #ifndef NO_AD_INITIALIZE
    pred_D_cHL_klb.initialize();
  #endif
  D_rHB_num.allocate(styr,endyr,1,nages,"D_rHB_num");
  #ifndef NO_AD_INITIALIZE
    D_rHB_num.initialize();
  #endif
  D_rHB_klb.allocate(styr,endyr,1,nages,"D_rHB_klb");
  #ifndef NO_AD_INITIALIZE
    D_rHB_klb.initialize();
  #endif
  obs_D_rHB.allocate(styr_D_rHB,endyr_D_rHB,"obs_D_rHB");
  #ifndef NO_AD_INITIALIZE
    obs_D_rHB.initialize();
  #endif
  pred_D_rHB_knum.allocate(styr,endyr,"pred_D_rHB_knum");
  #ifndef NO_AD_INITIALIZE
    pred_D_rHB_knum.initialize();
  #endif
  pred_D_rHB_klb.allocate(styr,endyr,"pred_D_rHB_klb");
  #ifndef NO_AD_INITIALIZE
    pred_D_rHB_klb.initialize();
  #endif
  D_rGN_num.allocate(styr,endyr,1,nages,"D_rGN_num");
  #ifndef NO_AD_INITIALIZE
    D_rGN_num.initialize();
  #endif
  D_rGN_klb.allocate(styr,endyr,1,nages,"D_rGN_klb");
  #ifndef NO_AD_INITIALIZE
    D_rGN_klb.initialize();
  #endif
  obs_D_rGN.allocate(styr_D_rGN,endyr_D_rGN,"obs_D_rGN");
  #ifndef NO_AD_INITIALIZE
    obs_D_rGN.initialize();
  #endif
  pred_D_rGN_knum.allocate(styr,endyr,"pred_D_rGN_knum");
  #ifndef NO_AD_INITIALIZE
    pred_D_rGN_knum.initialize();
  #endif
  pred_D_rGN_klb.allocate(styr,endyr,"pred_D_rGN_klb");
  #ifndef NO_AD_INITIALIZE
    pred_D_rGN_klb.initialize();
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
  Dmort_cHL.allocate("Dmort_cHL");
  #ifndef NO_AD_INITIALIZE
  Dmort_cHL.initialize();
  #endif
  Dmort_rHB.allocate("Dmort_rHB");
  #ifndef NO_AD_INITIALIZE
  Dmort_rHB.initialize();
  #endif
  Dmort_rGN.allocate("Dmort_rGN");
  #ifndef NO_AD_INITIALIZE
  Dmort_rGN.initialize();
  #endif
  F_prop_L_cHL.allocate("F_prop_L_cHL");
  #ifndef NO_AD_INITIALIZE
  F_prop_L_cHL.initialize();
  #endif
  F_prop_L_cTW.allocate("F_prop_L_cTW");
  #ifndef NO_AD_INITIALIZE
  F_prop_L_cTW.initialize();
  #endif
  F_prop_L_rHB.allocate("F_prop_L_rHB");
  #ifndef NO_AD_INITIALIZE
  F_prop_L_rHB.initialize();
  #endif
  F_prop_L_rGN.allocate("F_prop_L_rGN");
  #ifndef NO_AD_INITIALIZE
  F_prop_L_rGN.initialize();
  #endif
  F_prop_D_cHL.allocate("F_prop_D_cHL");
  #ifndef NO_AD_INITIALIZE
  F_prop_D_cHL.initialize();
  #endif
  F_prop_D_rHB.allocate("F_prop_D_rHB");
  #ifndef NO_AD_INITIALIZE
  F_prop_D_rHB.initialize();
  #endif
  F_prop_D_rGN.allocate("F_prop_D_rGN");
  #ifndef NO_AD_INITIALIZE
  F_prop_D_rGN.initialize();
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
  L_age_F30.allocate(1,nages,"L_age_F30");
  #ifndef NO_AD_INITIALIZE
    L_age_F30.initialize();
  #endif
  D_age_F30.allocate(1,nages,"D_age_F30");
  #ifndef NO_AD_INITIALIZE
    D_age_F30.initialize();
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
  log_avg_F_L_cHL.allocate(log_avg_F_L_cHL_LO,log_avg_F_L_cHL_HI,log_avg_F_L_cHL_PH,"log_avg_F_L_cHL");
  log_avg_F_L_cHL_out.allocate(1,8,"log_avg_F_L_cHL_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_L_cHL_out.initialize();
  #endif
  log_dev_F_L_cHL.allocate(styr_L_cHL,endyr_L_cHL,log_F_dev_L_cHL_LO,log_F_dev_L_cHL_HI,log_F_dev_L_cHL_PH,"log_dev_F_L_cHL");
  log_F_dev_L_cHL_out.allocate(styr_L_cHL,endyr_L_cHL,"log_F_dev_L_cHL_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_L_cHL_out.initialize();
  #endif
  F_cHL.allocate(styr,endyr,1,nages,"F_cHL");
  #ifndef NO_AD_INITIALIZE
    F_cHL.initialize();
  #endif
  F_cHL_out.allocate(styr,endyr,"F_cHL_out");
  #ifndef NO_AD_INITIALIZE
    F_cHL_out.initialize();
  #endif
  log_F_dev_init_L_cHL.allocate("log_F_dev_init_L_cHL");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_init_L_cHL.initialize();
  #endif
  log_F_dev_end_L_cHL.allocate("log_F_dev_end_L_cHL");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_L_cHL.initialize();
  #endif
  log_avg_F_L_cTW.allocate(log_avg_F_L_cTW_LO,log_avg_F_L_cTW_HI,log_avg_F_L_cTW_PH,"log_avg_F_L_cTW");
  log_avg_F_L_cTW_out.allocate(1,8,"log_avg_F_L_cTW_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_L_cTW_out.initialize();
  #endif
  log_dev_F_L_cTW.allocate(styr_L_cTW,endyr_L_cTW,log_F_dev_L_cTW_LO,log_F_dev_L_cTW_HI,log_F_dev_L_cTW_PH,"log_dev_F_L_cTW");
  log_F_dev_L_cTW_out.allocate(styr_L_cTW,endyr_L_cTW,"log_F_dev_L_cTW_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_L_cTW_out.initialize();
  #endif
  F_cTW.allocate(styr,endyr,1,nages,"F_cTW");
  #ifndef NO_AD_INITIALIZE
    F_cTW.initialize();
  #endif
  F_cTW_out.allocate(styr,endyr,"F_cTW_out");
  #ifndef NO_AD_INITIALIZE
    F_cTW_out.initialize();
  #endif
  log_F_dev_init_L_cTW.allocate("log_F_dev_init_L_cTW");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_init_L_cTW.initialize();
  #endif
  log_F_dev_end_L_cTW.allocate("log_F_dev_end_L_cTW");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_L_cTW.initialize();
  #endif
  log_avg_F_L_rHB.allocate(log_avg_F_L_rHB_LO,log_avg_F_L_rHB_HI,log_avg_F_L_rHB_PH,"log_avg_F_L_rHB");
  log_avg_F_L_rHB_out.allocate(1,8,"log_avg_F_L_rHB_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_L_rHB_out.initialize();
  #endif
  log_dev_F_L_rHB.allocate(styr_L_rHB,endyr_L_rHB,log_F_dev_L_rHB_LO,log_F_dev_L_rHB_HI,log_F_dev_L_rHB_PH,"log_dev_F_L_rHB");
  log_F_dev_L_rHB_out.allocate(styr_L_rHB,endyr_L_rHB,"log_F_dev_L_rHB_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_L_rHB_out.initialize();
  #endif
  F_rHB.allocate(styr,endyr,1,nages,"F_rHB");
  #ifndef NO_AD_INITIALIZE
    F_rHB.initialize();
  #endif
  F_rHB_out.allocate(styr,endyr,"F_rHB_out");
  #ifndef NO_AD_INITIALIZE
    F_rHB_out.initialize();
  #endif
  log_F_dev_init_L_rHB.allocate("log_F_dev_init_L_rHB");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_init_L_rHB.initialize();
  #endif
  log_F_dev_end_L_rHB.allocate("log_F_dev_end_L_rHB");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_L_rHB.initialize();
  #endif
  log_avg_F_L_rGN.allocate(log_avg_F_L_rGN_LO,log_avg_F_L_rGN_HI,log_avg_F_L_rGN_PH,"log_avg_F_L_rGN");
  log_avg_F_L_rGN_out.allocate(1,8,"log_avg_F_L_rGN_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_L_rGN_out.initialize();
  #endif
  log_dev_F_L_rGN.allocate(styr_L_rGN,endyr_L_rGN,log_F_dev_L_rGN_LO,log_F_dev_L_rGN_HI,log_F_dev_L_rGN_PH,"log_dev_F_L_rGN");
  log_F_dev_L_rGN_out.allocate(styr_L_rGN,endyr_L_rGN,"log_F_dev_L_rGN_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_L_rGN_out.initialize();
  #endif
  F_rGN.allocate(styr,endyr,1,nages,"F_rGN");
  #ifndef NO_AD_INITIALIZE
    F_rGN.initialize();
  #endif
  F_rGN_out.allocate(styr,endyr,"F_rGN_out");
  #ifndef NO_AD_INITIALIZE
    F_rGN_out.initialize();
  #endif
  log_F_dev_init_L_rGN.allocate("log_F_dev_init_L_rGN");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_init_L_rGN.initialize();
  #endif
  log_F_dev_end_L_rGN.allocate("log_F_dev_end_L_rGN");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_L_rGN.initialize();
  #endif
  log_avg_F_D_cHL.allocate(log_avg_F_D_cHL_LO,log_avg_F_D_cHL_HI,log_avg_F_D_cHL_PH,"log_avg_F_D_cHL");
  log_avg_F_D_cHL_out.allocate(1,8,"log_avg_F_D_cHL_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_D_cHL_out.initialize();
  #endif
  log_dev_F_D_cHL.allocate(styr_D_cHL,endyr_D_cHL,log_F_dev_D_cHL_LO,log_F_dev_D_cHL_HI,log_F_dev_D_cHL_PH,"log_dev_F_D_cHL");
  log_F_dev_D_cHL_out.allocate(styr_D_cHL,endyr_D_cHL,"log_F_dev_D_cHL_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_D_cHL_out.initialize();
  #endif
  F_D_cHL.allocate(styr,endyr,1,nages,"F_D_cHL");
  #ifndef NO_AD_INITIALIZE
    F_D_cHL.initialize();
  #endif
  F_D_cHL_out.allocate(styr,endyr,"F_D_cHL_out");
  #ifndef NO_AD_INITIALIZE
    F_D_cHL_out.initialize();
  #endif
  log_F_dev_end_D_cHL.allocate("log_F_dev_end_D_cHL");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_D_cHL.initialize();
  #endif
  log_F_avgdev_D_cHL.allocate("log_F_avgdev_D_cHL");
  #ifndef NO_AD_INITIALIZE
  log_F_avgdev_D_cHL.initialize();
  #endif
  log_avg_F_D_rHB.allocate(log_avg_F_D_rHB_LO,log_avg_F_D_rHB_HI,log_avg_F_D_rHB_PH,"log_avg_F_D_rHB");
  log_avg_F_D_rHB_out.allocate(1,8,"log_avg_F_D_rHB_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_D_rHB_out.initialize();
  #endif
  log_dev_F_D_rHB.allocate(styr_D_rHB,endyr_D_rHB,log_F_dev_D_rHB_LO,log_F_dev_D_rHB_HI,log_F_dev_D_rHB_PH,"log_dev_F_D_rHB");
  log_F_dev_D_rHB_out.allocate(styr_D_rHB,endyr_D_rHB,"log_F_dev_D_rHB_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_D_rHB_out.initialize();
  #endif
  F_D_rHB.allocate(styr,endyr,1,nages,"F_D_rHB");
  #ifndef NO_AD_INITIALIZE
    F_D_rHB.initialize();
  #endif
  F_D_rHB_out.allocate(styr,endyr,"F_D_rHB_out");
  #ifndef NO_AD_INITIALIZE
    F_D_rHB_out.initialize();
  #endif
  log_F_dev_end_D_rHB.allocate("log_F_dev_end_D_rHB");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_D_rHB.initialize();
  #endif
  log_F_avgdev_D_rHB.allocate("log_F_avgdev_D_rHB");
  #ifndef NO_AD_INITIALIZE
  log_F_avgdev_D_rHB.initialize();
  #endif
  log_avg_F_D_rGN.allocate(log_avg_F_D_rGN_LO,log_avg_F_D_rGN_HI,log_avg_F_D_rGN_PH,"log_avg_F_D_rGN");
  log_avg_F_D_rGN_out.allocate(1,8,"log_avg_F_D_rGN_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_D_rGN_out.initialize();
  #endif
  log_dev_F_D_rGN.allocate(styr_D_rGN,endyr_D_rGN,log_F_dev_D_rGN_LO,log_F_dev_D_rGN_HI,log_F_dev_D_rGN_PH,"log_dev_F_D_rGN");
  log_F_dev_D_rGN_out.allocate(styr_D_rGN,endyr_D_rGN,"log_F_dev_D_rGN_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_D_rGN_out.initialize();
  #endif
  F_D_rGN.allocate(styr,endyr,1,nages,"F_D_rGN");
  #ifndef NO_AD_INITIALIZE
    F_D_rGN.initialize();
  #endif
  F_D_rGN_out.allocate(styr,endyr,"F_D_rGN_out");
  #ifndef NO_AD_INITIALIZE
    F_D_rGN_out.initialize();
  #endif
  log_F_dev_end_D_rGN.allocate("log_F_dev_end_D_rGN");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_D_rGN.initialize();
  #endif
  log_F_avgdev_D_rGN.allocate("log_F_avgdev_D_rGN");
  #ifndef NO_AD_INITIALIZE
  log_F_avgdev_D_rGN.initialize();
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
  F_D_age_spr.allocate(1,nages,"F_D_age_spr");
  #ifndef NO_AD_INITIALIZE
    F_D_age_spr.initialize();
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
  spr_F0.allocate(styr,endyr,"spr_F0");
  #ifndef NO_AD_INITIALIZE
    spr_F0.initialize();
  #endif
  bpr_F0.allocate(styr,endyr,"bpr_F0");
  #ifndef NO_AD_INITIALIZE
    bpr_F0.initialize();
  #endif
  iter_inc_spr.allocate("iter_inc_spr");
  #ifndef NO_AD_INITIALIZE
  iter_inc_spr.initialize();
  #endif
  sdnr_agec_cHL.allocate("sdnr_agec_cHL");
  #ifndef NO_AD_INITIALIZE
  sdnr_agec_cHL.initialize();
  #endif
  sdnr_agec_rHB.allocate("sdnr_agec_rHB");
  #ifndef NO_AD_INITIALIZE
  sdnr_agec_rHB.initialize();
  #endif
  sdnr_agec_sCT.allocate("sdnr_agec_sCT");
  #ifndef NO_AD_INITIALIZE
  sdnr_agec_sCT.initialize();
  #endif
  sdnr_cpue_rHB.allocate("sdnr_cpue_rHB");
  #ifndef NO_AD_INITIALIZE
  sdnr_cpue_rHB.initialize();
  #endif
  sdnr_cpue_sCT.allocate("sdnr_cpue_sCT");
  #ifndef NO_AD_INITIALIZE
  sdnr_cpue_sCT.initialize();
  #endif
  w_L.allocate("w_L");
  #ifndef NO_AD_INITIALIZE
  w_L.initialize();
  #endif
  w_D.allocate("w_D");
  #ifndef NO_AD_INITIALIZE
  w_D.initialize();
  #endif
  w_cpue_rHB.allocate("w_cpue_rHB");
  #ifndef NO_AD_INITIALIZE
  w_cpue_rHB.initialize();
  #endif
  w_cpue_sCT.allocate("w_cpue_sCT");
  #ifndef NO_AD_INITIALIZE
  w_cpue_sCT.initialize();
  #endif
  w_lenc_cTW.allocate("w_lenc_cTW");
  #ifndef NO_AD_INITIALIZE
  w_lenc_cTW.initialize();
  #endif
  w_agec_cHL.allocate("w_agec_cHL");
  #ifndef NO_AD_INITIALIZE
  w_agec_cHL.initialize();
  #endif
  w_agec_rHB.allocate("w_agec_rHB");
  #ifndef NO_AD_INITIALIZE
  w_agec_rHB.initialize();
  #endif
  w_agec_sCT.allocate("w_agec_sCT");
  #ifndef NO_AD_INITIALIZE
  w_agec_sCT.initialize();
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
  f_L_cHL.allocate("f_L_cHL");
  #ifndef NO_AD_INITIALIZE
  f_L_cHL.initialize();
  #endif
  f_L_cTW.allocate("f_L_cTW");
  #ifndef NO_AD_INITIALIZE
  f_L_cTW.initialize();
  #endif
  f_L_rHB.allocate("f_L_rHB");
  #ifndef NO_AD_INITIALIZE
  f_L_rHB.initialize();
  #endif
  f_L_rGN.allocate("f_L_rGN");
  #ifndef NO_AD_INITIALIZE
  f_L_rGN.initialize();
  #endif
  f_D_cHL.allocate("f_D_cHL");
  #ifndef NO_AD_INITIALIZE
  f_D_cHL.initialize();
  #endif
  f_D_rHB.allocate("f_D_rHB");
  #ifndef NO_AD_INITIALIZE
  f_D_rHB.initialize();
  #endif
  f_D_rGN.allocate("f_D_rGN");
  #ifndef NO_AD_INITIALIZE
  f_D_rGN.initialize();
  #endif
  f_cpue_rHB.allocate("f_cpue_rHB");
  #ifndef NO_AD_INITIALIZE
  f_cpue_rHB.initialize();
  #endif
  f_cpue_sCT.allocate("f_cpue_sCT");
  #ifndef NO_AD_INITIALIZE
  f_cpue_sCT.initialize();
  #endif
  f_RWq_cpue_rHB.allocate("f_RWq_cpue_rHB");
  #ifndef NO_AD_INITIALIZE
  f_RWq_cpue_rHB.initialize();
  #endif
  f_lenc_cTW.allocate("f_lenc_cTW");
  #ifndef NO_AD_INITIALIZE
  f_lenc_cTW.initialize();
  #endif
  f_agec_cHL.allocate("f_agec_cHL");
  #ifndef NO_AD_INITIALIZE
  f_agec_cHL.initialize();
  #endif
  f_agec_rHB.allocate("f_agec_rHB");
  #ifndef NO_AD_INITIALIZE
  f_agec_rHB.initialize();
  #endif
  f_agec_sCT.allocate("f_agec_sCT");
  #ifndef NO_AD_INITIALIZE
  f_agec_sCT.initialize();
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
  D_knum_proj.allocate(styr_proj,endyr_proj,"D_knum_proj");
  #ifndef NO_AD_INITIALIZE
    D_knum_proj.initialize();
  #endif
  D_klb_proj.allocate(styr_proj,endyr_proj,"D_klb_proj");
  #ifndef NO_AD_INITIALIZE
    D_klb_proj.initialize();
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
  FD_age_proj.allocate(1,nages,"FD_age_proj");
  #ifndef NO_AD_INITIALIZE
    FD_age_proj.initialize();
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
  D_age_proj.allocate(styr_proj,endyr_proj,1,nages,"D_age_proj");
  #ifndef NO_AD_INITIALIZE
    D_age_proj.initialize();
  #endif
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{1000, 2000,3000, 10000, 10000, 10000, 10000;}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-1, 1e-2,1e-3, 1e-4, 1e-4, 1e-4, 1e-4;}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  Dmort_cHL=set_Dmort_cHL; Dmort_rHB=set_Dmort_rHB; Dmort_rGN=set_Dmort_rGN;
  for(iyear=styr_D_cHL; iyear<=endyr_D_cHL; iyear++)
	{obs_D_cHL(iyear)=Dmort_cHL*obs_released_cHL(iyear);
	}
	
  for(iyear=styr_D_rHB; iyear<=endyr_D_rHB; iyear++)
	{obs_D_rHB(iyear)=Dmort_rHB*obs_released_rHB(iyear);
	}
	
  for(iyear=styr_D_rGN; iyear<=endyr_D_rGN; iyear++)
	{obs_D_rGN(iyear)=Dmort_rGN*obs_released_rGN(iyear);
	}
 
 //Population		
  Linf=set_Linf(1);
  K=set_K(1);
  t0=set_t0(1);
  len_cv_val=set_len_cv_val(1);
  
  M=set_M; 
  M_constant=set_M_constant(1);
  smsy2msstM=1.0-M_constant;
  smsy2msst75=0.75;  
  
  log_R0=set_log_R0(1);
  steep=set_steep(1);
  R_autocorr=set_R_autocorr(1);
  rec_sigma=set_rec_sigma(1);
  
  log_dm_lenc_cTW=set_log_dm_lenc_cTW(1);
  log_dm_agec_cHL=set_log_dm_agec_cHL(1);
  log_dm_agec_rHB=set_log_dm_agec_rHB(1);
  log_dm_agec_sCT=set_log_dm_agec_sCT(1);
  
  log_q_cpue_rHB=set_log_q_cpue_rHB(1);
  log_q_cpue_sCT=set_log_q_cpue_sCT(1);
  
  q_rate=set_q_rate;
  q_rate_fcn_cpue_rHB=1.0;   
  q_DD_beta=set_q_DD_beta;
  q_DD_fcn=1.0;
  q_RW_log_dev_cpue_rHB.initialize(); 
  q_RW_log_dev_sCT.initialize();
  
   if (set_q_rate_phase<0 & q_rate!=0.0)
  {
    for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++) //LINETAG: cpue_rHB
      {   if (iyear>styr_cpue_rHB & iyear <=2003) //LINETAG: cpue_rHB
          {//q_rate_fcn_cpue_rHB(iyear)=(1.0+q_rate)*q_rate_fcn_cpue_rHB(iyear-1); //compound //LINETAG: cpue_rHB
             q_rate_fcn_cpue_rHB(iyear)=(1.0+(iyear-styr_cpue_rHB)*q_rate)*q_rate_fcn_cpue_rHB(styr_cpue_rHB);  //linear //LINETAG: cpue_rHB
          } //LINETAG: cpue_rHB
          if (iyear>2003) {q_rate_fcn_cpue_rHB(iyear)=q_rate_fcn_cpue_rHB(iyear-1);} //LINETAG: cpue_rHB
      }   //LINETAG: cpue_rHB
	  
  } //end q_rate conditional      
  w_L=set_w_L;
  w_D=set_w_D;
  
  w_cpue_rHB=set_w_cpue_rHB;
  w_cpue_sCT=set_w_cpue_sCT;
  
  w_lenc_cTW=set_w_lenc_cTW;  
    
  w_agec_cHL=set_w_agec_cHL;
  w_agec_rHB=set_w_agec_rHB;
  w_agec_sCT=set_w_agec_sCT;  
 
  
  w_Nage_init=set_w_Nage_init;
  w_rec=set_w_rec;
  w_rec_early=set_w_rec_early;
  w_rec_end=set_w_rec_end;
  w_fullF=set_w_fullF;
  w_Ftune=set_w_Ftune;
  log_avg_F_L_cHL=set_log_avg_F_L_cHL(1); //LINETAG: L_cHL
  log_avg_F_L_cTW=set_log_avg_F_L_cTW(1); //LINETAG: L_cTW
  log_avg_F_L_rHB=set_log_avg_F_L_rHB(1); //LINETAG: L_rHB 
  log_avg_F_L_rGN=set_log_avg_F_L_rGN(1); //LINETAG: L_rGN 
  log_avg_F_D_cHL=set_log_avg_F_D_cHL(1); //LINETAG: D_cHL
  log_avg_F_D_rHB=set_log_avg_F_D_rHB(1); //LINETAG: D_rHB
  log_avg_F_D_rGN=set_log_avg_F_D_rGN(1); //LINETAG: D_rGN 
    
  log_dev_F_L_cHL=set_log_dev_vals_F_L_cHL; //LINETAG: L_cHL
  log_dev_F_L_cTW=set_log_dev_vals_F_L_cTW; //LINETAG: L_cTW
  log_dev_F_L_rHB=set_log_dev_vals_F_L_rHB; //LINETAG: L_rHB
  log_dev_F_L_rGN=set_log_dev_vals_F_L_rGN; //LINETAG: L_rGN
  log_dev_F_D_cHL=set_log_dev_vals_F_D_cHL; //LINETAG: D_cHL
  log_dev_F_D_rHB=set_log_dev_vals_F_D_rHB; //LINETAG: D_rHB
  log_dev_F_D_rGN=set_log_dev_vals_F_D_rGN; //LINETAG: D_rGN
 
  A50_sel_cHL2=set_A50_sel_cHL2(1);
  slope_sel_cHL2=set_slope_sel_cHL2(1);
  A50_sel_cHL3=set_A50_sel_cHL3(1);
  slope_sel_cHL3=set_slope_sel_cHL3(1);
  A50_sel_cTW=set_A50_sel_cTW(1);
  slope_sel_cTW=set_slope_sel_cTW(1);
  
  A50_sel_rHB1=set_A50_sel_rHB1(1);
  slope_sel_rHB1=set_slope_sel_rHB1(1);
  A50_sel_rHB2=set_A50_sel_rHB2(1);
  slope_sel_rHB2=set_slope_sel_rHB2(1);
  A50_sel_rHB3=set_A50_sel_rHB3(1);
  slope_sel_rHB3=set_slope_sel_rHB3(1);
  
    
  A50_sel_sCT=set_A50_sel_sCT(1);
  slope_sel_sCT=set_slope_sel_sCT(1);
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
 maturity_m=obs_maturity_m;
 
 prop_f=obs_prop_f;
 prop_m=1.0-obs_prop_f;
  
	  nsamp_lenc_cTW_allyr=missing;  
      nsamp_agec_cHL_allyr=missing;
      nsamp_agec_rHB_allyr=missing;
	  nsamp_agec_sCT_allyr=missing;  	  
      
	  nfish_lenc_cTW_allyr=missing;  
      nfish_agec_cHL_allyr=missing;
      nfish_agec_rHB_allyr=missing;
	  nfish_agec_sCT_allyr=missing;  	  
   
      for (iyear=1; iyear<=nyr_lenc_cTW; iyear++) //LINETAG: lenc_cTW
         {if (nsamp_lenc_cTW(iyear)>=minSS_lenc_cTW) //LINETAG: lenc_cTW
           {nsamp_lenc_cTW_allyr(yrs_lenc_cTW(iyear))=nsamp_lenc_cTW(iyear); //LINETAG: lenc_cTW
            nfish_lenc_cTW_allyr(yrs_lenc_cTW(iyear))=nfish_lenc_cTW(iyear);}} //LINETAG: lenc_cTW		 	
	  for (iyear=1; iyear<=nyr_agec_cHL; iyear++) //LINETAG: agec_cHL
         {if (nsamp_agec_cHL(iyear)>=minSS_agec_cHL) //LINETAG: agec_cHL
           {nsamp_agec_cHL_allyr(yrs_agec_cHL(iyear))=nsamp_agec_cHL(iyear); //LINETAG: agec_cHL
            nfish_agec_cHL_allyr(yrs_agec_cHL(iyear))=nfish_agec_cHL(iyear);}} //LINETAG: agec_cHL
      for (iyear=1; iyear<=nyr_agec_rHB; iyear++) //LINETAG: agec_rHB
          {if (nsamp_agec_rHB(iyear)>=minSS_agec_rHB) //LINETAG: agec_rHB
            {nsamp_agec_rHB_allyr(yrs_agec_rHB(iyear))=nsamp_agec_rHB(iyear); //LINETAG: agec_rHB
             nfish_agec_rHB_allyr(yrs_agec_rHB(iyear))=nfish_agec_rHB(iyear);}}  //LINETAG: agec_rHB
	  for (iyear=1; iyear<=nyr_agec_sCT; iyear++) //LINETAG: agec_sCT  
          {if (nsamp_agec_sCT(iyear)>=minSS_agec_sCT) //LINETAG: agec_sCT
            {nsamp_agec_sCT_allyr(yrs_agec_sCT(iyear))=nsamp_agec_sCT(iyear); //LINETAG: agec_sCT
             nfish_agec_sCT_allyr(yrs_agec_sCT(iyear))=nfish_agec_sCT(iyear);}}  //LINETAG: agec_sCT
             
  F_msy(1)=0.0;  
  for (ff=2;ff<=n_iter_msy;ff++) {F_msy(ff)=F_msy(ff-1)+iter_inc_msy;}
  F_spr(1)=0.0;  
  for (ff=2;ff<=n_iter_spr;ff++) {F_spr(ff)=F_spr(ff-1)+iter_inc_spr;}
  F_cHL.initialize(); L_cHL_num.initialize();
  F_cTW.initialize(); L_cTW_num.initialize();
  F_rHB.initialize(); L_rHB_num.initialize();
  F_rGN.initialize(); L_rGN_num.initialize();
  F_D_cHL.initialize(); D_cHL_num.initialize();
  F_D_rHB.initialize(); D_rHB_num.initialize();
  F_D_rGN.initialize(); D_rGN_num.initialize();
  F_cHL_out.initialize();
  F_cTW_out.initialize();
  F_rHB_out.initialize();
  F_rGN_out.initialize();
  F_D_cHL_out.initialize();
  F_D_rHB_out.initialize();
  F_D_rGN_out.initialize();
  sel_cHL.initialize();
  sel_cTW.initialize();  
  sel_rHB.initialize();
  sel_rGN.initialize();
  sel_D_cHL.initialize();
  sel_D_rHB.initialize();
  sel_D_rGN.initialize();
  sel_sCT.initialize();
  sel_cHL_block1.initialize(); 
  sel_cHL_block2.initialize();
  sel_cHL_block3.initialize();  
  sel_cTW_vec.initialize();
  sel_rHB_block1.initialize();
  sel_rHB_block2.initialize();
  sel_rHB_block3.initialize();
  sel_rGN_block1.initialize();
  sel_rGN_block2.initialize();
  sel_rGN_block3.initialize();
  sel_sCT_vec.initialize();
  prob_belowsizelim_block2.initialize();
  prob_belowsizelim_block3.initialize();
  
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
  get_dead_discards(); 
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
    meanlen_TL=Linf*(1.0-mfexp(-K*(agebins-t0+0.5)));     
    wgt_kg=wgtpar_a*pow(meanlen_TL,wgtpar_b);             //whole wgt in kg 
    wgt_g=wgt_kg/g2kg;                                    //convert wgt in kg to weight in g    
    wgt_mt=wgt_g*g2mt;                                    //convert weight in g to weight in mt
    wgt_klb=mt2klb*wgt_mt;                                //1000 lb of whole wgt
    wgt_lb=mt2lb*wgt_mt;                                  //lb of whole wgt
}

void model_parameters::get_reprod(void)
{
  for (iyear=styr;iyear<=endyr;iyear++)
  {
   reprod(iyear)=elem_prod((elem_prod(prop_f,maturity_f(iyear))+elem_prod(prop_m,maturity_m)),wgt_mt);
  }
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
    //first length bin
	//population
    zscore_len=((lenbins(1)+0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage);
    cprob_lenvec(1)=cumd_norm(zscore_len);          //includes any probability mass below zero
    lenprob(iage,1)=cprob_lenvec(1)-cprob_lzero;    //removes any probability mass below zero
	zscore_lsizelim1=(sizelim1-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim1=cumd_norm(zscore_lsizelim1);                 //includes any probability mass below zero
	prob_belowsizelim_block2(iage)=	cprob_lsizelim1-cprob_lzero; //removes any probability mass below zero
	zscore_lsizelim2=(sizelim2-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim2=cumd_norm(zscore_lsizelim2);                 //includes any probability mass below zero
	prob_belowsizelim_block3(iage)=	cprob_lsizelim2-cprob_lzero; //removes any probability mass below zero
    //most other length bins  
    //population
    for (ilen=2;ilen<nlenbins;ilen++)
      {
        zscore_len=((lenbins(ilen)+0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage); 
		cprob_lenvec(ilen)=cumd_norm(zscore_len);
        lenprob(iage,ilen)=cprob_lenvec(ilen)-cprob_lenvec(ilen-1);
      }
    //last length bin is a plus group
	//population
    zscore_len=((lenbins(nlenbins)-0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage); 
	lenprob(iage,nlenbins)=1.0-cumd_norm(zscore_len);
      lenprob(iage)=lenprob(iage)/(1.0-cprob_lzero);  //renormalize to account for any prob mass below size=0
   }
  //fleet and survey specific length probs, all assumed here to equal the popn
  lenprob_lenc_cTW=lenprob; 
}

void model_parameters::get_weight_at_age_landings(void)
{
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    len_cHL_mm(iyear)=meanlen_TL;  
    wholewgt_cHL_klb(iyear)=wgt_klb; 
    len_cTW_mm(iyear)=meanlen_TL;  
    wholewgt_cTW_klb(iyear)=wgt_klb; 	
    len_rHB_mm(iyear)=meanlen_TL;
    wholewgt_rHB_klb(iyear)=wgt_klb;
    len_rGN_mm(iyear)=meanlen_TL;
    wholewgt_rGN_klb(iyear)=wgt_klb;
    len_D_cHL_mm(iyear)=meanlen_TL;  
    wholewgt_D_cHL_klb(iyear)=wgt_klb;
    len_D_rHB_mm(iyear)=meanlen_TL;
    wholewgt_D_rHB_klb(iyear)=wgt_klb;
    len_D_rGN_mm(iyear)=meanlen_TL;
    wholewgt_D_rGN_klb(iyear)=wgt_klb;
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
  for(iyear=styr; iyear<=endyr; iyear++)
  {
    spr_F0(iyear)=sum(elem_prod(N_spr_F0,reprod(iyear)));
    bpr_F0(iyear)=sum(elem_prod(N_spr_F0,  wgt_mt));
  }
}

void model_parameters::get_selectivity(void)
{
  sel_cHL_block1=logistic(agebins, A50_sel_cHL2, slope_sel_cHL2);
  sel_cHL_block2=logistic(agebins, A50_sel_cHL2, slope_sel_cHL2);
  sel_cHL_block3=logistic(agebins, A50_sel_cHL3, slope_sel_cHL3);
  sel_cTW_vec=logistic(agebins, A50_sel_cTW, slope_sel_cTW);
  sel_rHB_block1=logistic(agebins, A50_sel_rHB1, slope_sel_rHB1);
  sel_rHB_block2=logistic(agebins, A50_sel_rHB2, slope_sel_rHB2);
  sel_rHB_block3=logistic(agebins, A50_sel_rHB3, slope_sel_rHB3);
  sel_rGN_block1=sel_rHB_block1;
  sel_rGN_block2=sel_rHB_block2;
  sel_rGN_block3=sel_rHB_block3;
  sel_sCT_vec=logistic(agebins, A50_sel_sCT, slope_sel_sCT);
  //BLOCK 1 for selex 
  for (iyear=styr; iyear<=endyr_selex_phase1; iyear++)
   {     
    sel_cHL(iyear)=sel_cHL_block1;
	sel_cTW(iyear)=sel_cTW_vec;
    sel_rHB(iyear)=sel_rHB_block1;
    sel_rGN(iyear)=sel_rGN_block1;
	sel_sCT(iyear)=sel_sCT_vec;
   }
  //BLOCK 2 for selex
  for (iyear=(endyr_selex_phase1+1); iyear<=endyr_selex_phase2; iyear++)
   {
    sel_cHL(iyear)=sel_cHL_block2;
	sel_cTW(iyear)=sel_cTW_vec;
    sel_rHB(iyear)=sel_rHB_block2;
    sel_rGN(iyear)=sel_rGN_block2;
	sel_sCT(iyear)=sel_sCT_vec;
   }
  //BLOCK 3 for selex
   for (iyear=(endyr_selex_phase2+1); iyear<=endyr; iyear++)
   {   
    sel_cHL(iyear)=sel_cHL_block3;
	sel_cTW(iyear)=sel_cTW_vec;
    sel_rHB(iyear)=sel_rHB_block3;
    sel_rGN(iyear)=sel_rGN_block3;
	sel_sCT(iyear)=sel_sCT_vec;
   }  
  //DISCARDS
	for (iyear=styr; iyear<=endyr; iyear++)
	{
		for (iage=1; iage<=nages; iage++)
		{
		 sel_D_cHL(iyear,iage)=max(column(sel_cHL, iage)); 
		 sel_D_rHB(iyear,iage)=max(column(sel_rHB, iage));
		}
	}
	sel_D_rGN=sel_D_rHB;
}

void model_parameters::get_mortality(void)
{
  Fsum.initialize();
  Fapex.initialize();
  F.initialize();
  //initialization F is avg from first 3 yrs of observed landings
  log_F_dev_init_L_cHL=sum(log_dev_F_L_cHL(styr_L_cHL,(styr_L_cHL+2)))/3.0; //LINETAG: L_cHL  
  log_F_dev_init_L_cTW=sum(log_dev_F_L_cTW(styr_L_cTW,(styr_L_cTW+2)))/3.0; //LINETAG: L_cTW
  log_F_dev_init_L_rHB=sum(log_dev_F_L_rHB(styr_L_rHB,(styr_L_rHB+2)))/3.0; //LINETAG: L_rHB         
  log_F_dev_init_L_rGN=sum(log_dev_F_L_rGN(styr_L_rGN,(styr_L_rGN+2)))/3.0; //LINETAG: L_rGN         
  for (iyear=styr; iyear<=endyr; iyear++) 
  {
    if(iyear>=styr_L_cHL & iyear<=endyr_L_cHL) //spans full time series  //LINETAG: L_cHL
		{F_cHL_out(iyear)=mfexp(log_avg_F_L_cHL+log_dev_F_L_cHL(iyear));} //LINETAG: L_cHL     
    F_cHL(iyear)=sel_cHL(iyear)*F_cHL_out(iyear); //LINETAG: L_cHL
    Fsum(iyear)+=F_cHL_out(iyear); //LINETAG: L_cHL
    if(iyear>=styr_L_cTW & iyear<=endyr_L_cTW) //spans full time series  //LINETAG: L_cTW
		{F_cTW_out(iyear)=mfexp(log_avg_F_L_cTW+log_dev_F_L_cTW(iyear));} //LINETAG: L_cTW
	F_cTW(iyear)=sel_cTW(iyear)*F_cTW_out(iyear); //LINETAG: L_cTW
	Fsum(iyear)+=F_cTW_out(iyear); //LINETAG: L_cTW
    if(iyear>=styr_L_rHB & iyear<=endyr_L_rHB) //spans full time series //LINETAG: L_rHB
		{F_rHB_out(iyear)=mfexp(log_avg_F_L_rHB+log_dev_F_L_rHB(iyear));} //LINETAG: L_rHB     
    F_rHB(iyear)=sel_rHB(iyear)*F_rHB_out(iyear); //LINETAG: L_rHB
    Fsum(iyear)+=F_rHB_out(iyear); //LINETAG: L_rHB
    if(iyear>=styr_L_rGN & iyear<=endyr_L_rGN) //spans full time series //LINETAG: L_rGN
		{F_rGN_out(iyear)=mfexp(log_avg_F_L_rGN+log_dev_F_L_rGN(iyear));} //LINETAG: L_rGN    
	F_rGN(iyear)=sel_rGN(iyear)*F_rGN_out(iyear); //LINETAG: L_rGN 
    Fsum(iyear)+=F_rGN_out(iyear); //LINETAG: L_rGN
    log_F_avgdev_D_cHL=sum(log_dev_F_D_cHL(styr_D_cHL,endyr_D_cHL))/(endyr_D_cHL-styr_D_cHL+1.0);
    log_F_avgdev_D_rHB=sum(log_dev_F_D_rHB(styr_D_rHB,endyr_D_rHB))/(endyr_D_rHB-styr_D_rHB+1.0);   
    log_F_avgdev_D_rGN=sum(log_dev_F_D_rGN(styr_D_rGN,endyr_D_rGN))/(endyr_D_rGN-styr_D_rGN+1.0);	
    if(iyear>=styr_D_cHL & iyear<=endyr_D_cHL) //LINETAG: D_cHL
		{F_D_cHL_out(iyear)=mfexp(log_avg_F_D_cHL+log_dev_F_D_cHL(iyear));} //LINETAG: D_cHL    
    if(iyear>endyr_selex_phase1 & iyear<styr_D_cHL) //LINETAG: D_cHL
		{F_D_cHL_out(iyear)=mfexp(log_avg_F_D_cHL+log_F_avgdev_D_cHL);} 	 //LINETAG: D_cHL
    F_D_cHL(iyear)=sel_D_cHL(iyear)*F_D_cHL_out(iyear); //LINETAG: D_cHL
    Fsum(iyear)+=F_D_cHL_out(iyear); //LINETAG: D_cHL
    if(iyear>=styr_D_rHB & iyear<=endyr_D_rHB) //LINETAG: D_rHB
		{F_D_rHB_out(iyear)=mfexp(log_avg_F_D_rHB+log_dev_F_D_rHB(iyear));} //LINETAG: D_rHB    
    if(iyear>endyr_selex_phase1 & iyear<styr_D_rHB) //LINETAG: D_rHB
      {F_D_rHB_out(iyear)=mfexp(log_avg_F_D_rHB+log_F_avgdev_D_rHB);} //LINETAG: D_rHB 	
    F_D_rHB(iyear)=sel_D_rHB(iyear)*F_D_rHB_out(iyear); //LINETAG: D_rHB
    Fsum(iyear)+=F_D_rHB_out(iyear); //LINETAG: D_rHB
    if(iyear>=styr_D_rGN & iyear<=endyr_D_rGN) //LINETAG: D_rGN
		{F_D_rGN_out(iyear)=mfexp(log_avg_F_D_rGN+log_dev_F_D_rGN(iyear));}   //LINETAG: D_rGN
	F_D_rGN(iyear)=sel_D_rGN(iyear)*F_D_rGN_out(iyear);  //LINETAG: D_rGN
	Fsum(iyear)+=F_D_rGN_out(iyear); //LINETAG: D_rGN
    //Total F at age
    F(iyear)=F_cHL(iyear);  //first in additive series (NO +=)
    F(iyear)+=F_cTW(iyear);
    F(iyear)+=F_rHB(iyear);
    F(iyear)+=F_rGN(iyear);
    F(iyear)+=F_D_cHL(iyear);
    F(iyear)+=F_D_rHB(iyear);
    F(iyear)+=F_D_rGN(iyear);
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
  S0=spr_F0(styr)*R0;
  R_virgin=SR_eq_func(R0, steep, spr_F0(styr), spr_F0(styr), BiasCor, SR_switch);
  B0=bpr_F0(styr)*R_virgin;   
  B0_q_DD=R_virgin*sum(elem_prod(N_bpr_F0(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages))); 
  F_initial=sel_cHL(styr)*mfexp(log_avg_F_L_cHL+log_F_dev_init_L_cHL)+
            sel_cTW(styr)*mfexp(log_avg_F_L_cTW+log_F_dev_init_L_cTW)+
            sel_rHB(styr)*mfexp(log_avg_F_L_rHB+log_F_dev_init_L_rHB)+
            sel_rGN(styr)*mfexp(log_avg_F_L_rGN+log_F_dev_init_L_rGN);			
  Z_initial=M+F_initial;
  N_spr_initial(1)=1.0*mfexp(-1.0*Z_initial(1)*spawn_time_frac); //at peak spawning time;
  for (iage=2; iage<=nages; iage++)
    {
      N_spr_initial(iage)=N_spr_initial(iage-1)*
                   mfexp(-1.0*(Z_initial(iage-1)*(1.0-spawn_time_frac) + Z_initial(iage)*spawn_time_frac)); 
    }
  N_spr_initial(nages)=N_spr_initial(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
  spr_initial=sum(elem_prod(N_spr_initial,reprod(styr)));
  if (styr==styr_rec_dev) {R1=SR_eq_func(R0, steep, spr_F0(styr), spr_initial, 1.0, SR_switch);} //without bias correction (deviation added later)
  else {R1=SR_eq_func(R0, steep, spr_F0(styr), spr_initial, BiasCor, SR_switch);} //with bias correction
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
  SSB(styr)=sum(elem_prod(N_spawn(styr),reprod(styr)));
  B_q_DD(styr)=sum(elem_prod(N(styr)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
  for (iyear=styr; iyear<endyr; iyear++)
  {
    if(iyear<(styr_rec_dev-1)||iyear>(endyr_rec_dev-1)) //recruitment follows S-R curve (with bias correction) exactly
    {
        N(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0(iyear), SSB(iyear),SR_switch);
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); //plus group
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod(iyear+1)));
		B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));       
    }
    else   //recruitment follows S-R curve with lognormal deviation
    {
        N(iyear+1,1)=SR_func(R0, steep, spr_F0(iyear), SSB(iyear),SR_switch)*mfexp(log_dev_rec(iyear+1));
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); //plus group
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod(iyear+1)));
        B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
    }
  }
  //last year (projection) has no recruitment variability
  N(endyr+1,1)=BiasCor*SR_func(R0, steep, spr_F0(endyr), SSB(endyr),SR_switch);
  N(endyr+1)(2,nages)=++elem_prod(N(endyr)(1,nages-1),(mfexp(-1.*Z(endyr)(1,nages-1))));
  N(endyr+1,nages)+=N(endyr,nages)*mfexp(-1.*Z(endyr,nages)); //plus group
}

void model_parameters::get_landings_numbers(void)
{
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      L_cHL_num(iyear,iage)=N(iyear,iage)*F_cHL(iyear,iage)*(1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_cTW_num(iyear,iage)=N(iyear,iage)*F_cTW(iyear,iage)*(1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);		
      L_rHB_num(iyear,iage)=N(iyear,iage)*F_rHB(iyear,iage)*(1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_rGN_num(iyear,iage)=N(iyear,iage)*F_rGN(iyear,iage)*(1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);        
    }          
    pred_L_cHL_knum(iyear)=sum(L_cHL_num(iyear))/1000.0;
    pred_L_cTW_knum(iyear)=sum(L_cTW_num(iyear))/1000.0;
    pred_L_rHB_knum(iyear)=sum(L_rHB_num(iyear))/1000.0;
    pred_L_rGN_knum(iyear)=sum(L_rGN_num(iyear))/1000.0;
  }
}

void model_parameters::get_landings_wgt(void)
{
  for (iyear=styr; iyear<=endyr; iyear++)
  {    
    L_cHL_klb(iyear)=elem_prod(L_cHL_num(iyear),wholewgt_cHL_klb(iyear));     //in 1000 lb whole weight
	L_cTW_klb(iyear)=elem_prod(L_cTW_num(iyear),wholewgt_cTW_klb(iyear));     //in 1000 lb whole weight
    L_rHB_klb(iyear)=elem_prod(L_rHB_num(iyear),wholewgt_rHB_klb(iyear));     //in 1000 lb whole weight
    L_rGN_klb(iyear)=elem_prod(L_rGN_num(iyear),wholewgt_rGN_klb(iyear));     //in 1000 lb whole weight
    pred_L_cHL_klb(iyear)=sum(L_cHL_klb(iyear));
    pred_L_cTW_klb(iyear)=sum(L_cTW_klb(iyear));
    pred_L_rHB_klb(iyear)=sum(L_rHB_klb(iyear));
    pred_L_rGN_klb(iyear)=sum(L_rGN_klb(iyear));    
  }
}

void model_parameters::get_dead_discards(void)
{
  //dead discards at age (number fish) 
  //for (iyear=styr_D_cHL; iyear<=endyr_D_cHL; iyear++) //LINETAG: D_cHL
  for (iyear=styr; iyear<=endyr; iyear++) //LINETAG: D_cHL	  
  { //LINETAG: D_cHL
    for (iage=1; iage<=nages; iage++) //LINETAG: D_cHL
    { //LINETAG: D_cHL
      D_cHL_num(iyear,iage)=N(iyear,iage)*F_D_cHL(iyear,iage)*(1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage); //LINETAG: D_cHL
    } //LINETAG: D_cHL
    pred_D_cHL_knum(iyear)=sum(D_cHL_num(iyear))/1000.0;            //pred annual dead discards in 1000s (for matching data) //LINETAG: D_cHL
    pred_D_cHL_klb(iyear)=sum(elem_prod(D_cHL_num(iyear),wholewgt_D_cHL_klb(iyear))); //annual dead discards in 1000 lb whole (for output only) //LINETAG: D_cHL
  } //LINETAG: D_cHL
  //for (iyear=styr_D_rHB; iyear<=endyr_D_rHB; iyear++) //LINETAG: D_rHB
  for (iyear=styr; iyear<=endyr; iyear++) //LINETAG: D_rHB
  { //LINETAG: D_rHB
    for (iage=1; iage<=nages; iage++) //LINETAG: D_rHB
    { //LINETAG: D_rHB
      D_rHB_num(iyear,iage)=N(iyear,iage)*F_D_rHB(iyear,iage)*(1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage); //LINETAG: D_rHB
    }
    pred_D_rHB_knum(iyear)=sum(D_rHB_num(iyear))/1000.0;            //pred annual dead discards in 1000s (for matching data) //LINETAG: D_rHB
    pred_D_rHB_klb(iyear)=sum(elem_prod(D_rHB_num(iyear),wholewgt_D_rHB_klb(iyear))); //annual dead discards in 1000 lb whole (for output only)
  }
  //for (iyear=styr_D_rGN; iyear<=endyr_D_rGN; iyear++) //LINETAG: D_rGN
  for (iyear=styr; iyear<=endyr; iyear++) //LINETAG: D_rGN
  { //LINETAG: D_rGN
    for (iage=1; iage<=nages; iage++) //LINETAG: D_rGN
    { //LINETAG: D_rGN
      D_rGN_num(iyear,iage)=N(iyear,iage)*F_D_rGN(iyear,iage)*(1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage); //LINETAG: D_rGN
    } //LINETAG: D_rGN
    pred_D_rGN_knum(iyear)=sum(D_rGN_num(iyear))/1000.0;            //pred annual dead discards in 1000s (for matching data) //LINETAG: D_rGN
    pred_D_rGN_klb(iyear)=sum(elem_prod(D_rGN_num(iyear),wholewgt_D_rGN_klb(iyear))); //annual dead discards in 1000 lb whole (for output only) //LINETAG: D_rGN
  } //LINETAG: D_rGN
}

void model_parameters::get_catchability_fcns(void)
{
 //Get rate increase if estimated, otherwise fixed above
  if (set_q_rate_phase>0.0)
  {
      for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++) //LINETAG: cpue_rHB
      {   if (iyear>styr_cpue_rHB & iyear <=2003)  //LINETAG: cpue_rHB
          {//q_rate_fcn_cpue_rHB(iyear)=(1.0+q_rate)*q_rate_fcn_cpue_rHB(iyear-1); //compound //LINETAG: cpue_rHB
             q_rate_fcn_cpue_rHB(iyear)=(1.0+(iyear-styr_cpue_rHB)*q_rate)*q_rate_fcn_cpue_rHB(styr_cpue_rHB);  //linear //LINETAG: cpue_rHB
          } //LINETAG: cpue_rHB
          if (iyear>2003) {q_rate_fcn_cpue_rHB(iyear)=q_rate_fcn_cpue_rHB(iyear-1);}  //LINETAG: cpue_rHB
      }    //LINETAG: cpue_rHB
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
  q_cpue_rHB(styr_cpue_rHB)=mfexp(log_q_cpue_rHB);  //LINETAG: cpue_rHB
  for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++) //LINETAG: cpue_rHB
  {    //LINETAG: cpue_rHB
      N_rHB(iyear)=elem_prod(N_mdyr(iyear),sel_rHB(iyear));  //LINETAG: cpue_rHB
      pred_cpue_rHB(iyear)=q_cpue_rHB(iyear)*q_rate_fcn_cpue_rHB(iyear)*q_DD_fcn(iyear)*sum(N_rHB(iyear)); //LINETAG: cpue_rHB
      if (iyear<endyr_cpue_rHB){q_cpue_rHB(iyear+1)=q_cpue_rHB(iyear)*mfexp(q_RW_log_dev_cpue_rHB(iyear));} //LINETAG: cpue_rHB
  } //LINETAG: cpue_rHB
  q_cpue_sCT(styr_cpue_sCT)=mfexp(log_q_cpue_sCT);  //LINETAG: cpue_sCT
  for (iyear=styr_cpue_sCT; iyear<=endyr_cpue_sCT; iyear++)  //LINETAG: cpue_sCT
  {     //LINETAG: cpue_sCT
    N_sCT(iyear)=elem_prod(N_mdyr(iyear),sel_sCT(iyear));   //LINETAG: cpue_sCT
	pred_cpue_sCT(iyear)=q_cpue_sCT(iyear)*q_DD_fcn(iyear)*sum(N_sCT(iyear));  //LINETAG: cpue_sCT
    if (iyear<endyr_cpue_sCT){q_cpue_sCT(iyear+1)=q_cpue_sCT(iyear)*mfexp(q_RW_log_dev_sCT(iyear));}  //LINETAG: cpue_sCT
  }   //LINETAG: cpue_sCT
}

void model_parameters::get_length_comps(void)
{
  for (iyear=1;iyear<=nyr_lenc_pool_cTW;iyear++)
    {pred_lenc_cTW_yr(iyear)=(L_cTW_num(yrs_lenc_pool_cTW(iyear))*lenprob_lenc_cTW)/sum(L_cTW_num(yrs_lenc_pool_cTW(iyear)));}
  pred_lenc_cTW.initialize();
  for (iyear=1;iyear<=nyr_lenc_pool_cTW;iyear++)
    {pred_lenc_cTW(1) += nfish_lenc_pool_cTW(iyear) * pred_lenc_cTW_yr(iyear);}
  pred_lenc_cTW(1)=pred_lenc_cTW(1)/sum(nfish_lenc_pool_cTW);
}

void model_parameters::get_age_comps(void)
{
  for (iyear=1;iyear<=nyr_agec_cHL;iyear++)  //LINETAG: agec_cHL
  {  //LINETAG: agec_cHL
    ErrorFree_agec_cHL(iyear)=L_cHL_num(yrs_agec_cHL(iyear))/sum(L_cHL_num(yrs_agec_cHL(iyear)));    //LINETAG: agec_cHL
    pred_agec_cHL_allages(iyear)=age_error*(ErrorFree_agec_cHL(iyear)/sum(ErrorFree_agec_cHL(iyear)));     //LINETAG: agec_cHL
    for (iage=1; iage<=nages_agec; iage++) {pred_agec_cHL(iyear,iage)=pred_agec_cHL_allages(iyear,iage);}   //LINETAG: agec_cHL
    //for (iage=(nages_agec+1); iage<=nages; iage++) {pred_agec_cHL(iyear,nages_agec)+=pred_agec_cHL_allages(iyear,iage);} //plus group                               //LINETAG: agec_cHL
  }  //LINETAG: agec_cHL
 for (iyear=1;iyear<=nyr_agec_rHB;iyear++)  //LINETAG: agec_rHB
  {  //LINETAG: agec_rHB
    ErrorFree_agec_rHB(iyear)=L_rHB_num(yrs_agec_rHB(iyear))/sum(L_rHB_num(yrs_agec_rHB(iyear)));  //LINETAG: agec_rHB
    pred_agec_rHB_allages(iyear)=age_error*ErrorFree_agec_rHB(iyear);   //LINETAG: agec_rHB
    for (iage=1; iage<=nages_agec; iage++) {pred_agec_rHB(iyear,iage)=pred_agec_rHB_allages(iyear,iage);}   //LINETAG: agec_rHB
    //for (iage=(nages_agec+1); iage<=nages; iage++) {pred_agec_rHB(iyear,nages_agec_rHB)+=pred_agec_rHB_allages(iyear,iage);} //plus group       //LINETAG: agec_rHB                   
  }  //LINETAG: agec_rHB
 for (iyear=1;iyear<=nyr_agec_sCT;iyear++)  //LINETAG: agec_sCT
  {  //LINETAG: agec_sCT
    ErrorFree_agec_sCT(iyear)=N_sCT(yrs_agec_sCT(iyear))/sum(N_sCT(yrs_agec_sCT(iyear)));  //LINETAG: agec_sCT
    pred_agec_sCT_allages(iyear)=age_error*ErrorFree_agec_sCT(iyear);   //LINETAG: agec_sCT
    for (iage=1; iage<=nages_agec; iage++) {pred_agec_sCT(iyear,iage)=pred_agec_sCT_allages(iyear,iage);}   //LINETAG: agec_sCT
    //for (iage=(nages_agec+1); iage<=nages; iage++) {pred_agec_sCT(iyear,nages_agec)+=pred_agec_sCT_allages(iyear,iage);} //plus group //LINETAG: agec_sCT
  } //LINETAG: agec_sCT
}

void model_parameters::get_weighted_current(void)
{
  F_temp_sum=0.0;
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cHL+sum(log_dev_F_L_cHL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  //LINETAG: L_cHL
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rHB+sum(log_dev_F_L_rHB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  //LINETAG: L_rHB
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rGN+sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  //LINETAG: L_rGN
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_D_cHL+sum(log_dev_F_D_cHL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  //LINETAG: D_cHL
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_D_rHB+sum(log_dev_F_D_rHB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  //LINETAG: D_rHB
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_D_rGN+sum(log_dev_F_D_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  //LINETAG: D_rGN
  F_prop_L_cHL=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cHL+sum(log_dev_F_L_cHL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;  //LINETAG: L_cHL
  F_prop_L_rHB=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rHB+sum(log_dev_F_L_rHB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;  //LINETAG: L_rHB
  F_prop_L_rGN=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rGN+sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;  //LINETAG: L_rGN
  F_prop_D_cHL=mfexp((selpar_n_yrs_wgted*log_avg_F_D_cHL+sum(log_dev_F_D_cHL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;  //LINETAG: D_cHL
  F_prop_D_rHB=mfexp((selpar_n_yrs_wgted*log_avg_F_D_rHB+sum(log_dev_F_D_rHB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;  //LINETAG: D_rHB
  F_prop_D_rGN=mfexp((selpar_n_yrs_wgted*log_avg_F_D_rGN+sum(log_dev_F_D_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;  //LINETAG: D_rGN
  log_F_dev_end_L_cHL=sum(log_dev_F_L_cHL((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  //LINETAG: L_cHL
  log_F_dev_end_L_rHB=sum(log_dev_F_L_rHB((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  //LINETAG: L_rHB  
  log_F_dev_end_L_rGN=sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  //LINETAG: L_rGN  
  log_F_dev_end_D_cHL=sum(log_dev_F_D_cHL((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  //LINETAG: D_cHL
  log_F_dev_end_D_rHB=sum(log_dev_F_D_rHB((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  //LINETAG: D_rHB
  log_F_dev_end_D_rGN=sum(log_dev_F_D_rGN((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  //LINETAG: D_rGN
  F_end_L=
  sel_cHL(endyr)*mfexp(log_avg_F_L_cHL+log_F_dev_end_L_cHL)+  //LINETAG: L_cHL
  sel_rHB(endyr)*mfexp(log_avg_F_L_rHB+log_F_dev_end_L_rHB)+  //LINETAG: L_rHB
  sel_rGN(endyr)*mfexp(log_avg_F_L_rGN+log_F_dev_end_L_rGN)   //LINETAG: L_rGN
  ;   
  F_end_D=
  sel_D_cHL(endyr)*mfexp(log_avg_F_D_cHL+log_F_dev_end_D_cHL)+  //LINETAG: D_cHL
  sel_D_rHB(endyr)*mfexp(log_avg_F_D_rHB+log_F_dev_end_D_rHB)+  //LINETAG: D_rHB
  sel_D_rGN(endyr)*mfexp(log_avg_F_D_rGN+log_F_dev_end_D_rGN)  //LINETAG: D_rGN
  ;   
  F_end=F_end_L+F_end_D;
  F_end_apex=max(F_end);
  sel_wgted_tot=F_end/F_end_apex;
  sel_wgted_L=elem_prod(sel_wgted_tot, elem_div(F_end_L,F_end));
  sel_wgted_D=elem_prod(sel_wgted_tot, elem_div(F_end_D,F_end));
  wgt_wgted_L_denom=
  F_prop_L_cHL+  //LINETAG: L_cHL
  F_prop_L_rHB+  //LINETAG: L_rHB
  F_prop_L_rGN   //LINETAG: L_rGN
  ;
  wgt_wgted_L_klb=
  F_prop_L_cHL/wgt_wgted_L_denom*wholewgt_cHL_klb(endyr)+  //LINETAG: L_cHL 
  F_prop_L_rHB/wgt_wgted_L_denom*wholewgt_rHB_klb(endyr)+  //LINETAG: L_rHB
  F_prop_L_rGN/wgt_wgted_L_denom*wholewgt_rGN_klb(endyr)   //LINETAG: L_rGN
  ;                          
  wgt_wgted_D_denom=
  F_prop_D_cHL+  //LINETAG: D_cHL
  F_prop_D_rHB+  //LINETAG: D_rHB
  F_prop_D_rGN  //LINETAG: D_rGN
  ;
  wgt_wgted_D_klb=
  F_prop_D_cHL/wgt_wgted_D_denom*wholewgt_D_cHL_klb(endyr)+  //LINETAG: D_cHL 
  F_prop_D_rHB/wgt_wgted_D_denom*wholewgt_D_rHB_klb(endyr)+  //LINETAG: D_rHB
  F_prop_D_rGN/wgt_wgted_D_denom*wholewgt_D_rGN_klb(endyr)  //LINETAG: D_rGN
  ;                
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
    spr_msy(ff)=sum(elem_prod(N_age_msy_spawn,reprod(endyr)));
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
    SSB_eq(ff)=sum(elem_prod(N_age_msy_spawn,reprod(endyr)));
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
    spr_static(iyear)=sum(elem_prod(N_age_spr_spawn,reprod(iyear)))/spr_F0(iyear);
  }
  //compute SSB/R and YPR as functions of F
  for(ff=1; ff<=n_iter_spr; ff++)
  {
    //uses fishery-weighted F's, same as in MSY calculations
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
    spr_spr(ff)=sum(elem_prod(N_age_spr_spawn,reprod(endyr)));
    L_spr(ff)=0.0;
    for (iage=1; iage<=nages; iage++)
    {
      L_age_spr(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));
      L_spr(ff)+=L_age_spr(iage)*wgt_wgted_L_klb(iage)*1000.0; //in lb whole wgt
    }   
  }
  spr_ratio=spr_spr/spr_F0(endyr);
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
  SSB_F30_out=sum(elem_prod(N_age_spr_spawn,reprod(endyr)));
  B_F30_out=sum(elem_prod(N_age_spr,wgt_mt));
  L_F30_klb_out=sum(elem_prod(L_age_F30,wgt_wgted_L_klb)); //in whole weight
  L_F30_knum_out=sum(L_age_F30)/1000.0;  
  D_F30_klb_out=sum(elem_prod(D_age_F30,wgt_wgted_D_klb)); //in whole weight   
  D_F30_knum_out=sum(D_age_F30)/1000.0;    
}

void model_parameters::get_miscellaneous_stuff(void)
{
  if(var_rec_dev>0.0)
   {sigma_rec_dev=sqrt(var_rec_dev);} //sample SD of predicted residuals (may not equal rec_sigma)  
   else{sigma_rec_dev=0.0;}
  len_cv=elem_div(len_sd,meanlen_TL);
  //compute total landings- and discards-at-age in 1000 fish and klb whole weight
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
		L_total_klb_yr(iyear)=
		pred_L_cTW_klb(iyear)+
		pred_L_cHL_klb(iyear)+
		pred_L_rHB_klb(iyear)+
		pred_L_rGN_klb(iyear)
		;
		L_total_knum_yr(iyear)=
		pred_L_cTW_knum(iyear)+
		pred_L_cHL_knum(iyear)+
		pred_L_rHB_knum(iyear)+
		pred_L_rGN_knum(iyear)
		;
        B(iyear)=elem_prod(N(iyear),wgt_mt);
        totN(iyear)=sum(N(iyear));
        totB(iyear)=sum(B(iyear));   
        //if (iyear>=styr_D_cHL && iyear<=endyr_D_cHL) //LINETAG: D_cHL 
	    if (iyear>endyr_selex_phase1 && iyear<=endyr) //LINETAG: D_cHL 			
        { //LINETAG: D_cHL 
         D_total_knum_yr(iyear)+=pred_D_cHL_knum(iyear); //LINETAG: D_cHL 
         D_total_klb_yr(iyear)+=pred_D_cHL_klb(iyear); //LINETAG: D_cHL 
         D_cHL_klb(iyear)=elem_prod(D_cHL_num(iyear),wholewgt_D_cHL_klb(iyear));     //in 1000 lb  //LINETAG: D_cHL 
        } //LINETAG: D_cHL 
        //if (iyear>=styr_D_rHB && iyear<=endyr_D_rHB) //LINETAG: D_rHB
		if (iyear>endyr_selex_phase1 && iyear<=endyr) //LINETAG: D_rHB
        { //LINETAG: D_cHlB
         D_total_knum_yr(iyear)+=pred_D_rHB_knum(iyear); //LINETAG: D_rHB
         D_total_klb_yr(iyear)+=pred_D_rHB_klb(iyear); //LINETAG: D_rHB
         D_rHB_klb(iyear)=elem_prod(D_rHB_num(iyear),wholewgt_D_rHB_klb(iyear));     //in 1000 lb  //LINETAG: D_rHB
        }     //LINETAG: D_cHlB
        //if (iyear>=styr_D_rGN && iyear<=endyr_D_rGN) //LINETAG: D_rGN
		if (iyear>=styr_D_rGN && iyear<=endyr) //LINETAG: D_rGN	
        { //LINETAG: D_rGN
         D_total_knum_yr(iyear)+=pred_D_rGN_knum(iyear); //LINETAG: D_rGN
         D_total_klb_yr(iyear)+=pred_D_rGN_klb(iyear); //LINETAG: D_rGN
         D_rGN_klb(iyear)=elem_prod(D_rGN_num(iyear),wholewgt_D_rGN_klb(iyear));     //in 1000 lb  //LINETAG: D_rGN
        }           //LINETAG: D_rGN                
  }
  //landings at age in number fish
  L_total_num=
  L_cHL_num+
  L_cTW_num+
  L_rHB_num+
  L_rGN_num
  ;
  //landings at age in klb whole weight
  L_total_klb=
  L_cHL_klb+
  L_cTW_klb+
  L_rHB_klb+
  L_rGN_klb
  ;   
  //discards at age in number fish
  D_total_num=
  D_cHL_num+
  D_rHB_num+
  D_rGN_num
  ;          
  //discards at age in klb whole weight
  D_total_klb=
  D_cHL_klb+
  D_rHB_klb+
  D_rGN_klb
  ;
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
		FD_age_proj=sel_wgted_D*F_proj(iyear);
        Z_proj(iyear)=M+FL_age_proj+FD_age_proj;
        N_spawn_proj(iyear)(1,nages)=elem_prod(N_proj(iyear)(1,nages),(mfexp(-1.*(Z_proj(iyear)(1,nages))*spawn_time_frac))); //peak spawning time
		SSB_proj(iyear)= sum(elem_prod(N_spawn_proj(iyear),reprod(endyr)));
        B_proj(iyear)=sum(elem_prod(N_proj(iyear),wgt_mt)); //uses spawning weight
		for (iage=1; iage<=nages; iage++)
			{L_age_proj(iyear,iage)=N_proj(iyear,iage)*FL_age_proj(iage)*(1.-mfexp(-1.*Z_proj(iyear,iage)))/Z_proj(iyear,iage);
		     D_age_proj(iyear,iage)=N_proj(iyear,iage)*FD_age_proj(iage)*(1.-mfexp(-1.*Z_proj(iyear,iage)))/Z_proj(iyear,iage);
			}          
        L_knum_proj(iyear)=sum(L_age_proj(iyear))/1000.0;
	    D_knum_proj(iyear)=sum(D_age_proj(iyear))/1000.0;
	    L_klb_proj(iyear)=sum(elem_prod(L_age_proj(iyear),wgt_wgted_L_klb));     //in 1000 lb
        D_klb_proj(iyear)=sum(elem_prod(D_age_proj(iyear),wgt_wgted_D_klb));     //in 1000 lb
		if (iyear<endyr_proj) {
			N_proj(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0(endyr), SSB_proj(iyear),SR_switch);
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
  f_cpue_rHB=0.0;
  f_cpue_rHB=lk_lognormal(pred_cpue_rHB, obs_cpue_rHB, obs_cv_cpue_rHB, w_cpue_rHB);
  fval+=f_cpue_rHB;
  fval_data+=f_cpue_rHB;  
  f_cpue_sCT=0.0;
  f_cpue_sCT=lk_lognormal(pred_cpue_sCT, obs_cpue_sCT, obs_cv_cpue_sCT, w_cpue_sCT);
  fval+=f_cpue_sCT;
  fval_data+=f_cpue_sCT;  
  //f_L_cHL in 1000 lb whole wgt  
  f_L_cHL=lk_lognormal(pred_L_cHL_klb(styr_L_cHL,endyr_L_cHL), obs_L_cHL(styr_L_cHL,endyr_L_cHL),
                      obs_cv_L_cHL(styr_L_cHL,endyr_L_cHL), w_L);
  fval+=f_L_cHL;
  fval_data+=f_L_cHL;
  //f_L_cTW in 1000 lb whole wgt  
  f_L_cTW=lk_lognormal(pred_L_cTW_klb(styr_L_cTW,endyr_L_cTW), obs_L_cTW(styr_L_cTW,endyr_L_cTW),
                      obs_cv_L_cTW(styr_L_cTW,endyr_L_cTW), w_L);
  fval+=f_L_cTW;
  fval_data+=f_L_cTW;
  //f_L_rHB in 1000 fish
  f_L_rHB=lk_lognormal(pred_L_rHB_knum(styr_L_rHB,endyr_L_rHB), obs_L_rHB(styr_L_rHB,endyr_L_rHB), 
                      obs_cv_L_rHB(styr_L_rHB,endyr_L_rHB), w_L);
  fval+=f_L_rHB;
  fval_data+=f_L_rHB;  
  //f_L_rGN in 1000 fish
  f_L_rGN=lk_lognormal(pred_L_rGN_knum(styr_L_rGN,endyr_L_rGN), obs_L_rGN(styr_L_rGN,endyr_L_rGN), 
                      obs_cv_L_rGN(styr_L_rGN,endyr_L_rGN), w_L);
  fval+=f_L_rGN;
  fval_data+=f_L_rGN;  
  //f_D_cHL in 1000 fish
  f_D_cHL=lk_lognormal(pred_D_cHL_knum(styr_D_cHL,endyr_D_cHL), obs_D_cHL(styr_D_cHL,endyr_D_cHL), 
                      obs_cv_D_cHL(styr_D_cHL,endyr_D_cHL), w_D);
  fval+=f_D_cHL;
  fval_data+=f_D_cHL;  
  //f_D_rHB in 1000 fish
  f_D_rHB=lk_lognormal(pred_D_rHB_knum(styr_D_rHB,endyr_D_rHB), obs_D_rHB(styr_D_rHB,endyr_D_rHB), 
                      obs_cv_D_rHB(styr_D_rHB,endyr_D_rHB), w_D);
  fval+=f_D_rHB;
  fval_data+=f_D_rHB;  
  //f_D_rGN in 1000 fish
  f_D_rGN=lk_lognormal(pred_D_rGN_knum(styr_D_rGN,endyr_D_rGN), obs_D_rGN(styr_D_rGN,endyr_D_rGN), 
                      obs_cv_D_rGN(styr_D_rGN,endyr_D_rGN), w_D);
  fval+=f_D_rGN;
  fval_data+=f_D_rGN;  
  //f_lenc_cTW
  //f_lenc_cTW=lk_robust_multinomial(nsamp_lenc_cTW, pred_lenc_cTW, obs_lenc_cTW, nyr_lenc_cTW, double(nlenbins), minSS_lenc_cTW, w_lenc_cTW);
  //f_lenc_cTW=lk_logistic_normal(nsamp_lenc_cTW, pred_lenc_cTW, obs_lenc_cTW, nyr_lenc_cTW, double(nlenbins), minSS_lenc_cTW);
  f_lenc_cTW=lk_dirichlet_multinomial(nsamp_lenc_cTW, pred_lenc_cTW, obs_lenc_cTW, nyr_lenc_cTW, double(nlenbins), minSS_lenc_cTW, log_dm_lenc_cTW);
  fval+=f_lenc_cTW;
  fval_data+=f_lenc_cTW;
  //f_agec_cHL
  //f_agec_cHL=lk_robust_multinomial(nsamp_agec_cHL, pred_agec_cHL, obs_agec_cHL, nyr_agec_cHL, double(nages_agec), minSS_agec_cHL, w_agec_cHL);
  //f_agec_cHL=lk_logistic_normal(nsamp_agec_cHL, pred_agec_cHL, obs_agec_cHL, nyr_agec_cHL, double(nages_agec), minSS_agec_cHL);
  f_agec_cHL=lk_dirichlet_multinomial(nsamp_agec_cHL, pred_agec_cHL, obs_agec_cHL, nyr_agec_cHL, double(nages_agec), minSS_agec_cHL, log_dm_agec_cHL);
  fval+=f_agec_cHL;
  fval_data+=f_agec_cHL;
  //f_agec_rHB
  //f_agec_rHB=lk_robust_multinomial(nsamp_agec_rHB, pred_agec_rHB, obs_agec_rHB, nyr_agec_rHB, double(nages_agec), minSS_agec_rHB, w_agec_rHB);
  //f_agec_rHB=lk_logistic_normal(nsamp_agec_rHB, pred_agec_rHB, obs_agec_rHB, nyr_agec_rHB, double(nages_agec), minSS_agec_rHB);
  f_agec_rHB=lk_dirichlet_multinomial(nsamp_agec_rHB, pred_agec_rHB, obs_agec_rHB, nyr_agec_rHB, double(nages_agec), minSS_agec_rHB, log_dm_agec_rHB);
  fval+=f_agec_rHB;
  fval_data+=f_agec_rHB;
  //f_agec_sCT
  //f_agec_sCT=lk_robust_multinomial(nsamp_agec_sCT, pred_agec_sCT, obs_agec_sCT, nyr_agec_sCT, double(nages_agec), minSS_agec_sCT, w_agec_sCT);
  //f_agec_sCT=lk_logistic_normal(nsamp_agec_sCT, pred_agec_sCT, obs_agec_sCT, nyr_agec_sCT, double(nages_agec), minSS_agec_sCT);
  f_agec_sCT=lk_dirichlet_multinomial(nsamp_agec_sCT, pred_agec_sCT, obs_agec_sCT, nyr_agec_sCT, double(nages_agec), minSS_agec_sCT, log_dm_agec_sCT);
  fval+=f_agec_sCT;
  fval_data+=f_agec_sCT;
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
 f_RWq_cpue_rHB=0.0; //LINETAG: RWq_cpue_rHB
 for (iyear=styr_cpue_rHB; iyear<endyr_cpue_rHB; iyear++) //LINETAG: RWq_cpue_rHB
     {f_RWq_cpue_rHB+=square(q_RW_log_dev_cpue_rHB(iyear))/(2.0*set_RWq_var);} //LINETAG: RWq_cpue_rHB
 fval+=f_RWq_cpue_rHB; //LINETAG: RWq_cpue_rHB
  f_priors=0.0; 
  f_priors+=neg_log_prior(len_cv_val,set_len_cv_val(5),set_len_cv_val(6),set_len_cv_val(7));
  f_priors+=neg_log_prior(steep,set_steep(5),set_steep(6),set_steep(7)); 
  f_priors+=neg_log_prior(log_R0,set_log_R0(5),set_log_R0(6),set_log_R0(7)); 
  f_priors+=neg_log_prior(R_autocorr,set_R_autocorr(5),set_R_autocorr(6),set_R_autocorr(7));
  f_priors+=neg_log_prior(rec_sigma,set_rec_sigma(5),set_rec_sigma(6),set_rec_sigma(7));
  f_priors+=neg_log_prior(A50_sel_cHL2,set_A50_sel_cHL2(5), set_A50_sel_cHL2(6), set_A50_sel_cHL2(7));
  f_priors+=neg_log_prior(slope_sel_cHL2,set_slope_sel_cHL2(5), set_slope_sel_cHL2(6), set_slope_sel_cHL2(7));
  f_priors+=neg_log_prior(A50_sel_cHL3,set_A50_sel_cHL3(5), set_A50_sel_cHL3(6), set_A50_sel_cHL3(7));
  f_priors+=neg_log_prior(slope_sel_cHL3,set_slope_sel_cHL3(5), set_slope_sel_cHL3(6), set_slope_sel_cHL3(7));
  f_priors+=neg_log_prior(A50_sel_cTW,set_A50_sel_cTW(5), set_A50_sel_cTW(6), set_A50_sel_cTW(7));
  f_priors+=neg_log_prior(slope_sel_cTW,set_slope_sel_cTW(5), set_slope_sel_cTW(6), set_slope_sel_cTW(7));
  f_priors+=neg_log_prior(A50_sel_rHB1,set_A50_sel_rHB1(5), set_A50_sel_rHB1(6), set_A50_sel_rHB1(7));
  f_priors+=neg_log_prior(slope_sel_rHB1,set_slope_sel_rHB1(5), set_slope_sel_rHB1(6), set_slope_sel_rHB1(7));
  f_priors+=neg_log_prior(A50_sel_rHB2,set_A50_sel_rHB2(5), set_A50_sel_rHB2(6), set_A50_sel_rHB2(7));
  f_priors+=neg_log_prior(slope_sel_rHB2,set_slope_sel_rHB2(5), set_slope_sel_rHB2(6), set_slope_sel_rHB2(7));
  f_priors+=neg_log_prior(A50_sel_rHB3,set_A50_sel_rHB3(5), set_A50_sel_rHB3(6), set_A50_sel_rHB3(7));
  f_priors+=neg_log_prior(slope_sel_rHB3,set_slope_sel_rHB3(5), set_slope_sel_rHB3(6), set_slope_sel_rHB3(7));
  f_priors+=neg_log_prior(A50_sel_sCT,set_A50_sel_sCT(5), set_A50_sel_sCT(6), set_A50_sel_sCT(7));
  f_priors+=neg_log_prior(slope_sel_sCT,set_slope_sel_sCT(5), set_slope_sel_sCT(6), set_slope_sel_sCT(7));
  f_priors+=neg_log_prior(log_q_cpue_rHB,set_log_q_cpue_rHB(5),set_log_q_cpue_rHB(6),set_log_q_cpue_rHB(7));
  f_priors+=neg_log_prior(log_q_cpue_sCT,set_log_q_cpue_sCT(5),set_log_q_cpue_sCT(6),set_log_q_cpue_sCT(7));
  f_priors+=neg_log_prior(log_dm_lenc_cTW,set_log_dm_lenc_cTW(5),set_log_dm_lenc_cTW(6),set_log_dm_lenc_cTW(7));
  f_priors+=neg_log_prior(log_dm_agec_cHL,set_log_dm_agec_cHL(5),set_log_dm_agec_cHL(6),set_log_dm_agec_cHL(7));
  f_priors+=neg_log_prior(log_dm_agec_rHB,set_log_dm_agec_rHB(5),set_log_dm_agec_rHB(6),set_log_dm_agec_rHB(7));
  f_priors+=neg_log_prior(log_dm_agec_sCT,set_log_dm_agec_sCT(5),set_log_dm_agec_sCT(6),set_log_dm_agec_sCT(7));
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
  LkvalTmp=0.0; 
  dvar_vector nsamp_adjust=nsamp*mfexp(log_dir_par);
  //dvar_vector nsamp_adjust=mfexp(log_dir_par);
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
          if(ab_iq<=0) {
			cout << "Parameter input error: For beta priors, mu*(1-mu)/var must be greater than one. Try decreasing var." << endl;
            exit(0);
		  }
		  if(pred>=0 && pred<=1) LkvalTmp= (1.0-alpha)*log(pred)+(1.0-beta)*log(1.0-pred)-gammln(alpha+beta)+gammln(alpha)+gammln(beta);
          else LkvalTmp=big_number;
          break;
        default: // no such prior pdf currently available
          cout << "Parameter input error: Prior must be either 1(none), 2(lognormal), 3(normal), or 4(beta)." << endl;
          cout << "Presently at least one is " << pdf << endl;
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
      // cout<<"xdum = "<<xdum<<endl;
	  get_weighted_current();
      // cout<<"got weighted"<<endl;
      get_msy();
      // cout<<"got msy"<<endl;
      get_per_recruit_stuff();
      // cout<<"got per recruit"<<endl;  
      get_miscellaneous_stuff();
      // cout<<"got misc stuff"<<endl;
	  get_projection();
	  // cout<<"got projection"<<endl;
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
      cout <<"F status (FdF_msy_end_mean)="<<FdF_msy_end_mean<<endl;
      cout <<"Pop status (SdSSB_msy_end/smsy2msstM)="<<SdSSB_msy_end/smsy2msstM<<endl;
      cout << "h="<<steep<<"   R0="<<R0<<endl;
      //cout << "xdum " << xdum << endl;
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;  
      report << "TotalLikelihood " << fval << endl;
      report << "N" << endl;
      report << N<<endl;
      report << "F" << endl;
      report << F <<endl;
      report << "prob_belowsizelim_block3" << endl;
	  report<<prob_belowsizelim_block3<<endl;
      sdnr_agec_cHL=sdnr_multinomial(nyr_agec_cHL, agebins_agec, nsamp_agec_cHL, pred_agec_cHL, obs_agec_cHL, w_agec_cHL);  
      sdnr_agec_rHB=sdnr_multinomial(nyr_agec_rHB, agebins_agec, nsamp_agec_rHB, pred_agec_rHB, obs_agec_rHB, w_agec_rHB);  
      sdnr_agec_sCT=sdnr_multinomial(nyr_agec_sCT, agebins_agec, nsamp_agec_sCT, pred_agec_sCT, obs_agec_sCT, w_agec_sCT);  
      sdnr_cpue_rHB=sdnr_lognormal(pred_cpue_rHB, obs_cpue_rHB, obs_cv_cpue_rHB, w_cpue_rHB);
      sdnr_cpue_sCT=sdnr_lognormal(pred_cpue_sCT, obs_cpue_sCT, obs_cv_cpue_sCT, w_cpue_sCT);  
      //#################################################################################################
      //##  Passing parameters to vector for bounds check plotting
      //################################################################################################# 
       Linf_out(8)=Linf; Linf_out(1,7)=set_Linf; 
       K_out(8)=K; K_out(1,7)=set_K;
       t0_out(8)=t0; t0_out(1,7)=set_t0;
       len_cv_val_out(8)=len_cv_val; len_cv_val_out(1,7)=set_len_cv_val;
       log_R0_out(8)=log_R0; log_R0_out(1,7)=set_log_R0;
       M_constant_out(8)=M_constant; M_constant_out(1,7)=set_M_constant;
       steep_out(8)=steep; steep_out(1,7)=set_steep;
       rec_sigma_out(8)=rec_sigma; rec_sigma_out(1,7)=set_rec_sigma;
       R_autocorr_out(8)=R_autocorr; R_autocorr_out(1,7)=set_R_autocorr;
	   log_dm_lenc_cTW_out(8)=log_dm_lenc_cTW; log_dm_lenc_cTW_out(1,7)=set_log_dm_lenc_cTW;
	   log_dm_agec_cHL_out(8)=log_dm_agec_cHL; log_dm_agec_cHL_out(1,7)=set_log_dm_agec_cHL;
	   log_dm_agec_rHB_out(8)=log_dm_agec_rHB; log_dm_agec_rHB_out(1,7)=set_log_dm_agec_rHB;
	   log_dm_agec_sCT_out(8)=log_dm_agec_sCT; log_dm_agec_sCT_out(1,7)=set_log_dm_agec_sCT;
       A50_sel_cHL2_out(8)=A50_sel_cHL2; A50_sel_cHL2_out(1,7)=set_A50_sel_cHL2;
       slope_sel_cHL2_out(8)=slope_sel_cHL2; slope_sel_cHL2_out(1,7)=set_slope_sel_cHL2;
       A50_sel_cHL3_out(8)=A50_sel_cHL3; A50_sel_cHL3_out(1,7)=set_A50_sel_cHL3;
       slope_sel_cHL3_out(8)=slope_sel_cHL3; slope_sel_cHL3_out(1,7)=set_slope_sel_cHL3;
	   A50_sel_cTW_out(8)=A50_sel_cTW; A50_sel_cTW_out(1,7)=set_A50_sel_cTW;
       slope_sel_cTW_out(8)=slope_sel_cTW; slope_sel_cTW_out(1,7)=set_slope_sel_cTW;
       A50_sel_rHB1_out(8)=A50_sel_rHB1; A50_sel_rHB1_out(1,7)=set_A50_sel_rHB1;
       slope_sel_rHB1_out(8)=slope_sel_rHB1; slope_sel_rHB1_out(1,7)=set_slope_sel_rHB1;
       A50_sel_rHB2_out(8)=A50_sel_rHB2; A50_sel_rHB2_out(1,7)=set_A50_sel_rHB2;
       slope_sel_rHB2_out(8)=slope_sel_rHB2; slope_sel_rHB2_out(1,7)=set_slope_sel_rHB2;
       A50_sel_rHB3_out(8)=A50_sel_rHB3; A50_sel_rHB3_out(1,7)=set_A50_sel_rHB3;
       slope_sel_rHB3_out(8)=slope_sel_rHB3; slope_sel_rHB3_out(1,7)=set_slope_sel_rHB3;
       A50_sel_sCT_out(8)=A50_sel_sCT; A50_sel_sCT_out(1,7)=set_A50_sel_sCT;
       slope_sel_sCT_out(8)=slope_sel_sCT; slope_sel_sCT_out(1,7)=set_slope_sel_sCT;
       log_q_cpue_rHB_out(8)=log_q_cpue_rHB; log_q_cpue_rHB_out(1,7)=set_log_q_cpue_rHB;
       log_q_cpue_sCT_out(8)=log_q_cpue_sCT; log_q_cpue_sCT_out(1,7)=set_log_q_cpue_sCT;
       log_avg_F_L_cHL_out(8)=log_avg_F_L_cHL; log_avg_F_L_cHL_out(1,7)=set_log_avg_F_L_cHL;  //LINETAG: L_cHL
	   log_avg_F_L_cTW_out(8)=log_avg_F_L_cTW; log_avg_F_L_cTW_out(1,7)=set_log_avg_F_L_cTW;  //LINETAG: L_cTW
       log_avg_F_L_rHB_out(8)=log_avg_F_L_rHB; log_avg_F_L_rHB_out(1,7)=set_log_avg_F_L_rHB;  //LINETAG: L_rHB
       log_avg_F_L_rGN_out(8)=log_avg_F_L_rGN; log_avg_F_L_rGN_out(1,7)=set_log_avg_F_L_rGN;  //LINETAG: L_rGN       
       log_avg_F_D_cHL_out(8)=log_avg_F_D_cHL; log_avg_F_D_cHL_out(1,7)=set_log_avg_F_D_cHL;  //LINETAG: D_cHL
       log_avg_F_D_rHB_out(8)=log_avg_F_D_rHB; log_avg_F_D_rHB_out(1,7)=set_log_avg_F_D_rHB;  //LINETAG: D_rHB
       log_avg_F_D_rGN_out(8)=log_avg_F_D_rGN; log_avg_F_D_rGN_out(1,7)=set_log_avg_F_D_rGN;  //LINETAG: D_rGN
       log_rec_dev_out(styr_rec_dev, endyr_rec_dev)=log_dev_rec;
       log_F_dev_L_cHL_out(styr_L_cHL,endyr_L_cHL)=log_dev_F_L_cHL;  //LINETAG: L_cHL
       log_F_dev_L_cTW_out(styr_L_cTW,endyr_L_cTW)=log_dev_F_L_cTW;  //LINETAG: L_cTW
       log_F_dev_L_rHB_out(styr_L_rHB,endyr_L_rHB)=log_dev_F_L_rHB;  //LINETAG: L_rHB
       log_F_dev_L_rGN_out(styr_L_rGN,endyr_L_rGN)=log_dev_F_L_rGN;  //LINETAG: L_rGN
       log_F_dev_D_cHL_out(styr_D_cHL,endyr_D_cHL)=log_dev_F_D_cHL;  //LINETAG: D_cHL
       log_F_dev_D_rHB_out(styr_D_rHB,endyr_D_rHB)=log_dev_F_D_rHB;  //LINETAG: D_rHB
       log_F_dev_D_rGN_out(styr_D_rGN,endyr_D_rGN)=log_dev_F_D_rGN;  //LINETAG: D_rGN
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
