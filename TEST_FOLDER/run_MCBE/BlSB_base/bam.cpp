#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
  #include "admodel.h"          // Include AD class definitions
  #include "admb2r.cpp"    // Include S-compatible output functions (needs preceding)
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
  endyr_rec_phase1.allocate("endyr_rec_phase1");
  endyr_rec_phase2.allocate("endyr_rec_phase2");
  endyr_period1.allocate("endyr_period1");
  endyr_period2.allocate("endyr_period2");
  endyr_recr_period3.allocate("endyr_recr_period3");
  endyr_period4.allocate("endyr_period4");
  styr_cGN_closed.allocate("styr_cGN_closed");
  sizelim1.allocate("sizelim1");
  sizelim2.allocate("sizelim2");
  sizelim3.allocate("sizelim3");
  sizelim4.allocate("sizelim4");
  sizelim5.allocate("sizelim5");
  limit_disc.allocate("limit_disc");
  nages.allocate("nages");
  agebins.allocate(1,nages,"agebins");
   nyrs=endyr-styr+1.;
   //nyrs_rec=endyr-styr_rec_dev+1.;
   nyrs_rec=endyr_rec_phase2-styr_rec_dev+1.;
  nlenbins.allocate("nlenbins");
  nlenbins_plus.allocate("nlenbins_plus");
  lenbins_width.allocate("lenbins_width");
  lenbins.allocate(1,nlenbins,"lenbins");
  lenbins_plus.allocate(1,nlenbins_plus,"lenbins_plus");
   nlenbins_all=nlenbins+nlenbins_plus;   
  max_F_spr_msy.allocate("max_F_spr_msy");
  n_iter_spr.allocate("n_iter_spr");
  n_iter_msy.allocate("n_iter_msy");
		n_iter_msy=n_iter_spr; 
  selpar_n_yrs_wgted.allocate("selpar_n_yrs_wgted");
  set_BiasCor.allocate("set_BiasCor");
  BiasCor_exclude_yrs.allocate("BiasCor_exclude_yrs");
  styr_cpue_sBT.allocate("styr_cpue_sBT");
  endyr_cpue_sBT.allocate("endyr_cpue_sBT");
  obs_cpue_sBT.allocate(styr_cpue_sBT,endyr_cpue_sBT,"obs_cpue_sBT");
  obs_cv_cpue_sBT.allocate(styr_cpue_sBT,endyr_cpue_sBT,"obs_cv_cpue_sBT");
  nyr_lenc_sBT.allocate("nyr_lenc_sBT");
  yrs_lenc_sBT.allocate(1,nyr_lenc_sBT,"yrs_lenc_sBT");
  nsamp_lenc_sBT.allocate(1,nyr_lenc_sBT,"nsamp_lenc_sBT");
  nfish_lenc_sBT.allocate(1,nyr_lenc_sBT,"nfish_lenc_sBT");
  obs_lenc_sBT.allocate(1,nyr_lenc_sBT,1,nlenbins,"obs_lenc_sBT");
  nyr_agec_sBT.allocate("nyr_agec_sBT");
  yrs_agec_sBT.allocate(1,nyr_agec_sBT,"yrs_agec_sBT");
  nsamp_agec_sBT.allocate(1,nyr_agec_sBT,"nsamp_agec_sBT");
  nfish_agec_sBT.allocate(1,nyr_agec_sBT,"nfish_agec_sBT");
  obs_agec_sBT.allocate(1,nyr_agec_sBT,1,nages,"obs_agec_sBT");
  styr_cpue_sTV.allocate("styr_cpue_sTV");
  endyr_cpue_sTV.allocate("endyr_cpue_sTV");
  obs_cpue_sTV.allocate(styr_cpue_sTV,endyr_cpue_sTV,"obs_cpue_sTV");
  obs_cv_cpue_sTV.allocate(styr_cpue_sTV,endyr_cpue_sTV,"obs_cv_cpue_sTV");
  nyr_agec_sTV.allocate("nyr_agec_sTV");
  yrs_agec_sTV.allocate(1,nyr_agec_sTV,"yrs_agec_sTV");
  nsamp_agec_sTV.allocate(1,nyr_agec_sTV,"nsamp_agec_sTV");
  nfish_agec_sTV.allocate(1,nyr_agec_sTV,"nfish_agec_sTV");
  obs_agec_sTV.allocate(1,nyr_agec_sTV,1,nages,"obs_agec_sTV");
  styr_cpue_cHL.allocate("styr_cpue_cHL");
  endyr_cpue_cHL.allocate("endyr_cpue_cHL");
  obs_cpue_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL,"obs_cpue_cHL");
  obs_cv_cpue_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL,"obs_cv_cpue_cHL");
  styr_L_cHL.allocate("styr_L_cHL");
  endyr_L_cHL.allocate("endyr_L_cHL");
  obs_L_cHL.allocate(styr_L_cHL,endyr_L_cHL,"obs_L_cHL");
  obs_cv_L_cHL.allocate(styr_L_cHL,endyr_L_cHL,"obs_cv_L_cHL");
  styr_D_cHL.allocate("styr_D_cHL");
  endyr_D_cHL.allocate("endyr_D_cHL");
  obs_released_cHL.allocate(styr_D_cHL,endyr_D_cHL,"obs_released_cHL");
  obs_cv_D_cHL.allocate(styr_D_cHL,endyr_D_cHL,"obs_cv_D_cHL");
  styr_D_cHL_closed.allocate("styr_D_cHL_closed");
  endyr_D_cHL_closed.allocate("endyr_D_cHL_closed");
  obs_released_cHL_closed.allocate(styr_D_cHL_closed,endyr_D_cHL_closed,"obs_released_cHL_closed");
  nyr_lenc_cHL.allocate("nyr_lenc_cHL");
  yrs_lenc_cHL.allocate(1,nyr_lenc_cHL,"yrs_lenc_cHL");
  nsamp_lenc_cHL.allocate(1,nyr_lenc_cHL,"nsamp_lenc_cHL");
  nfish_lenc_cHL.allocate(1,nyr_lenc_cHL,"nfish_lenc_cHL");
  obs_lenc_cHL.allocate(1,nyr_lenc_cHL,1,nlenbins,"obs_lenc_cHL");
  nyr_agec_cHL.allocate("nyr_agec_cHL");
  yrs_agec_cHL.allocate(1,nyr_agec_cHL,"yrs_agec_cHL");
  nsamp_agec_cHL.allocate(1,nyr_agec_cHL,"nsamp_agec_cHL");
  nfish_agec_cHL.allocate(1,nyr_agec_cHL,"nfish_agec_cHL");
  obs_agec_cHL.allocate(1,nyr_agec_cHL,1,nages,"obs_agec_cHL");
  styr_L_cPT.allocate("styr_L_cPT");
  endyr_L_cPT.allocate("endyr_L_cPT");
  obs_L_cPT.allocate(styr_L_cPT,endyr_L_cPT,"obs_L_cPT");
  obs_cv_L_cPT.allocate(styr_L_cPT,endyr_L_cPT,"obs_cv_L_cPT");
  styr_D_cPT.allocate("styr_D_cPT");
  endyr_D_cPT.allocate("endyr_D_cPT");
  obs_released_cPT.allocate(styr_D_cPT,endyr_D_cPT,"obs_released_cPT");
  obs_cv_D_cPT.allocate(styr_D_cPT,endyr_D_cPT,"obs_cv_D_cPT");
  styr_D_cPT_closed.allocate("styr_D_cPT_closed");
  endyr_D_cPT_closed.allocate("endyr_D_cPT_closed");
  obs_released_cPT_closed.allocate(styr_D_cPT_closed,endyr_D_cPT_closed,"obs_released_cPT_closed");
  nyr_lenc_cPT.allocate("nyr_lenc_cPT");
  yrs_lenc_cPT.allocate(1,nyr_lenc_cPT,"yrs_lenc_cPT");
  nsamp_lenc_cPT.allocate(1,nyr_lenc_cPT,"nsamp_lenc_cPT");
  nfish_lenc_cPT.allocate(1,nyr_lenc_cPT,"nfish_lenc_cPT");
  obs_lenc_cPT.allocate(1,nyr_lenc_cPT,1,nlenbins,"obs_lenc_cPT");
  nyr_lenc_pool_cPT.allocate("nyr_lenc_pool_cPT");
  yrs_lenc_pool_cPT.allocate(1,nyr_lenc_pool_cPT,"yrs_lenc_pool_cPT");
  nsamp_lenc_pool_cPT.allocate(1,nyr_lenc_pool_cPT,"nsamp_lenc_pool_cPT");
  nyr_agec_cPT.allocate("nyr_agec_cPT");
  yrs_agec_cPT.allocate(1,nyr_agec_cPT,"yrs_agec_cPT");
  nsamp_agec_cPT.allocate(1,nyr_agec_cPT,"nsamp_agec_cPT");
  nfish_agec_cPT.allocate(1,nyr_agec_cPT,"nfish_agec_cPT");
  obs_agec_cPT.allocate(1,nyr_agec_cPT,1,nages,"obs_agec_cPT");
  styr_L_cTW.allocate("styr_L_cTW");
  endyr_L_cTW.allocate("endyr_L_cTW");
  obs_L_cTW.allocate(styr_L_cTW,endyr_L_cTW,"obs_L_cTW");
  obs_cv_L_cTW.allocate(styr_L_cTW,endyr_L_cTW,"obs_cv_L_cTW");
  styr_cpue_rHB.allocate("styr_cpue_rHB");
  endyr_cpue_rHB.allocate("endyr_cpue_rHB");
  obs_cpue_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"obs_cpue_rHB");
  obs_cv_cpue_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"obs_cv_cpue_rHB");
  styr_L_rHB.allocate("styr_L_rHB");
  endyr_L_rHB.allocate("endyr_L_rHB");
  obs_L_rHB.allocate(styr_L_rHB,endyr_L_rHB,"obs_L_rHB");
  obs_cv_L_rHB.allocate(styr_L_rHB,endyr_L_rHB,"obs_cv_L_rHB");
  styr_D_rHB.allocate("styr_D_rHB");
  endyr_D_rHB.allocate("endyr_D_rHB");
  obs_released_rHB.allocate(styr_D_rHB,endyr_D_rHB,"obs_released_rHB");
  obs_cv_D_rHB.allocate(styr_D_rHB,endyr_D_rHB,"obs_cv_D_rHB");
  nyr_lenc_rHB.allocate("nyr_lenc_rHB");
  yrs_lenc_rHB.allocate(1,nyr_lenc_rHB,"yrs_lenc_rHB");
  nsamp_lenc_rHB.allocate(1,nyr_lenc_rHB,"nsamp_lenc_rHB");
  nfish_lenc_rHB.allocate(1,nyr_lenc_rHB,"nfish_lenc_rHB");
  obs_lenc_rHB.allocate(1,nyr_lenc_rHB,1,nlenbins,"obs_lenc_rHB");
  nyr_agec_rHB.allocate("nyr_agec_rHB");
  yrs_agec_rHB.allocate(1,nyr_agec_rHB,"yrs_agec_rHB");
  nsamp_agec_rHB.allocate(1,nyr_agec_rHB,"nsamp_agec_rHB");
  nfish_agec_rHB.allocate(1,nyr_agec_rHB,"nfish_agec_rHB");
  obs_agec_rHB.allocate(1,nyr_agec_rHB,1,nages,"obs_agec_rHB");
  nyr_lenc_rHB_D.allocate("nyr_lenc_rHB_D");
  yrs_lenc_rHB_D.allocate(1,nyr_lenc_rHB_D,"yrs_lenc_rHB_D");
  nsamp_lenc_rHB_D.allocate(1,nyr_lenc_rHB_D,"nsamp_lenc_rHB_D");
  nfish_lenc_rHB_D.allocate(1,nyr_lenc_rHB_D,"nfish_lenc_rHB_D");
  obs_lenc_rHB_D.allocate(1,nyr_lenc_rHB_D,1,nlenbins,"obs_lenc_rHB_D");
  styr_L_rGN.allocate("styr_L_rGN");
  endyr_L_rGN.allocate("endyr_L_rGN");
  obs_L_rGN.allocate(styr_L_rGN,endyr_L_rGN,"obs_L_rGN");
  obs_cv_L_rGN.allocate(styr_L_rGN,endyr_L_rGN,"obs_cv_L_rGN");
  styr_D_rGN.allocate("styr_D_rGN");
  endyr_D_rGN.allocate("endyr_D_rGN");
  obs_released_rGN.allocate(styr_D_rGN,endyr_D_rGN,"obs_released_rGN");
  obs_cv_D_rGN.allocate(styr_D_rGN,endyr_D_rGN,"obs_cv_D_rGN");
  nyr_lenc_rGN.allocate("nyr_lenc_rGN");
  yrs_lenc_rGN.allocate(1,nyr_lenc_rGN,"yrs_lenc_rGN");
  nsamp_lenc_rGN.allocate(1,nyr_lenc_rGN,"nsamp_lenc_rGN");
  nfish_lenc_rGN.allocate(1,nyr_lenc_rGN,"nfish_lenc_rGN");
  obs_lenc_rGN.allocate(1,nyr_lenc_rGN,1,nlenbins,"obs_lenc_rGN");
  set_Dmort_HL.allocate("set_Dmort_HL");
  set_Dmort_rHB_HL.allocate("set_Dmort_rHB_HL");
  set_Dmort_rGN_HL.allocate("set_Dmort_rGN_HL");
  set_Dmort_cPT1.allocate("set_Dmort_cPT1");
  set_Dmort_cPT2.allocate("set_Dmort_cPT2");
  set_Linf.allocate(1,7,"set_Linf");
  set_K.allocate(1,7,"set_K");
  set_t0.allocate(1,7,"set_t0");
  set_len_cv.allocate(1,7,"set_len_cv");
  set_M_constant.allocate(1,7,"set_M_constant");
  set_steep.allocate(1,7,"set_steep");
  set_log_R0.allocate(1,7,"set_log_R0");
  set_R_autocorr.allocate(1,7,"set_R_autocorr");
  set_rec_sigma.allocate(1,7,"set_rec_sigma");
  set_log_dm_lenc_sBT.allocate(1,7,"set_log_dm_lenc_sBT");
  set_log_dm_lenc_sTV.allocate(1,7,"set_log_dm_lenc_sTV");
  set_log_dm_lenc_cHL.allocate(1,7,"set_log_dm_lenc_cHL");
  set_log_dm_lenc_cPT.allocate(1,7,"set_log_dm_lenc_cPT");
  set_log_dm_lenc_rHB.allocate(1,7,"set_log_dm_lenc_rHB");
  set_log_dm_lenc_rHB_D.allocate(1,7,"set_log_dm_lenc_rHB_D");
  set_log_dm_lenc_rGN.allocate(1,7,"set_log_dm_lenc_rGN");
  set_log_dm_agec_sBT.allocate(1,7,"set_log_dm_agec_sBT");
  set_log_dm_agec_sTV.allocate(1,7,"set_log_dm_agec_sTV");
  set_log_dm_agec_cHL.allocate(1,7,"set_log_dm_agec_cHL");
  set_log_dm_agec_cPT.allocate(1,7,"set_log_dm_agec_cPT");
  set_log_dm_agec_rHB.allocate(1,7,"set_log_dm_agec_rHB");
  set_selpar_A50_sBT.allocate(1,7,"set_selpar_A50_sBT");
  set_selpar_slope_sBT.allocate(1,7,"set_selpar_slope_sBT");
  set_selpar_A50_sTV.allocate(1,7,"set_selpar_A50_sTV");
  set_selpar_slope_sTV.allocate(1,7,"set_selpar_slope_sTV");
  set_selpar_A50_Vid.allocate(1,7,"set_selpar_A50_Vid");
  set_selpar_slope_Vid.allocate(1,7,"set_selpar_slope_Vid");
  set_selpar_A50_cHL2.allocate(1,7,"set_selpar_A50_cHL2");
  set_selpar_slope_cHL2.allocate(1,7,"set_selpar_slope_cHL2");
  set_selpar_A50_cHL3.allocate(1,7,"set_selpar_A50_cHL3");
  set_selpar_slope_cHL3.allocate(1,7,"set_selpar_slope_cHL3");
  set_selpar_A50_cHL4.allocate(1,7,"set_selpar_A50_cHL4");
  set_selpar_slope_cHL4.allocate(1,7,"set_selpar_slope_cHL4");
  set_selpar_A50_cPT2.allocate(1,7,"set_selpar_A50_cPT2");
  set_selpar_slope_cPT2.allocate(1,7,"set_selpar_slope_cPT2");
  set_selpar_A50_cPT3.allocate(1,7,"set_selpar_A50_cPT3");
  set_selpar_slope_cPT3.allocate(1,7,"set_selpar_slope_cPT3");
  set_selpar_A50_cPT4.allocate(1,7,"set_selpar_A50_cPT4");
  set_selpar_slope_cPT4.allocate(1,7,"set_selpar_slope_cPT4");
  set_selpar_A50_rHB1.allocate(1,7,"set_selpar_A50_rHB1");
  set_selpar_slope_rHB1.allocate(1,7,"set_selpar_slope_rHB1");
  set_selpar_A50_rHB2.allocate(1,7,"set_selpar_A50_rHB2");
  set_selpar_slope_rHB2.allocate(1,7,"set_selpar_slope_rHB2");
  set_selpar_A50_rHB3.allocate(1,7,"set_selpar_A50_rHB3");
  set_selpar_slope_rHB3.allocate(1,7,"set_selpar_slope_rHB3");
  set_selpar_A50_rHB4.allocate(1,7,"set_selpar_A50_rHB4");
  set_selpar_slope_rHB4.allocate(1,7,"set_selpar_slope_rHB4");
  set_selpar_A50_rHB5.allocate(1,7,"set_selpar_A50_rHB5");
  set_selpar_slope_rHB5.allocate(1,7,"set_selpar_slope_rHB5");
  set_selpar_A50_rGN1.allocate(1,7,"set_selpar_A50_rGN1");
  set_selpar_slope_rGN1.allocate(1,7,"set_selpar_slope_rGN1");
  set_selpar_A50_rGN2.allocate(1,7,"set_selpar_A50_rGN2");
  set_selpar_slope_rGN2.allocate(1,7,"set_selpar_slope_rGN2");
  set_selpar_A50_rGN3.allocate(1,7,"set_selpar_A50_rGN3");
  set_selpar_slope_rGN3.allocate(1,7,"set_selpar_slope_rGN3");
  set_selpar_A50_rGN4.allocate(1,7,"set_selpar_A50_rGN4");
  set_selpar_slope_rGN4.allocate(1,7,"set_selpar_slope_rGN4");
  set_selpar_A50_rGN5.allocate(1,7,"set_selpar_A50_rGN5");
  set_selpar_slope_rGN5.allocate(1,7,"set_selpar_slope_rGN5");
  set_selpar_logit_Age0_rHB_D.allocate(1,7,"set_selpar_logit_Age0_rHB_D");
  set_selpar_logit_Age1_rHB_D.allocate(1,7,"set_selpar_logit_Age1_rHB_D");
  set_selpar_logit_Age2_rHB_D.allocate(1,7,"set_selpar_logit_Age2_rHB_D");
  set_selpar_A50_rHB_D4.allocate(1,7,"set_selpar_A50_rHB_D4");
  set_selpar_slope_rHB_D4.allocate(1,7,"set_selpar_slope_rHB_D4");
  set_selpar_A502_rHB_D4.allocate(1,7,"set_selpar_A502_rHB_D4");
  set_selpar_slope2_rHB_D4.allocate(1,7,"set_selpar_slope2_rHB_D4");
  set_selpar_A50_rHB_D5.allocate(1,7,"set_selpar_A50_rHB_D5");
  set_selpar_slope_rHB_D5.allocate(1,7,"set_selpar_slope_rHB_D5");
  set_selpar_A502_rHB_D5.allocate(1,7,"set_selpar_A502_rHB_D5");
  set_selpar_slope2_rHB_D5.allocate(1,7,"set_selpar_slope2_rHB_D5");
  set_log_q_cpue_sBT.allocate(1,7,"set_log_q_cpue_sBT");
  set_log_q_cpue_sTV.allocate(1,7,"set_log_q_cpue_sTV");
  set_log_q_cpue_cHL.allocate(1,7,"set_log_q_cpue_cHL");
  set_log_q_cpue_rHB.allocate(1,7,"set_log_q_cpue_rHB");
  set_log_avg_F_L_cHL.allocate(1,7,"set_log_avg_F_L_cHL");
  set_log_avg_F_L_cPT.allocate(1,7,"set_log_avg_F_L_cPT");
  set_log_avg_F_L_cTW.allocate(1,7,"set_log_avg_F_L_cTW");
  set_log_avg_F_L_rHB.allocate(1,7,"set_log_avg_F_L_rHB");
  set_log_avg_F_L_rGN.allocate(1,7,"set_log_avg_F_L_rGN");
  set_log_avg_F_D_cGN.allocate(1,7,"set_log_avg_F_D_cGN");
  set_log_avg_F_D_rHB.allocate(1,7,"set_log_avg_F_D_rHB");
  set_log_avg_F_D_rGN.allocate(1,7,"set_log_avg_F_D_rGN");
  set_log_dev_F_L_cHL.allocate(1,3,"set_log_dev_F_L_cHL");
  set_log_dev_F_L_cPT.allocate(1,3,"set_log_dev_F_L_cPT");
  set_log_dev_F_L_cTW.allocate(1,3,"set_log_dev_F_L_cTW");
  set_log_dev_F_L_rHB.allocate(1,3,"set_log_dev_F_L_rHB");
  set_log_dev_F_L_rGN.allocate(1,3,"set_log_dev_F_L_rGN");
  set_log_dev_F_D_cGN.allocate(1,3,"set_log_dev_F_D_cGN");
  set_log_dev_F_D_rHB.allocate(1,3,"set_log_dev_F_D_rHB");
  set_log_dev_F_D_rGN.allocate(1,3,"set_log_dev_F_D_rGN");
  set_log_dev_RWq.allocate(1,3,"set_log_dev_RWq");
  set_log_dev_rec.allocate(1,3,"set_log_dev_rec");
  set_log_dev_Nage.allocate(1,3,"set_log_dev_Nage");
  set_log_dev_vals_F_L__cHL.allocate(styr_L_cHL,endyr_L_cHL,"set_log_dev_vals_F_L__cHL");
  set_log_dev_vals_F_L__cPT.allocate(styr_L_cPT,endyr_L_cPT,"set_log_dev_vals_F_L__cPT");
  set_log_dev_vals_F_L__rHB.allocate(styr_L_rHB,endyr_L_rHB,"set_log_dev_vals_F_L__rHB");
  set_log_dev_vals_F_L__rGN.allocate(styr_L_rGN,endyr_L_rGN,"set_log_dev_vals_F_L__rGN");
  set_log_dev_vals_F_D__commvals.allocate(styr_D_cHL,endyr_D_cHL,"set_log_dev_vals_F_D__commvals");
  set_log_dev_vals_F_D__HBvals.allocate(styr_D_rHB,endyr_D_rHB,"set_log_dev_vals_F_D__HBvals");
  set_log_dev_vals_F_D__mripvals.allocate(styr_D_rGN,endyr_D_rGN,"set_log_dev_vals_F_D__mripvals");
  set_log_dev_vals_rec.allocate(styr_rec_dev,endyr_rec_phase2,"set_log_dev_vals_rec");
  set_log_dev_vals_Nage.allocate(2,nages,"set_log_dev_vals_Nage");
  set_w_L.allocate("set_w_L");
  set_w_D.allocate("set_w_D");
  set_w_lenc_sBT.allocate("set_w_lenc_sBT");
  set_w_lenc_cHL.allocate("set_w_lenc_cHL");
  set_w_lenc_cPT.allocate("set_w_lenc_cPT");
  set_w_lenc_rHB.allocate("set_w_lenc_rHB");
  set_w_lenc_rHB_D.allocate("set_w_lenc_rHB_D");
  set_w_lenc_rGN.allocate("set_w_lenc_rGN");
  set_w_agec_sBT.allocate("set_w_agec_sBT");
  set_w_agec_sTV.allocate("set_w_agec_sTV");
  set_w_agec_cHL.allocate("set_w_agec_cHL");
  set_w_agec_cPT.allocate("set_w_agec_cPT");
  set_w_agec_rHB.allocate("set_w_agec_rHB");
  set_w_cpue_sBT.allocate("set_w_cpue_sBT");
  set_w_cpue_sTV.allocate("set_w_cpue_sTV");
  set_w_cpue_cHL.allocate("set_w_cpue_cHL");
  set_w_cpue_rHB.allocate("set_w_cpue_rHB");
  set_w_rec.allocate("set_w_rec");
  set_w_rec_early.allocate("set_w_rec_early");
  set_w_rec_end.allocate("set_w_rec_end");
  set_w_fullF.allocate("set_w_fullF");
  set_w_Ftune.allocate("set_w_Ftune");
  wgtpar_a.allocate("wgtpar_a");
  wgtpar_b.allocate("wgtpar_b");
  fecpar_a.allocate("fecpar_a");
  fecpar_b.allocate("fecpar_b");
  fecpar_batches.allocate("fecpar_batches");
  fecpar_scale.allocate("fecpar_scale");
  obs_maturity_f.allocate(1,nages,"obs_maturity_f");
  obs_prop_f.allocate(1,nages,"obs_prop_f");
  spawn_time_frac.allocate("spawn_time_frac");
  set_M.allocate(1,nages,"set_M");
  max_obs_age.allocate("max_obs_age");
  set_L_rHB_bias.allocate("set_L_rHB_bias");
  set_L_rGN_bias.allocate("set_L_rGN_bias");
  set_L_cGN_bias.allocate("set_L_cGN_bias");
  endyr_L_rHB_bias.allocate("endyr_L_rHB_bias");
  endyr_L_rGN_bias.allocate("endyr_L_rGN_bias");
  endyr_L_cGN_bias.allocate("endyr_L_cGN_bias");
  set_q_rate_phase.allocate("set_q_rate_phase");
  set_q_rate.allocate("set_q_rate");
  set_q_DD_phase.allocate("set_q_DD_phase");
  set_q_DD_beta.allocate("set_q_DD_beta");
  set_q_DD_beta_se.allocate("set_q_DD_beta_se");
  set_q_DD_stage.allocate("set_q_DD_stage");
  set_q_RW_phase.allocate("set_q_RW_phase");
  set_q_RW_cHL_var.allocate("set_q_RW_cHL_var");
  set_q_RW_rHB_var.allocate("set_q_RW_rHB_var");
  set_F_init_ratio.allocate("set_F_init_ratio");
  set_Ftune.allocate("set_Ftune");
  set_Ftune_yr.allocate("set_Ftune_yr");
  minSS_lenc.allocate("minSS_lenc");
  minSS_agec.allocate("minSS_agec");
  maxSS_lenc.allocate("maxSS_lenc");
  maxSS_agec.allocate("maxSS_agec");
  age_error.allocate(1,nages,1,nages,"age_error");
  set_p_lenc_cHL2.allocate("set_p_lenc_cHL2");
  set_p_lenc_cHL3.allocate("set_p_lenc_cHL3");
  set_p_lenc_cPT2.allocate("set_p_lenc_cPT2");
  set_p_lenc_cPT3.allocate("set_p_lenc_cPT3");
  set_p_lenc_cTW2.allocate("set_p_lenc_cTW2");
  set_p_lenc_cTW3.allocate("set_p_lenc_cTW3");
  set_p_lenc_rHB2.allocate("set_p_lenc_rHB2");
  set_p_lenc_rHB3.allocate("set_p_lenc_rHB3");
  set_p_lenc_rHB4.allocate("set_p_lenc_rHB4");
  set_p_lenc_rHB5.allocate("set_p_lenc_rHB5");
  set_p_lenc_rGN2.allocate("set_p_lenc_rGN2");
  set_p_lenc_rGN3.allocate("set_p_lenc_rGN3");
  set_p_lenc_rGN4.allocate("set_p_lenc_rGN4");
  set_p_lenc_rGN5.allocate("set_p_lenc_rGN5");
  set_p_lenc_cGN_D2.allocate("set_p_lenc_cGN_D2");
  set_p_lenc_cGN_D3.allocate("set_p_lenc_cGN_D3");
  set_p_lenc_cGN_D4.allocate("set_p_lenc_cGN_D4");
  set_p_lenc_rHB_D2.allocate("set_p_lenc_rHB_D2");
  set_p_lenc_rHB_D3.allocate("set_p_lenc_rHB_D3");
  set_p_lenc_rHB_D4.allocate("set_p_lenc_rHB_D4");
  set_p_lenc_rHB_D5.allocate("set_p_lenc_rHB_D5");
  set_p_lenc_rGN_D1.allocate("set_p_lenc_rGN_D1");
  set_p_lenc_rGN_D2.allocate("set_p_lenc_rGN_D2");
  set_p_lenc_rGN_D3.allocate("set_p_lenc_rGN_D3");
  set_p_lenc_rGN_D4.allocate("set_p_lenc_rGN_D4");
  set_p_lenc_rGN_D5.allocate("set_p_lenc_rGN_D5");
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
     
  const double M_constant_LO=set_M_constant(2); const double M_constant_HI=set_M_constant(3); const double M_constant_PH=set_M_constant(4);        
  const double steep_LO=set_steep(2); const double steep_HI=set_steep(3); const double steep_PH=set_steep(4);
  const double log_R0_LO=set_log_R0(2); const double log_R0_HI=set_log_R0(3); const double log_R0_PH=set_log_R0(4);
  const double R_autocorr_LO=set_R_autocorr(2); const double R_autocorr_HI=set_R_autocorr(3); const double R_autocorr_PH=set_R_autocorr(4);
  const double rec_sigma_LO=set_rec_sigma(2); const double rec_sigma_HI=set_rec_sigma(3); const double rec_sigma_PH=set_rec_sigma(4);
  
  const double log_dm_sBT_lc_LO=set_log_dm_lenc_sBT(2); const double log_dm_sBT_lc_HI=set_log_dm_lenc_sBT(3); const double log_dm_sBT_lc_PH=set_log_dm_lenc_sBT(4);
  const double log_dm_sTV_lc_LO=set_log_dm_lenc_sTV(2); const double log_dm_sTV_lc_HI=set_log_dm_lenc_sTV(3); const double log_dm_sTV_lc_PH=set_log_dm_lenc_sTV(4);
  const double log_dm_rHB_lc_LO=set_log_dm_lenc_rHB(2); const double log_dm_rHB_lc_HI=set_log_dm_lenc_rHB(3); const double log_dm_rHB_lc_PH=set_log_dm_lenc_rHB(4);
  const double log_dm_rGN_lc_LO=set_log_dm_lenc_rGN(2); const double log_dm_rGN_lc_HI=set_log_dm_lenc_rGN(3); const double log_dm_rGN_lc_PH=set_log_dm_lenc_rGN(4);
  const double log_dm_rHB_D_lc_LO=set_log_dm_lenc_rHB_D(2); const double log_dm_rHB_D_lc_HI=set_log_dm_lenc_rHB_D(3); const double log_dm_rHB_D_lc_PH=set_log_dm_lenc_rHB_D(4);
  const double log_dm_cHL_lc_LO=set_log_dm_lenc_cHL(2); const double log_dm_cHL_lc_HI=set_log_dm_lenc_cHL(3); const double log_dm_cHL_lc_PH=set_log_dm_lenc_cHL(4);
  const double log_dm_cPT_lc_LO=set_log_dm_lenc_cPT(2); const double log_dm_cPT_lc_HI=set_log_dm_lenc_cPT(3); const double log_dm_cPT_lc_PH=set_log_dm_lenc_cPT(4);
  const double log_dm_sBT_ac_LO=set_log_dm_agec_sBT(2); const double log_dm_sBT_ac_HI=set_log_dm_agec_sBT(3); const double log_dm_sBT_ac_PH=set_log_dm_agec_sBT(4);
  const double log_dm_sTV_ac_LO=set_log_dm_agec_sTV(2); const double log_dm_sTV_ac_HI=set_log_dm_agec_sTV(3); const double log_dm_sTV_ac_PH=set_log_dm_agec_sTV(4);
  const double log_dm_rHB_ac_LO=set_log_dm_agec_rHB(2); const double log_dm_rHB_ac_HI=set_log_dm_agec_rHB(3); const double log_dm_rHB_ac_PH=set_log_dm_agec_rHB(4);
  //const double log_dm_rGN_ac_LO=set_log_dm_rGN_ac(2); const double log_dm_rGN_ac_HI=set_log_dm_rGN_ac(3); const double log_dm_rGN_ac_PH=set_log_dm_rGN_ac(4);
  const double log_dm_cHL_ac_LO=set_log_dm_agec_cHL(2); const double log_dm_cHL_ac_HI=set_log_dm_agec_cHL(3); const double log_dm_cHL_ac_PH=set_log_dm_agec_cHL(4);
  const double log_dm_cPT_ac_LO=set_log_dm_agec_cPT(2); const double log_dm_cPT_ac_HI=set_log_dm_agec_cPT(3); const double log_dm_cPT_ac_PH=set_log_dm_agec_cPT(4);
  
  const double selpar_A50_sBT_LO=set_selpar_A50_sBT(2); const double selpar_A50_sBT_HI=set_selpar_A50_sBT(3); const double selpar_A50_sBT_PH=set_selpar_A50_sBT(4);
  const double selpar_slope_sBT_LO=set_selpar_slope_sBT(2); const double selpar_slope_sBT_HI=set_selpar_slope_sBT(3); const double selpar_slope_sBT_PH=set_selpar_slope_sBT(4);
  
  const double selpar_A50_sTV_LO=set_selpar_A50_sTV(2); const double selpar_A50_sTV_HI=set_selpar_A50_sTV(3); const double selpar_A50_sTV_PH=set_selpar_A50_sTV(4);
  const double selpar_slope_sTV_LO=set_selpar_slope_sTV(2); const double selpar_slope_sTV_HI=set_selpar_slope_sTV(3); const double selpar_slope_sTV_PH=set_selpar_slope_sTV(4);
  const double selpar_A50_cHL2_LO=set_selpar_A50_cHL2(2); const double selpar_A50_cHL2_HI=set_selpar_A50_cHL2(3); const double selpar_A50_cHL2_PH=set_selpar_A50_cHL2(4);
  const double selpar_slope_cHL2_LO=set_selpar_slope_cHL2(2); const double selpar_slope_cHL2_HI=set_selpar_slope_cHL2(3); const double selpar_slope_cHL2_PH=set_selpar_slope_cHL2(4);
  const double selpar_A50_cHL3_LO=set_selpar_A50_cHL3(2); const double selpar_A50_cHL3_HI=set_selpar_A50_cHL3(3); const double selpar_A50_cHL3_PH=set_selpar_A50_cHL3(4);
  const double selpar_slope_cHL3_LO=set_selpar_slope_cHL3(2); const double selpar_slope_cHL3_HI=set_selpar_slope_cHL3(3); const double selpar_slope_cHL3_PH=set_selpar_slope_cHL3(4);
  const double selpar_A50_cHL4_LO=set_selpar_A50_cHL4(2); const double selpar_A50_cHL4_HI=set_selpar_A50_cHL4(3); const double selpar_A50_cHL4_PH=set_selpar_A50_cHL4(4);
  const double selpar_slope_cHL4_LO=set_selpar_slope_cHL4(2); const double selpar_slope_cHL4_HI=set_selpar_slope_cHL4(3); const double selpar_slope_cHL4_PH=set_selpar_slope_cHL4(4);
  
  const double selpar_A50_cPT2_LO=set_selpar_A50_cPT2(2); const double selpar_A50_cPT2_HI=set_selpar_A50_cPT2(3); const double selpar_A50_cPT2_PH=set_selpar_A50_cPT2(4);
  const double selpar_slope_cPT2_LO=set_selpar_slope_cPT2(2); const double selpar_slope_cPT2_HI=set_selpar_slope_cPT2(3); const double selpar_slope_cPT2_PH=set_selpar_slope_cPT2(4);
  const double selpar_A50_cPT3_LO=set_selpar_A50_cPT3(2); const double selpar_A50_cPT3_HI=set_selpar_A50_cPT3(3); const double selpar_A50_cPT3_PH=set_selpar_A50_cPT3(4);
  const double selpar_slope_cPT3_LO=set_selpar_slope_cPT3(2); const double selpar_slope_cPT3_HI=set_selpar_slope_cPT3(3); const double selpar_slope_cPT3_PH=set_selpar_slope_cPT3(4);
  const double selpar_A50_cPT4_LO=set_selpar_A50_cPT4(2); const double selpar_A50_cPT4_HI=set_selpar_A50_cPT4(3); const double selpar_A50_cPT4_PH=set_selpar_A50_cPT4(4);
  const double selpar_slope_cPT4_LO=set_selpar_slope_cPT4(2); const double selpar_slope_cPT4_HI=set_selpar_slope_cPT4(3); const double selpar_slope_cPT4_PH=set_selpar_slope_cPT4(4);
  
  const double selpar_A50_rHB1_LO=set_selpar_A50_rHB1(2); const double selpar_A50_rHB1_HI=set_selpar_A50_rHB1(3); const double selpar_A50_rHB1_PH=set_selpar_A50_rHB1(4);
  const double selpar_slope_rHB1_LO=set_selpar_slope_rHB1(2); const double selpar_slope_rHB1_HI=set_selpar_slope_rHB1(3); const double selpar_slope_rHB1_PH=set_selpar_slope_rHB1(4);
  const double selpar_A50_rHB2_LO=set_selpar_A50_rHB2(2); const double selpar_A50_rHB2_HI=set_selpar_A50_rHB2(3); const double selpar_A50_rHB2_PH=set_selpar_A50_rHB2(4);
  const double selpar_slope_rHB2_LO=set_selpar_slope_rHB2(2); const double selpar_slope_rHB2_HI=set_selpar_slope_rHB2(3); const double selpar_slope_rHB2_PH=set_selpar_slope_rHB2(4);
  const double selpar_A50_rHB3_LO=set_selpar_A50_rHB3(2); const double selpar_A50_rHB3_HI=set_selpar_A50_rHB3(3); const double selpar_A50_rHB3_PH=set_selpar_A50_rHB3(4);
  const double selpar_slope_rHB3_LO=set_selpar_slope_rHB3(2); const double selpar_slope_rHB3_HI=set_selpar_slope_rHB3(3); const double selpar_slope_rHB3_PH=set_selpar_slope_rHB3(4);
  const double selpar_A50_rHB4_LO=set_selpar_A50_rHB4(2); const double selpar_A50_rHB4_HI=set_selpar_A50_rHB4(3); const double selpar_A50_rHB4_PH=set_selpar_A50_rHB4(4);
  const double selpar_slope_rHB4_LO=set_selpar_slope_rHB4(2); const double selpar_slope_rHB4_HI=set_selpar_slope_rHB4(3); const double selpar_slope_rHB4_PH=set_selpar_slope_rHB4(4);
  const double selpar_A50_rHB5_LO=set_selpar_A50_rHB5(2); const double selpar_A50_rHB5_HI=set_selpar_A50_rHB5(3); const double selpar_A50_rHB5_PH=set_selpar_A50_rHB5(4);
  const double selpar_slope_rHB5_LO=set_selpar_slope_rHB5(2); const double selpar_slope_rHB5_HI=set_selpar_slope_rHB5(3); const double selpar_slope_rHB5_PH=set_selpar_slope_rHB5(4);
  
  const double selpar_A50_rGN1_LO=set_selpar_A50_rGN1(2); const double selpar_A50_rGN1_HI=set_selpar_A50_rGN1(3); const double selpar_A50_rGN1_PH=set_selpar_A50_rGN1(4);
  const double selpar_slope_rGN1_LO=set_selpar_slope_rGN1(2); const double selpar_slope_rGN1_HI=set_selpar_slope_rGN1(3); const double selpar_slope_rGN1_PH=set_selpar_slope_rGN1(4);
  const double selpar_A50_rGN2_LO=set_selpar_A50_rGN2(2); const double selpar_A50_rGN2_HI=set_selpar_A50_rGN2(3); const double selpar_A50_rGN2_PH=set_selpar_A50_rGN2(4);
  const double selpar_slope_rGN2_LO=set_selpar_slope_rGN2(2); const double selpar_slope_rGN2_HI=set_selpar_slope_rGN2(3); const double selpar_slope_rGN2_PH=set_selpar_slope_rGN2(4);
  const double selpar_A50_rGN3_LO=set_selpar_A50_rGN3(2); const double selpar_A50_rGN3_HI=set_selpar_A50_rGN3(3); const double selpar_A50_rGN3_PH=set_selpar_A50_rGN3(4);
  const double selpar_slope_rGN3_LO=set_selpar_slope_rGN3(2); const double selpar_slope_rGN3_HI=set_selpar_slope_rGN3(3); const double selpar_slope_rGN3_PH=set_selpar_slope_rGN3(4);
  const double selpar_A50_rGN4_LO=set_selpar_A50_rGN4(2); const double selpar_A50_rGN4_HI=set_selpar_A50_rGN4(3); const double selpar_A50_rGN4_PH=set_selpar_A50_rGN4(4);
  const double selpar_slope_rGN4_LO=set_selpar_slope_rGN4(2); const double selpar_slope_rGN4_HI=set_selpar_slope_rGN4(3); const double selpar_slope_rGN4_PH=set_selpar_slope_rGN4(4);
  const double selpar_A50_rGN5_LO=set_selpar_A50_rGN5(2); const double selpar_A50_rGN5_HI=set_selpar_A50_rGN5(3); const double selpar_A50_rGN5_PH=set_selpar_A50_rGN5(4);
  const double selpar_slope_rGN5_LO=set_selpar_slope_rGN5(2); const double selpar_slope_rGN5_HI=set_selpar_slope_rGN5(3); const double selpar_slope_rGN5_PH=set_selpar_slope_rGN5(4);
    
  const double selpar_Age0_rHB_D_logit_LO=set_selpar_logit_Age0_rHB_D(2); const double selpar_Age0_rHB_D_logit_HI=set_selpar_logit_Age0_rHB_D(3); const double selpar_Age0_rHB_D_logit_PH=set_selpar_logit_Age0_rHB_D(4);
  const double selpar_Age1_rHB_D_logit_LO=set_selpar_logit_Age1_rHB_D(2); const double selpar_Age1_rHB_D_logit_HI=set_selpar_logit_Age1_rHB_D(3); const double selpar_Age1_rHB_D_logit_PH=set_selpar_logit_Age1_rHB_D(4);
  const double selpar_Age2_rHB_D_logit_LO=set_selpar_logit_Age2_rHB_D(2); const double selpar_Age2_rHB_D_logit_HI=set_selpar_logit_Age2_rHB_D(3); const double selpar_Age2_rHB_D_logit_PH=set_selpar_logit_Age2_rHB_D(4);
  
  const double selpar_A50_rHB_D4_LO=set_selpar_A50_rHB_D4(2); 
  const double selpar_A50_rHB_D4_HI=set_selpar_A50_rHB_D4(3); 
  const double selpar_A50_rHB_D4_PH=set_selpar_A50_rHB_D4(4);
  const double selpar_slope_rHB_D4_LO=set_selpar_slope_rHB_D4(2); 
  const double selpar_slope_rHB_D4_HI=set_selpar_slope_rHB_D4(3); 
  const double selpar_slope_rHB_D4_PH=set_selpar_slope_rHB_D4(4);
  const double selpar_A502_rHB_D4_LO=set_selpar_A502_rHB_D4(2); 
  const double selpar_A502_rHB_D4_HI=set_selpar_A502_rHB_D4(3); 
  const double selpar_A502_rHB_D4_PH=set_selpar_A502_rHB_D4(4);
  const double selpar_slope2_rHB_D4_LO=set_selpar_slope2_rHB_D4(2); 
  const double selpar_slope2_rHB_D4_HI=set_selpar_slope2_rHB_D4(3); 
  const double selpar_slope2_rHB_D4_PH=set_selpar_slope2_rHB_D4(4);
  const double selpar_A50_rHB_D5_LO=set_selpar_A50_rHB_D5(2); 
  const double selpar_A50_rHB_D5_HI=set_selpar_A50_rHB_D5(3); 
  const double selpar_A50_rHB_D5_PH=set_selpar_A50_rHB_D5(4);
  const double selpar_slope_rHB_D5_LO=set_selpar_slope_rHB_D5(2); 
  const double selpar_slope_rHB_D5_HI=set_selpar_slope_rHB_D5(3); 
  const double selpar_slope_rHB_D5_PH=set_selpar_slope_rHB_D5(4);
  const double selpar_A502_rHB_D5_LO=set_selpar_A502_rHB_D5(2); 
  const double selpar_A502_rHB_D5_HI=set_selpar_A502_rHB_D5(3); 
  const double selpar_A502_rHB_D5_PH=set_selpar_A502_rHB_D5(4);
  const double selpar_slope2_rHB_D5_LO=set_selpar_slope2_rHB_D5(2); 
  const double selpar_slope2_rHB_D5_HI=set_selpar_slope2_rHB_D5(3); 
  const double selpar_slope2_rHB_D5_PH=set_selpar_slope2_rHB_D5(4);
  
  const double log_q_sBT_LO=set_log_q_cpue_sBT(2); const double log_q_sBT_HI=set_log_q_cpue_sBT(3); const double log_q_sBT_PH=set_log_q_cpue_sBT(4);
  const double log_q_rHB_LO=set_log_q_cpue_rHB(2); const double log_q_rHB_HI=set_log_q_cpue_rHB(3); const double log_q_rHB_PH=set_log_q_cpue_rHB(4);
  const double log_q_sTV_LO=set_log_q_cpue_sTV(2); const double log_q_sTV_HI=set_log_q_cpue_sTV(3); const double log_q_sTV_PH=set_log_q_cpue_sTV(4);
  //const double log_q_Vid_LO=set_logq_Vid(2); const double log_q_Vid_HI=set_logq_Vid(3); const double log_q_Vid_PH=set_logq_Vid(4);
  const double log_q_cHL_LO=set_log_q_cpue_cHL(2); const double log_q_cHL_HI=set_log_q_cpue_cHL(3); const double log_q_cHL_PH=set_log_q_cpue_cHL(4);
  //const double log_q_rHB_D_LO=set_logq_rHB_D(2); const double log_q_rHB_D_HI=set_logq_rHB_D(3); const double log_q_rHB_D_PH=set_logq_rHB_D(4);
  
  const double log_avg_F_cHL_LO=set_log_avg_F_L_cHL(2); const double log_avg_F_cHL_HI=set_log_avg_F_L_cHL(3); const double log_avg_F_cHL_PH=set_log_avg_F_L_cHL(4);
  const double log_avg_F_cPT_LO=set_log_avg_F_L_cPT(2); const double log_avg_F_cPT_HI=set_log_avg_F_L_cPT(3); const double log_avg_F_cPT_PH=set_log_avg_F_L_cPT(4);
  const double log_avg_F_cTW_LO=set_log_avg_F_L_cTW(2); const double log_avg_F_cTW_HI=set_log_avg_F_L_cTW(3); const double log_avg_F_cTW_PH=set_log_avg_F_L_cTW(4);
  const double log_avg_F_rHB_LO=set_log_avg_F_L_rHB(2); const double log_avg_F_rHB_HI=set_log_avg_F_L_rHB(3); const double log_avg_F_rHB_PH=set_log_avg_F_L_rHB(4); 
  const double log_avg_F_rGN_LO=set_log_avg_F_L_rGN(2); const double log_avg_F_rGN_HI=set_log_avg_F_L_rGN(3); const double log_avg_F_rGN_PH=set_log_avg_F_L_rGN(4); 
  
  const double log_avg_F_cGN_D_LO=set_log_avg_F_D_cGN(2); const double log_avg_F_cGN_D_HI=set_log_avg_F_D_cGN(3); const double log_avg_F_cGN_D_PH=set_log_avg_F_D_cGN(4);
  const double log_avg_F_rHB_D_LO=set_log_avg_F_D_rHB(2); const double log_avg_F_rHB_D_HI=set_log_avg_F_D_rHB(3); const double log_avg_F_rHB_D_PH=set_log_avg_F_D_rHB(4); 
  const double log_avg_F_rGN_D_LO=set_log_avg_F_D_rGN(2); const double log_avg_F_rGN_D_HI=set_log_avg_F_D_rGN(3); const double log_avg_F_rGN_D_PH=set_log_avg_F_D_rGN(4); 
  
  //-dev vectors-----------------------------------------------------------------------------------------------------------  
  const double log_F_dev_cHL_LO=set_log_dev_F_L_cHL(1); const double log_F_dev_cHL_HI=set_log_dev_F_L_cHL(2); const double log_F_dev_cHL_PH=set_log_dev_F_L_cHL(3);  
  const double log_F_dev_cPT_LO=set_log_dev_F_L_cPT(1); const double log_F_dev_cPT_HI=set_log_dev_F_L_cPT(2); const double log_F_dev_cPT_PH=set_log_dev_F_L_cPT(3);   
  const double log_F_dev_cTW_LO=set_log_dev_F_L_cTW(1); const double log_F_dev_cTW_HI=set_log_dev_F_L_cTW(2); const double log_F_dev_cTW_PH=set_log_dev_F_L_cTW(3);   
  const double log_F_dev_rHB_LO=set_log_dev_F_L_rHB(1); const double log_F_dev_rHB_HI=set_log_dev_F_L_rHB(2); const double log_F_dev_rHB_PH=set_log_dev_F_L_rHB(3);   
  const double log_F_dev_rGN_LO=set_log_dev_F_L_rGN(1); const double log_F_dev_rGN_HI=set_log_dev_F_L_rGN(2); const double log_F_dev_rGN_PH=set_log_dev_F_L_rGN(3);   
  
  const double log_F_dev_cGN_D_LO=set_log_dev_F_D_cGN(1); const double log_F_dev_cGN_D_HI=set_log_dev_F_D_cGN(2); const double log_F_dev_cGN_D_PH=set_log_dev_F_D_cGN(3);   
  const double log_F_dev_rHB_D_LO=set_log_dev_F_D_rHB(1); const double log_F_dev_rHB_D_HI=set_log_dev_F_D_rHB(2); const double log_F_dev_rHB_D_PH=set_log_dev_F_D_rHB(3);   
  const double log_F_dev_rGN_D_LO=set_log_dev_F_D_rGN(1); const double log_F_dev_rGN_D_HI=set_log_dev_F_D_rGN(2); const double log_F_dev_rGN_D_PH=set_log_dev_F_D_rGN(3);   
  //const double F_init_ratio_LO=set_F_init_ratio(1); const double F_init_ratio_HI=set_F_init_ratio(2); const double F_init_ratio_PH=set_F_init_ratio(3); 
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
  fecundity.allocate(1,nages,"fecundity");
  #ifndef NO_AD_INITIALIZE
    fecundity.initialize();
  #endif
  len_cHL_mm.allocate(styr,endyr,1,nages,"len_cHL_mm");
  #ifndef NO_AD_INITIALIZE
    len_cHL_mm.initialize();
  #endif
  wgt_cHL_klb.allocate(styr,endyr,1,nages,"wgt_cHL_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_cHL_klb.initialize();
  #endif
  len_cPT_mm.allocate(styr,endyr,1,nages,"len_cPT_mm");
  #ifndef NO_AD_INITIALIZE
    len_cPT_mm.initialize();
  #endif
  wgt_cPT_klb.allocate(styr,endyr,1,nages,"wgt_cPT_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_cPT_klb.initialize();
  #endif
  len_cTW_mm.allocate(styr,endyr,1,nages,"len_cTW_mm");
  #ifndef NO_AD_INITIALIZE
    len_cTW_mm.initialize();
  #endif
  wgt_cTW_klb.allocate(styr,endyr,1,nages,"wgt_cTW_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_cTW_klb.initialize();
  #endif
  len_rHB_mm.allocate(styr,endyr,1,nages,"len_rHB_mm");
  #ifndef NO_AD_INITIALIZE
    len_rHB_mm.initialize();
  #endif
  wgt_rHB_klb.allocate(styr,endyr,1,nages,"wgt_rHB_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_rHB_klb.initialize();
  #endif
  len_rGN_mm.allocate(styr,endyr,1,nages,"len_rGN_mm");
  #ifndef NO_AD_INITIALIZE
    len_rGN_mm.initialize();
  #endif
  wgt_rGN_klb.allocate(styr,endyr,1,nages,"wgt_rGN_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_rGN_klb.initialize();
  #endif
  len_cGN_D_mm.allocate(styr,endyr,1,nages,"len_cGN_D_mm");
  #ifndef NO_AD_INITIALIZE
    len_cGN_D_mm.initialize();
  #endif
  wgt_cGN_D_klb.allocate(styr,endyr,1,nages,"wgt_cGN_D_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_cGN_D_klb.initialize();
  #endif
  len_rHB_D_mm.allocate(styr,endyr,1,nages,"len_rHB_D_mm");
  #ifndef NO_AD_INITIALIZE
    len_rHB_D_mm.initialize();
  #endif
  wgt_rHB_D_klb.allocate(styr,endyr,1,nages,"wgt_rHB_D_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_rHB_D_klb.initialize();
  #endif
  len_rGN_D_mm.allocate(styr,endyr,1,nages,"len_rGN_D_mm");
  #ifndef NO_AD_INITIALIZE
    len_rGN_D_mm.initialize();
  #endif
  wgt_rGN_D_klb.allocate(styr,endyr,1,nages,"wgt_rGN_D_klb");
  #ifndef NO_AD_INITIALIZE
    wgt_rGN_D_klb.initialize();
  #endif
  lenprob.allocate(1,nages,1,nlenbins,"lenprob");
  #ifndef NO_AD_INITIALIZE
    lenprob.initialize();
  #endif
  lenprob_plus.allocate(1,nages,1,nlenbins_plus,"lenprob_plus");
  #ifndef NO_AD_INITIALIZE
    lenprob_plus.initialize();
  #endif
  lenprob_all.allocate(1,nages,1,nlenbins_all,"lenprob_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_all.initialize();
  #endif
  lenbins_all.allocate(1,nlenbins_all,"lenbins_all");
  #ifndef NO_AD_INITIALIZE
    lenbins_all.initialize();
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
  lenprob_sBT.allocate(1,nages,1,nlenbins,"lenprob_sBT");
  #ifndef NO_AD_INITIALIZE
    lenprob_sBT.initialize();
  #endif
  lenprob_cHL1.allocate(1,nages,1,nlenbins,"lenprob_cHL1");
  #ifndef NO_AD_INITIALIZE
    lenprob_cHL1.initialize();
  #endif
  lenprob_cHL2.allocate(1,nages,1,nlenbins,"lenprob_cHL2");
  #ifndef NO_AD_INITIALIZE
    lenprob_cHL2.initialize();
  #endif
  lenprob_cHL3.allocate(1,nages,1,nlenbins,"lenprob_cHL3");
  #ifndef NO_AD_INITIALIZE
    lenprob_cHL3.initialize();
  #endif
  lenprob_cPT1.allocate(1,nages,1,nlenbins,"lenprob_cPT1");
  #ifndef NO_AD_INITIALIZE
    lenprob_cPT1.initialize();
  #endif
  lenprob_cPT2.allocate(1,nages,1,nlenbins,"lenprob_cPT2");
  #ifndef NO_AD_INITIALIZE
    lenprob_cPT2.initialize();
  #endif
  lenprob_cPT3.allocate(1,nages,1,nlenbins,"lenprob_cPT3");
  #ifndef NO_AD_INITIALIZE
    lenprob_cPT3.initialize();
  #endif
  lenprob_cTW1.allocate(1,nages,1,nlenbins,"lenprob_cTW1");
  #ifndef NO_AD_INITIALIZE
    lenprob_cTW1.initialize();
  #endif
  lenprob_cTW2.allocate(1,nages,1,nlenbins,"lenprob_cTW2");
  #ifndef NO_AD_INITIALIZE
    lenprob_cTW2.initialize();
  #endif
  lenprob_rHB1.allocate(1,nages,1,nlenbins,"lenprob_rHB1");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB1.initialize();
  #endif
  lenprob_rHB2.allocate(1,nages,1,nlenbins,"lenprob_rHB2");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB2.initialize();
  #endif
  lenprob_rHB3.allocate(1,nages,1,nlenbins,"lenprob_rHB3");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB3.initialize();
  #endif
  lenprob_rHB4.allocate(1,nages,1,nlenbins,"lenprob_rHB4");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB4.initialize();
  #endif
  lenprob_rHB5.allocate(1,nages,1,nlenbins,"lenprob_rHB5");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB5.initialize();
  #endif
  lenprob_rGN1.allocate(1,nages,1,nlenbins,"lenprob_rGN1");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN1.initialize();
  #endif
  lenprob_rGN2.allocate(1,nages,1,nlenbins,"lenprob_rGN2");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN2.initialize();
  #endif
  lenprob_rGN3.allocate(1,nages,1,nlenbins,"lenprob_rGN3");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN3.initialize();
  #endif
  lenprob_rGN4.allocate(1,nages,1,nlenbins,"lenprob_rGN4");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN4.initialize();
  #endif
  lenprob_rGN5.allocate(1,nages,1,nlenbins,"lenprob_rGN5");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN5.initialize();
  #endif
  lenprob_cGN_D2.allocate(1,nages,1,nlenbins,"lenprob_cGN_D2");
  #ifndef NO_AD_INITIALIZE
    lenprob_cGN_D2.initialize();
  #endif
  lenprob_cGN_D3.allocate(1,nages,1,nlenbins,"lenprob_cGN_D3");
  #ifndef NO_AD_INITIALIZE
    lenprob_cGN_D3.initialize();
  #endif
  lenprob_cGN_D4.allocate(1,nages,1,nlenbins,"lenprob_cGN_D4");
  #ifndef NO_AD_INITIALIZE
    lenprob_cGN_D4.initialize();
  #endif
  lenprob_rHB_D2.allocate(1,nages,1,nlenbins,"lenprob_rHB_D2");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB_D2.initialize();
  #endif
  lenprob_rHB_D3.allocate(1,nages,1,nlenbins,"lenprob_rHB_D3");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB_D3.initialize();
  #endif
  lenprob_rHB_D4.allocate(1,nages,1,nlenbins,"lenprob_rHB_D4");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB_D4.initialize();
  #endif
  lenprob_rHB_D5.allocate(1,nages,1,nlenbins,"lenprob_rHB_D5");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB_D5.initialize();
  #endif
  lenprob_rGN_D1.allocate(1,nages,1,nlenbins,"lenprob_rGN_D1");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN_D1.initialize();
  #endif
  lenprob_rGN_D2.allocate(1,nages,1,nlenbins,"lenprob_rGN_D2");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN_D2.initialize();
  #endif
  lenprob_rGN_D3.allocate(1,nages,1,nlenbins,"lenprob_rGN_D3");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN_D3.initialize();
  #endif
  lenprob_rGN_D4.allocate(1,nages,1,nlenbins,"lenprob_rGN_D4");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN_D4.initialize();
  #endif
  lenprob_rGN_D5.allocate(1,nages,1,nlenbins,"lenprob_rGN_D5");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN_D5.initialize();
  #endif
  lenprob_cHL1_all.allocate(1,nages,1,nlenbins_all,"lenprob_cHL1_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_cHL1_all.initialize();
  #endif
  lenprob_cHL2_all.allocate(1,nages,1,nlenbins_all,"lenprob_cHL2_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_cHL2_all.initialize();
  #endif
  lenprob_cHL3_all.allocate(1,nages,1,nlenbins_all,"lenprob_cHL3_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_cHL3_all.initialize();
  #endif
  lenprob_cPT1_all.allocate(1,nages,1,nlenbins_all,"lenprob_cPT1_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_cPT1_all.initialize();
  #endif
  lenprob_cPT2_all.allocate(1,nages,1,nlenbins_all,"lenprob_cPT2_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_cPT2_all.initialize();
  #endif
  lenprob_cPT3_all.allocate(1,nages,1,nlenbins_all,"lenprob_cPT3_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_cPT3_all.initialize();
  #endif
  lenprob_cTW1_all.allocate(1,nages,1,nlenbins_all,"lenprob_cTW1_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_cTW1_all.initialize();
  #endif
  lenprob_cTW2_all.allocate(1,nages,1,nlenbins_all,"lenprob_cTW2_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_cTW2_all.initialize();
  #endif
  lenprob_rHB1_all.allocate(1,nages,1,nlenbins_all,"lenprob_rHB1_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB1_all.initialize();
  #endif
  lenprob_rHB2_all.allocate(1,nages,1,nlenbins_all,"lenprob_rHB2_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB2_all.initialize();
  #endif
  lenprob_rHB3_all.allocate(1,nages,1,nlenbins_all,"lenprob_rHB3_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB3_all.initialize();
  #endif
  lenprob_rHB4_all.allocate(1,nages,1,nlenbins_all,"lenprob_rHB4_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB4_all.initialize();
  #endif
  lenprob_rHB5_all.allocate(1,nages,1,nlenbins_all,"lenprob_rHB5_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB5_all.initialize();
  #endif
  lenprob_rGN1_all.allocate(1,nages,1,nlenbins_all,"lenprob_rGN1_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN1_all.initialize();
  #endif
  lenprob_rGN2_all.allocate(1,nages,1,nlenbins_all,"lenprob_rGN2_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN2_all.initialize();
  #endif
  lenprob_rGN3_all.allocate(1,nages,1,nlenbins_all,"lenprob_rGN3_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN3_all.initialize();
  #endif
  lenprob_rGN4_all.allocate(1,nages,1,nlenbins_all,"lenprob_rGN4_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN4_all.initialize();
  #endif
  lenprob_rGN5_all.allocate(1,nages,1,nlenbins_all,"lenprob_rGN5_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN5_all.initialize();
  #endif
  lenprob_cGN_D2_all.allocate(1,nages,1,nlenbins_all,"lenprob_cGN_D2_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_cGN_D2_all.initialize();
  #endif
  lenprob_cGN_D3_all.allocate(1,nages,1,nlenbins_all,"lenprob_cGN_D3_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_cGN_D3_all.initialize();
  #endif
  lenprob_cGN_D4_all.allocate(1,nages,1,nlenbins_all,"lenprob_cGN_D4_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_cGN_D4_all.initialize();
  #endif
  lenprob_rHB_D2_all.allocate(1,nages,1,nlenbins_all,"lenprob_rHB_D2_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB_D2_all.initialize();
  #endif
  lenprob_rHB_D3_all.allocate(1,nages,1,nlenbins_all,"lenprob_rHB_D3_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB_D3_all.initialize();
  #endif
  lenprob_rHB_D4_all.allocate(1,nages,1,nlenbins_all,"lenprob_rHB_D4_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB_D4_all.initialize();
  #endif
  lenprob_rHB_D5_all.allocate(1,nages,1,nlenbins_all,"lenprob_rHB_D5_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB_D5_all.initialize();
  #endif
  lenprob_rGN_D1_all.allocate(1,nages,1,nlenbins_all,"lenprob_rGN_D1_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN_D1_all.initialize();
  #endif
  lenprob_rGN_D2_all.allocate(1,nages,1,nlenbins_all,"lenprob_rGN_D2_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN_D2_all.initialize();
  #endif
  lenprob_rGN_D3_all.allocate(1,nages,1,nlenbins_all,"lenprob_rGN_D3_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN_D3_all.initialize();
  #endif
  lenprob_rGN_D4_all.allocate(1,nages,1,nlenbins_all,"lenprob_rGN_D4_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN_D4_all.initialize();
  #endif
  lenprob_rGN_D5_all.allocate(1,nages,1,nlenbins_all,"lenprob_rGN_D5_all");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN_D5_all.initialize();
  #endif
  len_sd.allocate(1,nages,"len_sd");
  #ifndef NO_AD_INITIALIZE
    len_sd.initialize();
  #endif
  len_cv.allocate(1,nages,"len_cv");
  #ifndef NO_AD_INITIALIZE
    len_cv.initialize();
  #endif
  pred_sBT_lenc.allocate(1,nyr_lenc_sBT,1,nlenbins,"pred_sBT_lenc");
  #ifndef NO_AD_INITIALIZE
    pred_sBT_lenc.initialize();
  #endif
  pred_cHL_lenc.allocate(1,nyr_lenc_cHL,1,nlenbins,"pred_cHL_lenc");
  #ifndef NO_AD_INITIALIZE
    pred_cHL_lenc.initialize();
  #endif
  pred_cPT_lenc.allocate(1,nyr_lenc_cPT,1,nlenbins,"pred_cPT_lenc");
  #ifndef NO_AD_INITIALIZE
    pred_cPT_lenc.initialize();
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
  L_cPT_num_pool.allocate(1,nyr_lenc_cPT,1,nages,"L_cPT_num_pool");
  #ifndef NO_AD_INITIALIZE
    L_cPT_num_pool.initialize();
  #endif
  L_cPT_num_pool_yr.allocate(1,nyr_lenc_pool_cPT,1,nages,"L_cPT_num_pool_yr");
  #ifndef NO_AD_INITIALIZE
    L_cPT_num_pool_yr.initialize();
  #endif
  p_lenc_cHL2.allocate("p_lenc_cHL2");
  #ifndef NO_AD_INITIALIZE
  p_lenc_cHL2.initialize();
  #endif
  p_lenc_cHL3.allocate("p_lenc_cHL3");
  #ifndef NO_AD_INITIALIZE
  p_lenc_cHL3.initialize();
  #endif
  p_lenc_cPT2.allocate("p_lenc_cPT2");
  #ifndef NO_AD_INITIALIZE
  p_lenc_cPT2.initialize();
  #endif
  p_lenc_cPT3.allocate("p_lenc_cPT3");
  #ifndef NO_AD_INITIALIZE
  p_lenc_cPT3.initialize();
  #endif
  p_lenc_cTW2.allocate("p_lenc_cTW2");
  #ifndef NO_AD_INITIALIZE
  p_lenc_cTW2.initialize();
  #endif
  p_lenc_cTW3.allocate("p_lenc_cTW3");
  #ifndef NO_AD_INITIALIZE
  p_lenc_cTW3.initialize();
  #endif
  p_lenc_rHB2.allocate("p_lenc_rHB2");
  #ifndef NO_AD_INITIALIZE
  p_lenc_rHB2.initialize();
  #endif
  p_lenc_rHB3.allocate("p_lenc_rHB3");
  #ifndef NO_AD_INITIALIZE
  p_lenc_rHB3.initialize();
  #endif
  p_lenc_rHB4.allocate("p_lenc_rHB4");
  #ifndef NO_AD_INITIALIZE
  p_lenc_rHB4.initialize();
  #endif
  p_lenc_rHB5.allocate("p_lenc_rHB5");
  #ifndef NO_AD_INITIALIZE
  p_lenc_rHB5.initialize();
  #endif
  p_lenc_rGN2.allocate("p_lenc_rGN2");
  #ifndef NO_AD_INITIALIZE
  p_lenc_rGN2.initialize();
  #endif
  p_lenc_rGN3.allocate("p_lenc_rGN3");
  #ifndef NO_AD_INITIALIZE
  p_lenc_rGN3.initialize();
  #endif
  p_lenc_rGN4.allocate("p_lenc_rGN4");
  #ifndef NO_AD_INITIALIZE
  p_lenc_rGN4.initialize();
  #endif
  p_lenc_rGN5.allocate("p_lenc_rGN5");
  #ifndef NO_AD_INITIALIZE
  p_lenc_rGN5.initialize();
  #endif
  p_lenc_cGN_D2.allocate("p_lenc_cGN_D2");
  #ifndef NO_AD_INITIALIZE
  p_lenc_cGN_D2.initialize();
  #endif
  p_lenc_cGN_D3.allocate("p_lenc_cGN_D3");
  #ifndef NO_AD_INITIALIZE
  p_lenc_cGN_D3.initialize();
  #endif
  p_lenc_cGN_D4.allocate("p_lenc_cGN_D4");
  #ifndef NO_AD_INITIALIZE
  p_lenc_cGN_D4.initialize();
  #endif
  p_lenc_rHB_D2.allocate("p_lenc_rHB_D2");
  #ifndef NO_AD_INITIALIZE
  p_lenc_rHB_D2.initialize();
  #endif
  p_lenc_rHB_D3.allocate("p_lenc_rHB_D3");
  #ifndef NO_AD_INITIALIZE
  p_lenc_rHB_D3.initialize();
  #endif
  p_lenc_rHB_D4.allocate("p_lenc_rHB_D4");
  #ifndef NO_AD_INITIALIZE
  p_lenc_rHB_D4.initialize();
  #endif
  p_lenc_rHB_D5.allocate("p_lenc_rHB_D5");
  #ifndef NO_AD_INITIALIZE
  p_lenc_rHB_D5.initialize();
  #endif
  p_lenc_rGN_D1.allocate("p_lenc_rGN_D1");
  #ifndef NO_AD_INITIALIZE
  p_lenc_rGN_D1.initialize();
  #endif
  p_lenc_rGN_D2.allocate("p_lenc_rGN_D2");
  #ifndef NO_AD_INITIALIZE
  p_lenc_rGN_D2.initialize();
  #endif
  p_lenc_rGN_D3.allocate("p_lenc_rGN_D3");
  #ifndef NO_AD_INITIALIZE
  p_lenc_rGN_D3.initialize();
  #endif
  p_lenc_rGN_D4.allocate("p_lenc_rGN_D4");
  #ifndef NO_AD_INITIALIZE
  p_lenc_rGN_D4.initialize();
  #endif
  p_lenc_rGN_D5.allocate("p_lenc_rGN_D5");
  #ifndef NO_AD_INITIALIZE
  p_lenc_rGN_D5.initialize();
  #endif
  pred_sBT_agec.allocate(1,nyr_agec_sBT,1,nages,"pred_sBT_agec");
  #ifndef NO_AD_INITIALIZE
    pred_sBT_agec.initialize();
  #endif
  ErrorFree_sBT_agec.allocate(1,nyr_agec_sBT,1,nages,"ErrorFree_sBT_agec");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_sBT_agec.initialize();
  #endif
  pred_sTV_agec.allocate(1,nyr_agec_sTV,1,nages,"pred_sTV_agec");
  #ifndef NO_AD_INITIALIZE
    pred_sTV_agec.initialize();
  #endif
  ErrorFree_sTV_agec.allocate(1,nyr_agec_sTV,1,nages,"ErrorFree_sTV_agec");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_sTV_agec.initialize();
  #endif
  pred_cHL_agec.allocate(1,nyr_agec_cHL,1,nages,"pred_cHL_agec");
  #ifndef NO_AD_INITIALIZE
    pred_cHL_agec.initialize();
  #endif
  ErrorFree_cHL_agec.allocate(1,nyr_agec_cHL,1,nages,"ErrorFree_cHL_agec");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_cHL_agec.initialize();
  #endif
  pred_cPT_agec.allocate(1,nyr_agec_cPT,1,nages,"pred_cPT_agec");
  #ifndef NO_AD_INITIALIZE
    pred_cPT_agec.initialize();
  #endif
  ErrorFree_cPT_agec.allocate(1,nyr_agec_cPT,1,nages,"ErrorFree_cPT_agec");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_cPT_agec.initialize();
  #endif
  pred_rHB_agec.allocate(1,nyr_agec_rHB,1,nages,"pred_rHB_agec");
  #ifndef NO_AD_INITIALIZE
    pred_rHB_agec.initialize();
  #endif
  ErrorFree_rHB_agec.allocate(1,nyr_agec_rHB,1,nages,"ErrorFree_rHB_agec");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_rHB_agec.initialize();
  #endif
  nsamp_sBT_lenc_allyr.allocate(styr,endyr,"nsamp_sBT_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_sBT_lenc_allyr.initialize();
  #endif
  nsamp_cHL_lenc_allyr.allocate(styr,endyr,"nsamp_cHL_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_cHL_lenc_allyr.initialize();
  #endif
  nsamp_cPT_lenc_allyr.allocate(styr,endyr,"nsamp_cPT_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_cPT_lenc_allyr.initialize();
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
  nsamp_sBT_agec_allyr.allocate(styr,endyr,"nsamp_sBT_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_sBT_agec_allyr.initialize();
  #endif
  nsamp_sTV_agec_allyr.allocate(styr,endyr,"nsamp_sTV_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_sTV_agec_allyr.initialize();
  #endif
  nsamp_cHL_agec_allyr.allocate(styr,endyr,"nsamp_cHL_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_cHL_agec_allyr.initialize();
  #endif
  nsamp_cPT_agec_allyr.allocate(styr,endyr,"nsamp_cPT_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_cPT_agec_allyr.initialize();
  #endif
  nsamp_rHB_agec_allyr.allocate(styr,endyr,"nsamp_rHB_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_rHB_agec_allyr.initialize();
  #endif
  nfish_sBT_lenc_allyr.allocate(styr,endyr,"nfish_sBT_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_sBT_lenc_allyr.initialize();
  #endif
  nfish_cHL_lenc_allyr.allocate(styr,endyr,"nfish_cHL_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_cHL_lenc_allyr.initialize();
  #endif
  nfish_cPT_lenc_allyr.allocate(styr,endyr,"nfish_cPT_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_cPT_lenc_allyr.initialize();
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
  nfish_sBT_agec_allyr.allocate(styr,endyr,"nfish_sBT_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_sBT_agec_allyr.initialize();
  #endif
  nfish_sTV_agec_allyr.allocate(styr,endyr,"nfish_sTV_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_sTV_agec_allyr.initialize();
  #endif
  nfish_cHL_agec_allyr.allocate(styr,endyr,"nfish_cHL_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_cHL_agec_allyr.initialize();
  #endif
  nfish_cPT_agec_allyr.allocate(styr,endyr,"nfish_cPT_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_cPT_agec_allyr.initialize();
  #endif
  nfish_rHB_agec_allyr.allocate(styr,endyr,"nfish_rHB_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_rHB_agec_allyr.initialize();
  #endif
  neff_sBT_lenc_allyr_out.allocate(styr,endyr,"neff_sBT_lenc_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_sBT_lenc_allyr_out.initialize();
  #endif
  neff_cHL_lenc_allyr_out.allocate(styr,endyr,"neff_cHL_lenc_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_cHL_lenc_allyr_out.initialize();
  #endif
  neff_cPT_lenc_allyr_out.allocate(styr,endyr,"neff_cPT_lenc_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_cPT_lenc_allyr_out.initialize();
  #endif
  neff_rHB_lenc_allyr_out.allocate(styr,endyr,"neff_rHB_lenc_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_rHB_lenc_allyr_out.initialize();
  #endif
  neff_rHB_D_lenc_allyr_out.allocate(styr,endyr,"neff_rHB_D_lenc_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_rHB_D_lenc_allyr_out.initialize();
  #endif
  neff_rGN_lenc_allyr_out.allocate(styr,endyr,"neff_rGN_lenc_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_rGN_lenc_allyr_out.initialize();
  #endif
  neff_sBT_agec_allyr_out.allocate(styr,endyr,"neff_sBT_agec_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_sBT_agec_allyr_out.initialize();
  #endif
  neff_sTV_agec_allyr_out.allocate(styr,endyr,"neff_sTV_agec_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_sTV_agec_allyr_out.initialize();
  #endif
  neff_cHL_agec_allyr_out.allocate(styr,endyr,"neff_cHL_agec_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_cHL_agec_allyr_out.initialize();
  #endif
  neff_cPT_agec_allyr_out.allocate(styr,endyr,"neff_cPT_agec_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_cPT_agec_allyr_out.initialize();
  #endif
  neff_rHB_agec_allyr_out.allocate(styr,endyr,"neff_rHB_agec_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_rHB_agec_allyr_out.initialize();
  #endif
  N.allocate(styr,endyr,1,nages,"N");
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
  B.allocate(styr,endyr,1,nages,"B");
  #ifndef NO_AD_INITIALIZE
    B.initialize();
  #endif
  totB.allocate(styr,endyr,"totB");
  #ifndef NO_AD_INITIALIZE
    totB.initialize();
  #endif
  totN.allocate(styr,endyr,"totN");
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
  rec.allocate(styr,endyr,"rec");
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
  rec_sigma_sqd2.allocate("rec_sigma_sqd2");
  #ifndef NO_AD_INITIALIZE
  rec_sigma_sqd2.initialize();
  #endif
  rec_logL_add.allocate("rec_logL_add");
  #ifndef NO_AD_INITIALIZE
  rec_logL_add.initialize();
  #endif
  log_dev_rec.allocate(styr_rec_dev,endyr_rec_phase2,log_rec_dev_LO,log_rec_dev_HI,log_rec_dev_PH,"log_dev_rec");
  log_rec_dev_output.allocate(styr_rec_dev,endyr,"log_rec_dev_output");
  #ifndef NO_AD_INITIALIZE
    log_rec_dev_output.initialize();
  #endif
  log_rec_dev_out.allocate(styr_rec_dev,endyr_rec_phase2,"log_rec_dev_out");
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
  log_dm_lenc_sBT.allocate(log_dm_sBT_lc_LO,log_dm_sBT_lc_HI,log_dm_sBT_lc_PH,"log_dm_lenc_sBT");
  log_dm_lenc_sTV.allocate(log_dm_sTV_lc_LO,log_dm_sTV_lc_HI,log_dm_sTV_lc_PH,"log_dm_lenc_sTV");
  log_dm_lenc_cHL.allocate(log_dm_cHL_lc_LO,log_dm_cHL_lc_HI,log_dm_cHL_lc_PH,"log_dm_lenc_cHL");
  log_dm_lenc_cPT.allocate(log_dm_cPT_lc_LO,log_dm_cPT_lc_HI,log_dm_cPT_lc_PH,"log_dm_lenc_cPT");
  log_dm_lenc_rHB.allocate(log_dm_rHB_lc_LO,log_dm_rHB_lc_HI,log_dm_rHB_lc_PH,"log_dm_lenc_rHB");
  log_dm_lenc_rHB_D.allocate(log_dm_rHB_D_lc_LO,log_dm_rHB_D_lc_HI,log_dm_rHB_D_lc_PH,"log_dm_lenc_rHB_D");
  log_dm_lenc_rGN.allocate(log_dm_rGN_lc_LO,log_dm_rGN_lc_HI,log_dm_rGN_lc_PH,"log_dm_lenc_rGN");
  log_dm_agec_sBT.allocate(log_dm_sBT_ac_LO,log_dm_sBT_ac_HI,log_dm_sBT_ac_PH,"log_dm_agec_sBT");
  log_dm_agec_sTV.allocate(log_dm_sTV_ac_LO,log_dm_sTV_ac_HI,log_dm_sTV_ac_PH,"log_dm_agec_sTV");
  log_dm_agec_cHL.allocate(log_dm_cHL_ac_LO,log_dm_cHL_ac_HI,log_dm_cHL_ac_PH,"log_dm_agec_cHL");
  log_dm_agec_cPT.allocate(log_dm_cPT_ac_LO,log_dm_cPT_ac_HI,log_dm_cPT_ac_PH,"log_dm_agec_cPT");
  log_dm_agec_rHB.allocate(log_dm_rHB_ac_LO,log_dm_rHB_ac_HI,log_dm_rHB_ac_PH,"log_dm_agec_rHB");
  log_dm_sBT_lc_out.allocate(1,8,"log_dm_sBT_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_sBT_lc_out.initialize();
  #endif
  log_dm_sTV_lc_out.allocate(1,8,"log_dm_sTV_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_sTV_lc_out.initialize();
  #endif
  log_dm_cHL_lc_out.allocate(1,8,"log_dm_cHL_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_cHL_lc_out.initialize();
  #endif
  log_dm_cPT_lc_out.allocate(1,8,"log_dm_cPT_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_cPT_lc_out.initialize();
  #endif
  log_dm_rHB_lc_out.allocate(1,8,"log_dm_rHB_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_rHB_lc_out.initialize();
  #endif
  log_dm_rGN_lc_out.allocate(1,8,"log_dm_rGN_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_rGN_lc_out.initialize();
  #endif
  log_dm_rHB_D_lc_out.allocate(1,8,"log_dm_rHB_D_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_rHB_D_lc_out.initialize();
  #endif
  log_dm_sBT_ac_out.allocate(1,8,"log_dm_sBT_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_sBT_ac_out.initialize();
  #endif
  log_dm_sTV_ac_out.allocate(1,8,"log_dm_sTV_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_sTV_ac_out.initialize();
  #endif
  log_dm_cHL_ac_out.allocate(1,8,"log_dm_cHL_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_cHL_ac_out.initialize();
  #endif
  log_dm_cPT_ac_out.allocate(1,8,"log_dm_cPT_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_cPT_ac_out.initialize();
  #endif
  log_dm_rHB_ac_out.allocate(1,8,"log_dm_rHB_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_rHB_ac_out.initialize();
  #endif
  sel_sBT.allocate(styr,endyr,1,nages,"sel_sBT");
  #ifndef NO_AD_INITIALIZE
    sel_sBT.initialize();
  #endif
  sel_sBT_vec.allocate(1,nages,"sel_sBT_vec");
  #ifndef NO_AD_INITIALIZE
    sel_sBT_vec.initialize();
  #endif
  selpar_A50_sBT.allocate(selpar_A50_sBT_LO,selpar_A50_sBT_HI,selpar_A50_sBT_PH,"selpar_A50_sBT");
  selpar_slope_sBT.allocate(selpar_slope_sBT_LO,selpar_slope_sBT_HI,selpar_slope_sBT_PH,"selpar_slope_sBT");
  selpar_A50_sBT_out.allocate(1,8,"selpar_A50_sBT_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_sBT_out.initialize();
  #endif
  selpar_slope_sBT_out.allocate(1,8,"selpar_slope_sBT_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_sBT_out.initialize();
  #endif
  sel_sTV.allocate(styr,endyr,1,nages,"sel_sTV");
  #ifndef NO_AD_INITIALIZE
    sel_sTV.initialize();
  #endif
  sel_sTV_vec.allocate(1,nages,"sel_sTV_vec");
  #ifndef NO_AD_INITIALIZE
    sel_sTV_vec.initialize();
  #endif
  selpar_A50_sTV.allocate(selpar_A50_sTV_LO,selpar_A50_sTV_HI,selpar_A50_sTV_PH,"selpar_A50_sTV");
  selpar_slope_sTV.allocate(selpar_slope_sTV_LO,selpar_slope_sTV_HI,selpar_slope_sTV_PH,"selpar_slope_sTV");
  selpar_A50_sTV_out.allocate(1,8,"selpar_A50_sTV_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_sTV_out.initialize();
  #endif
  selpar_slope_sTV_out.allocate(1,8,"selpar_slope_sTV_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_sTV_out.initialize();
  #endif
  sel_cHL.allocate(styr,endyr,1,nages,"sel_cHL");
  #ifndef NO_AD_INITIALIZE
    sel_cHL.initialize();
  #endif
  sel_cHL_2.allocate(1,nages,"sel_cHL_2");
  #ifndef NO_AD_INITIALIZE
    sel_cHL_2.initialize();
  #endif
  sel_cHL_3.allocate(1,nages,"sel_cHL_3");
  #ifndef NO_AD_INITIALIZE
    sel_cHL_3.initialize();
  #endif
  sel_cHL_4.allocate(1,nages,"sel_cHL_4");
  #ifndef NO_AD_INITIALIZE
    sel_cHL_4.initialize();
  #endif
  selpar_A50_cHL2.allocate(selpar_A50_cHL2_LO,selpar_A50_cHL2_HI,selpar_A50_cHL2_PH,"selpar_A50_cHL2");
  selpar_slope_cHL2.allocate(selpar_slope_cHL2_LO,selpar_slope_cHL2_HI,selpar_slope_cHL2_PH,"selpar_slope_cHL2");
  selpar_A50_cHL3.allocate(selpar_A50_cHL3_LO,selpar_A50_cHL3_HI,selpar_A50_cHL3_PH,"selpar_A50_cHL3");
  selpar_slope_cHL3.allocate(selpar_slope_cHL3_LO,selpar_slope_cHL3_HI,selpar_slope_cHL3_PH,"selpar_slope_cHL3");
  selpar_A50_cHL4.allocate(selpar_A50_cHL4_LO,selpar_A50_cHL4_HI,selpar_A50_cHL4_PH,"selpar_A50_cHL4");
  selpar_slope_cHL4.allocate(selpar_slope_cHL4_LO,selpar_slope_cHL4_HI,selpar_slope_cHL4_PH,"selpar_slope_cHL4");
  selpar_A50_cHL2_out.allocate(1,8,"selpar_A50_cHL2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cHL2_out.initialize();
  #endif
  selpar_slope_cHL2_out.allocate(1,8,"selpar_slope_cHL2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cHL2_out.initialize();
  #endif
  selpar_A50_cHL3_out.allocate(1,8,"selpar_A50_cHL3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cHL3_out.initialize();
  #endif
  selpar_slope_cHL3_out.allocate(1,8,"selpar_slope_cHL3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cHL3_out.initialize();
  #endif
  selpar_A50_cHL4_out.allocate(1,8,"selpar_A50_cHL4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cHL4_out.initialize();
  #endif
  selpar_slope_cHL4_out.allocate(1,8,"selpar_slope_cHL4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cHL4_out.initialize();
  #endif
  sel_cGN_D.allocate(styr,endyr,1,nages,"sel_cGN_D");
  #ifndef NO_AD_INITIALIZE
    sel_cGN_D.initialize();
  #endif
  sel_cGN_D_2.allocate(1,nages,"sel_cGN_D_2");
  #ifndef NO_AD_INITIALIZE
    sel_cGN_D_2.initialize();
  #endif
  sel_cGN_D_3.allocate(1,nages,"sel_cGN_D_3");
  #ifndef NO_AD_INITIALIZE
    sel_cGN_D_3.initialize();
  #endif
  sel_cGN_D_4.allocate(1,nages,"sel_cGN_D_4");
  #ifndef NO_AD_INITIALIZE
    sel_cGN_D_4.initialize();
  #endif
  sel_cGN_D_quota3.allocate(1,nages,"sel_cGN_D_quota3");
  #ifndef NO_AD_INITIALIZE
    sel_cGN_D_quota3.initialize();
  #endif
  Dopen_cHL.allocate("Dopen_cHL");
  #ifndef NO_AD_INITIALIZE
  Dopen_cHL.initialize();
  #endif
  Dclosed_cHL.allocate("Dclosed_cHL");
  #ifndef NO_AD_INITIALIZE
  Dclosed_cHL.initialize();
  #endif
  Lopen_cHL.allocate("Lopen_cHL");
  #ifndef NO_AD_INITIALIZE
  Lopen_cHL.initialize();
  #endif
  Dopen_cPT.allocate("Dopen_cPT");
  #ifndef NO_AD_INITIALIZE
  Dopen_cPT.initialize();
  #endif
  Dclosed_cPT.allocate("Dclosed_cPT");
  #ifndef NO_AD_INITIALIZE
  Dclosed_cPT.initialize();
  #endif
  Lopen_cPT.allocate("Lopen_cPT");
  #ifndef NO_AD_INITIALIZE
  Lopen_cPT.initialize();
  #endif
  D_sum_cHLcPT.allocate("D_sum_cHLcPT");
  #ifndef NO_AD_INITIALIZE
  D_sum_cHLcPT.initialize();
  #endif
  Dprop_cGN_sel_D.allocate("Dprop_cGN_sel_D");
  #ifndef NO_AD_INITIALIZE
  Dprop_cGN_sel_D.initialize();
  #endif
  Dprop_cGN_sel_cHL.allocate("Dprop_cGN_sel_cHL");
  #ifndef NO_AD_INITIALIZE
  Dprop_cGN_sel_cHL.initialize();
  #endif
  Dprop_cGN_sel_cPT.allocate("Dprop_cGN_sel_cPT");
  #ifndef NO_AD_INITIALIZE
  Dprop_cGN_sel_cPT.initialize();
  #endif
  sel_cPT.allocate(styr,endyr,1,nages,"sel_cPT");
  #ifndef NO_AD_INITIALIZE
    sel_cPT.initialize();
  #endif
  sel_cPT_2.allocate(1,nages,"sel_cPT_2");
  #ifndef NO_AD_INITIALIZE
    sel_cPT_2.initialize();
  #endif
  sel_cPT_3.allocate(1,nages,"sel_cPT_3");
  #ifndef NO_AD_INITIALIZE
    sel_cPT_3.initialize();
  #endif
  sel_cPT_4.allocate(1,nages,"sel_cPT_4");
  #ifndef NO_AD_INITIALIZE
    sel_cPT_4.initialize();
  #endif
  selpar_A50_cPT2.allocate(selpar_A50_cPT2_LO,selpar_A50_cPT2_HI,selpar_A50_cPT2_PH,"selpar_A50_cPT2");
  selpar_slope_cPT2.allocate(selpar_slope_cPT2_LO,selpar_slope_cPT2_HI,selpar_slope_cPT2_PH,"selpar_slope_cPT2");
  selpar_A50_cPT3.allocate(selpar_A50_cPT3_LO,selpar_A50_cPT3_HI,selpar_A50_cPT3_PH,"selpar_A50_cPT3");
  selpar_slope_cPT3.allocate(selpar_slope_cPT3_LO,selpar_slope_cPT3_HI,selpar_slope_cPT3_PH,"selpar_slope_cPT3");
  selpar_A50_cPT4.allocate(selpar_A50_cPT4_LO,selpar_A50_cPT4_HI,selpar_A50_cPT4_PH,"selpar_A50_cPT4");
  selpar_slope_cPT4.allocate(selpar_slope_cPT4_LO,selpar_slope_cPT4_HI,selpar_slope_cPT4_PH,"selpar_slope_cPT4");
  selpar_A50_cPT2_out.allocate(1,8,"selpar_A50_cPT2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cPT2_out.initialize();
  #endif
  selpar_slope_cPT2_out.allocate(1,8,"selpar_slope_cPT2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cPT2_out.initialize();
  #endif
  selpar_A50_cPT3_out.allocate(1,8,"selpar_A50_cPT3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cPT3_out.initialize();
  #endif
  selpar_slope_cPT3_out.allocate(1,8,"selpar_slope_cPT3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cPT3_out.initialize();
  #endif
  selpar_A50_cPT4_out.allocate(1,8,"selpar_A50_cPT4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cPT4_out.initialize();
  #endif
  selpar_slope_cPT4_out.allocate(1,8,"selpar_slope_cPT4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cPT4_out.initialize();
  #endif
  sel_cTW.allocate(styr,endyr,1,nages,"sel_cTW");
  #ifndef NO_AD_INITIALIZE
    sel_cTW.initialize();
  #endif
  sel_rHB.allocate(styr,endyr,1,nages,"sel_rHB");
  #ifndef NO_AD_INITIALIZE
    sel_rHB.initialize();
  #endif
  sel_rHB_1.allocate(1,nages,"sel_rHB_1");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_1.initialize();
  #endif
  sel_rHB_2.allocate(1,nages,"sel_rHB_2");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_2.initialize();
  #endif
  sel_rHB_3.allocate(1,nages,"sel_rHB_3");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_3.initialize();
  #endif
  sel_rHB_4.allocate(1,nages,"sel_rHB_4");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_4.initialize();
  #endif
  sel_rHB_5.allocate(1,nages,"sel_rHB_5");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_5.initialize();
  #endif
  selpar_A50_rHB1.allocate(selpar_A50_rHB1_LO,selpar_A50_rHB1_HI,selpar_A50_rHB1_PH,"selpar_A50_rHB1");
  selpar_slope_rHB1.allocate(selpar_slope_rHB1_LO,selpar_slope_rHB1_HI,selpar_slope_rHB1_PH,"selpar_slope_rHB1");
  selpar_A50_rHB2.allocate(selpar_A50_rHB2_LO,selpar_A50_rHB2_HI,selpar_A50_rHB2_PH,"selpar_A50_rHB2");
  selpar_slope_rHB2.allocate(selpar_slope_rHB2_LO,selpar_slope_rHB2_HI,selpar_slope_rHB2_PH,"selpar_slope_rHB2");
  selpar_A50_rHB3.allocate(selpar_A50_rHB3_LO,selpar_A50_rHB3_HI,selpar_A50_rHB3_PH,"selpar_A50_rHB3");
  selpar_slope_rHB3.allocate(selpar_slope_rHB3_LO,selpar_slope_rHB3_HI,selpar_slope_rHB3_PH,"selpar_slope_rHB3");
  selpar_A50_rHB4.allocate(selpar_A50_rHB4_LO,selpar_A50_rHB4_HI,selpar_A50_rHB4_PH,"selpar_A50_rHB4");
  selpar_slope_rHB4.allocate(selpar_slope_rHB4_LO,selpar_slope_rHB4_HI,selpar_slope_rHB4_PH,"selpar_slope_rHB4");
  selpar_A50_rHB5.allocate(selpar_A50_rHB5_LO,selpar_A50_rHB5_HI,selpar_A50_rHB5_PH,"selpar_A50_rHB5");
  selpar_slope_rHB5.allocate(selpar_slope_rHB5_LO,selpar_slope_rHB5_HI,selpar_slope_rHB5_PH,"selpar_slope_rHB5");
  selpar_A50_rHB1_out.allocate(1,8,"selpar_A50_rHB1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rHB1_out.initialize();
  #endif
  selpar_slope_rHB1_out.allocate(1,8,"selpar_slope_rHB1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rHB1_out.initialize();
  #endif
  selpar_A50_rHB2_out.allocate(1,8,"selpar_A50_rHB2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rHB2_out.initialize();
  #endif
  selpar_slope_rHB2_out.allocate(1,8,"selpar_slope_rHB2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rHB2_out.initialize();
  #endif
  selpar_A50_rHB3_out.allocate(1,8,"selpar_A50_rHB3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rHB3_out.initialize();
  #endif
  selpar_slope_rHB3_out.allocate(1,8,"selpar_slope_rHB3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rHB3_out.initialize();
  #endif
  selpar_A50_rHB4_out.allocate(1,8,"selpar_A50_rHB4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rHB4_out.initialize();
  #endif
  selpar_slope_rHB4_out.allocate(1,8,"selpar_slope_rHB4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rHB4_out.initialize();
  #endif
  selpar_A50_rHB5_out.allocate(1,8,"selpar_A50_rHB5_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rHB5_out.initialize();
  #endif
  selpar_slope_rHB5_out.allocate(1,8,"selpar_slope_rHB5_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rHB5_out.initialize();
  #endif
  selpar_Age0_rHB_D.allocate("selpar_Age0_rHB_D");
  #ifndef NO_AD_INITIALIZE
  selpar_Age0_rHB_D.initialize();
  #endif
  selpar_Age1_rHB_D.allocate("selpar_Age1_rHB_D");
  #ifndef NO_AD_INITIALIZE
  selpar_Age1_rHB_D.initialize();
  #endif
  selpar_Age2_rHB_D.allocate("selpar_Age2_rHB_D");
  #ifndef NO_AD_INITIALIZE
  selpar_Age2_rHB_D.initialize();
  #endif
  sel_rHB_D.allocate(styr,endyr,1,nages,"sel_rHB_D");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_D.initialize();
  #endif
  sel_rHB_D_1.allocate(1,nages,"sel_rHB_D_1");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_D_1.initialize();
  #endif
  sel_rHB_D_2.allocate(1,nages,"sel_rHB_D_2");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_D_2.initialize();
  #endif
  sel_rHB_D_3.allocate(1,nages,"sel_rHB_D_3");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_D_3.initialize();
  #endif
  sel_rHB_D_4.allocate(1,nages,"sel_rHB_D_4");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_D_4.initialize();
  #endif
  sel_rHB_D_5.allocate(1,nages,"sel_rHB_D_5");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_D_5.initialize();
  #endif
  selpar_A50_rHB_D4.allocate(selpar_A50_rHB_D4_LO,selpar_A50_rHB_D4_HI,selpar_A50_rHB_D4_PH,"selpar_A50_rHB_D4");
  selpar_slope_rHB_D4.allocate(selpar_slope_rHB_D4_LO,selpar_slope_rHB_D4_HI,selpar_slope_rHB_D4_PH,"selpar_slope_rHB_D4");
  selpar_A502_rHB_D4.allocate(selpar_A502_rHB_D4_LO,selpar_A502_rHB_D4_HI,selpar_A502_rHB_D4_PH,"selpar_A502_rHB_D4");
  selpar_slope2_rHB_D4.allocate(selpar_slope2_rHB_D4_LO,selpar_slope2_rHB_D4_HI,selpar_slope2_rHB_D4_PH,"selpar_slope2_rHB_D4");
  selpar_A50_rHB_D5.allocate(selpar_A50_rHB_D5_LO,selpar_A50_rHB_D5_HI,selpar_A50_rHB_D5_PH,"selpar_A50_rHB_D5");
  selpar_slope_rHB_D5.allocate(selpar_slope_rHB_D5_LO,selpar_slope_rHB_D5_HI,selpar_slope_rHB_D5_PH,"selpar_slope_rHB_D5");
  selpar_A502_rHB_D5.allocate(selpar_A502_rHB_D5_LO,selpar_A502_rHB_D5_HI,selpar_A502_rHB_D5_PH,"selpar_A502_rHB_D5");
  selpar_slope2_rHB_D5.allocate(selpar_slope2_rHB_D5_LO,selpar_slope2_rHB_D5_HI,selpar_slope2_rHB_D5_PH,"selpar_slope2_rHB_D5");
  selpar_A50_rHB_D4_out.allocate(1,8,"selpar_A50_rHB_D4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rHB_D4_out.initialize();
  #endif
  selpar_slope_rHB_D4_out.allocate(1,8,"selpar_slope_rHB_D4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rHB_D4_out.initialize();
  #endif
  selpar_A502_rHB_D4_out.allocate(1,8,"selpar_A502_rHB_D4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_rHB_D4_out.initialize();
  #endif
  selpar_slope2_rHB_D4_out.allocate(1,8,"selpar_slope2_rHB_D4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_rHB_D4_out.initialize();
  #endif
  selpar_A50_rHB_D5_out.allocate(1,8,"selpar_A50_rHB_D5_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rHB_D5_out.initialize();
  #endif
  selpar_slope_rHB_D5_out.allocate(1,8,"selpar_slope_rHB_D5_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rHB_D5_out.initialize();
  #endif
  selpar_A502_rHB_D5_out.allocate(1,8,"selpar_A502_rHB_D5_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_rHB_D5_out.initialize();
  #endif
  selpar_slope2_rHB_D5_out.allocate(1,8,"selpar_slope2_rHB_D5_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_rHB_D5_out.initialize();
  #endif
  vecprob_rHB_D2.allocate(4,nages,"vecprob_rHB_D2");
  #ifndef NO_AD_INITIALIZE
    vecprob_rHB_D2.initialize();
  #endif
  vecprob_rHB_D3.allocate(4,nages,"vecprob_rHB_D3");
  #ifndef NO_AD_INITIALIZE
    vecprob_rHB_D3.initialize();
  #endif
  vecprob_rHB_D4.allocate(4,nages,"vecprob_rHB_D4");
  #ifndef NO_AD_INITIALIZE
    vecprob_rHB_D4.initialize();
  #endif
  vecprob_rHB_D5.allocate(4,nages,"vecprob_rHB_D5");
  #ifndef NO_AD_INITIALIZE
    vecprob_rHB_D5.initialize();
  #endif
  selpar_logit_Age0_rHB_D.allocate(selpar_Age0_rHB_D_logit_LO,selpar_Age0_rHB_D_logit_HI,selpar_Age0_rHB_D_logit_PH,"selpar_logit_Age0_rHB_D");
  selpar_logit_Age1_rHB_D.allocate(selpar_Age1_rHB_D_logit_LO,selpar_Age1_rHB_D_logit_HI,selpar_Age1_rHB_D_logit_PH,"selpar_logit_Age1_rHB_D");
  selpar_logit_Age2_rHB_D.allocate(selpar_Age2_rHB_D_logit_LO,selpar_Age2_rHB_D_logit_HI,selpar_Age2_rHB_D_logit_PH,"selpar_logit_Age2_rHB_D");
  selpar_Age0_rHB_D_logit_out.allocate(1,8,"selpar_Age0_rHB_D_logit_out");
  #ifndef NO_AD_INITIALIZE
    selpar_Age0_rHB_D_logit_out.initialize();
  #endif
  selpar_Age1_rHB_D_logit_out.allocate(1,8,"selpar_Age1_rHB_D_logit_out");
  #ifndef NO_AD_INITIALIZE
    selpar_Age1_rHB_D_logit_out.initialize();
  #endif
  selpar_Age2_rHB_D_logit_out.allocate(1,8,"selpar_Age2_rHB_D_logit_out");
  #ifndef NO_AD_INITIALIZE
    selpar_Age2_rHB_D_logit_out.initialize();
  #endif
  prob_belowsizelim_block1.allocate(1,nages,"prob_belowsizelim_block1");
  #ifndef NO_AD_INITIALIZE
    prob_belowsizelim_block1.initialize();
  #endif
  prob_belowsizelim_block2.allocate(1,nages,"prob_belowsizelim_block2");
  #ifndef NO_AD_INITIALIZE
    prob_belowsizelim_block2.initialize();
  #endif
  prob_belowsizelim_block3.allocate(1,nages,"prob_belowsizelim_block3");
  #ifndef NO_AD_INITIALIZE
    prob_belowsizelim_block3.initialize();
  #endif
  prob_belowsizelim_block4.allocate(1,nages,"prob_belowsizelim_block4");
  #ifndef NO_AD_INITIALIZE
    prob_belowsizelim_block4.initialize();
  #endif
  prob_belowsizelim_block5.allocate(1,nages,"prob_belowsizelim_block5");
  #ifndef NO_AD_INITIALIZE
    prob_belowsizelim_block5.initialize();
  #endif
  zscore_lsizelim1.allocate("zscore_lsizelim1");
  #ifndef NO_AD_INITIALIZE
  zscore_lsizelim1.initialize();
  #endif
  zscore_lsizelim2.allocate("zscore_lsizelim2");
  #ifndef NO_AD_INITIALIZE
  zscore_lsizelim2.initialize();
  #endif
  zscore_lsizelim3.allocate("zscore_lsizelim3");
  #ifndef NO_AD_INITIALIZE
  zscore_lsizelim3.initialize();
  #endif
  zscore_lsizelim4.allocate("zscore_lsizelim4");
  #ifndef NO_AD_INITIALIZE
  zscore_lsizelim4.initialize();
  #endif
  zscore_lsizelim5.allocate("zscore_lsizelim5");
  #ifndef NO_AD_INITIALIZE
  zscore_lsizelim5.initialize();
  #endif
  cprob_lsizelim1.allocate("cprob_lsizelim1");
  #ifndef NO_AD_INITIALIZE
  cprob_lsizelim1.initialize();
  #endif
  cprob_lsizelim2.allocate("cprob_lsizelim2");
  #ifndef NO_AD_INITIALIZE
  cprob_lsizelim2.initialize();
  #endif
  cprob_lsizelim3.allocate("cprob_lsizelim3");
  #ifndef NO_AD_INITIALIZE
  cprob_lsizelim3.initialize();
  #endif
  cprob_lsizelim4.allocate("cprob_lsizelim4");
  #ifndef NO_AD_INITIALIZE
  cprob_lsizelim4.initialize();
  #endif
  cprob_lsizelim5.allocate("cprob_lsizelim5");
  #ifndef NO_AD_INITIALIZE
  cprob_lsizelim5.initialize();
  #endif
  sel_rGN.allocate(styr,endyr,1,nages,"sel_rGN");
  #ifndef NO_AD_INITIALIZE
    sel_rGN.initialize();
  #endif
  sel_rGN_D.allocate(styr,endyr,1,nages,"sel_rGN_D");
  #ifndef NO_AD_INITIALIZE
    sel_rGN_D.initialize();
  #endif
  sel_rGN1.allocate(1,nages,"sel_rGN1");
  #ifndef NO_AD_INITIALIZE
    sel_rGN1.initialize();
  #endif
  sel_rGN2.allocate(1,nages,"sel_rGN2");
  #ifndef NO_AD_INITIALIZE
    sel_rGN2.initialize();
  #endif
  sel_rGN3.allocate(1,nages,"sel_rGN3");
  #ifndef NO_AD_INITIALIZE
    sel_rGN3.initialize();
  #endif
  sel_rGN4.allocate(1,nages,"sel_rGN4");
  #ifndef NO_AD_INITIALIZE
    sel_rGN4.initialize();
  #endif
  sel_rGN5.allocate(1,nages,"sel_rGN5");
  #ifndef NO_AD_INITIALIZE
    sel_rGN5.initialize();
  #endif
  selpar_A50_rGN1.allocate(selpar_A50_rGN1_LO,selpar_A50_rGN1_HI,selpar_A50_rGN1_PH,"selpar_A50_rGN1");
  selpar_slope_rGN1.allocate(selpar_slope_rGN1_LO,selpar_slope_rGN1_HI,selpar_slope_rGN1_PH,"selpar_slope_rGN1");
  selpar_A50_rGN2.allocate(selpar_A50_rGN2_LO,selpar_A50_rGN2_HI,selpar_A50_rGN2_PH,"selpar_A50_rGN2");
  selpar_slope_rGN2.allocate(selpar_slope_rGN2_LO,selpar_slope_rGN2_HI,selpar_slope_rGN2_PH,"selpar_slope_rGN2");
  selpar_A50_rGN3.allocate(selpar_A50_rGN3_LO,selpar_A50_rGN3_HI,selpar_A50_rGN3_PH,"selpar_A50_rGN3");
  selpar_slope_rGN3.allocate(selpar_slope_rGN3_LO,selpar_slope_rGN3_HI,selpar_slope_rGN3_PH,"selpar_slope_rGN3");
  selpar_A50_rGN4.allocate(selpar_A50_rGN4_LO,selpar_A50_rGN4_HI,selpar_A50_rGN4_PH,"selpar_A50_rGN4");
  selpar_slope_rGN4.allocate(selpar_slope_rGN4_LO,selpar_slope_rGN4_HI,selpar_slope_rGN4_PH,"selpar_slope_rGN4");
  selpar_A50_rGN5.allocate(selpar_A50_rGN5_LO,selpar_A50_rGN5_HI,selpar_A50_rGN5_PH,"selpar_A50_rGN5");
  selpar_slope_rGN5.allocate(selpar_slope_rGN5_LO,selpar_slope_rGN5_HI,selpar_slope_rGN5_PH,"selpar_slope_rGN5");
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
  selpar_A50_rGN3_out.allocate(1,8,"selpar_A50_rGN3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rGN3_out.initialize();
  #endif
  selpar_slope_rGN3_out.allocate(1,8,"selpar_slope_rGN3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rGN3_out.initialize();
  #endif
  selpar_A50_rGN4_out.allocate(1,8,"selpar_A50_rGN4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rGN4_out.initialize();
  #endif
  selpar_slope_rGN4_out.allocate(1,8,"selpar_slope_rGN4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rGN4_out.initialize();
  #endif
  selpar_A50_rGN5_out.allocate(1,8,"selpar_A50_rGN5_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rGN5_out.initialize();
  #endif
  selpar_slope_rGN5_out.allocate(1,8,"selpar_slope_rGN5_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rGN5_out.initialize();
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
  pred_sBT_cpue.allocate(styr_cpue_sBT,endyr_cpue_sBT,"pred_sBT_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_sBT_cpue.initialize();
  #endif
  N_sBT.allocate(styr_cpue_sBT,endyr_cpue_sBT,1,nages,"N_sBT");
  #ifndef NO_AD_INITIALIZE
    N_sBT.initialize();
  #endif
  pred_sTV_cpue.allocate(styr_cpue_sTV,endyr_cpue_sTV,"pred_sTV_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_sTV_cpue.initialize();
  #endif
  N_sTV.allocate(styr_cpue_sTV,endyr_cpue_sTV,1,nages,"N_sTV");
  #ifndef NO_AD_INITIALIZE
    N_sTV.initialize();
  #endif
  pred_cHL_cpue.allocate(styr_cpue_cHL,endyr_cpue_cHL,"pred_cHL_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_cHL_cpue.initialize();
  #endif
  N_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL,1,nages,"N_cHL");
  #ifndef NO_AD_INITIALIZE
    N_cHL.initialize();
  #endif
  pred_rHB_cpue.allocate(styr_cpue_rHB,endyr_cpue_rHB,"pred_rHB_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_rHB_cpue.initialize();
  #endif
  N_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,1,nages,"N_rHB");
  #ifndef NO_AD_INITIALIZE
    N_rHB.initialize();
  #endif
  q_rate.allocate(0.001,0.1,set_q_rate_phase,"q_rate");
  log_q_sBT.allocate(log_q_sBT_LO,log_q_sBT_HI,log_q_sBT_PH,"log_q_sBT");
  log_q_sTV.allocate(log_q_sTV_LO,log_q_sTV_HI,log_q_sTV_PH,"log_q_sTV");
  log_q_cHL.allocate(log_q_cHL_LO,log_q_cHL_HI,log_q_cHL_PH,"log_q_cHL");
  log_q_rHB.allocate(log_q_rHB_LO,log_q_rHB_HI,log_q_rHB_PH,"log_q_rHB");
  log_q_sBT_out.allocate(1,8,"log_q_sBT_out");
  #ifndef NO_AD_INITIALIZE
    log_q_sBT_out.initialize();
  #endif
  log_q_sTV_out.allocate(1,8,"log_q_sTV_out");
  #ifndef NO_AD_INITIALIZE
    log_q_sTV_out.initialize();
  #endif
  log_q_cHL_out.allocate(1,8,"log_q_cHL_out");
  #ifndef NO_AD_INITIALIZE
    log_q_cHL_out.initialize();
  #endif
  log_q_rHB_out.allocate(1,8,"log_q_rHB_out");
  #ifndef NO_AD_INITIALIZE
    log_q_rHB_out.initialize();
  #endif
  q_rate_fcn_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL,"q_rate_fcn_cHL");
  #ifndef NO_AD_INITIALIZE
    q_rate_fcn_cHL.initialize();
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
  q_RW_log_dev_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL-1,log_RWq_LO,log_RWq_HI,log_RWq_PH,"q_RW_log_dev_cHL");
  q_RW_log_dev_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB-1,log_RWq_LO,log_RWq_HI,log_RWq_PH,"q_RW_log_dev_rHB");
  q_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL,"q_cHL");
  #ifndef NO_AD_INITIALIZE
    q_cHL.initialize();
  #endif
  q_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"q_rHB");
  #ifndef NO_AD_INITIALIZE
    q_rHB.initialize();
  #endif
  L_rHB_bias.allocate("L_rHB_bias");
  #ifndef NO_AD_INITIALIZE
  L_rHB_bias.initialize();
  #endif
  L_rGN_bias.allocate("L_rGN_bias");
  #ifndef NO_AD_INITIALIZE
  L_rGN_bias.initialize();
  #endif
  L_cGN_bias.allocate("L_cGN_bias");
  #ifndef NO_AD_INITIALIZE
  L_cGN_bias.initialize();
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
  L_cPT_num.allocate(styr,endyr,1,nages,"L_cPT_num");
  #ifndef NO_AD_INITIALIZE
    L_cPT_num.initialize();
  #endif
  L_cPT_klb.allocate(styr,endyr,1,nages,"L_cPT_klb");
  #ifndef NO_AD_INITIALIZE
    L_cPT_klb.initialize();
  #endif
  pred_cPT_L_knum.allocate(styr,endyr,"pred_cPT_L_knum");
  #ifndef NO_AD_INITIALIZE
    pred_cPT_L_knum.initialize();
  #endif
  pred_cPT_L_klb.allocate(styr,endyr,"pred_cPT_L_klb");
  #ifndef NO_AD_INITIALIZE
    pred_cPT_L_klb.initialize();
  #endif
  L_cTW_num.allocate(styr,endyr,1,nages,"L_cTW_num");
  #ifndef NO_AD_INITIALIZE
    L_cTW_num.initialize();
  #endif
  L_cTW_klb.allocate(styr,endyr,1,nages,"L_cTW_klb");
  #ifndef NO_AD_INITIALIZE
    L_cTW_klb.initialize();
  #endif
  pred_cTW_L_knum.allocate(styr,endyr,"pred_cTW_L_knum");
  #ifndef NO_AD_INITIALIZE
    pred_cTW_L_knum.initialize();
  #endif
  pred_cTW_L_klb.allocate(styr,endyr,"pred_cTW_L_klb");
  #ifndef NO_AD_INITIALIZE
    pred_cTW_L_klb.initialize();
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
  obs_rHB_L_wbias.allocate(styr,endyr,"obs_rHB_L_wbias");
  #ifndef NO_AD_INITIALIZE
    obs_rHB_L_wbias.initialize();
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
  obs_rGN_L_wbias.allocate(styr,endyr,"obs_rGN_L_wbias");
  #ifndef NO_AD_INITIALIZE
    obs_rGN_L_wbias.initialize();
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
  D_cGN_num.allocate(styr,endyr,1,nages,"D_cGN_num");
  #ifndef NO_AD_INITIALIZE
    D_cGN_num.initialize();
  #endif
  D_cGN_klb.allocate(styr,endyr,1,nages,"D_cGN_klb");
  #ifndef NO_AD_INITIALIZE
    D_cGN_klb.initialize();
  #endif
  pred_cGN_D_knum.allocate(styr,endyr,"pred_cGN_D_knum");
  #ifndef NO_AD_INITIALIZE
    pred_cGN_D_knum.initialize();
  #endif
  pred_cGN_D_klb.allocate(styr,endyr,"pred_cGN_D_klb");
  #ifndef NO_AD_INITIALIZE
    pred_cGN_D_klb.initialize();
  #endif
  obs_cGN_D.allocate(styr_D_cHL,endyr_D_cHL,"obs_cGN_D");
  #ifndef NO_AD_INITIALIZE
    obs_cGN_D.initialize();
  #endif
  obs_cHL_D.allocate(styr_D_cHL,endyr_D_cHL,"obs_cHL_D");
  #ifndef NO_AD_INITIALIZE
    obs_cHL_D.initialize();
  #endif
  obs_cPT_D.allocate(styr_D_cHL,endyr_D_cHL,"obs_cPT_D");
  #ifndef NO_AD_INITIALIZE
    obs_cPT_D.initialize();
  #endif
  cGN_D_cv.allocate(styr_D_cHL,endyr_D_cHL,"cGN_D_cv");
  #ifndef NO_AD_INITIALIZE
    cGN_D_cv.initialize();
  #endif
  D_rHB_num.allocate(styr,endyr,1,nages,"D_rHB_num");
  #ifndef NO_AD_INITIALIZE
    D_rHB_num.initialize();
  #endif
  D_rHB_klb.allocate(styr,endyr,1,nages,"D_rHB_klb");
  #ifndef NO_AD_INITIALIZE
    D_rHB_klb.initialize();
  #endif
  pred_rHB_D_knum.allocate(styr,endyr,"pred_rHB_D_knum");
  #ifndef NO_AD_INITIALIZE
    pred_rHB_D_knum.initialize();
  #endif
  pred_rHB_D_klb.allocate(styr,endyr,"pred_rHB_D_klb");
  #ifndef NO_AD_INITIALIZE
    pred_rHB_D_klb.initialize();
  #endif
  obs_rHB_D.allocate(styr_D_rHB,endyr_D_rHB,"obs_rHB_D");
  #ifndef NO_AD_INITIALIZE
    obs_rHB_D.initialize();
  #endif
  D_rGN_num.allocate(styr,endyr,1,nages,"D_rGN_num");
  #ifndef NO_AD_INITIALIZE
    D_rGN_num.initialize();
  #endif
  D_rGN_klb.allocate(styr,endyr,1,nages,"D_rGN_klb");
  #ifndef NO_AD_INITIALIZE
    D_rGN_klb.initialize();
  #endif
  pred_rGN_D_knum.allocate(styr,endyr,"pred_rGN_D_knum");
  #ifndef NO_AD_INITIALIZE
    pred_rGN_D_knum.initialize();
  #endif
  pred_rGN_D_klb.allocate(styr,endyr,"pred_rGN_D_klb");
  #ifndef NO_AD_INITIALIZE
    pred_rGN_D_klb.initialize();
  #endif
  obs_rGN_D.allocate(styr_D_rGN,endyr_D_rGN,"obs_rGN_D");
  #ifndef NO_AD_INITIALIZE
    obs_rGN_D.initialize();
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
  F_cHL_prop.allocate("F_cHL_prop");
  #ifndef NO_AD_INITIALIZE
  F_cHL_prop.initialize();
  #endif
  F_cPT_prop.allocate("F_cPT_prop");
  #ifndef NO_AD_INITIALIZE
  F_cPT_prop.initialize();
  #endif
  F_rHB_prop.allocate("F_rHB_prop");
  #ifndef NO_AD_INITIALIZE
  F_rHB_prop.initialize();
  #endif
  F_rGN_prop.allocate("F_rGN_prop");
  #ifndef NO_AD_INITIALIZE
  F_rGN_prop.initialize();
  #endif
  F_cGN_D_prop.allocate("F_cGN_D_prop");
  #ifndef NO_AD_INITIALIZE
  F_cGN_D_prop.initialize();
  #endif
  F_rHB_D_prop.allocate("F_rHB_D_prop");
  #ifndef NO_AD_INITIALIZE
  F_rHB_D_prop.initialize();
  #endif
  F_rGN_D_prop.allocate("F_rGN_D_prop");
  #ifndef NO_AD_INITIALIZE
  F_rGN_D_prop.initialize();
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
  B_msy_out.allocate("B_msy_out");
  #ifndef NO_AD_INITIALIZE
  B_msy_out.initialize();
  #endif
  R_msy_out.allocate("R_msy_out");
  #ifndef NO_AD_INITIALIZE
  R_msy_out.initialize();
  #endif
  D_msy_knum_out.allocate("D_msy_knum_out");
  #ifndef NO_AD_INITIALIZE
  D_msy_knum_out.initialize();
  #endif
  D_msy_klb_out.allocate("D_msy_klb_out");
  #ifndef NO_AD_INITIALIZE
  D_msy_klb_out.initialize();
  #endif
  spr_msy_out.allocate("spr_msy_out");
  #ifndef NO_AD_INITIALIZE
  spr_msy_out.initialize();
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
  D_age_msy.allocate(1,nages,"D_age_msy");
  #ifndef NO_AD_INITIALIZE
    D_age_msy.initialize();
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
  SSB_eq.allocate(1,n_iter_msy,"SSB_eq");
  #ifndef NO_AD_INITIALIZE
    SSB_eq.initialize();
  #endif
  B_eq.allocate(1,n_iter_msy,"B_eq");
  #ifndef NO_AD_INITIALIZE
    B_eq.initialize();
  #endif
  D_eq_klb.allocate(1,n_iter_msy,"D_eq_klb");
  #ifndef NO_AD_INITIALIZE
    D_eq_klb.initialize();
  #endif
  D_eq_knum.allocate(1,n_iter_msy,"D_eq_knum");
  #ifndef NO_AD_INITIALIZE
    D_eq_knum.initialize();
  #endif
  FdF_msy.allocate(styr,endyr,"FdF_msy");
  #ifndef NO_AD_INITIALIZE
    FdF_msy.initialize();
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
  log_avg_F_L_cPT.allocate(log_avg_F_cPT_LO,log_avg_F_cPT_HI,log_avg_F_cPT_PH,"log_avg_F_L_cPT");
  log_avg_F_cPT_out.allocate(1,8,"log_avg_F_cPT_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_cPT_out.initialize();
  #endif
  log_dev_F_L_cPT.allocate(styr_L_cPT,endyr_L_cPT,log_F_dev_cPT_LO,log_F_dev_cPT_HI,log_F_dev_cPT_PH,"log_dev_F_L_cPT");
  log_F_dev_cPT_out.allocate(styr_L_cPT,endyr_L_cPT,"log_F_dev_cPT_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_cPT_out.initialize();
  #endif
  F_cPT.allocate(styr,endyr,1,nages,"F_cPT");
  #ifndef NO_AD_INITIALIZE
    F_cPT.initialize();
  #endif
  F_cPT_out.allocate(styr,endyr,"F_cPT_out");
  #ifndef NO_AD_INITIALIZE
    F_cPT_out.initialize();
  #endif
  log_F_dev_init_cPT.allocate("log_F_dev_init_cPT");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_init_cPT.initialize();
  #endif
  log_F_dev_end_cPT.allocate("log_F_dev_end_cPT");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_cPT.initialize();
  #endif
  log_avg_F_L_cTW.allocate(log_avg_F_cTW_LO,log_avg_F_cTW_HI,log_avg_F_cTW_PH,"log_avg_F_L_cTW");
  log_avg_F_cTW_out.allocate(1,8,"log_avg_F_cTW_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_cTW_out.initialize();
  #endif
  log_dev_F_L_cTW.allocate(styr_L_cTW,endyr_L_cTW,log_F_dev_cTW_LO,log_F_dev_cTW_HI,log_F_dev_cTW_PH,"log_dev_F_L_cTW");
  log_F_dev_cTW_out.allocate(styr_L_cTW,endyr_L_cTW,"log_F_dev_cTW_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_cTW_out.initialize();
  #endif
  F_cTW.allocate(styr,endyr,1,nages,"F_cTW");
  #ifndef NO_AD_INITIALIZE
    F_cTW.initialize();
  #endif
  F_cTW_out.allocate(styr,endyr,"F_cTW_out");
  #ifndef NO_AD_INITIALIZE
    F_cTW_out.initialize();
  #endif
  log_F_dev_init_cTW.allocate("log_F_dev_init_cTW");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_init_cTW.initialize();
  #endif
  log_F_dev_end_cTW.allocate("log_F_dev_end_cTW");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_cTW.initialize();
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
  log_F_init_rHB.allocate("log_F_init_rHB");
  #ifndef NO_AD_INITIALIZE
  log_F_init_rHB.initialize();
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
  F_init_ratio.allocate(0.1,1.5,1,"F_init_ratio");
  log_avg_F_D_cGN.allocate(log_avg_F_cGN_D_LO,log_avg_F_cGN_D_HI,log_avg_F_cGN_D_PH,"log_avg_F_D_cGN");
  log_avg_F_cGN_D_out.allocate(1,8,"log_avg_F_cGN_D_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_cGN_D_out.initialize();
  #endif
  log_dev_F_D_cGN.allocate(styr_D_cHL,endyr_D_cHL,log_F_dev_cGN_D_LO,log_F_dev_cGN_D_HI,log_F_dev_cGN_D_PH,"log_dev_F_D_cGN");
  log_F_dev_cGN_D_out.allocate(styr_D_cHL,endyr_D_cHL,"log_F_dev_cGN_D_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_cGN_D_out.initialize();
  #endif
  F_cGN_D.allocate(styr,endyr,1,nages,"F_cGN_D");
  #ifndef NO_AD_INITIALIZE
    F_cGN_D.initialize();
  #endif
  F_cGN_D_out.allocate(styr,endyr,"F_cGN_D_out");
  #ifndef NO_AD_INITIALIZE
    F_cGN_D_out.initialize();
  #endif
  log_F_dev_cGN_D2.allocate("log_F_dev_cGN_D2");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_cGN_D2.initialize();
  #endif
  log_F_dev_end_cGN_D.allocate("log_F_dev_end_cGN_D");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_cGN_D.initialize();
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
  Dmort_HL.allocate("Dmort_HL");
  #ifndef NO_AD_INITIALIZE
  Dmort_HL.initialize();
  #endif
  Dmort_rHB_HL.allocate("Dmort_rHB_HL");
  #ifndef NO_AD_INITIALIZE
  Dmort_rHB_HL.initialize();
  #endif
  Dmort_rGN_HL.allocate("Dmort_rGN_HL");
  #ifndef NO_AD_INITIALIZE
  Dmort_rGN_HL.initialize();
  #endif
  Dmort_cPT1.allocate("Dmort_cPT1");
  #ifndef NO_AD_INITIALIZE
  Dmort_cPT1.initialize();
  #endif
  Dmort_cPT2.allocate("Dmort_cPT2");
  #ifndef NO_AD_INITIALIZE
  Dmort_cPT2.initialize();
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
  sdnr_lc_sBT.allocate("sdnr_lc_sBT");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_sBT.initialize();
  #endif
  sdnr_lc_sTV.allocate("sdnr_lc_sTV");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_sTV.initialize();
  #endif
  sdnr_lc_cHL.allocate("sdnr_lc_cHL");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_cHL.initialize();
  #endif
  sdnr_lc_cPT.allocate("sdnr_lc_cPT");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_cPT.initialize();
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
  sdnr_ac_sBT.allocate("sdnr_ac_sBT");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_sBT.initialize();
  #endif
  sdnr_ac_sTV.allocate("sdnr_ac_sTV");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_sTV.initialize();
  #endif
  sdnr_ac_cHL.allocate("sdnr_ac_cHL");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_cHL.initialize();
  #endif
  sdnr_ac_cPT.allocate("sdnr_ac_cPT");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_cPT.initialize();
  #endif
  sdnr_ac_rHB.allocate("sdnr_ac_rHB");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_rHB.initialize();
  #endif
  sdnr_I_sBT.allocate("sdnr_I_sBT");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_sBT.initialize();
  #endif
  sdnr_I_sTV.allocate("sdnr_I_sTV");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_sTV.initialize();
  #endif
  sdnr_I_cHL.allocate("sdnr_I_cHL");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_cHL.initialize();
  #endif
  sdnr_I_rHB.allocate("sdnr_I_rHB");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_rHB.initialize();
  #endif
  w_L.allocate("w_L");
  #ifndef NO_AD_INITIALIZE
  w_L.initialize();
  #endif
  w_D.allocate("w_D");
  #ifndef NO_AD_INITIALIZE
  w_D.initialize();
  #endif
  w_lenc_sBT.allocate("w_lenc_sBT");
  #ifndef NO_AD_INITIALIZE
  w_lenc_sBT.initialize();
  #endif
  w_lenc_cHL.allocate("w_lenc_cHL");
  #ifndef NO_AD_INITIALIZE
  w_lenc_cHL.initialize();
  #endif
  w_lenc_cPT.allocate("w_lenc_cPT");
  #ifndef NO_AD_INITIALIZE
  w_lenc_cPT.initialize();
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
  w_agec_sBT.allocate("w_agec_sBT");
  #ifndef NO_AD_INITIALIZE
  w_agec_sBT.initialize();
  #endif
  w_agec_sTV.allocate("w_agec_sTV");
  #ifndef NO_AD_INITIALIZE
  w_agec_sTV.initialize();
  #endif
  w_agec_cHL.allocate("w_agec_cHL");
  #ifndef NO_AD_INITIALIZE
  w_agec_cHL.initialize();
  #endif
  w_agec_cPT.allocate("w_agec_cPT");
  #ifndef NO_AD_INITIALIZE
  w_agec_cPT.initialize();
  #endif
  w_agec_rHB.allocate("w_agec_rHB");
  #ifndef NO_AD_INITIALIZE
  w_agec_rHB.initialize();
  #endif
  w_cpue_sBT.allocate("w_cpue_sBT");
  #ifndef NO_AD_INITIALIZE
  w_cpue_sBT.initialize();
  #endif
  w_cpue_sTV.allocate("w_cpue_sTV");
  #ifndef NO_AD_INITIALIZE
  w_cpue_sTV.initialize();
  #endif
  w_cpue_cHL.allocate("w_cpue_cHL");
  #ifndef NO_AD_INITIALIZE
  w_cpue_cHL.initialize();
  #endif
  w_cpue_rHB.allocate("w_cpue_rHB");
  #ifndef NO_AD_INITIALIZE
  w_cpue_rHB.initialize();
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
  f_sBT_cpue.allocate("f_sBT_cpue");
  #ifndef NO_AD_INITIALIZE
  f_sBT_cpue.initialize();
  #endif
  f_sTV_cpue.allocate("f_sTV_cpue");
  #ifndef NO_AD_INITIALIZE
  f_sTV_cpue.initialize();
  #endif
  f_Vid_cpue.allocate("f_Vid_cpue");
  #ifndef NO_AD_INITIALIZE
  f_Vid_cpue.initialize();
  #endif
  f_cHL_cpue.allocate("f_cHL_cpue");
  #ifndef NO_AD_INITIALIZE
  f_cHL_cpue.initialize();
  #endif
  f_rHB_cpue.allocate("f_rHB_cpue");
  #ifndef NO_AD_INITIALIZE
  f_rHB_cpue.initialize();
  #endif
  f_cHL_L.allocate("f_cHL_L");
  #ifndef NO_AD_INITIALIZE
  f_cHL_L.initialize();
  #endif
  f_cPT_L.allocate("f_cPT_L");
  #ifndef NO_AD_INITIALIZE
  f_cPT_L.initialize();
  #endif
  f_cTW_L.allocate("f_cTW_L");
  #ifndef NO_AD_INITIALIZE
  f_cTW_L.initialize();
  #endif
  f_rHB_L.allocate("f_rHB_L");
  #ifndef NO_AD_INITIALIZE
  f_rHB_L.initialize();
  #endif
  f_rGN_L.allocate("f_rGN_L");
  #ifndef NO_AD_INITIALIZE
  f_rGN_L.initialize();
  #endif
  f_cGN_D.allocate("f_cGN_D");
  #ifndef NO_AD_INITIALIZE
  f_cGN_D.initialize();
  #endif
  f_rHB_D.allocate("f_rHB_D");
  #ifndef NO_AD_INITIALIZE
  f_rHB_D.initialize();
  #endif
  f_rGN_D.allocate("f_rGN_D");
  #ifndef NO_AD_INITIALIZE
  f_rGN_D.initialize();
  #endif
  f_sBT_lenc.allocate("f_sBT_lenc");
  #ifndef NO_AD_INITIALIZE
  f_sBT_lenc.initialize();
  #endif
  f_cHL_lenc.allocate("f_cHL_lenc");
  #ifndef NO_AD_INITIALIZE
  f_cHL_lenc.initialize();
  #endif
  f_cPT_lenc.allocate("f_cPT_lenc");
  #ifndef NO_AD_INITIALIZE
  f_cPT_lenc.initialize();
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
  f_sBT_agec.allocate("f_sBT_agec");
  #ifndef NO_AD_INITIALIZE
  f_sBT_agec.initialize();
  #endif
  f_sTV_agec.allocate("f_sTV_agec");
  #ifndef NO_AD_INITIALIZE
  f_sTV_agec.initialize();
  #endif
  f_cHL_agec.allocate("f_cHL_agec");
  #ifndef NO_AD_INITIALIZE
  f_cHL_agec.initialize();
  #endif
  f_cPT_agec.allocate("f_cPT_agec");
  #ifndef NO_AD_INITIALIZE
  f_cPT_agec.initialize();
  #endif
  f_rHB_agec.allocate("f_rHB_agec");
  #ifndef NO_AD_INITIALIZE
  f_rHB_agec.initialize();
  #endif
  f_rGN_agec.allocate("f_rGN_agec");
  #ifndef NO_AD_INITIALIZE
  f_rGN_agec.initialize();
  #endif
  f_cHL_RW_cpue.allocate("f_cHL_RW_cpue");
  #ifndef NO_AD_INITIALIZE
  f_cHL_RW_cpue.initialize();
  #endif
  f_rHB_RW_cpue.allocate("f_rHB_RW_cpue");
  #ifndef NO_AD_INITIALIZE
  f_rHB_RW_cpue.initialize();
  #endif
  f_rHB_D_RW_cpue.allocate("f_rHB_D_RW_cpue");
  #ifndef NO_AD_INITIALIZE
  f_rHB_D_RW_cpue.initialize();
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
  f_Ftune.allocate("f_Ftune");
  #ifndef NO_AD_INITIALIZE
  f_Ftune.initialize();
  #endif
  f_fullF_constraint.allocate("f_fullF_constraint");
  #ifndef NO_AD_INITIALIZE
  f_fullF_constraint.initialize();
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

void model_parameters::initializationfunction(void)
{
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{10000, 20000, 30000, 50000, 100000, 100000, 100000;}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-2, 1e-2,1e-2, 1e-4;}");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
 sqrt2pi=sqrt(2.*3.14159265);
 g2mt=0.000001;         //conversion of grams to metric tons
 g2kg=0.001;            //conversion of grams to kg 
 mt2klb=2.20462;        //conversion of metric tons to 1000 lb 
 mt2lb=mt2klb*1000.0;   //conversion of metric tons to lb
 g2klb=g2mt*mt2klb;     //conversion of grams to 1000 lb 
 dzero=0.00001;         
 huge_number=1.0e+10;   
 onehalf=0.5;
  Dmort_HL=set_Dmort_HL;
  Dmort_rHB_HL=set_Dmort_rHB_HL;
  Dmort_rGN_HL=set_Dmort_rGN_HL;
  Dmort_cPT1=set_Dmort_cPT1;
  Dmort_cPT2=set_Dmort_cPT2;
  //values used for weighting selex and avg weight of cGN discards in yrs with quotas
  //geometric mean of last three yrs
  //avg weights of landings were near 1 lb, so those values are left in weight
  Dopen_cHL=Dmort_HL*pow((obs_released_cHL(2014)*obs_released_cHL(2015)*obs_released_cHL(endyr_D_cHL)),(1.0/3.0));
  Dclosed_cHL=Dmort_HL*((obs_released_cHL_closed(2013)+obs_released_cHL_closed(2014)+obs_released_cHL_closed(2015)+obs_released_cHL_closed(endyr_D_cHL))/4.0);
  //Dclosed_cHL=Dmort_HL*pow((obs_released_cHL_closed(2014)*obs_released_cHL_closed(2015)*obs_released_cHL_closed(endyr_D_cHL)),(1.0/3.0));
  Lopen_cHL=Dmort_HL*pow((obs_L_cHL(2014)*obs_L_cHL(2015)*obs_L_cHL(endyr_L_cHL)),(1.0/3.0));
  
  Dopen_cPT=Dmort_cPT2*pow((obs_released_cPT(2014)*obs_released_cPT(2015)*obs_released_cPT(endyr_D_cPT)),(1.0/3.0));
  Dclosed_cPT=Dmort_cPT2*(obs_released_cPT_closed(2013)+(obs_released_cPT_closed(2014)+obs_released_cPT_closed(2015)+obs_released_cPT_closed(endyr_D_cPT))/4.0);
  //Dclosed_cPT=Dmort_cPT2*pow((obs_released_cPT_closed(2014)*obs_released_cPT_closed(2015)*obs_released_cPT_closed(endyr_D_cPT)),(1.0/3.0));
  Lopen_cPT=Dmort_cPT2*pow((obs_L_cPT(2014)*obs_L_cPT(2015)*obs_L_cPT(endyr_L_cPT)),(1.0/3.0));
  
  D_sum_cHLcPT=Dopen_cHL+Dclosed_cHL+Dopen_cPT+Dclosed_cPT;
  Dprop_cGN_sel_D=(Dopen_cHL + Dopen_cPT + Dclosed_cHL*(Dopen_cHL/(Dopen_cHL+Lopen_cHL)) +
                    Dclosed_cPT*(Dopen_cPT/(Dopen_cPT+Lopen_cPT)))/D_sum_cHLcPT; 
  Dprop_cGN_sel_cHL=Dclosed_cHL*(Lopen_cHL/(Dopen_cHL+Lopen_cHL))/D_sum_cHLcPT; 
  Dprop_cGN_sel_cPT=Dclosed_cPT*(Lopen_cPT/(Dopen_cPT+Lopen_cPT))/D_sum_cHLcPT; 
  //discards values for fitting, include discard mortality
  
  obs_cHL_D=Dmort_HL*obs_released_cHL;
  obs_cHL_D(styr_D_cHL_closed,endyr_D_cHL_closed)+=Dmort_HL*obs_released_cHL_closed;
  
  obs_cPT_D(styr_D_cPT,2006)=Dmort_cPT1*obs_released_cPT(styr_D_cPT,2006);
  obs_cPT_D(2007,endyr_D_cPT)=Dmort_cPT2*obs_released_cPT(2007,endyr_D_cPT);
  obs_cPT_D(styr_D_cPT_closed,endyr_D_cPT_closed)+=Dmort_cPT2*obs_released_cPT_closed;
  
  obs_cGN_D=obs_cHL_D+obs_cPT_D;
  
  obs_rHB_D=Dmort_rHB_HL*obs_released_rHB;
  obs_rGN_D=Dmort_rGN_HL*obs_released_rGN;
  
  cGN_D_cv=obs_cv_D_cHL;
 
  Linf=set_Linf(1);
  K=set_K(1);
  t0=set_t0(1);
  M=set_M; 
  M_constant=set_M_constant(1);
  
  
  steep=set_steep(1);
  R_autocorr=set_R_autocorr(1);
  rec_sigma=set_rec_sigma(1);
  
  log_dm_lenc_sBT=set_log_dm_lenc_sBT(1);
  log_dm_lenc_sTV=set_log_dm_lenc_sTV(1);
  log_dm_lenc_cHL=set_log_dm_lenc_cHL(1);
  log_dm_lenc_cPT=set_log_dm_lenc_cPT(1);
  log_dm_lenc_rHB=set_log_dm_lenc_rHB(1);
  log_dm_lenc_rHB_D=set_log_dm_lenc_rHB_D(1);
  log_dm_lenc_rGN=set_log_dm_lenc_rGN(1);
  
  log_dm_agec_sBT=set_log_dm_agec_sBT(1);
  log_dm_agec_sTV=set_log_dm_agec_sTV(1);
  log_dm_agec_cHL=set_log_dm_agec_cHL(1);
  log_dm_agec_cPT=set_log_dm_agec_cPT(1);
  log_dm_agec_rHB=set_log_dm_agec_rHB(1);
  //log_dm_rGN_ac=set_log_dm_rGN_ac(1);
    
  log_q_sBT=set_log_q_cpue_sBT(1);  
  log_q_sTV=set_log_q_cpue_sTV(1);
  log_q_cHL=set_log_q_cpue_cHL(1);
  log_q_rHB=set_log_q_cpue_rHB(1);
  //log_q_rHB_D=set_logq_rHB_D(1);
  q_rate=set_q_rate;
  q_rate_fcn_cHL=1.0;
  q_rate_fcn_rHB=1.0;
  //q_rate_fcn_rHB_D=1.0;
  q_DD_beta=set_q_DD_beta;
  q_DD_fcn=1.0;
  q_RW_log_dev_cHL.initialize();
  q_RW_log_dev_rHB.initialize();
  //q_RW_log_dev_rHB_D.initialize();
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
      //for (iyear=styr_rHB_D_cpue; iyear<=endyr_rHB_D_cpue; iyear++)
      //{   if (iyear>styr_rHB_D_cpue & iyear <=2003) 
       //   {//q_rate_fcn_rHB_D(iyear)=(1.0+q_rate)*q_rate_fcn_rHB_D(iyear-1); //compound
       //      q_rate_fcn_rHB_D(iyear)=(1.0+(iyear-styr_rHB_D_cpue)*q_rate)*q_rate_fcn_rHB_D(styr_rHB_D_cpue);  //linear
       //   }
       //   if (iyear>2003) {q_rate_fcn_rHB_D(iyear)=q_rate_fcn_rHB_D(iyear-1);} 
      //}
  } //end q_rate conditional      
  L_rHB_bias=set_L_rHB_bias;
  L_rGN_bias=set_L_rGN_bias;
  L_cGN_bias=set_L_cGN_bias;
  w_L=set_w_L;
  w_D=set_w_D;
  
  w_lenc_sBT=set_w_lenc_sBT;
  w_lenc_cHL=set_w_lenc_cHL;
  w_lenc_cPT=set_w_lenc_cPT;
  w_lenc_rHB=set_w_lenc_rHB;
  w_lenc_rHB_D=set_w_lenc_rHB_D;
  w_lenc_rGN=set_w_lenc_rGN;
  
  w_agec_sBT=set_w_agec_sBT;
  w_agec_sTV=set_w_agec_sTV;
  w_agec_cHL=set_w_agec_cHL;
  w_agec_cPT=set_w_agec_cPT;
  w_agec_rHB=set_w_agec_rHB;
  //w_ac_rGN=set_w_ac_rGN;  
  w_cpue_sTV=set_w_cpue_sTV;
  w_cpue_sBT=set_w_cpue_sBT;
  //w_I_Vid=set_w_I_Vid;
  w_cpue_cHL=set_w_cpue_cHL;
  w_cpue_rHB=set_w_cpue_rHB;
  //w_I_rHB_D=set_w_I_rHB_D;
  w_rec=set_w_rec;
  w_fullF=set_w_fullF;
  w_rec_early=set_w_rec_early;
  w_rec_end=set_w_rec_end;
  w_Ftune=set_w_Ftune;
  //w_cvlen_dev=set_w_cvlen_dev;
  //w_cvlen_diff=set_w_cvlen_diff;
  log_avg_F_L_cHL=set_log_avg_F_L_cHL(1);  
  log_avg_F_L_cPT=set_log_avg_F_L_cPT(1);
  log_avg_F_L_cTW=set_log_avg_F_L_cTW(1);  
  log_avg_F_L_rHB=set_log_avg_F_L_rHB(1);
  log_avg_F_L_rGN=set_log_avg_F_L_rGN(1);  
  F_init_ratio=set_F_init_ratio;
  
  log_avg_F_D_cGN=set_log_avg_F_D_cGN(1); 
  log_avg_F_D_rHB=set_log_avg_F_D_rHB(1);
  log_avg_F_D_rGN=set_log_avg_F_D_rGN(1);  
 
  len_cv_val=set_len_cv(1);
  log_R0=set_log_R0(1);
  selpar_A50_sBT= set_selpar_A50_sBT(1);
  selpar_slope_sBT=set_selpar_slope_sBT(1); 
  selpar_A50_sTV=set_selpar_A50_sTV(1);
  selpar_slope_sTV=set_selpar_slope_sTV(1); 
  selpar_A50_cHL2=set_selpar_A50_cHL2(1);
  selpar_slope_cHL2=set_selpar_slope_cHL2(1); 
  selpar_A50_cHL3=set_selpar_A50_cHL3(1);
  selpar_slope_cHL3=set_selpar_slope_cHL3(1); 
  selpar_A50_cHL4=set_selpar_A50_cHL4(1);
  selpar_slope_cHL4=set_selpar_slope_cHL4(1); 
  selpar_A50_cPT2=set_selpar_A50_cPT2(1);
  selpar_slope_cPT2=set_selpar_slope_cPT2(1); 
  selpar_A50_cPT3=set_selpar_A50_cPT3(1);
  selpar_slope_cPT3=set_selpar_slope_cPT3(1); 
  selpar_A50_cPT4=set_selpar_A50_cPT4(1);
  selpar_slope_cPT4=set_selpar_slope_cPT4(1);
  selpar_A50_rHB1=set_selpar_A50_rHB1(1);
  selpar_slope_rHB1=set_selpar_slope_rHB1(1); 
  selpar_A50_rHB2=set_selpar_A50_rHB2(1);
  selpar_slope_rHB2=set_selpar_slope_rHB2(1); 
  selpar_A50_rHB3=set_selpar_A50_rHB3(1);
  selpar_slope_rHB3=set_selpar_slope_rHB3(1); 
  selpar_A50_rHB4=set_selpar_A50_rHB4(1);
  selpar_slope_rHB4=set_selpar_slope_rHB4(1); 
  selpar_A50_rHB5=set_selpar_A50_rHB5(1);
  selpar_slope_rHB5=set_selpar_slope_rHB5(1); 
  selpar_A50_rGN1=set_selpar_A50_rGN1(1);
  selpar_slope_rGN1=set_selpar_slope_rGN1(1); 
  selpar_A50_rGN2=set_selpar_A50_rGN2(1);
  selpar_slope_rGN2=set_selpar_slope_rGN2(1); 
  selpar_A50_rGN3=set_selpar_A50_rGN3(1);
  selpar_slope_rGN3=set_selpar_slope_rGN3(1); 
  selpar_A50_rGN4=set_selpar_A50_rGN4(1);
  selpar_slope_rGN4=set_selpar_slope_rGN4(1); 
  selpar_A50_rGN5=set_selpar_A50_rGN5(1);
  selpar_slope_rGN5=set_selpar_slope_rGN5(1); 
  selpar_logit_Age0_rHB_D=set_selpar_logit_Age0_rHB_D(1);  
  selpar_logit_Age1_rHB_D=set_selpar_logit_Age1_rHB_D(1); 
  selpar_logit_Age2_rHB_D=set_selpar_logit_Age2_rHB_D(1);
  selpar_A50_rHB_D4=set_selpar_A50_rHB_D4(1);
  selpar_slope_rHB_D4=set_selpar_slope_rHB_D4(1); 
  selpar_A502_rHB_D4=set_selpar_A502_rHB_D4(1);
  selpar_slope2_rHB_D4=set_selpar_slope2_rHB_D4(1); 
  selpar_A50_rHB_D5=set_selpar_A50_rHB_D5(1);
  selpar_slope_rHB_D5=set_selpar_slope_rHB_D5(1);
  selpar_A502_rHB_D5=set_selpar_A502_rHB_D5(1);
  selpar_slope2_rHB_D5=set_selpar_slope2_rHB_D5(1);
   
 SSB_msy_out=0.0;
 iter_inc_msy=max_F_spr_msy/(n_iter_msy-1);
 iter_inc_spr=max_F_spr_msy/(n_iter_spr-1); 
 maturity_f=obs_maturity_f;
 prop_f=obs_prop_f;
 p_lenc_cHL2=set_p_lenc_cHL2;
 p_lenc_cHL3=set_p_lenc_cHL3; 
 p_lenc_cPT2=set_p_lenc_cPT2;
 p_lenc_cPT3=set_p_lenc_cPT3; 
 p_lenc_cTW2=set_p_lenc_cTW2;
 p_lenc_cTW3=set_p_lenc_cTW3;  
 p_lenc_rHB2=set_p_lenc_rHB2;
 p_lenc_rHB3=set_p_lenc_rHB3;
 p_lenc_rHB4=set_p_lenc_rHB4; 
 p_lenc_rHB5=set_p_lenc_rHB5; 
 p_lenc_rGN2=set_p_lenc_rGN2;
 p_lenc_rGN3=set_p_lenc_rGN3;
 p_lenc_rGN4=set_p_lenc_rGN4;
 p_lenc_rGN5=set_p_lenc_rGN5;
  
 p_lenc_cGN_D2=set_p_lenc_cGN_D2; 
 p_lenc_cGN_D3=set_p_lenc_cGN_D3;
 p_lenc_cGN_D4=set_p_lenc_cGN_D4;
 p_lenc_rHB_D2=set_p_lenc_rHB_D2;
 p_lenc_rHB_D3=set_p_lenc_rHB_D3; 
 p_lenc_rHB_D4=set_p_lenc_rHB_D4;
 p_lenc_rHB_D5=set_p_lenc_rHB_D5; 
 p_lenc_rGN_D1=set_p_lenc_rGN_D1; 
 p_lenc_rGN_D2=set_p_lenc_rGN_D2;
 p_lenc_rGN_D3=set_p_lenc_rGN_D3; 
 p_lenc_rGN_D4=set_p_lenc_rGN_D4; 
 p_lenc_rGN_D5=set_p_lenc_rGN_D5;
 
 lenbins_all(1,nlenbins)=lenbins(1,nlenbins);  
 for (iyear=1;iyear<=nlenbins_plus; iyear++) {lenbins_all(nlenbins+iyear)=lenbins_plus(iyear);} 
 //multiplicative bias for early rec data
 obs_rHB_L_wbias(styr_L_rHB,endyr_L_rHB_bias)=L_rHB_bias*obs_L_rHB(styr_L_rHB,endyr_L_rHB_bias); 
 obs_rHB_L_wbias((endyr_L_rHB_bias+1),endyr_L_rHB)=obs_L_rHB((endyr_L_rHB_bias+1),endyr_L_rHB); 
 obs_rGN_L_wbias(styr_L_rGN,endyr_L_rGN_bias)=L_rHB_bias*obs_L_rGN(styr_L_rGN,endyr_L_rGN_bias); 
 obs_rGN_L_wbias((endyr_L_rGN_bias+1),endyr_L_rGN)=obs_L_rGN((endyr_L_rGN_bias+1),endyr_L_rGN); 
 
  nsamp_sBT_lenc_allyr=missing;//"missing" defined in admb2r.cpp
  nsamp_cHL_lenc_allyr=missing;  
  nsamp_cPT_lenc_allyr=missing;
  nsamp_rHB_lenc_allyr=missing;
  nsamp_rHB_D_lenc_allyr=missing;
  nsamp_rGN_lenc_allyr=missing;
  nsamp_sBT_agec_allyr=missing;
  nsamp_sTV_agec_allyr=missing;  
  nsamp_cHL_agec_allyr=missing;
  nsamp_cPT_agec_allyr=missing;      
  nsamp_rHB_agec_allyr=missing;
  //nsamp_rGN_agec_allyr=missing;
  nfish_sBT_lenc_allyr=missing;//"missing" defined in admb2r.cpp
  nfish_cHL_lenc_allyr=missing;  
  nfish_cPT_lenc_allyr=missing;
  nfish_rHB_lenc_allyr=missing;
  nfish_rHB_D_lenc_allyr=missing;
  nfish_rGN_lenc_allyr=missing;
  nfish_sBT_agec_allyr=missing;
  nfish_sTV_agec_allyr=missing;  
  nfish_cHL_agec_allyr=missing;
  nfish_cPT_agec_allyr=missing;      
  nfish_rHB_agec_allyr=missing;
  //nfish_rGN_agec_allyr=missing;
          
      for (iyear=1; iyear<=nyr_lenc_sBT; iyear++)
         {  if (nsamp_lenc_sBT(iyear)>maxSS_lenc)
             {nsamp_lenc_sBT(iyear)=maxSS_lenc;}
            if (nsamp_lenc_sBT(iyear)>=minSS_lenc)
             {nsamp_sBT_lenc_allyr(yrs_lenc_sBT(iyear))=nsamp_lenc_sBT(iyear);
              nfish_sBT_lenc_allyr(yrs_lenc_sBT(iyear))=nfish_lenc_sBT(iyear);}}
            
      for (iyear=1; iyear<=nyr_lenc_cHL; iyear++)
         {  if (nsamp_lenc_cHL(iyear)>maxSS_lenc)
             {nsamp_lenc_cHL(iyear)=maxSS_lenc;}
            if (nsamp_lenc_cHL(iyear)>=minSS_lenc)
             {nsamp_cHL_lenc_allyr(yrs_lenc_cHL(iyear))=nsamp_lenc_cHL(iyear);
              nfish_cHL_lenc_allyr(yrs_lenc_cHL(iyear))=nfish_lenc_cHL(iyear);}}
         
      for (iyear=1; iyear<=nyr_lenc_cPT; iyear++)
         {  if (nsamp_lenc_cPT(iyear)>maxSS_lenc)
             {nsamp_lenc_cPT(iyear)=maxSS_lenc;}
            if (nsamp_lenc_cPT(iyear)>=minSS_lenc)
             {nsamp_cPT_lenc_allyr(yrs_lenc_cPT(iyear))=nsamp_lenc_cPT(iyear);
              nfish_cPT_lenc_allyr(yrs_lenc_cPT(iyear))=nfish_lenc_cPT(iyear);}}
      
      for (iyear=1; iyear<=nyr_lenc_rHB; iyear++)
         {  if (nsamp_lenc_rHB(iyear)>maxSS_lenc)
             {nsamp_lenc_rHB(iyear)=maxSS_lenc;}
            if (nsamp_lenc_rHB(iyear)>=minSS_lenc)
             {nsamp_rHB_lenc_allyr(yrs_lenc_rHB(iyear))=nsamp_lenc_rHB(iyear);
              nfish_rHB_lenc_allyr(yrs_lenc_rHB(iyear))=nfish_lenc_rHB(iyear);}}
      for (iyear=1; iyear<=nyr_lenc_rHB_D; iyear++)
         {  if (nsamp_lenc_rHB_D(iyear)>maxSS_lenc)
             {nsamp_lenc_rHB_D(iyear)=maxSS_lenc;}
            if (nsamp_lenc_rHB_D(iyear)>=minSS_lenc)
             {nsamp_rHB_D_lenc_allyr(yrs_lenc_rHB_D(iyear))=nsamp_lenc_rHB_D(iyear);
              nfish_rHB_D_lenc_allyr(yrs_lenc_rHB_D(iyear))=nfish_lenc_rHB_D(iyear);}}
         
      for (iyear=1; iyear<=nyr_lenc_rGN; iyear++)
         {  if (nsamp_lenc_rGN(iyear)>maxSS_lenc)
             {nsamp_lenc_rGN(iyear)=maxSS_lenc;}
            if (nsamp_lenc_rGN(iyear)>=minSS_lenc)
             {nsamp_rGN_lenc_allyr(yrs_lenc_rGN(iyear))=nsamp_lenc_rGN(iyear);
              nfish_rGN_lenc_allyr(yrs_lenc_rGN(iyear))=nfish_lenc_rGN(iyear);}}
                  
      for (iyear=1; iyear<=nyr_agec_sBT; iyear++)
         {  if (nsamp_agec_sBT(iyear)>maxSS_agec)
             {nsamp_agec_sBT(iyear)=maxSS_agec;}
            if (nsamp_agec_sBT(iyear)>=minSS_agec)
             {nsamp_sBT_agec_allyr(yrs_agec_sBT(iyear))=nsamp_agec_sBT(iyear);
              nfish_sBT_agec_allyr(yrs_agec_sBT(iyear))=nfish_agec_sBT(iyear);}} 
      for (iyear=1; iyear<=nyr_agec_sTV; iyear++)
         {  if (nsamp_agec_sTV(iyear)>maxSS_agec)
             {nsamp_agec_sTV(iyear)=maxSS_agec;}
            if (nsamp_agec_sTV(iyear)>=minSS_agec)
             {nsamp_sTV_agec_allyr(yrs_agec_sTV(iyear))=nsamp_agec_sTV(iyear);
              nfish_sTV_agec_allyr(yrs_agec_sTV(iyear))=nfish_agec_sTV(iyear);}} 
      for (iyear=1; iyear<=nyr_agec_cHL; iyear++)
         {  if (nsamp_agec_cHL(iyear)>maxSS_agec)
             {nsamp_agec_cHL(iyear)=maxSS_agec;}
            if (nsamp_agec_cHL(iyear)>=minSS_agec)
             {nsamp_cHL_agec_allyr(yrs_agec_cHL(iyear))=nsamp_agec_cHL(iyear);
              nfish_cHL_agec_allyr(yrs_agec_cHL(iyear))=nfish_agec_cHL(iyear);}} 
            
      for (iyear=1; iyear<=nyr_agec_cPT; iyear++)
         {  if (nsamp_agec_cPT(iyear)>maxSS_agec)
             {nsamp_agec_cPT(iyear)=maxSS_agec;}
            if (nsamp_agec_cPT(iyear)>=minSS_agec)
             {nsamp_cPT_agec_allyr(yrs_agec_cPT(iyear))=nsamp_agec_cPT(iyear);
              nfish_cPT_agec_allyr(yrs_agec_cPT(iyear))=nfish_agec_cPT(iyear);}}
             
      for (iyear=1; iyear<=nyr_agec_rHB; iyear++)
         {  if (nsamp_agec_rHB(iyear)>maxSS_agec)
             {nsamp_agec_rHB(iyear)=maxSS_agec;}
            if (nsamp_agec_rHB(iyear)>=minSS_agec)
             {nsamp_rHB_agec_allyr(yrs_agec_rHB(iyear))=nsamp_agec_rHB(iyear);
              nfish_rHB_agec_allyr(yrs_agec_rHB(iyear))=nfish_agec_rHB(iyear);}}
             
      //for (iyear=1; iyear<=nyr_rGN_agec; iyear++)
      //   {  if (nsamp_rGN_agec(iyear)>maxSS_agec)
      //       {nsamp_rGN_agec(iyear)=maxSS_agec;}
      //      if (nsamp_rGN_agec(iyear)>=minSS_agec)
      //       {nsamp_rGN_agec_allyr(yrs_rGN_agec(iyear))=nsamp_rGN_agec(iyear);
      //        nfish_rGN_agec_allyr(yrs_rGN_agec(iyear))=nfish_rGN_agec(iyear);}}  
  F_msy(1)=0.0;  
    for (ff=2;ff<=n_iter_msy;ff++) {F_msy(ff)=F_msy(ff-1)+iter_inc_msy;}
  F_spr(1)=0.0;  
    for (ff=2;ff<=n_iter_spr;ff++){F_spr(ff)=F_spr(ff-1)+iter_inc_spr;}
  F_cHL.initialize(); 
  F_cPT.initialize(); 
  F_cTW.initialize(); 
  F_rHB.initialize(); 
  F_rGN.initialize();
  
  F_cGN_D.initialize(); 
  F_rHB_D.initialize(); 
  F_rGN_D.initialize();  
      
  L_cHL_num.initialize(); 
  L_cPT_num.initialize(); 
  L_cTW_num.initialize(); 
  L_rHB_num.initialize(); 
  L_rGN_num.initialize();
  
  D_cGN_num.initialize(); 
  D_rHB_num.initialize(); 
  D_rGN_num.initialize();
      
  F_cHL_out.initialize(); 
  F_cPT_out.initialize(); 
  F_cTW_out.initialize(); 
  F_rHB_out.initialize(); 
  F_rGN_out.initialize();
  F_cGN_D_out.initialize(); 
  F_rHB_D_out.initialize(); 
  F_rGN_D_out.initialize();
  
  pred_cHL_L_klb.initialize();
  pred_cPT_L_klb.initialize();
  pred_cTW_L_klb.initialize(); 
  pred_rHB_L_klb.initialize(); 
  pred_rGN_L_klb.initialize();
  pred_cHL_L_knum.initialize();
  pred_cPT_L_knum.initialize();
  pred_cTW_L_knum.initialize(); 
  pred_rHB_L_knum.initialize(); 
  pred_rGN_L_knum.initialize();
  
  pred_cGN_D_klb.initialize(); 
  pred_rHB_D_klb.initialize();
  pred_rGN_D_klb.initialize();
  pred_cGN_D_knum.initialize(); 
  pred_rHB_D_knum.initialize(); 
  pred_rGN_D_knum.initialize();
 
  sel_sTV.initialize(); 
  sel_sBT.initialize();  
  sel_cHL.initialize(); 
  sel_cPT.initialize(); 
  sel_cTW.initialize();  
  sel_rHB.initialize(); 
  sel_rGN.initialize();  
  sel_cGN_D.initialize(); 
  sel_rHB_D.initialize(); 
  sel_rGN_D.initialize();  
  
  prob_belowsizelim_block1.initialize();
  prob_belowsizelim_block2.initialize();
  prob_belowsizelim_block3.initialize();
  prob_belowsizelim_block4.initialize();
  prob_belowsizelim_block5.initialize();
  
  log_rec_dev_output.initialize();
  log_Nage_dev_output.initialize();
  log_dev_rec.initialize();
  log_dev_Nage.initialize();
  
 
}

void model_parameters::userfunction(void)
{
  fval =0.0;
 R0=mfexp(log_R0);
 //cout<<"start"<<endl;
 //get_M_at_age(); //Needed only if M is estimated
 get_length_weight_at_age(); 
 //cout << "got length, weight, fecundity transitions" <<endl;
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
 //cout << "got discards" << endl;
 get_catchability_fcns();
 //cout << "got catchability_fcns" << endl;
 get_indices();
 //cout << "got indices" << endl;
 get_length_comps();
 //cout<< "got length comps"<< endl;
 get_age_comps();
 evaluate_objective_function();
 //cout << "objective function calculations complete" << endl;
}

void model_parameters::get_length_weight_at_age(void)
{
  //compute mean length (mm) and weight (whole) at age
    meanlen_TL=Linf*(1.0-mfexp(-K*(agebins-t0+0.5)));    //total length in mm
    wgt_g=wgtpar_a*pow(meanlen_TL,wgtpar_b);             //wgt in grams 
    wgt_kg=g2kg*wgt_g;                                   //wgt in kilograms 
    wgt_mt=g2mt*wgt_g;                                   //mt of whole wgt: g2mt converts g to mt
    wgt_klb=mt2klb*wgt_mt;                               //1000 lb of whole wgt
    wgt_lb=mt2lb*wgt_mt;                                 //1000 lb of whole wgt
    fecundity=fecpar_batches*mfexp(fecpar_a+wgt_g*fecpar_b)/fecpar_scale;    //fecundity at age, scaled
}

void model_parameters::get_reprod(void)
{
   //reprod is product of stuff going into reproductive capacity calcs
   reprod=elem_prod(elem_prod(prop_f,maturity_f),fecundity);  
   reprod2=elem_prod(elem_prod(prop_f,maturity_f),wgt_mt);  
}

void model_parameters::get_length_at_age_dist(void)
{
  //len_cv=len_cv_val; 
  //len_sd=elem_prod(len_cv, meanlen_TL);
  for (iage=1;iage<=nages;iage++)
   {len_cv(iage)=len_cv_val;
    len_sd(iage)=meanlen_TL(iage)*len_cv(iage);
    zscore_lzero=(0.0-meanlen_TL(iage))/len_sd(iage); 
	cprob_lzero=cumd_norm(zscore_lzero);	 
	//population
    zscore_len=((lenbins(1)+0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage);
    cprob_lenvec(1)=cumd_norm(zscore_len);          //includes any probability mass below zero
    lenprob(iage,1)=cprob_lenvec(1)-cprob_lzero;    //removes any probability mass below zero	
	//First size limit 8" Period 1 for both commercial and recreational discards (1984 through 1998)
	zscore_lsizelim1=(sizelim1-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim1=cumd_norm(zscore_lsizelim1);   //includes any probability mass below zero
	prob_belowsizelim_block1(iage)=	cprob_lsizelim1-cprob_lzero; //removes any probability mass below zero
	//Second size limit 10" Period 2 for both commercial and recreational discards and Period 3 for recreational landings (1999 through 2012 for commercial) (1999 through 2006 for recreational)
	zscore_lsizelim2=(sizelim2-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim2=cumd_norm(zscore_lsizelim2);   //includes any probability mass below zero
	prob_belowsizelim_block2(iage)=	cprob_lsizelim2-cprob_lzero; //removes any probability mass below zero
	//Third size limit 11" Period 3 for commercial (2013 through terminal year)
	zscore_lsizelim3=(sizelim3-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim3=cumd_norm(zscore_lsizelim3);                 //includes any probability mass below zero
	prob_belowsizelim_block3(iage)=	cprob_lsizelim3-cprob_lzero; //removes any probability mass below zero
	//Fourth size limit 12" Period 3 for recreational discards and Period 4 for recreational landings (2007 through 2012)
	zscore_lsizelim4=(sizelim4-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim4=cumd_norm(zscore_lsizelim4);                 //includes any probability mass below zero
	prob_belowsizelim_block4(iage)=	cprob_lsizelim4-cprob_lzero; //removes any probability mass below zero
	//Fifth size limit 13" Period 4 for recreational discards and Period 5 for recreational landings (2013 through the terminal year)
	zscore_lsizelim5=(sizelim5-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim5=cumd_norm(zscore_lsizelim5);                 //includes any probability mass below zero
	prob_belowsizelim_block5(iage)=	cprob_lsizelim5-cprob_lzero; //removes any probability mass below zero
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
  //fishery/fleet specific length probs, assumed normal prior to size limits
  lenprob_sBT=lenprob;
  lenprob_cHL1=lenprob; 
  lenprob_cHL2=lenprob;//_all;    //values may be adjusted based on size limit
  lenprob_cHL3=lenprob;//_all;    //values may be adjusted based on size limit 
                         ;//
  lenprob_cPT1=lenprob;   ;//
  lenprob_cPT2=lenprob;//_all;    //values may be adjusted based on size limit
  lenprob_cPT3=lenprob;//_all;    //values may be adjusted based on size limit
                         ;//
  lenprob_cTW1=lenprob;   ;//
  lenprob_cTW2=lenprob;//_all;    //values may be adjusted based on size limit
                         ;//
  lenprob_rHB1=lenprob;   ;//
  lenprob_rHB2=lenprob;//_all;    //values may be adjusted based on size limit 
  lenprob_rHB3=lenprob;//_all;    //values may be adjusted based on size limit 
  lenprob_rHB4=lenprob;//_all;    //values may be adjusted based on size limit
  lenprob_rHB5=lenprob;//_all;    //values may be adjusted based on size limit  
  lenprob_rGN1=lenprob; 
  lenprob_rGN2=lenprob;//_all;   //values may be adjusted based on size limit 
  lenprob_rGN3=lenprob;//_all;   //values may be adjusted based on size limit 
  lenprob_rGN4=lenprob;//_all;   //values may be adjusted based on size limit 
  lenprob_rGN5=lenprob;//_all;   //values may be adjusted based on size limit 
  lenprob_cGN_D2=lenprob;//_all; //values may be adjusted based on size limit 
  lenprob_cGN_D3=lenprob;//_all; //values may be adjusted based on size limit 
  lenprob_cGN_D4=lenprob;//_all; //values may be adjusted based on size limit 
  lenprob_rHB_D2=lenprob;//_all;   //values may be adjusted based on size limit 
  lenprob_rHB_D3=lenprob;//_all;   //values may be adjusted based on size limit 
  lenprob_rHB_D4=lenprob;//_all;   //values may be adjusted based on size limit
  lenprob_rHB_D5=lenprob;//_all;   //values may be adjusted based on size limit   
  lenprob_rGN_D1=lenprob;//_all; //values may be adjusted based on size limit 
  lenprob_rGN_D2=lenprob;//_all; //values may be adjusted based on size limit 
  lenprob_rGN_D3=lenprob;//_all; //values may be adjusted based on size limit 
  lenprob_rGN_D4=lenprob;//_all; //values may be adjusted based on size limit
  lenprob_rGN_D5=lenprob;//_all; //values may be adjusted based on size limit  
}

void model_parameters::get_weight_at_age_landings(void)
{
  //fleets under identical size limits are set equal at end of fcn
  for (iyear=styr; iyear<=endyr_period1; iyear++)
  {
    len_cHL_mm(iyear)=meanlen_TL;
    wgt_cHL_klb(iyear)=wgt_klb;      
    //len_cPT_mm(iyear)=meanlen_TL;  
    //wgt_cPT_klb(iyear)=wgt_klb;  
    //len_cTW_mm(iyear)=meanlen_TL;  
    //wgt_cTW_klb(iyear)=wgt_klb;  
    //len_rHB_mm(iyear)=meanlen_TL;
    //wgt_rHB_klb(iyear)=wgt_klb;  
    len_rGN_mm(iyear)=meanlen_TL;
    wgt_rGN_klb(iyear)=wgt_klb;  
    for (iage=1;iage<=nages; iage++)
    {
    len_rGN_D_mm(iyear,iage)=sum(elem_prod(lenprob_rGN_D2(iage),lenbins)); //assumes same size distn in period 1 as in period 2       
    }
    wgt_rGN_D_klb(iyear)=g2klb*wgtpar_a*pow(len_rGN_D_mm(iyear),wgtpar_b);
  } // end iyear loop
  for (iyear=(endyr_period1+1); iyear<=endyr_period2; iyear++)
  {
    for (iage=1;iage<=nages; iage++)
    {
    len_cHL_mm(iyear,iage)=sum(elem_prod(lenprob_cHL2(iage),lenbins)); 
    //len_cPT_mm(iyear,iage)=sum(elem_prod(lenprob_cPT2_all(iage),lenbins_all)); 
    //len_cTW_mm(iyear,iage)=sum(elem_prod(lenprob_cTW2_all(iage),lenbins_all));     
    //len_rHB_mm(iyear,iage)=sum(elem_prod(lenprob_rHB2_all(iage),lenbins_all));     
    len_rGN_mm(iyear,iage)=sum(elem_prod(lenprob_rGN2(iage),lenbins)); 
    len_cGN_D_mm(iyear,iage)=sum(elem_prod(lenprob_cGN_D2(iage),lenbins));
    //len_rHB_D_mm(iyear,iage)=sum(elem_prod(lenprob_rHB_D2_all(iage),lenbins_all));        
    len_rGN_D_mm(iyear,iage)=sum(elem_prod(lenprob_rGN_D2(iage),lenbins));             
    }
    wgt_cHL_klb(iyear)=g2klb*wgtpar_a*pow(len_cHL_mm(iyear),wgtpar_b);
    //wgt_cPT_klb(iyear)=g2klb*wgtpar_a*pow(len_cPT_mm(iyear),wgtpar_b);
    //wgt_cTW_klb(iyear)=g2klb*wgtpar_a*pow(len_cTW_mm(iyear),wgtpar_b);    
    //wgt_rHB_klb(iyear)=g2klb*wgtpar_a*pow(len_rHB_mm(iyear),wgtpar_b);
    wgt_rGN_klb(iyear)=g2klb*wgtpar_a*pow(len_rGN_mm(iyear),wgtpar_b);
    wgt_cGN_D_klb(iyear)=g2klb*wgtpar_a*pow(len_cGN_D_mm(iyear),wgtpar_b);  
    //wgt_rHB_D_klb(iyear)=g2klb*wgtpar_a*pow(len_rHB_D_mm(iyear),wgtpar_b);
    wgt_rGN_D_klb(iyear)=g2klb*wgtpar_a*pow(len_rGN_D_mm(iyear),wgtpar_b);               
  }
  for (iyear=(endyr_period2+1); iyear<=endyr; iyear++) //cGN only
  {
    for (iage=1;iage<=nages; iage++)
    {
    len_cHL_mm(iyear,iage)=sum(elem_prod(lenprob_cHL3(iage),lenbins)); 
    //len_cPT_mm(iyear,iage)=sum(elem_prod(lenprob_cPT3_all(iage),lenbins_all));      
    len_cGN_D_mm(iyear,iage)=sum(elem_prod(lenprob_cGN_D3(iage),lenbins));
    }
    wgt_cHL_klb(iyear)=g2klb*wgtpar_a*pow(len_cHL_mm(iyear),wgtpar_b);
    //wgt_cPT_klb(iyear)=g2klb*wgtpar_a*pow(len_cPT_mm(iyear),wgtpar_b);   
    wgt_cGN_D_klb(iyear)=g2klb*wgtpar_a*pow(len_cGN_D_mm(iyear),wgtpar_b);  
  }    
  for (iyear=(endyr_period2+1); iyear<=endyr_recr_period3; iyear++) //rec only
  {
    for (iage=1;iage<=nages; iage++)
    {
    //len_rHB_mm(iyear,iage)=sum(elem_prod(lenprob_rHB3_all(iage),lenbins_all));     
    len_rGN_mm(iyear,iage)=sum(elem_prod(lenprob_rGN3(iage),lenbins)); 
    //len_rHB_D_mm(iyear,iage)=sum(elem_prod(lenprob_rHB_D3_all(iage),lenbins_all));        
    len_rGN_D_mm(iyear,iage)=sum(elem_prod(lenprob_rGN_D3(iage),lenbins));             
    }
    //wgt_rHB_klb(iyear)=g2klb*wgtpar_a*pow(len_rHB_mm(iyear),wgtpar_b);
    wgt_rGN_klb(iyear)=g2klb*wgtpar_a*pow(len_rGN_mm(iyear),wgtpar_b);
    //wgt_rHB_D_klb(iyear)=g2klb*wgtpar_a*pow(len_rHB_D_mm(iyear),wgtpar_b);
    wgt_rGN_D_klb(iyear)=g2klb*wgtpar_a*pow(len_rGN_D_mm(iyear),wgtpar_b);               
  }
  for (iyear=(endyr_recr_period3+1); iyear<=endyr; iyear++) //rec only
  {
    for (iage=1;iage<=nages; iage++)
    {
    //len_rHB_mm(iyear,iage)=sum(elem_prod(lenprob_rHB4_all(iage),lenbins_all));     
    len_rGN_mm(iyear,iage)=sum(elem_prod(lenprob_rGN4(iage),lenbins)); 
    //len_rHB_D_mm(iyear,iage)=sum(elem_prod(lenprob_rHB_D4_all(iage),lenbins_all));        
    len_rGN_D_mm(iyear,iage)=sum(elem_prod(lenprob_rGN_D4(iage),lenbins));             
    }
    //wgt_rHB_klb(iyear)=g2klb*wgtpar_a*pow(len_rHB_mm(iyear),wgtpar_b);
    wgt_rGN_klb(iyear)=g2klb*wgtpar_a*pow(len_rGN_mm(iyear),wgtpar_b);
    //wgt_rHB_D_klb(iyear)=g2klb*wgtpar_a*pow(len_rHB_D_mm(iyear),wgtpar_b);
    wgt_rGN_D_klb(iyear)=g2klb*wgtpar_a*pow(len_rGN_D_mm(iyear),wgtpar_b);               
  }  
   for (iyear=(endyr_period4+1); iyear<=endyr; iyear++) //rec only
  {
    for (iage=1;iage<=nages; iage++)
    {
    //len_rHB_mm(iyear,iage)=sum(elem_prod(lenprob_rHB4_all(iage),lenbins_all));     
    len_rGN_mm(iyear,iage)=sum(elem_prod(lenprob_rGN5(iage),lenbins)); 
    //len_rHB_D_mm(iyear,iage)=sum(elem_prod(lenprob_rHB_D4_all(iage),lenbins_all));        
    len_rGN_D_mm(iyear,iage)=sum(elem_prod(lenprob_rGN_D4(iage),lenbins));             
    }
    //wgt_rHB_klb(iyear)=g2klb*wgtpar_a*pow(len_rHB_mm(iyear),wgtpar_b);
    wgt_rGN_klb(iyear)=g2klb*wgtpar_a*pow(len_rGN_mm(iyear),wgtpar_b);
    //wgt_rHB_D_klb(iyear)=g2klb*wgtpar_a*pow(len_rHB_D_mm(iyear),wgtpar_b);
    wgt_rGN_D_klb(iyear)=g2klb*wgtpar_a*pow(len_rGN_D_mm(iyear),wgtpar_b);               
  }  
    //identical fleets set equal here (for speed)      
    len_cPT_mm=len_cHL_mm;  wgt_cPT_klb=wgt_cHL_klb;  
    len_cTW_mm=len_cHL_mm;  wgt_cTW_klb=wgt_cHL_klb;  
    len_rHB_mm=len_rGN_mm;  wgt_rHB_klb=wgt_rGN_klb; 
    len_rHB_D_mm=len_rGN_D_mm;  wgt_rHB_D_klb=wgt_rGN_D_klb; 
  for (iyear=styr_cGN_closed; iyear<=endyr; iyear++) //overwrite last two yrs cGN discards, accnt for quotas
  { len_cGN_D_mm(iyear)=Dprop_cGN_sel_D*len_cGN_D_mm(iyear)+Dprop_cGN_sel_cHL*len_cHL_mm(iyear)+
                         Dprop_cGN_sel_cPT*len_cPT_mm(iyear);
    wgt_cGN_D_klb(iyear)=g2klb*wgtpar_a*pow(len_cGN_D_mm(iyear),wgtpar_b);                      
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
    N_spr_F0(iage)=N_spr_F0(iage-1)*
                   mfexp(-1.0*(M(iage-1)*(1.0-spawn_time_frac) + M(iage)*spawn_time_frac)); 
    N_bpr_F0(iage)=N_bpr_F0(iage-1)*mfexp(-1.0*(M(iage-1)));    
  }
  N_spr_F0(nages)=N_spr_F0(nages)/(1.0-mfexp(-1.0*M(nages))); //plus group (sum of geometric series)
  N_bpr_F0(nages)=N_bpr_F0(nages)/(1.0-mfexp(-1.0*M(nages)));
  spr_F0=sum(elem_prod(N_spr_F0,reprod));
  bpr_F0=sum(elem_prod(N_bpr_F0,wgt_mt));    
}

void model_parameters::get_selectivity(void)
{
  //---flat-topped sels---------------------------
  sel_sBT_vec=logistic(agebins, selpar_A50_sBT, selpar_slope_sBT);
  sel_sTV_vec=logistic(agebins, selpar_A50_sTV, selpar_slope_sTV);
  sel_cHL_2=logistic(agebins, selpar_A50_cHL2, selpar_slope_cHL2);
  sel_cHL_3=logistic(agebins, selpar_A50_cHL3, selpar_slope_cHL3);
  sel_cHL_4=logistic(agebins, selpar_A50_cHL4, selpar_slope_cHL4);
  sel_cPT_2=logistic(agebins, selpar_A50_cPT2, selpar_slope_cPT2);
  sel_cPT_3=logistic(agebins, selpar_A50_cPT3, selpar_slope_cPT3);
  sel_cPT_4=logistic(agebins, selpar_A50_cPT4, selpar_slope_cPT4);
  sel_rHB_1=logistic(agebins, selpar_A50_rHB1, selpar_slope_rHB1);
  sel_rHB_2=logistic(agebins, selpar_A50_rHB2, selpar_slope_rHB2);
  sel_rHB_3=logistic(agebins, selpar_A50_rHB3, selpar_slope_rHB3);
  sel_rHB_4=logistic(agebins, selpar_A50_rHB4, selpar_slope_rHB4);
  sel_rHB_5=logistic(agebins, selpar_A50_rHB5, selpar_slope_rHB5);
 // selpar_A50_rGN1=selpar_A50_rHB1;
 // selpar_slope_rGN1=selpar_slope_rHB1;
 // selpar_A50_rGN2=selpar_A50_rHB2;
 // selpar_slope_rGN2=selpar_slope_rHB2;
 // selpar_A50_rGN3=selpar_A50_rHB3;
 // selpar_slope_rGN3=selpar_slope_rHB3;    
 // selpar_A50_rGN4=selpar_A50_rHB4;
  sel_rGN1=logistic(agebins, selpar_A50_rGN1, selpar_slope_rGN1);
  sel_rGN2=logistic(agebins, selpar_A50_rGN2, selpar_slope_rGN2);
  sel_rGN3=logistic(agebins, selpar_A50_rGN3, selpar_slope_rGN3);
  sel_rGN4=logistic(agebins, selpar_A50_rGN4, selpar_slope_rGN4);
  sel_rGN5=logistic(agebins, selpar_A50_rGN5, selpar_slope_rGN5);
  //Period 1:   
  for (iyear=styr; iyear<=endyr_period1; iyear++)
  {
     sel_sBT(iyear)=sel_sBT_vec;
     sel_sTV(iyear)=sel_sTV_vec;
     sel_cHL(iyear)=sel_cHL_2; //commercial handline sel mirrors period 2
     sel_cPT(iyear)=sel_cPT_2; //commercial handline sel mirrors period 2
     sel_rHB(iyear)=sel_rHB_1; 
     sel_rGN(iyear)=sel_rGN1;      
  }
  //Period 2: 
  for (iyear=endyr_period1+1; iyear<=endyr_period2; iyear++)
  {     
     sel_sBT(iyear)=sel_sBT_vec;
     sel_sTV(iyear)=sel_sTV_vec;
     sel_cHL(iyear)=sel_cHL_2;
     sel_cPT(iyear)=sel_cPT_2; 
     sel_rHB(iyear)=sel_rHB_2;
     sel_rGN(iyear)=sel_rGN2;     
  }
  //Period 3 
  for (iyear=endyr_period2+1; iyear<=endyr; iyear++)
  {
     sel_sBT(iyear)=sel_sBT_vec;
     sel_sTV(iyear)=sel_sTV_vec;
     sel_cHL(iyear)=sel_cHL_3; 
     sel_cPT(iyear)=sel_cPT_3;     
     sel_rHB(iyear)=sel_rHB_3;
     sel_rGN(iyear)=sel_rGN3;     
  }   
  //Period 4: rec only, overwrites yrs calculated for period 3
  for (iyear=endyr_recr_period3; iyear<endyr_period4; iyear++)
  {
     sel_rHB(iyear)=sel_rHB_4;
     sel_rGN(iyear)=sel_rGN4;     
  }   
   //Period 5: Comm and rec, overwrites yrs previously calculated
   for (iyear=endyr_period4; iyear<=endyr; iyear++)
  {
     sel_cHL(iyear)=sel_cHL_4; 
     sel_cPT(iyear)=sel_cPT_4;
	 sel_rHB(iyear)=sel_rHB_5;
     sel_rGN(iyear)=sel_rGN5;     
  }   
  //set selectivities that mirror others
  sel_cTW=sel_cPT;  
  //sel_rGN=sel_rHB;
  //rGN and cGN mirror headboat discard selectivity
  selpar_Age0_rHB_D=1.0/(1.0+mfexp(-selpar_logit_Age0_rHB_D));
  selpar_Age1_rHB_D=1.0/(1.0+mfexp(-selpar_logit_Age1_rHB_D)); 
  selpar_Age2_rHB_D=1.0/(1.0+mfexp(-selpar_logit_Age2_rHB_D));
  sel_rHB_D_4=logistic_exponential(agebins, selpar_A50_rHB_D4, selpar_slope_rHB_D4, selpar_A502_rHB_D4, selpar_slope2_rHB_D4);//logistic_double(agebins, selpar_A50_rHB_D4, selpar_slope_rHB_D4, selpar_A502_rHB_D4, selpar_slope2_rHB_D4); 
  sel_rHB_D_5=logistic_exponential(agebins, selpar_A50_rHB_D5, selpar_slope_rHB_D5, selpar_A502_rHB_D5, selpar_slope2_rHB_D5);//logistic_double(agebins, selpar_A50_rHB_D5, selpar_slope_rHB_D5, selpar_A502_rHB_D5, selpar_slope2_rHB_D5); 
 //Assume same sel of age 0's across periods 
  sel_rHB_D_2(1)=selpar_Age0_rHB_D; 
  sel_rHB_D_3(1)=selpar_Age0_rHB_D; 
  //sel_rHB_D_4(1)=selpar_Age0_rHB_D;
  //sel_rHB_D_5(1)=selpar_Age0_rHB_D;
  sel_cGN_D_3(1)=selpar_Age0_rHB_D;
 //Assume same sel of age 1's across periods 
  sel_rHB_D_2(2)=selpar_Age1_rHB_D; 
  sel_rHB_D_3(2)=selpar_Age1_rHB_D; 
  //sel_rHB_D_4(2)=selpar_Age1_rHB_D;
  //sel_rHB_D_5(2)=selpar_Age1_rHB_D;
  sel_cGN_D_3(2)=selpar_Age1_rHB_D;
 //Assume same sel of age 2's across periods 
  sel_rHB_D_2(3)=selpar_Age2_rHB_D; 
  sel_rHB_D_3(3)=selpar_Age2_rHB_D; 
  //sel_rHB_D_4(3)=selpar_Age2_rHB_D;
  //sel_rHB_D_5(3)=selpar_Age2_rHB_D;  
  sel_cGN_D_3(3)=selpar_Age2_rHB_D;
 //Assume full sel at age 3 across periods 
  sel_rHB_D_2(4)=1.0; 
  sel_rHB_D_3(4)=1.0; 
  //sel_rHB_D_4(4)=1.0;
  //sel_rHB_D_5(4)=1.0; 
  sel_cGN_D_3(4)=1.0;  
  for (iage=5; iage<=nages; iage++)
      {sel_rHB_D_2(iage)=prob_belowsizelim_block1(iage);
       sel_rHB_D_3(iage)=prob_belowsizelim_block2(iage);
	   sel_cGN_D_3(iage)=prob_belowsizelim_block3(iage);
       //sel_rHB_D_4(iage)=prob_belowsizelim_block4(iage);
	   //sel_rHB_D_4(iage)=prob_belowsizelim_block5(iage);
       }
  //Period 1: assumed same as in period 1, no commercial discards
  for (iyear=styr; iyear<=endyr_period1; iyear++)
      {sel_rHB_D(iyear)=sel_rHB_D_2;}
  //Period 2: 
  for (iyear=endyr_period1+1; iyear<=endyr_period2; iyear++)
      {sel_rHB_D(iyear)=sel_rHB_D_2;
       sel_cGN_D(iyear)=sel_rHB_D_2;
      }
  //Period 3: Starts in 1999
  for (iyear=endyr_period2+1; iyear<=endyr; iyear++)
  {sel_rHB_D(iyear)=sel_rHB_D_3;
   sel_cGN_D(iyear)=sel_rHB_D_3;}  
  //Period 4: Starts in 2007, rHB and rGN only, overwrites last few yrs calculated for period 3
  for (iyear=endyr_recr_period3+1; iyear<=endyr_period4; iyear++)
  {sel_rHB_D(iyear)=sel_rHB_D_4;
   sel_cGN_D(iyear)=sel_rHB_D_3;
  }
  //Discard quota: wghted average,, overwrites last few yrs calculated for period 3 styr cGN closed=2009
  for (iyear=styr_cGN_closed; iyear<=endyr_period4; iyear++)
  {sel_cGN_D(iyear)=Dprop_cGN_sel_D*sel_rHB_D_3 + Dprop_cGN_sel_cHL*sel_cHL_2 +
                     Dprop_cGN_sel_cPT*sel_cPT_2;
   sel_cGN_D(iyear)=sel_cGN_D(iyear)/max(sel_cGN_D(iyear));} 
  //cout<<"sel_rHB_D after period 4 loop"<<sel_rHB_D<<endl;
  //Period 5: Starts in 2013 
  for (iyear=endyr_period4+1; iyear<=endyr; iyear++)
  {sel_rHB_D(iyear)=sel_rHB_D_5;
   sel_cGN_D(iyear)=sel_cGN_D_3;
   //Discard quota: wghted average,, overwrites last few yrs calculated for period 3 styr cGN closed=2009
   sel_cGN_D(iyear)=Dprop_cGN_sel_D*sel_cGN_D_3 + Dprop_cGN_sel_cHL*sel_cHL_2 +
                     Dprop_cGN_sel_cPT*sel_cPT_2;
   sel_cGN_D(iyear)=sel_cGN_D(iyear)/max(sel_cGN_D(iyear));} 
   sel_rGN_D=sel_rHB_D;
}

void model_parameters::get_mortality(void)
{
  Fsum.initialize();
  Fapex.initialize();
  F.initialize();
  //initialization F is avg from first 3 yrs of observed landings 
  log_F_dev_init_cHL=sum(log_dev_F_L_cHL(styr_L_cHL,(styr_L_cHL+2)))/3.0;
  log_F_dev_init_cPT=sum(log_dev_F_L_cHL(styr_L_cPT,(styr_L_cPT+2)))/3.0;
  log_F_dev_init_cTW=sum(log_dev_F_L_cHL(styr_L_cTW,(styr_L_cTW+2)))/3.0;
  log_F_init_rHB=sum(log_dev_F_L_rHB(styr_L_rHB,(styr_L_rHB+2)))/3.0;
  log_F_dev_init_rGN=sum(log_dev_F_L_rGN(styr_L_rGN,(styr_L_rGN+2)))/3.0;
  log_F_dev_init_rGN_D=sum(log_dev_F_D_rGN(styr_D_rGN,(styr_D_rGN+2)))/3.0;
  log_F_dev_cGN_D2=sum(log_dev_F_D_cGN(styr_D_cHL,(styr_D_cHL+5)))/6.0; //for cGN D 1984-1992
  //cout<<styr<<endl;  
  for (iyear=styr; iyear<=endyr; iyear++) 
  {
    //-------------
    if(iyear>=styr_L_cHL & iyear<=endyr_L_cHL)
    {  F_cHL_out(iyear)=mfexp(log_avg_F_L_cHL+log_dev_F_L_cHL(iyear)); //}    
       //if (iyear<styr_L_cHL){F_cHL_out(iyear)=mfexp(log_avg_F_L_cHL+log_F_dev_init_cHL);}        
       F_cHL(iyear)=sel_cHL(iyear)*F_cHL_out(iyear);
       Fsum(iyear)+=F_cHL_out(iyear);
    }
    //-------------
    if(iyear>=styr_L_cPT & iyear<=endyr_L_cPT)
    {  F_cPT_out(iyear)=mfexp(log_avg_F_L_cPT+log_dev_F_L_cPT(iyear)); //}
       //if (iyear<styr_L_cPT) {F_cPT_out(iyear)=0.0;}
       F_cPT(iyear)=sel_cPT(iyear)*F_cPT_out(iyear);
       Fsum(iyear)+=F_cPT_out(iyear);
    }
    //-------------
    if(iyear>=styr_L_cTW & iyear<=endyr_L_cTW)
    {  F_cTW_out(iyear)=mfexp(log_avg_F_L_cTW+log_dev_F_L_cTW(iyear)); //}
      // if (iyear<styr_L_cTW) {F_cTW_out(iyear)=0.0;}
       F_cTW(iyear)=sel_cTW(iyear)*F_cTW_out(iyear);
       Fsum(iyear)+=F_cTW_out(iyear);
    } 
    //-------------
    if(iyear>=styr_L_rHB & iyear<=endyr_L_rHB)
    {  F_rHB_out(iyear)=mfexp(log_avg_F_L_rHB+log_dev_F_L_rHB(iyear));//}
    // if (iyear<styr_L_rHB){F_rHB_out(iyear)=mfexp(log_avg_F_L_rHB+log_F_init_rHB);}
       F_rHB(iyear)=sel_rHB(iyear)*F_rHB_out(iyear);    
       Fsum(iyear)+=F_rHB_out(iyear);
    }
    //-------------
    if(iyear>=styr_L_rGN & iyear<=endyr_L_rGN)
       {F_rGN_out(iyear)=mfexp(log_avg_F_L_rGN+log_dev_F_L_rGN(iyear));}
    if (iyear<styr_L_rGN){F_rGN_out(iyear)=mfexp(log_avg_F_L_rGN+log_F_dev_init_rGN);}
    F_rGN(iyear)=sel_rGN(iyear)*F_rGN_out(iyear);    
    Fsum(iyear)+=F_rGN_out(iyear);
    //discards-------------
    if(iyear>=styr_D_cHL & iyear<=endyr_D_cHL)
      {F_cGN_D_out(iyear)=mfexp(log_avg_F_D_cGN+log_dev_F_D_cGN(iyear));}
    if(iyear > endyr_period1 & iyear < styr_D_cHL)
      {F_cGN_D_out(iyear)=mfexp(log_avg_F_D_cGN+log_F_dev_cGN_D2);}
    F_cGN_D(iyear)=sel_cGN_D(iyear)*F_cGN_D_out(iyear);  
    Fsum(iyear)+=F_cGN_D_out(iyear);    
    if(iyear>=styr_D_rHB & iyear<=endyr_D_rHB)
    {  F_rHB_D_out(iyear)=mfexp(log_avg_F_D_rHB+log_dev_F_D_rHB(iyear)); 
       F_rHB_D(iyear)=sel_rHB_D(iyear)*F_rHB_D_out(iyear); 
       Fsum(iyear)+=F_rHB_D_out(iyear); 
    }
    if(iyear<styr_D_rGN)
       {F_rGN_D_out(iyear)=mfexp(log_avg_F_D_rGN+log_F_dev_init_rGN_D);}
    if(iyear>=styr_D_rGN& iyear<=endyr_D_rGN)
       { F_rGN_D_out(iyear)=mfexp(log_avg_F_D_rGN+log_dev_F_D_rGN(iyear));}
    F_rGN_D(iyear)=sel_rGN_D(iyear)*F_rGN_D_out(iyear);  
    Fsum(iyear)+=F_rGN_D_out(iyear);
    //Total F at age
    F(iyear)=F_cHL(iyear); //first in additive series (NO +=)
    F(iyear)+=F_cPT(iyear);
    F(iyear)+=F_cTW(iyear);    
    F(iyear)+=F_rHB(iyear);
    F(iyear)+=F_rGN(iyear);
    F(iyear)+=F_cGN_D(iyear);
    F(iyear)+=F_rHB_D(iyear);
    F(iyear)+=F_rGN_D(iyear);
    Fapex(iyear)=max(F(iyear));
   // cout <<"Fapex" << Fapex <<endl;
    Z(iyear)=M+F(iyear);
  }  //end iyear
}

void model_parameters::get_bias_corr(void)
{
  //may exclude last BiasCor_exclude_yrs yrs bc constrained or lack info to estimate
  var_rec_dev=norm2(log_dev_rec(styr_rec_dev,(endyr-BiasCor_exclude_yrs))-
              sum(log_dev_rec(styr_rec_dev,(endyr-BiasCor_exclude_yrs)))
              /(nyrs_rec-BiasCor_exclude_yrs))/(nyrs_rec-BiasCor_exclude_yrs-1.0); 
			  //cout<<"var rec dev"<<var_rec_dev<<endl;
  //if (set_BiasCor <= 0.0) {BiasCor=mfexp(var_rec_dev/2.0);}   //bias correction
  rec_sigma_sq=square(rec_sigma);
  rec_sigma_sqd2=rec_sigma_sq/2.0;
  if (set_BiasCor <= 0.0) {BiasCor=mfexp(rec_sigma_sqd2);}   //bias correction               
  else {BiasCor=set_BiasCor;}
}

void model_parameters::get_numbers_at_age(void)
{
  S0=spr_F0*R0;
  R_virgin=(R0/((5.0*steep-1.0)*spr_F0))*(BiasCor*4.0*steep*spr_F0-spr_F0*(1.0-steep));
  B0=bpr_F0*R_virgin;   
  B0_q_DD=R_virgin*sum(elem_prod(N_bpr_F0(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages))); 
  F_initial=sel_cHL(styr)*mfexp(log_avg_F_L_cHL+log_F_dev_init_cHL)+
            sel_cPT(styr)*mfexp(log_avg_F_L_cPT+log_F_dev_init_cPT)+
            sel_cTW(styr)*mfexp(log_avg_F_L_cTW+log_F_dev_init_cTW)+
            sel_rHB(styr)*mfexp(log_avg_F_L_rHB+log_F_init_rHB)+
            sel_rGN(styr)*mfexp(log_avg_F_L_rGN+log_F_dev_init_rGN)+
            sel_rGN_D(styr)*mfexp(log_avg_F_D_rGN+log_F_dev_init_rGN_D);
  Z_initial=M+F_init_ratio*F_initial;
  N_spr_initial(1)=1.0*mfexp(-1.0*Z_initial(1)*spawn_time_frac); //at peak spawning time;
  for (iage=2; iage<=nages; iage++)
    {
      N_spr_initial(iage)=N_spr_initial(iage-1)*
                   mfexp(-1.0*(Z_initial(iage-1)*(1.0-spawn_time_frac) + Z_initial(iage)*spawn_time_frac)); 
    }
  N_spr_initial(nages)=N_spr_initial(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
  spr_initial=sum(elem_prod(N_spr_initial,reprod));
  if (styr==styr_rec_dev) {R1=(R0/((5.0*steep-1.0)*spr_initial))*
                 (4.0*steep*spr_initial-spr_F0*(1.0-steep));} //without bias correction (deviation added later)
  else {R1=(R0/((5.0*steep-1.0)*spr_initial))*
                 (BiasCor*4.0*steep*spr_initial-spr_F0*(1.0-steep));} //with bias correction                 
  if(R1<10.0) {R1=10.0;} //Avoid negative (or unreasonably low) popn sizes during search algorithm
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
  MatFemB(styr)=sum(elem_prod(N_spawn(styr),reprod2));
  B_q_DD(styr)=sum(elem_prod(N(styr)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
  for (iyear=styr; iyear<endyr; iyear++)
  {
    if(iyear<(styr_rec_dev-1)) //recruitment follows S-R curve exactly
    {
        N(iyear+1,1)=0.0;   //no age 0's mature in SSB calculations, value replaced below for abundance calcs
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages));//plus group
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod)); 
        MatFemB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod2));   
        B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));   
        N(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(iyear+1));
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year                    
    }
    else   //recruitment follows S-R curve with lognormal deviation
    {       
        N(iyear+1,1)=0.0;   //no age 0's mature in SSB calculations, value replaced below for abundance calcs
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages));//plus group
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));
        MatFemB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod2));
        B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
        N(iyear+1,1)=SR_func(R0, steep, spr_F0, SSB(iyear+1))*mfexp(log_dev_rec(iyear+1));
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
    }
  } //end iyear
    ////last year (projection) cannot compute recuitment of age 0's past terminal yr bc spawning occurs at spawn_time_frac
    //N(endyr+1,1)=SR_func(R0, steep, spr_F0, SSB(endyr));  
    //N(endyr+1)(2,nages)=++elem_prod(N(endyr)(1,nages-1),(mfexp(-1.*Z(endyr)(1,nages-1))));
    //N(endyr+1,nages)+=N(endyr,nages)*mfexp(-1.*Z(endyr,nages));//plus group
    //SSB(endyr+1)=sum(elem_prod(N(endyr+1),reprod));
  rec=column(N,1);
  SdS0=SSB/S0;
}

void model_parameters::get_landings_numbers(void)
{
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      L_cHL_num(iyear,iage)=N(iyear,iage)*F_cHL(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_cPT_num(iyear,iage)=N(iyear,iage)*F_cPT(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_cTW_num(iyear,iage)=N(iyear,iage)*F_cTW(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_rHB_num(iyear,iage)=N(iyear,iage)*F_rHB(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_rGN_num(iyear,iage)=N(iyear,iage)*F_rGN(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
    }
    pred_cHL_L_knum(iyear)=sum(L_cHL_num(iyear))/1000.0;
    pred_cPT_L_knum(iyear)=sum(L_cPT_num(iyear))/1000.0;
    pred_cTW_L_knum(iyear)=sum(L_cTW_num(iyear))/1000.0;    
    pred_rHB_L_knum(iyear)=sum(L_rHB_num(iyear))/1000.0;
    pred_rGN_L_knum(iyear)=sum(L_rGN_num(iyear))/1000.0;
  }
}

void model_parameters::get_landings_wgt(void)
{
  for (iyear=styr; iyear<=endyr; iyear++)
  {    
    L_cHL_klb(iyear)=elem_prod(L_cHL_num(iyear),wgt_cHL_klb(iyear));     //in 1000 lb
    L_cPT_klb(iyear)=elem_prod(L_cPT_num(iyear),wgt_cPT_klb(iyear));     //in 1000 lb
    L_cTW_klb(iyear)=elem_prod(L_cTW_num(iyear),wgt_cTW_klb(iyear));     //in 1000 lb    
    L_rHB_klb(iyear)=elem_prod(L_rHB_num(iyear),wgt_rHB_klb(iyear));     //in 1000 lb
    L_rGN_klb(iyear)=elem_prod(L_rGN_num(iyear),wgt_rGN_klb(iyear));  //in 1000 lb    
    pred_cHL_L_klb(iyear)=sum(L_cHL_klb(iyear));
    pred_cPT_L_klb(iyear)=sum(L_cPT_klb(iyear));
    pred_cTW_L_klb(iyear)=sum(L_cTW_klb(iyear));
    pred_rHB_L_klb(iyear)=sum(L_rHB_klb(iyear));
    pred_rGN_L_klb(iyear)=sum(L_rGN_klb(iyear));
  }
}

void model_parameters::get_dead_discards(void)
{
  //dead discards at age (number fish) 
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      D_cGN_num(iyear,iage)=N(iyear,iage)*F_cGN_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      D_rHB_num(iyear,iage)=N(iyear,iage)*F_rHB_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      D_rGN_num(iyear,iage)=N(iyear,iage)*F_rGN_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);        
    }
    pred_cGN_D_knum(iyear)=sum(D_cGN_num(iyear))/1000.0;         //pred annual dead discards in 1000s (for matching data)
    pred_cGN_D_klb(iyear)=sum(elem_prod(D_cGN_num(iyear),wgt_cGN_D_klb(iyear)));  //annual dead discards in 1000 lb (for output only)
    pred_rHB_D_knum(iyear)=sum(D_rHB_num(iyear))/1000.0;             //pred annual dead discards in 1000s (for matching data)
    pred_rHB_D_klb(iyear)=sum(elem_prod(D_rHB_num(iyear),wgt_rHB_D_klb(iyear)));        //annual dead discards in 1000 lb (for output only)
    pred_rGN_D_knum(iyear)=sum(D_rGN_num(iyear))/1000.0;         //pred annual dead discards in 1000s (for matching data)
    pred_rGN_D_klb(iyear)=sum(elem_prod(D_rGN_num(iyear),wgt_rGN_D_klb(iyear)));  //annual dead discards in 1000 lb (for output only)
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
      for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++)
      {   if (iyear>styr_cpue_rHB & iyear <=2003) 
          {//q_rate_fcn_rHB(iyear)=(1.0+q_rate)*q_rate_fcn_rHB(iyear-1); //compound
             q_rate_fcn_rHB(iyear)=(1.0+(iyear-styr_cpue_rHB)*q_rate)*q_rate_fcn_rHB(styr_cpue_rHB);  //linear
          }
          if (iyear>2003) {q_rate_fcn_rHB(iyear)=q_rate_fcn_rHB(iyear-1);} 
      }  
      //for (iyear=styr_rHB_D_cpue; iyear<=endyr_rHB_D_cpue; iyear++)
      //{   if (iyear>styr_rHB_D_cpue & iyear <=2003) 
      //    {//q_rate_fcn_rHB_D(iyear)=(1.0+q_rate)*q_rate_fcn_rHB_D(iyear-1); //compound
      //       q_rate_fcn_rHB_D(iyear)=(1.0+(iyear-styr_rHB_D_cpue)*q_rate)*q_rate_fcn_rHB_D(styr_rHB_D_cpue);  //linear
      //    }
      //    if (iyear>2003) {q_rate_fcn_rHB_D(iyear)=q_rate_fcn_rHB_D(iyear-1);} 
      //}
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
 //Survey 1: sBT
  for (iyear=styr_cpue_sBT; iyear<=endyr_cpue_sBT; iyear++)
  {   //index in number units
      N_sBT(iyear)=elem_prod(N_mdyr(iyear),sel_sBT(iyear)); 
      pred_sBT_cpue(iyear)=mfexp(log_q_sBT)*sum(N_sBT(iyear));
  }
 //Survey 2: sTV
  for (iyear=styr_cpue_sTV; iyear<=endyr_cpue_sTV; iyear++)
  {   //index in number units
      N_sTV(iyear)=elem_prod(N_mdyr(iyear),sel_sTV(iyear)); 
      pred_sTV_cpue(iyear)=mfexp(log_q_sTV)*sum(N_sTV(iyear));
  }
  //Survey 3: Vid
  //for (iyear=styr_Vid_cpue; iyear<=endyr_Vid_cpue; iyear++)
  //{   //index in number units
  //    N_Vid(iyear)=elem_prod(N_mdyr(iyear),sel_sTV(iyear)); 
  //    pred_Vid_cpue(iyear)=mfexp(log_q_Vid)*sum(N_Vid(iyear));
  //}
 //Commercial handline cpue
  q_cHL(styr_cpue_cHL)=mfexp(log_q_cHL); 
  for (iyear=styr_cpue_cHL; iyear<=endyr_cpue_cHL; iyear++)
  {   //index in weight units. original index in lb and re-scaled. predicted in klb, but difference is absorbed by q
      N_cHL(iyear)=elem_prod(elem_prod(N_mdyr(iyear),sel_cHL(iyear)),wgt_cHL_klb(iyear)); 
      pred_cHL_cpue(iyear)=q_cHL(iyear)*q_rate_fcn_cHL(iyear)*q_DD_fcn(iyear)*sum(N_cHL(iyear));
      if (iyear<endyr_cpue_cHL){q_cHL(iyear+1)=q_cHL(iyear)*mfexp(q_RW_log_dev_cHL(iyear));}
  }
 //Headboat cpue
  q_rHB(styr_cpue_rHB)=mfexp(log_q_rHB);
  for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++)
  {   //index in weight units. original index in lb and re-scaled. predicted in klb, but difference is absorbed by q
      N_rHB(iyear)=elem_prod(elem_prod(N_mdyr(iyear),sel_rHB(iyear)),wgt_rHB_klb(iyear)); 
      pred_rHB_cpue(iyear)=q_rHB(iyear)*q_rate_fcn_rHB(iyear)*q_DD_fcn(iyear)*sum(N_rHB(iyear));
      if (iyear<endyr_cpue_rHB){q_rHB(iyear+1)=q_rHB(iyear)*mfexp(q_RW_log_dev_rHB(iyear));}
  }
  //rHB_D cpue
  //q_rHB_D(styr_rHB_D_cpue)=mfexp(log_q_rHB_D);
  //for (iyear=styr_rHB_D_cpue; iyear<=endyr_rHB_D_cpue; iyear++)
  //{   //index in number units
  //    N_rHB_D(iyear)=elem_prod(N_mdyr(iyear),sel_rHB_D(iyear)); 
  //    pred_rHB_D_cpue(iyear)=q_rHB_D(iyear)*q_rate_fcn_rHB_D(iyear)*q_DD_fcn(iyear)*sum(N_rHB_D(iyear));
  //    if (iyear<endyr_rHB_D_cpue){q_rHB_D(iyear+1)=q_rHB_D(iyear)*mfexp(q_RW_log_dev_rHB_D(iyear));}
  //}
}

void model_parameters::get_length_comps(void)
{
  //sBT
  for (iyear=1;iyear<=nyr_lenc_sBT;iyear++)
      {pred_sBT_lenc(iyear)=(N_sBT(yrs_lenc_sBT(iyear))*lenprob)/sum(N_sBT(yrs_lenc_sBT(iyear)));} 
  //Commercial lines
  for (iyear=1;iyear<=nyr_lenc_cHL;iyear++) //all yrs within periods 2,3
  {  if (yrs_lenc_cHL(iyear)<=endyr_period2)
     {pred_cHL_lenc(iyear)=(L_cHL_num(yrs_lenc_cHL(iyear))*lenprob_cHL2)
                          /sum(L_cHL_num(yrs_lenc_cHL(iyear)));       
     } 
     if (yrs_lenc_cHL(iyear)>endyr_period2)
     {pred_cHL_lenc(iyear)=(L_cHL_num(yrs_lenc_cHL(iyear))*lenprob_cHL3)
                          /sum(L_cHL_num(yrs_lenc_cHL(iyear)));       
     }  
  }
  //Commercial pots: pooled all from period 2
  L_cPT_num_pool.initialize();
  for (iyear=1;iyear<=nyr_lenc_pool_cPT;iyear++)
  {  L_cPT_num_pool_yr(iyear)=nsamp_lenc_pool_cPT(iyear)*L_cPT_num(yrs_lenc_pool_cPT(iyear))
                            /sum(L_cPT_num(yrs_lenc_pool_cPT(iyear)));                             
     if (yrs_lenc_pool_cPT(iyear)<=endyr_period2) {L_cPT_num_pool(1)+=L_cPT_num_pool_yr(iyear);}                             
  } 
  for (iyear=1;iyear<=nyr_lenc_cPT;iyear++) //all yrs within periods 2
  {  if (yrs_lenc_cPT(iyear)<=endyr_period2)
       {pred_cPT_lenc(iyear)=(L_cPT_num_pool(iyear)*lenprob_cPT2)/sum(L_cPT_num_pool(iyear)); } 
	   pred_cPT_lenc(iyear)=(L_cPT_num(yrs_lenc_cPT(iyear))*lenprob_cPT3)  //added to calculate comps for last period
						  /sum(L_cPT_num(yrs_lenc_cPT(iyear)));
  }  
 //Headboat 
  for (iyear=1;iyear<=nyr_lenc_rHB;iyear++)  //all in periods 1,2,3
  {  if (yrs_lenc_rHB(iyear)<=endyr_period1)
     {pred_rHB_lenc(iyear)=(L_rHB_num(yrs_lenc_rHB(iyear))*lenprob_rHB1)
                          /sum(L_rHB_num(yrs_lenc_rHB(iyear)));       
     } 
     if (yrs_lenc_rHB(iyear)>endyr_period1 & yrs_lenc_rHB(iyear)<=endyr_period2)
     {pred_rHB_lenc(iyear)=(L_rHB_num(yrs_lenc_rHB(iyear))*lenprob_rHB2)
                          /sum(L_rHB_num(yrs_lenc_rHB(iyear)));       
     }  
     if (yrs_lenc_rHB(iyear)>endyr_period2)
     {pred_rHB_lenc(iyear)=(L_rHB_num(yrs_lenc_rHB(iyear))*lenprob_rHB3)
                          /sum(L_rHB_num(yrs_lenc_rHB(iyear)));       
     }  
  }
 //rHB discards 
  for (iyear=1;iyear<=nyr_lenc_rHB_D;iyear++) //all yrs within period 3,4
  {  if (yrs_lenc_rHB_D(iyear)<=endyr_recr_period3)
     {pred_rHB_D_lenc(iyear)=(D_rHB_num(yrs_lenc_rHB_D(iyear))*lenprob_rHB_D3)
                        /sum(D_rHB_num(yrs_lenc_rHB_D(iyear)));}      
    if (yrs_lenc_rHB_D(iyear)>endyr_recr_period3)
     {pred_rHB_D_lenc(iyear)=(D_rHB_num(yrs_lenc_rHB_D(iyear))*lenprob_rHB_D4)
                        /sum(D_rHB_num(yrs_lenc_rHB_D(iyear)));}                      
  }
 //MRIP
  for (iyear=1;iyear<=nyr_lenc_rGN;iyear++)  //all in periods 1,2,3
  {  if (yrs_lenc_rGN(iyear)<=endyr_period1)
     {pred_rGN_lenc(iyear)=(L_rGN_num(yrs_lenc_rGN(iyear))*lenprob_rGN1)
                          /sum(L_rGN_num(yrs_lenc_rGN(iyear)));       
     } 
     if (yrs_lenc_rGN(iyear)>endyr_period1 & yrs_lenc_rGN(iyear)<=endyr_period2)
     {pred_rGN_lenc(iyear)=(L_rGN_num(yrs_lenc_rGN(iyear))*lenprob_rGN2)
                          /sum(L_rGN_num(yrs_lenc_rGN(iyear)));       
     }  
     if (yrs_lenc_rGN(iyear)>endyr_period2 & yrs_lenc_rGN(iyear)<=endyr_recr_period3)
     {pred_rGN_lenc(iyear)=(L_rGN_num(yrs_lenc_rGN(iyear))*lenprob_rGN3)
                          /sum(L_rGN_num(yrs_lenc_rGN(iyear)));       
     }  
     if (yrs_lenc_rGN(iyear)>endyr_recr_period3 &
	 yrs_lenc_rGN(iyear)<=endyr_period4)
     {pred_rGN_lenc(iyear)=(L_rGN_num(yrs_lenc_rGN(iyear))*lenprob_rGN4)
                          /sum(L_rGN_num(yrs_lenc_rGN(iyear)));       
     }  
	 if (yrs_lenc_rGN(iyear)>endyr_period4)
     {pred_rGN_lenc(iyear)=(L_rGN_num(yrs_lenc_rGN(iyear))*lenprob_rGN5)
                          /sum(L_rGN_num(yrs_lenc_rGN(iyear)));       
     }  
  } 
}

void model_parameters::get_age_comps(void)
{
 //MARMAP bft
  for (iyear=1;iyear<=nyr_agec_sBT;iyear++)
  {
    ErrorFree_sBT_agec(iyear)=N_sBT(yrs_agec_sBT(iyear))/sum(N_sBT(yrs_agec_sBT(iyear)));
    pred_sBT_agec(iyear)=age_error*ErrorFree_sBT_agec(iyear);                      
  }
 //MARMAP cvt
  for (iyear=1;iyear<=nyr_agec_sTV;iyear++)
  {
    ErrorFree_sTV_agec(iyear)=N_sTV(yrs_agec_sTV(iyear))/sum(N_sTV(yrs_agec_sTV(iyear)));
    pred_sTV_agec(iyear)=age_error*ErrorFree_sTV_agec(iyear);                      
  }
 //Commercial lines
  for (iyear=1;iyear<=nyr_agec_cHL;iyear++)
  {
    ErrorFree_cHL_agec(iyear)=L_cHL_num(yrs_agec_cHL(iyear))/sum(L_cHL_num(yrs_agec_cHL(iyear)));
    pred_cHL_agec(iyear)=age_error*ErrorFree_cHL_agec(iyear);                      
  }
 //Commercial pots
  for (iyear=1;iyear<=nyr_agec_cPT;iyear++)
  {
    ErrorFree_cPT_agec(iyear)=L_cPT_num(yrs_agec_cPT(iyear))/sum(L_cPT_num(yrs_agec_cPT(iyear)));
    pred_cPT_agec(iyear)=age_error*ErrorFree_cPT_agec(iyear);                      
  }
 //Headboat
  for (iyear=1;iyear<=nyr_agec_rHB;iyear++)
  {
    ErrorFree_rHB_agec(iyear)=L_rHB_num(yrs_agec_rHB(iyear))/sum(L_rHB_num(yrs_agec_rHB(iyear)));
    pred_rHB_agec(iyear)=age_error*ErrorFree_rHB_agec(iyear);                    
  }
 ////rGN
 // for (iyear=1;iyear<=nyr_rGN_agec;iyear++)
 // {
 //   ErrorFree_rGN_agec(iyear)=L_rGN_num(yrs_rGN_agec(iyear))/sum(L_rGN_num(yrs_rGN_agec(iyear)));
 //   pred_rGN_agec(iyear)=age_error*ErrorFree_rGN_agec(iyear);                    
 // }
}

void model_parameters::get_weighted_current(void)
{
  F_temp_sum=0.0;
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cHL+
        sum(log_dev_F_L_cHL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cPT+
        sum(log_dev_F_L_cPT((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rHB+
        sum(log_dev_F_L_rHB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rGN+
        sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_D_cGN+
        sum(log_dev_F_D_cGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_D_rHB+
        sum(log_dev_F_D_rHB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_D_rGN+
        sum(log_dev_F_D_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
  F_cHL_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cHL+
        sum(log_dev_F_L_cHL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_cPT_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cPT+
        sum(log_dev_F_L_cPT((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_rHB_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rHB+
        sum(log_dev_F_L_rHB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_rGN_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rGN+
        sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_cGN_D_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_D_cGN+
        sum(log_dev_F_D_cGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_rHB_D_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_D_rHB+
        sum(log_dev_F_D_rHB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_rGN_D_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_D_rGN+
        sum(log_dev_F_D_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  log_F_dev_end_cHL=sum(log_dev_F_L_cHL((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_cPT=sum(log_dev_F_L_cPT((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  
  log_F_dev_end_rHB=sum(log_dev_F_L_rHB((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_rGN=sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_cGN_D=sum(log_dev_F_D_cGN((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_rHB_D=sum(log_dev_F_D_rHB((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_rGN_D=sum(log_dev_F_D_rGN((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  F_end_L=sel_cHL(endyr)*mfexp(log_avg_F_L_cHL+log_F_dev_end_cHL)+
          sel_cPT(endyr)*mfexp(log_avg_F_L_cPT+log_F_dev_end_cPT)+
          sel_rHB(endyr)*mfexp(log_avg_F_L_rHB+log_F_dev_end_rHB)+
          sel_rGN(endyr)*mfexp(log_avg_F_L_rGN+log_F_dev_end_rGN);
  F_end_D=sel_cGN_D(endyr)*mfexp(log_avg_F_D_cGN+log_F_dev_end_cGN_D)+
          sel_rHB_D(endyr)*mfexp(log_avg_F_D_rHB+log_F_dev_end_rHB_D)+
          sel_rGN_D(endyr)*mfexp(log_avg_F_D_rGN+log_F_dev_end_rGN_D);    
  F_end=F_end_L+F_end_D;
  F_end_apex=max(F_end);
  sel_wgted_tot=F_end/F_end_apex;
  sel_wgted_L=elem_prod(sel_wgted_tot, elem_div(F_end_L,F_end));
  sel_wgted_D=elem_prod(sel_wgted_tot, elem_div(F_end_D,F_end));
  wgt_wgted_L_denom=F_cHL_prop+F_cPT_prop+F_rHB_prop+F_rGN_prop;
  wgt_wgted_L_klb=F_cHL_prop/wgt_wgted_L_denom*wgt_cHL_klb(endyr)+
              F_cPT_prop/wgt_wgted_L_denom*wgt_cPT_klb(endyr)+
              F_rHB_prop/wgt_wgted_L_denom*wgt_rHB_klb(endyr)+
              F_rGN_prop/wgt_wgted_L_denom*wgt_rGN_klb(endyr);                
  wgt_wgted_D_denom=F_cGN_D_prop+F_rHB_D_prop+F_rGN_D_prop;
  wgt_wgted_D_klb=F_cGN_D_prop/wgt_wgted_D_denom*wgt_cGN_D_klb(endyr)+
              F_rHB_D_prop/wgt_wgted_D_denom*wgt_rHB_D_klb(endyr)+
              F_rGN_D_prop/wgt_wgted_D_denom*wgt_rGN_D_klb(endyr);                
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
    {
      N_age_msy(iage)=N_age_msy(iage-1)*mfexp(-1.*Z_age_msy(iage-1));
    }
    N_age_msy(nages)=N_age_msy(nages)/(1.0-mfexp(-1.*Z_age_msy(nages)));
    N_age_msy_spawn(1,(nages-1))=elem_prod(N_age_msy(1,(nages-1)),
                                   mfexp((-1.*Z_age_msy(1,(nages-1)))*spawn_time_frac));                 
    N_age_msy_spawn(nages)=(N_age_msy_spawn(nages-1)*
                          (mfexp(-1.*(Z_age_msy(nages-1)*(1.0-spawn_time_frac) + 
                                 Z_age_msy(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_msy(nages)));
    spr_msy(ff)=sum(elem_prod(N_age_msy_spawn,reprod));
    //Compute equilibrium values of R (including bias correction), SSB and Yield at each F
    R_eq(ff)=(R0/((5.0*steep-1.0)*spr_msy(ff)))*
                 (BiasCor*4.0*steep*spr_msy(ff)-spr_F0*(1.0-steep));
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
  }  //end ff loop
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
}

void model_parameters::get_per_recruit_stuff(void)
{
  //static per-recruit stuff
  for(iyear=styr; iyear<=endyr; iyear++)
  {
    N_age_spr(1)=1.0;
    for(iage=2; iage<=nages; iage++)
    {
      N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z(iyear,iage-1));
    }
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
    {
      N_age_spr(iage)=N_age_spr(iage-1)*mfexp(-1.*Z_age_spr(iage-1));
    }
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
      L_spr(ff)+=L_age_spr(iage)*wgt_wgted_L_klb(iage)*1000.0; //in lb
    }   
  }
  //--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
}

void model_parameters::get_miscellaneous_stuff(void)
{
  sigma_rec_dev=sqrt(var_rec_dev); //sample SD of predicted residuals (may not equal rec_sigma)  
  //len_cv=elem_div(len_sd,meanlen_TL);
  //compute total landings- and discards-at-age in 1000 fish and klb
  L_total_num.initialize();
  L_total_klb.initialize();
  D_total_num.initialize();
  D_total_klb.initialize();
  L_total_knum_yr.initialize();
  L_total_klb_yr.initialize();  
  D_total_knum_yr.initialize();
  D_total_klb_yr.initialize();
  for(iyear=styr; iyear<=endyr; iyear++)
  {
        L_total_klb_yr(iyear)= pred_cHL_L_klb(iyear)+pred_cPT_L_klb(iyear)+pred_cTW_L_klb(iyear)+
                               pred_rHB_L_klb(iyear)+pred_rGN_L_klb(iyear);
        L_total_knum_yr(iyear)=pred_cHL_L_knum(iyear)+pred_cPT_L_knum(iyear)+pred_cTW_L_knum(iyear)+
                               pred_rHB_L_knum(iyear)+pred_rGN_L_knum(iyear);
        D_total_knum_yr(iyear)+=pred_cGN_D_knum(iyear);
        D_total_klb_yr(iyear)+=pred_cGN_D_klb(iyear);
        D_total_knum_yr(iyear)+=pred_rHB_D_knum(iyear);
        D_total_klb_yr(iyear)+=pred_rHB_D_klb(iyear);
        D_total_knum_yr(iyear)+=pred_rGN_D_knum(iyear);
        D_total_klb_yr(iyear)+=pred_rGN_D_klb(iyear);
        D_cGN_klb(iyear)=elem_prod(D_cGN_num(iyear),wgt_cGN_D_klb(iyear));     //in 1000 lb
        D_rHB_klb(iyear)=elem_prod(D_rHB_num(iyear),wgt_rHB_D_klb(iyear));         //in 1000 lb
        D_rGN_klb(iyear)=elem_prod(D_rGN_num(iyear),wgt_rGN_D_klb(iyear));   //in 1000 lb 
        B(iyear)=elem_prod(N(iyear),wgt_mt);
        totN(iyear)=sum(N(iyear));
        totB(iyear)=sum(B(iyear));              
  }
  L_total_num=(L_cHL_num+L_cPT_num+L_cTW_num+L_rHB_num+L_rGN_num); //landings at age in number fish
  L_total_klb=L_cHL_klb+L_cPT_klb+L_cTW_klb+L_rHB_klb+L_rGN_klb;   //landings at age in klb whole weight
  D_total_num=(D_cGN_num+D_rHB_num+D_rGN_num);          //discards at age in number fish
  D_total_klb=D_cGN_klb+D_rHB_klb+D_rGN_klb;            //discards at age in klb whole weight
  //B(endyr+1)=elem_prod(N(endyr+1),wgt_mt);
  //totN(endyr+1)=sum(N(endyr+1));
  //totB(endyr+1)=sum(B(endyr+1));  
  if(F_msy_out>0)
    {
      FdF_msy=Fapex/F_msy_out;
      FdF_msy_end=FdF_msy(endyr_L_cHL);
      //cout<<"Fstatus2012"<< FdF_msy_end << endl;
      FdF_msy_end_mean=pow((FdF_msy(endyr_L_cHL)*FdF_msy(endyr_L_cHL-1)),(1.0/2.0));
      //FdF_msy_end_mean=pow((FdF_msy(endyr)*FdF_msy(endyr-1)*FdF_msy(endyr-2)),(1.0/3.0));
    }
  if(SSB_msy_out>0)
    {
      SdSSB_msy=SSB/SSB_msy_out;
      SdSSB_msy_end=SdSSB_msy(endyr);
    }  
   //fill in log recruitment deviations for yrs they are nonzero
   for(iyear=styr_rec_dev; iyear<=endyr; iyear++)
     {log_rec_dev_output(iyear)=log_dev_rec(iyear);}
   //fill in log Nage deviations for ages they are nonzero (ages2+)
   for(iage=2; iage<=nages; iage++)
   { log_Nage_dev_output(iage)=log_dev_Nage(iage);}
}

void model_parameters::get_effective_sample_sizes(void)
{
  neff_sBT_lenc_allyr_out=missing;//"missing" defined in admb2r.cpp
  neff_cHL_lenc_allyr_out=missing; 
  neff_cPT_lenc_allyr_out=missing;
  neff_rHB_lenc_allyr_out=missing;
  neff_rHB_D_lenc_allyr_out=missing;
  neff_rGN_lenc_allyr_out=missing;
  neff_sBT_agec_allyr_out=missing;
  neff_sTV_agec_allyr_out=missing;
  neff_cHL_agec_allyr_out=missing;
  neff_cPT_agec_allyr_out=missing;      
  neff_rHB_agec_allyr_out=missing;
  //neff_rGN_agec_allyr_out=missing;
      for (iyear=1; iyear<=nyr_lenc_sBT; iyear++)
         {if (nsamp_lenc_sBT(iyear)>=minSS_lenc)
            {neff_sBT_lenc_allyr_out(yrs_lenc_sBT(iyear))=multinom_eff_N(pred_sBT_lenc(iyear),obs_lenc_sBT(iyear));} 
          else {neff_sBT_lenc_allyr_out(yrs_lenc_sBT(iyear))=-99;}
         }                  
      for (iyear=1; iyear<=nyr_lenc_cHL; iyear++)
         {if (nsamp_lenc_cHL(iyear)>=minSS_lenc)
            {neff_cHL_lenc_allyr_out(yrs_lenc_cHL(iyear))=multinom_eff_N(pred_cHL_lenc(iyear),obs_lenc_cHL(iyear));} 
          else {neff_cHL_lenc_allyr_out(yrs_lenc_cHL(iyear))=-99;}
         }                  
      for (iyear=1; iyear<=nyr_lenc_cPT; iyear++)
         {if (nsamp_lenc_cPT(iyear)>=minSS_lenc)
             {neff_cPT_lenc_allyr_out(yrs_lenc_cPT(iyear))=multinom_eff_N(pred_cPT_lenc(iyear),obs_lenc_cPT(iyear));}                            
          else {neff_cPT_lenc_allyr_out(yrs_lenc_cPT(iyear))=-99;}
         }
      for (iyear=1; iyear<=nyr_lenc_rHB; iyear++)
         {if (nsamp_lenc_rHB(iyear)>=minSS_lenc)
             {neff_rHB_lenc_allyr_out(yrs_lenc_rHB(iyear))=multinom_eff_N(pred_rHB_lenc(iyear),obs_lenc_rHB(iyear));}                            
          else {neff_rHB_lenc_allyr_out(yrs_lenc_rHB(iyear))=-99;}
         }
      for (iyear=1; iyear<=nyr_lenc_rHB_D; iyear++)
         {if (nsamp_lenc_rHB_D(iyear)>=minSS_lenc)
            {neff_rHB_D_lenc_allyr_out(yrs_lenc_rHB_D(iyear))=multinom_eff_N(pred_rHB_D_lenc(iyear),obs_lenc_rHB_D(iyear));}                            
           else {neff_rHB_D_lenc_allyr_out(yrs_lenc_rHB_D(iyear))=-99;}
         }
      for (iyear=1; iyear<=nyr_lenc_rGN; iyear++)
         {if (nsamp_lenc_rGN(iyear)>=minSS_lenc)
            {neff_rGN_lenc_allyr_out(yrs_lenc_rGN(iyear))=multinom_eff_N(pred_rGN_lenc(iyear),obs_lenc_rGN(iyear)); }                           
           else {neff_rGN_lenc_allyr_out(yrs_lenc_rGN(iyear))=-99;}
         }
      for (iyear=1; iyear<=nyr_agec_sBT; iyear++)
         {if (nsamp_agec_sBT(iyear)>=minSS_agec)
            {neff_sBT_agec_allyr_out(yrs_agec_sBT(iyear))=multinom_eff_N(pred_sBT_agec(iyear),obs_agec_sBT(iyear));}                            
          else {neff_sBT_agec_allyr_out(yrs_agec_sBT(iyear))=-99;}
         }
      for (iyear=1; iyear<=nyr_agec_sTV; iyear++)
         {if (nsamp_agec_sTV(iyear)>=minSS_agec)
            {neff_sTV_agec_allyr_out(yrs_agec_sTV(iyear))=multinom_eff_N(pred_sTV_agec(iyear),obs_agec_sTV(iyear));}                            
          else {neff_sTV_agec_allyr_out(yrs_agec_sTV(iyear))=-99;}
         }
      for (iyear=1; iyear<=nyr_agec_cHL; iyear++)
         {if (nsamp_agec_cHL(iyear)>=minSS_agec)
            {neff_cHL_agec_allyr_out(yrs_agec_cHL(iyear))=multinom_eff_N(pred_cHL_agec(iyear),obs_agec_cHL(iyear));}                            
          else {neff_cHL_agec_allyr_out(yrs_agec_cHL(iyear))=-99;}
         }
      for (iyear=1; iyear<=nyr_agec_cPT; iyear++)
         {if (nsamp_agec_cPT(iyear)>=minSS_agec)
            {neff_cPT_agec_allyr_out(yrs_agec_cPT(iyear))=multinom_eff_N(pred_cPT_agec(iyear),obs_agec_cPT(iyear));}                            
          else {neff_cPT_agec_allyr_out(yrs_agec_cPT(iyear))=-99;}
         }
      for (iyear=1; iyear<=nyr_agec_rHB; iyear++)
         {if (nsamp_agec_rHB(iyear)>=minSS_agec)
            {neff_rHB_agec_allyr_out(yrs_agec_rHB(iyear))=multinom_eff_N(pred_rHB_agec(iyear),obs_agec_rHB(iyear));}                            
          else {neff_rHB_agec_allyr_out(yrs_agec_rHB(iyear))=-99;}
         }
      //for (iyear=1; iyear<=nyr_rGN_agec; iyear++)
      //   {if (nsamp_rGN_agec(iyear)>=minSS_agec)
      //      {neff_rGN_agec_allyr_out(yrs_rGN_agec(iyear))=multinom_eff_N(pred_rGN_agec(iyear),obs_rGN_agec(iyear));}                            
      //    else {neff_rGN_agec_allyr_out(yrs_rGN_agec(iyear))=-99;}
      //   }
}

void model_parameters::evaluate_objective_function(void)
{
  fval=0.0;
  fval_data=0.0;
  f_sBT_cpue=0.0;
  f_sBT_cpue=lk_lognormal(pred_sBT_cpue, obs_cpue_sBT, obs_cv_cpue_sBT, w_cpue_sBT);
  fval+=f_sBT_cpue;
  fval_data+=f_sBT_cpue;
  f_sTV_cpue=0.0;
  f_sTV_cpue=lk_lognormal(pred_sTV_cpue, obs_cpue_sTV, obs_cv_cpue_sTV, w_cpue_sTV);
  fval+=f_sTV_cpue;
  fval_data+=f_sTV_cpue;
  //f_Vid_cpue=0.0;
  //f_Vid_cpue=lk_lognormal(pred_Vid_cpue, obs_Vid_cpue, Vid_cpue_cv, w_I_Vid);
  //fval+=f_Vid_cpue;
  //fval_data+=f_Vid_cpue;
  f_cHL_cpue=0.0;
  f_cHL_cpue=lk_lognormal(pred_cHL_cpue, obs_cpue_cHL, obs_cv_cpue_cHL, w_cpue_cHL);
  fval+=f_cHL_cpue;
  fval_data+=f_cHL_cpue;
  f_rHB_cpue=0.0;
  f_rHB_cpue=lk_lognormal(pred_rHB_cpue, obs_cpue_rHB, obs_cv_cpue_rHB, w_cpue_rHB);
  fval+=f_rHB_cpue;
  fval_data+=f_rHB_cpue;
  //f_rHB_D_cpue=0.0;
  //f_rHB_D_cpue=lk_lognormal(pred_rHB_D_cpue, obs_rHB_D_cpue, rHB_D_cpue_cv, w_I_rHB_D);
  //fval+=f_rHB_D_cpue;
  //fval_data+=f_rHB_D_cpue;
  //f_cHL_L in 1000 lb ww
  f_cHL_L=lk_lognormal(pred_cHL_L_klb(styr_L_cHL,endyr_L_cHL), obs_L_cHL(styr_L_cHL,endyr_L_cHL), 
                      obs_cv_L_cHL(styr_L_cHL,endyr_L_cHL), w_L);
  fval+=f_cHL_L;
  fval_data+=f_cHL_L;  
  //f_cPT_L in 1000 lb ww
  f_cPT_L=lk_lognormal(pred_cPT_L_klb(styr_L_cPT,endyr_L_cPT), obs_L_cPT(styr_L_cPT,endyr_L_cPT), 
                      obs_cv_L_cPT(styr_L_cPT,endyr_L_cPT), w_L);
  fval+=f_cPT_L;
  fval_data+=f_cPT_L; 
  //f_cTW_L in 1000 lb ww
  f_cTW_L=lk_lognormal(pred_cTW_L_klb(styr_L_cTW,endyr_L_cTW), obs_L_cTW(styr_L_cTW,endyr_L_cTW), 
                      obs_cv_L_cTW(styr_L_cTW,endyr_L_cTW), w_L);
  fval+=f_cTW_L;
  fval_data+=f_cTW_L;
  //f_rHB_L in 1000 lb
  f_rHB_L=lk_lognormal(pred_rHB_L_klb(styr_L_rHB,endyr_L_rHB), obs_L_rHB(styr_L_rHB,endyr_L_rHB),
                      obs_cv_L_rHB(styr_L_rHB,endyr_L_rHB), w_L);
  fval+=f_rHB_L;
  fval_data+=f_rHB_L;
  //f_rGN_L in 1000 lb
  f_rGN_L=lk_lognormal(pred_rGN_L_klb(styr_L_rGN,endyr_L_rGN), obs_L_rGN(styr_L_rGN,endyr_L_rGN), 
                       obs_cv_L_rGN(styr_L_rGN,endyr_L_rGN), w_L);
  fval+=f_rGN_L;
  fval_data+=f_rGN_L;
  //f_cGN_D in 1000 fish
  f_cGN_D=lk_lognormal(pred_cGN_D_knum(styr_D_cHL,endyr_D_cHL), obs_cGN_D(styr_D_cHL,endyr_D_cHL), 
                      cGN_D_cv(styr_D_cHL,endyr_D_cHL), w_D);
  fval+=f_cGN_D;
  fval_data+=f_cGN_D;
  //f_rHB_D in 1000 fish
  f_rHB_D=lk_lognormal(pred_rHB_D_knum(styr_D_rHB,endyr_D_rHB), obs_rHB_D(styr_D_rHB,endyr_D_rHB), 
                         obs_cv_D_rHB(styr_D_rHB,endyr_D_rHB), w_D);
  fval+=f_rHB_D;
  fval_data+=f_rHB_D;
  //f_rGN_D in 1000 fish
  f_rGN_D=lk_lognormal(pred_rGN_D_knum(styr_D_rGN,endyr_D_rGN), obs_rGN_D(styr_D_rGN,endyr_D_rGN), 
                       obs_cv_D_rGN(styr_D_rGN,endyr_D_rGN), w_D);
  fval+=f_rGN_D;
  fval_data+=f_rGN_D;
  //f_sBT_lenc
  f_sBT_lenc=lk_dirichlet_multinomial(nsamp_lenc_sBT, pred_sBT_lenc, obs_lenc_sBT, nyr_lenc_sBT, double(nlenbins), minSS_lenc, log_dm_lenc_sBT); //w_lenc_sBT);
  fval+=f_sBT_lenc;
  fval_data+=f_sBT_lenc;
  //f_cHL_lenc
  f_cHL_lenc=lk_dirichlet_multinomial(nsamp_lenc_cHL, pred_cHL_lenc, obs_lenc_cHL, nyr_lenc_cHL, double(nlenbins), minSS_lenc, log_dm_lenc_cHL); //w_lenc_cHL);
  fval+=f_cHL_lenc;
  fval_data+=f_cHL_lenc;
  //f_cPT_lenc
  f_cPT_lenc=lk_dirichlet_multinomial(nsamp_lenc_cPT, pred_cPT_lenc, obs_lenc_cPT, nyr_lenc_cPT, double(nlenbins), minSS_lenc, log_dm_lenc_cPT); //w_lenc_cPT);
  fval+=f_cPT_lenc;
  fval_data+=f_cPT_lenc;
  //f_rHB_lenc
  f_rHB_lenc=lk_dirichlet_multinomial(nsamp_lenc_rHB, pred_rHB_lenc, obs_lenc_rHB, nyr_lenc_rHB, double(nlenbins), minSS_lenc, log_dm_lenc_rHB); //w_lenc_rHB);
  fval+=f_rHB_lenc;
  fval_data+=f_rHB_lenc;  
  //f_rHB_D_lenc
  f_rHB_D_lenc=lk_dirichlet_multinomial(nsamp_lenc_rHB_D, pred_rHB_D_lenc, obs_lenc_rHB_D, nyr_lenc_rHB_D, double(nlenbins), minSS_lenc, log_dm_lenc_rHB_D); //w_lenc_rHB_D);
  fval+=f_rHB_D_lenc;
  fval_data+=f_rHB_D_lenc;
  //f_rGN_lenc
  f_rGN_lenc=lk_dirichlet_multinomial(nsamp_lenc_rGN, pred_rGN_lenc, obs_lenc_rGN, nyr_lenc_rGN, double(nlenbins), minSS_lenc, log_dm_lenc_rGN); //w_lenc_rGN);
  fval+=f_rGN_lenc;
  fval_data+=f_rGN_lenc;
  //f_sBT_agec
  f_sBT_agec=lk_dirichlet_multinomial(nsamp_agec_sBT, pred_sBT_agec, obs_agec_sBT, nyr_agec_sBT, double(nages), minSS_agec, log_dm_agec_sBT); //w_agec_sBT);
  fval+=f_sBT_agec;
  fval_data+=f_sBT_agec;
  //f_sTV_agec
  f_sTV_agec=lk_dirichlet_multinomial(nsamp_agec_sTV, pred_sTV_agec, obs_agec_sTV, nyr_agec_sTV, double(nages), minSS_agec, log_dm_agec_sTV); //w_agec_sTV);
  fval+=f_sTV_agec;
  fval_data+=f_sTV_agec;
  //f_cHL_agec
  f_cHL_agec=lk_dirichlet_multinomial(nsamp_agec_cHL, pred_cHL_agec, obs_agec_cHL, nyr_agec_cHL, double(nages), minSS_agec, log_dm_agec_cHL); //w_agec_cHL);
  fval+=f_cHL_agec;
  fval_data+=f_cHL_agec;
  //f_cPT_agec
  f_cPT_agec=lk_dirichlet_multinomial(nsamp_agec_cPT, pred_cPT_agec, obs_agec_cPT, nyr_agec_cPT, double(nages), minSS_agec, log_dm_agec_cPT); //w_agec_cPT);
  fval+=f_cPT_agec;
  fval_data+=f_cPT_agec;
  //f_rHB_agec
  f_rHB_agec=lk_dirichlet_multinomial(nsamp_agec_rHB, pred_rHB_agec, obs_agec_rHB, nyr_agec_rHB, double(nages), minSS_agec, log_dm_agec_rHB); //w_agec_rHB);
  fval+=f_rHB_agec;
  fval_data+=f_rHB_agec;  
  ////f_rGN_agec
  //f_rGN_agec=lk_dirichlet_multinomial(nsamp_rGN_agec, pred_rGN_agec, obs_rGN_agec, nyr_rGN_agec, double(nages), minSS_agec, log_dm_rGN_ac); //w_ac_rGN);
  //fval+=f_rGN_agec;
  //fval_data+=f_rGN_agec;    
  ////f_rGN_agec=lk_multinomial(nsamp_rGN_agec, pred_rGN_agec, obs_rGN_agec, nyr_rGN_agec, minSS_agec, w_ac_rGN);
  ////fval+=f_rGN_agec;
  ////fval_data+=f_rGN_agec;  
  f_rec_dev=0.0;
  rec_logL_add=nyrs_rec*log(rec_sigma);
  f_rec_dev=(square(log_dev_rec(styr_rec_dev) + rec_sigma_sqd2)/(2.0*rec_sigma_sq));
  for(iyear=(styr_rec_dev+1); iyear<=endyr; iyear++)
  {f_rec_dev+=(square(log_dev_rec(iyear)-R_autocorr*log_dev_rec(iyear-1) + rec_sigma_sqd2)/
               (2.0*rec_sigma_sq));}
  f_rec_dev+=rec_logL_add;            
  fval+=w_rec*f_rec_dev;
  f_rec_dev_early=0.0; //possible extra constraint on early rec deviations
  if (w_rec_early>0.0)
    { if (styr_rec_dev<endyr_rec_phase1)
        {  
          f_rec_dev_early=(square(log_dev_rec(styr_rec_dev) + rec_sigma_sq/2.0)/(2.0*rec_sigma_sq)) + rec_logL_add;
          for(iyear=(styr_rec_dev+1); iyear<=endyr_rec_phase1; iyear++)
          {f_rec_dev_early+=(square(log_dev_rec(iyear)-R_autocorr*log_dev_rec(iyear-1) + rec_sigma_sqd2)/
                            (2.0*rec_sigma_sq)) + rec_logL_add;}
        }
  fval+=w_rec_early*f_rec_dev_early;
  }
  f_rec_dev_end=0.0; //possible extra constraint on ending rec deviations
  if (w_rec_end>0.0)
  { if (endyr_rec_phase2<endyr)
        {  
          for(iyear=(endyr_rec_phase2+1); iyear<=endyr; iyear++)
          {f_rec_dev_end+=(square(log_dev_rec(iyear)-R_autocorr*log_dev_rec(iyear-1) + rec_sigma_sqd2)/
                            (2.0*rec_sigma_sq)) + rec_logL_add;}
        }
      fval+=w_rec_end*f_rec_dev_end;
   }
  //applies if initial age structure is estimated
  fval+=norm2(log_dev_Nage); 
  //Random walk components of fishery dependent indices: these components equal zero if RW turned off
  f_cHL_RW_cpue=0.0;
  for (iyear=styr_cpue_cHL; iyear<endyr_cpue_cHL; iyear++)
      {f_cHL_RW_cpue+=square(q_RW_log_dev_cHL(iyear))/(2.0*set_q_RW_cHL_var);}
  fval+=f_cHL_RW_cpue;  
  f_rHB_RW_cpue=0.0;
  for (iyear=styr_cpue_rHB; iyear<endyr_cpue_rHB; iyear++)
      {f_rHB_RW_cpue+=square(q_RW_log_dev_rHB(iyear))/(2.0*set_q_RW_rHB_var);}      
  fval+=f_rHB_RW_cpue;   
  //f_rHB_D_RW_cpue=0.0;
  //for (iyear=styr_rHB_D_cpue; iyear<endyr_rHB_D_cpue; iyear++)
  //    {f_rHB_D_RW_cpue+=square(q_RW_log_dev_rHB_D(iyear))/(2.0*set_q_RW_rHB_D_var);}      
  //fval+=f_rHB_D_RW_cpue;  
  f_priors=0.0; 
  f_priors+=neg_log_prior(steep, set_steep(5), set_steep(6), set_steep(7)); 
  f_priors+=neg_log_prior(rec_sigma,set_rec_sigma(5),set_rec_sigma(6),set_rec_sigma(7));
  f_priors+=neg_log_prior(R_autocorr,set_R_autocorr(5),set_R_autocorr(6),set_R_autocorr(7));
  //f_priors+=neg_log_prior(q_DD_beta, set_q_DD_beta, square(set_q_DD_beta_se), 3);	
  //f_priors+=neg_log_prior(q_rate, set_q_rate, dzero+square(set_q_rate), 3);  
  f_priors+=neg_log_prior(F_init_ratio, set_F_init_ratio, -1.0 , 3);
  //f_priors+=neg_log_prior(M_constant, set_M_constant, square(set_M_constant_se), 3);	
  //f_priors+=neg_log_prior(Linf,set_Linf,square(set_Linf_se),3);
  //f_priors+=neg_log_prior(K,set_K,square(set_K_se),3);
  //f_priors+=neg_log_prior(t0,set_K,square(set_t0_se),3);
  f_priors+=neg_log_prior(len_cv_val,set_len_cv(5),set_len_cv(6),set_len_cv(7));
  f_priors+=neg_log_prior(selpar_A50_sBT, set_selpar_A50_sBT(5),set_selpar_A50_sBT(6),set_selpar_A50_sBT(7));  
  f_priors+=neg_log_prior(selpar_slope_sBT, set_selpar_slope_sBT(5),set_selpar_slope_sBT(6),set_selpar_slope_sBT(7));                  
  f_priors+=neg_log_prior(selpar_A50_sTV, set_selpar_A50_sTV(5),set_selpar_A50_sTV(6),set_selpar_A50_sTV(7));  
  f_priors+=neg_log_prior(selpar_slope_sTV, set_selpar_slope_sTV(5),set_selpar_slope_sTV(6),set_selpar_slope_sTV(7));                  
  f_priors+=neg_log_prior(selpar_A50_cHL2, set_selpar_A50_cHL2(5),set_selpar_A50_cHL2(6),set_selpar_A50_cHL2(7));  
  f_priors+=neg_log_prior(selpar_slope_cHL2, set_selpar_slope_cHL2(5),set_selpar_slope_cHL2(6),set_selpar_slope_cHL2(7));                  
  f_priors+=neg_log_prior(selpar_A50_cHL3, set_selpar_A50_cHL3(5),set_selpar_A50_cHL3(6),set_selpar_A50_cHL3(7));  
  f_priors+=neg_log_prior(selpar_slope_cHL3, set_selpar_slope_cHL3(5),set_selpar_slope_cHL3(6),set_selpar_slope_cHL3(7)); 
  f_priors+=neg_log_prior(selpar_A50_cHL4, set_selpar_A50_cHL4(5),set_selpar_A50_cHL4(6),set_selpar_A50_cHL4(7));  
  f_priors+=neg_log_prior(selpar_slope_cHL4, set_selpar_slope_cHL4(5),set_selpar_slope_cHL4(6),set_selpar_slope_cHL4(7));   
  f_priors+=neg_log_prior(selpar_A50_cPT2, set_selpar_A50_cPT2(5),set_selpar_A50_cPT2(6),set_selpar_A50_cPT2(7));  
  f_priors+=neg_log_prior(selpar_slope_cPT2, set_selpar_slope_cPT2(5),set_selpar_slope_cPT2(6),set_selpar_slope_cPT2(7));                  
  f_priors+=neg_log_prior(selpar_A50_cPT3, set_selpar_A50_cPT3(5),set_selpar_A50_cPT3(6),set_selpar_A50_cPT3(7));  
  f_priors+=neg_log_prior(selpar_slope_cPT3, set_selpar_slope_cPT3(5),set_selpar_slope_cPT3(6),set_selpar_slope_cPT3(7)); 
  f_priors+=neg_log_prior(selpar_A50_cPT4, set_selpar_A50_cPT4(5),set_selpar_A50_cPT4(6),set_selpar_A50_cPT4(7));  
  f_priors+=neg_log_prior(selpar_slope_cPT4, set_selpar_slope_cPT4(5),set_selpar_slope_cPT4(6),set_selpar_slope_cPT4(7));   
  f_priors+=neg_log_prior(selpar_A50_rHB1, set_selpar_A50_rHB1(5),set_selpar_A50_rHB1(6),set_selpar_A50_rHB1(7));  
  f_priors+=neg_log_prior(selpar_slope_rHB1, set_selpar_slope_rHB1(5),set_selpar_slope_rHB1(6),set_selpar_slope_rHB1(7));                  
  f_priors+=neg_log_prior(selpar_A50_rHB2, set_selpar_A50_rHB2(5),set_selpar_A50_rHB2(6),set_selpar_A50_rHB2(7));  
  f_priors+=neg_log_prior(selpar_slope_rHB2, set_selpar_slope_rHB2(5),set_selpar_slope_rHB2(6),set_selpar_slope_rHB2(7));                  
  f_priors+=neg_log_prior(selpar_A50_rHB3, set_selpar_A50_rHB3(5),set_selpar_A50_rHB3(6),set_selpar_A50_rHB3(7));  
  f_priors+=neg_log_prior(selpar_slope_rHB3, set_selpar_slope_rHB3(5),set_selpar_slope_rHB3(6),set_selpar_slope_rHB3(7));      
  f_priors+=neg_log_prior(selpar_A50_rHB4, set_selpar_A50_rHB4(5),set_selpar_A50_rHB4(6),set_selpar_A50_rHB4(7)); 
  f_priors+=neg_log_prior(selpar_slope_rHB4, set_selpar_slope_rHB4(5),set_selpar_slope_rHB4(6),set_selpar_slope_rHB4(7)); 
  f_priors+=neg_log_prior(selpar_A50_rHB5, set_selpar_A50_rHB5(5),set_selpar_A50_rHB5(6),set_selpar_A50_rHB5(7));   
  f_priors+=neg_log_prior(selpar_slope_rHB5, set_selpar_slope_rHB5(5),set_selpar_slope_rHB5(6),set_selpar_slope_rHB5(7));   
  f_priors+=neg_log_prior(selpar_A50_rGN1, set_selpar_A50_rGN1(5),set_selpar_A50_rGN1(6),set_selpar_A50_rGN1(7));  
  f_priors+=neg_log_prior(selpar_slope_rGN1, set_selpar_slope_rGN1(5),set_selpar_slope_rGN1(6),set_selpar_slope_rGN1(7));                  
  f_priors+=neg_log_prior(selpar_A50_rGN2, set_selpar_A50_rGN2(5),set_selpar_A50_rGN2(6),set_selpar_A50_rGN2(7));  
  f_priors+=neg_log_prior(selpar_slope_rGN2, set_selpar_slope_rGN2(5),set_selpar_slope_rGN2(6),set_selpar_slope_rGN2(7));                  
  f_priors+=neg_log_prior(selpar_A50_rGN3, set_selpar_A50_rGN3(5),set_selpar_A50_rGN3(6),set_selpar_A50_rGN3(7));  
  f_priors+=neg_log_prior(selpar_slope_rGN3, set_selpar_slope_rGN3(5),set_selpar_slope_rGN3(6),set_selpar_slope_rGN3(7));      
  f_priors+=neg_log_prior(selpar_A50_rGN4, set_selpar_A50_rGN4(5),set_selpar_A50_rGN4(6),set_selpar_A50_rGN4(7)); 
  f_priors+=neg_log_prior(selpar_slope_rGN4, set_selpar_slope_rGN4(5),set_selpar_slope_rGN4(6),set_selpar_slope_rGN4(7)); 
  f_priors+=neg_log_prior(selpar_A50_rGN5, set_selpar_A50_rGN5(5),set_selpar_A50_rGN5(6),set_selpar_A50_rGN5(7));   
  f_priors+=neg_log_prior(selpar_slope_rGN5, set_selpar_slope_rGN5(5),set_selpar_slope_rGN5(6),set_selpar_slope_rGN5(7));  
  f_priors+=neg_log_prior(selpar_logit_Age0_rHB_D, set_selpar_logit_Age0_rHB_D(5),set_selpar_logit_Age0_rHB_D(6),set_selpar_logit_Age0_rHB_D(7));
  f_priors+=neg_log_prior(selpar_logit_Age1_rHB_D, set_selpar_logit_Age1_rHB_D(5),set_selpar_logit_Age1_rHB_D(6),set_selpar_logit_Age1_rHB_D(7));
  f_priors+=neg_log_prior(selpar_logit_Age2_rHB_D, set_selpar_logit_Age2_rHB_D(5),set_selpar_logit_Age2_rHB_D(6), set_selpar_logit_Age2_rHB_D(7));
  f_priors+=neg_log_prior(selpar_A50_rHB_D4, set_selpar_A50_rHB_D4(5),set_selpar_A50_rHB_D4(6),set_selpar_A50_rHB_D4(7)); 
  f_priors+=neg_log_prior(selpar_slope_rHB_D4, set_selpar_slope_rHB_D4(5),set_selpar_slope_rHB_D4(6),set_selpar_slope_rHB_D4(7)); 
  f_priors+=neg_log_prior(selpar_A502_rHB_D4, set_selpar_A502_rHB_D4(5),set_selpar_A502_rHB_D4(6),set_selpar_A502_rHB_D4(7)); f_priors+=neg_log_prior(selpar_A50_rHB_D5, set_selpar_A50_rHB_D5(5),set_selpar_A50_rHB_D5(6),set_selpar_A50_rHB_D5(7));   
  f_priors+=neg_log_prior(selpar_slope_rHB_D5, set_selpar_slope_rHB_D5(5),set_selpar_slope_rHB_D5(6),set_selpar_slope_rHB_D5(7)); 
  f_priors+=neg_log_prior(selpar_A502_rHB_D5, set_selpar_A502_rHB_D5(5),set_selpar_A502_rHB_D5(6),set_selpar_A502_rHB_D5(7));
  fval+=f_priors;
  //cout << "fval = " << fval << "  fval_data = " << fval_data << endl;
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

dvar_vector model_parameters::logistic_double(const dvar_vector& ages, const dvariable& A501, const dvariable& slope1, const dvariable& A502, const dvariable& slope2)
{
  //ages=vector of ages, A50=age at 50% selectivity, slope=rate of increase, A502=age at 50% decrease additive to A501, slope2=slope of decrease
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=elem_prod( (1./(1.+mfexp(-1.*slope1*(ages-A501)))),(1.-(1./(1.+mfexp(-1.*slope2*(ages-(A501+A502)))))) );     
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;
  //Logistic-exponential: 4 parameters (but 1 is fixed)
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

dvariable model_parameters::SR_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& SSB)
{
  //R0=virgin recruitment, h=steepness, spr_F0=spawners per recruit @ F=0, SSB=spawning biomass
  RETURN_ARRAYS_INCREMENT();
  dvariable Recruits_Tmp;
  Recruits_Tmp=((0.8*R0*h*SSB)/(0.2*R0*spr_F0*(1.0-h)+(h-0.2)*SSB));
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
  //dzero is small value to avoid log(0) during search
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvar_vector var(cv.indexmin(),cv.indexmax()); //variance in log space
  var=log(1.0+square(cv/wgt_dat));   // convert cv in arithmetic space to variance in log space
  LkvalTmp=sum(0.5*elem_div(square(log(elem_div((pred+dzero),(obs+dzero)))),var) );
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;
}

dvariable model_parameters::lk_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const double& minSS, const dvariable& wgt_dat)
{
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, minSS=min N threshold, wgt_dat=scaling of N's
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  LkvalTmp=0.0;
  for (int ii=1; ii<=ncomp; ii++)
  {if (nsamp(ii)>=minSS)
    {LkvalTmp-=wgt_dat*nsamp(ii)*sum(elem_prod((obs_comp(ii)+dzero),
               log(elem_div((pred_comp(ii)+dzero), (obs_comp(ii)+dzero)))));
    }
  }  
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;
  //-----------------------------------------------------------------------------------
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

dvariable model_parameters::neg_log_prior(dvariable pred, const double& prior, dvariable var, int pdf)
{
  //prior=prior point estimate, var=variance (if negative, treated as CV in arithmetic space), pred=predicted value, pdf=prior type (1=none, 2=lognormal, 3=normal, 4=beta)
    dvariable LkvalTmp;
    dvariable alpha, beta, ab_iq;
    LkvalTmp=0.0;
    // compute generic pdf's
    switch(pdf) {
        case 1: //option to turn off prior
          LkvalTmp=0.0;
          break;
        case 2: // lognormal 
          if(prior<=0.0) cout << "YIKES: Don't use a lognormal distn for a negative prior" << endl;
          else if(pred<=0) LkvalTmp=huge_number;
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
          else LkvalTmp=huge_number;
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
      get_miscellaneous_stuff();
      cout<<"got misc stuff"<<endl;
      get_per_recruit_stuff();
      cout<<"got per recruit"<<endl;  
      get_effective_sample_sizes();
	  cout<<"got effective sample sizes"<<endl; 
      //get_projection();
	  //cout<<"got projection"<<endl;
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
      cout << "BC Fmsy=" << F_msy_out<< "   BC SSBmsy=" << SSB_msy_out <<endl;
      cout <<"F status="<<FdF_msy_end<<endl;
      cout <<"Pop status="<<SdSSB_msy_end<<endl;
      cout << "h="<<steep<<"   R0="<<R0<<endl;
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;  
      report << "TotalLikelihood " << fval << endl;
      report << "N" << endl;
      report << N<<endl;      
      report << "SSB" << endl;       
      report << SSB << endl;
	  //sdnr_lc_sBT=sdnr_multinomial(nyr_lenc_sBT, lenbins, nsamp_lenc_sBT, pred_sBT_lenc, obs_lenc_sBT, w_lenc_sBT); 
      //sdnr_lc_sTV=sdnr_multinomial(nyr_sTV_lenc, lenbins, nsamp_sTV_lenc, pred_sTV_lenc, obs_sTV_lenc, w_lc_sTV); 
      //sdnr_lc_cHL=sdnr_multinomial(nyr_lenc_cHL, lenbins, nsamp_lenc_cHL, pred_cHL_lenc, obs_lenc_cHL, w_lenc_cHL);  
      //sdnr_lc_cPT=sdnr_multinomial(nyr_lenc_cPT, lenbins, nsamp_lenc_cPT, pred_cPT_lenc, obs_lenc_cPT, w_lenc_cPT);  	  
      //sdnr_lc_rHB=sdnr_multinomial(nyr_lenc_rHB, lenbins, nsamp_lenc_rHB, pred_rHB_lenc, obs_lenc_rHB, w_lenc_rHB); 
	  //sdnr_lc_rHB_D=sdnr_multinomial(nyr_lenc_rHB_D, lenbins, nsamp_lenc_rHB_D, pred_rHB_D_lenc, obs_lenc_rHB_D, w_lenc_rHB_D); 
	  //sdnr_lc_rGN=sdnr_multinomial(nyr_lenc_rGN, lenbins, nsamp_lenc_rGN, pred_rGN_lenc, obs_lenc_rGN, w_lenc_rGN); 
      //sdnr_ac_sBT=sdnr_multinomial(nyr_agec_sBT, agebins_agec, nsamp_agec_sBT, pred_sBT_agec, obs_agec_sBT, w_agec_sBT);  
      //sdnr_ac_sTV=sdnr_multinomial(nyr_agec_sTV, agebins_agec, nsamp_agec_sTV, pred_sTV_agec, obs_agec_sTV, w_agec_sTV);  
      //sdnr_ac_cHL=sdnr_multinomial(nyr_agec_cHL, agebins_agec, nsamp_agec_cHL, pred_cHL_agec, obs_agec_cHL, w_agec_cHL);  
      //sdnr_ac_cPT=sdnr_multinomial(nyr_agec_cPT, agebins_agec, nsamp_agec_cPT, pred_cPT_agec, obs_agec_cPT, w_agec_cPT);  
      //sdnr_ac_rHB=sdnr_multinomial(nyr_agec_rHB, agebins_agec, nsamp_agec_rHB, pred_rHB_agec, obs_agec_rHB, w_agec_rHB);  
	  //sdnr_ac_rGN=sdnr_multinomial(nyr_rGN_agec, agebins_agec, nsamp_rGN_agec, pred_rGN_agec, obs_rGN_agec, w_ac_rGN);  
      sdnr_I_sBT=sdnr_lognormal(pred_sBT_cpue, obs_cpue_sBT, obs_cv_cpue_sBT, w_cpue_sBT);  
      sdnr_I_sTV=sdnr_lognormal(pred_sTV_cpue, obs_cpue_sTV, obs_cv_cpue_sTV, w_cpue_sTV);  
      sdnr_I_cHL=sdnr_lognormal(pred_cHL_cpue, obs_cpue_cHL, obs_cv_cpue_cHL, w_cpue_cHL);
      sdnr_I_rHB=sdnr_lognormal(pred_rHB_cpue, obs_cpue_rHB, obs_cv_cpue_rHB, w_cpue_rHB);
      //sdnr_I_Vid=sdnr_lognormal(pred_Vid_cpue, obs_Vid_cpue, Vid_cpue_cv, w_I_Vid);
	   Linf_out(8)=Linf; Linf_out(1,7)=set_Linf; 
       K_out(8)=K; K_out(1,7)=set_K;
       t0_out(8)=t0; t0_out(1,7)=set_t0;
       len_cv_val_out(8)=len_cv_val; len_cv_val_out(1,7)=set_len_cv;
       log_R0_out(8)=log_R0; log_R0_out(1,7)=set_log_R0;
       M_constant_out(8)=M_constant; M_constant_out(1,7)=set_M_constant;
       steep_out(8)=steep; steep_out(1,7)=set_steep;
       rec_sigma_out(8)=rec_sigma; rec_sigma_out(1,7)=set_rec_sigma;
       R_autocorr_out(8)=R_autocorr; R_autocorr_out(1,7)=set_R_autocorr;
	   log_dm_sBT_lc_out(8)=log_dm_lenc_sBT; log_dm_sBT_lc_out(1,7)=set_log_dm_lenc_sBT;
	   log_dm_sTV_lc_out(8)=log_dm_lenc_sTV; log_dm_sTV_lc_out(1,7)=set_log_dm_lenc_sTV;
	   log_dm_cHL_lc_out(8)=log_dm_lenc_cHL; log_dm_cHL_lc_out(1,7)=set_log_dm_lenc_cHL;
	   log_dm_cPT_lc_out(8)=log_dm_lenc_cPT; log_dm_cPT_lc_out(1,7)=set_log_dm_lenc_cPT;
	   log_dm_rHB_lc_out(8)=log_dm_lenc_rHB; log_dm_rHB_lc_out(1,7)=set_log_dm_lenc_rHB;
	   log_dm_rHB_D_lc_out(8)=log_dm_lenc_rHB_D; log_dm_rHB_D_lc_out(1,7)=set_log_dm_lenc_rHB_D;
	   log_dm_rGN_lc_out(8)=log_dm_lenc_rGN; log_dm_rGN_lc_out(1,7)=set_log_dm_lenc_rGN;
	   log_dm_sBT_ac_out(8)=log_dm_agec_sBT; log_dm_sBT_ac_out(1,7)=set_log_dm_agec_sBT;
       log_dm_sTV_ac_out(8)=log_dm_agec_sTV; log_dm_sTV_ac_out(1,7)=set_log_dm_agec_sTV;
       log_dm_cHL_ac_out(8)=log_dm_agec_cHL; log_dm_cHL_ac_out(1,7)=set_log_dm_agec_cHL;
	   log_dm_cPT_ac_out(8)=log_dm_agec_cPT; log_dm_cPT_ac_out(1,7)=set_log_dm_agec_cPT;
	   log_dm_rHB_ac_out(8)=log_dm_agec_rHB; log_dm_rHB_ac_out(1,7)=set_log_dm_agec_rHB;
	   //log_dm_rGN_ac_out(8)=log_dm_rGN_ac; log_dm_rGN_ac_out(1,7)=set_log_dm_rGN_ac;
       selpar_A50_sBT_out(8)=selpar_A50_sBT; selpar_A50_sBT_out(1,7)=set_selpar_A50_sBT;
       selpar_slope_sBT_out(8)=selpar_slope_sBT; selpar_slope_sBT_out(1,7)=set_selpar_slope_sBT;
	   selpar_A50_sTV_out(8)=selpar_A50_sTV; selpar_A50_sTV_out(1,7)=set_selpar_A50_sTV;
       selpar_slope_sTV_out(8)=selpar_slope_sTV; selpar_slope_sTV_out(1,7)=set_selpar_slope_sTV;
	   //selpar_A50_cHL1_out(8)=selpar_A50_cHL1; selpar_A50_cHL1_out(1,7)=set_selpar_A50_cHL1;
       //selpar_slope_cHL1_out(8)=selpar_slope_cHL1; selpar_slope_cHL1_out(1,7)=set_selpar_slope_cHL1;
       selpar_A50_cHL2_out(8)=selpar_A50_cHL2; selpar_A50_cHL2_out(1,7)=set_selpar_A50_cHL2;
       selpar_slope_cHL2_out(8)=selpar_slope_cHL2; selpar_slope_cHL2_out(1,7)=set_selpar_slope_cHL2;
       selpar_A50_cHL3_out(8)=selpar_A50_cHL3; selpar_A50_cHL3_out(1,7)=set_selpar_A50_cHL3;
       selpar_slope_cHL3_out(8)=selpar_slope_cHL3; selpar_slope_cHL3_out(1,7)=set_selpar_slope_cHL3;
	   selpar_A50_cHL4_out(8)=selpar_A50_cHL4; selpar_A50_cHL4_out(1,7)=set_selpar_A50_cHL4;
       selpar_slope_cHL4_out(8)=selpar_slope_cHL4; selpar_slope_cHL4_out(1,7)=set_selpar_slope_cHL4;
	   selpar_A50_cPT2_out(8)=selpar_A50_cPT2; selpar_A50_cPT2_out(1,7)=set_selpar_A50_cPT2;
       selpar_slope_cPT2_out(8)=selpar_slope_cPT2; selpar_slope_cPT2_out(1,7)=set_selpar_slope_cPT2;
       selpar_A50_cPT3_out(8)=selpar_A50_cPT3; selpar_A50_cPT3_out(1,7)=set_selpar_A50_cPT3;
       selpar_slope_cPT3_out(8)=selpar_slope_cPT3; selpar_slope_cPT3_out(1,7)=set_selpar_slope_cPT3;
	   selpar_A50_cPT4_out(8)=selpar_A50_cPT4; selpar_A50_cPT4_out(1,7)=set_selpar_A50_cPT4;
       selpar_slope_cPT4_out(8)=selpar_slope_cPT4; selpar_slope_cPT4_out(1,7)=set_selpar_slope_cPT4;
       selpar_A50_rHB1_out(8)=selpar_A50_rHB1; selpar_A50_rHB1_out(1,7)=set_selpar_A50_rHB1;
       selpar_slope_rHB1_out(8)=selpar_slope_rHB1; selpar_slope_rHB1_out(1,7)=set_selpar_slope_rHB1;
       selpar_A50_rHB2_out(8)=selpar_A50_rHB2; selpar_A50_rHB2_out(1,7)=set_selpar_A50_rHB2;
       selpar_slope_rHB2_out(8)=selpar_slope_rHB2; selpar_slope_rHB2_out(1,7)=set_selpar_slope_rHB2;
       selpar_A50_rHB3_out(8)=selpar_A50_rHB3; selpar_A50_rHB3_out(1,7)=set_selpar_A50_rHB3;
       selpar_slope_rHB3_out(8)=selpar_slope_rHB3; selpar_slope_rHB3_out(1,7)=set_selpar_slope_rHB3;
	   selpar_A50_rHB4_out(8)=selpar_A50_rHB4; selpar_A50_rHB4_out(1,7)=set_selpar_A50_rHB4;
       selpar_slope_rHB4_out(8)=selpar_slope_rHB4; selpar_slope_rHB4_out(1,7)=set_selpar_slope_rHB4;
	   selpar_A50_rHB5_out(8)=selpar_A50_rHB5; selpar_A50_rHB5_out(1,7)=set_selpar_A50_rHB5;
       selpar_slope_rHB5_out(8)=selpar_slope_rHB5; selpar_slope_rHB5_out(1,7)=set_selpar_slope_rHB5;
	   selpar_A50_rGN1_out(8)=selpar_A50_rGN1; selpar_A50_rGN1_out(1,7)=set_selpar_A50_rGN1;
       selpar_slope_rGN1_out(8)=selpar_slope_rGN1; selpar_slope_rGN1_out(1,7)=set_selpar_slope_rGN1;
       selpar_A50_rGN2_out(8)=selpar_A50_rGN2; selpar_A50_rGN2_out(1,7)=set_selpar_A50_rGN2;
       selpar_slope_rGN2_out(8)=selpar_slope_rGN2; selpar_slope_rGN2_out(1,7)=set_selpar_slope_rGN2;
       selpar_A50_rGN3_out(8)=selpar_A50_rGN3; selpar_A50_rGN3_out(1,7)=set_selpar_A50_rGN3;
       selpar_slope_rGN3_out(8)=selpar_slope_rGN3; selpar_slope_rGN3_out(1,7)=set_selpar_slope_rGN3;
	   selpar_A50_rGN4_out(8)=selpar_A50_rGN4; selpar_A50_rGN4_out(1,7)=set_selpar_A50_rGN4;
       selpar_slope_rGN4_out(8)=selpar_slope_rGN4; selpar_slope_rGN4_out(1,7)=set_selpar_slope_rGN4;
	   selpar_A50_rGN5_out(8)=selpar_A50_rGN5; selpar_A50_rGN5_out(1,7)=set_selpar_A50_rGN5;
       selpar_slope_rGN5_out(8)=selpar_slope_rGN5; selpar_slope_rGN5_out(1,7)=set_selpar_slope_rGN5;
	   selpar_Age0_rHB_D_logit_out(8)=selpar_logit_Age0_rHB_D; selpar_Age0_rHB_D_logit_out(1,7)=set_selpar_logit_Age0_rHB_D;
       selpar_Age1_rHB_D_logit_out(8)=selpar_logit_Age1_rHB_D; selpar_Age1_rHB_D_logit_out(1,7)=set_selpar_logit_Age1_rHB_D;
	   selpar_Age2_rHB_D_logit_out(8)=selpar_logit_Age2_rHB_D; selpar_Age2_rHB_D_logit_out(1,7)=set_selpar_logit_Age2_rHB_D;
	   selpar_A50_rHB_D4_out(8)=selpar_A50_rHB_D4; selpar_A50_rHB_D4_out(1,7)=set_selpar_A50_rHB_D4;
       selpar_slope_rHB_D4_out(8)=selpar_slope_rHB_D4; selpar_slope_rHB_D4_out(1,7)=set_selpar_slope_rHB_D4;
	   selpar_A502_rHB_D4_out(8)=selpar_A502_rHB_D4; selpar_A502_rHB_D4_out(1,7)=set_selpar_A502_rHB_D4;
       selpar_slope2_rHB_D4_out(8)=selpar_slope2_rHB_D4; selpar_slope2_rHB_D4_out(1,7)=set_selpar_slope2_rHB_D4;
	   selpar_A50_rHB_D5_out(8)=selpar_A50_rHB_D5; selpar_A50_rHB_D5_out(1,7)=set_selpar_A50_rHB_D5;
       selpar_slope_rHB_D5_out(8)=selpar_slope_rHB_D5; selpar_slope_rHB_D5_out(1,7)=set_selpar_slope_rHB_D5;
	   selpar_A502_rHB_D5_out(8)=selpar_A502_rHB_D5; selpar_A502_rHB_D5_out(1,7)=set_selpar_A502_rHB_D5;
       selpar_slope2_rHB_D5_out(8)=selpar_slope2_rHB_D5; selpar_slope2_rHB_D5_out(1,7)=set_selpar_slope2_rHB_D5;
       log_q_sBT_out(8)=log_q_sBT; log_q_sBT_out(1,7)=set_log_q_cpue_sBT;
	   log_q_sTV_out(8)=log_q_sTV; log_q_sTV_out(1,7)=set_log_q_cpue_sTV;
	   //log_q_Vid_out(8)=log_q_Vid; log_q_Vid_out(1,7)=set_logq_Vid;
       log_q_cHL_out(8)=log_q_cHL; log_q_cHL_out(1,7)=set_log_q_cpue_cHL;
	   log_q_rHB_out(8)=log_q_rHB; log_q_rHB_out(1,7)=set_log_q_cpue_rHB;
       log_avg_F_cHL_out(8)=log_avg_F_L_cHL; log_avg_F_cHL_out(1,7)=set_log_avg_F_L_cHL;
	   log_avg_F_cPT_out(8)=log_avg_F_L_cPT; log_avg_F_cPT_out(1,7)=set_log_avg_F_L_cPT;
       log_avg_F_rHB_out(8)=log_avg_F_L_rHB; log_avg_F_rHB_out(1,7)=set_log_avg_F_L_rHB;
       log_avg_F_cTW_out(8)=log_avg_F_L_cTW; log_avg_F_cTW_out(1,7)=set_log_avg_F_L_cTW;       
       log_avg_F_rGN_out(8)=log_avg_F_L_rGN; log_avg_F_rGN_out(1,7)=set_log_avg_F_L_rGN;
       log_avg_F_cGN_D_out(8)=log_avg_F_D_cGN; log_avg_F_cGN_D_out(1,7)=set_log_avg_F_D_cGN;
       log_avg_F_rHB_D_out(8)=log_avg_F_D_rHB; log_avg_F_rHB_D_out(1,7)=set_log_avg_F_D_rHB;
       log_avg_F_rGN_D_out(8)=log_avg_F_D_rGN; log_avg_F_rGN_D_out(1,7)=set_log_avg_F_D_rGN;
       log_rec_dev_out(styr_rec_dev, endyr_rec_phase2)=log_dev_rec;
       log_F_dev_cHL_out(styr_L_cHL,endyr_L_cHL)=log_dev_F_L_cHL;
	   log_F_dev_cPT_out(styr_L_cPT,endyr_L_cPT)=log_dev_F_L_cPT;
       log_F_dev_cTW_out(styr_L_cTW,endyr_L_cTW)=log_dev_F_L_cTW;
       log_F_dev_rHB_out(styr_L_rHB,endyr_L_rHB)=log_dev_F_L_rHB;
       log_F_dev_rGN_out(styr_L_rGN,endyr_L_rGN)=log_dev_F_L_rGN;
       log_F_dev_cGN_D_out(styr_D_cHL,endyr_D_cHL)=log_dev_F_D_cGN;
       log_F_dev_rHB_D_out(styr_D_rHB,endyr_D_rHB)=log_dev_F_D_rHB;
       log_F_dev_rGN_D_out(styr_D_rGN,endyr_D_rGN)=log_dev_F_D_rGN;
     #include "bam.cxx"   // write the S-compatible report
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
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(500);
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
