#ifdef DEBUG
  #ifndef __SUNPRO_C
    #include <cfenv>
    #include <cstdlib>
  #endif
#endif
  #include "admodel.h"          // Include AD class definitions
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
  endyr_selex_phase2.allocate("endyr_selex_phase2");
  endyr_selex_phase3.allocate("endyr_selex_phase3");
   nyrs=endyr-styr+1.;
   nyrs_rec=endyr_rec_dev-styr_rec_dev+1.;
  nages.allocate("nages");
  agebins.allocate(1,nages,"agebins");
  nages_agec.allocate("nages_agec");
  nages_agec_10p.allocate("nages_agec_10p");
  agebins_agec.allocate(1,nages_agec,"agebins_agec");
  agebins_agec_10p.allocate(1,nages_agec_10p,"agebins_agec_10p");
  nlenbins.allocate("nlenbins");
  lenbins_width.allocate("lenbins_width");
  lenbins.allocate(1,nlenbins,"lenbins");
  max_F_spr_msy.allocate("max_F_spr_msy");
  n_iter_spr.allocate("n_iter_spr");
		n_iter_msy=n_iter_spr; 
  styr_rec_spr.allocate("styr_rec_spr");
  endyr_rec_spr.allocate("endyr_rec_spr");
   nyrs_rec_spr=endyr_rec_spr-styr_rec_spr+1.;
  SPR_rec_switch.allocate("SPR_rec_switch");
  SR_switch.allocate("SR_switch");
  selpar_n_yrs_wgted.allocate("selpar_n_yrs_wgted");
  set_BiasCor.allocate("set_BiasCor");
  E_age_min.allocate("E_age_min");
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
  styr_D_cHL.allocate("styr_D_cHL");
  endyr_D_cHL.allocate("endyr_D_cHL");
  obs_released_cHL.allocate(styr_D_cHL,endyr_D_cHL,"obs_released_cHL");
  obs_cv_D_cHL.allocate(styr_D_cHL,endyr_D_cHL,"obs_cv_D_cHL");
  nyr_cHL_D_lenc_pool1.allocate("nyr_cHL_D_lenc_pool1");
  yrs_cHL_D_lenc_pool1.allocate(1,nyr_cHL_D_lenc_pool1,"yrs_cHL_D_lenc_pool1");
  nsamp_cHL_D_lenc_pool1.allocate(1,nyr_cHL_D_lenc_pool1,"nsamp_cHL_D_lenc_pool1");
  nyr_cHL_D_lenc_pool2.allocate("nyr_cHL_D_lenc_pool2");
  yrs_cHL_D_lenc_pool2.allocate(1,nyr_cHL_D_lenc_pool2,"yrs_cHL_D_lenc_pool2");
  nsamp_cHL_D_lenc_pool2.allocate(1,nyr_cHL_D_lenc_pool2,"nsamp_cHL_D_lenc_pool2");
  nyr_lenc_cHL_D.allocate("nyr_lenc_cHL_D");
  yrs_lenc_cHL_D.allocate(1,nyr_lenc_cHL_D,"yrs_lenc_cHL_D");
  nsamp_lenc_cHL_D.allocate(1,nyr_lenc_cHL_D,"nsamp_lenc_cHL_D");
  nfish_lenc_cHL_D.allocate(1,nyr_lenc_cHL_D,"nfish_lenc_cHL_D");
  obs_lenc_cHL_D.allocate(1,nyr_lenc_cHL_D,1,nlenbins,"obs_lenc_cHL_D");
  styr_cpue_rHB.allocate("styr_cpue_rHB");
  endyr_cpue_rHB.allocate("endyr_cpue_rHB");
  obs_cpue_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"obs_cpue_rHB");
  obs_cv_cpue_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"obs_cv_cpue_rHB");
  styr_L_rHB.allocate("styr_L_rHB");
  endyr_L_rHB.allocate("endyr_L_rHB");
  obs_L_rHB.allocate(styr_L_rHB,endyr_L_rHB,"obs_L_rHB");
  obs_cv_L_rHB.allocate(styr_L_rHB,endyr_L_rHB,"obs_cv_L_rHB");
  nyr_agec_rHB.allocate("nyr_agec_rHB");
  yrs_agec_rHB.allocate(1,nyr_agec_rHB,"yrs_agec_rHB");
  nsamp_agec_rHB.allocate(1,nyr_agec_rHB,"nsamp_agec_rHB");
  nfish_agec_rHB.allocate(1,nyr_agec_rHB,"nfish_agec_rHB");
  obs_agec_rHB.allocate(1,nyr_agec_rHB,1,nages_agec_10p,"obs_agec_rHB");
  styr_D_cpue_rHB.allocate("styr_D_cpue_rHB");
  endyr_D_cpue_rHB.allocate("endyr_D_cpue_rHB");
  obs_cpue_rHB_D.allocate(styr_D_cpue_rHB,endyr_D_cpue_rHB,"obs_cpue_rHB_D");
  obs_cv_cpue_rHB_D.allocate(styr_D_cpue_rHB,endyr_D_cpue_rHB,"obs_cv_cpue_rHB_D");
  styr_D_rHB.allocate("styr_D_rHB");
  endyr_D_rHB.allocate("endyr_D_rHB");
  obs_released_rHB.allocate(styr_D_rHB,endyr_D_rHB,"obs_released_rHB");
  obs_cv_D_rHB.allocate(styr_D_rHB,endyr_D_rHB,"obs_cv_D_rHB");
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
  nyr_agec_rGN.allocate("nyr_agec_rGN");
  yrs_agec_rGN.allocate(1,nyr_agec_rGN,"yrs_agec_rGN");
  nsamp_agec_rGN.allocate(1,nyr_agec_rGN,"nsamp_agec_rGN");
  nfish_agec_rGN.allocate(1,nyr_agec_rGN,"nfish_agec_rGN");
  obs_agec_rGN.allocate(1,nyr_agec_rGN,1,nages_agec_10p,"obs_agec_rGN");
  nyr_lenc_rGN_D.allocate("nyr_lenc_rGN_D");
  yrs_lenc_rGN_D.allocate(1,nyr_lenc_rGN_D,"yrs_lenc_rGN_D");
  nsamp_lenc_rGN_D.allocate(1,nyr_lenc_rGN_D,"nsamp_lenc_rGN_D");
  nfish_lenc_rGN_D.allocate(1,nyr_lenc_rGN_D,"nfish_lenc_rGN_D");
  obs_lenc_rGN_D.allocate(1,nyr_lenc_rGN_D,1,nlenbins,"obs_lenc_rGN_D");
  styr_cpue_sCT.allocate("styr_cpue_sCT");
  endyr_cpue_sCT.allocate("endyr_cpue_sCT");
  obs_cpue_sCT.allocate(styr_cpue_sCT,endyr_cpue_sCT,"obs_cpue_sCT");
  obs_cv_cpue_sCT.allocate(styr_cpue_sCT,endyr_cpue_sCT,"obs_cv_cpue_sCT");
  nyr_agec_sCT.allocate("nyr_agec_sCT");
  yrs_agec_sCT.allocate(1,nyr_agec_sCT,"yrs_agec_sCT");
  nsamp_agec_sCT.allocate(1,nyr_agec_sCT,"nsamp_agec_sCT");
  nfish_agec_sCT.allocate(1,nyr_agec_sCT,"nfish_agec_sCT");
  obs_agec_sCT.allocate(1,nyr_agec_sCT,1,nages_agec_10p,"obs_agec_sCT");
  styr_cpue_sVD.allocate("styr_cpue_sVD");
  endyr_cpue_sVD.allocate("endyr_cpue_sVD");
  obs_cpue_sVD.allocate(styr_cpue_sVD,endyr_cpue_sVD,"obs_cpue_sVD");
  obs_cv_cpue_sVD.allocate(styr_cpue_sVD,endyr_cpue_sVD,"obs_cv_cpue_sVD");
  set_Linf.allocate(1,7,"set_Linf");
  set_K.allocate(1,7,"set_K");
  set_t0.allocate(1,7,"set_t0");
  set_len_cv.allocate(1,7,"set_len_cv");
  set_Linf_L.allocate(1,7,"set_Linf_L");
  set_K_L.allocate(1,7,"set_K_L");
  set_t0_L.allocate(1,7,"set_t0_L");
  set_len_cv_L.allocate(1,7,"set_len_cv_L");
  set_Linf_20.allocate(1,7,"set_Linf_20");
  set_K_20.allocate(1,7,"set_K_20");
  set_t0_20.allocate(1,7,"set_t0_20");
  set_len_cv_20.allocate(1,7,"set_len_cv_20");
  set_M_constant.allocate(1,7,"set_M_constant");
  set_steep.allocate(1,7,"set_steep");
  set_log_R0.allocate(1,7,"set_log_R0");
  set_R_autocorr.allocate(1,7,"set_R_autocorr");
  set_rec_sigma.allocate(1,7,"set_rec_sigma");
  set_log_dm_lenc_cHL.allocate(1,7,"set_log_dm_lenc_cHL");
  set_log_dm_lenc_cHL_D.allocate(1,7,"set_log_dm_lenc_cHL_D");
  set_log_dm_lenc_rHB_D.allocate(1,7,"set_log_dm_lenc_rHB_D");
  set_log_dm_lenc_rGN_D.allocate(1,7,"set_log_dm_lenc_rGN_D");
  set_log_dm_agec_cHL.allocate(1,7,"set_log_dm_agec_cHL");
  set_log_dm_agec_rHB.allocate(1,7,"set_log_dm_agec_rHB");
  set_log_dm_agec_sCT.allocate(1,7,"set_log_dm_agec_sCT");
  set_log_dm_agec_rGN.allocate(1,7,"set_log_dm_agec_rGN");
  set_selpar_A50_cHL1.allocate(1,7,"set_selpar_A50_cHL1");
  set_selpar_slope_cHL1.allocate(1,7,"set_selpar_slope_cHL1");
  set_selpar_A50_cHL2.allocate(1,7,"set_selpar_A50_cHL2");
  set_selpar_slope_cHL2.allocate(1,7,"set_selpar_slope_cHL2");
  set_selpar_A50_cHL3.allocate(1,7,"set_selpar_A50_cHL3");
  set_selpar_slope_cHL3.allocate(1,7,"set_selpar_slope_cHL3");
  set_selpar_A50_rHB1.allocate(1,7,"set_selpar_A50_rHB1");
  set_selpar_slope_rHB1.allocate(1,7,"set_selpar_slope_rHB1");
  set_selpar_A502_rHB1.allocate(1,7,"set_selpar_A502_rHB1");
  set_selpar_slope2_rHB1.allocate(1,7,"set_selpar_slope2_rHB1");
  set_selpar_A50_rHB2.allocate(1,7,"set_selpar_A50_rHB2");
  set_selpar_slope_rHB2.allocate(1,7,"set_selpar_slope_rHB2");
  set_selpar_A502_rHB2.allocate(1,7,"set_selpar_A502_rHB2");
  set_selpar_slope2_rHB2.allocate(1,7,"set_selpar_slope2_rHB2");
  set_selpar_A50_rHB3.allocate(1,7,"set_selpar_A50_rHB3");
  set_selpar_slope_rHB3.allocate(1,7,"set_selpar_slope_rHB3");
  set_selpar_A502_rHB3.allocate(1,7,"set_selpar_A502_rHB3");
  set_selpar_slope2_rHB3.allocate(1,7,"set_selpar_slope2_rHB3");
  set_selpar_A50_rGN2.allocate(1,7,"set_selpar_A50_rGN2");
  set_selpar_slope_rGN2.allocate(1,7,"set_selpar_slope_rGN2");
  set_selpar_A502_rGN2.allocate(1,7,"set_selpar_A502_rGN2");
  set_selpar_slope2_rGN2.allocate(1,7,"set_selpar_slope2_rGN2");
  set_selpar_A50_rGN3.allocate(1,7,"set_selpar_A50_rGN3");
  set_selpar_slope_rGN3.allocate(1,7,"set_selpar_slope_rGN3");
  set_selpar_A50_sCT.allocate(1,7,"set_selpar_A50_sCT");
  set_selpar_slope_sCT.allocate(1,7,"set_selpar_slope_sCT");
  set_selpar_A502_sCT.allocate(1,7,"set_selpar_A502_sCT");
  set_selpar_slope2_sCT.allocate(1,7,"set_selpar_slope2_sCT");
  set_selpar_A50_cHL2_D.allocate(1,7,"set_selpar_A50_cHL2_D");
  set_selpar_slope_cHL2_D.allocate(1,7,"set_selpar_slope_cHL2_D");
  set_selpar_A502_cHL2_D.allocate(1,7,"set_selpar_A502_cHL2_D");
  set_selpar_slope2_cHL2_D.allocate(1,7,"set_selpar_slope2_cHL2_D");
  set_selpar_A50_cHL3_D.allocate(1,7,"set_selpar_A50_cHL3_D");
  set_selpar_slope_cHL3_D.allocate(1,7,"set_selpar_slope_cHL3_D");
  set_selpar_A50_rHB2_D.allocate(1,7,"set_selpar_A50_rHB2_D");
  set_selpar_slope_rHB2_D.allocate(1,7,"set_selpar_slope_rHB2_D");
  set_selpar_A502_rHB2_D.allocate(1,7,"set_selpar_A502_rHB2_D");
  set_selpar_slope2_rHB2_D.allocate(1,7,"set_selpar_slope2_rHB2_D");
  set_selpar_A50_rHB3_D.allocate(1,7,"set_selpar_A50_rHB3_D");
  set_selpar_slope_rHB3_D.allocate(1,7,"set_selpar_slope_rHB3_D");
  set_selpar_A502_rHB3_D.allocate(1,7,"set_selpar_A502_rHB3_D");
  set_selpar_slope2_rHB3_D.allocate(1,7,"set_selpar_slope2_rHB3_D");
  set_selpar_A50_rGN3_D.allocate(1,7,"set_selpar_A50_rGN3_D");
  set_selpar_slope_rGN3_D.allocate(1,7,"set_selpar_slope_rGN3_D");
  set_selpar_A502_rGN3_D.allocate(1,7,"set_selpar_A502_rGN3_D");
  set_selpar_slope2_rGN3_D.allocate(1,7,"set_selpar_slope2_rGN3_D");
  set_log_q_cpue_cHL.allocate(1,7,"set_log_q_cpue_cHL");
  set_log_q_cpue_rHB.allocate(1,7,"set_log_q_cpue_rHB");
  set_log_q_cpue_rHB_D.allocate(1,7,"set_log_q_cpue_rHB_D");
  set_log_q_cpue_sCT.allocate(1,7,"set_log_q_cpue_sCT");
  set_log_q_cpue_sVD.allocate(1,7,"set_log_q_cpue_sVD");
  set_F_init.allocate(1,7,"set_F_init");
  set_log_avg_F_L_cHL.allocate(1,7,"set_log_avg_F_L_cHL");
  set_log_avg_F_L_rHB.allocate(1,7,"set_log_avg_F_L_rHB");
  set_log_avg_F_L_rGN.allocate(1,7,"set_log_avg_F_L_rGN");
  set_log_avg_F_D_cHL.allocate(1,7,"set_log_avg_F_D_cHL");
  set_log_avg_F_D_rHB.allocate(1,7,"set_log_avg_F_D_rHB");
  set_log_avg_F_D_rGN.allocate(1,7,"set_log_avg_F_D_rGN");
  set_log_dev_F_L_cHL.allocate(1,3,"set_log_dev_F_L_cHL");
  set_log_dev_F_L_rHB.allocate(1,3,"set_log_dev_F_L_rHB");
  set_log_dev_F_L_rGN.allocate(1,3,"set_log_dev_F_L_rGN");
  set_log_dev_F_D_cHL.allocate(1,3,"set_log_dev_F_D_cHL");
  set_log_dev_F_D_rHB.allocate(1,3,"set_log_dev_F_D_rHB");
  set_log_dev_F_D_rGN.allocate(1,3,"set_log_dev_F_D_rGN");
  set_log_dev_rec.allocate(1,3,"set_log_dev_rec");
  set_log_dev_Nage.allocate(1,3,"set_log_dev_Nage");
  set_log_dev_vals_F_L__cHL.allocate(styr_L_cHL,endyr_L_cHL,"set_log_dev_vals_F_L__cHL");
  set_log_dev_vals_F_L__rHB.allocate(styr_L_rHB,endyr_L_rHB,"set_log_dev_vals_F_L__rHB");
  set_log_dev_vals_F_L__rGN.allocate(styr_L_rGN,endyr_L_rGN,"set_log_dev_vals_F_L__rGN");
  set_log_dev_vals_F_D__cHvals.allocate(styr_D_cHL,endyr_D_cHL,"set_log_dev_vals_F_D__cHvals");
  set_log_dev_vals_F_D__HBvals.allocate(styr_D_rHB,endyr_D_rHB,"set_log_dev_vals_F_D__HBvals");
  set_log_dev_vals_F_D__GRvals.allocate(styr_D_rGN,endyr_D_rGN,"set_log_dev_vals_F_D__GRvals");
  set_log_dev_vals_rec.allocate(styr_rec_dev,endyr_rec_dev,"set_log_dev_vals_rec");
  set_log_dev_vals_Nage.allocate(2,nages,"set_log_dev_vals_Nage");
  set_w_L.allocate("set_w_L");
  set_w_D.allocate("set_w_D");
  set_w_cpue_cHL.allocate("set_w_cpue_cHL");
  set_w_cpue_rHB.allocate("set_w_cpue_rHB");
  set_w_cpue_rHB_D.allocate("set_w_cpue_rHB_D");
  set_w_cpue_sCT.allocate("set_w_cpue_sCT");
  set_w_cpue_sVD.allocate("set_w_cpue_sVD");
  set_w_cpue_sCT_mult.allocate("set_w_cpue_sCT_mult");
  set_w_cpue_sVD_mult.allocate("set_w_cpue_sVD_mult");
  set_w_lenc_cHL.allocate("set_w_lenc_cHL");
  set_w_lenc_cHL_D.allocate("set_w_lenc_cHL_D");
  set_w_lenc_rHB_D.allocate("set_w_lenc_rHB_D");
  set_w_lenc_rGN_D.allocate("set_w_lenc_rGN_D");
  set_w_agec_cHL.allocate("set_w_agec_cHL");
  set_w_agec_rHB.allocate("set_w_agec_rHB");
  set_w_agec_sCT.allocate("set_w_agec_sCT");
  set_w_agec_rGN.allocate("set_w_agec_rGN");
  set_w_Nage_init.allocate("set_w_Nage_init");
  set_w_rec.allocate("set_w_rec");
  set_w_rec_early.allocate("set_w_rec_early");
  set_w_rec_end.allocate("set_w_rec_end");
  set_w_fullF.allocate("set_w_fullF");
  set_w_Ftune.allocate("set_w_Ftune");
  wgtpar_a.allocate("wgtpar_a");
  wgtpar_b.allocate("wgtpar_b");
  ww2gw.allocate("ww2gw");
  obs_maturity_f.allocate(1,nages,"obs_maturity_f");
  obs_prop_f.allocate(1,nages,"obs_prop_f");
  set_fecpar_batches.allocate(1,nages,"set_fecpar_batches");
  set_fecpar_a.allocate("set_fecpar_a");
  set_fecpar_b.allocate("set_fecpar_b");
  set_fecpar_c.allocate("set_fecpar_c");
  set_fecpar_thresh.allocate("set_fecpar_thresh");
  set_fecpar_min.allocate("set_fecpar_min");
  fecpar_scale.allocate("fecpar_scale");
  spawn_time_frac.allocate("spawn_time_frac");
  set_M.allocate(1,nages,"set_M");
  max_obs_age.allocate("max_obs_age");
  endyr_cHL_Dmort1.allocate("endyr_cHL_Dmort1");
  endyr_rHB_Dmort1.allocate("endyr_rHB_Dmort1");
  endyr_rGN_Dmort1.allocate("endyr_rGN_Dmort1");
  endyr_cHL_Dmort2.allocate("endyr_cHL_Dmort2");
  endyr_rHB_Dmort2.allocate("endyr_rHB_Dmort2");
  endyr_rGN_Dmort2.allocate("endyr_rGN_Dmort2");
  set_Dmort_cHL1.allocate("set_Dmort_cHL1");
  set_Dmort_rHB1.allocate("set_Dmort_rHB1");
  set_Dmort_rGN1.allocate("set_Dmort_rGN1");
  set_Dmort_cHL2.allocate("set_Dmort_cHL2");
  set_Dmort_rHB2.allocate("set_Dmort_rHB2");
  set_Dmort_rGN2.allocate("set_Dmort_rGN2");
  set_Dmort_cHL3.allocate("set_Dmort_cHL3");
  set_Dmort_rHB3.allocate("set_Dmort_rHB3");
  set_Dmort_rGN3.allocate("set_Dmort_rGN3");
  historic_Lrec_scale.allocate("historic_Lrec_scale");
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
  minSS_lenc_cHL_D.allocate("minSS_lenc_cHL_D");
  minSS_lenc_rHB_D.allocate("minSS_lenc_rHB_D");
  minSS_lenc_rGN_D.allocate("minSS_lenc_rGN_D");
  minSS_agec_cHL.allocate("minSS_agec_cHL");
  minSS_agec_rHB.allocate("minSS_agec_rHB");
  minSS_agec_sCT.allocate("minSS_agec_sCT");
  minSS_agec_rGN.allocate("minSS_agec_rGN");
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
  
  const double Linf_20_LO=set_Linf_20(2); const double Linf_20_HI=set_Linf_20(3); const double Linf_20_PH=set_Linf_20(4);
  const double K_20_LO=set_K_20(2); const double K_20_HI=set_K_20(3); const double K_20_PH=set_K_20(4);
  const double t0_20_LO=set_t0_20(2); const double t0_20_HI=set_t0_20(3); const double t0_20_PH=set_t0_20(4);  
  const double len_cv_20_LO=set_len_cv_20(2); const double len_cv_20_HI=set_len_cv_20(3); const double len_cv_20_PH=set_len_cv_20(4);
   
  const double M_constant_LO=set_M_constant(2); const double M_constant_HI=set_M_constant(3); const double M_constant_PH=set_M_constant(4);        
  const double steep_LO=set_steep(2); const double steep_HI=set_steep(3); const double steep_PH=set_steep(4);
  const double log_R0_LO=set_log_R0(2); const double log_R0_HI=set_log_R0(3); const double log_R0_PH=set_log_R0(4);
  const double R_autocorr_LO=set_R_autocorr(2); const double R_autocorr_HI=set_R_autocorr(3); const double R_autocorr_PH=set_R_autocorr(4);
  const double rec_sigma_LO=set_rec_sigma(2); const double rec_sigma_HI=set_rec_sigma(3); const double rec_sigma_PH=set_rec_sigma(4);
  const double log_dm_cHL_lc_LO=set_log_dm_lenc_cHL(2); const double log_dm_cHL_lc_HI=set_log_dm_lenc_cHL(3); const double log_dm_cHL_lc_PH=set_log_dm_lenc_cHL(4);
  const double log_dm_cHL_D_lc_LO=set_log_dm_lenc_cHL_D(2); const double log_dm_cHL_D_lc_HI=set_log_dm_lenc_cHL_D(3); const double log_dm_cHL_D_lc_PH=set_log_dm_lenc_cHL_D(4);
  const double log_dm_rHB_D_lc_LO=set_log_dm_lenc_rHB_D(2); const double log_dm_rHB_D_lc_HI=set_log_dm_lenc_rHB_D(3); const double log_dm_rHB_D_lc_PH=set_log_dm_lenc_rHB_D(4);
  const double log_dm_rGN_D_lc_LO=set_log_dm_lenc_rGN_D(2); const double log_dm_rGN_D_lc_HI=set_log_dm_lenc_rGN_D(3); const double log_dm_rGN_D_lc_PH=set_log_dm_lenc_rGN_D(4);
  const double log_dm_cHL_ac_LO=set_log_dm_agec_cHL(2); const double log_dm_cHL_ac_HI=set_log_dm_agec_cHL(3); const double log_dm_cHL_ac_PH=set_log_dm_agec_cHL(4);
  const double log_dm_rHB_ac_LO=set_log_dm_agec_rHB(2); const double log_dm_rHB_ac_HI=set_log_dm_agec_rHB(3); const double log_dm_rHB_ac_PH=set_log_dm_agec_rHB(4);
  const double log_dm_sCT_ac_LO=set_log_dm_agec_sCT(2); const double log_dm_sCT_ac_HI=set_log_dm_agec_sCT(3); const double log_dm_sCT_ac_PH=set_log_dm_agec_sCT(4);
  const double log_dm_rGN_ac_LO=set_log_dm_agec_rGN(2); const double log_dm_rGN_ac_HI=set_log_dm_agec_rGN(3); const double log_dm_rGN_ac_PH=set_log_dm_agec_rGN(4);
  
  const double selpar_A50_cHL1_LO=set_selpar_A50_cHL1(2); const double selpar_A50_cHL1_HI=set_selpar_A50_cHL1(3); const double selpar_A50_cHL1_PH=set_selpar_A50_cHL1(4);
  const double selpar_slope_cHL1_LO=set_selpar_slope_cHL1(2); const double selpar_slope_cHL1_HI=set_selpar_slope_cHL1(3); const double selpar_slope_cHL1_PH=set_selpar_slope_cHL1(4);
  const double selpar_A50_cHL2_LO=set_selpar_A50_cHL2(2); const double selpar_A50_cHL2_HI=set_selpar_A50_cHL2(3); const double selpar_A50_cHL2_PH=set_selpar_A50_cHL2(4);
  const double selpar_slope_cHL2_LO=set_selpar_slope_cHL2(2); const double selpar_slope_cHL2_HI=set_selpar_slope_cHL2(3); const double selpar_slope_cHL2_PH=set_selpar_slope_cHL2(4);
  const double selpar_A50_cHL3_LO=set_selpar_A50_cHL3(2); const double selpar_A50_cHL3_HI=set_selpar_A50_cHL3(3); const double selpar_A50_cHL3_PH=set_selpar_A50_cHL3(4);
  const double selpar_slope_cHL3_LO=set_selpar_slope_cHL3(2); const double selpar_slope_cHL3_HI=set_selpar_slope_cHL3(3); const double selpar_slope_cHL3_PH=set_selpar_slope_cHL3(4);
  
  const double selpar_A50_rHB1_LO=set_selpar_A50_rHB1(2); const double selpar_A50_rHB1_HI=set_selpar_A50_rHB1(3); const double selpar_A50_rHB1_PH=set_selpar_A50_rHB1(4);
  const double selpar_slope_rHB1_LO=set_selpar_slope_rHB1(2); const double selpar_slope_rHB1_HI=set_selpar_slope_rHB1(3); const double selpar_slope_rHB1_PH=set_selpar_slope_rHB1(4);
  const double selpar_A502_rHB1_LO=set_selpar_A502_rHB1(2); const double selpar_A502_rHB1_HI=set_selpar_A502_rHB1(3); const double selpar_A502_rHB1_PH=set_selpar_A502_rHB1(4);
  const double selpar_slope2_rHB1_LO=set_selpar_slope2_rHB1(2); const double selpar_slope2_rHB1_HI=set_selpar_slope2_rHB1(3); const double selpar_slope2_rHB1_PH=set_selpar_slope2_rHB1(4);
  const double selpar_A50_rHB2_LO=set_selpar_A50_rHB2(2); const double selpar_A50_rHB2_HI=set_selpar_A50_rHB2(3); const double selpar_A50_rHB2_PH=set_selpar_A50_rHB2(4);
  const double selpar_slope_rHB2_LO=set_selpar_slope_rHB2(2); const double selpar_slope_rHB2_HI=set_selpar_slope_rHB2(3); const double selpar_slope_rHB2_PH=set_selpar_slope_rHB2(4);
  const double selpar_A502_rHB2_LO=set_selpar_A502_rHB2(2); const double selpar_A502_rHB2_HI=set_selpar_A502_rHB2(3); const double selpar_A502_rHB2_PH=set_selpar_A502_rHB2(4);
  const double selpar_slope2_rHB2_LO=set_selpar_slope2_rHB2(2); const double selpar_slope2_rHB2_HI=set_selpar_slope2_rHB2(3); const double selpar_slope2_rHB2_PH=set_selpar_slope2_rHB2(4);
  const double selpar_A50_rHB3_LO=set_selpar_A50_rHB3(2); const double selpar_A50_rHB3_HI=set_selpar_A50_rHB3(3); const double selpar_A50_rHB3_PH=set_selpar_A50_rHB3(4);
  const double selpar_slope_rHB3_LO=set_selpar_slope_rHB3(2); const double selpar_slope_rHB3_HI=set_selpar_slope_rHB3(3); const double selpar_slope_rHB3_PH=set_selpar_slope_rHB3(4);
  const double selpar_A502_rHB3_LO=set_selpar_A502_rHB3(2); const double selpar_A502_rHB3_HI=set_selpar_A502_rHB3(3); const double selpar_A502_rHB3_PH=set_selpar_A502_rHB3(4);
  const double selpar_slope2_rHB3_LO=set_selpar_slope2_rHB3(2); const double selpar_slope2_rHB3_HI=set_selpar_slope2_rHB3(3); const double selpar_slope2_rHB3_PH=set_selpar_slope2_rHB3(4);
  
  const double selpar_A50_rGN2_LO=set_selpar_A50_rGN2(2); const double selpar_A50_rGN2_HI=set_selpar_A50_rGN2(3); const double selpar_A50_rGN2_PH=set_selpar_A50_rGN2(4);
  const double selpar_slope_rGN2_LO=set_selpar_slope_rGN2(2); const double selpar_slope_rGN2_HI=set_selpar_slope_rGN2(3); const double selpar_slope_rGN2_PH=set_selpar_slope_rGN2(4);
  const double selpar_A502_rGN2_LO=set_selpar_A502_rGN2(2); const double selpar_A502_rGN2_HI=set_selpar_A502_rGN2(3); const double selpar_A502_rGN2_PH=set_selpar_A502_rGN2(4);
  const double selpar_slope2_rGN2_LO=set_selpar_slope2_rGN2(2); const double selpar_slope2_rGN2_HI=set_selpar_slope2_rGN2(3); const double selpar_slope2_rGN2_PH=set_selpar_slope2_rGN2(4);
  const double selpar_A50_rGN3_LO=set_selpar_A50_rGN3(2); const double selpar_A50_rGN3_HI=set_selpar_A50_rGN3(3); const double selpar_A50_rGN3_PH=set_selpar_A50_rGN3(4);
  const double selpar_slope_rGN3_LO=set_selpar_slope_rGN3(2); const double selpar_slope_rGN3_HI=set_selpar_slope_rGN3(3); const double selpar_slope_rGN3_PH=set_selpar_slope_rGN3(4);
  
  const double selpar_A50_sCT_LO=set_selpar_A50_sCT(2); const double selpar_A50_sCT_HI=set_selpar_A50_sCT(3); const double selpar_A50_sCT_PH=set_selpar_A50_sCT(4);
  const double selpar_slope_sCT_LO=set_selpar_slope_sCT(2); const double selpar_slope_sCT_HI=set_selpar_slope_sCT(3); const double selpar_slope_sCT_PH=set_selpar_slope_sCT(4);
  const double selpar_A502_sCT_LO=set_selpar_A502_sCT(2); const double selpar_A502_sCT_HI=set_selpar_A502_sCT(3); const double selpar_A502_sCT_PH=set_selpar_A502_sCT(4);
  const double selpar_slope2_sCT_LO=set_selpar_slope2_sCT(2); const double selpar_slope2_sCT_HI=set_selpar_slope2_sCT(3); const double selpar_slope2_sCT_PH=set_selpar_slope2_sCT(4);
  
  const double selpar_A50_rHB2_D_LO=set_selpar_A50_rHB2_D(2); const double selpar_A50_rHB2_D_HI=set_selpar_A50_rHB2_D(3); const double selpar_A50_rHB2_D_PH=set_selpar_A50_rHB2_D(4);
  const double selpar_slope_rHB2_D_LO=set_selpar_slope_rHB2_D(2); const double selpar_slope_rHB2_D_HI=set_selpar_slope_rHB2_D(3); const double selpar_slope_rHB2_D_PH=set_selpar_slope_rHB2_D(4);
  const double selpar_A502_rHB2_D_LO=set_selpar_A502_rHB2_D(2); const double selpar_A502_rHB2_D_HI=set_selpar_A502_rHB2_D(3); const double selpar_A502_rHB2_D_PH=set_selpar_A502_rHB2_D(4);
  const double selpar_slope2_rHB2_D_LO=set_selpar_slope2_rHB2_D(2); const double selpar_slope2_rHB2_D_HI=set_selpar_slope2_rHB2_D(3); const double selpar_slope2_rHB2_D_PH=set_selpar_slope2_rHB2_D(4);
  
  const double selpar_A50_rHB3_D_LO=set_selpar_A50_rHB3_D(2); const double selpar_A50_rHB3_D_HI=set_selpar_A50_rHB3_D(3); const double selpar_A50_rHB3_D_PH=set_selpar_A50_rHB3_D(4);
  const double selpar_slope_rHB3_D_LO=set_selpar_slope_rHB3_D(2); const double selpar_slope_rHB3_D_HI=set_selpar_slope_rHB3_D(3); const double selpar_slope_rHB3_D_PH=set_selpar_slope_rHB3_D(4);
  const double selpar_A502_rHB3_D_LO=set_selpar_A502_rHB3_D(2); const double selpar_A502_rHB3_D_HI=set_selpar_A502_rHB3_D(3); const double selpar_A502_rHB3_D_PH=set_selpar_A502_rHB3_D(4);
  const double selpar_slope2_rHB3_D_LO=set_selpar_slope2_rHB3_D(2); const double selpar_slope2_rHB3_D_HI=set_selpar_slope2_rHB3_D(3); const double selpar_slope2_rHB3_D_PH=set_selpar_slope2_rHB3_D(4);
  const double selpar_A50_rGN3_D_LO=set_selpar_A50_rGN3_D(2); const double selpar_A50_rGN3_D_HI=set_selpar_A50_rGN3_D(3); const double selpar_A50_rGN3_D_PH=set_selpar_A50_rGN3_D(4);
  const double selpar_slope_rGN3_D_LO=set_selpar_slope_rGN3_D(2); const double selpar_slope_rGN3_D_HI=set_selpar_slope_rGN3_D(3); const double selpar_slope_rGN3_D_PH=set_selpar_slope_rGN3_D(4);
  const double selpar_A502_rGN3_D_LO=set_selpar_A502_rGN3_D(2); const double selpar_A502_rGN3_D_HI=set_selpar_A502_rGN3_D(3); const double selpar_A502_rGN3_D_PH=set_selpar_A502_rGN3_D(4);
  const double selpar_slope2_rGN3_D_LO=set_selpar_slope2_rGN3_D(2); const double selpar_slope2_rGN3_D_HI=set_selpar_slope2_rGN3_D(3); const double selpar_slope2_rGN3_D_PH=set_selpar_slope2_rGN3_D(4);
  
  const double selpar_A50_cHL2_D_LO=set_selpar_A50_cHL2_D(2); const double selpar_A50_cHL2_D_HI=set_selpar_A50_cHL2_D(3); const double selpar_A50_cHL2_D_PH=set_selpar_A50_cHL2_D(4);
  const double selpar_slope_cHL2_D_LO=set_selpar_slope_cHL2_D(2); const double selpar_slope_cHL2_D_HI=set_selpar_slope_cHL2_D(3); const double selpar_slope_cHL2_D_PH=set_selpar_slope_cHL2_D(4);
  const double selpar_A502_cHL2_D_LO=set_selpar_A502_cHL2_D(2); const double selpar_A502_cHL2_D_HI=set_selpar_A502_cHL2_D(3); const double selpar_A502_cHL2_D_PH=set_selpar_A502_cHL2_D(4);
  const double selpar_slope2_cHL2_D_LO=set_selpar_slope2_cHL2_D(2); const double selpar_slope2_cHL2_D_HI=set_selpar_slope2_cHL2_D(3); const double selpar_slope2_cHL2_D_PH=set_selpar_slope2_cHL2_D(4);
  
  const double selpar_A50_cHL3_D_LO=set_selpar_A50_cHL3_D(2); const double selpar_A50_cHL3_D_HI=set_selpar_A50_cHL3_D(3); const double selpar_A50_cHL3_D_PH=set_selpar_A50_cHL3_D(4);
  const double selpar_slope_cHL3_D_LO=set_selpar_slope_cHL3_D(2); const double selpar_slope_cHL3_D_HI=set_selpar_slope_cHL3_D(3); const double selpar_slope_cHL3_D_PH=set_selpar_slope_cHL3_D(4);
  
  const double log_q_cHL_LO=set_log_q_cpue_cHL(2); const double log_q_cHL_HI=set_log_q_cpue_cHL(3); const double log_q_cHL_PH=set_log_q_cpue_cHL(4);
  const double log_q_rHB_LO=set_log_q_cpue_rHB(2); const double log_q_rHB_HI=set_log_q_cpue_rHB(3); const double log_q_rHB_PH=set_log_q_cpue_rHB(4);
  const double log_q_rHB_D_LO=set_log_q_cpue_rHB_D(2); const double log_q_rHB_D_HI=set_log_q_cpue_rHB_D(3); const double log_q_rHB_D_PH=set_log_q_cpue_rHB_D(4);
  const double log_q_sCT_LO=set_log_q_cpue_sCT(2); const double log_q_sCT_HI=set_log_q_cpue_sCT(3); const double log_q_sCT_PH=set_log_q_cpue_sCT(4);
  const double log_q_sVD_LO=set_log_q_cpue_sVD(2); const double log_q_sVD_HI=set_log_q_cpue_sVD(3); const double log_q_sVD_PH=set_log_q_cpue_sVD(4);
  
  const double F_init_LO=set_F_init(2); const double F_init_HI=set_F_init(3); const double F_init_PH=set_F_init(4);
  const double log_avg_F_cHL_LO=set_log_avg_F_L_cHL(2); const double log_avg_F_cHL_HI=set_log_avg_F_L_cHL(3); const double log_avg_F_cHL_PH=set_log_avg_F_L_cHL(4);
  const double log_avg_F_rHB_LO=set_log_avg_F_L_rHB(2); const double log_avg_F_rHB_HI=set_log_avg_F_L_rHB(3); const double log_avg_F_rHB_PH=set_log_avg_F_L_rHB(4); 
  const double log_avg_F_rGN_LO=set_log_avg_F_L_rGN(2); const double log_avg_F_rGN_HI=set_log_avg_F_L_rGN(3); const double log_avg_F_rGN_PH=set_log_avg_F_L_rGN(4); 
  const double log_avg_F_cHL_D_LO=set_log_avg_F_D_cHL(2); const double log_avg_F_cHL_D_HI=set_log_avg_F_D_cHL(3); const double log_avg_F_cHL_D_PH=set_log_avg_F_D_cHL(4);
  const double log_avg_F_rHB_D_LO=set_log_avg_F_D_rHB(2); const double log_avg_F_rHB_D_HI=set_log_avg_F_D_rHB(3); const double log_avg_F_rHB_D_PH=set_log_avg_F_D_rHB(4); 
  const double log_avg_F_rGN_D_LO=set_log_avg_F_D_rGN(2); const double log_avg_F_rGN_D_HI=set_log_avg_F_D_rGN(3); const double log_avg_F_rGN_D_PH=set_log_avg_F_D_rGN(4); 
  
  //-dev vectors-----------------------------------------------------------------------------------------------------------  
  const double log_F_dev_cHL_LO=set_log_dev_F_L_cHL(1); const double log_F_dev_cHL_HI=set_log_dev_F_L_cHL(2); const double log_F_dev_cHL_PH=set_log_dev_F_L_cHL(3);   
  const double log_F_dev_rHB_LO=set_log_dev_F_L_rHB(1); const double log_F_dev_rHB_HI=set_log_dev_F_L_rHB(2); const double log_F_dev_rHB_PH=set_log_dev_F_L_rHB(3);   
  const double log_F_dev_rGN_LO=set_log_dev_F_L_rGN(1); const double log_F_dev_rGN_HI=set_log_dev_F_L_rGN(2); const double log_F_dev_rGN_PH=set_log_dev_F_L_rGN(3);   
  
  const double log_F_dev_cHL_D_LO=set_log_dev_F_D_cHL(1); const double log_F_dev_cHL_D_HI=set_log_dev_F_D_cHL(2); const double log_F_dev_cHL_D_PH=set_log_dev_F_D_cHL(3);   
  const double log_F_dev_rHB_D_LO=set_log_dev_F_D_rHB(1); const double log_F_dev_rHB_D_HI=set_log_dev_F_D_rHB(2); const double log_F_dev_rHB_D_PH=set_log_dev_F_D_rHB(3);   
  const double log_F_dev_rGN_D_LO=set_log_dev_F_D_rGN(1); const double log_F_dev_rGN_D_HI=set_log_dev_F_D_rGN(2); const double log_F_dev_rGN_D_PH=set_log_dev_F_D_rGN(3);   
  
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
  fecundity.allocate(1,nages,"fecundity");
  #ifndef NO_AD_INITIALIZE
    fecundity.initialize();
  #endif
  fecpar_batches.allocate(1,nages,"fecpar_batches");
  #ifndef NO_AD_INITIALIZE
    fecpar_batches.initialize();
  #endif
  fecpar_a.allocate("fecpar_a");
  #ifndef NO_AD_INITIALIZE
  fecpar_a.initialize();
  #endif
  fecpar_b.allocate("fecpar_b");
  #ifndef NO_AD_INITIALIZE
  fecpar_b.initialize();
  #endif
  fecpar_c.allocate("fecpar_c");
  #ifndef NO_AD_INITIALIZE
  fecpar_c.initialize();
  #endif
  fecpar_thresh.allocate("fecpar_thresh");
  #ifndef NO_AD_INITIALIZE
  fecpar_thresh.initialize();
  #endif
  fecpar_min.allocate("fecpar_min");
  #ifndef NO_AD_INITIALIZE
  fecpar_min.initialize();
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
  Linf_20.allocate(Linf_20_LO,Linf_20_HI,Linf_20_PH,"Linf_20");
  K_20.allocate(K_20_LO,K_20_HI,K_20_PH,"K_20");
  t0_20.allocate(t0_20_LO,t0_20_HI,t0_20_PH,"t0_20");
  len_cv_val_20.allocate(len_cv_20_LO,len_cv_20_HI,len_cv_20_PH,"len_cv_val_20");
  Linf_20_out.allocate(1,8,"Linf_20_out");
  #ifndef NO_AD_INITIALIZE
    Linf_20_out.initialize();
  #endif
  K_20_out.allocate(1,8,"K_20_out");
  #ifndef NO_AD_INITIALIZE
    K_20_out.initialize();
  #endif
  t0_20_out.allocate(1,8,"t0_20_out");
  #ifndef NO_AD_INITIALIZE
    t0_20_out.initialize();
  #endif
  len_cv_val_20_out.allocate(1,8,"len_cv_val_20_out");
  #ifndef NO_AD_INITIALIZE
    len_cv_val_20_out.initialize();
  #endif
  meanlen_TL_20.allocate(1,nages,"meanlen_TL_20");
  #ifndef NO_AD_INITIALIZE
    meanlen_TL_20.initialize();
  #endif
  wgt_g_20.allocate(1,nages,"wgt_g_20");
  #ifndef NO_AD_INITIALIZE
    wgt_g_20.initialize();
  #endif
  wgt_kg_20.allocate(1,nages,"wgt_kg_20");
  #ifndef NO_AD_INITIALIZE
    wgt_kg_20.initialize();
  #endif
  wgt_mt_20.allocate(1,nages,"wgt_mt_20");
  #ifndef NO_AD_INITIALIZE
    wgt_mt_20.initialize();
  #endif
  wgt_klb_20.allocate(1,nages,"wgt_klb_20");
  #ifndef NO_AD_INITIALIZE
    wgt_klb_20.initialize();
  #endif
  wgt_lb_20.allocate(1,nages,"wgt_lb_20");
  #ifndef NO_AD_INITIALIZE
    wgt_lb_20.initialize();
  #endif
  len_cHL_mm.allocate(styr,endyr,1,nages,"len_cHL_mm");
  #ifndef NO_AD_INITIALIZE
    len_cHL_mm.initialize();
  #endif
  wholewgt_cHL_klb.allocate(styr,endyr,1,nages,"wholewgt_cHL_klb");
  #ifndef NO_AD_INITIALIZE
    wholewgt_cHL_klb.initialize();
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
  len_cHL_D_mm.allocate(styr,endyr,1,nages,"len_cHL_D_mm");
  #ifndef NO_AD_INITIALIZE
    len_cHL_D_mm.initialize();
  #endif
  wholewgt_cHL_D_klb.allocate(styr,endyr,1,nages,"wholewgt_cHL_D_klb");
  #ifndef NO_AD_INITIALIZE
    wholewgt_cHL_D_klb.initialize();
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
  zscore_len_20.allocate("zscore_len_20");
  #ifndef NO_AD_INITIALIZE
  zscore_len_20.initialize();
  #endif
  cprob_lenvec_20.allocate(1,nlenbins,"cprob_lenvec_20");
  #ifndef NO_AD_INITIALIZE
    cprob_lenvec_20.initialize();
  #endif
  zscore_lzero_20.allocate("zscore_lzero_20");
  #ifndef NO_AD_INITIALIZE
  zscore_lzero_20.initialize();
  #endif
  cprob_lzero_20.allocate("cprob_lzero_20");
  #ifndef NO_AD_INITIALIZE
  cprob_lzero_20.initialize();
  #endif
  lenprob_cHL.allocate(1,nages,1,nlenbins,"lenprob_cHL");
  #ifndef NO_AD_INITIALIZE
    lenprob_cHL.initialize();
  #endif
  lenprob_cHL_D.allocate(1,nages,1,nlenbins,"lenprob_cHL_D");
  #ifndef NO_AD_INITIALIZE
    lenprob_cHL_D.initialize();
  #endif
  lenprob_rHB.allocate(1,nages,1,nlenbins,"lenprob_rHB");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB.initialize();
  #endif
  lenprob_rHB_D.allocate(1,nages,1,nlenbins,"lenprob_rHB_D");
  #ifndef NO_AD_INITIALIZE
    lenprob_rHB_D.initialize();
  #endif
  lenprob_rGN_D.allocate(1,nages,1,nlenbins,"lenprob_rGN_D");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN_D.initialize();
  #endif
  lenprob_sCT.allocate(1,nages,1,nlenbins,"lenprob_sCT");
  #ifndef NO_AD_INITIALIZE
    lenprob_sCT.initialize();
  #endif
  lenprob_rGN.allocate(1,nages,1,nlenbins,"lenprob_rGN");
  #ifndef NO_AD_INITIALIZE
    lenprob_rGN.initialize();
  #endif
  lenprob_20.allocate(1,nages,1,nlenbins,"lenprob_20");
  #ifndef NO_AD_INITIALIZE
    lenprob_20.initialize();
  #endif
  lenprob_L.allocate(1,nages,1,nlenbins,"lenprob_L");
  #ifndef NO_AD_INITIALIZE
    lenprob_L.initialize();
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
  len_sd_20.allocate(1,nages,"len_sd_20");
  #ifndef NO_AD_INITIALIZE
    len_sd_20.initialize();
  #endif
  len_cv_20.allocate(1,nages,"len_cv_20");
  #ifndef NO_AD_INITIALIZE
    len_cv_20.initialize();
  #endif
  pred_cHL_lenc.allocate(1,nyr_lenc_cHL,1,nlenbins,"pred_cHL_lenc");
  #ifndef NO_AD_INITIALIZE
    pred_cHL_lenc.initialize();
  #endif
  pred_cHL_D_lenc_yr1.allocate(1,nyr_cHL_D_lenc_pool1,1,nlenbins,"pred_cHL_D_lenc_yr1");
  #ifndef NO_AD_INITIALIZE
    pred_cHL_D_lenc_yr1.initialize();
  #endif
  pred_cHL_D_lenc_yr2.allocate(1,nyr_cHL_D_lenc_pool2,1,nlenbins,"pred_cHL_D_lenc_yr2");
  #ifndef NO_AD_INITIALIZE
    pred_cHL_D_lenc_yr2.initialize();
  #endif
  pred_cHL_D_lenc.allocate(1,nyr_lenc_cHL_D,1,nlenbins,"pred_cHL_D_lenc");
  #ifndef NO_AD_INITIALIZE
    pred_cHL_D_lenc.initialize();
  #endif
  pred_rHB_D_lenc.allocate(1,nyr_lenc_rHB_D,1,nlenbins,"pred_rHB_D_lenc");
  #ifndef NO_AD_INITIALIZE
    pred_rHB_D_lenc.initialize();
  #endif
  pred_rGN_D_lenc.allocate(1,nyr_lenc_rGN_D,1,nlenbins,"pred_rGN_D_lenc");
  #ifndef NO_AD_INITIALIZE
    pred_rGN_D_lenc.initialize();
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
  pred_rHB_agec.allocate(1,nyr_agec_rHB,1,nages_agec_10p,"pred_rHB_agec");
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
  pred_sCT_agec.allocate(1,nyr_agec_sCT,1,nages_agec_10p,"pred_sCT_agec");
  #ifndef NO_AD_INITIALIZE
    pred_sCT_agec.initialize();
  #endif
  pred_sCT_agec_allages.allocate(1,nyr_agec_sCT,1,nages,"pred_sCT_agec_allages");
  #ifndef NO_AD_INITIALIZE
    pred_sCT_agec_allages.initialize();
  #endif
  ErrorFree_sCT_agec.allocate(1,nyr_agec_sCT,1,nages,"ErrorFree_sCT_agec");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_sCT_agec.initialize();
  #endif
  pred_rGN_agec.allocate(1,nyr_agec_rGN,1,nages_agec_10p,"pred_rGN_agec");
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
  nsamp_cHL_lenc_allyr.allocate(styr,endyr,"nsamp_cHL_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_cHL_lenc_allyr.initialize();
  #endif
  nsamp_cHL_D_lenc_allyr.allocate(styr,endyr,"nsamp_cHL_D_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_cHL_D_lenc_allyr.initialize();
  #endif
  nsamp_rHB_D_lenc_allyr.allocate(styr,endyr,"nsamp_rHB_D_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_rHB_D_lenc_allyr.initialize();
  #endif
  nsamp_rGN_D_lenc_allyr.allocate(styr,endyr,"nsamp_rGN_D_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_rGN_D_lenc_allyr.initialize();
  #endif
  nsamp_cHL_agec_allyr.allocate(styr,endyr,"nsamp_cHL_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_cHL_agec_allyr.initialize();
  #endif
  nsamp_rHB_agec_allyr.allocate(styr,endyr,"nsamp_rHB_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_rHB_agec_allyr.initialize();
  #endif
  nsamp_sCT_agec_allyr.allocate(styr,endyr,"nsamp_sCT_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_sCT_agec_allyr.initialize();
  #endif
  nsamp_rGN_agec_allyr.allocate(styr,endyr,"nsamp_rGN_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_rGN_agec_allyr.initialize();
  #endif
  nfish_cHL_lenc_allyr.allocate(styr,endyr,"nfish_cHL_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_cHL_lenc_allyr.initialize();
  #endif
  nfish_cHL_D_lenc_allyr.allocate(styr,endyr,"nfish_cHL_D_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_cHL_D_lenc_allyr.initialize();
  #endif
  nfish_rHB_D_lenc_allyr.allocate(styr,endyr,"nfish_rHB_D_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_rHB_D_lenc_allyr.initialize();
  #endif
  nfish_rGN_D_lenc_allyr.allocate(styr,endyr,"nfish_rGN_D_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_rGN_D_lenc_allyr.initialize();
  #endif
  nfish_cHL_agec_allyr.allocate(styr,endyr,"nfish_cHL_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_cHL_agec_allyr.initialize();
  #endif
  nfish_rHB_agec_allyr.allocate(styr,endyr,"nfish_rHB_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_rHB_agec_allyr.initialize();
  #endif
  nfish_sCT_agec_allyr.allocate(styr,endyr,"nfish_sCT_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_sCT_agec_allyr.initialize();
  #endif
  nfish_rGN_agec_allyr.allocate(styr,endyr,"nfish_rGN_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_rGN_agec_allyr.initialize();
  #endif
  neff_cHL_lenc_allyr_out.allocate(styr,endyr,"neff_cHL_lenc_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_cHL_lenc_allyr_out.initialize();
  #endif
  neff_rHB_D_lenc_allyr_out.allocate(styr,endyr,"neff_rHB_D_lenc_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_rHB_D_lenc_allyr_out.initialize();
  #endif
  neff_rGN_D_lenc_allyr_out.allocate(styr,endyr,"neff_rGN_D_lenc_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_rGN_D_lenc_allyr_out.initialize();
  #endif
  neff_cHL_agec_allyr_out.allocate(styr,endyr,"neff_cHL_agec_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_cHL_agec_allyr_out.initialize();
  #endif
  neff_rHB_agec_allyr_out.allocate(styr,endyr,"neff_rHB_agec_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_rHB_agec_allyr_out.initialize();
  #endif
  neff_sCT_agec_allyr_out.allocate(styr,endyr,"neff_sCT_agec_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_sCT_agec_allyr_out.initialize();
  #endif
  neff_rGN_agec_allyr_out.allocate(styr,endyr,"neff_rGN_agec_allyr_out");
  #ifndef NO_AD_INITIALIZE
    neff_rGN_agec_allyr_out.initialize();
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
  maturity_f.allocate(1,nages,"maturity_f");
  #ifndef NO_AD_INITIALIZE
    maturity_f.initialize();
  #endif
  reprod.allocate(1,nages,"reprod");
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
  log_dm_lenc_cHL.allocate(log_dm_cHL_lc_LO,log_dm_cHL_lc_HI,log_dm_cHL_lc_PH,"log_dm_lenc_cHL");
  log_dm_lenc_cHL_D.allocate(log_dm_cHL_D_lc_LO,log_dm_cHL_D_lc_HI,log_dm_cHL_D_lc_PH,"log_dm_lenc_cHL_D");
  log_dm_lenc_rHB_D.allocate(log_dm_rHB_D_lc_LO,log_dm_rHB_D_lc_HI,log_dm_rHB_D_lc_PH,"log_dm_lenc_rHB_D");
  log_dm_lenc_rGN_D.allocate(log_dm_rGN_D_lc_LO,log_dm_rGN_D_lc_HI,log_dm_rGN_D_lc_PH,"log_dm_lenc_rGN_D");
  log_dm_agec_cHL.allocate(log_dm_cHL_ac_LO,log_dm_cHL_ac_HI,log_dm_cHL_ac_PH,"log_dm_agec_cHL");
  log_dm_agec_rHB.allocate(log_dm_rHB_ac_LO,log_dm_rHB_ac_HI,log_dm_rHB_ac_PH,"log_dm_agec_rHB");
  log_dm_agec_sCT.allocate(log_dm_sCT_ac_LO,log_dm_sCT_ac_HI,log_dm_sCT_ac_PH,"log_dm_agec_sCT");
  log_dm_agec_rGN.allocate(log_dm_rGN_ac_LO,log_dm_rGN_ac_HI,log_dm_rGN_ac_PH,"log_dm_agec_rGN");
  log_dm_cHL_lc_out.allocate(1,8,"log_dm_cHL_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_cHL_lc_out.initialize();
  #endif
  log_dm_cHL_D_lc_out.allocate(1,8,"log_dm_cHL_D_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_cHL_D_lc_out.initialize();
  #endif
  log_dm_rHB_D_lc_out.allocate(1,8,"log_dm_rHB_D_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_rHB_D_lc_out.initialize();
  #endif
  log_dm_rGN_D_lc_out.allocate(1,8,"log_dm_rGN_D_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_rGN_D_lc_out.initialize();
  #endif
  log_dm_cHL_ac_out.allocate(1,8,"log_dm_cHL_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_cHL_ac_out.initialize();
  #endif
  log_dm_rHB_ac_out.allocate(1,8,"log_dm_rHB_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_rHB_ac_out.initialize();
  #endif
  log_dm_sCT_ac_out.allocate(1,8,"log_dm_sCT_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_sCT_ac_out.initialize();
  #endif
  log_dm_rGN_ac_out.allocate(1,8,"log_dm_rGN_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_rGN_ac_out.initialize();
  #endif
  sel_cHL.allocate(styr,endyr,1,nages,"sel_cHL");
  #ifndef NO_AD_INITIALIZE
    sel_cHL.initialize();
  #endif
  selpar_A50_cHL1.allocate(selpar_A50_cHL1_LO,selpar_A50_cHL1_HI,selpar_A50_cHL1_PH,"selpar_A50_cHL1");
  selpar_slope_cHL1.allocate(selpar_slope_cHL1_LO,selpar_slope_cHL1_HI,selpar_slope_cHL1_PH,"selpar_slope_cHL1");
  selpar_A50_cHL2.allocate(selpar_A50_cHL2_LO,selpar_A50_cHL2_HI,selpar_A50_cHL2_PH,"selpar_A50_cHL2");
  selpar_slope_cHL2.allocate(selpar_slope_cHL2_LO,selpar_slope_cHL2_HI,selpar_slope_cHL2_PH,"selpar_slope_cHL2");
  selpar_A50_cHL3.allocate(selpar_A50_cHL3_LO,selpar_A50_cHL3_HI,selpar_A50_cHL3_PH,"selpar_A50_cHL3");
  selpar_slope_cHL3.allocate(selpar_slope_cHL3_LO,selpar_slope_cHL3_HI,selpar_slope_cHL3_PH,"selpar_slope_cHL3");
  selpar_A50_cHL1_out.allocate(1,8,"selpar_A50_cHL1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cHL1_out.initialize();
  #endif
  selpar_slope_cHL1_out.allocate(1,8,"selpar_slope_cHL1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cHL1_out.initialize();
  #endif
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
  selpar_A50_rHB1.allocate(selpar_A50_rHB1_LO,selpar_A50_rHB1_HI,selpar_A50_rHB1_PH,"selpar_A50_rHB1");
  selpar_slope_rHB1.allocate(selpar_slope_rHB1_LO,selpar_slope_rHB1_HI,selpar_slope_rHB1_PH,"selpar_slope_rHB1");
  selpar_A502_rHB1.allocate(selpar_A502_rHB1_LO,selpar_A502_rHB1_HI,selpar_A502_rHB1_PH,"selpar_A502_rHB1");
  selpar_slope2_rHB1.allocate(selpar_slope2_rHB1_LO,selpar_slope2_rHB1_HI,selpar_slope2_rHB1_PH,"selpar_slope2_rHB1");
  selpar_A50_rHB2.allocate(selpar_A50_rHB2_LO,selpar_A50_rHB2_HI,selpar_A50_rHB2_PH,"selpar_A50_rHB2");
  selpar_slope_rHB2.allocate(selpar_slope_rHB2_LO,selpar_slope_rHB2_HI,selpar_slope_rHB2_PH,"selpar_slope_rHB2");
  selpar_A502_rHB2.allocate(selpar_A502_rHB2_LO,selpar_A502_rHB2_HI,selpar_A502_rHB2_PH,"selpar_A502_rHB2");
  selpar_slope2_rHB2.allocate(selpar_slope2_rHB2_LO,selpar_slope2_rHB2_HI,selpar_slope2_rHB2_PH,"selpar_slope2_rHB2");
  selpar_A50_rHB3.allocate(selpar_A50_rHB3_LO,selpar_A50_rHB3_HI,selpar_A50_rHB3_PH,"selpar_A50_rHB3");
  selpar_slope_rHB3.allocate(selpar_slope_rHB3_LO,selpar_slope_rHB3_HI,selpar_slope_rHB3_PH,"selpar_slope_rHB3");
  selpar_A502_rHB3.allocate(selpar_A502_rHB3_LO,selpar_A502_rHB3_HI,selpar_A502_rHB3_PH,"selpar_A502_rHB3");
  selpar_slope2_rHB3.allocate(selpar_slope2_rHB3_LO,selpar_slope2_rHB3_HI,selpar_slope2_rHB3_PH,"selpar_slope2_rHB3");
  selpar_A50_rHB1_out.allocate(1,8,"selpar_A50_rHB1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rHB1_out.initialize();
  #endif
  selpar_slope_rHB1_out.allocate(1,8,"selpar_slope_rHB1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rHB1_out.initialize();
  #endif
  selpar_A502_rHB1_out.allocate(1,8,"selpar_A502_rHB1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_rHB1_out.initialize();
  #endif
  selpar_slope2_rHB1_out.allocate(1,8,"selpar_slope2_rHB1_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_rHB1_out.initialize();
  #endif
  selpar_A50_rHB2_out.allocate(1,8,"selpar_A50_rHB2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rHB2_out.initialize();
  #endif
  selpar_slope_rHB2_out.allocate(1,8,"selpar_slope_rHB2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rHB2_out.initialize();
  #endif
  selpar_A502_rHB2_out.allocate(1,8,"selpar_A502_rHB2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_rHB2_out.initialize();
  #endif
  selpar_slope2_rHB2_out.allocate(1,8,"selpar_slope2_rHB2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_rHB2_out.initialize();
  #endif
  selpar_A50_rHB3_out.allocate(1,8,"selpar_A50_rHB3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rHB3_out.initialize();
  #endif
  selpar_slope_rHB3_out.allocate(1,8,"selpar_slope_rHB3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rHB3_out.initialize();
  #endif
  selpar_A502_rHB3_out.allocate(1,8,"selpar_A502_rHB3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_rHB3_out.initialize();
  #endif
  selpar_slope2_rHB3_out.allocate(1,8,"selpar_slope2_rHB3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_rHB3_out.initialize();
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
  selpar_A50_rGN2.allocate(selpar_A50_rGN2_LO,selpar_A50_rGN2_HI,selpar_A50_rGN2_PH,"selpar_A50_rGN2");
  selpar_slope_rGN2.allocate(selpar_slope_rGN2_LO,selpar_slope_rGN2_HI,selpar_slope_rGN2_PH,"selpar_slope_rGN2");
  selpar_A502_rGN2.allocate(selpar_A502_rGN2_LO,selpar_A502_rGN2_HI,selpar_A502_rGN2_PH,"selpar_A502_rGN2");
  selpar_slope2_rGN2.allocate(selpar_slope2_rGN2_LO,selpar_slope2_rGN2_HI,selpar_slope2_rGN2_PH,"selpar_slope2_rGN2");
  selpar_A50_rGN3.allocate(selpar_A50_rGN3_LO,selpar_A50_rGN3_HI,selpar_A50_rGN3_PH,"selpar_A50_rGN3");
  selpar_slope_rGN3.allocate(selpar_slope_rGN3_LO,selpar_slope_rGN3_HI,selpar_slope_rGN3_PH,"selpar_slope_rGN3");
  selpar_A50_rGN2_out.allocate(1,8,"selpar_A50_rGN2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rGN2_out.initialize();
  #endif
  selpar_slope_rGN2_out.allocate(1,8,"selpar_slope_rGN2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rGN2_out.initialize();
  #endif
  selpar_A502_rGN2_out.allocate(1,8,"selpar_A502_rGN2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_rGN2_out.initialize();
  #endif
  selpar_slope2_rGN2_out.allocate(1,8,"selpar_slope2_rGN2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_rGN2_out.initialize();
  #endif
  selpar_A50_rGN3_out.allocate(1,8,"selpar_A50_rGN3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rGN3_out.initialize();
  #endif
  selpar_slope_rGN3_out.allocate(1,8,"selpar_slope_rGN3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rGN3_out.initialize();
  #endif
  sel_cHL_D.allocate(styr,endyr,1,nages,"sel_cHL_D");
  #ifndef NO_AD_INITIALIZE
    sel_cHL_D.initialize();
  #endif
  sel_rHB_D.allocate(styr,endyr,1,nages,"sel_rHB_D");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_D.initialize();
  #endif
  sel_rGN_D.allocate(styr,endyr,1,nages,"sel_rGN_D");
  #ifndef NO_AD_INITIALIZE
    sel_rGN_D.initialize();
  #endif
  sel_rHB_D_block3.allocate(1,nages,"sel_rHB_D_block3");
  #ifndef NO_AD_INITIALIZE
    sel_rHB_D_block3.initialize();
  #endif
  sel_rGN_D_block3.allocate(1,nages,"sel_rGN_D_block3");
  #ifndef NO_AD_INITIALIZE
    sel_rGN_D_block3.initialize();
  #endif
  selpar_A50_rHB2_D.allocate(selpar_A50_rHB2_D_LO,selpar_A50_rHB2_D_HI,selpar_A50_rHB2_D_PH,"selpar_A50_rHB2_D");
  selpar_slope_rHB2_D.allocate(selpar_slope_rHB2_D_LO,selpar_slope_rHB2_D_HI,selpar_slope_rHB2_D_PH,"selpar_slope_rHB2_D");
  selpar_A502_rHB2_D.allocate(selpar_A502_rHB2_D_LO,selpar_A502_rHB2_D_HI,selpar_A502_rHB2_D_PH,"selpar_A502_rHB2_D");
  selpar_slope2_rHB2_D.allocate(selpar_slope2_rHB2_D_LO,selpar_slope2_rHB2_D_HI,selpar_slope2_rHB2_D_PH,"selpar_slope2_rHB2_D");
  selpar_A50_rHB3_D.allocate(selpar_A50_rHB3_D_LO,selpar_A50_rHB3_D_HI,selpar_A50_rHB3_D_PH,"selpar_A50_rHB3_D");
  selpar_slope_rHB3_D.allocate(selpar_slope_rHB3_D_LO,selpar_slope_rHB3_D_HI,selpar_slope_rHB3_D_PH,"selpar_slope_rHB3_D");
  selpar_A502_rHB3_D.allocate(selpar_A502_rHB3_D_LO,selpar_A502_rHB3_D_HI,selpar_A502_rHB3_D_PH,"selpar_A502_rHB3_D");
  selpar_slope2_rHB3_D.allocate(selpar_slope2_rHB3_D_LO,selpar_slope2_rHB3_D_HI,selpar_slope2_rHB3_D_PH,"selpar_slope2_rHB3_D");
  selpar_A50_rGN3_D.allocate(selpar_A50_rGN3_D_LO,selpar_A50_rGN3_D_HI,selpar_A50_rGN3_D_PH,"selpar_A50_rGN3_D");
  selpar_slope_rGN3_D.allocate(selpar_slope_rGN3_D_LO,selpar_slope_rGN3_D_HI,selpar_slope_rGN3_D_PH,"selpar_slope_rGN3_D");
  selpar_A502_rGN3_D.allocate(selpar_A502_rGN3_D_LO,selpar_A502_rGN3_D_HI,selpar_A502_rGN3_D_PH,"selpar_A502_rGN3_D");
  selpar_slope2_rGN3_D.allocate(selpar_slope2_rGN3_D_LO,selpar_slope2_rGN3_D_HI,selpar_slope2_rGN3_D_PH,"selpar_slope2_rGN3_D");
  selpar_A50_cHL2_D.allocate(selpar_A50_cHL2_D_LO,selpar_A50_cHL2_D_HI,selpar_A50_cHL2_D_PH,"selpar_A50_cHL2_D");
  selpar_slope_cHL2_D.allocate(selpar_slope_cHL2_D_LO,selpar_slope_cHL2_D_HI,selpar_slope_cHL2_D_PH,"selpar_slope_cHL2_D");
  selpar_A502_cHL2_D.allocate(selpar_A502_cHL2_D_LO,selpar_A502_cHL2_D_HI,selpar_A502_cHL2_D_PH,"selpar_A502_cHL2_D");
  selpar_slope2_cHL2_D.allocate(selpar_slope2_cHL2_D_LO,selpar_slope2_cHL2_D_HI,selpar_slope2_cHL2_D_PH,"selpar_slope2_cHL2_D");
  selpar_A50_cHL3_D.allocate(selpar_A50_cHL3_D_LO,selpar_A50_cHL3_D_HI,selpar_A50_cHL3_D_PH,"selpar_A50_cHL3_D");
  selpar_slope_cHL3_D.allocate(selpar_slope_cHL3_D_LO,selpar_slope_cHL3_D_HI,selpar_slope_cHL3_D_PH,"selpar_slope_cHL3_D");
  selpar_A50_rHB2_D_out.allocate(1,8,"selpar_A50_rHB2_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rHB2_D_out.initialize();
  #endif
  selpar_slope_rHB2_D_out.allocate(1,8,"selpar_slope_rHB2_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rHB2_D_out.initialize();
  #endif
  selpar_A502_rHB2_D_out.allocate(1,8,"selpar_A502_rHB2_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_rHB2_D_out.initialize();
  #endif
  selpar_slope2_rHB2_D_out.allocate(1,8,"selpar_slope2_rHB2_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_rHB2_D_out.initialize();
  #endif
  selpar_A50_rHB3_D_out.allocate(1,8,"selpar_A50_rHB3_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rHB3_D_out.initialize();
  #endif
  selpar_slope_rHB3_D_out.allocate(1,8,"selpar_slope_rHB3_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rHB3_D_out.initialize();
  #endif
  selpar_A502_rHB3_D_out.allocate(1,8,"selpar_A502_rHB3_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_rHB3_D_out.initialize();
  #endif
  selpar_slope2_rHB3_D_out.allocate(1,8,"selpar_slope2_rHB3_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_rHB3_D_out.initialize();
  #endif
  selpar_A50_rGN3_D_out.allocate(1,8,"selpar_A50_rGN3_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_rGN3_D_out.initialize();
  #endif
  selpar_slope_rGN3_D_out.allocate(1,8,"selpar_slope_rGN3_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_rGN3_D_out.initialize();
  #endif
  selpar_A502_rGN3_D_out.allocate(1,8,"selpar_A502_rGN3_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_rGN3_D_out.initialize();
  #endif
  selpar_slope2_rGN3_D_out.allocate(1,8,"selpar_slope2_rGN3_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_rGN3_D_out.initialize();
  #endif
  selpar_A50_cHL2_D_out.allocate(1,8,"selpar_A50_cHL2_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cHL2_D_out.initialize();
  #endif
  selpar_slope_cHL2_D_out.allocate(1,8,"selpar_slope_cHL2_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cHL2_D_out.initialize();
  #endif
  selpar_A502_cHL2_D_out.allocate(1,8,"selpar_A502_cHL2_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_cHL2_D_out.initialize();
  #endif
  selpar_slope2_cHL2_D_out.allocate(1,8,"selpar_slope2_cHL2_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_cHL2_D_out.initialize();
  #endif
  selpar_A50_cHL3_D_out.allocate(1,8,"selpar_A50_cHL3_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cHL3_D_out.initialize();
  #endif
  selpar_slope_cHL3_D_out.allocate(1,8,"selpar_slope_cHL3_D_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cHL3_D_out.initialize();
  #endif
  sel_sCT.allocate(styr,endyr,1,nages,"sel_sCT");
  #ifndef NO_AD_INITIALIZE
    sel_sCT.initialize();
  #endif
  sel_sVD.allocate(styr,endyr,1,nages,"sel_sVD");
  #ifndef NO_AD_INITIALIZE
    sel_sVD.initialize();
  #endif
  sel_sCT_vec.allocate(1,nages,"sel_sCT_vec");
  #ifndef NO_AD_INITIALIZE
    sel_sCT_vec.initialize();
  #endif
  sel_sVD_vec.allocate(1,nages,"sel_sVD_vec");
  #ifndef NO_AD_INITIALIZE
    sel_sVD_vec.initialize();
  #endif
  selpar_A50_sCT.allocate(selpar_A50_sCT_LO,selpar_A50_sCT_HI,selpar_A50_sCT_PH,"selpar_A50_sCT");
  selpar_slope_sCT.allocate(selpar_slope_sCT_LO,selpar_slope_sCT_HI,selpar_slope_sCT_PH,"selpar_slope_sCT");
  selpar_A502_sCT.allocate(selpar_A502_sCT_LO,selpar_A502_sCT_HI,selpar_A502_sCT_PH,"selpar_A502_sCT");
  selpar_slope2_sCT.allocate(selpar_slope2_sCT_LO,selpar_slope2_sCT_HI,selpar_slope2_sCT_PH,"selpar_slope2_sCT");
  selpar_A50_sCT_out.allocate(1,8,"selpar_A50_sCT_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_sCT_out.initialize();
  #endif
  selpar_slope_sCT_out.allocate(1,8,"selpar_slope_sCT_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_sCT_out.initialize();
  #endif
  selpar_A502_sCT_out.allocate(1,8,"selpar_A502_sCT_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_sCT_out.initialize();
  #endif
  selpar_slope2_sCT_out.allocate(1,8,"selpar_slope2_sCT_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_sCT_out.initialize();
  #endif
  selpar_switch_sVD.allocate("selpar_switch_sVD");
  #ifndef NO_AD_INITIALIZE
  selpar_switch_sVD.initialize();
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
  pred_rHB_cpue.allocate(styr_cpue_rHB,endyr_cpue_rHB,"pred_rHB_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_rHB_cpue.initialize();
  #endif
  N_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,1,nages,"N_rHB");
  #ifndef NO_AD_INITIALIZE
    N_rHB.initialize();
  #endif
  pred_rHB_D_cpue.allocate(styr_D_cpue_rHB,endyr_D_cpue_rHB,"pred_rHB_D_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_rHB_D_cpue.initialize();
  #endif
  N_rHB_D.allocate(styr_D_cpue_rHB,endyr_D_cpue_rHB,1,nages,"N_rHB_D");
  #ifndef NO_AD_INITIALIZE
    N_rHB_D.initialize();
  #endif
  pred_sCT_cpue.allocate(styr_cpue_sCT,endyr_cpue_sCT,"pred_sCT_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_sCT_cpue.initialize();
  #endif
  N_sCT.allocate(styr_cpue_sCT,endyr_cpue_sCT,1,nages,"N_sCT");
  #ifndef NO_AD_INITIALIZE
    N_sCT.initialize();
  #endif
  pred_sVD_cpue.allocate(styr_cpue_sVD,endyr_cpue_sVD,"pred_sVD_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_sVD_cpue.initialize();
  #endif
  N_sVD.allocate(styr_cpue_sVD,endyr_cpue_sVD,1,nages,"N_sVD");
  #ifndef NO_AD_INITIALIZE
    N_sVD.initialize();
  #endif
  log_q_cpue_cHL.allocate(log_q_cHL_LO,log_q_cHL_HI,log_q_cHL_PH,"log_q_cpue_cHL");
  log_q_cpue_rHB.allocate(log_q_rHB_LO,log_q_rHB_HI,log_q_rHB_PH,"log_q_cpue_rHB");
  log_q_cpue_rHB_D.allocate(log_q_rHB_D_LO,log_q_rHB_D_HI,log_q_rHB_D_PH,"log_q_cpue_rHB_D");
  log_q_cpue_sCT.allocate(log_q_sCT_LO,log_q_sCT_HI,log_q_sCT_PH,"log_q_cpue_sCT");
  log_q_cpue_sVD.allocate(log_q_sVD_LO,log_q_sVD_HI,log_q_sVD_PH,"log_q_cpue_sVD");
  log_q_cHL_out.allocate(1,8,"log_q_cHL_out");
  #ifndef NO_AD_INITIALIZE
    log_q_cHL_out.initialize();
  #endif
  log_q_rHB_out.allocate(1,8,"log_q_rHB_out");
  #ifndef NO_AD_INITIALIZE
    log_q_rHB_out.initialize();
  #endif
  log_q_rHB_D_out.allocate(1,8,"log_q_rHB_D_out");
  #ifndef NO_AD_INITIALIZE
    log_q_rHB_D_out.initialize();
  #endif
  log_q_sCT_out.allocate(1,8,"log_q_sCT_out");
  #ifndef NO_AD_INITIALIZE
    log_q_sCT_out.initialize();
  #endif
  log_q_sVD_out.allocate(1,8,"log_q_sVD_out");
  #ifndef NO_AD_INITIALIZE
    log_q_sVD_out.initialize();
  #endif
  q_rate.allocate("q_rate");
  #ifndef NO_AD_INITIALIZE
  q_rate.initialize();
  #endif
  q_rate_fcn_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL,"q_rate_fcn_cHL");
  #ifndef NO_AD_INITIALIZE
    q_rate_fcn_cHL.initialize();
  #endif
  q_rate_fcn_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"q_rate_fcn_rHB");
  #ifndef NO_AD_INITIALIZE
    q_rate_fcn_rHB.initialize();
  #endif
  q_rate_fcn_rHB_D.allocate(styr_D_cpue_rHB,endyr_D_cpue_rHB,"q_rate_fcn_rHB_D");
  #ifndef NO_AD_INITIALIZE
    q_rate_fcn_rHB_D.initialize();
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
  q_RW_log_dev_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB-1,"q_RW_log_dev_rHB");
  #ifndef NO_AD_INITIALIZE
    q_RW_log_dev_rHB.initialize();
  #endif
  q_RW_log_dev_rHB_D.allocate(styr_D_cpue_rHB,endyr_D_cpue_rHB-1,"q_RW_log_dev_rHB_D");
  #ifndef NO_AD_INITIALIZE
    q_RW_log_dev_rHB_D.initialize();
  #endif
  q_RW_log_dev_sCT.allocate(styr_cpue_sCT,endyr_cpue_sCT-1,"q_RW_log_dev_sCT");
  #ifndef NO_AD_INITIALIZE
    q_RW_log_dev_sCT.initialize();
  #endif
  q_RW_log_dev_sVD.allocate(styr_cpue_sVD,endyr_cpue_sVD-1,"q_RW_log_dev_sVD");
  #ifndef NO_AD_INITIALIZE
    q_RW_log_dev_sVD.initialize();
  #endif
  q_cHL.allocate(styr_cpue_cHL,endyr_cpue_cHL,"q_cHL");
  #ifndef NO_AD_INITIALIZE
    q_cHL.initialize();
  #endif
  q_rHB.allocate(styr_cpue_rHB,endyr_cpue_rHB,"q_rHB");
  #ifndef NO_AD_INITIALIZE
    q_rHB.initialize();
  #endif
  q_rHB_D.allocate(styr_D_cpue_rHB,endyr_D_cpue_rHB,"q_rHB_D");
  #ifndef NO_AD_INITIALIZE
    q_rHB_D.initialize();
  #endif
  q_sCT.allocate(styr_cpue_sCT,endyr_cpue_sCT,"q_sCT");
  #ifndef NO_AD_INITIALIZE
    q_sCT.initialize();
  #endif
  q_sVD.allocate(styr_cpue_sVD,endyr_cpue_sVD,"q_sVD");
  #ifndef NO_AD_INITIALIZE
    q_sVD.initialize();
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
  D_cHL_num.allocate(styr,endyr,1,nages,"D_cHL_num");
  #ifndef NO_AD_INITIALIZE
    D_cHL_num.initialize();
  #endif
  D_cHL_klb.allocate(styr,endyr,1,nages,"D_cHL_klb");
  #ifndef NO_AD_INITIALIZE
    D_cHL_klb.initialize();
  #endif
  pred_cHL_D_knum.allocate(styr_D_cHL,endyr_D_cHL,"pred_cHL_D_knum");
  #ifndef NO_AD_INITIALIZE
    pred_cHL_D_knum.initialize();
  #endif
  obs_cHL_D.allocate(styr_D_cHL,endyr_D_cHL,"obs_cHL_D");
  #ifndef NO_AD_INITIALIZE
    obs_cHL_D.initialize();
  #endif
  pred_cHL_D_klb.allocate(styr_D_cHL,endyr_D_cHL,"pred_cHL_D_klb");
  #ifndef NO_AD_INITIALIZE
    pred_cHL_D_klb.initialize();
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
  Dmort_cHL1.allocate("Dmort_cHL1");
  #ifndef NO_AD_INITIALIZE
  Dmort_cHL1.initialize();
  #endif
  Dmort_rHB1.allocate("Dmort_rHB1");
  #ifndef NO_AD_INITIALIZE
  Dmort_rHB1.initialize();
  #endif
  Dmort_rGN1.allocate("Dmort_rGN1");
  #ifndef NO_AD_INITIALIZE
  Dmort_rGN1.initialize();
  #endif
  Dmort_cHL2.allocate("Dmort_cHL2");
  #ifndef NO_AD_INITIALIZE
  Dmort_cHL2.initialize();
  #endif
  Dmort_rHB2.allocate("Dmort_rHB2");
  #ifndef NO_AD_INITIALIZE
  Dmort_rHB2.initialize();
  #endif
  Dmort_rGN2.allocate("Dmort_rGN2");
  #ifndef NO_AD_INITIALIZE
  Dmort_rGN2.initialize();
  #endif
  Dmort_cHL3.allocate("Dmort_cHL3");
  #ifndef NO_AD_INITIALIZE
  Dmort_cHL3.initialize();
  #endif
  Dmort_rHB3.allocate("Dmort_rHB3");
  #ifndef NO_AD_INITIALIZE
  Dmort_rHB3.initialize();
  #endif
  Dmort_rGN3.allocate("Dmort_rGN3");
  #ifndef NO_AD_INITIALIZE
  Dmort_rGN3.initialize();
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
  F_cHL_D_prop.allocate("F_cHL_D_prop");
  #ifndef NO_AD_INITIALIZE
  F_cHL_D_prop.initialize();
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
  E_total_wgt.allocate(styr,endyr,"E_total_wgt");
  #ifndef NO_AD_INITIALIZE
    E_total_wgt.initialize();
  #endif
  E_L_wgt.allocate(styr,endyr,"E_L_wgt");
  #ifndef NO_AD_INITIALIZE
    E_L_wgt.initialize();
  #endif
  E_D_wgt.allocate(styr,endyr,"E_D_wgt");
  #ifndef NO_AD_INITIALIZE
    E_D_wgt.initialize();
  #endif
  E_total_num.allocate(styr,endyr,"E_total_num");
  #ifndef NO_AD_INITIALIZE
    E_total_num.initialize();
  #endif
  E_L_num.allocate(styr,endyr,"E_L_num");
  #ifndef NO_AD_INITIALIZE
    E_L_num.initialize();
  #endif
  E_D_num.allocate(styr,endyr,"E_D_num");
  #ifndef NO_AD_INITIALIZE
    E_D_num.initialize();
  #endif
  E_eq_wgt.allocate(1,n_iter_msy,"E_eq_wgt");
  #ifndef NO_AD_INITIALIZE
    E_eq_wgt.initialize();
  #endif
  E_eq_num.allocate(1,n_iter_msy,"E_eq_num");
  #ifndef NO_AD_INITIALIZE
    E_eq_num.initialize();
  #endif
  E_msy_wgt.allocate("E_msy_wgt");
  #ifndef NO_AD_INITIALIZE
  E_msy_wgt.initialize();
  #endif
  E_msy_num.allocate("E_msy_num");
  #ifndef NO_AD_INITIALIZE
  E_msy_num.initialize();
  #endif
  E_F30_wgt.allocate("E_F30_wgt");
  #ifndef NO_AD_INITIALIZE
  E_F30_wgt.initialize();
  #endif
  E_F30_num.allocate("E_F30_num");
  #ifndef NO_AD_INITIALIZE
  E_F30_num.initialize();
  #endif
  EdEmsy_wgt.allocate(styr,endyr,"EdEmsy_wgt");
  #ifndef NO_AD_INITIALIZE
    EdEmsy_wgt.initialize();
  #endif
  EdEmsy_num.allocate(styr,endyr,"EdEmsy_num");
  #ifndef NO_AD_INITIALIZE
    EdEmsy_num.initialize();
  #endif
  Eend_mean_temp_wgt.allocate("Eend_mean_temp_wgt");
  #ifndef NO_AD_INITIALIZE
  Eend_mean_temp_wgt.initialize();
  #endif
  E_wgt_end_mean.allocate("E_wgt_end_mean");
  #ifndef NO_AD_INITIALIZE
  E_wgt_end_mean.initialize();
  #endif
  Eend_mean_temp_num.allocate("Eend_mean_temp_num");
  #ifndef NO_AD_INITIALIZE
  Eend_mean_temp_num.initialize();
  #endif
  E_num_end_mean.allocate("E_num_end_mean");
  #ifndef NO_AD_INITIALIZE
  E_num_end_mean.initialize();
  #endif
  EdEmsy_wgt_end_mean.allocate("EdEmsy_wgt_end_mean");
  #ifndef NO_AD_INITIALIZE
  EdEmsy_wgt_end_mean.initialize();
  #endif
  EdEmsy_num_end_mean.allocate("EdEmsy_num_end_mean");
  #ifndef NO_AD_INITIALIZE
  EdEmsy_num_end_mean.initialize();
  #endif
  EdEF30_wgt.allocate(styr,endyr,"EdEF30_wgt");
  #ifndef NO_AD_INITIALIZE
    EdEF30_wgt.initialize();
  #endif
  EdEF30_num.allocate(styr,endyr,"EdEF30_num");
  #ifndef NO_AD_INITIALIZE
    EdEF30_num.initialize();
  #endif
  EdEF30_wgt_end_mean.allocate("EdEF30_wgt_end_mean");
  #ifndef NO_AD_INITIALIZE
  EdEF30_wgt_end_mean.initialize();
  #endif
  EdEF30_num_end_mean.allocate("EdEF30_num_end_mean");
  #ifndef NO_AD_INITIALIZE
  EdEF30_num_end_mean.initialize();
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
  log_avg_F_D_cHL.allocate(log_avg_F_cHL_D_LO,log_avg_F_cHL_D_HI,log_avg_F_cHL_D_PH,"log_avg_F_D_cHL");
  log_avg_F_cHL_D_out.allocate(1,8,"log_avg_F_cHL_D_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_cHL_D_out.initialize();
  #endif
  log_dev_F_D_cHL.allocate(styr_D_cHL,endyr_D_cHL,log_F_dev_cHL_D_LO,log_F_dev_cHL_D_HI,log_F_dev_cHL_D_PH,"log_dev_F_D_cHL");
  log_F_dev_cHL_D_out.allocate(styr_D_cHL,endyr_D_cHL,"log_F_dev_cHL_D_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_cHL_D_out.initialize();
  #endif
  F_cHL_D.allocate(styr,endyr,1,nages,"F_cHL_D");
  #ifndef NO_AD_INITIALIZE
    F_cHL_D.initialize();
  #endif
  F_cHL_D_out.allocate(styr,endyr,"F_cHL_D_out");
  #ifndef NO_AD_INITIALIZE
    F_cHL_D_out.initialize();
  #endif
  log_F_dev_end_cHL_D.allocate("log_F_dev_end_cHL_D");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_cHL_D.initialize();
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
  sdnr_lc_cHL_D.allocate("sdnr_lc_cHL_D");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_cHL_D.initialize();
  #endif
  sdnr_lc_rHB.allocate("sdnr_lc_rHB");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_rHB.initialize();
  #endif
  sdnr_lc_rHB_D.allocate("sdnr_lc_rHB_D");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_rHB_D.initialize();
  #endif
  sdnr_lc_sCT.allocate("sdnr_lc_sCT");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_sCT.initialize();
  #endif
  sdnr_lc_rGN.allocate("sdnr_lc_rGN");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_rGN.initialize();
  #endif
  sdnr_lc_rGN_D.allocate("sdnr_lc_rGN_D");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_rGN_D.initialize();
  #endif
  sdnr_ac_cHL.allocate("sdnr_ac_cHL");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_cHL.initialize();
  #endif
  sdnr_ac_rHB.allocate("sdnr_ac_rHB");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_rHB.initialize();
  #endif
  sdnr_ac_sCT.allocate("sdnr_ac_sCT");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_sCT.initialize();
  #endif
  sdnr_ac_rGN.allocate("sdnr_ac_rGN");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_rGN.initialize();
  #endif
  sdnr_I_cHL.allocate("sdnr_I_cHL");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_cHL.initialize();
  #endif
  sdnr_I_rHB.allocate("sdnr_I_rHB");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_rHB.initialize();
  #endif
  sdnr_I_rHB_D.allocate("sdnr_I_rHB_D");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_rHB_D.initialize();
  #endif
  sdnr_I_sCT.allocate("sdnr_I_sCT");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_sCT.initialize();
  #endif
  sdnr_I_sVD.allocate("sdnr_I_sVD");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_sVD.initialize();
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
  w_cpue_rHB.allocate("w_cpue_rHB");
  #ifndef NO_AD_INITIALIZE
  w_cpue_rHB.initialize();
  #endif
  w_cpue_rHB_D.allocate("w_cpue_rHB_D");
  #ifndef NO_AD_INITIALIZE
  w_cpue_rHB_D.initialize();
  #endif
  w_cpue_sCT.allocate("w_cpue_sCT");
  #ifndef NO_AD_INITIALIZE
  w_cpue_sCT.initialize();
  #endif
  w_cpue_sVD.allocate("w_cpue_sVD");
  #ifndef NO_AD_INITIALIZE
  w_cpue_sVD.initialize();
  #endif
  w_cpue_sCT_mult.allocate("w_cpue_sCT_mult");
  #ifndef NO_AD_INITIALIZE
  w_cpue_sCT_mult.initialize();
  #endif
  w_cpue_sVD_mult.allocate("w_cpue_sVD_mult");
  #ifndef NO_AD_INITIALIZE
  w_cpue_sVD_mult.initialize();
  #endif
  w_lenc_cHL.allocate("w_lenc_cHL");
  #ifndef NO_AD_INITIALIZE
  w_lenc_cHL.initialize();
  #endif
  w_lenc_cHL_D.allocate("w_lenc_cHL_D");
  #ifndef NO_AD_INITIALIZE
  w_lenc_cHL_D.initialize();
  #endif
  w_lenc_rHB_D.allocate("w_lenc_rHB_D");
  #ifndef NO_AD_INITIALIZE
  w_lenc_rHB_D.initialize();
  #endif
  w_lenc_rGN_D.allocate("w_lenc_rGN_D");
  #ifndef NO_AD_INITIALIZE
  w_lenc_rGN_D.initialize();
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
  f_cHL_D.allocate("f_cHL_D");
  #ifndef NO_AD_INITIALIZE
  f_cHL_D.initialize();
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
  f_rHB_cpue.allocate("f_rHB_cpue");
  #ifndef NO_AD_INITIALIZE
  f_rHB_cpue.initialize();
  #endif
  f_rHB_D_cpue.allocate("f_rHB_D_cpue");
  #ifndef NO_AD_INITIALIZE
  f_rHB_D_cpue.initialize();
  #endif
  f_sCT_cpue.allocate("f_sCT_cpue");
  #ifndef NO_AD_INITIALIZE
  f_sCT_cpue.initialize();
  #endif
  f_sVD_cpue.allocate("f_sVD_cpue");
  #ifndef NO_AD_INITIALIZE
  f_sVD_cpue.initialize();
  #endif
  f_cHL_lenc.allocate("f_cHL_lenc");
  #ifndef NO_AD_INITIALIZE
  f_cHL_lenc.initialize();
  #endif
  f_cHL_D_lenc.allocate("f_cHL_D_lenc");
  #ifndef NO_AD_INITIALIZE
  f_cHL_D_lenc.initialize();
  #endif
  f_rHB_lenc.allocate("f_rHB_lenc");
  #ifndef NO_AD_INITIALIZE
  f_rHB_lenc.initialize();
  #endif
  f_rHB_D_lenc.allocate("f_rHB_D_lenc");
  #ifndef NO_AD_INITIALIZE
  f_rHB_D_lenc.initialize();
  #endif
  f_sCT_lenc.allocate("f_sCT_lenc");
  #ifndef NO_AD_INITIALIZE
  f_sCT_lenc.initialize();
  #endif
  f_rGN_lenc.allocate("f_rGN_lenc");
  #ifndef NO_AD_INITIALIZE
  f_rGN_lenc.initialize();
  #endif
  f_rGN_D_lenc.allocate("f_rGN_D_lenc");
  #ifndef NO_AD_INITIALIZE
  f_rGN_D_lenc.initialize();
  #endif
  f_cHL_agec.allocate("f_cHL_agec");
  #ifndef NO_AD_INITIALIZE
  f_cHL_agec.initialize();
  #endif
  f_rHB_agec.allocate("f_rHB_agec");
  #ifndef NO_AD_INITIALIZE
  f_rHB_agec.initialize();
  #endif
  f_sCT_agec.allocate("f_sCT_agec");
  #ifndef NO_AD_INITIALIZE
  f_sCT_agec.initialize();
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
}

void model_parameters::set_runtime(void)
{
  dvector temp1("{1000, 2000, 3000, 5000; 5000; 5000; 10000; }");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-2, 1e-2, 1e-3, 1e-3; 1e-4; 1e-4; 1e-4; }");
  convergence_criteria.allocate(temp.indexmin(),temp.indexmax());
  convergence_criteria=temp;
}

void model_parameters::preliminary_calculations(void)
{

#if defined(USE_ADPVM)

  admaster_slave_variable_interface(*this);

#endif
  Dmort_cHL1=set_Dmort_cHL1; Dmort_rHB1=set_Dmort_rHB1; Dmort_rGN1=set_Dmort_rGN1;
  Dmort_cHL2=set_Dmort_cHL2; Dmort_rHB2=set_Dmort_rHB2; Dmort_rGN2=set_Dmort_rGN2;
  Dmort_cHL3=set_Dmort_cHL3; Dmort_rHB3=set_Dmort_rHB3; Dmort_rGN3=set_Dmort_rGN3;
  
  for(iyear=styr_D_cHL; iyear<=endyr_cHL_Dmort1; iyear++)
	{obs_cHL_D(iyear)=Dmort_cHL1*obs_released_cHL(iyear);
	}
  for(iyear=endyr_cHL_Dmort1+1; iyear<=endyr_cHL_Dmort2; iyear++)
	{obs_cHL_D(iyear)=Dmort_cHL2*obs_released_cHL(iyear);
	}
  for(iyear=endyr_cHL_Dmort2+1; iyear<=endyr_D_cHL; iyear++)
	{obs_cHL_D(iyear)=Dmort_cHL3*obs_released_cHL(iyear);
	}
	
  for(iyear=styr_D_rHB; iyear<=endyr_rHB_Dmort1; iyear++)
	{obs_rHB_D(iyear)=Dmort_rHB1*obs_released_rHB(iyear);
	}
  for(iyear=endyr_rHB_Dmort1+1; iyear<=endyr_rHB_Dmort2; iyear++)
	{obs_rHB_D(iyear)=Dmort_rHB2*obs_released_rHB(iyear);
	}
  for(iyear=endyr_rHB_Dmort2+1; iyear<=endyr_D_rHB; iyear++)
	{obs_rHB_D(iyear)=Dmort_rHB3*obs_released_rHB(iyear);
	}
	
  for(iyear=styr_D_rGN; iyear<=endyr_rGN_Dmort1; iyear++)
	{obs_rGN_D(iyear)=Dmort_rGN1*obs_released_rGN(iyear);
	}
  for(iyear=endyr_rGN_Dmort1+1; iyear<=endyr_rGN_Dmort2; iyear++)
	{obs_rGN_D(iyear)=Dmort_rGN2*obs_released_rGN(iyear);
	}
  for(iyear=endyr_rGN_Dmort2+1; iyear<=endyr_D_rGN; iyear++)
	{obs_rGN_D(iyear)=Dmort_rGN3*obs_released_rGN(iyear);
	}
	
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
  //20" size limit
  Linf_20=set_Linf_20(1);  
  K_20=set_K_20(1);
  t0_20=set_t0_20(1);
  len_cv_val_20=set_len_cv_20(1);
  
  M=set_M; 
  M_constant=set_M_constant(1);
  smsy2msst=1.0-M_constant;
  smsy2msst75=0.75;  
  
  log_R0=set_log_R0(1);
  steep=set_steep(1);
  R_autocorr=set_R_autocorr(1);
  rec_sigma=set_rec_sigma(1);
  
  log_dm_lenc_cHL=set_log_dm_lenc_cHL(1);
  log_dm_lenc_cHL_D=set_log_dm_lenc_cHL_D(1);
  log_dm_lenc_rHB_D=set_log_dm_lenc_rHB_D(1);
  log_dm_lenc_rGN_D=set_log_dm_lenc_rGN_D(1);
  log_dm_agec_cHL=set_log_dm_agec_cHL(1);
  log_dm_agec_rHB=set_log_dm_agec_rHB(1);
  log_dm_agec_sCT=set_log_dm_agec_sCT(1);
  log_dm_agec_rGN=set_log_dm_agec_rGN(1);  
  
  log_q_cpue_cHL=set_log_q_cpue_cHL(1);
  log_q_cpue_rHB=set_log_q_cpue_rHB(1);
  log_q_cpue_rHB_D=set_log_q_cpue_rHB_D(1);
  log_q_cpue_sCT=set_log_q_cpue_sCT(1);
  log_q_cpue_sVD=set_log_q_cpue_sVD(1);
  
  q_rate=set_q_rate;
  q_rate_fcn_cHL=1.0;   
  q_rate_fcn_rHB=1.0;   
  q_rate_fcn_rHB_D=1.0; 
  q_DD_beta=set_q_DD_beta;
  q_DD_fcn=1.0;
  q_RW_log_dev_cHL.initialize(); 
  q_RW_log_dev_rHB.initialize(); 
  q_RW_log_dev_rHB_D.initialize();
  q_RW_log_dev_sCT.initialize();
  q_RW_log_dev_sVD.initialize();
  
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
    for (iyear=styr_D_cpue_rHB; iyear<=endyr_D_cpue_rHB; iyear++)
      {   if (iyear>styr_D_cpue_rHB & iyear <=2003) 
          {//q_rate_fcn_rHB_D(iyear)=(1.0+q_rate)*q_rate_fcn_rHB_D(iyear-1); //compound
             q_rate_fcn_rHB_D(iyear)=(1.0+(iyear-styr_D_cpue_rHB)*q_rate)*q_rate_fcn_rHB_D(styr_D_cpue_rHB);  //linear
          }
          if (iyear>2003) {q_rate_fcn_rHB_D(iyear)=q_rate_fcn_rHB_D(iyear-1);} 
      }   
  } //end q_rate conditional      
  w_L=set_w_L;
  w_D=set_w_D;
  
  w_cpue_cHL=set_w_cpue_cHL;
  w_cpue_rHB=set_w_cpue_rHB;
  w_cpue_rHB_D=set_w_cpue_rHB_D;
  w_cpue_sCT=set_w_cpue_sCT;
  w_cpue_sVD=set_w_cpue_sVD;
  w_cpue_sCT_mult=set_w_cpue_sCT_mult;
  w_cpue_sVD_mult=set_w_cpue_sVD_mult;
  
  w_lenc_cHL=set_w_lenc_cHL;
  w_lenc_cHL_D=set_w_lenc_cHL_D;  
  w_lenc_rHB_D=set_w_lenc_rHB_D;
  w_lenc_rGN_D=set_w_lenc_rGN_D;  
    
  w_agec_cHL=set_w_agec_cHL;
  w_agec_rHB=set_w_agec_rHB;
  w_agec_sCT=set_w_agec_sCT;  
  w_agec_rGN=set_w_agec_rGN;  
  
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
  log_avg_F_D_cHL=set_log_avg_F_D_cHL(1);
  log_avg_F_D_rHB=set_log_avg_F_D_rHB(1); 
  log_avg_F_D_rGN=set_log_avg_F_D_rGN(1); 
    
  log_dev_F_L_cHL=set_log_dev_vals_F_L__cHL;
  log_dev_F_L_rHB=set_log_dev_vals_F_L__rHB;
  log_dev_F_L_rGN=set_log_dev_vals_F_L__rGN;
  log_dev_F_D_cHL=set_log_dev_vals_F_D__cHvals;
  log_dev_F_D_rHB=set_log_dev_vals_F_D__HBvals;
  log_dev_F_D_rGN=set_log_dev_vals_F_D__GRvals;
 
  selpar_A50_cHL1=set_selpar_A50_cHL1(1);
  selpar_slope_cHL1=set_selpar_slope_cHL1(1);
  selpar_A50_cHL2=set_selpar_A50_cHL2(1);
  selpar_slope_cHL2=set_selpar_slope_cHL2(1);
  selpar_A50_cHL3=set_selpar_A50_cHL3(1);
  selpar_slope_cHL3=set_selpar_slope_cHL3(1);
  selpar_A50_rHB1=set_selpar_A50_rHB1(1);
  selpar_slope_rHB1=set_selpar_slope_rHB1(1);
  selpar_A502_rHB1=set_selpar_A502_rHB1(1);
  selpar_slope2_rHB1=set_selpar_slope2_rHB1(1);
  selpar_A50_rHB2=set_selpar_A50_rHB2(1);
  selpar_slope_rHB2=set_selpar_slope_rHB2(1);
  selpar_A502_rHB2=set_selpar_A502_rHB2(1);
  selpar_slope2_rHB2=set_selpar_slope2_rHB2(1);
  selpar_A50_rHB3=set_selpar_A50_rHB3(1);
  selpar_slope_rHB3=set_selpar_slope_rHB3(1);
  selpar_A502_rHB3=set_selpar_A502_rHB3(1);
  selpar_slope2_rHB3=set_selpar_slope2_rHB3(1);
  
  selpar_A50_rGN2=set_selpar_A50_rGN2(1);
  selpar_slope_rGN2=set_selpar_slope_rGN2(1);
  selpar_A502_rGN2=set_selpar_A502_rGN2(1);
  selpar_slope2_rGN2=set_selpar_slope2_rGN2(1);
  selpar_A50_rGN3=set_selpar_A50_rGN3(1);
  selpar_slope_rGN3=set_selpar_slope_rGN3(1);
  
  selpar_A50_rHB2_D=set_selpar_A50_rHB2_D(1);
  selpar_slope_rHB2_D=set_selpar_slope_rHB2_D(1);
  selpar_A502_rHB2_D=set_selpar_A502_rHB2_D(1);
  selpar_slope2_rHB2_D=set_selpar_slope2_rHB2_D(1);
  
  selpar_A50_rHB3_D=set_selpar_A50_rHB3_D(1);
  selpar_slope_rHB3_D=set_selpar_slope_rHB3_D(1);
  selpar_A502_rHB3_D=set_selpar_A502_rHB3_D(1);
  selpar_slope2_rHB3_D=set_selpar_slope2_rHB3_D(1);
  selpar_A50_rGN3_D=set_selpar_A50_rGN3_D(1);
  selpar_slope_rGN3_D=set_selpar_slope_rGN3_D(1);
  selpar_A502_rGN3_D=set_selpar_A502_rGN3_D(1);
  selpar_slope2_rGN3_D=set_selpar_slope2_rGN3_D(1);  
 
  selpar_A50_cHL2_D=set_selpar_A50_cHL2_D(1);
  selpar_slope_cHL2_D=set_selpar_slope_cHL2_D(1);
  selpar_A502_cHL2_D=set_selpar_A502_cHL2_D(1);
  selpar_slope2_cHL2_D=set_selpar_slope2_cHL2_D(1);
  
  selpar_A50_cHL3_D=set_selpar_A50_cHL3_D(1);
  selpar_slope_cHL3_D=set_selpar_slope_cHL3_D(1);
    
  selpar_A50_sCT=set_selpar_A50_sCT(1);
  selpar_slope_sCT=set_selpar_slope_sCT(1);
  selpar_A502_sCT=set_selpar_A502_sCT(1);
  selpar_slope2_sCT=set_selpar_slope2_sCT(1);
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
 fecpar_a=set_fecpar_a;
 fecpar_b=set_fecpar_b;
 fecpar_c=set_fecpar_c;
 fecpar_thresh=set_fecpar_thresh;
 fecpar_min=set_fecpar_min;
 fecpar_batches=set_fecpar_batches;
 
      nsamp_cHL_lenc_allyr=missing;
	  nsamp_cHL_D_lenc_allyr=missing;  
      nsamp_rHB_D_lenc_allyr=missing;  
	  nsamp_rGN_D_lenc_allyr=missing; 
      nsamp_cHL_agec_allyr=missing;
      nsamp_rHB_agec_allyr=missing;
	  nsamp_sCT_agec_allyr=missing; 
	  nsamp_rGN_agec_allyr=missing;  
      
      nfish_cHL_lenc_allyr=missing;
	  nfish_cHL_D_lenc_allyr=missing;  
      nfish_rHB_D_lenc_allyr=missing; 
	  nfish_rGN_D_lenc_allyr=missing;
      nfish_cHL_agec_allyr=missing;
      nfish_rHB_agec_allyr=missing;
	  nfish_sCT_agec_allyr=missing;  
	  nfish_rGN_agec_allyr=missing;  
   
      for (iyear=1; iyear<=nyr_lenc_cHL; iyear++)
         {if (nsamp_lenc_cHL(iyear)>=minSS_lenc_cHL)
           {nsamp_cHL_lenc_allyr(yrs_lenc_cHL(iyear))=nsamp_lenc_cHL(iyear);
            nfish_cHL_lenc_allyr(yrs_lenc_cHL(iyear))=nfish_lenc_cHL(iyear);}}
	  for (iyear=1; iyear<=nyr_lenc_cHL_D; iyear++)                           
         {if (nsamp_lenc_cHL_D(iyear)>=minSS_lenc_cHL_D)
            {nsamp_cHL_D_lenc_allyr(yrs_lenc_cHL_D(iyear))=nsamp_lenc_cHL_D(iyear);
             nfish_cHL_D_lenc_allyr(yrs_lenc_cHL_D(iyear))=nfish_lenc_cHL_D(iyear);}}
      for (iyear=1; iyear<=nyr_lenc_rHB_D; iyear++)                          
         {if (nsamp_lenc_rHB_D(iyear)>=minSS_lenc_rHB_D)
            {nsamp_rHB_D_lenc_allyr(yrs_lenc_rHB_D(iyear))=nsamp_lenc_rHB_D(iyear);
             nfish_rHB_D_lenc_allyr(yrs_lenc_rHB_D(iyear))=nfish_lenc_rHB_D(iyear);}}
	  for (iyear=1; iyear<=nyr_lenc_rGN_D; iyear++)                          
         {if (nsamp_lenc_rGN_D(iyear)>=minSS_lenc_rGN_D)
            {nsamp_rGN_D_lenc_allyr(yrs_lenc_rGN_D(iyear))=nsamp_lenc_rGN_D(iyear);
             nfish_rGN_D_lenc_allyr(yrs_lenc_rGN_D(iyear))=nfish_lenc_rGN_D(iyear);}}
	  for (iyear=1; iyear<=nyr_agec_cHL; iyear++)
         {if (nsamp_agec_cHL(iyear)>=minSS_agec_cHL)
           {nsamp_cHL_agec_allyr(yrs_agec_cHL(iyear))=nsamp_agec_cHL(iyear);
            nfish_cHL_agec_allyr(yrs_agec_cHL(iyear))=nfish_agec_cHL(iyear);}}
      for (iyear=1; iyear<=nyr_agec_rHB; iyear++)
          {if (nsamp_agec_rHB(iyear)>=minSS_agec_rHB)
            {nsamp_rHB_agec_allyr(yrs_agec_rHB(iyear))=nsamp_agec_rHB(iyear);
             nfish_rHB_agec_allyr(yrs_agec_rHB(iyear))=nfish_agec_rHB(iyear);}} 
	  for (iyear=1; iyear<=nyr_agec_sCT; iyear++)  
          {if (nsamp_agec_sCT(iyear)>=minSS_agec_sCT)
            {nsamp_sCT_agec_allyr(yrs_agec_sCT(iyear))=nsamp_agec_sCT(iyear);
             nfish_sCT_agec_allyr(yrs_agec_sCT(iyear))=nfish_agec_sCT(iyear);}} 
      for (iyear=1; iyear<=nyr_agec_rGN; iyear++)  
         {if (nsamp_agec_rGN(iyear)>=minSS_agec_rGN)
           {nsamp_rGN_agec_allyr(yrs_agec_rGN(iyear))=nsamp_agec_rGN(iyear);
             nfish_rGN_agec_allyr(yrs_agec_rGN(iyear))=nfish_agec_rGN(iyear);}}  
             
  F_msy(1)=0.0;  
  for (ff=2;ff<=n_iter_msy;ff++) {F_msy(ff)=F_msy(ff-1)+iter_inc_msy;}
  F_spr(1)=0.0;  
  for (ff=2;ff<=n_iter_spr;ff++) {F_spr(ff)=F_spr(ff-1)+iter_inc_spr;}
  F_cHL.initialize(); L_cHL_num.initialize();
  F_rHB.initialize(); L_rHB_num.initialize();
  F_rGN.initialize(); L_rGN_num.initialize();
  F_cHL_D.initialize(); D_cHL_num.initialize();
  F_rHB_D.initialize(); D_rHB_num.initialize();
  F_rGN_D.initialize(); D_rGN_num.initialize();
  F_cHL_out.initialize();
  F_rHB_out.initialize();
  F_rGN_out.initialize();
  F_cHL_D_out.initialize();
  F_rHB_D_out.initialize();
  F_rGN_D_out.initialize();
      
  sel_cHL.initialize();
  sel_rHB.initialize();
  sel_rGN.initialize();
  sel_cHL_D.initialize();
  sel_rHB_D.initialize();
  sel_rGN_D.initialize();
  sel_sCT.initialize();
  sel_sVD.initialize();
  
  sel_rHB_block1.initialize();
  sel_rHB_block2.initialize();
  sel_rHB_block3.initialize();
  sel_rGN_block1.initialize();
  sel_rGN_block2.initialize();
  sel_rGN_block3.initialize();
  sel_rHB_D_block3.initialize();
  sel_rGN_D_block3.initialize();  
  sel_sCT_vec.initialize();
  sel_sVD_vec.initialize();
 
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
  //cout << "got repro stuff" << endl;
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
  //cout<< "got length comps"<< endl;
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
    fecundity=fecpar_a+fecpar_b*pow(meanlen_TL,fecpar_c); //relationship provided by Wyanski et al.  
    for (iage=1; iage<=nages; iage++)
		{ if(meanlen_TL(iage)<fecpar_thresh){fecundity(iage)=fecpar_min;}	
          else {break;}			
		}
	fecundity=fecundity/fecpar_scale; //batch fecundity of a mature female at age in units of fecpar_scale
	//20 in size limit
	meanlen_TL_20=Linf_20*(1.0-mfexp(-K_20*(agebins-t0_20+0.5)));    //Landings total length in mm
    wgt_kg_20=wgtpar_a*pow(meanlen_TL_20,wgtpar_b);             //whole wgt in kg 
    wgt_g_20=wgt_kg_20/g2kg;                                    //convert wgt in kg to weight in g    
    wgt_mt_20=wgt_g_20*g2mt;                                    //convert weight in g to weight in mt
    wgt_klb_20=mt2klb*wgt_mt_20;                               //1000 lb of whole wgt
    wgt_lb_20=mt2lb*wgt_mt_20;                                 //1000 lb of whole wgt
	//All fisheries
	meanlen_TL_L=Linf_L*(1.0-mfexp(-K_L*(agebins-t0_L+0.5)));    //Landings total length in mm
    wgt_kg_L=wgtpar_a*pow(meanlen_TL_L,wgtpar_b);             //whole wgt in kg 
    wgt_g_L=wgt_kg_L/g2kg;                                    //convert wgt in kg to weight in g    
    wgt_mt_L=wgt_g_L*g2mt;                                    //convert weight in g to weight in mt
    wgt_klb_L=mt2klb*wgt_mt_L;                               //1000 lb of whole wgt
    wgt_lb_L=mt2lb*wgt_mt_L;                                 //1000 lb of whole wgt
}

void model_parameters::get_reprod(void)
{
	reprod=elem_prod(elem_prod(elem_prod(prop_f,maturity_f),fecundity),fecpar_batches);
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
	//20 in size limit
	//len_cv_20(iage)=mfexp(log_len_cv_20+log_len_cv_dev_20(iage));
    len_cv_20(iage)=len_cv_val_20;
    len_sd_20(iage)=meanlen_TL_20(iage)*len_cv_20(iage);
	zscore_lzero_20=(0.0-meanlen_TL_20(iage))/len_sd_20(iage); 
	cprob_lzero_20=cumd_norm(zscore_lzero_20);
	//All fishery dependent
	//len_cv_L(iage)=mfexp(log_len_cv_L+log_len_cv_dev_L(iage));
    len_cv_L(iage)=len_cv_val_L;
    len_sd_L(iage)=meanlen_TL_L(iage)*len_cv_L(iage);
	zscore_lzero_L=(0.0-meanlen_TL_L(iage))/len_sd_L(iage); 
	cprob_lzero_L=cumd_norm(zscore_lzero_L);
    //first length bin
	//population
     zscore_len=((lenbins(1)+0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage);
    cprob_lenvec(1)=cumd_norm(zscore_len);          //includes any probability mass below zero
    lenprob(iage,1)=cprob_lenvec(1)-cprob_lzero;    //removes any probability mass below zero
	//20 in size limit
	 zscore_len_20=((lenbins(1)+0.5*lenbins_width)-meanlen_TL_20(iage)) / len_sd_20(iage);
    cprob_lenvec_20(1)=cumd_norm(zscore_len_20);          //includes any probability mass below zero
    lenprob_20(iage,1)=cprob_lenvec_20(1)-cprob_lzero_20;    //removes any probability mass below zero
	//All fishery dependent
    zscore_len_L=((lenbins(1)+0.5*lenbins_width)-meanlen_TL_L(iage)) / len_sd_L(iage);
    cprob_lenvec_L(1)=cumd_norm(zscore_len_L);          //includes any probability mass below zero
    lenprob_L(iage,1)=cprob_lenvec_L(1)-cprob_lzero_L;    //removes any probability mass below zero
    //most other length bins  
    //population
    for (ilen=2;ilen<nlenbins;ilen++)
      {
        zscore_len=((lenbins(ilen)+0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage); 
		cprob_lenvec(ilen)=cumd_norm(zscore_len);
        lenprob(iage,ilen)=cprob_lenvec(ilen)-cprob_lenvec(ilen-1);
      }
	//20 in size limit
	for (ilen=2;ilen<nlenbins;ilen++)
      {
        zscore_len_20=((lenbins(ilen)+0.5*lenbins_width)-meanlen_TL_20(iage)) / len_sd_20(iage); 
		cprob_lenvec_20(ilen)=cumd_norm(zscore_len_20);
        lenprob_20(iage,ilen)=cprob_lenvec_20(ilen)-cprob_lenvec_20(ilen-1);
      }
	//All fishery dependent	
    for (ilen=2;ilen<nlenbins;ilen++)
      {
        zscore_len_L=((lenbins(ilen)+0.5*lenbins_width)-meanlen_TL_L(iage)) / len_sd_L(iage); 
		cprob_lenvec_L(ilen)=cumd_norm(zscore_len_L);
        lenprob_L(iage,ilen)=cprob_lenvec_L(ilen)-cprob_lenvec_L(ilen-1);
      }
    //last length bin is a plus group
	//population
    zscore_len=((lenbins(nlenbins)-0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage); 
	lenprob(iage,nlenbins)=1.0-cumd_norm(zscore_len);
    lenprob(iage)=lenprob(iage)/(1.0-cprob_lzero);  //renormalize to account for any prob mass below size=0
	//20 in size limit
	zscore_len_20=((lenbins(nlenbins)-0.5*lenbins_width)-meanlen_TL_20(iage)) / len_sd_20(iage); 
	lenprob_20(iage,nlenbins)=1.0-cumd_norm(zscore_len_20);
    lenprob_20(iage)=lenprob_20(iage)/(1.0-cprob_lzero_20);  //renormalize to account for any prob mass below size=0
	//All fishery dependent
    zscore_len_L=((lenbins(nlenbins)-0.5*lenbins_width)-meanlen_TL_L(iage)) / len_sd_L(iage); 
	lenprob_L(iage,nlenbins)=1.0-cumd_norm(zscore_len_L);
    lenprob_L(iage)=lenprob_L(iage)/(1.0-cprob_lzero_L);  //renormalize to account for any prob mass below size=0
  }
  //fleet and survey specific length probs, all assumed here to equal the popn or to the landings
  lenprob_cHL=lenprob_L;
  lenprob_cHL_D=lenprob; //population
  lenprob_rHB=lenprob_L;
  lenprob_rHB_D=lenprob; //population
  lenprob_rGN_D=lenprob; //population
  lenprob_sCT=lenprob; //population
  lenprob_rGN=lenprob_L;  
}

void model_parameters::get_weight_at_age_landings(void)
{
  //for (iyear=styr; iyear<=endyr; iyear++)
  for (iyear=styr; iyear<=endyr_selex_phase2; iyear++) //start through 1991
  {
    len_cHL_mm(iyear)=meanlen_TL_L;  //using all fishery dependent growth for fisheries until 20 in size limit
    wholewgt_cHL_klb(iyear)=wgt_klb_L; 
    len_rHB_mm(iyear)=meanlen_TL_L;
    wholewgt_rHB_klb(iyear)=wgt_klb_L;
    len_rGN_mm(iyear)=meanlen_TL_L;
    wholewgt_rGN_klb(iyear)=wgt_klb_L;
    len_cHL_D_mm(iyear)=meanlen_TL;  //using the population growth to define discards
    wholewgt_cHL_D_klb(iyear)=wgt_klb;
    len_rHB_D_mm(iyear)=meanlen_TL;
    wholewgt_rHB_D_klb(iyear)=wgt_klb;
    len_rGN_D_mm(iyear)=meanlen_TL;
    wholewgt_rGN_D_klb(iyear)=wgt_klb;
  }  
  for (iyear=endyr_selex_phase2+1; iyear<=endyr_selex_phase3; iyear++)  //1992 through 2009
  {
    len_cHL_mm(iyear)=meanlen_TL_20;
    wholewgt_cHL_klb(iyear)=wgt_klb_20; 
    len_rHB_mm(iyear)=meanlen_TL_20;
    wholewgt_rHB_klb(iyear)=wgt_klb_20;
    len_rGN_mm(iyear)=meanlen_TL_20;
    wholewgt_rGN_klb(iyear)=wgt_klb_20;
    len_cHL_D_mm(iyear)=meanlen_TL;
    wholewgt_cHL_D_klb(iyear)=wgt_klb;   //using population growth curve for the discard weight
    len_rHB_D_mm(iyear)=meanlen_TL;
    wholewgt_rHB_D_klb(iyear)=wgt_klb;
    len_rGN_D_mm(iyear)=meanlen_TL;
    wholewgt_rGN_D_klb(iyear)=wgt_klb;
  }  
  for (iyear=endyr_selex_phase3+1; iyear<=endyr; iyear++)  //2010 through 2014 save in case we end up treating as a different growth curve
  {
    len_cHL_mm(iyear)=meanlen_TL_L;
    wholewgt_cHL_klb(iyear)=wgt_klb_L; 
    len_rHB_mm(iyear)=meanlen_TL_L;
    wholewgt_rHB_klb(iyear)=wgt_klb_L;
    len_rGN_mm(iyear)=meanlen_TL_L;
    wholewgt_rGN_klb(iyear)=wgt_klb_L;
    len_cHL_D_mm(iyear)=meanlen_TL;
    wholewgt_cHL_D_klb(iyear)=wgt_klb;
    len_rHB_D_mm(iyear)=meanlen_TL;
    wholewgt_rHB_D_klb(iyear)=wgt_klb;
    len_rGN_D_mm(iyear)=meanlen_TL;
    wholewgt_rGN_D_klb(iyear)=wgt_klb;
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
  sel_rHB_block1=logistic_double(agebins, selpar_A50_rHB1, selpar_slope_rHB1, selpar_A502_rHB1, selpar_slope2_rHB1);
  sel_rHB_block2=logistic_double(agebins, selpar_A50_rHB2, selpar_slope_rHB2, selpar_A502_rHB2,selpar_slope2_rHB2);
  sel_rHB_block3=logistic_double(agebins, selpar_A50_rHB3, selpar_slope_rHB3, selpar_A502_rHB3, selpar_slope2_rHB3);
  sel_rGN_block1=sel_rHB_block1;
  sel_rGN_block2=logistic_double(agebins, selpar_A50_rGN2, selpar_slope_rGN2, selpar_A502_rGN2,selpar_slope2_rGN2);
  sel_rGN_block3=logistic(agebins, selpar_A50_rGN3, selpar_slope_rGN3);
  sel_rHB_D_block3=logistic_double(agebins, selpar_A50_rHB3_D, selpar_slope_rHB3_D, selpar_A502_rHB3_D, selpar_slope2_rHB3_D);
  sel_rGN_D_block3=logistic_double(agebins, selpar_A50_rGN3_D, selpar_slope_rGN3_D, selpar_A502_rGN3_D, selpar_slope2_rGN3_D); 
  sel_sCT_vec=logistic_double(agebins, selpar_A50_sCT, selpar_slope_sCT, selpar_A502_sCT, selpar_slope2_sCT);
  sel_sCT_vec((nages_agec_10p+1),nages)=sel_sCT_vec(nages_agec_10p);
  selpar_switch_sVD=0.0;
  for (iage=1; iage<=nages; iage++)
  {if (selpar_switch_sVD==1.0 | sel_sCT_vec(iage)==max(sel_sCT_vec))
    {sel_sVD_vec(iage)=1.0; 
	 selpar_switch_sVD=1.0;
	} else {sel_sVD_vec(iage)=sel_sCT_vec(iage);}
  }
  //BLOCK 1 for selex. No limit through to 20" size limit  
  for (iyear=styr; iyear<=endyr_selex_phase2; iyear++)
   {     
    sel_cHL(iyear)=logistic(agebins, selpar_A50_cHL1, selpar_slope_cHL1);
    sel_rHB(iyear)(1,nages_agec_10p)=sel_rHB_block1(1,nages_agec_10p);
	sel_rHB(iyear)((nages_agec_10p+1),nages)=sel_rHB_block1(nages_agec_10p);
    sel_rGN(iyear)(1,nages_agec_10p)=sel_rGN_block1(1,nages_agec_10p); //NOTE: uses nages_agec_10p bc this selex mirrors rHB
	sel_rGN(iyear)((nages_agec_10p+1),nages)=sel_rGN_block1(nages_agec_10p);
	sel_cHL_D(iyear)=logistic_double(agebins, selpar_A50_cHL2_D, selpar_slope_cHL2_D, selpar_A502_cHL2_D, selpar_slope2_cHL2_D);
	sel_rHB_D(iyear)=logistic_double(agebins, selpar_A50_rHB2_D, selpar_slope_rHB2_D, selpar_A502_rHB2_D, selpar_slope2_rHB2_D);
	sel_rGN_D(iyear)=sel_rHB_D(iyear);
   }
  //BLOCK 2 for selex. 20" size limit until 2010
  for (iyear=(endyr_selex_phase2+1); iyear<=endyr_selex_phase3; iyear++)
   {
    sel_cHL(iyear)=logistic(agebins, selpar_A50_cHL2, selpar_slope_cHL2);
    sel_rHB(iyear)(1,nages_agec_10p)=sel_rHB_block2(1,nages_agec_10p);
	sel_rHB(iyear)((nages_agec_10p+1),nages)=sel_rHB_block2(nages_agec_10p);
    sel_rGN(iyear)(1,nages_agec_10p)=sel_rGN_block2(1,nages_agec_10p); 
	sel_rGN(iyear)((nages_agec_10p+1),nages)=sel_rGN_block2(nages_agec_10p);
    sel_cHL_D(iyear)=logistic_double(agebins, selpar_A50_cHL2_D, selpar_slope_cHL2_D, selpar_A502_cHL2_D, selpar_slope2_cHL2_D);  
	sel_rHB_D(iyear)=logistic_double(agebins, selpar_A50_rHB2_D, selpar_slope_rHB2_D, selpar_A502_rHB2_D, selpar_slope2_rHB2_D);  
	sel_rGN_D(iyear)=sel_rHB_D(iyear);
   }
  //BLOCK 3 for selex.  Moratorium and mini-seasons with no size limit during the mini-season
   for (iyear=(endyr_selex_phase3+1); iyear<=endyr; iyear++)
   {   
    sel_cHL(iyear)=logistic(agebins, selpar_A50_cHL3, selpar_slope_cHL3);
	sel_rHB(iyear)(1,nages_agec_10p)=sel_rHB_block3(1,nages_agec_10p);
	sel_rHB(iyear)((nages_agec_10p+1),nages)=sel_rHB_block3(nages_agec_10p);
    sel_rGN(iyear)(1,nages_agec)=sel_rGN_block3(1,nages_agec);
	sel_rGN(iyear)((nages_agec+1),nages)=sel_rGN_block3(nages_agec);
    sel_cHL_D(iyear)=logistic(agebins, selpar_A50_cHL3_D, selpar_slope_cHL3_D);
	sel_rHB_D(iyear)(1,nages_agec_10p)=sel_rHB_D_block3(1,nages_agec_10p);
	sel_rHB_D(iyear)((nages_agec_10p+1),nages)=sel_rHB_D_block3(nages_agec_10p);
	sel_rGN_D(iyear)(1,nages_agec_10p)=sel_rGN_D_block3(1,nages_agec_10p);
	sel_rGN_D(iyear)((nages_agec_10p+1),nages)=sel_rGN_D_block3(nages_agec_10p);
   }  
   //sel_rGN_D=sel_rHB_D;
   for (iyear=styr; iyear<=endyr; iyear++)
     {     
       sel_sCT(iyear)=sel_sCT_vec;
	   sel_sVD(iyear)=sel_sVD_vec;
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
       F_cHL(iyear)=sel_cHL(iyear)*F_cHL_out(iyear);
       Fsum(iyear)+=F_cHL_out(iyear);
    }
    if(iyear>=styr_L_rHB & iyear<=endyr_L_rHB)
    {  F_rHB_out(iyear)=mfexp(log_avg_F_L_rHB+log_dev_F_L_rHB(iyear)); //}    
       F_rHB(iyear)=sel_rHB(iyear)*F_rHB_out(iyear);
       Fsum(iyear)+=F_rHB_out(iyear);
    }
    if(iyear>=styr_L_rGN & iyear<=endyr_L_rGN) 
    {  F_rGN_out(iyear)=mfexp(log_avg_F_L_rGN+log_dev_F_L_rGN(iyear)); //}    
       F_rGN(iyear)=sel_rGN(iyear)*F_rGN_out(iyear); //(general rec shares headboat selex) fixed 
       Fsum(iyear)+=F_rGN_out(iyear);
    }
    if(iyear>=styr_D_cHL & iyear<=endyr_D_cHL)
    {  F_cHL_D_out(iyear)=mfexp(log_avg_F_D_cHL+log_dev_F_D_cHL(iyear)); //}    
       F_cHL_D(iyear)=sel_cHL_D(iyear)*F_cHL_D_out(iyear);
       Fsum(iyear)+=F_cHL_D_out(iyear);
    }
    if(iyear>=styr_D_rHB & iyear<=endyr_D_rHB)
    {  F_rHB_D_out(iyear)=mfexp(log_avg_F_D_rHB+log_dev_F_D_rHB(iyear)); //}    
       F_rHB_D(iyear)=sel_rHB_D(iyear)*F_rHB_D_out(iyear);
       Fsum(iyear)+=F_rHB_D_out(iyear);
    }
    if(iyear>=styr_D_rGN & iyear<=endyr_D_rGN)
    {  F_rGN_D_out(iyear)=mfexp(log_avg_F_D_rGN+log_dev_F_D_rGN(iyear)); //}    
       F_rGN_D(iyear)=sel_rGN_D(iyear)*F_rGN_D_out(iyear); //general rec shares headboat selex
       Fsum(iyear)+=F_rGN_D_out(iyear);
    }
    //Total F at age
    F(iyear)=F_cHL(iyear);  //first in additive series (NO +=)
    F(iyear)+=F_rHB(iyear);
    F(iyear)+=F_rGN(iyear);
    F(iyear)+=F_cHL_D(iyear);
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
  R_virgin=SR_eq_func(R0, steep, spr_F0, spr_F0, BiasCor, SR_switch);
  B0=bpr_F0*R_virgin;   
  B0_q_DD=R_virgin*sum(elem_prod(N_bpr_F0(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages))); 
  F_init_denom=mfexp(log_avg_F_L_cHL+log_F_dev_init_cHL)+mfexp(log_avg_F_L_rHB+log_F_dev_init_rHB)+mfexp(log_avg_F_L_rGN+log_F_dev_init_rGN);
  F_init_cHL_prop= 1.0; //mfexp(log_avg_F_L_cHL+log_F_dev_init_cHL)/F_init_denom;  //Leave this machinery in, in case a sensitivity is run with a later start year.
  F_init_rHB_prop= 0.0; //mfexp(log_avg_F_L_rHB+log_F_dev_init_rHB)/F_init_denom;
  F_init_rGN_prop= 0.0; //mfexp(log_avg_F_L_rGN+log_F_dev_init_rGN)/F_init_denom;
  F_initial=sel_cHL(styr)*F_init*F_init_cHL_prop+
            sel_rHB(styr)*F_init*F_init_rHB_prop+
            sel_rGN(styr)*F_init*F_init_rGN_prop; //(rGN uses rHB selex)  
  Z_initial=M+F_initial;
  N_spr_initial(1)=1.0*mfexp(-1.0*Z_initial(1)*spawn_time_frac); //at peak spawning time;
  for (iage=2; iage<=nages; iage++)
    {
      N_spr_initial(iage)=N_spr_initial(iage-1)*
                   mfexp(-1.0*(Z_initial(iage-1)*(1.0-spawn_time_frac) + Z_initial(iage)*spawn_time_frac)); 
    }
  N_spr_initial(nages)=N_spr_initial(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
  spr_initial=sum(elem_prod(N_spr_initial,reprod));
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
  SSB(styr)=sum(elem_prod(N_spawn(styr),reprod));
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
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));
		B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));       
    }
    else   //recruitment follows S-R curve with lognormal deviation
    {
        N(iyear+1,1)=SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch)*mfexp(log_dev_rec(iyear+1));
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); //plus group
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod));
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
  for (iyear=styr_D_cHL; iyear<=endyr_D_cHL; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      D_cHL_num(iyear,iage)=N(iyear,iage)*F_cHL_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
    }
    pred_cHL_D_knum(iyear)=sum(D_cHL_num(iyear))/1000.0;            //pred annual dead discards in 1000s (for matching data)
    pred_cHL_D_klb(iyear)=sum(elem_prod(D_cHL_num(iyear),wholewgt_cHL_D_klb(iyear))); //annual dead discards in 1000 lb whole (for output only)
  }
  for (iyear=styr_D_rHB; iyear<=endyr_D_rHB; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      D_rHB_num(iyear,iage)=N(iyear,iage)*F_rHB_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
    }
    pred_rHB_D_knum(iyear)=sum(D_rHB_num(iyear))/1000.0;            //pred annual dead discards in 1000s (for matHBing data)
    pred_rHB_D_klb(iyear)=sum(elem_prod(D_rHB_num(iyear),wholewgt_rHB_D_klb(iyear))); //annual dead discards in 1000 lb whole (for output only)
  }
  for (iyear=styr_D_rGN; iyear<=endyr_D_rGN; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      D_rGN_num(iyear,iage)=N(iyear,iage)*F_rGN_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
    }
    pred_rGN_D_knum(iyear)=sum(D_rGN_num(iyear))/1000.0;            //pred annual dead discards in 1000s (for matGRing data)
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
      for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++)
      {   if (iyear>styr_cpue_rHB & iyear <=2003) 
          {//q_rate_fcn_rHB(iyear)=(1.0+q_rate)*q_rate_fcn_rHB(iyear-1); //compound
             q_rate_fcn_rHB(iyear)=(1.0+(iyear-styr_cpue_rHB)*q_rate)*q_rate_fcn_rHB(styr_cpue_rHB);  //linear
          }
          if (iyear>2003) {q_rate_fcn_rHB(iyear)=q_rate_fcn_rHB(iyear-1);} 
      }   
      for (iyear=styr_D_cpue_rHB; iyear<=endyr_D_cpue_rHB; iyear++)
      {   if (iyear>styr_D_cpue_rHB & iyear <=2003) 
          {//q_rate_fcn_rHB_D(iyear)=(1.0+q_rate)*q_rate_fcn_rHB_D(iyear-1); //compound
             q_rate_fcn_rHB_D(iyear)=(1.0+(iyear-styr_D_cpue_rHB)*q_rate)*q_rate_fcn_rHB_D(styr_D_cpue_rHB);  //linear
          }
          if (iyear>2003) {q_rate_fcn_rHB_D(iyear)=q_rate_fcn_rHB_D(iyear-1);} 
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
      pred_cHL_cpue(iyear)=q_cHL(iyear)*q_rate_fcn_cHL(iyear)*q_DD_fcn(iyear)*sum(N_cHL(iyear));
      if (iyear<endyr_cpue_cHL){q_cHL(iyear+1)=q_cHL(iyear)*mfexp(q_RW_log_dev_cHL(iyear));}
  }
 //rHB  cpue
  q_rHB(styr_cpue_rHB)=mfexp(log_q_cpue_rHB); 
  for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++)
  {   
      N_rHB(iyear)=elem_prod(N_mdyr(iyear),sel_rHB(iyear)); 
      pred_rHB_cpue(iyear)=q_rHB(iyear)*q_rate_fcn_rHB(iyear)*q_DD_fcn(iyear)*sum(N_rHB(iyear));
      if (iyear<endyr_cpue_rHB){q_rHB(iyear+1)=q_rHB(iyear)*mfexp(q_RW_log_dev_rHB(iyear));}
  }
 //rHB disc cpue         
  q_rHB_D(styr_D_cpue_rHB)=mfexp(log_q_cpue_rHB_D); 
  for (iyear=styr_D_cpue_rHB; iyear<=endyr_D_cpue_rHB; iyear++)
  {   
      N_rHB_D(iyear)=elem_prod(N_mdyr(iyear),sel_rHB_D(endyr_selex_phase3));  //selectivity from block with 20inch size limit, as in the index
      pred_rHB_D_cpue(iyear)=q_rHB_D(iyear)*q_rate_fcn_rHB_D(iyear)*q_DD_fcn(iyear)*sum(N_rHB_D(iyear));
      if (iyear<endyr_D_cpue_rHB){q_rHB_D(iyear+1)=q_rHB_D(iyear)*mfexp(q_RW_log_dev_rHB_D(iyear));}
  } 
 //SERFS sCT cpue
  q_sCT(styr_cpue_sCT)=mfexp(log_q_cpue_sCT); 
  for (iyear=styr_cpue_sCT; iyear<=endyr_cpue_sCT; iyear++)
  {   
    N_sCT(iyear)=elem_prod(N_mdyr(iyear),sel_sCT(iyear)); 
	pred_sCT_cpue(iyear)=q_sCT(iyear)*q_DD_fcn(iyear)*sum(N_sCT(iyear));
    if (iyear<endyr_cpue_sCT){q_sCT(iyear+1)=q_sCT(iyear)*mfexp(q_RW_log_dev_sCT(iyear));}
  } 
 //SERFS sCT cpue
  q_sVD(styr_cpue_sVD)=mfexp(log_q_cpue_sVD); 
  for (iyear=styr_cpue_sVD; iyear<=endyr_cpue_sVD; iyear++)
  {   
    N_sVD(iyear)=elem_prod(N_mdyr(iyear),sel_sVD(iyear)); 
	pred_sVD_cpue(iyear)=q_sVD(iyear)*q_DD_fcn(iyear)*sum(N_sVD(iyear));
    if (iyear<endyr_cpue_sVD){q_sVD(iyear+1)=q_sVD(iyear)*mfexp(q_RW_log_dev_sVD(iyear));}
  } 
}

void model_parameters::get_length_comps(void)
{
  //------cGN handline 
  for (iyear=1;iyear<=nyr_lenc_cHL;iyear++)
  {pred_cHL_lenc(iyear)=(L_cHL_num(yrs_lenc_cHL(iyear))*lenprob_cHL)/sum(L_cHL_num(yrs_lenc_cHL(iyear)));}
  //------cGN discards  
  //for (iyear=1;iyear<=nyr_lenc_cHL_D;iyear++) 
  //{pred_cHL_D_lenc(iyear)=(D_cHL_num(yrs_lenc_cHL_D(iyear))*lenprob_cHL_D)/sum(D_cHL_num(yrs_lenc_cHL_D(iyear)));}
   for (iyear=1;iyear<=nyr_cHL_D_lenc_pool1;iyear++)
	{pred_cHL_D_lenc_yr1(iyear)=(D_cHL_num(yrs_cHL_D_lenc_pool1(iyear))*lenprob_cHL_D)/sum(D_cHL_num(yrs_cHL_D_lenc_pool1(iyear)));}
   for (iyear=1;iyear<=nyr_cHL_D_lenc_pool2;iyear++)
	{pred_cHL_D_lenc_yr2(iyear)=(D_cHL_num(yrs_cHL_D_lenc_pool2(iyear))*lenprob_cHL_D)/sum(D_cHL_num(yrs_cHL_D_lenc_pool2(iyear)));}
   pred_cHL_D_lenc.initialize();
   for (iyear=1;iyear<=nyr_cHL_D_lenc_pool1;iyear++)
     {pred_cHL_D_lenc(1) += nsamp_cHL_D_lenc_pool1(iyear) * pred_cHL_D_lenc_yr1(iyear);}
   pred_cHL_D_lenc(1)=pred_cHL_D_lenc(1)/sum(nsamp_cHL_D_lenc_pool1);
   for (iyear=1;iyear<=nyr_cHL_D_lenc_pool2;iyear++)
     {pred_cHL_D_lenc(2) += nsamp_cHL_D_lenc_pool2(iyear) * pred_cHL_D_lenc_yr2(iyear);}
   pred_cHL_D_lenc(2)=pred_cHL_D_lenc(2)/sum(nsamp_cHL_D_lenc_pool2);
  //------headboat discards  
  for (iyear=1;iyear<=nyr_lenc_rHB_D;iyear++) 
  {pred_rHB_D_lenc(iyear)=(D_rHB_num(yrs_lenc_rHB_D(iyear))*lenprob_rHB_D)/sum(D_rHB_num(yrs_lenc_rHB_D(iyear)));}
  //------gen rec discards  
  for (iyear=1;iyear<=nyr_lenc_rGN_D;iyear++) 
  {pred_rGN_D_lenc(iyear)=(D_rGN_num(yrs_lenc_rGN_D(iyear))*lenprob_rGN_D)/sum(D_rGN_num(yrs_lenc_rGN_D(iyear)));}
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
  //Headboat
 for (iyear=1;iyear<=nyr_agec_rHB;iyear++)
  {
    ErrorFree_rHB_agec(iyear)=L_rHB_num(yrs_agec_rHB(iyear))/sum(L_rHB_num(yrs_agec_rHB(iyear)));
    pred_rHB_agec_allages(iyear)=age_error*ErrorFree_rHB_agec(iyear); 
    for (iage=1; iage<=nages_agec_10p; iage++) {pred_rHB_agec(iyear,iage)=pred_rHB_agec_allages(iyear,iage);} 
    for (iage=(nages_agec_10p+1); iage<=nages; iage++) {pred_rHB_agec(iyear,nages_agec_10p)+=pred_rHB_agec_allages(iyear,iage);} //plus group                        
  }
   //General Recreational  
 for (iyear=1;iyear<=nyr_agec_rGN;iyear++)
  {
    ErrorFree_rGN_agec(iyear)=L_rGN_num(yrs_agec_rGN(iyear))/sum(L_rGN_num(yrs_agec_rGN(iyear)));
    pred_rGN_agec_allages(iyear)=age_error*ErrorFree_rGN_agec(iyear); 
    for (iage=1; iage<=nages_agec_10p; iage++) {pred_rGN_agec(iyear,iage)=pred_rGN_agec_allages(iyear,iage);} 
    for (iage=(nages_agec_10p+1); iage<=nages; iage++) {pred_rGN_agec(iyear,nages_agec_10p)+=pred_rGN_agec_allages(iyear,iage);} //plus group                        
  }
    //SERFS sCT
 for (iyear=1;iyear<=nyr_agec_sCT;iyear++)
  {
    ErrorFree_sCT_agec(iyear)=N_sCT(yrs_agec_sCT(iyear))/sum(N_sCT(yrs_agec_sCT(iyear)));
    pred_sCT_agec_allages(iyear)=age_error*ErrorFree_sCT_agec(iyear); 
    for (iage=1; iage<=nages_agec_10p; iage++) {pred_sCT_agec(iyear,iage)=pred_sCT_agec_allages(iyear,iage);} 
    for (iage=(nages_agec_10p+1); iage<=nages; iage++) {pred_sCT_agec(iyear,nages_agec_10p)+=pred_sCT_agec_allages(iyear,iage);} //plus group                        
  }
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
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_D_cHL+
        sum(log_dev_F_D_cHL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
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
  F_cHL_D_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_D_cHL+
        sum(log_dev_F_D_cHL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_rHB_D_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_D_rHB+
        sum(log_dev_F_D_rHB((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_rGN_D_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_D_rGN+
        sum(log_dev_F_D_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  log_F_dev_end_cHL=sum(log_dev_F_L_cHL((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_rHB=sum(log_dev_F_L_rHB((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  
  log_F_dev_end_rGN=sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  
  log_F_dev_end_cHL_D=sum(log_dev_F_D_cHL((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_rHB_D=sum(log_dev_F_D_rHB((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_rGN_D=sum(log_dev_F_D_rGN((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  F_end_L=sel_cHL(endyr)*mfexp(log_avg_F_L_cHL+log_F_dev_end_cHL)+
          sel_rHB(endyr)*mfexp(log_avg_F_L_rHB+log_F_dev_end_rHB)+
          sel_rGN(endyr)*mfexp(log_avg_F_L_rGN+log_F_dev_end_rGN);   // (rGN uses rHB selex ) 
  F_end_D=sel_cHL_D(endyr)*mfexp(log_avg_F_D_cHL+log_F_dev_end_cHL_D)+
          sel_rHB_D(endyr)*mfexp(log_avg_F_D_rHB+log_F_dev_end_rHB_D)+
          sel_rGN_D(endyr)*mfexp(log_avg_F_D_rGN+log_F_dev_end_rGN_D);  // (rGN uses rHB selex )  
  F_end=F_end_L+F_end_D;
  F_end_apex=max(F_end);
  sel_wgted_tot=F_end/F_end_apex;
  sel_wgted_L=elem_prod(sel_wgted_tot, elem_div(F_end_L,F_end));
  sel_wgted_D=elem_prod(sel_wgted_tot, elem_div(F_end_D,F_end));
  wgt_wgted_L_denom=F_cHL_prop+F_rHB_prop+F_rGN_prop;  
  wgt_wgted_L_klb=F_cHL_prop/wgt_wgted_L_denom*wholewgt_cHL_klb(endyr)+  
                  F_rHB_prop/wgt_wgted_L_denom*wholewgt_rHB_klb(endyr)+
                  F_rGN_prop/wgt_wgted_L_denom*wholewgt_rGN_klb(endyr);                          
  wgt_wgted_D_denom=F_cHL_D_prop+F_rHB_D_prop+F_rGN_D_prop;  
  wgt_wgted_D_klb=F_cHL_D_prop/wgt_wgted_D_denom*wholewgt_cHL_D_klb(endyr)+ 
                  F_rHB_D_prop/wgt_wgted_D_denom*wholewgt_rHB_D_klb(endyr)+
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
	E_eq_wgt(ff)=(sum(elem_prod(L_age_msy(E_age_min, nages),wgt_wgted_L_klb(E_age_min, nages)))+
				 sum(elem_prod(D_age_msy(E_age_min, nages),wgt_wgted_D_klb(E_age_min, nages))) )/
				 (mt2klb*sum(elem_prod(N_age_msy(E_age_min, nages),wgt_mt(E_age_min, nages))));
	E_eq_num(ff)=(sum(L_age_msy(E_age_min, nages))+sum(D_age_msy(E_age_min, nages)) )/
				 (sum(N_age_msy(E_age_min, nages)));
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
		E_msy_wgt=E_eq_wgt(ff);
		E_msy_num=E_eq_num(ff);
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
    F_D_age_spr=0.0;
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
    switch(SPR_rec_switch){
       case 1: //Observed mean over years specified
          rec_mean=sum(rec(styr_rec_spr, endyr_rec_spr))/nyrs_rec_spr;
          break;
       case 2: //Expected recruitment is R0 with bias correction (if that feature is turned on)
          rec_mean=BiasCor*R0;
          break;
       default: // no such switch available
          cout << "Error in input: SPR_rec_switch must be set to 1 or 2." << endl;
          cout << "Presently it is set to " << SPR_rec_switch <<"."<< endl;
          exit(0);          
	}	  
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
  L_F30_klb_out=sum(elem_prod(L_age_F30,wgt_wgted_L_klb)); //in whole weight
  L_F30_knum_out=sum(L_age_F30)/1000.0;  
  D_F30_klb_out=sum(elem_prod(D_age_F30,wgt_wgted_D_klb)); //in whole weight   
  D_F30_knum_out=sum(D_age_F30)/1000.0;    
  E_F30_wgt=(sum(elem_prod(L_age_F30(E_age_min, nages),wgt_wgted_L_klb(E_age_min, nages)))+
				 sum(elem_prod(D_age_F30(E_age_min, nages),wgt_wgted_D_klb(E_age_min, nages))) )/
				 (mt2klb*sum(elem_prod(N_age_spr(E_age_min, nages),wgt_mt(E_age_min, nages))));
  E_F30_num=(sum(L_age_F30(E_age_min, nages))+sum(D_age_F30(E_age_min, nages)) )/
				 (sum(N_age_spr(E_age_min, nages))); 
}

void model_parameters::get_miscellaneous_stuff(void)
{
  if(var_rec_dev>0.0)
   {sigma_rec_dev=sqrt(var_rec_dev);} //pow(var_rec_dev,0.5);  //sample SD of predicted residuals (may not equal rec_sigma)  
   else{sigma_rec_dev=0.0;}
  len_cv=elem_div(len_sd,meanlen_TL);
  len_cv_L=elem_div(len_sd_L,meanlen_TL_L);
  len_cv_20=elem_div(len_sd_20,meanlen_TL_20);
 //Time series of interest  
  B(endyr+1)=elem_prod(N(endyr+1),wgt_mt);
  totN(endyr+1)=sum(N(endyr+1));
  totB(endyr+1)=sum(B(endyr+1));  
  SdS0=SSB/S0;
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
  E_total_wgt.initialize();
  E_L_wgt.initialize();
  E_D_wgt.initialize();
  E_total_num.initialize();
  E_L_num.initialize();
  E_D_num.initialize();
  for(iyear=styr; iyear<=endyr; iyear++)
  {
        L_total_klb_yr(iyear)=pred_cHL_L_klb(iyear)+pred_rHB_L_klb(iyear)+pred_rGN_L_klb(iyear);
        L_total_knum_yr(iyear)=pred_cHL_L_knum(iyear)+pred_rHB_L_knum(iyear)+pred_rGN_L_knum(iyear);
        B(iyear)=elem_prod(N(iyear),wgt_mt);
        totN(iyear)=sum(N(iyear));
        totB(iyear)=sum(B(iyear));   
        if (iyear>=styr_D_cHL && iyear<=endyr_D_cHL)
        {
         D_total_knum_yr(iyear)+=pred_cHL_D_knum(iyear);
         D_total_klb_yr(iyear)+=pred_cHL_D_klb(iyear);
         D_cHL_klb(iyear)=elem_prod(D_cHL_num(iyear),wholewgt_cHL_D_klb(iyear));     //in 1000 lb 
        }
        if (iyear>=styr_D_rHB && iyear<=endyr_D_rHB)
        {
         D_total_knum_yr(iyear)+=pred_rHB_D_knum(iyear);
         D_total_klb_yr(iyear)+=pred_rHB_D_klb(iyear);
         D_rHB_klb(iyear)=elem_prod(D_rHB_num(iyear),wholewgt_rHB_D_klb(iyear));     //in 1000 lb 
        }    
        if (iyear>=styr_D_rGN && iyear<=endyr_D_rGN)
        {
         D_total_knum_yr(iyear)+=pred_rGN_D_knum(iyear);
         D_total_klb_yr(iyear)+=pred_rGN_D_klb(iyear);
         D_rGN_klb(iyear)=elem_prod(D_rGN_num(iyear),wholewgt_rGN_D_klb(iyear));     //in 1000 lb 
        }                          
	E_L_num(iyear)= (sum(L_cHL_num(iyear)(E_age_min, nages)) + 
	                 sum(L_rHB_num(iyear)(E_age_min, nages)) + 
					 sum(L_rGN_num(iyear)(E_age_min, nages)))/
					 sum(N(iyear)(E_age_min, nages));
	E_D_num(iyear)= (sum(D_cHL_num(iyear)(E_age_min, nages)) + 
					 sum(D_rHB_num(iyear)(E_age_min, nages)) + 
					 sum(D_rGN_num(iyear)(E_age_min, nages)))/
					 sum(N(iyear)(E_age_min, nages));
    E_total_num=E_L_num+E_D_num;	
	E_L_wgt(iyear)= (sum(L_cHL_klb(iyear)(E_age_min, nages)) + 
	                 sum(L_rHB_klb(iyear)(E_age_min, nages)) + 
					 sum(L_rGN_klb(iyear)(E_age_min, nages)))/
					 (mt2klb*sum(B(iyear)(E_age_min, nages)));		
	E_D_wgt(iyear)= (sum(D_cHL_klb(iyear)(E_age_min, nages)) + 
	                 sum(D_rHB_klb(iyear)(E_age_min, nages)) + 
					 sum(D_rGN_klb(iyear)(E_age_min, nages)))/
					 (mt2klb*sum(B(iyear)(E_age_min, nages)));		
    E_total_wgt=E_L_wgt+E_D_wgt;	
  }
  L_total_num=L_cHL_num+L_rHB_num+L_rGN_num;   //landings at age in number fish
  L_total_klb=L_cHL_klb+L_rHB_klb+L_rGN_klb;   //landings at age in klb whole weight
  D_total_num=(D_cHL_num+D_rHB_num+D_rGN_num);          //discards at age in number fish
  D_total_klb=D_cHL_klb+D_rHB_klb+D_rGN_klb;            //discards at age in klb whole weight
  Fend_mean_temp=1.0;
  for (iyear=1; iyear<=selpar_n_yrs_wgted; iyear++) {Fend_mean_temp*=Fapex(endyr-iyear+1);}
  Fend_mean=pow(Fend_mean_temp,(1.0/selpar_n_yrs_wgted));	  
  Eend_mean_temp_wgt=1.0;
  for (iyear=1; iyear<=selpar_n_yrs_wgted; iyear++) {Eend_mean_temp_wgt*=E_total_wgt(endyr-iyear+1);}
  E_wgt_end_mean=pow(Eend_mean_temp_wgt,(1.0/selpar_n_yrs_wgted));	  
  Eend_mean_temp_num=1.0;
  for (iyear=1; iyear<=selpar_n_yrs_wgted; iyear++) {Eend_mean_temp_num*=E_total_num(endyr-iyear+1);}
  E_num_end_mean=pow(Eend_mean_temp_num,(1.0/selpar_n_yrs_wgted));	  
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
   if (E_msy_wgt>0)
    {
	  EdEmsy_wgt=E_total_wgt/E_msy_wgt;
	  EdEF30_wgt=E_total_wgt/E_F30_wgt;
	  EdEmsy_wgt_end_mean=E_wgt_end_mean/E_msy_wgt;
	  EdEF30_wgt_end_mean=E_wgt_end_mean/E_F30_wgt;
    }	   
   if (E_msy_num>0)
    {
	  EdEmsy_num=E_total_num/E_msy_num;
	  EdEF30_num=E_total_wgt/E_F30_num;
	  EdEmsy_num_end_mean=E_num_end_mean/E_msy_num;
	  EdEF30_num_end_mean=E_num_end_mean/E_F30_num;
	}	   
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
	  neff_rHB_D_lenc_allyr_out=missing; 
	  neff_rGN_D_lenc_allyr_out=missing; 
      neff_cHL_agec_allyr_out=missing;
      neff_rHB_agec_allyr_out=missing;
	  neff_sCT_agec_allyr_out=missing;
	  neff_rGN_agec_allyr_out=missing;
      for (iyear=1; iyear<=nyr_lenc_cHL; iyear++)
         {if (nsamp_lenc_cHL(iyear)>=minSS_lenc_cHL)
            {neff_cHL_lenc_allyr_out(yrs_lenc_cHL(iyear))=multinom_eff_N(pred_cHL_lenc(iyear),obs_lenc_cHL(iyear));} 
          else {neff_cHL_lenc_allyr_out(yrs_lenc_cHL(iyear))=-99;}
         }
      for (iyear=1; iyear<=nyr_lenc_rHB_D; iyear++)
         {if (nsamp_lenc_rHB_D(iyear)>=minSS_lenc_rHB_D)
            {neff_rHB_D_lenc_allyr_out(yrs_lenc_rHB_D(iyear))=multinom_eff_N(pred_rHB_D_lenc(iyear),obs_lenc_rHB_D(iyear));} 
          else {neff_rHB_D_lenc_allyr_out(yrs_lenc_rHB_D(iyear))=-99;}
         }
     for (iyear=1; iyear<=nyr_lenc_rGN_D; iyear++)
         {if (nsamp_lenc_rGN_D(iyear)>=minSS_lenc_rGN_D)
            {neff_rGN_D_lenc_allyr_out(yrs_lenc_rGN_D(iyear))=multinom_eff_N(pred_rGN_D_lenc(iyear),obs_lenc_rGN_D(iyear));} 
          else {neff_rGN_D_lenc_allyr_out(yrs_lenc_rGN_D(iyear))=-99;}
         }
      for (iyear=1; iyear<=nyr_agec_cHL; iyear++)
         {if (nsamp_agec_cHL(iyear)>=minSS_agec_cHL)
            {neff_cHL_agec_allyr_out(yrs_agec_cHL(iyear))=multinom_eff_N(pred_cHL_agec(iyear),obs_agec_cHL(iyear));}                            
          else {neff_cHL_agec_allyr_out(yrs_agec_cHL(iyear))=-99;}
         }    
      for (iyear=1; iyear<=nyr_agec_rHB; iyear++)
         {if (nsamp_agec_rHB(iyear)>=minSS_agec_rHB)
            {neff_rHB_agec_allyr_out(yrs_agec_rHB(iyear))=multinom_eff_N(pred_rHB_agec(iyear),obs_agec_rHB(iyear));}                            
          else {neff_rHB_agec_allyr_out(yrs_agec_rHB(iyear))=-99;}
         }
	  for (iyear=1; iyear<=nyr_agec_sCT; iyear++)  
           {if (nsamp_agec_sCT(iyear)>=minSS_agec_sCT)
              {neff_sCT_agec_allyr_out(yrs_agec_sCT(iyear))=multinom_eff_N(pred_sCT_agec(iyear),obs_agec_sCT(iyear));}                            
            else {neff_sCT_agec_allyr_out(yrs_agec_sCT(iyear))=-99;}
           }
	  for (iyear=1; iyear<=nyr_agec_rGN; iyear++)  
           {if (nsamp_agec_rGN(iyear)>=minSS_agec_rGN)
              {neff_rGN_agec_allyr_out(yrs_agec_rGN(iyear))=multinom_eff_N(pred_rGN_agec(iyear),obs_agec_rGN(iyear));}                            
            else {neff_rGN_agec_allyr_out(yrs_agec_rGN(iyear))=-99;}
           }
}

void model_parameters::evaluate_objective_function(void)
{
  //fval=square(xdum-9.0);
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
  f_rHB_D_cpue=0.0;
  f_rHB_D_cpue=lk_lognormal(pred_rHB_D_cpue, obs_cpue_rHB_D, obs_cv_cpue_rHB_D, w_cpue_rHB_D);
  fval+=f_rHB_D_cpue;
  fval_data+=f_rHB_D_cpue;  
  f_sCT_cpue=0.0;
  f_sCT_cpue=lk_lognormal(pred_sCT_cpue, obs_cpue_sCT, obs_cv_cpue_sCT, w_cpue_sCT);
  fval+=w_cpue_sCT_mult*f_sCT_cpue;
  fval_data+=f_sCT_cpue;  
  f_sVD_cpue=0.0;
  f_sVD_cpue=lk_lognormal(pred_sVD_cpue, obs_cpue_sVD, obs_cv_cpue_sVD, w_cpue_sVD);
  fval+=w_cpue_sVD_mult*f_sVD_cpue;
  fval_data+=f_sVD_cpue;  
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
  //f_cHL_D in 1000 fish
  f_cHL_D=lk_lognormal(pred_cHL_D_knum(styr_D_cHL,endyr_D_cHL), obs_cHL_D(styr_D_cHL,endyr_D_cHL), 
                      obs_cv_D_cHL(styr_D_cHL,endyr_D_cHL), w_D);
  fval+=f_cHL_D;
  fval_data+=f_cHL_D;  
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
  //f_cHL_lenc
  //f_cHL_lenc=lk_robust_multinomial(nsamp_lenc_cHL, pred_cHL_lenc, obs_lenc_cHL, nyr_lenc_cHL, double(nlenbins), minSS_lenc_cHL, w_lenc_cHL);
  f_cHL_lenc=lk_dirichlet_multinomial(nsamp_lenc_cHL, pred_cHL_lenc, obs_lenc_cHL, nyr_lenc_cHL, double(nlenbins), minSS_lenc_cHL, log_dm_lenc_cHL);
  fval+=f_cHL_lenc;
  fval_data+=f_cHL_lenc;
  //f_cHL_D_lenc
  //f_cHL_D_lenc=lk_robust_multinomial(nsamp_lenc_cHL_D, pred_cHL_D_lenc, obs_lenc_cHL_D, nyr_lenc_cHL_D, double(nlenbins), minSS_lenc_cHL_D, w_lenc_cHL_D);
  f_cHL_D_lenc=lk_dirichlet_multinomial(nsamp_lenc_cHL_D, pred_cHL_D_lenc, obs_lenc_cHL_D, nyr_lenc_cHL_D, double(nlenbins), minSS_lenc_cHL_D, log_dm_lenc_cHL_D);
  fval+=f_cHL_D_lenc;
  fval_data+=f_cHL_D_lenc;
  //f_rHB_D_lenc
  //f_rHB_D_lenc=lk_robust_multinomial(nsamp_lenc_rHB_D, pred_rHB_D_lenc, obs_lenc_rHB_D, nyr_lenc_rHB_D, double(nlenbins), minSS_lenc_rHB_D, w_lenc_rHB_D);
  f_rHB_D_lenc=lk_dirichlet_multinomial(nsamp_lenc_rHB_D, pred_rHB_D_lenc, obs_lenc_rHB_D, nyr_lenc_rHB_D, double(nlenbins), minSS_lenc_rHB_D, log_dm_lenc_rHB_D); 
  fval+=f_rHB_D_lenc;
  fval_data+=f_rHB_D_lenc;
  //f_rGN_D_lenc
  //f_rGN_D_lenc=lk_robust_multinomial(nsamp_lenc_rGN_D, pred_rGN_D_lenc, obs_lenc_rGN_D, nyr_lenc_rGN_D, double(nlenbins), minSS_lenc_rGN_D, w_lenc_rGN_D);
  f_rGN_D_lenc=lk_dirichlet_multinomial(nsamp_lenc_rGN_D, pred_rGN_D_lenc, obs_lenc_rGN_D, nyr_lenc_rGN_D, double(nlenbins), minSS_lenc_rGN_D, log_dm_lenc_rGN_D); 
  fval+=f_rGN_D_lenc;
  fval_data+=f_rGN_D_lenc;
  //f_cHL_agec
  //f_cHL_agec=lk_robust_multinomial(nsamp_agec_cHL, pred_cHL_agec, obs_agec_cHL, nyr_agec_cHL, double(nages_agec), minSS_agec_cHL, w_agec_cHL);
  f_cHL_agec=lk_dirichlet_multinomial(nsamp_agec_cHL, pred_cHL_agec, obs_agec_cHL, nyr_agec_cHL, double(nages_agec), minSS_agec_cHL, log_dm_agec_cHL);
  fval+=f_cHL_agec;
  fval_data+=f_cHL_agec;
  //f_rHB_agec
  //f_rHB_agec=lk_robust_multinomial(nsamp_agec_rHB, pred_rHB_agec, obs_agec_rHB, nyr_agec_rHB, double(nages_agec_10p), minSS_agec_rHB, w_agec_rHB);
  f_rHB_agec=lk_dirichlet_multinomial(nsamp_agec_rHB, pred_rHB_agec, obs_agec_rHB, nyr_agec_rHB, double(nages_agec_10p), minSS_agec_rHB, log_dm_agec_rHB);
  fval+=f_rHB_agec;
  fval_data+=f_rHB_agec;
  //f_sCT_agec
  //f_sCT_agec=lk_robust_multinomial(nsamp_agec_sCT, pred_sCT_agec, obs_agec_sCT, nyr_agec_sCT, double(nages_agec_10p), minSS_agec_sCT, w_agec_sCT);
  f_sCT_agec=lk_dirichlet_multinomial(nsamp_agec_sCT, pred_sCT_agec, obs_agec_sCT, nyr_agec_sCT, double(nages_agec_10p), minSS_agec_sCT, log_dm_agec_sCT);  
  fval+=f_sCT_agec;
  fval_data+=f_sCT_agec;
  //f_rGN_agec
  //f_rGN_agec=lk_robust_multinomial(nsamp_agec_rGN, pred_rGN_agec, obs_agec_rGN, nyr_agec_rGN, double(nages_agec_10p), minSS_agec_rGN, w_agec_rGN);
  f_rGN_agec=lk_dirichlet_multinomial(nsamp_agec_rGN, pred_rGN_agec, obs_agec_rGN, nyr_agec_rGN, double(nages_agec_10p), minSS_agec_rGN, log_dm_agec_rGN);  
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
  f_priors=0.0; 
  f_priors+=neg_log_prior(len_cv_val,set_len_cv(5),set_len_cv(6),set_len_cv(7));
  f_priors+=neg_log_prior(len_cv_val_20,set_len_cv_20(5),set_len_cv_20(6),set_len_cv_20(7));
  f_priors+=neg_log_prior(len_cv_val_L,set_len_cv_L(5),set_len_cv_L(6),set_len_cv_L(7));
  //f_priors+=neg_log_prior(steep,set_steep(5),set_log_R0(6),set_log_R0(7)); 
  f_priors+=neg_log_prior(log_R0,set_log_R0(5),set_log_R0(6),set_log_R0(7)); 
  f_priors+=neg_log_prior(R_autocorr,set_R_autocorr(5),set_R_autocorr(6),set_R_autocorr(7));
  f_priors+=neg_log_prior(rec_sigma,set_rec_sigma(5),set_rec_sigma(6),set_rec_sigma(7));
  f_priors+=neg_log_prior(selpar_A50_cHL1,set_selpar_A50_cHL1(5), set_selpar_A50_cHL1(6), set_selpar_A50_cHL1(7));
  f_priors+=neg_log_prior(selpar_slope_cHL1,set_selpar_slope_cHL1(5), set_selpar_slope_cHL1(6), set_selpar_slope_cHL1(7));
  f_priors+=neg_log_prior(selpar_A50_cHL2,set_selpar_A50_cHL2(5), set_selpar_A50_cHL2(6), set_selpar_A50_cHL2(7));
  f_priors+=neg_log_prior(selpar_slope_cHL2,set_selpar_slope_cHL2(5), set_selpar_slope_cHL2(6), set_selpar_slope_cHL2(7));
  f_priors+=neg_log_prior(selpar_A50_cHL3,set_selpar_A50_cHL3(5), set_selpar_A50_cHL3(6), set_selpar_A50_cHL3(7));
  f_priors+=neg_log_prior(selpar_slope_cHL3,set_selpar_slope_cHL3(5), set_selpar_slope_cHL3(6), set_selpar_slope_cHL3(7));
  f_priors+=neg_log_prior(selpar_A50_rHB1,set_selpar_A50_rHB1(5), set_selpar_A50_rHB1(6), set_selpar_A50_rHB1(7));
  f_priors+=neg_log_prior(selpar_slope_rHB1,set_selpar_slope_rHB1(5), set_selpar_slope_rHB1(6), set_selpar_slope_rHB1(7));
  f_priors+=neg_log_prior(selpar_A502_rHB1,set_selpar_A502_rHB1(5), set_selpar_A502_rHB1(6), set_selpar_A502_rHB1(7));
  f_priors+=neg_log_prior(selpar_slope2_rHB1,set_selpar_slope2_rHB1(5), set_selpar_slope2_rHB1(6), set_selpar_slope2_rHB1(7));
  f_priors+=neg_log_prior(selpar_A50_rHB2,set_selpar_A50_rHB2(5), set_selpar_A50_rHB2(6), set_selpar_A50_rHB2(7));
  f_priors+=neg_log_prior(selpar_slope_rHB2,set_selpar_slope_rHB2(5), set_selpar_slope_rHB2(6), set_selpar_slope_rHB2(7));
  f_priors+=neg_log_prior(selpar_A502_rHB2,set_selpar_A502_rHB2(5), set_selpar_A502_rHB2(6), set_selpar_A502_rHB2(7));
  f_priors+=neg_log_prior(selpar_slope2_rHB2,set_selpar_slope2_rHB2(5), set_selpar_slope2_rHB2(6), set_selpar_slope2_rHB2(7));
  f_priors+=neg_log_prior(selpar_A50_rHB3,set_selpar_A50_rHB3(5), set_selpar_A50_rHB3(6), set_selpar_A50_rHB3(7));
  f_priors+=neg_log_prior(selpar_slope_rHB3,set_selpar_slope_rHB3(5), set_selpar_slope_rHB3(6), set_selpar_slope_rHB3(7));
  f_priors+=neg_log_prior(selpar_A502_rHB3,set_selpar_A502_rHB3(5), set_selpar_A502_rHB3(6), set_selpar_A502_rHB3(7));
  f_priors+=neg_log_prior(selpar_slope2_rHB3,set_selpar_slope2_rHB3(5), set_selpar_slope2_rHB3(6), set_selpar_slope2_rHB3(7));
  f_priors+=neg_log_prior(selpar_A50_rGN2,set_selpar_A50_rGN2(5), set_selpar_A50_rGN2(6), set_selpar_A50_rGN2(7));
  f_priors+=neg_log_prior(selpar_slope_rGN2,set_selpar_slope_rGN2(5), set_selpar_slope_rGN2(6), set_selpar_slope_rGN2(7));
  f_priors+=neg_log_prior(selpar_A502_rGN2,set_selpar_A502_rGN2(5), set_selpar_A502_rGN2(6), set_selpar_A502_rGN2(7));
  f_priors+=neg_log_prior(selpar_slope2_rGN2,set_selpar_slope2_rGN2(5), set_selpar_slope2_rGN2(6), set_selpar_slope2_rGN2(7));
  f_priors+=neg_log_prior(selpar_A50_rGN3,set_selpar_A50_rGN3(5), set_selpar_A50_rGN3(6), set_selpar_A50_rGN3(7));
  f_priors+=neg_log_prior(selpar_slope_rGN3,set_selpar_slope_rGN3(5), set_selpar_slope_rGN3(6), set_selpar_slope_rGN3(7));
  f_priors+=neg_log_prior(selpar_A50_cHL2_D,set_selpar_A50_cHL2_D(5), set_selpar_A50_cHL2_D(6), set_selpar_A50_cHL2_D(7));
  f_priors+=neg_log_prior(selpar_slope_cHL2_D,set_selpar_slope_cHL2_D(5), set_selpar_slope_cHL2_D(6), set_selpar_slope_cHL2_D(7));
  f_priors+=neg_log_prior(selpar_A502_cHL2_D,set_selpar_A502_cHL2_D(5), set_selpar_A502_cHL2_D(6), set_selpar_A502_cHL2_D(7));
  f_priors+=neg_log_prior(selpar_slope2_cHL2_D,set_selpar_slope2_cHL2_D(5), set_selpar_slope2_cHL2_D(6), set_selpar_slope2_cHL2_D(7));
  f_priors+=neg_log_prior(selpar_A50_cHL3_D,set_selpar_A50_cHL3_D(5), set_selpar_A50_cHL3_D(6), set_selpar_A50_cHL3_D(7));
  f_priors+=neg_log_prior(selpar_slope_cHL3_D,set_selpar_slope_cHL3_D(5), set_selpar_slope_cHL3_D(6), set_selpar_slope_cHL3_D(7));
  f_priors+=neg_log_prior(selpar_A50_rHB2_D,set_selpar_A50_rHB2_D(5), set_selpar_A50_rHB2_D(6), set_selpar_A50_rHB2_D(7));
  f_priors+=neg_log_prior(selpar_slope_rHB2_D,set_selpar_slope_rHB2_D(5), set_selpar_slope_rHB2_D(6), set_selpar_slope_rHB2_D(7));
  f_priors+=neg_log_prior(selpar_A502_rHB2_D,set_selpar_A502_rHB2_D(5), set_selpar_A502_rHB2_D(6), set_selpar_A502_rHB2_D(7));
  f_priors+=neg_log_prior(selpar_slope2_rHB2_D,set_selpar_slope2_rHB2_D(5), set_selpar_slope2_rHB2_D(6), set_selpar_slope2_rHB2_D(7));
  f_priors+=neg_log_prior(selpar_A50_rHB3_D,set_selpar_A50_rHB3_D(5), set_selpar_A50_rHB3_D(6), set_selpar_A50_rHB3_D(7));
  f_priors+=neg_log_prior(selpar_slope_rHB3_D,set_selpar_slope_rHB3_D(5), set_selpar_slope_rHB3_D(6), set_selpar_slope_rHB3_D(7));
  f_priors+=neg_log_prior(selpar_A502_rHB3_D,set_selpar_A502_rHB3_D(5), set_selpar_A502_rHB3_D(6), set_selpar_A502_rHB3_D(7));
  f_priors+=neg_log_prior(selpar_slope2_rHB3_D,set_selpar_slope2_rHB3_D(5), set_selpar_slope2_rHB3_D(6), set_selpar_slope2_rHB3_D(7));
  f_priors+=neg_log_prior(selpar_A50_rGN3_D,set_selpar_A50_rGN3_D(5), set_selpar_A50_rGN3_D(6), set_selpar_A50_rGN3_D(7));
  f_priors+=neg_log_prior(selpar_slope_rGN3_D,set_selpar_slope_rGN3_D(5), set_selpar_slope_rGN3_D(6), set_selpar_slope_rGN3_D(7));
  f_priors+=neg_log_prior(selpar_A502_rGN3_D,set_selpar_A502_rGN3_D(5), set_selpar_A502_rGN3_D(6), set_selpar_A502_rGN3_D(7));
  f_priors+=neg_log_prior(selpar_slope2_rGN3_D,set_selpar_slope2_rGN3_D(5), set_selpar_slope2_rGN3_D(6), set_selpar_slope2_rGN3_D(7));
  f_priors+=neg_log_prior(selpar_A50_sCT,set_selpar_A50_sCT(5), set_selpar_A50_sCT(6), set_selpar_A50_sCT(7));
  f_priors+=neg_log_prior(selpar_slope_sCT,set_selpar_slope_sCT(5), set_selpar_slope_sCT(6), set_selpar_slope_sCT(7));
  f_priors+=neg_log_prior(selpar_A502_sCT,set_selpar_A502_sCT(5), set_selpar_A502_sCT(6), set_selpar_A502_sCT(7));
  f_priors+=neg_log_prior(selpar_slope2_sCT,set_selpar_slope2_sCT(5), set_selpar_slope2_sCT(6), set_selpar_slope2_sCT(7));
  f_priors+=neg_log_prior(log_q_cpue_cHL,set_log_q_cpue_cHL(5),set_log_q_cpue_cHL(6),set_log_q_cpue_cHL(7));
  f_priors+=neg_log_prior(log_q_cpue_rHB,set_log_q_cpue_rHB(5),set_log_q_cpue_rHB(6),set_log_q_cpue_rHB(7));
  f_priors+=neg_log_prior(log_q_cpue_rHB_D,set_log_q_cpue_rHB_D(5),set_log_q_cpue_rHB_D(6),set_log_q_cpue_rHB_D(7));
  f_priors+=neg_log_prior(log_q_cpue_sCT,set_log_q_cpue_sCT(5),set_log_q_cpue_sCT(6),set_log_q_cpue_sCT(7));
  f_priors+=neg_log_prior(log_q_cpue_sVD,set_log_q_cpue_sVD(5),set_log_q_cpue_sVD(6),set_log_q_cpue_sVD(7));
  f_priors+=neg_log_prior(F_init,set_F_init(5),set_F_init(6),set_F_init(7));
  //f_priors+=neg_log_prior(log_avg_F_L_cHL,set_log_avg_F_L_cHL(5),set_log_avg_F_L_cHL(6),set_log_avg_F_L_cHL(7));
  //f_priors+=neg_log_prior(log_avg_F_cLL,set_log_avg_F_cLL(5),set_log_avg_F_cLL(6),set_log_avg_F_cLL(7));
  //f_priors+=neg_log_prior(log_avg_F_L_rHB,set_log_avg_F_L_rHB(5),set_log_avg_F_L_rHB(6),set_log_avg_F_L_rHB(7));
  //f_priors+=neg_log_prior(log_avg_F_L_rGN,set_log_avg_F_L_rGN(5),set_log_avg_F_L_rGN(6),set_log_avg_F_L_rGN(7));
  f_priors+=neg_log_prior(log_dm_lenc_cHL,set_log_dm_lenc_cHL(5),set_log_dm_lenc_cHL(6),set_log_dm_lenc_cHL(7));
  f_priors+=neg_log_prior(log_dm_lenc_cHL_D,set_log_dm_lenc_cHL_D(5),set_log_dm_lenc_cHL_D(6),set_log_dm_lenc_cHL_D(7));
  f_priors+=neg_log_prior(log_dm_lenc_rHB_D,set_log_dm_lenc_rHB_D(5),set_log_dm_lenc_rHB_D(6),set_log_dm_lenc_rHB_D(7));
  f_priors+=neg_log_prior(log_dm_lenc_rGN_D,set_log_dm_lenc_rGN_D(5),set_log_dm_lenc_rGN_D(6),set_log_dm_lenc_rGN_D(7));
  f_priors+=neg_log_prior(log_dm_agec_cHL,set_log_dm_agec_cHL(5),set_log_dm_agec_cHL(6),set_log_dm_agec_cHL(7));
  f_priors+=neg_log_prior(log_dm_agec_rHB,set_log_dm_agec_rHB(5),set_log_dm_agec_rHB(6),set_log_dm_agec_rHB(7));
  f_priors+=neg_log_prior(log_dm_agec_sCT,set_log_dm_agec_sCT(5),set_log_dm_agec_sCT(6),set_log_dm_agec_sCT(7));
  f_priors+=neg_log_prior(log_dm_agec_rGN,set_log_dm_agec_rGN(5),set_log_dm_agec_rGN(6),set_log_dm_agec_rGN(7));
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
  dvariable small_number=0.00001;
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
  dvariable small_number=0.00001;
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
  //nsamp=vector of N's, pred_comp=matrix of predicted comps, obs_comp=matrix of observed comps, ncomp = number of yrs in matrix, mbin=number of bins, minSS=min N threshold, log_dir_par=variance inflation factor
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
      get_weighted_current();
      cout<<"got weighted"<<endl;
      get_msy();
      cout<<"got msy"<<endl;
      get_per_recruit_stuff();
      cout<<"got per recruit"<<endl;  
      get_miscellaneous_stuff();
      cout<<"got misc stuff"<<endl;
	  get_effective_sample_sizes();
      cout<<"got effective samples sizes"<<endl;
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
      cout << "R0="<<R0<<endl;
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
      sdnr_lc_cHL_D=sdnr_multinomial(nyr_lenc_cHL_D, lenbins, nsamp_lenc_cHL_D, pred_cHL_D_lenc, obs_lenc_cHL_D, w_lenc_cHL_D); 
      sdnr_lc_rHB_D=sdnr_multinomial(nyr_lenc_rHB_D, lenbins, nsamp_lenc_rHB_D, pred_rHB_D_lenc, obs_lenc_rHB_D, w_lenc_rHB_D); 
	  sdnr_lc_rGN_D=sdnr_multinomial(nyr_lenc_rGN_D, lenbins, nsamp_lenc_rGN_D, pred_rGN_D_lenc, obs_lenc_rGN_D, w_lenc_rGN_D); 
      sdnr_ac_cHL=sdnr_multinomial(nyr_agec_cHL, agebins_agec, nsamp_agec_cHL, pred_cHL_agec, obs_agec_cHL, w_agec_cHL);  
      sdnr_ac_rHB=sdnr_multinomial(nyr_agec_rHB, agebins_agec_10p, nsamp_agec_rHB, pred_rHB_agec, obs_agec_rHB, w_agec_rHB);  
      sdnr_ac_sCT=sdnr_multinomial(nyr_agec_sCT, agebins_agec_10p, nsamp_agec_sCT, pred_sCT_agec, obs_agec_sCT, w_agec_sCT);  
      sdnr_ac_rGN=sdnr_multinomial(nyr_agec_rGN, agebins_agec_10p, nsamp_agec_rGN, pred_rGN_agec, obs_agec_rGN, w_agec_rGN);  
      sdnr_I_cHL=sdnr_lognormal(pred_cHL_cpue, obs_cpue_cHL, obs_cv_cpue_cHL, w_cpue_cHL);
      sdnr_I_rHB=sdnr_lognormal(pred_rHB_cpue, obs_cpue_rHB, obs_cv_cpue_rHB, w_cpue_rHB);
      sdnr_I_rHB_D=sdnr_lognormal(pred_rHB_D_cpue, obs_cpue_rHB_D, obs_cv_cpue_rHB_D, w_cpue_rHB_D);
      sdnr_I_sCT=sdnr_lognormal(pred_sCT_cpue, obs_cpue_sCT, obs_cv_cpue_sCT, w_cpue_sCT);  
      sdnr_I_sVD=sdnr_lognormal(pred_sVD_cpue, obs_cpue_sVD, obs_cv_cpue_sVD, w_cpue_sVD);   
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
	   Linf_20_out(8)=Linf_20; Linf_20_out(1,7)=set_Linf_20; 
       K_20_out(8)=K_20; K_20_out(1,7)=set_K_20;
       t0_20_out(8)=t0_20; t0_20_out(1,7)=set_t0_20;
       len_cv_val_20_out(8)=len_cv_val_20; len_cv_val_20_out(1,7)=set_len_cv_20;
       log_R0_out(8)=log_R0; log_R0_out(1,7)=set_log_R0;
       M_constant_out(8)=M_constant; M_constant_out(1,7)=set_M_constant;
       steep_out(8)=steep; steep_out(1,7)=set_steep;
       rec_sigma_out(8)=rec_sigma; rec_sigma_out(1,7)=set_rec_sigma;
       R_autocorr_out(8)=R_autocorr; R_autocorr_out(1,7)=set_R_autocorr;
       log_dm_cHL_lc_out(8)=log_dm_lenc_cHL; log_dm_cHL_lc_out(1,7)=set_log_dm_lenc_cHL;
	   log_dm_cHL_D_lc_out(8)=log_dm_lenc_cHL_D; log_dm_cHL_D_lc_out(1,7)=set_log_dm_lenc_cHL_D;
	   log_dm_rHB_D_lc_out(8)=log_dm_lenc_rHB_D; log_dm_rHB_D_lc_out(1,7)=set_log_dm_lenc_rHB_D;
	   log_dm_rGN_D_lc_out(8)=log_dm_lenc_rGN_D; log_dm_rGN_D_lc_out(1,7)=set_log_dm_lenc_rGN_D;
	   log_dm_cHL_ac_out(8)=log_dm_agec_cHL; log_dm_cHL_ac_out(1,7)=set_log_dm_agec_cHL;  
	   log_dm_rHB_ac_out(8)=log_dm_agec_rHB; log_dm_rHB_ac_out(1,7)=set_log_dm_agec_rHB;  
	   log_dm_sCT_ac_out(8)=log_dm_agec_sCT; log_dm_sCT_ac_out(1,7)=set_log_dm_agec_sCT;  
	   log_dm_rGN_ac_out(8)=log_dm_agec_rGN; log_dm_rGN_ac_out(1,7)=set_log_dm_agec_rGN;  
       selpar_A50_cHL1_out(8)=selpar_A50_cHL1; selpar_A50_cHL1_out(1,7)=set_selpar_A50_cHL1;
       selpar_slope_cHL1_out(8)=selpar_slope_cHL1; selpar_slope_cHL1_out(1,7)=set_selpar_slope_cHL1;
       selpar_A50_cHL2_out(8)=selpar_A50_cHL2; selpar_A50_cHL2_out(1,7)=set_selpar_A50_cHL2;
       selpar_slope_cHL2_out(8)=selpar_slope_cHL2; selpar_slope_cHL2_out(1,7)=set_selpar_slope_cHL2;
       selpar_A50_cHL3_out(8)=selpar_A50_cHL3; selpar_A50_cHL3_out(1,7)=set_selpar_A50_cHL3;
       selpar_slope_cHL3_out(8)=selpar_slope_cHL3; selpar_slope_cHL3_out(1,7)=set_selpar_slope_cHL3;
       selpar_A50_rHB1_out(8)=selpar_A50_rHB1; selpar_A50_rHB1_out(1,7)=set_selpar_A50_rHB1;
       selpar_slope_rHB1_out(8)=selpar_slope_rHB1; selpar_slope_rHB1_out(1,7)=set_selpar_slope_rHB1;
       selpar_A502_rHB1_out(8)=selpar_A502_rHB1; selpar_A502_rHB1_out(1,7)=set_selpar_A502_rHB1;
       selpar_slope2_rHB1_out(8)=selpar_slope2_rHB1; selpar_slope2_rHB1_out(1,7)=set_selpar_slope2_rHB1;
       selpar_A50_rHB2_out(8)=selpar_A50_rHB2; selpar_A50_rHB2_out(1,7)=set_selpar_A50_rHB2;
       selpar_slope_rHB2_out(8)=selpar_slope_rHB2; selpar_slope_rHB2_out(1,7)=set_selpar_slope_rHB2;
       selpar_A502_rHB2_out(8)=selpar_A502_rHB2; selpar_A502_rHB2_out(1,7)=set_selpar_A502_rHB2;
       selpar_slope2_rHB2_out(8)=selpar_slope2_rHB2; selpar_slope2_rHB2_out(1,7)=set_selpar_slope2_rHB2;
       selpar_A50_rHB3_out(8)=selpar_A50_rHB3; selpar_A50_rHB3_out(1,7)=set_selpar_A50_rHB3;
       selpar_slope_rHB3_out(8)=selpar_slope_rHB3; selpar_slope_rHB3_out(1,7)=set_selpar_slope_rHB3;
       selpar_A502_rHB3_out(8)=selpar_A502_rHB3; selpar_A502_rHB3_out(1,7)=set_selpar_A502_rHB3;
       selpar_slope2_rHB3_out(8)=selpar_slope2_rHB3; selpar_slope2_rHB3_out(1,7)=set_selpar_slope2_rHB3;
       selpar_A50_rHB2_D_out(8)=selpar_A50_rHB2_D; selpar_A50_rHB2_D_out(1,7)=set_selpar_A50_rHB2_D;
       selpar_slope_rHB2_D_out(8)=selpar_slope_rHB2_D; selpar_slope_rHB2_D_out(1,7)=set_selpar_slope_rHB2_D;
       selpar_A502_rHB2_D_out(8)=selpar_A502_rHB2_D; selpar_A502_rHB2_D_out(1,7)=set_selpar_A502_rHB2_D;
       selpar_slope2_rHB2_D_out(8)=selpar_slope2_rHB2_D; selpar_slope2_rHB2_D_out(1,7)=set_selpar_slope2_rHB2_D;
	   selpar_A50_rHB3_D_out(8)=selpar_A50_rHB3_D; selpar_A50_rHB3_D_out(1,7)=set_selpar_A50_rHB3_D;
	   selpar_slope_rHB3_D_out(8)=selpar_slope_rHB3_D; selpar_slope_rHB3_D_out(1,7)=set_selpar_slope_rHB3_D;
	   selpar_A502_rHB3_D_out(8)=selpar_A502_rHB3_D; selpar_A502_rHB3_D_out(1,7)=set_selpar_A502_rHB3_D;
	   selpar_slope2_rHB3_D_out(8)=selpar_slope2_rHB3_D; selpar_slope2_rHB3_D_out(1,7)=set_selpar_slope2_rHB3_D;
	   selpar_A50_rGN3_D_out(8)=selpar_A50_rGN3_D; selpar_A50_rGN3_D_out(1,7)=set_selpar_A50_rGN3_D;
	   selpar_slope_rGN3_D_out(8)=selpar_slope_rGN3_D; selpar_slope_rGN3_D_out(1,7)=set_selpar_slope_rGN3_D;
	   selpar_A502_rGN3_D_out(8)=selpar_A502_rGN3_D; selpar_A502_rGN3_D_out(1,7)=set_selpar_A502_rGN3_D;
	   selpar_slope2_rGN3_D_out(8)=selpar_slope2_rGN3_D; selpar_slope2_rGN3_D_out(1,7)=set_selpar_slope2_rGN3_D;	   
	   selpar_A50_cHL2_D_out(8)=selpar_A50_cHL2_D; selpar_A50_cHL2_D_out(1,7)=set_selpar_A50_cHL2_D;
	   selpar_slope_cHL2_D_out(8)=selpar_slope_cHL2_D; selpar_slope_cHL2_D_out(1,7)=set_selpar_slope_cHL2_D;
	   selpar_A502_cHL2_D_out(8)=selpar_A502_cHL2_D; selpar_A502_cHL2_D_out(1,7)=set_selpar_A502_cHL2_D;
	   selpar_slope2_cHL2_D_out(8)=selpar_slope2_cHL2_D; selpar_slope2_cHL2_D_out(1,7)=set_selpar_slope2_cHL2_D;
	   selpar_A50_cHL3_D_out(8)=selpar_A50_cHL3_D; selpar_A50_cHL3_D_out(1,7)=set_selpar_A50_cHL3_D;
	   selpar_slope_cHL3_D_out(8)=selpar_slope_cHL3_D; selpar_slope_cHL3_D_out(1,7)=set_selpar_slope_cHL3_D;
	   selpar_A50_rGN2_out(8)=selpar_A50_rGN2; selpar_A50_rGN2_out(1,7)=set_selpar_A50_rGN2;
       selpar_slope_rGN2_out(8)=selpar_slope_rGN2; selpar_slope_rGN2_out(1,7)=set_selpar_slope_rGN2;
       selpar_A502_rGN2_out(8)=selpar_A502_rGN2; selpar_A502_rGN2_out(1,7)=set_selpar_A502_rGN2;
       selpar_slope2_rGN2_out(8)=selpar_slope2_rGN2; selpar_slope2_rGN2_out(1,7)=set_selpar_slope2_rGN2;
       selpar_A50_rGN3_out(8)=selpar_A50_rGN3; selpar_A50_rGN3_out(1,7)=set_selpar_A50_rGN3;
       selpar_slope_rGN3_out(8)=selpar_slope_rGN3; selpar_slope_rGN3_out(1,7)=set_selpar_slope_rGN3;
       //selpar_A502_rGN3_out(8)=selpar_A502_rGN3; selpar_A502_rGN3_out(1,7)=set_selpar_A502_rGN3;
       //selpar_slope2_rGN3_out(8)=selpar_slope2_rGN3; selpar_slope2_rGN3_out(1,7)=set_selpar_slope2_rGN3;
       selpar_A50_sCT_out(8)=selpar_A50_sCT; selpar_A50_sCT_out(1,7)=set_selpar_A50_sCT;
       selpar_slope_sCT_out(8)=selpar_slope_sCT; selpar_slope_sCT_out(1,7)=set_selpar_slope_sCT;
	   selpar_A502_sCT_out(8)=selpar_A502_sCT; selpar_A502_sCT_out(1,7)=set_selpar_A502_sCT;
       selpar_slope2_sCT_out(8)=selpar_slope2_sCT; selpar_slope2_sCT_out(1,7)=set_selpar_slope2_sCT;
       log_q_cHL_out(8)=log_q_cpue_cHL; log_q_cHL_out(1,7)=set_log_q_cpue_cHL;
       log_q_rHB_out(8)=log_q_cpue_rHB; log_q_rHB_out(1,7)=set_log_q_cpue_rHB;
       log_q_rHB_D_out(8)=log_q_cpue_rHB_D; log_q_rHB_D_out(1,7)=set_log_q_cpue_rHB_D;
       log_q_sCT_out(8)=log_q_cpue_sCT; log_q_sCT_out(1,7)=set_log_q_cpue_sCT;
	   log_q_sVD_out(8)=log_q_cpue_sVD; log_q_sVD_out(1,7)=set_log_q_cpue_sVD;
       log_avg_F_cHL_out(8)=log_avg_F_L_cHL; log_avg_F_cHL_out(1,7)=set_log_avg_F_L_cHL;
       log_avg_F_rHB_out(8)=log_avg_F_L_rHB; log_avg_F_rHB_out(1,7)=set_log_avg_F_L_rHB;
       log_avg_F_rGN_out(8)=log_avg_F_L_rGN; log_avg_F_rGN_out(1,7)=set_log_avg_F_L_rGN;       
       log_avg_F_cHL_D_out(8)=log_avg_F_D_cHL; log_avg_F_cHL_D_out(1,7)=set_log_avg_F_D_cHL;
       log_avg_F_rHB_D_out(8)=log_avg_F_D_rHB; log_avg_F_rHB_D_out(1,7)=set_log_avg_F_D_rHB;
       log_avg_F_rGN_D_out(8)=log_avg_F_D_rGN; log_avg_F_rGN_D_out(1,7)=set_log_avg_F_D_rGN;
       F_init_out(8)=F_init; F_init_out(1,7)=set_F_init;
       log_rec_dev_out(styr_rec_dev, endyr_rec_dev)=log_dev_rec;
       log_F_dev_cHL_out(styr_L_cHL,endyr_L_cHL)=log_dev_F_L_cHL;
       log_F_dev_rHB_out(styr_L_rHB,endyr_L_rHB)=log_dev_F_L_rHB;
       log_F_dev_rGN_out(styr_L_rGN,endyr_L_rGN)=log_dev_F_L_rGN;
       log_F_dev_cHL_D_out(styr_D_cHL,endyr_D_cHL)=log_dev_F_D_cHL;
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
