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
  endyr_selex_phase3.allocate("endyr_selex_phase3");
  endyr_selex_phase4.allocate("endyr_selex_phase4");
  endyr_selex_phase5.allocate("endyr_selex_phase5");
  endyr_selex_phase6.allocate("endyr_selex_phase6");
  endyr_selex_phase7.allocate("endyr_selex_phase7");
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
  styr_cpue_nad.allocate("styr_cpue_nad");
  endyr_cpue_nad.allocate("endyr_cpue_nad");
  obs_cpue_nad.allocate(styr_cpue_nad,endyr_cpue_nad,"obs_cpue_nad");
  obs_cv_cpue_nad.allocate(styr_cpue_nad,endyr_cpue_nad,"obs_cv_cpue_nad");
  nyr_lenc_nad.allocate("nyr_lenc_nad");
  yrs_lenc_nad.allocate(1,nyr_lenc_nad,"yrs_lenc_nad");
  nsamp_lenc_nad.allocate(1,nyr_lenc_nad,"nsamp_lenc_nad");
  nfish_lenc_nad.allocate(1,nyr_lenc_nad,"nfish_lenc_nad");
  obs_lenc_nad.allocate(1,nyr_lenc_nad,1,nlenbins,"obs_lenc_nad");
  styr_cpue_mad.allocate("styr_cpue_mad");
  endyr_cpue_mad.allocate("endyr_cpue_mad");
  obs_cpue_mad.allocate(styr_cpue_mad,endyr_cpue_mad,"obs_cpue_mad");
  obs_cv_cpue_mad.allocate(styr_cpue_mad,endyr_cpue_mad,"obs_cv_cpue_mad");
  nyr_lenc_mad.allocate("nyr_lenc_mad");
  yrs_lenc_mad.allocate(1,nyr_lenc_mad,"yrs_lenc_mad");
  nsamp_lenc_mad.allocate(1,nyr_lenc_mad,"nsamp_lenc_mad");
  nfish_lenc_mad.allocate(1,nyr_lenc_mad,"nfish_lenc_mad");
  obs_lenc_mad.allocate(1,nyr_lenc_mad,1,nlenbins,"obs_lenc_mad");
  styr_cpue_sad.allocate("styr_cpue_sad");
  endyr_cpue_sad.allocate("endyr_cpue_sad");
  obs_cpue_sad.allocate(styr_cpue_sad,endyr_cpue_sad,"obs_cpue_sad");
  obs_cv_cpue_sad.allocate(styr_cpue_sad,endyr_cpue_sad,"obs_cv_cpue_sad");
  styr_cpue_jai.allocate("styr_cpue_jai");
  endyr_cpue_jai.allocate("endyr_cpue_jai");
  obs_cpue_jai.allocate(styr_cpue_jai,endyr_cpue_jai,"obs_cpue_jai");
  obs_cv_cpue_jai.allocate(styr_cpue_jai,endyr_cpue_jai,"obs_cv_cpue_jai");
  yr_q_change.allocate("yr_q_change");
  nyr_cpue_mareco.allocate("nyr_cpue_mareco");
  yrs_cpue_mareco.allocate(1,nyr_cpue_mareco,"yrs_cpue_mareco");
  obs_cpue_mareco.allocate(1,nyr_cpue_mareco,"obs_cpue_mareco");
  obs_cv_cpue_mareco.allocate(1,nyr_cpue_mareco,"obs_cv_cpue_mareco");
  styr_L_cRn.allocate("styr_L_cRn");
  endyr_L_cRn.allocate("endyr_L_cRn");
  obs_L_cRn.allocate(styr_L_cRn,endyr_L_cRn,"obs_L_cRn");
  obs_cv_L_cRn.allocate(styr_L_cRn,endyr_L_cRn,"obs_cv_L_cRn");
  nyr_agec_cRn.allocate("nyr_agec_cRn");
  yrs_agec_cRn.allocate(1,nyr_agec_cRn,"yrs_agec_cRn");
  nsamp_agec_cRn.allocate(1,nyr_agec_cRn,"nsamp_agec_cRn");
  nfish_agec_cRn.allocate(1,nyr_agec_cRn,"nfish_agec_cRn");
  obs_agec_cRn.allocate(1,nyr_agec_cRn,1,nages_agec,"obs_agec_cRn");
  styr_L_cRs.allocate("styr_L_cRs");
  endyr_L_cRs.allocate("endyr_L_cRs");
  obs_L_cRs.allocate(styr_L_cRs,endyr_L_cRs,"obs_L_cRs");
  obs_cv_L_cRs.allocate(styr_L_cRs,endyr_L_cRs,"obs_cv_L_cRs");
  nyr_agec_cRs.allocate("nyr_agec_cRs");
  yrs_agec_cRs.allocate(1,nyr_agec_cRs,"yrs_agec_cRs");
  nsamp_agec_cRs.allocate(1,nyr_agec_cRs,"nsamp_agec_cRs");
  nfish_agec_cRs.allocate(1,nyr_agec_cRs,"nfish_agec_cRs");
  obs_agec_cRs.allocate(1,nyr_agec_cRs,1,nages_agec,"obs_agec_cRs");
  styr_L_cBn.allocate("styr_L_cBn");
  endyr_L_cBn.allocate("endyr_L_cBn");
  obs_L_cBn.allocate(styr_L_cBn,endyr_L_cBn,"obs_L_cBn");
  obs_cv_L_cBn.allocate(styr_L_cBn,endyr_L_cBn,"obs_cv_L_cBn");
  nyr_agec_cBn.allocate("nyr_agec_cBn");
  yrs_agec_cBn.allocate(1,nyr_agec_cBn,"yrs_agec_cBn");
  nsamp_agec_cBn.allocate(1,nyr_agec_cBn,"nsamp_agec_cBn");
  nfish_agec_cBn.allocate(1,nyr_agec_cBn,"nfish_agec_cBn");
  obs_agec_cBn.allocate(1,nyr_agec_cBn,1,nages_agec,"obs_agec_cBn");
  styr_L_cBs.allocate("styr_L_cBs");
  endyr_L_cBs.allocate("endyr_L_cBs");
  obs_L_cBs.allocate(styr_L_cBs,endyr_L_cBs,"obs_L_cBs");
  obs_cv_L_cBs.allocate(styr_L_cBs,endyr_L_cBs,"obs_cv_L_cBs");
  nyr_agec_cBs.allocate("nyr_agec_cBs");
  yrs_agec_cBs.allocate(1,nyr_agec_cBs,"yrs_agec_cBs");
  nsamp_agec_cBs.allocate(1,nyr_agec_cBs,"nsamp_agec_cBs");
  nfish_agec_cBs.allocate(1,nyr_agec_cBs,"nfish_agec_cBs");
  obs_agec_cBs.allocate(1,nyr_agec_cBs,1,nages_agec,"obs_agec_cBs");
  set_Linf.allocate(1,7,"set_Linf");
  set_K.allocate(1,7,"set_K");
  set_t0.allocate(1,7,"set_t0");
  set_len_cv_nad.allocate(1,7,"set_len_cv_nad");
  set_len_cv_mad.allocate(1,7,"set_len_cv_mad");
  set_steep.allocate(1,7,"set_steep");
  set_log_R0.allocate(1,7,"set_log_R0");
  set_R_autocorr.allocate(1,7,"set_R_autocorr");
  set_rec_sigma.allocate(1,7,"set_rec_sigma");
  set_log_dm_lenc_nad.allocate(1,7,"set_log_dm_lenc_nad");
  set_log_dm_lenc_mad.allocate(1,7,"set_log_dm_lenc_mad");
  set_log_dm_agec_cRn.allocate(1,7,"set_log_dm_agec_cRn");
  set_log_dm_agec_cRs.allocate(1,7,"set_log_dm_agec_cRs");
  set_log_dm_agec_cBn.allocate(1,7,"set_log_dm_agec_cBn");
  set_log_dm_agec_cBs.allocate(1,7,"set_log_dm_agec_cBs");
  set_selpar_A50_cRn.allocate(1,7,"set_selpar_A50_cRn");
  set_selpar_slope_cRn.allocate(1,7,"set_selpar_slope_cRn");
  set_selpar_A502_cRn.allocate(1,7,"set_selpar_A502_cRn");
  set_selpar_slope2_cRn.allocate(1,7,"set_selpar_slope2_cRn");
  set_selpar_A50_cRn2.allocate(1,7,"set_selpar_A50_cRn2");
  set_selpar_slope_cRn2.allocate(1,7,"set_selpar_slope_cRn2");
  set_selpar_A502_cRn2.allocate(1,7,"set_selpar_A502_cRn2");
  set_selpar_slope2_cRn2.allocate(1,7,"set_selpar_slope2_cRn2");
  set_selpar_A50_cRn3.allocate(1,7,"set_selpar_A50_cRn3");
  set_selpar_slope_cRn3.allocate(1,7,"set_selpar_slope_cRn3");
  set_selpar_A502_cRn3.allocate(1,7,"set_selpar_A502_cRn3");
  set_selpar_slope2_cRn3.allocate(1,7,"set_selpar_slope2_cRn3");
  set_selpar_A50_cRn4.allocate(1,7,"set_selpar_A50_cRn4");
  set_selpar_slope_cRn4.allocate(1,7,"set_selpar_slope_cRn4");
  set_selpar_A502_cRn4.allocate(1,7,"set_selpar_A502_cRn4");
  set_selpar_slope2_cRn4.allocate(1,7,"set_selpar_slope2_cRn4");
  set_selpar_A50_cRs.allocate(1,7,"set_selpar_A50_cRs");
  set_selpar_slope_cRs.allocate(1,7,"set_selpar_slope_cRs");
  set_selpar_A502_cRs.allocate(1,7,"set_selpar_A502_cRs");
  set_selpar_slope2_cRs.allocate(1,7,"set_selpar_slope2_cRs");
  set_selpar_A50_cRs2.allocate(1,7,"set_selpar_A50_cRs2");
  set_selpar_slope_cRs2.allocate(1,7,"set_selpar_slope_cRs2");
  set_selpar_A502_cRs2.allocate(1,7,"set_selpar_A502_cRs2");
  set_selpar_slope2_cRs2.allocate(1,7,"set_selpar_slope2_cRs2");
  set_selpar_A50_cRs3.allocate(1,7,"set_selpar_A50_cRs3");
  set_selpar_slope_cRs3.allocate(1,7,"set_selpar_slope_cRs3");
  set_selpar_A502_cRs3.allocate(1,7,"set_selpar_A502_cRs3");
  set_selpar_slope2_cRs3.allocate(1,7,"set_selpar_slope2_cRs3");
  set_selpar_A50_cRs4.allocate(1,7,"set_selpar_A50_cRs4");
  set_selpar_slope_cRs4.allocate(1,7,"set_selpar_slope_cRs4");
  set_selpar_A502_cRs4.allocate(1,7,"set_selpar_A502_cRs4");
  set_selpar_slope2_cRs4.allocate(1,7,"set_selpar_slope2_cRs4");
  set_selpar_A50_cBn.allocate(1,7,"set_selpar_A50_cBn");
  set_selpar_slope_cBn.allocate(1,7,"set_selpar_slope_cBn");
  set_selpar_A502_cBn.allocate(1,7,"set_selpar_A502_cBn");
  set_selpar_slope2_cBn.allocate(1,7,"set_selpar_slope2_cBn");
  set_selpar_A50_cBn2.allocate(1,7,"set_selpar_A50_cBn2");
  set_selpar_slope_cBn2.allocate(1,7,"set_selpar_slope_cBn2");
  set_selpar_A502_cBn2.allocate(1,7,"set_selpar_A502_cBn2");
  set_selpar_slope2_cBn2.allocate(1,7,"set_selpar_slope2_cBn2");
  set_selpar_A50_cBs.allocate(1,7,"set_selpar_A50_cBs");
  set_selpar_slope_cBs.allocate(1,7,"set_selpar_slope_cBs");
  set_selpar_A502_cBs.allocate(1,7,"set_selpar_A502_cBs");
  set_selpar_slope2_cBs.allocate(1,7,"set_selpar_slope2_cBs");
  set_selpar_A50_cBs2.allocate(1,7,"set_selpar_A50_cBs2");
  set_selpar_slope_cBs2.allocate(1,7,"set_selpar_slope_cBs2");
  set_selpar_A502_cBs2.allocate(1,7,"set_selpar_A502_cBs2");
  set_selpar_slope2_cBs2.allocate(1,7,"set_selpar_slope2_cBs2");
  set_sel_age0_cRn.allocate(1,7,"set_sel_age0_cRn");
  set_sel_age1_cRn.allocate(1,7,"set_sel_age1_cRn");
  set_sel_age2_cRn.allocate(1,7,"set_sel_age2_cRn");
  set_sel_age3_cRn.allocate(1,7,"set_sel_age3_cRn");
  set_sel_age4_cRn.allocate(1,7,"set_sel_age4_cRn");
  set_sel_age5_cRn.allocate(1,7,"set_sel_age5_cRn");
  set_sel_age6_cRn.allocate(1,7,"set_sel_age6_cRn");
  set_sel_age0_cRn2.allocate(1,7,"set_sel_age0_cRn2");
  set_sel_age1_cRn2.allocate(1,7,"set_sel_age1_cRn2");
  set_sel_age2_cRn2.allocate(1,7,"set_sel_age2_cRn2");
  set_sel_age3_cRn2.allocate(1,7,"set_sel_age3_cRn2");
  set_sel_age4_cRn2.allocate(1,7,"set_sel_age4_cRn2");
  set_sel_age5_cRn2.allocate(1,7,"set_sel_age5_cRn2");
  set_sel_age6_cRn2.allocate(1,7,"set_sel_age6_cRn2");
  set_sel_age0_cRn3.allocate(1,7,"set_sel_age0_cRn3");
  set_sel_age1_cRn3.allocate(1,7,"set_sel_age1_cRn3");
  set_sel_age2_cRn3.allocate(1,7,"set_sel_age2_cRn3");
  set_sel_age3_cRn3.allocate(1,7,"set_sel_age3_cRn3");
  set_sel_age4_cRn3.allocate(1,7,"set_sel_age4_cRn3");
  set_sel_age5_cRn3.allocate(1,7,"set_sel_age5_cRn3");
  set_sel_age6_cRn3.allocate(1,7,"set_sel_age6_cRn3");
  set_sel_age0_cRn4.allocate(1,7,"set_sel_age0_cRn4");
  set_sel_age1_cRn4.allocate(1,7,"set_sel_age1_cRn4");
  set_sel_age2_cRn4.allocate(1,7,"set_sel_age2_cRn4");
  set_sel_age3_cRn4.allocate(1,7,"set_sel_age3_cRn4");
  set_sel_age4_cRn4.allocate(1,7,"set_sel_age4_cRn4");
  set_sel_age5_cRn4.allocate(1,7,"set_sel_age5_cRn4");
  set_sel_age6_cRn4.allocate(1,7,"set_sel_age6_cRn4");
  set_sel_age0_cRs.allocate(1,7,"set_sel_age0_cRs");
  set_sel_age1_cRs.allocate(1,7,"set_sel_age1_cRs");
  set_sel_age2_cRs.allocate(1,7,"set_sel_age2_cRs");
  set_sel_age3_cRs.allocate(1,7,"set_sel_age3_cRs");
  set_sel_age4_cRs.allocate(1,7,"set_sel_age4_cRs");
  set_sel_age5_cRs.allocate(1,7,"set_sel_age5_cRs");
  set_sel_age6_cRs.allocate(1,7,"set_sel_age6_cRs");
  set_sel_age0_cRs2.allocate(1,7,"set_sel_age0_cRs2");
  set_sel_age1_cRs2.allocate(1,7,"set_sel_age1_cRs2");
  set_sel_age2_cRs2.allocate(1,7,"set_sel_age2_cRs2");
  set_sel_age3_cRs2.allocate(1,7,"set_sel_age3_cRs2");
  set_sel_age4_cRs2.allocate(1,7,"set_sel_age4_cRs2");
  set_sel_age5_cRs2.allocate(1,7,"set_sel_age5_cRs2");
  set_sel_age6_cRs2.allocate(1,7,"set_sel_age6_cRs2");
  set_sel_age0_cRs3.allocate(1,7,"set_sel_age0_cRs3");
  set_sel_age1_cRs3.allocate(1,7,"set_sel_age1_cRs3");
  set_sel_age2_cRs3.allocate(1,7,"set_sel_age2_cRs3");
  set_sel_age3_cRs3.allocate(1,7,"set_sel_age3_cRs3");
  set_sel_age4_cRs3.allocate(1,7,"set_sel_age4_cRs3");
  set_sel_age5_cRs3.allocate(1,7,"set_sel_age5_cRs3");
  set_sel_age6_cRs3.allocate(1,7,"set_sel_age6_cRs3");
  set_sel_age0_cRs4.allocate(1,7,"set_sel_age0_cRs4");
  set_sel_age1_cRs4.allocate(1,7,"set_sel_age1_cRs4");
  set_sel_age2_cRs4.allocate(1,7,"set_sel_age2_cRs4");
  set_sel_age3_cRs4.allocate(1,7,"set_sel_age3_cRs4");
  set_sel_age4_cRs4.allocate(1,7,"set_sel_age4_cRs4");
  set_sel_age5_cRs4.allocate(1,7,"set_sel_age5_cRs4");
  set_sel_age6_cRs4.allocate(1,7,"set_sel_age6_cRs4");
  set_sel_age0_cBn.allocate(1,7,"set_sel_age0_cBn");
  set_sel_age1_cBn.allocate(1,7,"set_sel_age1_cBn");
  set_sel_age2_cBn.allocate(1,7,"set_sel_age2_cBn");
  set_sel_age3_cBn.allocate(1,7,"set_sel_age3_cBn");
  set_sel_age4_cBn.allocate(1,7,"set_sel_age4_cBn");
  set_sel_age5_cBn.allocate(1,7,"set_sel_age5_cBn");
  set_sel_age6_cBn.allocate(1,7,"set_sel_age6_cBn");
  set_sel_age0_cBn2.allocate(1,7,"set_sel_age0_cBn2");
  set_sel_age1_cBn2.allocate(1,7,"set_sel_age1_cBn2");
  set_sel_age2_cBn2.allocate(1,7,"set_sel_age2_cBn2");
  set_sel_age3_cBn2.allocate(1,7,"set_sel_age3_cBn2");
  set_sel_age4_cBn2.allocate(1,7,"set_sel_age4_cBn2");
  set_sel_age5_cBn2.allocate(1,7,"set_sel_age5_cBn2");
  set_sel_age6_cBn2.allocate(1,7,"set_sel_age6_cBn2");
  set_sel_age0_cBs.allocate(1,7,"set_sel_age0_cBs");
  set_sel_age1_cBs.allocate(1,7,"set_sel_age1_cBs");
  set_sel_age2_cBs.allocate(1,7,"set_sel_age2_cBs");
  set_sel_age3_cBs.allocate(1,7,"set_sel_age3_cBs");
  set_sel_age4_cBs.allocate(1,7,"set_sel_age4_cBs");
  set_sel_age5_cBs.allocate(1,7,"set_sel_age5_cBs");
  set_sel_age6_cBs.allocate(1,7,"set_sel_age6_cBs");
  set_sel_age0_cBs2.allocate(1,7,"set_sel_age0_cBs2");
  set_sel_age1_cBs2.allocate(1,7,"set_sel_age1_cBs2");
  set_sel_age2_cBs2.allocate(1,7,"set_sel_age2_cBs2");
  set_sel_age3_cBs2.allocate(1,7,"set_sel_age3_cBs2");
  set_sel_age4_cBs2.allocate(1,7,"set_sel_age4_cBs2");
  set_sel_age5_cBs2.allocate(1,7,"set_sel_age5_cBs2");
  set_sel_age6_cBs2.allocate(1,7,"set_sel_age6_cBs2");
  set_selpar_A50_nad.allocate(1,7,"set_selpar_A50_nad");
  set_selpar_slope_nad.allocate(1,7,"set_selpar_slope_nad");
  set_selpar_A502_nad.allocate(1,7,"set_selpar_A502_nad");
  set_selpar_slope2_nad.allocate(1,7,"set_selpar_slope2_nad");
  set_selpar_A50_mad.allocate(1,7,"set_selpar_A50_mad");
  set_selpar_slope_mad.allocate(1,7,"set_selpar_slope_mad");
  set_selpar_A502_mad.allocate(1,7,"set_selpar_A502_mad");
  set_selpar_slope2_mad.allocate(1,7,"set_selpar_slope2_mad");
  set_selpar_A50_sad.allocate(1,7,"set_selpar_A50_sad");
  set_selpar_slope_sad.allocate(1,7,"set_selpar_slope_sad");
  set_selpar_A502_sad.allocate(1,7,"set_selpar_A502_sad");
  set_selpar_slope2_sad.allocate(1,7,"set_selpar_slope2_sad");
  set_sel_age0_nad.allocate(1,7,"set_sel_age0_nad");
  set_sel_age1_nad.allocate(1,7,"set_sel_age1_nad");
  set_sel_age2_nad.allocate(1,7,"set_sel_age2_nad");
  set_sel_age3_nad.allocate(1,7,"set_sel_age3_nad");
  set_sel_age4_nad.allocate(1,7,"set_sel_age4_nad");
  set_sel_age5_nad.allocate(1,7,"set_sel_age5_nad");
  set_sel_age6_nad.allocate(1,7,"set_sel_age6_nad");
  set_sel_age0_mad.allocate(1,7,"set_sel_age0_mad");
  set_sel_age1_mad.allocate(1,7,"set_sel_age1_mad");
  set_sel_age2_mad.allocate(1,7,"set_sel_age2_mad");
  set_sel_age3_mad.allocate(1,7,"set_sel_age3_mad");
  set_sel_age4_mad.allocate(1,7,"set_sel_age4_mad");
  set_sel_age5_mad.allocate(1,7,"set_sel_age5_mad");
  set_sel_age6_mad.allocate(1,7,"set_sel_age6_mad");
  set_sel_age0_sad.allocate(1,7,"set_sel_age0_sad");
  set_sel_age1_sad.allocate(1,7,"set_sel_age1_sad");
  set_sel_age2_sad.allocate(1,7,"set_sel_age2_sad");
  set_sel_age3_sad.allocate(1,7,"set_sel_age3_sad");
  set_sel_age4_sad.allocate(1,7,"set_sel_age4_sad");
  set_sel_age5_sad.allocate(1,7,"set_sel_age5_sad");
  set_sel_age6_sad.allocate(1,7,"set_sel_age6_sad");
  set_log_q_cpue_nad.allocate(1,7,"set_log_q_cpue_nad");
  set_log_q_cpue_mad.allocate(1,7,"set_log_q_cpue_mad");
  set_log_q_cpue_sad.allocate(1,7,"set_log_q_cpue_sad");
  set_log_q_cpue_jai.allocate(1,7,"set_log_q_cpue_jai");
  set_log_q2_jai.allocate(1,7,"set_log_q2_jai");
  set_log_q_cpue_mar.allocate(1,7,"set_log_q_cpue_mar");
  set_log_q_cpue_eco.allocate(1,7,"set_log_q_cpue_eco");
  set_log_avg_F_L_cRn.allocate(1,7,"set_log_avg_F_L_cRn");
  set_log_avg_F_L_cRs.allocate(1,7,"set_log_avg_F_L_cRs");
  set_log_avg_F_L_cBn.allocate(1,7,"set_log_avg_F_L_cBn");
  set_log_avg_F_L_cBs.allocate(1,7,"set_log_avg_F_L_cBs");
  set_log_dev_F_L_cRn.allocate(1,3,"set_log_dev_F_L_cRn");
  set_log_dev_F_L_cRs.allocate(1,3,"set_log_dev_F_L_cRs");
  set_log_dev_F_L_cBn.allocate(1,3,"set_log_dev_F_L_cBn");
  set_log_dev_F_L_cBs.allocate(1,3,"set_log_dev_F_L_cBs");
  set_log_dev_RWq.allocate(1,3,"set_log_dev_RWq");
  set_log_dev_rec.allocate(1,3,"set_log_dev_rec");
  set_log_dev_Nage.allocate(1,3,"set_log_dev_Nage");
  set_log_dev_vals_F_L__cRn.allocate(styr_L_cRn,endyr_L_cRn,"set_log_dev_vals_F_L__cRn");
  set_log_dev_vals_F_L__cRs.allocate(styr_L_cRs,endyr_L_cRs,"set_log_dev_vals_F_L__cRs");
  set_log_dev_vals_F_L__cBn.allocate(styr_L_cBn,endyr_L_cBn,"set_log_dev_vals_F_L__cBn");
  set_log_dev_vals_F_L__cBs.allocate(styr_L_cBs,endyr_L_cBs,"set_log_dev_vals_F_L__cBs");
  set_log_dev_vals_rec.allocate(styr_rec_dev,endyr_rec_dev,"set_log_dev_vals_rec");
  set_log_dev_vals_Nage.allocate(2,nages,"set_log_dev_vals_Nage");
  set_w_L.allocate("set_w_L");
  set_w_cpue_nad.allocate("set_w_cpue_nad");
  set_w_cpue_mad.allocate("set_w_cpue_mad");
  set_w_cpue_sad.allocate("set_w_cpue_sad");
  set_w_cpue_jai.allocate("set_w_cpue_jai");
  set_w_cpue_mareco.allocate("set_w_cpue_mareco");
  set_w_lenc_nad.allocate("set_w_lenc_nad");
  set_w_lenc_mad.allocate("set_w_lenc_mad");
  set_w_agec_cRn.allocate("set_w_agec_cRn");
  set_w_agec_cRs.allocate("set_w_agec_cRs");
  set_w_agec_cBn.allocate("set_w_agec_cBn");
  set_w_agec_cBs.allocate("set_w_agec_cBs");
  set_w_Nage_init.allocate("set_w_Nage_init");
  set_w_rec.allocate("set_w_rec");
  set_w_rec_early.allocate("set_w_rec_early");
  set_w_rec_end.allocate("set_w_rec_end");
  set_w_fullF.allocate("set_w_fullF");
  set_w_Ftune.allocate("set_w_Ftune");
  wgtpar_a.allocate("wgtpar_a");
  wgtpar_b.allocate("wgtpar_b");
  obs_maturity_f.allocate(1,nages,"obs_maturity_f");
  obs_maturity_m.allocate(1,nages,"obs_maturity_m");
  maturity_obs_tv.allocate(styr,endyr,1,nages,"maturity_obs_tv");
  obs_prop_f.allocate(1,nages,"obs_prop_f");
  fecundity.allocate(1,nages,"fecundity");
  fecundity_tv.allocate(styr,endyr,1,nages,"fecundity_tv");
  wgt_spawn.allocate(1,nages,"wgt_spawn");
  wgt_spawn_tv.allocate(styr,endyr,1,nages,"wgt_spawn_tv");
  wgt_middle.allocate(1,nages,"wgt_middle");
  wgt_middle_tv.allocate(styr,endyr,1,nages,"wgt_middle_tv");
  len_apr15_tv.allocate(styr,endyr,1,nages,"len_apr15_tv");
  len_jun1_tv.allocate(styr,endyr,1,nages,"len_jun1_tv");
  len_oct15_tv.allocate(styr,endyr,1,nages,"len_oct15_tv");
  spawn_time_frac.allocate("spawn_time_frac");
  set_M.allocate(1,nages,"set_M");
  set_M_tv.allocate(styr,endyr,1,nages,"set_M_tv");
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
  minSS_lenc_nad.allocate("minSS_lenc_nad");
  minSS_lenc_mad.allocate("minSS_lenc_mad");
  minSS_agec_cRn.allocate("minSS_agec_cRn");
  minSS_agec_cRs.allocate("minSS_agec_cRs");
  minSS_agec_cBn.allocate("minSS_agec_cBn");
  minSS_agec_cBs.allocate("minSS_agec_cBs");
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
  const double len_cv_nad_LO=set_len_cv_nad(2); const double len_cv_nad_HI=set_len_cv_nad(3); const double len_cv_nad_PH=set_len_cv_nad(4); 
  const double len_cv_mad_LO=set_len_cv_mad(2); const double len_cv_mad_HI=set_len_cv_mad(3); const double len_cv_mad_PH=set_len_cv_mad(4); 
  //const double len_cv_sad_LO=set_len_cv_sad(2); const double len_cv_sad_HI=set_len_cv_sad(3); const double len_cv_sad_PH=set_len_cv_sad(4); 
  const double steep_LO=set_steep(2); const double steep_HI=set_steep(3); const double steep_PH=set_steep(4);
  const double log_R0_LO=set_log_R0(2); const double log_R0_HI=set_log_R0(3); const double log_R0_PH=set_log_R0(4);
  const double R_autocorr_LO=set_R_autocorr(2); const double R_autocorr_HI=set_R_autocorr(3); const double R_autocorr_PH=set_R_autocorr(4);
  const double rec_sigma_LO=set_rec_sigma(2); const double rec_sigma_HI=set_rec_sigma(3); const double rec_sigma_PH=set_rec_sigma(4);
  
  const double log_dm_cRn_ac_LO=set_log_dm_agec_cRn(2); const double log_dm_cRn_ac_HI=set_log_dm_agec_cRn(3); const double log_dm_cRn_ac_PH=set_log_dm_agec_cRn(4);
  const double log_dm_cRs_ac_LO=set_log_dm_agec_cRs(2); const double log_dm_cRs_ac_HI=set_log_dm_agec_cRs(3); const double log_dm_cRs_ac_PH=set_log_dm_agec_cRs(4);
  const double log_dm_cBn_ac_LO=set_log_dm_agec_cBn(2); const double log_dm_cBn_ac_HI=set_log_dm_agec_cBn(3); const double log_dm_cBn_ac_PH=set_log_dm_agec_cBn(4);
  const double log_dm_cBs_ac_LO=set_log_dm_agec_cBs(2); const double log_dm_cBs_ac_HI=set_log_dm_agec_cBs(3); const double log_dm_cBs_ac_PH=set_log_dm_agec_cBs(4);
  const double log_dm_nad_lc_LO=set_log_dm_lenc_nad(2); const double log_dm_nad_lc_HI=set_log_dm_lenc_nad(3); const double log_dm_nad_lc_PH=set_log_dm_lenc_nad(4);
  const double log_dm_mad_lc_LO=set_log_dm_lenc_mad(2); const double log_dm_mad_lc_HI=set_log_dm_lenc_mad(3); const double log_dm_mad_lc_PH=set_log_dm_lenc_mad(4);
  //const double log_dm_sad_lc_LO=set_log_dm_sad_lc(2); const double log_dm_sad_lc_HI=set_log_dm_sad_lc(3); const double log_dm_sad_lc_PH=set_log_dm_sad_lc(4);
  const double selpar_A50_cRn_LO=set_selpar_A50_cRn(2); const double selpar_A50_cRn_HI=set_selpar_A50_cRn(3); const double selpar_A50_cRn_PH=set_selpar_A50_cRn(4);
  const double selpar_slope_cRn_LO=set_selpar_slope_cRn(2); const double selpar_slope_cRn_HI=set_selpar_slope_cRn(3); const double selpar_slope_cRn_PH=set_selpar_slope_cRn(4);
  const double selpar_A502_cRn_LO=set_selpar_A502_cRn(2); const double selpar_A502_cRn_HI=set_selpar_A502_cRn(3); const double selpar_A502_cRn_PH=set_selpar_A502_cRn(4);
  const double selpar_slope2_cRn_LO=set_selpar_slope2_cRn(2); const double selpar_slope2_cRn_HI=set_selpar_slope2_cRn(3); const double selpar_slope2_cRn_PH=set_selpar_slope2_cRn(4);
  const double selpar_A50_cRn2_LO=set_selpar_A50_cRn2(2); const double selpar_A50_cRn2_HI=set_selpar_A50_cRn2(3); const double selpar_A50_cRn2_PH=set_selpar_A50_cRn2(4);
  const double selpar_slope_cRn2_LO=set_selpar_slope_cRn2(2); const double selpar_slope_cRn2_HI=set_selpar_slope_cRn2(3); const double selpar_slope_cRn2_PH=set_selpar_slope_cRn2(4);
  const double selpar_A502_cRn2_LO=set_selpar_A502_cRn2(2); const double selpar_A502_cRn2_HI=set_selpar_A502_cRn2(3); const double selpar_A502_cRn2_PH=set_selpar_A502_cRn2(4);
  const double selpar_slope2_cRn2_LO=set_selpar_slope2_cRn2(2); const double selpar_slope2_cRn2_HI=set_selpar_slope2_cRn2(3); const double selpar_slope2_cRn2_PH=set_selpar_slope2_cRn2(4);
  const double selpar_A50_cRn3_LO=set_selpar_A50_cRn3(2); const double selpar_A50_cRn3_HI=set_selpar_A50_cRn3(3); const double selpar_A50_cRn3_PH=set_selpar_A50_cRn3(4);
  const double selpar_slope_cRn3_LO=set_selpar_slope_cRn3(2); const double selpar_slope_cRn3_HI=set_selpar_slope_cRn3(3); const double selpar_slope_cRn3_PH=set_selpar_slope_cRn3(4);
  const double selpar_A502_cRn3_LO=set_selpar_A502_cRn3(2); const double selpar_A502_cRn3_HI=set_selpar_A502_cRn3(3); const double selpar_A502_cRn3_PH=set_selpar_A502_cRn3(4);
  const double selpar_slope2_cRn3_LO=set_selpar_slope2_cRn3(2); const double selpar_slope2_cRn3_HI=set_selpar_slope2_cRn3(3); const double selpar_slope2_cRn3_PH=set_selpar_slope2_cRn3(4);
  const double selpar_A50_cRn4_LO=set_selpar_A50_cRn4(2); const double selpar_A50_cRn4_HI=set_selpar_A50_cRn4(3); const double selpar_A50_cRn4_PH=set_selpar_A50_cRn4(4);
  const double selpar_slope_cRn4_LO=set_selpar_slope_cRn4(2); const double selpar_slope_cRn4_HI=set_selpar_slope_cRn4(3); const double selpar_slope_cRn4_PH=set_selpar_slope_cRn4(4);
  const double selpar_A502_cRn4_LO=set_selpar_A502_cRn4(2); const double selpar_A502_cRn4_HI=set_selpar_A502_cRn4(3); const double selpar_A502_cRn4_PH=set_selpar_A502_cRn4(4);
  const double selpar_slope2_cRn4_LO=set_selpar_slope2_cRn4(2); const double selpar_slope2_cRn4_HI=set_selpar_slope2_cRn4(3); const double selpar_slope2_cRn4_PH=set_selpar_slope2_cRn4(4);
  const double selpar_A50_cRs_LO=set_selpar_A50_cRs(2); const double selpar_A50_cRs_HI=set_selpar_A50_cRs(3); const double selpar_A50_cRs_PH=set_selpar_A50_cRs(4);
  const double selpar_slope_cRs_LO=set_selpar_slope_cRs(2); const double selpar_slope_cRs_HI=set_selpar_slope_cRs(3); const double selpar_slope_cRs_PH=set_selpar_slope_cRs(4);
  const double selpar_A502_cRs_LO=set_selpar_A502_cRs(2); const double selpar_A502_cRs_HI=set_selpar_A502_cRs(3); const double selpar_A502_cRs_PH=set_selpar_A502_cRs(4);
  const double selpar_slope2_cRs_LO=set_selpar_slope2_cRs(2); const double selpar_slope2_cRs_HI=set_selpar_slope2_cRs(3); const double selpar_slope2_cRs_PH=set_selpar_slope2_cRs(4);
  const double selpar_A50_cRs2_LO=set_selpar_A50_cRs2(2); const double selpar_A50_cRs2_HI=set_selpar_A50_cRs2(3); const double selpar_A50_cRs2_PH=set_selpar_A50_cRs2(4);
  const double selpar_slope_cRs2_LO=set_selpar_slope_cRs2(2); const double selpar_slope_cRs2_HI=set_selpar_slope_cRs2(3); const double selpar_slope_cRs2_PH=set_selpar_slope_cRs2(4);
  const double selpar_A502_cRs2_LO=set_selpar_A502_cRs2(2); const double selpar_A502_cRs2_HI=set_selpar_A502_cRs2(3); const double selpar_A502_cRs2_PH=set_selpar_A502_cRs2(4);
  const double selpar_slope2_cRs2_LO=set_selpar_slope2_cRs2(2); const double selpar_slope2_cRs2_HI=set_selpar_slope2_cRs2(3); const double selpar_slope2_cRs2_PH=set_selpar_slope2_cRs2(4);
  const double selpar_A50_cRs3_LO=set_selpar_A50_cRs3(2); const double selpar_A50_cRs3_HI=set_selpar_A50_cRs3(3); const double selpar_A50_cRs3_PH=set_selpar_A50_cRs3(4);
  const double selpar_slope_cRs3_LO=set_selpar_slope_cRs3(2); const double selpar_slope_cRs3_HI=set_selpar_slope_cRs3(3); const double selpar_slope_cRs3_PH=set_selpar_slope_cRs3(4);
  const double selpar_A502_cRs3_LO=set_selpar_A502_cRs3(2); const double selpar_A502_cRs3_HI=set_selpar_A502_cRs3(3); const double selpar_A502_cRs3_PH=set_selpar_A502_cRs3(4);
  const double selpar_slope2_cRs3_LO=set_selpar_slope2_cRs3(2); const double selpar_slope2_cRs3_HI=set_selpar_slope2_cRs3(3); const double selpar_slope2_cRs3_PH=set_selpar_slope2_cRs3(4);
  const double selpar_A50_cRs4_LO=set_selpar_A50_cRs4(2); const double selpar_A50_cRs4_HI=set_selpar_A50_cRs4(3); const double selpar_A50_cRs4_PH=set_selpar_A50_cRs4(4);
  const double selpar_slope_cRs4_LO=set_selpar_slope_cRs4(2); const double selpar_slope_cRs4_HI=set_selpar_slope_cRs4(3); const double selpar_slope_cRs4_PH=set_selpar_slope_cRs4(4);
  const double selpar_A502_cRs4_LO=set_selpar_A502_cRs4(2); const double selpar_A502_cRs4_HI=set_selpar_A502_cRs4(3); const double selpar_A502_cRs4_PH=set_selpar_A502_cRs4(4);
  const double selpar_slope2_cRs4_LO=set_selpar_slope2_cRs4(2); const double selpar_slope2_cRs4_HI=set_selpar_slope2_cRs4(3); const double selpar_slope2_cRs4_PH=set_selpar_slope2_cRs4(4);
  const double selpar_A50_cBn_LO=set_selpar_A50_cBn(2); const double selpar_A50_cBn_HI=set_selpar_A50_cBn(3); const double selpar_A50_cBn_PH=set_selpar_A50_cBn(4);
  const double selpar_slope_cBn_LO=set_selpar_slope_cBn(2); const double selpar_slope_cBn_HI=set_selpar_slope_cBn(3); const double selpar_slope_cBn_PH=set_selpar_slope_cBn(4);
  const double selpar_A502_cBn_LO=set_selpar_A502_cBn(2); const double selpar_A502_cBn_HI=set_selpar_A502_cBn(3); const double selpar_A502_cBn_PH=set_selpar_A502_cBn(4);
  const double selpar_slope2_cBn_LO=set_selpar_slope2_cBn(2); const double selpar_slope2_cBn_HI=set_selpar_slope2_cBn(3); const double selpar_slope2_cBn_PH=set_selpar_slope2_cBn(4);
  const double selpar_A50_cBn2_LO=set_selpar_A50_cBn2(2); const double selpar_A50_cBn2_HI=set_selpar_A50_cBn2(3); const double selpar_A50_cBn2_PH=set_selpar_A50_cBn2(4);
  const double selpar_slope_cBn2_LO=set_selpar_slope_cBn2(2); const double selpar_slope_cBn2_HI=set_selpar_slope_cBn2(3); const double selpar_slope_cBn2_PH=set_selpar_slope_cBn2(4);
  const double selpar_A502_cBn2_LO=set_selpar_A502_cBn2(2); const double selpar_A502_cBn2_HI=set_selpar_A502_cBn2(3); const double selpar_A502_cBn2_PH=set_selpar_A502_cBn2(4);
  const double selpar_slope2_cBn2_LO=set_selpar_slope2_cBn2(2); const double selpar_slope2_cBn2_HI=set_selpar_slope2_cBn2(3); const double selpar_slope2_cBn2_PH=set_selpar_slope2_cBn2(4);
  const double selpar_A50_cBs_LO=set_selpar_A50_cBs(2); const double selpar_A50_cBs_HI=set_selpar_A50_cBs(3); const double selpar_A50_cBs_PH=set_selpar_A50_cBs(4);
  const double selpar_slope_cBs_LO=set_selpar_slope_cBs(2); const double selpar_slope_cBs_HI=set_selpar_slope_cBs(3); const double selpar_slope_cBs_PH=set_selpar_slope_cBs(4);
  const double selpar_A502_cBs_LO=set_selpar_A502_cBs(2); const double selpar_A502_cBs_HI=set_selpar_A502_cBs(3); const double selpar_A502_cBs_PH=set_selpar_A502_cBs(4);
  const double selpar_slope2_cBs_LO=set_selpar_slope2_cBs(2); const double selpar_slope2_cBs_HI=set_selpar_slope2_cBs(3); const double selpar_slope2_cBs_PH=set_selpar_slope2_cBs(4);
  const double selpar_A50_cBs2_LO=set_selpar_A50_cBs2(2); const double selpar_A50_cBs2_HI=set_selpar_A50_cBs2(3); const double selpar_A50_cBs2_PH=set_selpar_A50_cBs2(4);
  const double selpar_slope_cBs2_LO=set_selpar_slope_cBs2(2); const double selpar_slope_cBs2_HI=set_selpar_slope_cBs2(3); const double selpar_slope_cBs2_PH=set_selpar_slope_cBs2(4);
  const double selpar_A502_cBs2_LO=set_selpar_A502_cBs2(2); const double selpar_A502_cBs2_HI=set_selpar_A502_cBs2(3); const double selpar_A502_cBs2_PH=set_selpar_A502_cBs2(4);
  const double selpar_slope2_cBs2_LO=set_selpar_slope2_cBs2(2); const double selpar_slope2_cBs2_HI=set_selpar_slope2_cBs2(3); const double selpar_slope2_cBs2_PH=set_selpar_slope2_cBs2(4);
  const double selpar_age0_cRn_LO=set_sel_age0_cRn(2); const double selpar_age0_cRn_HI=set_sel_age0_cRn(3); const double selpar_age0_cRn_PH=set_sel_age0_cRn(4);
  const double selpar_age1_cRn_LO=set_sel_age1_cRn(2); const double selpar_age1_cRn_HI=set_sel_age1_cRn(3); const double selpar_age1_cRn_PH=set_sel_age1_cRn(4);
  const double selpar_age2_cRn_LO=set_sel_age2_cRn(2); const double selpar_age2_cRn_HI=set_sel_age2_cRn(3); const double selpar_age2_cRn_PH=set_sel_age2_cRn(4);
  const double selpar_age3_cRn_LO=set_sel_age3_cRn(2); const double selpar_age3_cRn_HI=set_sel_age3_cRn(3); const double selpar_age3_cRn_PH=set_sel_age3_cRn(4);
  const double selpar_age4_cRn_LO=set_sel_age4_cRn(2); const double selpar_age4_cRn_HI=set_sel_age4_cRn(3); const double selpar_age4_cRn_PH=set_sel_age4_cRn(4);
  const double selpar_age5_cRn_LO=set_sel_age5_cRn(2); const double selpar_age5_cRn_HI=set_sel_age5_cRn(3); const double selpar_age5_cRn_PH=set_sel_age5_cRn(4);
  const double selpar_age6_cRn_LO=set_sel_age6_cRn(2); const double selpar_age6_cRn_HI=set_sel_age6_cRn(3); const double selpar_age6_cRn_PH=set_sel_age6_cRn(4);
  const double selpar_age0_cRn2_LO=set_sel_age0_cRn2(2); const double selpar_age0_cRn2_HI=set_sel_age0_cRn2(3); const double selpar_age0_cRn2_PH=set_sel_age0_cRn2(4);
  const double selpar_age1_cRn2_LO=set_sel_age1_cRn2(2); const double selpar_age1_cRn2_HI=set_sel_age1_cRn2(3); const double selpar_age1_cRn2_PH=set_sel_age1_cRn2(4);
  const double selpar_age2_cRn2_LO=set_sel_age2_cRn2(2); const double selpar_age2_cRn2_HI=set_sel_age2_cRn2(3); const double selpar_age2_cRn2_PH=set_sel_age2_cRn2(4);
  const double selpar_age3_cRn2_LO=set_sel_age3_cRn2(2); const double selpar_age3_cRn2_HI=set_sel_age3_cRn2(3); const double selpar_age3_cRn2_PH=set_sel_age3_cRn2(4);
  const double selpar_age4_cRn2_LO=set_sel_age4_cRn2(2); const double selpar_age4_cRn2_HI=set_sel_age4_cRn2(3); const double selpar_age4_cRn2_PH=set_sel_age4_cRn2(4);
  const double selpar_age5_cRn2_LO=set_sel_age5_cRn2(2); const double selpar_age5_cRn2_HI=set_sel_age5_cRn2(3); const double selpar_age5_cRn2_PH=set_sel_age5_cRn2(4);
  const double selpar_age6_cRn2_LO=set_sel_age6_cRn2(2); const double selpar_age6_cRn2_HI=set_sel_age6_cRn2(3); const double selpar_age6_cRn2_PH=set_sel_age6_cRn2(4);
  const double selpar_age0_cRn3_LO=set_sel_age0_cRn3(2); const double selpar_age0_cRn3_HI=set_sel_age0_cRn3(3); const double selpar_age0_cRn3_PH=set_sel_age0_cRn3(4);
  const double selpar_age1_cRn3_LO=set_sel_age1_cRn3(2); const double selpar_age1_cRn3_HI=set_sel_age1_cRn3(3); const double selpar_age1_cRn3_PH=set_sel_age1_cRn3(4);
  const double selpar_age2_cRn3_LO=set_sel_age2_cRn3(2); const double selpar_age2_cRn3_HI=set_sel_age2_cRn3(3); const double selpar_age2_cRn3_PH=set_sel_age2_cRn3(4);
  const double selpar_age3_cRn3_LO=set_sel_age3_cRn3(2); const double selpar_age3_cRn3_HI=set_sel_age3_cRn3(3); const double selpar_age3_cRn3_PH=set_sel_age3_cRn3(4);
  const double selpar_age4_cRn3_LO=set_sel_age4_cRn3(2); const double selpar_age4_cRn3_HI=set_sel_age4_cRn3(3); const double selpar_age4_cRn3_PH=set_sel_age4_cRn3(4);
  const double selpar_age5_cRn3_LO=set_sel_age5_cRn3(2); const double selpar_age5_cRn3_HI=set_sel_age5_cRn3(3); const double selpar_age5_cRn3_PH=set_sel_age5_cRn3(4);
  const double selpar_age6_cRn3_LO=set_sel_age6_cRn3(2); const double selpar_age6_cRn3_HI=set_sel_age6_cRn3(3); const double selpar_age6_cRn3_PH=set_sel_age6_cRn3(4);
  const double selpar_age0_cRn4_LO=set_sel_age0_cRn4(2); const double selpar_age0_cRn4_HI=set_sel_age0_cRn4(3); const double selpar_age0_cRn4_PH=set_sel_age0_cRn4(4);
  const double selpar_age1_cRn4_LO=set_sel_age1_cRn4(2); const double selpar_age1_cRn4_HI=set_sel_age1_cRn4(3); const double selpar_age1_cRn4_PH=set_sel_age1_cRn4(4);
  const double selpar_age2_cRn4_LO=set_sel_age2_cRn4(2); const double selpar_age2_cRn4_HI=set_sel_age2_cRn4(3); const double selpar_age2_cRn4_PH=set_sel_age2_cRn4(4);
  const double selpar_age3_cRn4_LO=set_sel_age3_cRn4(2); const double selpar_age3_cRn4_HI=set_sel_age3_cRn4(3); const double selpar_age3_cRn4_PH=set_sel_age3_cRn4(4);
  const double selpar_age4_cRn4_LO=set_sel_age4_cRn4(2); const double selpar_age4_cRn4_HI=set_sel_age4_cRn4(3); const double selpar_age4_cRn4_PH=set_sel_age4_cRn4(4);
  const double selpar_age5_cRn4_LO=set_sel_age5_cRn4(2); const double selpar_age5_cRn4_HI=set_sel_age5_cRn4(3); const double selpar_age5_cRn4_PH=set_sel_age5_cRn4(4);
  const double selpar_age6_cRn4_LO=set_sel_age6_cRn4(2); const double selpar_age6_cRn4_HI=set_sel_age6_cRn4(3); const double selpar_age6_cRn4_PH=set_sel_age6_cRn4(4);
  const double selpar_age0_cRs_LO=set_sel_age0_cRs(2); const double selpar_age0_cRs_HI=set_sel_age0_cRs(3); const double selpar_age0_cRs_PH=set_sel_age0_cRs(4);
  const double selpar_age1_cRs_LO=set_sel_age1_cRs(2); const double selpar_age1_cRs_HI=set_sel_age1_cRs(3); const double selpar_age1_cRs_PH=set_sel_age1_cRs(4);
  const double selpar_age2_cRs_LO=set_sel_age2_cRs(2); const double selpar_age2_cRs_HI=set_sel_age2_cRs(3); const double selpar_age2_cRs_PH=set_sel_age2_cRs(4);
  const double selpar_age3_cRs_LO=set_sel_age3_cRs(2); const double selpar_age3_cRs_HI=set_sel_age3_cRs(3); const double selpar_age3_cRs_PH=set_sel_age3_cRs(4);
  const double selpar_age4_cRs_LO=set_sel_age4_cRs(2); const double selpar_age4_cRs_HI=set_sel_age4_cRs(3); const double selpar_age4_cRs_PH=set_sel_age4_cRs(4);
  const double selpar_age5_cRs_LO=set_sel_age5_cRs(2); const double selpar_age5_cRs_HI=set_sel_age5_cRs(3); const double selpar_age5_cRs_PH=set_sel_age5_cRs(4);
  const double selpar_age6_cRs_LO=set_sel_age6_cRs(2); const double selpar_age6_cRs_HI=set_sel_age6_cRs(3); const double selpar_age6_cRs_PH=set_sel_age6_cRs(4);
  const double selpar_age0_cRs2_LO=set_sel_age0_cRs2(2); const double selpar_age0_cRs2_HI=set_sel_age0_cRs2(3); const double selpar_age0_cRs2_PH=set_sel_age0_cRs2(4);
  const double selpar_age1_cRs2_LO=set_sel_age1_cRs2(2); const double selpar_age1_cRs2_HI=set_sel_age1_cRs2(3); const double selpar_age1_cRs2_PH=set_sel_age1_cRs2(4);
  const double selpar_age2_cRs2_LO=set_sel_age2_cRs2(2); const double selpar_age2_cRs2_HI=set_sel_age2_cRs2(3); const double selpar_age2_cRs2_PH=set_sel_age2_cRs2(4);
  const double selpar_age3_cRs2_LO=set_sel_age3_cRs2(2); const double selpar_age3_cRs2_HI=set_sel_age3_cRs2(3); const double selpar_age3_cRs2_PH=set_sel_age3_cRs2(4);
  const double selpar_age4_cRs2_LO=set_sel_age4_cRs2(2); const double selpar_age4_cRs2_HI=set_sel_age4_cRs2(3); const double selpar_age4_cRs2_PH=set_sel_age4_cRs2(4);
  const double selpar_age5_cRs2_LO=set_sel_age5_cRs2(2); const double selpar_age5_cRs2_HI=set_sel_age5_cRs2(3); const double selpar_age5_cRs2_PH=set_sel_age5_cRs2(4);
  const double selpar_age6_cRs2_LO=set_sel_age6_cRs2(2); const double selpar_age6_cRs2_HI=set_sel_age6_cRs2(3); const double selpar_age6_cRs2_PH=set_sel_age6_cRs2(4);
  const double selpar_age0_cRs3_LO=set_sel_age0_cRs3(2); const double selpar_age0_cRs3_HI=set_sel_age0_cRs3(3); const double selpar_age0_cRs3_PH=set_sel_age0_cRs3(4);
  const double selpar_age1_cRs3_LO=set_sel_age1_cRs3(2); const double selpar_age1_cRs3_HI=set_sel_age1_cRs3(3); const double selpar_age1_cRs3_PH=set_sel_age1_cRs3(4);
  const double selpar_age2_cRs3_LO=set_sel_age2_cRs3(2); const double selpar_age2_cRs3_HI=set_sel_age2_cRs3(3); const double selpar_age2_cRs3_PH=set_sel_age2_cRs3(4);
  const double selpar_age3_cRs3_LO=set_sel_age3_cRs3(2); const double selpar_age3_cRs3_HI=set_sel_age3_cRs3(3); const double selpar_age3_cRs3_PH=set_sel_age3_cRs3(4);
  const double selpar_age4_cRs3_LO=set_sel_age4_cRs3(2); const double selpar_age4_cRs3_HI=set_sel_age4_cRs3(3); const double selpar_age4_cRs3_PH=set_sel_age4_cRs3(4);
  const double selpar_age5_cRs3_LO=set_sel_age5_cRs3(2); const double selpar_age5_cRs3_HI=set_sel_age5_cRs3(3); const double selpar_age5_cRs3_PH=set_sel_age5_cRs3(4);
  const double selpar_age6_cRs3_LO=set_sel_age6_cRs3(2); const double selpar_age6_cRs3_HI=set_sel_age6_cRs3(3); const double selpar_age6_cRs3_PH=set_sel_age6_cRs3(4);
  const double selpar_age0_cRs4_LO=set_sel_age0_cRs4(2); const double selpar_age0_cRs4_HI=set_sel_age0_cRs4(3); const double selpar_age0_cRs4_PH=set_sel_age0_cRs4(4);
  const double selpar_age1_cRs4_LO=set_sel_age1_cRs4(2); const double selpar_age1_cRs4_HI=set_sel_age1_cRs4(3); const double selpar_age1_cRs4_PH=set_sel_age1_cRs4(4);
  const double selpar_age2_cRs4_LO=set_sel_age2_cRs4(2); const double selpar_age2_cRs4_HI=set_sel_age2_cRs4(3); const double selpar_age2_cRs4_PH=set_sel_age2_cRs4(4);
  const double selpar_age3_cRs4_LO=set_sel_age3_cRs4(2); const double selpar_age3_cRs4_HI=set_sel_age3_cRs4(3); const double selpar_age3_cRs4_PH=set_sel_age3_cRs4(4);
  const double selpar_age4_cRs4_LO=set_sel_age4_cRs4(2); const double selpar_age4_cRs4_HI=set_sel_age4_cRs4(3); const double selpar_age4_cRs4_PH=set_sel_age4_cRs4(4);
  const double selpar_age5_cRs4_LO=set_sel_age5_cRs4(2); const double selpar_age5_cRs4_HI=set_sel_age5_cRs4(3); const double selpar_age5_cRs4_PH=set_sel_age5_cRs4(4);
  const double selpar_age6_cRs4_LO=set_sel_age6_cRs4(2); const double selpar_age6_cRs4_HI=set_sel_age6_cRs4(3); const double selpar_age6_cRs4_PH=set_sel_age6_cRs4(4);
  const double selpar_age0_cBn_LO=set_sel_age0_cBn(2); const double selpar_age0_cBn_HI=set_sel_age0_cBn(3); const double selpar_age0_cBn_PH=set_sel_age0_cBn(4);
  const double selpar_age1_cBn_LO=set_sel_age1_cBn(2); const double selpar_age1_cBn_HI=set_sel_age1_cBn(3); const double selpar_age1_cBn_PH=set_sel_age1_cBn(4);
  const double selpar_age2_cBn_LO=set_sel_age2_cBn(2); const double selpar_age2_cBn_HI=set_sel_age2_cBn(3); const double selpar_age2_cBn_PH=set_sel_age2_cBn(4);
  const double selpar_age3_cBn_LO=set_sel_age3_cBn(2); const double selpar_age3_cBn_HI=set_sel_age3_cBn(3); const double selpar_age3_cBn_PH=set_sel_age3_cBn(4);
  const double selpar_age4_cBn_LO=set_sel_age4_cBn(2); const double selpar_age4_cBn_HI=set_sel_age4_cBn(3); const double selpar_age4_cBn_PH=set_sel_age4_cBn(4);
  const double selpar_age5_cBn_LO=set_sel_age5_cBn(2); const double selpar_age5_cBn_HI=set_sel_age5_cBn(3); const double selpar_age5_cBn_PH=set_sel_age5_cBn(4);
  const double selpar_age6_cBn_LO=set_sel_age6_cBn(2); const double selpar_age6_cBn_HI=set_sel_age6_cBn(3); const double selpar_age6_cBn_PH=set_sel_age6_cBn(4);
  const double selpar_age0_cBn2_LO=set_sel_age0_cBn2(2); const double selpar_age0_cBn2_HI=set_sel_age0_cBn2(3); const double selpar_age0_cBn2_PH=set_sel_age0_cBn2(4);
  const double selpar_age1_cBn2_LO=set_sel_age1_cBn2(2); const double selpar_age1_cBn2_HI=set_sel_age1_cBn2(3); const double selpar_age1_cBn2_PH=set_sel_age1_cBn2(4);
  const double selpar_age2_cBn2_LO=set_sel_age2_cBn2(2); const double selpar_age2_cBn2_HI=set_sel_age2_cBn2(3); const double selpar_age2_cBn2_PH=set_sel_age2_cBn2(4);
  const double selpar_age3_cBn2_LO=set_sel_age3_cBn2(2); const double selpar_age3_cBn2_HI=set_sel_age3_cBn2(3); const double selpar_age3_cBn2_PH=set_sel_age3_cBn2(4);
  const double selpar_age4_cBn2_LO=set_sel_age4_cBn2(2); const double selpar_age4_cBn2_HI=set_sel_age4_cBn2(3); const double selpar_age4_cBn2_PH=set_sel_age4_cBn2(4);
  const double selpar_age5_cBn2_LO=set_sel_age5_cBn2(2); const double selpar_age5_cBn2_HI=set_sel_age5_cBn2(3); const double selpar_age5_cBn2_PH=set_sel_age5_cBn2(4);
  const double selpar_age6_cBn2_LO=set_sel_age6_cBn2(2); const double selpar_age6_cBn2_HI=set_sel_age6_cBn2(3); const double selpar_age6_cBn2_PH=set_sel_age6_cBn2(4);
  const double selpar_age0_cBs_LO=set_sel_age0_cBs(2); const double selpar_age0_cBs_HI=set_sel_age0_cBs(3); const double selpar_age0_cBs_PH=set_sel_age0_cBs(4);
  const double selpar_age1_cBs_LO=set_sel_age1_cBs(2); const double selpar_age1_cBs_HI=set_sel_age1_cBs(3); const double selpar_age1_cBs_PH=set_sel_age1_cBs(4);
  const double selpar_age2_cBs_LO=set_sel_age2_cBs(2); const double selpar_age2_cBs_HI=set_sel_age2_cBs(3); const double selpar_age2_cBs_PH=set_sel_age2_cBs(4);
  const double selpar_age3_cBs_LO=set_sel_age3_cBs(2); const double selpar_age3_cBs_HI=set_sel_age3_cBs(3); const double selpar_age3_cBs_PH=set_sel_age3_cBs(4);
  const double selpar_age4_cBs_LO=set_sel_age4_cBs(2); const double selpar_age4_cBs_HI=set_sel_age4_cBs(3); const double selpar_age4_cBs_PH=set_sel_age4_cBs(4);
  const double selpar_age5_cBs_LO=set_sel_age5_cBs(2); const double selpar_age5_cBs_HI=set_sel_age5_cBs(3); const double selpar_age5_cBs_PH=set_sel_age5_cBs(4);
  const double selpar_age6_cBs_LO=set_sel_age6_cBs(2); const double selpar_age6_cBs_HI=set_sel_age6_cBs(3); const double selpar_age6_cBs_PH=set_sel_age6_cBs(4);
  const double selpar_age0_cBs2_LO=set_sel_age0_cBs2(2); const double selpar_age0_cBs2_HI=set_sel_age0_cBs2(3); const double selpar_age0_cBs2_PH=set_sel_age0_cBs2(4);
  const double selpar_age1_cBs2_LO=set_sel_age1_cBs2(2); const double selpar_age1_cBs2_HI=set_sel_age1_cBs2(3); const double selpar_age1_cBs2_PH=set_sel_age1_cBs2(4);
  const double selpar_age2_cBs2_LO=set_sel_age2_cBs2(2); const double selpar_age2_cBs2_HI=set_sel_age2_cBs2(3); const double selpar_age2_cBs2_PH=set_sel_age2_cBs2(4);
  const double selpar_age3_cBs2_LO=set_sel_age3_cBs2(2); const double selpar_age3_cBs2_HI=set_sel_age3_cBs2(3); const double selpar_age3_cBs2_PH=set_sel_age3_cBs2(4);
  const double selpar_age4_cBs2_LO=set_sel_age4_cBs2(2); const double selpar_age4_cBs2_HI=set_sel_age4_cBs2(3); const double selpar_age4_cBs2_PH=set_sel_age4_cBs2(4);
  const double selpar_age5_cBs2_LO=set_sel_age5_cBs2(2); const double selpar_age5_cBs2_HI=set_sel_age5_cBs2(3); const double selpar_age5_cBs2_PH=set_sel_age5_cBs2(4);
  const double selpar_age6_cBs2_LO=set_sel_age6_cBs2(2); const double selpar_age6_cBs2_HI=set_sel_age6_cBs2(3); const double selpar_age6_cBs2_PH=set_sel_age6_cBs2(4);
  const double selpar_A50_nad_LO=set_selpar_A50_nad(2); const double selpar_A50_nad_HI=set_selpar_A50_nad(3); const double selpar_A50_nad_PH=set_selpar_A50_nad(4);
  const double selpar_slope_nad_LO=set_selpar_slope_nad(2); const double selpar_slope_nad_HI=set_selpar_slope_nad(3); const double selpar_slope_nad_PH=set_selpar_slope_nad(4);
  const double selpar_A502_nad_LO=set_selpar_A502_nad(2); const double selpar_A502_nad_HI=set_selpar_A502_nad(3); const double selpar_A502_nad_PH=set_selpar_A502_nad(4);
  const double selpar_slope2_nad_LO=set_selpar_slope2_nad(2); const double selpar_slope2_nad_HI=set_selpar_slope2_nad(3); const double selpar_slope2_nad_PH=set_selpar_slope2_nad(4);
  const double selpar_A50_mad_LO=set_selpar_A50_mad(2); const double selpar_A50_mad_HI=set_selpar_A50_mad(3); const double selpar_A50_mad_PH=set_selpar_A50_mad(4);
  const double selpar_slope_mad_LO=set_selpar_slope_mad(2); const double selpar_slope_mad_HI=set_selpar_slope_mad(3); const double selpar_slope_mad_PH=set_selpar_slope_mad(4);
  const double selpar_A502_mad_LO=set_selpar_A502_mad(2); const double selpar_A502_mad_HI=set_selpar_A502_mad(3); const double selpar_A502_mad_PH=set_selpar_A502_mad(4);
  const double selpar_slope2_mad_LO=set_selpar_slope2_mad(2); const double selpar_slope2_mad_HI=set_selpar_slope2_mad(3); const double selpar_slope2_mad_PH=set_selpar_slope2_mad(4);
  const double selpar_A50_sad_LO=set_selpar_A50_sad(2); const double selpar_A50_sad_HI=set_selpar_A50_sad(3); const double selpar_A50_sad_PH=set_selpar_A50_sad(4);
  const double selpar_slope_sad_LO=set_selpar_slope_sad(2); const double selpar_slope_sad_HI=set_selpar_slope_sad(3); const double selpar_slope_sad_PH=set_selpar_slope_sad(4);
  const double selpar_A502_sad_LO=set_selpar_A502_sad(2); const double selpar_A502_sad_HI=set_selpar_A502_sad(3); const double selpar_A502_sad_PH=set_selpar_A502_sad(4);
  const double selpar_slope2_sad_LO=set_selpar_slope2_sad(2); const double selpar_slope2_sad_HI=set_selpar_slope2_sad(3); const double selpar_slope2_sad_PH=set_selpar_slope2_sad(4);
  const double selpar_age0_nad_LO=set_sel_age0_nad(2); const double selpar_age0_nad_HI=set_sel_age0_nad(3); const double selpar_age0_nad_PH=set_sel_age0_nad(4);
  const double selpar_age1_nad_LO=set_sel_age1_nad(2); const double selpar_age1_nad_HI=set_sel_age1_nad(3); const double selpar_age1_nad_PH=set_sel_age1_nad(4);
  const double selpar_age2_nad_LO=set_sel_age2_nad(2); const double selpar_age2_nad_HI=set_sel_age2_nad(3); const double selpar_age2_nad_PH=set_sel_age2_nad(4);
  const double selpar_age3_nad_LO=set_sel_age3_nad(2); const double selpar_age3_nad_HI=set_sel_age3_nad(3); const double selpar_age3_nad_PH=set_sel_age3_nad(4);
  const double selpar_age4_nad_LO=set_sel_age4_nad(2); const double selpar_age4_nad_HI=set_sel_age4_nad(3); const double selpar_age4_nad_PH=set_sel_age4_nad(4);
  const double selpar_age5_nad_LO=set_sel_age5_nad(2); const double selpar_age5_nad_HI=set_sel_age5_nad(3); const double selpar_age5_nad_PH=set_sel_age5_nad(4);
  const double selpar_age6_nad_LO=set_sel_age6_nad(2); const double selpar_age6_nad_HI=set_sel_age6_nad(3); const double selpar_age6_nad_PH=set_sel_age6_nad(4);
  const double selpar_age0_mad_LO=set_sel_age0_mad(2); const double selpar_age0_mad_HI=set_sel_age0_mad(3); const double selpar_age0_mad_PH=set_sel_age0_mad(4);
  const double selpar_age1_mad_LO=set_sel_age1_mad(2); const double selpar_age1_mad_HI=set_sel_age1_mad(3); const double selpar_age1_mad_PH=set_sel_age1_mad(4);
  const double selpar_age2_mad_LO=set_sel_age2_mad(2); const double selpar_age2_mad_HI=set_sel_age2_mad(3); const double selpar_age2_mad_PH=set_sel_age2_mad(4);
  const double selpar_age3_mad_LO=set_sel_age3_mad(2); const double selpar_age3_mad_HI=set_sel_age3_mad(3); const double selpar_age3_mad_PH=set_sel_age3_mad(4);
  const double selpar_age4_mad_LO=set_sel_age4_mad(2); const double selpar_age4_mad_HI=set_sel_age4_mad(3); const double selpar_age4_mad_PH=set_sel_age4_mad(4);
  const double selpar_age5_mad_LO=set_sel_age5_mad(2); const double selpar_age5_mad_HI=set_sel_age5_mad(3); const double selpar_age5_mad_PH=set_sel_age5_mad(4);
  const double selpar_age6_mad_LO=set_sel_age6_mad(2); const double selpar_age6_mad_HI=set_sel_age6_mad(3); const double selpar_age6_mad_PH=set_sel_age6_mad(4);
  const double selpar_age0_sad_LO=set_sel_age0_sad(2); const double selpar_age0_sad_HI=set_sel_age0_sad(3); const double selpar_age0_sad_PH=set_sel_age0_sad(4);
  const double selpar_age1_sad_LO=set_sel_age1_sad(2); const double selpar_age1_sad_HI=set_sel_age1_sad(3); const double selpar_age1_sad_PH=set_sel_age1_sad(4);
  const double selpar_age2_sad_LO=set_sel_age2_sad(2); const double selpar_age2_sad_HI=set_sel_age2_sad(3); const double selpar_age2_sad_PH=set_sel_age2_sad(4);
  const double selpar_age3_sad_LO=set_sel_age3_sad(2); const double selpar_age3_sad_HI=set_sel_age3_sad(3); const double selpar_age3_sad_PH=set_sel_age3_sad(4);
  const double selpar_age4_sad_LO=set_sel_age4_sad(2); const double selpar_age4_sad_HI=set_sel_age4_sad(3); const double selpar_age4_sad_PH=set_sel_age4_sad(4);
  const double selpar_age5_sad_LO=set_sel_age5_sad(2); const double selpar_age5_sad_HI=set_sel_age5_sad(3); const double selpar_age5_sad_PH=set_sel_age5_sad(4);
  const double selpar_age6_sad_LO=set_sel_age6_sad(2); const double selpar_age6_sad_HI=set_sel_age6_sad(3); const double selpar_age6_sad_PH=set_sel_age6_sad(4);
  const double log_q_nad_LO=set_log_q_cpue_nad(2); const double log_q_nad_HI=set_log_q_cpue_nad(3); const double log_q_nad_PH=set_log_q_cpue_nad(4);
  const double log_q_mad_LO=set_log_q_cpue_mad(2); const double log_q_mad_HI=set_log_q_cpue_mad(3); const double log_q_mad_PH=set_log_q_cpue_mad(4);
  const double log_q_sad_LO=set_log_q_cpue_sad(2); const double log_q_sad_HI=set_log_q_cpue_sad(3); const double log_q_sad_PH=set_log_q_cpue_sad(4);
  const double log_q_jai_LO=set_log_q_cpue_jai(2); const double log_q_jai_HI=set_log_q_cpue_jai(3); const double log_q_jai_PH=set_log_q_cpue_jai(4);
  const double log_q2_jai_LO=set_log_q2_jai(2); const double log_q2_jai_HI=set_log_q2_jai(3); const double log_q2_jai_PH=set_log_q2_jai(4);
  const double log_q_mar_LO=set_log_q_cpue_mar(2); const double log_q_mar_HI=set_log_q_cpue_mar(3); const double log_q_mar_PH=set_log_q_cpue_mar(4);
  const double log_q_eco_LO=set_log_q_cpue_eco(2); const double log_q_eco_HI=set_log_q_cpue_eco(3); const double log_q_eco_PH=set_log_q_cpue_eco(4);
  
  const double log_avg_F_cRn_LO=set_log_avg_F_L_cRn(2); const double log_avg_F_cRn_HI=set_log_avg_F_L_cRn(3); const double log_avg_F_cRn_PH=set_log_avg_F_L_cRn(4);
  const double log_avg_F_cRs_LO=set_log_avg_F_L_cRs(2); const double log_avg_F_cRs_HI=set_log_avg_F_L_cRs(3); const double log_avg_F_cRs_PH=set_log_avg_F_L_cRs(4);
  const double log_avg_F_cBn_LO=set_log_avg_F_L_cBn(2); const double log_avg_F_cBn_HI=set_log_avg_F_L_cBn(3); const double log_avg_F_cBn_PH=set_log_avg_F_L_cBn(4);
  const double log_avg_F_cBs_LO=set_log_avg_F_L_cBs(2); const double log_avg_F_cBs_HI=set_log_avg_F_L_cBs(3); const double log_avg_F_cBs_PH=set_log_avg_F_L_cBs(4);
  //-dev vectors-----------------------------------------------------------------------------------------------------------  
  const double log_F_dev_cRn_LO=set_log_dev_F_L_cRn(1); const double log_F_dev_cRn_HI=set_log_dev_F_L_cRn(2); const double log_F_dev_cRn_PH=set_log_dev_F_L_cRn(3);  
  const double log_F_dev_cRs_LO=set_log_dev_F_L_cRs(1); const double log_F_dev_cRs_HI=set_log_dev_F_L_cRs(2); const double log_F_dev_cRs_PH=set_log_dev_F_L_cRs(3);  
  const double log_F_dev_cBn_LO=set_log_dev_F_L_cBn(1); const double log_F_dev_cBn_HI=set_log_dev_F_L_cBn(2); const double log_F_dev_cBn_PH=set_log_dev_F_L_cBn(3);  
  const double log_F_dev_cBs_LO=set_log_dev_F_L_cBs(1); const double log_F_dev_cBs_HI=set_log_dev_F_L_cBs(2); const double log_F_dev_cBs_PH=set_log_dev_F_L_cBs(3);  
  const double log_RWq_LO=set_log_dev_RWq(1); const double log_RWq_HI=set_log_dev_RWq(2); const double log_RWq_PH=set_log_dev_RWq(3);  
  const double log_rec_dev_LO=set_log_dev_rec(1); const double log_rec_dev_HI=set_log_dev_rec(2); const double log_rec_dev_PH=set_log_dev_rec(3);          
  const double log_Nage_dev_LO=set_log_dev_Nage(1); const double log_Nage_dev_HI=set_log_dev_Nage(2); const double log_Nage_dev_PH=set_log_dev_Nage(3);          
  Linf.allocate(Linf_LO,Linf_HI,Linf_PH,"Linf");
  K.allocate(K_LO,K_HI,K_PH,"K");
  t0.allocate(t0_LO,t0_HI,t0_PH,"t0");
  len_cv_nad_val.allocate(len_cv_nad_LO,len_cv_nad_HI,len_cv_nad_PH,"len_cv_nad_val");
  len_cv_mad_val.allocate(len_cv_mad_LO,len_cv_mad_HI,len_cv_mad_PH,"len_cv_mad_val");
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
  len_cv_nad_val_out.allocate(1,8,"len_cv_nad_val_out");
  #ifndef NO_AD_INITIALIZE
    len_cv_nad_val_out.initialize();
  #endif
  len_cv_mad_val_out.allocate(1,8,"len_cv_mad_val_out");
  #ifndef NO_AD_INITIALIZE
    len_cv_mad_val_out.initialize();
  #endif
  meanlen_FL.allocate(1,nages,"meanlen_FL");
  #ifndef NO_AD_INITIALIZE
    meanlen_FL.initialize();
  #endif
  meanlen_FL_apr15.allocate(styr,endyr,1,nages,"meanlen_FL_apr15");
  #ifndef NO_AD_INITIALIZE
    meanlen_FL_apr15.initialize();
  #endif
  meanlen_FL_jun1.allocate(styr,endyr,1,nages,"meanlen_FL_jun1");
  #ifndef NO_AD_INITIALIZE
    meanlen_FL_jun1.initialize();
  #endif
  meanlen_FL_oct15.allocate(styr,endyr,1,nages,"meanlen_FL_oct15");
  #ifndef NO_AD_INITIALIZE
    meanlen_FL_oct15.initialize();
  #endif
  wgt_fish_mt.allocate(1,nages,"wgt_fish_mt");
  #ifndef NO_AD_INITIALIZE
    wgt_fish_mt.initialize();
  #endif
  wgt_spawn_mt.allocate(1,nages,"wgt_spawn_mt");
  #ifndef NO_AD_INITIALIZE
    wgt_spawn_mt.initialize();
  #endif
  wgt_spawn_mt_tv.allocate(styr,endyr,1,nages,"wgt_spawn_mt_tv");
  #ifndef NO_AD_INITIALIZE
    wgt_spawn_mt_tv.initialize();
  #endif
  tv_wgt_middle_mt.allocate(styr,endyr,1,nages,"tv_wgt_middle_mt");
  #ifndef NO_AD_INITIALIZE
    tv_wgt_middle_mt.initialize();
  #endif
  tv_wgt_spawn_mt.allocate(styr,endyr,1,nages,"tv_wgt_spawn_mt");
  #ifndef NO_AD_INITIALIZE
    tv_wgt_spawn_mt.initialize();
  #endif
  len_cR_mm.allocate(styr,endyr,1,nages,"len_cR_mm");
  #ifndef NO_AD_INITIALIZE
    len_cR_mm.initialize();
  #endif
  wholewgt_cR_mt.allocate(styr,endyr,1,nages,"wholewgt_cR_mt");
  #ifndef NO_AD_INITIALIZE
    wholewgt_cR_mt.initialize();
  #endif
  lenprob_apr15.allocate(styr,endyr,1,nages,1,nlenbins,"lenprob_apr15");
  #ifndef NO_AD_INITIALIZE
    lenprob_apr15.initialize();
  #endif
  lenprob_jun1.allocate(styr,endyr,1,nages,1,nlenbins,"lenprob_jun1");
  #ifndef NO_AD_INITIALIZE
    lenprob_jun1.initialize();
  #endif
  lenprob_oct15.allocate(styr,endyr,1,nages,1,nlenbins,"lenprob_oct15");
  #ifndef NO_AD_INITIALIZE
    lenprob_oct15.initialize();
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
  lenprob_nad.allocate(styr,endyr,1,nages,1,nlenbins,"lenprob_nad");
  #ifndef NO_AD_INITIALIZE
    lenprob_nad.initialize();
  #endif
  lenprob_mad.allocate(styr,endyr,1,nages,1,nlenbins,"lenprob_mad");
  #ifndef NO_AD_INITIALIZE
    lenprob_mad.initialize();
  #endif
  len_sd_nad.allocate(styr,endyr,1,nages,"len_sd_nad");
  #ifndef NO_AD_INITIALIZE
    len_sd_nad.initialize();
  #endif
  len_sd_mad.allocate(styr,endyr,1,nages,"len_sd_mad");
  #ifndef NO_AD_INITIALIZE
    len_sd_mad.initialize();
  #endif
  len_cv_nad.allocate(1,nages,"len_cv_nad");
  #ifndef NO_AD_INITIALIZE
    len_cv_nad.initialize();
  #endif
  len_cv_mad.allocate(1,nages,"len_cv_mad");
  #ifndef NO_AD_INITIALIZE
    len_cv_mad.initialize();
  #endif
  len_cv_apr15.allocate(1,nages,"len_cv_apr15");
  #ifndef NO_AD_INITIALIZE
    len_cv_apr15.initialize();
  #endif
  len_cv_jun1.allocate(1,nages,"len_cv_jun1");
  #ifndef NO_AD_INITIALIZE
    len_cv_jun1.initialize();
  #endif
  len_cv_oct15.allocate(1,nages,"len_cv_oct15");
  #ifndef NO_AD_INITIALIZE
    len_cv_oct15.initialize();
  #endif
  pred_nad_lenc.allocate(1,nyr_lenc_nad,1,nlenbins,"pred_nad_lenc");
  #ifndef NO_AD_INITIALIZE
    pred_nad_lenc.initialize();
  #endif
  pred_mad_lenc.allocate(1,nyr_lenc_mad,1,nlenbins,"pred_mad_lenc");
  #ifndef NO_AD_INITIALIZE
    pred_mad_lenc.initialize();
  #endif
  pred_cRn_agec.allocate(1,nyr_agec_cRn,1,nages_agec,"pred_cRn_agec");
  #ifndef NO_AD_INITIALIZE
    pred_cRn_agec.initialize();
  #endif
  pred_cRn_agec_allages.allocate(1,nyr_agec_cRn,1,nages,"pred_cRn_agec_allages");
  #ifndef NO_AD_INITIALIZE
    pred_cRn_agec_allages.initialize();
  #endif
  ErrorFree_cRn_agec.allocate(1,nyr_agec_cRn,1,nages,"ErrorFree_cRn_agec");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_cRn_agec.initialize();
  #endif
  pred_cRs_agec.allocate(1,nyr_agec_cRs,1,nages_agec,"pred_cRs_agec");
  #ifndef NO_AD_INITIALIZE
    pred_cRs_agec.initialize();
  #endif
  pred_cRs_agec_allages.allocate(1,nyr_agec_cRs,1,nages,"pred_cRs_agec_allages");
  #ifndef NO_AD_INITIALIZE
    pred_cRs_agec_allages.initialize();
  #endif
  ErrorFree_cRs_agec.allocate(1,nyr_agec_cRs,1,nages,"ErrorFree_cRs_agec");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_cRs_agec.initialize();
  #endif
  pred_cBn_agec.allocate(1,nyr_agec_cBn,1,nages_agec,"pred_cBn_agec");
  #ifndef NO_AD_INITIALIZE
    pred_cBn_agec.initialize();
  #endif
  pred_cBn_agec_allages.allocate(1,nyr_agec_cBn,1,nages,"pred_cBn_agec_allages");
  #ifndef NO_AD_INITIALIZE
    pred_cBn_agec_allages.initialize();
  #endif
  ErrorFree_cBn_agec.allocate(1,nyr_agec_cBn,1,nages,"ErrorFree_cBn_agec");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_cBn_agec.initialize();
  #endif
  pred_cBs_agec.allocate(1,nyr_agec_cBs,1,nages_agec,"pred_cBs_agec");
  #ifndef NO_AD_INITIALIZE
    pred_cBs_agec.initialize();
  #endif
  pred_cBs_agec_allages.allocate(1,nyr_agec_cBs,1,nages,"pred_cBs_agec_allages");
  #ifndef NO_AD_INITIALIZE
    pred_cBs_agec_allages.initialize();
  #endif
  ErrorFree_cBs_agec.allocate(1,nyr_agec_cBs,1,nages,"ErrorFree_cBs_agec");
  #ifndef NO_AD_INITIALIZE
    ErrorFree_cBs_agec.initialize();
  #endif
  nsamp_nad_lenc_allyr.allocate(styr,endyr,"nsamp_nad_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_nad_lenc_allyr.initialize();
  #endif
  nsamp_mad_lenc_allyr.allocate(styr,endyr,"nsamp_mad_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_mad_lenc_allyr.initialize();
  #endif
  nsamp_cRn_agec_allyr.allocate(styr,endyr,"nsamp_cRn_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_cRn_agec_allyr.initialize();
  #endif
  nsamp_cRs_agec_allyr.allocate(styr,endyr,"nsamp_cRs_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_cRs_agec_allyr.initialize();
  #endif
  nsamp_cBn_agec_allyr.allocate(styr,endyr,"nsamp_cBn_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_cBn_agec_allyr.initialize();
  #endif
  nsamp_cBs_agec_allyr.allocate(styr,endyr,"nsamp_cBs_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nsamp_cBs_agec_allyr.initialize();
  #endif
  nfish_nad_lenc_allyr.allocate(styr,endyr,"nfish_nad_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_nad_lenc_allyr.initialize();
  #endif
  nfish_mad_lenc_allyr.allocate(styr,endyr,"nfish_mad_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_mad_lenc_allyr.initialize();
  #endif
  nfish_cRs_agec_allyr.allocate(styr,endyr,"nfish_cRs_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_cRs_agec_allyr.initialize();
  #endif
  nfish_cRn_agec_allyr.allocate(styr,endyr,"nfish_cRn_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_cRn_agec_allyr.initialize();
  #endif
  nfish_cBn_agec_allyr.allocate(styr,endyr,"nfish_cBn_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_cBn_agec_allyr.initialize();
  #endif
  nfish_cBs_agec_allyr.allocate(styr,endyr,"nfish_cBs_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    nfish_cBs_agec_allyr.initialize();
  #endif
  neff_nad_lenc_allyr.allocate(styr,endyr,"neff_nad_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    neff_nad_lenc_allyr.initialize();
  #endif
  neff_mad_lenc_allyr.allocate(styr,endyr,"neff_mad_lenc_allyr");
  #ifndef NO_AD_INITIALIZE
    neff_mad_lenc_allyr.initialize();
  #endif
  neff_cRn_agec_allyr.allocate(styr,endyr,"neff_cRn_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    neff_cRn_agec_allyr.initialize();
  #endif
  neff_cRs_agec_allyr.allocate(styr,endyr,"neff_cRs_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    neff_cRs_agec_allyr.initialize();
  #endif
  neff_cBn_agec_allyr.allocate(styr,endyr,"neff_cBn_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    neff_cBn_agec_allyr.initialize();
  #endif
  neff_cBs_agec_allyr.allocate(styr,endyr,"neff_cBs_agec_allyr");
  #ifndef NO_AD_INITIALIZE
    neff_cBs_agec_allyr.initialize();
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
  N_nad.allocate(styr,endyr,1,nages,"N_nad");
  #ifndef NO_AD_INITIALIZE
    N_nad.initialize();
  #endif
  N_mad.allocate(styr,endyr,1,nages,"N_mad");
  #ifndef NO_AD_INITIALIZE
    N_mad.initialize();
  #endif
  N_sad.allocate(styr,endyr,1,nages,"N_sad");
  #ifndef NO_AD_INITIALIZE
    N_sad.initialize();
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
  B_mdyr.allocate(styr,endyr+1,1,nages,"B_mdyr");
  #ifndef NO_AD_INITIALIZE
    B_mdyr.initialize();
  #endif
  totB.allocate(styr,endyr+1,"totB");
  #ifndef NO_AD_INITIALIZE
    totB.initialize();
  #endif
  totN.allocate(styr,endyr+1,"totN");
  #ifndef NO_AD_INITIALIZE
    totN.initialize();
  #endif
  SSB.allocate(styr,endyr+1,"SSB");
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
  maturity_f.allocate(1,nages,"maturity_f");
  #ifndef NO_AD_INITIALIZE
    maturity_f.initialize();
  #endif
  maturity_m.allocate(1,nages,"maturity_m");
  #ifndef NO_AD_INITIALIZE
    maturity_m.initialize();
  #endif
  tv_maturity.allocate(styr,endyr,1,nages,"tv_maturity");
  #ifndef NO_AD_INITIALIZE
    tv_maturity.initialize();
  #endif
  reprod.allocate(1,nages,"reprod");
  #ifndef NO_AD_INITIALIZE
    reprod.initialize();
  #endif
  reprod_tv.allocate(styr,endyr,1,nages,"reprod_tv");
  #ifndef NO_AD_INITIALIZE
    reprod_tv.initialize();
  #endif
  SSBatage.allocate(styr,endyr,1,nages,"SSBatage");
  #ifndef NO_AD_INITIALIZE
    SSBatage.initialize();
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
  SdS0.allocate(styr,endyr+1,"SdS0");
  #ifndef NO_AD_INITIALIZE
    SdS0.initialize();
  #endif
  log_dm_lenc_nad.allocate(log_dm_nad_lc_LO,log_dm_nad_lc_HI,log_dm_nad_lc_PH,"log_dm_lenc_nad");
  log_dm_lenc_mad.allocate(log_dm_mad_lc_LO,log_dm_mad_lc_HI,log_dm_mad_lc_PH,"log_dm_lenc_mad");
  log_dm_agec_cRn.allocate(log_dm_cRn_ac_LO,log_dm_cRn_ac_HI,log_dm_cRn_ac_PH,"log_dm_agec_cRn");
  log_dm_agec_cRs.allocate(log_dm_cRs_ac_LO,log_dm_cRs_ac_HI,log_dm_cRs_ac_PH,"log_dm_agec_cRs");
  log_dm_agec_cBn.allocate(log_dm_cBn_ac_LO,log_dm_cBn_ac_HI,log_dm_cBn_ac_PH,"log_dm_agec_cBn");
  log_dm_agec_cBs.allocate(log_dm_cBs_ac_LO,log_dm_cBs_ac_HI,log_dm_cBs_ac_PH,"log_dm_agec_cBs");
  log_dm_nad_lc_out.allocate(1,8,"log_dm_nad_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_nad_lc_out.initialize();
  #endif
  log_dm_mad_lc_out.allocate(1,8,"log_dm_mad_lc_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_mad_lc_out.initialize();
  #endif
  log_dm_cRn_ac_out.allocate(1,8,"log_dm_cRn_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_cRn_ac_out.initialize();
  #endif
  log_dm_cRs_ac_out.allocate(1,8,"log_dm_cRs_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_cRs_ac_out.initialize();
  #endif
  log_dm_cBn_ac_out.allocate(1,8,"log_dm_cBn_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_cBn_ac_out.initialize();
  #endif
  log_dm_cBs_ac_out.allocate(1,8,"log_dm_cBs_ac_out");
  #ifndef NO_AD_INITIALIZE
    log_dm_cBs_ac_out.initialize();
  #endif
  sel_cRn.allocate(styr,endyr,1,nages,"sel_cRn");
  #ifndef NO_AD_INITIALIZE
    sel_cRn.initialize();
  #endif
  sel_cRs.allocate(styr,endyr,1,nages,"sel_cRs");
  #ifndef NO_AD_INITIALIZE
    sel_cRs.initialize();
  #endif
  sel_cRn_block1.allocate(1,nages,"sel_cRn_block1");
  #ifndef NO_AD_INITIALIZE
    sel_cRn_block1.initialize();
  #endif
  sel_cRn_block2.allocate(1,nages,"sel_cRn_block2");
  #ifndef NO_AD_INITIALIZE
    sel_cRn_block2.initialize();
  #endif
  sel_cRn_block3.allocate(1,nages,"sel_cRn_block3");
  #ifndef NO_AD_INITIALIZE
    sel_cRn_block3.initialize();
  #endif
  sel_cRn_block4.allocate(1,nages,"sel_cRn_block4");
  #ifndef NO_AD_INITIALIZE
    sel_cRn_block4.initialize();
  #endif
  sel_cRs_block1.allocate(1,nages,"sel_cRs_block1");
  #ifndef NO_AD_INITIALIZE
    sel_cRs_block1.initialize();
  #endif
  sel_cRs_block2.allocate(1,nages,"sel_cRs_block2");
  #ifndef NO_AD_INITIALIZE
    sel_cRs_block2.initialize();
  #endif
  sel_cRs_block3.allocate(1,nages,"sel_cRs_block3");
  #ifndef NO_AD_INITIALIZE
    sel_cRs_block3.initialize();
  #endif
  sel_cRs_block4.allocate(1,nages,"sel_cRs_block4");
  #ifndef NO_AD_INITIALIZE
    sel_cRs_block4.initialize();
  #endif
  selpar_A50_cRn.allocate(selpar_A50_cRn_LO,selpar_A50_cRn_HI,selpar_A50_cRn_PH,"selpar_A50_cRn");
  selpar_slope_cRn.allocate(selpar_slope_cRn_LO,selpar_slope_cRn_HI,selpar_slope_cRn_PH,"selpar_slope_cRn");
  selpar_A502_cRn.allocate(selpar_A502_cRn_LO,selpar_A502_cRn_HI,selpar_A502_cRn_PH,"selpar_A502_cRn");
  selpar_slope2_cRn.allocate(selpar_slope2_cRn_LO,selpar_slope2_cRn_HI,selpar_slope2_cRn_PH,"selpar_slope2_cRn");
  selpar_A50_cRn_out.allocate(1,8,"selpar_A50_cRn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cRn_out.initialize();
  #endif
  selpar_slope_cRn_out.allocate(1,8,"selpar_slope_cRn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cRn_out.initialize();
  #endif
  selpar_A502_cRn_out.allocate(1,8,"selpar_A502_cRn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_cRn_out.initialize();
  #endif
  selpar_slope2_cRn_out.allocate(1,8,"selpar_slope2_cRn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_cRn_out.initialize();
  #endif
  selpar_A50_cRn2.allocate(selpar_A50_cRn2_LO,selpar_A50_cRn2_HI,selpar_A50_cRn2_PH,"selpar_A50_cRn2");
  selpar_slope_cRn2.allocate(selpar_slope_cRn2_LO,selpar_slope_cRn2_HI,selpar_slope_cRn2_PH,"selpar_slope_cRn2");
  selpar_A502_cRn2.allocate(selpar_A502_cRn2_LO,selpar_A502_cRn2_HI,selpar_A502_cRn2_PH,"selpar_A502_cRn2");
  selpar_slope2_cRn2.allocate(selpar_slope2_cRn2_LO,selpar_slope2_cRn2_HI,selpar_slope2_cRn2_PH,"selpar_slope2_cRn2");
  selpar_A50_cRn2_out.allocate(1,8,"selpar_A50_cRn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cRn2_out.initialize();
  #endif
  selpar_slope_cRn2_out.allocate(1,8,"selpar_slope_cRn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cRn2_out.initialize();
  #endif
  selpar_A502_cRn2_out.allocate(1,8,"selpar_A502_cRn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_cRn2_out.initialize();
  #endif
  selpar_slope2_cRn2_out.allocate(1,8,"selpar_slope2_cRn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_cRn2_out.initialize();
  #endif
  selpar_A50_cRn3.allocate(selpar_A50_cRn3_LO,selpar_A50_cRn3_HI,selpar_A50_cRn3_PH,"selpar_A50_cRn3");
  selpar_slope_cRn3.allocate(selpar_slope_cRn3_LO,selpar_slope_cRn3_HI,selpar_slope_cRn3_PH,"selpar_slope_cRn3");
  selpar_A502_cRn3.allocate(selpar_A502_cRn3_LO,selpar_A502_cRn3_HI,selpar_A502_cRn3_PH,"selpar_A502_cRn3");
  selpar_slope2_cRn3.allocate(selpar_slope2_cRn3_LO,selpar_slope2_cRn3_HI,selpar_slope2_cRn3_PH,"selpar_slope2_cRn3");
  selpar_A50_cRn3_out.allocate(1,8,"selpar_A50_cRn3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cRn3_out.initialize();
  #endif
  selpar_slope_cRn3_out.allocate(1,8,"selpar_slope_cRn3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cRn3_out.initialize();
  #endif
  selpar_A502_cRn3_out.allocate(1,8,"selpar_A502_cRn3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_cRn3_out.initialize();
  #endif
  selpar_slope2_cRn3_out.allocate(1,8,"selpar_slope2_cRn3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_cRn3_out.initialize();
  #endif
  selpar_A50_cRn4.allocate(selpar_A50_cRn4_LO,selpar_A50_cRn4_HI,selpar_A50_cRn4_PH,"selpar_A50_cRn4");
  selpar_slope_cRn4.allocate(selpar_slope_cRn4_LO,selpar_slope_cRn4_HI,selpar_slope_cRn4_PH,"selpar_slope_cRn4");
  selpar_A502_cRn4.allocate(selpar_A502_cRn4_LO,selpar_A502_cRn4_HI,selpar_A502_cRn4_PH,"selpar_A502_cRn4");
  selpar_slope2_cRn4.allocate(selpar_slope2_cRn4_LO,selpar_slope2_cRn4_HI,selpar_slope2_cRn4_PH,"selpar_slope2_cRn4");
  selpar_A50_cRn4_out.allocate(1,8,"selpar_A50_cRn4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cRn4_out.initialize();
  #endif
  selpar_slope_cRn4_out.allocate(1,8,"selpar_slope_cRn4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cRn4_out.initialize();
  #endif
  selpar_A502_cRn4_out.allocate(1,8,"selpar_A502_cRn4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_cRn4_out.initialize();
  #endif
  selpar_slope2_cRn4_out.allocate(1,8,"selpar_slope2_cRn4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_cRn4_out.initialize();
  #endif
  selpar_A50_cRs.allocate(selpar_A50_cRs_LO,selpar_A50_cRs_HI,selpar_A50_cRs_PH,"selpar_A50_cRs");
  selpar_slope_cRs.allocate(selpar_slope_cRs_LO,selpar_slope_cRs_HI,selpar_slope_cRs_PH,"selpar_slope_cRs");
  selpar_A502_cRs.allocate(selpar_A502_cRs_LO,selpar_A502_cRs_HI,selpar_A502_cRs_PH,"selpar_A502_cRs");
  selpar_slope2_cRs.allocate(selpar_slope2_cRs_LO,selpar_slope2_cRs_HI,selpar_slope2_cRs_PH,"selpar_slope2_cRs");
  selpar_A50_cRs_out.allocate(1,8,"selpar_A50_cRs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cRs_out.initialize();
  #endif
  selpar_slope_cRs_out.allocate(1,8,"selpar_slope_cRs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cRs_out.initialize();
  #endif
  selpar_A502_cRs_out.allocate(1,8,"selpar_A502_cRs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_cRs_out.initialize();
  #endif
  selpar_slope2_cRs_out.allocate(1,8,"selpar_slope2_cRs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_cRs_out.initialize();
  #endif
  selpar_A50_cRs2.allocate(selpar_A50_cRs2_LO,selpar_A50_cRs2_HI,selpar_A50_cRs2_PH,"selpar_A50_cRs2");
  selpar_slope_cRs2.allocate(selpar_slope_cRs2_LO,selpar_slope_cRs2_HI,selpar_slope_cRs2_PH,"selpar_slope_cRs2");
  selpar_A502_cRs2.allocate(selpar_A502_cRs2_LO,selpar_A502_cRs2_HI,selpar_A502_cRs2_PH,"selpar_A502_cRs2");
  selpar_slope2_cRs2.allocate(selpar_slope2_cRs2_LO,selpar_slope2_cRs2_HI,selpar_slope2_cRs2_PH,"selpar_slope2_cRs2");
  selpar_A50_cRs2_out.allocate(1,8,"selpar_A50_cRs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cRs2_out.initialize();
  #endif
  selpar_slope_cRs2_out.allocate(1,8,"selpar_slope_cRs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cRs2_out.initialize();
  #endif
  selpar_A502_cRs2_out.allocate(1,8,"selpar_A502_cRs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_cRs2_out.initialize();
  #endif
  selpar_slope2_cRs2_out.allocate(1,8,"selpar_slope2_cRs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_cRs2_out.initialize();
  #endif
  selpar_A50_cRs3.allocate(selpar_A50_cRs3_LO,selpar_A50_cRs3_HI,selpar_A50_cRs3_PH,"selpar_A50_cRs3");
  selpar_slope_cRs3.allocate(selpar_slope_cRs3_LO,selpar_slope_cRs3_HI,selpar_slope_cRs3_PH,"selpar_slope_cRs3");
  selpar_A502_cRs3.allocate(selpar_A502_cRs3_LO,selpar_A502_cRs3_HI,selpar_A502_cRs3_PH,"selpar_A502_cRs3");
  selpar_slope2_cRs3.allocate(selpar_slope2_cRs3_LO,selpar_slope2_cRs3_HI,selpar_slope2_cRs3_PH,"selpar_slope2_cRs3");
  selpar_A50_cRs3_out.allocate(1,8,"selpar_A50_cRs3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cRs3_out.initialize();
  #endif
  selpar_slope_cRs3_out.allocate(1,8,"selpar_slope_cRs3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cRs3_out.initialize();
  #endif
  selpar_A502_cRs3_out.allocate(1,8,"selpar_A502_cRs3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_cRs3_out.initialize();
  #endif
  selpar_slope2_cRs3_out.allocate(1,8,"selpar_slope2_cRs3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_cRs3_out.initialize();
  #endif
  selpar_A50_cRs4.allocate(selpar_A50_cRs4_LO,selpar_A50_cRs4_HI,selpar_A50_cRs4_PH,"selpar_A50_cRs4");
  selpar_slope_cRs4.allocate(selpar_slope_cRs4_LO,selpar_slope_cRs4_HI,selpar_slope_cRs4_PH,"selpar_slope_cRs4");
  selpar_A502_cRs4.allocate(selpar_A502_cRs4_LO,selpar_A502_cRs4_HI,selpar_A502_cRs4_PH,"selpar_A502_cRs4");
  selpar_slope2_cRs4.allocate(selpar_slope2_cRs4_LO,selpar_slope2_cRs4_HI,selpar_slope2_cRs4_PH,"selpar_slope2_cRs4");
  selpar_A50_cRs4_out.allocate(1,8,"selpar_A50_cRs4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cRs4_out.initialize();
  #endif
  selpar_slope_cRs4_out.allocate(1,8,"selpar_slope_cRs4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cRs4_out.initialize();
  #endif
  selpar_A502_cRs4_out.allocate(1,8,"selpar_A502_cRs4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_cRs4_out.initialize();
  #endif
  selpar_slope2_cRs4_out.allocate(1,8,"selpar_slope2_cRs4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_cRs4_out.initialize();
  #endif
  sel_age0_cRn_logit.allocate(selpar_age0_cRn_LO,selpar_age0_cRn_HI,selpar_age0_cRn_PH,"sel_age0_cRn_logit");
  sel_age1_cRn_logit.allocate(selpar_age1_cRn_LO,selpar_age1_cRn_HI,selpar_age1_cRn_PH,"sel_age1_cRn_logit");
  sel_age2_cRn_logit.allocate(selpar_age2_cRn_LO,selpar_age2_cRn_HI,selpar_age2_cRn_PH,"sel_age2_cRn_logit");
  sel_age3_cRn_logit.allocate(selpar_age3_cRn_LO,selpar_age3_cRn_HI,selpar_age3_cRn_PH,"sel_age3_cRn_logit");
  sel_age4_cRn_logit.allocate(selpar_age4_cRn_LO,selpar_age4_cRn_HI,selpar_age4_cRn_PH,"sel_age4_cRn_logit");
  sel_age5_cRn_logit.allocate(selpar_age5_cRn_LO,selpar_age5_cRn_HI,selpar_age5_cRn_PH,"sel_age5_cRn_logit");
  sel_age6_cRn_logit.allocate(selpar_age6_cRn_LO,selpar_age6_cRn_HI,selpar_age6_cRn_PH,"sel_age6_cRn_logit");
  sel_age_cRn_vec.allocate(1,nages,"sel_age_cRn_vec");
  #ifndef NO_AD_INITIALIZE
    sel_age_cRn_vec.initialize();
  #endif
  selpar_age0_cRn.allocate("selpar_age0_cRn");
  #ifndef NO_AD_INITIALIZE
  selpar_age0_cRn.initialize();
  #endif
  selpar_age1_cRn.allocate("selpar_age1_cRn");
  #ifndef NO_AD_INITIALIZE
  selpar_age1_cRn.initialize();
  #endif
  selpar_age2_cRn.allocate("selpar_age2_cRn");
  #ifndef NO_AD_INITIALIZE
  selpar_age2_cRn.initialize();
  #endif
  selpar_age3_cRn.allocate("selpar_age3_cRn");
  #ifndef NO_AD_INITIALIZE
  selpar_age3_cRn.initialize();
  #endif
  selpar_age4_cRn.allocate("selpar_age4_cRn");
  #ifndef NO_AD_INITIALIZE
  selpar_age4_cRn.initialize();
  #endif
  selpar_age5_cRn.allocate("selpar_age5_cRn");
  #ifndef NO_AD_INITIALIZE
  selpar_age5_cRn.initialize();
  #endif
  selpar_age6_cRn.allocate("selpar_age6_cRn");
  #ifndef NO_AD_INITIALIZE
  selpar_age6_cRn.initialize();
  #endif
  selpar_age0_cRn_out.allocate(1,8,"selpar_age0_cRn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age0_cRn_out.initialize();
  #endif
  selpar_age1_cRn_out.allocate(1,8,"selpar_age1_cRn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age1_cRn_out.initialize();
  #endif
  selpar_age2_cRn_out.allocate(1,8,"selpar_age2_cRn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age2_cRn_out.initialize();
  #endif
  selpar_age3_cRn_out.allocate(1,8,"selpar_age3_cRn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age3_cRn_out.initialize();
  #endif
  selpar_age4_cRn_out.allocate(1,8,"selpar_age4_cRn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age4_cRn_out.initialize();
  #endif
  selpar_age5_cRn_out.allocate(1,8,"selpar_age5_cRn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age5_cRn_out.initialize();
  #endif
  selpar_age6_cRn_out.allocate(1,8,"selpar_age6_cRn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age6_cRn_out.initialize();
  #endif
  sel_age0_cRn2_logit.allocate(selpar_age0_cRn2_LO,selpar_age0_cRn2_HI,selpar_age0_cRn2_PH,"sel_age0_cRn2_logit");
  sel_age1_cRn2_logit.allocate(selpar_age1_cRn2_LO,selpar_age1_cRn2_HI,selpar_age1_cRn2_PH,"sel_age1_cRn2_logit");
  sel_age2_cRn2_logit.allocate(selpar_age2_cRn2_LO,selpar_age2_cRn2_HI,selpar_age2_cRn2_PH,"sel_age2_cRn2_logit");
  sel_age3_cRn2_logit.allocate(selpar_age3_cRn2_LO,selpar_age3_cRn2_HI,selpar_age3_cRn2_PH,"sel_age3_cRn2_logit");
  sel_age4_cRn2_logit.allocate(selpar_age4_cRn2_LO,selpar_age4_cRn2_HI,selpar_age4_cRn2_PH,"sel_age4_cRn2_logit");
  sel_age5_cRn2_logit.allocate(selpar_age5_cRn2_LO,selpar_age5_cRn2_HI,selpar_age5_cRn2_PH,"sel_age5_cRn2_logit");
  sel_age6_cRn2_logit.allocate(selpar_age6_cRn2_LO,selpar_age6_cRn2_HI,selpar_age6_cRn2_PH,"sel_age6_cRn2_logit");
  sel_age_cRn2_vec.allocate(1,nages,"sel_age_cRn2_vec");
  #ifndef NO_AD_INITIALIZE
    sel_age_cRn2_vec.initialize();
  #endif
  selpar_age0_cRn2.allocate("selpar_age0_cRn2");
  #ifndef NO_AD_INITIALIZE
  selpar_age0_cRn2.initialize();
  #endif
  selpar_age1_cRn2.allocate("selpar_age1_cRn2");
  #ifndef NO_AD_INITIALIZE
  selpar_age1_cRn2.initialize();
  #endif
  selpar_age2_cRn2.allocate("selpar_age2_cRn2");
  #ifndef NO_AD_INITIALIZE
  selpar_age2_cRn2.initialize();
  #endif
  selpar_age3_cRn2.allocate("selpar_age3_cRn2");
  #ifndef NO_AD_INITIALIZE
  selpar_age3_cRn2.initialize();
  #endif
  selpar_age4_cRn2.allocate("selpar_age4_cRn2");
  #ifndef NO_AD_INITIALIZE
  selpar_age4_cRn2.initialize();
  #endif
  selpar_age5_cRn2.allocate("selpar_age5_cRn2");
  #ifndef NO_AD_INITIALIZE
  selpar_age5_cRn2.initialize();
  #endif
  selpar_age6_cRn2.allocate("selpar_age6_cRn2");
  #ifndef NO_AD_INITIALIZE
  selpar_age6_cRn2.initialize();
  #endif
  selpar_age0_cRn2_out.allocate(1,8,"selpar_age0_cRn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age0_cRn2_out.initialize();
  #endif
  selpar_age1_cRn2_out.allocate(1,8,"selpar_age1_cRn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age1_cRn2_out.initialize();
  #endif
  selpar_age2_cRn2_out.allocate(1,8,"selpar_age2_cRn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age2_cRn2_out.initialize();
  #endif
  selpar_age3_cRn2_out.allocate(1,8,"selpar_age3_cRn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age3_cRn2_out.initialize();
  #endif
  selpar_age4_cRn2_out.allocate(1,8,"selpar_age4_cRn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age4_cRn2_out.initialize();
  #endif
  selpar_age5_cRn2_out.allocate(1,8,"selpar_age5_cRn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age5_cRn2_out.initialize();
  #endif
  selpar_age6_cRn2_out.allocate(1,8,"selpar_age6_cRn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age6_cRn2_out.initialize();
  #endif
  sel_age0_cRn3_logit.allocate(selpar_age0_cRn3_LO,selpar_age0_cRn3_HI,selpar_age0_cRn3_PH,"sel_age0_cRn3_logit");
  sel_age1_cRn3_logit.allocate(selpar_age1_cRn3_LO,selpar_age1_cRn3_HI,selpar_age1_cRn3_PH,"sel_age1_cRn3_logit");
  sel_age2_cRn3_logit.allocate(selpar_age2_cRn3_LO,selpar_age2_cRn3_HI,selpar_age2_cRn3_PH,"sel_age2_cRn3_logit");
  sel_age3_cRn3_logit.allocate(selpar_age3_cRn3_LO,selpar_age3_cRn3_HI,selpar_age3_cRn3_PH,"sel_age3_cRn3_logit");
  sel_age4_cRn3_logit.allocate(selpar_age4_cRn3_LO,selpar_age4_cRn3_HI,selpar_age4_cRn3_PH,"sel_age4_cRn3_logit");
  sel_age5_cRn3_logit.allocate(selpar_age5_cRn3_LO,selpar_age5_cRn3_HI,selpar_age5_cRn3_PH,"sel_age5_cRn3_logit");
  sel_age6_cRn3_logit.allocate(selpar_age6_cRn3_LO,selpar_age6_cRn3_HI,selpar_age6_cRn3_PH,"sel_age6_cRn3_logit");
  sel_age_cRn3_vec.allocate(1,nages,"sel_age_cRn3_vec");
  #ifndef NO_AD_INITIALIZE
    sel_age_cRn3_vec.initialize();
  #endif
  selpar_age0_cRn3.allocate("selpar_age0_cRn3");
  #ifndef NO_AD_INITIALIZE
  selpar_age0_cRn3.initialize();
  #endif
  selpar_age1_cRn3.allocate("selpar_age1_cRn3");
  #ifndef NO_AD_INITIALIZE
  selpar_age1_cRn3.initialize();
  #endif
  selpar_age2_cRn3.allocate("selpar_age2_cRn3");
  #ifndef NO_AD_INITIALIZE
  selpar_age2_cRn3.initialize();
  #endif
  selpar_age3_cRn3.allocate("selpar_age3_cRn3");
  #ifndef NO_AD_INITIALIZE
  selpar_age3_cRn3.initialize();
  #endif
  selpar_age4_cRn3.allocate("selpar_age4_cRn3");
  #ifndef NO_AD_INITIALIZE
  selpar_age4_cRn3.initialize();
  #endif
  selpar_age5_cRn3.allocate("selpar_age5_cRn3");
  #ifndef NO_AD_INITIALIZE
  selpar_age5_cRn3.initialize();
  #endif
  selpar_age6_cRn3.allocate("selpar_age6_cRn3");
  #ifndef NO_AD_INITIALIZE
  selpar_age6_cRn3.initialize();
  #endif
  selpar_age0_cRn3_out.allocate(1,8,"selpar_age0_cRn3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age0_cRn3_out.initialize();
  #endif
  selpar_age1_cRn3_out.allocate(1,8,"selpar_age1_cRn3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age1_cRn3_out.initialize();
  #endif
  selpar_age2_cRn3_out.allocate(1,8,"selpar_age2_cRn3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age2_cRn3_out.initialize();
  #endif
  selpar_age3_cRn3_out.allocate(1,8,"selpar_age3_cRn3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age3_cRn3_out.initialize();
  #endif
  selpar_age4_cRn3_out.allocate(1,8,"selpar_age4_cRn3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age4_cRn3_out.initialize();
  #endif
  selpar_age5_cRn3_out.allocate(1,8,"selpar_age5_cRn3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age5_cRn3_out.initialize();
  #endif
  selpar_age6_cRn3_out.allocate(1,8,"selpar_age6_cRn3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age6_cRn3_out.initialize();
  #endif
  sel_age0_cRn4_logit.allocate(selpar_age0_cRn4_LO,selpar_age0_cRn4_HI,selpar_age0_cRn4_PH,"sel_age0_cRn4_logit");
  sel_age1_cRn4_logit.allocate(selpar_age1_cRn4_LO,selpar_age1_cRn4_HI,selpar_age1_cRn4_PH,"sel_age1_cRn4_logit");
  sel_age2_cRn4_logit.allocate(selpar_age2_cRn4_LO,selpar_age2_cRn4_HI,selpar_age2_cRn4_PH,"sel_age2_cRn4_logit");
  sel_age3_cRn4_logit.allocate(selpar_age3_cRn4_LO,selpar_age3_cRn4_HI,selpar_age3_cRn4_PH,"sel_age3_cRn4_logit");
  sel_age4_cRn4_logit.allocate(selpar_age4_cRn4_LO,selpar_age4_cRn4_HI,selpar_age4_cRn4_PH,"sel_age4_cRn4_logit");
  sel_age5_cRn4_logit.allocate(selpar_age5_cRn4_LO,selpar_age5_cRn4_HI,selpar_age5_cRn4_PH,"sel_age5_cRn4_logit");
  sel_age6_cRn4_logit.allocate(selpar_age6_cRn4_LO,selpar_age6_cRn4_HI,selpar_age6_cRn4_PH,"sel_age6_cRn4_logit");
  sel_age_cRn4_vec.allocate(1,nages,"sel_age_cRn4_vec");
  #ifndef NO_AD_INITIALIZE
    sel_age_cRn4_vec.initialize();
  #endif
  selpar_age0_cRn4.allocate("selpar_age0_cRn4");
  #ifndef NO_AD_INITIALIZE
  selpar_age0_cRn4.initialize();
  #endif
  selpar_age1_cRn4.allocate("selpar_age1_cRn4");
  #ifndef NO_AD_INITIALIZE
  selpar_age1_cRn4.initialize();
  #endif
  selpar_age2_cRn4.allocate("selpar_age2_cRn4");
  #ifndef NO_AD_INITIALIZE
  selpar_age2_cRn4.initialize();
  #endif
  selpar_age3_cRn4.allocate("selpar_age3_cRn4");
  #ifndef NO_AD_INITIALIZE
  selpar_age3_cRn4.initialize();
  #endif
  selpar_age4_cRn4.allocate("selpar_age4_cRn4");
  #ifndef NO_AD_INITIALIZE
  selpar_age4_cRn4.initialize();
  #endif
  selpar_age5_cRn4.allocate("selpar_age5_cRn4");
  #ifndef NO_AD_INITIALIZE
  selpar_age5_cRn4.initialize();
  #endif
  selpar_age6_cRn4.allocate("selpar_age6_cRn4");
  #ifndef NO_AD_INITIALIZE
  selpar_age6_cRn4.initialize();
  #endif
  selpar_age0_cRn4_out.allocate(1,8,"selpar_age0_cRn4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age0_cRn4_out.initialize();
  #endif
  selpar_age1_cRn4_out.allocate(1,8,"selpar_age1_cRn4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age1_cRn4_out.initialize();
  #endif
  selpar_age2_cRn4_out.allocate(1,8,"selpar_age2_cRn4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age2_cRn4_out.initialize();
  #endif
  selpar_age3_cRn4_out.allocate(1,8,"selpar_age3_cRn4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age3_cRn4_out.initialize();
  #endif
  selpar_age4_cRn4_out.allocate(1,8,"selpar_age4_cRn4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age4_cRn4_out.initialize();
  #endif
  selpar_age5_cRn4_out.allocate(1,8,"selpar_age5_cRn4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age5_cRn4_out.initialize();
  #endif
  selpar_age6_cRn4_out.allocate(1,8,"selpar_age6_cRn4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age6_cRn4_out.initialize();
  #endif
  sel_age0_cRs_logit.allocate(selpar_age0_cRs_LO,selpar_age0_cRs_HI,selpar_age0_cRs_PH,"sel_age0_cRs_logit");
  sel_age1_cRs_logit.allocate(selpar_age1_cRs_LO,selpar_age1_cRs_HI,selpar_age1_cRs_PH,"sel_age1_cRs_logit");
  sel_age2_cRs_logit.allocate(selpar_age2_cRs_LO,selpar_age2_cRs_HI,selpar_age2_cRs_PH,"sel_age2_cRs_logit");
  sel_age3_cRs_logit.allocate(selpar_age3_cRs_LO,selpar_age3_cRs_HI,selpar_age3_cRs_PH,"sel_age3_cRs_logit");
  sel_age4_cRs_logit.allocate(selpar_age4_cRs_LO,selpar_age4_cRs_HI,selpar_age4_cRs_PH,"sel_age4_cRs_logit");
  sel_age5_cRs_logit.allocate(selpar_age5_cRs_LO,selpar_age5_cRs_HI,selpar_age5_cRs_PH,"sel_age5_cRs_logit");
  sel_age6_cRs_logit.allocate(selpar_age6_cRs_LO,selpar_age6_cRs_HI,selpar_age6_cRs_PH,"sel_age6_cRs_logit");
  sel_age_cRs_vec.allocate(1,nages,"sel_age_cRs_vec");
  #ifndef NO_AD_INITIALIZE
    sel_age_cRs_vec.initialize();
  #endif
  selpar_age0_cRs.allocate("selpar_age0_cRs");
  #ifndef NO_AD_INITIALIZE
  selpar_age0_cRs.initialize();
  #endif
  selpar_age1_cRs.allocate("selpar_age1_cRs");
  #ifndef NO_AD_INITIALIZE
  selpar_age1_cRs.initialize();
  #endif
  selpar_age2_cRs.allocate("selpar_age2_cRs");
  #ifndef NO_AD_INITIALIZE
  selpar_age2_cRs.initialize();
  #endif
  selpar_age3_cRs.allocate("selpar_age3_cRs");
  #ifndef NO_AD_INITIALIZE
  selpar_age3_cRs.initialize();
  #endif
  selpar_age4_cRs.allocate("selpar_age4_cRs");
  #ifndef NO_AD_INITIALIZE
  selpar_age4_cRs.initialize();
  #endif
  selpar_age5_cRs.allocate("selpar_age5_cRs");
  #ifndef NO_AD_INITIALIZE
  selpar_age5_cRs.initialize();
  #endif
  selpar_age6_cRs.allocate("selpar_age6_cRs");
  #ifndef NO_AD_INITIALIZE
  selpar_age6_cRs.initialize();
  #endif
  selpar_age0_cRs_out.allocate(1,8,"selpar_age0_cRs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age0_cRs_out.initialize();
  #endif
  selpar_age1_cRs_out.allocate(1,8,"selpar_age1_cRs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age1_cRs_out.initialize();
  #endif
  selpar_age2_cRs_out.allocate(1,8,"selpar_age2_cRs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age2_cRs_out.initialize();
  #endif
  selpar_age3_cRs_out.allocate(1,8,"selpar_age3_cRs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age3_cRs_out.initialize();
  #endif
  selpar_age4_cRs_out.allocate(1,8,"selpar_age4_cRs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age4_cRs_out.initialize();
  #endif
  selpar_age5_cRs_out.allocate(1,8,"selpar_age5_cRs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age5_cRs_out.initialize();
  #endif
  selpar_age6_cRs_out.allocate(1,8,"selpar_age6_cRs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age6_cRs_out.initialize();
  #endif
  sel_age0_cRs2_logit.allocate(selpar_age0_cRs2_LO,selpar_age0_cRs2_HI,selpar_age0_cRs2_PH,"sel_age0_cRs2_logit");
  sel_age1_cRs2_logit.allocate(selpar_age1_cRs2_LO,selpar_age1_cRs2_HI,selpar_age1_cRs2_PH,"sel_age1_cRs2_logit");
  sel_age2_cRs2_logit.allocate(selpar_age2_cRs2_LO,selpar_age2_cRs2_HI,selpar_age2_cRs2_PH,"sel_age2_cRs2_logit");
  sel_age3_cRs2_logit.allocate(selpar_age3_cRs2_LO,selpar_age3_cRs2_HI,selpar_age3_cRs2_PH,"sel_age3_cRs2_logit");
  sel_age4_cRs2_logit.allocate(selpar_age4_cRs2_LO,selpar_age4_cRs2_HI,selpar_age4_cRs2_PH,"sel_age4_cRs2_logit");
  sel_age5_cRs2_logit.allocate(selpar_age5_cRs2_LO,selpar_age5_cRs2_HI,selpar_age5_cRs2_PH,"sel_age5_cRs2_logit");
  sel_age6_cRs2_logit.allocate(selpar_age6_cRs2_LO,selpar_age6_cRs2_HI,selpar_age6_cRs2_PH,"sel_age6_cRs2_logit");
  sel_age_cRs2_vec.allocate(1,nages,"sel_age_cRs2_vec");
  #ifndef NO_AD_INITIALIZE
    sel_age_cRs2_vec.initialize();
  #endif
  selpar_age0_cRs2.allocate("selpar_age0_cRs2");
  #ifndef NO_AD_INITIALIZE
  selpar_age0_cRs2.initialize();
  #endif
  selpar_age1_cRs2.allocate("selpar_age1_cRs2");
  #ifndef NO_AD_INITIALIZE
  selpar_age1_cRs2.initialize();
  #endif
  selpar_age2_cRs2.allocate("selpar_age2_cRs2");
  #ifndef NO_AD_INITIALIZE
  selpar_age2_cRs2.initialize();
  #endif
  selpar_age3_cRs2.allocate("selpar_age3_cRs2");
  #ifndef NO_AD_INITIALIZE
  selpar_age3_cRs2.initialize();
  #endif
  selpar_age4_cRs2.allocate("selpar_age4_cRs2");
  #ifndef NO_AD_INITIALIZE
  selpar_age4_cRs2.initialize();
  #endif
  selpar_age5_cRs2.allocate("selpar_age5_cRs2");
  #ifndef NO_AD_INITIALIZE
  selpar_age5_cRs2.initialize();
  #endif
  selpar_age6_cRs2.allocate("selpar_age6_cRs2");
  #ifndef NO_AD_INITIALIZE
  selpar_age6_cRs2.initialize();
  #endif
  selpar_age0_cRs2_out.allocate(1,8,"selpar_age0_cRs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age0_cRs2_out.initialize();
  #endif
  selpar_age1_cRs2_out.allocate(1,8,"selpar_age1_cRs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age1_cRs2_out.initialize();
  #endif
  selpar_age2_cRs2_out.allocate(1,8,"selpar_age2_cRs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age2_cRs2_out.initialize();
  #endif
  selpar_age3_cRs2_out.allocate(1,8,"selpar_age3_cRs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age3_cRs2_out.initialize();
  #endif
  selpar_age4_cRs2_out.allocate(1,8,"selpar_age4_cRs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age4_cRs2_out.initialize();
  #endif
  selpar_age5_cRs2_out.allocate(1,8,"selpar_age5_cRs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age5_cRs2_out.initialize();
  #endif
  selpar_age6_cRs2_out.allocate(1,8,"selpar_age6_cRs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age6_cRs2_out.initialize();
  #endif
  sel_age0_cRs3_logit.allocate(selpar_age0_cRs3_LO,selpar_age0_cRs3_HI,selpar_age0_cRs3_PH,"sel_age0_cRs3_logit");
  sel_age1_cRs3_logit.allocate(selpar_age1_cRs3_LO,selpar_age1_cRs3_HI,selpar_age1_cRs3_PH,"sel_age1_cRs3_logit");
  sel_age2_cRs3_logit.allocate(selpar_age2_cRs3_LO,selpar_age2_cRs3_HI,selpar_age2_cRs3_PH,"sel_age2_cRs3_logit");
  sel_age3_cRs3_logit.allocate(selpar_age3_cRs3_LO,selpar_age3_cRs3_HI,selpar_age3_cRs3_PH,"sel_age3_cRs3_logit");
  sel_age4_cRs3_logit.allocate(selpar_age4_cRs3_LO,selpar_age4_cRs3_HI,selpar_age4_cRs3_PH,"sel_age4_cRs3_logit");
  sel_age5_cRs3_logit.allocate(selpar_age5_cRs3_LO,selpar_age5_cRs3_HI,selpar_age5_cRs3_PH,"sel_age5_cRs3_logit");
  sel_age6_cRs3_logit.allocate(selpar_age6_cRs3_LO,selpar_age6_cRs3_HI,selpar_age6_cRs3_PH,"sel_age6_cRs3_logit");
  sel_age_cRs3_vec.allocate(1,nages,"sel_age_cRs3_vec");
  #ifndef NO_AD_INITIALIZE
    sel_age_cRs3_vec.initialize();
  #endif
  selpar_age0_cRs3.allocate("selpar_age0_cRs3");
  #ifndef NO_AD_INITIALIZE
  selpar_age0_cRs3.initialize();
  #endif
  selpar_age1_cRs3.allocate("selpar_age1_cRs3");
  #ifndef NO_AD_INITIALIZE
  selpar_age1_cRs3.initialize();
  #endif
  selpar_age2_cRs3.allocate("selpar_age2_cRs3");
  #ifndef NO_AD_INITIALIZE
  selpar_age2_cRs3.initialize();
  #endif
  selpar_age3_cRs3.allocate("selpar_age3_cRs3");
  #ifndef NO_AD_INITIALIZE
  selpar_age3_cRs3.initialize();
  #endif
  selpar_age4_cRs3.allocate("selpar_age4_cRs3");
  #ifndef NO_AD_INITIALIZE
  selpar_age4_cRs3.initialize();
  #endif
  selpar_age5_cRs3.allocate("selpar_age5_cRs3");
  #ifndef NO_AD_INITIALIZE
  selpar_age5_cRs3.initialize();
  #endif
  selpar_age6_cRs3.allocate("selpar_age6_cRs3");
  #ifndef NO_AD_INITIALIZE
  selpar_age6_cRs3.initialize();
  #endif
  selpar_age0_cRs3_out.allocate(1,8,"selpar_age0_cRs3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age0_cRs3_out.initialize();
  #endif
  selpar_age1_cRs3_out.allocate(1,8,"selpar_age1_cRs3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age1_cRs3_out.initialize();
  #endif
  selpar_age2_cRs3_out.allocate(1,8,"selpar_age2_cRs3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age2_cRs3_out.initialize();
  #endif
  selpar_age3_cRs3_out.allocate(1,8,"selpar_age3_cRs3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age3_cRs3_out.initialize();
  #endif
  selpar_age4_cRs3_out.allocate(1,8,"selpar_age4_cRs3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age4_cRs3_out.initialize();
  #endif
  selpar_age5_cRs3_out.allocate(1,8,"selpar_age5_cRs3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age5_cRs3_out.initialize();
  #endif
  selpar_age6_cRs3_out.allocate(1,8,"selpar_age6_cRs3_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age6_cRs3_out.initialize();
  #endif
  sel_age0_cRs4_logit.allocate(selpar_age0_cRs4_LO,selpar_age0_cRs4_HI,selpar_age0_cRs4_PH,"sel_age0_cRs4_logit");
  sel_age1_cRs4_logit.allocate(selpar_age1_cRs4_LO,selpar_age1_cRs4_HI,selpar_age1_cRs4_PH,"sel_age1_cRs4_logit");
  sel_age2_cRs4_logit.allocate(selpar_age2_cRs4_LO,selpar_age2_cRs4_HI,selpar_age2_cRs4_PH,"sel_age2_cRs4_logit");
  sel_age3_cRs4_logit.allocate(selpar_age3_cRs4_LO,selpar_age3_cRs4_HI,selpar_age3_cRs4_PH,"sel_age3_cRs4_logit");
  sel_age4_cRs4_logit.allocate(selpar_age4_cRs4_LO,selpar_age4_cRs4_HI,selpar_age4_cRs4_PH,"sel_age4_cRs4_logit");
  sel_age5_cRs4_logit.allocate(selpar_age5_cRs4_LO,selpar_age5_cRs4_HI,selpar_age5_cRs4_PH,"sel_age5_cRs4_logit");
  sel_age6_cRs4_logit.allocate(selpar_age6_cRs4_LO,selpar_age6_cRs4_HI,selpar_age6_cRs4_PH,"sel_age6_cRs4_logit");
  sel_age_cRs4_vec.allocate(1,nages,"sel_age_cRs4_vec");
  #ifndef NO_AD_INITIALIZE
    sel_age_cRs4_vec.initialize();
  #endif
  selpar_age0_cRs4.allocate("selpar_age0_cRs4");
  #ifndef NO_AD_INITIALIZE
  selpar_age0_cRs4.initialize();
  #endif
  selpar_age1_cRs4.allocate("selpar_age1_cRs4");
  #ifndef NO_AD_INITIALIZE
  selpar_age1_cRs4.initialize();
  #endif
  selpar_age2_cRs4.allocate("selpar_age2_cRs4");
  #ifndef NO_AD_INITIALIZE
  selpar_age2_cRs4.initialize();
  #endif
  selpar_age3_cRs4.allocate("selpar_age3_cRs4");
  #ifndef NO_AD_INITIALIZE
  selpar_age3_cRs4.initialize();
  #endif
  selpar_age4_cRs4.allocate("selpar_age4_cRs4");
  #ifndef NO_AD_INITIALIZE
  selpar_age4_cRs4.initialize();
  #endif
  selpar_age5_cRs4.allocate("selpar_age5_cRs4");
  #ifndef NO_AD_INITIALIZE
  selpar_age5_cRs4.initialize();
  #endif
  selpar_age6_cRs4.allocate("selpar_age6_cRs4");
  #ifndef NO_AD_INITIALIZE
  selpar_age6_cRs4.initialize();
  #endif
  selpar_age0_cRs4_out.allocate(1,8,"selpar_age0_cRs4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age0_cRs4_out.initialize();
  #endif
  selpar_age1_cRs4_out.allocate(1,8,"selpar_age1_cRs4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age1_cRs4_out.initialize();
  #endif
  selpar_age2_cRs4_out.allocate(1,8,"selpar_age2_cRs4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age2_cRs4_out.initialize();
  #endif
  selpar_age3_cRs4_out.allocate(1,8,"selpar_age3_cRs4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age3_cRs4_out.initialize();
  #endif
  selpar_age4_cRs4_out.allocate(1,8,"selpar_age4_cRs4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age4_cRs4_out.initialize();
  #endif
  selpar_age5_cRs4_out.allocate(1,8,"selpar_age5_cRs4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age5_cRs4_out.initialize();
  #endif
  selpar_age6_cRs4_out.allocate(1,8,"selpar_age6_cRs4_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age6_cRs4_out.initialize();
  #endif
  sel_cBn.allocate(styr,endyr,1,nages,"sel_cBn");
  #ifndef NO_AD_INITIALIZE
    sel_cBn.initialize();
  #endif
  sel_cBs.allocate(styr,endyr,1,nages,"sel_cBs");
  #ifndef NO_AD_INITIALIZE
    sel_cBs.initialize();
  #endif
  sel_cBn_block1.allocate(1,nages,"sel_cBn_block1");
  #ifndef NO_AD_INITIALIZE
    sel_cBn_block1.initialize();
  #endif
  sel_cBn_block2.allocate(1,nages,"sel_cBn_block2");
  #ifndef NO_AD_INITIALIZE
    sel_cBn_block2.initialize();
  #endif
  sel_cBs_block1.allocate(1,nages,"sel_cBs_block1");
  #ifndef NO_AD_INITIALIZE
    sel_cBs_block1.initialize();
  #endif
  sel_cBs_block2.allocate(1,nages,"sel_cBs_block2");
  #ifndef NO_AD_INITIALIZE
    sel_cBs_block2.initialize();
  #endif
  selpar_A50_cBn.allocate(selpar_A50_cBn_LO,selpar_A50_cBn_HI,selpar_A50_cBn_PH,"selpar_A50_cBn");
  selpar_slope_cBn.allocate(selpar_slope_cBn_LO,selpar_slope_cBn_HI,selpar_slope_cBn_PH,"selpar_slope_cBn");
  selpar_A502_cBn.allocate(selpar_A502_cBn_LO,selpar_A502_cBn_HI,selpar_A502_cBn_PH,"selpar_A502_cBn");
  selpar_slope2_cBn.allocate(selpar_slope2_cBn_LO,selpar_slope2_cBn_HI,selpar_slope2_cBn_PH,"selpar_slope2_cBn");
  selpar_A50_cBn_out.allocate(1,8,"selpar_A50_cBn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cBn_out.initialize();
  #endif
  selpar_slope_cBn_out.allocate(1,8,"selpar_slope_cBn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cBn_out.initialize();
  #endif
  selpar_A502_cBn_out.allocate(1,8,"selpar_A502_cBn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_cBn_out.initialize();
  #endif
  selpar_slope2_cBn_out.allocate(1,8,"selpar_slope2_cBn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_cBn_out.initialize();
  #endif
  selpar_A50_cBn2.allocate(selpar_A50_cBn2_LO,selpar_A50_cBn2_HI,selpar_A50_cBn2_PH,"selpar_A50_cBn2");
  selpar_slope_cBn2.allocate(selpar_slope_cBn2_LO,selpar_slope_cBn2_HI,selpar_slope_cBn2_PH,"selpar_slope_cBn2");
  selpar_A502_cBn2.allocate(selpar_A502_cBn2_LO,selpar_A502_cBn2_HI,selpar_A502_cBn2_PH,"selpar_A502_cBn2");
  selpar_slope2_cBn2.allocate(selpar_slope2_cBn2_LO,selpar_slope2_cBn2_HI,selpar_slope2_cBn2_PH,"selpar_slope2_cBn2");
  selpar_A50_cBn2_out.allocate(1,8,"selpar_A50_cBn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cBn2_out.initialize();
  #endif
  selpar_slope_cBn2_out.allocate(1,8,"selpar_slope_cBn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cBn2_out.initialize();
  #endif
  selpar_A502_cBn2_out.allocate(1,8,"selpar_A502_cBn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_cBn2_out.initialize();
  #endif
  selpar_slope2_cBn2_out.allocate(1,8,"selpar_slope2_cBn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_cBn2_out.initialize();
  #endif
  selpar_A50_cBs.allocate(selpar_A50_cBs_LO,selpar_A50_cBs_HI,selpar_A50_cBs_PH,"selpar_A50_cBs");
  selpar_slope_cBs.allocate(selpar_slope_cBs_LO,selpar_slope_cBs_HI,selpar_slope_cBs_PH,"selpar_slope_cBs");
  selpar_A502_cBs.allocate(selpar_A502_cBs_LO,selpar_A502_cBs_HI,selpar_A502_cBs_PH,"selpar_A502_cBs");
  selpar_slope2_cBs.allocate(selpar_slope2_cBs_LO,selpar_slope2_cBs_HI,selpar_slope2_cBs_PH,"selpar_slope2_cBs");
  selpar_A50_cBs_out.allocate(1,8,"selpar_A50_cBs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cBs_out.initialize();
  #endif
  selpar_slope_cBs_out.allocate(1,8,"selpar_slope_cBs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cBs_out.initialize();
  #endif
  selpar_A502_cBs_out.allocate(1,8,"selpar_A502_cBs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_cBs_out.initialize();
  #endif
  selpar_slope2_cBs_out.allocate(1,8,"selpar_slope2_cBs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_cBs_out.initialize();
  #endif
  selpar_A50_cBs2.allocate(selpar_A50_cBs2_LO,selpar_A50_cBs2_HI,selpar_A50_cBs2_PH,"selpar_A50_cBs2");
  selpar_slope_cBs2.allocate(selpar_slope_cBs2_LO,selpar_slope_cBs2_HI,selpar_slope_cBs2_PH,"selpar_slope_cBs2");
  selpar_A502_cBs2.allocate(selpar_A502_cBs2_LO,selpar_A502_cBs2_HI,selpar_A502_cBs2_PH,"selpar_A502_cBs2");
  selpar_slope2_cBs2.allocate(selpar_slope2_cBs2_LO,selpar_slope2_cBs2_HI,selpar_slope2_cBs2_PH,"selpar_slope2_cBs2");
  selpar_A50_cBs2_out.allocate(1,8,"selpar_A50_cBs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_cBs2_out.initialize();
  #endif
  selpar_slope_cBs2_out.allocate(1,8,"selpar_slope_cBs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_cBs2_out.initialize();
  #endif
  selpar_A502_cBs2_out.allocate(1,8,"selpar_A502_cBs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_cBs2_out.initialize();
  #endif
  selpar_slope2_cBs2_out.allocate(1,8,"selpar_slope2_cBs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_cBs2_out.initialize();
  #endif
  sel_age0_cBn_logit.allocate(selpar_age0_cBn_LO,selpar_age0_cBn_HI,selpar_age0_cBn_PH,"sel_age0_cBn_logit");
  sel_age1_cBn_logit.allocate(selpar_age1_cBn_LO,selpar_age1_cBn_HI,selpar_age1_cBn_PH,"sel_age1_cBn_logit");
  sel_age2_cBn_logit.allocate(selpar_age2_cBn_LO,selpar_age2_cBn_HI,selpar_age2_cBn_PH,"sel_age2_cBn_logit");
  sel_age3_cBn_logit.allocate(selpar_age3_cBn_LO,selpar_age3_cBn_HI,selpar_age3_cBn_PH,"sel_age3_cBn_logit");
  sel_age4_cBn_logit.allocate(selpar_age4_cBn_LO,selpar_age4_cBn_HI,selpar_age4_cBn_PH,"sel_age4_cBn_logit");
  sel_age5_cBn_logit.allocate(selpar_age5_cBn_LO,selpar_age5_cBn_HI,selpar_age5_cBn_PH,"sel_age5_cBn_logit");
  sel_age6_cBn_logit.allocate(selpar_age6_cBn_LO,selpar_age6_cBn_HI,selpar_age6_cBn_PH,"sel_age6_cBn_logit");
  sel_age_cBn_vec.allocate(1,nages,"sel_age_cBn_vec");
  #ifndef NO_AD_INITIALIZE
    sel_age_cBn_vec.initialize();
  #endif
  selpar_age0_cBn.allocate("selpar_age0_cBn");
  #ifndef NO_AD_INITIALIZE
  selpar_age0_cBn.initialize();
  #endif
  selpar_age1_cBn.allocate("selpar_age1_cBn");
  #ifndef NO_AD_INITIALIZE
  selpar_age1_cBn.initialize();
  #endif
  selpar_age2_cBn.allocate("selpar_age2_cBn");
  #ifndef NO_AD_INITIALIZE
  selpar_age2_cBn.initialize();
  #endif
  selpar_age3_cBn.allocate("selpar_age3_cBn");
  #ifndef NO_AD_INITIALIZE
  selpar_age3_cBn.initialize();
  #endif
  selpar_age4_cBn.allocate("selpar_age4_cBn");
  #ifndef NO_AD_INITIALIZE
  selpar_age4_cBn.initialize();
  #endif
  selpar_age5_cBn.allocate("selpar_age5_cBn");
  #ifndef NO_AD_INITIALIZE
  selpar_age5_cBn.initialize();
  #endif
  selpar_age6_cBn.allocate("selpar_age6_cBn");
  #ifndef NO_AD_INITIALIZE
  selpar_age6_cBn.initialize();
  #endif
  selpar_age0_cBn_out.allocate(1,8,"selpar_age0_cBn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age0_cBn_out.initialize();
  #endif
  selpar_age1_cBn_out.allocate(1,8,"selpar_age1_cBn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age1_cBn_out.initialize();
  #endif
  selpar_age2_cBn_out.allocate(1,8,"selpar_age2_cBn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age2_cBn_out.initialize();
  #endif
  selpar_age3_cBn_out.allocate(1,8,"selpar_age3_cBn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age3_cBn_out.initialize();
  #endif
  selpar_age4_cBn_out.allocate(1,8,"selpar_age4_cBn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age4_cBn_out.initialize();
  #endif
  selpar_age5_cBn_out.allocate(1,8,"selpar_age5_cBn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age5_cBn_out.initialize();
  #endif
  selpar_age6_cBn_out.allocate(1,8,"selpar_age6_cBn_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age6_cBn_out.initialize();
  #endif
  sel_age0_cBn2_logit.allocate(selpar_age0_cBn2_LO,selpar_age0_cBn2_HI,selpar_age0_cBn2_PH,"sel_age0_cBn2_logit");
  sel_age1_cBn2_logit.allocate(selpar_age1_cBn2_LO,selpar_age1_cBn2_HI,selpar_age1_cBn2_PH,"sel_age1_cBn2_logit");
  sel_age2_cBn2_logit.allocate(selpar_age2_cBn2_LO,selpar_age2_cBn2_HI,selpar_age2_cBn2_PH,"sel_age2_cBn2_logit");
  sel_age3_cBn2_logit.allocate(selpar_age3_cBn2_LO,selpar_age3_cBn2_HI,selpar_age3_cBn2_PH,"sel_age3_cBn2_logit");
  sel_age4_cBn2_logit.allocate(selpar_age4_cBn2_LO,selpar_age4_cBn2_HI,selpar_age4_cBn2_PH,"sel_age4_cBn2_logit");
  sel_age5_cBn2_logit.allocate(selpar_age5_cBn2_LO,selpar_age5_cBn2_HI,selpar_age5_cBn2_PH,"sel_age5_cBn2_logit");
  sel_age6_cBn2_logit.allocate(selpar_age6_cBn2_LO,selpar_age6_cBn2_HI,selpar_age6_cBn2_PH,"sel_age6_cBn2_logit");
  sel_age_cBn2_vec.allocate(1,nages,"sel_age_cBn2_vec");
  #ifndef NO_AD_INITIALIZE
    sel_age_cBn2_vec.initialize();
  #endif
  selpar_age0_cBn2.allocate("selpar_age0_cBn2");
  #ifndef NO_AD_INITIALIZE
  selpar_age0_cBn2.initialize();
  #endif
  selpar_age1_cBn2.allocate("selpar_age1_cBn2");
  #ifndef NO_AD_INITIALIZE
  selpar_age1_cBn2.initialize();
  #endif
  selpar_age2_cBn2.allocate("selpar_age2_cBn2");
  #ifndef NO_AD_INITIALIZE
  selpar_age2_cBn2.initialize();
  #endif
  selpar_age3_cBn2.allocate("selpar_age3_cBn2");
  #ifndef NO_AD_INITIALIZE
  selpar_age3_cBn2.initialize();
  #endif
  selpar_age4_cBn2.allocate("selpar_age4_cBn2");
  #ifndef NO_AD_INITIALIZE
  selpar_age4_cBn2.initialize();
  #endif
  selpar_age5_cBn2.allocate("selpar_age5_cBn2");
  #ifndef NO_AD_INITIALIZE
  selpar_age5_cBn2.initialize();
  #endif
  selpar_age6_cBn2.allocate("selpar_age6_cBn2");
  #ifndef NO_AD_INITIALIZE
  selpar_age6_cBn2.initialize();
  #endif
  selpar_age0_cBn2_out.allocate(1,8,"selpar_age0_cBn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age0_cBn2_out.initialize();
  #endif
  selpar_age1_cBn2_out.allocate(1,8,"selpar_age1_cBn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age1_cBn2_out.initialize();
  #endif
  selpar_age2_cBn2_out.allocate(1,8,"selpar_age2_cBn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age2_cBn2_out.initialize();
  #endif
  selpar_age3_cBn2_out.allocate(1,8,"selpar_age3_cBn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age3_cBn2_out.initialize();
  #endif
  selpar_age4_cBn2_out.allocate(1,8,"selpar_age4_cBn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age4_cBn2_out.initialize();
  #endif
  selpar_age5_cBn2_out.allocate(1,8,"selpar_age5_cBn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age5_cBn2_out.initialize();
  #endif
  selpar_age6_cBn2_out.allocate(1,8,"selpar_age6_cBn2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age6_cBn2_out.initialize();
  #endif
  sel_age0_cBs_logit.allocate(selpar_age0_cBs_LO,selpar_age0_cBs_HI,selpar_age0_cBs_PH,"sel_age0_cBs_logit");
  sel_age1_cBs_logit.allocate(selpar_age1_cBs_LO,selpar_age1_cBs_HI,selpar_age1_cBs_PH,"sel_age1_cBs_logit");
  sel_age2_cBs_logit.allocate(selpar_age2_cBs_LO,selpar_age2_cBs_HI,selpar_age2_cBs_PH,"sel_age2_cBs_logit");
  sel_age3_cBs_logit.allocate(selpar_age3_cBs_LO,selpar_age3_cBs_HI,selpar_age3_cBs_PH,"sel_age3_cBs_logit");
  sel_age4_cBs_logit.allocate(selpar_age4_cBs_LO,selpar_age4_cBs_HI,selpar_age4_cBs_PH,"sel_age4_cBs_logit");
  sel_age5_cBs_logit.allocate(selpar_age5_cBs_LO,selpar_age5_cBs_HI,selpar_age5_cBs_PH,"sel_age5_cBs_logit");
  sel_age6_cBs_logit.allocate(selpar_age6_cBs_LO,selpar_age6_cBs_HI,selpar_age6_cBs_PH,"sel_age6_cBs_logit");
  sel_age_cBs_vec.allocate(1,nages,"sel_age_cBs_vec");
  #ifndef NO_AD_INITIALIZE
    sel_age_cBs_vec.initialize();
  #endif
  selpar_age0_cBs.allocate("selpar_age0_cBs");
  #ifndef NO_AD_INITIALIZE
  selpar_age0_cBs.initialize();
  #endif
  selpar_age1_cBs.allocate("selpar_age1_cBs");
  #ifndef NO_AD_INITIALIZE
  selpar_age1_cBs.initialize();
  #endif
  selpar_age2_cBs.allocate("selpar_age2_cBs");
  #ifndef NO_AD_INITIALIZE
  selpar_age2_cBs.initialize();
  #endif
  selpar_age3_cBs.allocate("selpar_age3_cBs");
  #ifndef NO_AD_INITIALIZE
  selpar_age3_cBs.initialize();
  #endif
  selpar_age4_cBs.allocate("selpar_age4_cBs");
  #ifndef NO_AD_INITIALIZE
  selpar_age4_cBs.initialize();
  #endif
  selpar_age5_cBs.allocate("selpar_age5_cBs");
  #ifndef NO_AD_INITIALIZE
  selpar_age5_cBs.initialize();
  #endif
  selpar_age6_cBs.allocate("selpar_age6_cBs");
  #ifndef NO_AD_INITIALIZE
  selpar_age6_cBs.initialize();
  #endif
  selpar_age0_cBs_out.allocate(1,8,"selpar_age0_cBs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age0_cBs_out.initialize();
  #endif
  selpar_age1_cBs_out.allocate(1,8,"selpar_age1_cBs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age1_cBs_out.initialize();
  #endif
  selpar_age2_cBs_out.allocate(1,8,"selpar_age2_cBs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age2_cBs_out.initialize();
  #endif
  selpar_age3_cBs_out.allocate(1,8,"selpar_age3_cBs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age3_cBs_out.initialize();
  #endif
  selpar_age4_cBs_out.allocate(1,8,"selpar_age4_cBs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age4_cBs_out.initialize();
  #endif
  selpar_age5_cBs_out.allocate(1,8,"selpar_age5_cBs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age5_cBs_out.initialize();
  #endif
  selpar_age6_cBs_out.allocate(1,8,"selpar_age6_cBs_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age6_cBs_out.initialize();
  #endif
  sel_age0_cBs2_logit.allocate(selpar_age0_cBs2_LO,selpar_age0_cBs2_HI,selpar_age0_cBs2_PH,"sel_age0_cBs2_logit");
  sel_age1_cBs2_logit.allocate(selpar_age1_cBs2_LO,selpar_age1_cBs2_HI,selpar_age1_cBs2_PH,"sel_age1_cBs2_logit");
  sel_age2_cBs2_logit.allocate(selpar_age2_cBs2_LO,selpar_age2_cBs2_HI,selpar_age2_cBs2_PH,"sel_age2_cBs2_logit");
  sel_age3_cBs2_logit.allocate(selpar_age3_cBs2_LO,selpar_age3_cBs2_HI,selpar_age3_cBs2_PH,"sel_age3_cBs2_logit");
  sel_age4_cBs2_logit.allocate(selpar_age4_cBs2_LO,selpar_age4_cBs2_HI,selpar_age4_cBs2_PH,"sel_age4_cBs2_logit");
  sel_age5_cBs2_logit.allocate(selpar_age5_cBs2_LO,selpar_age5_cBs2_HI,selpar_age5_cBs2_PH,"sel_age5_cBs2_logit");
  sel_age6_cBs2_logit.allocate(selpar_age6_cBs2_LO,selpar_age6_cBs2_HI,selpar_age6_cBs2_PH,"sel_age6_cBs2_logit");
  sel_age_cBs2_vec.allocate(1,nages,"sel_age_cBs2_vec");
  #ifndef NO_AD_INITIALIZE
    sel_age_cBs2_vec.initialize();
  #endif
  selpar_age0_cBs2.allocate("selpar_age0_cBs2");
  #ifndef NO_AD_INITIALIZE
  selpar_age0_cBs2.initialize();
  #endif
  selpar_age1_cBs2.allocate("selpar_age1_cBs2");
  #ifndef NO_AD_INITIALIZE
  selpar_age1_cBs2.initialize();
  #endif
  selpar_age2_cBs2.allocate("selpar_age2_cBs2");
  #ifndef NO_AD_INITIALIZE
  selpar_age2_cBs2.initialize();
  #endif
  selpar_age3_cBs2.allocate("selpar_age3_cBs2");
  #ifndef NO_AD_INITIALIZE
  selpar_age3_cBs2.initialize();
  #endif
  selpar_age4_cBs2.allocate("selpar_age4_cBs2");
  #ifndef NO_AD_INITIALIZE
  selpar_age4_cBs2.initialize();
  #endif
  selpar_age5_cBs2.allocate("selpar_age5_cBs2");
  #ifndef NO_AD_INITIALIZE
  selpar_age5_cBs2.initialize();
  #endif
  selpar_age6_cBs2.allocate("selpar_age6_cBs2");
  #ifndef NO_AD_INITIALIZE
  selpar_age6_cBs2.initialize();
  #endif
  selpar_age0_cBs2_out.allocate(1,8,"selpar_age0_cBs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age0_cBs2_out.initialize();
  #endif
  selpar_age1_cBs2_out.allocate(1,8,"selpar_age1_cBs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age1_cBs2_out.initialize();
  #endif
  selpar_age2_cBs2_out.allocate(1,8,"selpar_age2_cBs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age2_cBs2_out.initialize();
  #endif
  selpar_age3_cBs2_out.allocate(1,8,"selpar_age3_cBs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age3_cBs2_out.initialize();
  #endif
  selpar_age4_cBs2_out.allocate(1,8,"selpar_age4_cBs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age4_cBs2_out.initialize();
  #endif
  selpar_age5_cBs2_out.allocate(1,8,"selpar_age5_cBs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age5_cBs2_out.initialize();
  #endif
  selpar_age6_cBs2_out.allocate(1,8,"selpar_age6_cBs2_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age6_cBs2_out.initialize();
  #endif
  sel_nad.allocate(styr_cpue_nad,endyr_cpue_nad,1,nages,"sel_nad");
  #ifndef NO_AD_INITIALIZE
    sel_nad.initialize();
  #endif
  sel_nad_block1.allocate(1,nages,"sel_nad_block1");
  #ifndef NO_AD_INITIALIZE
    sel_nad_block1.initialize();
  #endif
  selpar_A50_nad.allocate(selpar_A50_nad_LO,selpar_A50_nad_HI,selpar_A50_nad_PH,"selpar_A50_nad");
  selpar_slope_nad.allocate(selpar_slope_nad_LO,selpar_slope_nad_HI,selpar_slope_nad_PH,"selpar_slope_nad");
  selpar_A502_nad.allocate(selpar_A502_nad_LO,selpar_A502_nad_HI,selpar_A502_nad_PH,"selpar_A502_nad");
  selpar_slope2_nad.allocate(selpar_slope2_nad_LO,selpar_slope2_nad_HI,selpar_slope2_nad_PH,"selpar_slope2_nad");
  selpar_A50_nad_out.allocate(1,8,"selpar_A50_nad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_nad_out.initialize();
  #endif
  selpar_slope_nad_out.allocate(1,8,"selpar_slope_nad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_nad_out.initialize();
  #endif
  selpar_A502_nad_out.allocate(1,8,"selpar_A502_nad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_nad_out.initialize();
  #endif
  selpar_slope2_nad_out.allocate(1,8,"selpar_slope2_nad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_nad_out.initialize();
  #endif
  sel_age0_nad_logit.allocate(selpar_age0_nad_LO,selpar_age0_nad_HI,selpar_age0_nad_PH,"sel_age0_nad_logit");
  sel_age1_nad_logit.allocate(selpar_age1_nad_LO,selpar_age1_nad_HI,selpar_age1_nad_PH,"sel_age1_nad_logit");
  sel_age2_nad_logit.allocate(selpar_age2_nad_LO,selpar_age2_nad_HI,selpar_age2_nad_PH,"sel_age2_nad_logit");
  sel_age3_nad_logit.allocate(selpar_age3_nad_LO,selpar_age3_nad_HI,selpar_age3_nad_PH,"sel_age3_nad_logit");
  sel_age4_nad_logit.allocate(selpar_age4_nad_LO,selpar_age4_nad_HI,selpar_age4_nad_PH,"sel_age4_nad_logit");
  sel_age5_nad_logit.allocate(selpar_age5_nad_LO,selpar_age5_nad_HI,selpar_age5_nad_PH,"sel_age5_nad_logit");
  sel_age6_nad_logit.allocate(selpar_age6_nad_LO,selpar_age6_nad_HI,selpar_age6_nad_PH,"sel_age6_nad_logit");
  sel_age_nad_vec.allocate(1,nages,"sel_age_nad_vec");
  #ifndef NO_AD_INITIALIZE
    sel_age_nad_vec.initialize();
  #endif
  selpar_age0_nad.allocate("selpar_age0_nad");
  #ifndef NO_AD_INITIALIZE
  selpar_age0_nad.initialize();
  #endif
  selpar_age1_nad.allocate("selpar_age1_nad");
  #ifndef NO_AD_INITIALIZE
  selpar_age1_nad.initialize();
  #endif
  selpar_age2_nad.allocate("selpar_age2_nad");
  #ifndef NO_AD_INITIALIZE
  selpar_age2_nad.initialize();
  #endif
  selpar_age3_nad.allocate("selpar_age3_nad");
  #ifndef NO_AD_INITIALIZE
  selpar_age3_nad.initialize();
  #endif
  selpar_age4_nad.allocate("selpar_age4_nad");
  #ifndef NO_AD_INITIALIZE
  selpar_age4_nad.initialize();
  #endif
  selpar_age5_nad.allocate("selpar_age5_nad");
  #ifndef NO_AD_INITIALIZE
  selpar_age5_nad.initialize();
  #endif
  selpar_age6_nad.allocate("selpar_age6_nad");
  #ifndef NO_AD_INITIALIZE
  selpar_age6_nad.initialize();
  #endif
  selpar_age0_nad_out.allocate(1,8,"selpar_age0_nad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age0_nad_out.initialize();
  #endif
  selpar_age1_nad_out.allocate(1,8,"selpar_age1_nad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age1_nad_out.initialize();
  #endif
  selpar_age2_nad_out.allocate(1,8,"selpar_age2_nad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age2_nad_out.initialize();
  #endif
  selpar_age3_nad_out.allocate(1,8,"selpar_age3_nad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age3_nad_out.initialize();
  #endif
  selpar_age4_nad_out.allocate(1,8,"selpar_age4_nad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age4_nad_out.initialize();
  #endif
  selpar_age5_nad_out.allocate(1,8,"selpar_age5_nad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age5_nad_out.initialize();
  #endif
  selpar_age6_nad_out.allocate(1,8,"selpar_age6_nad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age6_nad_out.initialize();
  #endif
  sel_mad.allocate(styr_cpue_mad,endyr_cpue_mad,1,nages,"sel_mad");
  #ifndef NO_AD_INITIALIZE
    sel_mad.initialize();
  #endif
  sel_mad_block1.allocate(1,nages,"sel_mad_block1");
  #ifndef NO_AD_INITIALIZE
    sel_mad_block1.initialize();
  #endif
  selpar_A50_mad.allocate(selpar_A50_mad_LO,selpar_A50_mad_HI,selpar_A50_mad_PH,"selpar_A50_mad");
  selpar_slope_mad.allocate(selpar_slope_mad_LO,selpar_slope_mad_HI,selpar_slope_mad_PH,"selpar_slope_mad");
  selpar_A502_mad.allocate(selpar_A502_mad_LO,selpar_A502_mad_HI,selpar_A502_mad_PH,"selpar_A502_mad");
  selpar_slope2_mad.allocate(selpar_slope2_mad_LO,selpar_slope2_mad_HI,selpar_slope2_mad_PH,"selpar_slope2_mad");
  selpar_A50_mad_out.allocate(1,8,"selpar_A50_mad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_mad_out.initialize();
  #endif
  selpar_slope_mad_out.allocate(1,8,"selpar_slope_mad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_mad_out.initialize();
  #endif
  selpar_A502_mad_out.allocate(1,8,"selpar_A502_mad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_mad_out.initialize();
  #endif
  selpar_slope2_mad_out.allocate(1,8,"selpar_slope2_mad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_mad_out.initialize();
  #endif
  sel_age0_mad_logit.allocate(selpar_age0_mad_LO,selpar_age0_mad_HI,selpar_age0_mad_PH,"sel_age0_mad_logit");
  sel_age1_mad_logit.allocate(selpar_age1_mad_LO,selpar_age1_mad_HI,selpar_age1_mad_PH,"sel_age1_mad_logit");
  sel_age2_mad_logit.allocate(selpar_age2_mad_LO,selpar_age2_mad_HI,selpar_age2_mad_PH,"sel_age2_mad_logit");
  sel_age3_mad_logit.allocate(selpar_age3_mad_LO,selpar_age3_mad_HI,selpar_age3_mad_PH,"sel_age3_mad_logit");
  sel_age4_mad_logit.allocate(selpar_age4_mad_LO,selpar_age4_mad_HI,selpar_age4_mad_PH,"sel_age4_mad_logit");
  sel_age5_mad_logit.allocate(selpar_age5_mad_LO,selpar_age5_mad_HI,selpar_age5_mad_PH,"sel_age5_mad_logit");
  sel_age6_mad_logit.allocate(selpar_age6_mad_LO,selpar_age6_mad_HI,selpar_age6_mad_PH,"sel_age6_mad_logit");
  sel_age_mad_vec.allocate(1,nages,"sel_age_mad_vec");
  #ifndef NO_AD_INITIALIZE
    sel_age_mad_vec.initialize();
  #endif
  selpar_age0_mad.allocate("selpar_age0_mad");
  #ifndef NO_AD_INITIALIZE
  selpar_age0_mad.initialize();
  #endif
  selpar_age1_mad.allocate("selpar_age1_mad");
  #ifndef NO_AD_INITIALIZE
  selpar_age1_mad.initialize();
  #endif
  selpar_age2_mad.allocate("selpar_age2_mad");
  #ifndef NO_AD_INITIALIZE
  selpar_age2_mad.initialize();
  #endif
  selpar_age3_mad.allocate("selpar_age3_mad");
  #ifndef NO_AD_INITIALIZE
  selpar_age3_mad.initialize();
  #endif
  selpar_age4_mad.allocate("selpar_age4_mad");
  #ifndef NO_AD_INITIALIZE
  selpar_age4_mad.initialize();
  #endif
  selpar_age5_mad.allocate("selpar_age5_mad");
  #ifndef NO_AD_INITIALIZE
  selpar_age5_mad.initialize();
  #endif
  selpar_age6_mad.allocate("selpar_age6_mad");
  #ifndef NO_AD_INITIALIZE
  selpar_age6_mad.initialize();
  #endif
  selpar_age0_mad_out.allocate(1,8,"selpar_age0_mad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age0_mad_out.initialize();
  #endif
  selpar_age1_mad_out.allocate(1,8,"selpar_age1_mad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age1_mad_out.initialize();
  #endif
  selpar_age2_mad_out.allocate(1,8,"selpar_age2_mad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age2_mad_out.initialize();
  #endif
  selpar_age3_mad_out.allocate(1,8,"selpar_age3_mad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age3_mad_out.initialize();
  #endif
  selpar_age4_mad_out.allocate(1,8,"selpar_age4_mad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age4_mad_out.initialize();
  #endif
  selpar_age5_mad_out.allocate(1,8,"selpar_age5_mad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age5_mad_out.initialize();
  #endif
  selpar_age6_mad_out.allocate(1,8,"selpar_age6_mad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age6_mad_out.initialize();
  #endif
  sel_sad.allocate(styr_cpue_sad,endyr_cpue_sad,1,nages,"sel_sad");
  #ifndef NO_AD_INITIALIZE
    sel_sad.initialize();
  #endif
  sel_sad_block1.allocate(1,nages,"sel_sad_block1");
  #ifndef NO_AD_INITIALIZE
    sel_sad_block1.initialize();
  #endif
  selpar_A50_sad.allocate(selpar_A50_sad_LO,selpar_A50_sad_HI,selpar_A50_sad_PH,"selpar_A50_sad");
  selpar_slope_sad.allocate(selpar_slope_sad_LO,selpar_slope_sad_HI,selpar_slope_sad_PH,"selpar_slope_sad");
  selpar_A502_sad.allocate(selpar_A502_sad_LO,selpar_A502_sad_HI,selpar_A502_sad_PH,"selpar_A502_sad");
  selpar_slope2_sad.allocate(selpar_slope2_sad_LO,selpar_slope2_sad_HI,selpar_slope2_sad_PH,"selpar_slope2_sad");
  selpar_A50_sad_out.allocate(1,8,"selpar_A50_sad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A50_sad_out.initialize();
  #endif
  selpar_slope_sad_out.allocate(1,8,"selpar_slope_sad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope_sad_out.initialize();
  #endif
  selpar_A502_sad_out.allocate(1,8,"selpar_A502_sad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_A502_sad_out.initialize();
  #endif
  selpar_slope2_sad_out.allocate(1,8,"selpar_slope2_sad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_slope2_sad_out.initialize();
  #endif
  sel_age0_sad_logit.allocate(selpar_age0_sad_LO,selpar_age0_sad_HI,selpar_age0_sad_PH,"sel_age0_sad_logit");
  sel_age1_sad_logit.allocate(selpar_age1_sad_LO,selpar_age1_sad_HI,selpar_age1_sad_PH,"sel_age1_sad_logit");
  sel_age2_sad_logit.allocate(selpar_age2_sad_LO,selpar_age2_sad_HI,selpar_age2_sad_PH,"sel_age2_sad_logit");
  sel_age3_sad_logit.allocate(selpar_age3_sad_LO,selpar_age3_sad_HI,selpar_age3_sad_PH,"sel_age3_sad_logit");
  sel_age4_sad_logit.allocate(selpar_age4_sad_LO,selpar_age4_sad_HI,selpar_age4_sad_PH,"sel_age4_sad_logit");
  sel_age5_sad_logit.allocate(selpar_age5_sad_LO,selpar_age5_sad_HI,selpar_age5_sad_PH,"sel_age5_sad_logit");
  sel_age6_sad_logit.allocate(selpar_age6_sad_LO,selpar_age6_sad_HI,selpar_age6_sad_PH,"sel_age6_sad_logit");
  sel_age_sad_vec.allocate(1,nages,"sel_age_sad_vec");
  #ifndef NO_AD_INITIALIZE
    sel_age_sad_vec.initialize();
  #endif
  selpar_age0_sad.allocate("selpar_age0_sad");
  #ifndef NO_AD_INITIALIZE
  selpar_age0_sad.initialize();
  #endif
  selpar_age1_sad.allocate("selpar_age1_sad");
  #ifndef NO_AD_INITIALIZE
  selpar_age1_sad.initialize();
  #endif
  selpar_age2_sad.allocate("selpar_age2_sad");
  #ifndef NO_AD_INITIALIZE
  selpar_age2_sad.initialize();
  #endif
  selpar_age3_sad.allocate("selpar_age3_sad");
  #ifndef NO_AD_INITIALIZE
  selpar_age3_sad.initialize();
  #endif
  selpar_age4_sad.allocate("selpar_age4_sad");
  #ifndef NO_AD_INITIALIZE
  selpar_age4_sad.initialize();
  #endif
  selpar_age5_sad.allocate("selpar_age5_sad");
  #ifndef NO_AD_INITIALIZE
  selpar_age5_sad.initialize();
  #endif
  selpar_age6_sad.allocate("selpar_age6_sad");
  #ifndef NO_AD_INITIALIZE
  selpar_age6_sad.initialize();
  #endif
  selpar_age0_sad_out.allocate(1,8,"selpar_age0_sad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age0_sad_out.initialize();
  #endif
  selpar_age1_sad_out.allocate(1,8,"selpar_age1_sad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age1_sad_out.initialize();
  #endif
  selpar_age2_sad_out.allocate(1,8,"selpar_age2_sad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age2_sad_out.initialize();
  #endif
  selpar_age3_sad_out.allocate(1,8,"selpar_age3_sad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age3_sad_out.initialize();
  #endif
  selpar_age4_sad_out.allocate(1,8,"selpar_age4_sad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age4_sad_out.initialize();
  #endif
  selpar_age5_sad_out.allocate(1,8,"selpar_age5_sad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age5_sad_out.initialize();
  #endif
  selpar_age6_sad_out.allocate(1,8,"selpar_age6_sad_out");
  #ifndef NO_AD_INITIALIZE
    selpar_age6_sad_out.initialize();
  #endif
  sel_wgted_L.allocate(1,nages,"sel_wgted_L");
  #ifndef NO_AD_INITIALIZE
    sel_wgted_L.initialize();
  #endif
  sel_wgted_tot.allocate(1,nages,"sel_wgted_tot");
  #ifndef NO_AD_INITIALIZE
    sel_wgted_tot.initialize();
  #endif
  pred_nad_cpue.allocate(styr_cpue_nad,endyr_cpue_nad,"pred_nad_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_nad_cpue.initialize();
  #endif
  N_nad_cpue.allocate(styr_cpue_nad,endyr_cpue_nad,1,nages,"N_nad_cpue");
  #ifndef NO_AD_INITIALIZE
    N_nad_cpue.initialize();
  #endif
  pred_mad_cpue.allocate(styr_cpue_mad,endyr_cpue_mad,"pred_mad_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_mad_cpue.initialize();
  #endif
  N_mad_cpue.allocate(styr_cpue_mad,endyr_cpue_mad,1,nages,"N_mad_cpue");
  #ifndef NO_AD_INITIALIZE
    N_mad_cpue.initialize();
  #endif
  pred_sad_cpue.allocate(styr_cpue_sad,endyr_cpue_sad,"pred_sad_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_sad_cpue.initialize();
  #endif
  N_sad_cpue.allocate(styr_cpue_sad,endyr_cpue_sad,"N_sad_cpue");
  #ifndef NO_AD_INITIALIZE
    N_sad_cpue.initialize();
  #endif
  pred_jai_cpue.allocate(styr_cpue_jai,endyr_cpue_jai,"pred_jai_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_jai_cpue.initialize();
  #endif
  N_jai_cpue.allocate(styr_cpue_jai,endyr_cpue_jai,"N_jai_cpue");
  #ifndef NO_AD_INITIALIZE
    N_jai_cpue.initialize();
  #endif
  pred_mareco_cpue.allocate(1,nyr_cpue_mareco,"pred_mareco_cpue");
  #ifndef NO_AD_INITIALIZE
    pred_mareco_cpue.initialize();
  #endif
  SSB_mareco_cpue.allocate(1,nyr_cpue_mareco,"SSB_mareco_cpue");
  #ifndef NO_AD_INITIALIZE
    SSB_mareco_cpue.initialize();
  #endif
  log_q_cpue_nad.allocate(log_q_nad_LO,log_q_nad_HI,log_q_nad_PH,"log_q_cpue_nad");
  log_q_cpue_mad.allocate(log_q_mad_LO,log_q_mad_HI,log_q_mad_PH,"log_q_cpue_mad");
  log_q_cpue_sad.allocate(log_q_sad_LO,log_q_sad_HI,log_q_sad_PH,"log_q_cpue_sad");
  log_q_cpue_jai.allocate(log_q_jai_LO,log_q_jai_HI,log_q_jai_PH,"log_q_cpue_jai");
  log_q2_jai.allocate(log_q2_jai_LO,log_q2_jai_HI,log_q2_jai_PH,"log_q2_jai");
  log_q_cpue_mar.allocate(log_q_mar_LO,log_q_mar_HI,log_q_mar_PH,"log_q_cpue_mar");
  log_q_cpue_eco.allocate(log_q_eco_LO,log_q_eco_HI,log_q_eco_PH,"log_q_cpue_eco");
  log_q_nad_out.allocate(1,8,"log_q_nad_out");
  #ifndef NO_AD_INITIALIZE
    log_q_nad_out.initialize();
  #endif
  log_q_mad_out.allocate(1,8,"log_q_mad_out");
  #ifndef NO_AD_INITIALIZE
    log_q_mad_out.initialize();
  #endif
  log_q_sad_out.allocate(1,8,"log_q_sad_out");
  #ifndef NO_AD_INITIALIZE
    log_q_sad_out.initialize();
  #endif
  log_q_jai_out.allocate(1,8,"log_q_jai_out");
  #ifndef NO_AD_INITIALIZE
    log_q_jai_out.initialize();
  #endif
  log_q2_jai_out.allocate(1,8,"log_q2_jai_out");
  #ifndef NO_AD_INITIALIZE
    log_q2_jai_out.initialize();
  #endif
  log_q_mar_out.allocate(1,8,"log_q_mar_out");
  #ifndef NO_AD_INITIALIZE
    log_q_mar_out.initialize();
  #endif
  log_q_eco_out.allocate(1,8,"log_q_eco_out");
  #ifndef NO_AD_INITIALIZE
    log_q_eco_out.initialize();
  #endif
  q_rate.allocate("q_rate");
  #ifndef NO_AD_INITIALIZE
  q_rate.initialize();
  #endif
  q_rate_fcn_nad.allocate(styr_cpue_nad,endyr_cpue_nad,"q_rate_fcn_nad");
  #ifndef NO_AD_INITIALIZE
    q_rate_fcn_nad.initialize();
  #endif
  q_rate_fcn_mad.allocate(styr_cpue_mad,endyr_cpue_mad,"q_rate_fcn_mad");
  #ifndef NO_AD_INITIALIZE
    q_rate_fcn_mad.initialize();
  #endif
  q_rate_fcn_sad.allocate(styr_cpue_sad,endyr_cpue_sad,"q_rate_fcn_sad");
  #ifndef NO_AD_INITIALIZE
    q_rate_fcn_sad.initialize();
  #endif
  q_rate_fcn_jai.allocate(styr_cpue_jai,endyr_cpue_jai,"q_rate_fcn_jai");
  #ifndef NO_AD_INITIALIZE
    q_rate_fcn_jai.initialize();
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
  q_RW_log_dev_nad.allocate(styr_cpue_nad,endyr_cpue_nad-1,log_RWq_LO,log_RWq_HI,log_RWq_PH,"q_RW_log_dev_nad");
  q_RW_log_dev_mad.allocate(styr_cpue_mad,endyr_cpue_mad-1,log_RWq_LO,log_RWq_HI,log_RWq_PH,"q_RW_log_dev_mad");
  q_RW_log_dev_sad.allocate(styr_cpue_sad,endyr_cpue_sad-1,log_RWq_LO,log_RWq_HI,log_RWq_PH,"q_RW_log_dev_sad");
  q_RW_log_dev_jai.allocate(styr_cpue_jai,endyr_cpue_jai-1,log_RWq_LO,log_RWq_HI,log_RWq_PH,"q_RW_log_dev_jai");
  q_nad.allocate(styr_cpue_nad,endyr_cpue_nad,"q_nad");
  #ifndef NO_AD_INITIALIZE
    q_nad.initialize();
  #endif
  q_mad.allocate(styr_cpue_mad,endyr_cpue_mad,"q_mad");
  #ifndef NO_AD_INITIALIZE
    q_mad.initialize();
  #endif
  q_sad.allocate(styr_cpue_sad,endyr_cpue_sad,"q_sad");
  #ifndef NO_AD_INITIALIZE
    q_sad.initialize();
  #endif
  q_jai.allocate(styr_cpue_jai,endyr_cpue_jai,"q_jai");
  #ifndef NO_AD_INITIALIZE
    q_jai.initialize();
  #endif
  q2_jai.allocate(styr_cpue_jai,endyr_cpue_jai,"q2_jai");
  #ifndef NO_AD_INITIALIZE
    q2_jai.initialize();
  #endif
  q_mar.allocate(1,8,"q_mar");
  #ifndef NO_AD_INITIALIZE
    q_mar.initialize();
  #endif
  q_eco.allocate(9,26,"q_eco");
  #ifndef NO_AD_INITIALIZE
    q_eco.initialize();
  #endif
  L_cRn_num.allocate(styr,endyr,1,nages,"L_cRn_num");
  #ifndef NO_AD_INITIALIZE
    L_cRn_num.initialize();
  #endif
  L_cRn_mt.allocate(styr,endyr,1,nages,"L_cRn_mt");
  #ifndef NO_AD_INITIALIZE
    L_cRn_mt.initialize();
  #endif
  pred_cRn_L_mt.allocate(styr,endyr,"pred_cRn_L_mt");
  #ifndef NO_AD_INITIALIZE
    pred_cRn_L_mt.initialize();
  #endif
  L_cRs_num.allocate(styr,endyr,1,nages,"L_cRs_num");
  #ifndef NO_AD_INITIALIZE
    L_cRs_num.initialize();
  #endif
  L_cRs_mt.allocate(styr,endyr,1,nages,"L_cRs_mt");
  #ifndef NO_AD_INITIALIZE
    L_cRs_mt.initialize();
  #endif
  pred_cRs_L_mt.allocate(styr,endyr,"pred_cRs_L_mt");
  #ifndef NO_AD_INITIALIZE
    pred_cRs_L_mt.initialize();
  #endif
  L_cBn_num.allocate(styr,endyr,1,nages,"L_cBn_num");
  #ifndef NO_AD_INITIALIZE
    L_cBn_num.initialize();
  #endif
  L_cBn_mt.allocate(styr,endyr,1,nages,"L_cBn_mt");
  #ifndef NO_AD_INITIALIZE
    L_cBn_mt.initialize();
  #endif
  pred_cBn_L_mt.allocate(styr,endyr,"pred_cBn_L_mt");
  #ifndef NO_AD_INITIALIZE
    pred_cBn_L_mt.initialize();
  #endif
  L_cBs_num.allocate(styr,endyr,1,nages,"L_cBs_num");
  #ifndef NO_AD_INITIALIZE
    L_cBs_num.initialize();
  #endif
  L_cBs_mt.allocate(styr,endyr,1,nages,"L_cBs_mt");
  #ifndef NO_AD_INITIALIZE
    L_cBs_mt.initialize();
  #endif
  pred_cBs_L_mt.allocate(styr,endyr,"pred_cBs_L_mt");
  #ifndef NO_AD_INITIALIZE
    pred_cBs_L_mt.initialize();
  #endif
  L_total_num.allocate(styr,endyr,1,nages,"L_total_num");
  #ifndef NO_AD_INITIALIZE
    L_total_num.initialize();
  #endif
  L_total_mt.allocate(styr,endyr,1,nages,"L_total_mt");
  #ifndef NO_AD_INITIALIZE
    L_total_mt.initialize();
  #endif
  L_total_mt_yr.allocate(styr,endyr,"L_total_mt_yr");
  #ifndef NO_AD_INITIALIZE
    L_total_mt_yr.initialize();
  #endif
  F_cRn_prop.allocate("F_cRn_prop");
  #ifndef NO_AD_INITIALIZE
  F_cRn_prop.initialize();
  #endif
  F_cRs_prop.allocate("F_cRs_prop");
  #ifndef NO_AD_INITIALIZE
  F_cRs_prop.initialize();
  #endif
  F_cBn_prop.allocate("F_cBn_prop");
  #ifndef NO_AD_INITIALIZE
  F_cBn_prop.initialize();
  #endif
  F_cBs_prop.allocate("F_cBs_prop");
  #ifndef NO_AD_INITIALIZE
  F_cBs_prop.initialize();
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
  msy_mt_out.allocate("msy_mt_out");
  #ifndef NO_AD_INITIALIZE
  msy_mt_out.initialize();
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
  F35_dum.allocate("F35_dum");
  #ifndef NO_AD_INITIALIZE
  F35_dum.initialize();
  #endif
  F30_dum.allocate("F30_dum");
  #ifndef NO_AD_INITIALIZE
  F30_dum.initialize();
  #endif
  F40_dum.allocate("F40_dum");
  #ifndef NO_AD_INITIALIZE
  F40_dum.initialize();
  #endif
  F35_out.allocate("F35_out");
  #ifndef NO_AD_INITIALIZE
  F35_out.initialize();
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
  L_F30_mt_out.allocate("L_F30_mt_out");
  #ifndef NO_AD_INITIALIZE
  L_F30_mt_out.initialize();
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
  L_eq_mt.allocate(1,n_iter_msy,"L_eq_mt");
  #ifndef NO_AD_INITIALIZE
    L_eq_mt.initialize();
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
  SdSSB_msy.allocate(styr,endyr+1,"SdSSB_msy");
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
  SdSSB_F30.allocate(styr,endyr+1,"SdSSB_F30");
  #ifndef NO_AD_INITIALIZE
    SdSSB_F30.initialize();
  #endif
  SdSSB_F30_end.allocate("SdSSB_F30_end");
  #ifndef NO_AD_INITIALIZE
  SdSSB_F30_end.initialize();
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
  wgt_wgted_L_mt.allocate(1,nages,"wgt_wgted_L_mt");
  #ifndef NO_AD_INITIALIZE
    wgt_wgted_L_mt.initialize();
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
  M_tv.allocate(styr,endyr,1,nages,"M_tv");
  #ifndef NO_AD_INITIALIZE
    M_tv.initialize();
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
  log_avg_F_L_cRn.allocate(log_avg_F_cRn_LO,log_avg_F_cRn_HI,log_avg_F_cRn_PH,"log_avg_F_L_cRn");
  log_avg_F_cRn_out.allocate(1,8,"log_avg_F_cRn_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_cRn_out.initialize();
  #endif
  log_dev_F_L_cRn.allocate(styr_L_cRn,endyr_L_cRn,log_F_dev_cRn_LO,log_F_dev_cRn_HI,log_F_dev_cRn_PH,"log_dev_F_L_cRn");
  log_F_dev_cRn_out.allocate(styr_L_cRn,endyr_L_cRn,"log_F_dev_cRn_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_cRn_out.initialize();
  #endif
  F_cRn.allocate(styr,endyr,1,nages,"F_cRn");
  #ifndef NO_AD_INITIALIZE
    F_cRn.initialize();
  #endif
  F_cRn_out.allocate(styr,endyr,"F_cRn_out");
  #ifndef NO_AD_INITIALIZE
    F_cRn_out.initialize();
  #endif
  log_F_dev_init_cRn.allocate("log_F_dev_init_cRn");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_init_cRn.initialize();
  #endif
  log_F_dev_end_cRn.allocate("log_F_dev_end_cRn");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_cRn.initialize();
  #endif
  log_avg_F_L_cRs.allocate(log_avg_F_cRs_LO,log_avg_F_cRs_HI,log_avg_F_cRs_PH,"log_avg_F_L_cRs");
  log_avg_F_cRs_out.allocate(1,8,"log_avg_F_cRs_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_cRs_out.initialize();
  #endif
  log_dev_F_L_cRs.allocate(styr_L_cRs,endyr_L_cRs,log_F_dev_cRs_LO,log_F_dev_cRs_HI,log_F_dev_cRs_PH,"log_dev_F_L_cRs");
  log_F_dev_cRs_out.allocate(styr_L_cRs,endyr_L_cRs,"log_F_dev_cRs_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_cRs_out.initialize();
  #endif
  F_cRs.allocate(styr,endyr,1,nages,"F_cRs");
  #ifndef NO_AD_INITIALIZE
    F_cRs.initialize();
  #endif
  F_cRs_out.allocate(styr,endyr,"F_cRs_out");
  #ifndef NO_AD_INITIALIZE
    F_cRs_out.initialize();
  #endif
  log_F_dev_init_cRs.allocate("log_F_dev_init_cRs");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_init_cRs.initialize();
  #endif
  log_F_dev_end_cRs.allocate("log_F_dev_end_cRs");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_cRs.initialize();
  #endif
  log_avg_F_L_cBn.allocate(log_avg_F_cBn_LO,log_avg_F_cBn_HI,log_avg_F_cBn_PH,"log_avg_F_L_cBn");
  log_avg_F_cBn_out.allocate(1,8,"log_avg_F_cBn_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_cBn_out.initialize();
  #endif
  log_dev_F_L_cBn.allocate(styr_L_cBn,endyr_L_cBn,log_F_dev_cBn_LO,log_F_dev_cBn_HI,log_F_dev_cBn_PH,"log_dev_F_L_cBn");
  log_F_dev_cBn_out.allocate(styr_L_cBn,endyr_L_cBn,"log_F_dev_cBn_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_cBn_out.initialize();
  #endif
  F_cBn.allocate(styr,endyr,1,nages,"F_cBn");
  #ifndef NO_AD_INITIALIZE
    F_cBn.initialize();
  #endif
  F_cBn_out.allocate(styr,endyr,"F_cBn_out");
  #ifndef NO_AD_INITIALIZE
    F_cBn_out.initialize();
  #endif
  log_F_dev_init_cBn.allocate("log_F_dev_init_cBn");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_init_cBn.initialize();
  #endif
  log_F_dev_end_cBn.allocate("log_F_dev_end_cBn");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_cBn.initialize();
  #endif
  log_avg_F_L_cBs.allocate(log_avg_F_cBs_LO,log_avg_F_cBs_HI,log_avg_F_cBs_PH,"log_avg_F_L_cBs");
  log_avg_F_cBs_out.allocate(1,8,"log_avg_F_cBs_out");
  #ifndef NO_AD_INITIALIZE
    log_avg_F_cBs_out.initialize();
  #endif
  log_dev_F_L_cBs.allocate(styr_L_cBs,endyr_L_cBs,log_F_dev_cBs_LO,log_F_dev_cBs_HI,log_F_dev_cBs_PH,"log_dev_F_L_cBs");
  log_F_dev_cBs_out.allocate(styr_L_cBs,endyr_L_cBs,"log_F_dev_cBs_out");
  #ifndef NO_AD_INITIALIZE
    log_F_dev_cBs_out.initialize();
  #endif
  F_cBs.allocate(styr,endyr,1,nages,"F_cBs");
  #ifndef NO_AD_INITIALIZE
    F_cBs.initialize();
  #endif
  F_cBs_out.allocate(styr,endyr,"F_cBs_out");
  #ifndef NO_AD_INITIALIZE
    F_cBs_out.initialize();
  #endif
  log_F_dev_init_cBs.allocate("log_F_dev_init_cBs");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_init_cBs.initialize();
  #endif
  log_F_dev_end_cBs.allocate("log_F_dev_end_cBs");
  #ifndef NO_AD_INITIALIZE
  log_F_dev_end_cBs.initialize();
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
  sdnr_lc_nad.allocate("sdnr_lc_nad");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_nad.initialize();
  #endif
  sdnr_lc_mad.allocate("sdnr_lc_mad");
  #ifndef NO_AD_INITIALIZE
  sdnr_lc_mad.initialize();
  #endif
  sdnr_ac_cRn.allocate("sdnr_ac_cRn");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_cRn.initialize();
  #endif
  sdnr_ac_cRs.allocate("sdnr_ac_cRs");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_cRs.initialize();
  #endif
  sdnr_ac_cBn.allocate("sdnr_ac_cBn");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_cBn.initialize();
  #endif
  sdnr_ac_cBs.allocate("sdnr_ac_cBs");
  #ifndef NO_AD_INITIALIZE
  sdnr_ac_cBs.initialize();
  #endif
  sdnr_I_nad.allocate("sdnr_I_nad");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_nad.initialize();
  #endif
  sdnr_I_mad.allocate("sdnr_I_mad");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_mad.initialize();
  #endif
  sdnr_I_sad.allocate("sdnr_I_sad");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_sad.initialize();
  #endif
  sdnr_I_jai.allocate("sdnr_I_jai");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_jai.initialize();
  #endif
  sdnr_I_mareco.allocate("sdnr_I_mareco");
  #ifndef NO_AD_INITIALIZE
  sdnr_I_mareco.initialize();
  #endif
  w_L.allocate("w_L");
  #ifndef NO_AD_INITIALIZE
  w_L.initialize();
  #endif
  w_cpue_nad.allocate("w_cpue_nad");
  #ifndef NO_AD_INITIALIZE
  w_cpue_nad.initialize();
  #endif
  w_cpue_mad.allocate("w_cpue_mad");
  #ifndef NO_AD_INITIALIZE
  w_cpue_mad.initialize();
  #endif
  w_cpue_sad.allocate("w_cpue_sad");
  #ifndef NO_AD_INITIALIZE
  w_cpue_sad.initialize();
  #endif
  w_cpue_jai.allocate("w_cpue_jai");
  #ifndef NO_AD_INITIALIZE
  w_cpue_jai.initialize();
  #endif
  w_cpue_mareco.allocate("w_cpue_mareco");
  #ifndef NO_AD_INITIALIZE
  w_cpue_mareco.initialize();
  #endif
  w_lenc_nad.allocate("w_lenc_nad");
  #ifndef NO_AD_INITIALIZE
  w_lenc_nad.initialize();
  #endif
  w_lenc_mad.allocate("w_lenc_mad");
  #ifndef NO_AD_INITIALIZE
  w_lenc_mad.initialize();
  #endif
  w_agec_cRn.allocate("w_agec_cRn");
  #ifndef NO_AD_INITIALIZE
  w_agec_cRn.initialize();
  #endif
  w_agec_cRs.allocate("w_agec_cRs");
  #ifndef NO_AD_INITIALIZE
  w_agec_cRs.initialize();
  #endif
  w_agec_cBn.allocate("w_agec_cBn");
  #ifndef NO_AD_INITIALIZE
  w_agec_cBn.initialize();
  #endif
  w_agec_cBs.allocate("w_agec_cBs");
  #ifndef NO_AD_INITIALIZE
  w_agec_cBs.initialize();
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
  f_cRn_L.allocate("f_cRn_L");
  #ifndef NO_AD_INITIALIZE
  f_cRn_L.initialize();
  #endif
  f_cRs_L.allocate("f_cRs_L");
  #ifndef NO_AD_INITIALIZE
  f_cRs_L.initialize();
  #endif
  f_cBn_L.allocate("f_cBn_L");
  #ifndef NO_AD_INITIALIZE
  f_cBn_L.initialize();
  #endif
  f_cBs_L.allocate("f_cBs_L");
  #ifndef NO_AD_INITIALIZE
  f_cBs_L.initialize();
  #endif
  f_nad_cpue.allocate("f_nad_cpue");
  #ifndef NO_AD_INITIALIZE
  f_nad_cpue.initialize();
  #endif
  f_mad_cpue.allocate("f_mad_cpue");
  #ifndef NO_AD_INITIALIZE
  f_mad_cpue.initialize();
  #endif
  f_sad_cpue.allocate("f_sad_cpue");
  #ifndef NO_AD_INITIALIZE
  f_sad_cpue.initialize();
  #endif
  f_jai_cpue.allocate("f_jai_cpue");
  #ifndef NO_AD_INITIALIZE
  f_jai_cpue.initialize();
  #endif
  f_mareco_cpue.allocate("f_mareco_cpue");
  #ifndef NO_AD_INITIALIZE
  f_mareco_cpue.initialize();
  #endif
  f_nad_lenc.allocate("f_nad_lenc");
  #ifndef NO_AD_INITIALIZE
  f_nad_lenc.initialize();
  #endif
  f_mad_lenc.allocate("f_mad_lenc");
  #ifndef NO_AD_INITIALIZE
  f_mad_lenc.initialize();
  #endif
  f_cRn_agec.allocate("f_cRn_agec");
  #ifndef NO_AD_INITIALIZE
  f_cRn_agec.initialize();
  #endif
  f_cRs_agec.allocate("f_cRs_agec");
  #ifndef NO_AD_INITIALIZE
  f_cRs_agec.initialize();
  #endif
  f_cBn_agec.allocate("f_cBn_agec");
  #ifndef NO_AD_INITIALIZE
  f_cBn_agec.initialize();
  #endif
  f_cBs_agec.allocate("f_cBs_agec");
  #ifndef NO_AD_INITIALIZE
  f_cBs_agec.initialize();
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
  L_mt_proj.allocate(styr_proj,endyr_proj,"L_mt_proj");
  #ifndef NO_AD_INITIALIZE
    L_mt_proj.initialize();
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
  dvector temp1("{1000, 2000,3000, 5000, 10000, 10000, 10000;}");
  maximum_function_evaluations.allocate(temp1.indexmin(),temp1.indexmax());
  maximum_function_evaluations=temp1;
  dvector temp("{1e-2, 1e-2,1e-3, 1e-3, 1e-4, 1e-4, 1e-4;}");
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
  len_cv_nad_val=set_len_cv_nad(1);
  len_cv_mad_val=set_len_cv_mad(1);
  //len_cv_sad_val=set_len_cv_sad(1);
  
  M=set_M;
  M_tv=set_M_tv;
    
  log_R0=set_log_R0(1);
  steep=set_steep(1);
  R_autocorr=set_R_autocorr(1);
  rec_sigma=set_rec_sigma(1);
  
  log_dm_lenc_nad=set_log_dm_lenc_nad(1);
  log_dm_lenc_mad=set_log_dm_lenc_mad(1);
  //log_dm_sad_lc=set_log_dm_sad_lc(1);
  log_dm_agec_cRn=set_log_dm_agec_cRn(1);
  log_dm_agec_cRs=set_log_dm_agec_cRs(1);
  log_dm_agec_cBn=set_log_dm_agec_cBn(1);
  log_dm_agec_cBs=set_log_dm_agec_cBs(1);
  log_q_cpue_nad=set_log_q_cpue_nad(1);
  log_q_cpue_mad=set_log_q_cpue_mad(1);
  log_q_cpue_sad=set_log_q_cpue_sad(1);
  log_q_cpue_jai=set_log_q_cpue_jai(1);
  log_q2_jai=set_log_q2_jai(1);
  log_q_cpue_mar=set_log_q_cpue_mar(1);
  log_q_cpue_eco=set_log_q_cpue_eco(1);
 
  q_rate=set_q_rate;
  q_rate_fcn_nad=1.0;
  q_rate_fcn_mad=1.0;
  q_rate_fcn_sad=1.0;
  q_rate_fcn_jai=1.0;
  //q_rate_fcn_mar=1.0;
  //q_rate_fcn_eco=1.0;
  q_DD_beta=set_q_DD_beta;
  q_DD_fcn=1.0;
  q_RW_log_dev_nad.initialize(); 
  q_RW_log_dev_mad.initialize();
  q_RW_log_dev_sad.initialize();
  q_RW_log_dev_jai.initialize(); 
  //q_RW_log_dev_mar.initialize(); 
  //q_RW_log_dev_eco.initialize(); 
   
   if (set_q_rate_phase<0 & q_rate!=0.0)
  {
    for (iyear=styr_cpue_nad; iyear<=endyr_cpue_nad; iyear++)
      {   if (iyear>styr_cpue_nad & iyear <=2003) 
          {//q_rate_fcn_nad(iyear)=(1.0+q_rate)*q_rate_fcn_nad(iyear-1); //compound
             q_rate_fcn_nad(iyear)=(1.0+(iyear-styr_cpue_nad)*q_rate)*q_rate_fcn_nad(styr_cpue_nad);  //linear
          }
          if (iyear>2003) {q_rate_fcn_nad(iyear)=q_rate_fcn_nad(iyear-1);} 
      }
    for (iyear=styr_cpue_mad; iyear<=endyr_cpue_mad; iyear++)
      {   if (iyear>styr_cpue_mad & iyear <=2003) 
          {//q_rate_fcn_mad(iyear)=(1.0+q_rate)*q_rate_fcn_mad(iyear-1); //compound
             q_rate_fcn_mad(iyear)=(1.0+(iyear-styr_cpue_mad)*q_rate)*q_rate_fcn_mad(styr_cpue_mad);  //linear
          }
          if (iyear>2003) {q_rate_fcn_mad(iyear)=q_rate_fcn_mad(iyear-1);} 
      }
    for (iyear=styr_cpue_sad; iyear<=endyr_cpue_sad; iyear++)
      {   if (iyear>styr_cpue_sad & iyear <=2003) 
          {//q_rate_fcn_sad(iyear)=(1.0+q_rate)*q_rate_fcn_sad(iyear-1); //compound
             q_rate_fcn_sad(iyear)=(1.0+(iyear-styr_cpue_sad)*q_rate)*q_rate_fcn_sad(styr_cpue_sad);  //linear
          }
          if (iyear>2003) {q_rate_fcn_sad(iyear)=q_rate_fcn_sad(iyear-1);} 
      }
    for (iyear=styr_cpue_jai; iyear<=endyr_cpue_jai; iyear++)
      {   if (iyear>styr_cpue_jai & iyear <=2003) 
          {//q_rate_fcn_jai(iyear)=(1.0+q_rate)*q_rate_fcn_jai(iyear-1); //compound
             q_rate_fcn_jai(iyear)=(1.0+(iyear-styr_cpue_jai)*q_rate)*q_rate_fcn_jai(styr_cpue_jai);  //linear
          }
          if (iyear>2003) {q_rate_fcn_jai(iyear)=q_rate_fcn_jai(iyear-1);} 
      }
    //for (iyear=styr_mar_cpue; iyear<=endyr_mar_cpue; iyear++)
    // {   if (iyear>styr_mar_cpue & iyear <=2003) 
    //      {//q_rate_fcn_mar(iyear)=(1.0+q_rate)*q_rate_fcn_mar(iyear-1); //compound
    //         q_rate_fcn_mar(iyear)=(1.0+(iyear-styr_mar_cpue)*q_rate)*q_rate_fcn_mar(styr_mar_cpue);  //linear
    //      }
    //      if (iyear>2003) {q_rate_fcn_mar(iyear)=q_rate_fcn_mar(iyear-1);} 
    //  }
    //for (iyear=styr_eco_cpue; iyear<=endyr_eco_cpue; iyear++)
    //  {   if (iyear>styr_eco_cpue & iyear <=2003) 
    //      {//q_rate_fcn_eco(iyear)=(1.0+q_rate)*q_rate_fcn_eco(iyear-1); //compound
    //         q_rate_fcn_eco(iyear)=(1.0+(iyear-styr_eco_cpue)*q_rate)*q_rate_fcn_eco(styr_eco_cpue);  //linear
    //      }
    //      if (iyear>2003) {q_rate_fcn_eco(iyear)=q_rate_fcn_eco(iyear-1);} 
    //  }       
  } //end q_rate conditional      
  w_L=set_w_L;
  
  w_cpue_nad=set_w_cpue_nad;
  w_cpue_mad=set_w_cpue_mad;
  w_cpue_sad=set_w_cpue_sad;
  w_cpue_jai=set_w_cpue_jai;
  w_cpue_mareco=set_w_cpue_mareco;
  
  w_lenc_nad=set_w_lenc_nad;
  w_lenc_mad=set_w_lenc_mad;
  //w_lc_sad=set_w_lc_sad;
  
  w_agec_cRn=set_w_agec_cRn;
  w_agec_cRs=set_w_agec_cRs;
  w_agec_cBn=set_w_agec_cBn;
  w_agec_cBs=set_w_agec_cBs; 
  
  w_Nage_init=set_w_Nage_init;
  w_rec=set_w_rec;
  w_rec_early=set_w_rec_early;
  w_rec_end=set_w_rec_end;
  w_fullF=set_w_fullF;
  w_Ftune=set_w_Ftune;
  log_avg_F_L_cRn=set_log_avg_F_L_cRn(1);
  log_avg_F_L_cRs=set_log_avg_F_L_cRs(1);
  log_avg_F_L_cBn=set_log_avg_F_L_cBn(1);
  log_avg_F_L_cBs=set_log_avg_F_L_cBs(1);
  log_dev_F_L_cRn=set_log_dev_vals_F_L__cRn;
  log_dev_F_L_cRs=set_log_dev_vals_F_L__cRs;
  log_dev_F_L_cBn=set_log_dev_vals_F_L__cBn;
  log_dev_F_L_cBs=set_log_dev_vals_F_L__cBs;
  selpar_A50_cRn=set_selpar_A50_cRn(1);         //Block 1
  selpar_slope_cRn=set_selpar_slope_cRn(1);
  selpar_A502_cRn=set_selpar_A502_cRn(1);
  selpar_slope2_cRn=set_selpar_slope2_cRn(1);
  selpar_A50_cRn2=set_selpar_A50_cRn2(1);        //Block 2
  selpar_slope_cRn2=set_selpar_slope_cRn2(1);
  selpar_A502_cRn2=set_selpar_A502_cRn2(1);
  selpar_slope2_cRn2=set_selpar_slope2_cRn2(1);
  selpar_A50_cRn3=set_selpar_A50_cRn3(1);         //Block 3 
  selpar_slope_cRn3=set_selpar_slope_cRn3(1);
  selpar_A502_cRn3=set_selpar_A502_cRn3(1);
  selpar_slope2_cRn3=set_selpar_slope2_cRn3(1);
  selpar_A50_cRn4=set_selpar_A50_cRn4(1);         //Block 4 
  selpar_slope_cRn4=set_selpar_slope_cRn4(1);
  selpar_A502_cRn4=set_selpar_A502_cRn4(1);
  selpar_slope2_cRn4=set_selpar_slope2_cRn4(1);
  selpar_A50_cRs=set_selpar_A50_cRs(1);            //Block 1
  selpar_slope_cRs=set_selpar_slope_cRs(1);
  selpar_A502_cRs=set_selpar_A502_cRs(1);
  selpar_slope2_cRs=set_selpar_slope2_cRs(1);
  selpar_A50_cRs2=set_selpar_A50_cRs2(1);           //Block 2
  selpar_slope_cRs2=set_selpar_slope_cRs2(1);
  selpar_A502_cRs2=set_selpar_A502_cRs2(1);
  selpar_slope2_cRs2=set_selpar_slope2_cRs2(1);
  selpar_A50_cRs3=set_selpar_A50_cRs3(1);           //Block 3
  selpar_slope_cRs3=set_selpar_slope_cRs3(1);
  selpar_A502_cRs3=set_selpar_A502_cRs3(1);
  selpar_slope2_cRs3=set_selpar_slope2_cRs3(1);
  
  selpar_A50_cRs4=set_selpar_A50_cRs4(1);           //Block 4
  selpar_slope_cRs4=set_selpar_slope_cRs4(1);
  selpar_A502_cRs4=set_selpar_A502_cRs4(1);
  selpar_slope2_cRs4=set_selpar_slope2_cRs4(1);
  sel_age0_cRn_logit=set_sel_age0_cRn(1);  //setting cR selectivity at age in logit space - block 1
  sel_age1_cRn_logit=set_sel_age1_cRn(1);
  sel_age2_cRn_logit=set_sel_age2_cRn(1);
  sel_age3_cRn_logit=set_sel_age3_cRn(1);
  sel_age4_cRn_logit=set_sel_age4_cRn(1);
  sel_age5_cRn_logit=set_sel_age5_cRn(1);
  sel_age6_cRn_logit=set_sel_age6_cRn(1);
  sel_age0_cRn2_logit=set_sel_age0_cRn2(1);  //setting cR selectivity at age in logit space - block 2
  sel_age1_cRn2_logit=set_sel_age1_cRn2(1);
  sel_age2_cRn2_logit=set_sel_age2_cRn2(1);
  sel_age3_cRn2_logit=set_sel_age3_cRn2(1);
  sel_age4_cRn2_logit=set_sel_age4_cRn2(1);
  sel_age5_cRn2_logit=set_sel_age5_cRn2(1);
  sel_age6_cRn2_logit=set_sel_age6_cRn2(1);
  sel_age0_cRn3_logit=set_sel_age0_cRn3(1);  //setting cR selectivity at age in logit space - block 3
  sel_age1_cRn3_logit=set_sel_age1_cRn3(1);
  sel_age2_cRn3_logit=set_sel_age2_cRn3(1);
  sel_age3_cRn3_logit=set_sel_age3_cRn3(1);
  sel_age4_cRn3_logit=set_sel_age4_cRn3(1);
  sel_age5_cRn3_logit=set_sel_age5_cRn3(1);
  sel_age6_cRn3_logit=set_sel_age6_cRn3(1);
  sel_age0_cRn4_logit=set_sel_age0_cRn4(1);  //setting cR selectivity at age in logit space - block 4
  sel_age1_cRn4_logit=set_sel_age1_cRn4(1);
  sel_age2_cRn4_logit=set_sel_age2_cRn4(1);
  sel_age3_cRn4_logit=set_sel_age3_cRn4(1);
  sel_age4_cRn4_logit=set_sel_age4_cRn4(1);
  sel_age5_cRn4_logit=set_sel_age5_cRn4(1);
  sel_age6_cRn4_logit=set_sel_age6_cRn4(1);
  sel_age0_cRs_logit=set_sel_age0_cRs(1);  //setting cR selectivity at age in logit space - block 1
  sel_age1_cRs_logit=set_sel_age1_cRs(1);
  sel_age2_cRs_logit=set_sel_age2_cRs(1);
  sel_age3_cRs_logit=set_sel_age3_cRs(1);
  sel_age4_cRs_logit=set_sel_age4_cRs(1);
  sel_age5_cRs_logit=set_sel_age5_cRs(1);
  sel_age6_cRs_logit=set_sel_age6_cRs(1);
  sel_age0_cRs2_logit=set_sel_age0_cRs2(1);  //setting cR selectivity at age in logit space - block 2
  sel_age1_cRs2_logit=set_sel_age1_cRs2(1);
  sel_age2_cRs2_logit=set_sel_age2_cRs2(1);
  sel_age3_cRs2_logit=set_sel_age3_cRs2(1);
  sel_age4_cRs2_logit=set_sel_age4_cRs2(1);
  sel_age5_cRs2_logit=set_sel_age5_cRs2(1);
  sel_age6_cRs2_logit=set_sel_age6_cRs2(1);
  sel_age0_cRs3_logit=set_sel_age0_cRs3(1);  //setting cR selectivity at age in logit space - block 3
  sel_age1_cRs3_logit=set_sel_age1_cRs3(1);
  sel_age2_cRs3_logit=set_sel_age2_cRs3(1);
  sel_age3_cRs3_logit=set_sel_age3_cRs3(1);
  sel_age4_cRs3_logit=set_sel_age4_cRs3(1);
  sel_age5_cRs3_logit=set_sel_age5_cRs3(1);
  sel_age6_cRs3_logit=set_sel_age6_cRs3(1);
  sel_age0_cRs4_logit=set_sel_age0_cRs4(1);  //setting cR selectivity at age in logit space - block 4
  sel_age1_cRs4_logit=set_sel_age1_cRs4(1);
  sel_age2_cRs4_logit=set_sel_age2_cRs4(1);
  sel_age3_cRs4_logit=set_sel_age3_cRs4(1);
  sel_age4_cRs4_logit=set_sel_age4_cRs4(1);
  sel_age5_cRs4_logit=set_sel_age5_cRs4(1);
  sel_age6_cRs4_logit=set_sel_age6_cRs4(1);
  selpar_A50_cBn=set_selpar_A50_cBn(1);               //Block 1
  selpar_slope_cBn=set_selpar_slope_cBn(1);
  selpar_A502_cBn=set_selpar_A502_cBn(1);
  selpar_slope2_cBn=set_selpar_slope2_cBn(1);
  selpar_A50_cBn2=set_selpar_A50_cBn2(1);             //Block 2
  selpar_slope_cBn2=set_selpar_slope_cBn2(1);
  selpar_A502_cBn2=set_selpar_A502_cBn2(1);
  selpar_slope2_cBn2=set_selpar_slope2_cBn2(1);
  selpar_A50_cBs=set_selpar_A50_cBs(1);              //Block 1
  selpar_slope_cBs=set_selpar_slope_cBs(1);
  selpar_A502_cBs=set_selpar_A502_cBs(1);
  selpar_slope2_cBs=set_selpar_slope2_cBs(1);
  selpar_A50_cBs2=set_selpar_A50_cBs2(1);            //Block 2
  selpar_slope_cBs2=set_selpar_slope_cBs2(1);
  selpar_A502_cBs2=set_selpar_A502_cBs2(1);
  selpar_slope2_cBs2=set_selpar_slope2_cBs2(1);
  sel_age0_cBn_logit=set_sel_age0_cBn(1);  //setting cB selectivity at age in logit space - block 1
  sel_age1_cBn_logit=set_sel_age1_cBn(1);
  sel_age2_cBn_logit=set_sel_age2_cBn(1);
  sel_age3_cBn_logit=set_sel_age3_cBn(1);
  sel_age4_cBn_logit=set_sel_age4_cBn(1);
  sel_age5_cBn_logit=set_sel_age5_cBn(1);
  sel_age6_cBn_logit=set_sel_age6_cBn(1);
  sel_age0_cBn2_logit=set_sel_age0_cBn2(1);  //setting cB selectivity at age in logit space - block 2
  sel_age1_cBn2_logit=set_sel_age1_cBn2(1);
  sel_age2_cBn2_logit=set_sel_age2_cBn2(1);
  sel_age3_cBn2_logit=set_sel_age3_cBn2(1);
  sel_age4_cBn2_logit=set_sel_age4_cBn2(1);
  sel_age5_cBn2_logit=set_sel_age5_cBn2(1);
  sel_age6_cBn2_logit=set_sel_age6_cBn2(1);
  sel_age0_cBs_logit=set_sel_age0_cBs(1);  //setting cB selectivity at age in logit space - block 1
  sel_age1_cBs_logit=set_sel_age1_cBs(1);
  sel_age2_cBs_logit=set_sel_age2_cBs(1);
  sel_age3_cBs_logit=set_sel_age3_cBs(1);
  sel_age4_cBs_logit=set_sel_age4_cBs(1);
  sel_age5_cBs_logit=set_sel_age5_cBs(1);
  sel_age6_cBs_logit=set_sel_age6_cBs(1);
  
  sel_age0_cBs2_logit=set_sel_age0_cBs2(1);  //setting cB selectivity at age in logit space - block 2
  sel_age1_cBs2_logit=set_sel_age1_cBs2(1);
  sel_age2_cBs2_logit=set_sel_age2_cBs2(1);
  sel_age3_cBs2_logit=set_sel_age3_cBs2(1);
  sel_age4_cBs2_logit=set_sel_age4_cBs2(1);
  sel_age5_cBs2_logit=set_sel_age5_cBs2(1);
  sel_age6_cBs2_logit=set_sel_age6_cBs2(1);
  
  selpar_A50_nad=set_selpar_A50_nad(1);
  selpar_slope_nad=set_selpar_slope_nad(1);
  selpar_A502_nad=set_selpar_A502_nad(1);
  selpar_slope2_nad=set_selpar_slope2_nad(1);
  sel_age0_nad_logit=set_sel_age0_nad(1);  //setting NAD selectivity at age in logit space
  sel_age1_nad_logit=set_sel_age1_nad(1);
  sel_age2_nad_logit=set_sel_age2_nad(1);
  sel_age3_nad_logit=set_sel_age3_nad(1);
  sel_age4_nad_logit=set_sel_age4_nad(1);
  sel_age5_nad_logit=set_sel_age5_nad(1);
  sel_age6_nad_logit=set_sel_age6_nad(1);
  selpar_A50_mad=set_selpar_A50_mad(1);
  selpar_slope_mad=set_selpar_slope_mad(1);
  selpar_A502_mad=set_selpar_A502_mad(1);
  selpar_slope2_mad=set_selpar_slope2_mad(1);
  sel_age0_mad_logit=set_sel_age0_mad(1);  //setting MAD selectivity at age in logit space
  sel_age1_mad_logit=set_sel_age1_mad(1);
  sel_age2_mad_logit=set_sel_age2_mad(1);
  sel_age3_mad_logit=set_sel_age3_mad(1);
  sel_age4_mad_logit=set_sel_age4_mad(1);
  sel_age5_mad_logit=set_sel_age5_mad(1);
  sel_age6_mad_logit=set_sel_age6_mad(1);
  selpar_A50_sad=set_selpar_A50_sad(1);
  selpar_slope_sad=set_selpar_slope_sad(1);
  selpar_A502_sad=set_selpar_A502_sad(1);
  selpar_slope2_sad=set_selpar_slope2_sad(1);
  sel_age0_sad_logit=set_sel_age0_sad(1);  //setting SAD selectivity at age in logit space
  sel_age1_sad_logit=set_sel_age1_sad(1);
  sel_age2_sad_logit=set_sel_age2_sad(1);
  sel_age3_sad_logit=set_sel_age3_sad(1);
  sel_age4_sad_logit=set_sel_age4_sad(1);
  sel_age5_sad_logit=set_sel_age5_sad(1);
  sel_age6_sad_logit=set_sel_age6_sad(1);
   
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
 tv_maturity=maturity_obs_tv;
 prop_f=obs_prop_f;
 prop_m=1.0-obs_prop_f;
  
      nsamp_nad_lenc_allyr=missing;
      nsamp_mad_lenc_allyr=missing;
      //nsamp_sad_lenc_allyr=missing;
      
      nsamp_cRn_agec_allyr=missing;
      nsamp_cRs_agec_allyr=missing;
      nsamp_cBn_agec_allyr=missing;
      nsamp_cBs_agec_allyr=missing; 
      
      nfish_nad_lenc_allyr=missing;
      nfish_mad_lenc_allyr=missing;
      //nfish_sad_lenc_allyr=missing;
      
      nfish_cRn_agec_allyr=missing;
      nfish_cRs_agec_allyr=missing;
      nfish_cBn_agec_allyr=missing;
      nfish_cBs_agec_allyr=missing;
   
      for (iyear=1; iyear<=nyr_lenc_nad; iyear++)
         {if (nsamp_lenc_nad(iyear)>=minSS_lenc_nad)
           {nsamp_nad_lenc_allyr(yrs_lenc_nad(iyear))=nsamp_lenc_nad(iyear);
            nfish_nad_lenc_allyr(yrs_lenc_nad(iyear))=nfish_lenc_nad(iyear);}}
      for (iyear=1; iyear<=nyr_lenc_mad; iyear++)
         {if (nsamp_lenc_mad(iyear)>=minSS_lenc_mad)
           {nsamp_mad_lenc_allyr(yrs_lenc_mad(iyear))=nsamp_lenc_mad(iyear);
            nfish_mad_lenc_allyr(yrs_lenc_mad(iyear))=nfish_lenc_mad(iyear);}}
      //for (iyear=1; iyear<=nyr_sad_lenc; iyear++)
      //   {if (nsamp_sad_lenc(iyear)>=minSS_sad_lenc)
      //     {nsamp_sad_lenc_allyr(yrs_sad_lenc(iyear))=nsamp_sad_lenc(iyear);
      //      nfish_sad_lenc_allyr(yrs_sad_lenc(iyear))=nfish_sad_lenc(iyear);}}
      for (iyear=1; iyear<=nyr_agec_cRn; iyear++)
         {if (nsamp_agec_cRn(iyear)>=minSS_agec_cRn)
           {nsamp_cRn_agec_allyr(yrs_agec_cRn(iyear))=nsamp_agec_cRn(iyear);
            nfish_cRn_agec_allyr(yrs_agec_cRn(iyear))=nfish_agec_cRn(iyear);}}
      for (iyear=1; iyear<=nyr_agec_cRs; iyear++)
         {if (nsamp_agec_cRs(iyear)>=minSS_agec_cRs)
           {nsamp_cRs_agec_allyr(yrs_agec_cRs(iyear))=nsamp_agec_cRs(iyear);
            nfish_cRs_agec_allyr(yrs_agec_cRs(iyear))=nfish_agec_cRs(iyear);}}
            
      for (iyear=1; iyear<=nyr_agec_cBn; iyear++)
         {if (nsamp_agec_cBn(iyear)>=minSS_agec_cBn)
           {nsamp_cBn_agec_allyr(yrs_agec_cBn(iyear))=nsamp_agec_cBn(iyear);
            nfish_cBn_agec_allyr(yrs_agec_cBn(iyear))=nfish_agec_cBn(iyear);}}
            
      for (iyear=1; iyear<=nyr_agec_cBs; iyear++)
         {if (nsamp_agec_cBs(iyear)>=minSS_agec_cBs)
           {nsamp_cBs_agec_allyr(yrs_agec_cBs(iyear))=nsamp_agec_cBs(iyear);
            nfish_cBs_agec_allyr(yrs_agec_cBs(iyear))=nfish_agec_cBs(iyear);}}
  F_msy(1)=0.0;  
  for (ff=2;ff<=n_iter_msy;ff++) {F_msy(ff)=F_msy(ff-1)+iter_inc_msy;}
  F_spr(1)=0.0;  
  for (ff=2;ff<=n_iter_spr;ff++) {F_spr(ff)=F_spr(ff-1)+iter_inc_spr;}
  F_cRn.initialize(); L_cRn_num.initialize();
  F_cRs.initialize(); L_cRs_num.initialize();
  F_cBn.initialize(); L_cBn_num.initialize();
  F_cBs.initialize(); L_cBs_num.initialize();
  
  F_cRn_out.initialize();
  F_cRs_out.initialize();
  F_cBn_out.initialize();
  F_cBs_out.initialize();
  
  sel_cRn.initialize();
  sel_cRs.initialize();
  sel_cBn.initialize();
  sel_cBs.initialize();
  
  sel_nad.initialize();
  sel_mad.initialize();
  sel_sad.initialize();  
  
  sel_cRn_block1.initialize();
  sel_cRn_block2.initialize();
  sel_cRn_block3.initialize();
  sel_cRn_block4.initialize();
  sel_cRs_block1.initialize();
  sel_cRs_block2.initialize();
  sel_cRs_block3.initialize();
  sel_cRs_block4.initialize();
  sel_cBn_block1.initialize();
  sel_cBn_block2.initialize();
  sel_cBs_block1.initialize();
  sel_cBs_block2.initialize();
  
  sel_nad_block1.initialize();
  sel_mad_block1.initialize();
  sel_sad_block1.initialize();
    
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
  //cout << "got mortalities" << endl;
  get_bias_corr(); 
  //cout<< "got recruitment bias correction" << endl;
  get_numbers_at_age(); 
  //cout << "got numbers at age" << endl;
  get_landings_numbers();
  //cout << "got landings in numbers" << endl;
  get_landings_wgt();
  //cout << "got landings in wgt" << endl;
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
	//population fork length in mm
    //compute mean length (mm FL) and weight (whole) at age
    meanlen_FL=Linf*(1.0-mfexp(-K*(agebins-t0)));
    meanlen_FL_apr15=len_apr15_tv;
    meanlen_FL_jun1=len_jun1_tv;
    meanlen_FL_oct15=len_oct15_tv;
    wgt_fish_mt=g2mt*wgt_middle;  //wgt in mt - middle of year
    wgt_spawn_mt=g2mt*wgt_spawn;  //wgt in mt - spawning (Mar 1)
    tv_wgt_middle_mt=g2mt*wgt_middle_tv;
    tv_wgt_spawn_mt=g2mt*wgt_spawn_tv;
}

void model_parameters::get_reprod(void)
{
    reprod=elem_prod(elem_prod(prop_f,maturity_f),fecundity);
    for (iyear=styr; iyear<=endyr; iyear++)
    {
      reprod_tv(iyear)=elem_prod(elem_prod(prop_f,tv_maturity(iyear)),fecundity_tv(iyear));
    }
}

void model_parameters::get_length_at_age_dist(void)
{
   //compute matrix of length at age, based on the normal distribution
    //population
  for (iyear=styr; iyear<=endyr; iyear++)
  {
  for (iage=1;iage<=nages;iage++)
   {
    len_cv_mad(iage)=len_cv_mad_val;
    len_sd_mad(iyear,iage)=meanlen_FL_apr15(iyear,iage)*len_cv_mad(iage);
	zscore_lzero=(0.0-meanlen_FL_apr15(iyear,iage))/len_sd_mad(iyear,iage); 
	cprob_lzero=cumd_norm(zscore_lzero);
    //first length bin
	//population
    zscore_len=((lenbins(1)+0.5*lenbins_width)-meanlen_FL_apr15(iyear,iage)) / len_sd_mad(iyear,iage);
    cprob_lenvec(1)=cumd_norm(zscore_len);          //includes any probability mass below zero
    lenprob_apr15(iyear,iage,1)=cprob_lenvec(1)-cprob_lzero;    //removes any probability mass below zero
    //most other length bins  
    //population
    for (ilen=2;ilen<nlenbins;ilen++)
      {
        zscore_len=((lenbins(ilen)+0.5*lenbins_width)-meanlen_FL_apr15(iyear,iage)) / len_sd_mad(iyear,iage); 
		cprob_lenvec(ilen)=cumd_norm(zscore_len);
        lenprob_apr15(iyear,iage,ilen)=cprob_lenvec(ilen)-cprob_lenvec(ilen-1);
      }
    //last length bin is a plus group
	//population
    zscore_len=((lenbins(nlenbins)-0.5*lenbins_width)-meanlen_FL_apr15(iyear,iage)) / len_sd_mad(iyear,iage); 
	lenprob_apr15(iyear,iage,nlenbins)=1.0-cumd_norm(zscore_len);
      lenprob_apr15(iyear,iage)=lenprob_apr15(iyear,iage)/(1.0-cprob_lzero);  //renormalize to account for any prob mass below size=0
   }
  //fleet and survey specific length probs
  //lenprob_mad(iyear)=lenprob_apr15(iyear);
  lenprob_mad=lenprob_apr15;
  }
  //cout << "lenprob mad " << lenprob_mad << endl;
  //for (iyear=styr; iyear<=endyr; iyear++)
  //{
  //for (iage=1;iage<=nages;iage++)
  // {
  //  len_cv_sad(iage)=len_cv_sad_val;
  //  len_sd_sad(iyear,iage)=meanlen_FL_apr15(iyear,iage)*len_cv_sad(iage);
  //	zscore_lzero=(0.0-meanlen_FL_apr15(iyear,iage))/len_sd_sad(iyear,iage); 
  //	cprob_lzero=cumd_norm(zscore_lzero);
    //first length bin
	//population
  //  zscore_len=((lenbins(1)+0.5*lenbins_width)-meanlen_FL_apr15(iyear,iage)) / len_sd_sad(iyear,iage);
  //  cprob_lenvec(1)=cumd_norm(zscore_len);          //includes any probability mass below zero
  //  lenprob_apr15(iyear,iage,1)=cprob_lenvec(1)-cprob_lzero;    //removes any probability mass below zero
    //most other length bins  
    //population
  //  for (ilen=2;ilen<nlenbins;ilen++)
  //    {
  //      zscore_len=((lenbins(ilen)+0.5*lenbins_width)-meanlen_FL_apr15(iyear,iage)) / len_sd_sad(iyear,iage); 
  //		cprob_lenvec(ilen)=cumd_norm(zscore_len);
  //      lenprob_apr15(iyear,iage,ilen)=cprob_lenvec(ilen)-cprob_lenvec(ilen-1);
  //    }
    //last length bin is a plus group
	//population
  //  zscore_len=((lenbins(nlenbins)-0.5*lenbins_width)-meanlen_FL_apr15(iyear,iage)) / len_sd_sad(iyear,iage); 
  //	lenprob_apr15(iyear,iage,nlenbins)=1.0-cumd_norm(zscore_len);
  //    lenprob_apr15(iyear,iage)=lenprob_apr15(iyear,iage)/(1.0-cprob_lzero);  //renormalize to account for any prob mass below size=0
  // }
  //fleet and survey specific length probs
  //lenprob_sad(iyear)=lenprob_jun1(iyear);
  //lenprob_sad=lenprob_apr15;
  //}
  for (iyear=styr; iyear<=endyr; iyear++)
  {
  for (iage=1;iage<=nages;iage++)
   {
    len_cv_nad(iage)=len_cv_nad_val;
    len_sd_nad(iyear,iage)=meanlen_FL_oct15(iyear,iage)*len_cv_nad(iage);
	zscore_lzero=(0.0-meanlen_FL_oct15(iyear,iage))/len_sd_nad(iyear,iage); 
	cprob_lzero=cumd_norm(zscore_lzero);
    //first length bin
	//population
    zscore_len=((lenbins(1)+0.5*lenbins_width)-meanlen_FL_oct15(iyear,iage)) / len_sd_nad(iyear,iage);
    cprob_lenvec(1)=cumd_norm(zscore_len);          //includes any probability mass below zero
    lenprob_oct15(iyear,iage,1)=cprob_lenvec(1)-cprob_lzero;    //removes any probability mass below zero
    //most other length bins  
    //population
    for (ilen=2;ilen<nlenbins;ilen++)
      {
        zscore_len=((lenbins(ilen)+0.5*lenbins_width)-meanlen_FL_oct15(iyear,iage)) / len_sd_nad(iyear,iage); 
		cprob_lenvec(ilen)=cumd_norm(zscore_len);
        lenprob_oct15(iyear,iage,ilen)=cprob_lenvec(ilen)-cprob_lenvec(ilen-1);
      }
    //last length bin is a plus group
	//population
    zscore_len=((lenbins(nlenbins)-0.5*lenbins_width)-meanlen_FL_oct15(iyear,iage)) / len_sd_nad(iyear,iage); 
	lenprob_oct15(iyear,iage,nlenbins)=1.0-cumd_norm(zscore_len);
      lenprob_oct15(iyear,iage)=lenprob_oct15(iyear,iage)/(1.0-cprob_lzero);  //renormalize to account for any prob mass below size=0
   }
  //fleet and survey specific length probs
  //lenprob_nad(iyear)=lenprob_oct15(iyear);
  lenprob_nad=lenprob_oct15;
  }
}

void model_parameters::get_weight_at_age_landings(void)
{
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    wholewgt_cR_mt(iyear)=tv_wgt_middle_mt(iyear);
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
  bpr_F0=sum(elem_prod(N_bpr_F0,wgt_spawn_mt));    
}

void model_parameters::get_selectivity(void)
{
  sel_cRn_block1=logistic_double(agebins, selpar_A50_cRn, selpar_slope_cRn, selpar_A502_cRn,selpar_slope2_cRn);
  sel_cRn_block2=logistic_double(agebins, selpar_A50_cRn2, selpar_slope_cRn2, selpar_A502_cRn2,selpar_slope2_cRn2);
  sel_cRn_block3=logistic_double(agebins, selpar_A50_cRn3, selpar_slope_cRn3, selpar_A502_cRn3,selpar_slope2_cRn3);
  sel_cRn_block4=logistic_double(agebins, selpar_A50_cRn4, selpar_slope_cRn4, selpar_A502_cRn4,selpar_slope2_cRn4);
  sel_cRs_block1=logistic_double(agebins, selpar_A50_cRs, selpar_slope_cRs, selpar_A502_cRs,selpar_slope2_cRs);
  sel_cRs_block2=logistic_double(agebins, selpar_A50_cRs2, selpar_slope_cRs2, selpar_A502_cRs2,selpar_slope2_cRs2);
  sel_cRs_block3=logistic_double(agebins, selpar_A50_cRs3, selpar_slope_cRs3, selpar_A502_cRs3,selpar_slope2_cRs3);
  sel_cRs_block4=logistic_double(agebins, selpar_A50_cRs4, selpar_slope_cRs4, selpar_A502_cRs4,selpar_slope2_cRs4);
  sel_cBn_block1=logistic_double(agebins, selpar_A50_cBn, selpar_slope_cBn, selpar_A502_cBn,selpar_slope2_cBn);
  sel_cBn_block2=logistic_double(agebins, selpar_A50_cBn2, selpar_slope_cBn2, selpar_A502_cBn2,selpar_slope2_cBn2);
  sel_cBs_block1=logistic_double(agebins, selpar_A50_cBs, selpar_slope_cBs, selpar_A502_cBs,selpar_slope2_cBs);
  //sel_cBs_block2=logistic_double(agebins, selpar_A50_cBs2, selpar_slope_cBs2, selpar_A502_cBs2,selpar_slope2_cBs2);
  sel_nad_block1=logistic(agebins, selpar_A50_nad, selpar_slope_nad);
  //sel_nad_block1=logistic_double(agebins, selpar_A50_nad, selpar_slope_nad, selpar_A502_nad,selpar_slope2_nad);
  sel_mad_block1=logistic_double(agebins, selpar_A50_mad, selpar_slope_mad, selpar_A502_mad,selpar_slope2_mad);
  //sel_sad_block1=logistic_double(agebins, selpar_A50_sad, selpar_slope_sad, selpar_A502_sad,selpar_slope2_sad);
  //selpar_age0_cRn=1.0/(1.0+mfexp(-sel_age0_cRn_logit));   //block 1
  //selpar_age1_cRn=1.0/(1.0+mfexp(-sel_age1_cRn_logit));
  //selpar_age2_cRn=1.0/(1.0+mfexp(-sel_age2_cRn_logit));
  //selpar_age2_cRn=1.0;
  //selpar_age3_cRn=1.0/(1.0+mfexp(-sel_age3_cRn_logit));
  //selpar_age4_cRn=1.0/(1.0+mfexp(-sel_age4_cRn_logit));
  //selpar_age5_cRn=1.0/(1.0+mfexp(-sel_age5_cRn_logit));
  //selpar_age6_cRn=1.0/(1.0+mfexp(-sel_age6_cRn_logit));
  //sel_age_cRn_vec(1)=selpar_age0_cRn;
  //sel_age_cRn_vec(2)=selpar_age1_cRn;
  //sel_age_cRn_vec(3)=selpar_age2_cRn;
  //sel_age_cRn_vec(4)=selpar_age3_cRn;
  //sel_age_cRn_vec(5)=selpar_age4_cRn;
  //sel_age_cRn_vec(6)=selpar_age5_cRn;
  //sel_age_cRn_vec(7)=selpar_age6_cRn;
  //sel_cRn_block1=sel_age_cRn_vec;
  //selpar_age0_cRn2=1.0/(1.0+mfexp(-sel_age0_cRn2_logit));   //block 2
  //selpar_age1_cRn2=1.0/(1.0+mfexp(-sel_age1_cRn2_logit));
  //selpar_age2_cRn2=1.0/(1.0+mfexp(-sel_age2_cRn2_logit));
  //selpar_age3_cRn2=1.0;
  //selpar_age3_cRn2=1.0/(1.0+mfexp(-sel_age3_cRn2_logit));
  //selpar_age4_cRn2=1.0/(1.0+mfexp(-sel_age4_cRn2_logit));
  //selpar_age5_cRn2=1.0/(1.0+mfexp(-sel_age5_cRn2_logit));
  //selpar_age6_cRn2=1.0/(1.0+mfexp(-sel_age6_cRn2_logit));
  //sel_age_cRn2_vec(1)=selpar_age0_cRn2;
  //sel_age_cRn2_vec(2)=selpar_age1_cRn2;
  //sel_age_cRn2_vec(3)=selpar_age2_cRn2;
  //sel_age_cRn2_vec(4)=selpar_age3_cRn2;
  //sel_age_cRn2_vec(5)=selpar_age4_cRn2;
  //sel_age_cRn2_vec(6)=selpar_age5_cRn2;
  //sel_age_cRn2_vec(7)=selpar_age6_cRn2;
  //sel_cRn_block2=sel_age_cRn2_vec;
  //selpar_age0_cRn3=1.0/(1.0+mfexp(-sel_age0_cRn3_logit));   //block 3
  //selpar_age1_cRn3=1.0/(1.0+mfexp(-sel_age1_cRn3_logit));
  //selpar_age2_cRn3=1.0/(1.0+mfexp(-sel_age2_cRn3_logit));
  //selpar_age3_cRn3=1.0;
  //selpar_age3_cRn3=1.0/(1.0+mfexp(-sel_age3_cRn3_logit));
  //selpar_age4_cRn3=1.0/(1.0+mfexp(-sel_age4_cRn3_logit));
  //selpar_age5_cRn3=1.0/(1.0+mfexp(-sel_age5_cRn3_logit));
  //selpar_age6_cRn3=1.0/(1.0+mfexp(-sel_age6_cRn3_logit));
  //sel_age_cRn3_vec(1)=selpar_age0_cRn3;
  //sel_age_cRn3_vec(2)=selpar_age1_cRn3;
  //sel_age_cRn3_vec(3)=selpar_age2_cRn3;
  //sel_age_cRn3_vec(4)=selpar_age3_cRn3;
  //sel_age_cRn3_vec(5)=selpar_age4_cRn3;
  //sel_age_cRn3_vec(6)=selpar_age5_cRn3;
  //sel_age_cRn3_vec(7)=selpar_age6_cRn3;
  //sel_cRn_block3=sel_age_cRn3_vec;
  //selpar_age0_cRn4=1.0/(1.0+mfexp(-sel_age0_cRn4_logit));   //block 4
  //selpar_age1_cRn4=1.0/(1.0+mfexp(-sel_age1_cRn4_logit));
  //selpar_age2_cRn4=1.0/(1.0+mfexp(-sel_age2_cRn4_logit));
  //selpar_age3_cRn4=1.0;
  //selpar_age3_cRn4=1.0/(1.0+mfexp(-sel_age3_cRn4_logit));
  //selpar_age4_cRn4=1.0/(1.0+mfexp(-sel_age4_cRn4_logit));
  //selpar_age5_cRn4=1.0/(1.0+mfexp(-sel_age5_cRn4_logit));
  //selpar_age6_cRn4=1.0/(1.0+mfexp(-sel_age6_cRn4_logit));
  //sel_age_cRn4_vec(1)=selpar_age0_cRn4;
  //sel_age_cRn4_vec(2)=selpar_age1_cRn4;
  //sel_age_cRn4_vec(3)=selpar_age2_cRn4;
  //sel_age_cRn4_vec(4)=selpar_age3_cRn4;
  //sel_age_cRn4_vec(5)=selpar_age4_cRn4;
  //sel_age_cRn4_vec(6)=selpar_age5_cRn4;
  //sel_age_cRn4_vec(7)=selpar_age6_cRn4;
  //sel_cRn_block4=sel_age_cRn4_vec;
  //selpar_age0_cRs=1.0/(1.0+mfexp(-sel_age0_cRs_logit));     //block 1
  //selpar_age1_cRs=1.0/(1.0+mfexp(-sel_age1_cRs_logit));
  //selpar_age2_cRs=1.0/(1.0+mfexp(-sel_age2_cRs_logit));
  //selpar_age2_cRs=1.0;
  //selpar_age3_cRs=1.0/(1.0+mfexp(-sel_age3_cRs_logit));
  //selpar_age4_cRs=1.0/(1.0+mfexp(-sel_age4_cRs_logit));
  //selpar_age5_cRs=1.0/(1.0+mfexp(-sel_age5_cRs_logit));
  //selpar_age6_cRs=1.0/(1.0+mfexp(-sel_age6_cRs_logit));
  //sel_age_cRs_vec(1)=selpar_age0_cRs;
  //sel_age_cRs_vec(2)=selpar_age1_cRs;
  //sel_age_cRs_vec(3)=selpar_age2_cRs;
  //sel_age_cRs_vec(4)=selpar_age3_cRs;
  //sel_age_cRs_vec(5)=selpar_age4_cRs;
  //sel_age_cRs_vec(6)=selpar_age5_cRs;
  //sel_age_cRs_vec(7)=selpar_age6_cRs;
  //sel_cRs_block1=sel_age_cRs_vec;
  //selpar_age0_cRs2=1.0/(1.0+mfexp(-sel_age0_cRs2_logit));     //block 2
  //selpar_age1_cRs2=1.0/(1.0+mfexp(-sel_age1_cRs2_logit));
  //selpar_age2_cRs2=1.0/(1.0+mfexp(-sel_age2_cRs2_logit));
  //selpar_age2_cRs2=1.0;
  //selpar_age3_cRs2=1.0/(1.0+mfexp(-sel_age3_cRs2_logit));
  //selpar_age4_cRs2=1.0/(1.0+mfexp(-sel_age4_cRs2_logit));
  //selpar_age5_cRs2=1.0/(1.0+mfexp(-sel_age5_cRs2_logit));
  //selpar_age6_cRs2=1.0/(1.0+mfexp(-sel_age6_cRs2_logit));
  //sel_age_cRs2_vec(1)=selpar_age0_cRs2;
  //sel_age_cRs2_vec(2)=selpar_age1_cRs2;
  //sel_age_cRs2_vec(3)=selpar_age2_cRs2;
  //sel_age_cRs2_vec(4)=selpar_age3_cRs2;
  //sel_age_cRs2_vec(5)=selpar_age4_cRs2;
  //sel_age_cRs2_vec(6)=selpar_age5_cRs2;
  //sel_age_cRs2_vec(7)=selpar_age6_cRs2;
  //sel_cRs_block2=sel_age_cRs2_vec;
  //selpar_age0_cRs3=1.0/(1.0+mfexp(-sel_age0_cRs3_logit));     //block 3
  //selpar_age1_cRs3=1.0/(1.0+mfexp(-sel_age1_cRs3_logit));
  //selpar_age2_cRs3=1.0/(1.0+mfexp(-sel_age2_cRs3_logit));
  //selpar_age2_cRs3=1.0;
  //selpar_age3_cRs3=1.0/(1.0+mfexp(-sel_age3_cRs3_logit));
  //selpar_age4_cRs3=1.0/(1.0+mfexp(-sel_age4_cRs3_logit));
  //selpar_age5_cRs3=1.0/(1.0+mfexp(-sel_age5_cRs3_logit));
  //selpar_age6_cRs3=1.0/(1.0+mfexp(-sel_age6_cRs3_logit));
  //sel_age_cRs3_vec(1)=selpar_age0_cRs3;
  //sel_age_cRs3_vec(2)=selpar_age1_cRs3;
  //sel_age_cRs3_vec(3)=selpar_age2_cRs3;
  //sel_age_cRs3_vec(4)=selpar_age3_cRs3;
  //sel_age_cRs3_vec(5)=selpar_age4_cRs3;
  //sel_age_cRs3_vec(6)=selpar_age5_cRs3;
  //sel_age_cRs3_vec(7)=selpar_age6_cRs3;
  //sel_cRs_block3=sel_age_cRs3_vec;
  //selpar_age0_cRs4=1.0/(1.0+mfexp(-sel_age0_cRs4_logit));     //block 4
  //selpar_age1_cRs4=1.0/(1.0+mfexp(-sel_age1_cRs4_logit));
  //selpar_age2_cRs4=1.0/(1.0+mfexp(-sel_age2_cRs4_logit));
  //selpar_age2_cRs4=1.0;
  //selpar_age3_cRs4=1.0/(1.0+mfexp(-sel_age3_cRs4_logit));
  //selpar_age4_cRs4=1.0/(1.0+mfexp(-sel_age4_cRs4_logit));
  //selpar_age5_cRs4=1.0/(1.0+mfexp(-sel_age5_cRs4_logit));
  //selpar_age6_cRs4=1.0/(1.0+mfexp(-sel_age6_cRs4_logit));
  //sel_age_cRs4_vec(1)=selpar_age0_cRs4;
  //sel_age_cRs4_vec(2)=selpar_age1_cRs4;
  //sel_age_cRs4_vec(3)=selpar_age2_cRs4;
  //sel_age_cRs4_vec(4)=selpar_age3_cRs4;
  //sel_age_cRs4_vec(5)=selpar_age4_cRs4;
  //sel_age_cRs4_vec(6)=selpar_age5_cRs4;
  //sel_age_cRs4_vec(7)=selpar_age6_cRs4;
  //sel_cRs_block4=sel_age_cRs4_vec;
  //selpar_age0_cBn=1.0/(1.0+mfexp(-sel_age0_cBn_logit));      //block 1
  //selpar_age1_cBn=1.0/(1.0+mfexp(-sel_age1_cBn_logit));
  //selpar_age2_cBn=1.0/(1.0+mfexp(-sel_age2_cBn_logit));
  //selpar_age3_cBn=1.0;
  //selpar_age3_cBn=1.0/(1.0+mfexp(-sel_age3_cBn_logit));
  //selpar_age4_cBn=1.0/(1.0+mfexp(-sel_age4_cBn_logit));
  //selpar_age5_cBn=1.0/(1.0+mfexp(-sel_age5_cBn_logit));
  //selpar_age6_cBn=1.0/(1.0+mfexp(-sel_age6_cBn_logit));
  //sel_age_cBn_vec(1)=selpar_age0_cBn;
  //sel_age_cBn_vec(2)=selpar_age1_cBn;
  //sel_age_cBn_vec(3)=selpar_age2_cBn;
  //sel_age_cBn_vec(4)=selpar_age3_cBn;
  //sel_age_cBn_vec(5)=selpar_age4_cBn;
  //sel_age_cBn_vec(6)=selpar_age5_cBn;
  //sel_age_cBn_vec(7)=selpar_age6_cBn;
  //sel_cBn_block1=sel_age_cBn_vec;
  //selpar_age0_cBn2=1.0/(1.0+mfexp(-sel_age0_cBn2_logit));      //block 2
  //selpar_age1_cBn2=1.0/(1.0+mfexp(-sel_age1_cBn2_logit));
  //selpar_age2_cBn2=1.0/(1.0+mfexp(-sel_age2_cBn2_logit));
  //selpar_age3_cBn2=1.0;
  //selpar_age3_cBn2=1.0/(1.0+mfexp(-sel_age3_cBn2_logit));
  //selpar_age4_cBn2=1.0/(1.0+mfexp(-sel_age4_cBn2_logit));
  //selpar_age5_cBn2=1.0/(1.0+mfexp(-sel_age5_cBn2_logit));
  //selpar_age6_cBn2=1.0/(1.0+mfexp(-sel_age6_cBn2_logit));
  //sel_age_cBn2_vec(1)=selpar_age0_cBn2;
  //sel_age_cBn2_vec(2)=selpar_age1_cBn2;
  //sel_age_cBn2_vec(3)=selpar_age2_cBn2;
  //sel_age_cBn2_vec(4)=selpar_age3_cBn2;
  //sel_age_cBn2_vec(5)=selpar_age4_cBn2;
  //sel_age_cBn2_vec(6)=selpar_age5_cBn2;
  //sel_age_cBn2_vec(7)=selpar_age6_cBn2;
  //sel_cBn_block2=sel_age_cBn2_vec;
  //selpar_age0_cBs=1.0/(1.0+mfexp(-sel_age0_cBs_logit));       //block 1
  //selpar_age1_cBs=1.0/(1.0+mfexp(-sel_age1_cBs_logit));
  //selpar_age2_cBs=1.0/(1.0+mfexp(-sel_age2_cBs_logit));
  //selpar_age2_cBs=1.0;
  //selpar_age3_cBs=1.0/(1.0+mfexp(-sel_age3_cBs_logit));
  //selpar_age4_cBs=1.0/(1.0+mfexp(-sel_age4_cBs_logit));
  //selpar_age5_cBs=1.0/(1.0+mfexp(-sel_age5_cBs_logit));
  //selpar_age6_cBs=1.0/(1.0+mfexp(-sel_age6_cBs_logit));
  //sel_age_cBs_vec(1)=selpar_age0_cBs;
  //sel_age_cBs_vec(2)=selpar_age1_cBs;
  //sel_age_cBs_vec(3)=selpar_age2_cBs;
  //sel_age_cBs_vec(4)=selpar_age3_cBs;
  //sel_age_cBs_vec(5)=selpar_age4_cBs;
  //sel_age_cBs_vec(6)=selpar_age5_cBs;
  //sel_age_cBs_vec(7)=selpar_age6_cBs;
  //sel_cBs_block1=sel_age_cBs_vec;
  //selpar_age0_cBs2=1.0/(1.0+mfexp(-sel_age0_cBs2_logit));       //block 2
  //selpar_age1_cBs2=1.0/(1.0+mfexp(-sel_age1_cBs2_logit));
  //selpar_age2_cBs2=1.0/(1.0+mfexp(-sel_age2_cBs2_logit));
  //selpar_age2_cBs2=1.0;
  //selpar_age3_cBs2=1.0/(1.0+mfexp(-sel_age3_cBs2_logit));
  //selpar_age4_cBs2=1.0/(1.0+mfexp(-sel_age4_cBs2_logit));
  //selpar_age5_cBs2=1.0/(1.0+mfexp(-sel_age5_cBs2_logit));
  //selpar_age6_cBs2=1.0/(1.0+mfexp(-sel_age6_cBs2_logit));
  //sel_age_cBs2_vec(1)=selpar_age0_cBs2;
  //sel_age_cBs2_vec(2)=selpar_age1_cBs2;
  //sel_age_cBs2_vec(3)=selpar_age2_cBs2;
  //sel_age_cBs2_vec(4)=selpar_age3_cBs2;
  //sel_age_cBs2_vec(5)=selpar_age4_cBs2;
  //sel_age_cBs2_vec(6)=selpar_age5_cBs2;
  //sel_age_cBs2_vec(7)=selpar_age6_cBs2;
  //sel_cBs_block2=sel_age_cBs2_vec;
  //selpar_age0_nad=1.0/(1.0+mfexp(-sel_age0_nad_logit));
  //selpar_age1_nad=1.0/(1.0+mfexp(-sel_age1_nad_logit));
  //selpar_age2_nad=1.0/(1.0+mfexp(-sel_age2_nad_logit));
  //selpar_age2_nad=1.0;
  //selpar_age3_nad=1.0/(1.0+mfexp(-sel_age3_nad_logit));
  //selpar_age4_nad=1.0/(1.0+mfexp(-sel_age4_nad_logit));
  //selpar_age5_nad=1.0/(1.0+mfexp(-sel_age5_nad_logit));
  //selpar_age6_nad=1.0/(1.0+mfexp(-sel_age6_nad_logit));
  //sel_age_nad_vec(1)=selpar_age0_nad;
  //sel_age_nad_vec(2)=selpar_age1_nad;
  //sel_age_nad_vec(3)=selpar_age2_nad;
  //sel_age_nad_vec(4)=selpar_age3_nad;
  //sel_age_nad_vec(5)=selpar_age4_nad;
  //sel_age_nad_vec(6)=selpar_age5_nad;
  //sel_age_nad_vec(7)=selpar_age6_nad;
  //sel_nad_block1=sel_age_nad_vec;
  //selpar_age0_mad=1.0/(1.0+mfexp(-sel_age0_mad_logit));
  //selpar_age1_mad=1.0/(1.0+mfexp(-sel_age1_mad_logit));
  //selpar_age2_mad=1.0/(1.0+mfexp(-sel_age2_mad_logit));
  //selpar_age2_mad=1.0;
  //selpar_age3_mad=1.0/(1.0+mfexp(-sel_age3_mad_logit));
  //selpar_age4_mad=1.0/(1.0+mfexp(-sel_age4_mad_logit));
  //selpar_age5_mad=1.0/(1.0+mfexp(-sel_age5_mad_logit));
  //selpar_age6_mad=1.0/(1.0+mfexp(-sel_age6_mad_logit));
  //sel_age_mad_vec(1)=selpar_age0_nad;
  //sel_age_mad_vec(2)=selpar_age1_mad;
  //sel_age_mad_vec(3)=selpar_age2_mad;
  //sel_age_mad_vec(4)=selpar_age3_mad;
  //sel_age_mad_vec(5)=selpar_age4_mad;
  //sel_age_mad_vec(6)=selpar_age5_mad;
  //sel_age_mad_vec(7)=selpar_age6_mad;
  //sel_mad_block1=sel_age_mad_vec;
  //selpar_age0_sad=1.0/(1.0+mfexp(-sel_age0_sad_logit));
  //selpar_age1_sad=1.0/(1.0+mfexp(-sel_age1_sad_logit));
  //selpar_age2_sad=1.0/(1.0+mfexp(-sel_age2_sad_logit));
  //selpar_age2_sad=1.0;
  //selpar_age3_sad=1.0/(1.0+mfexp(-sel_age3_sad_logit));
  //selpar_age4_sad=1.0/(1.0+mfexp(-sel_age4_sad_logit));
  //selpar_age5_sad=1.0/(1.0+mfexp(-sel_age5_sad_logit));
  //selpar_age6_sad=1.0/(1.0+mfexp(-sel_age6_sad_logit));
  //sel_age_sad_vec(1)=selpar_age0_sad;
  //sel_age_sad_vec(2)=selpar_age1_sad;
  //sel_age_sad_vec(3)=selpar_age2_sad;
  //sel_age_sad_vec(4)=selpar_age3_sad;
  //sel_age_sad_vec(5)=selpar_age4_sad;
  //sel_age_sad_vec(6)=selpar_age5_sad;
  //sel_age_sad_vec(7)=selpar_age6_sad;
  //sel_sad_block1=sel_age_sad_vec;
  //BLOCK 1 for selex - 1955 to 1969  
  for (iyear=styr; iyear<=endyr_selex_phase1; iyear++)
   {     
    sel_cRn(iyear)=sel_cRn_block1;
    sel_cRs(iyear)=sel_cRs_block1;
    sel_cBn(iyear)=sel_cBn_block1;
    sel_cBs(iyear)=sel_cBs_block1;
   }
  //BLOCK 2 for selex - 1970 to 1971
  for (iyear=(endyr_selex_phase1+1); iyear<=endyr_selex_phase2; iyear++)
   {
    //sel_cRn(iyear)=sel_cRn_block1;
    sel_cRn(iyear)=sel_cRn_block2;
    sel_cRs(iyear)=sel_cRs_block1;
    sel_cBn(iyear)=sel_cBn_block1;
    sel_cBs(iyear)=sel_cBs_block1; 
   }
  //BLOCK 3 for selex - 1972 to 1984
  for (iyear=(endyr_selex_phase2+1); iyear<=endyr_selex_phase3; iyear++)
   {
    //sel_cRn(iyear)=sel_cRn_block1;
    sel_cRn(iyear)=sel_cRn_block2;
    //sel_cRs(iyear)=sel_cRs_block1;
    sel_cRs(iyear)=sel_cRs_block2;
    sel_cBn(iyear)=sel_cBn_block1;
    sel_cBs(iyear)=sel_cBs_block1;
   }
  //BLOCK 4 for selex - 1985 to 1989
  for (iyear=(endyr_selex_phase3+1); iyear<=endyr_selex_phase4; iyear++)
   {
    //sel_cRn(iyear)=sel_cRn_block1;
    sel_cRn(iyear)=sel_cRn_block2;
    //sel_cRs(iyear)=sel_cRs_block1;
    sel_cRs(iyear)=sel_cRs_block2;
    sel_cBn(iyear)=sel_cBn_block1;
    sel_cBs(iyear)=sel_cBs_block1;
    sel_mad(iyear)=sel_mad_block1; 
   }
  //BLOCK 5 for selex - 1990 to 1993
  for (iyear=(endyr_selex_phase4+1); iyear<=endyr_selex_phase5; iyear++)
   {   
    //sel_cRn(iyear)=sel_cRn_block1;
    sel_cRn(iyear)=sel_cRn_block2;
    //sel_cRs(iyear)=sel_cRs_block1;
    sel_cRs(iyear)=sel_cRs_block2;
    sel_cBn(iyear)=sel_cBn_block1;
    sel_cBs(iyear)=sel_cBs_block1;
    sel_nad(iyear)=sel_nad_block1;
    sel_mad(iyear)=sel_mad_block1;
    //sel_sad(iyear)=sel_sad_block1; 
   }  
  //BLOCK 6 for selex - 1994 to 2004
  for (iyear=(endyr_selex_phase5+1); iyear<=endyr_selex_phase6; iyear++)
   {   
    //sel_cRn(iyear)=sel_cRn_block1;
    sel_cRn(iyear)=sel_cRn_block3;
    //sel_cRs(iyear)=sel_cRs_block1;
    sel_cRs(iyear)=sel_cRs_block2;
    sel_cBn(iyear)=sel_cBn_block1;
    sel_cBs(iyear)=sel_cBs_block1;
    sel_nad(iyear)=sel_nad_block1;
    sel_mad(iyear)=sel_mad_block1;
    //sel_sad(iyear)=sel_sad_block1; 
   }  
  //BLOCK 7 for selex - 2005 to 2012
  for (iyear=(endyr_selex_phase6+1); iyear<=endyr_selex_phase7; iyear++)
   {   
    //sel_cRn(iyear)=sel_cRn_block1;
    sel_cRn(iyear)=sel_cRn_block3;
    //sel_cRs(iyear)=sel_cRs_block1;
    sel_cRs(iyear)=sel_cRs_block3;
    sel_cBn(iyear)=sel_cBn_block1;
    sel_cBs(iyear)=sel_cBs_block1;
    sel_nad(iyear)=sel_nad_block1;
    sel_mad(iyear)=sel_mad_block1;
    //sel_sad(iyear)=sel_sad_block1; 
   }  
  //BLOCK 8 for selex - 2012 to 2017
  for (iyear=(endyr_selex_phase7+1); iyear<=endyr; iyear++)
   {   
    //sel_cRn(iyear)=sel_cRn_block1;
    //sel_cRn(iyear)=sel_cRn_block3;
    sel_cRn(iyear)=sel_cRn_block4;
    //sel_cRs(iyear)=sel_cRs_block1;
    //sel_cRs(iyear)=sel_cRs_block3;
    sel_cRs(iyear)=sel_cRs_block4;
    //sel_cBn(iyear)=sel_cBn_block1;
    sel_cBn(iyear)=sel_cBn_block2;
    sel_cBs(iyear)=sel_cBs_block1;
    //sel_cBs(iyear)=sel_cBs_block2;
    sel_nad(iyear)=sel_nad_block1;
    sel_mad(iyear)=sel_mad_block1;
    //sel_sad(iyear)=sel_sad_block1; 
   }  
}

void model_parameters::get_mortality(void)
{
  Fsum.initialize();
  Fapex.initialize();
  F.initialize();
  //initialization F is avg from first 3 yrs of observed landings
  log_F_dev_init_cRn=sum(log_dev_F_L_cRn(styr_L_cRn,(styr_L_cRn+2)))/3.0;
  log_F_dev_init_cRs=sum(log_dev_F_L_cRs(styr_L_cRs,(styr_L_cRs+2)))/3.0;
  log_F_dev_init_cBn=sum(log_dev_F_L_cBn(styr_L_cBn,(styr_L_cBn+2)))/3.0;
  log_F_dev_init_cBs=sum(log_dev_F_L_cBs(styr_L_cBs,(styr_L_cBs+2)))/3.0; 
  for (iyear=styr; iyear<=endyr; iyear++) 
  {
    if(iyear>=styr_L_cRn & iyear<=endyr_L_cRn) //spans full time series
    {F_cRn_out(iyear)=mfexp(log_avg_F_L_cRn+log_dev_F_L_cRn(iyear));}     
     F_cRn(iyear)=sel_cRn(iyear)*F_cRn_out(iyear);
     Fsum(iyear)+=F_cRn_out(iyear);
    if(iyear>=styr_L_cRs & iyear<=endyr_L_cRs) //spans full time series
    {F_cRs_out(iyear)=mfexp(log_avg_F_L_cRs+log_dev_F_L_cRs(iyear));}     
     F_cRs(iyear)=sel_cRs(iyear)*F_cRs_out(iyear);
     Fsum(iyear)+=F_cRs_out(iyear);
    if(iyear>=styr_L_cBn & iyear<=endyr_L_cBn) //spans full time series
    {F_cBn_out(iyear)=mfexp(log_avg_F_L_cBn+log_dev_F_L_cBn(iyear));}     
     F_cBn(iyear)=sel_cBn(iyear)*F_cBn_out(iyear);
     Fsum(iyear)+=F_cBn_out(iyear);
    if(iyear>=styr_L_cBs & iyear<=endyr_L_cBs) //spans full time series
    {F_cBs_out(iyear)=mfexp(log_avg_F_L_cBs+log_dev_F_L_cBs(iyear));}     
     F_cBs(iyear)=sel_cBs(iyear)*F_cBs_out(iyear);
     Fsum(iyear)+=F_cBs_out(iyear);  
    //Total F at age
    F(iyear)=F_cRn(iyear);  //first in additive series (NO +=) 
    F(iyear)+=F_cRs(iyear);
    F(iyear)+=F_cBn(iyear);
    F(iyear)+=F_cBs(iyear);
    Fapex(iyear)=max(F(iyear));
    Z(iyear)=M_tv(iyear)+F(iyear);
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
  B0=bpr_F0*R_virgin*1000000;  //virgin biomass
  B0_q_DD=R_virgin*sum(elem_prod(N_bpr_F0(set_q_DD_stage,nages),wgt_spawn_mt(set_q_DD_stage,nages))); 
  F_initial=sel_cRn(styr)*mfexp(log_avg_F_L_cRn+log_F_dev_init_cRn)
            +sel_cRs(styr)*mfexp(log_avg_F_L_cRs+log_F_dev_init_cRs)
            +sel_cBn(styr)*mfexp(log_avg_F_L_cBn+log_F_dev_init_cBn)
            +sel_cBs(styr)*mfexp(log_avg_F_L_cBs+log_F_dev_init_cBs);			
  Z_initial=M_tv(styr)+F_initial;
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
  R1=SR_eq_func(R0, steep, spr_F0, spr_initial, BiasCor, SR_switch);
  if(R1<0.0) {R1=10.0;} //Avoid unrealistically low popn sizes during search algorithm
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
  //N(styr,1)=N_initial_eq(1)*mfexp(log_dev_rec(styr_rec_dev));
  //cout << "N initialization " << N(styr) << endl;
  N_mdyr(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*0.5))); //mid year 
  N_spawn(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*spawn_time_frac))); //peak spawning time 
  N_nad(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*0.63))); //October 15
  N_mad(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*0.13))); //April 15 
  N_sad(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*0.13))); //April 15
  SSB(styr)=sum(elem_prod(N_spawn(styr),reprod_tv(styr)));
  //B_q_DD(styr)=sum(elem_prod(N(styr)(set_q_DD_stage,nages),wgt_spawn_mt(set_q_DD_stage,nages)));
  for (iyear=styr; iyear<endyr; iyear++)
  {
    if(iyear<(styr_rec_dev-1)||iyear>(endyr_rec_dev-1)) //recruitment follows S-R curve (with bias correction) exactly
    {
        //N(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch);
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); //plus group
        //N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod_tv(iyear+1)));
	//B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_spawn_mt(set_q_DD_stage,nages)));
        N(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(iyear+1),SR_switch);
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year        
        N_nad(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.63))); //October 15        
        N_mad(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.13))); //April 15      
        N_sad(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.13))); //April 15        
    }
    else   //recruitment follows S-R curve with lognormal deviation
    {
        //N(iyear+1,1)=SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch)*mfexp(log_dev_rec(iyear+1));
        N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); //plus group
        //N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod_tv(iyear+1)));
        N(iyear+1,1)=SR_func(R0, steep, spr_F0, SSB(iyear+1),SR_switch)*mfexp(log_dev_rec(iyear+1));
        //cout << "R0 " << R0 << endl;
        //cout << "steep " << steep << endl;
        //cout << "SSB " << SSB(iyear+1) << endl;
        //cout << "SR_switch " << SR_switch << endl;
        //cout << "rec dev " << mfexp(log_dev_rec(iyear+1)) << endl;
        //cout << "iyear " << iyear << endl;
        //cout << "N " << N(iyear+1) << endl;
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year
        N_nad(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.63))); //October 15
        N_mad(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.13))); //April 15
        N_sad(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.13))); //April 15        
    }
  }
   //last year (projection) has no recruitment variability
  N(endyr+1)(2,nages)=++elem_prod(N(endyr)(1,nages-1),(mfexp(-1.*Z(endyr)(1,nages-1))));
  N(endyr+1,nages)+=N(endyr,nages)*mfexp(-1.*Z(endyr,nages)); //plus group
  SSB(endyr+1)=sum(elem_prod(N(endyr+1),reprod));
  N(endyr+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB(endyr+1),SR_switch);
   //Time series of interest
  rec=column(N,1);
  SdS0=SSB/S0;
}

void model_parameters::get_landings_numbers(void)
{
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      L_cRn_num(iyear,iage)=N(iyear,iage)*F_cRn(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_cRs_num(iyear,iage)=N(iyear,iage)*F_cRs(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_cBn_num(iyear,iage)=N(iyear,iage)*F_cBn(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_cBs_num(iyear,iage)=N(iyear,iage)*F_cBs(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
    }          
  }
}

void model_parameters::get_landings_wgt(void)
{
  for (iyear=styr; iyear<=endyr; iyear++)
  {    
    L_cRn_mt(iyear)=elem_prod(L_cRn_num(iyear),wholewgt_cR_mt(iyear))*1000000; //in 1000 mt whole weight
    L_cRs_mt(iyear)=elem_prod(L_cRs_num(iyear),wholewgt_cR_mt(iyear))*1000000; //in 1000 mt whole weight
    L_cBn_mt(iyear)=elem_prod(L_cBn_num(iyear),wholewgt_cR_mt(iyear))*1000000; //in 1000 mt whole weight
    L_cBs_mt(iyear)=elem_prod(L_cBs_num(iyear),wholewgt_cR_mt(iyear))*1000000; //in 1000 mt whole weight
    pred_cRn_L_mt(iyear)=sum(L_cRn_mt(iyear));
    pred_cRs_L_mt(iyear)=sum(L_cRs_mt(iyear));
    pred_cBn_L_mt(iyear)=sum(L_cBn_mt(iyear));
    pred_cBs_L_mt(iyear)=sum(L_cBs_mt(iyear));
  }
}

void model_parameters::get_catchability_fcns(void)
{
 //Get rate increase if estimated, otherwise fixed above
  if (set_q_rate_phase>0.0)
  {
      for (iyear=styr_cpue_nad; iyear<=endyr_cpue_nad; iyear++)
      {   if (iyear>styr_cpue_nad & iyear <=2003) 
          {//q_rate_fcn_nad(iyear)=(1.0+q_rate)*q_rate_fcn_nad(iyear-1); //compound
             q_rate_fcn_nad(iyear)=(1.0+(iyear-styr_cpue_nad)*q_rate)*q_rate_fcn_nad(styr_cpue_nad);  //linear
          }
          if (iyear>2003) {q_rate_fcn_nad(iyear)=q_rate_fcn_nad(iyear-1);} 
      }
      for (iyear=styr_cpue_mad; iyear<=endyr_cpue_mad; iyear++)
      {   if (iyear>styr_cpue_mad & iyear <=2003) 
          {//q_rate_fcn_mad(iyear)=(1.0+q_rate)*q_rate_fcn_mad(iyear-1); //compound
             q_rate_fcn_mad(iyear)=(1.0+(iyear-styr_cpue_mad)*q_rate)*q_rate_fcn_mad(styr_cpue_mad);  //linear
          }
          if (iyear>2003) {q_rate_fcn_mad(iyear)=q_rate_fcn_mad(iyear-1);} 
      }
      for (iyear=styr_cpue_sad; iyear<=endyr_cpue_sad; iyear++)
      {   if (iyear>styr_cpue_sad & iyear <=2003) 
          {//q_rate_fcn_sad(iyear)=(1.0+q_rate)*q_rate_fcn_sad(iyear-1); //compound
             q_rate_fcn_sad(iyear)=(1.0+(iyear-styr_cpue_sad)*q_rate)*q_rate_fcn_sad(styr_cpue_sad);  //linear
          }
          if (iyear>2003) {q_rate_fcn_sad(iyear)=q_rate_fcn_sad(iyear-1);} 
      }
      for (iyear=styr_cpue_jai; iyear<=endyr_cpue_jai; iyear++)
      {   if (iyear>styr_cpue_jai & iyear <=2003) 
          {//q_rate_fcn_jai(iyear)=(1.0+q_rate)*q_rate_fcn_jai(iyear-1); //compound
             q_rate_fcn_jai(iyear)=(1.0+(iyear-styr_cpue_jai)*q_rate)*q_rate_fcn_jai(styr_cpue_jai);  //linear
          }
          if (iyear>2003) {q_rate_fcn_jai(iyear)=q_rate_fcn_jai(iyear-1);} 
      }
      //for (iyear=styr_mar_cpue; iyear<=endyr_mar_cpue; iyear++)
      //{   if (iyear>styr_mar_cpue & iyear <=2003) 
      //    {//q_rate_fcn_mar(iyear)=(1.0+q_rate)*q_rate_fcn_mar(iyear-1); //compound
      //       q_rate_fcn_mar(iyear)=(1.0+(iyear-styr_mar_cpue)*q_rate)*q_rate_fcn_mar(styr_mar_cpue);  //linear
      //    }
      //    if (iyear>2003) {q_rate_fcn_mar(iyear)=q_rate_fcn_mar(iyear-1);} 
      //}
      //for (iyear=styr_eco_cpue; iyear<=endyr_eco_cpue; iyear++)
      //{   if (iyear>styr_eco_cpue & iyear <=2003) 
      //    {//q_rate_fcn_eco(iyear)=(1.0+q_rate)*q_rate_fcn_eco(iyear-1); //compound
      //       q_rate_fcn_eco(iyear)=(1.0+(iyear-styr_eco_cpue)*q_rate)*q_rate_fcn_eco(styr_eco_cpue);  //linear
      //    }
      //    if (iyear>2003) {q_rate_fcn_eco(iyear)=q_rate_fcn_eco(iyear-1);} 
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
  //Northern Adult index
  q_nad(styr_cpue_nad)=mfexp(log_q_cpue_nad); 
  for (iyear=styr_cpue_nad; iyear<=endyr_cpue_nad; iyear++)
  {   
      q_nad(iyear)=q_nad(styr_cpue_nad);
      N_nad_cpue(iyear)=elem_prod(N_nad(iyear),sel_nad(iyear));   //matched with October 15
      pred_nad_cpue(iyear)=q_nad(iyear)*sum(N_nad_cpue(iyear));
      //pred_nad_cpue(iyear)=q_nad(iyear)*q_rate_fcn_nad(iyear)*q_DD_fcn(iyear)*sum(N_nad_cpue(iyear));
      //if (iyear<endyr_cpue_nad){q_nad(iyear+1)=q_nad(iyear)*mfexp(q_RW_log_dev_nad(iyear));}
  }
  //Middle Adult index
  q_mad(styr_cpue_mad)=mfexp(log_q_cpue_mad); 
  for (iyear=styr_cpue_mad; iyear<=endyr_cpue_mad; iyear++)
  {   
      q_mad(iyear)=q_mad(styr_cpue_mad);
      N_mad_cpue(iyear)=elem_prod(N_mad(iyear),sel_mad(iyear));   //matched with April 15
      pred_mad_cpue(iyear)=q_mad(iyear)*sum(N_mad_cpue(iyear));
      //pred_mad_cpue(iyear)=q_mad(iyear)*q_rate_fcn_mad(iyear)*q_DD_fcn(iyear)*sum(N_mad_cpue(iyear));
      //if (iyear<endyr_cpue_mad){q_mad(iyear+1)=q_mad(iyear)*mfexp(q_RW_log_dev_mad(iyear));}
  }
  //Southern Adult index
  q_sad(styr_cpue_sad)=mfexp(log_q_cpue_sad); 
  for (iyear=styr_cpue_sad; iyear<=endyr_cpue_sad; iyear++)
  {   
      q_sad(iyear)=q_sad(styr_cpue_sad);
      N_sad_cpue(iyear)=N_sad(iyear,2);   //matched with April 15
      pred_sad_cpue(iyear)=q_sad(iyear)*N_sad_cpue(iyear);
      //pred_sad_cpue(iyear)=q_sad(iyear)*q_rate_fcn_sad(iyear)*q_DD_fcn(iyear)*sum(N_sad_cpue(iyear));
      //if (iyear<endyr_cpue_sad){q_sad(iyear+1)=q_sad(iyear)*mfexp(q_RW_log_dev_sad(iyear));}
  }
 //Juvenile abundance index
  q_jai(styr_cpue_jai)=mfexp(log_q_cpue_jai);
  q2_jai(styr_cpue_jai)=mfexp(log_q2_jai);
  for (iyear=styr_cpue_jai; iyear<=endyr_cpue_jai; iyear++)
  {   
      q_jai(iyear)=q_jai(styr_cpue_jai);
      N_jai_cpue(iyear)=N(iyear,1)*mfexp(-1.*(Z(iyear)(1)*0.25));//matching seine index with June 1  
      pred_jai_cpue(iyear)=q_jai(iyear)*N_jai_cpue(iyear);
        if(iyear>yr_q_change)
        {
          q2_jai(iyear)=q2_jai(styr_cpue_jai);
          pred_jai_cpue(iyear)=q2_jai(iyear)*N_jai_cpue(iyear);
        }
      //pred_jai_cpue(iyear)=q_jai(iyear)*q_rate_fcn_jai(iyear)*q_DD_fcn(iyear)*N_jai_cpue(iyear);
      //if (iyear<endyr_cpue_jai){q_jai(iyear+1)=q_jai(iyear)*mfexp(q_RW_log_dev_jai(iyear));}
  }
  //MARMAP and ECOMON
  q_mar(1)=mfexp(log_q_cpue_mar);
  //cout << "q_mar" << q_mar << endl;
  q_eco(9)=mfexp(log_q_cpue_eco);
  //cout << "q_eco" << q_eco << endl;
  for (iyear=1; iyear<=8; iyear++)
  {   
      q_mar(iyear)=q_mar(1);
      //cout << "q_mar" << q_mar << endl;
      SSB_mareco_cpue(iyear)=SSB(1980+iyear);   //matched with SSB
      //cout << "SSB_mareco_cpue" << SSB_mareco_cpue << endl;
      pred_mareco_cpue(iyear)=q_mar(iyear)*SSB_mareco_cpue(iyear);
      //cout << "pred_mareco_cpue" << pred_mareco_cpue << endl;     
  }
  for (iyear=9; iyear<=nyr_cpue_mareco; iyear++)
  {
      q_eco(iyear)=q_eco(9);
      //cout << "q_eco" << q_eco << endl;
      SSB_mareco_cpue(iyear)=SSB(1991+iyear);
      //cout << "SSB_mareco_cpue2" << SSB_mareco_cpue << endl;
      pred_mareco_cpue(iyear)=q_eco(iyear)*SSB_mareco_cpue(iyear);
      //cout << "pred_mareco_cpue2" << pred_mareco_cpue << endl;   
  } 
}

void model_parameters::get_length_comps(void)
{
 //Northern adult index
  for (iyear=1;iyear<=nyr_lenc_nad;iyear++)
  {
    pred_nad_lenc(iyear)=(N_nad_cpue(yrs_lenc_nad(iyear))
                         *lenprob_nad(yrs_lenc_nad(iyear)))
                          /sum(N_nad_cpue(yrs_lenc_nad(iyear)));
  }
  //Middle adult index
  for (iyear=1;iyear<=nyr_lenc_mad;iyear++)
  {
    pred_mad_lenc(iyear)=(N_mad_cpue(yrs_lenc_mad(iyear))
                         *lenprob_mad(yrs_lenc_mad(iyear)))
                          /sum(N_mad_cpue(yrs_lenc_mad(iyear)));
  }
  //Southern adult index
  //for (iyear=1;iyear<=nyr_sad_lenc;iyear++)
  //{
  //  pred_sad_lenc(iyear)=(N_sad_cpue(yrs_sad_lenc(iyear))
  //                       *lenprob_sad(yrs_sad_lenc(iyear)))
  //                        /sum(N_sad_cpue(yrs_sad_lenc(iyear)));
  //}
}

void model_parameters::get_age_comps(void)
{
  //Commerical reduction - north
  for (iyear=1;iyear<=nyr_agec_cRn;iyear++) 
  {
    ErrorFree_cRn_agec(iyear)=L_cRn_num(yrs_agec_cRn(iyear))/sum(L_cRn_num(yrs_agec_cRn(iyear)));  
    pred_cRn_agec_allages(iyear)=age_error*(ErrorFree_cRn_agec(iyear)/sum(ErrorFree_cRn_agec(iyear)));   
    for (iage=1; iage<=nages_agec; iage++) {pred_cRn_agec(iyear,iage)=pred_cRn_agec_allages(iyear,iage);} 
    //for (iage=(nages_agec+1); iage<=nages; iage++) {pred_cRn_agec(iyear,nages_agec)+=pred_cRn_agec_allages(iyear,iage);} //plus group                             
  }
  //Commerical reduction - south
  for (iyear=1;iyear<=nyr_agec_cRs;iyear++) 
  {
    ErrorFree_cRs_agec(iyear)=L_cRs_num(yrs_agec_cRs(iyear))/sum(L_cRs_num(yrs_agec_cRs(iyear)));  
    pred_cRs_agec_allages(iyear)=age_error*(ErrorFree_cRs_agec(iyear)/sum(ErrorFree_cRs_agec(iyear)));   
    for (iage=1; iage<=nages_agec; iage++) {pred_cRs_agec(iyear,iage)=pred_cRs_agec_allages(iyear,iage);} 
    //for (iage=(nages_agec+1); iage<=nages; iage++) {pred_cRs_agec(iyear,nages_agec)+=pred_cRs_agec_allages(iyear,iage);} //plus group                             
  }
  //Commerical bait - north
  for (iyear=1;iyear<=nyr_agec_cBn;iyear++) 
  {
    ErrorFree_cBn_agec(iyear)=L_cBn_num(yrs_agec_cBn(iyear))/sum(L_cBn_num(yrs_agec_cBn(iyear)));  
    pred_cBn_agec_allages(iyear)=age_error*(ErrorFree_cBn_agec(iyear)/sum(ErrorFree_cBn_agec(iyear)));   
    for (iage=1; iage<=nages_agec; iage++) {pred_cBn_agec(iyear,iage)=pred_cBn_agec_allages(iyear,iage);} 
    //for (iage=(nages_agec+1); iage<=nages; iage++) {pred_cBn_agec(iyear,nages_agec)+=pred_cBn_agec_allages(iyear,iage);} //plus group                             
  }
  //Commerical bait - south
  for (iyear=1;iyear<=nyr_agec_cBs;iyear++) 
  {
    ErrorFree_cBs_agec(iyear)=L_cBs_num(yrs_agec_cBs(iyear))/sum(L_cBs_num(yrs_agec_cBs(iyear)));  
    pred_cBs_agec_allages(iyear)=age_error*(ErrorFree_cBs_agec(iyear)/sum(ErrorFree_cBs_agec(iyear)));   
    for (iage=1; iage<=nages_agec; iage++) {pred_cBs_agec(iyear,iage)=pred_cBs_agec_allages(iyear,iage);} 
    //for (iage=(nages_agec+1); iage<=nages; iage++) {pred_cBs_agec(iyear,nages_agec)+=pred_cBs_agec_allages(iyear,iage);} //plus group                             
  }
}

void model_parameters::get_weighted_current(void)
{
  F_temp_sum=0.0;
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cRn+
        sum(log_dev_F_L_cRn((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cRs+
        sum(log_dev_F_L_cRs((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cBn+
        sum(log_dev_F_L_cBn((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cBs+
        sum(log_dev_F_L_cBs((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
  F_cRn_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cRn+
        sum(log_dev_F_L_cRn((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_cRs_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cRs+
        sum(log_dev_F_L_cRs((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_cBn_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cBn+
        sum(log_dev_F_L_cBn((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_cBs_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cBs+
        sum(log_dev_F_L_cBs((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  log_F_dev_end_cRn=sum(log_dev_F_L_cRn((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_cRs=sum(log_dev_F_L_cRs((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_cBn=sum(log_dev_F_L_cBn((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_cBs=sum(log_dev_F_L_cBs((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  F_end_L=sel_cRn(endyr)*mfexp(log_avg_F_L_cRn+log_F_dev_end_cRn)+
          sel_cRs(endyr)*mfexp(log_avg_F_L_cRs+log_F_dev_end_cRs)+
          sel_cBn(endyr)*mfexp(log_avg_F_L_cBn+log_F_dev_end_cBn)+
          sel_cBs(endyr)*mfexp(log_avg_F_L_cBs+log_F_dev_end_cBs);    
  F_end=F_end_L;
  F_end_apex=max(F_end);
  sel_wgted_tot=F_end/F_end_apex;
  sel_wgted_L=elem_prod(sel_wgted_tot, elem_div(F_end_L,F_end));
  wgt_wgted_L_denom=F_cRn_prop+F_cRs_prop+F_cBn_prop+F_cBs_prop;  
  wgt_wgted_L_mt=F_cRn_prop/wgt_wgted_L_denom*wholewgt_cR_mt(endyr)*1000+ //to scale to 1000s mt
                 F_cRs_prop/wgt_wgted_L_denom*wholewgt_cR_mt(endyr)*1000+ //to scale to 1000s mt
                 F_cBn_prop/wgt_wgted_L_denom*wholewgt_cR_mt(endyr)*1000+ //to scale to 1000s mt
                 F_cBs_prop/wgt_wgted_L_denom*wholewgt_cR_mt(endyr)*1000; //to scale to 1000s mt 
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
    R_eq(ff)=SR_eq_func(R0, steep, spr_msy(1), spr_msy(ff), BiasCor, SR_switch);
    if (R_eq(ff)<dzero) {R_eq(ff)=dzero;}    
    N_age_msy*=R_eq(ff);
    N_age_msy_spawn*=R_eq(ff);
    for (iage=1; iage<=nages; iage++)
    {
      L_age_msy(iage)=N_age_msy(iage)*(F_L_age_msy(iage)/Z_age_msy(iage))*
                      (1.-mfexp(-1.*Z_age_msy(iage)));                     
    }
    SSB_eq(ff)=sum(elem_prod(N_age_msy_spawn,reprod));
    B_eq(ff)=sum(elem_prod(N_age_msy,wgt_spawn_mt))*1000000; //to scale to 1000s mt and catch in 1000s
    L_eq_mt(ff)=sum(elem_prod(L_age_msy,wgt_wgted_L_mt))*1000; //to scale to catch in 1000s, wgt_wgted_L_mt is already scaled to 1000s mt
    //L_eq_knum(ff)=sum(L_age_msy)/1000.0;  
  }  
  msy_mt_out=max(L_eq_mt); //msy in whole weight 
  for(ff=1; ff<=n_iter_msy; ff++)
  {
   if(L_eq_mt(ff) == msy_mt_out) 
      {    
        SSB_msy_out=SSB_eq(ff);
        B_msy_out=B_eq(ff);
        R_msy_out=R_eq(ff);
        //msy_knum_out=L_eq_knum(ff); 
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
                          (mfexp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac))))
                          /(1.0-mfexp(-1.*Z_age_spr(nages)));
    spr_spr(ff)=sum(elem_prod(N_age_spr_spawn,reprod));
    L_spr(ff)=0.0;
    for (iage=1; iage<=nages; iage++)
    {
      L_age_spr(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));
      L_spr(ff)+=L_age_spr(iage)*wgt_wgted_L_mt(iage)*1000.0; //already scaled to 1000s mt, but need to scale to 1000s fish
    }   
  }
  spr_ratio=spr_spr/spr_F0;
  F35_dum=min(fabs(spr_ratio-0.35));
  F30_dum=min(fabs(spr_ratio-0.3));
  F40_dum=min(fabs(spr_ratio-0.4));
  for(ff=1; ff<=n_iter_spr; ff++)
  {   
      if (fabs(spr_ratio(ff)-0.35)==F35_dum) {F35_out=F_spr(ff);}	  
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
  B_F30_out=sum(elem_prod(N_age_spr,wgt_spawn))*1000000;
  L_F30_mt_out=sum(elem_prod(L_age_F30,wgt_wgted_L_mt))*1000; //in whole weight
  L_F30_knum_out=sum(L_age_F30)/1000.0;  
}

void model_parameters::get_miscellaneous_stuff(void)
{
  if(var_rec_dev>0.0)
   {sigma_rec_dev=sqrt(var_rec_dev);} //sample SD of predicted residuals (may not equal rec_sigma)  
   else{sigma_rec_dev=0.0;}
  //len_cv=elem_div(len_sd,meanlen_FL);
  for (iyear=styr; iyear<=endyr; iyear++)
  {
   len_cv_apr15=mean(elem_div(len_sd_mad(iyear),meanlen_FL_apr15(iyear)));
   //len_cv_jun1=mean(elem_div(len_sd_sad(iyear),meanlen_FL_jun1(iyear)));
   len_cv_oct15=mean(elem_div(len_sd_nad(iyear),meanlen_FL_oct15(iyear)));
  }
  //len_cv=(len_cv_apr15+len_cv_jun1+len_cv_oct15)/3;
  len_cv_mad=len_cv_apr15;
  //len_cv_sad=len_cv_jun1;
  len_cv_nad=len_cv_oct15;
  //compute total landings-at-age in 1000 fish and 1000s mt whole weight
  L_total_num.initialize();
  L_total_mt.initialize();
  //L_total_knum_yr.initialize();
  L_total_mt_yr.initialize();  
  for(iyear=styr; iyear<=endyr; iyear++)
  {
        L_total_mt_yr(iyear)=pred_cRn_L_mt(iyear)+pred_cRs_L_mt(iyear)+
                             pred_cBn_L_mt(iyear)+pred_cBs_L_mt(iyear);
        //L_total_knum_yr(iyear)=pred_cR_L_knum(iyear);
        B(iyear)=elem_prod(N(iyear),tv_wgt_spawn_mt(iyear))*1000000;  //scale to 1000s mt and 1000s fish landed
        B_mdyr(iyear)=elem_prod(N_mdyr(iyear),tv_wgt_middle_mt(iyear))*1000000;
        totN(iyear)=sum(N(iyear));   //in 1000s of fish
        totB(iyear)=sum(B(iyear));   //in 1000s of mt
        SSBatage(iyear)=elem_prod(N(iyear),reprod_tv(iyear));
  }
  L_total_num=L_cRn_num+L_cRs_num+L_cBn_num+L_cBs_num;   //landings at age in 1000s fish
  L_total_mt=L_cRn_mt+L_cRs_mt+L_cBn_mt+L_cBs_mt;   //landings at age in 1000s  mt whole weight
  //Time series of interest
  B(endyr+1)=elem_prod(N(endyr+1),wgt_spawn_mt)*1000000;  //scale to 1000s mt and 1000s fish
  //B_mdyr(endyr+1)=elem_prod(N_mdyr(endyr+1),wgt_fish_mt)*1000000;  //scale to 1000s mt and 1000s fish
  totN(endyr+1)=sum(N(endyr+1));  //in 1000s of fish
  totB(endyr+1)=sum(B(endyr+1));  //in 1000s of mt
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
      SdSSB_F30_end=SdSSB_F30(endyr);
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
        Z_proj(iyear)=M+FL_age_proj;
        N_spawn_proj(iyear)(1,nages)=elem_prod(N_proj(iyear)(1,nages),(mfexp(-1.*(Z_proj(iyear)(1,nages))*spawn_time_frac))); //peak spawning time
		SSB_proj(iyear)= sum(elem_prod(N_spawn_proj(iyear),reprod));
        B_proj(iyear)=sum(elem_prod(N_proj(iyear),wgt_spawn_mt))*1000000; //uses spawning weight
		for (iage=1; iage<=nages; iage++)
			{
                          L_age_proj(iyear,iage)=N_proj(iyear,iage)*FL_age_proj(iage)*(1.-mfexp(-1.*Z_proj(iyear,iage)))/Z_proj(iyear,iage);
			}          
        L_knum_proj(iyear)=sum(L_age_proj(iyear))/1000.0;
	L_mt_proj(iyear)=sum(elem_prod(L_age_proj(iyear),wgt_wgted_L_mt));     //in 1000 mt
		if (iyear<endyr_proj) {
			N_proj(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB_proj(iyear),SR_switch);
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
  f_nad_cpue=0.0;
  f_nad_cpue=lk_lognormal(pred_nad_cpue, obs_cpue_nad, obs_cv_cpue_nad, w_cpue_nad);
  fval+=f_nad_cpue;
  fval_data+=f_nad_cpue;
  f_mad_cpue=0.0;
  f_mad_cpue=lk_lognormal(pred_mad_cpue, obs_cpue_mad, obs_cv_cpue_mad, w_cpue_mad);
  fval+=f_mad_cpue;
  fval_data+=f_mad_cpue;
  f_sad_cpue=0.0;
  f_sad_cpue=lk_lognormal(pred_sad_cpue, obs_cpue_sad, obs_cv_cpue_sad, w_cpue_sad);
  fval+=f_sad_cpue;
  fval_data+=f_sad_cpue;  
  f_jai_cpue=0.0;
  f_jai_cpue=lk_lognormal(pred_jai_cpue, obs_cpue_jai, obs_cv_cpue_jai, w_cpue_jai);
  fval+=f_jai_cpue;
  fval_data+=f_jai_cpue;  
  f_mareco_cpue=0.0;
  f_mareco_cpue=lk_lognormal(pred_mareco_cpue, obs_cpue_mareco, obs_cv_cpue_mareco, w_cpue_mareco);
  fval+=f_mareco_cpue;
  fval_data+=f_mareco_cpue;  
  //---Landings-------------------------------
  //f_cRn_L in 1000 mt whole wgt  
  f_cRn_L=lk_lognormal(pred_cRn_L_mt(styr_L_cRn,endyr_L_cRn), obs_L_cRn(styr_L_cRn,endyr_L_cRn),
                      obs_cv_L_cRn(styr_L_cRn,endyr_L_cRn), w_L);
  fval+=f_cRn_L;
  fval_data+=f_cRn_L;
  //f_cRs_L in 1000 mt whole wgt  
  f_cRs_L=lk_lognormal(pred_cRs_L_mt(styr_L_cRs,endyr_L_cRs), obs_L_cRs(styr_L_cRs,endyr_L_cRs),
                      obs_cv_L_cRs(styr_L_cRs,endyr_L_cRs), w_L);
  fval+=f_cRs_L;
  fval_data+=f_cRs_L;
  //f_cBn_L in 1000 mt whole wgt  
  f_cBn_L=lk_lognormal(pred_cBn_L_mt(styr_L_cBn,endyr_L_cBn), obs_L_cBn(styr_L_cBn,endyr_L_cBn),
                      obs_cv_L_cBn(styr_L_cBn,endyr_L_cBn), w_L);
  fval+=f_cBn_L;
  fval_data+=f_cBn_L;
  //f_cBs_L in 1000 mt whole wgt  
  f_cBs_L=lk_lognormal(pred_cBs_L_mt(styr_L_cBs,endyr_L_cBs), obs_L_cBs(styr_L_cBs,endyr_L_cBs),
                      obs_cv_L_cBs(styr_L_cBs,endyr_L_cBs), w_L);
  fval+=f_cBs_L;
  fval_data+=f_cBs_L;
  //f_nad_lenc
  //f_nad_lenc=lk_robust_multinomial(nsamp_lenc_nad, pred_nad_lenc, obs_lenc_nad, nyr_lenc_nad, double(nlenbins), minSS_lenc_nad, w_lenc_nad);
  //f_nad_lenc=lk_logistic_normal(nsamp_lenc_nad, pred_nad_lenc, obs_lenc_nad, nyr_lenc_nad, double(nlenbins), minSS_lenc_nad);
  f_nad_lenc=lk_dirichlet_multinomial(nsamp_lenc_nad, pred_nad_lenc, obs_lenc_nad, nyr_lenc_nad, double(nlenbins), minSS_lenc_nad, log_dm_lenc_nad);
  fval+=f_nad_lenc;
  fval_data+=f_nad_lenc;
  //f_mad_lenc
  //f_mad_lenc=lk_robust_multinomial(nsamp_lenc_mad, pred_mad_lenc, obs_lenc_mad, nyr_lenc_mad, double(nlenbins), minSS_lenc_mad, w_lenc_mad);
  //f_mad_lenc=lk_logistic_normal(nsamp_lenc_mad, pred_mad_lenc, obs_lenc_mad, nyr_lenc_mad, double(nlenbins), minSS_lenc_mad);
  f_mad_lenc=lk_dirichlet_multinomial(nsamp_lenc_mad, pred_mad_lenc, obs_lenc_mad, nyr_lenc_mad, double(nlenbins), minSS_lenc_mad, log_dm_lenc_mad);
  fval+=f_mad_lenc;
  fval_data+=f_mad_lenc;
  //f_sad_lenc
  //f_sad_lenc=lk_robust_multinomial(nsamp_sad_lenc, pred_sad_lenc, obs_sad_lenc, nyr_sad_lenc, double(nlenbins), minSS_sad_lenc, w_lc_sad);
  //f_sad_lenc=lk_logistic_normal(nsamp_sad_lenc, pred_sad_lenc, obs_sad_lenc, nyr_sad_lenc, double(nlenbins), minSS_sad_lenc);
  //f_sad_lenc=lk_dirichlet_multinomial(nsamp_sad_lenc, pred_sad_lenc, obs_sad_lenc, nyr_sad_lenc, double(nlenbins), minSS_sad_lenc, log_dm_sad_lc);
  //fval+=f_sad_lenc;
  //fval_data+=f_sad_lenc;
  //f_cRn_agec
  //f_cRn_agec=lk_robust_multinomial(nsamp_agec_cRn, pred_cRn_agec, obs_agec_cRn, nyr_agec_cRn, double(nages_agec), minSS_agec_cRn, w_agec_cRn);
  //f_cRn_agec=lk_logistic_normal(nsamp_agec_cRn, pred_cRn_agec, obs_agec_cRn, nyr_agec_cRn, double(nages_agec), minSS_agec_cRn);
  f_cRn_agec=lk_dirichlet_multinomial(nsamp_agec_cRn, pred_cRn_agec, obs_agec_cRn, nyr_agec_cRn, double(nages_agec), minSS_agec_cRn, log_dm_agec_cRn);
  fval+=f_cRn_agec;
  fval_data+=f_cRn_agec;
  //f_cRs_agec
  //f_cRs_agec=lk_robust_multinomial(nsamp_agec_cRs, pred_cRs_agec, obs_agec_cRs, nyr_agec_cRs, double(nages_agec), minSS_agec_cRs, w_agec_cRs);
  //f_cRs_agec=lk_logistic_normal(nsamp_agec_cRs, pred_cRs_agec, obs_agec_cRs, nyr_agec_cRs, double(nages_agec), minSS_agec_cRs);
  f_cRs_agec=lk_dirichlet_multinomial(nsamp_agec_cRs, pred_cRs_agec, obs_agec_cRs, nyr_agec_cRs, double(nages_agec), minSS_agec_cRs, log_dm_agec_cRs);
  fval+=f_cRs_agec;
  fval_data+=f_cRs_agec;
  //f_cBn_agec
  //f_cBn_agec=lk_robust_multinomial(nsamp_agec_cBn, pred_cBn_agec, obs_agec_cBn, nyr_agec_cBn, double(nages_agec), minSS_agec_cBn, w_agec_cBn);
  //f_cBn_agec=lk_logistic_normal(nsamp_agec_cBn, pred_cBn_agec, obs_agec_cBn, nyr_agec_cBn, double(nages_agec), minSS_agec_cBn);
  f_cBn_agec=lk_dirichlet_multinomial(nsamp_agec_cBn, pred_cBn_agec, obs_agec_cBn, nyr_agec_cBn, double(nages_agec), minSS_agec_cBn, log_dm_agec_cBn);
  fval+=f_cBn_agec;
  fval_data+=f_cBn_agec;
  //f_cBs_agec
  //f_cBs_agec=lk_robust_multinomial(nsamp_agec_cBs, pred_cBs_agec, obs_agec_cBs, nyr_agec_cBs, double(nages_agec), minSS_agec_cBs, w_agec_cBs);
  //f_cBs_agec=lk_logistic_normal(nsamp_agec_cBs, pred_cBs_agec, obs_agec_cBs, nyr_agec_cBs, double(nages_agec), minSS_agec_cBs);
  f_cBs_agec=lk_dirichlet_multinomial(nsamp_agec_cBs, pred_cBs_agec, obs_agec_cBs, nyr_agec_cBs, double(nages_agec), minSS_agec_cBs, log_dm_agec_cBs);
  fval+=f_cBs_agec;
  fval_data+=f_cBs_agec;
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
 //Random walk components of fishery dependent indices - we don't have any fishery dependent indices
 //f_cHL_RWq_cpue=0.0;
 //for (iyear=styr_cHL_cpue; iyear<endyr_cHL_cpue; iyear++)
 //    {f_cHL_RWq_cpue+=square(q_RW_log_dev_cHL(iyear))/(2.0*set_RWq_var);}
 //fval+=f_cHL_RWq_cpue;   
  f_priors=0.0; 
  f_priors+=neg_log_prior(len_cv_nad_val,set_len_cv_nad(5),set_len_cv_nad(6),set_len_cv_nad(7));
  f_priors+=neg_log_prior(len_cv_mad_val,set_len_cv_mad(5),set_len_cv_mad(6),set_len_cv_mad(7));
  //f_priors+=neg_log_prior(len_cv_sad_val,set_len_cv_sad(5),set_len_cv_sad(6),set_len_cv_sad(7));
  f_priors+=neg_log_prior(steep,set_steep(5),set_log_R0(6),set_log_R0(7)); 
  f_priors+=neg_log_prior(log_R0,set_log_R0(5),set_log_R0(6),set_log_R0(7)); 
  f_priors+=neg_log_prior(R_autocorr,set_R_autocorr(5),set_R_autocorr(6),set_R_autocorr(7));
  f_priors+=neg_log_prior(rec_sigma,set_rec_sigma(5),set_rec_sigma(6),set_rec_sigma(7));
  f_priors+=neg_log_prior(selpar_A50_cRn,set_selpar_A50_cRn(5), set_selpar_A50_cRn(6), set_selpar_A50_cRn(7));
  f_priors+=neg_log_prior(selpar_slope_cRn,set_selpar_slope_cRn(5), set_selpar_slope_cRn(6), set_selpar_slope_cRn(7));
  f_priors+=neg_log_prior(selpar_A502_cRn,set_selpar_A502_cRn(5), set_selpar_A502_cRn(6), set_selpar_A502_cRn(7));
  f_priors+=neg_log_prior(selpar_slope2_cRn,set_selpar_slope2_cRn(5), set_selpar_slope2_cRn(6), set_selpar_slope2_cRn(7));
  f_priors+=neg_log_prior(selpar_A50_cRn2,set_selpar_A50_cRn2(5), set_selpar_A50_cRn2(6), set_selpar_A50_cRn2(7));
  f_priors+=neg_log_prior(selpar_slope_cRn2,set_selpar_slope_cRn2(5), set_selpar_slope_cRn2(6), set_selpar_slope_cRn2(7));
  f_priors+=neg_log_prior(selpar_A502_cRn2,set_selpar_A502_cRn2(5), set_selpar_A502_cRn2(6), set_selpar_A502_cRn2(7));
  f_priors+=neg_log_prior(selpar_slope2_cRn2,set_selpar_slope2_cRn2(5), set_selpar_slope2_cRn2(6), set_selpar_slope2_cRn2(7));
  f_priors+=neg_log_prior(selpar_A50_cRn3,set_selpar_A50_cRn3(5), set_selpar_A50_cRn3(6), set_selpar_A50_cRn3(7));
  f_priors+=neg_log_prior(selpar_slope_cRn3,set_selpar_slope_cRn3(5), set_selpar_slope_cRn3(6), set_selpar_slope_cRn3(7));
  f_priors+=neg_log_prior(selpar_A502_cRn3,set_selpar_A502_cRn3(5), set_selpar_A502_cRn3(6), set_selpar_A502_cRn3(7));
  f_priors+=neg_log_prior(selpar_slope2_cRn3,set_selpar_slope2_cRn3(5), set_selpar_slope2_cRn3(6), set_selpar_slope2_cRn3(7));
  f_priors+=neg_log_prior(selpar_A50_cRn4,set_selpar_A50_cRn4(5), set_selpar_A50_cRn4(6), set_selpar_A50_cRn4(7));
  f_priors+=neg_log_prior(selpar_slope_cRn4,set_selpar_slope_cRn4(5), set_selpar_slope_cRn4(6), set_selpar_slope_cRn4(7));
  f_priors+=neg_log_prior(selpar_A502_cRn4,set_selpar_A502_cRn4(5), set_selpar_A502_cRn4(6), set_selpar_A502_cRn4(7));
  f_priors+=neg_log_prior(selpar_slope2_cRn4,set_selpar_slope2_cRn4(5), set_selpar_slope2_cRn4(6), set_selpar_slope2_cRn4(7));
  f_priors+=neg_log_prior(selpar_A50_cRs,set_selpar_A50_cRs(5), set_selpar_A50_cRs(6), set_selpar_A50_cRs(7));
  f_priors+=neg_log_prior(selpar_slope_cRs,set_selpar_slope_cRs(5), set_selpar_slope_cRs(6), set_selpar_slope_cRs(7));
  f_priors+=neg_log_prior(selpar_A502_cRs,set_selpar_A502_cRs(5), set_selpar_A502_cRs(6), set_selpar_A502_cRs(7));
  f_priors+=neg_log_prior(selpar_slope2_cRs,set_selpar_slope2_cRs(5), set_selpar_slope2_cRs(6), set_selpar_slope2_cRs(7));
  f_priors+=neg_log_prior(selpar_A50_cRs2,set_selpar_A50_cRs2(5), set_selpar_A50_cRs2(6), set_selpar_A50_cRs2(7));
  f_priors+=neg_log_prior(selpar_slope_cRs2,set_selpar_slope_cRs2(5), set_selpar_slope_cRs2(6), set_selpar_slope_cRs2(7));
  f_priors+=neg_log_prior(selpar_A502_cRs2,set_selpar_A502_cRs2(5), set_selpar_A502_cRs2(6), set_selpar_A502_cRs2(7));
  f_priors+=neg_log_prior(selpar_slope2_cRs2,set_selpar_slope2_cRs2(5), set_selpar_slope2_cRs2(6), set_selpar_slope2_cRs2(7));
  f_priors+=neg_log_prior(selpar_A50_cRs3,set_selpar_A50_cRs3(5), set_selpar_A50_cRs3(6), set_selpar_A50_cRs3(7));
  f_priors+=neg_log_prior(selpar_slope_cRs3,set_selpar_slope_cRs3(5), set_selpar_slope_cRs3(6), set_selpar_slope_cRs3(7));
  f_priors+=neg_log_prior(selpar_A502_cRs3,set_selpar_A502_cRs3(5), set_selpar_A502_cRs3(6), set_selpar_A502_cRs3(7));
  f_priors+=neg_log_prior(selpar_slope2_cRs3,set_selpar_slope2_cRs3(5), set_selpar_slope2_cRs3(6), set_selpar_slope2_cRs3(7));
  f_priors+=neg_log_prior(selpar_A50_cRs4,set_selpar_A50_cRs4(5), set_selpar_A50_cRs4(6), set_selpar_A50_cRs4(7));
  f_priors+=neg_log_prior(selpar_slope_cRs4,set_selpar_slope_cRs4(5), set_selpar_slope_cRs4(6), set_selpar_slope_cRs4(7));
  f_priors+=neg_log_prior(selpar_A502_cRs4,set_selpar_A502_cRs4(5), set_selpar_A502_cRs4(6), set_selpar_A502_cRs4(7));
  f_priors+=neg_log_prior(selpar_slope2_cRs4,set_selpar_slope2_cRs4(5), set_selpar_slope2_cRs4(6), set_selpar_slope2_cRs4(7));
  f_priors+=neg_log_prior(selpar_A50_cBn,set_selpar_A50_cBn(5), set_selpar_A50_cBn(6), set_selpar_A50_cBn(7));
  f_priors+=neg_log_prior(selpar_slope_cBn,set_selpar_slope_cBn(5), set_selpar_slope_cBn(6), set_selpar_slope_cBn(7));
  f_priors+=neg_log_prior(selpar_A502_cBn,set_selpar_A502_cBn(5), set_selpar_A502_cBn(6), set_selpar_A502_cBn(7));
  f_priors+=neg_log_prior(selpar_slope2_cBn,set_selpar_slope2_cBn(5), set_selpar_slope2_cBn(6), set_selpar_slope2_cBn(7));
  f_priors+=neg_log_prior(selpar_A50_cBn2,set_selpar_A50_cBn2(5), set_selpar_A50_cBn2(6), set_selpar_A50_cBn2(7));
  f_priors+=neg_log_prior(selpar_slope_cBn2,set_selpar_slope_cBn2(5), set_selpar_slope_cBn2(6), set_selpar_slope_cBn2(7));
  f_priors+=neg_log_prior(selpar_A502_cBn2,set_selpar_A502_cBn2(5), set_selpar_A502_cBn2(6), set_selpar_A502_cBn2(7));
  f_priors+=neg_log_prior(selpar_slope2_cBn2,set_selpar_slope2_cBn2(5), set_selpar_slope2_cBn2(6), set_selpar_slope2_cBn2(7));
  f_priors+=neg_log_prior(selpar_A50_cBs,set_selpar_A50_cBs(5), set_selpar_A50_cBs(6), set_selpar_A50_cBs(7));
  f_priors+=neg_log_prior(selpar_slope_cBs,set_selpar_slope_cBs(5), set_selpar_slope_cBs(6), set_selpar_slope_cBs(7));
  f_priors+=neg_log_prior(selpar_A502_cBs,set_selpar_A502_cBs(5), set_selpar_A502_cBs(6), set_selpar_A502_cBs(7));
  f_priors+=neg_log_prior(selpar_slope2_cBs,set_selpar_slope2_cBs(5), set_selpar_slope2_cBs(6), set_selpar_slope2_cBs(7));
  f_priors+=neg_log_prior(selpar_A50_cBs2,set_selpar_A50_cBs2(5), set_selpar_A50_cBs2(6), set_selpar_A50_cBs2(7));
  f_priors+=neg_log_prior(selpar_slope_cBs2,set_selpar_slope_cBs2(5), set_selpar_slope_cBs2(6), set_selpar_slope_cBs2(7));
  f_priors+=neg_log_prior(selpar_A502_cBs2,set_selpar_A502_cBs2(5), set_selpar_A502_cBs2(6), set_selpar_A502_cBs2(7));
  f_priors+=neg_log_prior(selpar_slope2_cBs2,set_selpar_slope2_cBs2(5), set_selpar_slope2_cBs2(6), set_selpar_slope2_cBs2(7));
  f_priors+=neg_log_prior(sel_age0_cRn_logit,set_sel_age0_cRn(5),set_sel_age0_cRn(6), set_sel_age0_cRn(7));
  f_priors+=neg_log_prior(sel_age1_cRn_logit,set_sel_age1_cRn(5),set_sel_age1_cRn(6), set_sel_age1_cRn(7));
  f_priors+=neg_log_prior(sel_age2_cRn_logit,set_sel_age2_cRn(5),set_sel_age2_cRn(6), set_sel_age2_cRn(7));
  f_priors+=neg_log_prior(sel_age3_cRn_logit,set_sel_age3_cRn(5),set_sel_age3_cRn(6), set_sel_age3_cRn(7));
  f_priors+=neg_log_prior(sel_age4_cRn_logit,set_sel_age4_cRn(5),set_sel_age4_cRn(6), set_sel_age4_cRn(7));
  f_priors+=neg_log_prior(sel_age5_cRn_logit,set_sel_age5_cRn(5),set_sel_age5_cRn(6), set_sel_age5_cRn(7));
  f_priors+=neg_log_prior(sel_age6_cRn_logit,set_sel_age6_cRn(5),set_sel_age6_cRn(6), set_sel_age6_cRn(7));
  f_priors+=neg_log_prior(sel_age0_cRn2_logit,set_sel_age0_cRn2(5),set_sel_age0_cRn2(6), set_sel_age0_cRn2(7));
  f_priors+=neg_log_prior(sel_age1_cRn2_logit,set_sel_age1_cRn2(5),set_sel_age1_cRn2(6), set_sel_age1_cRn2(7));
  f_priors+=neg_log_prior(sel_age2_cRn2_logit,set_sel_age2_cRn2(5),set_sel_age2_cRn2(6), set_sel_age2_cRn2(7));
  f_priors+=neg_log_prior(sel_age3_cRn2_logit,set_sel_age3_cRn2(5),set_sel_age3_cRn2(6), set_sel_age3_cRn2(7));
  f_priors+=neg_log_prior(sel_age4_cRn2_logit,set_sel_age4_cRn2(5),set_sel_age4_cRn2(6), set_sel_age4_cRn2(7));
  f_priors+=neg_log_prior(sel_age5_cRn2_logit,set_sel_age5_cRn2(5),set_sel_age5_cRn2(6), set_sel_age5_cRn2(7));
  f_priors+=neg_log_prior(sel_age6_cRn2_logit,set_sel_age6_cRn2(5),set_sel_age6_cRn2(6), set_sel_age6_cRn2(7));
  f_priors+=neg_log_prior(sel_age0_cRn3_logit,set_sel_age0_cRn3(5),set_sel_age0_cRn3(6), set_sel_age0_cRn3(7));
  f_priors+=neg_log_prior(sel_age1_cRn3_logit,set_sel_age1_cRn3(5),set_sel_age1_cRn3(6), set_sel_age1_cRn3(7));
  f_priors+=neg_log_prior(sel_age2_cRn3_logit,set_sel_age2_cRn3(5),set_sel_age2_cRn3(6), set_sel_age2_cRn3(7));
  f_priors+=neg_log_prior(sel_age3_cRn3_logit,set_sel_age3_cRn3(5),set_sel_age3_cRn3(6), set_sel_age3_cRn3(7));
  f_priors+=neg_log_prior(sel_age4_cRn3_logit,set_sel_age4_cRn3(5),set_sel_age4_cRn3(6), set_sel_age4_cRn3(7));
  f_priors+=neg_log_prior(sel_age5_cRn3_logit,set_sel_age5_cRn3(5),set_sel_age5_cRn3(6), set_sel_age5_cRn3(7));
  f_priors+=neg_log_prior(sel_age6_cRn3_logit,set_sel_age6_cRn3(5),set_sel_age6_cRn3(6), set_sel_age6_cRn3(7));
  f_priors+=neg_log_prior(sel_age0_cRn4_logit,set_sel_age0_cRn4(5),set_sel_age0_cRn4(6), set_sel_age0_cRn4(7));
  f_priors+=neg_log_prior(sel_age1_cRn4_logit,set_sel_age1_cRn4(5),set_sel_age1_cRn4(6), set_sel_age1_cRn4(7));
  f_priors+=neg_log_prior(sel_age2_cRn4_logit,set_sel_age2_cRn4(5),set_sel_age2_cRn4(6), set_sel_age2_cRn4(7));
  f_priors+=neg_log_prior(sel_age3_cRn4_logit,set_sel_age3_cRn4(5),set_sel_age3_cRn4(6), set_sel_age3_cRn4(7));
  f_priors+=neg_log_prior(sel_age4_cRn4_logit,set_sel_age4_cRn4(5),set_sel_age4_cRn4(6), set_sel_age4_cRn4(7));
  f_priors+=neg_log_prior(sel_age5_cRn4_logit,set_sel_age5_cRn4(5),set_sel_age5_cRn4(6), set_sel_age5_cRn4(7));
  f_priors+=neg_log_prior(sel_age6_cRn4_logit,set_sel_age6_cRn4(5),set_sel_age6_cRn4(6), set_sel_age6_cRn4(7));
  f_priors+=neg_log_prior(sel_age0_cRs_logit,set_sel_age0_cRs(5),set_sel_age0_cRs(6), set_sel_age0_cRs(7));
  f_priors+=neg_log_prior(sel_age1_cRs_logit,set_sel_age1_cRs(5),set_sel_age1_cRs(6), set_sel_age1_cRs(7));
  f_priors+=neg_log_prior(sel_age2_cRs_logit,set_sel_age2_cRs(5),set_sel_age2_cRs(6), set_sel_age2_cRs(7));
  f_priors+=neg_log_prior(sel_age3_cRs_logit,set_sel_age3_cRs(5),set_sel_age3_cRs(6), set_sel_age3_cRs(7));
  f_priors+=neg_log_prior(sel_age4_cRs_logit,set_sel_age4_cRs(5),set_sel_age4_cRs(6), set_sel_age4_cRs(7));
  f_priors+=neg_log_prior(sel_age5_cRs_logit,set_sel_age5_cRs(5),set_sel_age5_cRs(6), set_sel_age5_cRs(7));
  f_priors+=neg_log_prior(sel_age6_cRs_logit,set_sel_age6_cRs(5),set_sel_age6_cRs(6), set_sel_age6_cRs(7)); 
  f_priors+=neg_log_prior(sel_age0_cRs2_logit,set_sel_age0_cRs2(5),set_sel_age0_cRs2(6), set_sel_age0_cRs2(7));
  f_priors+=neg_log_prior(sel_age1_cRs2_logit,set_sel_age1_cRs2(5),set_sel_age1_cRs2(6), set_sel_age1_cRs2(7));
  f_priors+=neg_log_prior(sel_age2_cRs2_logit,set_sel_age2_cRs2(5),set_sel_age2_cRs2(6), set_sel_age2_cRs2(7));
  f_priors+=neg_log_prior(sel_age3_cRs2_logit,set_sel_age3_cRs2(5),set_sel_age3_cRs2(6), set_sel_age3_cRs2(7));
  f_priors+=neg_log_prior(sel_age4_cRs2_logit,set_sel_age4_cRs2(5),set_sel_age4_cRs2(6), set_sel_age4_cRs2(7));
  f_priors+=neg_log_prior(sel_age5_cRs2_logit,set_sel_age5_cRs2(5),set_sel_age5_cRs2(6), set_sel_age5_cRs2(7));
  f_priors+=neg_log_prior(sel_age6_cRs2_logit,set_sel_age6_cRs2(5),set_sel_age6_cRs2(6), set_sel_age6_cRs2(7)); 
  f_priors+=neg_log_prior(sel_age0_cRs3_logit,set_sel_age0_cRs3(5),set_sel_age0_cRs3(6), set_sel_age0_cRs3(7));
  f_priors+=neg_log_prior(sel_age1_cRs3_logit,set_sel_age1_cRs3(5),set_sel_age1_cRs3(6), set_sel_age1_cRs3(7));
  f_priors+=neg_log_prior(sel_age2_cRs3_logit,set_sel_age2_cRs3(5),set_sel_age2_cRs3(6), set_sel_age2_cRs3(7));
  f_priors+=neg_log_prior(sel_age3_cRs3_logit,set_sel_age3_cRs3(5),set_sel_age3_cRs3(6), set_sel_age3_cRs3(7));
  f_priors+=neg_log_prior(sel_age4_cRs3_logit,set_sel_age4_cRs3(5),set_sel_age4_cRs3(6), set_sel_age4_cRs3(7));
  f_priors+=neg_log_prior(sel_age5_cRs3_logit,set_sel_age5_cRs3(5),set_sel_age5_cRs3(6), set_sel_age5_cRs3(7));
  f_priors+=neg_log_prior(sel_age6_cRs3_logit,set_sel_age6_cRs3(5),set_sel_age6_cRs3(6), set_sel_age6_cRs3(7)); 
  f_priors+=neg_log_prior(sel_age0_cRs4_logit,set_sel_age0_cRs4(5),set_sel_age0_cRs4(6), set_sel_age0_cRs4(7));
  f_priors+=neg_log_prior(sel_age1_cRs4_logit,set_sel_age1_cRs4(5),set_sel_age1_cRs4(6), set_sel_age1_cRs4(7));
  f_priors+=neg_log_prior(sel_age2_cRs4_logit,set_sel_age2_cRs4(5),set_sel_age2_cRs4(6), set_sel_age2_cRs4(7));
  f_priors+=neg_log_prior(sel_age3_cRs4_logit,set_sel_age3_cRs4(5),set_sel_age3_cRs4(6), set_sel_age3_cRs4(7));
  f_priors+=neg_log_prior(sel_age4_cRs4_logit,set_sel_age4_cRs4(5),set_sel_age4_cRs4(6), set_sel_age4_cRs4(7));
  f_priors+=neg_log_prior(sel_age5_cRs4_logit,set_sel_age5_cRs4(5),set_sel_age5_cRs4(6), set_sel_age5_cRs4(7));
  f_priors+=neg_log_prior(sel_age6_cRs4_logit,set_sel_age6_cRs4(5),set_sel_age6_cRs4(6), set_sel_age6_cRs4(7)); 
  f_priors+=neg_log_prior(sel_age0_cBn_logit,set_sel_age0_cBn(5),set_sel_age0_cBn(6), set_sel_age0_cBn(7));
  f_priors+=neg_log_prior(sel_age1_cBn_logit,set_sel_age1_cBn(5),set_sel_age1_cBn(6), set_sel_age1_cBn(7));
  f_priors+=neg_log_prior(sel_age2_cBn_logit,set_sel_age2_cBn(5),set_sel_age2_cBn(6), set_sel_age2_cBn(7));
  f_priors+=neg_log_prior(sel_age3_cBn_logit,set_sel_age3_cBn(5),set_sel_age3_cBn(6), set_sel_age3_cBn(7));
  f_priors+=neg_log_prior(sel_age4_cBn_logit,set_sel_age4_cBn(5),set_sel_age4_cBn(6), set_sel_age4_cBn(7));
  f_priors+=neg_log_prior(sel_age5_cBn_logit,set_sel_age5_cBn(5),set_sel_age5_cBn(6), set_sel_age5_cBn(7));
  f_priors+=neg_log_prior(sel_age6_cBn_logit,set_sel_age6_cBn(5),set_sel_age6_cBn(6), set_sel_age6_cBn(7));
  f_priors+=neg_log_prior(sel_age0_cBn2_logit,set_sel_age0_cBn2(5),set_sel_age0_cBn2(6), set_sel_age0_cBn2(7));
  f_priors+=neg_log_prior(sel_age1_cBn2_logit,set_sel_age1_cBn2(5),set_sel_age1_cBn2(6), set_sel_age1_cBn2(7));
  f_priors+=neg_log_prior(sel_age2_cBn2_logit,set_sel_age2_cBn2(5),set_sel_age2_cBn2(6), set_sel_age2_cBn2(7));
  f_priors+=neg_log_prior(sel_age3_cBn2_logit,set_sel_age3_cBn2(5),set_sel_age3_cBn2(6), set_sel_age3_cBn2(7));
  f_priors+=neg_log_prior(sel_age4_cBn2_logit,set_sel_age4_cBn2(5),set_sel_age4_cBn2(6), set_sel_age4_cBn2(7));
  f_priors+=neg_log_prior(sel_age5_cBn2_logit,set_sel_age5_cBn2(5),set_sel_age5_cBn2(6), set_sel_age5_cBn2(7));
  f_priors+=neg_log_prior(sel_age6_cBn2_logit,set_sel_age6_cBn2(5),set_sel_age6_cBn2(6), set_sel_age6_cBn2(7));
  f_priors+=neg_log_prior(sel_age0_cBs_logit,set_sel_age0_cBs(5),set_sel_age0_cBs(6), set_sel_age0_cBs(7));
  f_priors+=neg_log_prior(sel_age1_cBs_logit,set_sel_age1_cBs(5),set_sel_age1_cBs(6), set_sel_age1_cBs(7));
  f_priors+=neg_log_prior(sel_age2_cBs_logit,set_sel_age2_cBs(5),set_sel_age2_cBs(6), set_sel_age2_cBs(7));
  f_priors+=neg_log_prior(sel_age3_cBs_logit,set_sel_age3_cBs(5),set_sel_age3_cBs(6), set_sel_age3_cBs(7));
  f_priors+=neg_log_prior(sel_age4_cBs_logit,set_sel_age4_cBs(5),set_sel_age4_cBs(6), set_sel_age4_cBs(7));
  f_priors+=neg_log_prior(sel_age5_cBs_logit,set_sel_age5_cBs(5),set_sel_age5_cBs(6), set_sel_age5_cBs(7));
  f_priors+=neg_log_prior(sel_age6_cBs_logit,set_sel_age6_cBs(5),set_sel_age6_cBs(6), set_sel_age6_cBs(7)); 
  f_priors+=neg_log_prior(sel_age0_cBs2_logit,set_sel_age0_cBs2(5),set_sel_age0_cBs2(6), set_sel_age0_cBs2(7));
  f_priors+=neg_log_prior(sel_age1_cBs2_logit,set_sel_age1_cBs2(5),set_sel_age1_cBs2(6), set_sel_age1_cBs2(7));
  f_priors+=neg_log_prior(sel_age2_cBs2_logit,set_sel_age2_cBs2(5),set_sel_age2_cBs2(6), set_sel_age2_cBs2(7));
  f_priors+=neg_log_prior(sel_age3_cBs2_logit,set_sel_age3_cBs2(5),set_sel_age3_cBs2(6), set_sel_age3_cBs2(7));
  f_priors+=neg_log_prior(sel_age4_cBs2_logit,set_sel_age4_cBs2(5),set_sel_age4_cBs2(6), set_sel_age4_cBs2(7));
  f_priors+=neg_log_prior(sel_age5_cBs2_logit,set_sel_age5_cBs2(5),set_sel_age5_cBs2(6), set_sel_age5_cBs2(7));
  f_priors+=neg_log_prior(sel_age6_cBs2_logit,set_sel_age6_cBs2(5),set_sel_age6_cBs2(6), set_sel_age6_cBs2(7)); 
  f_priors+=neg_log_prior(selpar_A50_nad,set_selpar_A50_nad(5), set_selpar_A50_nad(6), set_selpar_A50_nad(7));
  f_priors+=neg_log_prior(selpar_slope_nad,set_selpar_slope_nad(5), set_selpar_slope_nad(6), set_selpar_slope_nad(7));
  f_priors+=neg_log_prior(selpar_A502_nad,set_selpar_A502_nad(5), set_selpar_A502_nad(6), set_selpar_A502_nad(7));
  f_priors+=neg_log_prior(selpar_slope2_nad,set_selpar_slope2_nad(5), set_selpar_slope2_nad(6), set_selpar_slope2_nad(7));
  f_priors+=neg_log_prior(sel_age0_nad_logit,set_sel_age0_nad(5),set_sel_age0_nad(6), set_sel_age0_nad(7));
  f_priors+=neg_log_prior(sel_age1_nad_logit,set_sel_age1_nad(5),set_sel_age1_nad(6), set_sel_age1_nad(7));
  f_priors+=neg_log_prior(sel_age2_nad_logit,set_sel_age2_nad(5),set_sel_age2_nad(6), set_sel_age2_nad(7));
  f_priors+=neg_log_prior(sel_age3_nad_logit,set_sel_age3_nad(5),set_sel_age3_nad(6), set_sel_age3_nad(7));
  f_priors+=neg_log_prior(sel_age4_nad_logit,set_sel_age4_nad(5),set_sel_age4_nad(6), set_sel_age4_nad(7));
  f_priors+=neg_log_prior(sel_age5_nad_logit,set_sel_age5_nad(5),set_sel_age5_nad(6), set_sel_age5_nad(7));
  f_priors+=neg_log_prior(sel_age6_nad_logit,set_sel_age6_nad(5),set_sel_age6_nad(6), set_sel_age6_nad(7)); 
  f_priors+=neg_log_prior(selpar_A50_mad,set_selpar_A50_mad(5), set_selpar_A50_mad(6), set_selpar_A50_mad(7));
  f_priors+=neg_log_prior(selpar_slope_mad,set_selpar_slope_mad(5), set_selpar_slope_mad(6), set_selpar_slope_mad(7));
  f_priors+=neg_log_prior(selpar_A502_mad,set_selpar_A502_mad(5), set_selpar_A502_mad(6), set_selpar_A502_mad(7));
  f_priors+=neg_log_prior(selpar_slope2_mad,set_selpar_slope2_mad(5), set_selpar_slope2_mad(6), set_selpar_slope2_mad(7));
  f_priors+=neg_log_prior(sel_age0_mad_logit,set_sel_age0_mad(5),set_sel_age0_mad(6), set_sel_age0_mad(7));
  f_priors+=neg_log_prior(sel_age1_mad_logit,set_sel_age1_mad(5),set_sel_age1_mad(6), set_sel_age1_mad(7));
  f_priors+=neg_log_prior(sel_age2_mad_logit,set_sel_age2_mad(5),set_sel_age2_mad(6), set_sel_age2_mad(7));
  f_priors+=neg_log_prior(sel_age3_mad_logit,set_sel_age3_mad(5),set_sel_age3_mad(6), set_sel_age3_mad(7));
  f_priors+=neg_log_prior(sel_age4_mad_logit,set_sel_age4_mad(5),set_sel_age4_mad(6), set_sel_age4_mad(7));
  f_priors+=neg_log_prior(sel_age5_mad_logit,set_sel_age5_mad(5),set_sel_age5_mad(6), set_sel_age5_mad(7));
  f_priors+=neg_log_prior(sel_age6_mad_logit,set_sel_age6_mad(5),set_sel_age6_mad(6), set_sel_age6_mad(7)); 
  f_priors+=neg_log_prior(selpar_A50_sad,set_selpar_A50_sad(5), set_selpar_A50_sad(6), set_selpar_A50_sad(7));
  f_priors+=neg_log_prior(selpar_slope_sad,set_selpar_slope_sad(5), set_selpar_slope_sad(6), set_selpar_slope_sad(7));
  f_priors+=neg_log_prior(selpar_A502_sad,set_selpar_A502_sad(5), set_selpar_A502_sad(6), set_selpar_A502_sad(7));
  f_priors+=neg_log_prior(selpar_slope2_sad,set_selpar_slope2_sad(5), set_selpar_slope2_sad(6), set_selpar_slope2_sad(7));
  f_priors+=neg_log_prior(sel_age0_sad_logit,set_sel_age0_sad(5),set_sel_age0_sad(6), set_sel_age0_sad(7));
  f_priors+=neg_log_prior(sel_age1_sad_logit,set_sel_age1_sad(5),set_sel_age1_sad(6), set_sel_age1_sad(7));
  f_priors+=neg_log_prior(sel_age2_sad_logit,set_sel_age2_sad(5),set_sel_age2_sad(6), set_sel_age2_sad(7));
  f_priors+=neg_log_prior(sel_age3_sad_logit,set_sel_age3_sad(5),set_sel_age3_sad(6), set_sel_age3_sad(7));
  f_priors+=neg_log_prior(sel_age4_sad_logit,set_sel_age4_sad(5),set_sel_age4_sad(6), set_sel_age4_sad(7));
  f_priors+=neg_log_prior(sel_age5_sad_logit,set_sel_age5_sad(5),set_sel_age5_sad(6), set_sel_age5_sad(7));
  f_priors+=neg_log_prior(sel_age6_sad_logit,set_sel_age6_sad(5),set_sel_age6_sad(6), set_sel_age6_sad(7)); 
  f_priors+=neg_log_prior(log_q_cpue_nad,set_log_q_cpue_nad(5),set_log_q_cpue_nad(6),set_log_q_cpue_nad(7));
  f_priors+=neg_log_prior(log_q_cpue_mad,set_log_q_cpue_mad(5),set_log_q_cpue_mad(6),set_log_q_cpue_mad(7));
  f_priors+=neg_log_prior(log_q_cpue_sad,set_log_q_cpue_sad(5),set_log_q_cpue_sad(6),set_log_q_cpue_sad(7));
  f_priors+=neg_log_prior(log_q_cpue_jai,set_log_q_cpue_jai(5),set_log_q_cpue_jai(6),set_log_q_cpue_jai(7));
  f_priors+=neg_log_prior(log_q2_jai,set_log_q2_jai(5),set_log_q2_jai(6),set_log_q2_jai(7));
  f_priors+=neg_log_prior(log_q_cpue_mar,set_log_q_cpue_mar(5),set_log_q_cpue_mar(6),set_log_q_cpue_mar(7));
  f_priors+=neg_log_prior(log_q_cpue_eco,set_log_q_cpue_eco(5),set_log_q_cpue_eco(6),set_log_q_cpue_eco(7));
  f_priors+=neg_log_prior(log_dm_lenc_nad,set_log_dm_lenc_nad(5),set_log_dm_lenc_nad(6),set_log_dm_lenc_nad(7));
  f_priors+=neg_log_prior(log_dm_lenc_mad,set_log_dm_lenc_mad(5),set_log_dm_lenc_mad(6),set_log_dm_lenc_mad(7));
  //f_priors+=neg_log_prior(log_dm_sad_lc,set_log_dm_sad_lc(5),set_log_dm_sad_lc(6),set_log_dm_sad_lc(7));
  f_priors+=neg_log_prior(log_dm_agec_cRn,set_log_dm_agec_cRn(5),set_log_dm_agec_cRn(6),set_log_dm_agec_cRn(7));
  f_priors+=neg_log_prior(log_dm_agec_cRs,set_log_dm_agec_cRs(5),set_log_dm_agec_cRs(6),set_log_dm_agec_cRs(7));
  f_priors+=neg_log_prior(log_dm_agec_cBn,set_log_dm_agec_cBn(5),set_log_dm_agec_cBn(6),set_log_dm_agec_cBn(7));
  f_priors+=neg_log_prior(log_dm_agec_cBs,set_log_dm_agec_cBs(5),set_log_dm_agec_cBs(6),set_log_dm_agec_cBs(7));
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
       //cout<<"got weighted"<<endl;
      get_msy();
       //cout<<"got msy"<<endl;
      get_per_recruit_stuff();
       //cout<<"got per recruit"<<endl;  
      get_miscellaneous_stuff();
       //cout<<"got misc stuff"<<endl;
      get_projection();
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
      //cout << "BC Fmsy=" << F_msy_out<< "   BC SSBmsy=" << SSB_msy_out <<endl;
      cout <<"F status="<<FdF_msy_end<<endl;
      cout <<"Pop status="<<SdSSB_msy_end<<endl;
      cout << "h="<<steep<<"   R0="<<R0<<endl;
      //cout << "xdum " << xdum << endl;
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;  
      report << "TotalLikelihood " << fval << endl;
      report << "N" << endl;
      report << N<<endl;
      report << "F" << endl;
      report << F <<endl;
      sdnr_lc_nad=sdnr_multinomial(nyr_lenc_nad, lenbins, nsamp_lenc_nad, pred_nad_lenc, obs_lenc_nad, w_lenc_nad);
      sdnr_lc_mad=sdnr_multinomial(nyr_lenc_mad, lenbins, nsamp_lenc_mad, pred_mad_lenc, obs_lenc_mad, w_lenc_mad);
      //sdnr_lc_sad=sdnr_multinomial(nyr_sad_lenc, lenbins, nsamp_sad_lenc, pred_sad_lenc, obs_sad_lenc, w_lc_sad);
      sdnr_ac_cRn=sdnr_multinomial(nyr_agec_cRn, agebins_agec, nsamp_agec_cRn, pred_cRn_agec, obs_agec_cRn, w_agec_cRn);
      sdnr_ac_cRs=sdnr_multinomial(nyr_agec_cRs, agebins_agec, nsamp_agec_cRs, pred_cRs_agec, obs_agec_cRs, w_agec_cRs);
      sdnr_ac_cBn=sdnr_multinomial(nyr_agec_cBn, agebins_agec, nsamp_agec_cBn, pred_cBn_agec, obs_agec_cBn, w_agec_cBn);
      sdnr_ac_cBs=sdnr_multinomial(nyr_agec_cBs, agebins_agec, nsamp_agec_cBs, pred_cBs_agec, obs_agec_cBs, w_agec_cBs); 
      sdnr_I_nad=sdnr_lognormal(pred_nad_cpue, obs_cpue_nad, obs_cv_cpue_nad, w_cpue_nad);
      sdnr_I_mad=sdnr_lognormal(pred_mad_cpue, obs_cpue_mad, obs_cv_cpue_mad, w_cpue_mad);
      sdnr_I_sad=sdnr_lognormal(pred_sad_cpue, obs_cpue_sad, obs_cv_cpue_sad, w_cpue_sad);
      sdnr_I_jai=sdnr_lognormal(pred_jai_cpue, obs_cpue_jai, obs_cv_cpue_jai, w_cpue_jai);
      sdnr_I_mareco=sdnr_lognormal(pred_mareco_cpue, obs_cpue_mareco, obs_cv_cpue_mareco, w_cpue_mareco);
      //#################################################################################################
      //##  Passing parameters to vector for bounds check plotting
      //################################################################################################# 
       Linf_out(8)=Linf; Linf_out(1,7)=set_Linf; 
       K_out(8)=K; K_out(1,7)=set_K;
       t0_out(8)=t0; t0_out(1,7)=set_t0;
       len_cv_nad_val_out(8)=len_cv_nad_val; len_cv_nad_val_out(1,7)=set_len_cv_nad;
       len_cv_mad_val_out(8)=len_cv_mad_val; len_cv_mad_val_out(1,7)=set_len_cv_mad;
       //len_cv_sad_val_out(8)=len_cv_sad_val; len_cv_sad_val_out(1,7)=set_len_cv_sad;
       log_R0_out(8)=log_R0; log_R0_out(1,7)=set_log_R0;
       steep_out(8)=steep; steep_out(1,7)=set_steep;
       rec_sigma_out(8)=rec_sigma; rec_sigma_out(1,7)=set_rec_sigma;
       R_autocorr_out(8)=R_autocorr; R_autocorr_out(1,7)=set_R_autocorr;
       log_dm_nad_lc_out(8)=log_dm_lenc_nad; log_dm_nad_lc_out(1,7)=set_log_dm_lenc_nad;
       log_dm_mad_lc_out(8)=log_dm_lenc_mad; log_dm_mad_lc_out(1,7)=set_log_dm_lenc_mad;
       //log_dm_sad_lc_out(8)=log_dm_sad_lc; log_dm_sad_lc_out(1,7)=set_log_dm_sad_lc;
       log_dm_cRn_ac_out(8)=log_dm_agec_cRn; log_dm_cRn_ac_out(1,7)=set_log_dm_agec_cRn;
       log_dm_cRs_ac_out(8)=log_dm_agec_cRs; log_dm_cRs_ac_out(1,7)=set_log_dm_agec_cRs;
       log_dm_cBn_ac_out(8)=log_dm_agec_cBn; log_dm_cBn_ac_out(1,7)=set_log_dm_agec_cBn;
       log_dm_cBs_ac_out(8)=log_dm_agec_cBs; log_dm_cBs_ac_out(1,7)=set_log_dm_agec_cBs;
       selpar_A50_cRn_out(8)=selpar_A50_cRn; selpar_A50_cRn_out(1,7)=set_selpar_A50_cRn;
       selpar_slope_cRn_out(8)=selpar_slope_cRn; selpar_slope_cRn_out(1,7)=set_selpar_slope_cRn;
       selpar_A502_cRn_out(8)=selpar_A502_cRn; selpar_A502_cRn_out(1,7)=set_selpar_A502_cRn;
       selpar_slope2_cRn_out(8)=selpar_slope2_cRn; selpar_slope2_cRn_out(1,7)=set_selpar_slope2_cRn;
       selpar_A50_cRn2_out(8)=selpar_A50_cRn2; selpar_A50_cRn2_out(1,7)=set_selpar_A50_cRn2;
       selpar_slope_cRn2_out(8)=selpar_slope_cRn2; selpar_slope_cRn2_out(1,7)=set_selpar_slope_cRn2;
       selpar_A502_cRn2_out(8)=selpar_A502_cRn2; selpar_A502_cRn2_out(1,7)=set_selpar_A502_cRn2;
       selpar_slope2_cRn2_out(8)=selpar_slope2_cRn2; selpar_slope2_cRn2_out(1,7)=set_selpar_slope2_cRn2;
       selpar_A50_cRn3_out(8)=selpar_A50_cRn3; selpar_A50_cRn3_out(1,7)=set_selpar_A50_cRn3;
       selpar_slope_cRn3_out(8)=selpar_slope_cRn3; selpar_slope_cRn3_out(1,7)=set_selpar_slope_cRn3;
       selpar_A502_cRn3_out(8)=selpar_A502_cRn3; selpar_A502_cRn3_out(1,7)=set_selpar_A502_cRn3;
       selpar_slope2_cRn3_out(8)=selpar_slope2_cRn3; selpar_slope2_cRn3_out(1,7)=set_selpar_slope2_cRn3;
       selpar_A50_cRn4_out(8)=selpar_A50_cRn4; selpar_A50_cRn4_out(1,7)=set_selpar_A50_cRn4;
       selpar_slope_cRn4_out(8)=selpar_slope_cRn4; selpar_slope_cRn4_out(1,7)=set_selpar_slope_cRn4;
       selpar_A502_cRn4_out(8)=selpar_A502_cRn4; selpar_A502_cRn4_out(1,7)=set_selpar_A502_cRn4;
       selpar_slope2_cRn4_out(8)=selpar_slope2_cRn4; selpar_slope2_cRn4_out(1,7)=set_selpar_slope2_cRn4;
       selpar_A50_cRs_out(8)=selpar_A50_cRs; selpar_A50_cRs_out(1,7)=set_selpar_A50_cRs;
       selpar_slope_cRs_out(8)=selpar_slope_cRs; selpar_slope_cRs_out(1,7)=set_selpar_slope_cRs;
       selpar_A502_cRs_out(8)=selpar_A502_cRs; selpar_A502_cRs_out(1,7)=set_selpar_A502_cRs;
       selpar_slope2_cRs_out(8)=selpar_slope2_cRs; selpar_slope2_cRs_out(1,7)=set_selpar_slope2_cRs;
       selpar_A50_cRs2_out(8)=selpar_A50_cRs2; selpar_A50_cRs2_out(1,7)=set_selpar_A50_cRs2;
       selpar_slope_cRs2_out(8)=selpar_slope_cRs2; selpar_slope_cRs2_out(1,7)=set_selpar_slope_cRs2;
       selpar_A502_cRs2_out(8)=selpar_A502_cRs2; selpar_A502_cRs2_out(1,7)=set_selpar_A502_cRs2;
       selpar_slope2_cRs2_out(8)=selpar_slope2_cRs2; selpar_slope2_cRs2_out(1,7)=set_selpar_slope2_cRs2;
       selpar_A50_cRs3_out(8)=selpar_A50_cRs3; selpar_A50_cRs3_out(1,7)=set_selpar_A50_cRs3;
       selpar_slope_cRs3_out(8)=selpar_slope_cRs3; selpar_slope_cRs3_out(1,7)=set_selpar_slope_cRs3;
       selpar_A502_cRs3_out(8)=selpar_A502_cRs3; selpar_A502_cRs3_out(1,7)=set_selpar_A502_cRs3;
       selpar_slope2_cRs3_out(8)=selpar_slope2_cRs3; selpar_slope2_cRs3_out(1,7)=set_selpar_slope2_cRs3;
       selpar_A50_cRs4_out(8)=selpar_A50_cRs4; selpar_A50_cRs4_out(1,7)=set_selpar_A50_cRs4;
       selpar_slope_cRs4_out(8)=selpar_slope_cRs4; selpar_slope_cRs4_out(1,7)=set_selpar_slope_cRs4;
       selpar_A502_cRs4_out(8)=selpar_A502_cRs4; selpar_A502_cRs4_out(1,7)=set_selpar_A502_cRs4;
       selpar_slope2_cRs4_out(8)=selpar_slope2_cRs4; selpar_slope2_cRs4_out(1,7)=set_selpar_slope2_cRs4;
       selpar_A50_cBn_out(8)=selpar_A50_cBn; selpar_A50_cBn_out(1,7)=set_selpar_A50_cBn;
       selpar_slope_cBn_out(8)=selpar_slope_cBn; selpar_slope_cBn_out(1,7)=set_selpar_slope_cBn;
       selpar_A502_cBn_out(8)=selpar_A502_cBn; selpar_A502_cBn_out(1,7)=set_selpar_A502_cBn;
       selpar_slope2_cBn_out(8)=selpar_slope2_cBn; selpar_slope2_cBn_out(1,7)=set_selpar_slope2_cBn;
       selpar_A50_cBn2_out(8)=selpar_A50_cBn2; selpar_A50_cBn2_out(1,7)=set_selpar_A50_cBn2;
       selpar_slope_cBn2_out(8)=selpar_slope_cBn2; selpar_slope_cBn2_out(1,7)=set_selpar_slope_cBn2;
       selpar_A502_cBn2_out(8)=selpar_A502_cBn2; selpar_A502_cBn2_out(1,7)=set_selpar_A502_cBn2;
       selpar_slope2_cBn2_out(8)=selpar_slope2_cBn2; selpar_slope2_cBn2_out(1,7)=set_selpar_slope2_cBn2;
       selpar_A50_cBs_out(8)=selpar_A50_cBs; selpar_A50_cBs_out(1,7)=set_selpar_A50_cBs;
       selpar_slope_cBs_out(8)=selpar_slope_cBs; selpar_slope_cBs_out(1,7)=set_selpar_slope_cBs;
       selpar_A502_cBs_out(8)=selpar_A502_cBs; selpar_A502_cBs_out(1,7)=set_selpar_A502_cBs;
       selpar_slope2_cBs_out(8)=selpar_slope2_cBs; selpar_slope2_cBs_out(1,7)=set_selpar_slope2_cBs;
       selpar_A50_cBs2_out(8)=selpar_A50_cBs2; selpar_A50_cBs2_out(1,7)=set_selpar_A50_cBs2;
       selpar_slope_cBs2_out(8)=selpar_slope_cBs2; selpar_slope_cBs2_out(1,7)=set_selpar_slope_cBs2;
       selpar_A502_cBs2_out(8)=selpar_A502_cBs2; selpar_A502_cBs2_out(1,7)=set_selpar_A502_cBs2;
       selpar_slope2_cBs2_out(8)=selpar_slope2_cBs2; selpar_slope2_cBs2_out(1,7)=set_selpar_slope2_cBs2;
       selpar_age0_cRn_out(8)=sel_age0_cRn_logit; selpar_age0_cRn_out(1,7)=set_sel_age0_cRn;
       selpar_age1_cRn_out(8)=sel_age1_cRn_logit; selpar_age1_cRn_out(1,7)=set_sel_age1_cRn;
       selpar_age2_cRn_out(8)=sel_age2_cRn_logit; selpar_age2_cRn_out(1,7)=set_sel_age2_cRn;
       selpar_age3_cRn_out(8)=sel_age3_cRn_logit; selpar_age3_cRn_out(1,7)=set_sel_age3_cRn;
       selpar_age4_cRn_out(8)=sel_age4_cRn_logit; selpar_age4_cRn_out(1,7)=set_sel_age4_cRn;
       selpar_age5_cRn_out(8)=sel_age5_cRn_logit; selpar_age5_cRn_out(1,7)=set_sel_age5_cRn;
       selpar_age6_cRn_out(8)=sel_age6_cRn_logit; selpar_age6_cRn_out(1,7)=set_sel_age6_cRn;
       selpar_age0_cRn2_out(8)=sel_age0_cRn2_logit; selpar_age0_cRn2_out(1,7)=set_sel_age0_cRn2;
       selpar_age1_cRn2_out(8)=sel_age1_cRn2_logit; selpar_age1_cRn2_out(1,7)=set_sel_age1_cRn2;
       selpar_age2_cRn2_out(8)=sel_age2_cRn2_logit; selpar_age2_cRn2_out(1,7)=set_sel_age2_cRn2;
       selpar_age3_cRn2_out(8)=sel_age3_cRn2_logit; selpar_age3_cRn2_out(1,7)=set_sel_age3_cRn2;
       selpar_age4_cRn2_out(8)=sel_age4_cRn2_logit; selpar_age4_cRn2_out(1,7)=set_sel_age4_cRn2;
       selpar_age5_cRn2_out(8)=sel_age5_cRn2_logit; selpar_age5_cRn2_out(1,7)=set_sel_age5_cRn2;
       selpar_age6_cRn2_out(8)=sel_age6_cRn2_logit; selpar_age6_cRn2_out(1,7)=set_sel_age6_cRn2;
       selpar_age0_cRn3_out(8)=sel_age0_cRn3_logit; selpar_age0_cRn3_out(1,7)=set_sel_age0_cRn3;
       selpar_age1_cRn3_out(8)=sel_age1_cRn3_logit; selpar_age1_cRn3_out(1,7)=set_sel_age1_cRn3;
       selpar_age2_cRn3_out(8)=sel_age2_cRn3_logit; selpar_age2_cRn3_out(1,7)=set_sel_age2_cRn3;
       selpar_age3_cRn3_out(8)=sel_age3_cRn3_logit; selpar_age3_cRn3_out(1,7)=set_sel_age3_cRn3;
       selpar_age4_cRn3_out(8)=sel_age4_cRn3_logit; selpar_age4_cRn3_out(1,7)=set_sel_age4_cRn3;
       selpar_age5_cRn3_out(8)=sel_age5_cRn3_logit; selpar_age5_cRn3_out(1,7)=set_sel_age5_cRn3;
       selpar_age6_cRn3_out(8)=sel_age6_cRn3_logit; selpar_age6_cRn3_out(1,7)=set_sel_age6_cRn3;
       selpar_age0_cRn4_out(8)=sel_age0_cRn4_logit; selpar_age0_cRn4_out(1,7)=set_sel_age0_cRn4;
       selpar_age1_cRn4_out(8)=sel_age1_cRn4_logit; selpar_age1_cRn4_out(1,7)=set_sel_age1_cRn4;
       selpar_age2_cRn4_out(8)=sel_age2_cRn4_logit; selpar_age2_cRn4_out(1,7)=set_sel_age2_cRn4;
       selpar_age3_cRn4_out(8)=sel_age3_cRn4_logit; selpar_age3_cRn4_out(1,7)=set_sel_age3_cRn4;
       selpar_age4_cRn4_out(8)=sel_age4_cRn4_logit; selpar_age4_cRn4_out(1,7)=set_sel_age4_cRn4;
       selpar_age5_cRn4_out(8)=sel_age5_cRn4_logit; selpar_age5_cRn4_out(1,7)=set_sel_age5_cRn4;
       selpar_age6_cRn4_out(8)=sel_age6_cRn4_logit; selpar_age6_cRn4_out(1,7)=set_sel_age6_cRn4;
       selpar_age0_cRs_out(8)=sel_age0_cRs_logit; selpar_age0_cRs_out(1,7)=set_sel_age0_cRs;
       selpar_age1_cRs_out(8)=sel_age1_cRs_logit; selpar_age1_cRs_out(1,7)=set_sel_age1_cRs;
       selpar_age2_cRs_out(8)=sel_age2_cRs_logit; selpar_age2_cRs_out(1,7)=set_sel_age2_cRs;
       selpar_age3_cRs_out(8)=sel_age3_cRs_logit; selpar_age3_cRs_out(1,7)=set_sel_age3_cRs;
       selpar_age4_cRs_out(8)=sel_age4_cRs_logit; selpar_age4_cRs_out(1,7)=set_sel_age4_cRs;
       selpar_age5_cRs_out(8)=sel_age5_cRs_logit; selpar_age5_cRs_out(1,7)=set_sel_age5_cRs;
       selpar_age6_cRs_out(8)=sel_age6_cRs_logit; selpar_age6_cRs_out(1,7)=set_sel_age6_cRs;
       selpar_age0_cRs2_out(8)=sel_age0_cRs2_logit; selpar_age0_cRs2_out(1,7)=set_sel_age0_cRs2;
       selpar_age1_cRs2_out(8)=sel_age1_cRs2_logit; selpar_age1_cRs2_out(1,7)=set_sel_age1_cRs2;
       selpar_age2_cRs2_out(8)=sel_age2_cRs2_logit; selpar_age2_cRs2_out(1,7)=set_sel_age2_cRs2;
       selpar_age3_cRs2_out(8)=sel_age3_cRs2_logit; selpar_age3_cRs2_out(1,7)=set_sel_age3_cRs2;
       selpar_age4_cRs2_out(8)=sel_age4_cRs2_logit; selpar_age4_cRs2_out(1,7)=set_sel_age4_cRs2;
       selpar_age5_cRs2_out(8)=sel_age5_cRs2_logit; selpar_age5_cRs2_out(1,7)=set_sel_age5_cRs2;
       selpar_age6_cRs2_out(8)=sel_age6_cRs2_logit; selpar_age6_cRs2_out(1,7)=set_sel_age6_cRs2;
       selpar_age0_cRs3_out(8)=sel_age0_cRs3_logit; selpar_age0_cRs3_out(1,7)=set_sel_age0_cRs3;
       selpar_age1_cRs3_out(8)=sel_age1_cRs3_logit; selpar_age1_cRs3_out(1,7)=set_sel_age1_cRs3;
       selpar_age2_cRs3_out(8)=sel_age2_cRs3_logit; selpar_age2_cRs3_out(1,7)=set_sel_age2_cRs3;
       selpar_age3_cRs3_out(8)=sel_age3_cRs3_logit; selpar_age3_cRs3_out(1,7)=set_sel_age3_cRs3;
       selpar_age4_cRs3_out(8)=sel_age4_cRs3_logit; selpar_age4_cRs3_out(1,7)=set_sel_age4_cRs3;
       selpar_age5_cRs3_out(8)=sel_age5_cRs3_logit; selpar_age5_cRs3_out(1,7)=set_sel_age5_cRs3;
       selpar_age6_cRs3_out(8)=sel_age6_cRs3_logit; selpar_age6_cRs3_out(1,7)=set_sel_age6_cRs3;
       selpar_age0_cRs4_out(8)=sel_age0_cRs4_logit; selpar_age0_cRs4_out(1,7)=set_sel_age0_cRs4;
       selpar_age1_cRs4_out(8)=sel_age1_cRs4_logit; selpar_age1_cRs4_out(1,7)=set_sel_age1_cRs4;
       selpar_age2_cRs4_out(8)=sel_age2_cRs4_logit; selpar_age2_cRs4_out(1,7)=set_sel_age2_cRs4;
       selpar_age3_cRs4_out(8)=sel_age3_cRs4_logit; selpar_age3_cRs4_out(1,7)=set_sel_age3_cRs4;
       selpar_age4_cRs4_out(8)=sel_age4_cRs4_logit; selpar_age4_cRs4_out(1,7)=set_sel_age4_cRs4;
       selpar_age5_cRs4_out(8)=sel_age5_cRs4_logit; selpar_age5_cRs4_out(1,7)=set_sel_age5_cRs4;
       selpar_age6_cRs4_out(8)=sel_age6_cRs4_logit; selpar_age6_cRs4_out(1,7)=set_sel_age6_cRs4;
       selpar_age0_cBn_out(8)=sel_age0_cBn_logit; selpar_age0_cBn_out(1,7)=set_sel_age0_cBn;
       selpar_age1_cBn_out(8)=sel_age1_cBn_logit; selpar_age1_cBn_out(1,7)=set_sel_age1_cBn;
       selpar_age2_cBn_out(8)=sel_age2_cBn_logit; selpar_age2_cBn_out(1,7)=set_sel_age2_cBn;
       selpar_age3_cBn_out(8)=sel_age3_cBn_logit; selpar_age3_cBn_out(1,7)=set_sel_age3_cBn;
       selpar_age4_cBn_out(8)=sel_age4_cBn_logit; selpar_age4_cBn_out(1,7)=set_sel_age4_cBn;
       selpar_age5_cBn_out(8)=sel_age5_cBn_logit; selpar_age5_cBn_out(1,7)=set_sel_age5_cBn;
       selpar_age6_cBn_out(8)=sel_age6_cBn_logit; selpar_age6_cBn_out(1,7)=set_sel_age6_cBn;
       selpar_age0_cBn2_out(8)=sel_age0_cBn2_logit; selpar_age0_cBn2_out(1,7)=set_sel_age0_cBn2;
       selpar_age1_cBn2_out(8)=sel_age1_cBn2_logit; selpar_age1_cBn2_out(1,7)=set_sel_age1_cBn2;
       selpar_age2_cBn2_out(8)=sel_age2_cBn2_logit; selpar_age2_cBn2_out(1,7)=set_sel_age2_cBn2;
       selpar_age3_cBn2_out(8)=sel_age3_cBn2_logit; selpar_age3_cBn2_out(1,7)=set_sel_age3_cBn2;
       selpar_age4_cBn2_out(8)=sel_age4_cBn2_logit; selpar_age4_cBn2_out(1,7)=set_sel_age4_cBn2;
       selpar_age5_cBn2_out(8)=sel_age5_cBn2_logit; selpar_age5_cBn2_out(1,7)=set_sel_age5_cBn2;
       selpar_age6_cBn2_out(8)=sel_age6_cBn2_logit; selpar_age6_cBn2_out(1,7)=set_sel_age6_cBn2;
       selpar_age0_cBs_out(8)=sel_age0_cBs_logit; selpar_age0_cBs_out(1,7)=set_sel_age0_cBs;
       selpar_age1_cBs_out(8)=sel_age1_cBs_logit; selpar_age1_cBs_out(1,7)=set_sel_age1_cBs;
       selpar_age2_cBs_out(8)=sel_age2_cBs_logit; selpar_age2_cBs_out(1,7)=set_sel_age2_cBs;
       selpar_age3_cBs_out(8)=sel_age3_cBs_logit; selpar_age3_cBs_out(1,7)=set_sel_age3_cBs;
       selpar_age4_cBs_out(8)=sel_age4_cBs_logit; selpar_age4_cBs_out(1,7)=set_sel_age4_cBs;
       selpar_age5_cBs_out(8)=sel_age5_cBs_logit; selpar_age5_cBs_out(1,7)=set_sel_age5_cBs;
       selpar_age6_cBs_out(8)=sel_age6_cBs_logit; selpar_age6_cBs_out(1,7)=set_sel_age6_cBs;
       selpar_age0_cBs2_out(8)=sel_age0_cBs2_logit; selpar_age0_cBs2_out(1,7)=set_sel_age0_cBs2;
       selpar_age1_cBs2_out(8)=sel_age1_cBs2_logit; selpar_age1_cBs2_out(1,7)=set_sel_age1_cBs2;
       selpar_age2_cBs2_out(8)=sel_age2_cBs2_logit; selpar_age2_cBs2_out(1,7)=set_sel_age2_cBs2;
       selpar_age3_cBs2_out(8)=sel_age3_cBs2_logit; selpar_age3_cBs2_out(1,7)=set_sel_age3_cBs2;
       selpar_age4_cBs2_out(8)=sel_age4_cBs2_logit; selpar_age4_cBs2_out(1,7)=set_sel_age4_cBs2;
       selpar_age5_cBs2_out(8)=sel_age5_cBs2_logit; selpar_age5_cBs2_out(1,7)=set_sel_age5_cBs2;
       selpar_age6_cBs2_out(8)=sel_age6_cBs2_logit; selpar_age6_cBs2_out(1,7)=set_sel_age6_cBs2;
       selpar_A50_nad_out(8)=selpar_A50_nad; selpar_A50_nad_out(1,7)=set_selpar_A50_nad;
       selpar_slope_nad_out(8)=selpar_slope_nad; selpar_slope_nad_out(1,7)=set_selpar_slope_nad;
       selpar_A502_nad_out(8)=selpar_A502_nad; selpar_A502_nad_out(1,7)=set_selpar_A502_nad;
       selpar_slope2_nad_out(8)=selpar_slope2_nad; selpar_slope2_nad_out(1,7)=set_selpar_slope2_nad;
       selpar_age0_nad_out(8)=sel_age0_nad_logit; selpar_age0_nad_out(1,7)=set_sel_age0_nad;
       selpar_age1_nad_out(8)=sel_age1_nad_logit; selpar_age1_nad_out(1,7)=set_sel_age1_nad;
       selpar_age2_nad_out(8)=sel_age2_nad_logit; selpar_age2_nad_out(1,7)=set_sel_age2_nad;
       selpar_age3_nad_out(8)=sel_age3_nad_logit; selpar_age3_nad_out(1,7)=set_sel_age3_nad;
       selpar_age4_nad_out(8)=sel_age4_nad_logit; selpar_age4_nad_out(1,7)=set_sel_age4_nad;
       selpar_age5_nad_out(8)=sel_age5_nad_logit; selpar_age5_nad_out(1,7)=set_sel_age5_nad;
       selpar_age6_nad_out(8)=sel_age6_nad_logit; selpar_age6_nad_out(1,7)=set_sel_age6_nad;
       selpar_A50_mad_out(8)=selpar_A50_mad; selpar_A50_mad_out(1,7)=set_selpar_A50_mad;
       selpar_slope_mad_out(8)=selpar_slope_mad; selpar_slope_mad_out(1,7)=set_selpar_slope_mad;
       selpar_A502_mad_out(8)=selpar_A502_mad; selpar_A502_mad_out(1,7)=set_selpar_A502_mad;
       selpar_slope2_mad_out(8)=selpar_slope2_mad; selpar_slope2_mad_out(1,7)=set_selpar_slope2_mad;
       selpar_age0_mad_out(8)=sel_age0_mad_logit; selpar_age0_mad_out(1,7)=set_sel_age0_mad;
       selpar_age1_mad_out(8)=sel_age1_mad_logit; selpar_age1_mad_out(1,7)=set_sel_age1_mad;
       selpar_age2_mad_out(8)=sel_age2_mad_logit; selpar_age2_mad_out(1,7)=set_sel_age2_mad;
       selpar_age3_mad_out(8)=sel_age3_mad_logit; selpar_age3_mad_out(1,7)=set_sel_age3_mad;
       selpar_age4_mad_out(8)=sel_age4_mad_logit; selpar_age4_mad_out(1,7)=set_sel_age4_mad;
       selpar_age5_mad_out(8)=sel_age5_mad_logit; selpar_age5_mad_out(1,7)=set_sel_age5_mad;
       selpar_age6_mad_out(8)=sel_age6_mad_logit; selpar_age6_mad_out(1,7)=set_sel_age6_mad;
       selpar_A50_sad_out(8)=selpar_A50_sad; selpar_A50_sad_out(1,7)=set_selpar_A50_sad;
       selpar_slope_sad_out(8)=selpar_slope_sad; selpar_slope_sad_out(1,7)=set_selpar_slope_sad;
       selpar_A502_sad_out(8)=selpar_A502_sad; selpar_A502_sad_out(1,7)=set_selpar_A502_sad;
       selpar_slope2_sad_out(8)=selpar_slope2_sad; selpar_slope2_sad_out(1,7)=set_selpar_slope2_sad;
       selpar_age0_sad_out(8)=sel_age0_sad_logit; selpar_age0_sad_out(1,7)=set_sel_age0_sad;
       selpar_age1_sad_out(8)=sel_age1_sad_logit; selpar_age1_sad_out(1,7)=set_sel_age1_sad;
       selpar_age2_sad_out(8)=sel_age2_sad_logit; selpar_age2_sad_out(1,7)=set_sel_age2_sad;
       selpar_age3_sad_out(8)=sel_age3_sad_logit; selpar_age3_sad_out(1,7)=set_sel_age3_sad;
       selpar_age4_sad_out(8)=sel_age4_sad_logit; selpar_age4_sad_out(1,7)=set_sel_age4_sad;
       selpar_age5_sad_out(8)=sel_age5_sad_logit; selpar_age5_sad_out(1,7)=set_sel_age5_sad;
       selpar_age6_sad_out(8)=sel_age6_sad_logit; selpar_age6_sad_out(1,7)=set_sel_age6_sad;
       log_q_nad_out(8)=log_q_cpue_nad; log_q_nad_out(1,7)=set_log_q_cpue_nad;
       log_q_mad_out(8)=log_q_cpue_mad; log_q_mad_out(1,7)=set_log_q_cpue_mad;
       log_q_sad_out(8)=log_q_cpue_sad; log_q_sad_out(1,7)=set_log_q_cpue_sad;
       log_q_jai_out(8)=log_q_cpue_jai; log_q_jai_out(1,7)=set_log_q_cpue_jai;
       log_q2_jai_out(8)=log_q2_jai; log_q2_jai_out(1,7)=set_log_q2_jai;
       log_q_mar_out(8)=log_q_cpue_mar; log_q_mar_out(1,7)=set_log_q_cpue_mar;
       log_q_eco_out(8)=log_q_cpue_eco; log_q_eco_out(1,7)=set_log_q_cpue_eco;
       log_avg_F_cRn_out(8)=log_avg_F_L_cRn; log_avg_F_cRn_out(1,7)=set_log_avg_F_L_cRn;
       log_avg_F_cRs_out(8)=log_avg_F_L_cRs; log_avg_F_cRs_out(1,7)=set_log_avg_F_L_cRs;
       log_avg_F_cBn_out(8)=log_avg_F_L_cBn; log_avg_F_cBn_out(1,7)=set_log_avg_F_L_cBn;
       log_avg_F_cBs_out(8)=log_avg_F_L_cBs; log_avg_F_cBs_out(1,7)=set_log_avg_F_L_cBs;
       log_rec_dev_out(styr_rec_dev, endyr_rec_dev)=log_dev_rec;
       log_F_dev_cRn_out(styr_L_cRn,endyr_L_cRn)=log_dev_F_L_cRn;
       log_F_dev_cRs_out(styr_L_cRs,endyr_L_cRs)=log_dev_F_L_cRs;
       log_F_dev_cBn_out(styr_L_cBn,endyr_L_cBn)=log_dev_F_L_cBn;
       log_F_dev_cBs_out(styr_L_cBs,endyr_L_cBs)=log_dev_F_L_cBs;
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
