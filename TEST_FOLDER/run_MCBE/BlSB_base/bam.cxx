//##  Author: NMFS, Beaufort Lab, Sustainable Fisheries Branch
//##  Analyst: Kate Siegfried
//##  Species: Black Sea Bass
//##  Region: US South Atlantic
//##  SEDAR: 56
//##  Date: 2022-11-18 22:17:15


// Create a file with an R object from AD Model Builder

// open the file using the default AD Model Builder file name, and
// 6 digits of precision
open_r_file(adprogram_name + ".rdat", 6);

// Example of an INFO object
open_r_info_list("info", true);
	wrt_r_item("title", "SEDAR 56 Standard Assessment");
	wrt_r_item("species", "Black sea bass");
	wrt_r_item("model", "Statistical Catch at Age");
	wrt_r_item("rec.model", "BH-steep");
	wrt_r_item("base.run", "bsbxxx.tpl");
	wrt_r_item("units.length", "mm");
	wrt_r_item("units.weight", "lb");
	wrt_r_item("units.biomass", "metric tons");
	wrt_r_item("units.ssb", "1E10 eggs");
	wrt_r_item("units.ypr", "lb");
//	wrt_r_item("units.landings", "1000 lb");
	wrt_r_item("units.discards", "1000 dead fish");
    wrt_r_item("units.numbers", "number fish");
    wrt_r_item("units.naa", "number fish");
	wrt_r_item("units.rec", "number fish");
close_r_info_list();


// VECTOR object of parameters and estimated quantities
open_r_info_list("parms", false);
	wrt_r_item("styr", styr);
	wrt_r_item("endyr", endyr);
	wrt_r_item("styrR", styr_rec_dev);
	wrt_r_item("M.msst", M_constant);
	wrt_r_item("Linf", Linf);
	wrt_r_item("K", K);
	wrt_r_item("t0", t0);
	wrt_r_item("wgt.a", wgtpar_a);
	wrt_r_item("wgt.b", wgtpar_b);
	wrt_r_item("fec.a", fecpar_a);
	wrt_r_item("fec.b", fecpar_b);
	wrt_r_item("fec.batches", fecpar_batches);
	wrt_r_item("fec.scale", fecpar_scale);
	    wrt_r_item("spawn.time", spawn_time_frac);
	wrt_r_item("D.mort.lines", Dmort_HL);
	wrt_r_item("D.mort.rHB.lines", Dmort_rHB_HL);
	wrt_r_item("D.mort.rGN.lines", Dmort_rGN_HL);
	wrt_r_item("D.mort.cPT1", Dmort_cPT1);
	wrt_r_item("D.mort.cPT2", Dmort_cPT2);
	wrt_r_item("q.sBT", mfexp(log_q_sBT));
	wrt_r_item("q.sTV", mfexp(log_q_sTV));
	wrt_r_item("q.cHL", mfexp(log_q_cHL));
	wrt_r_item("q.rHB", mfexp(log_q_rHB));
    //wrt_r_item("q.rHB.D", mfexp(log_q_rHD));
	wrt_r_item("q.rate",q_rate);
	wrt_r_item("q.beta",q_DD_beta);
	wrt_r_item("q.DD.B0.4plus", B0_q_DD);
	wrt_r_item("F.prop.cHL", F_cHL_prop);
	wrt_r_item("F.prop.cPT", F_cPT_prop);
	wrt_r_item("F.prop.rHB", F_rHB_prop);
	wrt_r_item("F.prop.rGN", F_rGN_prop);
	wrt_r_item("F.prop.cGN.D", F_cGN_D_prop);
	wrt_r_item("F.prop.rHB.D", F_rHB_D_prop);
	wrt_r_item("F.prop.rGN.D", F_rGN_D_prop);
	wrt_r_item("F.init.ratio", F_init_ratio);
	wrt_r_item("L.rGN.bias", L_rGN_bias);
	wrt_r_item("L.rHB.bias", L_rHB_bias);
	wrt_r_item("B0", B0);
	wrt_r_item("Bstyr.B0", totB(styr)/B0);
	wrt_r_item("SSB0", S0);
	wrt_r_item("SSBstyr.SSB0", SSB(styr)/S0);
	wrt_r_item("Rstyr.R0", rec(styr)/R0);
	wrt_r_item("BH.biascorr",BiasCor);
	wrt_r_item("BH.Phi0", spr_F0);
	wrt_r_item("BH.R0", R0);
	wrt_r_item("BH.steep", steep);
	wrt_r_item("R.sigma.logdevs", sigma_rec_dev);
	wrt_r_item("R.sigma.par",rec_sigma);
	wrt_r_item("R.autocorr",R_autocorr);
	wrt_r_item("R0", R0); //same as BH.R0, but used in BSR.time.plots
	wrt_r_item("R.virgin.bc", R_virgin); //bias-corrected virgin recruitment
	wrt_r_item("rec.lag", 0.0);
	wrt_r_item("msy.klb", msy_klb_out);
	wrt_r_item("msy.knum", msy_knum_out);
	wrt_r_item("Fmsy", F_msy_out);
	wrt_r_item("SSBmsy", SSB_msy_out);
	wrt_r_item("msst", (1.0-M_constant)*SSB_msy_out);
	wrt_r_item("Bmsy", B_msy_out);
	wrt_r_item("Rmsy", R_msy_out);
	wrt_r_item("sprmsy",spr_msy_out);
	wrt_r_item("Dmsy.knum", D_msy_knum_out);
	wrt_r_item("Dmsy.klb", D_msy_klb_out);
	wrt_r_item("Fend.Fmsy", FdF_msy_end);
	wrt_r_item("Fend.Fmsy.mean", FdF_msy_end_mean);
	wrt_r_item("SSBend.SSBmsy", SdSSB_msy_end);
	wrt_r_item("SSBend.MSST", SdSSB_msy_end/(1.0-M_constant));
close_r_info_list();

// MATRIX object of parameter constraints
open_r_df("parm.cons",1,8,2);
    wrt_r_namevector(1,8);
    wrt_r_df_col("Linf",Linf_out);
    wrt_r_df_col("K",K_out);
    wrt_r_df_col("t0",t0_out);
    wrt_r_df_col("len.cv.val",len_cv_val_out);
    wrt_r_df_col("log.R0",log_R0_out);
    wrt_r_df_col("steep",steep_out);
    wrt_r_df_col("rec.sigma",rec_sigma_out);
    wrt_r_df_col("R.autocorr",R_autocorr_out);
    wrt_r_df_col("log.dm.lenc.sTV",log_dm_sTV_lc_out);
	wrt_r_df_col("log.dm.lenc.sBT",log_dm_sBT_lc_out);
	wrt_r_df_col("log.dm.lenc.cHL",log_dm_cHL_lc_out);
	wrt_r_df_col("log.dm.lenc.cPT",log_dm_cPT_lc_out);
	//wrt_r_df_col("log.dm.cTW.lc",log_dm_cTW_lc_out);
	wrt_r_df_col("log.dm.lenc.rHB",log_dm_rHB_lc_out);
	wrt_r_df_col("log.dm.lenc.rGN",log_dm_rGN_lc_out);
	wrt_r_df_col("log.dm.lenc.rHB.D",log_dm_rHB_D_lc_out);
	wrt_r_df_col("log.dm.agec.sTV",log_dm_sTV_ac_out);
	wrt_r_df_col("log.dm.agec.sBT",log_dm_sBT_ac_out);
	wrt_r_df_col("log.dm.agec.cHL",log_dm_cHL_ac_out);
	wrt_r_df_col("log.dm.agec.cPT",log_dm_cPT_ac_out);
	wrt_r_df_col("log.dm.agec.rHB",log_dm_rHB_ac_out);
	//wrt_r_df_col("log.dm.rGN.ac",log_dm_rGN_ac_out);
			
			
	wrt_r_df_col("selpar.A50.sBT", selpar_A50_sBT_out);
	wrt_r_df_col("selpar.slope.sBT", selpar_slope_sBT_out);
   
	wrt_r_df_col("selpar.A50.sTV", selpar_A50_sTV_out);
	wrt_r_df_col("selpar.slope.sTV", selpar_slope_sTV_out);
	
	wrt_r_df_col("selpar.A50.cHL2", selpar_A50_cHL2_out);
	wrt_r_df_col("selpar.slope.cHL2", selpar_slope_cHL2_out);
	wrt_r_df_col("selpar.A50.cHL3", selpar_A50_cHL3_out);
	wrt_r_df_col("selpar.slope.cHL3", selpar_slope_cHL3_out);
	wrt_r_df_col("selpar.A50.cHL4", selpar_A50_cHL4_out);
	wrt_r_df_col("selpar.slope.cHL4", selpar_slope_cHL4_out);
    
	wrt_r_df_col("selpar.A50.cPT2", selpar_A50_cPT2_out);
	wrt_r_df_col("selpar.slope.cPT2", selpar_slope_cPT2_out);
	wrt_r_df_col("selpar.A50.cPT3", selpar_A50_cPT3_out);
	wrt_r_df_col("selpar.slope.cPT3", selpar_slope_cPT3_out);
	wrt_r_df_col("selpar.A50.cPT4", selpar_A50_cPT4_out);
	wrt_r_df_col("selpar.slope.cPT4", selpar_slope_cPT4_out);
    
	wrt_r_df_col("selpar.A50.rHB1", selpar_A50_rHB1_out);
	wrt_r_df_col("selpar.slope.rHB1", selpar_slope_rHB1_out);
	wrt_r_df_col("selpar.A50.rHB2", selpar_A50_rHB2_out);
	wrt_r_df_col("selpar.slope.rHB2", selpar_slope_rHB2_out);
	wrt_r_df_col("selpar.A50.rHB3", selpar_A50_rHB3_out);
	wrt_r_df_col("selpar.slope.rHB3", selpar_slope_rHB3_out);
	wrt_r_df_col("selpar.A50.rHB4", selpar_A50_rHB4_out);
	wrt_r_df_col("selpar.slope.rHB4", selpar_slope_rHB4_out);
	wrt_r_df_col("selpar.A50.rHB5", selpar_A50_rHB5_out);
	wrt_r_df_col("selpar.slope.rHB5", selpar_slope_rHB5_out);
   
	wrt_r_df_col("selpar.A50.rGN1", selpar_A50_rGN1_out);
	wrt_r_df_col("selpar.slope.rGN1", selpar_slope_rGN1_out);
	wrt_r_df_col("selpar.A50.rGN2", selpar_A50_rGN2_out);
	wrt_r_df_col("selpar.slope.rGN2", selpar_slope_rGN2_out);
	wrt_r_df_col("selpar.A50.rGN3", selpar_A50_rGN3_out);
	wrt_r_df_col("selpar.slope.rGN3", selpar_slope_rGN3_out);
	wrt_r_df_col("selpar.A50.rGN4", selpar_A50_rGN4_out);
	wrt_r_df_col("selpar.slope.rGN4", selpar_slope_rGN4_out);
	wrt_r_df_col("selpar.A50.rGN5", selpar_A50_rGN5_out);
	wrt_r_df_col("selpar.slope.rGN5", selpar_slope_rGN5_out);

    wrt_r_df_col("selpar.age0.rHB.D.logit",selpar_Age0_rHB_D_logit_out);
    wrt_r_df_col("selpar.age1.rHB.D.logit",selpar_Age1_rHB_D_logit_out);
    wrt_r_df_col("selpar.age2.rHB.D.logit",selpar_Age2_rHB_D_logit_out);
    
	wrt_r_df_col("selpar.A50.rHB.D4", selpar_A50_rHD4_out);
	wrt_r_df_col("selpar.slope.rHB.D4", selpar_slope_rHD4_out);
	wrt_r_df_col("selpar.A502.rHB.D4", selpar_A502_rHD4_out);
	wrt_r_df_col("selpar.slope2.rHB.D4", selpar_slope2_rHD4_out);
	wrt_r_df_col("selpar.A50.rHB.D5", selpar_A50_rHD5_out);
	wrt_r_df_col("selpar.slope.rHB.D5", selpar_slope_rHD5_out);
	wrt_r_df_col("selpar.A502.rHB.D5", selpar_A502_rHD5_out);
	wrt_r_df_col("selpar.slope2.rHB.D5", selpar_slope2_rHD5_out);
	
	wrt_r_df_col("log.q.cHL",log_q_cHL_out);
    wrt_r_df_col("log.q.rHB",log_q_rHB_out);
    wrt_r_df_col("log.q.sTV",log_q_sTV_out);
	wrt_r_df_col("log.q.sBT",log_q_sBT_out);
    //wrt_r_df_col("log.q.Vid",log_q_Vid_out);
    wrt_r_df_col("M.constant",M_constant_out);
    wrt_r_df_col("log.avg.F.L.cHL",log_avg_F_cHL_out);
	wrt_r_df_col("log.avg.F.L.cPT",log_avg_F_cPT_out);
	wrt_r_df_col("log.avg.F.L.cTW",log_avg_F_cTW_out);
    wrt_r_df_col("log.avg.F.L.rHB",log_avg_F_rHB_out);
    wrt_r_df_col("log.avg.F.L.rGN",log_avg_F_rGN_out);
    wrt_r_df_col("log.avg.F.D.cGN",log_avg_F_cGN_D_out);
    wrt_r_df_col("log.avg.F.D.rHB",log_avg_F_rHB_D_out);
    wrt_r_df_col("log.avg.F.D.rGN",log_avg_F_rGN_D_out);
    //wrt_r_df_col("F.init",F_init_out);
close_r_df();


// VECTOR object of likelihood contributions
open_r_vector("like");
	wrt_r_item("lk.total", fval);            //weighted likelihood
    wrt_r_item("lk.unwgt.data", fval_data);  //likelihood of just data components
    wrt_r_item("lk.U.sBT", f_sBT_cpue);
    wrt_r_item("lk.U.sTV", f_sTV_cpue);
    wrt_r_item("lk.U.cHL", f_cHL_cpue);
    wrt_r_item("lk.U.rHB", f_rHB_cpue);
	//wrt_r_item("lk.U.Vid", f_Vid_cpue);
    //wrt_r_item("lk.U.rHB.D", f_rHD_cpue);
    wrt_r_item("lk.L.cHL", f_cHL_L);
    wrt_r_item("lk.L.cPT", f_cPT_L);
    wrt_r_item("lk.L.cTW", f_cTW_L);
    wrt_r_item("lk.L.rHB", f_rHB_L);
    wrt_r_item("lk.L.rGN", f_rGN_L);
    wrt_r_item("lk.D.cGN", f_cGN_D);
    wrt_r_item("lk.D.rHB", f_rHB_D);
    wrt_r_item("lk.D.rGN", f_rGN_D);
    wrt_r_item("lk.lenc.sBT", f_sBT_lenc);
    wrt_r_item("lk.lenc.cHL", f_cHL_lenc);
    wrt_r_item("lk.lenc.cPT", f_cPT_lenc);
    wrt_r_item("lk.lenc.rHB", f_rHB_lenc);
    wrt_r_item("lk.lenc.rGN", f_rGN_lenc);
    wrt_r_item("lk.lenc.rHB.D", f_rHB_D_lenc);
    wrt_r_item("lk.agec.sBT", f_sBT_agec);
    wrt_r_item("lk.agec.sTV", f_sTV_agec);
    wrt_r_item("lk.agec.cHL", f_cHL_agec);
    wrt_r_item("lk.agec.cPT", f_cPT_agec);
    wrt_r_item("lk.agec.rHB", f_rHB_agec);
    //wrt_r_item("lk.agec.rGN", f_rGN_agec);
    wrt_r_item("lk.priors",f_priors);
    wrt_r_item("lk.U.RW.cHL",f_cHL_RW_cpue);
    wrt_r_item("lk.U.RW.rHB",f_rHB_RW_cpue);
    //wrt_r_item("lk.U.RW.rHB.D",f_rHD_RW_cpue);


    wrt_r_item("lk.SRfit", f_rec_dev);
    wrt_r_item("lk.SRearly", f_rec_dev_early);
    wrt_r_item("lk.SRend", f_rec_dev_end);
    //wrt_r_item("lk.Ftune", f_Ftune);
    //wrt_r_item("lk.fullF", f_fullF_constraint);
    //wrt_r_item("lk.cvlen.dev", f_cvlen_dev_constraint);
    //wrt_r_item("lk.cvlen.diff", f_cvlen_diff_constraint);

	wrt_r_item("w.L", w_L);
	wrt_r_item("w.D", w_D);
	wrt_r_item("w.lc.sBT", w_lenc_sBT);
	wrt_r_item("w.lc.cHL", w_lenc_cHL);
	wrt_r_item("w.lc.cPT", w_lenc_cPT);
	wrt_r_item("w.lc.rHB", w_lenc_rHB);
	wrt_r_item("w.lc.rHB.D", w_lenc_rHB_D);
	wrt_r_item("w.lc.rGN", w_lenc_rGN);
	wrt_r_item("w.ac.sBT", w_agec_sBT);
	wrt_r_item("w.ac.sTV", w_agec_sTV);
	wrt_r_item("w.ac.cHL", w_agec_cHL);
	wrt_r_item("w.ac.cPT", w_agec_cPT);
	wrt_r_item("w.ac.rHB", w_agec_rHB);
	//wrt_r_item("w.ac.rGN", w_ac_rGN);
	wrt_r_item("w.U.sBT", w_cpue_sBT);
	wrt_r_item("w.U.sTV", w_cpue_sTV);
	//wrt_r_item("w.U.Vid", w_I_Vid);
	wrt_r_item("w.U.cHL", w_cpue_cHL);
	wrt_r_item("w.U.rHB", w_cpue_rHB);
	//wrt_r_item("w.U.rHB.D", w_I_rHD);
	wrt_r_item("w.R", w_rec);
	wrt_r_item("w.R.init", w_rec_early);
	wrt_r_item("w.R.end", w_rec_end);
	//wrt_r_item("w.Ftune.early.phases", w_Ftune);
	//wrt_r_item("w.F.early.phases", w_fullF);
	//wrt_r_item("w.cvlen.dev", w_cvlen_dev);
	//wrt_r_item("w.cvlen.diff", w_cvlen_diff);
close_r_vector();

// VECTOR object of parameters and estimated quantities
open_r_info_list("sel.parms",false);

	wrt_r_item("selpar.A50.sBT", selpar_A50_sBT);
	wrt_r_item("selpar.slope.sBT", selpar_slope_sBT);

	wrt_r_item("selpar.A50.sTV", selpar_A50_sTV);
	wrt_r_item("selpar.slope.sTV", selpar_slope_sTV);
	
	//wrt_r_item("selpar.A50.Vid", selpar_A50_Vid);
	//wrt_r_item("selpar.slope.Vid", selpar_slope_Vid);

	wrt_r_item("selpar.A50.cHL2", selpar_A50_cHL2);
	wrt_r_item("selpar.slope.cHL2", selpar_slope_cHL2);
	wrt_r_item("selpar.A50.cHL3", selpar_A50_cHL3);
	wrt_r_item("selpar.slope.cHL3", selpar_slope_cHL3);
	wrt_r_item("selpar.A50.cHL4", selpar_A50_cHL4);
	wrt_r_item("selpar.slope.cHL4", selpar_slope_cHL4);

	wrt_r_item("selpar.A50.cPT2", selpar_A50_cPT2);
	wrt_r_item("selpar.slope.cPT2", selpar_slope_cPT2);
	wrt_r_item("selpar.A50.cPT3", selpar_A50_cPT3);
	wrt_r_item("selpar.slope.cPT3", selpar_slope_cPT3);
	wrt_r_item("selpar.A50.cPT4", selpar_A50_cPT4);
	wrt_r_item("selpar.slope.cPT4", selpar_slope_cPT4);

	wrt_r_item("selpar.A50.rHB1", selpar_A50_rHB1);
	wrt_r_item("selpar.slope.rHB1", selpar_slope_rHB1);
	wrt_r_item("selpar.A50.rHB2", selpar_A50_rHB2);
	wrt_r_item("selpar.slope.rHB2", selpar_slope_rHB2);
	wrt_r_item("selpar.A50.rHB3", selpar_A50_rHB3);
	wrt_r_item("selpar.slope.rHB3", selpar_slope_rHB3);
	wrt_r_item("selpar.A50.rHB4", selpar_A50_rHB4);
	wrt_r_item("selpar.slope.rHB4", selpar_slope_rHB4);
	wrt_r_item("selpar.A50.rHB5", selpar_A50_rHB5);
	wrt_r_item("selpar.slope.rHB5", selpar_slope_rHB5);

	wrt_r_item("selpar.A50.rGN1", selpar_A50_rGN1);
	wrt_r_item("selpar.slope.rGN1", selpar_slope_rGN1);
	wrt_r_item("selpar.A50.rGN2", selpar_A50_rGN2);
	wrt_r_item("selpar.slope.rGN2", selpar_slope_rGN2);
	wrt_r_item("selpar.A50.rGN3", selpar_A50_rGN3);
	wrt_r_item("selpar.slope.rGN3", selpar_slope_rGN3);
	wrt_r_item("selpar.A50.rGN4", selpar_A50_rGN4);
	wrt_r_item("selpar.slope.rGN4", selpar_slope_rGN4);
	wrt_r_item("selpar.A50.rGN5", selpar_A50_rGN5);
	wrt_r_item("selpar.slope.rGN5", selpar_slope_rGN5);

    wrt_r_item("selpar.Age0.rHB.D", selpar_Age0_rHB_D);
    wrt_r_item("selpar.Age1.rHB.D", selpar_Age1_rHB_D);
    wrt_r_item("selpar.Age2.rHB.D", selpar_Age2_rHB_D);

	wrt_r_df_col("selpar.A50.rHB.D4", selpar_A50_rHD4_out);
	wrt_r_df_col("selpar.slope.rHB.D4", selpar_slope_rHD4_out);
	wrt_r_df_col("selpar.A502.rHB.D4", selpar_A502_rHD4_out);
	wrt_r_df_col("selpar.slope2.rHB.D4", selpar_slope2_rHD4_out);
	wrt_r_df_col("selpar.A50.rHB.D5", selpar_A50_rHD5_out);
	wrt_r_df_col("selpar.slope.rHB.D5", selpar_slope_rHD5_out);
	wrt_r_df_col("selpar.A502.rHB.D5", selpar_A502_rHD5_out);
	wrt_r_df_col("selpar.slope2.rHB.D5", selpar_slope2_rHD5_out);
	
close_r_info_list();

// DATA FRAME of time series deviation vector estimates
// names used in this object must match the names used in the "parm.tvec.cons" object
open_r_df("parm.tvec", styr, endyr, 2);
	wrt_r_namevector(styr,endyr);
	wrt_r_df_col("year", styr,endyr);
    wrt_r_df_col("log.rec.dev", log_rec_dev_out);
    wrt_r_df_col("log.F.dev.cHL", log_F_dev_cHL_out);
	wrt_r_df_col("log.F.dev.cPT", log_F_dev_cPT_out);
	wrt_r_df_col("log.F.dev.cTW", log_F_dev_cTW_out);
    wrt_r_df_col("log.F.dev.rHB", log_F_dev_rHB_out);
    wrt_r_df_col("log.F.dev.rGN", log_F_dev_rGN_out);
    wrt_r_df_col("log.F.dev.cGN.D", log_F_dev_cGN_D_out);
    wrt_r_df_col("log.F.dev.rHB.D", log_F_dev_rHB_D_out);
    wrt_r_df_col("log.F.dev.rGN.D", log_F_dev_rGN_D_out);
	//wrt_r_df_col("log.q.dev.rHB.RWq", q_RW_log_dev_rHB);
	//wrt_r_df_col("log.q.dev.cHL.RWq", q_RW_log_dev_cHL);
	
close_r_df();

// MATRIX object of deviation vector constraints
// names used in this object must match the names used in the "parm.tvec" object
open_r_df("parm.tvec.cons",1,3,2);
    wrt_r_namevector(1,3);
    wrt_r_df_col("log.rec.dev",set_log_dev_rec);
    wrt_r_df_col("log.F.dev.cHL",set_log_dev_F_L_cHL);
	wrt_r_df_col("log.F.dev.cPT",set_log_dev_F_L_cPT);
	wrt_r_df_col("log.F.dev.cTW",set_log_dev_F_L_cTW);
    wrt_r_df_col("log.F.dev.rHB",set_log_dev_F_L_rHB);
    wrt_r_df_col("log.F.dev.rGN",set_log_dev_F_L_rGN);
    wrt_r_df_col("log.F.dev.cGN.D",set_log_dev_F_D_cGN);
    wrt_r_df_col("log.F.dev.rHB.D",set_log_dev_F_D_rHB);
    wrt_r_df_col("log.F.dev.rGN.D",set_log_dev_F_D_rGN);
	//wrt_r_df_col("log.q.dev.rHB.RWq", set_log_dev_RWq);
	//wrt_r_df_col("log.q.dev.cHL.RWq", set_log_dev_RWq);
close_r_df();

// DATA FRAME of age vector deviation estimates
// names used in this object must match the names used in the "parm.avec.cons" object
open_r_df("parm.avec", 1, nages, 2);
	wrt_r_namevector(1,nages);
	wrt_r_df_col("age", agebins); //deviations for first age not estimated
    wrt_r_df_col("log.Nage.dev", log_Nage_dev_output);
close_r_df();

// MATRIX object of age vector deviation constraints
// names used in this object must match the names used in the "parm.avec" object
open_r_df("parm.avec.cons",1,3,2);
    wrt_r_namevector(1,3);
    wrt_r_df_col("log.Nage.dev",set_log_dev_Nage);
close_r_df();


    open_r_matrix("N.age");
    wrt_r_matrix(N, 2, 2);
    wrt_r_namevector(styr, endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("N.age.mdyr");
    wrt_r_matrix(N_mdyr, 2, 2);
    wrt_r_namevector(styr, endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("N.age.spawn");
    wrt_r_matrix(N_spawn, 2, 2);
    wrt_r_namevector(styr, endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("B.age");
    wrt_r_matrix(B, 2, 2);
    wrt_r_namevector(styr, endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("F.age");
    wrt_r_matrix(F, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("Z.age");
    wrt_r_matrix(Z, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("L.age.pred.knum");
    wrt_r_matrix(L_total_num/1000.0, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("L.age.pred.klb");
    wrt_r_matrix(L_total_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();


// LIST object with annual selectivity at age by fishery

open_r_list("size.age.fishery");

    open_r_matrix("len.cHL.mm");
    wrt_r_matrix(len_cHL_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("len.cPT.mm");
    wrt_r_matrix(len_cPT_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("len.cTW.mm");
    wrt_r_matrix(len_cTW_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("len.rHB.mm");
    wrt_r_matrix(len_rHB_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("len.rGN.mm");
    wrt_r_matrix(len_rGN_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("len.cGN.D.mm");
    wrt_r_matrix(len_cGN_D_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("len.rHB.D.mm");
    wrt_r_matrix(len_rHB_D_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("len.rGN.D.mm");
    wrt_r_matrix(len_rGN_D_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();


    open_r_matrix("wgt.cHL.klb");
    wrt_r_matrix(wgt_cHL_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("wgt.cPT.klb");
    wrt_r_matrix(wgt_cPT_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("wgt.cTW.klb");
    wrt_r_matrix(wgt_cTW_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("wgt.rHB.klb");
    wrt_r_matrix(wgt_rHB_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("wgt.rGN.klb");
    wrt_r_matrix(wgt_rGN_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("wgt.cGN.D.klb");
    wrt_r_matrix(wgt_cGN_D_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("wgt.rHB.D.klb");
    wrt_r_matrix(wgt_rHB_D_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("wgt.rGN.D.klb");
    wrt_r_matrix(wgt_rGN_D_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

close_r_list();

open_r_list("sel.age");

    wrt_r_complete_vector("sel.v.wgted.L",sel_wgted_L, agebins);
    wrt_r_complete_vector("sel.v.wgted.D",sel_wgted_D, agebins);
    wrt_r_complete_vector("sel.v.wgted.tot",sel_wgted_tot, agebins);

    wrt_r_complete_vector("sel.v.sBT",sel_sBT_vec, agebins);
    wrt_r_complete_vector("sel.v.sTV",sel_sTV_vec, agebins);


    open_r_matrix("sel.m.cHL");
    wrt_r_matrix(sel_cHL, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("sel.m.cPT");
    wrt_r_matrix(sel_cPT, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("sel.m.cTW");
    wrt_r_matrix(sel_cTW, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("sel.m.rHB");
    wrt_r_matrix(sel_rHB, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("sel.m.rGN");
    wrt_r_matrix(sel_rGN, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("sel.m.cGN.D");
    wrt_r_matrix(sel_cGN_D, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("sel.m.rHB.D");
    wrt_r_matrix(sel_rHB_D, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("sel.m.rGN.D");
    wrt_r_matrix(sel_rGN_D, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

close_r_list();



//LIST object with predicted and observed composition data
open_r_list("comp.mats");
    open_r_matrix("lcomp.sBT.ob");  //Excluded because the blackfish trap comp is a single year
    wrt_r_matrix(obs_lenc_sBT, 2, 2);
    wrt_r_namevector(yrs_lenc_sBT);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.sBT.pr");
    wrt_r_matrix(pred_sBT_lenc, 2, 2);
    wrt_r_namevector(yrs_lenc_sBT);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.cHL.ob");
    wrt_r_matrix(obs_lenc_cHL, 2, 2);
    wrt_r_namevector(yrs_lenc_cHL);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.cHL.pr");
    wrt_r_matrix(pred_cHL_lenc, 2, 2);
    wrt_r_namevector(yrs_lenc_cHL);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.cPT.ob");
    wrt_r_matrix(obs_lenc_cPT, 2, 2);
    wrt_r_namevector(yrs_lenc_cPT);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.cPT.pr");
    wrt_r_matrix(pred_cPT_lenc, 2, 2);
    wrt_r_namevector(yrs_lenc_cPT);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.rHB.ob");
    wrt_r_matrix(obs_lenc_rHB, 2, 2);
    wrt_r_namevector(yrs_lenc_rHB);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.rHB.pr");
    wrt_r_matrix(pred_rHB_lenc, 2, 2);
    wrt_r_namevector(yrs_lenc_rHB);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.rGN.ob");
    wrt_r_matrix(obs_lenc_rGN, 2, 2);
    wrt_r_namevector(yrs_lenc_rGN);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.rGN.pr");
    wrt_r_matrix(pred_rGN_lenc, 2, 2);
    wrt_r_namevector(yrs_lenc_rGN);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.rHB.D.ob");
    wrt_r_matrix(obs_lenc_rHB_D, 2, 2);
    wrt_r_namevector(yrs_lenc_rHB_D);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.rHB.D.pr");
    wrt_r_matrix(pred_rHB_D_lenc, 2, 2);
    wrt_r_namevector(yrs_lenc_rHB_D);
    wrt_r_namevector(lenbins);
    close_r_matrix();
    //-----------------------------------

    open_r_matrix("acomp.sBT.ob");
    wrt_r_matrix(obs_agec_sBT, 2, 2);
    wrt_r_namevector(yrs_agec_sBT);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("acomp.sBT.pr");
    wrt_r_matrix(pred_sBT_agec, 2, 2);
    wrt_r_namevector(yrs_agec_sBT);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("acomp.sTV.ob");
    wrt_r_matrix(obs_agec_sTV, 2, 2);
    wrt_r_namevector(yrs_agec_sTV);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("acomp.sTV.pr");
    wrt_r_matrix(pred_sTV_agec, 2, 2);
    wrt_r_namevector(yrs_agec_sTV);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("acomp.cHL.ob");
    wrt_r_matrix(obs_agec_cHL, 2, 2);
    wrt_r_namevector(yrs_agec_cHL);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("acomp.cHL.pr");
    wrt_r_matrix(pred_cHL_agec, 2, 2);
    wrt_r_namevector(yrs_agec_cHL);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("acomp.cPT.ob");
    wrt_r_matrix(obs_agec_cPT, 2, 2);
    wrt_r_namevector(yrs_agec_cPT);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("acomp.cPT.pr");
    wrt_r_matrix(pred_cPT_agec, 2, 2);
    wrt_r_namevector(yrs_agec_cPT);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("acomp.rHB.ob");
    wrt_r_matrix(obs_agec_rHB, 2, 2);
    wrt_r_namevector(yrs_agec_rHB);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("acomp.rHB.pr");
    wrt_r_matrix(pred_rHB_agec, 2, 2);
    wrt_r_namevector(yrs_agec_rHB);
    wrt_r_namevector(agebins);
    close_r_matrix();

    //open_r_matrix("acomp.rGN.ob");
    //wrt_r_matrix(obs_rGN_agec, 2, 2);
    //wrt_r_namevector(yrs_rGN_agec);
    //wrt_r_namevector(agebins);
    //close_r_matrix();
    //
    //open_r_matrix("acomp.rGN.pr");
    //wrt_r_matrix(pred_rGN_agec, 2, 2);
    //wrt_r_namevector(yrs_rGN_agec);
    //wrt_r_namevector(agebins);
    //close_r_matrix();
close_r_list();

// VECTOR object of SDNR calculations
open_r_vector("sdnr");
    wrt_r_item("sdnr.U.cHL", sdnr_I_cHL);
	//wrt_r_item("sdnr.U.Vid", sdnr_I_Vid);
    wrt_r_item("sdnr.U.rHB", sdnr_I_rHB);
    wrt_r_item("sdnr.U.sTV", sdnr_I_sTV);
	wrt_r_item("sdnr.U.sBT", sdnr_I_sBT);
	 
    //wrt_r_item("sdnr.lc.sBT", sdnr_lc_sBT);
    //wrt_r_item("sdnr.lc.sTV", sdnr_lc_sTV);	
	//wrt_r_item("sdnr.lc.cHL", sdnr_lc_cHL);
	//wrt_r_item("sdnr.lc.cPT", sdnr_lc_cPT);	
    //wrt_r_item("sdnr.lc.rHB", sdnr_lc_rHB);
	//wrt_r_item("sdnr.lc.rHB.D", sdnr_lc_rHB_D);
	//wrt_r_item("sdnr.lc.rGN", sdnr_lc_rGN);
	//
    //wrt_r_item("sdnr.ac.sBT", sdnr_ac_sBT);
    //wrt_r_item("sdnr.ac.sTV", sdnr_ac_sTV);
	//wrt_r_item("sdnr.ac.cHL", sdnr_ac_cHL);
	//wrt_r_item("sdnr.ac.rHB", sdnr_ac_rHB);
	//wrt_r_item("sdnr.ac.Vid", sdnr_ac_Vid);
	
 close_r_vector();

// DATA FRAME of time series
open_r_df("t.series", styr, endyr, 2);
	wrt_r_namevector(styr,endyr);
	wrt_r_df_col("year", styr,endyr);
	wrt_r_df_col("F.Fmsy", FdF_msy);
	wrt_r_df_col("F.full", Fapex);
	wrt_r_df_col("F.cHL", F_cHL_out);
	wrt_r_df_col("F.cPT", F_cPT_out);
	wrt_r_df_col("F.cTW", F_cTW_out);
	wrt_r_df_col("F.rHB", F_rHB_out);
	wrt_r_df_col("F.rGN", F_rGN_out);
	wrt_r_df_col("F.cGN.D", F_cGN_D_out);
	wrt_r_df_col("F.rHB.D", F_rHB_D_out);
	wrt_r_df_col("F.rGN.D", F_rGN_D_out);
	wrt_r_df_col("Fsum", Fsum);
    wrt_r_df_col("N", totN);   //abundance at start of year
    wrt_r_df_col("recruits", rec);
    wrt_r_df_col("logR.dev", log_rec_dev_output); //places zeros in yrs deviations not estimated
    wrt_r_df_col("SSB", SSB);
    wrt_r_df_col("SSB.SSBmsy", SdSSB_msy);
    wrt_r_df_col("SSB.msst", SdSSB_msy/(1.0-M_constant));
    wrt_r_df_col("B", totB);
    wrt_r_df_col("B.B0", totB/B0);
    wrt_r_df_col("MatFemB", MatFemB);
    wrt_r_df_col("SPR.static", spr_static);

    wrt_r_df_col("total.L.klb", L_total_klb_yr);
    wrt_r_df_col("total.L.knum", L_total_knum_yr);
    wrt_r_df_col("total.D.klb",D_total_klb_yr);
    wrt_r_df_col("total.D.knum",D_total_knum_yr);

    wrt_r_df_col("U.sBT.ob", obs_cpue_sBT);
    wrt_r_df_col("U.sBT.pr", pred_sBT_cpue);
    wrt_r_df_col("cv.U.sBT", obs_cv_cpue_sBT);
    wrt_r_df_col("U.sTV.ob", obs_cpue_sTV);
    wrt_r_df_col("U.sTV.pr", pred_sTV_cpue);
    wrt_r_df_col("cv.U.sTV", obs_cv_cpue_sTV);
	//wrt_r_df_col("U.Vid.ob", obs_Vid_cpue);
    //wrt_r_df_col("U.Vid.pr", pred_Vid_cpue);
    //wrt_r_df_col("cv.U.Vid", Vid_cpue_cv);
    wrt_r_df_col("U.cHL.ob", obs_cpue_cHL);
    wrt_r_df_col("U.cHL.pr", pred_cHL_cpue);
    wrt_r_df_col("cv.U.cHL", obs_cv_cpue_cHL);
    wrt_r_df_col("U.rHB.ob", obs_cpue_rHB);
    wrt_r_df_col("U.rHB.pr", pred_rHB_cpue);
    wrt_r_df_col("cv.U.rHB", obs_cv_cpue_rHB);
    //wrt_r_df_col("U.rHB.D.ob", obs_rHD_cpue);
    //wrt_r_df_col("U.rHB.D.pr", pred_rHD_cpue);
    //wrt_r_df_col("cv.U.rHB.D", rHD_cpue_cv);

    wrt_r_df_col("q.cHL", q_cHL);
    wrt_r_df_col("q.rHB", q_rHB);
    
    wrt_r_df_col("q.cHL.rate.mult",q_rate_fcn_cHL);
    wrt_r_df_col("q.rHB.rate.mult",q_rate_fcn_rHB);
    //wrt_r_df_col("q.rHB.D.rate.mult",q_rate_fcn_rHD);
    wrt_r_df_col("q.DD.mult", q_DD_fcn);
    wrt_r_df_col("q.DD.B.4plus", B_q_DD);
    wrt_r_df_col("q.cHL.RW.log.dev",q_RW_log_dev_cHL);
    wrt_r_df_col("q.rHB.RW.log.dev",q_RW_log_dev_rHB);
    //wrt_r_df_col("q.rHB.D.RW.log.dev",q_RW_log_dev_rHD);

    wrt_r_df_col("L.cHL.ob", obs_L_cHL);
    wrt_r_df_col("L.cHL.pr", pred_cHL_L_klb);
    wrt_r_df_col("cv.L.cHL", obs_cv_L_cHL);
    wrt_r_df_col("L.cPT.ob", obs_L_cPT);
    wrt_r_df_col("L.cPT.pr", pred_cPT_L_klb);
    wrt_r_df_col("cv.L.cPT", obs_cv_L_cPT);
    wrt_r_df_col("L.cTW.ob", obs_L_cTW);
    wrt_r_df_col("L.cTW.pr", pred_cTW_L_klb);
    wrt_r_df_col("cv.L.cTW", obs_cv_L_cTW);
    wrt_r_df_col("L.rHB.ob", obs_L_rHB);
    wrt_r_df_col("L.rHB.pr", pred_rHB_L_klb);
    wrt_r_df_col("cv.L.rHB", obs_cv_L_rHB);
    wrt_r_df_col("L.rGN.ob", obs_L_rGN);
    wrt_r_df_col("L.rGN.pr", pred_rGN_L_klb);
    wrt_r_df_col("cv.L.rGN", obs_cv_L_rGN);

    wrt_r_df_col("D.cGN.ob", obs_cGN_D);
    wrt_r_df_col("D.cGN.pr", pred_cGN_D_knum);
    wrt_r_df_col("cv.D.cGN", cGN_D_cv);
    wrt_r_df_col("D.rHB.ob", obs_rHB_D);
    wrt_r_df_col("D.rHB.pr", pred_rHB_D_knum);
    wrt_r_df_col("cv.D.rHB", obs_cv_D_rHB);
    wrt_r_df_col("D.rGN.ob", obs_rGN_D);
    wrt_r_df_col("D.rGN.pr", pred_rGN_D_knum);
    wrt_r_df_col("cv.D.rGN", obs_cv_D_rGN);

    wrt_r_df_col("cHL.D.release", obs_released_cHL);
    wrt_r_df_col("cHL.D.closed.release", obs_released_cHL_closed);
    wrt_r_df_col("cPT.D.release", obs_released_cPT);
    wrt_r_df_col("cPT.D.closed.release", obs_released_cPT_closed);
    wrt_r_df_col("rHB.D.release", obs_released_rHB);
    wrt_r_df_col("rGN.D.release", obs_released_rGN);


    //comp sample sizes
    wrt_r_df_col("lcomp.sBT.n", nsamp_sBT_lenc_allyr);
    wrt_r_df_col("lcomp.cHL.n", nsamp_cHL_lenc_allyr);
    wrt_r_df_col("lcomp.cPT.n", nsamp_cPT_lenc_allyr);
    wrt_r_df_col("lcomp.rHB.n", nsamp_rHB_lenc_allyr);
    wrt_r_df_col("lcomp.rGN.n", nsamp_rGN_lenc_allyr);
    wrt_r_df_col("lcomp.rHB.D.n", nsamp_rHB_D_lenc_allyr);

    wrt_r_df_col("acomp.sBT.n", nsamp_sBT_agec_allyr);
    wrt_r_df_col("acomp.sTV.n", nsamp_sTV_agec_allyr);
    wrt_r_df_col("acomp.cHL.n", nsamp_cHL_agec_allyr);
    wrt_r_df_col("acomp.cPT.n", nsamp_cPT_agec_allyr);
    wrt_r_df_col("acomp.rHB.n", nsamp_rHB_agec_allyr);
    //wrt_r_df_col("acomp.rGN.n", nsamp_rGN_agec_allyr);

    wrt_r_df_col("lcomp.sBT.nfish", nfish_sBT_lenc_allyr);
    wrt_r_df_col("lcomp.cHL.nfish", nfish_cHL_lenc_allyr);
    wrt_r_df_col("lcomp.cPT.nfish", nfish_cPT_lenc_allyr);
    wrt_r_df_col("lcomp.rHB.nfish", nfish_rHB_lenc_allyr);
    wrt_r_df_col("lcomp.rGN.nfish", nfish_rGN_lenc_allyr);
    wrt_r_df_col("lcomp.rHB.D.nfish", nfish_rHB_D_lenc_allyr);

    wrt_r_df_col("acomp.sBT.nfish", nfish_sBT_agec_allyr);
    wrt_r_df_col("acomp.sTV.nfish", nfish_sTV_agec_allyr);
    wrt_r_df_col("acomp.cHL.nfish", nfish_cHL_agec_allyr);
    wrt_r_df_col("acomp.cPT.nfish", nfish_cPT_agec_allyr);
    wrt_r_df_col("acomp.rHB.nfish", nfish_rHB_agec_allyr);
    //wrt_r_df_col("acomp.rGN.nfish", nfish_rGN_agec_allyr);

close_r_df();

// DATA FRAME of L and D time series by fishery
open_r_df("LD.pr.tseries", styr, (endyr+1), 2);
	wrt_r_namevector(styr,(endyr+1));
	wrt_r_df_col("year", styr,(endyr+1));

    wrt_r_df_col("L.cHL.klb", pred_cHL_L_klb);
    wrt_r_df_col("L.cHL.knum", pred_cHL_L_knum);
    wrt_r_df_col("L.cPT.klb", pred_cPT_L_klb);
    wrt_r_df_col("L.cPT.knum", pred_cPT_L_knum);
    wrt_r_df_col("L.cTW.klb", pred_cTW_L_klb);
    wrt_r_df_col("L.cTW.knum", pred_cTW_L_knum);
    wrt_r_df_col("L.rHB.klb", pred_rHB_L_klb);
    wrt_r_df_col("L.rHB.knum", pred_rHB_L_knum);
    wrt_r_df_col("L.rGN.klb", pred_rGN_L_klb);
    wrt_r_df_col("L.rGN.knum", pred_rGN_L_knum);

    wrt_r_df_col("D.cGN.klb", pred_cGN_D_klb);
    wrt_r_df_col("D.cGN.knum", pred_cGN_D_knum);
    wrt_r_df_col("D.rHB.klb", pred_rHB_D_klb);
    wrt_r_df_col("D.rHB.knum", pred_rHB_D_knum);
    wrt_r_df_col("D.rGN.klb", pred_rGN_D_klb);
    wrt_r_df_col("D.rGN.knum", pred_rGN_D_knum);
close_r_df();


open_r_df("a.series", 1, nages, 2);
	wrt_r_namevector(agebins);
	wrt_r_df_col("age", agebins);
	wrt_r_df_col("length", meanlen_TL);
	wrt_r_df_col("length.cv", len_cv);
	wrt_r_df_col("length.sd", len_sd);
	wrt_r_df_col("weight", wgt_lb);     //for FishGraph
	wrt_r_df_col("wgt.klb", wgt_klb);
	wrt_r_df_col("wgt.mt", wgt_mt);
	wrt_r_df_col("wgt.wgted.L.klb", wgt_wgted_L_klb);
	wrt_r_df_col("wgt.wgted.D.klb", wgt_wgted_D_klb);
	wrt_r_df_col("prop.female", prop_f);
	wrt_r_df_col("mat.female", maturity_f);
	wrt_r_df_col("fecundity.scaled", fecundity);
	wrt_r_df_col("reprod", reprod);
	wrt_r_df_col("reprod.MatFemB", reprod2);
	wrt_r_df_col("M", M);
	wrt_r_df_col("F.initial", F_initial);
	wrt_r_df_col("Z.initial", Z_initial);
	wrt_r_df_col("Nage.eq.init",N_initial_eq);
	wrt_r_df_col("log.Nage2plus.init.dev",log_Nage_dev_output);
close_r_df();

open_r_df("eq.series", 1, n_iter_msy, 2);
	wrt_r_namevector(1,n_iter_msy);
	wrt_r_df_col("F.eq", F_msy);
	wrt_r_df_col("spr.eq", spr_msy);
	wrt_r_df_col("R.eq", R_eq);
	wrt_r_df_col("SSB.eq", SSB_eq);
	wrt_r_df_col("B.eq", B_eq);
	wrt_r_df_col("L.eq.klb", L_eq_klb);
	wrt_r_df_col("D.eq.knum", D_eq_knum);
close_r_df();

open_r_df("pr.series", 1, n_iter_spr, 2);
	wrt_r_namevector(1,n_iter_spr);
	wrt_r_df_col("F.spr", F_spr);
	wrt_r_df_col("spr", spr_spr);
	wrt_r_df_col("SPR", spr_spr/spr_F0);
	wrt_r_df_col("ypr.lb", L_spr);
close_r_df();


open_r_list("CLD.est.mats");

    open_r_matrix("Lw.cHL");
        wrt_r_matrix(L_cHL_klb, 1,1);
    close_r_matrix();

    open_r_matrix("Lw.cPT");
        wrt_r_matrix(L_cPT_klb, 1,1);
    close_r_matrix();

    open_r_matrix("Lw.cTW");
        wrt_r_matrix(L_cTW_klb, 1,1);
    close_r_matrix();

    open_r_matrix("Lw.rHB");
        wrt_r_matrix(L_rHB_klb, 1,1);
    close_r_matrix();

    open_r_matrix("Lw.rGN");
        wrt_r_matrix(L_rGN_klb, 1,1);
    close_r_matrix();

    open_r_matrix("Lw.total");
        wrt_r_matrix(L_total_klb, 1,1);
    close_r_matrix();

    open_r_matrix("Dw.cGN");
        wrt_r_matrix(D_cGN_klb, 1,1);
    close_r_matrix();

    open_r_matrix("Dw.rHB");
        wrt_r_matrix(D_rHB_klb, 1,1);
    close_r_matrix();

    open_r_matrix("Dw.rGN");
        wrt_r_matrix(D_rGN_klb, 1,1);
    close_r_matrix();

    open_r_matrix("Dw.total");
        wrt_r_matrix(D_total_klb, 1,1);
    close_r_matrix();


    open_r_matrix("Ln.cHL");
        wrt_r_matrix(L_cHL_num, 1,1);
    close_r_matrix();

    open_r_matrix("Ln.cPT");
        wrt_r_matrix(L_cPT_num, 1,1);
    close_r_matrix();

    open_r_matrix("Ln.cTW");
        wrt_r_matrix(L_cTW_num, 1,1);
    close_r_matrix();

    open_r_matrix("Ln.rHB");
        wrt_r_matrix(L_rHB_num, 1,1);
    close_r_matrix();

    open_r_matrix("Ln.rGN");
        wrt_r_matrix(L_rGN_num, 1,1);
    close_r_matrix();

    open_r_matrix("Ln.total");
        wrt_r_matrix(L_total_num, 1,1);
    close_r_matrix();

    open_r_matrix("Dn.cGN");
        wrt_r_matrix(D_cGN_num,1,1);
    close_r_matrix();

    open_r_matrix("Dn.rHB");
        wrt_r_matrix(D_rHB_num, 1,1);
    close_r_matrix();

    open_r_matrix("Dn.rGN");
        wrt_r_matrix(D_rGN_num, 1,1);
    close_r_matrix();

    open_r_matrix("Dn.total");
        wrt_r_matrix(D_total_num, 1,1);
    close_r_matrix();

close_r_list();



close_r_file();
