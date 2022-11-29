// Create a file with an R object from AD Model Builder

// open the file using the default AD Model Builder file name, and
// 6 digits of precision
open_r_file(adprogram_name + ".rdat", 6);

// Example of an INFO object
open_r_info_list("info", true);
	wrt_r_item("title", "SEDAR 58 Assessment");
	wrt_r_item("species", "SA Cobia");
	wrt_r_item("model", "Statistical Catch at Age");
	if(SR_switch==1)
	    {wrt_r_item("rec.model", "BH-steep");}
    if(SR_switch==2)
	    {wrt_r_item("rec.model", "Ricker-steep");}
	wrt_r_item("base.run", "");
	wrt_r_item("units.length", "mm");
	wrt_r_item("units.weight", "lb");
	wrt_r_item("units.biomass", "metric tons");
	wrt_r_item("units.ssb", "metric tons");
	wrt_r_item("units.ypr", "lb.whole");
	wrt_r_item("units.landings", "1000 lb whole");
	wrt_r_item("units.discards", "NA");
    wrt_r_item("units.numbers", "number fish");
    wrt_r_item("units.naa", "number fish");
	wrt_r_item("units.rec", "number fish");
	if(SR_switch==3)
   {wrt_r_item("rec.model", "Mean");}
close_r_info_list();


// VECTOR object of parameters and estimated quantities
open_r_info_list("parms", false);
	wrt_r_item("styr", styr);
	wrt_r_item("endyr", endyr);
	wrt_r_item("styrR", styr_rec_dev);
	wrt_r_item("M.constant", M_constant);
	wrt_r_item("smsy2msstM", smsy2msstM);
    wrt_r_item("smsy2msst75", smsy2msst75);
	wrt_r_item("Linf", Linf);
	wrt_r_item("K", K);
	wrt_r_item("t0", t0);
	wrt_r_item("len_cv_val",len_cv_val);
	wrt_r_item("Linf_L", Linf_L);
	wrt_r_item("K_L", K_L);
	wrt_r_item("t0_L", t0_L);
	wrt_r_item("len_cv_val_L",len_cv_val_L);
	wrt_r_item("Linf_F", Linf_F);
	wrt_r_item("K_F", K_F);
	wrt_r_item("t0_F", t0_F);
	wrt_r_item("len_cv_val_F",len_cv_val_F);
	wrt_r_item("wgt.a", wgtpar_a);
	wrt_r_item("wgt.b", wgtpar_b);
	wrt_r_item("spawn.time", spawn_time_frac);
	//wrt_r_item("q.cGN", mfexp(log_q_cGN));
	//wrt_r_item("q.cLL", mfexp(log_q_cLL));	
	wrt_r_item("q.rHB", mfexp(log_q_cpue_rHB));

    wrt_r_item("q.rate",q_rate);
	wrt_r_item("q.beta",q_DD_beta);
	wrt_r_item("q.DD.B0.exploitable", B0_q_DD);
	wrt_r_item("F.init", F_init);
	wrt_r_item("F.prop.cGN", F_cGN_prop);
	//wrt_r_item("F.prop.cLL", F_cLL_prop);	
	wrt_r_item("F.prop.rGN", F_rGN_prop);
	wrt_r_item("Fend.mean", Fend_mean);
	
	if(SR_switch==3)
    {
      wrt_r_item("Mean.biascorr",BiasCor);
      wrt_r_item("Mean.R0", R0);
    }

	wrt_r_item("B0", B0);
	wrt_r_item("Bstyr.B0", totB(styr)/B0);
	wrt_r_item("SSB0", S0);
	wrt_r_item("SSBstyr.SSB0", SSB(styr)/S0);
	wrt_r_item("Rstyr.R0", rec(styr)/R0);
	if(SR_switch==1)
	{
      wrt_r_item("BH.biascorr",BiasCor);
	  wrt_r_item("BH.Phi0", spr_F0);
	  wrt_r_item("BH.R0", R0);
	  wrt_r_item("BH.steep", steep);
    }
    if(SR_switch==2)
	{
      wrt_r_item("Ricker.biascorr",BiasCor);
	  wrt_r_item("Ricker.Phi0", spr_F0);
	  wrt_r_item("Ricker.R0", R0);
	  wrt_r_item("Ricker.steep", steep);
    }
	wrt_r_item("R.sigma.logdevs", sigma_rec_dev);
	wrt_r_item("R.sigma.par",rec_sigma);
	wrt_r_item("R.autocorr",R_autocorr);
	wrt_r_item("R0", R0); //same as BH.R0, but used in BSR.time.plots
	wrt_r_item("R.virgin.bc", R_virgin); //bias-corrected virgin recruitment
	wrt_r_item("rec.lag", 1.0);
	wrt_r_item("msy.klb", msy_klb_out);
	wrt_r_item("msy.knum", msy_knum_out);
	wrt_r_item("Fmsy", F_msy_out);
	wrt_r_item("SSBmsy", SSB_msy_out);
	//wrt_r_item("msst.msy", smsy2msst75*SSB_msy_out);
	wrt_r_item("msst", smsy2msst75*SSB_F40_out);
	wrt_r_item("Bmsy", B_msy_out);
	wrt_r_item("Rmsy", R_msy_out);
	wrt_r_item("sprmsy",spr_msy_out);
	wrt_r_item("Fend.Fmsy", FdF_msy_end);
	wrt_r_item("Fend.Fmsy.mean", FdF_msy_end_mean);
	wrt_r_item("SSBend.SSBmsy", SdSSB_msy_end);
	wrt_r_item("SSBend.MSST", SdSSB_msy_end/smsy2msst75);
	wrt_r_item("F30", F30_out); //For FishGraph
	wrt_r_item("SSB.F30", SSB_F30_out); //For FishGraph
	wrt_r_item("msst.F30", smsy2msst75*SSB_F30_out);
	//wrt_r_item("msst", smsy2msst75*SSB_F30_out); //For FishGraph
	 wrt_r_item("B.F30", B_F30_out);
	wrt_r_item("R.F30", R_F30_out);
    wrt_r_item("L.F30.knum", L_F30_knum_out);
    wrt_r_item("L.F30.klb", L_F30_klb_out);
	wrt_r_item("F40", F40_out); //For FishGraph
	wrt_r_item("SSB.F40", SSB_F40_out); //For FishGraph
	wrt_r_item("msst.F40", smsy2msst75*SSB_F40_out);
	//wrt_r_item("msst", smsy2msst75*SSB_F30_out); //For FishGraph
	 wrt_r_item("B.F40", B_F40_out);
	wrt_r_item("R.F40", R_F40_out);
    wrt_r_item("L.F40.knum", L_F40_knum_out);
    wrt_r_item("L.F40.klb", L_F40_klb_out);
 close_r_info_list();
	
 // VECTOR object of parameters and estimated quantities
open_r_info_list("spr.brps", false);
	wrt_r_item("F20", F20_out);
	wrt_r_item("F30", F30_out);
	wrt_r_item("F40", F40_out);
	wrt_r_item("SSB.F30", SSB_F30_out);
    wrt_r_item("msst.F30", smsy2msst75*SSB_F30_out);   
    wrt_r_item("B.F30", B_F30_out);
	wrt_r_item("R.F30", R_F30_out);
    wrt_r_item("L.F30.knum", L_F30_knum_out);
    wrt_r_item("L.F30.klb", L_F30_klb_out);       	 	
	wrt_r_item("Fend.F30.mean", FdF30_end_mean);
    wrt_r_item("SSBend.SSBF30", SdSSB_F30_end);
	wrt_r_item("SSBend.MSSTF30", Sdmsst_F30_end);	
	
	wrt_r_item("SSB.F40", SSB_F40_out);
    wrt_r_item("msst.F40", smsy2msst75*SSB_F40_out);   
    wrt_r_item("B.F40", B_F40_out);
	wrt_r_item("R.F40", R_F40_out);
    wrt_r_item("L.F40.knum", L_F40_knum_out);
    wrt_r_item("L.F40.klb", L_F40_klb_out);       	 	
	wrt_r_item("Fend.F40.mean", FdF40_end_mean);
    wrt_r_item("SSBend.SSBF40", SdSSB_F40_end);
	wrt_r_item("SSBend.MSSTF40", Sdmsst_F40_end);
 close_r_info_list();

// MATRIX object of parameter constraints
open_r_df("parm.cons",1,8,2);
    wrt_r_namevector(1,8);
    wrt_r_df_col("Linf",Linf_out);
    wrt_r_df_col("K",K_out);
    wrt_r_df_col("t0",t0_out);
    wrt_r_df_col("len_cv_val",len_cv_val_out);
	wrt_r_df_col("Linf_L", Linf_L_out);
	wrt_r_df_col("K_L", K_L_out);
	wrt_r_df_col("t0_L", t0_L_out);
	wrt_r_df_col("len_cv_val_L",len_cv_val_L_out);
	wrt_r_df_col("Linf_F", Linf_F_out);
	wrt_r_df_col("K_F", K_F_out);
	wrt_r_df_col("t0_F", t0_F_out);
	wrt_r_df_col("len_cv_val_F",len_cv_val_F_out);
    wrt_r_df_col("log_R0",log_R0_out);
    wrt_r_df_col("steep",steep_out);
    wrt_r_df_col("rec_sigma",rec_sigma_out);
    wrt_r_df_col("R_autocorr",R_autocorr_out);
    wrt_r_df_col("log_dm_lenc_cGN",log_dm_cGN_lc_out);
	wrt_r_df_col("log_dm_agec_rGN",log_dm_rGN_ac_out);

    wrt_r_df_col("selpar_A50_cGN1",selpar_A50_cGN1_out);
    wrt_r_df_col("selpar_slope_cGN1",selpar_slope_cGN1_out);
	
	wrt_r_df_col("selpar_A50_rGN1",selpar_A50_rGN1_out);
    wrt_r_df_col("selpar_slope_rGN1",selpar_slope_rGN1_out);
	wrt_r_df_col("selpar_A50_rGN2",selpar_A50_rGN2_out);
    wrt_r_df_col("selpar_slope_rGN2",selpar_slope_rGN2_out);
	
	wrt_r_df_col("log_q_cpue_rHB",log_q_rHB_out);
    wrt_r_df_col("M_constant",M_constant_out);
    wrt_r_df_col("log_avg_F_L_cGN",log_avg_F_cGN_out);
	wrt_r_df_col("log_avg_F_L_rGN",log_avg_F_rGN_out);
	
    wrt_r_df_col("F_init",F_init_out);
close_r_df();

// DATA FRAME of time series deviation vector estimates
// names used in this object must match the names used in the "parm.tvec.cons" object
open_r_df("parm.tvec", styr, endyr, 2);
	wrt_r_namevector(styr,endyr);
	wrt_r_df_col("year", styr,endyr);
    wrt_r_df_col("log.rec.dev", log_rec_dev_out);
	
    wrt_r_df_col("log.F.dev.cGN", log_F_dev_cGN_out);
	wrt_r_df_col("log.F.dev.rGN", log_F_dev_rGN_out);

	wrt_r_df_col("log.q.dev.rHB.RWq", q_RW_log_dev_rHB);
	
close_r_df();

// MATRIX object of deviation vector constraints
// names used in this object must match the names used in the "parm.tvec" object
open_r_df("parm.tvec.cons",1,3,2);
    wrt_r_namevector(1,3);
    wrt_r_df_col("log.rec.dev",set_log_dev_rec);
    wrt_r_df_col("log.F.dev.cGN",set_log_dev_F_L_cGN);
	wrt_r_df_col("log.F.dev.rGN",set_log_dev_F_L_rGN);

	wrt_r_df_col("log.q.dev.rHB.RWq", set_log_dev_RWq);
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

// VECTOR object of SDNR calculations
open_r_vector("sdnr");
    wrt_r_item("sdnr.U.rHB", sdnr_I_rHB);
	    
 close_r_vector();

// VECTOR object of likelihood contributions
open_r_vector("like");
    wrt_r_item("lk.total", fval); //weighted likelihood
    wrt_r_item("lk.unwgt.data", fval_data);  //likelihood of just data components

    wrt_r_item("lk.L.cGN", f_cGN_L);
    wrt_r_item("lk.L.rGN", f_rGN_L);

    wrt_r_item("lk.U.rHB", f_rHB_cpue);
	
    wrt_r_item("lk.lenc.cGN", f_cGN_lenc);
    wrt_r_item("lk.agec.rGN", f_rGN_agec);
	
    wrt_r_item("lk.Nage.init", f_Nage_init);
    wrt_r_item("lk.SRfit", f_rec_dev);
    wrt_r_item("lk.SRearly", f_rec_dev_early);
    wrt_r_item("lk.SRend", f_rec_dev_end);
    wrt_r_item("lk.fullF", f_fullF_constraint);
    wrt_r_item("lk.Ftune", f_Ftune);
    wrt_r_item("lk.U.rHB.RWq", f_rHB_RWq_cpue);
    wrt_r_item("lk.priors",f_priors);
    wrt_r_item("gradient.max",grad_max);

    wrt_r_item("w.L", w_L);

   // wrt_r_item("w.U.cGN", w_I_cGN);
	//wrt_r_item("w.U.cLL", w_I_cLL);
    wrt_r_item("w.U.rHB", w_cpue_rHB);
	
    wrt_r_item("w.lc.cGN", w_lenc_cGN);
   // wrt_r_item("w.lc.cLL", w_lc_cLL);	
	//wrt_r_item("w.lc.rGN", w_lc_rGN);
	wrt_r_item("w.ac.rGN", w_agec_rGN);
	
	wrt_r_item("w.Nage.init", w_Nage_init);
    wrt_r_item("w.R", w_rec);
    wrt_r_item("w.R.init", w_rec_early);
    wrt_r_item("w.R.end", w_rec_end);
    wrt_r_item("w.F.early.phases", w_fullF);
    wrt_r_item("w.Ftune.early.phases", w_Ftune);
	wrt_r_item("var.RWq", set_RWq_var);

 close_r_vector();

    open_r_matrix("N.age");
    wrt_r_matrix(N, 2, 2);
    wrt_r_namevector(styr, (endyr+1));
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
    wrt_r_namevector(styr, (endyr+1));
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

    open_r_matrix("len.cGN.mm");
    wrt_r_matrix(len_cGN_mm, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    //open_r_matrix("len.cLL.mm");
    //wrt_r_matrix(len_cLL_mm, 2, 2);
    //wrt_r_namevector(styr,endyr);
    //wrt_r_namevector(agebins);
    //close_r_matrix();
	
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

    open_r_matrix("wholewgt.cGN.lb");
    wrt_r_matrix(wholewgt_cGN_klb*1000.0, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

	//open_r_matrix("wholewgt.cLL.lb");
    //wrt_r_matrix(wholewgt_cLL_klb*1000.0, 2, 2);
    //wrt_r_namevector(styr,endyr);
    //wrt_r_namevector(agebins);
    //close_r_matrix();

    open_r_matrix("wholewgt.rHB.lb");
    wrt_r_matrix(wholewgt_rHB_klb*1000.0, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

    open_r_matrix("wholewgt.rGN.lb");
    wrt_r_matrix(wholewgt_rGN_klb*1000.0, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

 close_r_list();

open_r_list("sel.age");

    wrt_r_complete_vector("sel.v.wgted.L",sel_wgted_L, agebins);

    wrt_r_complete_vector("sel.v.wgted.tot",sel_wgted_tot, agebins);

    open_r_matrix("sel.m.cGN");
    wrt_r_matrix(sel_cGN, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);
    close_r_matrix();

	//open_r_matrix("sel.m.cLL");
    //wrt_r_matrix(sel_cLL, 2, 2);
    //wrt_r_namevector(styr,endyr);
    //wrt_r_namevector(agebins);
    //close_r_matrix();

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
	
 close_r_list();


//LIST object with predicted and observed composition data
open_r_list("comp.mats");

    open_r_matrix("lcomp.cGN.ob");
    wrt_r_matrix(obs_lenc_cGN, 2, 2);
    wrt_r_namevector(yrs_lenc_cGN);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    open_r_matrix("lcomp.cGN.pr");
    wrt_r_matrix(pred_cGN_lenc, 2, 2);
    wrt_r_namevector(yrs_lenc_cGN);
    wrt_r_namevector(lenbins);
    close_r_matrix();

    //open_r_matrix("lcomp.cLL.ob");
    //wrt_r_matrix(obs_cLL_lenc, 2, 2);
    //wrt_r_namevector(yrs_cLL_lenc);
    //wrt_r_namevector(lenbins);
    //close_r_matrix();
    //
    //open_r_matrix("lcomp.cLL.pr");
    //wrt_r_matrix(pred_cLL_lenc, 2, 2);
    //wrt_r_namevector(yrs_cLL_lenc);
    //wrt_r_namevector(lenbins);
    //close_r_matrix();

	//open_r_matrix("lcomp.rGN.ob");
    //wrt_r_matrix(obs_rGN_lenc, 2, 2);
    //wrt_r_namevector(yrs_rGN_lenc);
    //wrt_r_namevector(lenbins);
    //close_r_matrix();

    //open_r_matrix("lcomp.rGN.pr");
    //wrt_r_matrix(pred_rGN_lenc, 2, 2);
    //wrt_r_namevector(yrs_rGN_lenc);
    //wrt_r_namevector(lenbins);
    //close_r_matrix();

	open_r_matrix("acomp.rGN.ob");
    wrt_r_matrix(obs_agec_rGN, 2, 2);
    wrt_r_namevector(yrs_agec_rGN);
    wrt_r_namevector(agebins_agec);
    close_r_matrix();

    open_r_matrix("acomp.rGN.pr");
    wrt_r_matrix(pred_rGN_agec, 2, 2);
    wrt_r_namevector(yrs_agec_rGN);
    wrt_r_namevector(agebins_agec);
    close_r_matrix();
	
 close_r_list();

// DATA FRAME of time series
open_r_df("t.series", styr, (endyr+1), 2);
    wrt_r_namevector(styr,(endyr+1));
    wrt_r_df_col("year", styr,(endyr+1));
    wrt_r_df_col("F.Fmsy", FdF_msy);
	wrt_r_df_col("F.F30.ratio", FdF30);	//*.ratio extension is for FishGraph
	wrt_r_df_col("F.F40.ratio", FdF40);	//*.ratio extension is for FishGraph
    wrt_r_df_col("F.full", Fapex);
    wrt_r_df_col("F.cGN", F_cGN_out);
    //wrt_r_df_col("F.cLL", F_cLL_out);	;
    wrt_r_df_col("F.rGN", F_rGN_out);

    wrt_r_df_col("Fsum", Fsum);
    wrt_r_df_col("N", totN); //abundance at start of year
    wrt_r_df_col("recruits", rec);
    wrt_r_df_col("logR.dev", log_rec_dev_output); //places zeros in yrs deviations not estimated  //KWS
    wrt_r_df_col("SSB", SSB);
	wrt_r_df_col("SSB_knum", SSB_knum);
    wrt_r_df_col("SSB.SSBmsy", SdSSB_msy);
    //wrt_r_df_col("SSB.msst.M", SdSSB_msy/smsy2msstM);
    wrt_r_df_col("SSB.msst", SdSSB_msy/smsy2msst75);
	wrt_r_df_col("SSB.SSBF30", SdSSB_F30);
    wrt_r_df_col("SSB.msstF30", Sdmsst_F30);
	wrt_r_df_col("SSB.SSBF40", SdSSB_F40);
    wrt_r_df_col("SSB.msstF40", Sdmsst_F40);
    wrt_r_df_col("B", totB);
    wrt_r_df_col("B.B0", totB/B0);
    wrt_r_df_col("SPR.static", spr_static);

    wrt_r_df_col("total.L.klb", L_total_klb_yr);
    wrt_r_df_col("total.L.knum", L_total_knum_yr);

    //wrt_r_df_col("U.cGN.ob", obs_cGN_cpue);
    //wrt_r_df_col("U.cGN.pr", pred_cGN_cpue);
    //wrt_r_df_col("cv.U.cGN", cGN_cpue_cv/w_I_cGN);  //applied CV after weighting
    //wrt_r_df_col("cv.unwgted.U.cGN", cGN_cpue_cv); //CV before weighting
	
   // wrt_r_df_col("U.cLL.ob", obs_cLL_cpue);
   // wrt_r_df_col("U.cLL.pr", pred_cLL_cpue);
   // wrt_r_df_col("cv.U.cLL", cLL_cpue_cv/w_I_cLL);  //applied CV after weighting
   // wrt_r_df_col("cv.unwgted.U.cLL", cLL_cpue_cv); //CV before weighting	

    wrt_r_df_col("U.rHB.ob", obs_cpue_rHB);
    wrt_r_df_col("U.rHB.pr", pred_rHB_cpue);
    wrt_r_df_col("cv.U.rHB", obs_cv_cpue_rHB/w_cpue_rHB);
	wrt_r_df_col("cv.unwgted.U.rHB", obs_cv_cpue_rHB);
	
	//wrt_r_df_col("q.cGN", q_cGN);
    //wrt_r_df_col("q.cGN.rate.mult",q_rate_fcn_cGN);
    //wrt_r_df_col("q.cGN.RW.log.dev",q_RW_log_dev_cGN);
	
	//wrt_r_df_col("q.cLL", q_cLL);
    //wrt_r_df_col("q.cLL.rate.mult",q_rate_fcn_cLL);
    //wrt_r_df_col("q.cLL.RW.log.dev",q_RW_log_dev_cLL);	

    wrt_r_df_col("q.rHB", q_rHB);
    wrt_r_df_col("q.rHB.rate.mult",q_rate_fcn_rHB);
    wrt_r_df_col("q.rHB.RW.log.dev",q_RW_log_dev_rHB);

	
	wrt_r_df_col("q.DD.mult", q_DD_fcn);
    wrt_r_df_col("q.DD.B.exploitable", B_q_DD);

    wrt_r_df_col("L.cGN.ob", obs_L_cGN);
    wrt_r_df_col("L.cGN.pr", pred_cGN_L_klb);
    wrt_r_df_col("cv.L.cGN", obs_cv_L_cGN);

	//wrt_r_df_col("L.cLL.ob", obs_cLL_L);
    //wrt_r_df_col("L.cLL.pr", pred_cLL_L_klb);
    //wrt_r_df_col("cv.L.cLL", cLL_L_cv);

    wrt_r_df_col("L.rGN.ob", obs_L_rGN);
    wrt_r_df_col("L.rGN.pr", pred_rGN_L_knum);
    wrt_r_df_col("cv.L.rGN", obs_cv_L_rGN);


    //comp sample sizes
	wrt_r_df_col("lcomp.cGN.n", nsamp_cGN_lenc_allyr);
	//wrt_r_df_col("lcomp.cLL.n", nsamp_cLL_lenc_allyr);
	//wrt_r_df_col("lcomp.rGN.n", nsamp_rGN_lenc_allyr);

	wrt_r_df_col("acomp.rGN.n", nsamp_rGN_agec_allyr);
	
	wrt_r_df_col("lcomp.cGN.nfish", nfish_cGN_lenc_allyr);
	//wrt_r_df_col("lcomp.cLL.nfish", nfish_cLL_lenc_allyr);
	//wrt_r_df_col("lcomp.rGN.nfish", nfish_rGN_lenc_allyr);
	
	wrt_r_df_col("acomp.rGN.nfish", nfish_rGN_agec_allyr);
	
	wrt_r_df_col("lcomp.cGN.neff", (1+nsamp_cGN_lenc_allyr*exp(log_dm_cGN_lc_out(8)))/(1+exp(log_dm_cGN_lc_out(8))) );
	//wrt_r_df_col("lcomp.cLL.neff", (1+nsamp_cLL_lenc_allyr*exp(log_dm_cLL_lc_out(8)))/(1+exp(log_dm_cLL_lc_out(8))) );
	//wrt_r_df_col("lcomp.rGN.neff", (1+nsamp_rGN_lenc_allyr*exp(log_dm_rGN_lc_out(8)))/(1+exp(log_dm_rGN_lc_out(8))) );   
	
	wrt_r_df_col("acomp.rGN.neff", (1+nsamp_rGN_agec_allyr*exp(log_dm_rGN_ac_out(8)))/(1+exp(log_dm_rGN_ac_out(8))) );
	
 close_r_df();

// DATA FRAME of L and D time series by fishery
open_r_df("LD.pr.tseries", styr, endyr, 2);
	wrt_r_namevector(styr,endyr);
	wrt_r_df_col("year", styr, endyr);

    wrt_r_df_col("L.cGN.klb", pred_cGN_L_klb);
    wrt_r_df_col("L.cGN.knum", pred_cGN_L_knum);
    //wrt_r_df_col("L.cLL.klb", pred_cLL_L_klb);
    //wrt_r_df_col("L.cLL.knum", pred_cLL_L_knum);	
    wrt_r_df_col("L.rGN.klb", pred_rGN_L_klb);
    wrt_r_df_col("L.rGN.knum", pred_rGN_L_knum);

close_r_df();

 open_r_df("a.series", 1, nages, 2);
 	wrt_r_namevector(1,nages);
 	wrt_r_df_col("age", agebins);
 	wrt_r_df_col("length", meanlen_TL);
 	wrt_r_df_col("length.cv", len_cv);
 	wrt_r_df_col("length.sd", len_sd);
 	wrt_r_df_col("weight", wgt_lb);     //for FishGraph
 	wrt_r_df_col("wgt.klb", wgt_klb);
 	wrt_r_df_col("wgt.mt", wgt_mt);
 	wrt_r_df_col("wholewgt.wgted.L.klb", wgt_wgted_L_klb);
    wrt_r_df_col("prop.female", prop_f);	
   	wrt_r_df_col("mat.female", maturity_f);
 	wrt_r_df_col("reprod", reprod);
	wrt_r_df_col("reprodknum", reprodknum);
 	
 	wrt_r_df_col("M", M);
 	wrt_r_df_col("F.initial", F_initial);
 	wrt_r_df_col("Z.initial", Z_initial);
 	wrt_r_df_col("Nage.eq.init",N_initial_eq);
 	wrt_r_df_col("log.Nage.init.dev",log_Nage_dev_output);
 close_r_df();


 open_r_df("eq.series", 1, n_iter_msy, 2);
 	wrt_r_namevector(1,n_iter_msy);
 	wrt_r_df_col("F.eq", F_msy);
 	wrt_r_df_col("spr.eq", spr_msy);
 	wrt_r_df_col("R.eq", R_eq);
 	wrt_r_df_col("SSB.eq", SSB_eq);
 	wrt_r_df_col("B.eq", B_eq);
 	wrt_r_df_col("L.eq.wholeklb", L_eq_klb);
 	wrt_r_df_col("L.eq.knum", L_eq_knum);
 close_r_df();

 open_r_df("pr.series", 1, n_iter_spr, 2);
 	wrt_r_namevector(1,n_iter_spr);
 	wrt_r_df_col("F.spr", F_spr);
 	wrt_r_df_col("spr", spr_spr);
 	wrt_r_df_col("SPR", spr_ratio);
 	wrt_r_df_col("ypr", L_spr); //whole weight
 close_r_df();


 open_r_list("CLD.est.mats");

     open_r_matrix("Lw.cGN");
         wrt_r_matrix(L_cGN_klb, 1,1);
     close_r_matrix();

	// open_r_matrix("Lw.cLL");
    //     wrt_r_matrix(L_cLL_klb, 1,1);
    // close_r_matrix();

     open_r_matrix("Lw.rGN");
         wrt_r_matrix(L_rGN_klb, 1,1);
     close_r_matrix();

     open_r_matrix("Lw.total");
         wrt_r_matrix(L_total_klb, 1,1);
     close_r_matrix();


     open_r_matrix("Ln.cGN");
         wrt_r_matrix(L_cGN_num, 1,1);
     close_r_matrix();

	// open_r_matrix("Ln.cLL");
    //     wrt_r_matrix(L_cLL_num, 1,1);
    // close_r_matrix();
	 
     open_r_matrix("Ln.rGN");
         wrt_r_matrix(L_rGN_num, 1,1);
     close_r_matrix();

     open_r_matrix("Ln.total");
         wrt_r_matrix(L_total_num, 1,1);
     close_r_matrix();

  close_r_list();

  //LIST object of age error matrix
open_r_list("age.error");

    open_r_matrix("error.mat");		
    wrt_r_matrix(age_error, 2, 2);
    wrt_r_namevector(agebins);
    wrt_r_namevector(agebins);
    close_r_matrix();
  
close_r_list();

// DATA FRAME of projection time series
open_r_df("projection", styr_proj, endyr_proj, 2);
    wrt_r_namevector(styr_proj,endyr_proj);
    wrt_r_df_col("year", styr_proj,endyr_proj);
    wrt_r_df_col("F.proj", F_proj);
    wrt_r_df_col("SSB.proj", SSB_proj);
    wrt_r_df_col("B.proj", B_proj);
    wrt_r_df_col("R.proj", R_proj);
	wrt_r_df_col("L.knum.proj", L_knum_proj);
	wrt_r_df_col("L.klb.proj", L_klb_proj);	
 close_r_df();

close_r_file();
