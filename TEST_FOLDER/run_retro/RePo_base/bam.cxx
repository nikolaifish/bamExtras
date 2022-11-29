//##  Author: NMFS, Beaufort Lab, Sustainable Fisheries Branch
//##  Analyst: Nikolai Klibansky
//##  Species: Red Porgy
//##  Region: US South Atlantic
//##  SEDAR: 60
//##  Date: 2022-11-18 18:05:32


//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//###################################################################
//##  BAM R Object Creator File #####################################
//###################################################################
//##  Author: NMFS, Beaufort Lab, Sustainable Fisheries Branch
//##  Analyst: Nikolai Klibansky
//##  Species: Red Porgy
//##  Region: US South Atlantic
//##  SEDAR: 60
//##  Date: 2020-05-29 10:55:05
//###################################################################
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

// Create a file with an R object from AD Model Builder

// open the file using the default AD Model Builder file name, and
// 6 digits of precision
open_r_file(adprogram_name + ".rdat", 6);

// Example of an INFO object
open_r_info_list("info", true);
	wrt_r_item("title", "SEDAR 60");
	wrt_r_item("species", "Red Porgy");
	wrt_r_item("model", "Statistical Catch at Age");
	if(SR_switch==1)
	    {wrt_r_item("rec.model", "BH-steep");}
    if(SR_switch==2)
	    {wrt_r_item("rec.model", "Ricker-steep");}
	if(SR_switch==3)
	    {wrt_r_item("rec.model", "Mean");}	
	wrt_r_item("base.run", "RGX.tpl");
	wrt_r_item("units.length", "mm");
	wrt_r_item("units.weight", "lb");
	wrt_r_item("units.biomass", "metric tons");
	wrt_r_item("units.ssb", "metric tons");
	wrt_r_item("units.ypr", "lb whole");
	wrt_r_item("units.landings", "1000 lb whole");
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
	wrt_r_item("M.constant", M_constant);
	wrt_r_item("smsy2msstM", smsy2msstM);
    wrt_r_item("smsy2msst75", smsy2msst75);
	wrt_r_item("Linf", Linf);
	wrt_r_item("K", K);
	wrt_r_item("t0", t0);
	wrt_r_item("len_cv_val",len_cv_val);
	wrt_r_item("wgt.a", wgtpar_a);
	wrt_r_item("wgt.b", wgtpar_b);
    wrt_r_item("spawn.time", spawn_time_frac);
	wrt_r_item("D.mort.cHL", Dmort_cHL);
	wrt_r_item("D.mort.rHB", Dmort_rHB);
	wrt_r_item("D.mort.rGN", Dmort_rGN);
	wrt_r_item("q.rHB", mfexp(log_q_cpue_rHB));
	wrt_r_item("q.sCT", mfexp(log_q_cpue_sCT));

    wrt_r_item("q.rate",q_rate);
	wrt_r_item("q.beta",q_DD_beta);
	wrt_r_item("q.DD.B0.exploitable", B0_q_DD);
	wrt_r_item("F.prop.L.cHL", F_prop_L_cHL);
	wrt_r_item("F.prop.L.cTW", F_prop_L_cTW);	
	wrt_r_item("F.prop.L.rHB", F_prop_L_rHB);
	wrt_r_item("F.prop.L.rGN", F_prop_L_rGN);
	wrt_r_item("F.prop.D.cHL", F_prop_D_cHL);
	wrt_r_item("F.prop.D.rHB", F_prop_D_rHB);
	wrt_r_item("F.prop.D.rGN", F_prop_D_rGN);
	wrt_r_item("Fend.mean", Fend_mean);

	wrt_r_item("B0", B0);
	wrt_r_item("Bstyr.B0", totB(styr)/B0);
	wrt_r_item("SSB0", S0);
	wrt_r_item("SSBstyr.SSB0", SSB(styr)/S0);
	wrt_r_item("Rstyr.R0", rec(styr)/R0);
	if(SR_switch==1)
	{
      wrt_r_item("BH.biascorr",BiasCor);
	  wrt_r_item("BH.Phi0", spr_F0(endyr));
	  wrt_r_item("BH.R0", R0);
	  wrt_r_item("BH.steep", steep);
    }
    if(SR_switch==2)
	{
      wrt_r_item("Ricker.biascorr",BiasCor);
	  wrt_r_item("Ricker.Phi0", spr_F0(endyr));
	  wrt_r_item("Ricker.R0", R0);
	  wrt_r_item("Ricker.steep", steep);
    }
    if(SR_switch==3)
	{
      wrt_r_item("Mean.biascorr",BiasCor);
	  wrt_r_item("Mean.R0", R0);
    }		
	wrt_r_item("R.sigma.logdevs", sigma_rec_dev);
	wrt_r_item("R.sigma.par",rec_sigma);
	wrt_r_item("R.autocorr",R_autocorr);
	wrt_r_item("R0", R0); //same as BH.R0, but used in BSR.time.plots
	wrt_r_item("R.virgin.bc", R_virgin); //bias-corrected virgin recruitment
	wrt_r_item("rec.lag", 1.0);
	wrt_r_item("msy.klb", msy_klb_out);
	wrt_r_item("msy.knum", msy_knum_out);
	wrt_r_item("Dmsy.klb", D_msy_klb_out);
	wrt_r_item("Dmsy.knum", D_msy_knum_out);	
	wrt_r_item("Fmsy", F_msy_out);
	wrt_r_item("SSBmsy", SSB_msy_out);
	wrt_r_item("msst", smsy2msstM*SSB_msy_out);
	wrt_r_item("Bmsy", B_msy_out);
	wrt_r_item("Rmsy", R_msy_out);
	wrt_r_item("sprmsy",spr_msy_out);
	wrt_r_item("Fend.Fmsy", FdF_msy_end);
	wrt_r_item("Fend.Fmsy.mean", FdF_msy_end_mean);
	wrt_r_item("SSBend.SSBmsy", SdSSB_msy_end);
	wrt_r_item("SSBend.MSST", SdSSB_msy_end/smsy2msstM);
	wrt_r_item("F30", F30_out); //For FishGraph
	wrt_r_item("SSB.F30", SSB_F30_out); //For FishGraph
	wrt_r_item("msst.F30", smsy2msstM*SSB_F30_out);
	 wrt_r_item("B.F30", B_F30_out);
	wrt_r_item("R.F30", R_F30_out);
    wrt_r_item("L.F30.knum", L_F30_knum_out);
    wrt_r_item("L.F30.klb", L_F30_klb_out);
    wrt_r_item("D.F30.knum", D_F30_knum_out);
    wrt_r_item("D.F30.klb", D_F30_klb_out);
 close_r_info_list();
	
 // VECTOR object of parameters and estimated quantities
open_r_info_list("spr.brps", false);
	wrt_r_item("F20", F20_out);
	wrt_r_item("F30", F30_out);
	wrt_r_item("F40", F40_out);
	wrt_r_item("SSB.F30", SSB_F30_out);
    wrt_r_item("msst.F30", smsy2msstM*SSB_F30_out);   
    wrt_r_item("B.F30", B_F30_out);
	wrt_r_item("R.F30", R_F30_out);
    wrt_r_item("L.F30.knum", L_F30_knum_out);
    wrt_r_item("L.F30.klb", L_F30_klb_out);
    wrt_r_item("D.F30.knum", D_F30_knum_out);
    wrt_r_item("D.F30.klb", D_F30_klb_out);        	 	
	wrt_r_item("Fend.F30.mean", FdF30_end_mean);
    wrt_r_item("SSBend.SSBF30", SdSSB_F30_end);
	wrt_r_item("SSBend.MSSTF30", Sdmsst_F30_end);	
 close_r_info_list();

// MATRIX object of parameter constraints
open_r_df("parm.cons",1,8,2);
    wrt_r_namevector(1,8);
    wrt_r_df_col("Linf",Linf_out);
    wrt_r_df_col("K",K_out);
    wrt_r_df_col("t0",t0_out);
    wrt_r_df_col("len_cv_val",len_cv_val_out);
    wrt_r_df_col("log_R0",log_R0_out);
    wrt_r_df_col("steep",steep_out);
    wrt_r_df_col("rec_sigma",rec_sigma_out);
    wrt_r_df_col("R_autocorr",R_autocorr_out);
	wrt_r_df_col("log_dm_lenc_cTW",log_dm_lenc_cTW_out);
	wrt_r_df_col("log_dm_agec_cHL",log_dm_agec_cHL_out);
	wrt_r_df_col("log_dm_agec_rHB",log_dm_agec_rHB_out);
	wrt_r_df_col("log_dm_agec_sCT",log_dm_agec_sCT_out);
	
    wrt_r_df_col("A50_sel_cHL2",A50_sel_cHL2_out);
    wrt_r_df_col("slope_sel_cHL2",slope_sel_cHL2_out);
    wrt_r_df_col("A50_sel_cHL3",A50_sel_cHL3_out);
    wrt_r_df_col("slope_sel_cHL3",slope_sel_cHL3_out);
	

    wrt_r_df_col("A50_sel_cTW",A50_sel_cTW_out);
    wrt_r_df_col("slope_sel_cTW",slope_sel_cTW_out);
	
    wrt_r_df_col("A50_sel_rHB1",A50_sel_rHB1_out);
    wrt_r_df_col("slope_sel_rHB1",slope_sel_rHB1_out);
 	wrt_r_df_col("A50_sel_rHB2",A50_sel_rHB2_out);
    wrt_r_df_col("slope_sel_rHB2",slope_sel_rHB2_out);
	wrt_r_df_col("A50_sel_rHB3",A50_sel_rHB3_out);
    wrt_r_df_col("slope_sel_rHB3",slope_sel_rHB3_out);

		
    wrt_r_df_col("A50_sel_sCT",A50_sel_sCT_out);
    wrt_r_df_col("slope_sel_sCT",slope_sel_sCT_out);

	
    
    wrt_r_df_col("log_q_cpue_rHB",log_q_cpue_rHB_out);
	wrt_r_df_col("log_q_cpue_sCT",log_q_cpue_sCT_out);
    wrt_r_df_col("M_constant",M_constant_out);
    wrt_r_df_col("log_avg_F_L_cHL",log_avg_F_L_cHL_out);
	wrt_r_df_col("log_avg_F_L_cTW",log_avg_F_L_cTW_out);
    wrt_r_df_col("log_avg_F_L_rHB",log_avg_F_L_rHB_out);
    wrt_r_df_col("log_avg_F_L_rGN",log_avg_F_L_rGN_out);
    wrt_r_df_col("log_avg_F_D_cHL",log_avg_F_D_cHL_out);
    wrt_r_df_col("log_avg_F_D_rHB",log_avg_F_D_rHB_out);
    wrt_r_df_col("log_avg_F_D_rGN",log_avg_F_D_rGN_out);
    //wrt_r_df_col("F_init",F_init_out);
close_r_df();

// DATA FRAME of time series deviation vector estimates
// names used in this object must match the names used in the "parm.tvec.cons" object
open_r_df("parm.tvec", styr, endyr, 2);
	wrt_r_namevector(styr,endyr);
	wrt_r_df_col("year", styr,endyr);
    wrt_r_df_col("log.rec.dev", log_rec_dev_out);
    wrt_r_df_col("log.F.dev.L.cHL", log_F_dev_L_cHL_out);
	wrt_r_df_col("log.F.dev.L.cTW", log_F_dev_L_cTW_out);
    wrt_r_df_col("log.F.dev.L.rHB", log_F_dev_L_rHB_out);
    wrt_r_df_col("log.F.dev.L.rGN", log_F_dev_L_rGN_out);
    wrt_r_df_col("log.F.dev.D.cHL", log_F_dev_D_cHL_out);
    wrt_r_df_col("log.F.dev.D.rHB", log_F_dev_D_rHB_out);
    wrt_r_df_col("log.F.dev.D.rGN", log_F_dev_D_rGN_out);
	wrt_r_df_col("log.q.dev.rHB.RWq", q_RW_log_dev_cpue_rHB);
	
close_r_df();

// MATRIX object of deviation vector constraints
// names used in this object must match the names used in the "parm.tvec" object
open_r_df("parm.tvec.cons",1,3,2);
    wrt_r_namevector(1,3);
    wrt_r_df_col("log.rec.dev",set_log_dev_rec);
    wrt_r_df_col("log.F.dev.L.cHL",set_log_dev_F_L_cHL);
	wrt_r_df_col("log.F.dev.L.cTW",set_log_dev_F_L_cTW);
    wrt_r_df_col("log.F.dev.L.rHB",set_log_dev_F_L_rHB);
    wrt_r_df_col("log.F.dev.L.rGN",set_log_dev_F_L_rGN);
    wrt_r_df_col("log.F.dev.D.cHL",set_log_dev_F_D_cHL);
    wrt_r_df_col("log.F.dev.D.rHB",set_log_dev_F_D_rHB);
    wrt_r_df_col("log.F.dev.D.rGN",set_log_dev_F_D_rGN);
	wrt_r_df_col("log.q.dev.rHB.RWq", set_log_dev_RWq);
	wrt_r_df_col("log.q.dev.cHL.RWq", set_log_dev_RWq);
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
    wrt_r_item("sdnr.U.rHB", sdnr_cpue_rHB);
	wrt_r_item("sdnr.U.sCT", sdnr_cpue_sCT);
	
//     wrt_r_item("sdnr.lc.cTW", sdnr_lenc_cTW);	
	
    wrt_r_item("sdnr.ac.cHL", sdnr_agec_cHL);
    wrt_r_item("sdnr.ac.rHB", sdnr_agec_rHB);
	wrt_r_item("sdnr.ac.sCT", sdnr_agec_sCT);	
	
 close_r_vector();

// VECTOR object of likelihood contributions
open_r_vector("like");
    wrt_r_item("lk.total", fval); //weighted likelihood
    wrt_r_item("lk.unwgt.data", fval_data);  //likelihood of just data components

    wrt_r_item("lk.L.cHL", f_L_cHL);
    wrt_r_item("lk.L.cTW", f_L_cTW);	
    wrt_r_item("lk.L.rHB", f_L_rHB);
    wrt_r_item("lk.L.rGN", f_L_rGN);
    wrt_r_item("lk.D.cHL", f_D_cHL);
    wrt_r_item("lk.D.rHB", f_D_rHB);
    wrt_r_item("lk.D.rGN", f_D_rGN);

    wrt_r_item("lk.U.rHB", f_cpue_rHB);
	wrt_r_item("lk.U.sCT", f_cpue_sCT);
    wrt_r_item("lk.lenc.cTW", f_lenc_cTW);	
    wrt_r_item("lk.agec.cHL", f_agec_cHL);
    wrt_r_item("lk.agec.rHB", f_agec_rHB);
	wrt_r_item("lk.agec.sCT", f_agec_sCT);	
    wrt_r_item("lk.Nage.init", f_Nage_init);
    wrt_r_item("lk.SRfit", f_rec_dev);
    wrt_r_item("lk.SRearly", f_rec_dev_early);
    wrt_r_item("lk.SRend", f_rec_dev_end);
    wrt_r_item("lk.fullF", f_fullF_constraint);
    wrt_r_item("lk.Ftune", f_Ftune);
	wrt_r_item("lk.U.rHB.RWq", f_RWq_cpue_rHB);
    wrt_r_item("lk.priors",f_priors);
    wrt_r_item("gradient.max",grad_max);

    wrt_r_item("w.L", w_L);
    wrt_r_item("w.D", w_D);
    wrt_r_item("w.U.rHB", w_cpue_rHB);
	wrt_r_item("w.U.sCT", w_cpue_sCT);
	
    wrt_r_item("w.lc.cTW", w_lenc_cTW);	
	
    wrt_r_item("w.ac.cHL", w_agec_cHL);
    wrt_r_item("w.ac.rHB", w_agec_rHB);
	wrt_r_item("w.ac.sCT", w_agec_sCT);

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
	
	open_r_matrix("D.age.pred.knum");
    wrt_r_matrix(D_total_num/1000.0, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);   
    close_r_matrix();
    
    open_r_matrix("D.age.pred.klb");
    wrt_r_matrix(D_total_klb, 2, 2);
    wrt_r_namevector(styr,endyr);
    wrt_r_namevector(agebins);   
    close_r_matrix();
 
// LIST object with annual selectivity at age by fishery

open_r_list("size.age.fishery");

    open_r_matrix("len.cHL.mm"); //LINETAG: len.cHL.mm
    wrt_r_matrix(len_cHL_mm, 2, 2); //LINETAG: len.cHL.mm
    wrt_r_namevector(styr,endyr); //LINETAG: len.cHL.mm
    wrt_r_namevector(agebins); //LINETAG: len.cHL.mm
    close_r_matrix(); //LINETAG: len.cHL.mm


    open_r_matrix("len.cTW.mm"); //LINETAG: len.cTW.mm
    wrt_r_matrix(len_cTW_mm, 2, 2); //LINETAG: len.cTW.mm
    wrt_r_namevector(styr,endyr); //LINETAG: len.cTW.mm
    wrt_r_namevector(agebins); //LINETAG: len.cTW.mm
    close_r_matrix(); //LINETAG: len.cTW.mm
	
    open_r_matrix("len.rHB.mm"); //LINETAG: len.rHB.mm
    wrt_r_matrix(len_rHB_mm, 2, 2); //LINETAG: len.rHB.mm
    wrt_r_namevector(styr,endyr); //LINETAG: len.rHB.mm
    wrt_r_namevector(agebins); //LINETAG: len.rHB.mm
    close_r_matrix(); //LINETAG: len.rHB.mm

    open_r_matrix("len.rGN.mm"); //LINETAG: len.rGN.mm
    wrt_r_matrix(len_rGN_mm, 2, 2); //LINETAG: len.rGN.mm
    wrt_r_namevector(styr,endyr); //LINETAG: len.rGN.mm
    wrt_r_namevector(agebins); //LINETAG: len.rGN.mm
    close_r_matrix(); //LINETAG: len.rGN.mm

    open_r_matrix("len.cHL.D.mm"); //LINETAG: len.cHL.D.mm
    wrt_r_matrix(len_D_cHL_mm, 2, 2);  //LINETAG: len.cHL.D.mm
    wrt_r_namevector(styr,endyr);  //LINETAG: len.cHL.D.mm
    wrt_r_namevector(agebins);  //LINETAG: len.cHL.D.mm
    close_r_matrix();  //LINETAG: len.cHL.D.mm

    open_r_matrix("len.rHB.D.mm");  //LINETAG: len.rHB.D.mm
    wrt_r_matrix(len_D_rHB_mm, 2, 2);  //LINETAG: len.rHB.D.mm
    wrt_r_namevector(styr,endyr);  //LINETAG: len.rHB.D.mm
    wrt_r_namevector(agebins);  //LINETAG: len.rHB.D.mm
    close_r_matrix();  //LINETAG: len.rHB.D.mm

    open_r_matrix("len.rGN.D.mm");  //LINETAG: len.rGN.D.mm
    wrt_r_matrix(len_D_rGN_mm, 2, 2);  //LINETAG: len.rGN.D.mm
    wrt_r_namevector(styr,endyr);  //LINETAG: len.rGN.D.mm
    wrt_r_namevector(agebins);  //LINETAG: len.rGN.D.mm
    close_r_matrix();  //LINETAG: len.rGN.D.mm



    open_r_matrix("wholewgt.cHL.lb");  //LINETAG: wholewgt.cHL.lb
    wrt_r_matrix(wholewgt_cHL_klb*1000.0, 2, 2);  //LINETAG: wholewgt.cHL.lb
    wrt_r_namevector(styr,endyr);  //LINETAG: wholewgt.cHL.lb
    wrt_r_namevector(agebins);  //LINETAG: wholewgt.cHL.lb
    close_r_matrix();  //LINETAG: wholewgt.cHL.lb


	open_r_matrix("wholewgt.cTW.lb");  //LINETAG: wholewgt.cTW.lb
    wrt_r_matrix(wholewgt_cTW_klb*1000.0, 2, 2);  //LINETAG: wholewgt.cTW.lb
    wrt_r_namevector(styr,endyr);  //LINETAG: wholewgt.cTW.lb
    wrt_r_namevector(agebins);  //LINETAG: wholewgt.cTW.lb
    close_r_matrix();  //LINETAG: wholewgt.cTW.lb

    open_r_matrix("wholewgt.rHB.lb");  //LINETAG: wholewgt.rHB.lb
    wrt_r_matrix(wholewgt_rHB_klb*1000.0, 2, 2);  //LINETAG: wholewgt.rHB.lb
    wrt_r_namevector(styr,endyr);  //LINETAG: wholewgt.rHB.lb
    wrt_r_namevector(agebins);  //LINETAG: wholewgt.rHB.lb
    close_r_matrix();  //LINETAG: wholewgt.rHB.lb

    open_r_matrix("wholewgt.rGN.lb");  //LINETAG: wholewgt.rGN.lb
    wrt_r_matrix(wholewgt_rGN_klb*1000.0, 2, 2);  //LINETAG: wholewgt.rGN.lb
    wrt_r_namevector(styr,endyr);  //LINETAG: wholewgt.rGN.lb
    wrt_r_namevector(agebins);  //LINETAG: wholewgt.rGN.lb
    close_r_matrix();  //LINETAG: wholewgt.rGN.lb

    open_r_matrix("wholewgt.cHL.D.lb");  //LINETAG: wholewgt.cHL.D.lb
    wrt_r_matrix(wholewgt_D_cHL_klb*1000.0, 2, 2);  //LINETAG: wholewgt.cHL.D.lb
    wrt_r_namevector(styr,endyr);  //LINETAG: wholewgt.cHL.D.lb
    wrt_r_namevector(agebins);  //LINETAG: wholewgt.cHL.D.lb
    close_r_matrix();  //LINETAG: wholewgt.cHL.D.lb

    open_r_matrix("wholewgt.rHB.D.lb");  //LINETAG: wholewgt.rHB.D.lb
    wrt_r_matrix(wholewgt_D_rHB_klb*1000.0, 2, 2);  //LINETAG: wholewgt.rHB.D.lb
    wrt_r_namevector(styr,endyr);  //LINETAG: wholewgt.rHB.D.lb
    wrt_r_namevector(agebins);  //LINETAG: wholewgt.rHB.D.lb
    close_r_matrix();  //LINETAG: wholewgt.rHB.D.lb

    open_r_matrix("wholewgt.rGN.D.lb");  //LINETAG: wholewgt.rGN.D.lb
    wrt_r_matrix(wholewgt_D_rGN_klb*1000.0, 2, 2);  //LINETAG: wholewgt.rGN.D.lb
    wrt_r_namevector(styr,endyr);  //LINETAG: wholewgt.rGN.D.lb
    wrt_r_namevector(agebins);  //LINETAG: wholewgt.rGN.D.lb
    close_r_matrix();  //LINETAG: wholewgt.rGN.D.lb

 close_r_list();

open_r_list("sel.age");

    wrt_r_complete_vector("sel.v.wgted.L",sel_wgted_L, agebins);
    wrt_r_complete_vector("sel.v.wgted.D",sel_wgted_D, agebins);
    wrt_r_complete_vector("sel.v.wgted.tot",sel_wgted_tot, agebins);

    open_r_matrix("sel.m.cHL");  //LINETAG: sel.m.cHL
    wrt_r_matrix(sel_cHL, 2, 2);  //LINETAG: sel.m.cHL
    wrt_r_namevector(styr,endyr);  //LINETAG: sel.m.cHL
    wrt_r_namevector(agebins);  //LINETAG: sel.m.cHL
    close_r_matrix();  //LINETAG: sel.m.cHL


	open_r_matrix("sel.m.cTW");  //LINETAG: sel.m.cTW
    wrt_r_matrix(sel_cTW, 2, 2);  //LINETAG: sel.m.cTW
    wrt_r_namevector(styr,endyr);  //LINETAG: sel.m.cTW
    wrt_r_namevector(agebins);  //LINETAG: sel.m.cTW
    close_r_matrix();  //LINETAG: sel.m.cTW

    open_r_matrix("sel.m.rHB");  //LINETAG: sel.m.rHB
    wrt_r_matrix(sel_rHB, 2, 2);  //LINETAG: sel.m.rHB
    wrt_r_namevector(styr,endyr);  //LINETAG: sel.m.rHB
    wrt_r_namevector(agebins);  //LINETAG: sel.m.rHB
    close_r_matrix();  //LINETAG: sel.m.rHB

    open_r_matrix("sel.m.rGN");  //LINETAG: sel.m.rGN
    wrt_r_matrix(sel_rGN, 2, 2);  //LINETAG: sel.m.rGN
    wrt_r_namevector(styr,endyr);  //LINETAG: sel.m.rGN
    wrt_r_namevector(agebins);  //LINETAG: sel.m.rGN
    close_r_matrix();  //LINETAG: sel.m.rGN

    open_r_matrix("sel.m.cHL.D");  //LINETAG: sel.m.cHL.D
    wrt_r_matrix(sel_D_cHL, 2, 2);  //LINETAG: sel.m.cHL.D
    wrt_r_namevector(styr,endyr);  //LINETAG: sel.m.cHL.D
    wrt_r_namevector(agebins);  //LINETAG: sel.m.cHL.D
    close_r_matrix();  //LINETAG: sel.m.cHL.D
    open_r_matrix("sel.m.rHB.D");  //LINETAG: sel.m.rHB.D
    wrt_r_matrix(sel_D_rHB, 2, 2);  //LINETAG: sel.m.rHB.D
    wrt_r_namevector(styr,endyr);  //LINETAG: sel.m.rHB.D
    wrt_r_namevector(agebins);  //LINETAG: sel.m.rHB.D
    close_r_matrix();  //LINETAG: sel.m.rHB.D
    open_r_matrix("sel.m.rGN.D");  //LINETAG: sel.m.rGN.D
    wrt_r_matrix(sel_D_rGN, 2, 2);  //LINETAG: sel.m.rGN.D
    wrt_r_namevector(styr,endyr);  //LINETAG: sel.m.rGN.D
    wrt_r_namevector(agebins);  //LINETAG: sel.m.rGN.D
    close_r_matrix();  //LINETAG: sel.m.rGN.D
	
	open_r_matrix("sel.m.sCT");  //LINETAG: sel.m.sCT
    wrt_r_matrix(sel_sCT, 2, 2);  //LINETAG: sel.m.sCT
    wrt_r_namevector(styr,endyr);  //LINETAG: sel.m.sCT
    wrt_r_namevector(agebins);  //LINETAG: sel.m.sCT
    close_r_matrix();  //LINETAG: sel.m.sCT

	
 close_r_list();


//LIST object with predicted and observed composition data
open_r_list("comp.mats");




    open_r_matrix("lcomp.cTW.ob");  //LINETAG: lcomp.cTW.ob
    wrt_r_matrix(obs_lenc_cTW, 2, 2);  //LINETAG: lcomp.cTW.ob
    wrt_r_namevector(yrs_lenc_cTW);  //LINETAG: lcomp.cTW.ob
    wrt_r_namevector(lenbins);  //LINETAG: lcomp.cTW.ob
    close_r_matrix();  //LINETAG: lcomp.cTW.ob


    open_r_matrix("lcomp.cTW.pr");  //LINETAG: lcomp.cTW.pr
    wrt_r_matrix(pred_lenc_cTW, 2, 2);  //LINETAG: lcomp.cTW.pr
    wrt_r_namevector(yrs_lenc_cTW);  //LINETAG: lcomp.cTW.pr
    wrt_r_namevector(lenbins);  //LINETAG: lcomp.cTW.pr
    close_r_matrix();  //LINETAG: lcomp.cTW.pr




	

	



	
    open_r_matrix("acomp.cHL.ob");  //LINETAG: acomp.cHL.ob
    wrt_r_matrix(obs_agec_cHL, 2, 2);  //LINETAG: acomp.cHL.ob
    wrt_r_namevector(yrs_agec_cHL);  //LINETAG: acomp.cHL.ob
    wrt_r_namevector(agebins_agec);  //LINETAG: acomp.cHL.ob
    close_r_matrix();  //LINETAG: acomp.cHL.ob

    open_r_matrix("acomp.cHL.pr");  //LINETAG: acomp.cHL.pr
    wrt_r_matrix(pred_agec_cHL, 2, 2);  //LINETAG: acomp.cHL.pr
    wrt_r_namevector(yrs_agec_cHL);  //LINETAG: acomp.cHL.pr
    wrt_r_namevector(agebins_agec);  //LINETAG: acomp.cHL.pr
    close_r_matrix();  //LINETAG: acomp.cHL.pr

    open_r_matrix("acomp.rHB.ob");  //LINETAG: acomp.rHB.ob
    wrt_r_matrix(obs_agec_rHB, 2, 2);  //LINETAG: acomp.rHB.ob
    wrt_r_namevector(yrs_agec_rHB);  //LINETAG: acomp.rHB.ob
    wrt_r_namevector(agebins_agec);  //LINETAG: acomp.rHB.ob
    close_r_matrix();  //LINETAG: acomp.rHB.ob

    open_r_matrix("acomp.rHB.pr");  //LINETAG: acomp.rHB.pr
    wrt_r_matrix(pred_agec_rHB, 2, 2);  //LINETAG: acomp.rHB.pr
    wrt_r_namevector(yrs_agec_rHB);  //LINETAG: acomp.rHB.pr
    wrt_r_namevector(agebins_agec);  //LINETAG: acomp.rHB.pr
    close_r_matrix();  //LINETAG: acomp.rHB.pr
		


	open_r_matrix("acomp.sCT.ob");  //LINETAG: acomp.sCT.ob
    wrt_r_matrix(obs_agec_sCT, 2, 2);  //LINETAG: acomp.sCT.ob
    wrt_r_namevector(yrs_agec_sCT);  //LINETAG: acomp.sCT.ob
    wrt_r_namevector(agebins_agec);  //LINETAG: acomp.sCT.ob
    close_r_matrix();  //LINETAG: acomp.sCT.ob


    open_r_matrix("acomp.sCT.pr");  //LINETAG: acomp.sCT.pr
    wrt_r_matrix(pred_agec_sCT, 2, 2);  //LINETAG: acomp.sCT.pr
    wrt_r_namevector(yrs_agec_sCT);  //LINETAG: acomp.sCT.pr
    wrt_r_namevector(agebins_agec);  //LINETAG: acomp.sCT.pr
    close_r_matrix();  //LINETAG: acomp.sCT.pr


 close_r_list();

// DATA FRAME of time series
open_r_df("t.series", styr, (endyr+1), 2);
    wrt_r_namevector(styr,(endyr+1));
    wrt_r_df_col("year", styr,(endyr+1));
    wrt_r_df_col("F.Fmsy", FdF_msy);
	wrt_r_df_col("F.F30.ratio", FdF30);	//*.ratio extension is for FishGraph
    wrt_r_df_col("F.full", Fapex);
    wrt_r_df_col("F.cHL", F_cHL_out);
    wrt_r_df_col("F.cTW", F_cTW_out);	
    wrt_r_df_col("F.rHB", F_rHB_out);
    wrt_r_df_col("F.rGN", F_rGN_out);
    wrt_r_df_col("F.cHL.D", F_D_cHL_out);
    wrt_r_df_col("F.rHB.D", F_D_rHB_out);
    wrt_r_df_col("F.rGN.D", F_D_rGN_out);
    wrt_r_df_col("Fsum", Fsum);
    wrt_r_df_col("N", totN); //abundance at start of year
    wrt_r_df_col("recruits", rec);
    wrt_r_df_col("logR.dev", log_rec_dev_output); //places zeros in yrs deviations not estimated  //KWS
    wrt_r_df_col("SSB", SSB);
    wrt_r_df_col("SSB.SSBmsy", SdSSB_msy);
    //wrt_r_df_col("SSB.msst.M", SdSSB_msy/smsy2msstM);
    wrt_r_df_col("SSB.msst", SdSSB_msy/smsy2msstM);
	wrt_r_df_col("SSB.SSBF30", SdSSB_F30);
    wrt_r_df_col("SSB.msstF30", Sdmsst_F30);
    wrt_r_df_col("B", totB);
    wrt_r_df_col("B.B0", totB/B0);
    wrt_r_df_col("SPR.static", spr_static);

    wrt_r_df_col("total.L.klb", L_total_klb_yr);
    wrt_r_df_col("total.L.knum", L_total_knum_yr);
    wrt_r_df_col("total.D.klb", D_total_klb_yr);
    wrt_r_df_col("total.D.knum", D_total_knum_yr);


    wrt_r_df_col("U.rHB.ob", obs_cpue_rHB);
    wrt_r_df_col("U.rHB.pr", pred_cpue_rHB);
    wrt_r_df_col("cv.U.rHB", obs_cv_cpue_rHB/w_cpue_rHB);
	wrt_r_df_col("cv.unwgted.U.rHB", obs_cv_cpue_rHB);

	
	wrt_r_df_col("U.sCT.ob", obs_cpue_sCT);
    wrt_r_df_col("U.sCT.pr", pred_cpue_sCT);
   	wrt_r_df_col("cv.U.sCT", obs_cv_cpue_sCT/w_cpue_sCT);
	wrt_r_df_col("cv.unwgted.U.sCT", obs_cv_cpue_sCT);

	

    wrt_r_df_col("q.rHB", q_cpue_rHB);
    wrt_r_df_col("q.rHB.rate.mult",q_rate_fcn_cpue_rHB);
    wrt_r_df_col("q.rHB.RW.log.dev",q_RW_log_dev_cpue_rHB);


	wrt_r_df_col("q.sTV", q_cpue_sCT);
    wrt_r_df_col("q.sCT.RW.log.dev",q_RW_log_dev_sCT);
	
	wrt_r_df_col("q.DD.mult", q_DD_fcn);
    wrt_r_df_col("q.DD.B.exploitable", B_q_DD);

    wrt_r_df_col("L.cHL.ob", obs_L_cHL);
    wrt_r_df_col("L.cHL.pr", pred_L_cHL_klb);
    wrt_r_df_col("cv.L.cHL", obs_cv_L_cHL);


	wrt_r_df_col("L.cTW.ob", obs_L_cTW);
    wrt_r_df_col("L.cTW.pr", pred_L_cTW_klb);
    wrt_r_df_col("cv.L.cTW", obs_cv_L_cTW);

    wrt_r_df_col("L.rHB.ob", obs_L_rHB);
    wrt_r_df_col("L.rHB.pr", pred_L_rHB_knum);
    wrt_r_df_col("cv.L.rHB", obs_cv_L_rHB);

    wrt_r_df_col("L.rGN.ob", obs_L_rGN);
    wrt_r_df_col("L.rGN.pr", pred_L_rGN_knum);
    wrt_r_df_col("cv.L.rGN", obs_cv_L_rGN);

    wrt_r_df_col("D.cHL.ob", obs_D_cHL);
    wrt_r_df_col("D.cHL.pr", pred_D_cHL_knum);
    wrt_r_df_col("cv.D.cHL", obs_cv_D_cHL);

    wrt_r_df_col("D.rHB.ob", obs_D_rHB);
    wrt_r_df_col("D.rHB.pr", pred_D_rHB_knum);
    wrt_r_df_col("cv.D.rHB", obs_cv_D_rHB);

    wrt_r_df_col("D.rGN.ob", obs_D_rGN);
    wrt_r_df_col("D.rGN.pr", pred_D_rGN_knum);
    wrt_r_df_col("cv.D.rGN", obs_cv_D_rGN);

    //comp sample sizes
	wrt_r_df_col("lcomp.cTW.n", nsamp_lenc_cTW_allyr);
	
    wrt_r_df_col("acomp.cHL.n", nsamp_agec_cHL_allyr);
   	wrt_r_df_col("acomp.rHB.n", nsamp_agec_rHB_allyr);
	wrt_r_df_col("acomp.sCT.n", nsamp_agec_sCT_allyr);

	wrt_r_df_col("lcomp.cTW.nfish", nfish_lenc_cTW_allyr);
	
    wrt_r_df_col("acomp.cHL.nfish", nfish_agec_cHL_allyr);
   	wrt_r_df_col("acomp.rHB.nfish", nfish_agec_rHB_allyr);
	wrt_r_df_col("acomp.sCT.nfish", nfish_agec_sCT_allyr);
	
	wrt_r_df_col("lcomp.cTW.neff", (1+nsamp_lenc_cTW_allyr*exp(log_dm_lenc_cTW_out(8)))/(1+exp(log_dm_lenc_cTW_out(8))) );
	
	wrt_r_df_col("acomp.cHL.neff", (1+nsamp_agec_cHL_allyr*exp(log_dm_agec_cHL_out(8)))/(1+exp(log_dm_agec_cHL_out(8))) );
	wrt_r_df_col("acomp.rHB.neff", (1+nsamp_agec_rHB_allyr*exp(log_dm_agec_rHB_out(8)))/(1+exp(log_dm_agec_rHB_out(8))) );
	wrt_r_df_col("acomp.sCT.neff", (1+nsamp_agec_sCT_allyr*exp(log_dm_agec_sCT_out(8)))/(1+exp(log_dm_agec_sCT_out(8))) );	
 	
	
 close_r_df();

// DATA FRAME of L and D time series by fishery
open_r_df("LD.pr.tseries", styr, endyr, 2);
	wrt_r_namevector(styr,endyr);
	wrt_r_df_col("year", styr, endyr);

    wrt_r_df_col("L.cHL.klb", pred_L_cHL_klb);
    wrt_r_df_col("L.cHL.knum", pred_L_cHL_knum);

    wrt_r_df_col("L.cTW.klb", pred_L_cTW_klb);
    wrt_r_df_col("L.cTW.knum", pred_L_cTW_knum);	
    wrt_r_df_col("L.rHB.klb", pred_L_rHB_klb);
    wrt_r_df_col("L.rHB.knum", pred_L_rHB_knum);
    wrt_r_df_col("L.rGN.klb", pred_L_rGN_klb);
    wrt_r_df_col("L.rGN.knum", pred_L_rGN_knum);

    wrt_r_df_col("D.cHL.klb", pred_D_cHL_klb);
    wrt_r_df_col("D.cHL.knum", pred_D_cHL_knum);
    wrt_r_df_col("D.rHB.klb", pred_D_rHB_klb);
    wrt_r_df_col("D.rHB.knum", pred_D_rHB_knum);
    wrt_r_df_col("D.rGN.klb", pred_D_rGN_klb);
    wrt_r_df_col("D.rGN.knum", pred_D_rGN_knum);

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
 	wrt_r_df_col("wholewgt.wgted.D.klb", wgt_wgted_D_klb);
    wrt_r_df_col("prop.female", prop_f);
    wrt_r_df_col("prop.male", prop_m);	
   	wrt_r_df_col("mat.fem.endyr", maturity_f(endyr));
   	wrt_r_df_col("mat.male", maturity_m);	
 	wrt_r_df_col("reprod.endyr", reprod(endyr));
 	
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
 	wrt_r_df_col("D.eq.wholeklb", D_eq_klb);
 	wrt_r_df_col("D.eq.knum", D_eq_knum);
 close_r_df();

 open_r_df("pr.series", 1, n_iter_spr, 2);
 	wrt_r_namevector(1,n_iter_spr);
 	wrt_r_df_col("F.spr", F_spr);
 	wrt_r_df_col("spr", spr_spr);
 	wrt_r_df_col("SPR", spr_ratio);
 	wrt_r_df_col("ypr.lb.whole", L_spr); //whole weight
 close_r_df();


 open_r_list("CLD.est.mats");

     open_r_matrix("Lw.cHL");   //LINETAG: Lw.cHL
         wrt_r_matrix(L_cHL_klb, 1,1);   //LINETAG: Lw.cHL
     close_r_matrix();   //LINETAG: Lw.cHL


	 open_r_matrix("Lw.cTW");   //LINETAG: Lw.cTW
         wrt_r_matrix(L_cTW_klb, 1,1);   //LINETAG: Lw.cTW
     close_r_matrix();   //LINETAG: Lw.cTW
	 
     open_r_matrix("Lw.rHB");   //LINETAG: Lw.rHB
         wrt_r_matrix(L_rHB_klb, 1,1);   //LINETAG: Lw.rHB
     close_r_matrix();   //LINETAG: Lw.rHB

     open_r_matrix("Lw.rGN");   //LINETAG: Lw.rGN
         wrt_r_matrix(L_rGN_klb, 1,1);   //LINETAG: Lw.rGN
     close_r_matrix();   //LINETAG: Lw.rGN

     open_r_matrix("Lw.total");   //LINETAG: Lw.total
         wrt_r_matrix(L_total_klb, 1,1);   //LINETAG: Lw.total
     close_r_matrix();   //LINETAG: Lw.total


     open_r_matrix("Ln.cHL");   //LINETAG: Ln.cHL
         wrt_r_matrix(L_cHL_num, 1,1);   //LINETAG: Ln.cHL
     close_r_matrix();   //LINETAG: Ln.cHL


	 open_r_matrix("Ln.cTW");   //LINETAG: Ln.cTW
         wrt_r_matrix(L_cTW_num, 1,1);   //LINETAG: Ln.cTW
     close_r_matrix();   //LINETAG: Ln.cTW
	 
     open_r_matrix("Ln.rHB");   //LINETAG: Ln.rHB
         wrt_r_matrix(L_rHB_num, 1,1);   //LINETAG: Ln.rHB
     close_r_matrix();   //LINETAG: Ln.rHB

     open_r_matrix("Ln.rGN");   //LINETAG: Ln.rGN
         wrt_r_matrix(L_rGN_num, 1,1);   //LINETAG: Ln.rGN
     close_r_matrix();   //LINETAG: Ln.rGN

     open_r_matrix("Ln.total");   //LINETAG: Ln.total
         wrt_r_matrix(L_total_num, 1,1);   //LINETAG: Ln.total
     close_r_matrix();   //LINETAG: Ln.total


     open_r_matrix("Dw.cHL");   //LINETAG: Dw.cHL
         wrt_r_matrix(D_cHL_klb, 1,1);   //LINETAG: Dw.cHL
     close_r_matrix();   //LINETAG: Dw.cHL

     open_r_matrix("Dw.rHB");   //LINETAG: Dw.rHB
         wrt_r_matrix(D_rHB_klb, 1,1);   //LINETAG: Dw.rHB
     close_r_matrix();   //LINETAG: Dw.rHB

     open_r_matrix("Dw.rGN");   //LINETAG: Dw.rGN
         wrt_r_matrix(D_rGN_klb, 1,1);   //LINETAG: Dw.rGN
     close_r_matrix();   //LINETAG: Dw.rGN

     open_r_matrix("Dw.total");   //LINETAG: Dw.total
         wrt_r_matrix(D_total_klb, 1,1);   //LINETAG: Dw.total
     close_r_matrix();   //LINETAG: Dw.total


     open_r_matrix("Dn.cHL");   //LINETAG: Dn.cHL
         wrt_r_matrix(D_cHL_num, 1,1);   //LINETAG: Dn.cHL
     close_r_matrix();   //LINETAG: Dn.cHL

     open_r_matrix("Dn.rHB");   //LINETAG: Dn.rHB
         wrt_r_matrix(D_rHB_num, 1,1);   //LINETAG: Dn.rHB
     close_r_matrix();   //LINETAG: Dn.rHB

     open_r_matrix("Dn.rGN");   //LINETAG: Dn.rGN
         wrt_r_matrix(D_rGN_num, 1,1);   //LINETAG: Dn.rGN
     close_r_matrix();   //LINETAG: Dn.rGN

     open_r_matrix("Dn.total");   //LINETAG: Dn.total
         wrt_r_matrix(D_total_num, 1,1);   //LINETAG: Dn.total
     close_r_matrix();   //LINETAG: Dn.total

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
	wrt_r_df_col("D.knum.proj", D_knum_proj);
	wrt_r_df_col("D.klb.proj", D_klb_proj);	
 close_r_df();

close_r_file();
