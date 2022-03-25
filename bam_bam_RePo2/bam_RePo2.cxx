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
	wrt_r_item("D.mort.cHl", Dmort_cHl);
	wrt_r_item("D.mort.rHb", Dmort_rHb);
	wrt_r_item("D.mort.rGe", Dmort_rGe);
	wrt_r_item("q.rHb", mfexp(log_q_cpue_rHb));
	wrt_r_item("q.sCT", mfexp(log_q_cpue_sCT));

    wrt_r_item("q.rate",q_rate);
	wrt_r_item("q.beta",q_DD_beta);
	wrt_r_item("q.DD.B0.exploitable", B0_q_DD);
	wrt_r_item("F.prop.L.cHl", F_prop_L_cHl);
	wrt_r_item("F.prop.L.cTw", F_prop_L_cTw);	
	wrt_r_item("F.prop.L.rHb", F_prop_L_rHb);
	wrt_r_item("F.prop.L.rGe", F_prop_L_rGe);
	wrt_r_item("F.prop.D.cHl", F_prop_D_cHl);
	wrt_r_item("F.prop.D.rHb", F_prop_D_rHb);
	wrt_r_item("F.prop.D.rGe", F_prop_D_rGe);
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
	wrt_r_df_col("log_dm_lenc_cTw",log_dm_lenc_cTw_out);
	wrt_r_df_col("log_dm_agec_cHl",log_dm_agec_cHl_out);
	wrt_r_df_col("log_dm_agec_rHb",log_dm_agec_rHb_out);
	wrt_r_df_col("log_dm_agec_sCT",log_dm_agec_sCT_out);
	
    wrt_r_df_col("A50_sel_cHl2",A50_sel_cHl2_out);
    wrt_r_df_col("slope_sel_cHl2",slope_sel_cHl2_out);
    wrt_r_df_col("A50_sel_cHl3",A50_sel_cHl3_out);
    wrt_r_df_col("slope_sel_cHl3",slope_sel_cHl3_out);
	

    wrt_r_df_col("A50_sel_cTw",A50_sel_cTw_out);
    wrt_r_df_col("slope_sel_cTw",slope_sel_cTw_out);
	
    wrt_r_df_col("A50_sel_rHb1",A50_sel_rHb1_out);
    wrt_r_df_col("slope_sel_rHb1",slope_sel_rHb1_out);
 	wrt_r_df_col("A50_sel_rHb2",A50_sel_rHb2_out);
    wrt_r_df_col("slope_sel_rHb2",slope_sel_rHb2_out);
	wrt_r_df_col("A50_sel_rHb3",A50_sel_rHb3_out);
    wrt_r_df_col("slope_sel_rHb3",slope_sel_rHb3_out);

		
    wrt_r_df_col("A50_sel_sCT",A50_sel_sCT_out);
    wrt_r_df_col("slope_sel_sCT",slope_sel_sCT_out);

	
    
    wrt_r_df_col("log_q_cpue_rHb",log_q_cpue_rHb_out);
	wrt_r_df_col("log_q_cpue_sCT",log_q_cpue_sCT_out);
    wrt_r_df_col("M_constant",M_constant_out);
    wrt_r_df_col("log_avg_F_L_cHl",log_avg_F_L_cHl_out);
	wrt_r_df_col("log_avg_F_L_cTw",log_avg_F_L_cTw_out);
    wrt_r_df_col("log_avg_F_L_rHb",log_avg_F_L_rHb_out);
    wrt_r_df_col("log_avg_F_L_rGe",log_avg_F_L_rGe_out);
    wrt_r_df_col("log_avg_F_D_cHl",log_avg_F_D_cHl_out);
    wrt_r_df_col("log_avg_F_D_rHb",log_avg_F_D_rHb_out);
    wrt_r_df_col("log_avg_F_D_rGe",log_avg_F_D_rGe_out);
    //wrt_r_df_col("F_init",F_init_out);
close_r_df();

// DATA FRAME of time series deviation vector estimates
// names used in this object must match the names used in the "parm.tvec.cons" object
open_r_df("parm.tvec", styr, endyr, 2);
	wrt_r_namevector(styr,endyr);
	wrt_r_df_col("year", styr,endyr);
    wrt_r_df_col("log.rec.dev", log_rec_dev_out);
    wrt_r_df_col("log.F.dev.L.cHl", log_F_dev_L_cHl_out);
	wrt_r_df_col("log.F.dev.L.cTw", log_F_dev_L_cTw_out);
    wrt_r_df_col("log.F.dev.L.rHb", log_F_dev_L_rHb_out);
    wrt_r_df_col("log.F.dev.L.rGe", log_F_dev_L_rGe_out);
    wrt_r_df_col("log.F.dev.D.cHl", log_F_dev_D_cHl_out);
    wrt_r_df_col("log.F.dev.D.rHb", log_F_dev_D_rHb_out);
    wrt_r_df_col("log.F.dev.D.rGe", log_F_dev_D_rGe_out);
	wrt_r_df_col("log.q.dev.rHb.RWq", q_RW_log_dev_cpue_rHb);
	
close_r_df();

// MATRIX object of deviation vector constraints
// names used in this object must match the names used in the "parm.tvec" object
open_r_df("parm.tvec.cons",1,3,2);
    wrt_r_namevector(1,3);
    wrt_r_df_col("log.rec.dev",set_log_rec_dev);
    wrt_r_df_col("log.F.dev.L.cHl",set_log_F_dev_L_cHl);
	wrt_r_df_col("log.F.dev.L.cTw",set_log_F_dev_L_cTw);
    wrt_r_df_col("log.F.dev.L.rHb",set_log_F_dev_L_rHb);
    wrt_r_df_col("log.F.dev.L.rGe",set_log_F_dev_L_rGe);
    wrt_r_df_col("log.F.dev.D.cHl",set_log_F_dev_D_cHl);
    wrt_r_df_col("log.F.dev.D.rHb",set_log_F_dev_D_rHb);
    wrt_r_df_col("log.F.dev.D.rGe",set_log_F_dev_D_rGe);
	wrt_r_df_col("log.q.dev.rHb.RWq", set_log_RWq_dev);
	wrt_r_df_col("log.q.dev.cHl.RWq", set_log_RWq_dev);
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
    wrt_r_df_col("log.Nage.dev",set_log_Nage_dev);
close_r_df();

// VECTOR object of SDNR calculations
open_r_vector("sdnr");
    wrt_r_item("sdnr.U.rHb", sdnr_cpue_rHb);
	wrt_r_item("sdnr.U.sCT", sdnr_cpue_sCT);
	
//     wrt_r_item("sdnr.lc.cTw", sdnr_lenc_cTw);	
	
    wrt_r_item("sdnr.ac.cHl", sdnr_agec_cHl);
    wrt_r_item("sdnr.ac.rHb", sdnr_agec_rHb);
	wrt_r_item("sdnr.ac.sCT", sdnr_agec_sCT);	
	
 close_r_vector();

// VECTOR object of likelihood contributions
open_r_vector("like");
    wrt_r_item("lk.total", fval); //weighted likelihood
    wrt_r_item("lk.unwgt.data", fval_data);  //likelihood of just data components

    wrt_r_item("lk.L.cHl", f_L_cHl);
    wrt_r_item("lk.L.cTw", f_L_cTw);	
    wrt_r_item("lk.L.rHb", f_L_rHb);
    wrt_r_item("lk.L.rGe", f_L_rGe);
    wrt_r_item("lk.D.cHl", f_D_cHl);
    wrt_r_item("lk.D.rHb", f_D_rHb);
    wrt_r_item("lk.D.rGe", f_D_rGe);

    wrt_r_item("lk.U.rHb", f_cpue_rHb);
	wrt_r_item("lk.U.sCT", f_cpue_sCT);
    wrt_r_item("lk.lenc.cTw", f_lenc_cTw);	
    wrt_r_item("lk.agec.cHl", f_agec_cHl);
    wrt_r_item("lk.agec.rHb", f_agec_rHb);
	wrt_r_item("lk.agec.sCT", f_agec_sCT);	
    wrt_r_item("lk.Nage.init", f_Nage_init);
    wrt_r_item("lk.SRfit", f_rec_dev);
    wrt_r_item("lk.SRearly", f_rec_dev_early);
    wrt_r_item("lk.SRend", f_rec_dev_end);
    wrt_r_item("lk.fullF", f_fullF_constraint);
    wrt_r_item("lk.Ftune", f_Ftune);
	wrt_r_item("lk.U.rHb.RWq", f_RWq_cpue_rHb);
    wrt_r_item("lk.priors",f_priors);
    wrt_r_item("gradient.max",grad_max);

    wrt_r_item("w.L", w_L);
    wrt_r_item("w.D", w_D);
    wrt_r_item("w.U.rHb", w_cpue_rHb);
	wrt_r_item("w.U.sCT", w_cpue_sCT);
	
    wrt_r_item("w.lc.cTw", w_lenc_cTw);	
	
    wrt_r_item("w.ac.cHl", w_agec_cHl);
    wrt_r_item("w.ac.rHb", w_agec_rHb);
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

    open_r_matrix("len.cHl.mm"); //LINETAG: len.cHl.mm
    wrt_r_matrix(len_cHl_mm, 2, 2); //LINETAG: len.cHl.mm
    wrt_r_namevector(styr,endyr); //LINETAG: len.cHl.mm
    wrt_r_namevector(agebins); //LINETAG: len.cHl.mm
    close_r_matrix(); //LINETAG: len.cHl.mm


    open_r_matrix("len.cTw.mm"); //LINETAG: len.cTw.mm
    wrt_r_matrix(len_cTw_mm, 2, 2); //LINETAG: len.cTw.mm
    wrt_r_namevector(styr,endyr); //LINETAG: len.cTw.mm
    wrt_r_namevector(agebins); //LINETAG: len.cTw.mm
    close_r_matrix(); //LINETAG: len.cTw.mm
	
    open_r_matrix("len.rHb.mm"); //LINETAG: len.rHb.mm
    wrt_r_matrix(len_rHb_mm, 2, 2); //LINETAG: len.rHb.mm
    wrt_r_namevector(styr,endyr); //LINETAG: len.rHb.mm
    wrt_r_namevector(agebins); //LINETAG: len.rHb.mm
    close_r_matrix(); //LINETAG: len.rHb.mm

    open_r_matrix("len.rGe.mm"); //LINETAG: len.rGe.mm
    wrt_r_matrix(len_rGe_mm, 2, 2); //LINETAG: len.rGe.mm
    wrt_r_namevector(styr,endyr); //LINETAG: len.rGe.mm
    wrt_r_namevector(agebins); //LINETAG: len.rGe.mm
    close_r_matrix(); //LINETAG: len.rGe.mm

    open_r_matrix("len.cHl.D.mm"); //LINETAG: len.cHl.D.mm
    wrt_r_matrix(len_D_cHl_mm, 2, 2);  //LINETAG: len.cHl.D.mm
    wrt_r_namevector(styr,endyr);  //LINETAG: len.cHl.D.mm
    wrt_r_namevector(agebins);  //LINETAG: len.cHl.D.mm
    close_r_matrix();  //LINETAG: len.cHl.D.mm

    open_r_matrix("len.rHb.D.mm");  //LINETAG: len.rHb.D.mm
    wrt_r_matrix(len_D_rHb_mm, 2, 2);  //LINETAG: len.rHb.D.mm
    wrt_r_namevector(styr,endyr);  //LINETAG: len.rHb.D.mm
    wrt_r_namevector(agebins);  //LINETAG: len.rHb.D.mm
    close_r_matrix();  //LINETAG: len.rHb.D.mm

    open_r_matrix("len.rGe.D.mm");  //LINETAG: len.rGe.D.mm
    wrt_r_matrix(len_D_rGe_mm, 2, 2);  //LINETAG: len.rGe.D.mm
    wrt_r_namevector(styr,endyr);  //LINETAG: len.rGe.D.mm
    wrt_r_namevector(agebins);  //LINETAG: len.rGe.D.mm
    close_r_matrix();  //LINETAG: len.rGe.D.mm



    open_r_matrix("wholewgt.cHl.lb");  //LINETAG: wholewgt.cHl.lb
    wrt_r_matrix(wholewgt_cHl_klb*1000.0, 2, 2);  //LINETAG: wholewgt.cHl.lb
    wrt_r_namevector(styr,endyr);  //LINETAG: wholewgt.cHl.lb
    wrt_r_namevector(agebins);  //LINETAG: wholewgt.cHl.lb
    close_r_matrix();  //LINETAG: wholewgt.cHl.lb


	open_r_matrix("wholewgt.cTw.lb");  //LINETAG: wholewgt.cTw.lb
    wrt_r_matrix(wholewgt_cTw_klb*1000.0, 2, 2);  //LINETAG: wholewgt.cTw.lb
    wrt_r_namevector(styr,endyr);  //LINETAG: wholewgt.cTw.lb
    wrt_r_namevector(agebins);  //LINETAG: wholewgt.cTw.lb
    close_r_matrix();  //LINETAG: wholewgt.cTw.lb

    open_r_matrix("wholewgt.rHb.lb");  //LINETAG: wholewgt.rHb.lb
    wrt_r_matrix(wholewgt_rHb_klb*1000.0, 2, 2);  //LINETAG: wholewgt.rHb.lb
    wrt_r_namevector(styr,endyr);  //LINETAG: wholewgt.rHb.lb
    wrt_r_namevector(agebins);  //LINETAG: wholewgt.rHb.lb
    close_r_matrix();  //LINETAG: wholewgt.rHb.lb

    open_r_matrix("wholewgt.rGe.lb");  //LINETAG: wholewgt.rGe.lb
    wrt_r_matrix(wholewgt_rGe_klb*1000.0, 2, 2);  //LINETAG: wholewgt.rGe.lb
    wrt_r_namevector(styr,endyr);  //LINETAG: wholewgt.rGe.lb
    wrt_r_namevector(agebins);  //LINETAG: wholewgt.rGe.lb
    close_r_matrix();  //LINETAG: wholewgt.rGe.lb

    open_r_matrix("wholewgt.cHl.D.lb");  //LINETAG: wholewgt.cHl.D.lb
    wrt_r_matrix(wholewgt_D_cHl_klb*1000.0, 2, 2);  //LINETAG: wholewgt.cHl.D.lb
    wrt_r_namevector(styr,endyr);  //LINETAG: wholewgt.cHl.D.lb
    wrt_r_namevector(agebins);  //LINETAG: wholewgt.cHl.D.lb
    close_r_matrix();  //LINETAG: wholewgt.cHl.D.lb

    open_r_matrix("wholewgt.rHb.D.lb");  //LINETAG: wholewgt.rHb.D.lb
    wrt_r_matrix(wholewgt_D_rHb_klb*1000.0, 2, 2);  //LINETAG: wholewgt.rHb.D.lb
    wrt_r_namevector(styr,endyr);  //LINETAG: wholewgt.rHb.D.lb
    wrt_r_namevector(agebins);  //LINETAG: wholewgt.rHb.D.lb
    close_r_matrix();  //LINETAG: wholewgt.rHb.D.lb

    open_r_matrix("wholewgt.rGe.D.lb");  //LINETAG: wholewgt.rGe.D.lb
    wrt_r_matrix(wholewgt_D_rGe_klb*1000.0, 2, 2);  //LINETAG: wholewgt.rGe.D.lb
    wrt_r_namevector(styr,endyr);  //LINETAG: wholewgt.rGe.D.lb
    wrt_r_namevector(agebins);  //LINETAG: wholewgt.rGe.D.lb
    close_r_matrix();  //LINETAG: wholewgt.rGe.D.lb

 close_r_list();

open_r_list("sel.age");

    wrt_r_complete_vector("sel.v.wgted.L",sel_wgted_L, agebins);
    wrt_r_complete_vector("sel.v.wgted.D",sel_wgted_D, agebins);
    wrt_r_complete_vector("sel.v.wgted.tot",sel_wgted_tot, agebins);

    open_r_matrix("sel.m.cHl");  //LINETAG: sel.m.cHl
    wrt_r_matrix(sel_cHl, 2, 2);  //LINETAG: sel.m.cHl
    wrt_r_namevector(styr,endyr);  //LINETAG: sel.m.cHl
    wrt_r_namevector(agebins);  //LINETAG: sel.m.cHl
    close_r_matrix();  //LINETAG: sel.m.cHl


	open_r_matrix("sel.m.cTw");  //LINETAG: sel.m.cTw
    wrt_r_matrix(sel_cTw, 2, 2);  //LINETAG: sel.m.cTw
    wrt_r_namevector(styr,endyr);  //LINETAG: sel.m.cTw
    wrt_r_namevector(agebins);  //LINETAG: sel.m.cTw
    close_r_matrix();  //LINETAG: sel.m.cTw

    open_r_matrix("sel.m.rHb");  //LINETAG: sel.m.rHb
    wrt_r_matrix(sel_rHb, 2, 2);  //LINETAG: sel.m.rHb
    wrt_r_namevector(styr,endyr);  //LINETAG: sel.m.rHb
    wrt_r_namevector(agebins);  //LINETAG: sel.m.rHb
    close_r_matrix();  //LINETAG: sel.m.rHb

    open_r_matrix("sel.m.rGe");  //LINETAG: sel.m.rGe
    wrt_r_matrix(sel_rGe, 2, 2);  //LINETAG: sel.m.rGe
    wrt_r_namevector(styr,endyr);  //LINETAG: sel.m.rGe
    wrt_r_namevector(agebins);  //LINETAG: sel.m.rGe
    close_r_matrix();  //LINETAG: sel.m.rGe

    open_r_matrix("sel.m.cHl.D");  //LINETAG: sel.m.cHl.D
    wrt_r_matrix(sel_D_cHl, 2, 2);  //LINETAG: sel.m.cHl.D
    wrt_r_namevector(styr,endyr);  //LINETAG: sel.m.cHl.D
    wrt_r_namevector(agebins);  //LINETAG: sel.m.cHl.D
    close_r_matrix();  //LINETAG: sel.m.cHl.D
    open_r_matrix("sel.m.rHb.D");  //LINETAG: sel.m.rHb.D
    wrt_r_matrix(sel_D_rHb, 2, 2);  //LINETAG: sel.m.rHb.D
    wrt_r_namevector(styr,endyr);  //LINETAG: sel.m.rHb.D
    wrt_r_namevector(agebins);  //LINETAG: sel.m.rHb.D
    close_r_matrix();  //LINETAG: sel.m.rHb.D
    open_r_matrix("sel.m.rGe.D");  //LINETAG: sel.m.rGe.D
    wrt_r_matrix(sel_D_rGe, 2, 2);  //LINETAG: sel.m.rGe.D
    wrt_r_namevector(styr,endyr);  //LINETAG: sel.m.rGe.D
    wrt_r_namevector(agebins);  //LINETAG: sel.m.rGe.D
    close_r_matrix();  //LINETAG: sel.m.rGe.D
	
	open_r_matrix("sel.m.sCT");  //LINETAG: sel.m.sCT
    wrt_r_matrix(sel_sCT, 2, 2);  //LINETAG: sel.m.sCT
    wrt_r_namevector(styr,endyr);  //LINETAG: sel.m.sCT
    wrt_r_namevector(agebins);  //LINETAG: sel.m.sCT
    close_r_matrix();  //LINETAG: sel.m.sCT

	
 close_r_list();


//LIST object with predicted and observed composition data
open_r_list("comp.mats");




    open_r_matrix("lcomp.cTw.ob");  //LINETAG: lcomp.cTw.ob
    wrt_r_matrix(obs_lenc_cTw, 2, 2);  //LINETAG: lcomp.cTw.ob
    wrt_r_namevector(yrs_lenc_cTw);  //LINETAG: lcomp.cTw.ob
    wrt_r_namevector(lenbins);  //LINETAG: lcomp.cTw.ob
    close_r_matrix();  //LINETAG: lcomp.cTw.ob


    open_r_matrix("lcomp.cTw.pr");  //LINETAG: lcomp.cTw.pr
    wrt_r_matrix(pred_lenc_cTw, 2, 2);  //LINETAG: lcomp.cTw.pr
    wrt_r_namevector(yrs_lenc_cTw);  //LINETAG: lcomp.cTw.pr
    wrt_r_namevector(lenbins);  //LINETAG: lcomp.cTw.pr
    close_r_matrix();  //LINETAG: lcomp.cTw.pr




	

	



	
    open_r_matrix("acomp.cHl.ob");  //LINETAG: acomp.cHl.ob
    wrt_r_matrix(obs_agec_cHl, 2, 2);  //LINETAG: acomp.cHl.ob
    wrt_r_namevector(yrs_agec_cHl);  //LINETAG: acomp.cHl.ob
    wrt_r_namevector(agebins_agec);  //LINETAG: acomp.cHl.ob
    close_r_matrix();  //LINETAG: acomp.cHl.ob

    open_r_matrix("acomp.cHl.pr");  //LINETAG: acomp.cHl.pr
    wrt_r_matrix(pred_agec_cHl, 2, 2);  //LINETAG: acomp.cHl.pr
    wrt_r_namevector(yrs_agec_cHl);  //LINETAG: acomp.cHl.pr
    wrt_r_namevector(agebins_agec);  //LINETAG: acomp.cHl.pr
    close_r_matrix();  //LINETAG: acomp.cHl.pr

    open_r_matrix("acomp.rHb.ob");  //LINETAG: acomp.rHb.ob
    wrt_r_matrix(obs_agec_rHb, 2, 2);  //LINETAG: acomp.rHb.ob
    wrt_r_namevector(yrs_agec_rHb);  //LINETAG: acomp.rHb.ob
    wrt_r_namevector(agebins_agec);  //LINETAG: acomp.rHb.ob
    close_r_matrix();  //LINETAG: acomp.rHb.ob

    open_r_matrix("acomp.rHb.pr");  //LINETAG: acomp.rHb.pr
    wrt_r_matrix(pred_agec_rHb, 2, 2);  //LINETAG: acomp.rHb.pr
    wrt_r_namevector(yrs_agec_rHb);  //LINETAG: acomp.rHb.pr
    wrt_r_namevector(agebins_agec);  //LINETAG: acomp.rHb.pr
    close_r_matrix();  //LINETAG: acomp.rHb.pr
		


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
    wrt_r_df_col("F.cHl", F_cHl_out);
    wrt_r_df_col("F.cTw", F_cTw_out);	
    wrt_r_df_col("F.rHb", F_rHb_out);
    wrt_r_df_col("F.rGe", F_rGe_out);
    wrt_r_df_col("F.cHl.D", F_D_cHl_out);
    wrt_r_df_col("F.rHb.D", F_D_rHb_out);
    wrt_r_df_col("F.rGe.D", F_D_rGe_out);
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


    wrt_r_df_col("U.rHb.ob", obs_cpue_rHb);
    wrt_r_df_col("U.rHb.pr", pred_cpue_rHb);
    wrt_r_df_col("cv.U.rHb", cv_cpue_rHb/w_cpue_rHb);
	wrt_r_df_col("cv.unwgted.U.rHb", cv_cpue_rHb);

	
	wrt_r_df_col("U.sCT.ob", obs_cpue_sCT);
    wrt_r_df_col("U.sCT.pr", pred_cpue_sCT);
   	wrt_r_df_col("cv.U.sCT", cv_cpue_sCT/w_cpue_sCT);
	wrt_r_df_col("cv.unwgted.U.sCT", cv_cpue_sCT);

	

    wrt_r_df_col("q.rHb", q_cpue_rHb);
    wrt_r_df_col("q.rHb.rate.mult",q_rate_fcn_cpue_rHb);
    wrt_r_df_col("q.rHb.RW.log.dev",q_RW_log_dev_cpue_rHb);


	wrt_r_df_col("q.CVID", q_cpue_sCT);
    wrt_r_df_col("q.sCT.RW.log.dev",q_RW_log_dev_sCT);
	
	wrt_r_df_col("q.DD.mult", q_DD_fcn);
    wrt_r_df_col("q.DD.B.exploitable", B_q_DD);

    wrt_r_df_col("L.cHl.ob", obs_L_cHl);
    wrt_r_df_col("L.cHl.pr", pred_L_cHl_klb);
    wrt_r_df_col("cv.L.cHl", cv_L_cHl);


	wrt_r_df_col("L.cTw.ob", obs_L_cTw);
    wrt_r_df_col("L.cTw.pr", pred_L_cTw_klb);
    wrt_r_df_col("cv.L.cTw", cv_L_cTw);

    wrt_r_df_col("L.rHb.ob", obs_L_rHb);
    wrt_r_df_col("L.rHb.pr", pred_L_rHb_knum);
    wrt_r_df_col("cv.L.rHb", cv_L_rHb);

    wrt_r_df_col("L.rGe.ob", obs_L_rGe);
    wrt_r_df_col("L.rGe.pr", pred_L_rGe_knum);
    wrt_r_df_col("cv.L.rGe", cv_L_rGe);

    wrt_r_df_col("D.cHl.ob", obs_D_cHl);
    wrt_r_df_col("D.cHl.pr", pred_D_cHl_knum);
    wrt_r_df_col("cv.D.cHl", cv_released_cHl);

    wrt_r_df_col("D.rHb.ob", obs_D_rHb);
    wrt_r_df_col("D.rHb.pr", pred_D_rHb_knum);
    wrt_r_df_col("cv.D.rHb", cv_released_rHb);

    wrt_r_df_col("D.rGe.ob", obs_D_rGe);
    wrt_r_df_col("D.rGe.pr", pred_D_rGe_knum);
    wrt_r_df_col("cv.D.rGe", cv_released_rGe);

    //comp sample sizes
	wrt_r_df_col("lcomp.cTw.n", nsamp_lenc_cTw_allyr);
	
    wrt_r_df_col("acomp.cHl.n", nsamp_agec_cHl_allyr);
   	wrt_r_df_col("acomp.rHb.n", nsamp_agec_rHb_allyr);
	wrt_r_df_col("acomp.sCT.n", nsamp_agec_sCT_allyr);

	wrt_r_df_col("lcomp.cTw.nfish", nfish_lenc_cTw_allyr);
	
    wrt_r_df_col("acomp.cHl.nfish", nfish_agec_cHl_allyr);
   	wrt_r_df_col("acomp.rHb.nfish", nfish_agec_rHb_allyr);
	wrt_r_df_col("acomp.sCT.nfish", nfish_agec_sCT_allyr);
	
	wrt_r_df_col("lcomp.cTw.neff", (1+nsamp_lenc_cTw_allyr*exp(log_dm_lenc_cTw_out(8)))/(1+exp(log_dm_lenc_cTw_out(8))) );
	
	wrt_r_df_col("acomp.cHl.neff", (1+nsamp_agec_cHl_allyr*exp(log_dm_agec_cHl_out(8)))/(1+exp(log_dm_agec_cHl_out(8))) );
	wrt_r_df_col("acomp.rHb.neff", (1+nsamp_agec_rHb_allyr*exp(log_dm_agec_rHb_out(8)))/(1+exp(log_dm_agec_rHb_out(8))) );
	wrt_r_df_col("acomp.sCT.neff", (1+nsamp_agec_sCT_allyr*exp(log_dm_agec_sCT_out(8)))/(1+exp(log_dm_agec_sCT_out(8))) );	
 	
	
 close_r_df();

// DATA FRAME of L and D time series by fishery
open_r_df("LD.pr.tseries", styr, endyr, 2);
	wrt_r_namevector(styr,endyr);
	wrt_r_df_col("year", styr, endyr);

    wrt_r_df_col("L.cHl.klb", pred_L_cHl_klb);
    wrt_r_df_col("L.cHl.knum", pred_L_cHl_knum);

    wrt_r_df_col("L.cTw.klb", pred_L_cTw_klb);
    wrt_r_df_col("L.cTw.knum", pred_L_cTw_knum);	
    wrt_r_df_col("L.rHb.klb", pred_L_rHb_klb);
    wrt_r_df_col("L.rHb.knum", pred_L_rHb_knum);
    wrt_r_df_col("L.rGe.klb", pred_L_rGe_klb);
    wrt_r_df_col("L.rGe.knum", pred_L_rGe_knum);

    wrt_r_df_col("D.cHl.klb", pred_D_cHl_klb);
    wrt_r_df_col("D.cHl.knum", pred_D_cHl_knum);
    wrt_r_df_col("D.rHb.klb", pred_D_rHb_klb);
    wrt_r_df_col("D.rHb.knum", pred_D_rHb_knum);
    wrt_r_df_col("D.rGe.klb", pred_D_rGe_klb);
    wrt_r_df_col("D.rGe.knum", pred_D_rGe_knum);

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

     open_r_matrix("Lw.cHl");   //LINETAG: Lw.cHl
         wrt_r_matrix(L_cHl_klb, 1,1);   //LINETAG: Lw.cHl
     close_r_matrix();   //LINETAG: Lw.cHl


	 open_r_matrix("Lw.cTw");   //LINETAG: Lw.cTw
         wrt_r_matrix(L_cTw_klb, 1,1);   //LINETAG: Lw.cTw
     close_r_matrix();   //LINETAG: Lw.cTw
	 
     open_r_matrix("Lw.rHb");   //LINETAG: Lw.rHb
         wrt_r_matrix(L_rHb_klb, 1,1);   //LINETAG: Lw.rHb
     close_r_matrix();   //LINETAG: Lw.rHb

     open_r_matrix("Lw.rGe");   //LINETAG: Lw.rGe
         wrt_r_matrix(L_rGe_klb, 1,1);   //LINETAG: Lw.rGe
     close_r_matrix();   //LINETAG: Lw.rGe

     open_r_matrix("Lw.total");   //LINETAG: Lw.total
         wrt_r_matrix(L_total_klb, 1,1);   //LINETAG: Lw.total
     close_r_matrix();   //LINETAG: Lw.total


     open_r_matrix("Ln.cHl");   //LINETAG: Ln.cHl
         wrt_r_matrix(L_cHl_num, 1,1);   //LINETAG: Ln.cHl
     close_r_matrix();   //LINETAG: Ln.cHl


	 open_r_matrix("Ln.cTw");   //LINETAG: Ln.cTw
         wrt_r_matrix(L_cTw_num, 1,1);   //LINETAG: Ln.cTw
     close_r_matrix();   //LINETAG: Ln.cTw
	 
     open_r_matrix("Ln.rHb");   //LINETAG: Ln.rHb
         wrt_r_matrix(L_rHb_num, 1,1);   //LINETAG: Ln.rHb
     close_r_matrix();   //LINETAG: Ln.rHb

     open_r_matrix("Ln.rGe");   //LINETAG: Ln.rGe
         wrt_r_matrix(L_rGe_num, 1,1);   //LINETAG: Ln.rGe
     close_r_matrix();   //LINETAG: Ln.rGe

     open_r_matrix("Ln.total");   //LINETAG: Ln.total
         wrt_r_matrix(L_total_num, 1,1);   //LINETAG: Ln.total
     close_r_matrix();   //LINETAG: Ln.total


     open_r_matrix("Dw.cHl");   //LINETAG: Dw.cHl
         wrt_r_matrix(D_cHl_klb, 1,1);   //LINETAG: Dw.cHl
     close_r_matrix();   //LINETAG: Dw.cHl

     open_r_matrix("Dw.rHb");   //LINETAG: Dw.rHb
         wrt_r_matrix(D_rHb_klb, 1,1);   //LINETAG: Dw.rHb
     close_r_matrix();   //LINETAG: Dw.rHb

     open_r_matrix("Dw.rGe");   //LINETAG: Dw.rGe
         wrt_r_matrix(D_rGe_klb, 1,1);   //LINETAG: Dw.rGe
     close_r_matrix();   //LINETAG: Dw.rGe

     open_r_matrix("Dw.total");   //LINETAG: Dw.total
         wrt_r_matrix(D_total_klb, 1,1);   //LINETAG: Dw.total
     close_r_matrix();   //LINETAG: Dw.total


     open_r_matrix("Dn.cHl");   //LINETAG: Dn.cHl
         wrt_r_matrix(D_cHl_num, 1,1);   //LINETAG: Dn.cHl
     close_r_matrix();   //LINETAG: Dn.cHl

     open_r_matrix("Dn.rHb");   //LINETAG: Dn.rHb
         wrt_r_matrix(D_rHb_num, 1,1);   //LINETAG: Dn.rHb
     close_r_matrix();   //LINETAG: Dn.rHb

     open_r_matrix("Dn.rGe");   //LINETAG: Dn.rGe
         wrt_r_matrix(D_rGe_num, 1,1);   //LINETAG: Dn.rGe
     close_r_matrix();   //LINETAG: Dn.rGe

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
