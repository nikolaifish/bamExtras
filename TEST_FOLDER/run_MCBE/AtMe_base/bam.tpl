//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##
//##  SEDAR 69  Atlantic menhaden assessment model, 2019
//##  NMFS, Beaufort Lab, Sustainable Fisheries Branch
//##
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
DATA_SECTION

!!cout << "Starting Beaufort Assessment Model" << endl;
!!cout << endl;
!!cout << "                BAM!" << endl;
!!cout << endl;

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: set-up section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

// Starting and ending year of the model (year data starts)
init_int styr;
init_int endyr;

//Starting year to estimate recruitment deviation from S-R curve
init_int styr_rec_dev;
//Ending year to estimate recruitment deviation from S-R curve
init_int endyr_rec_dev;
//possible 3 phases of constraints on recruitment deviations
init_int endyr_rec_phase1;
init_int endyr_rec_phase2;

// ending years for selectivity blocks
init_int endyr_selex_phase1;
init_int endyr_selex_phase2;
init_int endyr_selex_phase3;
init_int endyr_selex_phase4;
init_int endyr_selex_phase5;
init_int endyr_selex_phase6;
init_int endyr_selex_phase7;

//number assessment years
number nyrs;
number nyrs_rec;

//this section MUST BE INDENTED!!!
 LOCAL_CALCS
   nyrs=endyr-styr+1.;
   nyrs_rec=endyr_rec_dev-styr_rec_dev+1.;
 END_CALCS
 
//Total number of ages in population model
init_int nages;
// Vector of ages for age bins in population model
init_vector agebins(1,nages);
 
//Total number of ages used to match age comps: plus group may differ from popn, first age must not
init_int nages_agec;

//Vector of ages for age bins in age comps
init_vector agebins_agec(1,nages_agec);

//Total number of length bins for each matrix and width of bins)
init_int nlenbins;          //used to match data
init_number lenbins_width;  //width of length bins (mm)

//Vector of lengths for length bins (mm)(midpoint) 
init_vector lenbins(1,nlenbins);
 
//Max F used in spr and msy calcs
init_number max_F_spr_msy;
//Total number of iterations for spr calcs
init_int n_iter_spr;
//Total number of iterations for msy calcs
int n_iter_msy;
 LOCAL_CALCS
		n_iter_msy=n_iter_spr; 
 END_CALCS

//Starting year to compute arithmetic average recruitment for SPR-related values
init_int styr_rec_spr;
//Ending year to compute arithmetic average recruitment for SPR-related values
init_int endyr_rec_spr;
number nyrs_rec_spr;
 LOCAL_CALCS
   nyrs_rec_spr=endyr_rec_spr-styr_rec_spr+1.;
 END_CALCS 
 
//Number years at end of time series over which to average sector F's, for weighted selectivities
init_int selpar_n_yrs_wgted;
//bias correction (set to 1.0 for no bias correction or a negative value to compute from rec variance)
init_number set_BiasCor;

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: observed data section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

//###Northern adult index (NAD)#######################################################
//CPUE
init_int styr_cpue_nad;
init_int endyr_cpue_nad;
init_vector obs_cpue_nad(styr_cpue_nad,endyr_cpue_nad);   //Observed CPUE
init_vector obs_cv_cpue_nad(styr_cpue_nad,endyr_cpue_nad);    //CV of cpue

// Length Compositions 
init_int nyr_lenc_nad;
init_ivector yrs_lenc_nad(1,nyr_lenc_nad);
init_vector nsamp_lenc_nad(1,nyr_lenc_nad);
init_vector nfish_lenc_nad(1,nyr_lenc_nad);
init_matrix obs_lenc_nad(1,nyr_lenc_nad,1,nlenbins);

//###Middle adult index (MAD)#######################################################
//CPUE
init_int styr_cpue_mad;
init_int endyr_cpue_mad;
init_vector obs_cpue_mad(styr_cpue_mad,endyr_cpue_mad);   //Observed CPUE
init_vector obs_cv_cpue_mad(styr_cpue_mad,endyr_cpue_mad);    //CV of cpue

// Length Compositions 
init_int nyr_lenc_mad;
init_ivector yrs_lenc_mad(1,nyr_lenc_mad);
init_vector nsamp_lenc_mad(1,nyr_lenc_mad);
init_vector nfish_lenc_mad(1,nyr_lenc_mad);
init_matrix obs_lenc_mad(1,nyr_lenc_mad,1,nlenbins);

//###Southern adult index (SAD)#######################################################
//CPUE
init_int styr_cpue_sad;
init_int endyr_cpue_sad;
init_vector obs_cpue_sad(styr_cpue_sad,endyr_cpue_sad);   //Observed CPUE
init_vector obs_cv_cpue_sad(styr_cpue_sad,endyr_cpue_sad);    //CV of cpue

// Length Compositions 
//init_int nyr_sad_lenc;
//init_ivector yrs_sad_lenc(1,nyr_sad_lenc);
//init_vector nsamp_sad_lenc(1,nyr_sad_lenc);
//init_vector nfish_sad_lenc(1,nyr_sad_lenc);
//init_matrix obs_sad_lenc(1,nyr_sad_lenc,1,nlenbins);

//###Young of the year (YOY)################################################################
//CPUE
init_int styr_cpue_jai;
init_int endyr_cpue_jai;
init_vector obs_cpue_jai(styr_cpue_jai,endyr_cpue_jai);   //Observed CPUE
init_vector obs_cv_cpue_jai(styr_cpue_jai,endyr_cpue_jai);    //CV of cpue
init_int yr_q_change;

//###MARMAP AND ECOMON  INDEX################################################################
//CPUE
init_int nyr_cpue_mareco;
init_ivector yrs_cpue_mareco(1,nyr_cpue_mareco);
init_vector obs_cpue_mareco(1,nyr_cpue_mareco);   //Observed CPUE
init_vector obs_cv_cpue_mareco(1,nyr_cpue_mareco);    //CV of cpue

//################Commercial reduction fishery fleet - north #######################################
// Comm reduction  landings (1000 mt)
init_int styr_L_cRn;
init_int endyr_L_cRn;
init_vector obs_L_cRn(styr_L_cRn,endyr_L_cRn);
init_vector obs_cv_L_cRn(styr_L_cRn,endyr_L_cRn);

// Comm reduction age compositions
init_int nyr_agec_cRn;
init_ivector yrs_agec_cRn(1,nyr_agec_cRn);
init_vector nsamp_agec_cRn(1,nyr_agec_cRn);
init_vector nfish_agec_cRn(1,nyr_agec_cRn);
init_matrix obs_agec_cRn(1,nyr_agec_cRn,1,nages_agec);

//################Commercial reduction fishery fleet - south #######################################
// Comm reduction  landings (1000 mt)
init_int styr_L_cRs;
init_int endyr_L_cRs;
init_vector obs_L_cRs(styr_L_cRs,endyr_L_cRs);
init_vector obs_cv_L_cRs(styr_L_cRs,endyr_L_cRs);

// Comm reduction age compositions
init_int nyr_agec_cRs;
init_ivector yrs_agec_cRs(1,nyr_agec_cRs);
init_vector nsamp_agec_cRs(1,nyr_agec_cRs);
init_vector nfish_agec_cRs(1,nyr_agec_cRs);
init_matrix obs_agec_cRs(1,nyr_agec_cRs,1,nages_agec);

//################Commercial bait fishery fleet - north #######################################
// Comm bait  landings (1000 mt) - includes rec
init_int styr_L_cBn;
init_int endyr_L_cBn;
init_vector obs_L_cBn(styr_L_cBn,endyr_L_cBn);
init_vector obs_cv_L_cBn(styr_L_cBn,endyr_L_cBn);

// Comm reduction age compositions
init_int nyr_agec_cBn;
init_ivector yrs_agec_cBn(1,nyr_agec_cBn);
init_vector nsamp_agec_cBn(1,nyr_agec_cBn);
init_vector nfish_agec_cBn(1,nyr_agec_cBn);
init_matrix obs_agec_cBn(1,nyr_agec_cBn,1,nages_agec);

//################Commercial bait fishery fleet - south #######################################
// Comm bait  landings (1000 mt) - includes rec
init_int styr_L_cBs;
init_int endyr_L_cBs;
init_vector obs_L_cBs(styr_L_cBs,endyr_L_cBs);
init_vector obs_cv_L_cBs(styr_L_cBs,endyr_L_cBs);

// Comm reduction age compositions
init_int nyr_agec_cBs;
init_ivector yrs_agec_cBs(1,nyr_agec_cBs);
init_vector nsamp_agec_cBs(1,nyr_agec_cBs);
init_vector nfish_agec_cBs(1,nyr_agec_cBs);
init_matrix obs_agec_cBs(1,nyr_agec_cBs,1,nages_agec);


//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: parameter section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##################Single Parameter values and initial guesses #################################
// von Bert parms in FL mm population
init_vector set_Linf(1,7);
init_vector set_K(1,7);
init_vector set_t0(1,7);
init_vector set_len_cv_nad(1,7);
init_vector set_len_cv_mad(1,7);
//init_vector set_len_cv_sad(1,7);

//Spawner-recruit parameters (Initial guesses or fixed values)
init_vector set_steep(1,7);         //recruitment steepness
init_vector set_log_R0(1,7);        //recruitment R0
init_vector set_R_autocorr(1,7);    //recruitment autocorrelation
init_vector set_rec_sigma(1,7);     //recruitment standard deviation in log space

//Dirichlet-multinomial overdispersion parameters
init_vector set_log_dm_lenc_nad(1,7);    //Dirichlet-multinomial overdispersion parameter-NAD
init_vector set_log_dm_lenc_mad(1,7);    //Dirichlet-multinomial overdispersion parameter-MAD
//init_vector set_log_dm_sad_lc(1,7);    //Dirichlet-multinomial overdispersion parameter-SAD

init_vector set_log_dm_agec_cRn(1,7);    //Dirichlet-multinomial overdispersion parameter-Comm reduction north
init_vector set_log_dm_agec_cRs(1,7);    //Dirichlet-multinomial overdispersion parameter-Comm reduction south
init_vector set_log_dm_agec_cBn(1,7);    //Dirichlet-multinomial overdispersion parameter-Comm bait north
init_vector set_log_dm_agec_cBs(1,7);    //Dirichlet-multinomial overdispersion parameter-Comm bait south

//Initial guesses or fixed values of estimated selectivity parameters
init_vector set_selpar_A50_cRn(1,7);     //paramaters for cGN reduction double logistic
init_vector set_selpar_slope_cRn(1,7);   //north - block 1
init_vector set_selpar_A502_cRn(1,7);
init_vector set_selpar_slope2_cRn(1,7);

init_vector set_selpar_A50_cRn2(1,7);     //paramaters for cGN reduction double logistic
init_vector set_selpar_slope_cRn2(1,7);   //north - block 2
init_vector set_selpar_A502_cRn2(1,7);
init_vector set_selpar_slope2_cRn2(1,7);

init_vector set_selpar_A50_cRn3(1,7);     //paramaters for cGN reduction double logistic
init_vector set_selpar_slope_cRn3(1,7);   //north - block 3
init_vector set_selpar_A502_cRn3(1,7);
init_vector set_selpar_slope2_cRn3(1,7);

init_vector set_selpar_A50_cRn4(1,7);     //paramaters for cGN reduction double logistic
init_vector set_selpar_slope_cRn4(1,7);   //north - block 4
init_vector set_selpar_A502_cRn4(1,7);
init_vector set_selpar_slope2_cRn4(1,7);

init_vector set_selpar_A50_cRs(1,7);     //paramaters for cGN reduction double logistic
init_vector set_selpar_slope_cRs(1,7);   //south - block 1
init_vector set_selpar_A502_cRs(1,7);
init_vector set_selpar_slope2_cRs(1,7);

init_vector set_selpar_A50_cRs2(1,7);     //paramaters for cGN reduction double logistic
init_vector set_selpar_slope_cRs2(1,7);   //south - block 2
init_vector set_selpar_A502_cRs2(1,7);
init_vector set_selpar_slope2_cRs2(1,7);

init_vector set_selpar_A50_cRs3(1,7);     //paramaters for cGN reduction double logistic
init_vector set_selpar_slope_cRs3(1,7);   //south - block 3
init_vector set_selpar_A502_cRs3(1,7);
init_vector set_selpar_slope2_cRs3(1,7);

init_vector set_selpar_A50_cRs4(1,7);     //paramaters for cGN reduction double logistic
init_vector set_selpar_slope_cRs4(1,7);   //south - block 4
init_vector set_selpar_A502_cRs4(1,7);
init_vector set_selpar_slope2_cRs4(1,7);

init_vector set_selpar_A50_cBn(1,7);     //paramaters for cGN bait double logistic
init_vector set_selpar_slope_cBn(1,7);   //north - block 1
init_vector set_selpar_A502_cBn(1,7);
init_vector set_selpar_slope2_cBn(1,7);

init_vector set_selpar_A50_cBn2(1,7);     //paramaters for cGN bait double logistic
init_vector set_selpar_slope_cBn2(1,7);   //north - block 2
init_vector set_selpar_A502_cBn2(1,7);
init_vector set_selpar_slope2_cBn2(1,7);

init_vector set_selpar_A50_cBs(1,7);     //paramaters for cGN bait double logistic
init_vector set_selpar_slope_cBs(1,7);   //south - block 1
init_vector set_selpar_A502_cBs(1,7);
init_vector set_selpar_slope2_cBs(1,7);

init_vector set_selpar_A50_cBs2(1,7);     //paramaters for cGN bait double logistic
init_vector set_selpar_slope_cBs2(1,7);   //south - block 2
init_vector set_selpar_A502_cBs2(1,7);
init_vector set_selpar_slope2_cBs2(1,7);

init_vector set_sel_age0_cRn(1,7); //input in logit space by age
init_vector set_sel_age1_cRn(1,7); //Commerical reduction fishery - north - block 1
init_vector set_sel_age2_cRn(1,7);
init_vector set_sel_age3_cRn(1,7);
init_vector set_sel_age4_cRn(1,7);
init_vector set_sel_age5_cRn(1,7);
init_vector set_sel_age6_cRn(1,7);

init_vector set_sel_age0_cRn2(1,7); //input in logit space by age
init_vector set_sel_age1_cRn2(1,7); //Commerical reduction fishery - north - block 2
init_vector set_sel_age2_cRn2(1,7);
init_vector set_sel_age3_cRn2(1,7);
init_vector set_sel_age4_cRn2(1,7);
init_vector set_sel_age5_cRn2(1,7);
init_vector set_sel_age6_cRn2(1,7);

init_vector set_sel_age0_cRn3(1,7); //input in logit space by age
init_vector set_sel_age1_cRn3(1,7); //Commerical reduction fishery - north - block 3
init_vector set_sel_age2_cRn3(1,7);
init_vector set_sel_age3_cRn3(1,7);
init_vector set_sel_age4_cRn3(1,7);
init_vector set_sel_age5_cRn3(1,7);
init_vector set_sel_age6_cRn3(1,7);

init_vector set_sel_age0_cRn4(1,7); //input in logit space by age
init_vector set_sel_age1_cRn4(1,7); //Commerical reduction fishery - north - block 4
init_vector set_sel_age2_cRn4(1,7);
init_vector set_sel_age3_cRn4(1,7);
init_vector set_sel_age4_cRn4(1,7);
init_vector set_sel_age5_cRn4(1,7);
init_vector set_sel_age6_cRn4(1,7);

init_vector set_sel_age0_cRs(1,7); //input in logit space by age
init_vector set_sel_age1_cRs(1,7); //Commerical reduction fishery - south - block 1
init_vector set_sel_age2_cRs(1,7);
init_vector set_sel_age3_cRs(1,7);
init_vector set_sel_age4_cRs(1,7);
init_vector set_sel_age5_cRs(1,7);
init_vector set_sel_age6_cRs(1,7);

init_vector set_sel_age0_cRs2(1,7); //input in logit space by age
init_vector set_sel_age1_cRs2(1,7); //Commerical reduction fishery - south - block 2
init_vector set_sel_age2_cRs2(1,7);
init_vector set_sel_age3_cRs2(1,7);
init_vector set_sel_age4_cRs2(1,7);
init_vector set_sel_age5_cRs2(1,7);
init_vector set_sel_age6_cRs2(1,7);

init_vector set_sel_age0_cRs3(1,7); //input in logit space by age
init_vector set_sel_age1_cRs3(1,7); //Commerical reduction fishery - south - block 3
init_vector set_sel_age2_cRs3(1,7);
init_vector set_sel_age3_cRs3(1,7);
init_vector set_sel_age4_cRs3(1,7);
init_vector set_sel_age5_cRs3(1,7);
init_vector set_sel_age6_cRs3(1,7);

init_vector set_sel_age0_cRs4(1,7); //input in logit space by age
init_vector set_sel_age1_cRs4(1,7); //Commerical reduction fishery - south - block 4
init_vector set_sel_age2_cRs4(1,7);
init_vector set_sel_age3_cRs4(1,7);
init_vector set_sel_age4_cRs4(1,7);
init_vector set_sel_age5_cRs4(1,7);
init_vector set_sel_age6_cRs4(1,7);

init_vector set_sel_age0_cBn(1,7); //input in logit space by age
init_vector set_sel_age1_cBn(1,7); //Commerical bait fishery - north - block 1
init_vector set_sel_age2_cBn(1,7);
init_vector set_sel_age3_cBn(1,7);
init_vector set_sel_age4_cBn(1,7);
init_vector set_sel_age5_cBn(1,7);
init_vector set_sel_age6_cBn(1,7);

init_vector set_sel_age0_cBn2(1,7); //input in logit space by age
init_vector set_sel_age1_cBn2(1,7); //Commerical bait fishery - north - block 2
init_vector set_sel_age2_cBn2(1,7);
init_vector set_sel_age3_cBn2(1,7);
init_vector set_sel_age4_cBn2(1,7);
init_vector set_sel_age5_cBn2(1,7);
init_vector set_sel_age6_cBn2(1,7);

init_vector set_sel_age0_cBs(1,7); //input in logit space by age
init_vector set_sel_age1_cBs(1,7); //Commerical bait fishery - south - block 1
init_vector set_sel_age2_cBs(1,7);
init_vector set_sel_age3_cBs(1,7);
init_vector set_sel_age4_cBs(1,7);
init_vector set_sel_age5_cBs(1,7);
init_vector set_sel_age6_cBs(1,7);

init_vector set_sel_age0_cBs2(1,7); //input in logit space by age
init_vector set_sel_age1_cBs2(1,7); //Commerical bait fishery - south - block 2
init_vector set_sel_age2_cBs2(1,7);
init_vector set_sel_age3_cBs2(1,7);
init_vector set_sel_age4_cBs2(1,7);
init_vector set_sel_age5_cBs2(1,7);
init_vector set_sel_age6_cBs2(1,7);


init_vector set_selpar_A50_nad(1,7);    //parameters for northern adult index (NAD)
init_vector set_selpar_slope_nad(1,7);  //double logistic or logistic
init_vector set_selpar_A502_nad(1,7);
init_vector set_selpar_slope2_nad(1,7);

init_vector set_selpar_A50_mad(1,7);    //parameters for middle adult index (MAD)
init_vector set_selpar_slope_mad(1,7);  //double logistic or logistic
init_vector set_selpar_A502_mad(1,7);
init_vector set_selpar_slope2_mad(1,7);

init_vector set_selpar_A50_sad(1,7);    //parameters for southern adult index (SAD)
init_vector set_selpar_slope_sad(1,7);  //double logistic or logistic
init_vector set_selpar_A502_sad(1,7);
init_vector set_selpar_slope2_sad(1,7);

init_vector set_sel_age0_nad(1,7); //input in logit space by age
init_vector set_sel_age1_nad(1,7); //NAD index
init_vector set_sel_age2_nad(1,7);
init_vector set_sel_age3_nad(1,7);
init_vector set_sel_age4_nad(1,7);
init_vector set_sel_age5_nad(1,7);
init_vector set_sel_age6_nad(1,7);

init_vector set_sel_age0_mad(1,7); //input in logit space by age
init_vector set_sel_age1_mad(1,7); //MAD index
init_vector set_sel_age2_mad(1,7);
init_vector set_sel_age3_mad(1,7);
init_vector set_sel_age4_mad(1,7);
init_vector set_sel_age5_mad(1,7);
init_vector set_sel_age6_mad(1,7);

init_vector set_sel_age0_sad(1,7); //input in logit space by age
init_vector set_sel_age1_sad(1,7); //SAD index
init_vector set_sel_age2_sad(1,7);
init_vector set_sel_age3_sad(1,7);
init_vector set_sel_age4_sad(1,7);
init_vector set_sel_age5_sad(1,7);
init_vector set_sel_age6_sad(1,7);

//--index catchability-----------------------------------------------------------------------------------
init_vector set_log_q_cpue_nad(1,7);    //catchability coefficient (log) for NAD index
init_vector set_log_q_cpue_mad(1,7);    //catchability coefficient (log) for MAD index
init_vector set_log_q_cpue_sad(1,7);    //catchability coefficient (log) for SAD index
init_vector set_log_q_cpue_jai(1,7);    //catchability coefficient (log) for JAI index
init_vector set_log_q2_jai(1,7);    //catchability coefficient (log) for JAI index
init_vector set_log_q_cpue_mar(1,7);    //catchability coefficient (log) for JAI index
init_vector set_log_q_cpue_eco(1,7);    //catchability coefficient (log) for JAI index

//--mean F's in log space --------------------------------
init_vector set_log_avg_F_L_cRn(1,7);    //commercial reduction fishery - north
init_vector set_log_avg_F_L_cRs(1,7);    //commercial reduction fishery - south
init_vector set_log_avg_F_L_cBn(1,7);    //commercial bait fishery - north
init_vector set_log_avg_F_L_cBs(1,7);    //commercial bait fishery - south

//##################Dev Vector Parameter values (vals) and bounds #################################
//--F vectors---------------------------
init_vector set_log_dev_F_L_cRn(1,3);   //commerical reduction F devs - north
init_vector set_log_dev_F_L_cRs(1,3);   //commerical reduction F devs - south
init_vector set_log_dev_F_L_cBn(1,3);   //commerical bait F devs - north
init_vector set_log_dev_F_L_cBs(1,3);   //commerical bait F devs - south

init_vector set_log_dev_RWq(1,3);    //rand walk on q devs
init_vector set_log_dev_rec(1,3);    //recruitment devs
init_vector set_log_dev_Nage(1,3);   //initial age structure devs

init_vector set_log_dev_vals_F_L__cRn(styr_L_cRn,endyr_L_cRn);
init_vector set_log_dev_vals_F_L__cRs(styr_L_cRs,endyr_L_cRs);
init_vector set_log_dev_vals_F_L__cBn(styr_L_cBn,endyr_L_cBn);
init_vector set_log_dev_vals_F_L__cBs(styr_L_cBs,endyr_L_cBs);

init_vector set_log_dev_vals_rec(styr_rec_dev,endyr_rec_dev);
init_vector set_log_dev_vals_Nage(2,nages);             


//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: likelihood weights section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
init_number set_w_L;            //weight for landings

init_number set_w_cpue_nad;        //weight for NAD index
init_number set_w_cpue_mad;        //weight for MAD index
init_number set_w_cpue_sad;        //weight for SAD index
init_number set_w_cpue_jai;        //weight for JAI index
init_number set_w_cpue_mareco;     //weight for MARMAP and ECOMON index

init_number set_w_lenc_nad;      //weight for NAD len comps
init_number set_w_lenc_mad;      //weight for MAD len comps
//init_number set_w_lc_sad;      //weight for SAD len comps

init_number set_w_agec_cRn;        //weight for cGN reduction age comps - north
init_number set_w_agec_cRs;        //weight for cGN reduction age comps - south
init_number set_w_agec_cBn;        //weight for cGN bait age comps - north
init_number set_w_agec_cBs;        //weight for cGN bait age comps - south

init_number set_w_Nage_init;    //for fitting initial abundance at age (excluding first age)
init_number set_w_rec;          //for fitting S-R curve
init_number set_w_rec_early;    //additional constraint on early years recruitment
init_number set_w_rec_end;      //additional constraint on ending years recruitment 
init_number set_w_fullF;        //penalty for any Fapex>3(removed in final phase of optimization)
init_number set_w_Ftune;        //weight applied to tuning F (removed in final phase of optimization)


//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: miscellaneous stuff section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//FL(mm)-weight(whole weight in g) relationship: W=aL^b
init_number wgtpar_a;
init_number wgtpar_b;

///Maturity and proportion female at age
init_vector obs_maturity_f(1,nages);            //proportion females mature at age
init_vector obs_maturity_m(1,nages);            //proportion males mature at age
init_matrix maturity_obs_tv(styr,endyr,1,nages);    //Proportion of females mature at age over time
init_vector obs_prop_f(1,nages);                //sex ratio; prop female at age
init_vector fecundity(1,nages);                 //fecundity at age
init_matrix fecundity_tv(styr,endyr,1,nages);   //time-varying fecundity at age
init_vector wgt_spawn(1,nages);                 //weight at age at spawning
init_matrix wgt_spawn_tv(styr,endyr,1,nages);   //time-varying weight at age at spawning
init_vector wgt_middle(1,nages);                //weight at age at middle of year
init_matrix wgt_middle_tv(styr,endyr,1,nages);  //time-varying weight at age at the middle of fishing year
init_matrix len_apr15_tv(styr,endyr,1,nages);   //time-varying length at age on May 15
init_matrix len_jun1_tv(styr,endyr,1,nages);    //time-varying length at age on June 1
init_matrix len_oct15_tv(styr,endyr,1,nages);   //time-varying length at age on Oct 15
init_number spawn_time_frac; //time of year of peak spawning, as a fraction of the year

// Natural mortality
init_vector set_M(1,nages);     //age-dependent: used in model
init_matrix set_M_tv(styr,endyr,1,nages);  //time-varying and age-varying  natural mortality

//Spawner-recruit parameters (Initial guesses or fixed values)
init_int SR_switch;

//rate of increase on q
init_int set_q_rate_phase;  //value sets estimation phase of rate increase, negative value turns it off
init_number set_q_rate;

//density dependence on fishery q's 
init_int set_q_DD_phase;      //value sets estimation phase of random walk, negative value turns it off
init_number set_q_DD_beta;    //value of 0.0 is density indepenent
init_number set_q_DD_beta_se;
init_int set_q_DD_stage;      //age to begin counting biomass, should be near full exploitation

//random walk on fishery q's 
init_number set_RWq_var;     //assumed variance of RW q

//Tune Fapex (tuning removed in final year of optimization)
init_number set_Ftune;
init_int set_Ftune_yr;

//threshold sample sizes for length comps 
init_number minSS_lenc_nad;
init_number minSS_lenc_mad;
//init_number minSS_sad_lenc;

//threshold sample sizes for age comps
init_number minSS_agec_cRn;
init_number minSS_agec_cRs;
init_number minSS_agec_cBn;
init_number minSS_agec_cBs;

//input for deterministic F-based projections
init_int endyr_proj; //last year of projections, by default, first year is endyr+1
init_int styr_regs;  //apply current F until styr_regs, then the projection F
init_int Fproj_switch; //Value to use for Fproj: 1=current F, 2=Fmsy, 3=F30, 4=F40
init_number Fproj_mult; //multiplier to the Fproj
int styr_proj;
 LOCAL_CALCS
   styr_proj=endyr+1;
 END_CALCS
 
//ageing error matrix (columns are true ages, rows are ages as read for age comps: columns should sum to one)
init_matrix age_error(1,nages,1,nages);

// #######Indexing integers for year(iyear), age(iage),length(ilen) ###############
int iyear;
int iage;
int ilen;
int ff;

number sqrt2pi;
number g2mt;                    //conversion of grams to metric tons 
number g2kg;                    //conversion of grams to kg   
number g2klb;                   //conversion of grams to 1000 lb   
number mt2klb;                  //conversion of metric tons to 1000 lb
number mt2lb;                   //conversion of metric tons to lb
number dzero;                   //small additive constant to prevent division by zero
number huge_number;             //huge number, to avoid irregular parameter space

init_number end_of_data_file;

//!!cout << "start year of reduction landings " << styr_L_cRn << endl;
//!!cout << "Start year of JAI " << styr_cpue_jai << endl;
//!!cout << "Linf set up " << set_Linf << endl;
//!!cout << "start year of cBs " << styr_L_cBs << endl;
//!!cout << "start year of cRs " << styr_L_cRs << endl;
//!!cout << "start year of cBn " << styr_L_cBn << endl;

//this section MUST BE INDENTED!!!
 LOCAL_CALCS
   if(end_of_data_file!=999)
   {
       cout << "*** WARNING: Data File NOT READ CORRECTLY ****" << endl;
       exit(0);  
   }
   else
   {cout << "Data File read correctly" << endl;} 
 END_CALCS   


PARAMETER_SECTION

 LOCAL_CALCS
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

 END_CALCS
 
////--------------Growth--------------------------------------------------------------------------- 
  //Population growth parms and conversions
  init_bounded_number Linf(Linf_LO,Linf_HI,Linf_PH);
  init_bounded_number K(K_LO,K_HI,K_PH);
  init_bounded_number t0(t0_LO,t0_HI,t0_PH);
  init_bounded_number len_cv_nad_val(len_cv_nad_LO,len_cv_nad_HI,len_cv_nad_PH);  
  init_bounded_number len_cv_mad_val(len_cv_mad_LO,len_cv_mad_HI,len_cv_mad_PH);
  //init_bounded_number len_cv_sad_val(len_cv_sad_LO,len_cv_sad_HI,len_cv_sad_PH);
  vector Linf_out(1,8);
  vector K_out(1,8);
  vector t0_out(1,8);
  vector len_cv_nad_val_out(1,8);
  vector len_cv_mad_val_out(1,8);
  //vector len_cv_sad_val_out(1,8);

  vector meanlen_FL(1,nages);   //mean fork length (mm) at age all fish
  matrix meanlen_FL_apr15(styr,endyr,1,nages);  //mean fork length April 15
  matrix meanlen_FL_jun1(styr,endyr,1,nages);   //mean fork length June 1
  matrix meanlen_FL_oct15(styr,endyr,1,nages);  //mean fork length October 15
  vector wgt_fish_mt(1,nages);
  vector wgt_spawn_mt(1,nages);
  matrix wgt_spawn_mt_tv(styr,endyr,1,nages);   //time varying weight at spawn
  matrix tv_wgt_middle_mt(styr,endyr,1,nages);
  matrix tv_wgt_spawn_mt(styr,endyr,1,nages);

  matrix len_cR_mm(styr,endyr,1,nages);          //mean length at age of commercial reduction landings in mm 
  matrix wholewgt_cR_mt(styr,endyr,1,nages);    //whole wgt of commercial reduction landings in 1000 mt   
   
  3darray lenprob_apr15(styr,endyr,1,nages,1,nlenbins);  //distn of size at age (age-length key, 10 mm bins) in population
  3darray lenprob_jun1(styr,endyr,1,nages,1,nlenbins);   //distn of size at age (age-length key, 10 mm bins) in population
  3darray lenprob_oct15(styr,endyr,1,nages,1,nlenbins);  //distn of size at age (age-length key, 10 mm bins) in population
  matrix lenprob(1,nages,1,nlenbins);           //distn of size at age (age-length key, 10 mm bins) in population
  number zscore_len;                            //standardized normal values used for computing lenprob
  vector cprob_lenvec(1,nlenbins);              //cumulative probabilities used for computing lenprob
  number zscore_lzero;                          //standardized normal values for length = 0
  number cprob_lzero;                           //length probability mass below zero, used for computing lenprob
    
  //matrices below are used to match length comps
  3darray lenprob_nad(styr,endyr,1,nages,1,nlenbins);    
  3darray lenprob_mad(styr,endyr,1,nages,1,nlenbins);    
  //3darray lenprob_sad(styr,endyr,1,nages,1,nlenbins);   
  
  //init_bounded_dev_vector log_len_cv_dev(1,nages,-2,2,3)
  matrix len_sd_nad(styr,endyr,1,nages);
  matrix len_sd_mad(styr,endyr,1,nages);
  //matrix len_sd_sad(styr,endyr,1,nages);
  vector len_cv_nad(1,nages); //for fishgraph
  vector len_cv_mad(1,nages); //for fishgraph
  //vector len_cv_sad(1,nages); //for fishgraph
  vector len_cv_apr15(1,nages);
  vector len_cv_jun1(1,nages);
  vector len_cv_oct15(1,nages);

//----Predicted length and age compositions
  matrix pred_nad_lenc(1,nyr_lenc_nad,1,nlenbins);
  matrix pred_mad_lenc(1,nyr_lenc_mad,1,nlenbins);
  //matrix pred_sad_lenc(1,nyr_sad_lenc,1,nlenbins);

  matrix pred_cRn_agec(1,nyr_agec_cRn,1,nages_agec);
  matrix pred_cRn_agec_allages(1,nyr_agec_cRn,1,nages);
  matrix ErrorFree_cRn_agec(1,nyr_agec_cRn,1,nages);    

  matrix pred_cRs_agec(1,nyr_agec_cRs,1,nages_agec);
  matrix pred_cRs_agec_allages(1,nyr_agec_cRs,1,nages);
  matrix ErrorFree_cRs_agec(1,nyr_agec_cRs,1,nages);    

  matrix pred_cBn_agec(1,nyr_agec_cBn,1,nages_agec);
  matrix pred_cBn_agec_allages(1,nyr_agec_cBn,1,nages);
  matrix ErrorFree_cBn_agec(1,nyr_agec_cBn,1,nages);
  
  matrix pred_cBs_agec(1,nyr_agec_cBs,1,nages_agec);
  matrix pred_cBs_agec_allages(1,nyr_agec_cBs,1,nages);
  matrix ErrorFree_cBs_agec(1,nyr_agec_cBs,1,nages);    
  
//Sample size (perhaps adjusted herein) used in fitting comp data
  vector nsamp_nad_lenc_allyr(styr,endyr);
  vector nsamp_mad_lenc_allyr(styr,endyr);
  //vector nsamp_sad_lenc_allyr(styr,endyr);

  vector nsamp_cRn_agec_allyr(styr,endyr);
  vector nsamp_cRs_agec_allyr(styr,endyr);
  vector nsamp_cBn_agec_allyr(styr,endyr);
  vector nsamp_cBs_agec_allyr(styr,endyr);

//Nfish used in MCB analysis (not used in fitting)
  vector nfish_nad_lenc_allyr(styr,endyr);
  vector nfish_mad_lenc_allyr(styr,endyr);
  //vector nfish_sad_lenc_allyr(styr,endyr);

  vector nfish_cRs_agec_allyr(styr,endyr);
  vector nfish_cRn_agec_allyr(styr,endyr);
  vector nfish_cBn_agec_allyr(styr,endyr);
  vector nfish_cBs_agec_allyr(styr,endyr);
  
//Computed effective sample size for output (not used in fitting)
  vector neff_nad_lenc_allyr(styr,endyr);
  vector neff_mad_lenc_allyr(styr,endyr);
  //vector neff_sad_lenc_allyr(styr,endyr);

  vector neff_cRn_agec_allyr(styr,endyr);
  vector neff_cRs_agec_allyr(styr,endyr);
  vector neff_cBn_agec_allyr(styr,endyr);
  vector neff_cBs_agec_allyr(styr,endyr);
  
//-----Population-----------------------------------------------------------------------------------
  matrix N(styr,endyr+1,1,nages);           //Population numbers by year and age at start of yr
  matrix N_mdyr(styr,endyr,1,nages);        //Population numbers by year and age at mdpt of yr: used for comps and cpue
  matrix N_spawn(styr,endyr,1,nages);       //Population numbers by year and age at peaking spawning: used for SSB
  matrix N_nad(styr,endyr,1,nages);
  matrix N_mad(styr,endyr,1,nages);
  matrix N_sad(styr,endyr,1,nages);
  init_bounded_vector log_dev_Nage(2,nages,log_Nage_dev_LO,log_Nage_dev_HI,log_Nage_dev_PH);
  vector log_Nage_dev_output(1,nages);      //used in output. equals zero for first age
  matrix B(styr,endyr+1,1,nages);           //Population biomass by year and age at start of yr
  matrix B_mdyr(styr,endyr+1,1,nages);      //population biomass at the middle of the yr
  vector totB(styr,endyr+1);                //Total biomass by year
  vector totN(styr,endyr+1);                //Total abundance by year
  vector SSB(styr,endyr+1);                   //Total spawning biomass by year (female + male mature biomass) 
  vector rec(styr,endyr+1);                 //Recruits by year
  vector prop_f(1,nages);
  vector prop_m(1,nages);
  vector maturity_f(1,nages);
  vector maturity_m(1,nages);
  matrix tv_maturity(styr,endyr,1,nages);
  vector reprod(1,nages);
  matrix reprod_tv(styr,endyr,1,nages);
  matrix SSBatage(styr,endyr,1,nages);
 
//---Stock-Recruit Function (Beverton-Holt, steepness parameterization)----------
  init_bounded_number log_R0(log_R0_LO,log_R0_HI,log_R0_PH);        //log(virgin Recruitment)
  vector log_R0_out(1,8);
  number R0;                                  //virgin recruitment
  init_bounded_number steep(steep_LO,steep_HI,steep_PH); //steepness
  vector steep_out(1,8);
  init_bounded_number rec_sigma(rec_sigma_LO,rec_sigma_HI,rec_sigma_PH);  //sd recruitment residuals  
  vector rec_sigma_out(1,8);
  init_bounded_number R_autocorr(R_autocorr_LO,R_autocorr_HI,R_autocorr_PH);  //autocorrelation in SR  
  vector R_autocorr_out(1,8);

  number rec_sigma_sq;                        //square of rec_sigma      
  number rec_logL_add;                        //additive term in -logL term   
 
  init_bounded_dev_vector log_dev_rec(styr_rec_dev,endyr_rec_dev,log_rec_dev_LO,log_rec_dev_HI,log_rec_dev_PH);
  vector log_rec_dev_output(styr,endyr+1);             //used in t.series output. equals zero except for yrs in log_dev_rec
  vector log_rec_dev_out(styr_rec_dev,endyr_rec_dev);  //used in output for bound checking
  
  number var_rec_dev;                           //variance of log recruitment deviations, from yrs with unconstrainted S-R(XXXX-XXXX)
  number sigma_rec_dev;                         //sample SD of log residuals (may not equal rec_sigma 
  number BiasCor;                               //Bias correction in equilibrium recruits
  number S0;                                    //equal to spr_F0*R0 = virgin SSB
  number B0;                                    //equal to bpr_F0*R0 = virgin B  
  number R1;                                    //Recruits in styr
  number R_virgin;                              //unfished recruitment with bias correction
  vector SdS0(styr,endyr+1);                      //Spawners relative to the unfished level              
  
  init_bounded_number log_dm_lenc_nad(log_dm_nad_lc_LO,log_dm_nad_lc_HI,log_dm_nad_lc_PH);
  init_bounded_number log_dm_lenc_mad(log_dm_mad_lc_LO,log_dm_mad_lc_HI,log_dm_mad_lc_PH);
  //init_bounded_number log_dm_sad_lc(log_dm_sad_lc_LO,log_dm_sad_lc_HI,log_dm_sad_lc_PH);
  
  init_bounded_number log_dm_agec_cRn(log_dm_cRn_ac_LO,log_dm_cRn_ac_HI,log_dm_cRn_ac_PH);
  init_bounded_number log_dm_agec_cRs(log_dm_cRs_ac_LO,log_dm_cRs_ac_HI,log_dm_cRs_ac_PH);
  init_bounded_number log_dm_agec_cBn(log_dm_cBn_ac_LO,log_dm_cBn_ac_HI,log_dm_cBn_ac_PH);
  init_bounded_number log_dm_agec_cBs(log_dm_cBs_ac_LO,log_dm_cBs_ac_HI,log_dm_cBs_ac_PH);
  
  vector log_dm_nad_lc_out(1,8);
  vector log_dm_mad_lc_out(1,8);
  //vector log_dm_sad_lc_out(1,8);

  vector log_dm_cRn_ac_out(1,8);
  vector log_dm_cRs_ac_out(1,8);
  vector log_dm_cBn_ac_out(1,8);
  vector log_dm_cBs_ac_out(1,8);
  
//-----------------------------------------------------------------------------------------------------------------------------------------------
////---Selectivity-------------------------------------------------------------------------

//Commercial reduction------------------------------------------------------------------
  matrix sel_cRn(styr,endyr,1,nages);  
  matrix sel_cRs(styr,endyr,1,nages);
  vector sel_cRn_block1(1,nages);
  vector sel_cRn_block2(1,nages);
  vector sel_cRn_block3(1,nages);
  vector sel_cRn_block4(1,nages);
  vector sel_cRs_block1(1,nages);
  vector sel_cRs_block2(1,nages);
  vector sel_cRs_block3(1,nages);
  vector sel_cRs_block4(1,nages);

  //north - block 1
  init_bounded_number selpar_A50_cRn(selpar_A50_cRn_LO,selpar_A50_cRn_HI,selpar_A50_cRn_PH);
  init_bounded_number selpar_slope_cRn(selpar_slope_cRn_LO,selpar_slope_cRn_HI,selpar_slope_cRn_PH);
  init_bounded_number selpar_A502_cRn(selpar_A502_cRn_LO,selpar_A502_cRn_HI,selpar_A502_cRn_PH);
  init_bounded_number selpar_slope2_cRn(selpar_slope2_cRn_LO,selpar_slope2_cRn_HI,selpar_slope2_cRn_PH);
  vector selpar_A50_cRn_out(1,8);
  vector selpar_slope_cRn_out(1,8);
  vector selpar_A502_cRn_out(1,8);
  vector selpar_slope2_cRn_out(1,8);

  //north - block 2
  init_bounded_number selpar_A50_cRn2(selpar_A50_cRn2_LO,selpar_A50_cRn2_HI,selpar_A50_cRn2_PH);
  init_bounded_number selpar_slope_cRn2(selpar_slope_cRn2_LO,selpar_slope_cRn2_HI,selpar_slope_cRn2_PH);
  init_bounded_number selpar_A502_cRn2(selpar_A502_cRn2_LO,selpar_A502_cRn2_HI,selpar_A502_cRn2_PH);
  init_bounded_number selpar_slope2_cRn2(selpar_slope2_cRn2_LO,selpar_slope2_cRn2_HI,selpar_slope2_cRn2_PH);
  vector selpar_A50_cRn2_out(1,8);
  vector selpar_slope_cRn2_out(1,8);
  vector selpar_A502_cRn2_out(1,8);
  vector selpar_slope2_cRn2_out(1,8);

  //north - block 3
  init_bounded_number selpar_A50_cRn3(selpar_A50_cRn3_LO,selpar_A50_cRn3_HI,selpar_A50_cRn3_PH);
  init_bounded_number selpar_slope_cRn3(selpar_slope_cRn3_LO,selpar_slope_cRn3_HI,selpar_slope_cRn3_PH);
  init_bounded_number selpar_A502_cRn3(selpar_A502_cRn3_LO,selpar_A502_cRn3_HI,selpar_A502_cRn3_PH);
  init_bounded_number selpar_slope2_cRn3(selpar_slope2_cRn3_LO,selpar_slope2_cRn3_HI,selpar_slope2_cRn3_PH);
  vector selpar_A50_cRn3_out(1,8);
  vector selpar_slope_cRn3_out(1,8);
  vector selpar_A502_cRn3_out(1,8);
  vector selpar_slope2_cRn3_out(1,8);

  //north - block 4
  init_bounded_number selpar_A50_cRn4(selpar_A50_cRn4_LO,selpar_A50_cRn4_HI,selpar_A50_cRn4_PH);
  init_bounded_number selpar_slope_cRn4(selpar_slope_cRn4_LO,selpar_slope_cRn4_HI,selpar_slope_cRn4_PH);
  init_bounded_number selpar_A502_cRn4(selpar_A502_cRn4_LO,selpar_A502_cRn4_HI,selpar_A502_cRn4_PH);
  init_bounded_number selpar_slope2_cRn4(selpar_slope2_cRn4_LO,selpar_slope2_cRn4_HI,selpar_slope2_cRn4_PH);
  vector selpar_A50_cRn4_out(1,8);
  vector selpar_slope_cRn4_out(1,8);
  vector selpar_A502_cRn4_out(1,8);
  vector selpar_slope2_cRn4_out(1,8);

  //south - block 1
  init_bounded_number selpar_A50_cRs(selpar_A50_cRs_LO,selpar_A50_cRs_HI,selpar_A50_cRs_PH);
  init_bounded_number selpar_slope_cRs(selpar_slope_cRs_LO,selpar_slope_cRs_HI,selpar_slope_cRs_PH);
  init_bounded_number selpar_A502_cRs(selpar_A502_cRs_LO,selpar_A502_cRs_HI,selpar_A502_cRs_PH);
  init_bounded_number selpar_slope2_cRs(selpar_slope2_cRs_LO,selpar_slope2_cRs_HI,selpar_slope2_cRs_PH);
  vector selpar_A50_cRs_out(1,8);
  vector selpar_slope_cRs_out(1,8);
  vector selpar_A502_cRs_out(1,8);
  vector selpar_slope2_cRs_out(1,8);

  //south - block 2
  init_bounded_number selpar_A50_cRs2(selpar_A50_cRs2_LO,selpar_A50_cRs2_HI,selpar_A50_cRs2_PH);
  init_bounded_number selpar_slope_cRs2(selpar_slope_cRs2_LO,selpar_slope_cRs2_HI,selpar_slope_cRs2_PH);
  init_bounded_number selpar_A502_cRs2(selpar_A502_cRs2_LO,selpar_A502_cRs2_HI,selpar_A502_cRs2_PH);
  init_bounded_number selpar_slope2_cRs2(selpar_slope2_cRs2_LO,selpar_slope2_cRs2_HI,selpar_slope2_cRs2_PH);
  vector selpar_A50_cRs2_out(1,8);
  vector selpar_slope_cRs2_out(1,8);
  vector selpar_A502_cRs2_out(1,8);
  vector selpar_slope2_cRs2_out(1,8);

  //south - block 3
  init_bounded_number selpar_A50_cRs3(selpar_A50_cRs3_LO,selpar_A50_cRs3_HI,selpar_A50_cRs3_PH);
  init_bounded_number selpar_slope_cRs3(selpar_slope_cRs3_LO,selpar_slope_cRs3_HI,selpar_slope_cRs3_PH);
  init_bounded_number selpar_A502_cRs3(selpar_A502_cRs3_LO,selpar_A502_cRs3_HI,selpar_A502_cRs3_PH);
  init_bounded_number selpar_slope2_cRs3(selpar_slope2_cRs3_LO,selpar_slope2_cRs3_HI,selpar_slope2_cRs3_PH);
  vector selpar_A50_cRs3_out(1,8);
  vector selpar_slope_cRs3_out(1,8);
  vector selpar_A502_cRs3_out(1,8);
  vector selpar_slope2_cRs3_out(1,8);

  //south - block 4
  init_bounded_number selpar_A50_cRs4(selpar_A50_cRs4_LO,selpar_A50_cRs4_HI,selpar_A50_cRs4_PH);
  init_bounded_number selpar_slope_cRs4(selpar_slope_cRs4_LO,selpar_slope_cRs4_HI,selpar_slope_cRs4_PH);
  init_bounded_number selpar_A502_cRs4(selpar_A502_cRs4_LO,selpar_A502_cRs4_HI,selpar_A502_cRs4_PH);
  init_bounded_number selpar_slope2_cRs4(selpar_slope2_cRs4_LO,selpar_slope2_cRs4_HI,selpar_slope2_cRs4_PH);
  vector selpar_A50_cRs4_out(1,8);
  vector selpar_slope_cRs4_out(1,8);
  vector selpar_A502_cRs4_out(1,8);
  vector selpar_slope2_cRs4_out(1,8);

  //logit based selectivity cR north - block 1
  init_bounded_number sel_age0_cRn_logit(selpar_age0_cRn_LO,selpar_age0_cRn_HI,selpar_age0_cRn_PH);  
  init_bounded_number sel_age1_cRn_logit(selpar_age1_cRn_LO,selpar_age1_cRn_HI,selpar_age1_cRn_PH);
  init_bounded_number sel_age2_cRn_logit(selpar_age2_cRn_LO,selpar_age2_cRn_HI,selpar_age2_cRn_PH);
  init_bounded_number sel_age3_cRn_logit(selpar_age3_cRn_LO,selpar_age3_cRn_HI,selpar_age3_cRn_PH);
  init_bounded_number sel_age4_cRn_logit(selpar_age4_cRn_LO,selpar_age4_cRn_HI,selpar_age4_cRn_PH);
  init_bounded_number sel_age5_cRn_logit(selpar_age5_cRn_LO,selpar_age5_cRn_HI,selpar_age5_cRn_PH);
  init_bounded_number sel_age6_cRn_logit(selpar_age6_cRn_LO,selpar_age6_cRn_HI,selpar_age6_cRn_PH);
  vector sel_age_cRn_vec(1,nages);
  number selpar_age0_cRn;
  number selpar_age1_cRn;
  number selpar_age2_cRn;
  number selpar_age3_cRn;
  number selpar_age4_cRn;
  number selpar_age5_cRn;
  number selpar_age6_cRn;
  vector selpar_age0_cRn_out(1,8);
  vector selpar_age1_cRn_out(1,8);
  vector selpar_age2_cRn_out(1,8);
  vector selpar_age3_cRn_out(1,8);
  vector selpar_age4_cRn_out(1,8);
  vector selpar_age5_cRn_out(1,8);
  vector selpar_age6_cRn_out(1,8);

  //logit based selectivity cR north - block 2
  init_bounded_number sel_age0_cRn2_logit(selpar_age0_cRn2_LO,selpar_age0_cRn2_HI,selpar_age0_cRn2_PH);  
  init_bounded_number sel_age1_cRn2_logit(selpar_age1_cRn2_LO,selpar_age1_cRn2_HI,selpar_age1_cRn2_PH);
  init_bounded_number sel_age2_cRn2_logit(selpar_age2_cRn2_LO,selpar_age2_cRn2_HI,selpar_age2_cRn2_PH);
  init_bounded_number sel_age3_cRn2_logit(selpar_age3_cRn2_LO,selpar_age3_cRn2_HI,selpar_age3_cRn2_PH);
  init_bounded_number sel_age4_cRn2_logit(selpar_age4_cRn2_LO,selpar_age4_cRn2_HI,selpar_age4_cRn2_PH);
  init_bounded_number sel_age5_cRn2_logit(selpar_age5_cRn2_LO,selpar_age5_cRn2_HI,selpar_age5_cRn2_PH);
  init_bounded_number sel_age6_cRn2_logit(selpar_age6_cRn2_LO,selpar_age6_cRn2_HI,selpar_age6_cRn2_PH);
  vector sel_age_cRn2_vec(1,nages);
  number selpar_age0_cRn2;
  number selpar_age1_cRn2;
  number selpar_age2_cRn2;
  number selpar_age3_cRn2;
  number selpar_age4_cRn2;
  number selpar_age5_cRn2;
  number selpar_age6_cRn2;
  vector selpar_age0_cRn2_out(1,8);
  vector selpar_age1_cRn2_out(1,8);
  vector selpar_age2_cRn2_out(1,8);
  vector selpar_age3_cRn2_out(1,8);
  vector selpar_age4_cRn2_out(1,8);
  vector selpar_age5_cRn2_out(1,8);
  vector selpar_age6_cRn2_out(1,8);

  //logit based selectivity cR north - block 3
  init_bounded_number sel_age0_cRn3_logit(selpar_age0_cRn3_LO,selpar_age0_cRn3_HI,selpar_age0_cRn3_PH);  
  init_bounded_number sel_age1_cRn3_logit(selpar_age1_cRn3_LO,selpar_age1_cRn3_HI,selpar_age1_cRn3_PH);
  init_bounded_number sel_age2_cRn3_logit(selpar_age2_cRn3_LO,selpar_age2_cRn3_HI,selpar_age2_cRn3_PH);
  init_bounded_number sel_age3_cRn3_logit(selpar_age3_cRn3_LO,selpar_age3_cRn3_HI,selpar_age3_cRn3_PH);
  init_bounded_number sel_age4_cRn3_logit(selpar_age4_cRn3_LO,selpar_age4_cRn3_HI,selpar_age4_cRn3_PH);
  init_bounded_number sel_age5_cRn3_logit(selpar_age5_cRn3_LO,selpar_age5_cRn3_HI,selpar_age5_cRn3_PH);
  init_bounded_number sel_age6_cRn3_logit(selpar_age6_cRn3_LO,selpar_age6_cRn3_HI,selpar_age6_cRn3_PH);
  vector sel_age_cRn3_vec(1,nages);
  number selpar_age0_cRn3;
  number selpar_age1_cRn3;
  number selpar_age2_cRn3;
  number selpar_age3_cRn3;
  number selpar_age4_cRn3;
  number selpar_age5_cRn3;
  number selpar_age6_cRn3;
  vector selpar_age0_cRn3_out(1,8);
  vector selpar_age1_cRn3_out(1,8);
  vector selpar_age2_cRn3_out(1,8);
  vector selpar_age3_cRn3_out(1,8);
  vector selpar_age4_cRn3_out(1,8);
  vector selpar_age5_cRn3_out(1,8);
  vector selpar_age6_cRn3_out(1,8);

  //logit based selectivity cR north - block 4
  init_bounded_number sel_age0_cRn4_logit(selpar_age0_cRn4_LO,selpar_age0_cRn4_HI,selpar_age0_cRn4_PH);  
  init_bounded_number sel_age1_cRn4_logit(selpar_age1_cRn4_LO,selpar_age1_cRn4_HI,selpar_age1_cRn4_PH);
  init_bounded_number sel_age2_cRn4_logit(selpar_age2_cRn4_LO,selpar_age2_cRn4_HI,selpar_age2_cRn4_PH);
  init_bounded_number sel_age3_cRn4_logit(selpar_age3_cRn4_LO,selpar_age3_cRn4_HI,selpar_age3_cRn4_PH);
  init_bounded_number sel_age4_cRn4_logit(selpar_age4_cRn4_LO,selpar_age4_cRn4_HI,selpar_age4_cRn4_PH);
  init_bounded_number sel_age5_cRn4_logit(selpar_age5_cRn4_LO,selpar_age5_cRn4_HI,selpar_age5_cRn4_PH);
  init_bounded_number sel_age6_cRn4_logit(selpar_age6_cRn4_LO,selpar_age6_cRn4_HI,selpar_age6_cRn4_PH);
  vector sel_age_cRn4_vec(1,nages);
  number selpar_age0_cRn4;
  number selpar_age1_cRn4;
  number selpar_age2_cRn4;
  number selpar_age3_cRn4;
  number selpar_age4_cRn4;
  number selpar_age5_cRn4;
  number selpar_age6_cRn4;
  vector selpar_age0_cRn4_out(1,8);
  vector selpar_age1_cRn4_out(1,8);
  vector selpar_age2_cRn4_out(1,8);
  vector selpar_age3_cRn4_out(1,8);
  vector selpar_age4_cRn4_out(1,8);
  vector selpar_age5_cRn4_out(1,8);
  vector selpar_age6_cRn4_out(1,8);
    
  //logit based selectivity cR south - block 1
  init_bounded_number sel_age0_cRs_logit(selpar_age0_cRs_LO,selpar_age0_cRs_HI,selpar_age0_cRs_PH);  
  init_bounded_number sel_age1_cRs_logit(selpar_age1_cRs_LO,selpar_age1_cRs_HI,selpar_age1_cRs_PH);
  init_bounded_number sel_age2_cRs_logit(selpar_age2_cRs_LO,selpar_age2_cRs_HI,selpar_age2_cRs_PH);
  init_bounded_number sel_age3_cRs_logit(selpar_age3_cRs_LO,selpar_age3_cRs_HI,selpar_age3_cRs_PH);
  init_bounded_number sel_age4_cRs_logit(selpar_age4_cRs_LO,selpar_age4_cRs_HI,selpar_age4_cRs_PH);
  init_bounded_number sel_age5_cRs_logit(selpar_age5_cRs_LO,selpar_age5_cRs_HI,selpar_age5_cRs_PH);
  init_bounded_number sel_age6_cRs_logit(selpar_age6_cRs_LO,selpar_age6_cRs_HI,selpar_age6_cRs_PH);
  vector sel_age_cRs_vec(1,nages);
  number selpar_age0_cRs;
  number selpar_age1_cRs;
  number selpar_age2_cRs;
  number selpar_age3_cRs;
  number selpar_age4_cRs;
  number selpar_age5_cRs;
  number selpar_age6_cRs;
  vector selpar_age0_cRs_out(1,8);
  vector selpar_age1_cRs_out(1,8);
  vector selpar_age2_cRs_out(1,8);
  vector selpar_age3_cRs_out(1,8);
  vector selpar_age4_cRs_out(1,8);
  vector selpar_age5_cRs_out(1,8);
  vector selpar_age6_cRs_out(1,8);

  //logit based selectivity cR south - block 2
  init_bounded_number sel_age0_cRs2_logit(selpar_age0_cRs2_LO,selpar_age0_cRs2_HI,selpar_age0_cRs2_PH);  
  init_bounded_number sel_age1_cRs2_logit(selpar_age1_cRs2_LO,selpar_age1_cRs2_HI,selpar_age1_cRs2_PH);
  init_bounded_number sel_age2_cRs2_logit(selpar_age2_cRs2_LO,selpar_age2_cRs2_HI,selpar_age2_cRs2_PH);
  init_bounded_number sel_age3_cRs2_logit(selpar_age3_cRs2_LO,selpar_age3_cRs2_HI,selpar_age3_cRs2_PH);
  init_bounded_number sel_age4_cRs2_logit(selpar_age4_cRs2_LO,selpar_age4_cRs2_HI,selpar_age4_cRs2_PH);
  init_bounded_number sel_age5_cRs2_logit(selpar_age5_cRs2_LO,selpar_age5_cRs2_HI,selpar_age5_cRs2_PH);
  init_bounded_number sel_age6_cRs2_logit(selpar_age6_cRs2_LO,selpar_age6_cRs2_HI,selpar_age6_cRs2_PH);
  vector sel_age_cRs2_vec(1,nages);
  number selpar_age0_cRs2;
  number selpar_age1_cRs2;
  number selpar_age2_cRs2;
  number selpar_age3_cRs2;
  number selpar_age4_cRs2;
  number selpar_age5_cRs2;
  number selpar_age6_cRs2;
  vector selpar_age0_cRs2_out(1,8);
  vector selpar_age1_cRs2_out(1,8);
  vector selpar_age2_cRs2_out(1,8);
  vector selpar_age3_cRs2_out(1,8);
  vector selpar_age4_cRs2_out(1,8);
  vector selpar_age5_cRs2_out(1,8);
  vector selpar_age6_cRs2_out(1,8);

  //logit based selectivity cR south - block 3
  init_bounded_number sel_age0_cRs3_logit(selpar_age0_cRs3_LO,selpar_age0_cRs3_HI,selpar_age0_cRs3_PH);  
  init_bounded_number sel_age1_cRs3_logit(selpar_age1_cRs3_LO,selpar_age1_cRs3_HI,selpar_age1_cRs3_PH);
  init_bounded_number sel_age2_cRs3_logit(selpar_age2_cRs3_LO,selpar_age2_cRs3_HI,selpar_age2_cRs3_PH);
  init_bounded_number sel_age3_cRs3_logit(selpar_age3_cRs3_LO,selpar_age3_cRs3_HI,selpar_age3_cRs3_PH);
  init_bounded_number sel_age4_cRs3_logit(selpar_age4_cRs3_LO,selpar_age4_cRs3_HI,selpar_age4_cRs3_PH);
  init_bounded_number sel_age5_cRs3_logit(selpar_age5_cRs3_LO,selpar_age5_cRs3_HI,selpar_age5_cRs3_PH);
  init_bounded_number sel_age6_cRs3_logit(selpar_age6_cRs3_LO,selpar_age6_cRs3_HI,selpar_age6_cRs3_PH);
  vector sel_age_cRs3_vec(1,nages);
  number selpar_age0_cRs3;
  number selpar_age1_cRs3;
  number selpar_age2_cRs3;
  number selpar_age3_cRs3;
  number selpar_age4_cRs3;
  number selpar_age5_cRs3;
  number selpar_age6_cRs3;
  vector selpar_age0_cRs3_out(1,8);
  vector selpar_age1_cRs3_out(1,8);
  vector selpar_age2_cRs3_out(1,8);
  vector selpar_age3_cRs3_out(1,8);
  vector selpar_age4_cRs3_out(1,8);
  vector selpar_age5_cRs3_out(1,8);
  vector selpar_age6_cRs3_out(1,8);

  //logit based selectivity cR south - block 4
  init_bounded_number sel_age0_cRs4_logit(selpar_age0_cRs4_LO,selpar_age0_cRs4_HI,selpar_age0_cRs4_PH);  
  init_bounded_number sel_age1_cRs4_logit(selpar_age1_cRs4_LO,selpar_age1_cRs4_HI,selpar_age1_cRs4_PH);
  init_bounded_number sel_age2_cRs4_logit(selpar_age2_cRs4_LO,selpar_age2_cRs4_HI,selpar_age2_cRs4_PH);
  init_bounded_number sel_age3_cRs4_logit(selpar_age3_cRs4_LO,selpar_age3_cRs4_HI,selpar_age3_cRs4_PH);
  init_bounded_number sel_age4_cRs4_logit(selpar_age4_cRs4_LO,selpar_age4_cRs4_HI,selpar_age4_cRs4_PH);
  init_bounded_number sel_age5_cRs4_logit(selpar_age5_cRs4_LO,selpar_age5_cRs4_HI,selpar_age5_cRs4_PH);
  init_bounded_number sel_age6_cRs4_logit(selpar_age6_cRs4_LO,selpar_age6_cRs4_HI,selpar_age6_cRs4_PH);
  vector sel_age_cRs4_vec(1,nages);
  number selpar_age0_cRs4;
  number selpar_age1_cRs4;
  number selpar_age2_cRs4;
  number selpar_age3_cRs4;
  number selpar_age4_cRs4;
  number selpar_age5_cRs4;
  number selpar_age6_cRs4;
  vector selpar_age0_cRs4_out(1,8);
  vector selpar_age1_cRs4_out(1,8);
  vector selpar_age2_cRs4_out(1,8);
  vector selpar_age3_cRs4_out(1,8);
  vector selpar_age4_cRs4_out(1,8);
  vector selpar_age5_cRs4_out(1,8);
  vector selpar_age6_cRs4_out(1,8);

//Commercial bait--------------------------------------------------------------------
  matrix sel_cBn(styr,endyr,1,nages);  
  matrix sel_cBs(styr,endyr,1,nages);
  vector sel_cBn_block1(1,nages);
  vector sel_cBn_block2(1,nages);
  vector sel_cBs_block1(1,nages);
  vector sel_cBs_block2(1,nages);

  //north - block 1
  init_bounded_number selpar_A50_cBn(selpar_A50_cBn_LO,selpar_A50_cBn_HI,selpar_A50_cBn_PH);
  init_bounded_number selpar_slope_cBn(selpar_slope_cBn_LO,selpar_slope_cBn_HI,selpar_slope_cBn_PH);
  init_bounded_number selpar_A502_cBn(selpar_A502_cBn_LO,selpar_A502_cBn_HI,selpar_A502_cBn_PH);
  init_bounded_number selpar_slope2_cBn(selpar_slope2_cBn_LO,selpar_slope2_cBn_HI,selpar_slope2_cBn_PH);
  vector selpar_A50_cBn_out(1,8);
  vector selpar_slope_cBn_out(1,8);
  vector selpar_A502_cBn_out(1,8);
  vector selpar_slope2_cBn_out(1,8);

  //north - block 2
  init_bounded_number selpar_A50_cBn2(selpar_A50_cBn2_LO,selpar_A50_cBn2_HI,selpar_A50_cBn2_PH);
  init_bounded_number selpar_slope_cBn2(selpar_slope_cBn2_LO,selpar_slope_cBn2_HI,selpar_slope_cBn2_PH);
  init_bounded_number selpar_A502_cBn2(selpar_A502_cBn2_LO,selpar_A502_cBn2_HI,selpar_A502_cBn2_PH);
  init_bounded_number selpar_slope2_cBn2(selpar_slope2_cBn2_LO,selpar_slope2_cBn2_HI,selpar_slope2_cBn2_PH);
  vector selpar_A50_cBn2_out(1,8);
  vector selpar_slope_cBn2_out(1,8);
  vector selpar_A502_cBn2_out(1,8);
  vector selpar_slope2_cBn2_out(1,8);

  //south - block 1
  init_bounded_number selpar_A50_cBs(selpar_A50_cBs_LO,selpar_A50_cBs_HI,selpar_A50_cBs_PH);
  init_bounded_number selpar_slope_cBs(selpar_slope_cBs_LO,selpar_slope_cBs_HI,selpar_slope_cBs_PH);
  init_bounded_number selpar_A502_cBs(selpar_A502_cBs_LO,selpar_A502_cBs_HI,selpar_A502_cBs_PH);
  init_bounded_number selpar_slope2_cBs(selpar_slope2_cBs_LO,selpar_slope2_cBs_HI,selpar_slope2_cBs_PH);
  vector selpar_A50_cBs_out(1,8);
  vector selpar_slope_cBs_out(1,8);
  vector selpar_A502_cBs_out(1,8);
  vector selpar_slope2_cBs_out(1,8);

  //south - block 2
  init_bounded_number selpar_A50_cBs2(selpar_A50_cBs2_LO,selpar_A50_cBs2_HI,selpar_A50_cBs2_PH);
  init_bounded_number selpar_slope_cBs2(selpar_slope_cBs2_LO,selpar_slope_cBs2_HI,selpar_slope_cBs2_PH);
  init_bounded_number selpar_A502_cBs2(selpar_A502_cBs2_LO,selpar_A502_cBs2_HI,selpar_A502_cBs2_PH);
  init_bounded_number selpar_slope2_cBs2(selpar_slope2_cBs2_LO,selpar_slope2_cBs2_HI,selpar_slope2_cBs2_PH);
  vector selpar_A50_cBs2_out(1,8);
  vector selpar_slope_cBs2_out(1,8);
  vector selpar_A502_cBs2_out(1,8);
  vector selpar_slope2_cBs2_out(1,8);

  //logit based selectivity cB north - block 1
  init_bounded_number sel_age0_cBn_logit(selpar_age0_cBn_LO,selpar_age0_cBn_HI,selpar_age0_cBn_PH);  
  init_bounded_number sel_age1_cBn_logit(selpar_age1_cBn_LO,selpar_age1_cBn_HI,selpar_age1_cBn_PH);
  init_bounded_number sel_age2_cBn_logit(selpar_age2_cBn_LO,selpar_age2_cBn_HI,selpar_age2_cBn_PH);
  init_bounded_number sel_age3_cBn_logit(selpar_age3_cBn_LO,selpar_age3_cBn_HI,selpar_age3_cBn_PH);
  init_bounded_number sel_age4_cBn_logit(selpar_age4_cBn_LO,selpar_age4_cBn_HI,selpar_age4_cBn_PH);
  init_bounded_number sel_age5_cBn_logit(selpar_age5_cBn_LO,selpar_age5_cBn_HI,selpar_age5_cBn_PH);
  init_bounded_number sel_age6_cBn_logit(selpar_age6_cBn_LO,selpar_age6_cBn_HI,selpar_age6_cBn_PH);
  vector sel_age_cBn_vec(1,nages);
  number selpar_age0_cBn;
  number selpar_age1_cBn;
  number selpar_age2_cBn;
  number selpar_age3_cBn;
  number selpar_age4_cBn;
  number selpar_age5_cBn;
  number selpar_age6_cBn;
  vector selpar_age0_cBn_out(1,8);
  vector selpar_age1_cBn_out(1,8);
  vector selpar_age2_cBn_out(1,8);
  vector selpar_age3_cBn_out(1,8);
  vector selpar_age4_cBn_out(1,8);
  vector selpar_age5_cBn_out(1,8);
  vector selpar_age6_cBn_out(1,8);

  //logit based selectivity cB north - block 2
  init_bounded_number sel_age0_cBn2_logit(selpar_age0_cBn2_LO,selpar_age0_cBn2_HI,selpar_age0_cBn2_PH);  
  init_bounded_number sel_age1_cBn2_logit(selpar_age1_cBn2_LO,selpar_age1_cBn2_HI,selpar_age1_cBn2_PH);
  init_bounded_number sel_age2_cBn2_logit(selpar_age2_cBn2_LO,selpar_age2_cBn2_HI,selpar_age2_cBn2_PH);
  init_bounded_number sel_age3_cBn2_logit(selpar_age3_cBn2_LO,selpar_age3_cBn2_HI,selpar_age3_cBn2_PH);
  init_bounded_number sel_age4_cBn2_logit(selpar_age4_cBn2_LO,selpar_age4_cBn2_HI,selpar_age4_cBn2_PH);
  init_bounded_number sel_age5_cBn2_logit(selpar_age5_cBn2_LO,selpar_age5_cBn2_HI,selpar_age5_cBn2_PH);
  init_bounded_number sel_age6_cBn2_logit(selpar_age6_cBn2_LO,selpar_age6_cBn2_HI,selpar_age6_cBn2_PH);
  vector sel_age_cBn2_vec(1,nages);
  number selpar_age0_cBn2;
  number selpar_age1_cBn2;
  number selpar_age2_cBn2;
  number selpar_age3_cBn2;
  number selpar_age4_cBn2;
  number selpar_age5_cBn2;
  number selpar_age6_cBn2;
  vector selpar_age0_cBn2_out(1,8);
  vector selpar_age1_cBn2_out(1,8);
  vector selpar_age2_cBn2_out(1,8);
  vector selpar_age3_cBn2_out(1,8);
  vector selpar_age4_cBn2_out(1,8);
  vector selpar_age5_cBn2_out(1,8);
  vector selpar_age6_cBn2_out(1,8);

  //logit based selectivity cB south - block 1
  init_bounded_number sel_age0_cBs_logit(selpar_age0_cBs_LO,selpar_age0_cBs_HI,selpar_age0_cBs_PH);  
  init_bounded_number sel_age1_cBs_logit(selpar_age1_cBs_LO,selpar_age1_cBs_HI,selpar_age1_cBs_PH);
  init_bounded_number sel_age2_cBs_logit(selpar_age2_cBs_LO,selpar_age2_cBs_HI,selpar_age2_cBs_PH);
  init_bounded_number sel_age3_cBs_logit(selpar_age3_cBs_LO,selpar_age3_cBs_HI,selpar_age3_cBs_PH);
  init_bounded_number sel_age4_cBs_logit(selpar_age4_cBs_LO,selpar_age4_cBs_HI,selpar_age4_cBs_PH);
  init_bounded_number sel_age5_cBs_logit(selpar_age5_cBs_LO,selpar_age5_cBs_HI,selpar_age5_cBs_PH);
  init_bounded_number sel_age6_cBs_logit(selpar_age6_cBs_LO,selpar_age6_cBs_HI,selpar_age6_cBs_PH);
  vector sel_age_cBs_vec(1,nages);
  number selpar_age0_cBs;
  number selpar_age1_cBs;
  number selpar_age2_cBs;
  number selpar_age3_cBs;
  number selpar_age4_cBs;
  number selpar_age5_cBs;
  number selpar_age6_cBs;
  vector selpar_age0_cBs_out(1,8);
  vector selpar_age1_cBs_out(1,8);
  vector selpar_age2_cBs_out(1,8);
  vector selpar_age3_cBs_out(1,8);
  vector selpar_age4_cBs_out(1,8);
  vector selpar_age5_cBs_out(1,8);
  vector selpar_age6_cBs_out(1,8);

  //logit based selectivity cB south - block 2
  init_bounded_number sel_age0_cBs2_logit(selpar_age0_cBs2_LO,selpar_age0_cBs2_HI,selpar_age0_cBs2_PH);  
  init_bounded_number sel_age1_cBs2_logit(selpar_age1_cBs2_LO,selpar_age1_cBs2_HI,selpar_age1_cBs2_PH);
  init_bounded_number sel_age2_cBs2_logit(selpar_age2_cBs2_LO,selpar_age2_cBs2_HI,selpar_age2_cBs2_PH);
  init_bounded_number sel_age3_cBs2_logit(selpar_age3_cBs2_LO,selpar_age3_cBs2_HI,selpar_age3_cBs2_PH);
  init_bounded_number sel_age4_cBs2_logit(selpar_age4_cBs2_LO,selpar_age4_cBs2_HI,selpar_age4_cBs2_PH);
  init_bounded_number sel_age5_cBs2_logit(selpar_age5_cBs2_LO,selpar_age5_cBs2_HI,selpar_age5_cBs2_PH);
  init_bounded_number sel_age6_cBs2_logit(selpar_age6_cBs2_LO,selpar_age6_cBs2_HI,selpar_age6_cBs2_PH);
  vector sel_age_cBs2_vec(1,nages);
  number selpar_age0_cBs2;
  number selpar_age1_cBs2;
  number selpar_age2_cBs2;
  number selpar_age3_cBs2;
  number selpar_age4_cBs2;
  number selpar_age5_cBs2;
  number selpar_age6_cBs2;
  vector selpar_age0_cBs2_out(1,8);
  vector selpar_age1_cBs2_out(1,8);
  vector selpar_age2_cBs2_out(1,8);
  vector selpar_age3_cBs2_out(1,8);
  vector selpar_age4_cBs2_out(1,8);
  vector selpar_age5_cBs2_out(1,8);
  vector selpar_age6_cBs2_out(1,8);

//Northern adult index selectivity---------------------------------------------------
  matrix sel_nad(styr_cpue_nad,endyr_cpue_nad,1,nages);
  vector sel_nad_block1(1,nages);
   
  init_bounded_number selpar_A50_nad(selpar_A50_nad_LO,selpar_A50_nad_HI,selpar_A50_nad_PH);
  init_bounded_number selpar_slope_nad(selpar_slope_nad_LO,selpar_slope_nad_HI,selpar_slope_nad_PH);
  init_bounded_number selpar_A502_nad(selpar_A502_nad_LO,selpar_A502_nad_HI,selpar_A502_nad_PH);
  init_bounded_number selpar_slope2_nad(selpar_slope2_nad_LO,selpar_slope2_nad_HI,selpar_slope2_nad_PH);
  vector selpar_A50_nad_out(1,8);
  vector selpar_slope_nad_out(1,8);
  vector selpar_A502_nad_out(1,8);
  vector selpar_slope2_nad_out(1,8);

  //logit based selectivity for NAD index
  init_bounded_number sel_age0_nad_logit(selpar_age0_nad_LO,selpar_age0_nad_HI,selpar_age0_nad_PH);  
  init_bounded_number sel_age1_nad_logit(selpar_age1_nad_LO,selpar_age1_nad_HI,selpar_age1_nad_PH);
  init_bounded_number sel_age2_nad_logit(selpar_age2_nad_LO,selpar_age2_nad_HI,selpar_age2_nad_PH);
  init_bounded_number sel_age3_nad_logit(selpar_age3_nad_LO,selpar_age3_nad_HI,selpar_age3_nad_PH);
  init_bounded_number sel_age4_nad_logit(selpar_age4_nad_LO,selpar_age4_nad_HI,selpar_age4_nad_PH);
  init_bounded_number sel_age5_nad_logit(selpar_age5_nad_LO,selpar_age5_nad_HI,selpar_age5_nad_PH);
  init_bounded_number sel_age6_nad_logit(selpar_age6_nad_LO,selpar_age6_nad_HI,selpar_age6_nad_PH);
  vector sel_age_nad_vec(1,nages);
  number selpar_age0_nad;
  number selpar_age1_nad;
  number selpar_age2_nad;
  number selpar_age3_nad;
  number selpar_age4_nad;
  number selpar_age5_nad;
  number selpar_age6_nad;
  vector selpar_age0_nad_out(1,8);
  vector selpar_age1_nad_out(1,8);
  vector selpar_age2_nad_out(1,8);
  vector selpar_age3_nad_out(1,8);
  vector selpar_age4_nad_out(1,8);
  vector selpar_age5_nad_out(1,8);
  vector selpar_age6_nad_out(1,8);

//Middle adult index selectivity---------------------------------------------------
  matrix sel_mad(styr_cpue_mad,endyr_cpue_mad,1,nages);
  vector sel_mad_block1(1,nages);
   
  init_bounded_number selpar_A50_mad(selpar_A50_mad_LO,selpar_A50_mad_HI,selpar_A50_mad_PH);
  init_bounded_number selpar_slope_mad(selpar_slope_mad_LO,selpar_slope_mad_HI,selpar_slope_mad_PH);
  init_bounded_number selpar_A502_mad(selpar_A502_mad_LO,selpar_A502_mad_HI,selpar_A502_mad_PH);
  init_bounded_number selpar_slope2_mad(selpar_slope2_mad_LO,selpar_slope2_mad_HI,selpar_slope2_mad_PH);
  vector selpar_A50_mad_out(1,8);
  vector selpar_slope_mad_out(1,8);
  vector selpar_A502_mad_out(1,8);
  vector selpar_slope2_mad_out(1,8);

  //logit based selectivity for MAD index
  init_bounded_number sel_age0_mad_logit(selpar_age0_mad_LO,selpar_age0_mad_HI,selpar_age0_mad_PH);  
  init_bounded_number sel_age1_mad_logit(selpar_age1_mad_LO,selpar_age1_mad_HI,selpar_age1_mad_PH);
  init_bounded_number sel_age2_mad_logit(selpar_age2_mad_LO,selpar_age2_mad_HI,selpar_age2_mad_PH);
  init_bounded_number sel_age3_mad_logit(selpar_age3_mad_LO,selpar_age3_mad_HI,selpar_age3_mad_PH);
  init_bounded_number sel_age4_mad_logit(selpar_age4_mad_LO,selpar_age4_mad_HI,selpar_age4_mad_PH);
  init_bounded_number sel_age5_mad_logit(selpar_age5_mad_LO,selpar_age5_mad_HI,selpar_age5_mad_PH);
  init_bounded_number sel_age6_mad_logit(selpar_age6_mad_LO,selpar_age6_mad_HI,selpar_age6_mad_PH);
  vector sel_age_mad_vec(1,nages);
  number selpar_age0_mad;
  number selpar_age1_mad;
  number selpar_age2_mad;
  number selpar_age3_mad;
  number selpar_age4_mad;
  number selpar_age5_mad;
  number selpar_age6_mad;
  vector selpar_age0_mad_out(1,8);
  vector selpar_age1_mad_out(1,8);
  vector selpar_age2_mad_out(1,8);
  vector selpar_age3_mad_out(1,8);
  vector selpar_age4_mad_out(1,8);
  vector selpar_age5_mad_out(1,8);
  vector selpar_age6_mad_out(1,8);

//Southern adult index selectivity---------------------------------------------------
  matrix sel_sad(styr_cpue_sad,endyr_cpue_sad,1,nages);
  vector sel_sad_block1(1,nages);
   
  init_bounded_number selpar_A50_sad(selpar_A50_sad_LO,selpar_A50_sad_HI,selpar_A50_sad_PH);
  init_bounded_number selpar_slope_sad(selpar_slope_sad_LO,selpar_slope_sad_HI,selpar_slope_sad_PH);
  init_bounded_number selpar_A502_sad(selpar_A502_sad_LO,selpar_A502_sad_HI,selpar_A502_sad_PH);
  init_bounded_number selpar_slope2_sad(selpar_slope2_sad_LO,selpar_slope2_sad_HI,selpar_slope2_sad_PH);
  vector selpar_A50_sad_out(1,8);
  vector selpar_slope_sad_out(1,8);
  vector selpar_A502_sad_out(1,8);
  vector selpar_slope2_sad_out(1,8);

  //logit based selectivity for SAD index
  init_bounded_number sel_age0_sad_logit(selpar_age0_sad_LO,selpar_age0_sad_HI,selpar_age0_sad_PH);  
  init_bounded_number sel_age1_sad_logit(selpar_age1_sad_LO,selpar_age1_sad_HI,selpar_age1_sad_PH);
  init_bounded_number sel_age2_sad_logit(selpar_age2_sad_LO,selpar_age2_sad_HI,selpar_age2_sad_PH);
  init_bounded_number sel_age3_sad_logit(selpar_age3_sad_LO,selpar_age3_sad_HI,selpar_age3_sad_PH);
  init_bounded_number sel_age4_sad_logit(selpar_age4_sad_LO,selpar_age4_sad_HI,selpar_age4_sad_PH);
  init_bounded_number sel_age5_sad_logit(selpar_age5_sad_LO,selpar_age5_sad_HI,selpar_age5_sad_PH);
  init_bounded_number sel_age6_sad_logit(selpar_age6_sad_LO,selpar_age6_sad_HI,selpar_age6_sad_PH);
  vector sel_age_sad_vec(1,nages);
  number selpar_age0_sad;
  number selpar_age1_sad;
  number selpar_age2_sad;
  number selpar_age3_sad;
  number selpar_age4_sad;
  number selpar_age5_sad;
  number selpar_age6_sad;
  vector selpar_age0_sad_out(1,8);
  vector selpar_age1_sad_out(1,8);
  vector selpar_age2_sad_out(1,8);
  vector selpar_age3_sad_out(1,8);
  vector selpar_age4_sad_out(1,8);
  vector selpar_age5_sad_out(1,8);
  vector selpar_age6_sad_out(1,8);

//Weighted total selectivity--------------------------------------------  
 //effort-weighted, recent selectivities
  vector sel_wgted_L(1,nages);  //toward landings 
  vector sel_wgted_tot(1,nages);//toward Z, landings plus deads discards

//-----------------------------------------------------------------------------------------------------------------------------------------------
//-------CPUE Predictions--------------------------------
  //Northern adult index
  vector pred_nad_cpue(styr_cpue_nad,endyr_cpue_nad);    //predicted NAD index (number fish per effort)
  matrix N_nad_cpue(styr_cpue_nad,endyr_cpue_nad,1,nages);    //used to compute NAD index

  //Middle adult index
  vector pred_mad_cpue(styr_cpue_mad,endyr_cpue_mad);    //predicted MAD index (number fish per effort)
  matrix N_mad_cpue(styr_cpue_mad,endyr_cpue_mad,1,nages);    //used to compute MAD index

  //Southern adult index
  vector pred_sad_cpue(styr_cpue_sad,endyr_cpue_sad);    //predicted SAD index (number fish per effort)
  vector N_sad_cpue(styr_cpue_sad,endyr_cpue_sad);    //used to compute SAD index

  //Juvenile abundance index
  vector pred_jai_cpue(styr_cpue_jai,endyr_cpue_jai);    //predicted JAI index (number fish per effort)
  vector N_jai_cpue(styr_cpue_jai,endyr_cpue_jai);    //used to compute JAI index

  //MARMAP and ECOMON index
  vector pred_mareco_cpue(1,nyr_cpue_mareco);
  vector SSB_mareco_cpue(1,nyr_cpue_mareco);

  
//---Catchability (CPUE q's)----------------------------------------------------------
  init_bounded_number log_q_cpue_nad(log_q_nad_LO,log_q_nad_HI,log_q_nad_PH);
  init_bounded_number log_q_cpue_mad(log_q_mad_LO,log_q_mad_HI,log_q_mad_PH);
  init_bounded_number log_q_cpue_sad(log_q_sad_LO,log_q_sad_HI,log_q_sad_PH);
  init_bounded_number log_q_cpue_jai(log_q_jai_LO,log_q_jai_HI,log_q_jai_PH);
  init_bounded_number log_q2_jai(log_q2_jai_LO,log_q2_jai_HI,log_q2_jai_PH);
  init_bounded_number log_q_cpue_mar(log_q_mar_LO,log_q_mar_HI,log_q_mar_PH);
  init_bounded_number log_q_cpue_eco(log_q_eco_LO,log_q_eco_HI,log_q_eco_PH);
  vector log_q_nad_out(1,8);
  vector log_q_mad_out(1,8);
  vector log_q_sad_out(1,8);
  vector log_q_jai_out(1,8);
  vector log_q2_jai_out(1,8);
  vector log_q_mar_out(1,8);
  vector log_q_eco_out(1,8);

  number q_rate;
  vector q_rate_fcn_nad(styr_cpue_nad,endyr_cpue_nad);    //increase due to technology creep (saturates in 2003)-shouldn't really be used since fishery-independent 
  vector q_rate_fcn_mad(styr_cpue_mad,endyr_cpue_mad);    //increase due to technology creep (saturates in 2003)-shouldn't really be used since fishery-independent 
  vector q_rate_fcn_sad(styr_cpue_sad,endyr_cpue_sad);    //increase due to technology creep (saturates in 2003)-shouldn't really be used since fishery-independent 
  vector q_rate_fcn_jai(styr_cpue_jai,endyr_cpue_jai);    //increase due to technology creep (saturates in 2003)-shouldn't really be used since fishery-independent
  //vector q_rate_fcn_mar(styr_mar_cpue,endyr_mar_cpue);    //increase due to technology creep (saturates in 2003)-shouldn't really be used since fishery-independent
  //vector q_rate_fcn_eco(styr_eco_cpue,endyr_eco_cpue);    //increase due to technology creep (saturates in 2003)-shouldn't really be used since fishery-independent
  
  //init_bounded_number q_DD_beta(0.1,0.9,set_q_DD_phase);    //not estimated so commented out and declared as number (below)
  number q_DD_beta;
  vector q_DD_fcn(styr,endyr);    //density dependent function as a multiple of q (scaled a la Katsukawa and Matsuda. 2003)
  number B0_q_DD;                 //B0 of ages q_DD_age plus
  vector B_q_DD(styr,endyr);      //annual biomass of ages q_DD_age plus

  //Random walk catchability - not using
  init_bounded_vector q_RW_log_dev_nad(styr_cpue_nad,endyr_cpue_nad-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
  init_bounded_vector q_RW_log_dev_mad(styr_cpue_mad,endyr_cpue_mad-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
  init_bounded_vector q_RW_log_dev_sad(styr_cpue_sad,endyr_cpue_sad-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
  init_bounded_vector q_RW_log_dev_jai(styr_cpue_jai,endyr_cpue_jai-1,log_RWq_LO,log_RWq_HI,log_RWq_PH); 
  //init_bounded_vector q_RW_log_dev_mar(styr_mar_cpue,endyr_mar_cpue-1,log_RWq_LO,log_RWq_HI,log_RWq_PH); 
  //init_bounded_vector q_RW_log_dev_eco(styr_eco_cpue,endyr_eco_cpue-1,log_RWq_LO,log_RWq_HI,log_RWq_PH); 
  
  //Fishery independent catchability over time, may be constant
  vector q_nad(styr_cpue_nad,endyr_cpue_nad);
  vector q_mad(styr_cpue_mad,endyr_cpue_mad);
  vector q_sad(styr_cpue_sad,endyr_cpue_sad);
  vector q_jai(styr_cpue_jai,endyr_cpue_jai); 
  vector q2_jai(styr_cpue_jai,endyr_cpue_jai); 
  vector q_mar(1,8); 
  vector q_eco(9,26); 
  
//----------------------------------------------------------------------------------------------------------------------------------------------- 
//---Landings in numbers (total or 1000 fish) and in wgt (1000s mt)--------------------------------------------------
  matrix L_cRn_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_cRn_mt(styr,endyr,1,nages);    //landings (1000 mt whole weight) at age    
  vector pred_cRn_L_mt(styr,endyr);       //yearly landings in 1000 mt whole summed over ages

  matrix L_cRs_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_cRs_mt(styr,endyr,1,nages);    //landings (1000 mt whole weight) at age    
  vector pred_cRs_L_mt(styr,endyr);       //yearly landings in 1000 mt whole summed over ages

  matrix L_cBn_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_cBn_mt(styr,endyr,1,nages);    //landings (1000 mt whole weight) at age    
  vector pred_cBn_L_mt(styr,endyr);       //yearly landings in 1000 mt whole summed over ages

  matrix L_cBs_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_cBs_mt(styr,endyr,1,nages);    //landings (1000 mt whole weight) at age    
  vector pred_cBs_L_mt(styr,endyr);       //yearly landings in 1000 mt whole summed over ages

  matrix L_total_num(styr,endyr,1,nages); //total landings in number at age
  matrix L_total_mt(styr,endyr,1,nages);  //landings in 1000 mt whole wgt at age 
  vector L_total_mt_yr(styr,endyr);       //total landings (1000 mt whole wgt) by yr summed over ages
  

////---MSY calcs----------------------------------------------------------------------------
  number F_cRn_prop;       //proportion of F_sum attributable to cRn, last X=selpar_n_yrs_wgted yrs  
  number F_cRs_prop;       //proportion of F_sum attributable to cRs, last X=selpar_n_yrs_wgted yrs  
  number F_cBn_prop;       //proportion of F_sum attributable to cBn, last X=selpar_n_yrs_wgted yrs  
  number F_cBs_prop;       //proportion of F_sum attributable to cBs, last X=selpar_n_yrs_wgted yrs  
  number F_temp_sum;      //sum of geom mean Fsum's in last X yrs, used to compute F_fishery_prop

  vector F_end(1,nages);
  vector F_end_L(1,nages);  
  number F_end_apex;
  
  number SSB_msy_out;           //SSB (total mature biomass) at msy
  number F_msy_out;             //F at msy
  number msy_mt_out;           //max sustainable yield (1000 mt whole wgt)
  number B_msy_out;             //total biomass at MSY 
  number R_msy_out;             //equilibrium recruitment at F=Fmsy
  number spr_msy_out;           //spr at F=Fmsy

  number F35_dum;				//intermediate calculation for F35
  number F30_dum;				//intermediate calculation for F30
  number F40_dum;				//intermediate calculation for F40
  number F35_out;              	//F35
  number F30_out;              	//F30
  number F40_out;              	//F40
  number SSB_F30_out;  
  number B_F30_out;
  number R_F30_out;
  number L_F30_knum_out;
  number L_F30_mt_out;  
  number rec_mean;  	//arithemetic average recruitment used in SPR-related quantities

  vector N_age_msy(1,nages);         //numbers at age for MSY calculations: beginning of yr
  vector N_age_msy_spawn(1,nages);   //numbers at age for MSY calculations: time of peak spawning  
  vector L_age_msy(1,nages);         //landings at age for MSY calculations
  vector Z_age_msy(1,nages);         //total mortality at age for MSY calculations
  vector F_L_age_msy(1,nages);       //fishing mortality landings (not discards) at age for MSY calculations
  vector F_msy(1,n_iter_msy);        //values of full F to be used in equilibrium calculations
  vector spr_msy(1,n_iter_msy);      //reproductive capacity-per-recruit values corresponding to F values in F_msy
  vector R_eq(1,n_iter_msy);         //equilibrium recruitment values corresponding to F values in F_msy
  vector L_eq_mt(1,n_iter_msy);     //equilibrium landings(1000 mt whole wgt) values corresponding to F values in F_msy
  vector SSB_eq(1,n_iter_msy);       //equilibrium reproductive capacity values corresponding to F values in F_msy
  vector B_eq(1,n_iter_msy);         //equilibrium biomass values corresponding to F values in F_msy
  
  vector FdF_msy(styr,endyr);
  vector FdF30(styr,endyr);
  vector SdSSB_msy(styr,endyr+1);	 
  number SdSSB_msy_end;
  number FdF_msy_end;
  number FdF_msy_end_mean;           //geometric mean of last X yrs  
  vector SdSSB_F30(styr,endyr+1);	  
  number SdSSB_F30_end;
  number FdF30_end_mean;             //geometric mean of last selpar_n_yrs_wgted yrs  
  number Fend_mean_temp;   	     //intermediate calc for geometric mean of last selpar_n_yrs_wgted yrs
  number Fend_mean;		      //geometric mean of last selpar_n_yrs_wgted yrs
  vector L_age_F30(1,nages);         //landings at age for F30 calculations
  
  vector wgt_wgted_L_mt(1,nages);   //fishery-weighted average weight at age of landings in whole weight
  number wgt_wgted_L_denom;          //used in intermediate calculations
  
  number iter_inc_msy;               //increments used to compute msy, equals 1/(n_iter_msy-1)
  
////--------Mortality------------------------------------------------------------------

  vector M(1,nages);                         //age-dependent natural mortality
  matrix M_tv(styr,endyr,1,nages);           //age and time varying natural mortality
 
  matrix F(styr,endyr,1,nages);
  vector Fsum(styr,endyr);                   //Full fishing mortality rate by year
  vector Fapex(styr,endyr);                  //Max across ages, fishing mortality rate by year (may differ from Fsum bc of dome-shaped sel 
  matrix Z(styr,endyr,1,nages);

  init_bounded_number log_avg_F_L_cRn(log_avg_F_cRn_LO,log_avg_F_cRn_HI,log_avg_F_cRn_PH);
  vector log_avg_F_cRn_out(1,8);  
  init_bounded_dev_vector log_dev_F_L_cRn(styr_L_cRn,endyr_L_cRn,log_F_dev_cRn_LO,log_F_dev_cRn_HI,log_F_dev_cRn_PH);
  vector log_F_dev_cRn_out(styr_L_cRn,endyr_L_cRn);
  matrix F_cRn(styr,endyr,1,nages);
  vector F_cRn_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_cRn;
  number log_F_dev_end_cRn; 

  init_bounded_number log_avg_F_L_cRs(log_avg_F_cRs_LO,log_avg_F_cRs_HI,log_avg_F_cRs_PH);
  vector log_avg_F_cRs_out(1,8);  
  init_bounded_dev_vector log_dev_F_L_cRs(styr_L_cRs,endyr_L_cRs,log_F_dev_cRs_LO,log_F_dev_cRs_HI,log_F_dev_cRs_PH);
  vector log_F_dev_cRs_out(styr_L_cRs,endyr_L_cRs);
  matrix F_cRs(styr,endyr,1,nages);
  vector F_cRs_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_cRs;
  number log_F_dev_end_cRs; 

  init_bounded_number log_avg_F_L_cBn(log_avg_F_cBn_LO,log_avg_F_cBn_HI,log_avg_F_cBn_PH);
  vector log_avg_F_cBn_out(1,8);  
  init_bounded_dev_vector log_dev_F_L_cBn(styr_L_cBn,endyr_L_cBn,log_F_dev_cBn_LO,log_F_dev_cBn_HI,log_F_dev_cBn_PH);
  vector log_F_dev_cBn_out(styr_L_cBn,endyr_L_cBn);
  matrix F_cBn(styr,endyr,1,nages);
  vector F_cBn_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_cBn;
  number log_F_dev_end_cBn; 

  init_bounded_number log_avg_F_L_cBs(log_avg_F_cBs_LO,log_avg_F_cBs_HI,log_avg_F_cBs_PH);
  vector log_avg_F_cBs_out(1,8);  
  init_bounded_dev_vector log_dev_F_L_cBs(styr_L_cBs,endyr_L_cBs,log_F_dev_cBs_LO,log_F_dev_cBs_HI,log_F_dev_cBs_PH);
  vector log_F_dev_cBs_out(styr_L_cBs,endyr_L_cBs);
  matrix F_cBs(styr,endyr,1,nages);
  vector F_cBs_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_cBs;
  number log_F_dev_end_cBs; 
 
//---Per-recruit stuff----------------------------------------------------------------------------------
  vector N_age_spr(1,nages);         //numbers at age for SPR calculations: beginning of year
  vector N_age_spr_spawn(1,nages);   //numbers at age for SPR calculations: time of peak spawning  
  vector L_age_spr(1,nages);         //catch at age for SPR calculations
  vector Z_age_spr(1,nages);         //total mortality at age for SPR calculations
  vector spr_static(styr,endyr);     //vector of static SPR values by year
  vector F_L_age_spr(1,nages);       //fishing mortality of landings (not discards) at age for SPR calculations
  vector F_spr(1,n_iter_spr);        //values of full F to be used in per-recruit calculations
  vector spr_spr(1,n_iter_spr);      //reproductive capacity-per-recruit values corresponding to F values in F_spr
  vector spr_ratio(1,n_iter_spr);    //reproductive capacity-per-recruit relative to spr_F0 values corresponding to F values in F_spr
  vector L_spr(1,n_iter_spr);        //landings(lb)-per-recruit (ypr) values corresponding to F values in F_spr

  vector N_spr_F0(1,nages);          //Used to compute spr at F=0: at time of peak spawning
  vector N_bpr_F0(1,nages);          //Used to compute bpr at F=0: at start of year  
  vector N_spr_initial(1,nages);     //Initial spawners per recruit at age given initial F
  vector N_initial_eq(1,nages);      //Initial equilibrium abundance at age
  vector F_initial(1,nages);         //initial F at age
  vector Z_initial(1,nages);         //initial Z at age
  number spr_initial;                //initial spawners per recruit
  number spr_F0;                     //Spawning biomass per recruit at F=0
  number bpr_F0;                     //Biomass per recruit at F=0

  number iter_inc_spr;               //increments used to compute msy, equals max_F_spr_msy/(n_iter_spr-1)

////-------SDNR output-----------------------------------------------------------------------------
 
  number sdnr_lc_nad;
  number sdnr_lc_mad;
  //number sdnr_lc_sad;
 
  number sdnr_ac_cRn;
  number sdnr_ac_cRs;
  number sdnr_ac_cBn;
  number sdnr_ac_cBs;
   
  number sdnr_I_nad;
  number sdnr_I_mad;
  number sdnr_I_sad;
  number sdnr_I_jai;
  number sdnr_I_mareco;

////-------Objective function components-----------------------------------------------------------------------------
  number w_L;
  
  number w_cpue_nad;
  number w_cpue_mad;
  number w_cpue_sad;
  number w_cpue_jai;
  number w_cpue_mareco;
    
  number w_lenc_nad;
  number w_lenc_mad;
  //number w_lc_sad;
  
  number w_agec_cRn;
  number w_agec_cRs;
  number w_agec_cBn;
  number w_agec_cBs;
  
  number w_Nage_init;  
  number w_rec;
  number w_rec_early;
  number w_rec_end;
  number w_fullF;  
  number w_Ftune;

  number f_cRn_L;
  number f_cRs_L;
  number f_cBn_L;
  number f_cBs_L; 
  
  number f_nad_cpue;
  number f_mad_cpue;
  number f_sad_cpue;
  number f_jai_cpue;
  number f_mareco_cpue;
  
  number f_nad_lenc;
  number f_mad_lenc;
  //number f_sad_lenc;
  
  number f_cRn_agec;
  number f_cRs_agec;
  number f_cBn_agec;
  number f_cBs_agec;
  
//  Penalties and constraints. Not all are used.
  number f_Nage_init;              //weight on log devs to estimate initial abundance (excluding first age)
  number f_rec_dev;                //weight on recruitment deviations to fit S-R curve
  number f_rec_dev_early;          //extra weight on deviations in first recruitment stanza
  number f_rec_dev_end;            //extra weight on deviations in ending recruitment stanza
  number f_fullF_constraint;       //penalty for Fapex>X
  number f_Ftune;                  //penalty for tuning F in Ftune yr.  Not applied in final optimization phase.
  number f_priors;                 //prior information on parameters
  
   //init_number xdum;
   objective_function_value fval;
   number fval_data;
   number grad_max;
  
 //--Dummy variables ----
   number denom;                   //denominator used in some calculations
   number numer;                   //numerator used in some calculations

  //---------- Projection quantities--------------------------------------------------------------
  number F_reg_proj;			       //value used to define the projections				
  vector F_proj(styr_proj,endyr_proj);         //F by yr for projections (=F_reg_proj after regulations start, current F till then)  
  vector L_knum_proj(styr_proj,endyr_proj);    //total landings in 1000 fish for projections  
  vector L_mt_proj(styr_proj,endyr_proj);      //total landings in weight (1000 mt) for projections
  
  vector B_proj(styr_proj,endyr_proj);         //Biomass for projections  
  vector SSB_proj(styr_proj,endyr_proj);       //SSB for projections  
  vector R_proj(styr_proj,endyr_proj);         //recruits for projections
  vector FL_age_proj(1,nages);      	       //F (landings) by age for projections    
  
  matrix N_proj(styr_proj,endyr_proj,1,nages);           //Population numbers by year and age at start of yr
  matrix N_spawn_proj(styr_proj,endyr_proj,1,nages);     //Population numbers by year and age at peaking spawning: used for SSB in projections 
  matrix Z_proj(styr_proj,endyr_proj,1,nages);           //Z by year and age for projections 
  matrix L_age_proj(styr_proj,endyr_proj,1,nages);       //Projected landings at age in numbers 
  
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//INITIALIZATION_SECTION
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
GLOBALS_SECTION
  #include "admodel.h"          // Include AD class definitions
  #include "admb2r.cpp"         // Include S-compatible output functions (needs preceding)
  #include <time.h>
	time_t start,finish;
	long hour,minute,second;	
	double elapsed_time;
	
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
RUNTIME_SECTION
 maximum_function_evaluations 1000, 2000,3000, 5000, 10000, 10000, 10000;
 convergence_criteria 1e-2, 1e-2,1e-3, 1e-3, 1e-4, 1e-4, 1e-4;
 
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PRELIMINARY_CALCS_SECTION

// Set values of fixed parameters or set initial guess of estimated parameters
  
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
  
//Fill in sample sizes of comps, possibly sampled in nonconsec yrs 
//Used primarily for output in R object   

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


//fill in Fs for msy and per-recruit analyses
  F_msy(1)=0.0;  
  for (ff=2;ff<=n_iter_msy;ff++) {F_msy(ff)=F_msy(ff-1)+iter_inc_msy;}
  F_spr(1)=0.0;  
  for (ff=2;ff<=n_iter_spr;ff++) {F_spr(ff)=F_spr(ff-1)+iter_inc_spr;}


//fill in F's, Catch matrices, and log rec dev with zero's
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
 
 
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
TOP_OF_MAIN_SECTION
  time(&start);
  arrmblsize=20000000;
  gradient_structure::set_MAX_NVAR_OFFSET(1600);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(2000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(10000);

//>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PROCEDURE_SECTION

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
  
 
FUNCTION get_length_weight_at_age
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

FUNCTION get_reprod 
 
    reprod=elem_prod(elem_prod(prop_f,maturity_f),fecundity);

    for (iyear=styr; iyear<=endyr; iyear++)
    {
      reprod_tv(iyear)=elem_prod(elem_prod(prop_f,tv_maturity(iyear)),fecundity_tv(iyear));
    }

 
FUNCTION get_length_at_age_dist
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
    
FUNCTION get_weight_at_age_landings  ///whole weight in mt
  
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    wholewgt_cR_mt(iyear)=tv_wgt_middle_mt(iyear);
  }  
 
 
FUNCTION get_spr_F0
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



FUNCTION get_selectivity

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

   
FUNCTION get_mortality
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

FUNCTION get_bias_corr
  var_rec_dev=norm2(log_dev_rec(styr_rec_dev,endyr_rec_dev)-
              sum(log_dev_rec(styr_rec_dev,endyr_rec_dev))/nyrs_rec)
              /(nyrs_rec-1.0);                           
  //if (set_BiasCor <= 0.0) {BiasCor=mfexp(var_rec_dev/2.0);}   //bias correction based on empirical residuals
  rec_sigma_sq=square(rec_sigma);
  if (set_BiasCor <= 0.0) {BiasCor=mfexp(rec_sigma_sq/2.0);}   //bias correction based on Rsigma               
  else {BiasCor=set_BiasCor;}

FUNCTION get_numbers_at_age
//Initialization
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

//Initial equilibrium age structure
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

//Compute equilibrium age structure for first year
  N_initial_eq(1)=R1;
  for (iage=2; iage<=nages; iage++)
  {
    N_initial_eq(iage)=N_initial_eq(iage-1)*
        mfexp(-1.0*(Z_initial(iage-1)));    
  }
  //plus group calculation
  N_initial_eq(nages)=N_initial_eq(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
  
//Add deviations to initial equilibrium N
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
  
//Rest of years 
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

   
FUNCTION get_landings_numbers //Baranov catch eqn
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

 
FUNCTION get_landings_wgt
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
 
   
FUNCTION get_catchability_fcns    
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

     
FUNCTION get_indices
//---Predicted CPUEs------------------------

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
FUNCTION get_length_comps

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


FUNCTION get_age_comps
   
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
 
  
////--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_weighted_current 
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

    
FUNCTION get_msy
  
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

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_per_recruit_stuff

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
  
	
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_miscellaneous_stuff

//switch here if var_rec_dev <=dzero 
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

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_projection
  
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
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   


FUNCTION evaluate_objective_function
  //fval=square(xdum-9.0);
  
  fval=0.0;
  fval_data=0.0;
  
//---likelihoods---------------------------

//---Indices-------------------------------

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
  
//---Length comps-------------------------------

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
  
//---Age comps-------------------------------

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
 
//-----------Constraints and penalties--------------------------------
  
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
  
//---Priors---------------------------------------------------
//neg_log_prior arguments: estimate, prior mean, prior var/-CV, pdf type
//Variance input as a negative value is considered to be CV in arithmetic space (CV=-1 implies loose prior) 
//pdf type 1=none, 2=lognormal, 3=normal, 4=beta 
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
  
//----------------------------------------------------------------------------------
//Logistic function: 2 parameters
FUNCTION dvar_vector logistic(const dvar_vector& ages, const dvariable& A50, const dvariable& slope)
  //ages=vector of ages, A50=age at 50% selectivity, slope=rate of increase
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1./(1.+mfexp(-1.*slope*(ages-A50))); //logistic;  
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;

//-----------------------------------------------------------------------------------                                                                           
//Logistic-exponential: 4 parameters (but 1 is fixed)
FUNCTION dvar_vector logistic_exponential(const dvar_vector& ages, const dvariable& A50, const dvariable& slope, const dvariable& sigma, const dvariable& joint)
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

//-----------------------------------------------------------------------------------
//Logistic function: 4 parameters
FUNCTION dvar_vector logistic_double(const dvar_vector& ages, const dvariable& A501, const dvariable& slope1, const dvariable& A502, const dvariable& slope2)
  //ages=vector of ages, A50=age at 50% selectivity, slope=rate of increase, A502=age at 50% decrease additive to A501, slope2=slope of decrease
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=elem_prod( (1./(1.+mfexp(-1.*slope1*(ages-A501)))),(1.-(1./(1.+mfexp(-1.*slope2*(ages-(A501+A502)))))) );     
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;

//-----------------------------------------------------------------------------------
//Jointed logistic function: 6 parameters (increasing and decreasing logistics joined at peak selectivity)
FUNCTION dvar_vector logistic_joint(const dvar_vector& ages, const dvariable& A501, const dvariable& slope1, const dvariable& A502, const dvariable& slope2, const dvariable& satval, const dvariable& joint)
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

//-----------------------------------------------------------------------------------  
//Double Gaussian function: 6 parameters (as in SS3)
FUNCTION dvar_vector gaussian_double(const dvar_vector& ages, const dvariable& peak, const dvariable& top, const dvariable& ascwid, const dvariable& deswid, const dvariable& init, const dvariable& final)
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
    
//-----------------------------------------------------------------------------------    
//Spawner-recruit function (Beverton-Holt or Ricker)
FUNCTION dvariable SR_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& SSB, int func)
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
  
//-----------------------------------------------------------------------------------    
//Spawner-recruit equilibrium function (Beverton-Holt or Ricker)
FUNCTION dvariable SR_eq_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& spr_F, const dvariable& BC, int func)
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

//-----------------------------------------------------------------------------------
//compute multinomial effective sample size for a single yr
FUNCTION dvariable multinom_eff_N(const dvar_vector& pred_comp, const dvar_vector& obs_comp)
  //pred_comp=vector of predicted comps, obscomp=vector of observed comps
  dvariable EffN_Tmp; dvariable numer; dvariable denom;
  RETURN_ARRAYS_INCREMENT();
  numer=sum( elem_prod(pred_comp,(1.0-pred_comp)) );
  denom=sum( square(obs_comp-pred_comp) );
  if (denom>0.0) {EffN_Tmp=numer/denom;}
  else {EffN_Tmp=-missing;}                            
  RETURN_ARRAYS_DECREMENT();
  return EffN_Tmp;

//-----------------------------------------------------------------------------------
//Likelihood contribution: lognormal
FUNCTION dvariable lk_lognormal(const dvar_vector& pred, const dvar_vector& obs, const dvar_vector& cv, const dvariable& wgt_dat)
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

//-----------------------------------------------------------------------------------
//Likelihood contribution: multinomial
FUNCTION dvariable lk_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const double& minSS, const dvariable& wgt_dat)
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

//-----------------------------------------------------------------------------------
//Likelihood contribution: robust multinomial
FUNCTION dvariable lk_robust_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS, const dvariable& wgt_dat)
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
  
//-----------------------------------------------------------------------------------
//Likelihood contribution: Dirichlet-multinomial
FUNCTION dvariable lk_dirichlet_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS, const dvariable& log_dir_par)
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

//-----------------------------------------------------------------------------------
//Likelihood contribution: logistic normal (aka multivariate logistic in iSCAM; logistic normal in Francis' terminology)
FUNCTION dvariable lk_logistic_normal(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS)
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
  
//-----------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------
//Likelihood contribution: priors
FUNCTION  dvariable neg_log_prior(dvariable pred, const double& prior, dvariable var, int pdf)
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

//-----------------------------------------------------------------------------------
//SDNR: age comp likelihood (assumes fits are done with the robust multinomial function)
FUNCTION dvariable sdnr_multinomial(const double& ncomp, const dvar_vector& ages, const dvar_vector& nsamp, 
                                    const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const dvariable& wgt_dat)
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

//-----------------------------------------------------------------------------------
//SDNR: lognormal likelihood
FUNCTION dvariable sdnr_lognormal(const dvar_vector& pred, const dvar_vector& obs, const dvar_vector& cv, const dvariable& wgt_dat)
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

//-----------------------------------------------------------------------------------
REPORT_SECTION

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
  
 
