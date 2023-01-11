//##  Author: NMFS, Beaufort Lab, Sustainable Fisheries Branch
//##  Analyst: Kate Siegfried
//##  Species: Black Sea Bass
//##  Region: US South Atlantic
//##  SEDAR: 56
//##  Date: 2022-11-18 22:17:15


//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##
//##  SEDAR56 Standard Assessment: Black Sea Bass, 2017 - 8
//##
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

//Reg blocks -- 3 possible phases of cGN size regs: styr-83 no restrictions, 1984-98 8-inch TL, 1999-2010 10-in TL
//           -- 4 possible phases of recr size regs: styr-83 no restrictions, 1984-98 8-inch TL, 1999-2006 10-in TL, 2007-2010 12-in TL 
init_int endyr_selex_phase1;
init_int endyr_selex_phase2;
init_int endyr_selex_phase3;
init_int endyr_selex_phase4;
//first yr commercial fisheries were closed due to quotas
init_int styr_cGN_closed;  
//size limits
init_number sizelim1;  //8 inch limit in mm
init_number sizelim2;  //10 inch limit in mm
init_number sizelim3;  //11 inch limit in mm
init_number sizelim4;  //12 inch limit in mm
init_number sizelim5;  //13 inch limit in mm
init_number limit_disc; //max size applied to discards in block one, prior to fed regs

//Total number of ages
init_int nages;

// Vector of ages for age bins
init_vector agebins(1,nages);

//number assessment years
number nyrs;
number nyrs_rec;
//this section MUST BE INDENTED!!!
 LOCAL_CALCS
   nyrs=endyr-styr+1.;
   //nyrs_rec=endyr-styr_rec_dev+1.;
   nyrs_rec=endyr_rec_dev-styr_rec_dev+1.;
 END_CALCS

//Total number of length bins for each matrix and length bins used to compute mass in largest bin (plus group)
init_int nlenbins;       //used to match data
init_int nlenbins_plus;  //used to compute density of largest bin (plus group)
init_number lenbins_width;  //width of length bins (mm)

//Vector of lengths for length bins (mm)(midpoint) and bins used in computation of plus group
init_ivector lenbins(1,nlenbins);
init_ivector lenbins_plus(1,nlenbins_plus);
int nlenbins_all;    //largest size class used to compute average lengths and weights

//this section MUST BE INDENTED!!!
 LOCAL_CALCS
   nlenbins_all=nlenbins+nlenbins_plus;   
 END_CALCS
  

//Max F used in spr and msy calcs
init_number max_F_spr_msy;
//Total number of iterations for spr calcs
init_int n_iter_spr;
//Total number of iterations for msy calcs
init_int n_iter_msy;
 
 LOCAL_CALCS
		n_iter_msy=n_iter_spr; 
 END_CALCS

//Number years at end of time series over which to average sector F's, for weighted selectivities
init_int selpar_n_yrs_wgted;
//bias correction (set to 1.0 for no bias correction or a negative value to compute from rec variance)
init_number set_BiasCor;
//exclude these years from end of time series for computing bias correction
init_number BiasCor_exclude_yrs;

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: observed data section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

//###MARMAP bft/fst#################
//CPUE
init_int styr_cpue_sBT;
init_int endyr_cpue_sBT;
init_vector obs_cpue_sBT(styr_cpue_sBT,endyr_cpue_sBT);   //Observed CPUE
init_vector obs_cv_cpue_sBT(styr_cpue_sBT,endyr_cpue_sBT);    //CV of cpue

// Length Compositions (1 cm bins)
init_int nyr_lenc_sBT;
init_ivector yrs_lenc_sBT(1,nyr_lenc_sBT);
init_vector nsamp_lenc_sBT(1,nyr_lenc_sBT);
init_vector nfish_lenc_sBT(1,nyr_lenc_sBT);
init_matrix obs_lenc_sBT(1,nyr_lenc_sBT,1,nlenbins);

// Age Compositions 
init_int nyr_agec_sBT;
init_ivector yrs_agec_sBT(1,nyr_agec_sBT);
init_vector nsamp_agec_sBT(1,nyr_agec_sBT);
init_vector nfish_agec_sBT(1,nyr_agec_sBT);
init_matrix obs_agec_sBT(1,nyr_agec_sBT,1,nages);


//###SERFS Video#################
//CPUE
//init_int styr_Vid_cpue;
//init_int endyr_Vid_cpue;
//init_vector obs_Vid_cpue(styr_Vid_cpue,endyr_Vid_cpue);   //Observed CPUE
//init_vector Vid_cpue_cv(styr_Vid_cpue,endyr_Vid_cpue);    //CV of cpue

//###MARMAP cvt/sTV#################
//CPUE
init_int styr_cpue_sTV;
init_int endyr_cpue_sTV;
init_vector obs_cpue_sTV(styr_cpue_sTV,endyr_cpue_sTV);   //Observed CPUE
init_vector obs_cv_cpue_sTV(styr_cpue_sTV,endyr_cpue_sTV);    //CV of cpue

// Age Compositions 
init_int nyr_agec_sTV;
init_ivector yrs_agec_sTV(1,nyr_agec_sTV);
init_vector nsamp_agec_sTV(1,nyr_agec_sTV);
init_vector nfish_agec_sTV(1,nyr_agec_sTV);
init_matrix obs_agec_sTV(1,nyr_agec_sTV,1,nages);

//###################Commercial Hook and Line fishery #########################
//CPUE
init_int styr_cpue_cHL;                                             
init_int endyr_cpue_cHL;                                            
init_vector obs_cpue_cHL(styr_cpue_cHL,endyr_cpue_cHL);//Observed CPUE
init_vector obs_cv_cpue_cHL(styr_cpue_cHL,endyr_cpue_cHL); //CV of cpue
 
// Landings  (1000 lb whole weight)
init_int styr_L_cHL;
init_int endyr_L_cHL;
init_vector obs_L_cHL(styr_L_cHL,endyr_L_cHL);   //vector of observed landings by year 
init_vector obs_cv_L_cHL(styr_L_cHL,endyr_L_cHL);    //vector of CV of landings by year

// Discards (1000 fish)
init_int styr_D_cHL;
init_int endyr_D_cHL;
init_vector obs_released_cHL(styr_D_cHL,endyr_D_cHL); //vector of observed releases by year, gets multiplied by discard mortality for fitting
init_vector obs_cv_D_cHL(styr_D_cHL,endyr_D_cHL);         //vector of CV of discards by year
// Discards (1000 fish) during closed season
init_int styr_D_cHL_closed;
init_int endyr_D_cHL_closed;
init_vector obs_released_cHL_closed(styr_D_cHL_closed,endyr_D_cHL_closed); //vector of observed releases by year, gets multiplied by discard mortality for fitting

// Length Compositions (1 cm bins)
init_int nyr_lenc_cHL;
init_ivector yrs_lenc_cHL(1,nyr_lenc_cHL);
init_vector nsamp_lenc_cHL(1,nyr_lenc_cHL);
init_vector nfish_lenc_cHL(1,nyr_lenc_cHL);
init_matrix obs_lenc_cHL(1,nyr_lenc_cHL,1,nlenbins);
// Age Compositions 
init_int nyr_agec_cHL;
init_ivector yrs_agec_cHL(1,nyr_agec_cHL);
init_vector nsamp_agec_cHL(1,nyr_agec_cHL);
init_vector nfish_agec_cHL(1,nyr_agec_cHL);
init_matrix obs_agec_cHL(1,nyr_agec_cHL,1,nages);


//#############################################################################
//##Commercial pot (+ other) fleet 
// Landings (1000 lb whole weight)
init_int styr_L_cPT;
init_int endyr_L_cPT;
init_vector obs_L_cPT(styr_L_cPT,endyr_L_cPT);
init_vector obs_cv_L_cPT(styr_L_cPT,endyr_L_cPT);    //vector of CV of landings by year

// Discards (1000 fish)
init_int styr_D_cPT;
init_int endyr_D_cPT;
init_vector obs_released_cPT(styr_D_cPT,endyr_D_cPT); //vector of observed releases by year, gets multiplied by discard mortality for fitting
init_vector obs_cv_D_cPT(styr_D_cPT,endyr_D_cPT);         //vector of CV of discards by year
// Discards (1000 fish) during closed season
init_int styr_D_cPT_closed;
init_int endyr_D_cPT_closed;
init_vector obs_released_cPT_closed(styr_D_cPT_closed,endyr_D_cPT_closed); //vector of observed releases by year, gets multiplied by discard mortality for fitting

// Length Compositions (1 cm bins)
init_int nyr_lenc_cPT;
init_ivector yrs_lenc_cPT(1,nyr_lenc_cPT);
init_vector nsamp_lenc_cPT(1,nyr_lenc_cPT);
init_vector nfish_lenc_cPT(1,nyr_lenc_cPT);
init_matrix obs_lenc_cPT(1,nyr_lenc_cPT,1,nlenbins);
init_int nyr_lenc_pool_cPT;     //years and weights to pool predicted cPT length comps to match pooled observations
init_ivector yrs_lenc_pool_cPT(1,nyr_lenc_pool_cPT);
init_vector nsamp_lenc_pool_cPT(1,nyr_lenc_pool_cPT);

// Age Compositions
init_int nyr_agec_cPT;
init_ivector yrs_agec_cPT(1,nyr_agec_cPT);
init_vector nsamp_agec_cPT(1,nyr_agec_cPT);
init_vector nfish_agec_cPT(1,nyr_agec_cPT);
init_matrix obs_agec_cPT(1,nyr_agec_cPT,1,nages);

//#############################################################################
//#############################################################################
//##Commercial Trawl fleet 
// Landings (1000 lb whole weight)
init_int styr_L_cTW;
init_int endyr_L_cTW;
init_vector obs_L_cTW(styr_L_cTW,endyr_L_cTW);
init_vector obs_cv_L_cTW(styr_L_cTW,endyr_L_cTW);    //vector of CV of landings by year


//#############################################################################
//################################Headboat fleet ########################################
//CPUE
init_int styr_cpue_rHB;
init_int endyr_cpue_rHB;
init_vector obs_cpue_rHB(styr_cpue_rHB,endyr_cpue_rHB);//Observed CPUE
init_vector obs_cv_cpue_rHB(styr_cpue_rHB,endyr_cpue_rHB); //CV of cpue
//###rHD index (headboat discards from at sea observer program#################
//init_int styr_rHD_cpue;
//init_int endyr_rHD_cpue;
//init_vector obs_rHD_cpue(styr_rHD_cpue,endyr_rHD_cpue);   //Observed CPUE
//init_vector rHD_cpue_cv(styr_rHD_cpue,endyr_rHD_cpue);    //CV of cpue
// Landings (1000 lb)
init_int styr_L_rHB;
init_int endyr_L_rHB;
init_vector obs_L_rHB(styr_L_rHB,endyr_L_rHB);
init_vector obs_cv_L_rHB(styr_L_rHB,endyr_L_rHB);
// Discards (1000s)
init_int styr_D_rHB;
init_int endyr_D_rHB;
init_vector obs_released_rHB(styr_D_rHB,endyr_D_rHB);  //vector of observed releases by year, multiplied by discard mortality for fitting  
init_vector obs_cv_D_rHB(styr_D_rHB,endyr_D_rHB);          //vector of CV of discards by year
// Length Compositions (1 cm bins) of landings
init_int nyr_lenc_rHB;
init_ivector yrs_lenc_rHB(1,nyr_lenc_rHB);
init_vector nsamp_lenc_rHB(1,nyr_lenc_rHB);
init_vector nfish_lenc_rHB(1,nyr_lenc_rHB);
init_matrix obs_lenc_rHB(1,nyr_lenc_rHB,1,nlenbins);
// Age compositions of landings
init_int nyr_agec_rHB;
init_ivector yrs_agec_rHB(1,nyr_agec_rHB);
init_vector nsamp_agec_rHB(1,nyr_agec_rHB);
init_vector nfish_agec_rHB(1,nyr_agec_rHB);
init_matrix obs_agec_rHB(1,nyr_agec_rHB,1,nages);
//Length Compositions (1 cm bins) of rHB discards
init_int nyr_lenc_rHB_D;
init_ivector yrs_lenc_rHB_D(1,nyr_lenc_rHB_D);
init_vector nsamp_lenc_rHB_D(1,nyr_lenc_rHB_D);
init_vector nfish_lenc_rHB_D(1,nyr_lenc_rHB_D);
init_matrix obs_lenc_rHB_D(1,nyr_lenc_rHB_D,1,nlenbins);

//#############################################################################
//############################rGN recreational fleet #################################
// Landings (1000 lb)
init_int styr_L_rGN;
init_int endyr_L_rGN;
init_vector obs_L_rGN(styr_L_rGN,endyr_L_rGN);
init_vector obs_cv_L_rGN(styr_L_rGN,endyr_L_rGN);
// Discards (1000s)
init_int styr_D_rGN;
init_int endyr_D_rGN;
init_vector obs_released_rGN(styr_D_rGN,endyr_D_rGN); //vector of observed releases by year, multiplied by discard mortality for fitting 
init_vector obs_cv_D_rGN(styr_D_rGN,endyr_D_rGN);         //vector of CV of discards by year
// Length Compositions (1 cm bins)
init_int nyr_lenc_rGN;
init_ivector yrs_lenc_rGN(1,nyr_lenc_rGN);
init_vector nsamp_lenc_rGN(1,nyr_lenc_rGN);
init_vector nfish_lenc_rGN(1,nyr_lenc_rGN);
init_matrix obs_lenc_rGN(1,nyr_lenc_rGN,1,nlenbins);
//init_int nyr_rGN_lenc_pool;     //years and weights to pool predicted rGN length comps to match pooled observations
//init_ivector yrs_rGN_lenc_pool(1,nyr_rGN_lenc_pool);
//init_vector nsamp_rGN_lenc_pool(1,nyr_rGN_lenc_pool);
// Age Compositions 
//init_int nyr_rGN_agec;
//init_ivector yrs_rGN_agec(1,nyr_rGN_agec);
//init_vector nsamp_rGN_agec(1,nyr_rGN_agec);
//init_vector nfish_rGN_agec(1,nyr_rGN_agec);
//init_matrix obs_rGN_agec(1,nyr_rGN_agec,1,nages);

//Discard mortality constants
init_number set_Dmort_HL;   //handline (commercial)
init_number set_Dmort_rHB_HL; //headboat-specific hook and line
init_number set_Dmort_rGN_HL; //charterboat and private hook and line
init_number set_Dmort_cPT1;  //pots 1.5 inch panel
init_number set_Dmort_cPT2;  //pots 2.0 inch panel
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: parameter section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##################Single Parameter values and initial guesses #################################

// Von Bert parameters in TL mm 
init_vector set_Linf(1,7);
init_vector set_K(1,7);
init_vector set_t0(1,7);
init_vector set_len_cv(1,7);
//init_number set_len_cv_se(1,7);
//Scalar used only for computing MSST.
init_vector set_M_constant(1,7);     //age-independent: used only for MSST and to scale age dependent M, prior if M is estimated

//Standard errors of von bert params
//init_number set_Linf_se;
//init_number set_K_se;
//init_number set_t0_se;
//init_number set_len_cv_se;

//Spawner-recruit parameters (Initial guesses or fixed values)
init_vector set_steep(1,7);         //recruitment steepness
//init_number set_steep_se;      //SE of recruitment steepness
//init_int steep_prior_pdf;      //(1=none, 2=lognormal, 3=normal, 4=beta)
init_vector set_log_R0(1,7);        //recruitment R0
init_vector set_R_autocorr(1,7);    //recruitment autocorrelation
init_vector set_rec_sigma(1,7);     //recruitment standard deviation in log space
//init_number set_rec_sigma_se;  //SE of recruitment standard deviation in log space
//init_int rec_sigma_prior_pdf;      //(1=none, 2=lognormal, 3=normal, 4=beta)

init_vector set_log_dm_lenc_sBT(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_lenc_sTV(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_lenc_cHL(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_lenc_cPT(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_lenc_rHB(1,7);  //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_lenc_rHB_D(1,7);   //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_lenc_rGN(1,7);   //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_agec_sBT(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_agec_sTV(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_agec_cHL(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_agec_cPT(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_agec_rHB(1,7);   //Dirichlet-multinomial overdispersion parameter
//init_vector set_log_dm_rGN_ac(1,7);   //Dirichlet-multinomial overdispersion parameter

//Initial guesses or fixed values of estimated selectivity parameters

init_vector set_selpar_A50_sBT(1,7);
init_vector set_selpar_slope_sBT(1,7);

init_vector set_selpar_A50_sTV(1,7);
init_vector set_selpar_slope_sTV(1,7);

init_vector set_selpar_A50_Vid(1,7);
init_vector set_selpar_slope_Vid(1,7);

init_vector set_selpar_A50_cHL2(1,7);
init_vector set_selpar_slope_cHL2(1,7);
init_vector set_selpar_A50_cHL3(1,7);
init_vector set_selpar_slope_cHL3(1,7);
init_vector set_selpar_A50_cHL4(1,7);
init_vector set_selpar_slope_cHL4(1,7);

init_vector set_selpar_A50_cPT2(1,7);
init_vector set_selpar_slope_cPT2(1,7);
init_vector set_selpar_A50_cPT3(1,7);
init_vector set_selpar_slope_cPT3(1,7);
init_vector set_selpar_A50_cPT4(1,7);
init_vector set_selpar_slope_cPT4(1,7);

init_vector set_selpar_A50_rHB1(1,7);
init_vector set_selpar_slope_rHB1(1,7);
init_vector set_selpar_A50_rHB2(1,7);
init_vector set_selpar_slope_rHB2(1,7);
init_vector set_selpar_A50_rHB3(1,7);
init_vector set_selpar_slope_rHB3(1,7);
init_vector set_selpar_A50_rHB4(1,7);
init_vector set_selpar_slope_rHB4(1,7);
init_vector set_selpar_A50_rHB5(1,7);
init_vector set_selpar_slope_rHB5(1,7);

init_vector set_selpar_A50_rGN1(1,7);
init_vector set_selpar_slope_rGN1(1,7);
init_vector set_selpar_A50_rGN2(1,7);
init_vector set_selpar_slope_rGN2(1,7);
init_vector set_selpar_A50_rGN3(1,7);
init_vector set_selpar_slope_rGN3(1,7);
init_vector set_selpar_A50_rGN4(1,7);
init_vector set_selpar_slope_rGN4(1,7);
init_vector set_selpar_A50_rGN5(1,7);
init_vector set_selpar_slope_rGN5(1,7);

init_vector set_selpar_logit_Age0_rHB_D(1,7);
init_vector set_selpar_logit_Age1_rHB_D(1,7);
init_vector set_selpar_logit_Age2_rHB_D(1,7);

init_vector set_selpar_A50_rHD4(1,7);
init_vector set_selpar_slope_rHD4(1,7);
init_vector set_selpar_A502_rHD4(1,7);
init_vector set_selpar_slope2_rHD4(1,7);
init_vector set_selpar_A50_rHD5(1,7);
init_vector set_selpar_slope_rHD5(1,7);
init_vector set_selpar_A502_rHD5(1,7);
init_vector set_selpar_slope2_rHD5(1,7);

//--index catchability------------------------------------------------------------------------------------------------------------
init_vector set_log_q_cpue_sBT(1,7);    //catchability coefficient (log) for Blackfish trap
init_vector set_log_q_cpue_sTV(1,7);    //catchability coefficient (log) for MARMAP chevron trap
//init_vector set_logq_Vid(1,7);    //catchability coefficient (log) for SERFS video
init_vector set_log_q_cpue_cHL(1,7);      //catchability coefficient (log) for commercial logbook index
init_vector set_log_q_cpue_rHB(1,7);      //catchability coefficient (log) for the headboat index
//init_vector set_logq_rHD(1,7);     //catchability coefficient (log) for rHD

////--mean F's in log space--------------------------------
init_vector set_log_avg_F_L_cHL(1,7);
init_vector set_log_avg_F_L_cPT(1,7);
init_vector set_log_avg_F_L_cTW(1,7);
init_vector set_log_avg_F_L_rHB(1,7);
init_vector set_log_avg_F_L_rGN(1,7);
////--discard F's-----------------------
init_vector set_log_avg_F_D_cGN(1,7);
init_vector set_log_avg_F_D_rHB(1,7);
init_vector set_log_avg_F_D_rGN(1,7);

////Dev's-------------------------------------
init_vector set_log_dev_F_L_cHL(1,3);
init_vector set_log_dev_F_L_cPT(1,3);
init_vector set_log_dev_F_L_cTW(1,3);
init_vector set_log_dev_F_L_rHB(1,3);
init_vector set_log_dev_F_L_rGN(1,3);

init_vector set_log_dev_F_D_cGN(1,3);
init_vector set_log_dev_F_D_rHB(1,3);
init_vector set_log_dev_F_D_rGN(1,3);

init_vector set_log_dev_RWq(1,3);
init_vector set_log_dev_rec(1,3);
init_vector set_log_dev_Nage(1,3);

init_vector set_log_dev_vals_F_L_cHL(styr_L_cHL,endyr_L_cHL);
init_vector set_log_dev_vals_F_L_cPT(styr_L_cPT,endyr_L_cPT);
//init_vector set_log_F_dev_cTW_vals(styr_L_cTW,endyr_L_cTW);
init_vector set_log_dev_vals_F_L_rHB(styr_L_rHB,endyr_L_rHB);
init_vector set_log_dev_vals_F_L_rGN(styr_L_rGN,endyr_L_rGN);
init_vector set_log_dev_vals_F_D_cGN(styr_D_cHL,endyr_D_cHL);
init_vector set_log_dev_vals_F_D_rHB(styr_D_rHB,endyr_D_rHB);
init_vector set_log_dev_vals_F_D_rGN(styr_D_rGN,endyr_D_rGN);
init_vector set_log_dev_vals_rec(styr_rec_dev,endyr_rec_dev);
init_vector set_log_dev_vals_Nage(2,nages);   

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: likelihood weights section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

//--weights for likelihood components-------------------------------------------------------------------------------
init_number set_w_L;
init_number set_w_D;

init_number set_w_lenc_sBT; //Won't need with the dirichlet
init_number set_w_lenc_cHL;
init_number set_w_lenc_cPT;
init_number set_w_lenc_rHB;
init_number set_w_lenc_rHB_D;
init_number set_w_lenc_rGN;

init_number set_w_agec_sBT;  //Won't need with the dirichlet
init_number set_w_agec_sTV;
init_number set_w_agec_cHL;
init_number set_w_agec_cPT;
init_number set_w_agec_rHB;
//init_number set_w_ac_rGN;

init_number set_w_cpue_sBT;
init_number set_w_cpue_sTV;
//init_number set_w_I_Vid;
init_number set_w_cpue_cHL;
init_number set_w_cpue_rHB;
//init_number set_w_I_rHD;

init_number set_w_rec;             //for fitting S-R curve
init_number set_w_rec_early;       //additional constraint on early years recruitment
init_number set_w_rec_end;         //additional constraint on ending years recruitment 
init_number set_w_fullF;           //penalty for any Fapex>3(removed in final phase of optimization)
init_number set_w_Ftune;           //weight applied to tuning F (removed in final phase of optimization)
//init_number set_w_cvlen_dev;         //penalty on cv deviations at age
//init_number set_w_cvlen_diff;       //penalty on first difference of cv deviations at age

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: miscellaneous stuff section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

//TL(mm)-weight(whole weight in kg) relationship: W=aL^b
init_number wgtpar_a;
init_number wgtpar_b;

//weight(whole weight)- fecundity (units=eggs/batch) relationship: log(y)=a+bW
init_number fecpar_a;
init_number fecpar_b;
init_number fecpar_batches;  //number of annual batches may be important if age-dependent, otherwise just a scalar
init_number fecpar_scale;    //used for scaling annual egg production (10^X eggs)

//Female maturity and proportion female at age
init_vector obs_maturity_f(1,nages);            //proportion females mature at age
init_vector obs_prop_f(1,nages);                //proportion female at age
init_number spawn_time_frac; //time of year of peak spawning, as a fraction of the year

init_vector set_M(1,nages);     //age-dependent: used in model
init_number max_obs_age;        //max observed age, used to scale M

//value or initial guess for recreational rHB and rGN historic landings multiplicative bias, last yr the bias applies  
//this feature is not currently implemented in the likelihood fcn
init_number set_L_rHB_bias;
init_number set_L_rGN_bias;
init_number set_L_cGN_bias;
init_number endyr_L_rHB_bias;
init_number endyr_L_rGN_bias;
init_number endyr_L_cGN_bias;

//rate of increase on q
init_int set_q_rate_phase;  //value sets estimation phase of rate increase, negative value turns it off
init_number set_q_rate;
//density dependence on fishery q's 
init_int set_q_DD_phase;      //value sets estimation phase of random walk, negative value turns it off
init_number set_q_DD_beta;    //value of 0.0 is density indepenent
init_number set_q_DD_beta_se;
init_int set_q_DD_stage;      //age to begin counting biomass, should be near full exploitation

//random walk on fishery q's 
init_int set_q_RW_phase;         //value sets estimation phase of random walk, negative value turns it off
init_number set_q_RW_cHL_var;     //assumed variance of RW q
init_number set_q_RW_rHB_var;     //assumed variance of RW q
//init_number set_q_RW_rHD_var;    //assumed variance of RW q

//init_vector set_F_init_ratio;  //defines initialization F as a ratio of that from first several yrs of assessment
init_number set_F_init_ratio

//Tune Fapex (tuning removed in final year of optimization)
init_number set_Ftune;
init_int set_Ftune_yr;

//threshold sample sizes for including length comps, age comps, respectively 
init_number minSS_lenc;
init_number minSS_agec;

//maximum allowable annual sample sizes for length comps, age comps, respectively
init_number maxSS_lenc;
init_number maxSS_agec;

//ageing error matrix (columns are true ages, rows are ages as read for age comps: columns should sum to one)
init_matrix age_error(1,nages,1,nages);

//proportion of length comp mass below size limit considered when matching length comp
//note: these need length comp and age comp data to be estimable
init_number set_p_lenc_cHL2; 
init_number set_p_lenc_cHL3; 
init_number set_p_lenc_cPT2;
init_number set_p_lenc_cPT3;
init_number set_p_lenc_cTW2;
init_number set_p_lenc_cTW3;
init_number set_p_lenc_rHB2;
init_number set_p_lenc_rHB3;
init_number set_p_lenc_rHB4;
init_number set_p_lenc_rHB5; 
init_number set_p_lenc_rGN2;  
init_number set_p_lenc_rGN3;
init_number set_p_lenc_rGN4;
init_number set_p_lenc_rGN5;

init_number set_p_lenc_cGN_D2;
init_number set_p_lenc_cGN_D3;
init_number set_p_lenc_cGN_D4;
init_number set_p_lenc_rHB_D2;
init_number set_p_lenc_rHB_D3;
init_number set_p_lenc_rHB_D4;
init_number set_p_lenc_rHB_D5;
init_number set_p_lenc_rGN_D1;
init_number set_p_lenc_rGN_D2;
init_number set_p_lenc_rGN_D3;
init_number set_p_lenc_rGN_D4;
init_number set_p_lenc_rGN_D5;

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
number onehalf;                 //0.5

init_number end_of_data_file;
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


//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PARAMETER_SECTION
LOCAL_CALCS
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
  
  const double selpar_A50_rHD4_LO=set_selpar_A50_rHD4(2); 
  const double selpar_A50_rHD4_HI=set_selpar_A50_rHD4(3); 
  const double selpar_A50_rHD4_PH=set_selpar_A50_rHD4(4);
  const double selpar_slope_rHD4_LO=set_selpar_slope_rHD4(2); 
  const double selpar_slope_rHD4_HI=set_selpar_slope_rHD4(3); 
  const double selpar_slope_rHD4_PH=set_selpar_slope_rHD4(4);
  const double selpar_A502_rHD4_LO=set_selpar_A502_rHD4(2); 
  const double selpar_A502_rHD4_HI=set_selpar_A502_rHD4(3); 
  const double selpar_A502_rHD4_PH=set_selpar_A502_rHD4(4);
  const double selpar_slope2_rHD4_LO=set_selpar_slope2_rHD4(2); 
  const double selpar_slope2_rHD4_HI=set_selpar_slope2_rHD4(3); 
  const double selpar_slope2_rHD4_PH=set_selpar_slope2_rHD4(4);
  const double selpar_A50_rHD5_LO=set_selpar_A50_rHD5(2); 
  const double selpar_A50_rHD5_HI=set_selpar_A50_rHD5(3); 
  const double selpar_A50_rHD5_PH=set_selpar_A50_rHD5(4);
  const double selpar_slope_rHD5_LO=set_selpar_slope_rHD5(2); 
  const double selpar_slope_rHD5_HI=set_selpar_slope_rHD5(3); 
  const double selpar_slope_rHD5_PH=set_selpar_slope_rHD5(4);
  const double selpar_A502_rHD5_LO=set_selpar_A502_rHD5(2); 
  const double selpar_A502_rHD5_HI=set_selpar_A502_rHD5(3); 
  const double selpar_A502_rHD5_PH=set_selpar_A502_rHD5(4);
  const double selpar_slope2_rHD5_LO=set_selpar_slope2_rHD5(2); 
  const double selpar_slope2_rHD5_HI=set_selpar_slope2_rHD5(3); 
  const double selpar_slope2_rHD5_PH=set_selpar_slope2_rHD5(4);
  
  const double log_q_sBT_LO=set_log_q_cpue_sBT(2); const double log_q_sBT_HI=set_log_q_cpue_sBT(3); const double log_q_sBT_PH=set_log_q_cpue_sBT(4);
  const double log_q_rHB_LO=set_log_q_cpue_rHB(2); const double log_q_rHB_HI=set_log_q_cpue_rHB(3); const double log_q_rHB_PH=set_log_q_cpue_rHB(4);
  const double log_q_sTV_LO=set_log_q_cpue_sTV(2); const double log_q_sTV_HI=set_log_q_cpue_sTV(3); const double log_q_sTV_PH=set_log_q_cpue_sTV(4);
  //const double log_q_Vid_LO=set_logq_Vid(2); const double log_q_Vid_HI=set_logq_Vid(3); const double log_q_Vid_PH=set_logq_Vid(4);
  const double log_q_cHL_LO=set_log_q_cpue_cHL(2); const double log_q_cHL_HI=set_log_q_cpue_cHL(3); const double log_q_cHL_PH=set_log_q_cpue_cHL(4);
  //const double log_q_rHD_LO=set_logq_rHD(2); const double log_q_rHD_HI=set_logq_rHD(3); const double log_q_rHD_PH=set_logq_rHD(4);
  
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
  
 END_CALCS
 
////--------------Growth---------------------------------------------------------------------------
 
  //Population growth parms and conversions
  init_bounded_number Linf(Linf_LO,Linf_HI,Linf_PH);
  init_bounded_number K(K_LO,K_HI,K_PH);
  init_bounded_number t0(t0_LO,t0_HI,t0_PH);
  init_bounded_number len_cv_val(len_cv_LO,len_cv_HI,len_cv_PH);  
  vector Linf_out(1,8);
  vector K_out(1,8);
  vector t0_out(1,8);
  vector len_cv_val_out(1,8);
  
  ////init_bounded_number Linf(300,800,3);
  ////init_bounded_number K(0.05,0.5,3);
  ////init_bounded_number t0(-1.5,-0.01,3);
  //number Linf;
  //number K;
  //number t0;  
  //init_bounded_number len_cv_val(0.07,0.3,4);
  ////  init_bounded_dev_vector log_len_cv_dev(1,nages,-2,2,-4)
  ////number len_cv_val;
  //vector len_sd(1,nages);
  //vector len_cv(1,nages);  
  
  vector meanlen_TL(1,nages);   //mean Total length (mm) at age
  
  vector wgt_g(1,nages);        //whole wgt in g
  vector wgt_kg(1,nages);       //whole wgt in kg
  vector wgt_mt(1,nages);       //whole wgt in mt
  vector wgt_klb(1,nages);      //whole wgt in 1000 lb
  vector wgt_lb(1,nages);       //whole wgt in lb
  vector fecundity(1,nages);    //fecundity at age, perhaps scaled  

  matrix len_cHL_mm(styr,endyr,1,nages);       //mean length at age of cHL landings in mm (may differ from popn mean)    
  matrix wgt_cHL_klb(styr,endyr,1,nages);      //whole wgt of cHL landings in 1000 lb  
  matrix len_cPT_mm(styr,endyr,1,nages);       //mean length at age of cPT landings in mm (may differ from popn mean)    
  matrix wgt_cPT_klb(styr,endyr,1,nages);      //whole wgt of cPT landings in 1000 lb  
  matrix len_cTW_mm(styr,endyr,1,nages);       //mean length at age of cTW landings in mm (may differ from popn mean)    
  matrix wgt_cTW_klb(styr,endyr,1,nages);      //whole wgt of cTW landings in 1000 lb  
  matrix len_rHB_mm(styr,endyr,1,nages);       //mean length at age of rHB landings in mm (may differ from popn mean)    
  matrix wgt_rHB_klb(styr,endyr,1,nages);      //whole wgt of rHB landings in 1000 lb  
  matrix len_rGN_mm(styr,endyr,1,nages);     //mean length at age of rGN landings in mm (may differ from popn mean)    
  matrix wgt_rGN_klb(styr,endyr,1,nages);    //whole wgt of rGN landings in 1000 lb  

  matrix len_cGN_D_mm(styr,endyr,1,nages);    //mean length at age of cHL discards in mm (may differ from popn mean)    
  matrix wgt_cGN_D_klb(styr,endyr,1,nages);   //whole wgt of cHL discards in 1000 lb  
  matrix len_rHB_D_mm(styr,endyr,1,nages);      //mean length at age of cHL discards in mm (may differ from popn mean)    
  matrix wgt_rHB_D_klb(styr,endyr,1,nages);     //whole wgt of cHL discards in 1000 lb  
  matrix len_rGN_D_mm(styr,endyr,1,nages);    //mean length at age of cHL discards in mm (may differ from popn mean)    
  matrix wgt_rGN_D_klb(styr,endyr,1,nages);   //whole wgt of cHL discards in 1000 lb  
 
  matrix lenprob(1,nages,1,nlenbins);           //distn of size at age (age-length key, 1 cm bins) in population
    
  matrix lenprob_plus(1,nages,1,nlenbins_plus); //used to compute mass in last length bin (a plus group)   
  matrix lenprob_all(1,nages,1,nlenbins_all);   //extended lenprob
  vector lenbins_all(1,nlenbins_all);
  
  number zscore_len;                            //standardized normal values used for computing lenprob
  vector cprob_lenvec(1,nlenbins);              //cumulative probabilities used for computing lenprob
  number zscore_lzero;                          //standardized normal values for length = 0
  number cprob_lzero;                           //length probability mass below zero, used for computing lenprob
  
  //matrices below are used to match length comps
  matrix lenprob_sBT(1,nages,1,nlenbins);    //distn of size at age of sBT
  matrix lenprob_cHL1(1,nages,1,nlenbins);     //distn of size at age in cHL block 1
  matrix lenprob_cHL2(1,nages,1,nlenbins);     //distn of size at age in cHL block 2
  matrix lenprob_cHL3(1,nages,1,nlenbins);     //distn of size at age in cHL block 3
  matrix lenprob_cPT1(1,nages,1,nlenbins);     //distn of size at age in cPT block 1
  matrix lenprob_cPT2(1,nages,1,nlenbins);     //distn of size at age in cPT block 2
  matrix lenprob_cPT3(1,nages,1,nlenbins);     //distn of size at age in cPT block 3
  matrix lenprob_cTW1(1,nages,1,nlenbins);     //distn of size at age in cTW block 1
  matrix lenprob_cTW2(1,nages,1,nlenbins);     //distn of size at age in cTW block 2
  matrix lenprob_rHB1(1,nages,1,nlenbins);     //distn of size at age in rHB block 1
  matrix lenprob_rHB2(1,nages,1,nlenbins);     //distn of size at age in rHB block 2
  matrix lenprob_rHB3(1,nages,1,nlenbins);     //distn of size at age in rHB block 3
  matrix lenprob_rHB4(1,nages,1,nlenbins);     //distn of size at age in rHB block 4 
  matrix lenprob_rHB5(1,nages,1,nlenbins);     //distn of size at age in rHB block 5    
  matrix lenprob_rGN1(1,nages,1,nlenbins);   //distn of size at age in rGN block 2
  matrix lenprob_rGN2(1,nages,1,nlenbins);   //distn of size at age in rGN block 2
  matrix lenprob_rGN3(1,nages,1,nlenbins);   //distn of size at age in rGN block 3
  matrix lenprob_rGN4(1,nages,1,nlenbins);   //distn of size at age in rGN block 4
  matrix lenprob_rGN5(1,nages,1,nlenbins);   //distn of size at age in rGN block 5

  matrix lenprob_cGN_D2(1,nages,1,nlenbins);    //distn of size at age in cGN discards cGN block 2   
  matrix lenprob_cGN_D3(1,nages,1,nlenbins);    //distn of size at age in cGN discards cGN block 3
  matrix lenprob_cGN_D4(1,nages,1,nlenbins);    //distn of size at age in cGN discards cGN block 4
  matrix lenprob_rHB_D2(1,nages,1,nlenbins);      //distn of size at age in rHB discards rec block 2 
  matrix lenprob_rHB_D3(1,nages,1,nlenbins);      //distn of size at age in rHB discards rec block 3
  matrix lenprob_rHB_D4(1,nages,1,nlenbins);      //distn of size at age in rHB discards rec block 4
  matrix lenprob_rHB_D5(1,nages,1,nlenbins);      //distn of size at age in rHB discards rec block 5
  matrix lenprob_rGN_D1(1,nages,1,nlenbins);    //distn of size at age in rGN discards rec block 1
  matrix lenprob_rGN_D2(1,nages,1,nlenbins);    //distn of size at age in rGN discards rec block 2
  matrix lenprob_rGN_D3(1,nages,1,nlenbins);    //distn of size at age in rGN discards rec block 3
  matrix lenprob_rGN_D4(1,nages,1,nlenbins);    //distn of size at age in rGN discards rec block 4
  matrix lenprob_rGN_D5(1,nages,1,nlenbins);    //distn of size at age in rGN discards rec block 5

  //matrices below used to compute mean weights
  matrix lenprob_cHL1_all(1,nages,1,nlenbins_all);     //distn of size at age in cHL block 1
  matrix lenprob_cHL2_all(1,nages,1,nlenbins_all);     //distn of size at age in cHL block 2
  matrix lenprob_cHL3_all(1,nages,1,nlenbins_all);     //distn of size at age in cHL block 3
  matrix lenprob_cPT1_all(1,nages,1,nlenbins_all);     //distn of size at age in cPT block 1  
  matrix lenprob_cPT2_all(1,nages,1,nlenbins_all);     //distn of size at age in cPT block 2
  matrix lenprob_cPT3_all(1,nages,1,nlenbins_all);     //distn of size at age in cPT block 3
  matrix lenprob_cTW1_all(1,nages,1,nlenbins_all);     //distn of size at age in cTW block 1  
  matrix lenprob_cTW2_all(1,nages,1,nlenbins_all);     //distn of size at age in cTW block 2
  matrix lenprob_rHB1_all(1,nages,1,nlenbins_all);     //distn of size at age in rHB block 1
  matrix lenprob_rHB2_all(1,nages,1,nlenbins_all);     //distn of size at age in rHB block 2
  matrix lenprob_rHB3_all(1,nages,1,nlenbins_all);     //distn of size at age in rHB block 3
  matrix lenprob_rHB4_all(1,nages,1,nlenbins_all);     //distn of size at age in rHB block 4
  matrix lenprob_rHB5_all(1,nages,1,nlenbins_all);     //distn of size at age in rHB block 5
  matrix lenprob_rGN1_all(1,nages,1,nlenbins_all);   //distn of size at age in rGN block 1
  matrix lenprob_rGN2_all(1,nages,1,nlenbins_all);   //distn of size at age in rGN block 2
  matrix lenprob_rGN3_all(1,nages,1,nlenbins_all);   //distn of size at age in rGN block 3
  matrix lenprob_rGN4_all(1,nages,1,nlenbins_all);   //distn of size at age in rGN block 4
  matrix lenprob_rGN5_all(1,nages,1,nlenbins_all);   //distn of size at age in rGN block 5

  matrix lenprob_cGN_D2_all(1,nages,1,nlenbins_all);   //distn of size at age in cHL discards cGN block 2   
  matrix lenprob_cGN_D3_all(1,nages,1,nlenbins_all);   //distn of size at age in cHL discards cGN block 3
  matrix lenprob_cGN_D4_all(1,nages,1,nlenbins_all);   //distn of size at age in cHL discards cGN block 4
  matrix lenprob_rHB_D2_all(1,nages,1,nlenbins_all);     //distn of size at age in rHB discards rec block 2 
  matrix lenprob_rHB_D3_all(1,nages,1,nlenbins_all);     //distn of size at age in rHB discards rec block 3
  matrix lenprob_rHB_D4_all(1,nages,1,nlenbins_all);     //distn of size at age in rHB discards rec block 4
  matrix lenprob_rHB_D5_all(1,nages,1,nlenbins_all);     //distn of size at age in rHB discards rec block 5
  matrix lenprob_rGN_D1_all(1,nages,1,nlenbins_all);   //distn of size at age in rGN discards rec block 1
  matrix lenprob_rGN_D2_all(1,nages,1,nlenbins_all);   //distn of size at age in rGN discards rec block 2
  matrix lenprob_rGN_D3_all(1,nages,1,nlenbins_all);   //distn of size at age in rGN discards rec block 3
  matrix lenprob_rGN_D4_all(1,nages,1,nlenbins_all);   //distn of size at age in rGN discards rec block 4
  matrix lenprob_rGN_D5_all(1,nages,1,nlenbins_all);   //distn of size at age in rGN discards rec block 5

  vector len_sd(1,nages);
  vector len_cv(1,nages); //for fishgraph 
  
//----Predicted length and age compositions
  matrix pred_sBT_lenc(1,nyr_lenc_sBT,1,nlenbins);
  matrix pred_cHL_lenc(1,nyr_lenc_cHL,1,nlenbins);
  matrix pred_cPT_lenc(1,nyr_lenc_cPT,1,nlenbins);
  matrix pred_rHB_lenc(1,nyr_lenc_rHB,1,nlenbins);
  matrix pred_rHB_D_lenc(1,nyr_lenc_rHB_D,1,nlenbins);
  matrix pred_rGN_lenc(1,nyr_lenc_rGN,1,nlenbins);  
  
  matrix L_cPT_num_pool(1,nyr_lenc_cPT,1,nages);          //landings (numbers) at age pooled for length comps
  matrix L_cPT_num_pool_yr(1,nyr_lenc_pool_cPT,1,nages);  //scaled and weighted landings (numbers) for pooling length comps
  //matrix L_rGN_num_pool(1,nyr_lenc_rGN,1,nages);          //landings (numbers) at age pooled for length comps
  //matrix L_rGN_num_pool_yr(1,nyr_rGN_lenc_pool,1,nages);  //scaled and weighted landings (numbers) for pooling length comps
  
//  //##p_lenc_fishery pars require age comp and length comp data for estimation
  number p_lenc_cHL2;
  number p_lenc_cHL3;  
  number p_lenc_cPT2;
  number p_lenc_cPT3;  
  number p_lenc_cTW2;
  number p_lenc_cTW3;  
  number p_lenc_rHB2;
  number p_lenc_rHB3;
  number p_lenc_rHB4;
  number p_lenc_rHB5;
  number p_lenc_rGN2;
  number p_lenc_rGN3;
  number p_lenc_rGN4;
  number p_lenc_rGN5;
    
  number p_lenc_cGN_D2;
  number p_lenc_cGN_D3;
  number p_lenc_cGN_D4;
  number p_lenc_rHB_D2; 
  number p_lenc_rHB_D3;
  number p_lenc_rHB_D4;
  number p_lenc_rHB_D5;  
  number p_lenc_rGN_D1; 
  number p_lenc_rGN_D2; 
  number p_lenc_rGN_D3;
  number p_lenc_rGN_D4;
  number p_lenc_rGN_D5;
  
  matrix pred_sBT_agec(1,nyr_agec_sBT,1,nages);
  matrix ErrorFree_sBT_agec(1,nyr_agec_sBT,1,nages);
  matrix pred_sTV_agec(1,nyr_agec_sTV,1,nages);
  matrix ErrorFree_sTV_agec(1,nyr_agec_sTV,1,nages);  
  matrix pred_cHL_agec(1,nyr_agec_cHL,1,nages);
  matrix ErrorFree_cHL_agec(1,nyr_agec_cHL,1,nages);
  matrix pred_cPT_agec(1,nyr_agec_cPT,1,nages);
  matrix ErrorFree_cPT_agec(1,nyr_agec_cPT,1,nages);  
  matrix pred_rHB_agec(1,nyr_agec_rHB,1,nages);
  matrix ErrorFree_rHB_agec(1,nyr_agec_rHB,1,nages);
  //matrix pred_rGN_agec(1,nyr_rGN_agec,1,nages);
  //matrix ErrorFree_rGN_agec(1,nyr_rGN_agec,1,nages);
 
//effective sample size applied in multinomial distributions
  vector nsamp_sBT_lenc_allyr(styr,endyr);
  vector nsamp_cHL_lenc_allyr(styr,endyr);
  vector nsamp_cPT_lenc_allyr(styr,endyr);
  vector nsamp_rHB_lenc_allyr(styr,endyr);
  vector nsamp_rHB_D_lenc_allyr(styr,endyr);
  vector nsamp_rGN_lenc_allyr(styr,endyr);
  vector nsamp_sBT_agec_allyr(styr,endyr);
  vector nsamp_sTV_agec_allyr(styr,endyr);
  vector nsamp_cHL_agec_allyr(styr,endyr);
  vector nsamp_cPT_agec_allyr(styr,endyr);  
  vector nsamp_rHB_agec_allyr(styr,endyr);
  //vector nsamp_rGN_agec_allyr(styr,endyr);

//Nfish used in MCB analysis (not used in fitting)
  vector nfish_sBT_lenc_allyr(styr,endyr);
  vector nfish_cHL_lenc_allyr(styr,endyr);
  vector nfish_cPT_lenc_allyr(styr,endyr);
  vector nfish_rHB_lenc_allyr(styr,endyr);
  vector nfish_rHB_D_lenc_allyr(styr,endyr);
  vector nfish_rGN_lenc_allyr(styr,endyr);
  vector nfish_sBT_agec_allyr(styr,endyr);
  vector nfish_sTV_agec_allyr(styr,endyr);
  vector nfish_cHL_agec_allyr(styr,endyr);
  vector nfish_cPT_agec_allyr(styr,endyr);  
  vector nfish_rHB_agec_allyr(styr,endyr);
  //vector nfish_rGN_agec_allyr(styr,endyr);
  
//Computed effective sample size for output (not used in fitting)
  vector neff_sBT_lenc_allyr_out(styr,endyr);
  vector neff_cHL_lenc_allyr_out(styr,endyr);
  vector neff_cPT_lenc_allyr_out(styr,endyr);
  vector neff_rHB_lenc_allyr_out(styr,endyr);
  vector neff_rHB_D_lenc_allyr_out(styr,endyr);
  vector neff_rGN_lenc_allyr_out(styr,endyr);
  vector neff_sBT_agec_allyr_out(styr,endyr);
  vector neff_sTV_agec_allyr_out(styr,endyr);
  vector neff_cHL_agec_allyr_out(styr,endyr);
  vector neff_cPT_agec_allyr_out(styr,endyr);  
  vector neff_rHB_agec_allyr_out(styr,endyr);
  //vector neff_rGN_agec_allyr_out(styr,endyr);


//-----Population-----------------------------------------------------------------------------------
  matrix N(styr,endyr,1,nages);             //Population numbers by year and age at start of yr
  matrix N_mdyr(styr,endyr,1,nages);        //Population numbers by year and age at mdpt of yr: used for comps and cpue
  matrix N_spawn(styr,endyr,1,nages);       //Population numbers by year and age at peaking spawning: used for SSB  
  //init_bounded_vector log_dev_Nage(2,nages,-5,3,1); //log deviations on initial abundance at age
  init_bounded_vector log_dev_Nage(2,nages,log_Nage_dev_LO,log_Nage_dev_HI,log_Nage_dev_PH);
  //vector log_dev_Nage(2,nages);
  vector log_Nage_dev_output(1,nages);      //used in output. equals zero for first age
  matrix B(styr,endyr,1,nages);             //Population biomass by year and age at start of yr
  vector totB(styr,endyr);                  //Total biomass by year
  vector totN(styr,endyr);                  //Total abundance by year
  vector SSB(styr,endyr);                   //Total spawning biomass by year (scaled popn fecundity)
  vector MatFemB(styr,endyr);               //Total spawning biomass by year (total mature female biomass)  
  vector rec(styr,endyr);                   //Recruits by year
  vector prop_f(1,nages);                   //Proportion female by age
  vector maturity_f(1,nages);               //Proportion of female mature at age
  vector reprod(1,nages);                   //vector used to compute spawning biomass (scaled popn fecundity)
  vector reprod2(1,nages);                  //vector used to compute mature female biomass 

////---Stock-Recruit Function (Beverton-Holt, steepness parameterization)----------
  //init_bounded_number log_R0(13,20,1);        //log(virgin Recruitment)
  //number log_R0;
  init_bounded_number log_R0(log_R0_LO,log_R0_HI,log_R0_PH);        //log(virgin Recruitment)
  vector log_R0_out(1,8);
  number R0;                                  //virgin recruitment
  //init_bounded_number steep(0.21,0.991,3);    //steepness
  init_bounded_number steep(steep_LO,steep_HI,steep_PH); //steepness
  vector steep_out(1,8);
  //number steep;  //uncomment to fix steepness, comment line directly above
  //init_bounded_number rec_sigma(0.1,1.5,4);  //sd recruitment residuals
  init_bounded_number rec_sigma(rec_sigma_LO,rec_sigma_HI,rec_sigma_PH);  //sd recruitment residuals  
  vector rec_sigma_out(1,8);
  init_bounded_number R_autocorr(R_autocorr_LO,R_autocorr_HI,R_autocorr_PH);  //autocorrelation in SR  
  vector R_autocorr_out(1,8);
  
  number rec_sigma_sq;                        //square of rec_sigma      
  number rec_sigma_sqd2;                      //square of rec_sigma divided by two
  number rec_logL_add;                        //additive term in -logL term   
    
  //init_bounded_dev_vector log_dev_rec(styr_rec_dev,endyr,-3,3,2);  //log recruitment deviations
  //vector log_dev_rec(styr_rec_dev,endyr);
  init_bounded_dev_vector log_dev_rec(styr_rec_dev,endyr_rec_dev,log_rec_dev_LO,log_rec_dev_HI,log_rec_dev_PH);
  vector log_rec_dev_output(styr_rec_dev,endyr);             //used in output. equals zero except for yrs in log_dev_rec
  vector log_rec_dev_out(styr_rec_dev,endyr_rec_dev);  //used in output for bound checking

  number var_rec_dev;                                //variance of log recruitment deviations, from yrs with unconstrainted S-R(XXXX-XXXX)
  number sigma_rec_dev;                              //sample SD of log residuals (may not equal rec_sigma 
  number BiasCor;                               //Bias correction in equilibrium recruits
  //number R_autocorr;
  number S0;                                    //equal to spr_F0*R0 = virgin SSB
  number B0;                                    //equal to bpr_F0*R0 = virgin B  
  number R1;                                    //Recruits in styr
  number R_virgin;                              //unfished recruitment with bias correction
  vector SdS0(styr,endyr);                      //SSB / virgin SSB

  init_bounded_number log_dm_lenc_sBT(log_dm_sBT_lc_LO,log_dm_sBT_lc_HI,log_dm_sBT_lc_PH);
  init_bounded_number log_dm_lenc_sTV(log_dm_sTV_lc_LO,log_dm_sTV_lc_HI,log_dm_sTV_lc_PH);
  init_bounded_number log_dm_lenc_cHL(log_dm_cHL_lc_LO,log_dm_cHL_lc_HI,log_dm_cHL_lc_PH);
  init_bounded_number log_dm_lenc_cPT(log_dm_cPT_lc_LO,log_dm_cPT_lc_HI,log_dm_cPT_lc_PH);
  init_bounded_number log_dm_lenc_rHB(log_dm_rHB_lc_LO,log_dm_rHB_lc_HI,log_dm_rHB_lc_PH);
  init_bounded_number log_dm_lenc_rHB_D(log_dm_rHB_D_lc_LO,log_dm_rHB_D_lc_HI,log_dm_rHB_D_lc_PH);
  init_bounded_number log_dm_lenc_rGN(log_dm_rGN_lc_LO,log_dm_rGN_lc_HI,log_dm_rGN_lc_PH);
  init_bounded_number log_dm_agec_sBT(log_dm_sBT_ac_LO,log_dm_sBT_ac_HI,log_dm_sBT_ac_PH);
  init_bounded_number log_dm_agec_sTV(log_dm_sTV_ac_LO,log_dm_sTV_ac_HI,log_dm_sTV_ac_PH);
  init_bounded_number log_dm_agec_cHL(log_dm_cHL_ac_LO,log_dm_cHL_ac_HI,log_dm_cHL_ac_PH);
  init_bounded_number log_dm_agec_cPT(log_dm_cPT_ac_LO,log_dm_cPT_ac_HI,log_dm_cPT_ac_PH);
  init_bounded_number log_dm_agec_rHB(log_dm_rHB_ac_LO,log_dm_rHB_ac_HI,log_dm_rHB_ac_PH);
  //init_bounded_number log_dm_rGN_ac(log_dm_rGN_ac_LO,log_dm_rGN_ac_HI,log_dm_rGN_ac_PH);
  vector log_dm_sBT_lc_out(1,8);
  vector log_dm_sTV_lc_out(1,8);
  vector log_dm_cHL_lc_out(1,8);
  vector log_dm_cPT_lc_out(1,8);
  vector log_dm_rHB_lc_out(1,8);
  vector log_dm_rGN_lc_out(1,8);
  vector log_dm_rHB_D_lc_out(1,8);
  vector log_dm_sBT_ac_out(1,8);
  vector log_dm_sTV_ac_out(1,8);
  vector log_dm_cHL_ac_out(1,8);
  vector log_dm_cPT_ac_out(1,8);
  vector log_dm_rHB_ac_out(1,8);
  //vector log_dm_rGN_ac_out(1,8);
  
//-----------------------------------------------------------------------------------------------------------------------------------------------
//---Selectivity-------------------------------------------------------------------------

//MARMAP sBT selectivity -------------------------------------------------------------------------
  matrix sel_sBT(styr,endyr,1,nages);
  vector sel_sBT_vec(1,nages);
  
  init_bounded_number selpar_A50_sBT(selpar_A50_sBT_LO,selpar_A50_sBT_HI,selpar_A50_sBT_PH);
  init_bounded_number selpar_slope_sBT(selpar_slope_sBT_LO,selpar_slope_sBT_HI,selpar_slope_sBT_PH);
   
  vector selpar_A50_sBT_out(1,8);
  vector selpar_slope_sBT_out(1,8);
  
//MARMAP sTV selectivity -------------------------------------------------------------------------
  matrix sel_sTV(styr,endyr,1,nages);
  vector sel_sTV_vec(1,nages);
  
  init_bounded_number selpar_A50_sTV(selpar_A50_sTV_LO,selpar_A50_sTV_HI,selpar_A50_sTV_PH);
  init_bounded_number selpar_slope_sTV(selpar_slope_sTV_LO,selpar_slope_sTV_HI,selpar_slope_sTV_PH);
   
  vector selpar_A50_sTV_out(1,8);
  vector selpar_slope_sTV_out(1,8);
  
//Commercial handline selectivity-------------------------------------------------
  matrix sel_cHL(styr,endyr,1,nages);  
  //vector sel_cHL_1(1,nages); //sel in selex_phase1 assumed equal to selex_phase2
  vector sel_cHL_2(1,nages); //sel in selex_phase2
  vector sel_cHL_3(1,nages); //sel in selex_phase3 
  vector sel_cHL_4(1,nages); //sel in selex_phase4 
  
  init_bounded_number selpar_A50_cHL2(selpar_A50_cHL2_LO,selpar_A50_cHL2_HI,selpar_A50_cHL2_PH);
  init_bounded_number selpar_slope_cHL2(selpar_slope_cHL2_LO,selpar_slope_cHL2_HI,selpar_slope_cHL2_PH);
  init_bounded_number selpar_A50_cHL3(selpar_A50_cHL3_LO,selpar_A50_cHL3_HI,selpar_A50_cHL3_PH);
  init_bounded_number selpar_slope_cHL3(selpar_slope_cHL3_LO,selpar_slope_cHL3_HI,selpar_slope_cHL3_PH);
  init_bounded_number selpar_A50_cHL4(selpar_A50_cHL4_LO,selpar_A50_cHL4_HI,selpar_A50_cHL4_PH);
  init_bounded_number selpar_slope_cHL4(selpar_slope_cHL4_LO,selpar_slope_cHL4_HI,selpar_slope_cHL4_PH);
 
  vector selpar_A50_cHL2_out(1,8);
  vector selpar_slope_cHL2_out(1,8);
  vector selpar_A50_cHL3_out(1,8);
  vector selpar_slope_cHL3_out(1,8);
  vector selpar_A50_cHL4_out(1,8);
  vector selpar_slope_cHL4_out(1,8);
     
//commercial discards (handline + pots)
  matrix sel_cGN_D(styr,endyr,1,nages); 
  vector sel_cGN_D_2(1,nages);         //sel in selex_phase2
  vector sel_cGN_D_3(1,nages);         //sel in selex_phase3 
  vector sel_cGN_D_4(1,nages);         //sel in selex_phase4  
  vector sel_cGN_D_quota3(1,nages);    //sel in selex_phase3 when quotas were in place (2009,2010) Also in 2011-2012   
 
//values used for weighting selex and avg weights of discards during yrs with quotas
  number Dopen_cHL; number Dclosed_cHL; number Lopen_cHL;  
  number Dopen_cPT; number Dclosed_cPT; number Lopen_cPT;
  number D_sum_cHLcPT; 
  number Dprop_cGN_sel_D; number Dprop_cGN_sel_cHL; number Dprop_cGN_sel_cPT; 
         
//Commercial pots selectivity -------------------------------------------------            
  matrix sel_cPT(styr,endyr,1,nages); 
  //  vector sel_cPT_1(1,nages); //sel vector in selex_phase1 assumed equal to selex_phase2 
  vector sel_cPT_2(1,nages);   //sel vector in selex_phase2  
  vector sel_cPT_3(1,nages);   //sel vector in selex_phase3  
  vector sel_cPT_4(1,nages);   //sel vector in selex_phase3    
  
  init_bounded_number selpar_A50_cPT2(selpar_A50_cPT2_LO,selpar_A50_cPT2_HI,selpar_A50_cPT2_PH);
  init_bounded_number selpar_slope_cPT2(selpar_slope_cPT2_LO,selpar_slope_cPT2_HI,selpar_slope_cPT2_PH);
  init_bounded_number selpar_A50_cPT3(selpar_A50_cPT3_LO,selpar_A50_cPT3_HI,selpar_A50_cPT3_PH);
  init_bounded_number selpar_slope_cPT3(selpar_slope_cPT3_LO,selpar_slope_cPT3_HI,selpar_slope_cPT3_PH);
  init_bounded_number selpar_A50_cPT4(selpar_A50_cPT4_LO,selpar_A50_cPT4_HI,selpar_A50_cPT4_PH);
  init_bounded_number selpar_slope_cPT4(selpar_slope_cPT4_LO,selpar_slope_cPT4_HI,selpar_slope_cPT4_PH);
 
  vector selpar_A50_cPT2_out(1,8);
  vector selpar_slope_cPT2_out(1,8);
  vector selpar_A50_cPT3_out(1,8);
  vector selpar_slope_cPT3_out(1,8);
  vector selpar_A50_cPT4_out(1,8);
  vector selpar_slope_cPT4_out(1,8);

//Commercial trawl selectivity -------------------------------------------------            
  matrix sel_cTW(styr,endyr,1,nages);  //mirrors cGN pot sel

//Headboat selectivity -------------------------------------------------
  matrix sel_rHB(styr,endyr,1,nages);  
  vector sel_rHB_1(1,nages); //sel in selex_phase1 
  vector sel_rHB_2(1,nages); //sel in selex_phase2
  vector sel_rHB_3(1,nages); //sel in selex_phase3 
  vector sel_rHB_4(1,nages); //sel in selex_phase4  
  vector sel_rHB_5(1,nages); //sel in selex_phase5   
  
  init_bounded_number selpar_A50_rHB1(selpar_A50_rHB1_LO,selpar_A50_rHB1_HI,selpar_A50_rHB1_PH);
  init_bounded_number selpar_slope_rHB1(selpar_slope_rHB1_LO,selpar_slope_rHB1_HI,selpar_slope_rHB1_PH);
  init_bounded_number selpar_A50_rHB2(selpar_A50_rHB2_LO,selpar_A50_rHB2_HI,selpar_A50_rHB2_PH);
  init_bounded_number selpar_slope_rHB2(selpar_slope_rHB2_LO,selpar_slope_rHB2_HI,selpar_slope_rHB2_PH);    
  init_bounded_number selpar_A50_rHB3(selpar_A50_rHB3_LO,selpar_A50_rHB3_HI,selpar_A50_rHB3_PH);
  init_bounded_number selpar_slope_rHB3(selpar_slope_rHB3_LO,selpar_slope_rHB3_HI,selpar_slope_rHB3_PH);
  init_bounded_number selpar_A50_rHB4(selpar_A50_rHB4_LO,selpar_A50_rHB4_HI,selpar_A50_rHB4_PH);
  init_bounded_number selpar_slope_rHB4(selpar_slope_rHB4_LO,selpar_slope_rHB4_HI,selpar_slope_rHB4_PH);
  init_bounded_number selpar_A50_rHB5(selpar_A50_rHB5_LO,selpar_A50_rHB5_HI,selpar_A50_rHB5_PH);
  init_bounded_number selpar_slope_rHB5(selpar_slope_rHB5_LO,selpar_slope_rHB5_HI,selpar_slope_rHB5_PH);
  
  vector selpar_A50_rHB1_out(1,8);
  vector selpar_slope_rHB1_out(1,8);
  vector selpar_A50_rHB2_out(1,8);
  vector selpar_slope_rHB2_out(1,8);
  vector selpar_A50_rHB3_out(1,8);
  vector selpar_slope_rHB3_out(1,8);
  vector selpar_A50_rHB4_out(1,8);
  vector selpar_slope_rHB4_out(1,8);
  vector selpar_A50_rHB5_out(1,8);
  vector selpar_slope_rHB5_out(1,8);
  
  number selpar_Age0_rHB_D;
  number selpar_Age1_rHB_D;
  number selpar_Age2_rHB_D;
  
  //---headboat discards--------------------------------- 
  matrix sel_rHB_D(styr,endyr,1,nages); 
  vector sel_rHB_D_1(1,nages); //sel in selex_phase1, assumed equal to selex_phase2
  vector sel_rHB_D_2(1,nages); //sel in selex_phase2
  vector sel_rHB_D_3(1,nages); //sel in selex_phase3 
  vector sel_rHB_D_4(1,nages); //sel in selex_phase4  
  vector sel_rHB_D_5(1,nages); //sel in selex_phase5 
  
  init_bounded_number selpar_A50_rHD4(selpar_A50_rHD4_LO,selpar_A50_rHD4_HI,selpar_A50_rHD4_PH);
  init_bounded_number selpar_slope_rHD4(selpar_slope_rHD4_LO,selpar_slope_rHD4_HI,selpar_slope_rHD4_PH);
  init_bounded_number selpar_A502_rHD4(selpar_A502_rHD4_LO,selpar_A502_rHD4_HI,selpar_A502_rHD4_PH);
  init_bounded_number selpar_slope2_rHD4(selpar_slope2_rHD4_LO,selpar_slope2_rHD4_HI,selpar_slope2_rHD4_PH);
  
  init_bounded_number selpar_A50_rHD5(selpar_A50_rHD5_LO,selpar_A50_rHD5_HI,selpar_A50_rHD5_PH);
  init_bounded_number selpar_slope_rHD5(selpar_slope_rHD5_LO,selpar_slope_rHD5_HI,selpar_slope_rHD5_PH);
  init_bounded_number selpar_A502_rHD5(selpar_A502_rHD5_LO,selpar_A502_rHD5_HI,selpar_A502_rHD5_PH);
  init_bounded_number selpar_slope2_rHD5(selpar_slope2_rHD5_LO,selpar_slope2_rHD5_HI,selpar_slope2_rHD5_PH);
  
  vector selpar_A50_rHD4_out(1,8);
  vector selpar_slope_rHD4_out(1,8);
  vector selpar_A502_rHD4_out(1,8);
  vector selpar_slope2_rHD4_out(1,8);
  vector selpar_A50_rHD5_out(1,8);
  vector selpar_slope_rHD5_out(1,8);
  vector selpar_A502_rHD5_out(1,8);
  vector selpar_slope2_rHD5_out(1,8);
  
  vector vecprob_rHB_D2(4,nages);     //prob of less than size limit
  vector vecprob_rHB_D3(4,nages);     //prob of less than size limit
  vector vecprob_rHB_D4(4,nages);     //prob of less than size limit
  vector vecprob_rHB_D5(4,nages);     //prob of less than size limit

  init_bounded_number selpar_logit_Age0_rHB_D(selpar_Age0_rHB_D_logit_LO,selpar_Age0_rHB_D_logit_HI,selpar_Age0_rHB_D_logit_PH);
  init_bounded_number selpar_logit_Age1_rHB_D(selpar_Age1_rHB_D_logit_LO,selpar_Age1_rHB_D_logit_HI,selpar_Age1_rHB_D_logit_PH);
  init_bounded_number selpar_logit_Age2_rHB_D(selpar_Age2_rHB_D_logit_LO,selpar_Age2_rHB_D_logit_HI,selpar_Age2_rHB_D_logit_PH);
  
  vector selpar_Age0_rHB_D_logit_out(1,8);                                  
  vector selpar_Age1_rHB_D_logit_out(1,8);                                  
  vector selpar_Age2_rHB_D_logit_out(1,8); 
  
  vector prob_belowsizelim_block1(1,nages);
  vector prob_belowsizelim_block2(1,nages);
  vector prob_belowsizelim_block3(1,nages);
  vector prob_belowsizelim_block4(1,nages);
  vector prob_belowsizelim_block5(1,nages);
  
  number zscore_lsizelim1;
  number zscore_lsizelim2;
  number zscore_lsizelim3;
  number zscore_lsizelim4;
  number zscore_lsizelim5;
  
  number cprob_lsizelim1;
  number cprob_lsizelim2;
  number cprob_lsizelim3;
  number cprob_lsizelim4;
  number cprob_lsizelim5;
  
    //rGN selectivity  -------------------------------------------------
  matrix sel_rGN(styr,endyr,1,nages);  
  matrix sel_rGN_D(styr,endyr,1,nages); 
  vector sel_rGN1(1,nages); //sel in selex_phase1 
  vector sel_rGN2(1,nages); //sel in selex_phase2
  vector sel_rGN3(1,nages); //sel in selex_phase3 
  vector sel_rGN4(1,nages); //sel in selex_phase4  
  vector sel_rGN5(1,nages); //sel in selex_phase5   
     
  init_bounded_number selpar_A50_rGN1(selpar_A50_rGN1_LO,selpar_A50_rGN1_HI,selpar_A50_rGN1_PH);
  init_bounded_number selpar_slope_rGN1(selpar_slope_rGN1_LO,selpar_slope_rGN1_HI,selpar_slope_rGN1_PH);
  init_bounded_number selpar_A50_rGN2(selpar_A50_rGN2_LO,selpar_A50_rGN2_HI,selpar_A50_rGN2_PH);
  init_bounded_number selpar_slope_rGN2(selpar_slope_rGN2_LO,selpar_slope_rGN2_HI,selpar_slope_rGN2_PH);    
  init_bounded_number selpar_A50_rGN3(selpar_A50_rGN3_LO,selpar_A50_rGN3_HI,selpar_A50_rGN3_PH);
  init_bounded_number selpar_slope_rGN3(selpar_slope_rGN3_LO,selpar_slope_rGN3_HI,selpar_slope_rGN3_PH);
  init_bounded_number selpar_A50_rGN4(selpar_A50_rGN4_LO,selpar_A50_rGN4_HI,selpar_A50_rGN4_PH);
  init_bounded_number selpar_slope_rGN4(selpar_slope_rGN4_LO,selpar_slope_rGN4_HI,selpar_slope_rGN4_PH);
  init_bounded_number selpar_A50_rGN5(selpar_A50_rGN5_LO,selpar_A50_rGN5_HI,selpar_A50_rGN5_PH);
  init_bounded_number selpar_slope_rGN5(selpar_slope_rGN5_LO,selpar_slope_rGN5_HI,selpar_slope_rGN5_PH);
  
  vector selpar_A50_rGN1_out(1,8);
  vector selpar_slope_rGN1_out(1,8);
  vector selpar_A50_rGN2_out(1,8);
  vector selpar_slope_rGN2_out(1,8);
  vector selpar_A50_rGN3_out(1,8);
  vector selpar_slope_rGN3_out(1,8);
  vector selpar_A50_rGN4_out(1,8);
  vector selpar_slope_rGN4_out(1,8);
  vector selpar_A50_rGN5_out(1,8);
  vector selpar_slope_rGN5_out(1,8);
  
//effort-weighted, recent selectivities
  vector sel_wgted_L(1,nages);  //toward landings
  vector sel_wgted_D(1,nages);  //toward discards  
  vector sel_wgted_tot(1,nages);//toward Z, landings plus deads discards

//-----------------------------------------------------------------------------------------------------------------------------------------------
//-------CPUE Predictions--------------------------------
  vector pred_sBT_cpue(styr_cpue_sBT,endyr_cpue_sBT);       //predicted sBT U (fish/trap-hour)
  matrix N_sBT(styr_cpue_sBT,endyr_cpue_sBT,1,nages);       //used to compute sBT index
  vector pred_sTV_cpue(styr_cpue_sTV,endyr_cpue_sTV);       //predicted sTV U (fish/trap-hour)
  matrix N_sTV(styr_cpue_sTV,endyr_cpue_sTV,1,nages);       //used to compute sTV index
  //vector pred_Vid_cpue(styr_Vid_cpue,endyr_Vid_cpue);       //predicted sTV U (fish/trap-hour)
  //matrix N_Vid(styr_Vid_cpue,endyr_Vid_cpue,1,nages);       //used to compute sTV index
  vector pred_cHL_cpue(styr_cpue_cHL,endyr_cpue_cHL);             //predicted cHL U (pounds/hook-hour)
  matrix N_cHL(styr_cpue_cHL,endyr_cpue_cHL,1,nages);             //used to compute cHL index
  vector pred_rHB_cpue(styr_cpue_rHB,endyr_cpue_rHB);             //predicted rHB U (pounds/hour)
  matrix N_rHB(styr_cpue_rHB,endyr_cpue_rHB,1,nages);             //used to compute rHB index
  //vector pred_rHD_cpue(styr_rHD_cpue,endyr_rHD_cpue);          //predicted rHD U (fish/angler-hour)
  //matrix N_rHD(styr_rHD_cpue,endyr_rHD_cpue,1,nages);          //used to compute rHD index

//---Catchability (CPUE q's)----------------------------------------------------------
  //init_bounded_number log_q_sBT(-20,-10,1);
  //init_bounded_number log_q_sTV(-20,-10,1);
  //init_bounded_number log_q_Vid(-20,-10,1);  
  //init_bounded_number log_q_cHL(-20,-5,1);
  //init_bounded_number log_q_rHB(-20,-5,1);
  //init_bounded_number log_q_rHD(-20,-10,1);
  init_bounded_number q_rate(0.001,0.1,set_q_rate_phase);
  
  init_bounded_number log_q_sBT(log_q_sBT_LO,log_q_sBT_HI,log_q_sBT_PH);
  init_bounded_number log_q_sTV(log_q_sTV_LO,log_q_sTV_HI,log_q_sTV_PH);
  init_bounded_number log_q_cHL(log_q_cHL_LO,log_q_cHL_HI,log_q_cHL_PH);
  init_bounded_number log_q_rHB(log_q_rHB_LO,log_q_rHB_HI,log_q_rHB_PH);
  //init_bounded_number log_q_rHD(log_q_rHD_LO,log_q_rHD_HI,log_q_rHD_PH);
  //init_bounded_number log_q_Vid(log_q_Vid_LO,log_q_Vid_HI,log_q_Vid_PH);
  
  vector log_q_sBT_out(1,8);
  vector log_q_sTV_out(1,8);
  vector log_q_cHL_out(1,8);
  vector log_q_rHB_out(1,8);
  //vector log_q_rHD_out(1,8);
  //vector log_q_Vid_out(1,8);
  
  //number q_rate;
  vector q_rate_fcn_cHL(styr_cpue_cHL,endyr_cpue_cHL);       //increase due to technology creep (saturates in 2003)
  vector q_rate_fcn_rHB(styr_cpue_rHB,endyr_cpue_rHB);         //increase due to technology creep (saturates in 2003)
  //vector q_rate_fcn_rHD(styr_rHD_cpue,endyr_rHD_cpue);      //increase due to technology creep (saturates in 2003)
  
  //init_bounded_number q_DD_beta(0.1,0.9,set_q_DD_phase);  
  number q_DD_beta;
  vector q_DD_fcn(styr,endyr);    //density dependent function as a multiple of q (scaled a la Katsukawa and Matsuda. 2003)
  number B0_q_DD;                 //B0 of ages q_DD_age plus
  vector B_q_DD(styr,endyr);      //annual biomass of ages q_DD_age plus

  //init_bounded_vector q_RW_log_dev_cHL(styr_cpue_cHL,endyr_cpue_cHL-1,-3.0,3.0,set_q_RW_phase);
  //init_bounded_vector q_RW_log_dev_rHB(styr_cpue_rHB,endyr_cpue_rHB-1,-3.0,3.0,set_q_RW_phase);
  //init_bounded_vector q_RW_log_dev_rHD(styr_rHD_cpue,endyr_rHD_cpue-1,-3.0,3.0,set_q_RW_phase);
 
  //init_bounded_vector q_RW_log_dev_cHL(styr_cpue_cHL,endyr_cpue_cHL-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
  //init_bounded_vector q_RW_log_dev_rHB(styr_cpue_rHB,endyr_cpue_rHB-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
  //init_bounded_vector q_RW_log_dev_rHD(styr_rHD_cpue,endyr_rHD_cpue-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
  init_bounded_vector q_RW_log_dev_cHL(styr_cpue_cHL,endyr_cpue_cHL-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
  init_bounded_vector q_RW_log_dev_rHB(styr_cpue_rHB,endyr_cpue_rHB-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
  //init_bounded_vector q_RW_log_dev_rHD(styr_rHD_cpue,endyr_rHD_cpue-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
 
  vector q_cHL(styr_cpue_cHL,endyr_cpue_cHL);
  vector q_rHB(styr_cpue_rHB,endyr_cpue_rHB);
  //vector q_rHD(styr_rHD_cpue,endyr_rHD_cpue);

//----------------------------------------------------------------------------------------------------------------------------------------------- 
//---Landings Bias for recreational landings------------------------------------------------------------------
  //init_bounded_number L_rGN_bias(0.1,10.0,3);  
  number L_rHB_bias;
  number L_rGN_bias;
  number L_cGN_bias;

//---Landings in numbers (total or 1000 fish) and in wgt (klb)--------------------------------------------------
  matrix L_cHL_num(styr,endyr,1,nages);  //landings (numbers) at age
  matrix L_cHL_klb(styr,endyr,1,nages);  //landings (1000 lb whole weight) at age
  vector pred_cHL_L_knum(styr,endyr);    //yearly landings in 1000 fish summed over ages
  vector pred_cHL_L_klb(styr,endyr);     //yearly landings in 1000 lb summed over ages

  matrix L_cPT_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_cPT_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age    
  vector pred_cPT_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages   
  vector pred_cPT_L_klb(styr,endyr);      //yearly landings in 1000 lb summed over ages 

  matrix L_cTW_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_cTW_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age    
  vector pred_cTW_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages   
  vector pred_cTW_L_klb(styr,endyr);      //yearly landings in 1000 lb summed over ages 

  matrix L_rHB_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_rHB_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age    
  vector pred_rHB_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages  
  vector pred_rHB_L_klb(styr,endyr);      //yearly landings in 1000 lb summed over ages
  vector obs_rHB_L_wbias(styr,endyr);     //yearly landings observed, perhaps adjusted for multiplicative bias

  matrix L_rGN_num(styr,endyr,1,nages);  //landings (numbers) at age
  matrix L_rGN_klb(styr,endyr,1,nages);  //landings (1000 lb whole weight) at age    
  vector pred_rGN_L_knum(styr,endyr);    //yearly landings in 1000 fish summed over ages  
  vector pred_rGN_L_klb(styr,endyr);     //yearly landings in 1000 lb summed over ages
  vector obs_rGN_L_wbias(styr,endyr);    //yearly landings observed, perhaps adjusted for multiplicative bias

  matrix L_total_num(styr,endyr,1,nages);//total landings in number at age
  matrix L_total_klb(styr,endyr,1,nages);//landings in klb at age 
  vector L_total_knum_yr(styr,endyr);    //total landings in 1000 fish by yr summed over ages  
  vector L_total_klb_yr(styr,endyr);     //total landings (klb) by yr summed over ages

//---Dead discards in numbers (total or 1000 fish) and in wgt (klb) --------------------------------------------------
  matrix D_cGN_num(styr,endyr,1,nages);   //discards (numbers) at age
  matrix D_cGN_klb(styr,endyr,1,nages);   //discards (1000 lb) at age
  vector pred_cGN_D_knum(styr,endyr);     //yearly discards summed over ages
  vector pred_cGN_D_klb(styr,endyr);      //yearly discards in klb summed over ages
  vector obs_cGN_D(styr_D_cHL,endyr_D_cHL); //observed releases multiplied by discard mortality
  vector obs_cHL_D(styr_D_cHL,endyr_D_cHL);   //observed releases multiplied by discard mortality
  vector obs_cPT_D(styr_D_cHL,endyr_D_cHL);   //observed releases multiplied by discard mortality
  vector cGN_D_cv(styr_D_cHL,endyr_D_cHL);  //CVs for fitting combined discards
  
  matrix D_rHB_num(styr,endyr,1,nages);     //discards (numbers) at age
  matrix D_rHB_klb(styr,endyr,1,nages);     //discards (1000 lb) at age
  vector pred_rHB_D_knum(styr,endyr);       //yearly discards summed over ages
  vector pred_rHB_D_klb(styr,endyr);        //yearly discards in klb summed over ages
  vector obs_rHB_D(styr_D_rHB,endyr_D_rHB);   //observed releases multiplied by discard mortality

  matrix D_rGN_num(styr,endyr,1,nages);   //discards (numbers) at age
  matrix D_rGN_klb(styr,endyr,1,nages);   //discards (1000 lb) at age
  vector pred_rGN_D_knum(styr,endyr);     //yearly discards summed over ages
  vector pred_rGN_D_klb(styr,endyr);      //yearly discards in klb summed over ages
  vector obs_rGN_D(styr_D_rGN,endyr_D_rGN); //observed releases multiplied by discard mortality
  
  matrix D_total_num(styr,endyr,1,nages);  //total discards in number at age
  matrix D_total_klb(styr,endyr,1,nages);  //discards in klb at age 
  vector D_total_knum_yr(styr,endyr);      //total discards in 1000 fish by yr summed over ages  
  vector D_total_klb_yr(styr,endyr);       //total discards (klb) by yr summed over ages

////---MSY calcs----------------------------------------------------------------------------
  number F_cHL_prop;       //proportion of F_sum attributable to hal, last X=selpar_n_yrs_wgted yrs, used for avg body weights
  number F_cPT_prop;       //proportion of F_sum attributable to pots, last X yrs
  number F_rHB_prop;       //proportion of F_sum attributable to headboat, last X yrs
  number F_rGN_prop;     //proportion of F_sum attributable to rGN, last X yrs
  number F_cGN_D_prop;   //proportion of F_sum attributable to cGN discards, last X yrs
  number F_rHB_D_prop;     //proportion of F_sum attributable to headboat discards, last X yrs
  number F_rGN_D_prop;   //proportion of F_sum attributable to rGN discards, last X yrs
  number F_temp_sum;      //sum of geom mean Fsum's in last X yrs, used to compute F_fishery_prop

  vector F_end(1,nages);
  vector F_end_L(1,nages);
  vector F_end_D(1,nages);    
  number F_end_apex;
  
  number SSB_msy_out;           //SSB (popn fecudity) at msy
  number F_msy_out;             //F at msy
  number msy_klb_out;           //max sustainable yield (1000 lb)
  number msy_knum_out;          //max sustainable yield (1000 fish)  
  number B_msy_out;             //total biomass at MSY 
  number R_msy_out;             //equilibrium recruitment at F=Fmsy
  number D_msy_knum_out;        //equilibrium dead discards (1000 fish) at F=Fmsy
  number D_msy_klb_out;         //equilibrium dead discards (1000 lb) at F=Fmsy  
  number spr_msy_out;           //spr at F=Fmsy

  vector N_age_msy(1,nages);         //numbers at age for MSY calculations: beginning of yr
  vector N_age_msy_spawn(1,nages);   //numbers at age for MSY calculations: time of peak spawning
  vector L_age_msy(1,nages);         //catch at age for MSY calculations
  vector Z_age_msy(1,nages);         //total mortality at age for MSY calculations
  vector D_age_msy(1,nages);         //discard mortality (dead discards) at age for MSY calculations
  vector F_L_age_msy(1,nages);       //fishing mortality landings (not discards) at age for MSY calculations
  vector F_D_age_msy(1,nages);       //fishing mortality of discards at age for MSY calculations
  vector F_msy(1,n_iter_msy);        //values of full F to be used in equilibrium calculations
  vector spr_msy(1,n_iter_msy);      //reproductive capacity-per-recruit values corresponding to F values in F_msy
  vector R_eq(1,n_iter_msy);         //equilibrium recruitment values corresponding to F values in F_msy
  vector L_eq_klb(1,n_iter_msy);     //equilibrium landings(klb) values corresponding to F values in F_msy
  vector L_eq_knum(1,n_iter_msy);    //equilibrium landings(1000 fish) values corresponding to F values in F_msy
  vector SSB_eq(1,n_iter_msy);       //equilibrium reproductive capacity values corresponding to F values in F_msy
  vector B_eq(1,n_iter_msy);         //equilibrium biomass values corresponding to F values in F_msy
  vector D_eq_klb(1,n_iter_msy);     //equilibrium discards (klb) corresponding to F values in F_msy
  vector D_eq_knum(1,n_iter_msy);    //equilibrium discards (1000s) corresponding to F values in F_msy

  vector FdF_msy(styr,endyr);
  vector SdSSB_msy(styr,endyr);
  number SdSSB_msy_end;
  number FdF_msy_end;
  number FdF_msy_end_mean;          //geometric mean of last 3 yrs  

  vector wgt_wgted_L_klb(1,nages);  //fishery-weighted average weight at age of landings
  vector wgt_wgted_D_klb(1,nages);  //fishery-weighted average weight at age of discards  
  number wgt_wgted_L_denom;         //used in intermediate calculations
  number wgt_wgted_D_denom;         //used in intermediate calculations
  
  number iter_inc_msy;               //increments used to compute msy, equals 1/(n_iter_msy-1)
  
////--------Mortality------------------------------------------------------------------

// Stuff immediately below used only if M is estimated
//  //init_bounded_number M_constant(0.1,0.2,1);  //age-indpendent: used only for MSST
//  vector Mscale_ages(1,max_obs_age);
//  vector Mscale_len(1,max_obs_age);
//  vector Mscale_wgt_g(1,max_obs_age); 
//  vector M_lorenzen(1,max_obs_age);   
//  number cum_surv_1plus;  

  vector M(1,nages);                         //age-dependent natural mortality
  //number M_constant;                         //age-indpendent: used only for MSST
  init_bounded_number M_constant(M_constant_LO,M_constant_HI,M_constant_PH);   //age-indpendent: used only for MSST
  vector M_constant_out(1,8);
 
  number smsy2msstM;                         //scales Smsy to get msst using (1-M). Used only in output.
  number smsy2msst75;                        //scales Smsy to get msst using 75%. Used only in output.  
  
  matrix F(styr,endyr,1,nages);
  vector Fsum(styr,endyr);                   //Full fishing mortality rate by year
  vector Fapex(styr,endyr);                  //Max across ages, fishing mortality rate by year (may differ from Fsum bc of dome-shaped sel 
//  sdreport_vector fullF_sd(styr,endyr);
  matrix Z(styr,endyr,1,nages);

  //init_bounded_number log_avg_F_L_cHL(-10,0.0,1); 
  //init_bounded_dev_vector log_dev_F_L_cHL(styr_L_cHL,endyr_L_cHL,-10.0,5.0,2);  
  
  init_bounded_number log_avg_F_L_cHL(log_avg_F_cHL_LO,log_avg_F_cHL_HI,log_avg_F_cHL_PH);
  vector log_avg_F_cHL_out(1,8);  
  init_bounded_dev_vector log_dev_F_L_cHL(styr_L_cHL,endyr_L_cHL,log_F_dev_cHL_LO,log_F_dev_cHL_HI,log_F_dev_cHL_PH);
  vector log_F_dev_cHL_out(styr_L_cHL,endyr_L_cHL);
  matrix F_cHL(styr,endyr,1,nages);
  vector F_cHL_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_cHL;
  number log_F_dev_end_cHL;  

  //init_bounded_number log_avg_F_L_cPT(-10,0.0,1);
  //init_bounded_dev_vector log_dev_F_L_cPT(styr_L_cPT,endyr_L_cPT,-10.0,5.0,2);
  init_bounded_number log_avg_F_L_cPT(log_avg_F_cPT_LO,log_avg_F_cPT_HI,log_avg_F_cPT_PH);
  vector log_avg_F_cPT_out(1,8);  
  init_bounded_dev_vector log_dev_F_L_cPT(styr_L_cPT,endyr_L_cPT,log_F_dev_cPT_LO,log_F_dev_cPT_HI,log_F_dev_cPT_PH);
  vector log_F_dev_cPT_out(styr_L_cPT,endyr_L_cPT);
  matrix F_cPT(styr,endyr,1,nages);
  vector F_cPT_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_cPT;
  number log_F_dev_end_cPT;  

  //init_bounded_number log_avg_F_L_cTW(-10,0.0,1);
  //init_bounded_dev_vector log_dev_F_L_cTW(styr_L_cTW,endyr_L_cTW,-10.0,5.0,2);
  init_bounded_number log_avg_F_L_cTW(log_avg_F_cTW_LO,log_avg_F_cTW_HI,log_avg_F_cTW_PH);
  vector log_avg_F_cTW_out(1,8);  
  init_bounded_dev_vector log_dev_F_L_cTW(styr_L_cTW,endyr_L_cTW,log_F_dev_cTW_LO,log_F_dev_cTW_HI,log_F_dev_cTW_PH);
  vector log_F_dev_cTW_out(styr_L_cTW,endyr_L_cTW);
  matrix F_cTW(styr,endyr,1,nages);
  vector F_cTW_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_cTW;
  number log_F_dev_end_cTW;  
 
  //init_bounded_number log_avg_F_L_rHB(-10.0,0.0,1);  
  //init_bounded_dev_vector log_dev_F_L_rHB(styr_L_rHB,endyr_L_rHB,-10.0,5.0,2);
  init_bounded_number log_avg_F_L_rHB(log_avg_F_rHB_LO,log_avg_F_rHB_HI,log_avg_F_rHB_PH);
  vector log_avg_F_rHB_out(1,8); 
  init_bounded_dev_vector log_dev_F_L_rHB(styr_L_rHB,endyr_L_rHB,log_F_dev_rHB_LO,log_F_dev_rHB_HI,log_F_dev_rHB_PH);    
  vector log_F_dev_rHB_out(styr_L_rHB,endyr_L_rHB);
  matrix F_rHB(styr,endyr,1,nages);
  vector F_rHB_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_init_rHB;
  number log_F_dev_end_rHB;    

  //init_bounded_number log_avg_F_L_rGN(-10,0.0,1);
  //init_bounded_dev_vector log_dev_F_L_rGN(styr_L_rGN,endyr_L_rGN,-10.0,5.0,2);
  init_bounded_number log_avg_F_L_rGN(log_avg_F_rGN_LO,log_avg_F_rGN_HI,log_avg_F_rGN_PH);
  vector log_avg_F_rGN_out(1,8); 
  init_bounded_dev_vector log_dev_F_L_rGN(styr_L_rGN,endyr_L_rGN,log_F_dev_rGN_LO,log_F_dev_rGN_HI,log_F_dev_rGN_PH);    
  vector log_F_dev_rGN_out(styr_L_rGN,endyr_L_rGN);
  matrix F_rGN(styr,endyr,1,nages);
  vector F_rGN_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_rGN;
  number log_F_dev_end_rGN;  
  
  init_bounded_number F_init_ratio(0.1,1.5,1); //scales initial F, which is geometric mean first three yrs
  //init_bounded_number F_init_ratio(F_init_ratio_LO,F_init_ratio_HI,F_init_ratio_PH); //number F_init_ratio;
  //vector F_init_ratio_out(1,8); 

  //--Discard mortality stuff------------------------------------------------------------------------------
  //init_bounded_number log_avg_F_D_cGN(-10.0,0.0,1);
  //init_bounded_dev_vector log_dev_F_D_cGN(styr_D_cHL,endyr_D_cHL,-10.0,5.0,2);
  init_bounded_number log_avg_F_D_cGN(log_avg_F_cGN_D_LO,log_avg_F_cGN_D_HI,log_avg_F_cGN_D_PH);
  vector log_avg_F_cGN_D_out(1,8);  
  init_bounded_dev_vector log_dev_F_D_cGN(styr_D_cHL,endyr_D_cHL,log_F_dev_cGN_D_LO,log_F_dev_cGN_D_HI,log_F_dev_cGN_D_PH);
  vector log_F_dev_cGN_D_out(styr_D_cHL,endyr_D_cHL);
  matrix F_cGN_D(styr,endyr,1,nages);
  vector F_cGN_D_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_cGN_D2;   //avg log deviations in reg selex_phase2 (for estimation 1984-1992, prior to data) 
  number log_F_dev_end_cGN_D;  
  
  //init_bounded_number log_avg_F_D_rHB(-10.0,0.0,1);
  //init_bounded_dev_vector log_dev_F_D_rHB(styr_D_rHB,endyr_D_rHB,-10.0,5.0,2);
  init_bounded_number log_avg_F_D_rHB(log_avg_F_rHB_D_LO,log_avg_F_rHB_D_HI,log_avg_F_rHB_D_PH);
  vector log_avg_F_rHB_D_out(1,8);  
  init_bounded_dev_vector log_dev_F_D_rHB(styr_D_rHB,endyr_D_rHB,log_F_dev_rHB_D_LO,log_F_dev_rHB_D_HI,log_F_dev_rHB_D_PH);
  vector log_F_dev_rHB_D_out(styr_D_rHB,endyr_D_rHB);
  matrix F_rHB_D(styr,endyr,1,nages);
  vector F_rHB_D_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_end_rHB_D;  
       
  //init_bounded_number log_avg_F_D_rGN(-10.0,0.0,1);
  //init_bounded_dev_vector log_dev_F_D_rGN(styr_D_rGN,endyr_D_rGN,-10.0,5.0,2);
  init_bounded_number log_avg_F_D_rGN(log_avg_F_rGN_D_LO,log_avg_F_rGN_D_HI,log_avg_F_rGN_D_PH);
  vector log_avg_F_rGN_D_out(1,8);  
  init_bounded_dev_vector log_dev_F_D_rGN(styr_D_rGN,endyr_D_rGN,log_F_dev_rGN_D_LO,log_F_dev_rGN_D_HI,log_F_dev_rGN_D_PH);
  vector log_F_dev_rGN_D_out(styr_D_rGN,endyr_D_rGN);
  matrix F_rGN_D(styr,endyr,1,nages);
  vector F_rGN_D_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_rGN_D;
  number log_F_dev_end_rGN_D;  

  number Dmort_HL;
  number Dmort_rHB_HL;
  number Dmort_rGN_HL;
  number Dmort_cPT1;
  number Dmort_cPT2;  

//---Per-recruit stuff----------------------------------------------------------------------------------
  vector N_age_spr(1,nages);         //numbers at age for SPR calculations: beginning of year
  vector N_age_spr_spawn(1,nages);   //numbers at age for SPR calculations: at time of peak spawning  
  vector L_age_spr(1,nages);         //catch at age for SPR calculations
  vector Z_age_spr(1,nages);         //total mortality at age for SPR calculations
  vector spr_static(styr,endyr);     //vector of static SPR values by year
  vector F_L_age_spr(1,nages);       //fishing mortality of landings (not discards) at age for SPR calculations
  vector F_spr(1,n_iter_spr);        //values of full F to be used in per-recruit calculations
  vector spr_spr(1,n_iter_spr);      //reproductive capacity-per-recruit values corresponding to F values in F_spr
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
 
  number sdnr_lc_sBT;
  number sdnr_lc_sTV;
  number sdnr_lc_cHL;
  number sdnr_lc_cPT;
  number sdnr_lc_rHB;  
  number sdnr_lc_rHB_D;  
  number sdnr_lc_rGN;  
 
  number sdnr_ac_sBT;
  number sdnr_ac_sTV;
  number sdnr_ac_cHL; 
  number sdnr_ac_cPT;  
  number sdnr_ac_rHB;
  //number sdnr_ac_rGN;
  
  number sdnr_I_sBT;
  number sdnr_I_sTV;
  number sdnr_I_cHL;
  number sdnr_I_rHB;  
  //number sdnr_I_Vid;    
  
//-------Objective function components-----------------------------------------------------------------------------
  number w_L;
  number w_D;

  number w_lenc_sBT;
  number w_lenc_cHL;
  number w_lenc_cPT;
  number w_lenc_rHB;
  number w_lenc_rHB_D;
  number w_lenc_rGN;

  number w_agec_sBT; 
  number w_agec_sTV; 
  number w_agec_cHL; 
  number w_agec_cPT;  
  number w_agec_rHB;
  //number w_ac_rGN;
  
  number w_cpue_sBT;  
  number w_cpue_sTV;
  //number w_I_Vid;
  number w_cpue_cHL;
  number w_cpue_rHB;
  //number w_I_rHD;
  
  number w_rec;
  number w_rec_early;
  number w_rec_end;
  number w_fullF;  
  number w_Ftune;
//  number w_cvlen_dev;
//  number w_cvlen_diff;
  
  number f_sBT_cpue;
  number f_sTV_cpue;
  number f_Vid_cpue;
  number f_cHL_cpue;
  number f_rHB_cpue;
  //number f_rHD_cpue;
 
  number f_cHL_L;  
  number f_cPT_L;
  number f_cTW_L;
  number f_rHB_L;
  number f_rGN_L;

  number f_cGN_D; 
  number f_rHB_D;
  number f_rGN_D;

  number f_sBT_lenc;   
  number f_cHL_lenc;
  number f_cPT_lenc;
  number f_rHB_lenc;
  number f_rHB_D_lenc;
  number f_rGN_lenc;

  number f_sBT_agec;
  number f_sTV_agec;  
  number f_cHL_agec;  
  number f_cPT_agec;    
  number f_rHB_agec; 
  number f_rGN_agec;
  
  number f_cHL_RW_cpue; //random walk component of indices
  number f_rHB_RW_cpue;
  number f_rHD_RW_cpue;   
  
//Penalties and constraints. Not all are used.
  number f_rec_dev;                //weight on recruitment deviations to fit S-R curve
  number f_rec_dev_early;          //extra weight on deviations in first recruitment stanza
  number f_rec_dev_end;            //extra weight on deviations in first recruitment stanza
  number f_Ftune;                  //penalty for tuning F in Ftune yr.  Not applied in final optimization phase.
  number f_fullF_constraint;       //penalty for Fapex>X
  number f_priors;                  //prior information on parameters
  //number f_cvlen_dev_constraint; //deviation penalty on cv's of length at age
  //number f_cvlen_diff_constraint;//first diff penalty on cv's of length at age
  
  objective_function_value fval;
  number fval_data;
  number grad_max;
  
//--Dummy variables ----
  number denom;                   //denominator used in some calculations
  number numer;                   //numerator used in some calculations
     
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
INITIALIZATION_SECTION


//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
GLOBALS_SECTION
  #include "admodel.h"          // Include AD class definitions
  #include "admb2r.cpp"    // Include S-compatible output functions (needs preceding)
	time_t start,finish;
	long hour,minute,second;	
	double elapsed_time;
	
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
RUNTIME_SECTION
 maximum_function_evaluations 10000, 20000, 30000, 50000, 100000, 100000, 100000;
 convergence_criteria 1e-2, 1e-2,1e-2, 1e-4;
 
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PRELIMINARY_CALCS_SECTION

//// Set values of fixed parameters or set initial guess of estimated parameters
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
  Dopen_cHL=Dmort_HL*pow((obs_released_cHL(endyr_D_cHL-2)*obs_released_cHL(endyr_D_cHL-1)*obs_released_cHL(endyr_D_cHL)),(1.0/3.0));
  Dclosed_cHL=Dmort_HL*((obs_released_cHL_closed(endyr_D_cHL-3)+obs_released_cHL_closed(endyr_D_cHL-2)+obs_released_cHL_closed(endyr_D_cHL-1)+obs_released_cHL_closed(endyr_D_cHL))/4.0);
  //Dclosed_cHL=Dmort_HL*pow((obs_released_cHL_closed(endyr_D_cHL-2)*obs_released_cHL_closed(endyr_D_cHL-1)*obs_released_cHL_closed(endyr_D_cHL)),(1.0/3.0));
  Lopen_cHL=Dmort_HL*pow((obs_L_cHL(endyr_L_cHL-2)*obs_L_cHL(endyr_L_cHL-1)*obs_L_cHL(endyr_L_cHL)),(1.0/3.0));
  
  Dopen_cPT=Dmort_cPT2*pow((obs_released_cPT(endyr_D_cPT-2)*obs_released_cPT(endyr_D_cPT-1)*obs_released_cPT(endyr_D_cPT)),(1.0/3.0));
  Dclosed_cPT=Dmort_cPT2*(obs_released_cPT_closed(endyr_D_cPT-3)+(obs_released_cPT_closed(endyr_D_cPT-2)+obs_released_cPT_closed(endyr_D_cPT-1)+obs_released_cPT_closed(endyr_D_cPT))/4.0);
  //Dclosed_cPT=Dmort_cPT2*pow((obs_released_cPT_closed(endyr_D_cPT-2)*obs_released_cPT_closed(endyr_D_cPT-1)*obs_released_cPT_closed(endyr_D_cPT)),(1.0/3.0));
  Lopen_cPT=Dmort_cPT2*pow((obs_L_cPT(endyr_L_cPT-2)*obs_L_cPT(endyr_L_cPT-1)*obs_L_cPT(endyr_L_cPT)),(1.0/3.0));
  
  D_sum_cHLcPT=Dopen_cHL+Dclosed_cHL+Dopen_cPT+Dclosed_cPT;

  Dprop_cGN_sel_D=(Dopen_cHL + Dopen_cPT + Dclosed_cHL*(Dopen_cHL/(Dopen_cHL+Lopen_cHL)) +
                    Dclosed_cPT*(Dopen_cPT/(Dopen_cPT+Lopen_cPT)))/D_sum_cHLcPT; 
  Dprop_cGN_sel_cHL=Dclosed_cHL*(Lopen_cHL/(Dopen_cHL+Lopen_cHL))/D_sum_cHLcPT; 
  Dprop_cGN_sel_cPT=Dclosed_cPT*(Lopen_cPT/(Dopen_cPT+Lopen_cPT))/D_sum_cHLcPT; 

  //discards values for fitting, include discard mortality
  
  obs_cHL_D=Dmort_HL*obs_released_cHL;
  obs_cHL_D(styr_D_cHL_closed,endyr_D_cHL_closed)+=Dmort_HL*obs_released_cHL_closed;
  
  obs_cPT_D(styr_D_cPT,endyr_selex_phase3)=Dmort_cPT1*obs_released_cPT(styr_D_cPT,endyr_selex_phase3);
  obs_cPT_D(endyr_selex_phase3+1,endyr_D_cPT)=Dmort_cPT2*obs_released_cPT(endyr_selex_phase3+1,endyr_D_cPT);
  obs_cPT_D(styr_D_cPT_closed,endyr_D_cPT_closed)+=Dmort_cPT2*obs_released_cPT_closed;
  
  obs_cGN_D=obs_cHL_D+obs_cPT_D;
  
  obs_rHB_D=Dmort_rHB_HL*obs_released_rHB;
  obs_rGN_D=Dmort_rGN_HL*obs_released_rGN;
  
  cGN_D_cv=obs_cv_D_cHL;
 
  Linf=set_Linf(1);
  K=set_K(1);
  t0=set_t0(1);

//  age_limit_8in=t0-log(1.0-limit_8in/Linf)/K; //age at size limit: 8" limit;
//  age_limit_10in=t0-log(1.0-limit_10in/Linf)/K; //age at size limit: 10" limit;
//  age_limit_12in=t0-log(1.0-limit_12in/Linf)/K; //age at size limit: 12" limit;

  M=set_M; 
  M_constant=set_M_constant(1);
//  for (iage=1;iage<=max_obs_age;iage++){Mscale_ages(iage)=iage;}
  
  
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
  //log_q_rHD=set_logq_rHD(1);
  q_rate=set_q_rate;
  q_rate_fcn_cHL=1.0;
  q_rate_fcn_rHB=1.0;
  //q_rate_fcn_rHD=1.0;
  q_DD_beta=set_q_DD_beta;
  q_DD_fcn=1.0;
  q_RW_log_dev_cHL.initialize();
  q_RW_log_dev_rHB.initialize();
  //q_RW_log_dev_rHD.initialize();

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
      //for (iyear=styr_rHD_cpue; iyear<=endyr_rHD_cpue; iyear++)
      //{   if (iyear>styr_rHD_cpue & iyear <=2003) 
       //   {//q_rate_fcn_rHD(iyear)=(1.0+q_rate)*q_rate_fcn_rHD(iyear-1); //compound
       //      q_rate_fcn_rHD(iyear)=(1.0+(iyear-styr_rHD_cpue)*q_rate)*q_rate_fcn_rHD(styr_rHD_cpue);  //linear
       //   }
       //   if (iyear>2003) {q_rate_fcn_rHD(iyear)=q_rate_fcn_rHD(iyear-1);} 
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
  //w_I_rHD=set_w_I_rHD;

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

  selpar_A50_rHD4=set_selpar_A50_rHD4(1);
  selpar_slope_rHD4=set_selpar_slope_rHD4(1); 
  selpar_A502_rHD4=set_selpar_A502_rHD4(1);
  selpar_slope2_rHD4=set_selpar_slope2_rHD4(1); 
  selpar_A50_rHD5=set_selpar_A50_rHD5(1);
  selpar_slope_rHD5=set_selpar_slope_rHD5(1);
  selpar_A502_rHD5=set_selpar_A502_rHD5(1);
  selpar_slope2_rHD5=set_selpar_slope2_rHD5(1);
   
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

 
//Fill in sample sizes of comps, possibly sampled in nonconsec yrs 
//Used primarily for output in R object   
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

//fill in Fs for msy and per-recruit analyses
  F_msy(1)=0.0;  
    for (ff=2;ff<=n_iter_msy;ff++) {F_msy(ff)=F_msy(ff-1)+iter_inc_msy;}
  F_spr(1)=0.0;  
    for (ff=2;ff<=n_iter_spr;ff++){F_spr(ff)=F_spr(ff-1)+iter_inc_spr;}


//fill in F's, Catch matrices, and log rec dev with zero's
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
  
 
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
TOP_OF_MAIN_SECTION
  time(&start);
  arrmblsize=20000000;
  gradient_structure::set_MAX_NVAR_OFFSET(1600);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(2000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(2000000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(500);
  

//>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PROCEDURE_SECTION

 R0=mfexp(log_R0);
 
 //cout<<"start"<<endl;

 //get_M_at_age(); //Needed only if M is estimated
 
 get_length_weight_at_age(); 
 //cout << "got length, weight, fecundity transitions" <<endl;
 get_reprod();
// cout << "got reprod" << endl;
 get_length_at_age_dist(); 
// cout<< "got predicted length at age distribution"<<endl;
 get_weight_at_age_landings();
// cout<< "got weight at age of landings"<<endl; 
 get_spr_F0();
// cout << "got F0 spr" << endl;
 get_selectivity(); 
// cout << "got selectivity" << endl;
 get_mortality(); 
// cout << "got mortalities" << endl;
 get_bias_corr(); 
// cout<< "got recruitment bias correction" << endl;
 get_numbers_at_age(); 
// cout << "got numbers at age" << endl;
 get_landings_numbers();
// cout << "got catch at age" << endl;
 get_landings_wgt();
// cout << "got landings" << endl;
 get_dead_discards();
 //cout << "got discards" << endl;
 get_catchability_fcns();
 //cout << "got catchability_fcns" << endl;
 get_indices();
 //cout << "got indices" << endl;
 get_length_comps();
 //cout<< "got length comps"<< endl;
 get_age_comps();
// cout<< "got age comps"<< endl;
 evaluate_objective_function();
 //cout << "objective function calculations complete" << endl;


//FUNCTION get_M_at_age
//    Mscale_len=Linf*(1.0-mfexp(-K*(Mscale_ages-t0+0.5))); 
//    Mscale_wgt_g=wgtpar_a*pow(Mscale_len,wgtpar_b); 
//    M_lorenzen=3.69*pow(Mscale_wgt_g,-0.305);
//    cum_surv_1plus=mfexp(-max_obs_age*M_constant);  
//    M=M_lorenzen(1,nages)*(-log(cum_surv_1plus)/sum(M_lorenzen(1,max_obs_age)));    


FUNCTION get_length_weight_at_age
  //compute mean length (mm) and weight (whole) at age
    meanlen_TL=Linf*(1.0-mfexp(-K*(agebins-t0+0.5)));    //total length in mm
    wgt_g=wgtpar_a*pow(meanlen_TL,wgtpar_b);             //wgt in grams 
    wgt_kg=g2kg*wgt_g;                                   //wgt in kilograms 
    wgt_mt=g2mt*wgt_g;                                   //mt of whole wgt: g2mt converts g to mt
    wgt_klb=mt2klb*wgt_mt;                               //1000 lb of whole wgt
    wgt_lb=mt2lb*wgt_mt;                                 //1000 lb of whole wgt
    fecundity=fecpar_batches*mfexp(fecpar_a+wgt_g*fecpar_b)/fecpar_scale;    //fecundity at age, scaled

FUNCTION get_reprod 
   //reprod is product of stuff going into reproductive capacity calcs
   reprod=elem_prod(elem_prod(prop_f,maturity_f),fecundity);  
   reprod2=elem_prod(elem_prod(prop_f,maturity_f),wgt_mt);  
 
FUNCTION get_length_at_age_dist
  
//compute matrix of length at age, based on the normal distribution
  //len_cv=len_cv_val; 
  //len_sd=elem_prod(len_cv, meanlen_TL);
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

	//First size limit 8" selex_phase1 for both commercial and recreational discards (1984 through 1998)
	zscore_lsizelim1=(sizelim1-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim1=cumd_norm(zscore_lsizelim1);   //includes any probability mass below zero
	prob_belowsizelim_block1(iage)=	cprob_lsizelim1-cprob_lzero; //removes any probability mass below zero

	//Second size limit 10" selex_phase2 for both commercial and recreational discards and selex_phase3 for recreational landings (1999 through 2012 for commercial) (1999 through 2006 for recreational)
	zscore_lsizelim2=(sizelim2-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim2=cumd_norm(zscore_lsizelim2);   //includes any probability mass below zero
	prob_belowsizelim_block2(iage)=	cprob_lsizelim2-cprob_lzero; //removes any probability mass below zero

	//Third size limit 11" selex_phase3 for commercial (2013 through terminal year)
	zscore_lsizelim3=(sizelim3-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim3=cumd_norm(zscore_lsizelim3);                 //includes any probability mass below zero
	prob_belowsizelim_block3(iage)=	cprob_lsizelim3-cprob_lzero; //removes any probability mass below zero

	//Fourth size limit 12" selex_phase3 for recreational discards and selex_phase4 for recreational landings (2007 through 2012)
	zscore_lsizelim4=(sizelim4-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim4=cumd_norm(zscore_lsizelim4);                 //includes any probability mass below zero
	prob_belowsizelim_block4(iage)=	cprob_lsizelim4-cprob_lzero; //removes any probability mass below zero

	//Fifth size limit 13" selex_phase4 for recreational discards and selex_phase5 for recreational landings (2013 through the terminal year)
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
   
FUNCTION get_weight_at_age_landings
  //fleets under identical size limits are set equal at end of fcn
  for (iyear=styr; iyear<=endyr_selex_phase1; iyear++)
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
    len_rGN_D_mm(iyear,iage)=sum(elem_prod(lenprob_rGN_D2(iage),lenbins)); //assumes same size distn in selex_phase1 as in selex_phase2       
    }
    wgt_rGN_D_klb(iyear)=g2klb*wgtpar_a*pow(len_rGN_D_mm(iyear),wgtpar_b);
  } // end iyear loop
  
  
  for (iyear=(endyr_selex_phase1+1); iyear<=endyr_selex_phase2; iyear++)
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

  for (iyear=(endyr_selex_phase2+1); iyear<=endyr; iyear++) //cGN only
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
  
  for (iyear=(endyr_selex_phase2+1); iyear<=endyr_selex_phase3; iyear++) //rec only
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
  
  for (iyear=(endyr_selex_phase3+1); iyear<=endyr; iyear++) //rec only
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
   
   for (iyear=(endyr_selex_phase4+1); iyear<=endyr; iyear++) //rec only
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
    
FUNCTION get_spr_F0
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
 

FUNCTION get_selectivity

// ------- compute landings selectivities by selex_phase

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
//  selpar_slope_rGN4=selpar_slope_rHB4;
  
  sel_rGN1=logistic(agebins, selpar_A50_rGN1, selpar_slope_rGN1);
  sel_rGN2=logistic(agebins, selpar_A50_rGN2, selpar_slope_rGN2);
  sel_rGN3=logistic(agebins, selpar_A50_rGN3, selpar_slope_rGN3);
  sel_rGN4=logistic(agebins, selpar_A50_rGN4, selpar_slope_rGN4);
  sel_rGN5=logistic(agebins, selpar_A50_rGN5, selpar_slope_rGN5);


//-----------fill in years--------------------------------------------
  
  //selex_phase1:   
  for (iyear=styr; iyear<=endyr_selex_phase1; iyear++)
  {
     sel_sBT(iyear)=sel_sBT_vec;
     sel_sTV(iyear)=sel_sTV_vec;
     sel_cHL(iyear)=sel_cHL_2; //commercial handline sel mirrors selex_phase2
     sel_cPT(iyear)=sel_cPT_2; //commercial handline sel mirrors selex_phase2
     sel_rHB(iyear)=sel_rHB_1; 
     sel_rGN(iyear)=sel_rGN1;      
  }

  //selex_phase2: 
  for (iyear=endyr_selex_phase1+1; iyear<=endyr_selex_phase2; iyear++)
  {     
     sel_sBT(iyear)=sel_sBT_vec;
     sel_sTV(iyear)=sel_sTV_vec;
     sel_cHL(iyear)=sel_cHL_2;
     sel_cPT(iyear)=sel_cPT_2; 
     sel_rHB(iyear)=sel_rHB_2;
     sel_rGN(iyear)=sel_rGN2;     
  }

  //selex_phase3 
  for (iyear=endyr_selex_phase2+1; iyear<=endyr; iyear++)
  {
     sel_sBT(iyear)=sel_sBT_vec;
     sel_sTV(iyear)=sel_sTV_vec;
     sel_cHL(iyear)=sel_cHL_3; 
     sel_cPT(iyear)=sel_cPT_3;     
     sel_rHB(iyear)=sel_rHB_3;
     sel_rGN(iyear)=sel_rGN3;     
  }   
  //selex_phase4: rec only, overwrites yrs calculated for selex_phase3
  for (iyear=endyr_selex_phase3; iyear<endyr_selex_phase4; iyear++)
  {
     sel_rHB(iyear)=sel_rHB_4;
     sel_rGN(iyear)=sel_rGN4;     
  }   
   //selex_phase5: Comm and rec, overwrites yrs previously calculated
   for (iyear=endyr_selex_phase4; iyear<=endyr; iyear++)
  {
     sel_cHL(iyear)=sel_cHL_4; 
     sel_cPT(iyear)=sel_cPT_4;
	 sel_rHB(iyear)=sel_rHB_5;
     sel_rGN(iyear)=sel_rGN5;     
  }   
  //set selectivities that mirror others
  sel_cTW=sel_cPT;  
  //sel_rGN=sel_rHB;
  
  
//---Discard selectivities---------------------------------    
//---------------------------------------------------------   
  //rGN and cGN mirror headboat discard selectivity

  selpar_Age0_rHB_D=1.0/(1.0+mfexp(-selpar_logit_Age0_rHB_D));
  selpar_Age1_rHB_D=1.0/(1.0+mfexp(-selpar_logit_Age1_rHB_D)); 
  selpar_Age2_rHB_D=1.0/(1.0+mfexp(-selpar_logit_Age2_rHB_D));

  sel_rHB_D_4=logistic_exponential(agebins, selpar_A50_rHD4, selpar_slope_rHD4, selpar_A502_rHD4, selpar_slope2_rHD4);//logistic_double(agebins, selpar_A50_rHD4, selpar_slope_rHD4, selpar_A502_rHD4, selpar_slope2_rHD4); 
  sel_rHB_D_5=logistic_exponential(agebins, selpar_A50_rHD5, selpar_slope_rHD5, selpar_A502_rHD5, selpar_slope2_rHD5);//logistic_double(agebins, selpar_A50_rHD5, selpar_slope_rHD5, selpar_A502_rHD5, selpar_slope2_rHD5); 
   
 //Assume same sel of age 0's across selex_phases 
  sel_rHB_D_2(1)=selpar_Age0_rHB_D; 
  sel_rHB_D_3(1)=selpar_Age0_rHB_D; 
  //sel_rHB_D_4(1)=selpar_Age0_rHB_D;
  //sel_rHB_D_5(1)=selpar_Age0_rHB_D;
  sel_cGN_D_3(1)=selpar_Age0_rHB_D;
 //Assume same sel of age 1's across selex_phases 
  sel_rHB_D_2(2)=selpar_Age1_rHB_D; 
  sel_rHB_D_3(2)=selpar_Age1_rHB_D; 
  //sel_rHB_D_4(2)=selpar_Age1_rHB_D;
  //sel_rHB_D_5(2)=selpar_Age1_rHB_D;
  sel_cGN_D_3(2)=selpar_Age1_rHB_D;
 //Assume same sel of age 2's across selex_phases 
  sel_rHB_D_2(3)=selpar_Age2_rHB_D; 
  sel_rHB_D_3(3)=selpar_Age2_rHB_D; 
  //sel_rHB_D_4(3)=selpar_Age2_rHB_D;
  //sel_rHB_D_5(3)=selpar_Age2_rHB_D;  
  sel_cGN_D_3(3)=selpar_Age2_rHB_D;
 //Assume full sel at age 3 across selex_phases 
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

  //selex_phase1: assumed same as in selex_phase1, no commercial discards
  for (iyear=styr; iyear<=endyr_selex_phase1; iyear++)
      {sel_rHB_D(iyear)=sel_rHB_D_2;}

  //selex_phase2: 
  for (iyear=endyr_selex_phase1+1; iyear<=endyr_selex_phase2; iyear++)
      {sel_rHB_D(iyear)=sel_rHB_D_2;
       sel_cGN_D(iyear)=sel_rHB_D_2;
      }
  //selex_phase3: Starts in 1999
  for (iyear=endyr_selex_phase2+1; iyear<=endyr; iyear++)
  {sel_rHB_D(iyear)=sel_rHB_D_3;
   sel_cGN_D(iyear)=sel_rHB_D_3;}  
  
  //selex_phase4: Starts in 2007, rHB and rGN only, overwrites last few yrs calculated for selex_phase3
  for (iyear=endyr_selex_phase3+1; iyear<=endyr_selex_phase4; iyear++)
  {sel_rHB_D(iyear)=sel_rHB_D_4;
   sel_cGN_D(iyear)=sel_rHB_D_3;
   
  }
  //Discard quota: wghted average,, overwrites last few yrs calculated for selex_phase3 styr cGN closed=2009
  for (iyear=styr_cGN_closed; iyear<=endyr_selex_phase4; iyear++)
  {sel_cGN_D(iyear)=Dprop_cGN_sel_D*sel_rHB_D_3 + Dprop_cGN_sel_cHL*sel_cHL_2 +
                     Dprop_cGN_sel_cPT*sel_cPT_2;
   sel_cGN_D(iyear)=sel_cGN_D(iyear)/max(sel_cGN_D(iyear));} 
   
  //cout<<"sel_rHB_D after selex_phase4 loop"<<sel_rHB_D<<endl;
  //selex_phase5: Starts in 2013 
  for (iyear=endyr_selex_phase4+1; iyear<=endyr; iyear++)
  {sel_rHB_D(iyear)=sel_rHB_D_5;
   sel_cGN_D(iyear)=sel_cGN_D_3;
   //Discard quota: wghted average,, overwrites last few yrs calculated for selex_phase3 styr cGN closed=2009
   sel_cGN_D(iyear)=Dprop_cGN_sel_D*sel_cGN_D_3 + Dprop_cGN_sel_cHL*sel_cHL_2 +
                     Dprop_cGN_sel_cPT*sel_cPT_2;
   sel_cGN_D(iyear)=sel_cGN_D(iyear)/max(sel_cGN_D(iyear));} 
  
   sel_rGN_D=sel_rHB_D;

FUNCTION get_mortality
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
    if(iyear > endyr_selex_phase1 & iyear < styr_D_cHL)
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
 
 
FUNCTION get_bias_corr
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

//  cout<< "get bias corr" << var_rec_dev <<endl;

FUNCTION get_numbers_at_age
//Initialization
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


//Initial equilibrium age structure
  N_spr_initial(1)=1.0*mfexp(-1.0*Z_initial(1)*spawn_time_frac); //at peak spawning time;
  for (iage=2; iage<=nages; iage++)
    {
      N_spr_initial(iage)=N_spr_initial(iage-1)*
                   mfexp(-1.0*(Z_initial(iage-1)*(1.0-spawn_time_frac) + Z_initial(iage)*spawn_time_frac)); 
    }
  N_spr_initial(nages)=N_spr_initial(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
//    N_spr_F_init_mdyr(1,(nages-1))=elem_prod(N_spr_initial(1,(nages-1)),
//                                   mfexp((-1.*(M(nages-1)+ F_initial))/2.0));   

  spr_initial=sum(elem_prod(N_spr_initial,reprod));

  if (styr==styr_rec_dev) {R1=(R0/((5.0*steep-1.0)*spr_initial))*
                 (4.0*steep*spr_initial-spr_F0*(1.0-steep));} //without bias correction (deviation added later)
  else {R1=(R0/((5.0*steep-1.0)*spr_initial))*
                 (BiasCor*4.0*steep*spr_initial-spr_F0*(1.0-steep));} //with bias correction                 

  if(R1<10.0) {R1=10.0;} //Avoid negative (or unreasonably low) popn sizes during search algorithm

  
   
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
  
  N_mdyr(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*0.5))); //mid year 
  N_spawn(styr)(1,nages)=elem_prod(N(styr)(1,nages),(mfexp(-1.*(Z_initial(1,nages))*spawn_time_frac))); //peak spawning time 

  SSB(styr)=sum(elem_prod(N_spawn(styr),reprod));
  MatFemB(styr)=sum(elem_prod(N_spawn(styr),reprod2));
  B_q_DD(styr)=sum(elem_prod(N(styr)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
    
//Rest of years 
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


//Time series of interest
  rec=column(N,1);
  SdS0=SSB/S0;
  
FUNCTION get_landings_numbers //Baranov catch eqn
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
 
 
FUNCTION get_landings_wgt

////---Predicted landings------------------------
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
 
   
    
FUNCTION get_dead_discards //Baranov catch eqn
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

FUNCTION get_catchability_fcns    
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
      //for (iyear=styr_rHD_cpue; iyear<=endyr_rHD_cpue; iyear++)
      //{   if (iyear>styr_rHD_cpue & iyear <=2003) 
      //    {//q_rate_fcn_rHD(iyear)=(1.0+q_rate)*q_rate_fcn_rHD(iyear-1); //compound
      //       q_rate_fcn_rHD(iyear)=(1.0+(iyear-styr_rHD_cpue)*q_rate)*q_rate_fcn_rHD(styr_rHD_cpue);  //linear
      //    }
      //    if (iyear>2003) {q_rate_fcn_rHD(iyear)=q_rate_fcn_rHD(iyear-1);} 
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
  //rHD cpue
  //q_rHD(styr_rHD_cpue)=mfexp(log_q_rHD);
  //for (iyear=styr_rHD_cpue; iyear<=endyr_rHD_cpue; iyear++)
  //{   //index in number units
  //    N_rHD(iyear)=elem_prod(N_mdyr(iyear),sel_rHB_D(iyear)); 
  //    pred_rHD_cpue(iyear)=q_rHD(iyear)*q_rate_fcn_rHD(iyear)*q_DD_fcn(iyear)*sum(N_rHD(iyear));
  //    if (iyear<endyr_rHD_cpue){q_rHD(iyear+1)=q_rHD(iyear)*mfexp(q_RW_log_dev_rHD(iyear));}
  //}
  
FUNCTION get_length_comps
  
  //sBT
  for (iyear=1;iyear<=nyr_lenc_sBT;iyear++)
      {pred_sBT_lenc(iyear)=(N_sBT(yrs_lenc_sBT(iyear))*lenprob)/sum(N_sBT(yrs_lenc_sBT(iyear)));} 

  //Commercial lines
  for (iyear=1;iyear<=nyr_lenc_cHL;iyear++) //all yrs within selex_phases 2,3
  {  if (yrs_lenc_cHL(iyear)<=endyr_selex_phase2)
     {pred_cHL_lenc(iyear)=(L_cHL_num(yrs_lenc_cHL(iyear))*lenprob_cHL2)
                          /sum(L_cHL_num(yrs_lenc_cHL(iyear)));       
     } 
     if (yrs_lenc_cHL(iyear)>endyr_selex_phase2)
     {pred_cHL_lenc(iyear)=(L_cHL_num(yrs_lenc_cHL(iyear))*lenprob_cHL3)
                          /sum(L_cHL_num(yrs_lenc_cHL(iyear)));       
     }  
  }

  //Commercial pots: pooled all from selex_phase2
  L_cPT_num_pool.initialize();
  for (iyear=1;iyear<=nyr_lenc_pool_cPT;iyear++)
  {  L_cPT_num_pool_yr(iyear)=nsamp_lenc_pool_cPT(iyear)*L_cPT_num(yrs_lenc_pool_cPT(iyear))
                            /sum(L_cPT_num(yrs_lenc_pool_cPT(iyear)));                             
     if (yrs_lenc_pool_cPT(iyear)<=endyr_selex_phase2) {L_cPT_num_pool(1)+=L_cPT_num_pool_yr(iyear);}                             
  } 
  for (iyear=1;iyear<=nyr_lenc_cPT;iyear++) //all yrs within selex_phase 2
  {  if (yrs_lenc_cPT(iyear)<=endyr_selex_phase2)
       {pred_cPT_lenc(iyear)=(L_cPT_num_pool(iyear)*lenprob_cPT2)/sum(L_cPT_num_pool(iyear)); } 
	   pred_cPT_lenc(iyear)=(L_cPT_num(yrs_lenc_cPT(iyear))*lenprob_cPT3)  //added to calculate comps for last selex_phase
						  /sum(L_cPT_num(yrs_lenc_cPT(iyear)));
  }  
   

 
 //Headboat 
  for (iyear=1;iyear<=nyr_lenc_rHB;iyear++)  //all in selex_phases 1,2,3
  {  if (yrs_lenc_rHB(iyear)<=endyr_selex_phase1)
     {pred_rHB_lenc(iyear)=(L_rHB_num(yrs_lenc_rHB(iyear))*lenprob_rHB1)
                          /sum(L_rHB_num(yrs_lenc_rHB(iyear)));       
     } 
     if (yrs_lenc_rHB(iyear)>endyr_selex_phase1 & yrs_lenc_rHB(iyear)<=endyr_selex_phase2)
     {pred_rHB_lenc(iyear)=(L_rHB_num(yrs_lenc_rHB(iyear))*lenprob_rHB2)
                          /sum(L_rHB_num(yrs_lenc_rHB(iyear)));       
     }  
     if (yrs_lenc_rHB(iyear)>endyr_selex_phase2)
     {pred_rHB_lenc(iyear)=(L_rHB_num(yrs_lenc_rHB(iyear))*lenprob_rHB3)
                          /sum(L_rHB_num(yrs_lenc_rHB(iyear)));       
     }  
  }
 //rHB discards 
  for (iyear=1;iyear<=nyr_lenc_rHB_D;iyear++) //all yrs within selex_phase3,4
  {  if (yrs_lenc_rHB_D(iyear)<=endyr_selex_phase3)
     {pred_rHB_D_lenc(iyear)=(D_rHB_num(yrs_lenc_rHB_D(iyear))*lenprob_rHB_D3)
                        /sum(D_rHB_num(yrs_lenc_rHB_D(iyear)));}      
    if (yrs_lenc_rHB_D(iyear)>endyr_selex_phase3)
     {pred_rHB_D_lenc(iyear)=(D_rHB_num(yrs_lenc_rHB_D(iyear))*lenprob_rHB_D4)
                        /sum(D_rHB_num(yrs_lenc_rHB_D(iyear)));}                      
  }

 
 //MRIP
  for (iyear=1;iyear<=nyr_lenc_rGN;iyear++)  //all in selex_phases 1,2,3
  {  if (yrs_lenc_rGN(iyear)<=endyr_selex_phase1)
     {pred_rGN_lenc(iyear)=(L_rGN_num(yrs_lenc_rGN(iyear))*lenprob_rGN1)
                          /sum(L_rGN_num(yrs_lenc_rGN(iyear)));       
     } 
     if (yrs_lenc_rGN(iyear)>endyr_selex_phase1 & yrs_lenc_rGN(iyear)<=endyr_selex_phase2)
     {pred_rGN_lenc(iyear)=(L_rGN_num(yrs_lenc_rGN(iyear))*lenprob_rGN2)
                          /sum(L_rGN_num(yrs_lenc_rGN(iyear)));       
     }  
     if (yrs_lenc_rGN(iyear)>endyr_selex_phase2 & yrs_lenc_rGN(iyear)<=endyr_selex_phase3)
     {pred_rGN_lenc(iyear)=(L_rGN_num(yrs_lenc_rGN(iyear))*lenprob_rGN3)
                          /sum(L_rGN_num(yrs_lenc_rGN(iyear)));       
     }  
     if (yrs_lenc_rGN(iyear)>endyr_selex_phase3 &
	 yrs_lenc_rGN(iyear)<=endyr_selex_phase4)
     {pred_rGN_lenc(iyear)=(L_rGN_num(yrs_lenc_rGN(iyear))*lenprob_rGN4)
                          /sum(L_rGN_num(yrs_lenc_rGN(iyear)));       
     }  
	 if (yrs_lenc_rGN(iyear)>endyr_selex_phase4)
     {pred_rGN_lenc(iyear)=(L_rGN_num(yrs_lenc_rGN(iyear))*lenprob_rGN5)
                          /sum(L_rGN_num(yrs_lenc_rGN(iyear)));       
     }  
  } 
 
  
FUNCTION get_age_comps

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
  
    
////--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_weighted_current 
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
  
FUNCTION get_msy
  
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


//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_per_recruit_stuff

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
FUNCTION get_miscellaneous_stuff
  
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
  
//  steep_sd=steep;
//  fullF_sd=Fsum;
  
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


FUNCTION get_effective_sample_sizes
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
//FUNCTION get_projection
//  
//    switch(Fproj_switch){
//       case 1: //F=Fcurrent
//          F_reg_proj=Fend_mean;
//          break;
//       case 2: //F=Fmsy
//          F_reg_proj=F_msy_out;
//          break;
//       case 3: //F=F30
//          F_reg_proj=F30_out;
//          break;     
//       case 4: //F=F40
//          F_reg_proj=F40_out;
//          break;          		  
//       default: // no such switch available
//          cout << "Error in input: Projection switch Fproj_switch must be set to 1, 2, 3, or 4." << endl;
//          cout << "Presently it is set to " << Fproj_switch <<"."<< endl;
//          exit(0);          
//   }
//
//  N_proj(styr_proj)=N(endyr+1); //initial conditions computed previously
// 
//  for (iyear=styr_proj; iyear<=endyr_proj; iyear++) //recruitment follows S-R curve (with bias correction) exactly
//  {     
//        if (iyear<styr_regs) {F_proj(iyear)=Fend_mean;}
//		else {F_proj(iyear)=Fproj_mult*F_reg_proj;}
//		
//		FL_age_proj=sel_wgted_L*F_proj(iyear);
//		FD_age_proj=sel_wgted_D*F_proj(iyear);
//		
//        Z_proj(iyear)=M+FL_age_proj+FD_age_proj;
//        N_spawn_proj(iyear)(1,nages)=elem_prod(N_proj(iyear)(1,nages),(mfexp(-1.*(Z_proj(iyear)(1,nages))*spawn_time_frac))); //peak spawning time
//		SSB_proj(iyear)= sum(elem_prod(N_spawn_proj(iyear),reprod));
//        B_proj(iyear)=sum(elem_prod(N_proj(iyear),wgt_mt)); //uses spawning weight
//	     
//		for (iage=1; iage<=nages; iage++)
//			{L_age_proj(iyear,iage)=N_proj(iyear,iage)*FL_age_proj(iage)*(1.-mfexp(-1.*Z_proj(iyear,iage)))/Z_proj(iyear,iage);
//		     D_age_proj(iyear,iage)=N_proj(iyear,iage)*FD_age_proj(iage)*(1.-mfexp(-1.*Z_proj(iyear,iage)))/Z_proj(iyear,iage);
//			}          
//        L_knum_proj(iyear)=sum(L_age_proj(iyear))/1000.0;
//	    D_knum_proj(iyear)=sum(D_age_proj(iyear))/1000.0;
//	    L_klb_proj(iyear)=sum(elem_prod(L_age_proj(iyear),wgt_wgted_L_klb));     //in 1000 lb
//        D_klb_proj(iyear)=sum(elem_prod(D_age_proj(iyear),wgt_wgted_D_klb));     //in 1000 lb
//		
//		if (iyear<endyr_proj) {
//			N_proj(iyear+1,1)=BiasCor*SR_func(R0, steep, spr_F0, SSB_proj(iyear),SR_switch);
//			N_proj(iyear+1)(2,nages)=++elem_prod(N_proj(iyear)(1,nages-1),(mfexp(-1.*Z_proj(iyear)(1,nages-1))));
//			N_proj(iyear+1,nages)+=N_proj(iyear,nages)*mfexp(-1.*Z_proj(iyear,nages)); //plus group		
//		}
//  }
//   R_proj=column(N_proj,1);                          
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   


FUNCTION evaluate_objective_function
  fval=0.0;
  fval_data=0.0;

//---likelihoods---------------------------

//---Indices-------------------------------

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
  
  //f_rHD_cpue=0.0;
  //f_rHD_cpue=lk_lognormal(pred_rHD_cpue, obs_rHD_cpue, rHD_cpue_cv, w_I_rHD);
  //fval+=f_rHD_cpue;
  //fval_data+=f_rHD_cpue;

////---Landings-------------------------------
 
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


//---Discards-------------------------------

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


//---Length comps-------------------------------

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


//---Age comps-------------------------------

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
//-----------Constraints and penalties--------------------------------

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
//	cout<< "log rec devs" << log_dev_rec << endl;

//  f_rec_dev_early=0.0; //possible extra constraint on early rec deviations
//  if (styr_rec_dev<endyr_rec_phase1)
//    {  
//      f_rec_dev_early=pow(log_dev_rec(styr_rec_dev),2);
//      for(iyear=(styr_rec_dev+1); iyear<=endyr_rec_phase1; iyear++)
//      {f_rec_dev_early+=pow((log_dev_rec(iyear)-R_autocorr*log_dev_rec(iyear-1)),2);}
//    }
//  fval+=w_rec_early*f_rec_dev_early;

//  f_rec_dev_end=0.0; //possible extra constraint on ending rec deviations
//  if (endyr_rec_phase2<endyr)
//    {  
//      for(iyear=(endyr_rec_phase2+1); iyear<=endyr; iyear++)
//      {f_rec_dev_end+=pow((log_dev_rec(iyear)-R_autocorr*log_dev_rec(iyear-1)),2);}
//    }
//  fval+=w_rec_end*f_rec_dev_end;

//  f_Ftune=0.0; 
//  if (set_Ftune>0.0 && !last_phase()) {f_Ftune=square(Fapex(set_Ftune_yr)-set_Ftune);}
//  fval+=w_Ftune*f_Ftune;
  
//  //code below contingent on four phases
//  f_fullF_constraint=0.0;
//  if (!last_phase())
//  {for (iyear=styr; iyear<=endyr; iyear++)
//       {if (Fapex(iyear)>3.0){f_fullF_constraint+=mfexp(Fapex(iyear)-3.0);}}
//   if (current_phase()==1) {w_fullF=set_w_fullF;}
//   if (current_phase()==2) {w_fullF=set_w_fullF/10.0;}  
//   if (current_phase()==3) {w_fullF=set_w_fullF/100.0;}
//  }

//  fval+=w_fullF*f_fullF_constraint;
//
//  f_fullF_constraint=0.0;
//  for (iyear=styr; iyear<=endyr; iyear++)
//       {if (Fapex(iyear)>3.0){f_fullF_constraint+=mfexp(Fapex(iyear)-3.0);}}      
//  fval+=w_fullF*f_fullF_constraint;  
    
//  f_cvlen_diff_constraint=0.0;
//    f_cvlen_diff_constraint=norm2(first_difference(log_len_cv_dev));
//  fval+=w_cvlen_diff*f_cvlen_diff_constraint;
//  
//  f_cvlen_dev_constraint=0.0;
//    f_cvlen_dev_constraint=norm2(log_len_cv_dev);  
//  fval+=w_cvlen_dev*f_cvlen_dev_constraint;
  
  
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
  
  //f_rHD_RW_cpue=0.0;
  //for (iyear=styr_rHD_cpue; iyear<endyr_rHD_cpue; iyear++)
  //    {f_rHD_RW_cpue+=square(q_RW_log_dev_rHD(iyear))/(2.0*set_q_RW_rHD_var);}      
  //fval+=f_rHD_RW_cpue;  
  
  
//---Priors---------------------------------------------------
//neg_log_prior arguments: estimate, prior, variance, pdf type
//Variance input as a negative value is considered to be CV in arithmetic space (CV=-1 implies loose prior) 
//pdf type(1=none, 2=lognormal, 3=normal, 4=beta)

  f_priors=0.0; 
  f_priors+=neg_log_prior(steep, set_steep(5), set_steep(6), set_steep(7)); 
//	cout<<f_priors<<endl;
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
    
  f_priors+=neg_log_prior(selpar_A50_rHD4, set_selpar_A50_rHD4(5),set_selpar_A50_rHD4(6),set_selpar_A50_rHD4(7)); 
  f_priors+=neg_log_prior(selpar_slope_rHD4, set_selpar_slope_rHD4(5),set_selpar_slope_rHD4(6),set_selpar_slope_rHD4(7)); 
  f_priors+=neg_log_prior(selpar_A502_rHD4, set_selpar_A502_rHD4(5),set_selpar_A502_rHD4(6),set_selpar_A502_rHD4(7)); f_priors+=neg_log_prior(selpar_A50_rHD5, set_selpar_A50_rHD5(5),set_selpar_A50_rHD5(6),set_selpar_A50_rHD5(7));   
  f_priors+=neg_log_prior(selpar_slope_rHD5, set_selpar_slope_rHD5(5),set_selpar_slope_rHD5(6),set_selpar_slope_rHD5(7)); 
  f_priors+=neg_log_prior(selpar_A502_rHD5, set_selpar_A502_rHD5(5),set_selpar_A502_rHD5(6),set_selpar_A502_rHD5(7));
  
  fval+=f_priors;

  //cout << "fval = " << fval << "  fval_data = " << fval_data << endl;


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
//Spawner-recruit function (Beverton-Holt)
FUNCTION dvariable SR_func(const dvariable& R0, const dvariable& h, const dvariable& spr_F0, const dvariable& SSB)
  //R0=virgin recruitment, h=steepness, spr_F0=spawners per recruit @ F=0, SSB=spawning biomass
  RETURN_ARRAYS_INCREMENT();
  dvariable Recruits_Tmp;
  Recruits_Tmp=((0.8*R0*h*SSB)/(0.2*R0*spr_F0*(1.0-h)+(h-0.2)*SSB));
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
  //dzero is small value to avoid log(0) during search
  RETURN_ARRAYS_INCREMENT();
  dvariable LkvalTmp;
  dvar_vector var(cv.indexmin(),cv.indexmax()); //variance in log space
  var=log(1.0+square(cv/wgt_dat));   // convert cv in arithmetic space to variance in log space
  LkvalTmp=sum(0.5*elem_div(square(log(elem_div((pred+dzero),(obs+dzero)))),var) );
  RETURN_ARRAYS_DECREMENT();
  return LkvalTmp;

//-----------------------------------------------------------------------------------
//Likelihood contribution: multinomial
FUNCTION dvariable lk_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const double& minSS, const dvariable& wgt_dat)
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
//Likelihood contribution: Dirichlet-multinomial
FUNCTION dvariable lk_dirichlet_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS, const dvariable& log_dir_par)
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

// //Likelihood contribution: Dirichlet-multinomial
// FUNCTION dvariable lk_dirichlet_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS, const dvariable& log_dir_par)
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
  

//-----------------------------------------------------------------------------------
//Likelihood contribution: priors
FUNCTION  dvariable neg_log_prior(dvariable pred, const double& prior, dvariable var, int pdf)
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

//-----------------------------------------------------------------------------------
REPORT_SECTION

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

	   selpar_A50_rHD4_out(8)=selpar_A50_rHD4; selpar_A50_rHD4_out(1,7)=set_selpar_A50_rHD4;
       selpar_slope_rHD4_out(8)=selpar_slope_rHD4; selpar_slope_rHD4_out(1,7)=set_selpar_slope_rHD4;
	   selpar_A502_rHD4_out(8)=selpar_A502_rHD4; selpar_A502_rHD4_out(1,7)=set_selpar_A502_rHD4;
       selpar_slope2_rHD4_out(8)=selpar_slope2_rHD4; selpar_slope2_rHD4_out(1,7)=set_selpar_slope2_rHD4;
	   selpar_A50_rHD5_out(8)=selpar_A50_rHD5; selpar_A50_rHD5_out(1,7)=set_selpar_A50_rHD5;
       selpar_slope_rHD5_out(8)=selpar_slope_rHD5; selpar_slope_rHD5_out(1,7)=set_selpar_slope_rHD5;
	   selpar_A502_rHD5_out(8)=selpar_A502_rHD5; selpar_A502_rHD5_out(1,7)=set_selpar_A502_rHD5;
       selpar_slope2_rHD5_out(8)=selpar_slope2_rHD5; selpar_slope2_rHD5_out(1,7)=set_selpar_slope2_rHD5;
	   
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
        
       log_rec_dev_out(styr_rec_dev, endyr_rec_dev)=log_dev_rec;
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
  
  
