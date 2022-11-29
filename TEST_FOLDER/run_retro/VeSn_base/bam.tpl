//##  Author: NMFS, Beaufort Lab, Sustainable Fisheries Branch
//##  Analyst: Kevin Craig 
//##  Species: Vermilion Snapper
//##  Region: US South Atlantic
//##  SEDAR: 55
//##  Date: 2022-11-10 17:28:29


//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##
//##  SEDAR 55  Vermilion Snapper  assessment model, Fall 2017
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
//init_int endyr_rec_phase3;

// ending years for selectivity blocks
init_int endyr_selex_phase1;
init_int endyr_selex_phase2;
init_int endyr_selex_phase3;
init_int endyr_2008; //commercial discards include legal sized fish in 2009-2011

//Size limits
//init_number limit_10in;  //10 inch limit in mm
//init_number limit_11in;  //11 inch limit in mm
//init_number limit_12in; //12 inch limit in mm

init_number sizelim1;  //10 inch limit in mm
init_number sizelim2;  //11 inch limit in mm
init_number sizelim3; //12 inch limit in mm

////4 periods of rec size regs: styr-91 no restrictions, 1992-98 10-inch TL, 1999-06 11-in TL, 2006-endyr 12-in TL
////2 periods of cGN size regs: styr-91 no restrictions, 1992-endyr 12-in TL
//init_int endyr_period1;
//init_int endyr_period2;
//init_int endyr_period3;
//init_int endyr_2008; //commercial discards include legal sized fish in 2009-2011
//
//init_number limit_10in;  //10 inch limit in mm
//init_number limit_11in;  //11 inch limit in mm
//init_number limit_12in; //12 inch limit in mm

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
//init_int nlenbins_plus;
init_number lenbins_width;  //width of length bins (mm)


//Vector of lengths for length bins (mm)(midpoint) 
init_vector lenbins(1,nlenbins);
//kc init_vector lenbins_plus(1,nlenbins_plus);

//Max F used in spr and msy calcs
init_number max_F_spr_msy;
//!!cout << max_F_spr_msy<< endl;

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

//###MARMAP Florida Snapper Trap combination#################
//CPUE
init_int styr_cpue_sFT;
init_int endyr_cpue_sFT;
init_vector obs_cpue_sFT(styr_cpue_sFT,endyr_cpue_sFT);   //Observed CPUE
init_vector obs_cv_cpue_sFT(styr_cpue_sFT,endyr_cpue_sFT);    //CV of cpue

// Length Compositions 
init_int nyr_lenc_sFT;
init_ivector yrs_lenc_sFT(1,nyr_lenc_sFT);
init_vector nsamp_lenc_sFT(1,nyr_lenc_sFT);
init_vector nfish_lenc_sFT(1,nyr_lenc_sFT);
init_matrix obs_lenc_sFT(1,nyr_lenc_sFT,1,nlenbins);
init_vector sel_sFT_vec(1,nages);
//!!cout<<sel_sFT_vec<<endl;
//###MARMAP sCT/sVD combination#################
//CPUE
init_int styr_cpue_sCT;
init_int endyr_cpue_sCT;
init_vector obs_cpue_sCT(styr_cpue_sCT,endyr_cpue_sCT);   //Observed CPUE
init_vector obs_cv_cpue_sCT(styr_cpue_sCT,endyr_cpue_sCT);    //CV of cpue

// Length Compositions 
//init_int nyr_sCT_lenc;
//init_ivector yrs_sCT_lenc(1,nyr_sCT_lenc);
//init_vector nsamp_sCT_lenc(1,nyr_sCT_lenc);
//init_vector nfish_sCT_lenc(1,nyr_sCT_lenc);
//init_matrix obs_sCT_lenc(1,nyr_sCT_lenc,1,nlenbins);

// Age Compositions 
init_int nyr_agec_sCT;
init_ivector yrs_agec_sCT(1,nyr_agec_sCT);
init_vector nsamp_agec_sCT(1,nyr_agec_sCT);
init_vector nfish_agec_sCT(1,nyr_agec_sCT);
init_matrix obs_agec_sCT(1,nyr_agec_sCT,1,nages_agec);


//################Commercial handline fleet #######################################
// Comm HL CPUE
init_int styr_cpue_cHL;                                             
init_int endyr_cpue_cHL;                                            
init_vector obs_cpue_cHL(styr_cpue_cHL,endyr_cpue_cHL); //Observed CPUE
init_vector obs_cv_cpue_cHL(styr_cpue_cHL,endyr_cpue_cHL);  //CV of cpue

// Comm HL  Landings (1000 lb whole weight)
init_int styr_L_cHL;
init_int endyr_L_cHL;
init_vector obs_L_cHL(styr_L_cHL,endyr_L_cHL);
init_vector obs_cv_L_cHL(styr_L_cHL,endyr_L_cHL);

// Comm HL length Compositions (3 cm bins)
init_int nyr_lenc_cHL;
init_ivector yrs_lenc_cHL(1,nyr_lenc_cHL);
init_vector nsamp_lenc_cHL(1,nyr_lenc_cHL);
init_vector nfish_lenc_cHL(1,nyr_lenc_cHL);
init_matrix obs_lenc_cHL(1,nyr_lenc_cHL,1,nlenbins);

// Comm HL age compositions
init_int nyr_agec_cHL;
init_ivector yrs_agec_cHL(1,nyr_agec_cHL);
init_vector nsamp_agec_cHL(1,nyr_agec_cHL);
init_vector nfish_agec_cHL(1,nyr_agec_cHL);
init_matrix obs_agec_cHL(1,nyr_agec_cHL,1,nages_agec);

// Comm HL  Discards (1000 fish)
init_int styr_D_cHL;
init_int endyr_D_cHL;
init_vector obs_released_cHL(styr_D_cHL,endyr_D_cHL);
init_vector obs_cv_D_cHL(styr_D_cHL,endyr_D_cHL);

// Comm HL Discard length Compositions (1 cm bins??)
init_int nyr_lenc_cHL_D;
init_ivector yrs_lenc_cHL_D(1,nyr_lenc_cHL_D);
init_vector nsamp_lenc_cHL_D(1,nyr_lenc_cHL_D);
init_vector nfish_lenc_cHL_D(1,nyr_lenc_cHL_D);
init_matrix obs_lenc_cHL_D(1,nyr_lenc_cHL_D,1,nlenbins);
//################Commercial handline fleet #######################################

// Comm Other  Landings (1000 lb whole weight)
init_int styr_L_cOT;
init_int endyr_L_cOT;
init_vector obs_L_cOT(styr_L_cOT,endyr_L_cOT);
init_vector obs_cv_L_cOT(styr_L_cOT,endyr_L_cOT);

// Comm Other length Compositions (3 cm bins)
init_int nyr_lenc_cOT;
init_ivector yrs_lenc_cOT(1,nyr_lenc_cOT);
init_vector nsamp_lenc_cOT(1,nyr_lenc_cOT);
init_vector nfish_lenc_cOT(1,nyr_lenc_cOT);
init_matrix obs_lenc_cOT(1,nyr_lenc_cOT,1,nlenbins);


//################Commercial historic trawl fleet #######################################

// Comm Other  Landings (1000 lb whole weight)
init_int styr_L_cTW;
init_int endyr_L_cTW;
init_vector obs_L_cTW(styr_L_cTW,endyr_L_cTW);
init_vector obs_cv_L_cTW(styr_L_cTW,endyr_L_cTW);

//###################Headboat fleet ##########################################
// rHB CPUE
init_int styr_cpue_rHB;                                             
init_int endyr_cpue_rHB;                                            
init_vector obs_cpue_rHB(styr_cpue_rHB,endyr_cpue_rHB);//Observed CPUE
init_vector obs_cv_cpue_rHB(styr_cpue_rHB,endyr_cpue_rHB); //CV of cpue

// rHB Landings (1000s fish)
init_int styr_L_rHB;
init_int endyr_L_rHB;
init_vector obs_L_rHB(styr_L_rHB,endyr_L_rHB);   //vector of observed landings by year 
init_vector obs_cv_L_rHB(styr_L_rHB,endyr_L_rHB);    //vector of CV of landings by year

//rHB length compositions
//init_int nyr_rHB_lenc;
//init_ivector yrs_rHB_lenc(1,nyr_rHB_lenc);
//init_vector nsamp_rHB_lenc(1,nyr_rHB_lenc);
//init_vector nfish_rHB_lenc(1,nyr_rHB_lenc);
//init_matrix obs_rHB_lenc(1,nyr_rHB_lenc,1,nlenbins);

// rHB age compositions
init_int nyr_agec_rHB;
init_ivector yrs_agec_rHB(1,nyr_agec_rHB);
init_vector nsamp_agec_rHB(1,nyr_agec_rHB);
init_vector nfish_agec_rHB(1,nyr_agec_rHB);
init_matrix obs_agec_rHB(1,nyr_agec_rHB,1,nages_agec);

// rHB  Discards (1000 fish)
init_int styr_D_rHB;
init_int endyr_D_rHB;
init_vector obs_released_rHB(styr_D_rHB,endyr_D_rHB);
init_vector obs_cv_D_rHB(styr_D_rHB,endyr_D_rHB);

//rHB discards length compositions
init_int nyr_lenc_rHB_D;
init_ivector yrs_lenc_rHB_D(1,nyr_lenc_rHB_D);
init_vector nsamp_lenc_rHB_D(1,nyr_lenc_rHB_D);
init_vector nfish_lenc_rHB_D(1,nyr_lenc_rHB_D);
init_matrix obs_lenc_rHB_D(1,nyr_lenc_rHB_D,1,nlenbins);


//###################General Recreational fleet ##########################################
// // rGN CPUE
init_int styr_cpue_rGN;                                             
init_int endyr_cpue_rGN;                                            
init_vector obs_cpue_rGN(styr_cpue_rGN,endyr_cpue_rGN);//Observed CPUE
init_vector obs_cv_cpue_rGN(styr_cpue_rGN,endyr_cpue_rGN); //CV of cpue

// Landings (1000s fish)
init_int styr_L_rGN;
init_int endyr_L_rGN;
init_vector obs_L_rGN(styr_L_rGN,endyr_L_rGN);   //vector of observed landings by year 
init_vector obs_cv_L_rGN(styr_L_rGN,endyr_L_rGN);    //vector of CV of landings by year

//  Discards (1000 fish)
init_int styr_D_rGN;
init_int endyr_D_rGN;
init_vector obs_released_rGN(styr_D_rGN,endyr_D_rGN);
init_vector obs_cv_D_rGN(styr_D_rGN,endyr_D_rGN);

//rGN length compositions
//init_int nyr_rGN_lenc;
//init_ivector yrs_rGN_lenc(1,nyr_rGN_lenc);
//init_vector nsamp_rGN_lenc(1,nyr_rGN_lenc);
//init_vector nfish_rGN_lenc(1,nyr_rGN_lenc);
//init_matrix obs_rGN_lenc(1,nyr_rGN_lenc,1,nlenbins);

//  Age compositions
init_int nyr_agec_rGN;
init_ivector yrs_agec_rGN(1,nyr_agec_rGN);
init_vector nsamp_agec_rGN(1,nyr_agec_rGN);
init_vector nfish_agec_rGN(1,nyr_agec_rGN);
init_matrix obs_agec_rGN(1,nyr_agec_rGN,1,nages_agec);
//!!cout << nyr_agec_rGN << endl;


//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: parameter section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##################Single Parameter values and initial guesses #################################
// von Bert parms in TL mm population
init_vector set_Linf(1,7);
init_vector set_K(1,7);
init_vector set_t0(1,7);
init_vector set_len_cv(1,7);

//Scalar used only for computing MSST. For RG, this scalar is MSST = 0.75 SSBmsy
init_vector set_M_constant(1,7);     
//Spawner-recruit parameters (Initial guesses or fixed values)
init_vector set_steep(1,7);         //recruitment steepness
init_vector set_log_R0(1,7);        //recruitment R0
init_vector set_R_autocorr(1,7);    //recruitment autocorrelation
init_vector set_rec_sigma(1,7);     //recruitment standard deviation in log space

init_vector set_log_dm_lenc_cHL(1,7);    //Dirichlet-multinomial overdispersion parameter
//init_vector set_log_dm_cOT_lc(1,7);    //Dirichlet-multinomial overdispersion parameter
//init_vector set_log_dm_rHB_lc(1,7);    //Dirichlet-multinomial overdispersion parameter
//init_vector set_log_dm_rGN_lc(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_lenc_rHB_D(1,7);  //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_lenc_cHL_D(1,7);  //Dirichlet-multinomial overdispersion parameter
//init_vector set_log_dm_sFT_lc(1,7);   //Dirichlet-multinomial overdispersion parameter
////init_vector set_log_dm_sCT_lc(1,7);   //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_agec_cHL(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_agec_rHB(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_agec_rGN(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_agec_sCT(1,7);   //Dirichlet-multinomial overdispersion parameter

//cout << set_log_dm_agec_sCT << endl;
// //Initial guesses or fixed values of estimated selectivity parameters

//init_number set_selpar_A50_sFT; 
//init_number set_selpar_slope_sFT;
//init_number set_selpar_A502_sFT; 
//init_number set_selpar_slope2_sFT;

init_vector set_selpar_A50_sCT(1,7); //was logistic exponential in update
init_vector set_selpar_slope_sCT(1,7);
init_vector set_selpar_A502_sCT(1,7); 
init_vector set_selpar_slope2_sCT(1,7);

init_vector set_selpar_logit_Age1_sCT(1,7);
init_vector set_selpar_logit_Age2_sCT(1,7);

init_vector set_selpar_A50_cHL1(1,7);
init_vector set_selpar_slope_cHL1(1,7);
init_vector set_selpar_A502_cHL1(1,7);
init_vector set_selpar_slope2_cHL1(1,7);
init_vector set_selpar_A50_cHL2(1,7);
init_vector set_selpar_slope_cHL2(1,7);
init_vector set_selpar_A502_cHL2(1,7); 
init_vector set_selpar_slope2_cHL2(1,7);

init_vector set_selpar_A50_cOT1(1,7);
init_vector set_selpar_slope_cOT1(1,7);
init_vector set_selpar_A502_cOT1(1,7); 
init_vector set_selpar_slope2_cOT1(1,7);

//init_vector set_selpar_A50_cTW(1,7); 
//init_vector set_selpar_slope_cTW(1,7);

init_vector set_selpar_A50_rHB1(1,7);
init_vector set_selpar_slope_rHB1(1,7);
init_vector set_selpar_A502_rHB1(1,7);
init_vector set_selpar_slope2_rHB1(1,7);
init_vector set_selpar_A50_rHB2(1,7);
init_vector set_selpar_slope_rHB2(1,7);
init_vector set_selpar_A502_rHB2(1,7);
init_vector set_selpar_slope2_rHB2(1,7);
init_vector set_selpar_A50_rHB3(1,7); 
init_vector set_selpar_slope_rHB3(1,7);
init_vector set_selpar_A502_rHB3(1,7); 
init_vector set_selpar_slope2_rHB3(1,7);
init_vector set_selpar_A50_rHB4(1,7); 
init_vector set_selpar_slope_rHB4(1,7);
init_vector set_selpar_A502_rHB4(1,7); 
init_vector set_selpar_slope2_rHB4(1,7);

init_vector set_selpar_A50_rGN2(1,7);  //KC 
init_vector set_selpar_slope_rGN2(1,7); //KC
init_vector set_selpar_A50_rGN3(1,7); 
init_vector set_selpar_slope_rGN3(1,7);
init_vector set_selpar_A502_rGN3(1,7); 
init_vector set_selpar_slope2_rGN3(1,7);
init_vector set_selpar_A50_rGN4(1,7); 
init_vector set_selpar_slope_rGN4(1,7);
init_vector set_selpar_A502_rGN4(1,7); 
init_vector set_selpar_slope2_rGN4(1,7);

//Discard selectivities
init_vector set_sel_cHL_D_2(1,nages);
init_vector set_selpar_logit_Age2_cHL_D2(1,7);
init_vector set_selpar_A502_cHL_D3(1,7);

init_number set_selpar_Age1_cHL_D2;
init_number set_selpar_A50_cHL_D3;  
init_number set_selpar_slope_cHL_D3;
init_number set_selpar_slope2_cHL_D3;

init_vector set_sel_rHB_D_2(1,nages);
init_vector set_sel_rHB_D_3(1,nages);
init_number set_selpar_Age1_rHB_D3;
init_vector set_sel_rHB_D_4(1,nages);
init_number set_selpar_Age1_rHB_D4;

init_vector set_selpar_logit_Age2_rHB_D3(1,7);
init_vector set_selpar_logit_Age2_rHB_D4(1,7);

//init_vector set_selpar_age1logit_D(1,7);
//!!cout << set_selpar_logit_Age2_rHB_D4 << endl;

//--index catchability-----------------------------------------------------------------------------------
init_vector set_log_q_cpue_cHL(1,7);      //catchability coefficient (log) for cGN handline index
init_vector set_log_q_cpue_rHB(1,7);      //catchability coefficient (log) for headboat index
init_vector set_log_q_cpue_rGN(1,7);      //catchability coefficient (log) for rGN index
init_vector set_log_q_cpue_sFT(1,7);     //catchability coefficient (log) for SERFS Florida-snapper trapindex  
init_vector set_log_q_cpue_sCT(1,7);     //catchability coefficient (log) for SERFS MARMAP trap index  

//--mean F's in log space --------------------------------
init_vector set_log_avg_F_L_cHL(1,7);
init_vector set_log_avg_F_L_cOT(1,7);
init_vector set_log_avg_F_L_cTW(1,7);
init_vector set_log_avg_F_L_rHB(1,7);
init_vector set_log_avg_F_L_rGN(1,7);
init_vector set_log_avg_F_D_cHL(1,7);
init_vector set_log_avg_F_D_rHB(1,7);
init_vector set_log_avg_F_D_rGN(1,7);


//##################Dev Vector Parameter values (vals) and bounds #################################
//--F vectors---------------------------
init_vector set_log_dev_F_L_cHL(1,3);
init_vector set_log_dev_F_L_cOT(1,3); 
init_vector set_log_dev_F_L_cTW(1,3);
init_vector set_log_dev_F_L_rHB(1,3);
init_vector set_log_dev_F_L_rGN(1,3);
init_vector set_log_dev_F_D_cHL(1,3); 
init_vector set_log_dev_F_D_rHB(1,3);
init_vector set_log_dev_F_D_rGN(1,3);
init_vector set_log_dev_RWq(1,3);
init_vector set_log_dev_rec(1,3);
init_vector set_log_dev_Nage(1,3);

init_vector set_log_dev_vals_F_L_cHL(styr_L_cHL,endyr_L_cHL);
init_vector set_log_dev_vals_F_L_cOT(styr_L_cOT,endyr_L_cOT);
init_vector set_log_dev_vals_F_L_cTW(styr_L_cTW,endyr_L_cTW);
init_vector set_log_dev_vals_F_L_rHB(styr_L_rHB,endyr_L_rHB);
init_vector set_log_dev_vals_F_L_rGN(styr_L_rGN,endyr_L_rGN);
init_vector set_log_dev_vals_F_D_cHL(styr_D_cHL,endyr_D_cHL);
init_vector set_log_dev_vals_F_D_rHB(styr_D_rHB,endyr_D_rHB);
init_vector set_log_dev_vals_F_D_rGN(styr_D_rGN,endyr_D_rGN);
init_vector set_log_dev_vals_rec(styr_rec_dev,endyr_rec_dev);
init_vector set_log_dev_vals_Nage(2,nages);             

//!!cout << set_log_dev_vals_Nage << endl;


//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: likelihood weights section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

init_number set_w_L;
init_number set_w_D;            //weight for discards
init_number set_w_cpue_cHL;         //weight for cGN handline index
init_number set_w_cpue_rHB;         //weight for headboat index
init_number set_w_cpue_rGN;         //weight for rGN index
init_number set_w_cpue_sFT;        //weight for sFT index
init_number set_w_cpue_sCT;        //weight for sCT index
init_number set_w_lenc_cHL;        //weight for cGN handline len comps
init_number set_w_lenc_cOT;        //weight for cGN other len comps
//init_number set_w_lc_rHB;        //weight for headboat len comps
//init_number set_w_lc_rGN;		//weight for the recreational length comps
init_number set_w_lenc_rHB_D;		//weight for the headboat discard length comps
init_number set_w_lenc_cHL_D;		//weight for the headboat discard length comps   
init_number set_w_lenc_sFT;		//weight for the Florida-snapper trap length comps
//init_number set_w_lc_sCT;		//weight for the SERFS length comps  
init_number set_w_agec_cHL;        //weight for cGN handline age comps
init_number set_w_agec_rHB;        //weight for headboat age comps
init_number set_w_agec_rGN;		//weight for the Recreational age comps  
init_number set_w_agec_sCT;		//weight for the SERFS age comps 
init_number set_w_Nage_init;    //for fitting initial abundance at age (excluding first age)
init_number set_w_rec;          //for fitting S-R curve
init_number set_w_rec_early;    //additional constraint on early years recruitment
init_number set_w_rec_end;      //additional constraint on ending years recruitment 
init_number set_w_fullF;        //penalty for any Fapex>3(removed in final phase of optimization)
init_number set_w_Ftune;        //weight applied to tuning F (removed in final phase of optimization)

//!!cout << set_w_agec_sCT << endl;
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: miscellaneous stuff section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

//TL(mm)-weight(whole weight in kg) relationship: W=aL^b
init_number wgtpar_a;
init_number wgtpar_b;

//TL(mm)-FL(mm) relationship: FL=a+b*TL
init_number lenpar_TL2FL_a;
init_number lenpar_TL2FL_b;

///Maturity and proportion female at age
init_vector obs_maturity_f(1,nages);            //proportion females mature at age
//init_vector maturity_m_obs(1,nages);            //proportion males mature at age
init_vector obs_prop_f(1,nages);
init_number spawn_time_frac; //time of year of peak spawning, as a fraction of the year

//Fecundity-length relationship: Eggs/batch=a + bFL (length in FL); batches in number batches per yr
//kc init_number fecpar_batches;
init_vector fecpar_batches(1,nages); //kc convert to vector for age-based batch number
init_number fecpar_a;
init_number fecpar_b;
init_number fecpar_scale;    //could be used for scaling annual egg production (10^X eggs)

// Natural mortality
init_vector set_M(1,nages);     //age-dependent: used in model

//discard mortality constants
init_number set_Dmort_cHL;
init_number set_Dmort_rHB;
init_number set_Dmort_rGN;

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
init_number minSS_lenc_cHL;
init_number minSS_lenc_cOT;
//init_number minSS_rHB_lenc;
init_number minSS_lenc_rHB_D;
init_number minSS_lenc_cHL_D;
//init_number minSS_rGN_lenc;
init_number minSS_lenc_sFT;
//init_number minSS_sCT_lenc;

//threshold sample sizes for age comps
init_number minSS_agec_cHL;
init_number minSS_agec_rHB;
init_number minSS_agec_rGN;
init_number minSS_agec_sFT;
init_number minSS_agec_sCT;


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

//proportion of length comp mass below size limit considered when matching length comp
//note: these need both length comp and age comp data to be estimable
//kcinit_number set_p_lenc_cHL; 
//kcinit_number set_p_lenc_cOT;
//kcinit_number set_p_lenc_rHB2;
//kcinit_number set_p_lenc_rHB3;
//kcinit_number set_p_lenc_rHB4;
//kcinit_number set_p_lenc_rGN2;  
//kcinit_number set_p_lenc_rGN3;
//kcinit_number set_p_lenc_rGN4;

//kcinit_number set_p_lenc_cHL_D;
//kcinit_number set_p_lenc_rHB_D2;
//kcinit_number set_p_lenc_rHB_D3;
//kcinit_number set_p_lenc_rHB_D4;
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
  const double Linf_LO=set_Linf(2); 
  const double Linf_HI=set_Linf(3); 
  const double Linf_PH=set_Linf(4);
  const double K_LO=set_K(2); 
  const double K_HI=set_K(3); 
  const double K_PH=set_K(4);
  const double t0_LO=set_t0(2); 
  const double t0_HI=set_t0(3); 
  const double t0_PH=set_t0(4);  
  const double len_cv_LO=set_len_cv(2); 
  const double len_cv_HI=set_len_cv(3); 
  const double len_cv_PH=set_len_cv(4); 
     
  const double M_constant_LO=set_M_constant(2); 
  const double M_constant_HI=set_M_constant(3); 
  const double M_constant_PH=set_M_constant(4);        
  const double steep_LO=set_steep(2); 
  const double steep_HI=set_steep(3); 
  const double steep_PH=set_steep(4);
  const double log_R0_LO=set_log_R0(2); 
  const double log_R0_HI=set_log_R0(3); 
  const double log_R0_PH=set_log_R0(4);
  const double R_autocorr_LO=set_R_autocorr(2); 
  const double R_autocorr_HI=set_R_autocorr(3); 
  const double R_autocorr_PH=set_R_autocorr(4);
  const double rec_sigma_LO=set_rec_sigma(2); 
  const double rec_sigma_HI=set_rec_sigma(3); 
  const double rec_sigma_PH=set_rec_sigma(4);
  
  const double log_dm_cHL_lc_LO=set_log_dm_lenc_cHL(2); 
  const double log_dm_cHL_lc_HI=set_log_dm_lenc_cHL(3); 
  const double log_dm_cHL_lc_PH=set_log_dm_lenc_cHL(4);
  //const double log_dm_cOT_lc_LO=set_log_dm_cOT_lc(2); 
  //const double log_dm_cOT_lc_HI=set_log_dm_cOT_lc(3); 
  //const double log_dm_cOT_lc_PH=set_log_dm_cOT_lc(4);
  //const double log_dm_rHB_lc_LO=set_log_dm_rHB_lc(2); 
  //const double log_dm_rHB_lc_HI=set_log_dm_rHB_lc(3); 
  //const double log_dm_rHB_lc_PH=set_log_dm_rHB_lc(4);
  //const double log_dm_rGN_lc_LO=set_log_dm_rGN_lc(2); 
  //const double log_dm_rGN_lc_HI=set_log_dm_rGN_lc(3); 
  //const double log_dm_rGN_lc_PH=set_log_dm_rGN_lc(4);
  const double log_dm_rHB_D_lc_LO=set_log_dm_lenc_rHB_D(2); 
  const double log_dm_rHB_D_lc_HI=set_log_dm_lenc_rHB_D(3); 
  const double log_dm_rHB_D_lc_PH=set_log_dm_lenc_rHB_D(4);
  const double log_dm_cHL_D_lc_LO=set_log_dm_lenc_cHL_D(2); 
  const double log_dm_cHL_D_lc_HI=set_log_dm_lenc_cHL_D(3); 
  const double log_dm_cHL_D_lc_PH=set_log_dm_lenc_cHL_D(4);
  //const double log_dm_sFT_lc_LO=set_log_dm_sFT_lc(2); 
  //const double log_dm_sFT_lc_HI=set_log_dm_sFT_lc(3); 
  //const double log_dm_sFT_lc_PH=set_log_dm_sFT_lc(4);
  ////const double log_dm_sCT_lc_LO=set_log_dm_sCT_lc(2); 
  ////const double log_dm_sCT_lc_HI=set_log_dm_sCT_lc(3); 
  ////const double log_dm_sCT_lc_PH=set_log_dm_sCT_lc(4);
  const double log_dm_cHL_ac_LO=set_log_dm_agec_cHL(2); 
  const double log_dm_cHL_ac_HI=set_log_dm_agec_cHL(3); 
  const double log_dm_cHL_ac_PH=set_log_dm_agec_cHL(4);
  const double log_dm_rHB_ac_LO=set_log_dm_agec_rHB(2); 
  const double log_dm_rHB_ac_HI=set_log_dm_agec_rHB(3); 
  const double log_dm_rHB_ac_PH=set_log_dm_agec_rHB(4);
  const double log_dm_rGN_ac_LO=set_log_dm_agec_rGN(2); 
  const double log_dm_rGN_ac_HI=set_log_dm_agec_rGN(3); 
  const double log_dm_rGN_ac_PH=set_log_dm_agec_rGN(4);
  const double log_dm_sCT_ac_LO=set_log_dm_agec_sCT(2); 
  const double log_dm_sCT_ac_HI=set_log_dm_agec_sCT(3); 
  const double log_dm_sCT_ac_PH=set_log_dm_agec_sCT(4);
  
  const double selpar_A50_cHL1_LO=set_selpar_A50_cHL1(2); 
  const double selpar_A50_cHL1_HI=set_selpar_A50_cHL1(3); 
  const double selpar_A50_cHL1_PH=set_selpar_A50_cHL1(4);
  const double selpar_slope_cHL1_LO=set_selpar_slope_cHL1(2); 
  const double selpar_slope_cHL1_HI=set_selpar_slope_cHL1(3); 
  const double selpar_slope_cHL1_PH=set_selpar_slope_cHL1(4);
  const double selpar_A502_cHL1_LO=set_selpar_A502_cHL1(2); 
  const double selpar_A502_cHL1_HI=set_selpar_A502_cHL1(3); 
  const double selpar_A502_cHL1_PH=set_selpar_A502_cHL1(4);
  const double selpar_slope2_cHL1_LO=set_selpar_slope2_cHL1(2); 
  const double selpar_slope2_cHL1_HI=set_selpar_slope2_cHL1(3); 
  const double selpar_slope2_cHL1_PH=set_selpar_slope2_cHL1(4);
  const double selpar_A50_cHL2_LO=set_selpar_A50_cHL2(2); 
  const double selpar_A50_cHL2_HI=set_selpar_A50_cHL2(3); 
  const double selpar_A50_cHL2_PH=set_selpar_A50_cHL2(4);
  const double selpar_slope_cHL2_LO=set_selpar_slope_cHL2(2); 
  const double selpar_slope_cHL2_HI=set_selpar_slope_cHL2(3); 
  const double selpar_slope_cHL2_PH=set_selpar_slope_cHL2(4);
  const double selpar_A502_cHL2_LO=set_selpar_A502_cHL2(2); 
  const double selpar_A502_cHL2_HI=set_selpar_A502_cHL2(3); 
  const double selpar_A502_cHL2_PH=set_selpar_A502_cHL2(4);
  const double selpar_slope2_cHL2_LO=set_selpar_slope2_cHL2(2); 
  const double selpar_slope2_cHL2_HI=set_selpar_slope2_cHL2(3); 
  const double selpar_slope2_cHL2_PH=set_selpar_slope2_cHL2(4);

  const double selpar_A50_cOT1_LO=set_selpar_A50_cOT1(2); 
  const double selpar_A50_cOT1_HI=set_selpar_A50_cOT1(3); 
  const double selpar_A50_cOT1_PH=set_selpar_A50_cOT1(4);
  const double selpar_slope_cOT1_LO=set_selpar_slope_cOT1(2); 
  const double selpar_slope_cOT1_HI=set_selpar_slope_cOT1(3); 
  const double selpar_slope_cOT1_PH=set_selpar_slope_cOT1(4);
  const double selpar_A502_cOT1_LO=set_selpar_A502_cOT1(2); 
  const double selpar_A502_cOT1_HI=set_selpar_A502_cOT1(3); 
  const double selpar_A502_cOT1_PH=set_selpar_A502_cOT1(4);
  const double selpar_slope2_cOT1_LO=set_selpar_slope2_cOT1(2); 
  const double selpar_slope2_cOT1_HI=set_selpar_slope2_cOT1(3); 
  const double selpar_slope2_cOT1_PH=set_selpar_slope2_cOT1(4);

  const double selpar_A50_rHB1_LO=set_selpar_A50_rHB1(2); 
  const double selpar_A50_rHB1_HI=set_selpar_A50_rHB1(3); 
  const double selpar_A50_rHB1_PH=set_selpar_A50_rHB1(4);
  const double selpar_slope_rHB1_LO=set_selpar_slope_rHB1(2); 
  const double selpar_slope_rHB1_HI=set_selpar_slope_rHB1(3); 
  const double selpar_slope_rHB1_PH=set_selpar_slope_rHB1(4);
  const double selpar_A502_rHB1_LO=set_selpar_A502_rHB1(2); 
  const double selpar_A502_rHB1_HI=set_selpar_A502_rHB1(3); 
  const double selpar_A502_rHB1_PH=set_selpar_A502_rHB1(4);
  const double selpar_slope2_rHB1_LO=set_selpar_slope2_rHB1(2); 
  const double selpar_slope2_rHB1_HI=set_selpar_slope2_rHB1(3); 
  const double selpar_slope2_rHB1_PH=set_selpar_slope2_rHB1(4);
  const double selpar_A50_rHB2_LO=set_selpar_A50_rHB2(2); 
  const double selpar_A50_rHB2_HI=set_selpar_A50_rHB2(3); 
  const double selpar_A50_rHB2_PH=set_selpar_A50_rHB2(4);
  const double selpar_slope_rHB2_LO=set_selpar_slope_rHB2(2); 
  const double selpar_slope_rHB2_HI=set_selpar_slope_rHB2(3); 
  const double selpar_slope_rHB2_PH=set_selpar_slope_rHB2(4);
  const double selpar_A502_rHB2_LO=set_selpar_A502_rHB2(2); 
  const double selpar_A502_rHB2_HI=set_selpar_A502_rHB2(3); 
  const double selpar_A502_rHB2_PH=set_selpar_A502_rHB2(4);
  const double selpar_slope2_rHB2_LO=set_selpar_slope2_rHB2(2); 
  const double selpar_slope2_rHB2_HI=set_selpar_slope2_rHB2(3); 
  const double selpar_slope2_rHB2_PH=set_selpar_slope2_rHB2(4);
  const double selpar_A50_rHB3_LO=set_selpar_A50_rHB3(2); 
  const double selpar_A50_rHB3_HI=set_selpar_A50_rHB3(3); 
  const double selpar_A50_rHB3_PH=set_selpar_A50_rHB3(4);
  const double selpar_slope_rHB3_LO=set_selpar_slope_rHB3(2); 
  const double selpar_slope_rHB3_HI=set_selpar_slope_rHB3(3); 
  const double selpar_slope_rHB3_PH=set_selpar_slope_rHB3(4);
  const double selpar_A502_rHB3_LO=set_selpar_A502_rHB3(2); 
  const double selpar_A502_rHB3_HI=set_selpar_A502_rHB3(3); 
  const double selpar_A502_rHB3_PH=set_selpar_A502_rHB3(4);
  const double selpar_slope2_rHB3_LO=set_selpar_slope2_rHB3(2); 
  const double selpar_slope2_rHB3_HI=set_selpar_slope2_rHB3(3); 
  const double selpar_slope2_rHB3_PH=set_selpar_slope2_rHB3(4);
  const double selpar_A50_rHB4_LO=set_selpar_A50_rHB4(2); 
  const double selpar_A50_rHB4_HI=set_selpar_A50_rHB4(3); 
  const double selpar_A50_rHB4_PH=set_selpar_A50_rHB4(4);
  const double selpar_slope_rHB4_LO=set_selpar_slope_rHB4(2); 
  const double selpar_slope_rHB4_HI=set_selpar_slope_rHB4(3); 
  const double selpar_slope_rHB4_PH=set_selpar_slope_rHB4(4);
  const double selpar_A502_rHB4_LO=set_selpar_A502_rHB4(2); 
  const double selpar_A502_rHB4_HI=set_selpar_A502_rHB4(3); 
  const double selpar_A502_rHB4_PH=set_selpar_A502_rHB4(4);
  const double selpar_slope2_rHB4_LO=set_selpar_slope2_rHB4(2); 
  const double selpar_slope2_rHB4_HI=set_selpar_slope2_rHB4(3); 
  const double selpar_slope2_rHB4_PH=set_selpar_slope2_rHB4(4);

  const double selpar_A50_rGN2_LO=set_selpar_A50_rGN2(2); 
  const double selpar_A50_rGN2_HI=set_selpar_A50_rGN2(3); 
  const double selpar_A50_rGN2_PH=set_selpar_A50_rGN2(4);
  const double selpar_slope_rGN2_LO=set_selpar_slope_rGN2(2); 
  const double selpar_slope_rGN2_HI=set_selpar_slope_rGN2(3); 
  const double selpar_slope_rGN2_PH=set_selpar_slope_rGN2(4);
  const double selpar_A50_rGN3_LO=set_selpar_A50_rGN3(2); 
  const double selpar_A50_rGN3_HI=set_selpar_A50_rGN3(3); 
  const double selpar_A50_rGN3_PH=set_selpar_A50_rGN3(4);
  const double selpar_slope_rGN3_LO=set_selpar_slope_rGN3(2); 
  const double selpar_slope_rGN3_HI=set_selpar_slope_rGN3(3); 
  const double selpar_slope_rGN3_PH=set_selpar_slope_rGN3(4);
  const double selpar_A502_rGN3_LO=set_selpar_A502_rGN3(2); 
  const double selpar_A502_rGN3_HI=set_selpar_A502_rGN3(3); 
  const double selpar_A502_rGN3_PH=set_selpar_A502_rGN3(4);
  const double selpar_slope2_rGN3_LO=set_selpar_slope2_rGN3(2); 
  const double selpar_slope2_rGN3_HI=set_selpar_slope2_rGN3(3); 
  const double selpar_slope2_rGN3_PH=set_selpar_slope2_rGN3(4);
  const double selpar_A50_rGN4_LO=set_selpar_A50_rGN4(2); 
  const double selpar_A50_rGN4_HI=set_selpar_A50_rGN4(3); 
  const double selpar_A50_rGN4_PH=set_selpar_A50_rGN4(4);
  const double selpar_slope_rGN4_LO=set_selpar_slope_rGN4(2); 
  const double selpar_slope_rGN4_HI=set_selpar_slope_rGN4(3); 
  const double selpar_slope_rGN4_PH=set_selpar_slope_rGN4(4);
  const double selpar_A502_rGN4_LO=set_selpar_A502_rGN4(2); 
  const double selpar_A502_rGN4_HI=set_selpar_A502_rGN4(3); 
  const double selpar_A502_rGN4_PH=set_selpar_A502_rGN4(4);
  const double selpar_slope2_rGN4_LO=set_selpar_slope2_rGN4(2); 
  const double selpar_slope2_rGN4_HI=set_selpar_slope2_rGN4(3); 
  const double selpar_slope2_rGN4_PH=set_selpar_slope2_rGN4(4);

  const double selpar_A50_sCT_LO=set_selpar_A50_sCT(2); 
  const double selpar_A50_sCT_HI=set_selpar_A50_sCT(3); 
  const double selpar_A50_sCT_PH=set_selpar_A50_sCT(4);
  const double selpar_slope_sCT_LO=set_selpar_slope_sCT(2); 
  const double selpar_slope_sCT_HI=set_selpar_slope_sCT(3); 
  const double selpar_slope_sCT_PH=set_selpar_slope_sCT(4);
  const double selpar_A502_sCT_LO=set_selpar_A502_sCT(2); 
  const double selpar_A502_sCT_HI=set_selpar_A502_sCT(3); 
  const double selpar_A502_sCT_PH=set_selpar_A502_sCT(4);
  const double selpar_slope2_sCT_LO=set_selpar_slope2_sCT(2); 
  const double selpar_slope2_sCT_HI=set_selpar_slope2_sCT(3); 
  const double selpar_slope2_sCT_PH=set_selpar_slope2_sCT(4);
  const double selpar_Age1_sCT_logit_LO=set_selpar_logit_Age1_sCT(2); 
  const double selpar_Age1_sCT_logit_HI=set_selpar_logit_Age1_sCT(3); 
  const double selpar_Age1_sCT_logit_PH=set_selpar_logit_Age1_sCT(4);
  const double selpar_Age2_sCT_logit_LO=set_selpar_logit_Age2_sCT(2); 
  const double selpar_Age2_sCT_logit_HI=set_selpar_logit_Age2_sCT(3); 
  const double selpar_Age2_sCT_logit_PH=set_selpar_logit_Age2_sCT(4);
  
  const double selpar_Age2_cHL_D2_logit_LO=set_selpar_logit_Age2_cHL_D2(2); 
  const double selpar_Age2_cHL_D2_logit_HI=set_selpar_logit_Age2_cHL_D2(3); 
  const double selpar_Age2_cHL_D2_logit_PH=set_selpar_logit_Age2_cHL_D2(4);
  
  const double selpar_A502_cHL_D3_LO=set_selpar_A502_cHL_D3(2); 
  const double selpar_A502_cHL_D3_HI=set_selpar_A502_cHL_D3(3); 
  const double selpar_A502_cHL_D3_PH=set_selpar_A502_cHL_D3(4);
  
  const double selpar_Age2_rHB_D3_logit_LO=set_selpar_logit_Age2_rHB_D3(2); 
  const double selpar_Age2_rHB_D3_logit_HI=set_selpar_logit_Age2_rHB_D3(3); 
  const double selpar_Age2_rHB_D3_logit_PH=set_selpar_logit_Age2_rHB_D3(4);
  
  const double selpar_Age2_rHB_D4_logit_LO=set_selpar_logit_Age2_rHB_D4(2); 
  const double selpar_Age2_rHB_D4_logit_HI=set_selpar_logit_Age2_rHB_D4(3); 
  const double selpar_Age2_rHB_D4_logit_PH=set_selpar_logit_Age2_rHB_D4(4);
  
  
  const double log_q_cHL_LO=set_log_q_cpue_cHL(2); const double log_q_cHL_HI=set_log_q_cpue_cHL(3); const double log_q_cHL_PH=set_log_q_cpue_cHL(4);
  const double log_q_rHB_LO=set_log_q_cpue_rHB(2); const double log_q_rHB_HI=set_log_q_cpue_rHB(3); const double log_q_rHB_PH=set_log_q_cpue_rHB(4);
  const double log_q_rGN_LO=set_log_q_cpue_rGN(2); const double log_q_rGN_HI=set_log_q_cpue_rGN(3); const double log_q_rGN_PH=set_log_q_cpue_rGN(4);
  const double log_q_sFT_LO=set_log_q_cpue_sFT(2); const double log_q_sFT_HI=set_log_q_cpue_sFT(3); const double log_q_sFT_PH=set_log_q_cpue_sFT(4);
  const double log_q_sCT_LO=set_log_q_cpue_sCT(2); const double log_q_sCT_HI=set_log_q_cpue_sCT(3); const double log_q_sCT_PH=set_log_q_cpue_sCT(4);
  const double log_avg_F_cHL_LO=set_log_avg_F_L_cHL(2); const double log_avg_F_cHL_HI=set_log_avg_F_L_cHL(3); const double log_avg_F_cHL_PH=set_log_avg_F_L_cHL(4);
  const double log_avg_F_cOT_LO=set_log_avg_F_L_cOT(2); const double log_avg_F_cOT_HI=set_log_avg_F_L_cOT(3); const double log_avg_F_cOT_PH=set_log_avg_F_L_cOT(4);
  const double log_avg_F_cTW_LO=set_log_avg_F_L_cTW(2); const double log_avg_F_cTW_HI=set_log_avg_F_L_cTW(3); const double log_avg_F_cTW_PH=set_log_avg_F_L_cTW(4);
  const double log_avg_F_rHB_LO=set_log_avg_F_L_rHB(2); const double log_avg_F_rHB_HI=set_log_avg_F_L_rHB(3); const double log_avg_F_rHB_PH=set_log_avg_F_L_rHB(4); 
  const double log_avg_F_rGN_LO=set_log_avg_F_L_rGN(2); const double log_avg_F_rGN_HI=set_log_avg_F_L_rGN(3); const double log_avg_F_rGN_PH=set_log_avg_F_L_rGN(4); 
  const double log_avg_F_cHL_D_LO=set_log_avg_F_D_cHL(2); const double log_avg_F_cHL_D_HI=set_log_avg_F_D_cHL(3); const double log_avg_F_cHL_D_PH=set_log_avg_F_D_cHL(4);
  const double log_avg_F_rHB_D_LO=set_log_avg_F_D_rHB(2); const double log_avg_F_rHB_D_HI=set_log_avg_F_D_rHB(3); const double log_avg_F_rHB_D_PH=set_log_avg_F_D_rHB(4); 
  const double log_avg_F_rGN_D_LO=set_log_avg_F_D_rGN(2); const double log_avg_F_rGN_D_HI=set_log_avg_F_D_rGN(3); const double log_avg_F_rGN_D_PH=set_log_avg_F_D_rGN(4); 
  
  //-dev vectors-----------------------------------------------------------------------------------------------------------  
  const double log_F_dev_cHL_LO=set_log_dev_F_L_cHL(1); const double log_F_dev_cHL_HI=set_log_dev_F_L_cHL(2); const double log_F_dev_cHL_PH=set_log_dev_F_L_cHL(3);  
  const double log_F_dev_cOT_LO=set_log_dev_F_L_cOT(1); const double log_F_dev_cOT_HI=set_log_dev_F_L_cOT(2); const double log_F_dev_cOT_PH=set_log_dev_F_L_cOT(3);
  const double log_F_dev_cTW_LO=set_log_dev_F_L_cTW(1); const double log_F_dev_cTW_HI=set_log_dev_F_L_cTW(2); const double log_F_dev_cTW_PH=set_log_dev_F_L_cTW(3);   
  const double log_F_dev_rHB_LO=set_log_dev_F_L_rHB(1); const double log_F_dev_rHB_HI=set_log_dev_F_L_rHB(2); const double log_F_dev_rHB_PH=set_log_dev_F_L_rHB(3);   
  const double log_F_dev_rGN_LO=set_log_dev_F_L_rGN(1); const double log_F_dev_rGN_HI=set_log_dev_F_L_rGN(2); const double log_F_dev_rGN_PH=set_log_dev_F_L_rGN(3);   
  
  const double log_F_dev_cHL_D_LO=set_log_dev_F_D_cHL(1); const double log_F_dev_cHL_D_HI=set_log_dev_F_D_cHL(2); const double log_F_dev_cHL_D_PH=set_log_dev_F_D_cHL(3);   
  const double log_F_dev_rHB_D_LO=set_log_dev_F_D_rHB(1); const double log_F_dev_rHB_D_HI=set_log_dev_F_D_rHB(2); const double log_F_dev_rHB_D_PH=set_log_dev_F_D_rHB(3);   
  const double log_F_dev_rGN_D_LO=set_log_dev_F_D_rGN(1); const double log_F_dev_rGN_D_HI=set_log_dev_F_D_rGN(2); const double log_F_dev_rGN_D_PH=set_log_dev_F_D_rGN(3);   
  
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

  vector meanlen_TL(1,nages);   //mean total length (mm) at age all fish
  vector meanlen_FL(1,nages);        //mean Fork length (mm) at age
  vector fecundity(1,nages);         // annual fec at age per mat female: fecpar_batches*(fecpar_a+fecpar_bFL)
  
  vector wgt_g(1,nages);        //whole wgt in g
  vector wgt_kg(1,nages);       //whole wgt in kg
  vector wgt_mt(1,nages);       //whole wgt in mt
  vector wgt_klb(1,nages);      //whole wgt in 1000 lb
  vector wgt_lb(1,nages);       //whole wgt in lb  
    
    
  matrix len_cHL_mm(styr,endyr,1,nages);          //mean length at age of commercial handline landings in mm 
  matrix wholewgt_cHL_klb(styr,endyr,1,nages);    //whole wgt of commercial handline landings in 1000 lb   
  matrix len_cOT_mm(styr,endyr,1,nages);          //mean length at age of commercial other landings in mm
  matrix wholewgt_cOT_klb(styr,endyr,1,nages);    //whole wgt of commercial other landings in 1000 lb
  matrix len_cTW_mm(styr,endyr,1,nages);          //mean length at age of commercial other landings in mm 
  matrix wholewgt_cTW_klb(styr,endyr,1,nages);    //whole wgt of commercial other landings in 1000 lb 
  matrix len_rHB_mm(styr,endyr,1,nages);          //mean length at age of rHB landings in mm     
  matrix wholewgt_rHB_klb(styr,endyr,1,nages);    //whole wgt of rHB landings in 1000 lb  
  matrix len_rGN_mm(styr,endyr,1,nages);          //mean length at age of rGN landings in mm     
  matrix wholewgt_rGN_klb(styr,endyr,1,nages);    //whole wgt of rGN landings in 1000 lb    

  matrix len_cHL_D_mm(styr,endyr,1,nages);          //mean length at age of commercial handline discards in mm 
  matrix wholewgt_cHL_D_klb(styr,endyr,1,nages);    //whole wgt of commercial handline discards in 1000 lb   
  matrix len_rHB_D_mm(styr,endyr,1,nages);          //mean length at age of rHB discards in mm     
  matrix wholewgt_rHB_D_klb(styr,endyr,1,nages);    //whole wgt of rHB discards in 1000 lb  
  matrix len_rGN_D_mm(styr,endyr,1,nages);          //mean length at age of rGN discards in mm     
  matrix wholewgt_rGN_D_klb(styr,endyr,1,nages);    //whole wgt of rGN discards in 1000 lb  
  
  matrix lenprob(1,nages,1,nlenbins);           //distn of size at age (age-length key, 3 cm bins) in population
  //kc matrix lenprob_plus(1,nages,1,nlenbins_plus); //used to compute mass in last length bin (a plus group)   
  number zscore_len;                            //standardized normal values used for computing lenprob
  vector cprob_lenvec(1,nlenbins);              //cumulative probabilities used for computing lenprob
  number zscore_lzero;                          //standardized normal values for length = 0
  number cprob_lzero;                           //length probability mass below zero, used for computing lenprob
    
  //matrices below are used to match length comps
  matrix lenprob_cHL(1,nages,1,nlenbins);   //distn of size at age in cHL   kc added
  matrix lenprob_cOT(1,nages,1,nlenbins);   //distn of size at age in cOT   kc added
  //matrix lenprob_rHB(1,nages,1,nlenbins);   //distn of size at age in rHB   kc added
  matrix lenprob_rHB_D(1,nages,1,nlenbins); //distn of size at age in rHB_D  kc added
 // matrix lenprob_rGN(1,nages,1,nlenbins);   //distn of size at age in rGN   kc added
  matrix lenprob_cHL_D(1,nages,1,nlenbins);   //distn of size at age in cHL_D   kc added
  matrix lenprob_sFT(1,nages,1,nlenbins); //distn of size at age in sFT   kc added
  //matrix lenprob_sCT(1,nages,1,nlenbins); //distn of size at age in sCT   kc added

  //matrix lenprob_cHL1(1,nages,1,nlenbins);     //distn of size at age in cHL  kc
  //matrix lenprob_cHL2(1,nages,1,nlenbins);     //distn of size at age in cHL  kc
  //matrix lenprob_cOT1(1,nages,1,nlenbins);     //distn of size at age in cOT  kc
  //matrix lenprob_cOT2(1,nages,1,nlenbins);     //distn of size at age in cOT  kc
  //matrix lenprob_rHB1(1,nages,1,nlenbins);     //distn of size at age in rHB   kc
  //matrix lenprob_rHB2(1,nages,1,nlenbins);     //distn of size at age in rHB   kc
  //matrix lenprob_rHB3(1,nages,1,nlenbins);     //distn of size at age in rHB   kc
  //matrix lenprob_rHB4(1,nages,1,nlenbins);     //distn of size at age in rHB   kc
  //matrix lenprob_rHB_D2(1,nages,1,nlenbins);   //distn of size at age in rHB discards kc  
  //matrix lenprob_rHB_D3(1,nages,1,nlenbins);   //distn of size at age in rHB discards kc 
  //matrix lenprob_rHB_D4(1,nages,1,nlenbins);   //distn of size at age in rHB discards kc 
  //matrix lenprob_rGN1(1,nages,1,nlenbins);     //distn of size at age in rGN  kc
  //matrix lenprob_rGN2(1,nages,1,nlenbins);     //distn of size at age in rGN  kc
  //matrix lenprob_rGN3(1,nages,1,nlenbins);     //distn of size at age in rGN  kc
  //matrix lenprob_rGN4(1,nages,1,nlenbins);     //distn of size at age in rGN  kc
  //matrix lenprob_cHL_D2(1,nages,1,nlenbins);   //distn of size at age in rHB discards kc  
  //matrix lenprob_cHL_D3(1,nages,1,nlenbins);   //distn of size at age in rHB discards  kc
  //matrix lenprob_sFT(1,nages,1,nlenbins);    //distn of size at age in sFT  kc
  //matrix lenprob_sCT(1,nages,1,nlenbins);    //distn of size at age in sCT  kc
  
  //init_bounded_dev_vector log_len_cv_dev(1,nages,-2,2,3)
  vector len_sd(1,nages);
  vector len_cv(1,nages); //for fishgraph 

//----Predicted length and age compositions
  matrix pred_cHL_lenc(1,nyr_lenc_cHL,1,nlenbins);
  matrix pred_cOT_lenc(1,nyr_lenc_cOT,1,nlenbins);
  //matrix pred_rHB_lenc(1,nyr_rHB_lenc,1,nlenbins);  
  matrix pred_rHB_D_lenc(1,nyr_lenc_rHB_D,1,nlenbins); 
  matrix pred_cHL_D_lenc(1,nyr_lenc_cHL_D,1,nlenbins);  
  //matrix pred_rGN_lenc(1,nyr_rGN_lenc,1,nlenbins);
  matrix pred_sFT_lenc(1,nyr_lenc_sFT,1,nlenbins); 
  //matrix pred_sCT_lenc(1,nyr_sCT_lenc,1,nlenbins);  
  
  matrix pred_cHL_agec(1,nyr_agec_cHL,1,nages_agec);
  matrix pred_cHL_agec_allages(1,nyr_agec_cHL,1,nages);
  matrix ErrorFree_cHL_agec(1,nyr_agec_cHL,1,nages);    
  matrix pred_rHB_agec(1,nyr_agec_rHB,1,nages_agec);
  matrix pred_rHB_agec_allages(1,nyr_agec_rHB,1,nages);  
  matrix ErrorFree_rHB_agec(1,nyr_agec_rHB,1,nages);
  matrix pred_rGN_agec(1,nyr_agec_rGN,1,nages_agec);  
  matrix pred_rGN_agec_allages(1,nyr_agec_rGN,1,nages);  
  matrix ErrorFree_rGN_agec(1,nyr_agec_rGN,1,nages);  
  matrix pred_sCT_agec(1,nyr_agec_sCT,1,nages_agec);  
  matrix pred_sCT_agec_allages(1,nyr_agec_sCT,1,nages);  
  matrix ErrorFree_sCT_agec(1,nyr_agec_sCT,1,nages);  
 
 //p_lenc_fishery pars require age comp and length comp data for estimation
  //kc number p_lenc_cHL;
  //kc number p_lenc_cOT;
  //kc number p_lenc_rHB2;
  //kc number p_lenc_rHB3;
  //kc number p_lenc_rHB4;
  //kc number p_lenc_rGN2;
  //kc number p_lenc_rGN3;
  //kc number p_lenc_rGN4;
    
  //kc number p_lenc_cHL_D;
  //kc number p_lenc_rHB_D2; //no comp data in this period, this par only used for avg weight
  //kc number p_lenc_rHB_D3;
  //kcnumber p_lenc_rHB_D4;
  
//Sample size (perhaps adjusted herein) used in fitting comp data
  vector nsamp_cHL_lenc_allyr(styr,endyr);
  vector nsamp_cOT_lenc_allyr(styr,endyr);
  //vector nsamp_rHB_lenc_allyr(styr,endyr);
  vector nsamp_rHB_D_lenc_allyr(styr,endyr);  
  vector nsamp_cHL_D_lenc_allyr(styr,endyr);
  //vector nsamp_rGN_lenc_allyr(styr,endyr);
  vector nsamp_sFT_lenc_allyr(styr,endyr);
  //vector nsamp_sCT_lenc_allyr(styr,endyr);
  
  vector nsamp_cHL_agec_allyr(styr,endyr);
  vector nsamp_rHB_agec_allyr(styr,endyr);
  vector nsamp_rGN_agec_allyr(styr,endyr);  
  vector nsamp_sCT_agec_allyr(styr,endyr);  
  
//Nfish used in MCB analysis (not used in fitting)
  vector nfish_cHL_lenc_allyr(styr,endyr);
  vector nfish_cOT_lenc_allyr(styr,endyr);
  //vector nfish_rHB_lenc_allyr(styr,endyr);
  vector nfish_rHB_D_lenc_allyr(styr,endyr);  
  vector nfish_cHL_D_lenc_allyr(styr,endyr);  
  //vector nfish_rGN_lenc_allyr(styr,endyr);
  vector nfish_sFT_lenc_allyr(styr,endyr);
 // vector nfish_sCT_lenc_allyr(styr,endyr);
  
  vector nfish_cHL_agec_allyr(styr,endyr);
  vector nfish_rHB_agec_allyr(styr,endyr);
  vector nfish_rGN_agec_allyr(styr,endyr);  
  vector nfish_sCT_agec_allyr(styr,endyr);  

//Computed effective sample size for output (not used in fitting)
  vector neff_cHL_lenc_allyr(styr,endyr);
  vector neff_cOT_lenc_allyr(styr,endyr);
  //vector neff_rHB_lenc_allyr(styr,endyr);
  vector neff_rHB_D_lenc_allyr(styr,endyr);  
  vector neff_cHL_D_lenc_allyr(styr,endyr);  
  //vector neff_rGN_lenc_allyr(styr,endyr);
  vector neff_sFT_lenc_allyr(styr,endyr);
  //vector neff_sCT_lenc_allyr(styr,endyr);
  
  vector neff_cHL_agec_allyr(styr,endyr);
  vector neff_rHB_agec_allyr(styr,endyr);
  vector neff_rGN_agec_allyr(styr,endyr);  
  vector neff_sCT_agec_allyr(styr,endyr);  

  //kc number age_limit_10in; //age corresponding to 10-inch size limit, given mean growth curve
  //kc number age_limit_11in; //age corresponding to 11-inch size limit, given mean growth curve
  //kc number age_limit_12in; //age corresponding to 12-inch size limit, given mean growth curve
  
//-----Population-----------------------------------------------------------------------------------
  matrix N(styr,endyr+1,1,nages);           //Population numbers by year and age at start of yr
  matrix N_mdyr(styr,endyr,1,nages);        //Population numbers by year and age at mdpt of yr: used for comps and cpue
  matrix N_spawn(styr,endyr,1,nages);       //Population numbers by year and age at peaking spawning: used for SSB  
  init_bounded_vector log_dev_Nage(2,nages,log_Nage_dev_LO,log_Nage_dev_HI,log_Nage_dev_PH);
  vector log_Nage_dev_output(1,nages);      //used in output. equals zero for first age
  matrix B(styr,endyr+1,1,nages);           //Population biomass by year and age at start of yr
  vector totB(styr,endyr+1);                //Total biomass by year
  vector totN(styr,endyr+1);                //Total abundance by year
  vector SSB(styr,endyr);                   //Total spawning biomass by year (female + male mature biomass) 
  vector rec(styr,endyr+1);                 //Recruits by year
  vector prop_f(1,nages);
  vector prop_m(1,nages);
  vector maturity_f(1,nages);
  vector maturity_m(1,nages);
  vector reprod(1,nages);
 
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
  vector SdS0(styr,endyr);                      //Spawners relative to the unfished level              

  
  init_bounded_number log_dm_lenc_cHL(log_dm_cHL_lc_LO,log_dm_cHL_lc_HI,log_dm_cHL_lc_PH);
  //init_bounded_number log_dm_cOT_lc(log_dm_cOT_lc_LO,log_dm_cOT_lc_HI,log_dm_cOT_lc_PH);
 // init_bounded_number log_dm_rHB_lc(log_dm_rHB_lc_LO,log_dm_rHB_lc_HI,log_dm_rHB_lc_PH);
 // init_bounded_number log_dm_rGN_lc(log_dm_rGN_lc_LO,log_dm_rGN_lc_HI,log_dm_rGN_lc_PH);
  init_bounded_number log_dm_lenc_rHB_D(log_dm_rHB_D_lc_LO,log_dm_rHB_D_lc_HI,log_dm_rHB_D_lc_PH);
  init_bounded_number log_dm_lenc_cHL_D(log_dm_cHL_D_lc_LO,log_dm_cHL_D_lc_HI,log_dm_cHL_D_lc_PH);
  //init_bounded_number log_dm_sFT_lc(log_dm_sFT_lc_LO,log_dm_sFT_lc_HI,log_dm_sFT_lc_PH);
  ////init_bounded_number log_dm_sCT_lc(log_dm_sCT_lc_LO,log_dm_sCT_lc_HI,log_dm_sCT_lc_PH);
  init_bounded_number log_dm_agec_cHL(log_dm_cHL_ac_LO,log_dm_cHL_ac_HI,log_dm_cHL_ac_PH);
  init_bounded_number log_dm_agec_rHB(log_dm_rHB_ac_LO,log_dm_rHB_ac_HI,log_dm_rHB_ac_PH);
  init_bounded_number log_dm_agec_rGN(log_dm_rGN_ac_LO,log_dm_rGN_ac_HI,log_dm_rGN_ac_PH);
  init_bounded_number log_dm_agec_sCT(log_dm_sCT_ac_LO,log_dm_sCT_ac_HI,log_dm_sCT_ac_PH);
  vector log_dm_cHL_lc_out(1,8);
  //vector log_dm_cOT_lc_out(1,8);
  //vector log_dm_rHB_lc_out(1,8);
  //vector log_dm_rGN_lc_out(1,8);
  vector log_dm_rHB_D_lc_out(1,8);
  vector log_dm_cHL_D_lc_out(1,8);
  //vector log_dm_sFT_lc_out(1,8);
  ////vector log_dm_sCT_lc_out(1,8);
  vector log_dm_cHL_ac_out(1,8);
  vector log_dm_rHB_ac_out(1,8);
  vector log_dm_rGN_ac_out(1,8);
  vector log_dm_sCT_ac_out(1,8);

//-----------------------------------------------------------------------------------------------------------------------------------------------
////---Selectivity-------------------------------------------------------------------------

//Florida Snapper Trap index selectivity  
  matrix sel_sFT(styr,endyr,1,nages);
  //vector sel_sFT_vec(1,nages);
  
  //number selpar_A50_sFT;
  //number selpar_slope_sFT;
  //number selpar_A502_sFT;
  //number selpar_slope2_sFT;
    
//SEFIS index selectivity  
  matrix sel_sCT(styr,endyr,1,nages);
  vector sel_sCT_vec(1,nages);
  
  init_bounded_number selpar_A50_sCT(selpar_A50_sCT_LO,selpar_A50_sCT_HI,selpar_A50_sCT_PH);
  init_bounded_number selpar_slope_sCT(selpar_slope_sCT_LO,selpar_slope_sCT_HI,selpar_slope_sCT_PH);
  init_bounded_number selpar_A502_sCT(selpar_A502_sCT_LO,selpar_A502_sCT_HI,selpar_A502_sCT_PH);
  init_bounded_number selpar_slope2_sCT(selpar_slope2_sCT_LO,selpar_slope2_sCT_HI,selpar_slope2_sCT_PH);
  init_bounded_number selpar_logit_Age1_sCT(selpar_Age1_sCT_logit_LO,selpar_Age1_sCT_logit_HI,selpar_Age1_sCT_logit_PH);
  init_bounded_number selpar_logit_Age2_sCT(selpar_Age2_sCT_logit_LO,selpar_Age2_sCT_logit_HI,selpar_Age2_sCT_logit_PH);
  
  number selpar_Age2_sCT; //age 1 value in arithmetic space, after logit backtransform 
  number selpar_Age1_sCT; //age 1 value in arithmetic space, after logit backtransform 
 
  vector selpar_A50_sCT_out(1,8);
  vector selpar_slope_sCT_out(1,8);
  vector selpar_A502_sCT_out(1,8);
  vector selpar_slope2_sCT_out(1,8);
  vector selpar_Age1_sCT_logit_out(1,8);
  vector selpar_Age2_sCT_logit_out(1,8);

  //Commercial handline-------------------------------------------------
  matrix sel_cHL(styr,endyr,1,nages);
  vector sel_cHL_block1(1,nages);
  vector sel_cHL_block2(1,nages);
     
  init_bounded_number selpar_A50_cHL1(selpar_A50_cHL1_LO,selpar_A50_cHL1_HI,selpar_A50_cHL1_PH);
  init_bounded_number selpar_slope_cHL1(selpar_slope_cHL1_LO,selpar_slope_cHL1_HI,selpar_slope_cHL1_PH);
  number selpar_A502_cHL1;
  number selpar_slope2_cHL1;
  
  //init_bounded_number selpar_A502_cHL1(selpar_A502_cHL1_LO,selpar_A502_cHL1_HI,selpar_A502_cHL1_PH);
  //init_bounded_number selpar_slope2_cHL1(selpar_slope2_cHL1_LO,selpar_slope2_cHL1_HI,selpar_slope2_cHL1_PH);
  
  init_bounded_number selpar_A50_cHL2(selpar_A50_cHL2_LO,selpar_A50_cHL2_HI,selpar_A50_cHL2_PH);
  init_bounded_number selpar_slope_cHL2(selpar_slope_cHL2_LO,selpar_slope_cHL2_HI,selpar_slope_cHL2_PH);
  number selpar_A502_cHL2;
  number selpar_slope2_cHL2;
  
  //init_bounded_number selpar_A502_cHL2(selpar_A502_cHL2_LO,selpar_A502_cHL2_HI,selpar_A502_cHL2_PH);
  //init_bounded_number selpar_slope2_cHL2(selpar_slope2_cHL2_LO,selpar_slope2_cHL2_HI,selpar_slope2_cHL2_PH);
  
  vector selpar_A50_cHL1_out(1,8);
  vector selpar_slope_cHL1_out(1,8);
  //vector selpar_A502_cHL1_out(1,8);
  //vector selpar_slope2_cHL1_out(1,8);
  vector selpar_A50_cHL2_out(1,8);
  vector selpar_slope_cHL2_out(1,8);
  //vector selpar_A502_cHL2_out(1,8);
  //vector selpar_slope2_cHL2_out(1,8);
  
 //Commercial handline discards
  matrix sel_cHL_D(styr,endyr,1,nages);  
  //!!cout<<sel_cHL_D<<endl;
  number selpar_Age1_cHL_D2;
  number selpar_Age2_cHL_D2;
  number selpar_A50_cHL_D3;  
  number selpar_slope_cHL_D3;
  number selpar_slope2_cHL_D3;
  vector sel_cHL_D_2(1,nages); //sel in period 2  
  vector sel_cHL_D_3(1,nages); //sel in period 3  
  
  init_bounded_number selpar_logit_Age2_cHL_D2(selpar_Age2_cHL_D2_logit_LO,selpar_Age2_cHL_D2_logit_HI,selpar_Age2_cHL_D2_logit_PH);
  init_bounded_number selpar_A502_cHL_D3(selpar_A502_cHL_D3_LO,selpar_A502_cHL_D3_HI,selpar_A502_cHL_D3_PH);
  
  vector selpar_Age2_cHL_D2_logit_out(1,8);
  vector selpar_A502_cHL_D3_out(1,8);
  
  //Commercial other-------------------------------------------------
  matrix sel_cOT(styr,endyr,1,nages);
  vector sel_cOT_block1(1,nages);
  vector sel_cOT_block2(1,nages);
  //vector sel_cOT_block3(1,nages);

  init_bounded_number selpar_A50_cOT1(selpar_A50_cOT1_LO,selpar_A50_cOT1_HI,selpar_A50_cOT1_PH);
  init_bounded_number selpar_slope_cOT1(selpar_slope_cOT1_LO,selpar_slope_cOT1_HI,selpar_slope_cOT1_PH);
  init_bounded_number selpar_A502_cOT1(selpar_A502_cOT1_LO,selpar_A502_cOT1_HI,selpar_A502_cOT1_PH);
  init_bounded_number selpar_slope2_cOT1(selpar_slope2_cOT1_LO,selpar_slope2_cOT1_HI,selpar_slope2_cOT1_PH);
  
  number selpar_A50_cOT2;
  number selpar_slope_cOT2;
  number selpar_A502_cOT2;
  number selpar_slope2_cOT2;
  
  vector selpar_A50_cOT1_out(1,8);
  vector selpar_slope_cOT1_out(1,8);
  vector selpar_A502_cOT1_out(1,8);
  vector selpar_slope2_cOT1_out(1,8);

  //Historic Trawl 
  matrix sel_cTW(styr,endyr,1,nages);

    //Headboat -------------------------------------------------
  matrix sel_rHB(styr,endyr,1,nages);
  vector sel_rHB_block1(1,nages);
  vector sel_rHB_block2(1,nages);
  vector sel_rHB_block3(1,nages);
  vector sel_rHB_block4(1,nages);
    
  init_bounded_number selpar_A50_rHB1(selpar_A50_rHB1_LO,selpar_A50_rHB1_HI,selpar_A50_rHB1_PH);
  init_bounded_number selpar_slope_rHB1(selpar_slope_rHB1_LO,selpar_slope_rHB1_HI,selpar_slope_rHB1_PH);
  init_bounded_number selpar_A502_rHB1(selpar_A502_rHB1_LO,selpar_A502_rHB1_HI,selpar_A502_rHB1_PH);
  init_bounded_number selpar_slope2_rHB1(selpar_slope2_rHB1_LO,selpar_slope2_rHB1_HI,selpar_slope2_rHB1_PH);
  init_bounded_number selpar_A50_rHB2(selpar_A50_rHB2_LO,selpar_A50_rHB2_HI,selpar_A50_rHB2_PH);
  init_bounded_number selpar_slope_rHB2(selpar_slope_rHB2_LO,selpar_slope_rHB2_HI,selpar_slope_rHB2_PH);    
  init_bounded_number selpar_A502_rHB2(selpar_A502_rHB2_LO,selpar_A502_rHB2_HI,selpar_A502_rHB2_PH);
  init_bounded_number selpar_slope2_rHB2(selpar_slope2_rHB2_LO,selpar_slope2_rHB2_HI,selpar_slope2_rHB2_PH);    
  init_bounded_number selpar_A50_rHB3(selpar_A50_rHB3_LO,selpar_A50_rHB3_HI,selpar_A50_rHB3_PH);
  init_bounded_number selpar_slope_rHB3(selpar_slope_rHB3_LO,selpar_slope_rHB3_HI,selpar_slope_rHB3_PH);
  init_bounded_number selpar_A502_rHB3(selpar_A502_rHB3_LO,selpar_A502_rHB3_HI,selpar_A502_rHB3_PH);
  init_bounded_number selpar_slope2_rHB3(selpar_slope2_rHB3_LO,selpar_slope2_rHB3_HI,selpar_slope2_rHB3_PH);
  init_bounded_number selpar_A50_rHB4(selpar_A50_rHB4_LO,selpar_A50_rHB4_HI,selpar_A50_rHB4_PH);
  init_bounded_number selpar_slope_rHB4(selpar_slope_rHB4_LO,selpar_slope_rHB4_HI,selpar_slope_rHB4_PH);
  init_bounded_number selpar_A502_rHB4(selpar_A502_rHB4_LO,selpar_A502_rHB4_HI,selpar_A502_rHB4_PH);
  init_bounded_number selpar_slope2_rHB4(selpar_slope2_rHB4_LO,selpar_slope2_rHB4_HI,selpar_slope2_rHB4_PH);
  
  vector selpar_A50_rHB1_out(1,8);
  vector selpar_slope_rHB1_out(1,8);
  vector selpar_A502_rHB1_out(1,8);
  vector selpar_slope2_rHB1_out(1,8);
  vector selpar_A50_rHB2_out(1,8);
  vector selpar_slope_rHB2_out(1,8);
  vector selpar_A502_rHB2_out(1,8);
  vector selpar_slope2_rHB2_out(1,8);
  vector selpar_A50_rHB3_out(1,8);
  vector selpar_slope_rHB3_out(1,8);
  vector selpar_A502_rHB3_out(1,8);
  vector selpar_slope2_rHB3_out(1,8);
  vector selpar_A50_rHB4_out(1,8);
  vector selpar_slope_rHB4_out(1,8);
  vector selpar_A502_rHB4_out(1,8);
  vector selpar_slope2_rHB4_out(1,8);

  //Headboat discards
  matrix sel_rHB_D(styr,endyr,1,nages);
  number selpar_Age1_rHB_D3; //period 1,2
  number selpar_Age1_rHB_D4; //period 1,2 
  number selpar_Age2_rHB_D3;
  number selpar_Age2_rHB_D4;
  
  init_bounded_number selpar_logit_Age2_rHB_D3(selpar_Age2_rHB_D3_logit_LO,selpar_Age2_rHB_D3_logit_HI,selpar_Age2_rHB_D3_logit_PH);
  init_bounded_number selpar_logit_Age2_rHB_D4(selpar_Age2_rHB_D4_logit_LO,selpar_Age2_rHB_D4_logit_HI,selpar_Age2_rHB_D4_logit_PH);
  
  //vector sel_rHB_D_1(1,nages); //sel in period 1
  vector sel_rHB_D_2(1,nages); //sel in period 2
  vector sel_rHB_D_3(1,nages); //sel in period 3  
  vector sel_rHB_D_4(1,nages); //sel in period 4 

  vector selpar_Age2_rHB_D3_logit_out(1,8);
  vector selpar_Age2_rHB_D4_logit_out(1,8);
  
  //General Rec 
  matrix sel_rGN(styr,endyr,1,nages);
  vector sel_rGN_block1(1,nages);
  vector sel_rGN_block2(1,nages);
  vector sel_rGN_block3(1,nages);
  vector sel_rGN_block4(1,nages);

  init_bounded_number selpar_A50_rGN2(selpar_A50_rGN2_LO,selpar_A50_rGN2_HI,selpar_A50_rGN2_PH);
  init_bounded_number selpar_slope_rGN2(selpar_slope_rGN2_LO,selpar_slope_rGN2_HI,selpar_slope_rGN2_PH); 
  init_bounded_number selpar_A50_rGN3(selpar_A50_rGN3_LO,selpar_A50_rGN3_HI,selpar_A50_rGN3_PH);
  init_bounded_number selpar_slope_rGN3(selpar_slope_rGN3_LO,selpar_slope_rGN3_HI,selpar_slope_rGN3_PH);
  init_bounded_number selpar_A502_rGN3(selpar_A502_rGN3_LO,selpar_A502_rGN3_HI,selpar_A502_rGN3_PH);
  init_bounded_number selpar_slope2_rGN3(selpar_slope2_rGN3_LO,selpar_slope2_rGN3_HI,selpar_slope2_rGN3_PH);
  init_bounded_number selpar_A50_rGN4(selpar_A50_rGN4_LO,selpar_A50_rGN4_HI,selpar_A50_rGN4_PH);
  init_bounded_number selpar_slope_rGN4(selpar_slope_rGN4_LO,selpar_slope_rGN4_HI,selpar_slope_rGN4_PH);
  init_bounded_number selpar_A502_rGN4(selpar_A502_rGN4_LO,selpar_A502_rGN4_HI,selpar_A502_rGN4_PH);
  init_bounded_number selpar_slope2_rGN4(selpar_slope2_rGN4_LO,selpar_slope2_rGN4_HI,selpar_slope2_rGN4_PH);

  vector selpar_A50_rGN2_out(1,8);
  vector selpar_slope_rGN2_out(1,8);
  vector selpar_A50_rGN3_out(1,8);
  vector selpar_slope_rGN3_out(1,8);
  vector selpar_A502_rGN3_out(1,8);
  vector selpar_slope2_rGN3_out(1,8);  
  vector selpar_A50_rGN4_out(1,8);
  vector selpar_slope_rGN4_out(1,8);  
  vector selpar_A502_rGN4_out(1,8);
  vector selpar_slope2_rGN4_out(1,8);
  
//General rec discards - get the same selex as rHB
  matrix sel_rGN_D(styr,endyr,1,nages);  
  
//Discard selectivities
  //matrix sel_D(styr,endyr,1,nages);
  //vector sel_D_block1(1,nages);
  //vector sel_D_block2(1,nages);
  //vector sel_D_block3(1,nages);
  vector prob_belowsizelim_block2(1,nages);  //kc uncommented
  vector prob_belowsizelim_block3(1,nages);  //kc uncommented
  vector prob_belowsizelim_block4(1,nages);  //kc added
  number zscore_lsizelim1;  //kc uncommented
  number zscore_lsizelim2;  //kc uncommented
  number zscore_lsizelim3;  //kc added
  number cprob_lsizelim1;  //kc uncommented
  number cprob_lsizelim2;  //kc uncommented
  number cprob_lsizelim3;  //kc adeed
	
  //init_bounded_number selpar_age1logit_D(selpar_age1logit_D_LO,selpar_age1logit_D_HI,selpar_age1logit_D_PH);
  //number selpar_age1_D; //age 1 value in arithmetic space, after logit backtransform 
 // vector selpar_age1logit_D_out(1,8);


  
//Weighted total selectivity--------------------------------------------  
  //effort-weighted, recent selectivities
  vector sel_wgted_L(1,nages);  //toward landings 
  vector sel_wgted_D(1,nages);  //toward discards    
  vector sel_wgted_tot(1,nages);//toward Z, landings plus deads discards

//-----------------------------------------------------------------------------------------------------------------------------------------------
//-------CPUE Predictions--------------------------------
  vector pred_cHL_cpue(styr_cpue_cHL,endyr_cpue_cHL);                   //predicted cHL index (weight fish per effort)
  matrix N_cHL(styr_cpue_cHL,endyr_cpue_cHL,1,nages);                   //used to compute cHL index
  vector pred_rHB_cpue(styr_cpue_rHB,endyr_cpue_rHB);                   //predicted rHB index (number fish per effort)
  matrix N_rHB(styr_cpue_rHB,endyr_cpue_rHB,1,nages);                   //used to compute rHB index
  vector pred_rGN_cpue(styr_cpue_rGN,endyr_cpue_rGN);                   //predicted rGN index (number fish per effort)
  matrix N_rGN(styr_cpue_rGN,endyr_cpue_rGN,1,nages);                   //used to compute rGN index
  vector pred_sFT_cpue(styr_cpue_sFT,endyr_cpue_sFT);                //predicted Florida snapper trap index (number fish per effort)
  matrix N_sFT(styr_cpue_sFT,endyr_cpue_sFT,1,nages);                //used to compute SERFS index
  vector pred_sCT_cpue(styr_cpue_sCT,endyr_cpue_sCT);                //predicted SERFS index (number fish per effort)
  matrix N_sCT(styr_cpue_sCT,endyr_cpue_sCT,1,nages);                //used to compute SERFS index
  //matrix sel_rGN_cpue(styr,endyr,1,nages);                   //includes discards 
  
//---Catchability (CPUE q's)----------------------------------------------------------
  init_bounded_number log_q_cpue_cHL(log_q_cHL_LO,log_q_cHL_HI,log_q_cHL_PH);
  init_bounded_number log_q_cpue_rHB(log_q_rHB_LO,log_q_rHB_HI,log_q_rHB_PH);
  init_bounded_number log_q_cpue_rGN(log_q_rGN_LO,log_q_rGN_HI,log_q_rGN_PH);
  init_bounded_number log_q_cpue_sFT(log_q_sFT_LO,log_q_sFT_HI,log_q_sFT_PH);
  init_bounded_number log_q_cpue_sCT(log_q_sCT_LO,log_q_sCT_HI,log_q_sCT_PH);
  
  vector log_q_cHL_out(1,8);
  vector log_q_rHB_out(1,8);
  vector log_q_rGN_out(1,8);
  vector log_q_sFT_out(1,8);
  vector log_q_sCT_out(1,8);  
  
  number q_rate;
  vector q_rate_fcn_cHL(styr_cpue_cHL,endyr_cpue_cHL);         //increase due to technology creep (saturates in 2003) 
  vector q_rate_fcn_rHB(styr_cpue_rHB,endyr_cpue_rHB);         //increase due to technology creep (saturates in 2003) 
  vector q_rate_fcn_rGN(styr_cpue_rGN,endyr_cpue_rGN);       //increase due to technology creep (saturates in 2003) 
  
//init_bounded_number q_DD_beta(0.1,0.9,set_q_DD_phase);    //not estimated so commented out and declared as number (below)
  number q_DD_beta;
  vector q_DD_fcn(styr,endyr);    //density dependent function as a multiple of q (scaled a la Katsukawa and Matsuda. 2003)
  number B0_q_DD;                 //B0 of ages q_DD_age plus
  vector B_q_DD(styr,endyr);      //annual biomass of ages q_DD_age plus

//Fishery dependent random walk catchability
//kc vector q_RW_log_dev_cHL(styr_cpue_cHL,endyr_cpue_cHL-1); 
//kc vector q_RW_log_dev_rHB(styr_cpue_rHB,endyr_cpue_rHB-1);
//kc vector q_RW_log_dev_rGN(styr_cpue_rGN,endyr_cpue_rGN-1);
 init_bounded_vector q_RW_log_dev_cHL(styr_cpue_cHL,endyr_cpue_cHL-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
 init_bounded_vector q_RW_log_dev_rHB(styr_cpue_rHB,endyr_cpue_rHB-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
 init_bounded_vector q_RW_log_dev_rGN(styr_cpue_rGN,endyr_cpue_rGN-1,log_RWq_LO,log_RWq_HI,log_RWq_PH);
 vector q_RW_log_dev_sFT(styr_cpue_sFT,endyr_cpue_sFT-1); 
 vector q_RW_log_dev_sCT(styr_cpue_sCT,endyr_cpue_sCT-1); 
 
//Fishery dependent catchability over time, may be constant
 vector q_cHL(styr_cpue_cHL,endyr_cpue_cHL); 
 vector q_rHB(styr_cpue_rHB,endyr_cpue_rHB); 
 vector q_rGN(styr_cpue_rGN,endyr_cpue_rGN);
 vector q_sFT(styr_cpue_sFT,endyr_cpue_sFT); 
 vector q_sCT(styr_cpue_sCT,endyr_cpue_sCT); 
 
//----------------------------------------------------------------------------------------------------------------------------------------------- 
//---Landings in numbers (total or 1000 fish) and in wgt (whole klb)--------------------------------------------------
  matrix L_cHL_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_cHL_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age    
  vector pred_cHL_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages  
  vector pred_cHL_L_klb(styr,endyr);      //yearly landings in 1000 lb whole summed over ages

  matrix L_cOT_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_cOT_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age    
  vector pred_cOT_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages  
  vector pred_cOT_L_klb(styr,endyr);      //yearly landings in 1000 lb whole summed over ages

  matrix L_cTW_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_cTW_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age    
  vector pred_cTW_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages  
  vector pred_cTW_L_klb(styr,endyr);      //yearly landings in 1000 lb whole summed over ages
  
  matrix L_rHB_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_rHB_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age
  vector pred_rHB_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages
  vector pred_rHB_L_klb(styr,endyr);      //yearly landings in 1000 lb whole summed over ages

  matrix L_rGN_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_rGN_klb(styr,endyr,1,nages);   //landings (1000 lb whole weight) at age
  vector pred_rGN_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages
  vector pred_rGN_L_klb(styr,endyr);      //yearly landings in 1000 lb whole summed over ages

  matrix L_total_num(styr,endyr,1,nages);//total landings in number at age
  matrix L_total_klb(styr,endyr,1,nages);//landings in klb whole wgt at age 
  vector L_total_knum_yr(styr,endyr);    //total landings in 1000 fish by yr summed over ages  
  vector L_total_klb_yr(styr,endyr);     //total landings (klb whole wgt) by yr summed over ages
  
//---Discards (number dead fish) --------------------------------------------------
  matrix D_cHL_num(styr,endyr,1,nages);            //discards (numbers) at age
  matrix D_cHL_klb(styr,endyr,1,nages);             //discards (1000 lb whole) at age
  //vector pred_cHL_D_knum(styr_D_cHL,endyr_D_cHL);     //yearly dead discards summed over ages
  vector obs_cHL_D(styr_D_cHL,endyr_D_cHL);           //observed releases multiplied by discard mortality
  //vector pred_cHL_D_klb(styr_D_cHL,endyr_D_cHL);      //yearly dead discards in klb whole wgt summed over ages
  vector pred_cHL_D_knum(styr,endyr);     //yearly dead discards summed over ages
  vector pred_cHL_D_klb(styr,endyr);      //yearly dead discards in klb whole wgt summed over ages
  
  matrix D_rHB_num(styr,endyr,1,nages);             //discards (numbers) at age
  matrix D_rHB_klb(styr,endyr,1,nages);             //discards (1000 lb) at age
  //vector pred_rHB_D_knum(styr_D_rHB,endyr_D_rHB);     //yearly dead discards summed over ages
  vector obs_rHB_D(styr_D_rHB,endyr_D_rHB);           //observed releases multiplied by discard mortality
  //vector pred_rHB_D_klb(styr_D_rHB,endyr_D_rHB);      //yearly dead discards in klb whole wgt summed over ages
  vector pred_rHB_D_knum(styr,endyr);     //yearly dead discards summed over ages
  vector pred_rHB_D_klb(styr,endyr);      //yearly dead discards in klb whole wgt summed over ages
  
  matrix D_rGN_num(styr,endyr,1,nages);             //discards (numbers) at age
  matrix D_rGN_klb(styr,endyr,1,nages);             //discards (1000 lb) at age
  //vector pred_rGN_D_knum(styr_D_rGN,endyr_D_rGN);     //yearly dead discards summed over ages
  vector obs_rGN_D(styr_D_rGN,endyr_D_rGN);           //observed releases multiplied by discard mortality
  //vector pred_rGN_D_klb(styr_D_rGN,endyr_D_rGN);      //yearly dead discards in klb whole wgt summed over ages
  vector pred_rGN_D_knum(styr,endyr);     //yearly dead discards summed over ages
  vector pred_rGN_D_klb(styr,endyr);      //yearly dead discards in klb whole wgt summed over ages
  
  matrix D_total_num(styr,endyr,1,nages);          //total discards in number at age
  matrix D_total_klb(styr,endyr,1,nages);          //discards in klb wgt at age 
  vector D_total_knum_yr(styr,endyr);              //total discards in 1000 fish by yr summed over ages  
  vector D_total_klb_yr(styr,endyr);               //total discards (klb whole wgt) by yr summed over ages

  number Dmort_cHL;
  number Dmort_rHB;
  number Dmort_rGN;  

////---MSY calcs----------------------------------------------------------------------------
  number F_cHL_prop;       //proportion of F_sum attributable to cHL, last X=selpar_n_yrs_wgted yrs
  number F_cOT_prop;       //proportion of F_sum attributable to cHL, last X=selpar_n_yrs_wgted yrs  
  number F_cTW_prop;      //proportion of F_sum attributable to HTR, last X=selpar_n_yrs_wgted yrs 
  number F_rHB_prop;       //proportion of F_sum attributable to rHB, last X=selpar_n_yrs_wgted yrs 
  number F_rGN_prop;       //proportion of F_sum attributable to rGN, last X=selpar_n_yrs_wgted yrs
  number F_cHL_D_prop;     //proportion of F_sum attributable to cHL discards, last X=selpar_n_yrs_wgted yrs
  number F_rHB_D_prop;     //proportion of F_sum attributable to rHB discards, last X=selpar_n_yrs_wgted yrs
  number F_rGN_D_prop;     //proportion of F_sum attributable to rGN discards, last X=selpar_n_yrs_wgted yrs
    
  number F_temp_sum;      //sum of geom mean Fsum's in last X yrs, used to compute F_fishery_prop

  vector F_end(1,nages);
  vector F_end_L(1,nages);  
  vector F_end_D(1,nages);    
  number F_end_apex;
  
  number SSB_msy_out;           //SSB (total mature biomass) at msy
  number F_msy_out;             //F at msy
  number msy_klb_out;           //max sustainable yield (1000 lb whole wgt)
  number msy_knum_out;          //max sustainable yield (1000 fish)  
  number D_msy_klb_out;         //discards associated with msy (1000 lb whole wgt)
  number D_msy_knum_out;        //discards associated with msy (1000 fish)  
  number B_msy_out;             //total biomass at MSY 
  number R_msy_out;             //equilibrium recruitment at F=Fmsy
  number spr_msy_out;           //spr at F=Fmsy

  number F20_dum;				//intermediate calculation for F20
  number F30_dum;				//intermediate calculation for F30
  number F40_dum;				//intermediate calculation for F40
  number F20_out;              	//F20
  number F30_out;              	//F30
  number F40_out;              	//F40
  number SSB_F30_out;  
  number B_F30_out;
  number R_F30_out;
  number L_F30_knum_out;
  number L_F30_klb_out;
  number D_F30_knum_out;
  number D_F30_klb_out;    
  number rec_mean;  			//arithemetic average recruitment used in SPR-related quantities

  vector N_age_msy(1,nages);         //numbers at age for MSY calculations: beginning of yr
  vector N_age_msy_spawn(1,nages);   //numbers at age for MSY calculations: time of peak spawning  
  vector L_age_msy(1,nages);         //landings at age for MSY calculations
  vector D_age_msy(1,nages);         //discard mortality (dead discards) at age for MSY calculations
  vector Z_age_msy(1,nages);         //total mortality at age for MSY calculations
  vector F_L_age_msy(1,nages);       //fishing mortality landings (not discards) at age for MSY calculations
  vector F_D_age_msy(1,nages);       //fishing mortality of discards at age for MSY calculations
  vector F_msy(1,n_iter_msy);        //values of full F to be used in equilibrium calculations
  vector spr_msy(1,n_iter_msy);      //reproductive capacity-per-recruit values corresponding to F values in F_msy
  vector R_eq(1,n_iter_msy);         //equilibrium recruitment values corresponding to F values in F_msy
  vector L_eq_klb(1,n_iter_msy);     //equilibrium landings(klb whole wgt) values corresponding to F values in F_msy
  vector L_eq_knum(1,n_iter_msy);    //equilibrium landings(1000 fish) values corresponding to F values in F_msy
  vector D_eq_klb(1,n_iter_msy);     //equilibrium discards(klb whole wgt) values corresponding to F values in F_msy
  vector D_eq_knum(1,n_iter_msy);    //equilibrium discards(1000 fish) values corresponding to F values in F_msy
  vector SSB_eq(1,n_iter_msy);       //equilibrium reproductive capacity values corresponding to F values in F_msy
  vector B_eq(1,n_iter_msy);         //equilibrium biomass values corresponding to F values in F_msy
  
  vector FdF_msy(styr,endyr);
  vector FdF30(styr,endyr);
  vector SdSSB_msy(styr,endyr);	 
  number SdSSB_msy_end;
  number FdF_msy_end;
  number FdF_msy_end_mean;           //geometric mean of last X yrs  
  vector SdSSB_F30(styr,endyr);	 
  vector Sdmsst_F30(styr,endyr);	 
  number SdSSB_F30_end;
  number Sdmsst_F30_end;
  number FdF30_end_mean;             //geometric mean of last selpar_n_yrs_wgted yrs  
  number Fend_mean_temp;			 //intermediate calc for geometric mean of last selpar_n_yrs_wgted yrs
  number Fend_mean;					 //geometric mean of last selpar_n_yrs_wgted yrs
  vector L_age_F30(1,nages);         //landings at age for F30 calculations
  vector D_age_F30(1,nages);         //discard mortality (dead discards) at age for F30 calculations
  
  vector wgt_wgted_L_klb(1,nages);   //fishery-weighted average weight at age of landings in whole weight
  vector wgt_wgted_D_klb(1,nages);   //fishery-weighted average weight at age of discards in whole weight  
  number wgt_wgted_L_denom;          //used in intermediate calculations
  number wgt_wgted_D_denom;          //used in intermediate calculations

  number iter_inc_msy;               //increments used to compute msy, equals 1/(n_iter_msy-1)
  
////--------Mortality------------------------------------------------------------------

  vector M(1,nages);                         //age-dependent natural mortality
  init_bounded_number M_constant(M_constant_LO,M_constant_HI,M_constant_PH);   //age-indpendent: used only for MSST
  vector M_constant_out(1,8);
  number smsy2msstM;                         //scales Smsy to get msst using (1-M). Used only in output.
  number smsy2msst75;                        //scales Smsy to get msst using 75%. Used only in output.  
  
  matrix F(styr,endyr,1,nages);
  vector Fsum(styr,endyr);                   //Full fishing mortality rate by year
  vector Fapex(styr,endyr);                  //Max across ages, fishing mortality rate by year (may differ from Fsum bc of dome-shaped sel 
  matrix Z(styr,endyr,1,nages);

  init_bounded_number log_avg_F_L_cHL(log_avg_F_cHL_LO,log_avg_F_cHL_HI,log_avg_F_cHL_PH);
  vector log_avg_F_cHL_out(1,8);  
  init_bounded_dev_vector log_dev_F_L_cHL(styr_L_cHL,endyr_L_cHL,log_F_dev_cHL_LO,log_F_dev_cHL_HI,log_F_dev_cHL_PH);
  vector log_F_dev_cHL_out(styr_L_cHL,endyr_L_cHL);
  matrix F_cHL(styr,endyr,1,nages);
  vector F_cHL_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_cHL;
  number log_F_dev_end_cHL; 

  init_bounded_number log_avg_F_L_cOT(log_avg_F_cOT_LO,log_avg_F_cOT_HI,log_avg_F_cOT_PH);
  vector log_avg_F_cOT_out(1,8);  
  init_bounded_dev_vector log_dev_F_L_cOT(styr_L_cOT,endyr_L_cOT,log_F_dev_cOT_LO,log_F_dev_cOT_HI,log_F_dev_cOT_PH);
  vector log_F_dev_cOT_out(styr_L_cOT,endyr_L_cOT);
  matrix F_cOT(styr,endyr,1,nages);
  vector F_cOT_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_cOT;
  number log_F_dev_end_cOT; 

  init_bounded_number log_avg_F_L_cTW(log_avg_F_cTW_LO,log_avg_F_cTW_HI,log_avg_F_cTW_PH);
  vector log_avg_F_cTW_out(1,8);  
  init_bounded_dev_vector log_dev_F_L_cTW(styr_L_cTW,endyr_L_cTW,log_F_dev_cTW_LO,log_F_dev_cTW_HI,log_F_dev_cTW_PH);
  vector log_F_dev_cTW_out(styr_L_cTW,endyr_L_cTW);
  matrix F_cTW(styr,endyr,1,nages);
  vector F_cTW_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_cTW;
  number log_F_dev_end_cTW;
  
  init_bounded_number log_avg_F_L_rHB(log_avg_F_rHB_LO,log_avg_F_rHB_HI,log_avg_F_rHB_PH);
  vector log_avg_F_rHB_out(1,8); 
  init_bounded_dev_vector log_dev_F_L_rHB(styr_L_rHB,endyr_L_rHB,log_F_dev_rHB_LO,log_F_dev_rHB_HI,log_F_dev_rHB_PH);    
  vector log_F_dev_rHB_out(styr_L_rHB,endyr_L_rHB);
  matrix F_rHB(styr,endyr,1,nages);
  vector F_rHB_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_rHB;    
  number log_F_dev_end_rHB;  

  init_bounded_number log_avg_F_L_rGN(log_avg_F_rGN_LO,log_avg_F_rGN_HI,log_avg_F_rGN_PH);
  vector log_avg_F_rGN_out(1,8); 
  init_bounded_dev_vector log_dev_F_L_rGN(styr_L_rGN,endyr_L_rGN,log_F_dev_rGN_LO,log_F_dev_rGN_HI,log_F_dev_rGN_PH);    
  vector log_F_dev_rGN_out(styr_L_rGN,endyr_L_rGN);
  matrix F_rGN(styr,endyr,1,nages);
  vector F_rGN_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_rGN;    
  number log_F_dev_end_rGN;  

  init_bounded_number log_avg_F_D_cHL(log_avg_F_cHL_D_LO,log_avg_F_cHL_D_HI,log_avg_F_cHL_D_PH);
  vector log_avg_F_cHL_D_out(1,8);  
  init_bounded_dev_vector log_dev_F_D_cHL(styr_D_cHL,endyr_D_cHL,log_F_dev_cHL_D_LO,log_F_dev_cHL_D_HI,log_F_dev_cHL_D_PH);
  vector log_F_dev_cHL_D_out(styr_D_cHL,endyr_D_cHL);
  matrix F_cHL_D(styr,endyr,1,nages);
  vector F_cHL_D_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_dev_end_cHL_D; 
  number log_F_avgdev_cHL_D;

  init_bounded_number log_avg_F_D_rHB(log_avg_F_rHB_D_LO,log_avg_F_rHB_D_HI,log_avg_F_rHB_D_PH);
  vector log_avg_F_rHB_D_out(1,8);  
  init_bounded_dev_vector log_dev_F_D_rHB(styr_D_rHB,endyr_D_rHB,log_F_dev_rHB_D_LO,log_F_dev_rHB_D_HI,log_F_dev_rHB_D_PH);
  vector log_F_dev_rHB_D_out(styr_D_rHB,endyr_D_rHB);
  matrix F_rHB_D(styr,endyr,1,nages);
  vector F_rHB_D_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_dev_end_rHB_D; 
  number log_F_avgdev_rHB_D;
  
  init_bounded_number log_avg_F_D_rGN(log_avg_F_rGN_D_LO,log_avg_F_rGN_D_HI,log_avg_F_rGN_D_PH);
  vector log_avg_F_rGN_D_out(1,8);  
  init_bounded_dev_vector log_dev_F_D_rGN(styr_D_rGN,endyr_D_rGN,log_F_dev_rGN_D_LO,log_F_dev_rGN_D_HI,log_F_dev_rGN_D_PH);
  vector log_F_dev_rGN_D_out(styr_D_rGN,endyr_D_rGN);
  matrix F_rGN_D(styr,endyr,1,nages);
  vector F_rGN_D_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_dev_end_rGN_D; 
  number log_F_avgdev_rGN_D;

  //number F_init_cHL_prop;  //proportion of F_init attributable to cHL, first X yrs, No discards in initial yrs
  //number F_init_cOT_prop;  //proportion of F_init attributable to cOT, first X yrs
  //number F_init_rHB_prop;  //proportion of F_init attributable to rHB, first X yrs
  //number F_init_rGN_prop;  //proportion of F_init attributable to rGN, first X yrs
  //number F_init_denom;    //intermediate calcs

//---Per-recruit stuff----------------------------------------------------------------------------------
  vector N_age_spr(1,nages);         //numbers at age for SPR calculations: beginning of year
  vector N_age_spr_spawn(1,nages);   //numbers at age for SPR calculations: time of peak spawning  
  vector L_age_spr(1,nages);         //catch at age for SPR calculations
  vector Z_age_spr(1,nages);         //total mortality at age for SPR calculations
  vector spr_static(styr,endyr);     //vector of static SPR values by year
  vector F_L_age_spr(1,nages);       //fishing mortality of landings (not discards) at age for SPR calculations
  vector F_D_age_spr(1,nages);       //fishing mortality of discards at age for SPR calculations
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
 
 number sdnr_lc_cHL;
 number sdnr_lc_cOT;
// number sdnr_lc_rHB;
 number sdnr_lc_rHB_D; 
 number sdnr_lc_cHL_D;  
// number sdnr_lc_rGN;
 number sdnr_lc_sFT; 
//number sdnr_lc_sCT; 

 number sdnr_ac_cHL;
 number sdnr_ac_rHB;
 number sdnr_ac_rGN;
 number sdnr_ac_sCT;  
  
  number sdnr_I_cHL;
  number sdnr_I_rHB;
  number sdnr_I_rGN;
  number sdnr_I_sFT;
  number sdnr_I_sCT;
    
////-------Objective function components-----------------------------------------------------------------------------
  number w_L;
  number w_D;
  
  number w_cpue_cHL;
  number w_cpue_rHB;
  number w_cpue_rGN;
  number w_cpue_sFT;
  number w_cpue_sCT;
   
 number w_lenc_cHL;
 number w_lenc_cOT; 
// number w_lc_rHB;
 number w_lenc_rHB_D; 
 number w_lenc_cHL_D;  
// number w_lc_rGN;
 number w_lenc_sFT;
 //number w_lc_sCT;  
 
 number w_agec_cHL; 
 number w_agec_rHB;
 number w_agec_rGN;  
 number w_agec_sCT; 
  
  number w_Nage_init;  
  number w_rec;
  number w_rec_early;
  number w_rec_end;
  number w_fullF;  
  number w_Ftune;

  number f_cHL_L; 
  number f_cOT_L;
  number f_cTW_L; 
  number f_rHB_L; 
  number f_rGN_L; 

  number f_cHL_D; 
  number f_rHB_D; 
  number f_rGN_D; 

  number f_cHL_cpue;
  number f_rHB_cpue;
  number f_rGN_cpue;
  number f_sFT_cpue; 
  number f_sCT_cpue;  
   
  number f_rHB_RWq_cpue;
  number f_cHL_RWq_cpue;
   number f_rGN_RWq_cpue;   //kc
   
  number f_cHL_lenc;
  number f_cOT_lenc;  
 // number f_rHB_lenc;
  number f_rHB_D_lenc;
  number f_cHL_D_lenc;  
//  number f_rGN_lenc;
  number f_sFT_lenc; 
  //number f_sCT_lenc;  
  
  number f_cHL_agec;
  number f_rHB_agec;   
  number f_rGN_agec; 
  number f_sCT_agec; 
  
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
  number F_reg_proj;						   //value used to define the projections				
  vector F_proj(styr_proj,endyr_proj);         //F by yr for projections (=F_reg_proj after regulations start, current F till then)  
  vector L_knum_proj(styr_proj,endyr_proj);    //total landings in 1000 fish for projections  
  vector L_klb_proj(styr_proj,endyr_proj);     //total landings in weight (1000 lb) for projections
  vector D_knum_proj(styr_proj,endyr_proj);    //total discards in 1000 fish for projections  
  vector D_klb_proj(styr_proj,endyr_proj);     //total landings in weight (1000 lb) for projections

  vector B_proj(styr_proj,endyr_proj);         //Biomass for projections  
  vector SSB_proj(styr_proj,endyr_proj);       //SSB for projections  
  vector R_proj(styr_proj,endyr_proj);     	   //recruits for projections
  vector FL_age_proj(1,nages);      		   //F (landings) by age for projections    
  vector FD_age_proj(1,nages);      		   //F (discards) by age for projections    
  
  matrix N_proj(styr_proj,endyr_proj,1,nages);           //Population numbers by year and age at start of yr
  matrix N_spawn_proj(styr_proj,endyr_proj,1,nages);     //Population numbers by year and age at peaking spawning: used for SSB in projections 
  matrix Z_proj(styr_proj,endyr_proj,1,nages);           //Z by year and age for projections 
  matrix L_age_proj(styr_proj,endyr_proj,1,nages);       //Projected landings at age in numbers 
  matrix D_age_proj(styr_proj,endyr_proj,1,nages);       //Projected discards at age in numbers    
  
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
 convergence_criteria 1e-2, 1e-2, 1e-3, 1e-3, 1e-3, 1e-4, 1e-4;
 
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PRELIMINARY_CALCS_SECTION

// Set values of fixed parameters or set initial guess of estimated parameters
  Dmort_cHL=set_Dmort_cHL; 
  Dmort_rHB=set_Dmort_rHB; 
  Dmort_rGN=set_Dmort_rGN;
  
  for(iyear=styr_D_cHL; iyear<=endyr_D_cHL; iyear++)
	{obs_cHL_D(iyear)=Dmort_cHL*obs_released_cHL(iyear);
	}
	
  for(iyear=styr_D_rHB; iyear<=endyr_D_rHB; iyear++)
	{obs_rHB_D(iyear)=Dmort_rHB*obs_released_rHB(iyear);
	}
	
  for(iyear=styr_D_rGN; iyear<=endyr_D_rGN; iyear++)
	{obs_rGN_D(iyear)=Dmort_rGN*obs_released_rGN(iyear);
	}
 
  //age_limit_10in=t0-log(1.0-sizelim1/Linf)/K; //age at size limit: 10" limit;
  //age_limit_11in=t0-log(1.0-sizelim2/Linf)/K; //age at size limit: 11" limit;
  //age_limit_12in=t0-log(1.0-sizelim3/Linf)/K; //age at size limit: 12" limit;
  
  
 //Population		
  Linf=set_Linf(1);
  K=set_K(1);
  t0=set_t0(1);
  len_cv_val=set_len_cv(1);
  
  
  M=set_M; 
  M_constant=set_M_constant(1);
  smsy2msstM=1.0-M_constant;
  smsy2msst75=0.75;  
  
  log_R0=set_log_R0(1);
  steep=set_steep(1);
  R_autocorr=set_R_autocorr(1);
  rec_sigma=set_rec_sigma(1);
  
  log_dm_lenc_cHL=set_log_dm_lenc_cHL(1);
  //log_dm_cOT_lc=set_log_dm_cOT_lc(1);
  //log_dm_rHB_lc=set_log_dm_rHB_lc(1);
 // log_dm_rGN_lc=set_log_dm_rGN_lc(1);
  log_dm_lenc_rHB_D=set_log_dm_lenc_rHB_D(1);
  log_dm_lenc_cHL_D=set_log_dm_lenc_cHL_D(1);
  //log_dm_sFT_lc=set_log_dm_sFT_lc(1);
  ////log_dm_sCT_lc=set_log_dm_sCT_lc(1);
  log_dm_agec_cHL=set_log_dm_agec_cHL(1);
  log_dm_agec_rHB=set_log_dm_agec_rHB(1);
  log_dm_agec_rGN=set_log_dm_agec_rGN(1);
  log_dm_agec_sCT=set_log_dm_agec_sCT(1);
  
  log_q_cpue_cHL=set_log_q_cpue_cHL(1);
  log_q_cpue_rHB=set_log_q_cpue_rHB(1);
  log_q_cpue_rGN=set_log_q_cpue_rGN(1);
  log_q_cpue_sFT=set_log_q_cpue_sFT(1);
  log_q_cpue_sCT=set_log_q_cpue_sCT(1);
  
  q_rate=set_q_rate;
  q_rate_fcn_cHL=1.0;
  q_rate_fcn_rHB=1.0;   
  q_rate_fcn_rGN=1.0; 
  q_DD_beta=set_q_DD_beta;
  q_DD_fcn=1.0;

  q_RW_log_dev_cHL.initialize(); 
  q_RW_log_dev_rHB.initialize(); 
  q_RW_log_dev_rGN.initialize();
  q_RW_log_dev_sFT.initialize();
  q_RW_log_dev_sCT.initialize();
  
   if (set_q_rate_phase<0 & q_rate!=0.0)
  {
    for (iyear=styr_cpue_cHL; iyear<=endyr_cpue_cHL; iyear++)
      {   if (iyear>styr_cpue_cHL & iyear <=2008)  //changed from 2003 KC 
          {//q_rate_fcn_cHL(iyear)=(1.0+q_rate)*q_rate_fcn_cHL(iyear-1); //compound
             q_rate_fcn_cHL(iyear)=(1.0+(iyear-styr_cpue_cHL)*q_rate)*q_rate_fcn_cHL(styr_cpue_cHL);  //linear
          }
          if (iyear>2008) {q_rate_fcn_cHL(iyear)=q_rate_fcn_cHL(iyear-1);} 
      }   
    for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++)
      {   if (iyear>styr_cpue_rHB & iyear <=2008)  //changed from 2003 KC 
          {//q_rate_fcn_rHB(iyear)=(1.0+q_rate)*q_rate_fcn_rHB(iyear-1); //compound
             q_rate_fcn_rHB(iyear)=(1.0+(iyear-styr_cpue_rHB)*q_rate)*q_rate_fcn_rHB(styr_cpue_rHB);  //linear
          }
          if (iyear>2008) {q_rate_fcn_rHB(iyear)=q_rate_fcn_rHB(iyear-1);} 
      }   
     for (iyear=styr_cpue_rGN; iyear<=endyr_cpue_rGN; iyear++)
       {   if (iyear>styr_cpue_rGN & iyear <=2008)  //changed from 2003 KC 
           {//q_rate_fcn_rGN(iyear)=(1.0+q_rate)*q_rate_fcn_rGN(iyear-1); //compound
              q_rate_fcn_rGN(iyear)=(1.0+(iyear-styr_cpue_rGN)*q_rate)*q_rate_fcn_rGN(styr_cpue_rGN);  //linear
           }
           if (iyear>2008) {q_rate_fcn_rGN(iyear)=q_rate_fcn_rGN(iyear-1);} 
       }   
	  
  } //end q_rate conditional      

  w_L=set_w_L;
  w_D=set_w_D;
  
  w_cpue_cHL=set_w_cpue_cHL;
  w_cpue_rHB=set_w_cpue_rHB;
  w_cpue_rGN=set_w_cpue_rGN;
  w_cpue_sFT=set_w_cpue_sFT;
  w_cpue_sCT=set_w_cpue_sCT;
  
 //w_lenc_cHL=set_w_lenc_cHL;
 w_lenc_cOT=set_w_lenc_cOT;  
 //w_lc_rHB=set_w_lc_rHB;
 //w_lenc_rHB_D=set_w_lenc_rHB_D; 
 //w_lenc_cHL_D=set_w_lenc_cHL_D;
 //w_lc_rGN=set_w_lc_rGN;
 w_lenc_sFT=set_w_lenc_sFT; 
 //w_lc_sCT=set_w_lc_sCT;      
 //w_agec_cHL=set_w_agec_cHL;
 //w_agec_rHB=set_w_agec_rHB;
 //w_agec_rGN=set_w_agec_rGN; 
 //w_agec_sCT=set_w_agec_sCT;  
 
  
  w_Nage_init=set_w_Nage_init;
  w_rec=set_w_rec;
  w_rec_early=set_w_rec_early;
  w_rec_end=set_w_rec_end;
  w_fullF=set_w_fullF;
  w_Ftune=set_w_Ftune;

  log_avg_F_L_cHL=set_log_avg_F_L_cHL(1);
  log_avg_F_L_cOT=set_log_avg_F_L_cOT(1);
  log_avg_F_L_cTW=set_log_avg_F_L_cTW(1);
  log_avg_F_L_rHB=set_log_avg_F_L_rHB(1); 
  log_avg_F_L_rGN=set_log_avg_F_L_rGN(1); 
  log_avg_F_D_cHL=set_log_avg_F_D_cHL(1);
  log_avg_F_D_rHB=set_log_avg_F_D_rHB(1); 
  log_avg_F_D_rGN=set_log_avg_F_D_rGN(1); 
    
  log_dev_F_L_cHL=set_log_dev_vals_F_L_cHL;
  log_dev_F_L_cOT=set_log_dev_vals_F_L_cOT;
  log_dev_F_L_cTW=set_log_dev_vals_F_L_cTW;
  log_dev_F_L_rHB=set_log_dev_vals_F_L_rHB;
  log_dev_F_L_rGN=set_log_dev_vals_F_L_rGN;
  log_dev_F_D_cHL=set_log_dev_vals_F_D_cHL;
  log_dev_F_D_rHB=set_log_dev_vals_F_D_rHB;
  log_dev_F_D_rGN=set_log_dev_vals_F_D_rGN;
 
  selpar_A50_cHL1=set_selpar_A50_cHL1(1);
  selpar_slope_cHL1=set_selpar_slope_cHL1(1);
  selpar_A502_cHL1=set_selpar_A502_cHL1(1);
  selpar_slope2_cHL1=set_selpar_slope2_cHL1(1);
  selpar_A50_cHL2=set_selpar_A50_cHL2(1);
  selpar_slope_cHL2=set_selpar_slope_cHL2(1);
  selpar_A502_cHL2=set_selpar_A50_cHL2(1);
  selpar_slope2_cHL2=set_selpar_slope2_cHL2(1);
  
  selpar_A50_cOT1=set_selpar_A50_cOT1(1);
  selpar_slope_cOT1=set_selpar_slope_cOT1(1);
  selpar_A502_cOT1=set_selpar_A502_cOT1(1);
  selpar_slope2_cOT1=set_selpar_slope2_cOT1(1);
  //sel_cOT_block1.initialize();
  sel_cOT_block1(1)=1.0;
  sel_cOT_block1(2)=0.5;
  
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
  selpar_A50_rHB4=set_selpar_A50_rHB4(1);
  selpar_slope_rHB4=set_selpar_slope_rHB4(1);
  selpar_A502_rHB4=set_selpar_A502_rHB4(1);
  selpar_slope2_rHB4=set_selpar_slope2_rHB4(1);

  selpar_A50_rGN2=set_selpar_A50_rGN2(1);
  selpar_slope_rGN2=set_selpar_slope_rGN2(1);
  selpar_A50_rGN3=set_selpar_A50_rGN3(1);
  selpar_slope_rGN3=set_selpar_slope_rGN3(1);
  selpar_A502_rGN3=set_selpar_A502_rGN3(1);
  selpar_slope2_rGN3=set_selpar_slope2_rGN3(1);
  selpar_A50_rGN4=set_selpar_A50_rGN4(1);
  selpar_slope_rGN4=set_selpar_slope_rGN4(1);
  selpar_A502_rGN4=set_selpar_A502_rGN4(1);
  selpar_slope2_rGN4=set_selpar_slope2_rGN4(1);
  
  selpar_A50_sCT=set_selpar_A50_sCT(1);
  selpar_slope_sCT=set_selpar_slope_sCT(1);
  selpar_A502_sCT=set_selpar_A502_sCT(1);
  selpar_slope2_sCT=set_selpar_slope2_sCT(1);
  selpar_logit_Age1_sCT=set_selpar_logit_Age1_sCT(1);  
  selpar_logit_Age2_sCT=set_selpar_logit_Age2_sCT(1);  
  
  //Discards
  //selpar_age1logit_D=set_selpar_age1logit_D(1);
  selpar_logit_Age2_cHL_D2=set_selpar_logit_Age2_cHL_D2(1);  
  selpar_A502_cHL_D3=set_selpar_A502_cHL_D3(1);
  sel_cHL_D_2=set_sel_cHL_D_2;
  selpar_A50_cHL_D3=set_selpar_A50_cHL_D3;  
  selpar_slope_cHL_D3=set_selpar_slope_cHL_D3;
  selpar_slope2_cHL_D3=set_selpar_slope2_cHL_D3;
  
  selpar_Age1_rHB_D3=set_selpar_Age1_rHB_D3;
  selpar_Age1_rHB_D4=set_selpar_Age1_rHB_D4;
  //selpar_Age2_rHB_D3=set_selpar_Age2_rHB_D3;
  //selpar_Age2_rHB_D4=set_selpar_Age2_rHB_D4;
  
  selpar_logit_Age2_rHB_D3=set_selpar_logit_Age2_rHB_D3(1);
  selpar_logit_Age2_rHB_D4=set_selpar_logit_Age2_rHB_D4(1);
  
  sel_rHB_D_2=set_sel_rHB_D_2;
  sel_rHB_D_3=set_sel_rHB_D_3;
  sel_rHB_D_4=set_sel_rHB_D_4;
  
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
 //maturity_m=maturity_m_obs;
 
 prop_f=obs_prop_f;
 prop_m=1.0-obs_prop_f;
 
 //kc p_lenc_cHL=set_p_lenc_cHL;
 //kc p_lenc_cOT=set_p_lenc_cOT;
 //kc p_lenc_rHB2=set_p_lenc_rHB2;
 //kc p_lenc_rHB3=set_p_lenc_rHB3;
 //kc p_lenc_rHB4=set_p_lenc_rHB4;  
 //kc p_lenc_rGN2=set_p_lenc_rGN2;
 //kc p_lenc_rGN3=set_p_lenc_rGN3;
 //kc p_lenc_rGN4=set_p_lenc_rGN4;
   
 //kc p_lenc_cHL_D=set_p_lenc_cHL_D;
 //kc p_lenc_rHB_D2=set_p_lenc_rHB_D2;
 //kc p_lenc_rHB_D3=set_p_lenc_rHB_D3; 
 //k p_lenc_rHB_D4=set_p_lenc_rHB_D4;  
 
//Fill in sample sizes of comps, possibly sampled in nonconsec yrs 
//Used primarily for output in R object   

      nsamp_cHL_lenc_allyr=missing;
      nsamp_cOT_lenc_allyr=missing;  
     //nsamp_rHB_lenc_allyr=missing;
      nsamp_rHB_D_lenc_allyr=missing;
      nsamp_cHL_D_lenc_allyr=missing;	  
     // nsamp_rGN_lenc_allyr=missing;
      nsamp_sFT_lenc_allyr=missing;
      //nsamp_sCT_lenc_allyr=missing;	  
      nsamp_cHL_agec_allyr=missing;
      nsamp_rHB_agec_allyr=missing;
      nsamp_rGN_agec_allyr=missing;
      nsamp_sCT_agec_allyr=missing;  	  
      
      nfish_cHL_lenc_allyr=missing;
      nfish_cOT_lenc_allyr=missing;  
      //nfish_rHB_lenc_allyr=missing;
      nfish_rHB_D_lenc_allyr=missing;  
      nfish_cHL_D_lenc_allyr=missing;
      //nfish_rGN_lenc_allyr=missing;
      nfish_sFT_lenc_allyr=missing;
      //nfish_sCT_lenc_allyr=missing;	  
      nfish_cHL_agec_allyr=missing;
      nfish_rHB_agec_allyr=missing;
      nfish_rGN_agec_allyr=missing;  
      nfish_sCT_agec_allyr=missing;  	  
   
      for (iyear=1; iyear<=nyr_lenc_cHL; iyear++)
         {if (nsamp_lenc_cHL(iyear)>=minSS_lenc_cHL)
           {nsamp_cHL_lenc_allyr(yrs_lenc_cHL(iyear))=nsamp_lenc_cHL(iyear);
            nfish_cHL_lenc_allyr(yrs_lenc_cHL(iyear))=nfish_lenc_cHL(iyear);}}
      for (iyear=1; iyear<=nyr_lenc_cOT; iyear++)
         {if (nsamp_lenc_cOT(iyear)>=minSS_lenc_cOT)
           {nsamp_cOT_lenc_allyr(yrs_lenc_cOT(iyear))=nsamp_lenc_cOT(iyear);
            nfish_cOT_lenc_allyr(yrs_lenc_cOT(iyear))=nfish_lenc_cOT(iyear);}}			
//	  for (iyear=1; iyear<=nyr_rHB_lenc; iyear++)                           
//         {if (nsamp_rHB_lenc(iyear)>=minSS_rHB_lenc)
//            {nsamp_rHB_lenc_allyr(yrs_rHB_lenc(iyear))=nsamp_rHB_lenc(iyear);
//             nfish_rHB_lenc_allyr(yrs_rHB_lenc(iyear))=nfish_rHB_lenc(iyear);}}
      for (iyear=1; iyear<=nyr_lenc_rHB_D; iyear++)                           
         {if (nsamp_lenc_rHB_D(iyear)>=minSS_lenc_rHB_D)
            {nsamp_rHB_D_lenc_allyr(yrs_lenc_rHB_D(iyear))=nsamp_lenc_rHB_D(iyear);
             nfish_rHB_D_lenc_allyr(yrs_lenc_rHB_D(iyear))=nfish_lenc_rHB_D(iyear);}}
			 
	  for (iyear=1; iyear<=nyr_lenc_cHL_D; iyear++)                           
         {if (nsamp_lenc_cHL_D(iyear)>=minSS_lenc_cHL_D)
            {nsamp_cHL_D_lenc_allyr(yrs_lenc_cHL_D(iyear))=nsamp_lenc_cHL_D(iyear);
             nfish_cHL_D_lenc_allyr(yrs_lenc_cHL_D(iyear))=nfish_lenc_cHL_D(iyear);}}
			 
//	  for (iyear=1; iyear<=nyr_rGN_lenc; iyear++)                           
//         {if (nsamp_rGN_lenc(iyear)>=minSS_rGN_lenc)
//            {nsamp_rGN_lenc_allyr(yrs_rGN_lenc(iyear))=nsamp_rGN_lenc(iyear);
//             nfish_rGN_lenc_allyr(yrs_rGN_lenc(iyear))=nfish_rGN_lenc(iyear);}}
      for (iyear=1; iyear<=nyr_lenc_sFT; iyear++)                           
         {if (nsamp_lenc_sFT(iyear)>=minSS_lenc_sFT)
            {nsamp_sFT_lenc_allyr(yrs_lenc_sFT(iyear))=nsamp_lenc_sFT(iyear);
             nfish_sFT_lenc_allyr(yrs_lenc_sFT(iyear))=nfish_lenc_sFT(iyear);}}
	  //for (iyear=1; iyear<=nyr_sCT_lenc; iyear++)                           
      //   {if (nsamp_sCT_lenc(iyear)>=minSS_sCT_lenc)
      //      {nsamp_sCT_lenc_allyr(yrs_sCT_lenc(iyear))=nsamp_sCT_lenc(iyear);
      //       nfish_sCT_lenc_allyr(yrs_sCT_lenc(iyear))=nfish_sCT_lenc(iyear);}}
      //
	  for (iyear=1; iyear<=nyr_agec_cHL; iyear++)
         {if (nsamp_agec_cHL(iyear)>=minSS_agec_cHL)
           {nsamp_cHL_agec_allyr(yrs_agec_cHL(iyear))=nsamp_agec_cHL(iyear);
            nfish_cHL_agec_allyr(yrs_agec_cHL(iyear))=nfish_agec_cHL(iyear);}}
      for (iyear=1; iyear<=nyr_agec_rHB; iyear++)
          {if (nsamp_agec_rHB(iyear)>=minSS_agec_rHB)
            {nsamp_rHB_agec_allyr(yrs_agec_rHB(iyear))=nsamp_agec_rHB(iyear);
             nfish_rHB_agec_allyr(yrs_agec_rHB(iyear))=nfish_agec_rHB(iyear);}} 
      for (iyear=1; iyear<=nyr_agec_rGN; iyear++)  
         {if (nsamp_agec_rGN(iyear)>=minSS_agec_rGN)
           {nsamp_rGN_agec_allyr(yrs_agec_rGN(iyear))=nsamp_agec_rGN(iyear);
             nfish_rGN_agec_allyr(yrs_agec_rGN(iyear))=nfish_agec_rGN(iyear);}}  
	  for (iyear=1; iyear<=nyr_agec_sCT; iyear++)  
          {if (nsamp_agec_sCT(iyear)>=minSS_agec_sCT)
            {nsamp_sCT_agec_allyr(yrs_agec_sCT(iyear))=nsamp_agec_sCT(iyear);
             nfish_sCT_agec_allyr(yrs_agec_sCT(iyear))=nfish_agec_sCT(iyear);}} 

             
//fill in Fs for msy and per-recruit analyses
  F_msy(1)=0.0;  
  for (ff=2;ff<=n_iter_msy;ff++) {F_msy(ff)=F_msy(ff-1)+iter_inc_msy;}
  F_spr(1)=0.0;  
  for (ff=2;ff<=n_iter_spr;ff++) {F_spr(ff)=F_spr(ff-1)+iter_inc_spr;}


//fill in F's, Catch matrices, and log rec dev with zero's
  F_cTW.initialize(); L_cTW_num.initialize();
  F_cHL.initialize(); L_cHL_num.initialize();
  F_cOT.initialize(); L_cOT_num.initialize();
  F_rHB.initialize(); L_rHB_num.initialize();
  F_rGN.initialize(); L_rGN_num.initialize();
  
  F_cHL_D.initialize(); D_cHL_num.initialize();
  F_rHB_D.initialize(); D_rHB_num.initialize();
  F_rGN_D.initialize(); D_rGN_num.initialize();

  F_cHL_out.initialize();
  F_cOT_out.initialize();
  F_cTW_out.initialize();
  F_rHB_out.initialize();
  F_rGN_out.initialize();
  F_cHL_D_out.initialize();
  F_rHB_D_out.initialize();
  F_rGN_D_out.initialize();

  sel_cHL.initialize();
  //sel_cOT.initialize();
  //sel_cTW.initialize();
  sel_rHB.initialize();
  sel_rGN.initialize();
  //sel_D.initialize();
  sel_sFT.initialize();
  sel_sCT.initialize();
  //sel_cHL_D.initialize();  //kc
  sel_rHB_D.initialize();  //kc

  sel_cHL_block1.initialize(); 
  sel_cHL_block2.initialize();
  //sel_cHL_block3.initialize();
  //sel_cOT_block1.initialize();
  //sel_cOT_block2.initialize();
  //sel_cOT_block3.initialize();    
  sel_rHB_block1.initialize();
  sel_rHB_block2.initialize();
  sel_rHB_block3.initialize();
  sel_rHB_block4.initialize();
  sel_rGN_block1.initialize();
  sel_rGN_block2.initialize();
  sel_rGN_block3.initialize();
  sel_rGN_block4.initialize();
  //sel_D_block1.initialize();   
  //sel_D_block2.initialize();  
  //sel_D_block3.initialize();
  //sel_sFT_vec.initialize();
  sel_sCT_vec.initialize();
  prob_belowsizelim_block2.initialize();//kc uncommented
  prob_belowsizelim_block3.initialize();//kc uncommented
  prob_belowsizelim_block4.initialize();//kc added
    
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
 
 
FUNCTION get_length_weight_at_age
	//population total length in mm
    //compute mean length (mm TL) and weight (whole) at age
    meanlen_TL=Linf*(1.0-mfexp(-K*(agebins-t0+0.5)));   //total length in mm 
    meanlen_FL=lenpar_TL2FL_a+lenpar_TL2FL_b*meanlen_TL; //fork length in mm
 
    wgt_g=wgtpar_a*pow(meanlen_TL,wgtpar_b);
    wgt_kg=g2kg*wgt_g;
//KC    wgt_g=wgt_kg/g2kg;                                    //convert wgt in kg to weight in g    
//KC    wgt_kg=wgtpar_a*pow(meanlen_TL,wgtpar_b);             //whole wgt in kg 
	wgt_mt=wgt_g*g2mt;                                    //convert weight in g to weight in mt
    wgt_klb=mt2klb*wgt_mt;                                //1000 lb of whole wgt
    wgt_lb=mt2lb*wgt_mt;                                  //lb of whole wgt

//KC	fecundity=fecpar_batches*fecpar_a*pow(meanlen_FL,fecpar_b)/fecpar_scale; //annual egg production of a mature female at age in units of fecpar_scale
//KC    fecundity=fecpar_batches*(fecpar_a+fecpar_b*meanlen_FL)/fecpar_scale;
    fecundity=(fecpar_a+fecpar_b*meanlen_FL)/fecpar_scale;  //kc
    
FUNCTION get_reprod 
 
	//reprod=elem_prod(elem_prod(prop_f,maturity_f),fecundity); //kc
        reprod=elem_prod(elem_prod(elem_prod(prop_f,maturity_f),fecundity),fecpar_batches);  //kc  
	//reprod=elem_prod((elem_prod(prop_f,maturity_f)+elem_prod(prop_m,maturity_m)),wgt_mt);
 
FUNCTION get_length_at_age_dist
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
	
	zscore_lsizelim3=(sizelim3-meanlen_TL(iage)) / len_sd(iage);
	cprob_lsizelim3=cumd_norm(zscore_lsizelim3);                 //includes any probability mass below zero
	prob_belowsizelim_block4(iage)=	cprob_lsizelim3-cprob_lzero; //removes any probability mass below zero
	
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
   lenprob_cHL=lenprob;
   lenprob_cOT=lenprob; 
   //lenprob_rHB=lenprob;
   lenprob_rHB_D=lenprob; 
   lenprob_cHL_D=lenprob;
   //lenprob_rGN=lenprob;
   lenprob_sFT=lenprob; 
   //lenprob_sCT=lenprob; 
  
//compute matrix of length at age, based on the normal distribution
//  for (iage=1;iage<=nages;iage++)
//  {
//    //len_cv(iage)=mfexp(log_len_cv+log_len_cv_dev(iage));
//    len_cv(iage)=len_cv_val;
//    for (ilen=1;ilen<=nlenbins;ilen++)
//    {
//      lenprob(iage,ilen)=(mfexp(-(square(lenbins(ilen)-meanlen_TL(iage))/
//      (2.*square(len_cv(iage)*meanlen_TL(iage)))))/(sqrt2pi*len_cv(iage)*meanlen_TL(iage)));
//    }
//    for (ilen=1;ilen<=nlenbins_plus;ilen++)
//    {
//      lenprob_plus(iage,ilen)=(mfexp(-(square(lenbins_plus(ilen)-meanlen_TL(iage))/
//      (2.*square(len_cv(iage)*meanlen_TL(iage)))))/(sqrt2pi*len_cv(iage)*meanlen_TL(iage)));
//    } 

//    lenprob(iage)(nlenbins)=lenprob(iage)(nlenbins)+sum(lenprob_plus(iage)); //add mass to plus group
//    lenprob(iage)/=sum(lenprob(iage)); //standardize to account for truncated normal (i.e., no sizes<smallest)
//  }
  
//  //fishery specific length probs, assumed normal prior to size limits
    
  //kc lenprob_cHL1=lenprob;
  //kc lenprob_cHL2=lenprob; //values will be adjusted based on size limit 
  //kc lenprob_cOT1=lenprob; 
  //kc lenprob_cOT2=lenprob; //values will be adjusted based on size limit 
  //kc lenprob_rHB1=lenprob; 
  //kc lenprob_rHB2=lenprob; //values will be adjusted based on size limit 
  //kc lenprob_rHB3=lenprob; //values will be adjusted based on size limit 
  //kc lenprob_rHB4=lenprob; //values will be adjusted based on size limit  
  //kc lenprob_rGN1=lenprob; 
  //kc lenprob_rGN2=lenprob; //values will be adjusted based on size limit 
  //kc lenprob_rGN3=lenprob; //values will be adjusted based on size limit 
  //kc lenprob_rGN4=lenprob; //values will be adjusted based on size limit 
  //kc lenprob_sFT=lenprob; 
  
  //kc lenprob_cHL_D2=lenprob; //values will be adjusted based on size limit 
  //kc lenprob_cHL_D3=lenprob;  
  //kc lenprob_rHB_D2=lenprob; //values will be adjusted based on size limit 
  //kc lenprob_rHB_D3=lenprob; //values will be adjusted based on size limit 
  //kc lenprob_rHB_D4=lenprob; //values will be adjusted based on size limit   
  
//  for (iage=1;iage<=nages;iage++)
//  {
//    for (ilen=1;ilen<=nlenbins;ilen++)
//    {
//        if (lenbins(ilen) < sizelim1)
//        {
//          lenprob_rHB2(iage,ilen)=p_lenc_rHB2*lenprob(iage,ilen);
//          lenprob_rGN2(iage,ilen)=p_lenc_rGN2*lenprob(iage,ilen);
          
//        }
//        if (lenbins(ilen) < sizelim2)
//        {
//          lenprob_rHB3(iage,ilen)=p_lenc_rHB3*lenprob(iage,ilen);
//          lenprob_rGN3(iage,ilen)=p_lenc_rGN3*lenprob(iage,ilen);
          
//        }
//        if (lenbins(ilen) < sizelim3)
//        {
//          lenprob_cHL2(iage,ilen)=p_lenc_cHL*lenprob(iage,ilen);
//          lenprob_cOT2(iage,ilen)=p_lenc_cOT*lenprob(iage,ilen);
//          lenprob_rHB4(iage,ilen)=p_lenc_rHB4*lenprob(iage,ilen);
//          lenprob_rGN4(iage,ilen)=p_lenc_rGN4*lenprob(iage,ilen);
          
//        }

//        if (lenbins(ilen) > sizelim1)
//        {
//          lenprob_rHB_D2(iage,ilen)=p_lenc_rHB_D2*lenprob(iage,ilen);
//        }


//        if (lenbins(ilen) > sizelim2)
//        {
//          lenprob_rHB_D3(iage,ilen)=p_lenc_rHB_D3*lenprob(iage,ilen);
//        }


//        if (lenbins(ilen) > sizelim3)
//        {
//          lenprob_cHL_D2(iage,ilen)=p_lenc_cHL_D*lenprob(iage,ilen);
//          lenprob_rHB_D4(iage,ilen)=p_lenc_rHB_D4*lenprob(iage,ilen);
//        }


//    }
    
    
//    lenprob_cHL2(iage)/=sum(lenprob_cHL2(iage));   //standardize
//    lenprob_cOT2(iage)/=sum(lenprob_cOT2(iage));   //standardize
//    lenprob_rHB2(iage)/=sum(lenprob_rHB2(iage));       //standardize
//    lenprob_rGN2(iage)/=sum(lenprob_rGN2(iage)); //standardize
//    lenprob_rHB3(iage)/=sum(lenprob_rHB3(iage));       //standardize
//    lenprob_rGN3(iage)/=sum(lenprob_rGN3(iage)); //standardize
//    lenprob_rHB4(iage)/=sum(lenprob_rHB4(iage));       //standardize
//    lenprob_rGN4(iage)/=sum(lenprob_rGN4(iage)); //standardize
    
//    lenprob_cHL_D2(iage)/=sum(lenprob_cHL_D2(iage)); //standardize
//    lenprob_rHB_D2(iage)/=sum(lenprob_rHB_D2(iage)); //standardize
//    lenprob_rHB_D3(iage)/=sum(lenprob_rHB_D3(iage)); //standardize
//    lenprob_rHB_D4(iage)/=sum(lenprob_rHB_D4(iage)); //standardize
        
//  }
  
FUNCTION get_weight_at_age_landings  ///***in whole weight
  
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    len_cHL_mm(iyear)=meanlen_TL;  
    wholewgt_cHL_klb(iyear)=wgt_klb; 
    len_cOT_mm(iyear)=meanlen_TL;  
    wholewgt_cOT_klb(iyear)=wgt_klb;
    len_cTW_mm(iyear)=meanlen_TL;  
    wholewgt_cTW_klb(iyear)=wgt_klb;
    len_rHB_mm(iyear)=meanlen_TL;
    wholewgt_rHB_klb(iyear)=wgt_klb;
    len_rGN_mm(iyear)=meanlen_TL;
    wholewgt_rGN_klb(iyear)=wgt_klb;

    len_cHL_D_mm(iyear)=meanlen_TL;  
    wholewgt_cHL_D_klb(iyear)=wgt_klb;
    len_rHB_D_mm(iyear)=meanlen_TL;
    wholewgt_rHB_D_klb(iyear)=wgt_klb;
    len_rGN_D_mm(iyear)=meanlen_TL;
    wholewgt_rGN_D_klb(iyear)=wgt_klb;
  
  //rHB and rec discards assumed same size as period one, consistent with selectivity assumptions
  //kc for (iage=1;iage<=nages; iage++)
  //kc      {len_rHB_D_mm(iyear,iage)=sum(elem_prod(lenprob_rHB_D2(iage),lenbins));}
  //kc  wholewgt_rHB_D_klb(iyear)=g2klb*wgtpar_a*pow(len_rHB_D_mm(iyear),wgtpar_b);
  //kc }  
 
 //kc for (iyear=(endyr_selex_phase1+1); iyear<=endyr_selex_phase2; iyear++)
 //kc  {
 //kc   for (iage=1;iage<=nages; iage++)
 //kc   {
 //kc   len_cHL_mm(iyear,iage)=sum(elem_prod(lenprob_cHL2(iage),lenbins)); 
 //kc   len_cOT_mm(iyear,iage)=sum(elem_prod(lenprob_cOT2(iage),lenbins)); 
 //kc   len_rHB_mm(iyear,iage)=sum(elem_prod(lenprob_rHB2(iage),lenbins));     
 //kc   len_rGN_mm(iyear,iage)=sum(elem_prod(lenprob_rGN2(iage),lenbins));     
 //kc   len_cHL_D_mm(iyear,iage)=sum(elem_prod(lenprob_cHL_D2(iage),lenbins)); 
 //kc   len_rHB_D_mm(iyear,iage)=sum(elem_prod(lenprob_rHB_D2(iage),lenbins));       
 //kc   }
 //kc   wholewgt_cHL_klb(iyear)=g2klb*wgtpar_a*pow(len_cHL_mm(iyear),wgtpar_b);
 //kc   wholewgt_cOT_klb(iyear)=g2klb*wgtpar_a*pow(len_cOT_mm(iyear),wgtpar_b);
 //kc   wholewgt_rHB_klb(iyear)=g2klb*wgtpar_a*pow(len_rHB_mm(iyear),wgtpar_b);
 //kc   wholewgt_rGN_klb(iyear)=g2klb*wgtpar_a*pow(len_rGN_mm(iyear),wgtpar_b);
 //kc   wholewgt_cHL_D_klb(iyear)=g2klb*wgtpar_a*pow(len_cHL_D_mm(iyear),wgtpar_b);  
 //kc   wholewgt_rHB_D_klb(iyear)=g2klb*wgtpar_a*pow(len_rHB_D_mm(iyear),wgtpar_b);           
 // }
  
 // for (iyear=(endyr_selex_phase2+1); iyear<=endyr_selex_phase3; iyear++)
  //{
 //   for (iage=1;iage<=nages; iage++)
 //   {
 //   len_cHL_mm(iyear,iage)=sum(elem_prod(lenprob_cHL2(iage),lenbins)); 
 //   len_cOT_mm(iyear,iage)=sum(elem_prod(lenprob_cOT2(iage),lenbins)); 
 //   len_rHB_mm(iyear,iage)=sum(elem_prod(lenprob_rHB3(iage),lenbins));     
 //   len_rGN_mm(iyear,iage)=sum(elem_prod(lenprob_rGN3(iage),lenbins));     
 //   len_cHL_D_mm(iyear,iage)=sum(elem_prod(lenprob_cHL_D2(iage),lenbins)); 
 //   len_rHB_D_mm(iyear,iage)=sum(elem_prod(lenprob_rHB_D3(iage),lenbins));
 //   }
 //   wholewgt_cHL_klb(iyear)=g2klb*wgtpar_a*pow(len_cHL_mm(iyear),wgtpar_b);
 //   wholewgt_cOT_klb(iyear)=g2klb*wgtpar_a*pow(len_cOT_mm(iyear),wgtpar_b); 
 //   wholewgt_rHB_klb(iyear)=g2klb*wgtpar_a*pow(len_rHB_mm(iyear),wgtpar_b);
 //   wholewgt_rGN_klb(iyear)=g2klb*wgtpar_a*pow(len_rGN_mm(iyear),wgtpar_b);       
 //   wholewgt_cHL_D_klb(iyear)=g2klb*wgtpar_a*pow(len_cHL_D_mm(iyear),wgtpar_b); 
 //   wholewgt_rHB_D_klb(iyear)=g2klb*wgtpar_a*pow(len_rHB_D_mm(iyear),wgtpar_b);                
 // }
  
//  for (iyear=(endyr_selex_phase3+1); iyear<=endyr; iyear++)
//  {
//    for (iage=1;iage<=nages; iage++)
//    {
//    len_cHL_mm(iyear,iage)=sum(elem_prod(lenprob_cHL2(iage),lenbins)); 
//    len_cOT_mm(iyear,iage)=sum(elem_prod(lenprob_cOT2(iage),lenbins)); 
//    len_rHB_mm(iyear,iage)=sum(elem_prod(lenprob_rHB4(iage),lenbins));     
//    len_rGN_mm(iyear,iage)=sum(elem_prod(lenprob_rGN4(iage),lenbins));         
//    len_cHL_D_mm(iyear,iage)=sum(elem_prod(lenprob_cHL_D2(iage),lenbins)); 
//    len_rHB_D_mm(iyear,iage)=sum(elem_prod(lenprob_rHB_D4(iage),lenbins));
//    }
//    wholewgt_cHL_klb(iyear)=g2klb*wgtpar_a*pow(len_cHL_mm(iyear),wgtpar_b);
//    wholewgt_cOT_klb(iyear)=g2klb*wgtpar_a*pow(len_cOT_mm(iyear),wgtpar_b); 
//    wholewgt_rHB_klb(iyear)=g2klb*wgtpar_a*pow(len_rHB_mm(iyear),wgtpar_b);
//   wholewgt_rGN_klb(iyear)=g2klb*wgtpar_a*pow(len_rGN_mm(iyear),wgtpar_b);       
//    wholewgt_cHL_D_klb(iyear)=g2klb*wgtpar_a*pow(len_cHL_D_mm(iyear),wgtpar_b);  
//    wholewgt_rHB_D_klb(iyear)=g2klb*wgtpar_a*pow(len_rHB_D_mm(iyear),wgtpar_b);           
//  }

//  //overwrite values for cGN discards
//  for (iyear=(endyr_2008+1); iyear<=endyr; iyear++)
//  {
//    for (iage=1;iage<=nages; iage++)
//    {len_cHL_D_mm(iyear,iage)=sum(elem_prod(lenprob_cHL_D3(iage),lenbins)); }
//    wholewgt_cHL_D_klb(iyear)=g2klb*wgtpar_a*pow(len_cHL_D_mm(iyear),wgtpar_b);  
  }


//  //rGN_D assumed same as rHB_D
//  len_rGN_D_mm=len_rHB_D_mm;
//  wholewgt_rGN_D_klb=wholewgt_rHB_D_klb;                  



 
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
  bpr_F0=sum(elem_prod(N_bpr_F0,wgt_mt));    



FUNCTION get_selectivity
    
	//sel_sFT_vec=logistic_double(agebins, selpar_A50_sFT, selpar_slope_sFT, selpar_A502_sFT,selpar_slope2_sFT);
	
  selpar_A50_cOT2=selpar_A50_cHL2; //setting all commercial other selectivity parameters to cHL
  selpar_slope_cOT2=selpar_slope_cHL2;
  selpar_A502_cOT2=selpar_A502_cHL2;
  selpar_slope2_cOT2=selpar_slope2_cHL2;
  
  //sel_sCT_vec=logistic_double(agebins, selpar_A50_sCT, selpar_slope_sCT, selpar_A502_sCT,selpar_slope2_sCT);
  sel_sCT_vec=logistic_exponential(agebins, selpar_A50_sCT, selpar_slope_sCT, selpar_A502_sCT,selpar_slope2_sCT);
  //KCselpar_Age1_sCT=1.0/(1.0+mfexp(-selpar_logit_Age1_sCT));
  //KCselpar_Age2_sCT=1.0/(1.0+mfexp(-selpar_logit_Age2_sCT));
  //KCsel_sCT_vec(1)=selpar_Age1_sCT;
  //KC sel_sCT_vec(2)=selpar_Age2_sCT;
  
  sel_cHL_block1=logistic(agebins, selpar_A50_cHL1, selpar_slope_cHL1);
  sel_cHL_block2=logistic(agebins, selpar_A50_cHL2, selpar_slope_cHL2);
  //sel_cHL_block1=logistic_double(agebins, selpar_A50_cHL1, selpar_slope_cHL1, selpar_A502_cHL1, selpar_slope2_cHL1);
  //sel_cHL_block2=logistic_double(agebins, selpar_A50_cHL2, selpar_slope_cHL2, selpar_A502_cHL2, selpar_slope2_cHL2);
  //sel_cHL_block1=sel_cHL_block2;
  
  sel_cOT_block2=logistic(agebins, selpar_A50_cOT2, selpar_slope_cOT2);
    
  sel_rHB_block1=logistic(agebins, selpar_A50_rHB1, selpar_slope_rHB1);
  sel_rHB_block2=logistic(agebins, selpar_A50_rHB2, selpar_slope_rHB2);
  sel_rHB_block3=logistic(agebins, selpar_A50_rHB3, selpar_slope_rHB3);
  sel_rHB_block4=logistic(agebins, selpar_A50_rHB4, selpar_slope_rHB4);
  
  //sel_rGN_block1=sel_rHB_block1;
  //sel_rGN_block2=sel_rHB_block2;
  sel_rGN_block2=logistic(agebins, selpar_A50_rGN2, selpar_slope_rGN2);//KC
  sel_rGN_block3=logistic(agebins, selpar_A50_rGN3, selpar_slope_rGN3);
  sel_rGN_block4=logistic(agebins, selpar_A50_rGN4, selpar_slope_rGN4);
  
  //-------------Discard selectivities------------------------------
  selpar_Age2_rHB_D3=1.0/(1.0+mfexp(-selpar_logit_Age2_rHB_D3));
  selpar_Age2_rHB_D4=1.0/(1.0+mfexp(-selpar_logit_Age2_rHB_D4));
  
  selpar_Age2_cHL_D2=1.0/(1.0+mfexp(-selpar_logit_Age2_cHL_D2));
  
  sel_cHL_D_2(1)=selpar_Age1_cHL_D2;
  sel_cHL_D_2(2)=selpar_Age2_cHL_D2; 
  //cout << sel_cHL_D_2 << endl;
  sel_rHB_D_2(1)=selpar_Age1_rHB_D3;
  sel_rHB_D_2(2)=selpar_Age2_rHB_D3;
  sel_rHB_D_3(1)=selpar_Age1_rHB_D3;
  sel_rHB_D_3(2)=selpar_Age2_rHB_D3;
  sel_rHB_D_4(1)=selpar_Age1_rHB_D4;
  sel_rHB_D_4(2)=selpar_Age2_rHB_D4;
  
  //sel_cHL_D_3=logistic_double(agebins, selpar_A50_cHL_D3, selpar_slope_cHL_D3, selpar_A502_cHL_D3, selpar_slope2_cHL_D3);
  sel_cHL_D_3=logistic_exponential(agebins, selpar_A50_cHL_D3, selpar_slope_cHL_D3, selpar_A502_cHL_D3, selpar_slope2_cHL_D3);
  sel_cHL_D_3(1)=set_sel_cHL_D_2(1);
  sel_cHL_D_3(2)=set_sel_cHL_D_2(2);
  //cout <<"sel_cHL_D_2(2)"<< sel_cHL_D_2(2)<<endl;
  //cout <<"sel_cHL_D_3(2)"<< sel_cHL_D_3(2)<<endl;
  
  //-----------fill in years--------------------------------------------
  //selpar_age1_D=mfexp(selpar_age1logit_D)/(1.0+mfexp(selpar_age1logit_D));
  //sel_D_block2(1)=selpar_age1_D; sel_D_block2(2)=1.0; //sel_D_block2(3,nages)=prob_belowsizelim_block2(3,nages);
  //sel_D_block3(1)=selpar_age1_D; sel_D_block3(2)=1.0; //sel_D_block3(3,nages)=prob_belowsizelim_block3(3,nages);
  //sel_D_block1=sel_D_block2;
    
  for (iyear=styr; iyear<=endyr; iyear++)
  {  //time-invariant selectivities
     sel_sFT(iyear)=sel_sFT_vec;
     sel_sCT(iyear)=sel_sCT_vec;
	 
  //Period 1: MRFSS same as rHB.  Historical trawl same as CMB.  
  if (iyear<=endyr_selex_phase1)
   {     
    sel_cTW(iyear)=sel_cOT_block1;
	sel_cHL(iyear)=sel_cHL_block1;
	sel_rHB(iyear)=sel_rHB_block1;
    sel_rGN(iyear)=sel_rHB_block1;
    sel_cOT(iyear)=sel_cOT_block1;
	sel_rHB_D(iyear)=sel_rHB_D_2;
    sel_rGN_D(iyear)=sel_rHB_D_2; 
		 //sel_D(iyear)=sel_D_block1;
   }
   
    //Period 2: rGN same as rHB
  
  if (iyear>endyr_selex_phase1 & iyear<=endyr_selex_phase2)
  
   {
    sel_cHL(iyear)=sel_cHL_block2;
	sel_rHB(iyear)=sel_rHB_block2;
        sel_rGN(iyear)=sel_rHB_block2;
    //KC	sel_rGN(iyear)=sel_rGN_block2;
 	sel_cHL_D(iyear)=sel_cHL_D_2;
	//cout<<sel_cHL_D<<endl;
    sel_rHB_D(iyear)=sel_rHB_D_2;
    sel_rGN_D(iyear)=sel_rHB_D_2;
		 //sel_D(iyear)=sel_D_block2;
   }
  if (iyear>1988 & iyear<=endyr_selex_phase2)
  {sel_cOT(iyear)=sel_cOT_block2;}

  //Period 3 (commercial doesn't change from period 2)
   if (iyear>endyr_selex_phase2 & iyear<=endyr_selex_phase3)
   {   
    sel_cHL(iyear)=sel_cHL_block2;
	sel_cOT(iyear)=sel_cOT_block2;
    sel_rHB(iyear)=sel_rHB_block3;
    //KC sel_rGN(iyear)=sel_rGN_block3;
     sel_rGN(iyear)=sel_rHB_block3;  //KC
     
 	sel_cHL_D(iyear)=sel_cHL_D_2;     
    sel_rHB_D(iyear)=sel_rHB_D_3;
    sel_rGN_D(iyear)=sel_rHB_D_3;
		 //sel_D(iyear)=sel_D_block3;
   }  
   if (iyear>endyr_selex_phase3)
   {
         sel_cHL(iyear)=sel_cHL_block2; 
         sel_cOT(iyear)=sel_cOT_block2;    
         sel_rHB(iyear)=sel_rHB_block4;
         //KC sel_rGN(iyear)=sel_rGN_block4;
         sel_rGN(iyear)=sel_rHB_block4; //KC

		 //sel_D(iyear)=sel_D_block3;
         sel_rHB_D(iyear)=sel_rHB_D_4;
         sel_rGN_D(iyear)=sel_rHB_D_4;
         
         if (iyear<=endyr_2008){sel_cHL_D(iyear)=sel_cHL_D_2;} 
		//cout<<sel_cHL_D_2<<endl;		 
         if (iyear>endyr_2008){sel_cHL_D(iyear)=sel_cHL_D_3;}
		//cout<<sel_cHL_D<<endl;
                 sel_cHL_D_3(1)=set_sel_cHL_D_2(1);
                 sel_cHL_D_3(2)=set_sel_cHL_D_2(2);
	}
  
  }   
	
	  
FUNCTION get_mortality
  Fsum.initialize();
  Fapex.initialize();
  F.initialize();
  //initialization F is avg from first 3 yrs of observed landings
  log_F_dev_init_cHL=sum(log_dev_F_L_cHL(styr_L_cHL,(styr_L_cHL+2)))/3.0;  
  log_F_dev_init_cOT=sum(log_dev_F_L_cOT(styr_L_cOT,(styr_L_cOT+2)))/3.0;
  log_F_dev_init_cTW=log_dev_F_L_cTW(styr_L_cTW); 
  log_F_dev_init_rHB=sum(log_dev_F_L_rHB(styr_L_rHB,(styr_L_rHB+2)))/3.0;         
  log_F_dev_init_rGN=sum(log_dev_F_L_rGN(styr_L_rGN,(styr_L_rGN+2)))/3.0;         
  
  for (iyear=styr; iyear<=endyr; iyear++) 
  {
    if(iyear>=styr_L_cHL & iyear<=endyr_L_cHL) //spans full time series
		{F_cHL_out(iyear)=mfexp(log_avg_F_L_cHL+log_dev_F_L_cHL(iyear));}     
    F_cHL(iyear)=sel_cHL(iyear)*F_cHL_out(iyear);
    Fsum(iyear)+=F_cHL_out(iyear);
    

	if(iyear>=styr_L_cOT & iyear<=endyr_L_cOT) //spans full time series
		{F_cOT_out(iyear)=mfexp(log_avg_F_L_cOT+log_dev_F_L_cOT(iyear));}     
    F_cOT(iyear)=sel_cOT(iyear)*F_cOT_out(iyear);
    Fsum(iyear)+=F_cOT_out(iyear);
    

	if(iyear>=styr_L_cTW & iyear<=endyr_L_cTW) //spans full time series
		{F_cTW_out(iyear)=mfexp(log_avg_F_L_cTW+log_dev_F_L_cTW(iyear));}    
    F_cTW(iyear)=sel_cTW(iyear)*F_cTW_out(iyear);
    Fsum(iyear)+=F_cTW_out(iyear);
    

    if(iyear>=styr_L_rHB & iyear<=endyr_L_rHB) //spans full time series
		{F_rHB_out(iyear)=mfexp(log_avg_F_L_rHB+log_dev_F_L_rHB(iyear));}     
    F_rHB(iyear)=sel_rHB(iyear)*F_rHB_out(iyear);
    Fsum(iyear)+=F_rHB_out(iyear);
    
   
    if(iyear>=styr_L_rGN & iyear<=endyr_L_rGN) //starts in 1981
		{F_rGN_out(iyear)=mfexp(log_avg_F_L_rGN+log_dev_F_L_rGN(iyear));}    
    if (iyear<styr_L_rGN)
		{F_rGN_out(iyear)=mfexp(log_avg_F_L_rGN+log_F_dev_init_rGN);}
	F_rGN(iyear)=sel_rGN(iyear)*F_rGN_out(iyear); 
    Fsum(iyear)+=F_rGN_out(iyear);

	
    log_F_avgdev_cHL_D=sum(log_dev_F_D_cHL(styr_D_cHL,endyr_D_cHL))/(endyr_D_cHL-styr_D_cHL+1.0);
    log_F_avgdev_rHB_D=sum(log_dev_F_D_rHB(styr_D_rHB,endyr_D_rHB))/(endyr_D_rHB-styr_D_rHB+1.0);   
    log_F_avgdev_rGN_D=sum(log_dev_F_D_rGN(styr_D_rGN,endyr_D_rGN))/(endyr_D_rGN-styr_D_rGN+1.0);	

    if(iyear>=styr_D_cHL & iyear<=endyr_D_cHL)
		{F_cHL_D_out(iyear)=mfexp(log_avg_F_D_cHL+log_dev_F_D_cHL(iyear));}    
    if(iyear>endyr_selex_phase1 & iyear<styr_D_cHL)
		{F_cHL_D_out(iyear)=mfexp(log_avg_F_D_cHL+log_F_avgdev_cHL_D);} 	
    F_cHL_D(iyear)=sel_cHL_D(iyear)*F_cHL_D_out(iyear);
    Fsum(iyear)+=F_cHL_D_out(iyear);
    
    if(iyear>=styr_D_rHB & iyear<=endyr_D_rHB)
		{F_rHB_D_out(iyear)=mfexp(log_avg_F_D_rHB+log_dev_F_D_rHB(iyear));}    
    if(iyear>endyr_selex_phase1 & iyear<styr_D_rHB)
      {F_rHB_D_out(iyear)=mfexp(log_avg_F_D_rHB+log_F_avgdev_rHB_D);} 	
    F_rHB_D(iyear)=sel_rHB_D(iyear)*F_rHB_D_out(iyear);
    Fsum(iyear)+=F_rHB_D_out(iyear);
 
    if(iyear>=styr_D_rGN & iyear<=endyr_D_rGN)
		{F_rGN_D_out(iyear)=mfexp(log_avg_F_D_rGN+log_dev_F_D_rGN(iyear));}  
	F_rGN_D(iyear)=sel_rHB_D(iyear)*F_rGN_D_out(iyear); 
	Fsum(iyear)+=F_rGN_D_out(iyear);
    
 
    //Total F at age
    F(iyear)=F_cHL(iyear);  //first in additive series (NO +=)
    F(iyear)+=F_cOT(iyear);
    F(iyear)+=F_cTW(iyear);
    F(iyear)+=F_rHB(iyear);
    F(iyear)+=F_rGN(iyear);
    F(iyear)+=F_cHL_D(iyear);
    F(iyear)+=F_rHB_D(iyear);
    F(iyear)+=F_rGN_D(iyear);
    
    Fapex(iyear)=max(F(iyear));
    Z(iyear)=M+F(iyear);
    
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
 
  B0=bpr_F0*R_virgin;   
  B0_q_DD=R_virgin*sum(elem_prod(N_bpr_F0(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages))); 
  			
  //KC F_initial=sel_cHL(styr)*mfexp(log_avg_F_L_cHL+log_F_dev_init_cHL)+
  //KC          sel_cOT(styr)*mfexp(log_avg_F_L_cOT+log_F_dev_init_cOT)+
  //KC          sel_cTW(styr)*mfexp(log_avg_F_L_cTW+log_F_dev_init_cTW)+
  //KC          sel_rHB(styr)*mfexp(log_avg_F_L_rHB+log_F_dev_init_rHB)+
  //KC          sel_rGN(styr)*mfexp(log_avg_F_L_rGN+log_F_dev_init_rGN);			
    
  Z_initial=M;    //KC+F_initial;

//Initial equilibrium age structure
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
  B_q_DD(styr)=sum(elem_prod(N(styr)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
    
//Rest of years 
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

  
FUNCTION get_landings_numbers //Baranov catch eqn
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      L_cHL_num(iyear,iage)=N(iyear,iage)*F_cHL(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_cOT_num(iyear,iage)=N(iyear,iage)*F_cOT(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_cTW_num(iyear,iage)=N(iyear,iage)*F_cTW(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);		
      L_rHB_num(iyear,iage)=N(iyear,iage)*F_rHB(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_rGN_num(iyear,iage)=N(iyear,iage)*F_rGN(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);        
    }          
        pred_cHL_L_knum(iyear)=sum(L_cHL_num(iyear))/1000.0;
	pred_cOT_L_knum(iyear)=sum(L_cOT_num(iyear))/1000.0;
	pred_cTW_L_knum(iyear)=sum(L_cTW_num(iyear))/1000.0;
        pred_rHB_L_knum(iyear)=sum(L_rHB_num(iyear))/1000.0;
        pred_rGN_L_knum(iyear)=sum(L_rGN_num(iyear))/1000.0;
  }

 
FUNCTION get_landings_wgt
  for (iyear=styr; iyear<=endyr; iyear++)
  {    
    L_cHL_klb(iyear)=elem_prod(L_cHL_num(iyear),wholewgt_cHL_klb(iyear));     //in 1000 lb whole weight
    L_cOT_klb(iyear)=elem_prod(L_cOT_num(iyear),wholewgt_cOT_klb(iyear));     //in 1000 lb whole weight
    L_cTW_klb(iyear)=elem_prod(L_cTW_num(iyear),wholewgt_cTW_klb(iyear));     //in 1000 lb whole weight
    L_rHB_klb(iyear)=elem_prod(L_rHB_num(iyear),wholewgt_rHB_klb(iyear));     //in 1000 lb whole weight
    L_rGN_klb(iyear)=elem_prod(L_rGN_num(iyear),wholewgt_rGN_klb(iyear));     //in 1000 lb whole weight
    
    pred_cHL_L_klb(iyear)=sum(L_cHL_klb(iyear));
    pred_cOT_L_klb(iyear)=sum(L_cOT_klb(iyear));
    pred_cTW_L_klb(iyear)=sum(L_cTW_klb(iyear));
    pred_rHB_L_klb(iyear)=sum(L_rHB_klb(iyear));
    pred_rGN_L_klb(iyear)=sum(L_rGN_klb(iyear));    
  }
 
FUNCTION get_dead_discards
  //dead discards at age (number fish) 

  //for (iyear=styr_D_cHL; iyear<=endyr_D_cHL; iyear++)
  for (iyear=styr; iyear<=endyr; iyear++)	  
  {
    for (iage=1; iage<=nages; iage++)
    {
      D_cHL_num(iyear,iage)=N(iyear,iage)*F_cHL_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
    }
    pred_cHL_D_knum(iyear)=sum(D_cHL_num(iyear))/1000.0;            //pred annual dead discards in 1000s (for matching data)
    pred_cHL_D_klb(iyear)=sum(elem_prod(D_cHL_num(iyear),wholewgt_cHL_D_klb(iyear))); //annual dead discards in 1000 lb whole (for output only)
  }

  //for (iyear=styr_D_rHB; iyear<=endyr_D_rHB; iyear++)
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      D_rHB_num(iyear,iage)=N(iyear,iage)*F_rHB_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
    }
    pred_rHB_D_knum(iyear)=sum(D_rHB_num(iyear))/1000.0;            //pred annual dead discards in 1000s (for matHBing data)
    pred_rHB_D_klb(iyear)=sum(elem_prod(D_rHB_num(iyear),wholewgt_rHB_D_klb(iyear))); //annual dead discards in 1000 lb whole (for output only)
  }

  //for (iyear=styr_D_rGN; iyear<=endyr_D_rGN; iyear++)
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      D_rGN_num(iyear,iage)=N(iyear,iage)*F_rGN_D(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
    }
    pred_rGN_D_knum(iyear)=sum(D_rGN_num(iyear))/1000.0;            //pred annual dead discards in 1000s (for matGRing data)
    pred_rGN_D_klb(iyear)=sum(elem_prod(D_rGN_num(iyear),wholewgt_rGN_D_klb(iyear))); //annual dead discards in 1000 lb whole (for output only)
  }
   
FUNCTION get_catchability_fcns    
 //Get rate increase if estimated, otherwise fixed above
  if (set_q_rate_phase>0.0)
  {

      for (iyear=styr_cpue_cHL; iyear<=endyr_cpue_cHL; iyear++)
      {   //KC if (iyear>styr_cpue_cHL & iyear <=2003)
          {//q_rate_fcn_cHL(iyear)=(1.0+q_rate)*q_rate_fcn_cHL(iyear-1); //compound
             q_rate_fcn_cHL(iyear)=(1.0+(iyear-styr_cpue_cHL)*q_rate)*q_rate_fcn_cHL(styr_cpue_cHL);  //linear
          }
//KC          if (iyear>2003) {q_rate_fcn_cHL(iyear)=q_rate_fcn_cHL(iyear-1);} 
      }   
      for (iyear=styr_cpue_rHB; iyear<=endyr_cpue_rHB; iyear++)
      {   //KC if (iyear>styr_cpue_rHB & iyear <=2003) 
          {//q_rate_fcn_rHB(iyear)=(1.0+q_rate)*q_rate_fcn_rHB(iyear-1); //compound
             q_rate_fcn_rHB(iyear)=(1.0+(iyear-styr_cpue_rHB)*q_rate)*q_rate_fcn_rHB(styr_cpue_rHB);  //linear
          }
//KC          if (iyear>2003) {q_rate_fcn_rHB(iyear)=q_rate_fcn_rHB(iyear-1);} 
      }   
       for (iyear=styr_cpue_rGN; iyear<=endyr_cpue_rGN; iyear++)
       {   //KC if (iyear>styr_cpue_rGN & iyear <=2003) 
           {//q_rate_fcn_rGN(iyear)=(1.0+q_rate)*q_rate_fcn_rGN(iyear-1); //compound
              q_rate_fcn_rGN(iyear)=(1.0+(iyear-styr_cpue_rGN)*q_rate)*q_rate_fcn_rGN(styr_cpue_rGN);  //linear
           }
//KC           if (iyear>2003) {q_rate_fcn_rGN(iyear)=q_rate_fcn_rGN(iyear-1);} 
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

     
FUNCTION get_indices
//---Predicted CPUEs------------------------

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

 //rGN cpue 
  // sel_rGN_cpue.initialize(); //approximate selectivity; includes discards (live + dead)
   q_rGN(styr_cpue_rGN)=mfexp(log_q_cpue_rGN); 
   for (iyear=styr_cpue_rGN; iyear<=endyr_cpue_rGN; iyear++)
   {   
      // sel_rGN_cpue(iyear)=F_rGN(iyear)+F_rGN_D(iyear)/Dmort_rGN;
      // sel_rGN_cpue(iyear)/=max(sel_rGN_cpue(iyear));
      N_rGN(iyear)=elem_prod(N_mdyr(iyear),sel_rGN(iyear)); 
      pred_rGN_cpue(iyear)=q_rGN(iyear)*q_rate_fcn_rGN(iyear)*q_DD_fcn(iyear)*sum(N_rGN(iyear));
      if (iyear<endyr_cpue_rGN){q_rGN(iyear+1)=q_rGN(iyear)*mfexp(q_RW_log_dev_rGN(iyear));}
   } 
 
 //MARMAP Snapper trap cpue
  q_sFT(styr_cpue_sFT)=mfexp(log_q_cpue_sFT); 
  for (iyear=styr_cpue_sFT; iyear<=endyr_cpue_sFT; iyear++)
  {   
    N_sFT(iyear)=elem_prod(N_mdyr(iyear),sel_sFT(iyear)); 
	pred_sFT_cpue(iyear)=q_sFT(iyear)*q_DD_fcn(iyear)*sum(N_sFT(iyear));
    if (iyear<endyr_cpue_sFT){q_sFT(iyear+1)=q_sFT(iyear)*mfexp(q_RW_log_dev_sFT(iyear));}
  }
  
//SEFIS sCT cpue
  q_sCT(styr_cpue_sCT)=mfexp(log_q_cpue_sCT); 
  for (iyear=styr_cpue_sCT; iyear<=endyr_cpue_sCT; iyear++)
  {   
    N_sCT(iyear)=elem_prod(N_mdyr(iyear),sel_sCT(iyear)); 
	pred_sCT_cpue(iyear)=q_sCT(iyear)*q_DD_fcn(iyear)*sum(N_sCT(iyear));
    if (iyear<endyr_cpue_sCT){q_sCT(iyear+1)=q_sCT(iyear)*mfexp(q_RW_log_dev_sCT(iyear));}
  } 

FUNCTION get_length_comps

//  //cGN lines 
  for (iyear=1;iyear<=nyr_lenc_cHL;iyear++)
  {pred_cHL_lenc(iyear)=(L_cHL_num(yrs_lenc_cHL(iyear))*lenprob_cHL)/sum(L_cHL_num(yrs_lenc_cHL(iyear)));}
 
  //cGN other 
  for (iyear=1;iyear<=nyr_lenc_cOT;iyear++)
  {pred_cOT_lenc(iyear)=(L_cOT_num(yrs_lenc_cOT(iyear))*lenprob_cOT)/sum(L_cOT_num(yrs_lenc_cOT(iyear)));}

  //headboat    
 // for (iyear=1;iyear<=nyr_rHB_lenc;iyear++) 
 // {pred_rHB_lenc(iyear)=(L_rHB_num(yrs_rHB_lenc(iyear))*lenprob_rHB)/sum(L_rHB_num(yrs_rHB_lenc(iyear)));}
  
  //headboat discards   
  for (iyear=1;iyear<=nyr_lenc_rHB_D;iyear++) 
  {pred_rHB_D_lenc(iyear)=(D_rHB_num(yrs_lenc_rHB_D(iyear))*lenprob_rHB_D)/sum(D_rHB_num(yrs_lenc_rHB_D(iyear)));}

  //commercial handline discards   
  for (iyear=1;iyear<=nyr_lenc_cHL_D;iyear++) 
  {pred_cHL_D_lenc(iyear)=(D_cHL_num(yrs_lenc_cHL_D(iyear))*lenprob_cHL_D)/sum(D_cHL_num(yrs_lenc_cHL_D(iyear)));}

 //general rec   
//  for (iyear=1;iyear<=nyr_rGN_lenc;iyear++) 
//  {pred_rGN_lenc(iyear)=(L_rGN_num(yrs_rGN_lenc(iyear))*lenprob_rGN)/sum(L_rGN_num(yrs_rGN_lenc(iyear)));}
 
  //sCT    
  //for (iyear=1;iyear<=nyr_sCT_lenc;iyear++) 
   //{pred_sCT_lenc(iyear)=(N_sCT(yrs_sCT_lenc(iyear))*lenprob_sCT)/sum(N_sCT(yrs_sCT_lenc(iyear)));}

////////below is old method
  
//   //MARMAP fishery independent
//  for (iyear=1;iyear<=nyr_lenc_sFT;iyear++)
//    {pred_sFT_lenc(iyear)=(N_sFT(yrs_lenc_sFT(iyear))*lenprob_sFT)/sum(N_sFT(yrs_lenc_sFT(iyear)));}
  
// //Commercial
//   //lines
//  for (iyear=1;iyear<=nyr_lenc_cHL;iyear++)
//   
//   {if (yrs_lenc_cHL(iyear)<=endyr_selex_phase1)
//		{pred_cHL_lenc(iyear)=(L_cHL_num(yrs_lenc_cHL(iyear))*lenprob_cHL1)/sum(L_cHL_num(yrs_lenc_cHL(iyear)));}
	
//	if (yrs_lenc_cHL(iyear)>endyr_selex_phase1)
//		{pred_cHL_lenc(iyear)=(L_cHL_num(yrs_lenc_cHL(iyear))*lenprob_cHL2)/sum(L_cHL_num(yrs_lenc_cHL(iyear)));}
//	}
//  //Other 
//  for (iyear=1;iyear<=nyr_lenc_cOT;iyear++)
//  {
//     if (yrs_lenc_cOT(iyear)<=endyr_selex_phase1)
//     {pred_cOT_lenc(iyear)=(L_cOT_num(yrs_lenc_cOT(iyear))*lenprob_cOT1)/sum(L_cOT_num(yrs_lenc_cOT(iyear)));}
     
//	 else
//     {pred_cOT_lenc(iyear)=(L_cOT_num(yrs_lenc_cOT(iyear))*lenprob_cOT2)/sum(L_cOT_num(yrs_lenc_cOT(iyear))); }      
//  }
//  //Discards 
// for (iyear=1;iyear<=nyr_lenc_cHL_D;iyear++)
// {
//    if (yrs_lenc_cHL_D(iyear)<=endyr_2008)
//    {pred_cHL_D_lenc(iyear)=(D_cHL_num(yrs_lenc_cHL_D(iyear))*lenprob_cHL_D2)
//                            /sum(D_cHL_num(yrs_lenc_cHL_D(iyear)));}
//    else
//    {pred_cHL_D_lenc(iyear)=(D_cHL_num(yrs_lenc_cHL_D(iyear))*lenprob_cHL_D3)
//                            /sum(D_cHL_num(yrs_lenc_cHL_D(iyear)));}
// }

//Headboat
// for (iyear=1;iyear<=nyr_rHB_lenc;iyear++)
// {	if (yrs_rHB_lenc(iyear)<=endyr_selex_phase1)
//		{pred_rHB_lenc(iyear)=(L_rHB_num(yrs_rHB_lenc(iyear))*lenprob_rHB1)/sum(L_rHB_num(yrs_rHB_lenc(iyear)));}
//   if (yrs_rHB_lenc(iyear)>endyr_selex_phase1 &    		
//	yrs_rHB_lenc(iyear)<=endyr_selex_phase2)
//		{pred_rHB_lenc(iyear)=(L_rHB_num(yrs_rHB_lenc(iyear))*lenprob_rHB2)/sum(L_rHB_num(yrs_rHB_lenc(iyear)));}
//	if (yrs_rHB_lenc(iyear)>endyr_selex_phase2 &  
//	yrs_rHB_lenc(iyear)<=endyr_selex_phase3)
//		{pred_rHB_lenc(iyear)=(L_rHB_num(yrs_rHB_lenc(iyear))*lenprob_rHB3)/sum(L_rHB_num(yrs_rHB_lenc(iyear)));}
//	if (yrs_rHB_lenc(iyear)>endyr_selex_phase3) 
//		{pred_rHB_lenc(iyear)=(L_rHB_num(yrs_rHB_lenc(iyear))*lenprob_rHB4)/sum(L_rHB_num(yrs_rHB_lenc(iyear)));}  
//  }   
	//Discards
//  for (iyear=1;iyear<=nyr_lenc_rHB_D;iyear++)
//  { 
//	if (yrs_lenc_rHB_D(iyear)<=endyr_selex_phase3) 			{pred_rHB_D_lenc(iyear)=(D_rHB_num(yrs_lenc_rHB_D(iyear))*lenprob_rHB_D3)/sum(D_rHB_num(yrs_lenc_rHB_D(iyear)));}
    
//	else    		{pred_rHB_D_lenc(iyear)=(D_rHB_num(yrs_lenc_rHB_D(iyear))*lenprob_rHB_D4)/sum(D_rHB_num(yrs_lenc_rHB_D(iyear)));}
//   }
 
//rGN
// for (iyear=1;iyear<=nyr_rGN_lenc;iyear++)
//  { if (yrs_rGN_lenc(iyear)<=endyr_selex_phase1) 
//		{pred_rGN_lenc(iyear)=(L_rGN_num(yrs_rGN_lenc(iyear))*lenprob_rGN1)/sum(L_rGN_num(yrs_rGN_lenc(iyear)));}
	
//	if (yrs_rGN_lenc(iyear)>endyr_selex_phase1 &
//	yrs_rGN_lenc(iyear)<=endyr_selex_phase2)
//		{pred_rGN_lenc(iyear)=(L_rGN_num(yrs_rGN_lenc(iyear))*lenprob_rGN2)/sum(L_rGN_num(yrs_rGN_lenc(iyear)));}
 
//    if (yrs_rGN_lenc(iyear)>endyr_selex_phase2 & 
//	yrs_rGN_lenc(iyear)<=endyr_selex_phase3)
//		{pred_rGN_lenc(iyear)=(L_rGN_num(yrs_rGN_lenc(iyear))*lenprob_rGN3)/sum(L_rGN_num(yrs_rGN_lenc(iyear)));}

//    if (yrs_rGN_lenc(iyear)>endyr_selex_phase3)
//		{pred_rGN_lenc(iyear)=(L_rGN_num(yrs_rGN_lenc(iyear))*lenprob_rGN4)/sum(L_rGN_num(yrs_rGN_lenc(iyear)));}
// }
  
FUNCTION get_age_comps
   
  //Commercial handline
  for (iyear=1;iyear<=nyr_agec_cHL;iyear++) 
  {
    ErrorFree_cHL_agec(iyear)=L_cHL_num(yrs_agec_cHL(iyear))/sum(L_cHL_num(yrs_agec_cHL(iyear)));  
    pred_cHL_agec_allages(iyear)=age_error*(ErrorFree_cHL_agec(iyear)/sum(ErrorFree_cHL_agec(iyear)));   
    for (iage=1; iage<=nages_agec; iage++) {pred_cHL_agec(iyear,iage)=pred_cHL_agec_allages(iyear,iage);} 
    //for (iage=(nages_agec+1); iage<=nages; iage++) {pred_cHL_agec(iyear,nages_agec)+=pred_cHL_agec_allages(iyear,iage);} //plus group                             
  }
 
  //Headboat
 for (iyear=1;iyear<=nyr_agec_rHB;iyear++)
  {
    ErrorFree_rHB_agec(iyear)=L_rHB_num(yrs_agec_rHB(iyear))/sum(L_rHB_num(yrs_agec_rHB(iyear)));
    pred_rHB_agec_allages(iyear)=age_error*ErrorFree_rHB_agec(iyear); 
    for (iage=1; iage<=nages_agec; iage++) {pred_rHB_agec(iyear,iage)=pred_rHB_agec_allages(iyear,iage);} 
    //for (iage=(nages_agec+1); iage<=nages; iage++) {pred_rHB_agec(iyear,nages_agec_rHB)+=pred_rHB_agec_allages(iyear,iage);} //plus group                        
  }
  
  
   //Recreational 
 for (iyear=1;iyear<=nyr_agec_rGN;iyear++)
  {
    ErrorFree_rGN_agec(iyear)=L_rGN_num(yrs_agec_rGN(iyear))/sum(L_rGN_num(yrs_agec_rGN(iyear)));
    pred_rGN_agec_allages(iyear)=age_error*ErrorFree_rGN_agec(iyear); 
    for (iage=1; iage<=nages_agec; iage++) {pred_rGN_agec(iyear,iage)=pred_rGN_agec_allages(iyear,iage);} 
    //for (iage=(nages_agec+1); iage<=nages; iage++) {pred_rGN_agec(iyear,nages_agec)+=pred_rGN_agec_allages(iyear,iage);} //plus group                        
  }


   //MARMAP sCT 
 for (iyear=1;iyear<=nyr_agec_sCT;iyear++)
  {
    ErrorFree_sCT_agec(iyear)=N_sCT(yrs_agec_sCT(iyear))/sum(N_sCT(yrs_agec_sCT(iyear)));
    pred_sCT_agec_allages(iyear)=age_error*ErrorFree_sCT_agec(iyear); 
    for (iage=1; iage<=nages_agec; iage++) {pred_sCT_agec(iyear,iage)=pred_sCT_agec_allages(iyear,iage);} 
    //for (iage=(nages_agec+1); iage<=nages; iage++) {pred_sCT_agec(iyear,nages_agec)+=pred_sCT_agec_allages(iyear,iage);} //plus group                        
  }
  
////--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_weighted_current //should not and does not include historic trawl
  F_temp_sum=0.0;
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cHL+
        sum(log_dev_F_L_cHL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cOT+
        sum(log_dev_F_L_cOT((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
  //F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cTW+
  //      sum(log_dev_F_L_cTW((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);
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
  F_cOT_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cOT+
        sum(log_dev_F_L_cOT((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  //F_cTW_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cTW+
  //      //sum(log_dev_F_L_cTW((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_te//mp_sum;
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
  log_F_dev_end_cOT=sum(log_dev_F_L_cOT((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  //log_F_dev_end_cTW=sum(log_dev_F_L_cTW((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_rHB=sum(log_dev_F_L_rHB((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  
  log_F_dev_end_rGN=sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  

  log_F_dev_end_cHL_D=sum(log_dev_F_D_cHL((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_rHB_D=sum(log_dev_F_D_rHB((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_rGN_D=sum(log_dev_F_D_rGN((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;

  F_end_L=sel_cHL(endyr)*mfexp(log_avg_F_L_cHL+log_F_dev_end_cHL)+
          sel_cOT(endyr)*mfexp(log_avg_F_L_cOT+log_F_dev_end_cOT)+
          //sel_cTW(endyr)*mfexp(log_avg_F_L_cTW+log_F_dev_end_cTW)+
          sel_rHB(endyr)*mfexp(log_avg_F_L_rHB+log_F_dev_end_rHB)+
          sel_rGN(endyr)*mfexp(log_avg_F_L_rGN+log_F_dev_end_rGN);   

 //cout << F_end_L << endl;

  F_end_D=sel_cHL_D(endyr)*mfexp(log_avg_F_D_cHL+log_F_dev_end_cHL_D)+
          sel_rHB_D(endyr)*mfexp(log_avg_F_D_rHB+log_F_dev_end_rHB_D)+
          sel_rHB_D(endyr)*mfexp(log_avg_F_D_rGN+log_F_dev_end_rGN_D);   
    
  F_end=F_end_L+F_end_D;
  F_end_apex=max(F_end);
  
  sel_wgted_tot=F_end/F_end_apex;
  sel_wgted_L=elem_prod(sel_wgted_tot, elem_div(F_end_L,F_end));
  sel_wgted_D=elem_prod(sel_wgted_tot, elem_div(F_end_D,F_end));
  
  wgt_wgted_L_denom=F_cHL_prop+F_cOT_prop+F_rHB_prop+F_rGN_prop;  
  wgt_wgted_L_klb=F_cHL_prop/wgt_wgted_L_denom*wholewgt_cHL_klb(endyr)+ 
		  F_cOT_prop/wgt_wgted_L_denom*wholewgt_cOT_klb(endyr)+
                  F_rHB_prop/wgt_wgted_L_denom*wholewgt_rHB_klb(endyr)+
                  F_rGN_prop/wgt_wgted_L_denom*wholewgt_rGN_klb(endyr);                          

  wgt_wgted_D_denom=F_cHL_D_prop+F_rHB_D_prop+F_rGN_D_prop;  
  wgt_wgted_D_klb=F_cHL_D_prop/wgt_wgted_D_denom*wholewgt_cHL_D_klb(endyr)+ 
                  F_rHB_D_prop/wgt_wgted_D_denom*wholewgt_rHB_D_klb(endyr)+
                  F_rGN_D_prop/wgt_wgted_D_denom*wholewgt_rGN_D_klb(endyr);                
  
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
  SSB_F30_out=sum(elem_prod(N_age_spr_spawn,reprod));
  B_F30_out=sum(elem_prod(N_age_spr,wgt_mt));
  L_F30_klb_out=sum(elem_prod(L_age_F30,wgt_wgted_L_klb)); //in whole weight
  L_F30_knum_out=sum(L_age_F30)/1000.0;  
  D_F30_klb_out=sum(elem_prod(D_age_F30,wgt_wgted_D_klb)); //in whole weight   
  D_F30_knum_out=sum(D_age_F30)/1000.0;    

	
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_miscellaneous_stuff

//switch here if var_rec_dev <=dzero 
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
        L_total_klb_yr(iyear)=pred_cHL_L_klb(iyear)+pred_cOT_L_klb(iyear)+pred_cTW_L_klb(iyear)+pred_rHB_L_klb(iyear)+pred_rGN_L_klb(iyear);
        L_total_knum_yr(iyear)=pred_cHL_L_knum(iyear)+pred_cOT_L_knum(iyear)+pred_cTW_L_knum(iyear)+pred_rHB_L_knum(iyear)+pred_rGN_L_knum(iyear);
                
        B(iyear)=elem_prod(N(iyear),wgt_mt);
        totN(iyear)=sum(N(iyear));
        totB(iyear)=sum(B(iyear));   
        
        //if (iyear>=styr_D_cHL && iyear<=endyr_D_cHL)
	    if (iyear>endyr_selex_phase1 && iyear<=endyr)			
        {
         D_total_knum_yr(iyear)+=pred_cHL_D_knum(iyear);
         D_total_klb_yr(iyear)+=pred_cHL_D_klb(iyear);
         D_cHL_klb(iyear)=elem_prod(D_cHL_num(iyear),wholewgt_cHL_D_klb(iyear));     //in 1000 lb 
        }
                
        //if (iyear>=styr_D_rHB && iyear<=endyr_D_rHB)
		if (iyear>endyr_selex_phase1 && iyear<=endyr)
        {
         D_total_knum_yr(iyear)+=pred_rHB_D_knum(iyear);
         D_total_klb_yr(iyear)+=pred_rHB_D_klb(iyear);
         D_rHB_klb(iyear)=elem_prod(D_rHB_num(iyear),wholewgt_rHB_D_klb(iyear));     //in 1000 lb 
        }    
    	
        //if (iyear>=styr_D_rGN && iyear<=endyr_D_rGN)
		if (iyear>=styr_D_rGN && iyear<=endyr)	
        {
         D_total_knum_yr(iyear)+=pred_rGN_D_knum(iyear);
         D_total_klb_yr(iyear)+=pred_rGN_D_klb(iyear);
         D_rGN_klb(iyear)=elem_prod(D_rGN_num(iyear),wholewgt_rGN_D_klb(iyear));     //in 1000 lb 
        }                          
  }
  
  L_total_num=L_cHL_num+L_cOT_num+L_rHB_num+L_cTW_num+L_rGN_num;   //landings at age in number fish
  L_total_klb=L_cHL_klb+L_cOT_klb+L_cTW_klb+L_rHB_klb+L_rGN_klb;   //landings at age in klb whole weight
 
  D_total_num=(D_cHL_num+D_rHB_num+D_rGN_num);          //discards at age in number fish
  D_total_klb=D_cHL_klb+D_rHB_klb+D_rGN_klb;            //discards at age in klb whole weight
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
		FD_age_proj=sel_wgted_D*F_proj(iyear);
		
        Z_proj(iyear)=M+FL_age_proj+FD_age_proj;
        N_spawn_proj(iyear)(1,nages)=elem_prod(N_proj(iyear)(1,nages),(mfexp(-1.*(Z_proj(iyear)(1,nages))*spawn_time_frac))); //peak spawning time
		SSB_proj(iyear)= sum(elem_prod(N_spawn_proj(iyear),reprod));
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

  f_cHL_cpue=0.0;
  f_cHL_cpue=lk_lognormal(pred_cHL_cpue, obs_cpue_cHL, obs_cv_cpue_cHL, w_cpue_cHL);
  fval+=f_cHL_cpue;
  fval_data+=f_cHL_cpue;  

  f_rHB_cpue=0.0;
  f_rHB_cpue=lk_lognormal(pred_rHB_cpue, obs_cpue_rHB, obs_cv_cpue_rHB, w_cpue_rHB);
  fval+=f_rHB_cpue;
  fval_data+=f_rHB_cpue;  

  f_rGN_cpue=0.0;
  f_rGN_cpue=lk_lognormal(pred_rGN_cpue, obs_cpue_rGN, obs_cv_cpue_rGN, w_cpue_rGN);
  fval+=f_rGN_cpue;
  fval_data+=f_rGN_cpue;  

  f_sFT_cpue=0.0;
  f_sFT_cpue=lk_lognormal(pred_sFT_cpue, obs_cpue_sFT, obs_cv_cpue_sFT, w_cpue_sFT);
  fval+=f_sFT_cpue;
  fval_data+=f_sFT_cpue;
  
  f_sCT_cpue=0.0;
  f_sCT_cpue=lk_lognormal(pred_sCT_cpue, obs_cpue_sCT, obs_cv_cpue_sCT, w_cpue_sCT);
  fval+=f_sCT_cpue;
  fval_data+=f_sCT_cpue;
  
//---Landings-------------------------------
  
  //f_cHL_L in 1000 lb whole wgt  
    f_cHL_L=lk_lognormal(pred_cHL_L_klb(styr_L_cHL,endyr_L_cHL), obs_L_cHL(styr_L_cHL,endyr_L_cHL),
                      obs_cv_L_cHL(styr_L_cHL,endyr_L_cHL), w_L);
    fval+=f_cHL_L;
    fval_data+=f_cHL_L;

  //f_cOT_L in 1000 lb whole wgt  
    f_cOT_L=lk_lognormal(pred_cOT_L_klb(styr_L_cOT,endyr_L_cOT), obs_L_cOT(styr_L_cOT,endyr_L_cOT),
                      obs_cv_L_cOT(styr_L_cOT,endyr_L_cOT), w_L);
                      
    fval+=f_cOT_L;
    fval_data+=f_cOT_L;

  //f_ccHTR_L in 1000 lb whole wgt  
  f_cTW_L=lk_lognormal(pred_cTW_L_klb(styr_L_cTW,endyr_L_cTW), obs_L_cTW(styr_L_cTW,endyr_L_cTW),
                      obs_cv_L_cTW(styr_L_cTW,endyr_L_cTW), w_L);
                      
  fval+=f_cTW_L;
  fval_data+=f_cTW_L;
  
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

//---Discards-------------------------------
  
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
  
//---Length comps-------------------------------

  //f_cHL_lenc
  //f_cHL_lenc=lk_multinomial(nsamp_lenc_cHL, pred_cHL_lenc, obs_lenc_cHL, nyr_lenc_cHL, minSS_lenc_cHL, w_lenc_cHL);
  //KC f_cHL_lenc=lk_robust_multinomial(nsamp_lenc_cHL, pred_cHL_lenc, obs_lenc_cHL, nyr_lenc_cHL, double(nlenbins), minSS_lenc_cHL, w_lenc_cHL);
  //f_cHL_lenc=lk_logistic_normal(nsamp_lenc_cHL, pred_cHL_lenc, obs_lenc_cHL, nyr_lenc_cHL, double(nlenbins), minSS_lenc_cHL);
  f_cHL_lenc=lk_dirichlet_multinomial(nsamp_lenc_cHL, pred_cHL_lenc, obs_lenc_cHL, nyr_lenc_cHL, double(nlenbins), minSS_lenc_cHL, log_dm_lenc_cHL);
  fval+=f_cHL_lenc;
  fval_data+=f_cHL_lenc;

  //f_cOT_lenc
  f_cOT_lenc=lk_multinomial(nsamp_lenc_cOT, pred_cOT_lenc, obs_lenc_cOT, nyr_lenc_cOT, minSS_lenc_cOT, w_lenc_cOT);
  //KC f_cOT_lenc=lk_robust_multinomial(nsamp_lenc_cOT, pred_cOT_lenc, obs_lenc_cOT, nyr_lenc_cOT, double(nlenbins), minSS_lenc_cOT, w_lenc_cOT);
  //f_cOT_lenc=lk_logistic_normal(nsamp_lenc_cOT, pred_cOT_lenc, obs_lenc_cOT, nyr_lenc_cOT, double(nlenbins), minSS_lenc_cOT);
  //f_cOT_lenc=lk_dirichlet_multinomial(nsamp_lenc_cOT, pred_cOT_lenc, obs_lenc_cOT, nyr_lenc_cOT, double(nlenbins), minSS_lenc_cOT, log_dm_cOT_lc);
  fval+=f_cOT_lenc;
  fval_data+=f_cOT_lenc;

  ////f_rHB_lenc
  ////f_rHB_lenc=lk_multinomial(nsamp_rHB_lenc, pred_rHB_lenc, obs_rHB_lenc, nyr_rHB_lenc, minSS_rHB_lenc, w_lc_rHB);
  ////KC f_rHB_lenc=lk_robust_multinomial(nsamp_rHB_lenc, pred_rHB_lenc, obs_rHB_lenc, nyr_rHB_lenc, double(nlenbins), minSS_rHB_lenc, w_lc_rHB);
  ////f_rHB_lenc=lk_logistic_normal(nsamp_rHB_lenc, pred_rHB_lenc, obs_rHB_lenc, nyr_rHB_lenc, double(nlenbins), minSS_rHB_lenc);
  // f_rHB_lenc=lk_dirichlet_multinomial(nsamp_rHB_lenc, pred_rHB_lenc, obs_rHB_lenc, nyr_rHB_lenc, double(nlenbins), minSS_rHB_lenc, log_dm_rHB_lc);
  // fval+=f_rHB_lenc;
  //fval_data+=f_rHB_lenc;
  
 // //f_rGN_lenc
 // //f_rGN_lenc=lk_multinomial(nsamp_rGN_lenc, pred_rGN_lenc, obs_rGN_lenc, nyr_rGN_lenc, minSS_rGN_lenc, w_lc_rGN);
 // //KC f_rGN_lenc=lk_robust_multinomial(nsamp_rGN_lenc, pred_rGN_lenc, obs_rGN_lenc, nyr_rGN_lenc, double(nlenbins), minSS_rGN_lenc, w_lc_rGN);
 // //f_rGN_lenc=lk_logistic_normal(nsamp_rGN_lenc, pred_rGN_lenc, obs_rGN_lenc, nyr_rGN_lenc, double(nlenbins), minSS_rGN_lenc);
 // f_rGN_lenc=lk_dirichlet_multinomial(nsamp_rGN_lenc, pred_rGN_lenc, obs_rGN_lenc, nyr_rGN_lenc, double(nlenbins), minSS_rGN_lenc, log_dm_rGN_lc);
 // fval+=f_rGN_lenc;
 // fval_data+=f_rGN_lenc;
  
  //f_rHB_D_lenc
  //f_rHB_D_lenc=lk_multinomial(nsamp_lenc_rHB_D, pred_rHB_D_lenc, obs_lenc_rHB_D, nyr_lenc_rHB_D, minSS_lenc_rHB_D, w_lenc_rHB_D);
  //KC f_rHB_D_lenc=lk_robust_multinomial(nsamp_lenc_rHB_D, pred_rHB_D_lenc, obs_lenc_rHB_D, nyr_lenc_rHB_D, double(nlenbins), minSS_lenc_rHB_D, w_lenc_rHB_D);
  //f_rHB_D_lenc=lk_logistic_normal(nsamp_lenc_rHB_D, pred_rHB_D_lenc, obs_lenc_rHB_D, nyr_lenc_rHB_D, double(nlenbins), minSS_lenc_rHB_D);
  f_rHB_D_lenc=lk_dirichlet_multinomial(nsamp_lenc_rHB_D, pred_rHB_D_lenc, obs_lenc_rHB_D, nyr_lenc_rHB_D, double(nlenbins), minSS_lenc_rHB_D, log_dm_lenc_rHB_D);
  fval+=f_rHB_D_lenc;
  fval_data+=f_rHB_D_lenc;

  //f_cHL_D_lenc
  //f_cHL_D_lenc=lk_multinomial(nsamp_lenc_cHL_D, pred_cHL_D_lenc, obs_lenc_cHL_D, nyr_lenc_cHL_D, minSS_lenc_cHL_D, w_lenc_cHL_D);
  //KC f_cHL_D_lenc=lk_robust_multinomial(nsamp_lenc_cHL_D, pred_cHL_D_lenc, obs_lenc_cHL_D, nyr_lenc_cHL_D, double(nlenbins), minSS_lenc_cHL_D, w_lenc_cHL_D);
  f_cHL_D_lenc=lk_dirichlet_multinomial(nsamp_lenc_cHL_D, pred_cHL_D_lenc, obs_lenc_cHL_D, nyr_lenc_cHL_D, double(nlenbins), minSS_lenc_cHL_D, log_dm_lenc_cHL_D);
  fval+=f_cHL_D_lenc;
  fval_data+=f_cHL_D_lenc;

  //f_sFT_lenc
  f_sFT_lenc=lk_multinomial(nsamp_lenc_sFT, pred_sFT_lenc, obs_lenc_sFT, nyr_lenc_sFT, minSS_lenc_sFT, w_lenc_sFT);
  //KC f_sFT_lenc=lk_robust_multinomial(nsamp_lenc_sFT, pred_sFT_lenc, obs_lenc_sFT, nyr_lenc_sFT, double(nlenbins), minSS_lenc_sFT, w_lenc_sFT);
  //f_sFT_lenc=lk_logistic_normal(nsamp_lenc_sFT, pred_sFT_lenc, obs_lenc_sFT, nyr_lenc_sFT, double(nlenbins), minSS_lenc_sFT);
  //f_sFT_lenc=lk_dirichlet_multinomial(nsamp_lenc_sFT, pred_sFT_lenc, obs_lenc_sFT, nyr_lenc_sFT, double(nlenbins), minSS_lenc_sFT, log_dm_sFT_lc);
  fval+=f_sFT_lenc;
  fval_data+=f_sFT_lenc;
  
  //f_sCT_lenc
  //f_sCT_lenc=lk_robust_multinomial(nsamp_sCT_lenc, pred_sCT_lenc, obs_sCT_lenc, nyr_sCT_lenc, double(nlenbins), minSS_sCT_lenc, w_lc_sCT);
  //f_sCT_lenc=lk_logistic_normal(nsamp_sCT_lenc, pred_sCT_lenc, obs_sCT_lenc, nyr_sCT_lenc, double(nlenbins), minSS_sCT_lenc);
  //f_sCT_lenc=lk_dirichlet_multinomial(nsamp_sCT_lenc, pred_sCT_lenc, obs_sCT_lenc, nyr_sCT_lenc, double(nlenbins), minSS_sCT_lenc, log_dm_sCT_lc);
  //fval+=f_sCT_lenc;
  //fval_data+=f_sCT_lenc;
   
//---Age comps-------------------------------

  //f_cHL_agec
  //f_cHL_agec=lk_multinomial(nsamp_agec_cHL, pred_cHL_agec, obs_agec_cHL, nyr_agec_cHL, minSS_agec_cHL, w_agec_cHL);
  //KC f_cHL_agec=lk_robust_multinomial(nsamp_agec_cHL, pred_cHL_agec, obs_agec_cHL, nyr_agec_cHL, double(nages_agec), minSS_agec_cHL, w_agec_cHL);
  //f_cHL_agec=lk_logistic_normal(nsamp_agec_cHL, pred_cHL_agec, obs_agec_cHL, nyr_agec_cHL, double(nages_agec), minSS_agec_cHL);
  f_cHL_agec=lk_dirichlet_multinomial(nsamp_agec_cHL, pred_cHL_agec, obs_agec_cHL, nyr_agec_cHL, double(nages_agec), minSS_agec_cHL, log_dm_agec_cHL);
  fval+=f_cHL_agec;
  fval_data+=f_cHL_agec;

  //f_rHB_agec
  //f_rHB_agec=lk_multinomial(nsamp_agec_rHB, pred_rHB_agec, obs_agec_rHB, nyr_agec_rHB, minSS_agec_rHB, w_agec_rHB);
  //KC f_rHB_agec=lk_robust_multinomial(nsamp_agec_rHB, pred_rHB_agec, obs_agec_rHB, nyr_agec_rHB, double(nages_agec), minSS_agec_rHB, w_agec_rHB);
  //f_rHB_agec=lk_logistic_normal(nsamp_agec_rHB, pred_rHB_agec, obs_agec_rHB, nyr_agec_rHB, double(nages_agec), minSS_agec_rHB);
  f_rHB_agec=lk_dirichlet_multinomial(nsamp_agec_rHB, pred_rHB_agec, obs_agec_rHB, nyr_agec_rHB, double(nages_agec), minSS_agec_rHB, log_dm_agec_rHB);
  fval+=f_rHB_agec;
  fval_data+=f_rHB_agec;

  //f_rGN_agec
  //f_rGN_agec=lk_multinomial(nsamp_agec_rGN, pred_rGN_agec, obs_agec_rGN, nyr_agec_rGN, minSS_agec_rGN, w_agec_rGN);
  //KC f_rGN_agec=lk_robust_multinomial(nsamp_agec_rGN, pred_rGN_agec, obs_agec_rGN, nyr_agec_rGN, double(nages_agec), minSS_agec_rGN, w_agec_rGN);
  //f_rGN_agec=lk_logistic_normal(nsamp_agec_rGN, pred_rGN_agec, obs_agec_rGN, nyr_agec_rGN, double(nages_agec), minSS_agec_rGN);
  f_rGN_agec=lk_dirichlet_multinomial(nsamp_agec_rGN, pred_rGN_agec, obs_agec_rGN, nyr_agec_rGN, double(nages_agec), minSS_agec_rGN, log_dm_agec_rGN);
  fval+=f_rGN_agec;
  fval_data+=f_rGN_agec;
  
  //f_sCT_agec
  //f_sCT_agec=lk_multinomial(nsamp_agec_sCT, pred_sCT_agec, obs_agec_sCT, nyr_agec_sCT, minSS_agec_sCT, w_agec_sCT);
  //KC f_sCT_agec=lk_robust_multinomial(nsamp_agec_sCT, pred_sCT_agec, obs_agec_sCT, nyr_agec_sCT, double(nages_agec), minSS_agec_sCT, w_agec_sCT);
  //f_sCT_agec=lk_logistic_normal(nsamp_agec_sCT, pred_sCT_agec, obs_agec_sCT, nyr_agec_sCT, double(nages_agec), minSS_agec_sCT);
  f_sCT_agec=lk_dirichlet_multinomial(nsamp_agec_sCT, pred_sCT_agec, obs_agec_sCT, nyr_agec_sCT, double(nages_agec), minSS_agec_sCT, log_dm_agec_sCT);
  fval+=f_sCT_agec;
  fval_data+=f_sCT_agec;
  


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
  
 //Random walk components of fishery dependent indices
 f_rHB_RWq_cpue=0.0;
 for (iyear=styr_cpue_rHB; iyear<endyr_cpue_rHB; iyear++)
     {f_rHB_RWq_cpue+=square(q_RW_log_dev_rHB(iyear))/(2.0*set_RWq_var);}
 fval+=f_rHB_RWq_cpue;   

 f_cHL_RWq_cpue=0.0;
 for (iyear=styr_cpue_cHL; iyear<endyr_cpue_cHL; iyear++)
     {f_cHL_RWq_cpue+=square(q_RW_log_dev_cHL(iyear))/(2.0*set_RWq_var);}
 fval+=f_cHL_RWq_cpue;

//kc added below
 f_rGN_RWq_cpue=0.0;
 for (iyear=styr_cpue_rGN; iyear<endyr_cpue_rGN; iyear++)
     {f_rGN_RWq_cpue+=square(q_RW_log_dev_rGN(iyear))/(2.0*set_RWq_var);}
 fval+=f_rGN_RWq_cpue; 
  
//---Priors---------------------------------------------------
//neg_log_prior arguments: estimate, prior mean, prior var/-CV, pdf type
//Variance input as a negative value is considered to be CV in arithmetic space (CV=-1 implies loose prior) 
//pdf type 1=none, 2=lognormal, 3=normal, 4=beta 
  f_priors=0.0; 
  f_priors+=neg_log_prior(len_cv_val,set_len_cv(5),set_len_cv(6),set_len_cv(7));
   
  f_priors+=neg_log_prior(steep,set_steep(5),set_log_R0(6),set_log_R0(7)); 
  f_priors+=neg_log_prior(log_R0,set_log_R0(5),set_log_R0(6),set_log_R0(7)); 
  f_priors+=neg_log_prior(R_autocorr,set_R_autocorr(5),set_R_autocorr(6),set_R_autocorr(7));
  f_priors+=neg_log_prior(rec_sigma,set_rec_sigma(5),set_rec_sigma(6),set_rec_sigma(7));
 
  f_priors+=neg_log_prior(selpar_A50_cHL1,set_selpar_A50_cHL1(5), set_selpar_A50_cHL1(6), set_selpar_A50_cHL1(7));
  f_priors+=neg_log_prior(selpar_slope_cHL1,set_selpar_slope_cHL1(5), set_selpar_slope_cHL1(6), set_selpar_slope_cHL1(7));
  f_priors+=neg_log_prior(selpar_A50_cHL2,set_selpar_A50_cHL2(5), set_selpar_A50_cHL2(6), set_selpar_A50_cHL2(7));
  f_priors+=neg_log_prior(selpar_slope_cHL2,set_selpar_slope_cHL2(5), set_selpar_slope_cHL2(6), set_selpar_slope_cHL2(7));
  f_priors+=neg_log_prior(selpar_A50_cHL1,set_selpar_A50_cHL1(5), set_selpar_A50_cHL1(6), set_selpar_A50_cHL1(7));
  f_priors+=neg_log_prior(selpar_slope_cHL1,set_selpar_slope_cHL1(5), set_selpar_slope_cHL1(6), set_selpar_slope_cHL1(7));

  f_priors+=neg_log_prior(selpar_A50_cOT1,set_selpar_A50_cOT1(5), set_selpar_A50_cOT1(6), set_selpar_A50_cOT1(7));
  f_priors+=neg_log_prior(selpar_slope_cOT1,set_selpar_slope_cOT1(5), set_selpar_slope_cOT1(6), set_selpar_slope_cOT1(7));
  f_priors+=neg_log_prior(selpar_A502_cOT1,set_selpar_A502_cOT1(5), set_selpar_A502_cOT1(6), set_selpar_A502_cOT1(7));
  f_priors+=neg_log_prior(selpar_slope2_cOT1,set_selpar_slope2_cOT1(5), set_selpar_slope2_cOT1(6), set_selpar_slope2_cOT1(7));

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
  f_priors+=neg_log_prior(selpar_A50_rHB4,set_selpar_A50_rHB4(5), set_selpar_A50_rHB4(6), set_selpar_A50_rHB4(7));
  f_priors+=neg_log_prior(selpar_slope_rHB4,set_selpar_slope_rHB4(5), set_selpar_slope_rHB4(6), set_selpar_slope_rHB4(7));
  f_priors+=neg_log_prior(selpar_A502_rHB4,set_selpar_A502_rHB4(5), set_selpar_A502_rHB4(6), set_selpar_A502_rHB4(7));
  f_priors+=neg_log_prior(selpar_slope2_rHB4,set_selpar_slope2_rHB4(5), set_selpar_slope2_rHB4(6), set_selpar_slope2_rHB4(7));

  f_priors+=neg_log_prior(selpar_A50_rGN2,set_selpar_A50_rGN2(5), set_selpar_A50_rGN2(6), set_selpar_A50_rGN2(7));
  f_priors+=neg_log_prior(selpar_slope_rGN2,set_selpar_slope_rGN2(5), set_selpar_slope_rGN2(6), set_selpar_slope_rGN2(7));
  f_priors+=neg_log_prior(selpar_A50_rGN3,set_selpar_A50_rGN3(5), set_selpar_A50_rGN3(6), set_selpar_A50_rGN3(7));
  f_priors+=neg_log_prior(selpar_slope_rGN3,set_selpar_slope_rGN3(5), set_selpar_slope_rGN3(6), set_selpar_slope_rGN3(7));
  f_priors+=neg_log_prior(selpar_A502_rGN3,set_selpar_A502_rGN3(5), set_selpar_A502_rGN3(6), set_selpar_A502_rGN3(7));
  f_priors+=neg_log_prior(selpar_slope2_rGN3,set_selpar_slope2_rGN3(5), set_selpar_slope2_rGN3(6), set_selpar_slope2_rGN3(7));
  f_priors+=neg_log_prior(selpar_A50_rGN4,set_selpar_A50_rGN4(5), set_selpar_A50_rGN4(6), set_selpar_A50_rGN4(7));
  f_priors+=neg_log_prior(selpar_slope_rGN4,set_selpar_slope_rGN4(5), set_selpar_slope_rGN4(6), set_selpar_slope_rGN4(7));
  f_priors+=neg_log_prior(selpar_A502_rGN4,set_selpar_A502_rGN4(5), set_selpar_A502_rGN4(6), set_selpar_A502_rGN4(7));
  f_priors+=neg_log_prior(selpar_slope2_rGN4,set_selpar_slope2_rGN4(5), set_selpar_slope2_rGN4(6), set_selpar_slope2_rGN4(7));

  f_priors+=neg_log_prior(selpar_A50_sCT,set_selpar_A50_sCT(5), set_selpar_A50_sCT(6), set_selpar_A50_sCT(7));
  f_priors+=neg_log_prior(selpar_slope_sCT,set_selpar_slope_sCT(5), set_selpar_slope_sCT(6), set_selpar_slope_sCT(7));
  f_priors+=neg_log_prior(selpar_A502_sCT,set_selpar_A502_sCT(5), set_selpar_A502_sCT(6), set_selpar_A502_sCT(7));
  f_priors+=neg_log_prior(selpar_slope2_sCT,set_selpar_slope2_sCT(5), set_selpar_slope2_sCT(6), set_selpar_slope2_sCT(7));
  //f_priors+=neg_log_prior(selpar_logit_Age1_sCT, set_selpar_logit_Age1_sCT(5),set_selpar_logit_Age1_sCT(6), set_selpar_logit_Age1_sCT(7));    
  //f_priors+=neg_log_prior(selpar_logit_Age2_sCT, set_selpar_logit_Age2_sCT(5),set_selpar_logit_Age2_sCT(6),set_selpar_logit_Age2_sCT(7)); 

  f_priors+=neg_log_prior(selpar_logit_Age2_cHL_D2, set_selpar_logit_Age2_cHL_D2(5),set_selpar_logit_Age2_cHL_D2(6),set_selpar_logit_Age2_cHL_D2(7));
  f_priors+=neg_log_prior(selpar_A502_cHL_D3, set_selpar_A502_cHL_D3(5),set_selpar_A502_cHL_D3(6),set_selpar_A502_cHL_D3(7));

  f_priors+=neg_log_prior(selpar_logit_Age2_rHB_D3, set_selpar_logit_Age2_rHB_D3(5),set_selpar_logit_Age2_rHB_D3(6),set_selpar_logit_Age2_rHB_D3(7));
  f_priors+=neg_log_prior(selpar_logit_Age2_rHB_D4, set_selpar_logit_Age2_rHB_D4(5),set_selpar_logit_Age2_rHB_D4(6),set_selpar_logit_Age2_rHB_D4(7));
  //f_priors+=neg_log_prior(selpar_age1logit_D,set_selpar_age1logit_D(5), set_selpar_age1logit_D(6), set_selpar_age1logit_D(7));

  f_priors+=neg_log_prior(log_q_cpue_cHL,set_log_q_cpue_cHL(5),set_log_q_cpue_cHL(6),set_log_q_cpue_cHL(7));
  f_priors+=neg_log_prior(log_q_cpue_rHB,set_log_q_cpue_rHB(5),set_log_q_cpue_rHB(6),set_log_q_cpue_rHB(7));
  f_priors+=neg_log_prior(log_q_cpue_rGN,set_log_q_cpue_rGN(5),set_log_q_cpue_rGN(6),set_log_q_cpue_rGN(7));
  f_priors+=neg_log_prior(log_q_cpue_sFT,set_log_q_cpue_sFT(5),set_log_q_cpue_sFT(6),set_log_q_cpue_sFT(7));
  f_priors+=neg_log_prior(log_q_cpue_sCT,set_log_q_cpue_sCT(5),set_log_q_cpue_sCT(6),set_log_q_cpue_sCT(7));
  
  f_priors+=neg_log_prior(log_dm_lenc_cHL,set_log_dm_lenc_cHL(5),set_log_dm_lenc_cHL(6),set_log_dm_lenc_cHL(7));
  //f_priors+=neg_log_prior(log_dm_cOT_lc,set_log_dm_cOT_lc(5),set_log_dm_cOT_lc(6),set_log_dm_cOT_lc(7));
  //f_priors+=neg_log_prior(log_dm_rHB_lc,set_log_dm_rHB_lc(5),set_log_dm_rHB_lc(6),set_log_dm_rHB_lc(7));
  //f_priors+=neg_log_prior(log_dm_rGN_lc,set_log_dm_rGN_lc(5),set_log_dm_rGN_lc(6),set_log_dm_rGN_lc(7));
  f_priors+=neg_log_prior(log_dm_lenc_rHB_D,set_log_dm_lenc_rHB_D(5),set_log_dm_lenc_rHB_D(6),set_log_dm_lenc_rHB_D(7));
  f_priors+=neg_log_prior(log_dm_lenc_cHL_D,set_log_dm_lenc_cHL_D(5),set_log_dm_lenc_cHL_D(6),set_log_dm_lenc_cHL_D(7));
  //f_priors+=neg_log_prior(log_dm_sFT_lc,set_log_dm_sFT_lc(5),set_log_dm_sFT_lc(6),set_log_dm_sFT_lc(7));
  ////f_priors+=neg_log_prior(log_dm_sCT_lc,set_log_dm_sCT_lc(5),set_log_dm_sCT_lc(6),set_log_dm_sCT_lc(7));
  f_priors+=neg_log_prior(log_dm_agec_cHL,set_log_dm_agec_cHL(5),set_log_dm_agec_cHL(6),set_log_dm_agec_cHL(7));
  f_priors+=neg_log_prior(log_dm_agec_rHB,set_log_dm_agec_rHB(5),set_log_dm_agec_rHB(6),set_log_dm_agec_rHB(7));
  f_priors+=neg_log_prior(log_dm_agec_rGN,set_log_dm_agec_rGN(5),set_log_dm_agec_rGN(6),set_log_dm_agec_rGN(7));
  f_priors+=neg_log_prior(log_dm_agec_sCT,set_log_dm_agec_sCT(5),set_log_dm_agec_sCT(6),set_log_dm_agec_sCT(7));
  
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
 //      cout<<"xdum = "<<xdum<<endl;
	  get_weighted_current();
       cout<<"got weighted"<<endl;
      get_msy();
       cout<<"got msy"<<endl;
      get_per_recruit_stuff();
       cout<<"got per recruit"<<endl;  
      get_miscellaneous_stuff();
      cout<<"got misc stuff"<<endl;
	  get_projection();
	   cout<<"got projection"<<endl;
	  
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
      
      //report << "TotalLikelihood " << fval << endl;
      //report << "N" << endl;
      //report << N<<endl;
      //report << "F" << endl;
      //report << F <<endl;
      //report << "prob_belowsizelim_block3" << endl;
	  //report << prob_belowsizelim_block3 <<endl;
	  //
      //sdnr_lc_cHL=sdnr_multinomial(nyr_lenc_cHL, lenbins, nsamp_lenc_cHL, pred_cHL_lenc, obs_lenc_cHL, w_lenc_cHL);  
      //sdnr_lc_cOT=sdnr_multinomial(nyr_lenc_cOT, lenbins, nsamp_lenc_cOT, pred_cOT_lenc, obs_lenc_cOT, w_lenc_cOT);  	  
      //sdnr_lc_rHB=sdnr_multinomial(nyr_rHB_lenc, lenbins, nsamp_rHB_lenc, pred_rHB_lenc, obs_rHB_lenc, w_lc_rHB); 
      //sdnr_lc_rGN=sdnr_multinomial(nyr_rGN_lenc, lenbins, nsamp_rGN_lenc, pred_rGN_lenc, obs_rGN_lenc, w_lc_rGN);
      //sdnr_lc_sFT=sdnr_multinomial(nyr_lenc_sFT, lenbins, nsamp_lenc_sFT, pred_sFT_lenc, obs_lenc_sFT, w_lenc_sFT);	  
      //sdnr_lc_rHB_D=sdnr_multinomial(nyr_lenc_rHB_D, lenbins, nsamp_lenc_rHB_D, pred_rHB_D_lenc, obs_lenc_rHB_D, w_lenc_rHB_D); 
      //sdnr_lc_cHL_D=sdnr_multinomial(nyr_lenc_cHL_D, lenbins, nsamp_lenc_cHL_D, pred_cHL_D_lenc, obs_lenc_cHL_D, w_lenc_cHL_D);
      //sdnr_lc_sCT=sdnr_multinomial(nyr_sCT_lenc, lenbins, nsamp_sCT_lenc, pred_sCT_lenc, obs_sCT_lenc, w_lc_sCT); 
      //sdnr_ac_cHL=sdnr_multinomial(nyr_agec_cHL, agebins_agec, nsamp_agec_cHL, pred_cHL_agec, obs_agec_cHL, w_agec_cHL);  
      //sdnr_ac_rHB=sdnr_multinomial(nyr_agec_rHB, agebins_agec, nsamp_agec_rHB, pred_rHB_agec, obs_agec_rHB, w_agec_rHB);  
      //sdnr_ac_rGN=sdnr_multinomial(nyr_agec_rGN, agebins_agec, nsamp_agec_rGN, pred_rGN_agec, obs_agec_rGN, w_agec_rGN);  
      //sdnr_ac_sCT=sdnr_multinomial(nyr_agec_sCT, agebins_agec, nsamp_agec_sCT, pred_sCT_agec, obs_agec_sCT, w_agec_sCT);  
      
      sdnr_I_cHL=sdnr_lognormal(pred_cHL_cpue, obs_cpue_cHL, obs_cv_cpue_cHL, w_cpue_cHL);
      sdnr_I_rHB=sdnr_lognormal(pred_rHB_cpue, obs_cpue_rHB, obs_cv_cpue_rHB, w_cpue_rHB);
      sdnr_I_rGN=sdnr_lognormal(pred_rGN_cpue, obs_cpue_rGN, obs_cv_cpue_rGN, w_cpue_rGN);
      sdnr_I_sFT=sdnr_lognormal(pred_sFT_cpue, obs_cpue_sFT, obs_cv_cpue_sFT, w_cpue_sFT);
      sdnr_I_sCT=sdnr_lognormal(pred_sCT_cpue, obs_cpue_sCT, obs_cv_cpue_sCT, w_cpue_sCT);  
      
      //#################################################################################################
      //##  Passing parameters to vector for bounds check plotting
      //################################################################################################# 
       Linf_out(8)=Linf; Linf_out(1,7)=set_Linf; 
       K_out(8)=K; K_out(1,7)=set_K;
       t0_out(8)=t0; t0_out(1,7)=set_t0;
       len_cv_val_out(8)=len_cv_val; len_cv_val_out(1,7)=set_len_cv;
	   //!!cout <<"len cv val"<<len_cv_val<<endl;	   
       log_R0_out(8)=log_R0; log_R0_out(1,7)=set_log_R0;
       M_constant_out(8)=M_constant; M_constant_out(1,7)=set_M_constant;
       steep_out(8)=steep; steep_out(1,7)=set_steep;
       rec_sigma_out(8)=rec_sigma; rec_sigma_out(1,7)=set_rec_sigma;
       R_autocorr_out(8)=R_autocorr; R_autocorr_out(1,7)=set_R_autocorr;
	   
	   log_dm_cHL_lc_out(8)=log_dm_lenc_cHL; log_dm_cHL_lc_out(1,7)=set_log_dm_lenc_cHL;
	   //log_dm_cOT_lc_out(8)=log_dm_cOT_lc; log_dm_cOT_lc_out(1,7)=set_log_dm_cOT_lc;
	   //log_dm_rHB_lc_out(8)=log_dm_rHB_lc; log_dm_rHB_lc_out(1,7)=set_log_dm_rHB_lc;
	   //log_dm_rGN_lc_out(8)=log_dm_rGN_lc; log_dm_rGN_lc_out(1,7)=set_log_dm_rGN_lc;
	   log_dm_rHB_D_lc_out(8)=log_dm_lenc_rHB_D; log_dm_rHB_D_lc_out(1,7)=set_log_dm_lenc_rHB_D;
       log_dm_cHL_D_lc_out(8)=log_dm_lenc_cHL_D; log_dm_cHL_D_lc_out(1,7)=set_log_dm_lenc_cHL_D;   //kc check this? was rHB_D (repeated)
       //log_dm_sFT_lc_out(8)=log_dm_sFT_lc; log_dm_sFT_lc_out(1,7)=set_log_dm_sFT_lc;
	   ////log_dm_sCT_lc_out(8)=log_dm_sCT_lc; log_dm_sCT_lc_out(1,7)=set_log_dm_sCT_lc;
	   log_dm_cHL_ac_out(8)=log_dm_agec_cHL; log_dm_cHL_ac_out(1,7)=set_log_dm_agec_cHL;
	   log_dm_rHB_ac_out(8)=log_dm_agec_rHB; log_dm_rHB_ac_out(1,7)=set_log_dm_agec_rHB;
	   log_dm_rGN_ac_out(8)=log_dm_agec_rGN; log_dm_rGN_ac_out(1,7)=set_log_dm_agec_rGN;
	   log_dm_sCT_ac_out(8)=log_dm_agec_sCT; log_dm_sCT_ac_out(1,7)=set_log_dm_agec_sCT;
       
       selpar_A50_cHL1_out(8)=selpar_A50_cHL1; selpar_A50_cHL1_out(1,7)=set_selpar_A50_cHL1;
       selpar_slope_cHL1_out(8)=selpar_slope_cHL1; selpar_slope_cHL1_out(1,7)=set_selpar_slope_cHL1;
       selpar_A50_cHL2_out(8)=selpar_A50_cHL2; selpar_A50_cHL2_out(1,7)=set_selpar_A50_cHL2;
       selpar_slope_cHL2_out(8)=selpar_slope_cHL2; selpar_slope_cHL2_out(1,7)=set_selpar_slope_cHL2;
       //selpar_A50_cHL3_out(8)=selpar_A50_cHL3; selpar_A50_cHL3_out(1,7)=set_selpar_A50_cHL3;
       //selpar_slope_cHL3_out(8)=selpar_slope_cHL3; selpar_slope_cHL3_out(1,7)=set_selpar_slope_cHL3;

	   selpar_A50_cOT1_out(8)=selpar_A50_cOT1; selpar_A50_cOT1_out(1,7)=set_selpar_A50_cOT1;
       selpar_slope_cOT1_out(8)=selpar_slope_cOT1; selpar_slope_cOT1_out(1,7)=set_selpar_slope_cOT1;
       selpar_A502_cOT1_out(8)=selpar_A502_cOT1; selpar_A502_cOT1_out(1,7)=set_selpar_A502_cOT1;
       selpar_slope2_cOT1_out(8)=selpar_slope2_cOT1; selpar_slope2_cOT1_out(1,7)=set_selpar_slope2_cOT1;	  

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
	   selpar_A50_rHB4_out(8)=selpar_A50_rHB4; selpar_A50_rHB4_out(1,7)=set_selpar_A50_rHB4;
       selpar_slope_rHB4_out(8)=selpar_slope_rHB4; selpar_slope_rHB4_out(1,7)=set_selpar_slope_rHB4;
	   selpar_A502_rHB4_out(8)=selpar_A502_rHB4; selpar_A502_rHB4_out(1,7)=set_selpar_A502_rHB4;
       selpar_slope2_rHB4_out(8)=selpar_slope2_rHB4; selpar_slope2_rHB4_out(1,7)=set_selpar_slope2_rHB4;

	   selpar_A50_rGN2_out(8)=selpar_A50_rGN2; selpar_A50_rGN2_out(1,7)=set_selpar_A50_rGN2;
       selpar_slope_rGN2_out(8)=selpar_slope_rGN2; selpar_slope_rGN2_out(1,7)=set_selpar_slope_rGN2;
	   selpar_A50_rGN3_out(8)=selpar_A50_rGN3; selpar_A50_rGN3_out(1,7)=set_selpar_A50_rGN3;
       selpar_slope_rGN3_out(8)=selpar_slope_rGN3; selpar_slope_rGN3_out(1,7)=set_selpar_slope_rGN3;
	   selpar_A502_rGN3_out(8)=selpar_A502_rGN3; selpar_A502_rGN3_out(1,7)=set_selpar_A502_rGN3;
       selpar_slope2_rGN3_out(8)=selpar_slope2_rGN3; selpar_slope2_rGN3_out(1,7)=set_selpar_slope2_rGN3;
	   selpar_A50_rGN4_out(8)=selpar_A50_rGN4; selpar_A50_rGN4_out(1,7)=set_selpar_A50_rGN4;
       selpar_slope_rGN4_out(8)=selpar_slope_rGN4; selpar_slope_rGN4_out(1,7)=set_selpar_slope_rGN4;
	   selpar_A502_rGN4_out(8)=selpar_A502_rGN4; selpar_A502_rGN4_out(1,7)=set_selpar_A502_rGN4;
       selpar_slope2_rGN4_out(8)=selpar_slope2_rGN4; selpar_slope2_rGN4_out(1,7)=set_selpar_slope2_rGN4;
       
       selpar_A50_sCT_out(8)=selpar_A50_sCT; selpar_A50_sCT_out(1,7)=set_selpar_A50_sCT;
       selpar_slope_sCT_out(8)=selpar_slope_sCT; selpar_slope_sCT_out(1,7)=set_selpar_slope_sCT;
       selpar_A502_sCT_out(8)=selpar_A502_sCT; selpar_A502_sCT_out(1,7)=set_selpar_A502_sCT;
       selpar_slope2_sCT_out(8)=selpar_slope2_sCT; selpar_slope2_sCT_out(1,7)=set_selpar_slope2_sCT;
       
	   //selpar_Age1_sCT_logit_out(8)=selpar_logit_Age1_sCT; selpar_Age1_sCT_logit_out(1,7)=set_selpar_logit_Age1_sCT;
	   //KC  selpar_Age2_sCT_logit_out(8)=selpar_logit_Age2_sCT; selpar_Age2_sCT_logit_out(1,7)=set_selpar_logit_Age2_sCT;
	   
	   //selpar_age1logit_D_out(8)=selpar_age1logit_D; selpar_age1logit_D_out(1,7)=set_selpar_age1logit_D;

       log_q_cHL_out(8)=log_q_cpue_cHL; log_q_cHL_out(1,7)=set_log_q_cpue_cHL;
       log_q_rHB_out(8)=log_q_cpue_rHB; log_q_rHB_out(1,7)=set_log_q_cpue_rHB;
       log_q_rGN_out(8)=log_q_cpue_rGN; log_q_rGN_out(1,7)=set_log_q_cpue_rGN;
       log_q_sFT_out(8)=log_q_cpue_sFT; log_q_sFT_out(1,7)=set_log_q_cpue_sFT;
       log_q_sCT_out(8)=log_q_cpue_sCT; log_q_sCT_out(1,7)=set_log_q_cpue_sCT;
                     
       log_avg_F_cHL_out(8)=log_avg_F_L_cHL; log_avg_F_cHL_out(1,7)=set_log_avg_F_L_cHL;
       log_avg_F_cOT_out(8)=log_avg_F_L_cOT; log_avg_F_cOT_out(1,7)=set_log_avg_F_L_cOT;
       log_avg_F_cTW_out(8)=log_avg_F_L_cTW; log_avg_F_cTW_out(1,7)=set_log_avg_F_L_cTW;
       log_avg_F_rHB_out(8)=log_avg_F_L_rHB; log_avg_F_rHB_out(1,7)=set_log_avg_F_L_rHB;
       log_avg_F_rGN_out(8)=log_avg_F_L_rGN; log_avg_F_rGN_out(1,7)=set_log_avg_F_L_rGN;       
       log_avg_F_cHL_D_out(8)=log_avg_F_D_cHL; log_avg_F_cHL_D_out(1,7)=set_log_avg_F_D_cHL;
       log_avg_F_rHB_D_out(8)=log_avg_F_D_rHB; log_avg_F_rHB_D_out(1,7)=set_log_avg_F_D_rHB;
       log_avg_F_rGN_D_out(8)=log_avg_F_D_rGN; log_avg_F_rGN_D_out(1,7)=set_log_avg_F_D_rGN;
        
       log_rec_dev_out(styr_rec_dev, endyr_rec_dev)=log_dev_rec;
       log_F_dev_cHL_out(styr_L_cHL,endyr_L_cHL)=log_dev_F_L_cHL;
       log_F_dev_cOT_out(styr_L_cOT,endyr_L_cOT)=log_dev_F_L_cOT;
       log_F_dev_cTW_out(styr_L_cTW,endyr_L_cTW)=log_dev_F_L_cTW;
       log_F_dev_rHB_out(styr_L_rHB,endyr_L_rHB)=log_dev_F_L_rHB;
       log_F_dev_rGN_out(styr_L_rGN,endyr_L_rGN)=log_dev_F_L_rGN;
       log_F_dev_cHL_D_out(styr_D_cHL,endyr_D_cHL)=log_dev_F_D_cHL;
       log_F_dev_rHB_D_out(styr_D_rHB,endyr_D_rHB)=log_dev_F_D_rHB;
       log_F_dev_rGN_D_out(styr_D_rGN,endyr_D_rGN)=log_dev_F_D_rGN;
           
   #include "bam.cxx"   // write the R-compatible report

  } //endl last phase loop     
  
 
