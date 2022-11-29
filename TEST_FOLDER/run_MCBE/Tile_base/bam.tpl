//##  Author: NMFS, Beaufort Lab, Sustainable Fisheries Branch
//##  Analyst: Nikolai Klibansky
//##  Species: Golden Tilefish
//##  Region: US South Atlantic
//##  SEDAR: 66
//##  Date: 2022-10-25 15:21:35


//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//###################################################################
//##  BAM Template File #############################################
//###################################################################
//##  Author: NMFS, Beaufort Lab, Sustainable Fisheries Branch
//##  Analyst: Nikolai Klibansky
//##  Species: Golden Tilefish
//##  Region: US South Atlantic
//##  SEDAR: 66
//##  Date: 2021-05-12 17:25:37
//###################################################################
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

// starting and ending years of the specify a range of years to average recruitment over for setting alternative recruitment for years at the end of the assessment when rec devs aren't estimated
init_int styr_rec_alt1
init_int endyr_rec_alt1

// 0 = don't use alternative recruitment values, 1 = use alternative recruitment values
init_int rec_alt1_switch

//number assessment years
number nyrs;
number nyrs_rec;
//this section MUST BE INDENTED!!!
 LOCAL_CALCS
   nyrs=endyr-styr+1.;
   nyrs_rec=endyr_rec_dev-styr_rec_dev+1.;
 END_CALCS

// ending years for selectivity blocks
init_int endyr_selex_phase1;

//Total number of ages in population model
init_int nages;
// Vector of ages for age bins in population model
init_vector agebins(1,nages);
 
//Total number of ages used to match age comps: plus group may differ from popn, first age must not
init_int nages_agec;
init_int nages_agec_sBL;
//Vector of ages for age bins in age comps
init_vector agebins_agec(1,nages_agec);
init_vector agebins_agec_sBL(1,nages_agec_sBL);

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
init_int n_iter_msy;
// styr_rec_spr: Start year to compute arithmetic average recruitment for SPR-related values
init_int styr_rec_spr;
// endyr_rec_spr: End year to compute arithmetic average recruitment for SPR-related values
init_int endyr_rec_spr;

// nyrs_rec_spr: Number of years to compute arithmetic average recruitment for SPR-related values
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

//################ Commercial handline fleet (cHL) #######################################
// Comm HL  Landings (1000 lb gutted weight)
init_int styr_L_cHL;
init_int endyr_L_cHL;
init_vector obs_L_cHL(styr_L_cHL,endyr_L_cHL);
init_vector obs_cv_L_cHL(styr_L_cHL,endyr_L_cHL);

// Comm HL length Compositions (3 cm bins)
 
// Comm HL age compositions
init_int nyr_agec_cHL;
init_ivector yrs_agec_cHL(1,nyr_agec_cHL);
init_vector nsamp_agec_cHL(1,nyr_agec_cHL);
init_vector nfish_agec_cHL(1,nyr_agec_cHL);
init_matrix obs_agec_cHL(1,nyr_agec_cHL,1,nages_agec);

//################ Commercial Longline fleet (cLL) #######################################
// cLL CPUE
init_int styr_cpue_cLL;                                             
init_int endyr_cpue_cLL;                                            
init_vector obs_cpue_cLL(styr_cpue_cLL,endyr_cpue_cLL);//Observed CPUE
init_vector obs_cv_cpue_cLL(styr_cpue_cLL,endyr_cpue_cLL); //CV of cpue

// Comm Longline  Landings (1000 lb gutted weight)
init_int styr_L_cLL;
init_int endyr_L_cLL;
init_vector obs_L_cLL(styr_L_cLL,endyr_L_cLL);
init_vector obs_cv_L_cLL(styr_L_cLL,endyr_L_cLL);

// Comm Longline length Compositions (3 cm bins)
 
// Comm Longline age compositions
init_int nyr_agec_cLL;
init_ivector yrs_agec_cLL(1,nyr_agec_cLL);
init_vector nsamp_agec_cLL(1,nyr_agec_cLL);
init_vector nfish_agec_cLL(1,nyr_agec_cLL);
init_matrix obs_agec_cLL(1,nyr_agec_cLL,1,nages_agec);

//################### All Recreational fleets (rGN) ##########################################
// rGN Landings (1000 lb gutted weight)
init_int styr_L_rGN;
init_int endyr_L_rGN;
init_vector obs_L_rGN(styr_L_rGN,endyr_L_rGN);   //vector of observed landings by year 
init_vector obs_cv_L_rGN(styr_L_rGN,endyr_L_rGN);    //vector of CV of landings by year

// rGN length Compositions (3 cm bins)
init_int nyr_lenc_rGN;
init_ivector yrs_lenc_rGN(1,nyr_lenc_rGN);
init_vector nsamp_lenc_rGN(1,nyr_lenc_rGN);
init_vector nfish_lenc_rGN(1,nyr_lenc_rGN);
init_matrix obs_lenc_rGN(1,nyr_lenc_rGN,1,nlenbins);

//################### MARMAP index (sBL) ##########################################
// MARMAP index
init_int nyr_cpue_sBL;                                             
init_ivector yrs_cpue_sBL(1,nyr_cpue_sBL);  //middle year in group of years to be combined                                           
init_vector obs_cpue_sBL(1,nyr_cpue_sBL);//Observed CPUE
init_vector obs_cv_cpue_sBL(1,nyr_cpue_sBL); //CV of cpue

// Age Compositions (data from MARMAP)
init_int nyr_agec_sBL;
init_ivector yrs_agec_sBL(1,nyr_agec_sBL);  //middle year in group of years to be combined
init_vector nsamp_agec_sBL(1,nyr_agec_sBL);
init_vector nfish_agec_sBL(1,nyr_agec_sBL); 
init_matrix obs_agec_sBL(1,nyr_agec_sBL,1,nages_agec_sBL);

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: parameter section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

//##################Single Parameter values and initial guesses #################################
// Von Bert parameters in TL mm all fish
init_vector set_Linf(1,7);
init_vector set_K(1,7);
init_vector set_t0(1,7);

// Von Bert parameters in TL mm females only
init_vector set_Linf_f(1,7);
init_vector set_K_f(1,7);
init_vector set_t0_f(1,7);

//CV of length at age and its standard error all fish
init_vector set_len_cv(1,7);

//Scalar used only for computing MSST. 
init_vector set_M_constant(1,7);     

//Spawner-recruit parameters (Initial guesses or fixed values)
init_vector set_steep(1,7);         //recruitment steepness
init_vector set_log_R0(1,7);        //recruitment R0
init_vector set_R_autocorr(1,7);    //recruitment autocorrelation
init_vector set_rec_sigma(1,7);     //recruitment standard deviation in log space

init_vector set_log_dm_lenc_rGN(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_agec_cHL(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_agec_cLL(1,7);    //Dirichlet-multinomial overdispersion parameter
init_vector set_log_dm_agec_sBL(1,7);    //Dirichlet-multinomial overdispersion parameter


//Initial guesses or fixed values of estimated selectivity parameters

init_vector set_selpar_L50_cHL(1,7);
init_vector set_selpar_L50_cHL2(1,7);
init_vector set_selpar_slope_cHL(1,7);
init_vector set_selpar_slope_cHL2(1,7);

init_vector set_selpar_L50_cLL(1,7);
init_vector set_selpar_L50_cLL2(1,7);
init_vector set_selpar_slope_cLL(1,7);
init_vector set_selpar_slope_cLL2(1,7);

init_vector set_selpar_L50_rGN(1,7);
init_vector set_selpar_slope_rGN(1,7);

init_vector set_selpar_L50_sBL(1,7);
init_vector set_selpar_slope_sBL(1,7);

//--index catchability-----------------------------------------------------------------------------------
init_vector set_log_q_cpue_cLL(1,7);      //catchability coefficient (log) for cGN longline index
init_vector set_log_q_cpue_sBL(1,7);      //catchability coefficient (log) for MARMAP index

//initial F
init_vector set_F_init(1,7);  //scales initial F 
//--mean F's in log space --------------------------------
init_vector set_log_avg_F_L_cHL(1,7);
init_vector set_log_avg_F_L_cLL(1,7);
init_vector set_log_avg_F_L_rGN(1,7);

//##################Dev Vector Parameter values (vals) and bounds #################################
//--F vectors---------------------------
init_vector set_log_dev_F_L_cHL(1,3);        
init_vector set_log_dev_F_L_cLL(1,3); 
init_vector set_log_dev_F_L_rGN(1,3);
init_vector set_log_dev_rec(1,3);
init_vector set_log_dev_Nage(1,3);

init_vector set_log_dev_vals_F_L_cHL(styr_L_cHL,endyr_L_cHL);
init_vector set_log_dev_vals_F_L_cLL(styr_L_cLL,endyr_L_cLL);
init_vector set_log_dev_vals_F_L_rGN(styr_L_rGN,endyr_L_rGN);
init_vector set_log_dev_vals_rec(styr_rec_dev,endyr_rec_dev);
init_vector set_log_dev_vals_Nage(2,nages);             


//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: likelihood weights section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

init_number set_w_L;            //weight for landings
init_number set_w_cpue_cLL;         //weight for cGN handline index
init_number set_w_cpue_sBL;         //weight for MARMAP index
init_number set_w_lenc_rGN;        //weight for rGN len comps
init_number set_w_agec_cHL;        //weight for cGN handline age comps
init_number set_w_agec_cLL;        //weight for cGN longline age comps
init_number set_w_agec_sBL;        //weight for MARMAP age comps
init_number set_w_Nage_init;    //for fitting initial abundance at age (excluding first age)
init_number set_w_rec;          //for fitting S-R curve
init_number set_w_rec_early;    //additional constraint on early years recruitment
init_number set_w_rec_end;      //additional constraint on ending years recruitment 
init_number set_w_fullF;        //penalty for any Fapex>3(removed in final phase of optimization)
init_number set_w_Ftune;        //weight applied to tuning F (removed in final phase of optimization)

//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//-- BAM DATA_SECTION: miscellaneous stuff section
//--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>

//TL(mm)-weight(whole weight in mt) relationship: W=aL^b
init_number wgtpar_a;
init_number wgtpar_b;
//weight(whole weight)-gonad weight (units=g) relationship: GW=a+b*W
init_number gwgtpar_a;
init_number gwgtpar_b;
//gutted weight to whole weight conversion WW=1.05893*GW
init_number gut2whole;


//Maturity and proportion female at age
init_vector obs_maturity_f(1,nages);            //proportion females mature at age
init_vector obs_maturity_m(1,nages);            //proportion males mature at age
init_vector obs_prop_f(1,nages);     			//proportion female at age

init_number spawn_time_frac; //time of year of peak spawning, as a fraction of the year

// Natural mortality
init_vector set_M(1,nages);     //age-dependent: used in model
init_number max_obs_age;        //max observed age, used to scale M, if estimated

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
init_int set_q_RW_phase;         //value sets estimation phase of random walk, negative value turns it off
init_number set_q_RW_rec_var;     //assumed variance of RW q

//Tune Fapex (tuning removed in final year of optimization)
init_number set_Ftune;
init_int set_Ftune_yr;

//threshold sample sizes for length comps 
init_number minSS_lenc_rGN;

//threshold sample sizes for age comps
init_number minSS_agec_cHL;
init_number minSS_agec_cLL;
init_number minSS_agec_sBL;

//ageing error matrix (columns are true ages, rows are ages as read for age comps: columns should sum to one)
init_matrix age_error(1,nages,1,nages);

// #######Indexing integers for year(iyear), age(iage),length(ilen) ###############
int iyear;
int jyear;
int iage;
int ilen;
int ff;

number sqrt2pi;
number g2mt;                    //conversion of grams to metric tons 
number g2kg;                    //conversion of grams to kg   
//number g2klb;                   //conversion of grams to 1000 lb   
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
  const double Linf_LO=set_Linf(2); const double Linf_HI=set_Linf(3); const double Linf_PH=set_Linf(4);
  const double K_LO=set_K(2); const double K_HI=set_K(3); const double K_PH=set_K(4);
  const double t0_LO=set_t0(2); const double t0_HI=set_t0(3); const double t0_PH=set_t0(4);
  const double Linf_f_LO=set_Linf_f(2); const double Linf_f_HI=set_Linf_f(3); const double Linf_f_PH=set_Linf_f(4);
  const double K_f_LO=set_K_f(2); const double K_f_HI=set_K_f(3); const double K_f_PH=set_K_f(4);
  const double t0_f_LO=set_t0_f(2); const double t0_f_HI=set_t0_f(3); const double t0_f_PH=set_t0_f(4);  
  const double len_cv_LO=set_len_cv(2); const double len_cv_HI=set_len_cv(3); const double len_cv_PH=set_len_cv(4); 
  const double M_constant_LO=set_M_constant(2); const double M_constant_HI=set_M_constant(3); const double M_constant_PH=set_M_constant(4);        
  const double steep_LO=set_steep(2); const double steep_HI=set_steep(3); const double steep_PH=set_steep(4);
  const double log_R0_LO=set_log_R0(2); const double log_R0_HI=set_log_R0(3); const double log_R0_PH=set_log_R0(4);
  const double R_autocorr_LO=set_R_autocorr(2); const double R_autocorr_HI=set_R_autocorr(3); const double R_autocorr_PH=set_R_autocorr(4);
  const double rec_sigma_LO=set_rec_sigma(2); const double rec_sigma_HI=set_rec_sigma(3); const double rec_sigma_PH=set_rec_sigma(4);

  const double log_dm_lenc_rGN_LO=set_log_dm_lenc_rGN(2); const double log_dm_lenc_rGN_HI=set_log_dm_lenc_rGN(3); const double log_dm_lenc_rGN_PH=set_log_dm_lenc_rGN(4);
  const double log_dm_agec_cHL_LO=set_log_dm_agec_cHL(2); const double log_dm_agec_cHL_HI=set_log_dm_agec_cHL(3); const double log_dm_agec_cHL_PH=set_log_dm_agec_cHL(4);
  const double log_dm_agec_cLL_LO=set_log_dm_agec_cLL(2); const double log_dm_agec_cLL_HI=set_log_dm_agec_cLL(3); const double log_dm_agec_cLL_PH=set_log_dm_agec_cLL(4);
  const double log_dm_agec_sBL_LO=set_log_dm_agec_sBL(2); const double log_dm_agec_sBL_HI=set_log_dm_agec_sBL(3); const double log_dm_agec_sBL_PH=set_log_dm_agec_sBL(4);
  
  
  const double selpar_L50_cHL_LO=set_selpar_L50_cHL(2); const double selpar_L50_cHL_HI=set_selpar_L50_cHL(3); const double selpar_L50_cHL_PH=set_selpar_L50_cHL(4);
  const double selpar_L50_cHL2_LO=set_selpar_L50_cHL2(2); const double selpar_L50_cHL2_HI=set_selpar_L50_cHL2(3); const double selpar_L50_cHL2_PH=set_selpar_L50_cHL2(4);
  const double selpar_slope_cHL_LO=set_selpar_slope_cHL(2); const double selpar_slope_cHL_HI=set_selpar_slope_cHL(3); const double selpar_slope_cHL_PH=set_selpar_slope_cHL(4);
  const double selpar_slope_cHL2_LO=set_selpar_slope_cHL2(2); const double selpar_slope_cHL2_HI=set_selpar_slope_cHL2(3); const double selpar_slope_cHL2_PH=set_selpar_slope_cHL2(4);
  const double selpar_L50_cLL_LO=set_selpar_L50_cLL(2); const double selpar_L50_cLL_HI=set_selpar_L50_cLL(3); const double selpar_L50_cLL_PH=set_selpar_L50_cLL(4);
  const double selpar_L50_cLL2_LO=set_selpar_L50_cLL2(2); const double selpar_L50_cLL2_HI=set_selpar_L50_cLL2(3); const double selpar_L50_cLL2_PH=set_selpar_L50_cLL2(4);
  const double selpar_slope_cLL_LO=set_selpar_slope_cLL(2); const double selpar_slope_cLL_HI=set_selpar_slope_cLL(3); const double selpar_slope_cLL_PH=set_selpar_slope_cLL(4);
  const double selpar_slope_cLL2_LO=set_selpar_slope_cLL2(2); const double selpar_slope_cLL2_HI=set_selpar_slope_cLL2(3); const double selpar_slope_cLL2_PH=set_selpar_slope_cLL2(4);

  //const double selpar_afull_cLL_LO=set_selpar_afull_cLL(2); const double selpar_afull_cLL_HI=set_selpar_afull_cLL(3); const double selpar_afull_cLL_PH=set_selpar_afull_cLL(4);
  //const double selpar_sigma_cLL_LO=set_selpar_afull_cLL(2); const double selpar_sigma_cLL_HI=set_selpar_afull_cLL(3); const double selpar_sigma_cLL_PH=set_selpar_afull_cLL(4);

  const double selpar_L50_rGN_LO=set_selpar_L50_rGN(2); const double selpar_L50_rGN_HI=set_selpar_L50_rGN(3); const double selpar_L50_rGN_PH=set_selpar_L50_rGN(4);
  const double selpar_slope_rGN_LO=set_selpar_slope_rGN(2); const double selpar_slope_rGN_HI=set_selpar_slope_rGN(3); const double selpar_slope_rGN_PH=set_selpar_slope_rGN(4);

  const double selpar_L50_sBL_LO=set_selpar_L50_sBL(2); const double selpar_L50_sBL_HI=set_selpar_L50_sBL(3); const double selpar_L50_sBL_PH=set_selpar_L50_sBL(4);
  const double selpar_slope_sBL_LO=set_selpar_slope_sBL(2); const double selpar_slope_sBL_HI=set_selpar_slope_sBL(3); const double selpar_slope_sBL_PH=set_selpar_slope_sBL(4);

  const double log_q_cpue_cLL_LO=set_log_q_cpue_cLL(2); const double log_q_cpue_cLL_HI=set_log_q_cpue_cLL(3); const double log_q_cpue_cLL_PH=set_log_q_cpue_cLL(4);
  const double log_q_cpue_sBL_LO=set_log_q_cpue_sBL(2); const double log_q_cpue_sBL_HI=set_log_q_cpue_sBL(3); const double log_q_cpue_sBL_PH=set_log_q_cpue_sBL(4);

  const double F_init_LO=set_F_init(2); const double F_init_HI=set_F_init(3); const double F_init_PH=set_F_init(4);
  const double log_avg_F_L_cHL_LO=set_log_avg_F_L_cHL(2); const double log_avg_F_L_cHL_HI=set_log_avg_F_L_cHL(3); const double log_avg_F_L_cHL_PH=set_log_avg_F_L_cHL(4);
  const double log_avg_F_L_cLL_LO=set_log_avg_F_L_cLL(2); const double log_avg_F_L_cLL_HI=set_log_avg_F_L_cLL(3); const double log_avg_F_L_cLL_PH=set_log_avg_F_L_cLL(4); 
  const double log_avg_F_L_rGN_LO=set_log_avg_F_L_rGN(2); const double log_avg_F_L_rGN_HI=set_log_avg_F_L_rGN(3); const double log_avg_F_L_rGN_PH=set_log_avg_F_L_rGN(4); 

  //-dev vectors-----------------------------------------------------------------------------------------------------------  
  const double log_F_dev_L_cHL_LO=set_log_dev_F_L_cHL(1); const double log_F_dev_L_cHL_HI=set_log_dev_F_L_cHL(2); const double log_F_dev_L_cHL_PH=set_log_dev_F_L_cHL(3);   
  const double log_F_dev_L_cLL_LO=set_log_dev_F_L_cLL(1); const double log_F_dev_L_cLL_HI=set_log_dev_F_L_cLL(2); const double log_F_dev_L_cLL_PH=set_log_dev_F_L_cLL(3);   
  const double log_F_dev_L_rGN_LO=set_log_dev_F_L_rGN(1); const double log_F_dev_L_rGN_HI=set_log_dev_F_L_rGN(2); const double log_F_dev_L_rGN_PH=set_log_dev_F_L_rGN(3);   

  const double log_rec_dev_LO=set_log_dev_rec(1); const double log_rec_dev_HI=set_log_dev_rec(2); const double log_rec_dev_PH=set_log_dev_rec(3);          
  const double log_Nage_dev_LO=set_log_dev_Nage(1); const double log_Nage_dev_HI=set_log_dev_Nage(2); const double log_Nage_dev_PH=set_log_dev_Nage(3);          

 END_CALCS
 
////--------------Growth--------------------------------------------------------------------------- 
  init_bounded_number Linf(Linf_LO,Linf_HI,Linf_PH);
  init_bounded_number K(K_LO,K_HI,K_PH);
  init_bounded_number t0(t0_LO,t0_HI,t0_PH);
  init_bounded_number Linf_f(Linf_f_LO,Linf_f_HI,Linf_f_PH);
  init_bounded_number K_f(K_f_LO,K_f_HI,K_f_PH);
  init_bounded_number t0_f(t0_f_LO,t0_f_HI,t0_f_PH);
  init_bounded_number len_cv_val(len_cv_LO,len_cv_HI,len_cv_PH);  
  vector Linf_out(1,8);
  vector K_out(1,8);
  vector t0_out(1,8);
  vector Linf_f_out(1,8);
  vector K_f_out(1,8);
  vector t0_f_out(1,8);
  vector len_cv_val_out(1,8);
  
  vector meanlen_TL(1,nages);   //mean total length (mm) at age all fish
  vector meanlen_TL_f(1,nages);   //mean total length (mm) at age females

  vector wgt_g(1,nages);        //whole wgt in g
  vector wgt_kg(1,nages);       //whole wgt in kg
  vector wgt_mt(1,nages);       //whole wgt in mt
  vector wgt_klb(1,nages);      //whole wgt in 1000 lb
  vector wgt_lb(1,nages);       //whole wgt in lb  
  vector wgt_klb_gut(1,nages);  //gutted wgt in 1000 lb  
  vector wgt_lb_gut(1,nages);   //gutted wgt in lb
 
  vector wgt_g_f(1,nages);        //whole wgt in g
  vector wgt_kg_f(1,nages);       //whole wgt in kg
  vector wgt_mt_f(1,nages);       //whole wgt in mt
  vector wgt_klb_f(1,nages);      //whole wgt in 1000 lb
  vector wgt_lb_f(1,nages);       //whole wgt in lb  
  vector gonad_wgt_mt(1,nages); //gonad wgt in mt  
 
  matrix len_cLL_mm(styr,endyr,1,nages);          //mean length at age of commercial longline landings in mm (may differ from popn mean)
  matrix wholewgt_cLL_klb(styr,endyr,1,nages);    //whole wgt of commercial longline landings in 1000 lb   
  matrix gutwgt_cLL_klb(styr,endyr,1,nages);      //gutted wgt of commercial longline landings in 1000 lb   

  matrix len_cHL_mm(styr,endyr,1,nages);          //mean length at age of commercial handline landings in mm (may differ from popn mean)
  matrix gutwgt_cHL_klb(styr,endyr,1,nages);      //gutted wgt of commercial handline landings in 1000 lb   
  
  matrix len_rGN_mm(styr,endyr,1,nages);          //mean length at age of rGN landings in mm (may differ from popn mean)    
  matrix gutwgt_rGN_klb(styr,endyr,1,nages);    //gutted wgt of rGN landings in 1000 lb
  
  matrix len_sBL_mm(styr,endyr,1,nages);      //mean length at age of MARMAP landings in mm (may differ from popn mean)    
  matrix wgt_sBL_klb(styr,endyr,1,nages);     //whole wgt of MARMAP landings in 1000 lb  
 
  matrix lenprob(1,nages,1,nlenbins);           //distn of size at age (age-length key, 3 cm bins) in population
  number zscore_len;                            //standardized normal values used for computing lenprob
  vector cprob_lenvec(1,nlenbins);              //cumulative probabilities used for computing lenprob
  number zscore_lzero;                          //standardized normal values for length = 0
  number cprob_lzero;                           //length probability mass below zero, used for computing lenprob

  //matrices below are used to match length comps
  matrix lenprob_cHL(1,nages,1,nlenbins);     //distn of size at age in cHL
  matrix lenprob_cLL(1,nages,1,nlenbins);     //distn of size at age in cLL
  matrix lenprob_rGN(1,nages,1,nlenbins);     //distn of size at age in rGN
 
//  //init_bounded_dev_vector log_len_cv_dev(1,nages,-2,2,3)
//  number log_len_cv
  vector len_sd(1,nages);
  vector len_cv(1,nages); //for fishgraph 

//----Predicted length compositions
  matrix pred_lenc_rGN(1,nyr_lenc_rGN,1,nlenbins); 

//----Predicted age compositions
  matrix pred_agec_cHL(1,nyr_agec_cHL,1,nages_agec);
  matrix pred_agec_cHL_allages(1,nyr_agec_cHL,1,nages);
  matrix ErrorFree_agec_cHL(1,nyr_agec_cHL,1,nages);  
  matrix pred_agec_cLL(1,nyr_agec_cLL,1,nages_agec);
  matrix pred_agec_cLL_allages(1,nyr_agec_cLL,1,nages);  
  matrix ErrorFree_agec_cLL(1,nyr_agec_cLL,1,nages);
  matrix pred_agec_sBL(1,nyr_agec_sBL,1,nages_agec_sBL);
  matrix pred_agec_sBL_allages(1,nyr_agec_sBL,1,nages);  
  matrix ErrorFree_agec_sBL(1,nyr_agec_sBL,1,nages);
  
//effective sample size applied in multinomial distributions
  vector nsamp_lenc_rGN_allyr(styr,endyr);

  vector nsamp_agec_cHL_allyr(styr,endyr);
  vector nsamp_agec_cLL_allyr(styr,endyr);
  vector nsamp_agec_sBL_allyr(styr,endyr);

//Nfish used in MCB analysis (not used in fitting)
  vector nfish_lenc_rGN_allyr(styr,endyr);
  vector nfish_agec_cHL_allyr(styr,endyr);
  vector nfish_agec_cLL_allyr(styr,endyr);
  vector nfish_agec_sBL_allyr(styr,endyr);

//Computed effective sample size for output (not used in fitting)
  vector neff_lenc_rGN_allyr_out(styr,endyr); 
  vector neff_agec_cHL_allyr_out(styr,endyr);
  vector neff_agec_cLL_allyr_out(styr,endyr);
  vector neff_agec_sBL_allyr_out(styr,endyr);

//-----Population-----------------------------------------------------------------------------------
  matrix N(styr,endyr+1,1,nages);           //Population numbers by year and age at start of yr
  matrix N_mdyr(styr,endyr,1,nages);        //Population numbers by year and age at mdpt of yr: used for comps and cpue
  matrix N_spawn(styr,endyr+1,1,nages);     //Population numbers by year and age at peaking spawning: used for SSB (proj yr ok bc of ssb on Jan1)  
  init_bounded_vector log_dev_Nage(2,nages,log_Nage_dev_LO,log_Nage_dev_HI,log_Nage_dev_PH);
  vector log_Nage_dev_output(1,nages);      //used in output. equals zero for first age
  matrix B(styr,endyr+1,1,nages);           //Population biomass by year and age at start of yr
  vector totB(styr,endyr+1);                //Total biomass by year
  vector totN(styr,endyr+1);                //Total abundance by year
  vector SSB(styr,endyr+1);                 //Total spawning biomass by year (female + male mature biomass) (proj yr ok bc of ssb on Jan1)
  vector MatFemB(styr,endyr+1);               //Total spawning biomass by year (mature female biomass) (proj yr ok bc of ssb on Jan1) 
  vector rec(styr,endyr+1);                 //Recruits by year
  vector prop_m(1,nages);      //Year-dependent proportion male by age  
  vector prop_f(1,nages);      //Proportion female by age
  vector maturity_f(1,nages);               //Proportion of female mature at age
  vector maturity_m(1,nages);               //Proportion of male mature at age
  vector reprod(1,nages);                   //vector used to compute spawning biomass (total mature biomass - males + females)
  vector reprod2(1,nages);                  //vector used to compute mature female biomass 

//---Stock-Recruit Function (Beverton-Holt, steepness parameterization)----------
  init_bounded_number log_R0(log_R0_LO,log_R0_HI,log_R0_PH);        //log(virgin Recruitment)
  vector log_R0_out(1,8);
  number R0;                                  //virgin recruitment
  init_bounded_number steep(steep_LO,steep_HI,steep_PH); //steepness
  vector steep_out(1,8);
  init_bounded_number rec_sigma(rec_sigma_LO,rec_sigma_HI,rec_sigma_PH);  //sd recruitment residuals  //KC comment out; this should fix it at the initial value
  vector rec_sigma_out(1,8);
  init_bounded_number R_autocorr(R_autocorr_LO,R_autocorr_HI,R_autocorr_PH);  //autocorrelation in SR  KC commented out since not estimated
  vector R_autocorr_out(1,8);

  number rec_sigma_sq;                        //square of rec_sigma      
  number rec_logL_add;                        //additive term in -logL term   
 
  init_bounded_dev_vector log_dev_rec(styr_rec_dev,endyr_rec_dev,log_rec_dev_LO,log_rec_dev_HI,log_rec_dev_PH);
  vector log_rec_dev_output(styr,endyr+1);             //used in t.series output. equals zero except for yrs in log_dev_rec
  vector log_rec_dev_out(styr_rec_dev,endyr_rec_dev);  //used in output for bound checking
  
  number var_rec_dev;                                //variance of log recruitment deviations, from yrs with unconstrainted S-R(XXXX-XXXX)
  number sigma_rec_dev;                              //sample SD of log residuals (may not equal rec_sigma 
  
  number rec_mean_alt1_temp          //intermediate calc for geometric mean of recruitment values during alternative period
  number rec_mean_alt1               //geometric mean of recruitment values during alternative period
  number nyrs_rec_alt1               //number of years during alternative period
  
  number BiasCor;                               //Bias correction in equilibrium recruits
  number S0;                                    //equal to spr_F0*R0 = virgin SSB
  number B0;                                    //equal to bpr_F0*R0 = virgin B  
  number R1;                                    //Recruits in styr
  number R_virgin;                              //unfished recruitment with bias correction
  vector SdS0(styr,endyr+1);                      //SSB / virgin SSB (projection yr possible bc of SSB on Jan1

  init_bounded_number log_dm_lenc_rGN(log_dm_lenc_rGN_LO,log_dm_lenc_rGN_HI,log_dm_lenc_rGN_PH);
  init_bounded_number log_dm_agec_cHL(log_dm_agec_cHL_LO,log_dm_agec_cHL_HI,log_dm_agec_cHL_PH);
  init_bounded_number log_dm_agec_cLL(log_dm_agec_cLL_LO,log_dm_agec_cLL_HI,log_dm_agec_cLL_PH);
  init_bounded_number log_dm_agec_sBL(log_dm_agec_sBL_LO,log_dm_agec_sBL_HI,log_dm_agec_sBL_PH);

  vector log_dm_lenc_rGN_out(1,8);
  vector log_dm_agec_cHL_out(1,8);
  vector log_dm_agec_cLL_out(1,8);
  vector log_dm_agec_sBL_out(1,8);


  
//-----------------------------------------------------------------------------------------------------------------------------------------------
////---Selectivity-------------------------------------------------------------------------

//Commercial handline-------------------------------------------------
  matrix sel_cHL(styr,endyr,1,nages);
  init_bounded_number selpar_L50_cHL(selpar_L50_cHL_LO,selpar_L50_cHL_HI,selpar_L50_cHL_PH);
  init_bounded_number selpar_L50_cHL2(selpar_L50_cHL2_LO,selpar_L50_cHL2_HI,selpar_L50_cHL2_PH);
  init_bounded_number selpar_slope_cHL(selpar_slope_cHL_LO,selpar_slope_cHL_HI,selpar_slope_cHL_PH);
  init_bounded_number selpar_slope_cHL2(selpar_slope_cHL2_LO,selpar_slope_cHL2_HI,selpar_slope_cHL2_PH);
  vector selpar_L50_cHL_out(1,8);
  vector selpar_L50_cHL2_out(1,8);
  vector selpar_slope_cHL_out(1,8);
  vector selpar_slope_cHL2_out(1,8);
  
//Commercial longline-------------------------------------------------
  matrix sel_cLL(styr,endyr,1,nages);
  init_bounded_number selpar_L50_cLL(selpar_L50_cLL_LO,selpar_L50_cLL_HI,selpar_L50_cLL_PH);
  init_bounded_number selpar_L50_cLL2(selpar_L50_cLL2_LO,selpar_L50_cLL2_HI,selpar_L50_cLL2_PH);
  init_bounded_number selpar_slope_cLL(selpar_slope_cLL_LO,selpar_slope_cLL_HI,selpar_slope_cLL_PH);
  init_bounded_number selpar_slope_cLL2(selpar_slope_cLL2_LO,selpar_slope_cLL2_HI,selpar_slope_cLL2_PH);
  //init_bounded_number selpar_afull_cLL(selpar_afull_cLL_LO,selpar_afull_cLL_HI,selpar_afull_cLL_PH);
  //init_bounded_number selpar_sigma_cLL(selpar_sigma_cLL_LO,selpar_sigma_cLL_HI,selpar_sigma_cLL_PH);
  vector selpar_L50_cLL_out(1,8);
  vector selpar_L50_cLL2_out(1,8);
  vector selpar_slope_cLL_out(1,8);
  vector selpar_slope_cLL2_out(1,8);
  //vector selpar_afull_cLL_out(1,8);
  //vector selpar_sigma_cLL_out(1,8);

//Recreational (rGN)-------------------------------------------------
  matrix sel_rGN(styr,endyr,1,nages);
  init_bounded_number selpar_L50_rGN(selpar_L50_rGN_LO,selpar_L50_rGN_HI,selpar_L50_rGN_PH);
  init_bounded_number selpar_slope_rGN(selpar_slope_rGN_LO,selpar_slope_rGN_HI,selpar_slope_rGN_PH);
  vector selpar_L50_rGN_out(1,8);
  vector selpar_slope_rGN_out(1,8);

//MARMAP (sBL)-------------------------------------------------
  matrix sel_sBL(styr,endyr,1,nages);
  init_bounded_number selpar_L50_sBL(selpar_L50_sBL_LO,selpar_L50_sBL_HI,selpar_L50_sBL_PH);
  init_bounded_number selpar_slope_sBL(selpar_slope_sBL_LO,selpar_slope_sBL_HI,selpar_slope_sBL_PH);
  vector selpar_L50_sBL_out(1,8);
  vector selpar_slope_sBL_out(1,8);


//Weighted total selectivity--------------------------------------------  
  //effort-weighted, recent selectivities
  vector sel_wgted_L(1,nages);  //toward landings 
  vector sel_wgted_tot(1,nages);//toward Z, landings plus dead discards (none in this assmt, but kept structure)

//-----------------------------------------------------------------------------------------------------------------------------------------------
//-------CPUE Predictions--------------------------------
  vector pred_cLL_cpue(styr_cpue_cLL,endyr_cpue_cLL);                   //predicted cLL index (weight fish per effort)
  matrix N_cLL(styr_cpue_cLL,endyr_cpue_cLL,1,nages);                   //used to compute cLL index
  
  vector pred_sBL_cpue(1,nyr_cpue_sBL);                   //predicted sBL index (number fish per effort)
  matrix N_sBL(1,nyr_cpue_sBL,1,nages);                   //used to compute sBL index
  vector pred_sBL_cpue_allyr(styr,endyr);  // used for graphing purposes, fills in consec years
  vector obs_cpue_sBL_allyr(styr,endyr);   // used for graphing purposes, fills in consec years
  vector obs_cv_cpue_sBL_allyr(styr,endyr);    // used for graphing purposes, fills in consec years


//---Catchability (CPUE q's)----------------------------------------------------------
  init_bounded_number log_q_cpue_cLL(log_q_cpue_cLL_LO,log_q_cpue_cLL_HI,log_q_cpue_cLL_PH);
  init_bounded_number log_q_cpue_sBL(log_q_cpue_sBL_LO,log_q_cpue_sBL_HI,log_q_cpue_sBL_PH);
  vector log_q_cpue_cLL_out(1,8);
  vector log_q_cpue_sBL_out(1,8);
  
  //init_bounded_number q_rate(0.001,0.1,set_q_rate_phase);  //not estimated so commented out, declared as number
  number q_rate;
  vector q_rate_fcn_cLL(styr_cpue_cLL,endyr_cpue_cLL);         //increase due to technology creep (saturates in 2003) 
  
//  init_bounded_number q_DD_beta(0.1,0.9,set_q_DD_phase);    //not estimated so commented out and declared as number (below)
  number q_DD_beta;
  vector q_DD_fcn(styr,endyr);    //density dependent function as a multiple of q (scaled a la Katsukawa and Matsuda. 2003)
  number B0_q_DD;                 //B0 of ages q_DD_age plus
  vector B_q_DD(styr,endyr);      //annual biomass of ages q_DD_age plus

//Fishery dependent random walk catchability
//  init_bounded_vector q_RW_log_dev_rHB(styr_rHB_cpue,endyr_rHB_cpue-1,-3.0,3.0,set_q_RW_phase); //NOT estimated in this model 
 vector q_RW_log_dev_cLL(styr_cpue_cLL,endyr_cpue_cLL-1); 

//Catchability vector over time, may be constant
  vector q_cLL(styr_cpue_cLL,endyr_cpue_cLL); 
  number q_sBL;
 

//----------------------------------------------------------------------------------------------------------------------------------------------- 
//---Landings in numbers (total or 1000 fish) and in wgt (gutted klb)--------------------------------------------------
  matrix L_cHL_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_cHL_klb(styr,endyr,1,nages);   //landings (1000 lb gutted weight) at age    
  vector pred_cHL_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages  
  vector pred_cHL_L_klb(styr,endyr);      //yearly landings in 1000 lb gutted summed over ages
  //vector obs_L_cHL_wbias(styr,endyr);     //yearly landings observed, perhaps adjusted for multiplicitive bias

  matrix L_cLL_num(styr,endyr,1,nages);   //landings (numbers) at age
  matrix L_cLL_klb(styr,endyr,1,nages);   //landings (1000 lb gutted weight) at age    
  vector pred_cLL_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages  
  vector pred_cLL_L_klb(styr,endyr);      //yearly landings in 1000 lb gutted summed over ages

  matrix L_rGN_num(styr,endyr,1,nages);   //landings (numbers) at age
matrix L_rGN_klb(styr,endyr,1,nages);   //landings (1000 lb gutted weight) at age
  vector pred_rGN_L_knum(styr,endyr);     //yearly landings in 1000 fish summed over ages
  vector pred_rGN_L_klb(styr,endyr);      //yearly landings in 1000 lb gutted summed over ages

  matrix L_total_num(styr,endyr,1,nages);//total landings in number at age
  matrix L_total_klb(styr,endyr,1,nages);//landings in klb gutted wgt at age 
  vector L_total_knum_yr(styr,endyr);    //total landings in 1000 fish by yr summed over ages  
  vector L_total_klb_yr(styr,endyr);     //total landings (klb gutted wgt) by yr summed over ages
  

////---MSY calcs----------------------------------------------------------------------------
  number F_cHL_prop;       //proportion of F_sum attributable to cHL, last X=selpar_n_yrs_wgted yrs
  number F_cLL_prop;       //proportion of F_sum attributable to cLL, last X=selpar_n_yrs_wgted yrs
  number F_rGN_prop;       //proportion of F_sum attributable to rGN, last X=selpar_n_yrs_wgted yrs
  
  number F_init_cHL_prop;  //proportion of F_init attributable to cHL, first X yrs, No diving or discards in initial yrs
  number F_init_cLL_prop;  //proportion of F_init attributable to cLL, first X yrs
  number F_init_rGN_prop;  //proportion of F_init attributable to rGN, first X yrs
  
  number F_temp_sum;      //sum of geom mean Fsum's in last X yrs, used to compute F_fishery_prop

  vector F_end(1,nages);
  vector F_end_L(1,nages);  
  vector F_end_D(1,nages);    
  number F_end_apex;
  
  number SSB_msy_out;           //SSB (total mature biomass) at msy
  number F_msy_out;             //F at msy
  number msy_klb_out;           //max sustainable yield (1000 lb gutted wgt)
  number msy_knum_out;          //max sustainable yield (1000 fish)  
  number B_msy_out;             //total biomass at MSY 
  number R_msy_out;             //equilibrium recruitment at F=Fmsy
  number spr_msy_out;           //spr at F=Fmsy
  
  // Stuff that goes into spr.brps matrix in cxx file
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
  number rec_mean;  			     //arithmetic average recruitment used in SPR-related quantities

  vector N_age_msy(1,nages);         //numbers at age for MSY calculations: beginning of yr
  vector N_age_msy_spawn(1,nages);   //numbers at age for MSY calculations: time of peak spawning  
  vector L_age_msy(1,nages);         //landings at age for MSY calculations
  vector Z_age_msy(1,nages);         //total mortality at age for MSY calculations
  vector F_L_age_msy(1,nages);       //fishing mortality landings (not discards) at age for MSY calculations
  vector F_msy(1,n_iter_msy);        //values of full F to be used in equilibrium calculations
  vector spr_msy(1,n_iter_msy);      //reproductive capacity-per-recruit values corresponding to F values in F_msy
  vector R_eq(1,n_iter_msy);         //equilibrium recruitment values corresponding to F values in F_msy
  vector L_eq_klb(1,n_iter_msy);     //equilibrium landings(klb gutted wgt) values corresponding to F values in F_msy
  vector L_eq_knum(1,n_iter_msy);    //equilibrium landings(1000 fish) values corresponding to F values in F_msy
  vector SSB_eq(1,n_iter_msy);       //equilibrium reproductive capacity values corresponding to F values in F_msy
  vector B_eq(1,n_iter_msy);         //equilibrium biomass values corresponding to F values in F_msy
  
  vector FdF_msy(styr,endyr);
  vector FdF30(styr,endyr);
  vector SdSSB_msy(styr,endyr+1);	 //(proj yr ok bc of ssb on Jan1)	 
  number SdSSB_msy_end;
  number FdF_msy_end;
  number FdF_msy_end_mean;           //geometric mean of last X yrs  
  vector SdSSB_F30(styr,endyr+1);	 //(proj yr ok bc of ssb on Jan1)	
  vector Sdmsst_F30(styr,endyr+1);	 //(proj yr ok bc of ssb on Jan1)		 
  number SdSSB_F30_end;
  number Sdmsst_F30_end;
  number FdF30_end_mean;             //geometric mean of last selpar_n_yrs_wgted yrs  
  number Fend_mean_temp;			 //intermediate calc for geometric mean of last selpar_n_yrs_wgted yrs
  number Fend_mean;					 //geometric mean of last selpar_n_yrs_wgted yrs
  vector L_age_F30(1,nages);         //landings at age for F30 calculations
  //vector D_age_F30(1,nages);         //discard mortality (dead discards) at age for F30 calculations
 
  
  vector wgt_wgted_L_klb(1,nages);   //fishery-weighted average weight at age of landings in gutted weight
  number wgt_wgted_L_denom;          //used in intermediate calculations

  number iter_inc_msy;               //increments used to compute msy, equals 1/(n_iter_msy-1)
  
////--------Mortality------------------------------------------------------------------

// Stuff immediately below used only if M is estimated
//  //init_bounded_number M_constant(0.1,0.2,1);                         //age-indpendent: used only for MSST
//  vector Mscale_ages(1,max_obs_age);
//  vector Mscale_len(1,max_obs_age);
//  vector Mscale_wgt_g(1,max_obs_age); 
//  vector M_lorenzen(1,max_obs_age);   
//  number cum_surv_1plus;  

  vector M(1,nages);                         //age-dependent natural mortality
  init_bounded_number M_constant(M_constant_LO,M_constant_HI,M_constant_PH);                         //age-indpendent: used only for MSST
  vector M_constant_out(1,8);
  number smsy2msst;                           //scales Smsy to get msst using (1-M). Used only in output.
  number smsy2msst75;                         //scales Smsy to get msst using 75%. Used only in output.  
  
  matrix F(styr,endyr,1,nages);
  vector Fsum(styr,endyr);                   //Full fishing mortality rate by year
  vector Fapex(styr,endyr);                  //Max across ages, fishing mortality rate by year (may differ from Fsum bc of dome-shaped sel 
  matrix Z(styr,endyr,1,nages);

  init_bounded_number log_avg_F_L_cHL(log_avg_F_L_cHL_LO,log_avg_F_L_cHL_HI,log_avg_F_L_cHL_PH);
  vector log_avg_F_L_cHL_out(1,8);  
  init_bounded_dev_vector log_dev_F_L_cHL(styr_L_cHL,endyr_L_cHL,log_F_dev_L_cHL_LO,log_F_dev_L_cHL_HI,log_F_dev_L_cHL_PH);
  vector log_F_dev_L_cHL_out(styr_L_cHL,endyr_L_cHL);
  matrix F_cHL(styr,endyr,1,nages);
  vector F_cHL_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_cHL;
  number log_F_dev_end_cHL; 

  init_bounded_number log_avg_F_L_cLL(log_avg_F_L_cLL_LO,log_avg_F_L_cLL_HI,log_avg_F_L_cLL_PH);
  vector log_avg_F_L_cLL_out(1,8);  
  init_bounded_dev_vector log_dev_F_L_cLL(styr_L_cLL,endyr_L_cLL,log_F_dev_L_cLL_LO,log_F_dev_L_cLL_HI,log_F_dev_L_cLL_PH);
  vector log_F_dev_L_cLL_out(styr_L_cLL,endyr_L_cLL);
  matrix F_cLL(styr,endyr,1,nages);
  vector F_cLL_out(styr,endyr);        //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_cLL;         //cLL landings do not extend to beginning of time series
  number log_F_dev_end_cLL; 

  init_bounded_number log_avg_F_L_rGN(log_avg_F_L_rGN_LO,log_avg_F_L_rGN_HI,log_avg_F_L_rGN_PH);
  vector log_avg_F_L_rGN_out(1,8); 
  init_bounded_dev_vector log_dev_F_L_rGN(styr_L_rGN,endyr_L_rGN,log_F_dev_L_rGN_LO,log_F_dev_L_rGN_HI,log_F_dev_L_rGN_PH);    
  vector log_F_dev_L_rGN_out(styr_L_rGN,endyr_L_rGN);
  matrix F_rGN(styr,endyr,1,nages);
  vector F_rGN_out(styr,endyr); //used for intermediate calculations in fcn get_mortality
  number log_F_dev_init_rGN;    
  number log_F_dev_end_rGN;  

  init_bounded_number F_init(F_init_LO,F_init_HI,F_init_PH); //scales early F for initialization
  vector F_init_out(1,8); 
  number F_init_denom;  //interim calculation
  vector sel_initial(1,nages);       //initial selectivity (a combination of for-hire and commercial selectivities)

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
  vector L_spr(1,n_iter_spr);        //landings(lb gutted)-per-recruit (ypr) values corresponding to F values in F_spr

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
  number sdnr_lc_cLL;
  number sdnr_lc_rGN;
  
  number sdnr_ac_cHL;
  number sdnr_ac_cLL;
  number sdnr_ac_sBL;
    
  number sdnr_I_cLL;
  number sdnr_I_sBL;
        
////-------Objective function components-----------------------------------------------------------------------------
  number w_L;
   
  number w_cpue_cLL;
  number w_cpue_sBL;
   
  number w_lenc_rGN;
  
  number w_agec_cHL; 
  number w_agec_cLL;
  number w_agec_sBL;
  
  number w_Nage_init;  
  number w_rec;
  number w_rec_early;
  number w_rec_end;
  number w_fullF;  
  number w_Ftune;

  number f_cHL_L; 
  number f_cLL_L; 
  number f_rGN_L; 

  number f_cLL_cpue;
  number f_sBL_cpue;
   
  number f_lenc_rGN;

  number f_agec_cHL;
  number f_agec_cLL;
  number f_agec_sBL;       

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
   
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//INITIALIZATION_SECTION


//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
GLOBALS_SECTION
  #include "admodel.h"          // Include AD class definitions
  #include "admb2r.cpp"    // Include S-compatible output functions (needs preceding)
  #include <time.h>
	time_t start,finish;
	long hour,minute,second;	
	double elapsed_time;
	
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
RUNTIME_SECTION
 maximum_function_evaluations 1000, 1000, 2000, 3000, 10000;
 convergence_criteria 1e-2, 1e-3, 1e-3, 1e-3, 1e-4;
 
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
//##--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>--><>
PRELIMINARY_CALCS_SECTION

// Set values of fixed parameters or set initial guess of estimated parameters
  Linf=set_Linf(1);
  K=set_K(1);
  t0=set_t0(1);
  Linf_f=set_Linf_f(1);
  K_f=set_K_f(1);
  t0_f=set_t0_f(1);
  len_cv_val=set_len_cv(1);

  M=set_M; 
  M_constant=set_M_constant(1);
  smsy2msst=1.0-M_constant;
  smsy2msst75=0.75;  
//  for (iage=1;iage<=max_obs_age;iage++){Mscale_ages(iage)=iage;}
  
  log_R0=set_log_R0(1);
  steep=set_steep(1);
  R_autocorr=set_R_autocorr(1);
  rec_sigma=set_rec_sigma(1);
 
  log_dm_lenc_rGN=set_log_dm_lenc_rGN(1);
  log_dm_agec_cHL=set_log_dm_agec_cHL(1);
  log_dm_agec_cLL=set_log_dm_agec_cLL(1);
  log_dm_agec_sBL=set_log_dm_agec_sBL(1);
 
  log_q_cpue_cLL=set_log_q_cpue_cLL(1);
  log_q_cpue_sBL=set_log_q_cpue_sBL(1);
    
  q_rate=set_q_rate;
  q_rate_fcn_cLL=1.0;   
  
  q_RW_log_dev_cLL.initialize(); 
  //q_RW_log_dev_sBL.initialize(); 
        
  if (set_q_rate_phase<0 & q_rate!=0.0)
  {
    for (iyear=styr_cpue_cLL; iyear<=endyr_cpue_cLL; iyear++)
      {   if (iyear>styr_cpue_cLL & iyear <=2003) 
          {//q_rate_fcn_cLL(iyear)=(1.0+q_rate)*q_rate_fcn_cLL(iyear-1); //compound
             q_rate_fcn_cLL(iyear)=(1.0+(iyear-styr_cpue_cLL)*q_rate)*q_rate_fcn_cLL(styr_cpue_cLL);  //linear
          }
          if (iyear>2003) {q_rate_fcn_cLL(iyear)=q_rate_fcn_cLL(iyear-1);} 
      }   
  } //end q_rate conditional      

  w_L=set_w_L;
  
  w_cpue_cLL=set_w_cpue_cLL;
  w_cpue_sBL=set_w_cpue_sBL;

  w_lenc_rGN=set_w_lenc_rGN;
    
  w_agec_cHL=set_w_agec_cHL;
  w_agec_cLL=set_w_agec_cLL;
  w_agec_sBL=set_w_agec_sBL;

  w_Nage_init=set_w_Nage_init;
  w_rec=set_w_rec;
  w_rec_early=set_w_rec_early;
  w_rec_end=set_w_rec_end;
  w_fullF=set_w_fullF;
  w_Ftune=set_w_Ftune;

  F_init=set_F_init(1);

  log_avg_F_L_cHL=set_log_avg_F_L_cHL(1);
  log_avg_F_L_cLL=set_log_avg_F_L_cLL(1);
  log_avg_F_L_rGN=set_log_avg_F_L_rGN(1); 
    
  log_dev_F_L_cHL=set_log_dev_vals_F_L_cHL;
  log_dev_F_L_cLL=set_log_dev_vals_F_L_cLL;
  log_dev_F_L_rGN=set_log_dev_vals_F_L_rGN;
 
  selpar_L50_cHL=set_selpar_L50_cHL(1);
  selpar_L50_cHL2=set_selpar_L50_cHL2(1);
  selpar_slope_cHL=set_selpar_slope_cHL(1);
  selpar_slope_cHL2=set_selpar_slope_cHL2(1);

  selpar_L50_cLL=set_selpar_L50_cLL(1);
  selpar_L50_cLL2=set_selpar_L50_cLL2(1);
  selpar_slope_cLL=set_selpar_slope_cLL(1);
  selpar_slope_cLL2=set_selpar_slope_cLL2(1);
  //selpar_afull_cLL=set_selpar_afull_cLL(1);
  //selpar_sigma_cLL=set_selpar_afull_cLL(1);

  selpar_L50_rGN=set_selpar_L50_rGN(1);
  selpar_slope_rGN=set_selpar_slope_rGN(1);

  selpar_L50_sBL=set_selpar_L50_sBL(1);
  selpar_slope_sBL=set_selpar_slope_sBL(1);


 sqrt2pi=sqrt(2.*3.14159265);
 g2mt=0.000001;         //conversion of grams to metric tons
 g2kg=0.001;            //conversion of grams to kg 
 mt2klb=2.20462;        //conversion of metric tons to 1000 lb 
 mt2lb=mt2klb*1000.0;   //conversion of metric tons to lb
 //g2klb=g2mt*mt2klb;     //conversion of grams to 1000 lb 
 dzero=0.00001;         
 huge_number=1.0e+10;   
 
 SSB_msy_out=0.0;

 iter_inc_msy=max_F_spr_msy/(n_iter_msy-1);
 iter_inc_spr=max_F_spr_msy/(n_iter_spr-1); 

 maturity_f=obs_maturity_f;
 maturity_m=obs_maturity_m; 
 prop_f=obs_prop_f;
 
 //lbins=lenbins; //NOT NEEDED
 
//Fill in sample sizes of comps, possibly sampled in nonconsec yrs 
//Used primarily for output in R object   

	  nsamp_lenc_rGN_allyr=missing;
      nsamp_agec_cHL_allyr=missing;
      nsamp_agec_cLL_allyr=missing;
      nsamp_agec_sBL_allyr=missing;

      nfish_lenc_rGN_allyr=missing;
      nfish_agec_cHL_allyr=missing;
      nfish_agec_cLL_allyr=missing;
      nfish_agec_sBL_allyr=missing;
	  
	  pred_sBL_cpue_allyr=missing;
      obs_cpue_sBL_allyr=missing;
      obs_cv_cpue_sBL_allyr=missing;

   


      for (iyear=1; iyear<=nyr_lenc_rGN; iyear++)
         {if (nsamp_lenc_rGN(iyear)>=minSS_lenc_rGN)
            {nsamp_lenc_rGN_allyr(yrs_lenc_rGN(iyear))=nsamp_lenc_rGN(iyear);
             nfish_lenc_rGN_allyr(yrs_lenc_rGN(iyear))=nfish_lenc_rGN(iyear);}}

      for (iyear=1; iyear<=nyr_agec_cHL; iyear++)
         {if (nsamp_agec_cHL(iyear)>=minSS_agec_cHL)
           {nsamp_agec_cHL_allyr(yrs_agec_cHL(iyear))=nsamp_agec_cHL(iyear);
            nfish_agec_cHL_allyr(yrs_agec_cHL(iyear))=nfish_agec_cHL(iyear);}}

      for (iyear=1; iyear<=nyr_agec_cLL; iyear++)
         {if (nsamp_agec_cLL(iyear)>=minSS_agec_cLL)
           {nsamp_agec_cLL_allyr(yrs_agec_cLL(iyear))=nsamp_agec_cLL(iyear);
            nfish_agec_cLL_allyr(yrs_agec_cLL(iyear))=nfish_agec_cLL(iyear);}}

      for (iyear=1; iyear<=nyr_agec_sBL; iyear++)
         {if (nsamp_agec_sBL(iyear)>=minSS_agec_sBL)
           {nsamp_agec_sBL_allyr(yrs_agec_sBL(iyear))=nsamp_agec_sBL(iyear);
            nfish_agec_sBL_allyr(yrs_agec_sBL(iyear))=nfish_agec_sBL(iyear);}} 
			
      for (iyear=1; iyear<=nyr_cpue_sBL; iyear++)
         {
           obs_cpue_sBL_allyr(yrs_cpue_sBL(iyear))=obs_cpue_sBL(iyear);
           obs_cv_cpue_sBL_allyr(yrs_cpue_sBL(iyear))=obs_cv_cpue_sBL(iyear);
         }		

 		 
             
//fill in Fs for msy and per-recruit analyses
  F_msy(1)=0.0;  
  for (ff=2;ff<=n_iter_msy;ff++) {F_msy(ff)=F_msy(ff-1)+iter_inc_msy;}
  F_spr(1)=0.0;  
  for (ff=2;ff<=n_iter_spr;ff++) {F_spr(ff)=F_spr(ff-1)+iter_inc_spr;}


//fill in F's, Catch matrices, and log rec dev with zero's
  F_cHL.initialize(); L_cHL_num.initialize();
  F_cLL.initialize(); L_cLL_num.initialize();
  F_rGN.initialize(); L_rGN_num.initialize();

  F_cHL_out.initialize();
  F_cLL_out.initialize();
  F_rGN_out.initialize();
      
  sel_cHL.initialize();
  sel_cLL.initialize();
  sel_rGN.initialize();
  sel_sBL.initialize();
  
  log_rec_dev_output.initialize();  
  log_dev_rec=set_log_dev_vals_rec;
  log_Nage_dev_output.initialize();
  log_dev_Nage=set_log_dev_vals_Nage;

  f_cHL_L.initialize(); 
  f_cLL_L.initialize(); 
  f_rGN_L.initialize(); 

  f_cLL_cpue.initialize();
  f_sBL_cpue.initialize();
   
  f_lenc_rGN.initialize();

  f_agec_cHL.initialize();
  f_agec_cLL.initialize();
  f_agec_sBL.initialize();       

  f_Nage_init.initialize();              
  f_rec_dev.initialize();                
  f_rec_dev_early.initialize();          
  f_rec_dev_end.initialize();            
  f_fullF_constraint.initialize();       
  f_Ftune.initialize();                  
  f_priors.initialize();                 
  
  fval.initialize();
  fval_data.initialize();  
  
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
    meanlen_TL=Linf*(1.0-mfexp(-K*(agebins-t0+0.5)));    //mean total length at age in mm
    wgt_mt=wgtpar_a*pow(meanlen_TL,wgtpar_b);             //mean whole wgt at age in mt
    wgt_g=wgt_mt/g2mt;                                    //mean whole wgt at age in grams
    wgt_kg=g2kg*wgt_g;                                   //mean whole wgt at age in kilograms
    wgt_klb=mt2klb*wgt_mt;                               //mean whole weight at age in 1000 lbs
    wgt_lb=mt2lb*wgt_mt;                                 //mean whole weight at age in lbs
	wgt_klb_gut=wgt_klb/gut2whole;                            //mean gutted weight at age in 1000 lbs
    wgt_lb_gut=wgt_lb/gut2whole;                              //mean gutted weight at age in lbs

    meanlen_TL_f=Linf_f*(1.0-mfexp(-K_f*(agebins-t0_f+0.5)));    //total length in mm
    wgt_mt_f=wgtpar_a*pow(meanlen_TL_f,wgtpar_b);             //wgt in mt
    wgt_g_f=wgt_mt_f/g2mt;                                    //wgt in grams
    wgt_kg_f=g2kg*wgt_g_f;                                   //wgt in kilograms 
    wgt_klb_f=mt2klb*wgt_mt_f;                               //1000 lb of whole wgt
    wgt_lb_f=mt2lb*wgt_mt_f;                                 //lb of whole wgt
    
	gonad_wgt_mt=g2mt*mfexp(gwgtpar_a+gwgtpar_b*log(wgt_g_f));    //gonad wgt in mt
    for (iage=1;iage<=nages;iage++)
    {
      if(gonad_wgt_mt(iage)<0){gonad_wgt_mt(iage)=0;}
    }


FUNCTION get_reprod 
   //reprod is product of stuff going into reproductive capacity calcs
  //for (iyear=styr; iyear<=endyr; iyear++)
  //{
   //reprod(iyear)=elem_prod((elem_prod(prop_f(iyear),maturity_f)+elem_prod((prop_m(iyear)),maturity_m)),wgt_mt);       
   //reprod2(iyear)=elem_prod(elem_prod(prop_f(iyear),maturity_f),wgt_mt); 
  //} 
   reprod=elem_prod(elem_prod(prop_f,maturity_f),gonad_wgt_mt);  
   reprod2=elem_prod(elem_prod(prop_f,maturity_f),wgt_mt_f); 

FUNCTION get_length_at_age_dist
  //compute matrix of length at age, based on the normal distribution
  
  for (iage=1;iage<=nages;iage++)
  {
    //len_cv(iage)=mfexp(log_len_cv+log_len_cv_dev(iage));
    len_cv(iage)=len_cv_val;
    len_sd(iage)=meanlen_TL(iage)*len_cv(iage);

    zscore_lzero=(0.0-meanlen_TL(iage))/len_sd(iage);
    cprob_lzero=cumd_norm(zscore_lzero);
    
    //first length bin
    zscore_len=((lenbins(1)+0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage);
    cprob_lenvec(1)=cumd_norm(zscore_len);          //includes any probability mass below zero
    lenprob(iage,1)=cprob_lenvec(1)-cprob_lzero;    //removes any probability mass below zero
         
    //most other length bins     
    for (ilen=2;ilen<nlenbins;ilen++)
      {
        zscore_len=((lenbins(ilen)+0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage);
        cprob_lenvec(ilen)=cumd_norm(zscore_len);
        lenprob(iage,ilen)=cprob_lenvec(ilen)-cprob_lenvec(ilen-1);
      }
    //last length bin is a plus group
    zscore_len=((lenbins(nlenbins)-0.5*lenbins_width)-meanlen_TL(iage)) / len_sd(iage);
    lenprob(iage,nlenbins)=1.0-cumd_norm(zscore_len);
  
    lenprob(iage)=lenprob(iage)/(1.0-cprob_lzero);  //renormalize to account for any prob mass below size=0
  }
  //fleet and survey specific length probs, all assumed here to equal the popn
  lenprob_cHL=lenprob;
  lenprob_cLL=lenprob;
  lenprob_rGN=lenprob;   
  

FUNCTION get_weight_at_age_landings
  
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    len_cHL_mm(iyear)=meanlen_TL;
    gutwgt_cHL_klb(iyear)=wgt_klb_gut; 
    
	len_cLL_mm(iyear)=meanlen_TL;
    gutwgt_cLL_klb(iyear)=wgt_klb_gut; 
	wholewgt_cLL_klb(iyear)=wgt_klb; //wholeweight used to match index
    
    len_rGN_mm(iyear)=meanlen_TL;
    gutwgt_rGN_klb(iyear)=wgt_klb_gut;
	
	len_sBL_mm(iyear)=meanlen_TL;
    wgt_sBL_klb(iyear)=wgt_klb;    

    
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
  bpr_F0=sum(elem_prod(N_bpr_F0,wgt_mt));    

FUNCTION get_selectivity
  for (iyear=styr; iyear<=endyr_selex_phase1; iyear++)
   {sel_cHL(iyear)=logistic(agebins, selpar_L50_cHL, selpar_slope_cHL);
   sel_cLL(iyear)=logistic(agebins, selpar_L50_cLL, selpar_slope_cLL);}

  for (iyear=(endyr_selex_phase1+1); iyear<=endyr; iyear++)
   {sel_cHL(iyear)=logistic(agebins, selpar_L50_cHL2, selpar_slope_cHL2);
   sel_cLL(iyear)=logistic(agebins, selpar_L50_cLL2, selpar_slope_cLL2);}

  for (iyear=styr; iyear<=endyr; iyear++)
   {sel_rGN(iyear)=logistic(agebins, selpar_L50_rGN, selpar_slope_rGN);
   sel_sBL(iyear)=logistic(agebins, selpar_L50_sBL, selpar_slope_sBL); }

   sel_initial=sel_cLL(styr);
  
FUNCTION get_mortality
  Fsum.initialize();
  Fapex.initialize();
  F.initialize();
  //initialization F is avg from first 3 yrs of observed landings
  log_F_dev_init_cHL=sum(log_dev_F_L_cHL(styr_L_cHL,(styr_L_cHL+2)))/3.0;         
  log_F_dev_init_cLL=sum(log_dev_F_L_cLL(styr_L_cLL,(styr_L_cLL+2)))/3.0;         
  log_F_dev_init_rGN=sum(log_dev_F_L_rGN(styr_L_rGN,(styr_L_rGN+2)))/3.0;   
  
  
  for (iyear=styr; iyear<=endyr; iyear++) 
  {
    if(iyear>=styr_L_cHL & iyear<=endyr_L_cHL)
    {  F_cHL_out(iyear)=mfexp(log_avg_F_L_cHL+log_dev_F_L_cHL(iyear)); //}    
    // if (iyear<styr_L_cHL){F_cHL_out(iyear)=mfexp(log_avg_F_L_cHL+log_F_dev_init_cHL);}        
       F_cHL(iyear)=sel_cHL(iyear)*F_cHL_out(iyear);
       Fsum(iyear)+=F_cHL_out(iyear);
    }

    if(iyear>=styr_L_cLL & iyear<=endyr_L_cLL)
    {  F_cLL_out(iyear)=mfexp(log_avg_F_L_cLL+log_dev_F_L_cLL(iyear)); //}    
    // if (iyear<styr_L_cLL){F_cLL_out(iyear)=mfexp(log_avg_F_L_cLL+log_F_dev_init_cLL);}        
       F_cLL(iyear)=sel_cLL(iyear)*F_cLL_out(iyear);
       Fsum(iyear)+=F_cLL_out(iyear);
    }

    if(iyear>=styr_L_rGN & iyear<=endyr_L_rGN)
    {  F_rGN_out(iyear)=mfexp(log_avg_F_L_rGN+log_dev_F_L_rGN(iyear)); //}    
    // if (iyear<styr_L_rGN){F_rGN_out(iyear)=mfexp(log_avg_F_L_rGN+log_F_dev_init_rGN);}        
       F_rGN(iyear)=sel_rGN(iyear)*F_rGN_out(iyear);
       Fsum(iyear)+=F_rGN_out(iyear);
    }
   
    //Total F at age
    F(iyear)=F_cHL(iyear);  //first in additive series (NO +=)
    F(iyear)+=F_cLL(iyear);
    F(iyear)+=F_rGN(iyear);
    
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
  //R_virgin=(R0/((5.0*steep-1.0)*spr_F0))*
  //               (BiasCor*4.0*steep*spr_F0-spr_F0*(1.0-steep));
  //R_virgin=R0/(spr_F0/spr_F0)*BiasCor*(1.0+log(spr_F0/spr_F0)/steep); //Ricker
 
  R_virgin=SR_eq_func(R0, steep, spr_F0, spr_F0, BiasCor, SR_switch);
 
  B0=bpr_F0*R_virgin;   
  B0_q_DD=R_virgin*sum(elem_prod(N_bpr_F0(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages))); 
  
  //F_init_denom=mfexp(log_avg_F_L_cHL+log_F_dev_init_cHL)+mfexp(log_avg_F_rHB+log_F_dev_init_rHB)+mfexp(log_avg_F_L_rGN+log_F_dev_init_rGN);
  //F_init_cHL_prop=mfexp(log_avg_F_L_cHL+log_F_dev_init_cHL)/F_init_denom;
  //F_init_rHB_prop=mfexp(log_avg_F_rHB+log_F_dev_init_rHB)/F_init_denom;
  //F_init_rGN_prop=mfexp(log_avg_F_L_rGN+log_F_dev_init_rGN)/F_init_denom;

  F_initial=sel_initial*F_init;

  //F_initial=sel_cHL(styr)*F_init*F_init_cHL_prop+
            //sel_rHB(styr)*F_init*F_init_rHB_prop+
            //sel_rHB(styr)*F_init*F_init_rGN_prop; //rGN uses rHB selex
  Z_initial=M+F_initial;

//Initial equilibrium age structure
  N_spr_initial(1)=1.0*mfexp(-1.0*Z_initial(1)*spawn_time_frac); //at peak spawning time;
  for (iage=2; iage<=nages; iage++)
    {
      N_spr_initial(iage)=N_spr_initial(iage-1)*
                   mfexp(-1.0*(Z_initial(iage-1)*(1.0-spawn_time_frac) + Z_initial(iage)*spawn_time_frac)); 
    }
  N_spr_initial(nages)=N_spr_initial(nages)/(1.0-mfexp(-1.0*Z_initial(nages))); //plus group
  spr_initial=sum(elem_prod(N_spr_initial,reprod));

  //if (styr==styr_rec_dev) {R1=(R0/((5.0*steep-1.0)*spr_initial))*
  //               (4.0*steep*spr_initial-spr_F0*(1.0-steep));} //without bias correction (deviation added later)
  //else {R1=(R0/((5.0*steep-1.0)*spr_initial))*
  //               (BiasCor*4.0*steep*spr_initial-spr_F0*(1.0-steep));} //with bias correction                 

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
  MatFemB(styr)=sum(elem_prod(N_spawn(styr),reprod2));
  B_q_DD(styr)=sum(elem_prod(N(styr)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages)));
    
//Rest of years
	nyrs_rec_alt1 = endyr_rec_alt1-styr_rec_alt1+1.0; // Number of years in time period for estimated alternative recruitment value
 
  for (iyear=styr; iyear<endyr; iyear++)
  {
    if(iyear<(styr_rec_dev-1)||iyear>(endyr_rec_dev-1)) //recruitment follows S-R curve (with bias correction) exactly
    {
        N(iyear+1,1)=SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch)*BiasCor;
		
		if(rec_alt1_switch==1 && iyear>(endyr_rec_dev-1)) // recruitment is estimated as geometric mean of a specified set of years
		{ 
			rec_mean_alt1_temp=1.0;
			for (jyear=1; jyear<=nyrs_rec_alt1; jyear++) {rec_mean_alt1_temp*=mfexp(log_dev_rec(styr_rec_alt1+jyear-1));}			
			rec_mean_alt1=log(pow(rec_mean_alt1_temp,(1.0/nyrs_rec_alt1)));
			log_rec_dev_output(iyear+1)=rec_mean_alt1;
			N(iyear+1,1)=SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch)*mfexp(rec_mean_alt1);
		}
			
    }
    else   //recruitment follows S-R curve with lognormal deviation
    {
        N(iyear+1,1)=SR_func(R0, steep, spr_F0, SSB(iyear),SR_switch)*mfexp(log_dev_rec(iyear+1));
    }
	    N(iyear+1)(2,nages)=++elem_prod(N(iyear)(1,nages-1),(mfexp(-1.*Z(iyear)(1,nages-1))));
        N(iyear+1,nages)+=N(iyear,nages)*mfexp(-1.*Z(iyear,nages)); //plus group
        N_mdyr(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*0.5))); //mid year 
        N_spawn(iyear+1)(1,nages)=elem_prod(N(iyear+1)(1,nages),(mfexp(-1.*(Z(iyear+1)(1,nages))*spawn_time_frac))); //peak spawning time 
        SSB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod)); 
        MatFemB(iyear+1)=sum(elem_prod(N_spawn(iyear+1),reprod2));   
        B_q_DD(iyear+1)=sum(elem_prod(N(iyear+1)(set_q_DD_stage,nages),wgt_mt(set_q_DD_stage,nages))); 

  }

    //last year (projection) has no recruitment variability
  if(rec_alt1_switch==1){
	  N(iyear+1,1)=SR_func(R0, steep, spr_F0, SSB(endyr),SR_switch)*mfexp(rec_mean_alt1);
  }
  else{
	  N(endyr+1,1)=SR_func(R0, steep, spr_F0, SSB(endyr),SR_switch)*BiasCor;	  
  }
  N(endyr+1)(2,nages)=++elem_prod(N(endyr)(1,nages-1),(mfexp(-1.*Z(endyr)(1,nages-1))));
  N(endyr+1,nages)+=N(endyr,nages)*mfexp(-1.*Z(endyr,nages)); //plus group

  
FUNCTION get_landings_numbers //Baranov catch eqn
  for (iyear=styr; iyear<=endyr; iyear++)
  {
    for (iage=1; iage<=nages; iage++)
    {
      L_cHL_num(iyear,iage)=N(iyear,iage)*F_cHL(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_cLL_num(iyear,iage)=N(iyear,iage)*F_cLL(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);
      L_rGN_num(iyear,iage)=N(iyear,iage)*F_rGN(iyear,iage)*
        (1.-mfexp(-1.*Z(iyear,iage)))/Z(iyear,iage);        
    }          
    pred_cHL_L_knum(iyear)=sum(L_cHL_num(iyear))/1000.0;
    pred_cLL_L_knum(iyear)=sum(L_cLL_num(iyear))/1000.0;    
    pred_rGN_L_knum(iyear)=sum(L_rGN_num(iyear))/1000.0;
  }

 
FUNCTION get_landings_wgt
  for (iyear=styr; iyear<=endyr; iyear++)
  {    
    L_cHL_klb(iyear)=elem_prod(L_cHL_num(iyear),gutwgt_cHL_klb(iyear));     //in 1000 lb gutted weight
    L_cLL_klb(iyear)=elem_prod(L_cLL_num(iyear),gutwgt_cLL_klb(iyear));     //in 1000 lb gutted weight
    L_rGN_klb(iyear)=elem_prod(L_rGN_num(iyear),gutwgt_rGN_klb(iyear));     //in 1000 lb gutted weight
    
    pred_cHL_L_klb(iyear)=sum(L_cHL_klb(iyear));
    pred_cLL_L_klb(iyear)=sum(L_cLL_klb(iyear));
    pred_rGN_L_klb(iyear)=sum(L_rGN_klb(iyear));    
  }
 
   
FUNCTION get_catchability_fcns    
 //Get rate increase if estimated, otherwise fixed above
  if (set_q_rate_phase>0.0)
  {

      for (iyear=styr_cpue_cLL; iyear<=endyr_cpue_cLL; iyear++)
      {   if (iyear>styr_cpue_cLL & iyear <=2003) 
          {//q_rate_fcn_cLL(iyear)=(1.0+q_rate)*q_rate_fcn_cLL(iyear-1); //compound
             q_rate_fcn_cLL(iyear)=(1.0+(iyear-styr_cpue_cLL)*q_rate)*q_rate_fcn_cLL(styr_cpue_cLL);  //linear
          }
          if (iyear>2003) {q_rate_fcn_cLL(iyear)=q_rate_fcn_cLL(iyear-1);} 
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

 //cLL  cpue
  q_cLL(styr_cpue_cLL)=mfexp(log_q_cpue_cLL); 
  for (iyear=styr_cpue_cLL; iyear<=endyr_cpue_cLL; iyear++)
  {//index in weight units. original index in lb and re-scaled. predicted in klb, but difference is absorbed by q
      N_cLL(iyear)=elem_prod(elem_prod(N_mdyr(iyear),sel_cLL(iyear)),wholewgt_cLL_klb(iyear));   
	  //N_cLL(iyear)=elem_prod(elem_prod(N_mdyr(iyear),sel_cLL(iyear)),gutwgt_cLL_klb(iyear));   
      //pred_cLL_cpue(iyear)=q_cLL(iyear)*q_rate_fcn_cLL(iyear)*q_DD_fcn(iyear)*sum(N_cLL(iyear));
	  pred_cLL_cpue(iyear)=q_cLL(iyear)*q_rate_fcn_cLL(iyear)*sum(N_cLL(iyear));
      if (iyear<endyr_cpue_cLL){q_cLL(iyear+1)=q_cLL(iyear)*mfexp(q_RW_log_dev_cLL(iyear));}}

  q_sBL=mfexp(log_q_cpue_sBL);
  for (iyear=1; iyear<=nyr_cpue_sBL; iyear++)
   {N_sBL(iyear)=elem_prod(elem_prod(N_mdyr(yrs_cpue_sBL(iyear)),sel_sBL(yrs_cpue_sBL(iyear))),wgt_sBL_klb(yrs_cpue_sBL(iyear)));
     pred_sBL_cpue(iyear)=q_sBL*sum(N_sBL(iyear));}
  for (iyear=1; iyear<=nyr_cpue_sBL; iyear++)
    {pred_sBL_cpue_allyr(yrs_cpue_sBL(iyear))=pred_sBL_cpue(iyear);}  // for graphing purposes only
  
  for (iyear=1; iyear<=nyr_cpue_sBL; iyear++)
  {pred_sBL_cpue_allyr(yrs_cpue_sBL(iyear))=pred_sBL_cpue(iyear);}  // for graphing purposes only
  
FUNCTION get_length_comps


 //cGN handline

 //cGN longline

 //rGN
  for (iyear=1;iyear<=nyr_lenc_rGN;iyear++) 
  {pred_lenc_rGN(iyear)=(L_rGN_num(yrs_lenc_rGN(iyear))*lenprob_rGN)/sum(L_rGN_num(yrs_lenc_rGN(iyear)));}
  

FUNCTION get_age_comps
   

  //Commercial handline
  for (iyear=1;iyear<=nyr_agec_cHL;iyear++) 
  {
    ErrorFree_agec_cHL(iyear)=L_cHL_num(yrs_agec_cHL(iyear))/sum(L_cHL_num(yrs_agec_cHL(iyear)));  
    //ErrorFree_agec_cHL(iyear)=elem_prod(N(yrs_agec_cHL(iyear)),sel_cHL(yrs_agec_cHL(iyear)));
    pred_agec_cHL_allages(iyear)=age_error*(ErrorFree_agec_cHL(iyear)/sum(ErrorFree_agec_cHL(iyear)));   
    for (iage=1; iage<=nages_agec; iage++) {pred_agec_cHL(iyear,iage)=pred_agec_cHL_allages(iyear,iage);} 
    for (iage=(nages_agec+1); iage<=nages; iage++) {pred_agec_cHL(iyear,nages_agec)+=pred_agec_cHL_allages(iyear,iage);} //plus group                             
  }
 
  //Commercial longline
  for (iyear=1;iyear<=nyr_agec_cLL;iyear++) 
  {
    ErrorFree_agec_cLL(iyear)=L_cLL_num(yrs_agec_cLL(iyear))/sum(L_cLL_num(yrs_agec_cLL(iyear)));
    //ErrorFree_agec_cLL(iyear)=elem_prod(N(yrs_agec_cLL(iyear)),sel_cLL(yrs_agec_cLL(iyear)));
    pred_agec_cLL_allages(iyear)=age_error*(ErrorFree_agec_cLL(iyear)/sum(ErrorFree_agec_cLL(iyear)));     
    for (iage=1; iage<=nages_agec; iage++) {pred_agec_cLL(iyear,iage)=pred_agec_cLL_allages(iyear,iage);} 
    for (iage=(nages_agec+1); iage<=nages; iage++) {pred_agec_cLL(iyear,nages_agec)+=pred_agec_cLL_allages(iyear,iage);} //plus group                           
  }

 // sBL
  for (iyear=1;iyear<=nyr_agec_sBL;iyear++)
   {ErrorFree_agec_sBL(iyear)=N_sBL(iyear)/sum(N_sBL(iyear));
   //ErrorFree_agec_sBL(iyear)=elem_prod(N(yrs_agec_sBL(iyear)),sel_sBL(yrs_agec_sBL(iyear)));
   pred_agec_sBL_allages(iyear)=age_error*(ErrorFree_agec_sBL(iyear)/sum(ErrorFree_agec_sBL(iyear)));
   for (iage=1; iage<=nages_agec_sBL; iage++) {pred_agec_sBL(iyear,iage)=pred_agec_sBL_allages(iyear,iage);}
   for (iage=(nages_agec_sBL+1); iage<=nages; iage++) {pred_agec_sBL(iyear,nages_agec_sBL)+=pred_agec_sBL_allages(iyear,iage);}} //plus group
 
  
////--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_weighted_current 
  F_temp_sum=0.0;
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cHL+
        sum(log_dev_F_L_cHL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cLL+
        sum(log_dev_F_L_cLL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);  
  F_temp_sum+=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rGN+
        sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted);

  F_cHL_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cHL+
        sum(log_dev_F_L_cHL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_cLL_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_cLL+
        sum(log_dev_F_L_cLL((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  F_rGN_prop=mfexp((selpar_n_yrs_wgted*log_avg_F_L_rGN+
        sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr)))/selpar_n_yrs_wgted)/F_temp_sum;
  
  log_F_dev_end_cHL=sum(log_dev_F_L_cHL((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_cLL=sum(log_dev_F_L_cLL((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;
  log_F_dev_end_rGN=sum(log_dev_F_L_rGN((endyr-selpar_n_yrs_wgted+1),endyr))/selpar_n_yrs_wgted;  

  F_end_L=sel_cHL(endyr)*mfexp(log_avg_F_L_cHL+log_F_dev_end_cHL)+
          sel_cLL(endyr)*mfexp(log_avg_F_L_cLL+log_F_dev_end_cLL)+
          sel_rGN(endyr)*mfexp(log_avg_F_L_rGN+log_F_dev_end_rGN);  
    
  F_end=F_end_L;
  F_end_apex=max(F_end);
  
  sel_wgted_tot=F_end/F_end_apex;
  sel_wgted_L=elem_prod(sel_wgted_tot, elem_div(F_end_L,F_end));
  
  wgt_wgted_L_denom=F_cHL_prop+F_cLL_prop+F_rGN_prop;
  wgt_wgted_L_klb=F_cHL_prop/wgt_wgted_L_denom*gutwgt_cHL_klb(endyr)+
                  F_cLL_prop/wgt_wgted_L_denom*gutwgt_cLL_klb(endyr)+
                  F_rGN_prop/wgt_wgted_L_denom*gutwgt_rGN_klb(endyr);                          

				  
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
    B_eq(ff)=sum(elem_prod(N_age_msy,wgt_mt));
    L_eq_klb(ff)=sum(elem_prod(L_age_msy,wgt_wgted_L_klb)); //in gutted weight
    L_eq_knum(ff)=sum(L_age_msy)/1000.0;  
  }  
  
  msy_klb_out=max(L_eq_klb); //msy in gutted weight
  
  for(ff=1; ff<=n_iter_msy; ff++)
  {
   if(L_eq_klb(ff) == msy_klb_out) 
      {    
        SSB_msy_out=SSB_eq(ff);
        B_msy_out=B_eq(ff);
        R_msy_out=R_eq(ff);
        msy_knum_out=L_eq_knum(ff);
        F_msy_out=F_msy(ff);  
        spr_msy_out=spr_msy(ff);      
      }
  }

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
FUNCTION get_miscellaneous_stuff

//switch here if var_rec_dev <=dzero 
  if(var_rec_dev>0.0)
   {sigma_rec_dev=sqrt(var_rec_dev);} //pow(var_rec_dev,0.5);  //sample SD of predicted residuals (may not equal rec_sigma)  
   else{sigma_rec_dev=0.0;}

  len_cv=elem_div(len_sd,meanlen_TL);
  
  //compute total landings in 1000 fish and klb gutted weight
  L_total_num.initialize();
  L_total_klb.initialize();
  L_total_knum_yr.initialize();
  L_total_klb_yr.initialize();  
  
  for(iyear=styr; iyear<=endyr; iyear++)
  {
        L_total_klb_yr(iyear)=pred_cHL_L_klb(iyear)+pred_cLL_L_klb(iyear)+pred_rGN_L_klb(iyear);
        L_total_knum_yr(iyear)=pred_cHL_L_knum(iyear)+pred_cLL_L_knum(iyear)+pred_rGN_L_knum(iyear);
                
        B(iyear)=elem_prod(N(iyear),wgt_mt);
        totN(iyear)=sum(N(iyear));
        totB(iyear)=sum(B(iyear));                  
  }
  
  L_total_num=L_cHL_num+L_cLL_num+L_rGN_num;   //landings at age in number fish
  L_total_klb=L_cHL_klb+L_cLL_klb+L_rGN_klb;   //landings at age in klb gutted weight
 
  //Time series of interest

  B(endyr+1)=elem_prod(N(endyr+1),wgt_mt);
  totN(endyr+1)=sum(N(endyr+1));
  totB(endyr+1)=sum(B(endyr+1));  
  N_spawn(endyr+1)=N(endyr+1);
  SSB(endyr+1)=sum(elem_prod(N_spawn(endyr+1),reprod));
 // SSB(endyr+1)=sum(elem_prod(N_spawn(endyr+1),reprod));
  MatFemB(endyr+1)=sum(elem_prod(N_spawn(endyr+1),reprod2));
  rec=column(N,1);
  SdS0=SSB/S0;

  Fend_mean_temp=1.0;
  for (iyear=1; iyear<=selpar_n_yrs_wgted; iyear++) {Fend_mean_temp*=Fapex(endyr-iyear+1);}
  Fend_mean=pow(Fend_mean_temp,(1.0/selpar_n_yrs_wgted));	    
  if(F_msy_out>0)
    {
      FdF_msy=Fapex/F_msy_out;
      FdF_msy_end=FdF_msy(endyr);
      FdF_msy_end_mean=pow((FdF_msy(endyr)*FdF_msy(endyr-1)*FdF_msy(endyr-2)),(1.0/3.0));
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
                          (mfexp(-1.*(Z_age_spr(nages-1)*(1.0-spawn_time_frac) + Z_age_spr(nages)*spawn_time_frac) )))
                          /(1.0-mfexp(-1.*Z_age_spr(nages)));
    spr_spr(ff)=sum(elem_prod(N_age_spr_spawn,reprod));
    L_spr(ff)=0.0;
    for (iage=1; iage<=nages; iage++)
    {
      L_age_spr(iage)=N_age_spr(iage)*(F_L_age_spr(iage)/Z_age_spr(iage))*
                      (1.-mfexp(-1.*Z_age_spr(iage)));
      L_spr(ff)+=L_age_spr(iage)*wgt_wgted_L_klb(iage)*1000.0; //in lb gutted wgt
    }   
  }
  // Compute stuff for spr.brps
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
  //F_D_age_spr=F30_out*sel_wgted_D;
  Z_age_spr=M+F_L_age_spr;//+F_D_age_spr;

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
    //  D_age_F30(iage)=N_age_spr(iage)*(F_D_age_spr(iage)/Z_age_spr(iage))*
    //                  (1.-mfexp(-1.0*Z_age_spr(iage)));                      
    }
  SSB_F30_out=sum(elem_prod(N_age_spr_spawn,reprod));
  B_F30_out=sum(elem_prod(N_age_spr,wgt_mt));
  L_F30_klb_out=sum(elem_prod(L_age_F30,wgt_wgted_L_klb)); //in gutted weight
  L_F30_knum_out=sum(L_age_F30)/1000.0;  
  //D_F30_klb_out=sum(elem_prod(D_age_F30,wgt_wgted_D_klb)); //in whole weight   
  //D_F30_knum_out=sum(D_age_F30)/1000.0;      

//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------     
FUNCTION get_effective_sample_sizes
      neff_lenc_rGN_allyr_out=missing;
      
      neff_agec_cHL_allyr_out=missing;
      neff_agec_cLL_allyr_out=missing;
      neff_agec_sBL_allyr_out=missing;
        
		


      for (iyear=1; iyear<=nyr_lenc_rGN; iyear++)
         {if (nsamp_lenc_rGN(iyear)>=minSS_lenc_rGN)
            {neff_lenc_rGN_allyr_out(yrs_lenc_rGN(iyear))=multinom_eff_N(pred_lenc_rGN(iyear),obs_lenc_rGN(iyear));} 
		 else {neff_lenc_rGN_allyr_out(yrs_lenc_rGN(iyear))=-99;}}


      for (iyear=1; iyear<=nyr_agec_cHL; iyear++)
         {if (nsamp_agec_cHL(iyear)>=minSS_agec_cHL)
            {neff_agec_cHL_allyr_out(yrs_agec_cHL(iyear))=multinom_eff_N(pred_agec_cHL(iyear),obs_agec_cHL(iyear));}                            
		 else {neff_agec_cHL_allyr_out(yrs_agec_cHL(iyear))=-99;}}
         
      for (iyear=1; iyear<=nyr_agec_cLL; iyear++)
         {if (nsamp_agec_cLL(iyear)>=minSS_agec_cLL)
            {neff_agec_cLL_allyr_out(yrs_agec_cLL(iyear))=multinom_eff_N(pred_agec_cLL(iyear),obs_agec_cLL(iyear));}                            
		 else {neff_agec_cLL_allyr_out(yrs_agec_cLL(iyear))=-99;}}

      for (iyear=1; iyear<=nyr_agec_sBL; iyear++)
         {if (nsamp_agec_sBL(iyear)>=minSS_agec_sBL)
            {neff_agec_sBL_allyr_out(yrs_agec_sBL(iyear))=multinom_eff_N(pred_agec_sBL(iyear),obs_agec_sBL(iyear));}                            
          else {neff_agec_sBL_allyr_out(yrs_agec_sBL(iyear))=-99;}}
    
                           
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   


FUNCTION evaluate_objective_function
  //fval=square(xdum-9.0);
  
  fval=0.0;
  fval_data=0.0;  
//---likelihoods---------------------------

//---Indices-------------------------------

  f_cLL_cpue=0.0;
  f_cLL_cpue=lk_lognormal(pred_cLL_cpue, obs_cpue_cLL, obs_cv_cpue_cLL, w_cpue_cLL);
  fval+=f_cLL_cpue;
  fval_data+=f_cLL_cpue;  

  f_sBL_cpue=0.0;
  f_sBL_cpue=lk_lognormal(pred_sBL_cpue, obs_cpue_sBL, obs_cv_cpue_sBL, w_cpue_sBL);
  fval+=f_sBL_cpue;
  fval_data+=f_sBL_cpue;  

//---Landings-------------------------------
  
  //f_cHL_L in 1000 lb gutted wgt
  f_cHL_L=lk_lognormal(pred_cHL_L_klb(styr_L_cHL,endyr_L_cHL), obs_L_cHL(styr_L_cHL,endyr_L_cHL),
                      obs_cv_L_cHL(styr_L_cHL,endyr_L_cHL), w_L);
  fval+=f_cHL_L;
  fval_data+=f_cHL_L;

  //f_cLL_L in 1000 lb gutted wgt
  f_cLL_L=lk_lognormal(pred_cLL_L_klb(styr_L_cLL,endyr_L_cLL), obs_L_cLL(styr_L_cLL,endyr_L_cLL),
                      obs_cv_L_cLL(styr_L_cLL,endyr_L_cLL), w_L);
  fval+=f_cLL_L;
  fval_data+=f_cLL_L;

//f_rGN_L in 1000 lb gutted wgt
  f_rGN_L=lk_lognormal(pred_rGN_L_klb(styr_L_rGN,endyr_L_rGN), obs_L_rGN(styr_L_rGN,endyr_L_rGN), 
                      obs_cv_L_rGN(styr_L_rGN,endyr_L_rGN), w_L);
  fval+=f_rGN_L;
  fval_data+=f_rGN_L;  

  
//---Length comps-------------------------------


 
  
  //f_lenc_rGN
  f_lenc_rGN=lk_dirichlet_multinomial(nsamp_lenc_rGN, pred_lenc_rGN, obs_lenc_rGN, nyr_lenc_rGN, double(nlenbins), minSS_lenc_rGN, log_dm_lenc_rGN);
  //f_lenc_rGN=lk_robust_multinomial(nsamp_lenc_rGN, pred_lenc_rGN, obs_lenc_rGN, nyr_lenc_rGN, double(nlenbins), minSS_lenc_rGN, w_lenc_rGN);
  //f_lenc_rGN=lk_multinomial(nsamp_lenc_rGN, pred_lenc_rGN, obs_lenc_rGN, nyr_lenc_rGN, minSS_lenc_rGN, w_lenc_rGN);
  fval+=set_w_lenc_rGN*f_lenc_rGN;
  fval_data+=f_lenc_rGN;
  
 
//---Age comps-------------------------------

  //f_agec_cHL
  f_agec_cHL=lk_dirichlet_multinomial(nsamp_agec_cHL, pred_agec_cHL, obs_agec_cHL, nyr_agec_cHL, double(nages_agec), minSS_agec_cHL, log_dm_agec_cHL);
  //f_agec_cHL=lk_robust_multinomial(nsamp_agec_cHL, pred_agec_cHL, obs_agec_cHL, nyr_agec_cHL, double(nages_agec), minSS_agec_cHL, w_agec_cHL);
  //f_agec_cHL=lk_multinomial(nsamp_agec_cHL, pred_agec_cHL, obs_agec_cHL, nyr_agec_cHL, minSS_agec_cHL, w_agec_cHL);
  fval+=set_w_agec_cHL*f_agec_cHL;
  fval_data+=f_agec_cHL;

  //f_agec_cLL
  f_agec_cLL=lk_dirichlet_multinomial(nsamp_agec_cLL, pred_agec_cLL, obs_agec_cLL, nyr_agec_cLL, double(nages_agec), minSS_agec_cLL, log_dm_agec_cLL);
  //f_agec_cLL=lk_robust_multinomial(nsamp_agec_cLL, pred_agec_cLL, obs_agec_cLL, nyr_agec_cLL, double(nages_agec), minSS_agec_cLL, w_agec_cLL);
  //f_agec_cLL=lk_multinomial(nsamp_agec_cLL, pred_agec_cLL, obs_agec_cLL, nyr_agec_cLL, minSS_agec_cLL, w_agec_cLL);
  fval+=set_w_agec_cLL*f_agec_cLL;
  fval_data+=f_agec_cLL;

  //f_agec_sBL
  f_agec_sBL=lk_dirichlet_multinomial(nsamp_agec_sBL, pred_agec_sBL, obs_agec_sBL, nyr_agec_sBL, double(nages_agec_sBL), minSS_agec_sBL, log_dm_agec_sBL);
  //f_agec_sBL=lk_robust_multinomial(nsamp_agec_sBL, pred_agec_sBL, obs_agec_sBL, nyr_agec_sBL, double(nages_agec_sBL), minSS_agec_sBL, w_agec_sBL);
  //f_agec_sBL=lk_multinomial(nsamp_agec_sBL, pred_agec_sBL, obs_agec_sBL, nyr_agec_sBL, minSS_agec_sBL, w_agec_sBL);
  fval+=set_w_agec_sBL*f_agec_sBL;
  fval_data+=f_agec_sBL;

//-----------Constraints and penalties--------------------------------
  
  //Light penalty applied to log_dev_Nage for deviation from zero. If not estimated, this penalty equals zero.
  f_Nage_init=0.0;
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
  
//  //Random walk components of fishery dependent indices
//  f_rHB_RW_cpue=0.0;
//  for (iyear=styr_rHB_cpue; iyear<endyr_rHB_cpue; iyear++)
//      {f_rHB_RW_cpue+=square(q_RW_log_dev_rHB(iyear))/(2.0*set_q_RW_rHB_var);}
//  fval+=f_rHB_RW_cpue;   

  
//---Priors---------------------------------------------------
//neg_log_prior arguments: estimate, prior mean, prior var/-CV, pdf type
//Variance input as a negative value is considered to be CV in arithmetic space (CV=-1 implies loose prior) 
//pdf type 1=none, 2=lognormal, 3=normal, 4=beta 
  f_priors=0.0; 
  //f_priors+=neg_log_prior(Linf,set_Linf(5),set_Linf(6),set_Linf(7));
  //f_priors+=neg_log_prior(K,set_K(5),set_K(6),set_K(7));
  //f_priors+=neg_log_prior(t0,set_t0(5),set_t0(6),set_t0(7));
  f_priors+=neg_log_prior(len_cv_val,set_len_cv(5),set_len_cv(6),set_len_cv(7));
  //f_priors+=neg_log_prior(M_constant,set_M_constant(5),set_M_constant(6),set_M_constant(7));
   
  //f_priors+=neg_log_prior(steep,set_steep(5),set_log_R0(6),set_log_R0(7)); 
  //f_priors+=neg_log_prior(log_R0,set_log_R0(5),set_log_R0(6),set_log_R0(7)); 
  //f_priors+=neg_log_prior(R_autocorr,set_R_autocorr(5),set_R_autocorr(6),set_R_autocorr(7));
  f_priors+=neg_log_prior(rec_sigma,set_rec_sigma(5),set_rec_sigma(6),set_rec_sigma(7));
 
  f_priors+=neg_log_prior(selpar_L50_cHL,set_selpar_L50_cHL(5), set_selpar_L50_cHL(6), set_selpar_L50_cHL(7));
  f_priors+=neg_log_prior(selpar_L50_cHL2,set_selpar_L50_cHL2(5), set_selpar_L50_cHL2(6), set_selpar_L50_cHL2(7));
  f_priors+=neg_log_prior(selpar_slope_cHL,set_selpar_slope_cHL(5), set_selpar_slope_cHL(6), set_selpar_slope_cHL(7));
  f_priors+=neg_log_prior(selpar_slope_cHL2,set_selpar_slope_cHL2(5), set_selpar_slope_cHL2(6), set_selpar_slope_cHL2(7));
  f_priors+=neg_log_prior(selpar_L50_cLL,set_selpar_L50_cLL(5), set_selpar_L50_cLL(6), set_selpar_L50_cLL(7));
  f_priors+=neg_log_prior(selpar_L50_cLL2,set_selpar_L50_cLL2(5), set_selpar_L50_cLL2(6), set_selpar_L50_cLL2(7));
  f_priors+=neg_log_prior(selpar_slope_cLL,set_selpar_slope_cLL(5), set_selpar_slope_cLL(6), set_selpar_slope_cLL(7));
  f_priors+=neg_log_prior(selpar_slope_cLL2,set_selpar_slope_cLL2(5), set_selpar_slope_cLL2(6), set_selpar_slope_cLL2(7));
  f_priors+=neg_log_prior(selpar_L50_rGN,set_selpar_L50_rGN(5), set_selpar_L50_rGN(6), set_selpar_L50_rGN(7));
  f_priors+=neg_log_prior(selpar_slope_rGN,set_selpar_slope_rGN(5), set_selpar_slope_rGN(6), set_selpar_slope_rGN(7));

  f_priors+=neg_log_prior(selpar_L50_sBL,set_selpar_L50_sBL(5), set_selpar_L50_sBL(6), set_selpar_L50_sBL(7));
  f_priors+=neg_log_prior(selpar_slope_sBL,set_selpar_slope_sBL(5), set_selpar_slope_sBL(6), set_selpar_slope_sBL(7));

  //f_priors+=neg_log_prior(selpar_afull_cLL,set_selpar_afull_cLL(5), set_selpar_afull_cLL(6), set_selpar_afull_cLL(7));
  //f_priors+=neg_log_prior(selpar_sigma_cLL,set_selpar_afull_cLL(5), set_selpar_afull_cLL(6), set_selpar_afull_cLL(7));
  
  //f_priors+=neg_log_prior(log_q_cpue_cLL,set_log_q_cpue_cLL(5),set_log_q_cpue_cLL(6),set_log_q_cpue_cLL(7));
  //f_priors+=neg_log_prior(log_q_cpue_sBL,set_log_q_cpue_sBL(5),set_log_q_cpue_sBL(6),set_log_q_cpue_sBL(7));
  
  //f_priors+=neg_log_prior(F_init,set_F_init(5),set_F_init(6),set_F_init(7));
//  f_priors+=neg_log_prior(log_avg_F_L_cHL,set_log_avg_F_L_cHL(5),set_log_avg_F_L_cHL(6),set_log_avg_F_L_cHL(7));
//  f_priors+=neg_log_prior(log_avg_F_L_cLL,set_log_avg_F_L_cLL(5),set_log_avg_F_L_cLL(6),set_log_avg_F_L_cLL(7));
//  f_priors+=neg_log_prior(log_avg_F_rHB,set_log_avg_F_rHB(5),set_log_avg_F_rHB(6),set_log_avg_F_rHB(7));
//  f_priors+=neg_log_prior(log_avg_F_L_rGN,set_log_avg_F_L_rGN(5),set_log_avg_F_L_rGN(6),set_log_avg_F_L_rGN(7));
 
  f_priors+=neg_log_prior(log_dm_lenc_rGN,set_log_dm_lenc_rGN(5),set_log_dm_lenc_rGN(6),set_log_dm_lenc_rGN(7));
  f_priors+=neg_log_prior(log_dm_agec_cHL,set_log_dm_agec_cHL(5),set_log_dm_agec_cHL(6),set_log_dm_agec_cHL(7));
  f_priors+=neg_log_prior(log_dm_agec_cLL,set_log_dm_agec_cLL(5),set_log_dm_agec_cLL(6),set_log_dm_agec_cLL(7));
  f_priors+=neg_log_prior(log_dm_agec_sBL,set_log_dm_agec_sBL(5),set_log_dm_agec_sBL(6),set_log_dm_agec_sBL(7));

 
  fval+=f_priors;
  fval=fval/1.0;
  
  //cout << "fval = " << fval << "  fval_data = " << fval_data << endl;
  //cout << endl;

//----------------------------------------------------------------------------------
//Logistic function: 2 parameters
FUNCTION dvar_vector logistic(const dvar_vector& ages, const dvariable& L50, const dvariable& slope)
  //ages=vector of ages, L50=age at 50% selectivity, slope=rate of increase
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1./(1.+mfexp(-1.*slope*(ages-L50))); //logistic;  
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;

//-----------------------------------------------------------------------------------                                                                           
//Logistic-exponential: 4 parameters (but 1 is fixed)
FUNCTION dvar_vector logistic_exponential(const dvar_vector& ages, const dvariable& L50, const dvariable& slope, const dvariable& sigma, const dvariable& joint)
  //ages=vector of ages, L50=age at 50% sel (ascending limb), slope=rate of increase, sigma=controls rate of descent (descending)                               
  //joint=age to join curves                                                                                                                                    
  RETURN_ARRAYS_INCREMENT();                                                                                                                                    
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());                                                                                                         
  Sel_Tmp=1.0;                                                                                                                                                  
  for (iage=1; iage<=nages; iage++)                                                                                                                             
  {                                                                                                                                                             
   if (ages(iage)<joint) {Sel_Tmp(iage)=1./(1.+mfexp(-1.*slope*(ages(iage)-L50)));}                                                                             
   if (ages(iage)>joint){Sel_Tmp(iage)=mfexp(-1.*square((ages(iage)-joint)/sigma));}                                                                            
  }                                                                                                                                                             
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);                                                                                                                                 
  RETURN_ARRAYS_DECREMENT();                                                                                                                                    
  return Sel_Tmp;   

//-----------------------------------------------------------------------------------
//Logistic function: 4 parameters
FUNCTION dvar_vector logistic_double(const dvar_vector& ages, const dvariable& L501, const dvariable& slope1, const dvariable& L502, const dvariable& slope2)
  //ages=vector of ages, L50=age at 50% selectivity, slope=rate of increase, L502=age at 50% decrease additive to L501, slope2=slope of decrease
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=elem_prod( (1./(1.+mfexp(-1.*slope1*(ages-L501)))),(1.-(1./(1.+mfexp(-1.*slope2*(ages-(L501+L502)))))) );     
  Sel_Tmp=Sel_Tmp/max(Sel_Tmp);
  RETURN_ARRAYS_DECREMENT();
  return Sel_Tmp;

//-----------------------------------------------------------------------------------
//Jointed logistic function: 6 parameters (increasing and decreasing logistics joined at peak selectivity)
FUNCTION dvar_vector logistic_joint(const dvar_vector& ages, const dvariable& L501, const dvariable& slope1, const dvariable& L502, const dvariable& slope2, const dvariable& satval, const dvariable& joint)
  //ages=vector of ages, L501=age at 50% sel (ascending limb), slope1=rate of increase,L502=age at 50% sel (descending), slope1=rate of increase (ascending), 
  //satval=saturation value of descending limb, joint=location in age vector to join curves (may equal age or age + 1 if age-0 is included)
  RETURN_ARRAYS_INCREMENT();
  dvar_vector Sel_Tmp(ages.indexmin(),ages.indexmax());
  Sel_Tmp=1.0; 
  for (iage=1; iage<=nages; iage++)
  {
   if (double(iage)<joint) {Sel_Tmp(iage)=1./(1.+mfexp(-1.*slope1*(ages(iage)-L501)));}  
   if (double(iage)>joint){Sel_Tmp(iage)=1.0-(1.0-satval)/(1.+mfexp(-1.*slope2*(ages(iage)-L502)));}  
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
  dvariable small_number=0.00001;
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


//-----------------------------------------------------------------------------------
//Likelihood contribution: multinomial
FUNCTION dvariable lk_robust_multinomial(const dvar_vector& nsamp, const dvar_matrix& pred_comp, const dvar_matrix& obs_comp, const double& ncomp, const dvariable& mbin, const double& minSS, const dvariable& wgt_dat)
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
      //cout<<"start report"<<endl;
      get_weighted_current();
      //cout<<"got weighted"<<endl;
      get_msy();
      //cout<<"got msy"<<endl;
      get_per_recruit_stuff();
      //cout<<"got per recruit"<<endl;  	  
      get_miscellaneous_stuff();
      //cout<<"got misc stuff"<<endl;
      get_effective_sample_sizes();
     
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
      //cout << "len_cv = "<<len_cv_val<<endl;
      //cout << "xdum " << xdum << endl;
      cout << "><>--><>--><>--><>--><>--><>--><>--><>--><>--><>"  <<endl;  
      // cout << F_initial << endl;
      
      report << "TotalLikelihood " << fval << endl;
      report << "N" << endl;
      report << N<<endl;
      report << "F" << endl;
      report << F <<endl;
	  report << obs_cpue_sBL_allyr << endl;
	  report << pred_sBL_cpue_allyr << endl;
	  report << obs_cv_cpue_sBL_allyr << endl;
	  report << obs_cpue_cLL << endl;
	  report << pred_cLL_cpue << endl;
	  report << obs_cv_cpue_cLL << endl;
	  
//        report <<"lenprob" <<endl;
//        report << lenprob<<endl;      

	  sdnr_lc_rGN=sdnr_multinomial(nyr_lenc_rGN, lenbins, nsamp_lenc_rGN, pred_lenc_rGN, obs_lenc_rGN, w_lenc_rGN); 
	  
      sdnr_ac_cHL=sdnr_multinomial(nyr_agec_cHL, agebins_agec, nsamp_agec_cHL, pred_agec_cHL, obs_agec_cHL, w_agec_cHL);  
      sdnr_ac_cLL=sdnr_multinomial(nyr_agec_cLL, agebins_agec, nsamp_agec_cLL, pred_agec_cLL, obs_agec_cLL, w_agec_cLL);  
      sdnr_ac_sBL=sdnr_multinomial(nyr_agec_sBL, agebins_agec_sBL, nsamp_agec_sBL, pred_agec_sBL, obs_agec_sBL, w_agec_sBL);  
       
      sdnr_I_cLL=sdnr_lognormal(pred_cLL_cpue, obs_cpue_cLL, obs_cv_cpue_cLL, w_cpue_cLL);
      sdnr_I_sBL=sdnr_lognormal(pred_sBL_cpue, obs_cpue_sBL, obs_cv_cpue_sBL, w_cpue_sBL);
      
	  	  

      
      //#################################################################################################
      //##  Passing parameters to vector for bounds check plotting
      //################################################################################################# 
       Linf_out(8)=Linf; Linf_out(1,7)=set_Linf; 
       K_out(8)=K; K_out(1,7)=set_K;
       t0_out(8)=t0; t0_out(1,7)=set_t0;
       Linf_f_out(8)=Linf_f; Linf_f_out(1,7)=set_Linf_f; 
       K_f_out(8)=K_f; K_f_out(1,7)=set_K_f;
       t0_f_out(8)=t0_f; t0_f_out(1,7)=set_t0_f;
       len_cv_val_out(8)=len_cv_val; len_cv_val_out(1,7)=set_len_cv;
       log_R0_out(8)=log_R0; log_R0_out(1,7)=set_log_R0;
       M_constant_out(8)=M_constant; M_constant_out(1,7)=set_M_constant;
       steep_out(8)=steep; steep_out(1,7)=set_steep;
       rec_sigma_out(8)=rec_sigma; rec_sigma_out(1,7)=set_rec_sigma;
       R_autocorr_out(8)=R_autocorr; R_autocorr_out(1,7)=set_R_autocorr;

	   log_dm_lenc_rGN_out(8)=log_dm_lenc_rGN; log_dm_lenc_rGN_out(1,7)=set_log_dm_lenc_rGN;
	   log_dm_agec_cHL_out(8)=log_dm_agec_cHL; log_dm_agec_cHL_out(1,7)=set_log_dm_agec_cHL;
	   log_dm_agec_cLL_out(8)=log_dm_agec_cLL; log_dm_agec_cLL_out(1,7)=set_log_dm_agec_cLL;
	   log_dm_agec_sBL_out(8)=log_dm_agec_sBL; log_dm_agec_sBL_out(1,7)=set_log_dm_agec_sBL;
       
       selpar_L50_cHL_out(8)=selpar_L50_cHL; selpar_L50_cHL_out(1,7)=set_selpar_L50_cHL;
       selpar_L50_cHL2_out(8)=selpar_L50_cHL2; selpar_L50_cHL2_out(1,7)=set_selpar_L50_cHL2;
       selpar_slope_cHL_out(8)=selpar_slope_cHL; selpar_slope_cHL_out(1,7)=set_selpar_slope_cHL;
       selpar_slope_cHL2_out(8)=selpar_slope_cHL2; selpar_slope_cHL2_out(1,7)=set_selpar_slope_cHL2;

       selpar_L50_cLL_out(8)=selpar_L50_cLL; selpar_L50_cLL_out(1,7)=set_selpar_L50_cLL;
       selpar_L50_cLL2_out(8)=selpar_L50_cLL2; selpar_L50_cLL2_out(1,7)=set_selpar_L50_cLL2;
       selpar_slope_cLL_out(8)=selpar_slope_cLL; selpar_slope_cLL_out(1,7)=set_selpar_slope_cLL;
       selpar_slope_cLL2_out(8)=selpar_slope_cLL2; selpar_slope_cLL2_out(1,7)=set_selpar_slope_cLL2;
       //selpar_afull_cLL_out(8)=selpar_afull_cLL; selpar_afull_cLL_out(1,7)=set_selpar_afull_cLL;
       //selpar_sigma_cLL_out(8)=selpar_sigma_cLL; selpar_sigma_cLL_out(1,7)=set_selpar_afull_cLL;

       selpar_L50_rGN_out(8)=selpar_L50_rGN; selpar_L50_rGN_out(1,7)=set_selpar_L50_rGN;
       selpar_slope_rGN_out(8)=selpar_slope_rGN; selpar_slope_rGN_out(1,7)=set_selpar_slope_rGN;
	   
	   selpar_L50_sBL_out(8)=selpar_L50_sBL; selpar_L50_sBL_out(1,7)=set_selpar_L50_sBL;
       selpar_slope_sBL_out(8)=selpar_slope_sBL; selpar_slope_sBL_out(1,7)=set_selpar_slope_sBL;

       log_q_cpue_cLL_out(8)=log_q_cpue_cLL; log_q_cpue_cLL_out(1,7)=set_log_q_cpue_cLL;
       log_q_cpue_sBL_out(8)=log_q_cpue_sBL; log_q_cpue_sBL_out(1,7)=set_log_q_cpue_sBL;
                     
       log_avg_F_L_cHL_out(8)=log_avg_F_L_cHL; log_avg_F_L_cHL_out(1,7)=set_log_avg_F_L_cHL;
       log_avg_F_L_cLL_out(8)=log_avg_F_L_cLL; log_avg_F_L_cLL_out(1,7)=set_log_avg_F_L_cLL;
       log_avg_F_L_rGN_out(8)=log_avg_F_L_rGN; log_avg_F_L_rGN_out(1,7)=set_log_avg_F_L_rGN;       
       F_init_out(8)=F_init; F_init_out(1,7)=set_F_init;
       
       log_rec_dev_out(styr_rec_dev, endyr_rec_dev)=log_dev_rec;
       log_F_dev_L_cHL_out(styr_L_cHL,endyr_L_cHL)=log_dev_F_L_cHL;
       log_F_dev_L_cLL_out(styr_L_cLL,endyr_L_cLL)=log_dev_F_L_cLL;
       log_F_dev_L_rGN_out(styr_L_rGN,endyr_L_rGN)=log_dev_F_L_rGN;
           
    #include "bam.cxx"   // write the R-compatible report
	 //cout<<"All done!"<<endl;
	 //system("type beep.txt"); 
  
  } //endl last phase loop     

