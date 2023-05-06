rm(list=ls())

CommonName = "ScampGrouper"
fileName="bam"
dir_bam_sim  = "ScmG_sim"
dir_bam_base = "ScmG_base"
bam=NULL
dat_file=NULL
tpl_file=NULL
cxx_file=NULL
dat_obj=NULL
tpl_obj=NULL
cxx_obj=NULL
data_sim =  list(cv_U=NULL,cv_L=NULL,cv_D=NULL)
par_default = list(cv_U=0.2,cv_L=0.2,cv_D=0.2)
standardize=TRUE
nsim=10
sclim_gen = c(0.9,1.1)
sclim = list()
data_type_resamp = c()#c("U","L","D","age","len")
fn_par = list(
  M = "runif(nsim,M_min,M_max)",
  K = "runif(nsim,K_min,K_max)",
  Linf = "runif(nsim,Linf_min,Linf_max)",
  t0 = "runif(nsim,t0_min,t0_max)",
  steep = "runif(nsim,steep_min,steep_max)",
  rec_sigma = "rtnorm(n=nsim,mean=0.6,sd=0.15,lower=0.3,upper=1.0)",
  Dmort = "apply(Dmort_lim,2,function(x){runif(nsim,min(x),max(x))})",
  Pfa=    "runif(nsim,Pfa_min,Pfa_max)",     # these are scalars, not actual parameter values
  Pfb=    "runif(nsim,Pfb_min,Pfb_max)",     # these are scalars, not actual parameter values
  Pfma=   "runif(nsim,Pfma_min,Pfma_max)",  # these are scalars, not actual parameter values
  Pfmb=   "runif(nsim,Pfmb_min,Pfmb_max)"   # these are scalars, not actual parameter values
)
fix_par = c()
# fn_par = list(
#   steep = expression(seq(0.21,0.99,length=nsim))
# )
# fix_par = "set_steep"
# parallel=TRUE, # Right now it has to be in parallel
coresUse=NULL
ndigits=4
unlink_dir_bam_base=FALSE
run_bam_base=TRUE
overwrite_bam_base=TRUE
admb_switch_base = '-nox'
run_sim=TRUE
admb_switch_sim = '-est -nox -ind'
prompt_me=FALSE
subset_rdat=list("eq.series"=101,"pr.series"=101)
random_seed=12345
