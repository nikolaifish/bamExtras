CommonName = "AtlanticMenhaden"
fileName="bam"
dir_bam_sim  = "AtMe_sim"
dir_bam_base = "AtMe_base"
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
# fn_par = list(
#   M = expression(runif(nsim,min(M_lim),max(M_lim))),
#   K = expression(runif(nsim,min(K_lim),max(K_lim))),
#   Linf = expression(runif(nsim,min(Linf_lim),max(Linf_lim))),
#   t0 = expression(runif(nsim,min(t0_lim),max(t0_lim))),
#   steep = expression(runif(nsim,min(steep_lim),max(steep_lim))),
#   rec_sigma = expression(rtnorm(n=nsim,mean=0.6,sd=0.15,lower=0.3,upper=1.0)),
#   Dmort = expression(apply(Dmort_lim,2,function(x){runif(nsim,min(x),max(x))}))
# )
# fix_par = c()
fn_par = list(
  steep = expression(seq(0.21,0.99,length=nsim))
)
fix_par = "set_steep"
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
