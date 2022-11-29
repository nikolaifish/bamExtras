CommonName = "ScampGrouper"
fileName="bam"
dir_bam_sim  = "ScGr_sim"
dir_bam_base = "ScGr_base"
bam=NULL
dat_file=NULL
tpl_file=NULL
cxx_file=NULL
dat_obj=NULL
tpl_obj=NULL
cxx_obj=NULL
standardize=TRUE
nsim=10
M_FUN="M_Lorenzen"
M_FUN_args=list(aP=0)
sc = 0.1
scLim = sc*c(-1,1)+1
Mc_scLim = 0.1*c(-1,1)+1
steep_scLim = c(1,1)
Linf_scLim = scLim
K_scLim = scLim
t0_scLim = scLim
Dmort_scLim = scLim
coresUse=NULL
ndigits=4 # number of digits to round simulated values to
unlink_dir_bam_base=FALSE
run_bam_base=TRUE # If FALSE, the function will look for an executable named fileName.exe in dir_bam_base and use it as the base model.
# If TRUE and overwrite_bam_base=TRUE, the function will call run_bam.
overwrite_bam_base=FALSE
admb_switch_sim = '-est -nox -ind'
prompt_me=FALSE
run_sim=TRUE # If FALSE, the simulated data will be generated but won't be used in new BAM runs
admb_switch_base = '-nox'
subset_rdat=list("eq.series"=101,"pr.series"=101)

data_sim =  list(cv_U=NULL,
                 cv_L=NULL,
                 cv_D=NULL)

par_default = list(cv_U=0.2,
                   cv_L=0.2,
                   cv_D=0.2)

# data_sim = list(cv_D = read.csv("cv_D_BlSB.csv",row.names = 1),
#                  cv_L = read.csv("cv_L_BlSB.csv",row.names = 1),
#                  cv_U = read.csv("cv_U_BlSB.csv",row.names = 1)
# )
