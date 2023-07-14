library(bamExtras)

CommonName = NULL
fileName     = "bam"
dir_bam_sim  = "sim"
dir_bam_base = "base"
bam=bam=bam2r_out
dat_file=NULL
tpl_file=NULL
cxx_file=NULL
dat_obj=NULL
tpl_obj=NULL
cxx_obj=NULL
standardize=TRUE
nyr_remove=1:5
# parallel=TRUE, # Right now it has to be in parallel
ncores=NULL
ndigits=4
unlink_dir_bam_base=FALSE
run_bam_base=TRUE
overwrite_bam_base=TRUE
admb_options_base = '-nox'
run_sim=TRUE
admb_options_sim = '-est -nox -ind'
prompt_me=FALSE
subset_rdat=list("eq.series"=101,"pr.series"=101)
random_seed=12345
admb2r_obj = admb2r.cpp
cleanup = list(del=c("*.r0*","*.p0*","*.b0*","*.log","*.rpt","*.obj",
                     "*.htp","*.eva","*.bar","*.tds","*.o","tmp_admb",
                     "variance","*.dep","*.hes","*.tmp"))

# writeLines(dat_RedGrouper,"RedGrouper.dat")
# writeLines(tpl_RedGrouper,"RedGrouper.tpl")
# writeLines(cxx_RedGrouper,"RedGrouper.cxx")
#
# redg2 <- bam2r(dat_file = "RedGrouper.dat",tpl_file = "RedGrouper.tpl",cxx_file = "RedGrouper.cxx")
# #
# # run_retro(bam=redg2)
