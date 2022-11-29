library(bamExtras)

CommonName = "BlackSeaBass";
fileName     = "bam";
dir_bam_sim  = "BlSB_sim";
dir_bam_base = "BlSB_base";
bam=NULL;
dat_file=NULL;
tpl_file=NULL;
cxx_file=NULL;
dat_obj=NULL;
tpl_obj=NULL;
cxx_obj=NULL;
standardize=FALSE;
nyr_remove=1:5;
# parallel=TRUE, # Right now it has to be in parallel
coresUse=NULL;
ndigits=4;
unlink_dir_bam_base=FALSE;
run_bam_base=TRUE;
overwrite_bam_base=TRUE;
admb_switch_base = '-nox';
run_sim=TRUE;
admb_switch_sim = '-est -nox -ind';
prompt_me=FALSE;
subset_rdat=list("eq.series"=101,"pr.series"=101);
random_seed=12345;
admb2r_obj = admb2r.cpp;
cleanup = list(del=c("*.r0*","*.p0*","*.b0*","*.log","*.rpt","*.obj",
                     "*.htp","*.eva","*.bar","*.tds","*.o","tmp_admb",
                     "variance","*.dep","*.hes","*.tmp"))

