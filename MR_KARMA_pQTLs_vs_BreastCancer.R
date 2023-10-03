#### all pQTLs - IVs against cancer
library(TwoSampleMR)
library(data.table)

pqtl_file=read.table('data_exposure_mr.txt',sep='\t',header=T)
outcome=read.table('data_outcome_mr.txt',sep='\t',header=T)
dir.create('resultMR')
tmpDir=getwd()

uProt=unique(pqtl_file$ProteinID)

for(i in 1:length(uProt)){
  pqtl_file_tmp=pqtl_file[which(pqtl_file$ProteinID == uProt[i]),]
  write.table(pqtl_file_tmp,file='pqtl_tmp.txt',sep='\t',row.names=F,quote=F)
  pqtl_in='pqtl_tmp.txt'
  pqtl_exp_data_in=read_exposure_data(pqtl_in,snp_col= 'SNP', beta_col='beta',se_col='se',
                                      eaf_col='eaf', effect_allele_col = "effect_allele",
                                      other_allele_col = "other_allele",
                                      pval_col = "pval",
                                      samplesize_col = "samplesize",
                                      gene_col = "ProteinID", clump=F, sep='\t')
  outcome_dat_in <- read_outcome_data(
    snps = pqtl_in_exp_data$ID,
    filename = "data_outcome_mr.txt",
    sep = "\t",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = 'eaf',
    pval_col = 'pval',
    ncase_col = "ncases",
    ncontrol_col="ncontrols"
  )
  dat <- harmonise_data(
    exposure_dat = pqtl_exp_data_in, 
    outcome_dat = outcome_dat_in
  )
  if(sum(dat$mr_keep == 'TRUE') < 1){
    next()
  }
  dat$exposure = paste(uProt[i],'_cis_IVs',sep='')
  dat$outcome = 'BREASTCANCER BCAC'
  res=mr(dat)
  res_single <- mr_singlesnp(dat)
  write.table(res,file=paste(tmpDir,'/resultMR/MR_',uProt[i],'_KARMA_3K_vs_BCAC.txt',sep=''),sep='\t',row.names=F)
  write.table(res_single,file=paste(tmpDir,'/resultMR/MR_',uProt[i],'_KARMA_3K_cis_vs_BCAC_singlesnps.txt',sep=''),sep='\t',row.names=F)
  rm(dat)
  rm(res)
  rm(res_single)
}