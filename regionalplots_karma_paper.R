library(RACER) 
pqtl_karma=as.data.frame(fread('KARMA_pQTL_DNPH1_cisregion.txt'))
bc_cisregion=as.data.frame(fread('BreastCancer_BCAC_DNPH1_cisregion_b37.txt'))

gene='DNPH1'
genechr=6
topsnp=pqtl_karma$SNP[pqtl_karma$pvalue == min(pqtl_karma$pvalue)]

### LD regional plots of pQTL and breast cancer gwas
pqtl_karma_f = RACER::formatRACER(assoc_data = pqtl_karma, chr_col = 2, pos_col = 3, p_col = 10)
bc_cisregion_f = RACER::formatRACER(assoc_data = bc_cisregion, chr_col = 2, pos_col = 3, p_col = 9)

pqtl_karma_f_ld = RACER::ldRACER(assoc_data = pqtl_karma_f, rs_col = 1, pops = "EUR", lead_snp = topsnp)
bc_cisregion_f_ld = RACER::ldRACER(assoc_data = bc_cisregion_f, rs_col = 1, pops = "EUR", lead_snp = topsnp)

p=mirrorPlotRACER(assoc_data1 = pqtl_karma_f_ld, assoc_data2 = bc_cisregion_f_ld, chr = genechr,  build = "hg19",plotby = "gene", gene_plot = gene,name1=paste(gene,'pQTL',sep=' '),name2="Breast cancer BCAC")

png(paste('mirrorplot_',gene,'_pqtl_vs_breastcancer.png',sep=''))
print(p)
dev.off()
pdf(paste('mirrorplot_',gene,'_pqtl_vs_breastcancer.pdf',sep=''))
print(p)
dev.off()