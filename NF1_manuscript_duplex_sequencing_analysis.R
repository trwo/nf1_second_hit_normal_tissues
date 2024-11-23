### Henry Lee-Six, 2024
# analyse targeted nanorate sequencing calls in the NF1 gene.

library(dndscv)
library(RColorBrewer)

# read in the data
tns = read.csv('neurofibromatosis_study_tns_output.txt', sep='\t', header=T, stringsAsFactors = F)

# read in the gene list. These are the genes in the baitset
gene_list = as.vector(read.csv('gene_list_for_dnds_analysis.txt', header=F)$V1)

# read in the duplex coverage statistics.
duplex_cov = read.csv('duplex_coverage_stats.txt', sep='\t', header = F, stringsAsFactors = F)
dc = duplex_cov$V2
names(dc) = duplex_cov$V1

# get a measure of how highly each gene is expressed relative to the median
duplex_cov$rel_to_median = duplex_cov$V2/median(duplex_cov$V2)
rownames(duplex_cov) = duplex_cov$V1

##################################
## analyse by germ layre, separating out spleen ##

# cohort 1 NF1 wt (PD50297, PD51123) vs NF1 mutant (PD51122), neuroectoderm vs mesoderm vs analysed separately. 
nf1_wt_neuroectoderm_data = tns[tns$patient %in% c('PD50297', 'PD51123') & tns$Embryology %in% c('CNS', 'PNS', 'Ectoderm / skin'), c('sampleID', 'chr', 'pos', 'ref', 'mut', 'gene', 'duplex_cov')] 
nf1_mt_neuroectoderm_data = tns[tns$patient %in% c('PD51122') & tns$Embryology %in% c('CNS', 'PNS', 'Ectoderm / skin'), c('sampleID', 'chr', 'pos', 'ref', 'mut', 'gene', 'duplex_cov')] 
nf1_wt_endomesoderm_data = tns[tns$patient %in% c('PD50297', 'PD51123') & tns$Embryology %in% c('Endoderm', 'Endoderm / mesoderm' , 'Mesoderm', 'Mesoderm / blood'), c('sampleID', 'chr', 'pos', 'ref', 'mut', 'gene', 'duplex_cov')] 
nf1_mt_endomesoderm_data = tns[tns$patient %in% c('PD51122') & tns$Embryology %in% c('Endoderm', 'Endoderm / mesoderm' , 'Mesoderm', 'Mesoderm / blood'), c('sampleID', 'chr', 'pos', 'ref', 'mut', 'gene', 'duplex_cov')] 
nf1_wt_spleen_data = tns[tns$patient %in% c('PD50297', 'PD51123') & tns$Embryology %in% c('Mesoderm / spleen'), c('sampleID', 'chr', 'pos', 'ref', 'mut', 'gene', 'duplex_cov')] 
nf1_mt_spleen_data = tns[tns$patient %in% c('PD51122') & tns$Embryology %in% c('Mesoderm / spleen'), c('sampleID', 'chr', 'pos', 'ref', 'mut', 'gene', 'duplex_cov')] 

# make a vector of duplex coverage. This should be the mean per gene. Should recalculate for every group. Where don't have a mutation called in the group and so can't calculate coverage, use:
# the mean for all samples combined x the ratio of coverage of this gene to all genes in the overall group. 
nf1_wt_neuroectoderm_dc = sapply(gene_list, function(gene) mean(nf1_wt_neuroectoderm_data$duplex_cov[nf1_wt_neuroectoderm_data$gene==gene & !is.na(nf1_wt_neuroectoderm_data$gene)], na.rm=T))
median_dc_all_genes_this_group = median(nf1_wt_neuroectoderm_data$duplex_cov, na.rm=T)
nf1_wt_neuroectoderm_dc[is.na(nf1_wt_neuroectoderm_dc)] <- duplex_cov[names(nf1_wt_neuroectoderm_dc[is.na(nf1_wt_neuroectoderm_dc)]),'rel_to_median'] * median_dc_all_genes_this_group

nf1_mt_neuroectoderm_dc = sapply(gene_list, function(gene) mean(nf1_mt_neuroectoderm_data$duplex_cov[nf1_mt_neuroectoderm_data$gene==gene & !is.na(nf1_mt_neuroectoderm_data$gene)], na.rm=T))
median_dc_all_genes_this_group = median(nf1_mt_neuroectoderm_data$duplex_cov, na.rm=T)
nf1_mt_neuroectoderm_dc[is.na(nf1_mt_neuroectoderm_dc)] <- duplex_cov[names(nf1_mt_neuroectoderm_dc[is.na(nf1_mt_neuroectoderm_dc)]),'rel_to_median'] * median_dc_all_genes_this_group

nf1_wt_endomesoderm_dc = sapply(gene_list, function(gene) mean(nf1_wt_endomesoderm_data$duplex_cov[nf1_wt_endomesoderm_data$gene==gene & !is.na(nf1_wt_endomesoderm_data$gene)], na.rm=T))
median_dc_all_genes_this_group = median(nf1_wt_endomesoderm_data$duplex_cov, na.rm=T)
nf1_wt_endomesoderm_dc[is.na(nf1_wt_endomesoderm_dc)] <- duplex_cov[names(nf1_wt_endomesoderm_dc[is.na(nf1_wt_endomesoderm_dc)]),'rel_to_median'] * median_dc_all_genes_this_group

nf1_mt_endomesoderm_dc = sapply(gene_list, function(gene) mean(nf1_mt_endomesoderm_data$duplex_cov[nf1_mt_endomesoderm_data$gene==gene & !is.na(nf1_mt_endomesoderm_data$gene)], na.rm=T))
median_dc_all_genes_this_group = median(nf1_mt_endomesoderm_data$duplex_cov, na.rm=T)
nf1_mt_endomesoderm_dc[is.na(nf1_mt_endomesoderm_dc)] <- duplex_cov[names(nf1_mt_endomesoderm_dc[is.na(nf1_mt_endomesoderm_dc)]),'rel_to_median'] * median_dc_all_genes_this_group

nf1_wt_spleen_dc = sapply(gene_list, function(gene) mean(nf1_wt_spleen_data$duplex_cov[nf1_wt_spleen_data$gene==gene & !is.na(nf1_wt_spleen_data$gene)], na.rm=T))
median_dc_all_genes_this_group = median(nf1_wt_spleen_data$duplex_cov, na.rm=T)
nf1_wt_spleen_dc[is.na(nf1_wt_spleen_dc)] <- duplex_cov[names(nf1_wt_spleen_dc[is.na(nf1_wt_spleen_dc)]),'rel_to_median'] * median_dc_all_genes_this_group

nf1_mt_spleen_dc = sapply(gene_list, function(gene) mean(nf1_mt_spleen_data$duplex_cov[nf1_mt_spleen_data$gene==gene & !is.na(nf1_mt_spleen_data$gene)], na.rm=T))
median_dc_all_genes_this_group = median(nf1_mt_spleen_data$duplex_cov, na.rm=T)
nf1_mt_spleen_dc[is.na(nf1_mt_spleen_dc)] <- duplex_cov[names(nf1_mt_spleen_dc[is.na(nf1_mt_spleen_dc)]),'rel_to_median'] * median_dc_all_genes_this_group

# now run dNdS on each of these.
dnds_nf1_wt_neuroectoderm = dndscv::dndscv(mutations=nf1_wt_neuroectoderm_data[,1:5], dc = nf1_wt_neuroectoderm_dc,
                                           refdb = '/Users/hl11/Documents/Research/Methods/reference_files/grch38_refcds.rda', outmats=T,
                                           max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, onesided = T, gene_list = gene_list)

dnds_nf1_mt_neuroectoderm = dndscv::dndscv(mutations=nf1_mt_neuroectoderm_data[,1:5], dc = nf1_mt_neuroectoderm_dc,
                                           refdb = '/Users/hl11/Documents/Research/Methods/reference_files/grch38_refcds.rda', outmats=T,
                                           max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, onesided = T, gene_list = gene_list)

dnds_nf1_wt_endomesoderm = dndscv::dndscv(mutations=nf1_wt_endomesoderm_data[,1:5], dc = nf1_wt_endomesoderm_dc,
                                          refdb = '/Users/hl11/Documents/Research/Methods/reference_files/grch38_refcds.rda', outmats=T,
                                          max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, onesided = T, gene_list = gene_list)

dnds_nf1_mt_endomesoderm = dndscv::dndscv(mutations=nf1_mt_endomesoderm_data[,1:5], dc = nf1_mt_endomesoderm_dc,
                                          refdb = '/Users/hl11/Documents/Research/Methods/reference_files/grch38_refcds.rda', outmats=T,
                                          max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, onesided = T, gene_list = gene_list)

dnds_nf1_wt_spleen = dndscv::dndscv(mutations=nf1_wt_spleen_data[,1:5], dc = nf1_wt_spleen_dc,
                                    refdb = '/Users/hl11/Documents/Research/Methods/reference_files/grch38_refcds.rda', outmats=T,
                                    max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, onesided = T, gene_list = gene_list)

dnds_nf1_mt_spleen = dndscv::dndscv(mutations=nf1_mt_spleen_data[,1:5], dc = nf1_mt_spleen_dc,
                                    refdb = '/Users/hl11/Documents/Research/Methods/reference_files/grch38_refcds.rda', outmats=T,
                                    max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, onesided = T, gene_list = gene_list)


# now extract the confidence intervals for NF1
nf1_ne_nf1wt <- geneci(dnds_nf1_wt_neuroectoderm, gene_list = 'NF1', level = 0.95)
nf1_ne_nf1mt <- geneci(dnds_nf1_mt_neuroectoderm, gene_list = 'NF1', level = 0.95)
nf1_em_nf1wt <- geneci(dnds_nf1_wt_endomesoderm, gene_list = 'NF1', level = 0.95)
nf1_em_nf1mt <- geneci(dnds_nf1_mt_endomesoderm, gene_list = 'NF1', level = 0.95)
nf1_spleen_nf1wt <- geneci(dnds_nf1_wt_spleen, gene_list = 'NF1', level = 0.95)
nf1_spleen_nf1mt <- geneci(dnds_nf1_mt_spleen, gene_list = 'NF1', level = 0.95)

nf1_emb_spl <- as.data.frame(matrix(c(as.vector(unlist(nf1_ne_nf1mt[c('tru_mle', 'tru_low', 'tru_high')])), 
                                      as.vector(unlist(nf1_ne_nf1wt[c('tru_mle', 'tru_low', 'tru_high')])), 
                                      as.vector(unlist(nf1_em_nf1mt[c('tru_mle', 'tru_low', 'tru_high')])), 
                                      as.vector(unlist(nf1_em_nf1wt[c('tru_mle', 'tru_low', 'tru_high')])),
                                      as.vector(unlist(nf1_spleen_nf1mt[c('tru_mle', 'tru_low', 'tru_high')])), 
                                      as.vector(unlist(nf1_spleen_nf1wt[c('tru_mle', 'tru_low', 'tru_high')]))), 
                                    nrow=6, byrow = T))
colnames(nf1_emb_spl) <- c('tru_mle', 'tru_low', 'tru_high')
rownames(nf1_emb_spl) <- c('ne_nf1mt', 'ne_nf1wt', 'em_nf1mt', 'em_nf1wt', 'spleen_nf1mt', 'spleen_nf1wt')

write.table(nf1_emb_spl, 'nf1_truncating_dnds_by_germ_layer_and_germline_status_spleen_separate.txt', sep='\t', col.names=T, row.names=T, quote=F)

# now plot figure 2f
nf1_emb_spl[nf1_emb_spl==0] <- 0.00001

xlim_vals = c(1,1.5,2.5,3,4,4.5)
embcols=brewer.pal('Pastel2', n=8)[1:3]

necol = adjustcolor(embcols[1],alpha.f = 0.4)
mecol = adjustcolor(embcols[2],alpha.f = 0.4)
spleencol = adjustcolor(embcols[3],alpha.f = 0.4)

names(pdfFonts())

pdf('nf1_dnds_values_spleen_separate.pdf', height=5, width=6)
par(family="ArialMT")
plot(nf1_emb_spl$tru_mle~xlim_vals, log='y', pch=16, cex=2, col=c('black', 'grey'),
     xlab='', ylab='dN/dS ratio', axes=F, ylim=c(0.1,600), xlim=c(0.5,5), type='n')
#box(col='black')
abline(h=1, col='red')
rect(xleft = 0.5, xright = 2, ybottom = 0.01, ytop = 1000, col=necol, xpd=NA, border=NA)
text(x=1.25, y=0.03, labels='NEURO\nECTODERM', xpd=NA)
rect(xleft = 2, xright = 3.5, ybottom = 0.01, ytop = 1000, col=mecol, xpd=NA, border=NA)
text(x=2.75, y=0.03, labels='MESO/\nENDODERM', xpd=NA)
rect(xleft = 3.5, xright = 5, ybottom = 0.01, ytop = 1000, col=spleencol, xpd=NA, border=NA)
text(x=4.25, y=0.03, labels='SPLEEN', xpd=NA)

points(nf1_emb_spl$tru_mle~xlim_vals, cex=2, col=c('black', 'grey'),pch=16, add=T)
for (i in 1:nrow(nf1_emb_spl)) {
  tcol=c('grey', 'black')[i%%2 + 1]
  lines(x=c(xlim_vals[i],xlim_vals[i]),y=c(nf1_emb_spl[i,'tru_low'],nf1_emb_spl[i,'tru_high']), lwd=2, col=c(tcol))
}
axis(side=2, at=c(0.1,1,10,100), labels = c(0.1,1,10,100), tick=F, las=2)

dev.off()


# now plot figure 2g.
# this is just a stacked barplot for each patient of silent vs missense/inframe vs truncating mutations
table(tns$consequence2)
constab = t(table(tns[tns$gene=='NF1' & tns$consequence2 != 'noncoding', c('patient', 'consequence2')]))

# re-order.
table(tns[!tns$patient %in% c('PD50297', 'PD51122', 'PD51123'),c('patient', 'Tissue')]) # only one tissue per patient in cohort 2. 

ptorder <- c('PD50297', 'PD51123', 'PD51122', unique(tns$patient[tns$Tissue=='Nerve']), unique(tns$patient[tns$Tissue=='Muscle']), unique(tns$patient[tns$Tissue=='Blood']))
# drop PD61048 as that did not have a convincing germline mutations.
ptorder = ptorder[ptorder != 'PD61048' & !is.na(ptorder)]

ctab = constab[c('truncating', 'missense_inframe', 'synonymous'),ptorder]

# also need to calculate the duplex coverage over NF1 per patient
# for every sample I need to calculate the average, and then sum up all the samples per patient.

patient = 'PD51122'
per_pt_nf1_cov = sapply(ptorder, function(patient) {
  tblocks = unique(tns[tns$patient==patient,'sampleID'])
  tblockmeans = sapply(tblocks, function(block) {
    mean(tns$duplex_cov[tns$gene == 'NF1' & !is.na(tns$gene) & tns$sampleID==block], na.rm=T)
  })
  return(sum(tblockmeans, na.rm=T))
})


# plot

nervepts = unique(tns$patient[tns$Tissue=='Nerve'])
nervepts = nervepts[!is.na(nervepts)]
musclepts = unique(tns$patient[tns$Tissue=='Muscle'])
musclepts = musclepts[!is.na(musclepts)]
bldpts = unique(tns$patient[tns$Tissue=='Blood'])
bldpts = bldpts[!is.na(bldpts) & bldpts!='PD61048']
nf1wtpts = c('PD50297', 'PD51123')

dtab=cbind(rowSums(ctab[,nf1wtpts]), ctab[,c('PD51122')], rowSums(ctab[,nervepts]), rowSums(ctab[,musclepts]), rowSums(ctab[,bldpts]))
colnames(dtab)=c('nf1wt','PD51122', 'nerve', 'muscle', 'blood')

per_tissue_nf1_cov = c(sum(per_pt_nf1_cov[nf1wtpts]), per_pt_nf1_cov["PD51122"], sum(per_pt_nf1_cov[nervepts]), sum(per_pt_nf1_cov[musclepts]), sum(per_pt_nf1_cov[bldpts]))
names(per_tissue_nf1_cov) = c('nf1wt','PD51122', 'nerve', 'muscle', 'blood')

# I should also calculate the dNdS for each group. 
# first analysis
# cohort 1 NF1 wt (PD50297, PD51123) vs NF1 mutant (PD51122), neuroectoderm and mesoderm analysed separately. 

nf1wt_data = tns[tns$patient %in% c('PD50297', 'PD51123'), c('sampleID', 'chr', 'pos', 'ref', 'mut', 'gene', 'duplex_cov')] 
pd51122_data = tns[tns$patient %in% c('PD51122'), c('sampleID', 'chr', 'pos', 'ref', 'mut', 'gene', 'duplex_cov')] 
nerve_data = tns[tns$patient %in% nervepts, c('sampleID', 'chr', 'pos', 'ref', 'mut', 'gene', 'duplex_cov')] 
muscle_data = tns[tns$patient %in% musclepts, c('sampleID', 'chr', 'pos', 'ref', 'mut', 'gene', 'duplex_cov')] 
blood_data = tns[tns$patient %in% bldpts, c('sampleID', 'chr', 'pos', 'ref', 'mut', 'gene', 'duplex_cov')] 

# make a vector of duplex coverage. This should be the mean per gene. Should recalculate for every group. Where don't have a mutation called in the group and so can't calculate coverage, use:
# the mean for all samples combined x the ratio of coverage of this gene to all genes in the overall group. 
nf1wt_dc = sapply(gene_list, function(gene) mean(nf1wt_data$duplex_cov[nf1wt_data$gene==gene & !is.na(nf1wt_data$gene)], na.rm=T))
median_dc_all_genes_this_group = median(nf1wt_data$duplex_cov, na.rm=T)
nf1wt_dc[is.na(nf1wt_dc)] <- duplex_cov[names(nf1wt_dc[is.na(nf1wt_dc)]),'rel_to_median'] * median_dc_all_genes_this_group

pd51122_dc = sapply(gene_list, function(gene) mean(pd51122_data$duplex_cov[pd51122_data$gene==gene & !is.na(pd51122_data$gene)], na.rm=T))
median_dc_all_genes_this_group = median(pd51122_data$duplex_cov, na.rm=T)
pd51122_dc[is.na(pd51122_dc)] <- duplex_cov[names(pd51122_dc[is.na(pd51122_dc)]),'rel_to_median'] * median_dc_all_genes_this_group

nerve_dc = sapply(gene_list, function(gene) mean(nerve_data$duplex_cov[nerve_data$gene==gene & !is.na(nerve_data$gene)], na.rm=T))
median_dc_all_genes_this_group = median(nerve_data$duplex_cov, na.rm=T)
nerve_dc[is.na(nerve_dc)] <- duplex_cov[names(nerve_dc[is.na(nerve_dc)]),'rel_to_median'] * median_dc_all_genes_this_group

muscle_dc = sapply(gene_list, function(gene) mean(muscle_data$duplex_cov[muscle_data$gene==gene & !is.na(muscle_data$gene)], na.rm=T))
median_dc_all_genes_this_group = median(muscle_data$duplex_cov, na.rm=T)
muscle_dc[is.na(muscle_dc)] <- duplex_cov[names(muscle_dc[is.na(muscle_dc)]),'rel_to_median'] * median_dc_all_genes_this_group

blood_dc = sapply(gene_list, function(gene) mean(blood_data$duplex_cov[blood_data$gene==gene & !is.na(blood_data$gene)], na.rm=T))
median_dc_all_genes_this_group = median(blood_data$duplex_cov, na.rm=T)
blood_dc[is.na(blood_dc)] <- duplex_cov[names(blood_dc[is.na(blood_dc)]),'rel_to_median'] * median_dc_all_genes_this_group

# now run dNdS on each of these.
dnds_nf1wt = dndscv::dndscv(mutations=nf1wt_data[,1:5], dc = nf1wt_dc,
                            refdb = '/Users/hl11/Documents/Research/Methods/reference_files/grch38_refcds.rda', outmats=T,
                            max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, onesided = T, gene_list = gene_list)
dnds_pd51122 = dndscv::dndscv(mutations=pd51122_data[,1:5], dc = pd51122_dc,
                              refdb = '/Users/hl11/Documents/Research/Methods/reference_files/grch38_refcds.rda', outmats=T,
                              max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, onesided = T, gene_list = gene_list)
dnds_nerve = dndscv::dndscv(mutations=nerve_data[,1:5], dc = nerve_dc,
                            refdb = '/Users/hl11/Documents/Research/Methods/reference_files/grch38_refcds.rda', outmats=T,
                            max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, onesided = T, gene_list = gene_list)
dnds_muscle = dndscv::dndscv(mutations=muscle_data[,1:5], dc = muscle_dc,
                             refdb = '/Users/hl11/Documents/Research/Methods/reference_files/grch38_refcds.rda', outmats=T,
                             max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, onesided = T, gene_list = gene_list)
dnds_blood = dndscv::dndscv(mutations=blood_data[,1:5], dc = blood_dc,
                            refdb = '/Users/hl11/Documents/Research/Methods/reference_files/grch38_refcds.rda', outmats=T,
                            max_muts_per_gene_per_sample = Inf, max_coding_muts_per_sample = Inf, onesided = T, gene_list = gene_list)


# now extract the confidence intervals for NF1
nf1_nf1wt <- geneci(dnds_nf1wt, gene_list = 'NF1', level = 0.95)
nf1_pd51122 <- geneci(dnds_pd51122, gene_list = 'NF1', level = 0.95)
nf1_nerve <- geneci(dnds_nerve, gene_list = 'NF1', level = 0.95)
nf1_muscle <- geneci(dnds_muscle, gene_list = 'NF1', level = 0.95)
nf1_blood <- geneci(dnds_blood, gene_list = 'NF1', level = 0.95)

nf1_by_tissue <- as.data.frame(matrix(c(
  as.vector(unlist(nf1_nf1wt[c('tru_mle', 'tru_low', 'tru_high')])), 
  as.vector(unlist(nf1_pd51122[c('tru_mle', 'tru_low', 'tru_high')])), 
  as.vector(unlist(nf1_nerve[c('tru_mle', 'tru_low', 'tru_high')])), 
  as.vector(unlist(nf1_muscle[c('tru_mle', 'tru_low', 'tru_high')])), 
  as.vector(unlist(nf1_blood[c('tru_mle', 'tru_low', 'tru_high')]))),
  nrow=5, byrow = T))
colnames(nf1_by_tissue) <- c('tru_mle', 'tru_low', 'tru_high')
rownames(nf1_by_tissue) <- c('nf1wt', 'pd51122', 'nerve', 'muscle', 'blood')

write.table(nf1_by_tissue, 'nf1_truncating_dnds_by_pt_and_tissue_for_cohort2_nf1wt_pts_merged.txt', sep='\t', col.names=T, row.names=T, quote=F)

nf1_by_tissue[nf1_by_tissue==0] <- 0.00001

pdf('nf1_mutation_counts_by_tissue_with_dnds_nf1wt_merged.pdf', height=6, width=4)
layout(matrix(c(1,2,2,3), ncol=1))
par(mar=c(2,4,1,2))

plot(nf1_by_tissue$tru_mle~c(1.5,2.5,3.5,4.5,5.5), log='y', pch=16, cex=2, col=c('grey', 'black', 'black', 'black', 'black'),
     xlab='', ylab='dN/dS ratio', axes=F, ylim=c(0.1,500), xlim=c(1,6))
for (i in 1:nrow(nf1_by_tissue)) {
  if (i == 1) {
    tcol='grey'
  } else {
    tcol='black'
  }
  lines(x=c(i+0.5,i+0.5),y=c(nf1_by_tissue[i,'tru_low'],nf1_by_tissue[i,'tru_high']), lwd=2, col=tcol)
}
axis(side=2, at=c(0.1,1,10,100), labels = c(0.1,1,10,100), tick=F, las=2)
box(col='black')
abline(h=1, col='red')

barplot(dtab, beside=F, col=c('darkblue', 'lightblue', 'white'), las=2, ylab='Variants in NF1', axes=F, names.arg=rep('', ncol(dtab)))
axis(side=2, at=c(0,15,30), labels=c(0,15,30),tick = F, las=2)
legend(box.lty=0, 'topright', fill=c('darkblue', 'lightblue', 'white'), legend=c('truncating', 'missense/inframe', 'synonymous'))
barplot(per_tissue_nf1_cov, las=2, ylab='duplex coverage', axes=F, names.arg='', border=NA)
axis(side=2, at=c(0,4000,8000), labels=c(0,4000,8000),tick = F, las=2)
dev.off()
