manifest_wgs = read.table('/lustre/scratch126/casm/team274sb/to3/paed_brain/autopsy_reanalysis/wgs_data/00_supplementary_data/wgs_paeds_pm_sample_metadata_20240225.txt', header = T, sep = '\t', stringsAsFactors = F, quote = '')
tum_root = read.table("/lustre/scratch126/casm/team274sb/to3/paed_brain/autopsy_reanalysis/wgs_data/03_DNA_analyses/tum_root_vars/PD51122_put_tum_root_snp_combined_method_reviewed.txt", header = T, sep = "\t", stringsAsFactors = F)
tum_root_vars = tum_root[tum_root$Keep == "Y", 1] #mutations shared between glioma and normal tissues prior to malignant transformation
#"chr15_59361799_T_A" "chr17_4949465_C_A"  "chr4_178386347_G_A" "chr8_133329433_A_G"

soi = manifest_wgs[manifest_wgs$nf1_ssm_pileup_call == "Y",]$sample #all samples that contain a second NF1 hit (sub or indel) after genotyping

#identify normal samples without second NF1 hit
norm_pd51122_cns_samples = manifest_wgs[manifest_wgs$case.id == "PD51122" & 
                                          manifest_wgs$CNS == "Y" & 
                                          manifest_wgs$purity_trunk < 0.01 & 
                                          manifest_wgs$nf1_ssm_pileup_call == "N" &
                                          is.na(manifest_wgs$nf1_dh_clone_size),]$sample
norm_pd51122_meso_samples = manifest_wgs[manifest_wgs$case.id == "PD51122" & 
                                           manifest_wgs$predominant_germ_layer == "Mesoderm" & 
                                           manifest_wgs$purity_trunk < 0.01,]$sample
norm_pd51122_cns_samples[norm_pd51122_cns_samples %in% norm_pd51122_meso_samples] #0

#define set of pure tumour biopsies
high_purity_tumour_tissues = manifest_wgs[manifest_wgs$purity_trunk >= 0.4 &
                                            manifest_wgs$case.id == "PD51122" & 
                                            manifest_wgs$CNS == "Y" &
                                            manifest_wgs$nf1_ssm_pileup_call == "N",]$sample

#get subs for relevant WGS samples
all_subs = read.table("/lustre/scratch126/casm/team274sb/to3/paed_brain/autopsy_reanalysis/wgs_data/02_DNA_processing/variant_calls/01_subs/09_swater/all_sample_subs_distribution_annotated_20240302.txt", header = T, sep = '\t', stringsAsFactors = F, quote = '')
all_subs_pure_tumour = unique(all_subs[all_subs$Sample %in% high_purity_tumour_tissues,]$mutID)

#get calls matrix - what passed error rate?
er_data = read.table("/lustre/scratch126/casm/team274sb/to3/paed_brain/autopsy_reanalysis/wgs_data/02_DNA_processing/shearwater_framework/PD51122/shearwater_combined_snv_qval_pass_mat.txt", header = T, sep = " ", stringsAsFactors = F) #matrix containing all subs and samples for case with binary mutation present (1) or absent (0) above the locus-specific error rate
er_data_norm = er_data[(!row.names(er_data) %in% all_subs_pure_tumour) |  row.names(er_data) %in% tum_root_vars,] #not in pure tumour unless root mutation, minimise risk of tumour contamination
er_data_norm_shared = er_data_norm[rowSums(er_data_norm) > 1,] #variants found in normal tissues that are shared between at least two samples

#list of second hits and which samples they are found in
nf1_wgs_hits = list(chr17_31200443_C_T = c("PD51122g_lo0003", "PD51122g_lo0005", "PD51122p"),
                    chr17_31214524_A_G = c("PD51122b_lo0018", "PD51122z_lo0006"),
                    chr17_31226459_A_AC = c("PD51122c", "PD51122d_lo0004", "PD51122d_lo0005", "PD51122d_lo0010", "PD51122d_lo0011", "PD51122d"))
manifest_wgs[manifest_wgs$sample %in% as.character(unlist(nf1_wgs_hits)), c("sample", "purity_trunk")]
#          sample purity_trunk
# PD51122b_lo0018  0.000000000
#        PD51122c  0.160300435
#        PD51122d  0.001972213
# PD51122d_lo0004  0.000000000
# PD51122d_lo0005  0.035280218
# PD51122d_lo0010  0.069459457
# PD51122d_lo0011  0.000000000
# PD51122g_lo0003  0.000000000
# PD51122g_lo0005  0.000000000
#        PD51122p  0.000000000
# PD51122z_lo0006  0.000000000

input = er_data_norm_shared

#PAIRWISE COMPARISONS
##Normal tissues harbours second NF1 mutations
nf1_norm_null = list()
nf1_norm_null_share_median = list()
for(mut in names(nf1_wgs_hits)){
  
  mut_tissues = nf1_wgs_hits[[mut]]
  combination_mat = combn(length(mut_tissues),2)
  nf1_comparisons = c()
  share_median = c()
  
  for(i in 1:ncol(combination_mat)){
    
    pair = mut_tissues[combination_mat[,i]]
    sample1 = input[,pair[1]]
    sample2 = input[,pair[2]]
    
    region_involved = unique(manifest_wgs[manifest_wgs$sample %in% c(pair[1], pair[2]),]$bulk_sample)
    er_data_mod = er_data[, !colnames(er_data) %in% manifest_wgs[manifest_wgs$bulk_sample %in% region_involved,]$sample]
    
    iter_muts = intersect(row.names(input)[sample1 == 1], row.names(input)[sample2 == 1])
    
    share_vec = c()
    
    if(length(iter_muts) == 0){
      
      share_vec = 0
      
    } else {
      for(j in iter_muts){

        share_vec = c(share_vec, rowSums(er_data_mod[j,]) / ncol(er_data_mod))
        
      }
    }
    
    iter_length = length(iter_muts)
    
    nf1_comparisons = c(nf1_comparisons, iter_length)
    share_median = c(share_median, median(share_vec))
    
  }
  
  nf1_norm_null[[mut]] = nf1_comparisons
  nf1_norm_null_share_median[[mut]] = share_median
  
}

###what is the average number of subs shared between normal tissues harbouring each of these NF1 variants?
lapply(nf1_norm_null, mean)
###$chr17_31200443_C_T
###5.666667

###$chr17_31214524_A_G
###2

###$chr17_31226459_A_AC
###5.8

##Normal CNS without a second NF1 hit vs normal CNS without a second NF1 hit
choose(length(norm_pd51122_cns_samples),2) #12090
combination_mat = combn(length(norm_pd51122_cns_samples),2) 

non_nf1_comparisons = c()
norm_share_median = c()
for(i in 1:ncol(combination_mat)){
  
  pair = norm_pd51122_cns_samples[combination_mat[,i]]
  sample1 = input[,pair[1]]
  sample2 = input[,pair[2]]
  
  region_involved = unique(manifest_wgs[manifest_wgs$sample %in% c(pair[1], pair[2]),]$bulk_sample)
  er_data_mod = er_data[, !colnames(er_data) %in% manifest_wgs[manifest_wgs$bulk_sample %in% region_involved,]$sample]
  
  iter_muts = intersect(row.names(input)[sample1 == 1], row.names(input)[sample2 == 1])
  
  share_vec = c()
  
  if(length(iter_muts) == 0){
    
    share_vec = 0
    
  } else {
    for(j in iter_muts){
      
      share_vec = c(share_vec, rowSums(er_data_mod[j,]) / ncol(er_data_mod))
      
    }
  }
  
  iter_length = length(iter_muts)
  
  non_nf1_comparisons = c(non_nf1_comparisons, iter_length)
  norm_share_median = c(norm_share_median, median(share_vec))
  
}

mean(non_nf1_comparisons)
###1.874359
median(non_nf1_comparisons)
###2

##Normal CNS without a second NF1 hit vs normal mesoderm-derived tissue without a second NF1 hit
cns_v_meso_comparisons = c()
cns_v_meso_share_median = c()
for(cns_sample in norm_pd51122_cns_samples){
  
  for(meso_sample in norm_pd51122_meso_samples){
    
    sample1 = input[,cns_sample]
    sample2 = input[,meso_sample]
    
    region_involved = unique(manifest_wgs[manifest_wgs$sample %in% c(cns_sample, meso_sample),]$bulk_sample)
    er_data_mod = er_data[, !colnames(er_data) %in% manifest_wgs[manifest_wgs$bulk_sample %in% region_involved,]$sample]
    
    iter_muts = intersect(row.names(input)[sample1 == 1], row.names(input)[sample2 == 1])
    
    share_vec = c()
    
    if(length(iter_muts) == 0){
      share_vec = 0
    } else {
      for(j in iter_muts){
        
        share_vec = c(share_vec, rowSums(er_data_mod[j,]) / ncol(er_data_mod))
        
      }
    }
    
    iter_length = length(iter_muts)
    
    cns_v_meso_comparisons = c(cns_v_meso_comparisons, iter_length)
    cns_v_meso_share_median = c(cns_v_meso_share_median, median(share_vec))
    
  }
  
}

mean(cns_v_meso_comparisons)
###2.362981
median(cns_v_meso_comparisons)
###2

cns_meso_pairs = c()
for(cns in norm_pd51122_cns_samples){
  for(meso in norm_pd51122_meso_samples){
    
    cns_meso_pairs = c(cns_meso_pairs, paste0(cns, "_", meso))
    
  }
}

length(cns_meso_pairs) #8739
length(cns_v_meso_comparisons) #8739

#MERGE DATA AND PLOT
df = rbind(data.frame(group = rep("chr17_31200443_C_T", length(nf1_norm_null[["chr17_31200443_C_T"]])), count = nf1_norm_null[["chr17_31200443_C_T"]], median_shared = nf1_norm_null_share_median[["chr17_31200443_C_T"]]),
           data.frame(group = rep("chr17_31214524_A_G", length(nf1_norm_null[["chr17_31214524_A_G"]])), count = nf1_norm_null[["chr17_31214524_A_G"]], median_shared = nf1_norm_null_share_median[["chr17_31214524_A_G"]]),
           data.frame(group = rep("chr17_31226459_A_AC", length(nf1_norm_null[["chr17_31226459_A_AC"]])), count = nf1_norm_null[["chr17_31226459_A_AC"]], median_shared = nf1_norm_null_share_median[["chr17_31226459_A_AC"]]),
           data.frame(group = rep("normal_cns", length(non_nf1_comparisons)), count = non_nf1_comparisons, median_shared = norm_share_median),
           data.frame(group = rep("cns_meso", length(cns_v_meso_comparisons)), count = cns_v_meso_comparisons, median_shared = cns_v_meso_share_median))
      
df$group = factor(df$group, levels = c("chr17_31200443_C_T", "chr17_31226459_A_AC", "normal_cns", "chr17_31214524_A_G", "cns_meso"))
      
hist(df$median_shared, breaks = 100)
      
table(df$group)
##chr17_31200443_C_T chr17_31226459_A_AC          normal_cns  chr17_31214524_A_G            cns_meso 
##                 3                  15               12090                   1                8736
aggregate(count ~ group, df, mean)
##              group    count
##  chr17_31200443_C_T 5.666667
## chr17_31226459_A_AC 5.800000
##          normal_cns 1.874359
##  chr17_31214524_A_G 2.000000
##            cns_meso 2.362981
      
##write.table(df, "/lustre/scratch126/casm/team274sb/to3/paed_brain/autopsy_reanalysis/wgs_data/03_DNA_analyses/draft_figures/nf1_sharedness_analysis_data_20240411.txt", col.names = T, row.names = F, sep = "\t", quote = F)
##df = read.table("/lustre/scratch126/casm/team274sb/to3/paed_brain/autopsy_reanalysis/wgs_data/03_DNA_analyses/draft_figures/nf1_sharedness_analysis_data_20240411.txt", header = T, sep = "\t", stringsAsFactors = F)
df$group = factor(df$group, levels = c("chr17_31200443_C_T", "chr17_31226459_A_AC", "normal_cns", "chr17_31214524_A_G", "cns_meso"))
      
##generate P values - permutation test
###chr17_31200443_C_T
nRand = 1000
set.seed(42)
chr17_31200443_C_T_mean_shared = c()
      
choose(nrow(df[df$group == "normal_cns",]), nrow(df[df$group == "chr17_31200443_C_T",])) #294455641480
      
for(i in 1:nRand){
  rand_sample = sample(df[df$group == "normal_cns",]$count, nrow(df[df$group == "chr17_31200443_C_T",]), replace = F)
        
  chr17_31200443_C_T_mean_shared = c(chr17_31200443_C_T_mean_shared, (mean(rand_sample)))
}
      
hist(chr17_31200443_C_T_mean_shared, breaks = 100)
length(chr17_31200443_C_T_mean_shared[chr17_31200443_C_T_mean_shared > mean(df[df$group == "chr17_31200443_C_T",]$count)]) / length(chr17_31200443_C_T_mean_shared) #0, <0.001
      
##chr17_31226459_A_AC
nRand = 1000
set.seed(42)
chr17_31226459_A_AC_mean_shared = c()
      
choose(nrow(df[df$group == "normal_cns",]), nrow(df[df$group == "chr17_31226459_A_AC",])) #1.306537e+49
      
for(i in 1:nRand){
  rand_sample = sample(df[df$group == "normal_cns",]$count, nrow(df[df$group == "chr17_31226459_A_AC",]), replace = F)
        
  chr17_31226459_A_AC_mean_shared = c(chr17_31226459_A_AC_mean_shared, (mean(rand_sample)))
}
      
hist(chr17_31226459_A_AC_mean_shared, breaks = 100)
length(chr17_31226459_A_AC_mean_shared[chr17_31226459_A_AC_mean_shared > mean(df[df$group == "chr17_31226459_A_AC",]$count)]) / length(chr17_31226459_A_AC_mean_shared) #0, <0.001
      
##chr17_31214524_A_G
nRand = 1000
set.seed(42)
chr17_31214524_A_G_mean_shared = c()
      
choose(nrow(df[df$group == "cns_meso",]), nrow(df[df$group == "chr17_31214524_A_G",])) #8736
      
for(i in 1:nRand){
  rand_sample = sample(df[df$group == "cns_meso",]$count, nrow(df[df$group == "chr17_31214524_A_G",]), replace = F)
        
  chr17_31214524_A_G_mean_shared = c(chr17_31214524_A_G_mean_shared, (mean(rand_sample)))
}
      
hist(chr17_31214524_A_G_mean_shared, breaks = 100)
      
length(chr17_31214524_A_G_mean_shared[chr17_31214524_A_G_mean_shared > mean(df[df$group == "chr17_31214524_A_G",]$count)]) / length(chr17_31214524_A_G_mean_shared) #0.373
      
length(df[df$group == "cns_meso",]$count[df[df$group == "cns_meso",]$count > mean(df[df$group == "chr17_31214524_A_G",]$count)]) / length(df[df$group == "cns_meso",]$count) #0.4078526 if you just compare the single pairwise comparison between the second hit tissues with all of those in the CNS vs MESO pairwise category

##Plot
library(ggplot2)
library(ggbeeswarm)
library(circlize)
col_fun = colorRamp2(c(0, 1), c("white", "black"))
      
pdf("/lustre/scratch126/casm/team274sb/to3/paed_brain/autopsy_reanalysis/wgs_data/03_DNA_analyses/draft_figures/pairwise_shared_mut_comparison_20240413.pdf", height = 3.5, width = 7, useDingbats = F)
ggplot(df) + 
  geom_quasirandom(mapping = aes(x = group, y = count, fill = group), pch = 21, stroke = NA, show.legend = F) +
  geom_pointrange(mapping = aes(x = group, y = count), stat = "summary", 
                  shape=19, 
                  fun.min = function(z) {quantile(z,0.25)},
                  fun.max = function(z) {quantile(z,0.75)},
                  fun = median, 
                  fill="black",
                  size = 1) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.x = element_blank()) +
  scale_fill_manual(values = c('#E83B3F', '#E83B3F', '#9A9999', '#E83B3F', '#9A9999')) +
  scale_x_discrete(labels = c("p.R304*", "p.I679fs*21", "CNS vs CNS", "p.Y489C", "CNS vs MES"))
dev.off()
      