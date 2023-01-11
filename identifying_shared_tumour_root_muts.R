
### PD50297 ###
setwd("/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/03_DNA_analyses/tum_root_vars")

patient = "PD50297"
truncal_cluster = 9

manifest = read.table('/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/00_supplementary_data/wgs_paeds_pm_sample_metadata_20221015.txt', header = T, sep = '\t', stringsAsFactors = F, quote = '')

tumour_bulk_samples = manifest[manifest$case.id == patient & manifest$keep == "Y" & manifest$purity_trunk >= 0.4 & manifest$experimental_arm == "WGS_bulk",]$sample
non_tumour_samples = manifest[manifest$case.id == patient & manifest$keep == "Y" & manifest$purity_trunk == 0 & manifest$picard_median_coverage >= 30,]$sample #high depth so that low depth with low level tumour purity not misidentified as shared

#get variants in trunk
vars = read.table(paste0('/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/03_DNA_analyses/nd_dpclust/01_Input/', patient, '/md_out/', patient, '_allClusterassignmentsFromParallelRuns.txt'), header = T, sep = '\t', stringsAsFactors = F)
cluster_vars = vars[vars$cluster.no == truncal_cluster,]
cluster_vars$ID = paste0(cluster_vars$chr, "_", cluster_vars$pos)

truncal_vars = paste0("chr", cluster_vars$ID)
length(truncal_vars) #929

#read in VAF dataframe
VAF = read.table(paste0("/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data_mk2/merged_final_subs/", patient, "_VAF_final_merged_snps.txt"), header = T, sep = "\t", stringsAsFactors = F)

VAF_filter = VAF
VAF_filter$ID = paste0(unlist(lapply(strsplit(row.names(VAF_filter), '_'), "[[", 1)), "_", unlist(lapply(strsplit(row.names(VAF_filter), '_'), "[[", 2)))

trunc_normal_VAF = VAF[VAF_filter$ID %in% truncal_vars, non_tumour_samples]
nrow(trunc_normal_VAF) #929

put_ee_vars = row.names(trunc_normal_VAF)[rowSums(trunc_normal_VAF > 0.1) > 0] #greater than

#intersection method
sub.files = list.files("/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/02_DNA_processing/variant_calls/subs/07_indel", ".standard.gnomAD.AF.0.001.PON.read.direction.low_cov.bq.mq.binom.bbinom.5bp.indel.filtered.txt", full.names = T)
patient.sub.files = sub.files[grepl(patient, sub.files)]

manifest_flt = manifest[manifest$keep == "Y" & manifest$case.id == patient,]

#read in patient subs
all_subs = data.frame()
for(file in patient.sub.files){
  
  data = read.table(file, header = T, sep = '\t', stringsAsFactors = F, colClasses = c("character", "character", "character", "numeric", "character", "character", "numeric", "numeric", "numeric", "character", "character", "character", "character", "character", "character"))
  all_subs = rbind(all_subs, data)
}

all_subs$mutID = paste(all_subs$Chr, all_subs$Position, all_subs$Wildtype, all_subs$Mutant, sep = "_")
all_subs_flt = all_subs[all_subs$Sample %in% manifest_flt$sample,]

intersect_trunk_vars = names(table(all_subs_flt[all_subs_flt$Sample %in% tumour_bulk_samples,]$mutID)[table(all_subs_flt[all_subs_flt$Sample %in% tumour_bulk_samples,]$mutID) == length(tumour_bulk_samples)])
intersect_root_vars = unique(all_subs_flt[all_subs_flt$mutID %in% intersect_trunk_vars & all_subs_flt$Sample %in% non_tumour_samples,]$mutID)

aggregate(VAF ~ mutID, all_subs_flt[all_subs_flt$mutID %in% intersect_root_vars & all_subs_flt$Sample %in% non_tumour_samples,], max) #3 are seen at a VAF less than 0.1 - where are these found?

#                 mutID        VAF
#  chr12_112505263_G_A 0.09000000
#   chr17_67362034_C_A 0.69696970
#   chr17_72915405_C_A 0.13953488
#   chr19_37176351_C_T 0.12195122
#    chr4_68023594_C_T 0.56000000
#   chr6_120782973_G_A 0.18367347
#    chr6_78642457_C_T 0.08823529
#   chr7_143282253_C_T 0.08163265
#   chr7_146198061_C_T 0.20588235
#   chr8_87155015_G_A 0.47058824
#   chrX_52010308_C_T 1.00000000

all_subs_flt[all_subs_flt$mutID %in% c("chr12_112505263_G_A", "chr6_78642457_C_T", "chr7_143282253_C_T") & all_subs_flt$Sample %in% non_tumour_samples,] #all have read depths much higher than 30x
min(all_subs_flt[all_subs_flt$mutID %in% c("chr12_112505263_G_A", "chr6_78642457_C_T", "chr7_143282253_C_T") & all_subs_flt$Sample %in% non_tumour_samples,]$Total_read_depth) #68

all_put_root_vars = sort(unique(c(put_ee_vars, intersect_root_vars)))

write.table(VAF[all_put_root_vars, ], paste0(patient, "_put_tum_root_snp_combined_method.txt"), col.names = T, row.names = T, sep = "\t", quote = F)

#### PD51122 ###
setwd("/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/03_DNA_analyses/tum_root_vars")

patient = "PD51122"
truncal_cluster = "6_9"

manifest = read.table('/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/00_supplementary_data/wgs_paeds_pm_sample_metadata_20221015.txt', header = T, sep = '\t', stringsAsFactors = F, quote = '')

tumour_bulk_samples = manifest[manifest$case.id == patient & manifest$keep == "Y" & manifest$purity_trunk >= 0.4 & manifest$experimental_arm == "WGS_bulk",]$sample
non_tumour_samples = manifest[manifest$case.id == patient & manifest$keep == "Y" & manifest$purity_trunk == 0 & manifest$picard_median_coverage >= 30,]$sample #high depth so that low depth with low level tumour purity not misidentified as shared

#get variants in trunk
vars = read.table(paste0('/lustre/scratch119/realdata/mdt1/team274/to3/paed_brain/autopsy/wgs_data/03_DNA_analyses/nd_dpclust/01_Input/', patient, '/md_out/', patient, '_allClusterassignmentsFromParallelRuns_clusters_reviewed.txt'), header = T, sep = '\t', stringsAsFactors = F) # truncal cluster different from the original output
cluster_vars = vars[vars$cluster.no == truncal_cluster,]
cluster_vars$ID = paste0(cluster_vars$chr, "_", cluster_vars$pos)

truncal_vars = paste0("chr", cluster_vars$ID)
length(truncal_vars) #972

#read in VAF dataframe
VAF = read.table(paste0("/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data_mk2/merged_final_subs/", patient, "_VAF_final_merged_snps.txt"), header = T, sep = "\t", stringsAsFactors = F)

VAF_filter = VAF
VAF_filter$ID = paste0(unlist(lapply(strsplit(row.names(VAF_filter), '_'), "[[", 1)), "_", unlist(lapply(strsplit(row.names(VAF_filter), '_'), "[[", 2)))

trunc_normal_VAF = VAF[VAF_filter$ID %in% truncal_vars, non_tumour_samples]
nrow(trunc_normal_VAF) #972
put_ee_vars = row.names(trunc_normal_VAF)[rowSums(trunc_normal_VAF > 0.1) > 0] #greater than

#intersection method
sub.files = list.files("/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/02_DNA_processing/variant_calls/subs/07_indel", ".standard.gnomAD.AF.0.001.PON.read.direction.low_cov.bq.mq.binom.bbinom.5bp.indel.filtered.txt", full.names = T)
patient.sub.files = sub.files[grepl(patient, sub.files)]

manifest_flt = manifest[manifest$keep == "Y" & manifest$case.id == patient,]

#read in patient subs
all_subs = data.frame()
for(file in patient.sub.files){
  
  data = read.table(file, header = T, sep = '\t', stringsAsFactors = F, colClasses = c("character", "character", "character", "numeric", "character", "character", "numeric", "numeric", "numeric", "character", "character", "character", "character", "character", "character"))
  all_subs = rbind(all_subs, data)
}

all_subs$mutID = paste(all_subs$Chr, all_subs$Position, all_subs$Wildtype, all_subs$Mutant, sep = "_")
all_subs_flt = all_subs[all_subs$Sample %in% manifest_flt$sample,]

intersect_trunk_vars = names(table(all_subs_flt[all_subs_flt$Sample %in% tumour_bulk_samples,]$mutID)[table(all_subs_flt[all_subs_flt$Sample %in% tumour_bulk_samples,]$mutID) == length(tumour_bulk_samples)])
intersect_root_vars = unique(all_subs_flt[all_subs_flt$mutID %in% intersect_trunk_vars & all_subs_flt$Sample %in% non_tumour_samples,]$mutID)

aggregate(VAF ~ mutID, all_subs_flt[all_subs_flt$mutID %in% intersect_root_vars & all_subs_flt$Sample %in% non_tumour_samples,], max) #1 is seen at a VAF less than 0.1 - where is this found?

#               mutID       VAF
#chr1_169603900_G_T 0.5909091
#chr15_59361799_T_A 0.4117647
# chr17_4949465_C_A 0.0776699
#chr8_133329433_A_G 0.3703704

all_subs_flt[all_subs_flt$mutID %in% c("chr17_4949465_C_A") & all_subs_flt$Sample %in% non_tumour_samples,] #read depths much higher than 30x
min(all_subs_flt[all_subs_flt$mutID %in% c("chr17_4949465_C_A") & all_subs_flt$Sample %in% non_tumour_samples,]$Total_read_depth) #103

all_put_root_vars = sort(unique(c(put_ee_vars, intersect_root_vars)))

write.table(VAF[all_put_root_vars, ], paste0(patient, "_put_tum_root_snp_combined_method.txt"), col.names = T, row.names = T, sep = "\t", quote = F)

### PD51123 ###
setwd("/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/03_DNA_analyses/tum_root_vars")

patient = "PD51123"
truncal_cluster = 2

manifest = read.table('/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/00_supplementary_data/wgs_paeds_pm_sample_metadata_20221015.txt', header = T, sep = '\t', stringsAsFactors = F, quote = '')

tumour_bulk_samples = manifest[manifest$case.id == patient & manifest$keep == "Y" & manifest$purity_trunk >= 0.4 & manifest$experimental_arm == "WGS_bulk",]$sample
non_tumour_samples = manifest[manifest$case.id == patient & manifest$keep == "Y" & manifest$purity_trunk == 0 & manifest$picard_median_coverage >= 30,]$sample #high depth so that low depth with low level tumour purity not misidentified as shared

#get variants in trunk
vars = read.table(paste0('/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/03_DNA_analyses/nd_dpclust/01_Input/', patient, '/md_out/', patient, '_allClusterassignmentsFromParallelRuns.txt'), header = T, sep = '\t', stringsAsFactors = F)
cluster_vars = vars[vars$cluster.no == truncal_cluster,]
cluster_vars$ID = paste0(cluster_vars$chr, "_", cluster_vars$pos)

truncal_vars = paste0("chr", cluster_vars$ID)
length(truncal_vars) #699

#read in VAF dataframe
VAF = read.table(paste0("/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data_mk2/merged_final_subs/", patient, "_VAF_final_merged_snps.txt"), header = T, sep = "\t", stringsAsFactors = F)

VAF_filter = VAF
VAF_filter$ID = paste0(unlist(lapply(strsplit(row.names(VAF_filter), '_'), "[[", 1)), "_", unlist(lapply(strsplit(row.names(VAF_filter), '_'), "[[", 2)))

trunc_normal_VAF = VAF[VAF_filter$ID %in% truncal_vars, non_tumour_samples]
nrow(trunc_normal_VAF) #699
put_ee_vars = row.names(trunc_normal_VAF)[rowSums(trunc_normal_VAF > 0.1) > 0] #greater than

#intersection method
sub.files = list.files("/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/02_DNA_processing/variant_calls/subs/07_indel", ".standard.gnomAD.AF.0.001.PON.read.direction.low_cov.bq.mq.binom.bbinom.5bp.indel.filtered.txt", full.names = T)
patient.sub.files = sub.files[grepl(patient, sub.files)]

manifest_flt = manifest[manifest$keep == "Y" & manifest$case.id == patient,]

#read in patient subs
all_subs = data.frame()
for(file in patient.sub.files){
  
  data = read.table(file, header = T, sep = '\t', stringsAsFactors = F, colClasses = c("character", "character", "character", "numeric", "character", "character", "numeric", "numeric", "numeric", "character", "character", "character", "character", "character", "character"))
  all_subs = rbind(all_subs, data)
}

all_subs$mutID = paste(all_subs$Chr, all_subs$Position, all_subs$Wildtype, all_subs$Mutant, sep = "_")
all_subs_flt = all_subs[all_subs$Sample %in% manifest_flt$sample,]

intersect_trunk_vars = names(table(all_subs_flt[all_subs_flt$Sample %in% tumour_bulk_samples,]$mutID)[table(all_subs_flt[all_subs_flt$Sample %in% tumour_bulk_samples,]$mutID) == length(tumour_bulk_samples)])
intersect_root_vars = unique(all_subs_flt[all_subs_flt$mutID %in% intersect_trunk_vars & all_subs_flt$Sample %in% non_tumour_samples,]$mutID)

aggregate(VAF ~ mutID, all_subs_flt[all_subs_flt$mutID %in% intersect_root_vars & all_subs_flt$Sample %in% non_tumour_samples,], max) #none seen at a max VAF <0.1
#mutID       VAF
#chr1_162729402_A_G 0.6000000
#chr6_14609040_A_G 0.5142857
#chr7_157223062_G_A 0.6296296
#chr7_158935187_G_A 0.1315789
#chr7_159159464_G_A 0.1282051
#chr9_132072514_C_A 0.6923077

all_put_root_vars = sort(unique(c(put_ee_vars, intersect_root_vars)))

write.table(VAF[all_put_root_vars, ], paste0(patient, "_put_tum_root_snp_combined_method.txt"), col.names = T, row.names = T, sep = "\t", quote = F)
