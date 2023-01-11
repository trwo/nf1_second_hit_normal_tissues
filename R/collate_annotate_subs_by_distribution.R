setwd("/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/02_DNA_processing/variant_calls/subs/07_indel")

manifest = read.table('/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/00_supplementary_data/wgs_paeds_pm_sample_metadata_20221015.txt', header = T, sep = '\t', stringsAsFactors = F, quote = '')

manifest_flt = manifest[manifest$keep == "Y",] #4 samples failed QC, not included in study

all_subs = data.frame()
for(sample in manifest_flt$sample){
  data = read.table(paste0(sample, ".standard.gnomAD.AF.0.001.PON.read.direction.low_cov.bq.mq.binom.bbinom.5bp.indel.filtered.tier.annot.txt"), header = T, sep = '\t', stringsAsFactors = F)
  print(unique(data$Sample))
  all_subs = rbind(all_subs, data)
}

#bulk biopsy from which multiple LCM and bulk sequences were derived
all_subs$bulk_sample = unlist(lapply(strsplit(all_subs$Sample, "_"), "[[", 1))
all_subs[all_subs$bulk_sample == "PD51123x3",]$bulk_sample = "PD51123x"
all_subs[all_subs$bulk_sample == "PD51123aa3",]$bulk_sample = "PD51123aa"
all_subs[all_subs$bulk_sample == "PD51123r3",]$bulk_sample = "PD51123r"
all_subs[all_subs$bulk_sample == "PD51122ae2",]$bulk_sample = "PD51122ae"

#add new columns
all_subs$mut_status = NA
all_subs$pure_tumour_mut = F

#single sample only called 
pd50297_private_muts = names(table(all_subs[all_subs$Case_ID == "PD50297",]$mutID))[table(all_subs[all_subs$Case_ID == "PD50297",]$mutID) == 1]
pd51122_private_muts = names(table(all_subs[all_subs$Case_ID == "PD51122",]$mutID))[table(all_subs[all_subs$Case_ID == "PD51122",]$mutID) == 1]
pd51123_private_muts = names(table(all_subs[all_subs$Case_ID == "PD51123",]$mutID))[table(all_subs[all_subs$Case_ID == "PD51123",]$mutID) == 1]

#bulk only called, anything taken from the original bulk sample
bulk_sample.list = list()
for(bulk in unique(all_subs$bulk_sample)){
  
  bulk_sample.list[[bulk]] = unique(all_subs[all_subs$bulk_sample == bulk,]$mutID)
  
}

bulk_sample.df = data.frame(bulk_sample = rep(names(bulk_sample.list), as.numeric(unlist(lapply(bulk_sample.list, length)))),
                            mutID = as.character(unlist(bulk_sample.list)))
bulk_sample.df$case.id = substr(bulk_sample.df$bulk_sample, 0, 7)

pd50297_single_bulk_muts = names(table(bulk_sample.df[bulk_sample.df$case.id == "PD50297",]$mutID))[table(bulk_sample.df[bulk_sample.df$case.id == "PD50297",]$mutID) == 1]
pd50297_multi_bulk_muts = names(table(bulk_sample.df[bulk_sample.df$case.id == "PD50297",]$mutID))[table(bulk_sample.df[bulk_sample.df$case.id == "PD50297",]$mutID) > 1]
pd51122_single_bulk_muts = names(table(bulk_sample.df[bulk_sample.df$case.id == "PD51122",]$mutID))[table(bulk_sample.df[bulk_sample.df$case.id == "PD51122",]$mutID) == 1]
pd51122_multi_bulk_muts = names(table(bulk_sample.df[bulk_sample.df$case.id == "PD51122",]$mutID))[table(bulk_sample.df[bulk_sample.df$case.id == "PD51122",]$mutID) > 1]
pd51123_single_bulk_muts = names(table(bulk_sample.df[bulk_sample.df$case.id == "PD51123",]$mutID))[table(bulk_sample.df[bulk_sample.df$case.id == "PD51123",]$mutID) == 1]
pd51123_multi_bulk_muts = names(table(bulk_sample.df[bulk_sample.df$case.id == "PD51123",]$mutID))[table(bulk_sample.df[bulk_sample.df$case.id == "PD51123",]$mutID) > 1]

#pure tumour variants
pd50297_tumour_samples = manifest[manifest$keep == "Y" & manifest$purity_trunk >= 0.4 & manifest$case.id == "PD50297",]$sample
pd51122_tumour_samples = manifest[manifest$keep == "Y" & manifest$purity_trunk >= 0.4 & manifest$case.id == "PD51122",]$sample
pd51123_tumour_samples = manifest[manifest$keep == "Y" & manifest$purity_trunk >= 0.4 & manifest$case.id == "PD51123",]$sample

pd50297_tumour_muts = unique(all_subs[all_subs$Sample %in% pd50297_tumour_samples,]$mutID)
pd51122_tumour_muts = unique(all_subs[all_subs$Sample %in% pd51122_tumour_samples,]$mutID)
pd51123_tumour_muts = unique(all_subs[all_subs$Sample %in% pd51123_tumour_samples,]$mutID)

#add annotations
all_subs[all_subs$mutID %in% pd50297_private_muts & all_subs$Case_ID == "PD50297",]$mut_status = "sample"
all_subs[all_subs$mutID %in% pd51122_private_muts & all_subs$Case_ID == "PD51122",]$mut_status = "sample"
all_subs[all_subs$mutID %in% pd51123_private_muts & all_subs$Case_ID == "PD51123",]$mut_status = "sample"

all_subs[all_subs$mutID %in% pd50297_single_bulk_muts & all_subs$Case_ID == "PD50297" & !all_subs$mutID %in% pd50297_private_muts,]$mut_status = "bulk"
all_subs[all_subs$mutID %in% pd51122_single_bulk_muts & all_subs$Case_ID == "PD51122" & !all_subs$mutID %in% pd51122_private_muts,]$mut_status = "bulk"
all_subs[all_subs$mutID %in% pd51123_single_bulk_muts & all_subs$Case_ID == "PD51123" & !all_subs$mutID %in% pd51123_private_muts,]$mut_status = "bulk"

all_subs[all_subs$mutID %in% pd50297_multi_bulk_muts & all_subs$Case_ID == "PD50297",]$mut_status = "multi"
all_subs[all_subs$mutID %in% pd51122_multi_bulk_muts & all_subs$Case_ID == "PD51122",]$mut_status = "multi"
all_subs[all_subs$mutID %in% pd51123_multi_bulk_muts & all_subs$Case_ID == "PD51123",]$mut_status = "multi"

#all_subs[is.na(all_subs$mut_status),] none
#table(all_subs$mut_status)

all_subs[all_subs$mutID %in% pd50297_tumour_muts & all_subs$Case_ID == "PD50297",]$pure_tumour_mut = T
all_subs[all_subs$mutID %in% pd51122_tumour_muts & all_subs$Case_ID == "PD51122",]$pure_tumour_mut = T
all_subs[all_subs$mutID %in% pd51123_tumour_muts & all_subs$Case_ID == "PD51123",]$pure_tumour_mut = T
#table(all_subs$pure_tumour_mut)

write.table(all_subs, "all_sample_subs_distribution_annotated_20221201.txt", col.names = T, row.names = F, sep = "\t", quote = F)
