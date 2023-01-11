setwd("/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/02_DNA_processing/variant_calls/indels/04_bbinom")

manifest = read.table('/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/00_supplementary_data/wgs_paeds_pm_sample_metadata_20221015.txt', header = T, sep = '\t', stringsAsFactors = F, quote = '')

manifest_flt = manifest[manifest$keep == "Y",] #4 samples failed QC, not included in study

indel_files = list.files(".", ".pindel.HP.9.read.direction.germline.beta.binom.unk.amb.min.depth.filtered.txt")

all_indels = data.frame()
for(file in indel_files){
  
  sample = unlist(strsplit(file, ".pindel.HP.9.read.direction.germline.beta.binom.unk.amb.min.depth.filtered.txt"))
  
  if(file.exists(file)){
    
    data = read.table(file, header = T, sep = '\t', stringsAsFactors = F, colClasses = c("character", "character", "character", "numeric", "character", "character", "numeric", "numeric", "numeric", "character", "character", "character", "character", "character", "character"))
    print(sample)
    all_indels = rbind(all_indels, data)
    
  } 
}

#bulk biopsy from which multiple LCM and bulk sequences were derived
all_indels$bulk_sample = unlist(lapply(strsplit(all_indels$Sample, "_"), "[[", 1))
all_indels = all_indels[all_indels$Sample %in% manifest_flt$sample,]

table(all_indels$bulk_sample)
#all_indels[all_indels$bulk_sample == "PD51123x3",]$bulk_sample = "PD51123x"
#all_indels[all_indels$bulk_sample == "PD51123aa3",]$bulk_sample = "PD51123aa"
#all_indels[all_indels$bulk_sample == "PD51123r3",]$bulk_sample = "PD51123r"
#all_indels[all_indels$bulk_sample == "PD51122ae2",]$bulk_sample = "PD51122ae"

#add new columns
all_indels$mut_status = NA
all_indels$pure_tumour_mut = F
all_indels$mutID = paste(all_indels$Chr, all_indels$Position, all_indels$Wildtype, all_indels$Mutant, sep = "_")

#single sample only called 
pd50297_private_muts = names(table(all_indels[all_indels$Case_ID == "PD50297",]$mutID))[table(all_indels[all_indels$Case_ID == "PD50297",]$mutID) == 1]
pd51122_private_muts = names(table(all_indels[all_indels$Case_ID == "PD51122",]$mutID))[table(all_indels[all_indels$Case_ID == "PD51122",]$mutID) == 1]
pd51123_private_muts = names(table(all_indels[all_indels$Case_ID == "PD51123",]$mutID))[table(all_indels[all_indels$Case_ID == "PD51123",]$mutID) == 1]

#bulk only called, anything taken from the original bulk sample
bulk_sample.list = list()
for(bulk in unique(all_indels$bulk_sample)){
  
  bulk_sample.list[[bulk]] = unique(all_indels[all_indels$bulk_sample == bulk,]$mutID)
  
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

pd50297_tumour_muts = unique(all_indels[all_indels$Sample %in% pd50297_tumour_samples,]$mutID)
pd51122_tumour_muts = unique(all_indels[all_indels$Sample %in% pd51122_tumour_samples,]$mutID)
pd51123_tumour_muts = unique(all_indels[all_indels$Sample %in% pd51123_tumour_samples,]$mutID)

#add annotations
all_indels[all_indels$mutID %in% pd50297_private_muts & all_indels$Case_ID == "PD50297",]$mut_status = "sample"
all_indels[all_indels$mutID %in% pd51122_private_muts & all_indels$Case_ID == "PD51122",]$mut_status = "sample"
all_indels[all_indels$mutID %in% pd51123_private_muts & all_indels$Case_ID == "PD51123",]$mut_status = "sample"

all_indels[all_indels$mutID %in% pd50297_single_bulk_muts & all_indels$Case_ID == "PD50297" & !all_indels$mutID %in% pd50297_private_muts,]$mut_status = "bulk"
all_indels[all_indels$mutID %in% pd51122_single_bulk_muts & all_indels$Case_ID == "PD51122" & !all_indels$mutID %in% pd51122_private_muts,]$mut_status = "bulk"
all_indels[all_indels$mutID %in% pd51123_single_bulk_muts & all_indels$Case_ID == "PD51123" & !all_indels$mutID %in% pd51123_private_muts,]$mut_status = "bulk"

all_indels[all_indels$mutID %in% pd50297_multi_bulk_muts & all_indels$Case_ID == "PD50297",]$mut_status = "multi"
all_indels[all_indels$mutID %in% pd51122_multi_bulk_muts & all_indels$Case_ID == "PD51122",]$mut_status = "multi"
all_indels[all_indels$mutID %in% pd51123_multi_bulk_muts & all_indels$Case_ID == "PD51123",]$mut_status = "multi"

#all_indels[is.na(all_indels$mut_status),] none
#table(all_indels$mut_status)

all_indels[all_indels$mutID %in% pd50297_tumour_muts & all_indels$Case_ID == "PD50297",]$pure_tumour_mut = T
all_indels[all_indels$mutID %in% pd51122_tumour_muts & all_indels$Case_ID == "PD51122",]$pure_tumour_mut = T
all_indels[all_indels$mutID %in% pd51123_tumour_muts & all_indels$Case_ID == "PD51123",]$pure_tumour_mut = T
#table(all_indels$pure_tumour_mut)

write.table(all_indels, "all_sample_indels_distribution_annotated_20221206.txt", col.names = T, row.names = F, sep = "\t", quote = F)
