
### FUNCTION ###

mean_no_outliers <- function(x, remove_na = TRUE) { #early embryonic variants might skew overall calculation
  qnt <- quantile(x, probs = c(.25, .75), na.rm = remove_na)
  val <- 1.5 * IQR(x, na.rm = remove_na)
  y <- x
  y[x < (qnt[1] - val)] <- NA
  y[x > (qnt[2] + val)] <- NA
  z = mean(y, na.rm = remove_na)
  return(z)
}

### PD51122 ###
patient = "PD51122"
truncal_cluster = "6_9" #merged

manifest = read.table('/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/00_supplementary_data/wgs_paeds_pm_sample_metadata_20220717.txt', header = T, sep = '\t', stringsAsFactors = F, quote = '')
samples = manifest[manifest$case.id == patient & manifest$battenberg_purity >= 0.4 & !is.na(manifest$battenberg_purity), ]$sample

#get variants in trunk
vars = read.table(paste0('/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/03_DNA_analyses/nd_dpclust/01_Input/', patient, '/md_out/', patient, '_allClusterassignmentsFromParallelRuns_clusters_reviewed.txt'), header = T, sep = '\t', stringsAsFactors = F) #reviewed and merged clusters

cluster_vars = vars[vars$cluster.no == truncal_cluster,]
cluster_vars$ID = paste0(cluster_vars$chr, "_", cluster_vars$pos)

#identify which truncal variants cover two chromomes copies in 2+0 sites
dp_files = list.files(paste0('/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/03_DNA_analyses/nd_dpclust/01_Input/', patient), "allDirichletProcessInfo_fromMultipleSamplesWithoutChromosomeLosses.txt", full.names = T)

trunc_pre_dup = list()
for(file in dp_files){
  data = read.table(file, header = T, sep = '\t', stringsAsFactors = F)
  sample = unlist(strsplit(basename(file), "_allDirichletProcessInfo_fromMultipleSamplesWithoutChromosomeLosses.txt"))
  data$ID = paste0(data$chr, "_", data$pos)
  trunc_pre_dup[[sample]] = data[data$ID %in% cluster_vars$ID & data$subclonal.CN == 2 & is.na(data$nMaj2) & data$nMin1 == 0 & data$no.chrs.bearing.mut == 2,]$ID
}

hi_conf_loci = Reduce(intersect, trunc_pre_dup)
hi_conf_loci = paste0("chr", hi_conf_loci)
length(hi_conf_loci) #60

#read in VAF dataframe containing every substitution (rows) and its VAF across all samples (columns) from the case
VAF = read.table(paste0("/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data_mk2/merged_final_subs/", patient, "_VAF_final_merged_snps.txt"), header = T, sep = "\t", stringsAsFactors = F)

VAF_filter = VAF
VAF_filter$ID = paste0(unlist(lapply(strsplit(row.names(VAF_filter), '_'), "[[", 1)), "_", unlist(lapply(strsplit(row.names(VAF_filter), '_'), "[[", 2)))

tum_infiltration_est = data.frame(apply(VAF[VAF_filter$ID %in% hi_conf_loci, ], 2, function(x) mean_no_outliers(x)))
tum_infiltration_est$PD = row.names(tum_infiltration_est)
colnames(tum_infiltration_est)[1] = c("purity_trunc")
tum_infiltration_est = tum_infiltration_est[, c(2,1)]

manifest$purity_trunk = NA

for(i in 1:nrow(manifest)){
  if(manifest$case.id[i] == patient & manifest$keep[i] == "Y"){
    manifest$purity_trunk[i] = tum_infiltration_est[tum_infiltration_est$PD == manifest$sample[i], ]$purity_trunc
  }
}

### PD50297 ###
patient = "PD50297"
truncal_cluster = 9

samples = manifest[manifest$case.id == patient & manifest$battenberg_purity >= 0.4 & !is.na(manifest$battenberg_purity), ]$sample

#get variants in trunk
vars = read.table(paste0('/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/03_DNA_analyses/nd_dpclust/01_Input/', patient, '/md_out/', patient, '_allClusterassignmentsFromParallelRuns.txt'), header = T, sep = '\t', stringsAsFactors = F) #from original as root unchanged

cluster_vars = vars[vars$cluster.no == truncal_cluster,]
cluster_vars$ID = paste0(cluster_vars$chr, "_", cluster_vars$pos)

#identify which truncal variants occur at 1+1 sites
dp_files = list.files(paste0('/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/03_DNA_analyses/nd_dpclust/01_Input/', patient), "allDirichletProcessInfo_fromMultipleSamplesWithoutChromosomeLosses.txt", full.names = T)

trunc_pre_dup = list()
for(file in dp_files){
  data = read.table(file, header = T, sep = '\t', stringsAsFactors = F)
  sample = unlist(strsplit(basename(file), "_allDirichletProcessInfo_fromMultipleSamplesWithoutChromosomeLosses.txt"))
  data$ID = paste0(data$chr, "_", data$pos)
  trunc_pre_dup[[sample]] = data[data$ID %in% cluster_vars$ID & data$subclonal.CN == 2 & is.na(data$nMaj2) & data$nMin1 == 1 & data$nMaj1 == 1,]$ID
}

hi_conf_loci = Reduce(intersect, trunc_pre_dup)
hi_conf_loci = paste0("chr", hi_conf_loci)
length(hi_conf_loci) #385

#read in VAF dataframe
VAF = read.table(paste0("/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data_mk2/merged_final_subs/", patient, "_VAF_final_merged_snps.txt"), header = T, sep = "\t", stringsAsFactors = F)

VAF_filter = VAF
VAF_filter$ID = paste0(unlist(lapply(strsplit(row.names(VAF_filter), '_'), "[[", 1)), "_", unlist(lapply(strsplit(row.names(VAF_filter), '_'), "[[", 2)))

tum_infiltration_est = data.frame(apply(VAF[VAF_filter$ID %in% hi_conf_loci, ], 2, function(x) 2 * mean_no_outliers(x))) #x2 at 1+1 sites for cellularity
tum_infiltration_est$PD = row.names(tum_infiltration_est)
colnames(tum_infiltration_est)[1] = c("purity_trunc")
tum_infiltration_est = tum_infiltration_est[, c(2,1)]

for(i in 1:nrow(manifest)){
  if(manifest$case.id[i] == patient & manifest$keep[i] == "Y"){
    manifest$purity_trunk[i] = tum_infiltration_est[tum_infiltration_est$PD == manifest$sample[i], ]$purity_trunc
  }
}

### PD51123 ###
patient = "PD51123"
truncal_cluster = 2

samples = manifest[manifest$case.id == patient & manifest$battenberg_purity >= 0.4 & !is.na(manifest$battenberg_purity), ]$sample

#get variants in trunk
vars = read.table(paste0('/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/03_DNA_analyses/nd_dpclust/01_Input/', patient, '/md_out/', patient, '_allClusterassignmentsFromParallelRuns.txt'), header = T, sep = '\t', stringsAsFactors = F) #from original as root unchanged
cluster_vars = vars[vars$cluster.no == truncal_cluster,]
cluster_vars$ID = paste0(cluster_vars$chr, "_", cluster_vars$pos)

#identify which truncal variants occur at 1+1 sites
dp_files = list.files(paste0('/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/03_DNA_analyses/nd_dpclust/01_Input/', patient), "allDirichletProcessInfo_fromMultipleSamplesWithoutChromosomeLosses.txt", full.names = T)

trunc_pre_dup = list()
for(file in dp_files){
  data = read.table(file, header = T, sep = '\t', stringsAsFactors = F)
  sample = unlist(strsplit(basename(file), "_allDirichletProcessInfo_fromMultipleSamplesWithoutChromosomeLosses.txt"))
  data$ID = paste0(data$chr, "_", data$pos)
  trunc_pre_dup[[sample]] = data[data$ID %in% cluster_vars$ID & data$subclonal.CN == 2 & is.na(data$nMaj2) & data$nMin1 == 1 & data$nMaj1 == 1,]$ID
}

hi_conf_loci = Reduce(intersect, trunc_pre_dup)
hi_conf_loci = paste0("chr", hi_conf_loci)
length(hi_conf_loci) #260

#read in VAF dataframe
VAF = read.table(paste0("/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data_mk2/merged_final_subs/", patient, "_VAF_final_merged_snps.txt"), header = T, sep = "\t", stringsAsFactors = F)

VAF_filter = VAF
VAF_filter$ID = paste0(unlist(lapply(strsplit(row.names(VAF_filter), '_'), "[[", 1)), "_", unlist(lapply(strsplit(row.names(VAF_filter), '_'), "[[", 2)))

tum_infiltration_est = data.frame(apply(VAF[VAF_filter$ID %in% hi_conf_loci, ], 2, function(x) 2 * mean_no_outliers(x))) #x2 at 1+1 sites for cellularity
tum_infiltration_est$PD = row.names(tum_infiltration_est)
colnames(tum_infiltration_est)[1] = c("purity_trunc")
tum_infiltration_est = tum_infiltration_est[, c(2,1)]

for(i in 1:nrow(manifest)){
  if(manifest$case.id[i] == patient & manifest$keep[i] == "Y"){
    manifest$purity_trunk[i] = tum_infiltration_est[tum_infiltration_est$PD == manifest$sample[i], ]$purity_trunc
  }
}

write.table(manifest, '/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/00_supplementary_data/wgs_paeds_pm_sample_metadata_20220717.txt', sep = '\t', quote = F, col.names = T, row.names = F)
