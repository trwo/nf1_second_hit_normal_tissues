
### CALCULATE ###
setwd("/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/02_DNA_processing/localised_copy_number_calling")

manifest = read.table('/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/00_supplementary_data/wgs_paeds_pm_sample_metadata_20221015.txt', header = T, sep = '\t', stringsAsFactors = F, quote = '')

#read in binned coverage files
tumour1_bins = read.table("PD51123o_lo0008.tumour.chr4.54229000.54298500.seg.bins_averaged.txt", header = F, sep = "\t", stringsAsFactors = F)
names(tumour1_bins) = c("chr", "start", "end", "winnum", "coverage")
tumour1_bins$plot_coord = apply(tumour1_bins[, 2:3], 1, mean)
tumour1_bins$sample = "PD51123o_lo0008"

tumour2_bins = read.table("PD51123o_lo0012.tumour.chr4.54229000.54298500.seg.bins_averaged.txt", header = F, sep = "\t", stringsAsFactors = F)
names(tumour2_bins) = c("chr", "start", "end", "winnum", "coverage")
tumour2_bins$plot_coord = apply(tumour2_bins[, 2:3], 1, mean)
tumour2_bins$sample = "PD51123o_lo0012"

normal_bins = read.table("PD51123o_lo0008.normal.chr4.54229000.54298500.seg.bins_averaged.txt", header = F, sep = "\t", stringsAsFactors = F)
names(normal_bins) = c("chr", "start", "end", "winnum", "coverage")
normal_bins$plot_coord = apply(normal_bins[, 2:3], 1, mean)
normal_bins$sample = "PD51123aa3"

all_cov = rbind(tumour1_bins, tumour2_bins, normal_bins)

all_cov$sample = factor(all_cov$sample, levels = c("PD51123o_lo0008", "PD51123o_lo0012", "PD51123aa3"))

#PDGFRA coordinates
min(elocs$start_position) #54229280
max(elocs$end_position) #54298245

#PDGFRA breakpoints
svs = read.table("/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/02_DNA_processing/variant_calls/structural_variants/paeds_autopsy_annotated_svs_uniq_svs_labelled_purity_annotated_20221206.txt", header = T, sep = "\t", stringsAsFactors = F)
svs[svs$gene == "PDGFRA" & svs$gene2 == "PDGFRA" & svs$sample %in% c("PD51123o_lo0008", "PD51123o_lo0012"), c(1:6, 10)]
#chro1   start1     end1      chro2   start2      end2          sample
#chr4    54253544   54253548  chr4    54268457    54268461      PD51123o_lo0008
#chr4    54267788   54267791  chr4    54272641    54272644      PD51123o_lo0012

#within deletion cn calculation
##PD51123o_lo0008
within_del_cov = mean(all_cov[all_cov$sample == "PD51123o_lo0008" & all_cov$start >= 54253544 & all_cov$end <= 54268461,]$coverage) #31.06579
outside_del_cov = mean(all_cov[all_cov$sample == "PD51123o_lo0008" & ((all_cov$start >= 54229280 & all_cov$end <= 54253544) | (all_cov$start >= 54268461 & all_cov$end <= 54298245)),]$coverage) #408.9732
tumour_purity = manifest[manifest$sample == "PD51123o_lo0008",]$purity_trunk
tumour_cn = 18 #outside deletion
normal_cn = 2

tum_cov = outside_del_cov * (tumour_purity * tumour_cn / (tumour_purity * tumour_cn + (1 - tumour_purity) * normal_cn)) #407.4018
norm_cov = outside_del_cov - tum_cov #1.571355, should be consistent across the autosomal genome
tum_cov_per_cn = tum_cov / tumour_cn
tum_cov_within = within_del_cov - norm_cov
tum_cov_within / tum_cov_per_cn #1.303135

##PD51123o_lo0012
within_del_cov = mean(all_cov[all_cov$sample == "PD51123o_lo0012" & all_cov$start >= 54267788 & all_cov$end <= 54272644,]$coverage)  #65.77267
outside_del_cov = mean(all_cov[all_cov$sample == "PD51123o_lo0012" & ((all_cov$start >= 54229280 & all_cov$end <= 54267788) | (all_cov$start >= 54272644 & all_cov$end <= 54298245)),]$coverage) #732.1705
tumour_purity = manifest[manifest$sample == "PD51123o_lo0012",]$purity_trunk
tumour_cn = 33 #outside deletion
normal_cn = 2

tum_cov = outside_del_cov * (tumour_purity * tumour_cn / (tumour_purity * tumour_cn + (1 - tumour_purity) * normal_cn)) #724.5423
norm_cov = outside_del_cov - tum_cov #7.62821, should be consistent across the autosomal genome
tum_cov_per_cn = tum_cov / tumour_cn
tum_cov_within = within_del_cov - norm_cov
tum_cov_within / tum_cov_per_cn #2.648247

### PLOT ###

#PDGFRA breakpoints
svs = read.table("/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/02_DNA_processing/variant_calls/structural_variants/paeds_autopsy_annotated_svs_uniq_svs_labelled_purity_annotated_20221206.txt", header = T, sep = "\t", stringsAsFactors = F)
svs[svs$gene == "PDGFRA" & svs$gene2 == "PDGFRA" & svs$sample %in% c("PD51123o_lo0008", "PD51123o_lo0012"), c(1:6, 10)]
#chro1   start1     end1      chro2   start2      end2          sample
#chr4    54253544   54253548  chr4    54268457    54268461      PD51123o_lo0008
#chr4    54267788   54267791  chr4    54272641    54272644      PD51123o_lo0012

#create df
est_cn_df = data.frame(sample = c("PD51123o_lo0008", "PD51123o_lo0008", "PD51123o_lo0008", "PD51123o_lo0012", "PD51123o_lo0012", "PD51123o_lo0012", "PD51123p", "PD51123aa3"), 
                       start = c(54228500, 54253544, 54268461, 54228500, 54267788, 54272644, 54228500, 54228500), 
                       end = c(54253543, 54268460, 54298501, 54267787, 54272643, 54298501, 54298501, 54298501), 
                       tot_cn = c(18, 1, 18, 33, 3, 33, 19, 2))
est_cn_df$sample = factor(est_cn_df$sample, levels = c("PD51123o_lo0008", "PD51123o_lo0012", "PD51123p", "PD51123aa3"))

p2 = ggplot(est_cn_df) +
  geom_segment(mapping = aes(x = start, xend = end, y = tot_cn, yend = tot_cn), size = 2) +
  facet_grid(sample ~ .) +
  theme_bw() +
  scale_x_continuous(expand = c(0,0), limits = c(54228500, 54298501)) +
  ylab("Total copy number") +
  xlab("Chromosome 4 position") +
  ylim(0, 40)

p1 = ggplot() +
  geom_line(mapping = aes(x = c(54229280, 54298245), y = c(0.5, 0.5)), size = 2) +
  geom_segment(elocs, mapping = aes(x = exon_chrom_start, xend = exon_chrom_end, y = 0.5, yend = 0.5), size = 10) +
  scale_x_continuous(expand = c(0,0), limits = c(54228500, 54298501)) +
  ylim(0.45, 0.55) +
  theme_void()

pdf("/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/03_DNA_analyses/draft_figures/pdgra_deletions_20221211.pdf", width = 4, height = 4, useDingbats = F)
plot_grid(p1, p2, ncol = 1, rel_heights = c(0.05, 0.95), align = "v", axis = "lr")
dev.off()

write.table(est_cn_df, "/lustre/scratch119/casm/team274sb/to3/paed_brain/autopsy/wgs_data/03_DNA_analyses/draft_figures/pdgra_deletions_data_20221222.txt", col.names = T, row.names = F, sep = "\t", quote = F)
