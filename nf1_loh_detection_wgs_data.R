library(ggplot2)
library(cowplot)
library(plyr)
library(zoo)

setwd("/lustre/scratch126/casm/team274sb/to3/paed_brain/autopsy/wgs_data/03_DNA_analyses/PD51122_phasing/chr17_allelecount")
options(stringsAsFactors = F)

allelecount_files = list.files(".", "_chr17_alleleFrequencies.txt") #per sample output from alleleCounter for SNPs on chromosome 17
phased_chr17 = read.table("PD51122m_chr17_hetSNPs_phased.txt", header = T, comment.char = "", sep = "\t") #identified heterozygous SNPs in tumour sample used for phasing other samples against
row.names(phased_chr17) = phased_chr17$coord

manifest_wgs = read.table('/lustre/scratch126/casm/team274sb/to3/paed_brain/autopsy_reanalysis/wgs_data/00_supplementary_data/wgs_paeds_pm_sample_metadata_20240225.txt', header = T, sep = '\t', stringsAsFactors = F, quote = '')
manifest_wgs$exp_tum_allele_baf = NA
manifest_wgs$obs_tum_allele_baf = NA
manifest_wgs$exp.tum.baf.pval = NA

for(file in allelecount_files){
  
  sample = unlist(strsplit(file, "_chr17_alleleFrequencies.txt"))
  print(sample)
  
  data = read.table(file, sep = "\t", comment.char = "", header = T)
  data$coord = paste0(data$X.CHR, "_", data$POS)
  data_flt = data[data$coord %in% phased_chr17$coord,]
  data_flt$tum_allele = paste0("Count_", phased_chr17[data_flt$coord, ]$tum_allele)
  data_flt$normal_allele = paste0("Count_", phased_chr17[data_flt$coord, ]$normal_allele)
  
  data_flt$tum_count = as.numeric(apply(data_flt, 1, function(x) x[x[9]]))
  data_flt$normal_count = as.numeric(apply(data_flt, 1, function(x) x[x[10]]))
  
  data_flt = data_flt[data_flt$Good_depth >= 10,] #missing regions for WES
  
  data_flt$tum_VAF = data_flt$tum_count / data_flt$Good_depth
  data_flt$normal_VAF = data_flt$normal_count / data_flt$Good_depth
  data_flt$tum_VAF_rollmean = rollmean(data_flt$tum_VAF, 50, fill = NA)
  data_flt$normal_VAF_rollmean = rollmean(data_flt$normal_VAF, 50, fill = NA)
  
  data_flt_nf1 = data_flt[data_flt$POS >= 31094977 - 1.5e5 & data_flt$POS <= 31377675 + 1.5e5, ]
  data_flt_nf1$tum_VAF_rollmean = rollmean(data_flt_nf1$tum_VAF, 10, fill = NA)
  data_flt_nf1$normal_VAF_rollmean = rollmean(data_flt_nf1$normal_VAF, 10, fill = NA)
  
  p1 = ggplot(data_flt) +
    geom_point(mapping = aes(x = POS, y = Good_depth), pch = ".") +
    theme_bw() +
    ylab("Good depth") +
    xlab("") +
    geom_vline(xintercept = c(31094977, 31377675)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm"),
          plot.title = element_text(hjust = 0.5)) +
    ylim(0, round_any(mean(data_flt$Good_depth) + 6 * sd(data_flt$Good_depth), 100, f = ceiling)) +
    ggtitle(paste0(sample))
  p2 = ggplot(data_flt) +
    geom_point(mapping = aes(x = POS, y = normal_VAF), col = "blue", pch = ".", alpha = 0.1) +
    geom_point(mapping = aes(x = POS, y = tum_VAF), col = "red", pch = ".", alpha = 0.1) +
    geom_vline(xintercept = c(31094977, 31377675)) +
    ylab("BAF") +
    xlab("Chromosome 17") +
    theme_bw() +
    theme(plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")) +
    ylim(c(0, 1))
  p2.5 = ggplot(data_flt) +
    geom_line(mapping = aes(x = POS, y = tum_VAF_rollmean), col = "red", alpha = 0.5) +
    geom_line(mapping = aes(x = POS, y = normal_VAF_rollmean), col = "blue", alpha = 0.5) +
    geom_vline(xintercept = c(31094977, 31377675)) +
    ylab("BAF") +
    xlab("Chromosome 17") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylim(c(0, 1)) +
    ggtitle(paste0(sample))
  
  png(paste0("/lustre/scratch126/casm/team274sb/to3/paed_brain/autopsy_reanalysis/wgs_data/03_DNA_analyses/PD51122_phasing/phased_plots/", sample, "_phased_chr17_hetSNPs.png"), width = 12, height = 8, units = "cm", res = 800)
  print(plot_grid(p1, p2, align = "v", ncol = 1))
  dev.off()
  
  png(paste0("/lustre/scratch126/casm/team274sb/to3/paed_brain/autopsy_reanalysis/wgs_data/03_DNA_analyses/PD51122_phasing/phased_plots/", sample, "_phased_chr17_smoothed_BAF_hetSNPs.png"), width = 12, height = 8, units = "cm", res = 800)
  print(p2.5)
  dev.off()
  
  #tumour has 2+0 at NF1
  if(sample %in% manifest_wgs[manifest_wgs$keep == "Y",]$sample){
    
    exp_tum_allele_baf = manifest_wgs[manifest_wgs$sample == sample,]$purity_trunk + (1 - manifest_wgs[manifest_wgs$sample == sample,]$purity_trunk) * 0.5
    
    nf1_tum_allele_count = sum(data_flt[data_flt$POS >= 31094977 & data_flt$POS <= 31377675, ]$tum_count)
    nf1_total_depth_count = sum(data_flt[data_flt$POS >= 31094977 & data_flt$POS <= 31377675, ]$Good_depth)
    
    obs_tum_allele_baf = nf1_tum_allele_count / nf1_total_depth_count
    
    manifest_wgs[manifest_wgs$sample == sample,]$exp_tum_allele_baf = exp_tum_allele_baf
    manifest_wgs[manifest_wgs$sample == sample,]$obs_tum_allele_baf = obs_tum_allele_baf
    manifest_wgs[manifest_wgs$sample == sample,]$exp.tum.baf.pval = binom.test(nf1_tum_allele_count, nf1_total_depth_count, p = exp_tum_allele_baf)$p.value
    
    p3 = ggplot(data_flt_nf1) +
      geom_point(mapping = aes(x = POS, y = Good_depth), pch = 19) +
      theme_bw() +
      ylab("Good depth") +
      xlab("") +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm"),
            plot.title = element_text(hjust = 0.5)) +
      ylim(0, round_any(mean(data_flt$Good_depth) + 6 * sd(data_flt$Good_depth), 100, f = ceiling)) +
      geom_vline(xintercept = c(31094977, 31377675)) +
      geom_vline(xintercept = 31230383, col = "orange") +
      ggtitle(paste0(sample))
    
    p4 = ggplot(data_flt_nf1) +
      geom_point(mapping = aes(x = POS, y = normal_VAF), col = "dodgerblue", pch = 20, alpha = 0.15) +
      geom_point(mapping = aes(x = POS, y = tum_VAF), col = "tomato", pch = 20, alpha = 0.15) +
      ylab("BAF") +
      xlab("NF1 genomic coordinate") +
      theme_bw() +
      theme(plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")) +
      ylim(c(0, 1)) +
      geom_line(mapping = aes(x = POS, y = tum_VAF_rollmean), col = "red") +
      geom_line(mapping = aes(x = POS, y = normal_VAF_rollmean), col = "blue") +
      geom_vline(xintercept = c(31094977, 31377675)) +
      geom_vline(xintercept = 31230383, col = "orange") +
      geom_hline(yintercept = exp_tum_allele_baf, lty = 2)
    
  } else {
    
    p3 = ggplot(data_flt_nf1) +
      geom_point(mapping = aes(x = POS, y = Good_depth), pch = 19) +
      theme_bw() +
      ylab("Good depth") +
      xlab("") +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.margin = unit(c(0.5, 0.5, 0, 0.5), "cm"),
            plot.title = element_text(hjust = 0.5)) +
      ylim(0, round_any(mean(data_flt$Good_depth) + 6 * sd(data_flt$Good_depth), 100, f = ceiling)) +
      geom_vline(xintercept = c(31094977, 31377675)) +
      geom_vline(xintercept = 31230383, col = "orange") +
      ggtitle(paste0(sample))
    
    p4 = ggplot(data_flt_nf1) +
      geom_point(mapping = aes(x = POS, y = normal_VAF), col = "dodgerblue", pch = 20, alpha = 0.15) +
      geom_point(mapping = aes(x = POS, y = tum_VAF), col = "tomato", pch = 20, alpha = 0.15) +
      ylab("BAF") +
      xlab("NF1 genomic coordinate") +
      theme_bw() +
      theme(plot.margin = unit(c(0, 0.5, 0.5, 0.5), "cm")) +
      ylim(c(0, 1)) +
      geom_line(mapping = aes(x = POS, y = tum_VAF_rollmean), col = "red") +
      geom_line(mapping = aes(x = POS, y = normal_VAF_rollmean), col = "blue") +
      geom_vline(xintercept = c(31094977, 31377675)) +
      geom_vline(xintercept = 31230383, col = "orange")
    
  }
  
  png(paste0("/lustre/scratch126/casm/team274sb/to3/paed_brain/autopsy_reanalysis/wgs_data/03_DNA_analyses/PD51122_phasing/phased_plots/", sample, "_phased_nf1_hetSNPs.png"), width = 12, height = 8, units = "cm", res = 800)
  print(plot_grid(p3, p4, align = "v", ncol = 1))
  dev.off()
  
}

manifest_wgs$exp.tum.baf.qval = NA
manifest_wgs[!is.na(manifest_wgs$exp_tum_allele_baf),]$exp.tum.baf.qval = p.adjust(manifest_wgs[!is.na(manifest_wgs$exp_tum_allele_baf),]$exp.tum.baf.pval, "BH")

write.table(manifest_wgs, '/lustre/scratch126/casm/team274sb/to3/paed_brain/autopsy_reanalysis/wgs_data/00_supplementary_data/wgs_paeds_pm_sample_metadata_20240225.txt', sep = '\t', quote = F, col.names = T, row.names = F)