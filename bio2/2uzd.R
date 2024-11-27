data <- readRDS("data.rds")
genomemap <- readRDS("genomemap.rds")
samplekey <- readRDS("samplekey.rds")
study_participants<-read.csv("1-s2.0-S1872497323001151-mmc5.csv",sep = ";",
                             row.names = NULL)
library(tidyverse)
library(ggrepel)
library(RColorBrewer)
library(qqman)
library(purrr)
library(reshape)
study_participants <- dplyr::rename(study_participants, age = Actual.age)
study_participants <- dplyr::rename(study_participants, disease = D_Underlying)

study_participants <- mutate(study_participants,
                             donor = paste(study_participants$Sample,
                                           "_",
                                           study_participants$age,
                                           sep = ""))

sampleykey_new_names <- paste(samplekey$donor,"_", samplekey$sex,"_", samplekey$tissue, sep ="")
rownames(samplekey) <- sampleykey_new_names
colnames(data) <- sampleykey_new_names

study_participants$possible_depressed <- ifelse(study_participants$Manner == "Suicide", 1, 2)
depreessed_participants <- study_participants[study_participants$Manner == "Suicide",]

possible_depressed_values <- numeric()
for (i in 1:nrow(samplekey)) {
  if (samplekey$donor[i] %in% depreessed_participants$donor) {
    possible_depressed_values[i] <- 1
  } else {
    possible_depressed_values[i] <- 2
  }
}

samplekey$possible_depressed <- possible_depressed_values

########################################

for (tissue in unique(samplekey$tissue)) {
  data <- readRDS("data.rds")
  colnames(data) <- sampleykey_new_names
  samplekey_index <- grep( tissue, rownames(samplekey))
  samplekey_data <- samplekey[samplekey_index,]
  
  depressed_index <- rownames(samplekey_data[
    which(samplekey_data$possible_depressed == 1), ])
  
  non_depressed_index <- rownames(samplekey_data[
    which(samplekey_data$possible_depressed == 2), ])
  
  depressed_data <- data[, depressed_index]
  non_depressed_data <- data[, non_depressed_index]
  p_values <- matrix(NA, nrow = nrow(data), ncol = 1,
                     dimnames = list(rownames(data), c("P_value")))
  
  effect_size <- matrix(NA, nrow = nrow(data), ncol = 1,
                        dimnames = list(rownames(data), c("Effect_size")))
  rm(data)
  depressed_data <- t(depressed_data)
  non_depressed_data <- t(non_depressed_data)
  
  p_values <- mapply(function(x, y)  t.test(x,y)$p.value,
                     as.data.frame(depressed_data),
                     as.data.frame(non_depressed_data))
  p_values <- as.data.frame(p_values)
  
  
  effect_size <- mapply(function(x, y) mean(x) - mean(y),
                        as.data.frame(depressed_data),
                        as.data.frame(non_depressed_data))
  effect_size <- as.data.frame(effect_size)
  file_path<-paste0("dep_not_dep/", tissue)
  p_values_file <- paste0(file_path,"/p_values_", tissue, ".rds")
  effect_size_file <- paste0(file_path,"/effect_size_", tissue, ".rds")
  
  saveRDS(p_values, file = p_values_file)
  saveRDS(effect_size, file = effect_size_file)
  
  p_df<-p_values
  
  p_df$fdr<-p.adjust(p_df$p_values, method = "fdr",
                     n = length(p_df$p_values))
  
  p_df$bonferroni<-p.adjust(p_df$p_values, method = "bonferroni",
                      n = length(p_df$p_values))
  
  no_cor <-ggplot(p_df, aes(x = p_values)) + 
    geom_histogram(binwidth = 0.025, color = "grey", fill = "salmon") +
    geom_vline(xintercept = 0.05, color = "red", size = 0.7,
               linetype = "longdash") +
    labs(x = "P value", y = "Count",
         title = paste0("P-value distribution histogram_",tissue)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.05))
  
  
  bonferroni_cor<- ggplot(p_df, aes(x = bonferroni)) + 
    geom_histogram(binwidth = 0.025, color = "grey", fill = "salmon") +
    geom_vline(xintercept = 0.05, color = "red", size = 0.7,
               linetype = "longdash") +
    labs(x = "P value", y = "Count",
         title = paste0("Bonferroni P-value distribution histogram_",tissue)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.05))
  
  
  fdr_cor<-ggplot(p_df, aes(x = fdr)) + 
    geom_histogram(binwidth = 0.025, color = "grey", fill = "salmon") +
    geom_vline(xintercept = 0.05, color = "red", size = 0.7,
               linetype = "longdash") +
    labs(x = "P value", y = "Count",
         title = paste0("FDR P-value distribution histogram_",tissue)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.05))
  
  ggsave(paste0(file_path, "/no_cor_histogram.png"),
         plot = no_cor, width = 8, height = 6, units = "in")
  ggsave(paste0(file_path, "/bonferroni_cor_histogram.png"),
         plot = bonferroni_cor, width = 8, height = 6, units = "in")
  ggsave(paste0(file_path, "/fdr_cor_histogram.png"),
         plot = fdr_cor, width = 8, height = 6, units = "in")
  
  p_lower_05 <- sum(p_df$p_values < 0.05)
  p_lower_05_fdr <- sum(p_df$fdr < 0.05)
  p_lower_05_b <- sum(p_df$bonferroni < 0.05)
  data_df <- data.frame(no_corr = p_lower_05,
                        fdr = p_lower_05_fdr, bon = p_lower_05_b)
  write.table(data_df, file =paste0(file_path,"/count.tsv"),
              sep = "\t", row.names = FALSE)
  
}



coronary <- c("2020_127_38_M","2020_14_48_M","2020_16_77_M",
               "2020_93_74_M","2020_25_56_M","2020_23_47_F",
               "2020_138_50_M","2020_7_43_M","2020_8_66_M",
               "2020_17_78_M")
non_coronary <- c("2020_109_52_F","2020_115_52_M","2020_118_77_M",
                   "2020_120_39_M","2020_131_18_M","2020_157_41_F",
                   "2020_165_35_M","2020_147_63_M","2020_13_46_M",
                   "2020_15_56_F")

coronary_data <- data[, grep(paste(coronary, collapse = "|"),
                             colnames(data))]

non_coronary_data <- data[, grep(paste(non_coronary, collapse = "|"),
                                 colnames(data))]

p_values <- readRDS("dep_not_dep/blood/p_values_blood.rds")

for (tissue in unique(samplekey$tissue)) {
  p_values_file <- paste0("coronary/", tissue,"/p_values_", tissue, ".rds")
  p_values <- readRDS(p_values_file)
  p_df<-as.data.frame(p_values)
  
  p_df$fdr<-p.adjust(p_df$P_value, method = "fdr",
                     n = length(p_df$P_value))
  
  p_df$bonferroni<-p.adjust(p_df$P_value, method = "bonferroni",
                            n = length(p_df$P_value))
  
  no_cor <-ggplot(p_df, aes(x = P_value)) + 
    geom_histogram(binwidth = 0.025, color = "grey", fill = "salmon") +
    geom_vline(xintercept = 0.05, color = "red", size = 0.7,
               linetype = "longdash") +
    labs(x = "P value", y = "Count",
         title = paste0("P-value distribution histogram_",tissue)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.05))
  
  
  bonferroni_cor<- ggplot(p_df, aes(x = bonferroni)) + 
    geom_histogram(binwidth = 0.025, color = "grey", fill = "salmon") +
    geom_vline(xintercept = 0.05, color = "red", size = 0.7,
               linetype = "longdash") +
    labs(x = "P value", y = "Count",
         title = paste0("Bonferroni P-value distribution histogram_",tissue)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.05))
  
  
  fdr_cor<-ggplot(p_df, aes(x = fdr)) + 
    geom_histogram(binwidth = 0.025, color = "grey", fill = "salmon") +
    geom_vline(xintercept = 0.05, color = "red", size = 0.7,
               linetype = "longdash") +
    labs(x = "P value", y = "Count",
         title = paste0("FDR P-value distribution histogram_",tissue)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.05))
  
  ggsave(paste0("coronary/", tissue, "/no_cor_histogram.png"),
         plot = no_cor, width = 8, height = 6, units = "in")
  ggsave(paste0("coronary/", tissue, "/bonferroni_cor_histogram.png"),
         plot = bonferroni_cor, width = 8, height = 6, units = "in")
  ggsave(paste0("coronary/", tissue, "/fdr_cor_histogram.png"),
         plot = fdr_cor, width = 8, height = 6, units = "in")
  
  p_lower_05 <- 0
  p_lower_05_fdr <- 0
  p_lower_05_b <- 0
  for (row in 1:nrow(p_df)) {
    p_value <- p_df[row, "P_value"]
    p_value_fdr <- p_df[row, "fdr"]
    p_value_b <- p_df[row, "bonferroni"]
    if (p_value < 0.05) {
      p_lower_05 <- p_lower_05 + 1
    }
    
    if (p_value_fdr < 0.05) {
      p_lower_05_fdr <- p_lower_05_fdr + 1
    }
    if (p_value_b < 0.05) {
      p_lower_05_b <- p_lower_05_b + 1
    }
  }
  data_df <- data.frame(no_corr = p_lower_05,
                        fdr = p_lower_05_fdr, bon = p_lower_05_b)
  write.table(data_df, file =paste0("coronary/", tissue,"/count.tsv"),
              sep = "\t", row.names = FALSE)
  
}



for (tissue in unique(samplekey$tissue)) {
  data <- readRDS("data.rds")
  colnames(data) <- sampleykey_new_names
  male_index <- grep( paste0("M_", tissue), colnames(data))
  male_data <- data[, male_index]
  
  female_index <- grep( paste0("F_", tissue), colnames(data))
  female_data <- data[, female_index]

  p_values <- matrix(NA, nrow = nrow(data), ncol = 1,
          dimnames = list(rownames(data), c("P_value")))
  
  effect_size <- matrix(NA, nrow = nrow(data), ncol = 1,
            dimnames = list(rownames(data), c("Effect_size")))
  rm(data)
  male_data <- t(male_data)
  female_data <- t(female_data)
  
  p_values <- mapply(function(x, y)  t.test(x,y)$p.value,as.data.frame(male_data),
                     as.data.frame(female_data))
  p_values <- as.data.frame(p_values)
  
  effect_size <- mapply(function(x, y) mean(x) - mean(y),
            as.data.frame(male_data), as.data.frame(female_data))
  effect_size <- as.data.frame(p_values)
  file_path<-paste0("male_female/", tissue)
  p_values_file <- paste0(file_path,"/p_values_", tissue, ".rds")
  effect_size_file <- paste0(file_path,"/effect_size_", tissue, ".rds")
  
  saveRDS(p_values, file = p_values_file)
  saveRDS(effect_size, file = effect_size_file)
}

for (tissue in unique(samplekey$tissue)) {
  file_path<-paste0("male_female/", tissue)
  p_values_file <- paste0(file_path,"/p_values_", tissue, ".rds")
  p_values <- readRDS(p_values_file)
  p_df<-as.data.frame(p_values)
  
  p_df$fdr<-p.adjust(p_df$P_value, method = "fdr",
                     n = length(p_df$P_value))
  
  p_df$bonferroni<-p.adjust(p_df$P_value, method = "bonferroni",
                            n = length(p_df$P_value))
  
  no_cor <-ggplot(p_df, aes(x = P_value)) + 
    geom_histogram(binwidth = 0.025, color = "grey", fill = "salmon") +
    geom_vline(xintercept = 0.05, color = "red", size = 0.7,
               linetype = "longdash") +
    labs(x = "P value", y = "Count",
         title = paste0("P-value distribution histogram_",tissue)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.05))
  
  
  bonferroni_cor<- ggplot(p_df, aes(x = bonferroni)) + 
    geom_histogram(binwidth = 0.025, color = "grey", fill = "salmon") +
    geom_vline(xintercept = 0.05, color = "red", size = 0.7,
               linetype = "longdash") +
    labs(x = "P value", y = "Count",
         title = paste0("Bonferroni P-value distribution histogram_",tissue)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.05))
  
  
  fdr_cor<-ggplot(p_df, aes(x = fdr)) + 
    geom_histogram(binwidth = 0.025, color = "grey", fill = "salmon") +
    geom_vline(xintercept = 0.05, color = "red", size = 0.7,
               linetype = "longdash") +
    labs(x = "P value", y = "Count",
         title = paste0("FDR P-value distribution histogram_",tissue)) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.05))
  
  ggsave(paste0(file_path, tissue, "/no_cor_histogram.png"),
         plot = no_cor, width = 8, height = 6, units = "in")
  ggsave(paste0(file_path, tissue, "/bonferroni_cor_histogram.png"),
         plot = bonferroni_cor, width = 8, height = 6, units = "in")
  ggsave(paste0(file_path, tissue, "/fdr_cor_histogram.png"),
         plot = fdr_cor, width = 8, height = 6, units = "in")
  
  p_lower_05 <- 0
  p_lower_05_fdr <- 0
  p_lower_05_b <- 0
  for (row in 1:nrow(p_df)) {
    p_value <- p_df[row, "P_value"]
    p_value_fdr <- p_df[row, "fdr"]
    p_value_b <- p_df[row, "bonferroni"]
    if (p_value < 0.05) {
      p_lower_05 <- p_lower_05 + 1
    }
    
    if (p_value_fdr < 0.05) {
      p_lower_05_fdr <- p_lower_05_fdr + 1
    }
    if (p_value_b < 0.05) {
      p_lower_05_b <- p_lower_05_b + 1
    }
  }
  data_df <- data.frame(no_corr = p_lower_05,
          fdr = p_lower_05_fdr, bon = p_lower_05_b)
  write.table(data_df, file =paste0(file_path, tissue,"/count.tsv"),
                                    sep = "\t", row.names = FALSE)
}

#Volcan plot############ 
p_value<- readRDS("dep_not_dep/kidney/p_values_kidney.rds")
effect_size<- readRDS("dep_not_dep/kidney/effect_size_kidney.rds")

p_value_df <- as.data.frame(p_value)
effect_df<-as.data.frame(effect_size)
p_value_df$effect_size<-effect_df$effect_size
p_value_df <- p_value_df %>%
  mutate(significance = case_when(effect_size > 0 & p_values <= 0.05 ~ "up",
                                  effect_size <= 0 & p_values <= 0.05 ~ "down",
                                  TRUE ~ "ns"))
sorted_df <- p_value_df[order(p_value_df$p_values), ]
sorted_index <- rownames(head(sorted_df, 3))
p_value_df$delabel <- ifelse(rownames(p_value_df) %in% sorted_index,
                      rownames(p_value_df), NA)
sorted_df <- p_value_df[order(p_value_df$effect_size), ]
sorted_index <- rownames(head(sorted_df, 1))
p_value_df[sorted_index,]$delabel <- rownames(p_value_df[sorted_index,])

ggplot(p_value_df, aes(x = effect_size, y = -log10(p_values),
                      col = significance, label= delabel)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
  theme( axis.title.y = element_text(face = "bold",
      margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
    axis.title.x = element_text(hjust = 0.5, face = "bold",
      margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
    plot.title = element_text(hjust = 0.5)
  ) +
  geom_point() +
  scale_color_manual(values = c("#00AFBB", "grey", "salmon"),
      labels = c("Downregulated", "Not significant", "Upregulated")) +
  labs(x = "Effect size", y = "-log10(p-value)",
       title = "Volcano plot, depressed/non-depressed group, kidney tissue") +
  geom_text_repel(max.overlaps = Inf)

# manhattan plot ##############################################################
chromosome_names <- split(rownames(genomemap), genomemap$chr)
count_list <- list()
total_list <- list()
percent_list <- list()
for (chromosome in names(chromosome_names)) {
  row_names <- chromosome_names[[chromosome]]
  p_values <- p_value_df[row_names, "p_values"]
  count <- sum(p_values < 0.05, na.rm = TRUE)
  total_positions <- length(row_names)
  count_list[[chromosome]] <- count
  total_list[[chromosome]] <- total_positions
  percent <- (count / total_positions) * 100
  percent_list[[chromosome]] <- percent
}
p_value_X_df <- p_value_df[chromosome_names[["chrX"]],]
genomemap_X <- genomemap[chromosome_names[["chrX"]],]
p_value_X_df$Pos <- genomemap_X$pos
p_value_X_df$gene <- genomemap_X$ucsc_refgene_name
p_value_X_df$gene <- ifelse(p_value_X_df$gene == "",
                             "No gene", p_value_X_df$gene)
p_value_X_df$gene <- sapply(strsplit(p_value_X_df$gene, ";"), `[`, 1)
ggplot(p_value_X_df, aes(x = Pos, y = -log10(p_values))) +
  geom_point(aes(color = -log10(p_values) > 6), show.legend = FALSE) +
  scale_color_manual(values = c("FALSE" = "grey",
                                "TRUE" = "skyblue")) +
  labs(x = "Position", y = "P-value", title = "Manhattan Plot") +
  theme_minimal() +
  scale_x_continuous(labels = scales::number_format()) +
  geom_text(data = subset(p_value_X_df, -log10(p_values) > 6),
            aes(label = gene),
            hjust = 1.2, vjust = 0.5, size = 3)

p_value_df$pos <- genomemap$pos
p_value_df$chr <- genomemap$chr
p_value_df$chr <- gsub("^chr", "", p_value_df$chr)
p_value_df$chr <- gsub("X", "23", p_value_df$chr)
p_value_df$chr <- gsub("Y", "24", p_value_df$chr)
p_value_df$chr <- as.numeric(p_value_df$chr)
p_value_df$snp = rownames(p_value_df)
manhattan(p_value_df, chr="chr", bp="pos",snp = "snp", p="p_values")
# Dnr profile ##################################################################
top_10_df <- head(p_value_df[order(p_value_df$p_values), ],10)
samplekey_b_li <- samplekey[rownames(samplekey$celltype=="bcell"), ]
depressed_index <- rownames(samplekey[
  which(samplekey$possible_depressed == 1), ])

non_depressed_index <- rownames(samplekey[
  which(samplekey$possible_depressed == 2), ])

depressed_data <- data[rownames(top_10_df), depressed_index]
non_depressed_data <- data[rownames(top_10_df), non_depressed_index]

depressed_data <- melt(depressed_data)
depressed_data <- mutate(depressed_data, X2 = "Depressed")

non_depressed_data <- melt(non_depressed_data)
non_depressed_data <- mutate(non_depressed_data, X2 = "Non_depressed")
mod_dnr <- rbind(depressed_data, non_depressed_data)
mod_dnr<-setNames(mod_dnr,
                  c("Position", "Group", "Modification_value"))
mod_dnr$Group <- factor(mod_dnr$Group)
pos_names <- levels(mod_dnr$Position)
genomemap_top_10 <- genomemap[pos_names,]

indices <- match(mod_dnr$Position, rownames(genomemap_top_10))
mod_dnr$genes <- genomemap_top_10[indices, ]$ucsc_refgene_name
mod_dnr$genes <- sapply(strsplit(mod_dnr$genes, ";"), `[`, 1)
mod_dnr$genes <- ifelse(is.na(mod_dnr$genes), "No gene", mod_dnr$genes)
mod_dnr$genes <- factor(mod_dnr$genes)

mod_duplicate <- which(mod_dnr$Position == "cg14520512")
mod_dnr$genes <- as.character(mod_dnr$genes)
mod_dnr[mod_duplicate, ]$genes <- paste0("MBNL3*")
mod_dnr$genes <- as.character(mod_dnr$genes)
ggplot(mod_dnr, aes(x = jitter(as.numeric(Group)), y = Modification_value,
                    color = genes)) +
  geom_point(alpha = 0.7) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Depressed", "Non-depressed")) +
  scale_color_manual(values=c("red", "purple", "black", "yellow",
                              "salmon", "burlywood", "blue", "#984EA3",
                              "#69b3a2","skyblue")) +
  labs(x = "Groups", y = "Modification value",
       title = "Top 10 modification value positions's DNR profile")
#foreground and background ####################################################

fdr<-p.adjust(p_value_df$p_values, method = "fdr",
                     n = length(p_value_df$p_values))
p_value_df$fdr<-fdr

significant_gene_df <- p_value_df[p_value_df$fdr < 0.05, ]

foreground_genomemap<-genomemap[rownames(significant_gene_df),]
significant_gene_df$genes <- foreground_genomemap$ucsc_refgene_name

foreground_df <- data.frame(gene =significant_gene_df$genes,
                    p_value =significant_gene_df$fdr)
foreground_df <- subset(foreground_df, gene != "")
foregound <- unique(foreground_df$gene)
foregound <- sapply(strsplit(foregound, ";"), `[`, 1)
foregound <- unique(foregound)
backgroundd_df <- data.frame(gene =genomemap$ucsc_refgene_name,
                                    p_value =p_value_df$fdr)
backgroundd_df <- subset(backgroundd_df, gene != "")
background <- unique(backgroundd_df$gene)
background <- sapply(strsplit(background, ";"), `[`, 1)
background <- unique(background)
write.table(foreground_df, "foreground.tsv", sep = "\t", row.names = FALSE)
write.table(backgroundd_df, "background.tsv", sep = "\t", row.names = FALSE)
writeLines(background, "background.txt")
writeLines(foregound, "foregound.txt")
