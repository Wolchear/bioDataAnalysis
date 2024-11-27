data <- readRDS("data2.rds")
genomemap <- readRDS("genomemap2.rds")
samplekey <- readRDS("samplekey2.rds")

library(matrixTests)
library(dplyr)
library(ggplot2)
library(remotes)
library(ggrepel)
library(reshape)
library(methylCIPHER)
#install_github("karoliskoncevicius/matrixTests@dev_lm")
############################## Hipotezes ######################################
mod_1  <- model.matrix(~ age + smoking, data = samplekey)
mod_0 <-  model.matrix(~ smoking, data = samplekey)
res<-row_lm_f(data, m = mod_1, null = mod_0)
p_value <- data.frame(res$pvalue)
rownames(p_value)<- rownames(res)
p_value<-dplyr::rename(p_value, p_value = "res.pvalue")

############################ Korekcija ########################################
p_value$fdr<-p.adjust(p_value$p_value, method = "fdr",
                   n = length(p_value$p_value))

p_value$bonferroni<-p.adjust(p_value$p_value, method = "bonferroni",
                          n = length(p_value$p_value))

################ Histogramos ###############################################
ggplot(p_value, aes(x = p_value)) + 
  geom_histogram(binwidth = 0.025, color = "grey", fill = "salmon") +
  geom_vline(xintercept = 0.05, color = "red", size = 0.7,
             linetype = "longdash") +
  labs(x = "P value", y = "Count",
       title = "P-value distribution histogram, no correction") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.05))

ggplot(p_value, aes(x = fdr)) + 
  geom_histogram(binwidth = 0.025, color = "grey", fill = "salmon") +
  geom_vline(xintercept = 0.05, color = "red", size = 0.7,
             linetype = "longdash") +
  labs(x = "P value", y = "Count",
       title = "P-value distribution histogram, fdr correcction") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.05))

ggplot(p_value, aes(x = bonferroni)) + 
  geom_histogram(binwidth = 0.025, color = "grey", fill = "salmon") +
  geom_vline(xintercept = 0.05, color = "red", size = 0.7,
             linetype = "longdash") +
  labs(x = "P value", y = "Count",
       title = "P-value distribution histogram, bonferoni correction") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.05))

sum(p_value$p_value < 0.05)
sum(p_value$fdr < 0.05)
sum(p_value$bonferroni < 0.05)
#################### Voolcano_plot ############################################
p_value$effect_size <- res$beta.2
p_value <- p_value %>%
  mutate(significance = case_when(effect_size > 0 & p_value <= 0.05 ~ "up",
                                  effect_size <= 0 & p_value <= 0.05 ~ "down",
                                  TRUE ~ "ns"))

sorted_df <- p_value[order(p_value$p_value), ]
sorted_index <- rownames(head(sorted_df, 3))
p_value$delabel <- ifelse(rownames(p_value) %in% sorted_index,
                             rownames(p_value), NA)

ggplot(p_value, aes(x = effect_size, y = -log10(p_value),
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
       title = "Volcano plot") +
  geom_text_repel(max.overlaps = Inf) +
  coord_cartesian(xlim = c(-0.3, 0.3)) 
############### Manhaten plot #################################################
p_value$chr <- genomemap$chr
chrTable <- table(p_value$chr[p_value$p_value < 0.05])
chrTable_all <- table(p_value$chr)
sumPos <- sum(table(genomemap$chr))
chrPercent <- (chrTable / chrTable_all) * 100
chrPercent_df <- as.data.frame(chrPercent)
colnames(chrPercent_df) <- c("Chromosome", "Percentage")

genomemap_chr_1 <- genomemap[genomemap$chr == "chr18",]
p_value_1 <- p_value[rownames(genomemap_chr_1),]
p_value_1$Pos <- genomemap_chr_1$pos
p_value_1$gene <- genomemap_chr_1$ucsc_refgene_name
p_value_1$gene <- ifelse(p_value_1$gene == "",
                            "No gene", p_value_1$gene)
p_value_1$gene <- sapply(strsplit(p_value_1$gene, ";"), `[`, 1)
ggplot(p_value_1, aes(x = Pos, y = -log10(p_value))) +
  geom_point(aes(color = -log10(p_value) > 10), show.legend = FALSE) +
  scale_color_manual(values = c("FALSE" = "grey",
                                "TRUE" = "skyblue")) +
  labs(x = "Position", y = "P-value", title = "Manhattan Plot") +
  theme_minimal() +
  scale_x_continuous(labels = scales::number_format()) +
  geom_text(data = subset(p_value_1, -log10(p_value) > 10),
            aes(label = gene),
            hjust = -0.4, vjust = -1, size = 2)

########################### DNA profiles #######################################
p_value_top_10 <- head(p_value[order(p_value$p_value), ],10)

samplekey_b_li <- subset(samplekey, celltype == 'bcell')

non_smoker_index<- rownames(samplekey_b_li[
        which(samplekey_b_li$smoking =="non-smoker"), ])
smoker_index<- rownames(samplekey_b_li[
  which(samplekey_b_li$smoking =="smoker"), ])

smoker_data <- data[rownames(p_value_top_10), smoker_index]
non_smoker_data <- data[rownames(p_value_top_10), non_smoker_index]

smoker_data <- melt(smoker_data)
smoker_data <- mutate(smoker_data, X2 = "Smoker")

non_smoker_data <- melt(non_smoker_data)
non_smoker_data <- mutate(non_smoker_data, X2 = "Non_smoker")
mod_dnr <- rbind(smoker_data, non_smoker_data)
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

ggplot(mod_dnr, aes(x = jitter(as.numeric(Group)), y = Modification_value,
                    color = genes)) +
  geom_point(alpha = 0.7) +
  scale_x_continuous(breaks = c(1, 2), labels = c("Smoker", "Non-smoker")) +
  scale_color_manual(values=c("red", "purple", "black", "yellow",
                              "salmon", "burlywood", "blue", "#984EA3",
                              "#69b3a2","skyblue")) +
  labs(x = "Groups", y = "Modification value",
       title = "Top 10 modification value positions's DNR profile")

############ horvarth ##########################################################
horvarth <- read.csv("horvarth.csv", sep = ",")
# data_horvarth <- horvarth[horvarth$CpGmarker %in% rownames(data) ,]
intercept <-horvarth[1,]$CoefficientTraining

common_cpg_sites <- intersect(rownames(data), horvarth$CpGmarker)
horvath_cpg <- horvarth[horvarth$CpGmarker %in% common_cpg_sites,]
methylation_data <- data[common_cpg_sites,]
weights <- horvath_cpg$CoefficientTraining
names(weights) <- horvath_cpg$CpGmarker
m_age <- t(methylation_data) %*% matrix(data=weights[rownames(methylation_data)])
m_age <- m_age + as.numeric(intercept)
anti.trafo <- function(x, adult.age = 20) { 
  ifelse( x < 0, (1 + adult.age) * exp(x) - 1, (1 + adult.age) * x + adult.age) 
}
result <- anti.trafo(m_age)

donors <- unique(samplekey$donor)
avg_age <- data.frame(donor = character(), avg.age = numeric())

for (donor in unique(donors)) {
  donor_indices <- grep(donor, rownames(result))
  avg_value <- mean(result[donor_indices])
  avg_age <- rbind(avg_age, data.frame(donor = donor, avg.age = avg_value))
}
unique_donors <- unique(samplekey[c("donor", "age")])
combined_data <- merge(unique_donors, avg_age, by = "donor")

ggplot(combined_data, aes(x = age, y = avg.age)) +
  geom_point(color = "blue") +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  labs(x = "Real_age", y = "Predicted_age",
       title = "Scatter plot for real and predicted age") + 
  theme_minimal()   # 
combined_data$age_diff<- combined_data$age - combined_data$avg.age
combined_data_mean<-mean(combined_data$age_diff)
combined_data_max<-max(combined_data$age_diff)
combined_data_min<-min(combined_data$age_diff)
######################### Other epigenetic clock ###############################

methylCIPHER_Horvath1 <- calcHorvath1(t(data), samplekey, imputation = F)
methylCIPHER_Horvath2 <- calcHorvath2(t(data), samplekey, imputation = F)
unique_donors_Horvath1 <- unique(methylCIPHER_Horvath1[c("donor", "age","Horvath1")])
avg_Horvath1 <- data.frame(donor = character(), Horvath1 = numeric())
for (donor in unique(donors)) {
  donor_indices <- grep(donor, rownames(methylCIPHER_Horvath1))
  avg_value <- mean(methylCIPHER_Horvath1$Horvath1[donor_indices])
  avg_Horvath1 <- rbind(avg_Horvath1, data.frame(donor = donor, Horvath1 = avg_value))
}
combined_data <- merge(combined_data, avg_Horvath1, by = "donor")
combined_data$age_diff_Horvath1<- combined_data$age - combined_data$Horvath1

unique_donors_Horvath1 <- unique(methylCIPHER_Horvath1[c("donor", "age","Horvath1")])
avg_Horvath2 <- data.frame(donor = character(), Horvath2 = numeric())
for (donor in unique(donors)) {
  donor_indices <- grep(donor, rownames(methylCIPHER_Horvath2))
  avg_value <- mean(methylCIPHER_Horvath2$Horvath2[donor_indices])
  avg_Horvath2 <- rbind(avg_Horvath2, data.frame(donor = donor, Horvath2 = avg_value))
}
combined_data <- merge(combined_data, avg_Horvath2, by = "donor")
combined_data$age_diff_Horvath2<- combined_data$age - combined_data$Horvath1
combined_data_mean_Horvath1<-mean(combined_data$age_diff_Horvath1)
combined_data_max_Horvath1<-max(combined_data$age_diff_Horvath1)
combined_data_min_Horvath1<-min(combined_data$age_diff_Horvath1)

combined_data_mean_Horvath2<-mean(combined_data$age_diff_Horvath2)
combined_data_max_Horvath2<-max(combined_data$age_diff_Horvath2)
combined_data_min_Horvath2<-min(combined_data$age_diff_Horvath2)


################# Papildoma ###################################################
genes <- c("ELOVL2", "FHL2", "EDARADD", "ASPA", "PDE4C", "PENK", "C1orf132", "TRIM59", "KLF14")

matching_rows <- genomemap[grepl(paste(genes, collapse="|"), genomemap$ucsc_refgene_name), ]
gene_predictors_data <- data[rownames(matching_rows),]
gene_predictors_Horvath1 <- calcHorvath1(t(gene_predictors_data), samplekey, imputation = F)

