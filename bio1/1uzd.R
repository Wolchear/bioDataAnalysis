data <- readRDS("data.rds")
genomemap <- readRDS("genomemap.rds")
samplekey <- readRDS("samplekey.rds")
study_participants<-read.csv("1-s2.0-S1872497323001151-mmc5.csv",sep = ";",
                              row.names = NULL)
library(WGCNA)
library(dplyr)
library(reshape)
library(ggplot2)
library(patchwork)
library(cluster)
library(gridExtra)
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


####### Kokybes kontrole
openSea <- genomemap[genomemap$relation_to_island == "OpenSea", ]
mean_openSea <- rowMeans(data[rownames(openSea), ])

island <- genomemap[genomemap$relation_to_island == "Island", ]
mean_island <- rowMeans(data[rownames(island), ])

n_shore <- genomemap[genomemap$relation_to_island == "N_Shore", ]
mean_n_shore <-  rowMeans(data[rownames(n_shore), ])

s_shore <- genomemap[genomemap$relation_to_island == "S_Shore", ]
mean_s_shore <- rowMeans(data[rownames(s_shore), ])

n_shelf <- genomemap[genomemap$relation_to_island == "N_Shelf", ]
mean_n_shelf <-  rowMeans(data[rownames(n_shelf), ])

s_shelf <- genomemap[genomemap$relation_to_island == "S_Shelf", ]
mean_s_shelf <-  rowMeans(data[rownames(s_shelf), ])

mean_by_group <- c(mean_island, mean_openSea, mean_n_shore,
                   mean_s_shore,mean_n_shelf, mean_s_shelf)

cpg_regions <- c(rep("Island", length(mean_island)),
                 rep("OpenSea", length(mean_openSea)),
                 rep("S_shore", length(mean_s_shore)),
                 rep("N_shore", length(mean_n_shore)),
                 rep("N_shelf", length(mean_n_shelf)),
                 rep("S_shelf", length(mean_s_shelf)))

cpg_regions<- data.frame(mean_by_group, cpg_regions)

ggplot(cpg_regions, aes(x = mean_by_group, fill = cpg_regions)) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot",
       x = "Mean Values",
       y = "Density") +
  theme_minimal()


brain_columns <- grep("brain", colnames(data), value = TRUE)
heart_columns <- grep("heart", colnames(data), value = TRUE)
dermis_columns <- grep("dermis", colnames(data), value = TRUE)
lung_columns <- grep("lung", colnames(data), value = TRUE)
epidermis_columns <- grep("epidermis", colnames(data), value = TRUE)
kidney_columns <- grep("kidney", colnames(data), value = TRUE)
blood_columns <- grep("blood", colnames(data), value = TRUE)
muscle_columns <- grep("muscle", colnames(data), value = TRUE)
liver_columns <- grep("liver", colnames(data), value = TRUE)

brain_columns_mean <- colMeans(data[,brain_columns ])
heart_columns_mean <- colMeans(data[,heart_columns ])
dermis_columns_mean <- colMeans(data[,dermis_columns ])
lung_columns_mean <- colMeans(data[,lung_columns ])
epidermis_columns_mean <- colMeans(data[,epidermis_columns ])
kidney_columns_mean <- colMeans(data[,kidney_columns ])
blood_columns_mean <- colMeans(data[,blood_columns ])
muscle_columns_mean <- colMeans(data[,muscle_columns ])
liver_columns_mean <- colMeans(data[,liver_columns ])

tissues <- c(rep("brain", length(brain_columns)),
                 rep("heart", length(heart_columns)),
                 rep("dermis", length(dermis_columns)),
                 rep("lung", length(lung_columns)),
                 rep("kidney", length(kidney_columns)),
                 rep("epidermis", length(epidermis_columns_mean)),
                 rep("blood", length(blood_columns)),
                 rep("muscle", length(muscle_columns)),
                 rep("liver", length(liver_columns))
                 )

tissue_by_group<-c(brain_columns_mean,heart_columns_mean,dermis_columns_mean,
                   lung_columns_mean,epidermis_columns_mean,kidney_columns_mean,
                   blood_columns_mean,muscle_columns_mean,liver_columns_mean)

tissues_df<- data.frame(tissue_by_group, tissues)

tissue_colors2<- labels2colors(tissues_df$tissues,
                        colorSeq =c("red","salmon", "green","blue","seagreen1",
                                  "honeydew2","black","white", "yellow"))



data_matrix <- matrix(tissue_by_group, nrow = length(tissues), ncol = 1)
hc <- hclust(dist(data_matrix), method = "complete")
plotDendroAndColors(hc,
                    cbind(
                      tissue_colors2
                    ),
                    hang=-1,
                    dendroLabels = FALSE,
                    main = "Cluster dendrogram")


#### Clustering
d <- dist(1 - cor(data))
clust <- hclust(d, method = "complete")
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

possible_depressed_colors <-labels2colors(samplekey$possible_depressed,
                                          colorSeq = c("grey","yellow"))

groups_colors<- labels2colors(samplekey$sex,
            colorSeq = c("pink","cyan"))

tissue_colors<- labels2colors(samplekey$tissue,
  colorSeq =c("red","black", "green","blue","seagreen1",
  "honeydew2","lightcoral","white", "yellow"))



plotDendroAndColors(clust,
                    cbind(
                            groups_colors,
                            tissue_colors,
                            possible_depressed_colors
                            ),
                    hang=-1,
                    dendroLabels = FALSE,
                    main = "Cluster dendrogram")

# data X chromosome

genomrows_X <- which(genomemap$chr == "chrX")
genomrows_open_sea <- which(genomemap$relation_to_island == "OpenSea")
genomrows_X_open_sea <- genomrows_open_sea[genomrows_open_sea %in% genomrows_X]
data_X<-data[rownames(genomemap)[genomrows_X_open_sea], ]

Coronary_patients <- study_participants[grepl("Coronary atherosclerosis",
                                        study_participants$disease), ]

samplekeys_coronary <- samplekey[samplekey$donor %in% Coronary_patients$donor, ]
samplekeys_heart_coronary2<-samplekeys_coronary[samplekeys_coronary$tissue == "heart",]
samplekeys_heart_coronary<-rownames(samplekeys_coronary[samplekeys_coronary$tissue == "heart",])
data_x_heart<-data_X[, colnames(data_X) %in% samplekeys_heart_coronary ]

data_x_heart <- as.data.frame(data_x_heart)
row_variability <- apply(data_x_heart, 1, sd)
top_10_variability <- data_x_heart[order(row_variability, decreasing = TRUE), ][1:10, ]


new_names <- c("CG1","CG2","CG3","CG4","CG5","CG6","CG7","CG8","CG9","CG10")
old_names <- rownames(top_10_variability)
rownames(top_10_variability) <- new_names
top_10_variability <- as.matrix(top_10_variability)

top_10_variability_melt <- melt(top_10_variability, as.is = TRUE)

heatmap<-ggplot(top_10_variability_melt, aes(X1, X2)) +
  labs(x = "Posions in X chromosome", y = "Patients") + 
  geom_tile(aes(fill = value), color = "black") +
  geom_text(aes(label = round(value, 2)),  size = 3) +
  scale_fill_gradient(low = "white", high = "red")

legend_text <- paste(old_names,
  "->", new_names, collapse = "\n")

legend <- ggplot() +
  geom_blank() +
  theme_void() +
  annotate("text", x = 0.5, y = 0.5, label = legend_text, size = 3)

new_heatmap<-heatmap + legend

print(new_heatmap)

#####
heatmap_list <- list()
for (tissue in unique(samplekey$tissue)) {
  tissue_index <- grep(paste0("_", tissue), colnames(data))
  tissue_data <- data[, tissue_index]
  tissue_correlation <- cor(tissue_data)
  colnames(tissue_correlation) <- paste0("p", seq_len(ncol(tissue_correlation)))
  rownames(tissue_correlation) <- paste0("p", seq_len(ncol(tissue_correlation)))
  tissue_correlation_melt <- melt(tissue_correlation, as.is = TRUE)
  
  heatmap<-ggplot(data = tissue_correlation_melt, aes(X2, X1, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "red") +
    geom_tile(aes(fill = value), color = "black") +
    labs(title = paste("Correlation Heatmap for", tissue)) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          plot.title = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  heatmap_list[[tissue]] <- heatmap
  ggsave(file.path("corelation_heatmaps_by_tissues",
                   paste0("Correlation_Heatmap_", tissue, ".png")), heatmap)
}
do.call(grid.arrange, c(heatmap_list, ncol = 3))

non_and_depressed_women<-c("2020_15_56_F_brain", "2020_109_52_F_brain",
                      "2020_157_41_F_brain", "2020_23_47_F_brain")
data_box_plot <- data[, non_and_depressed_women]
data_box_plot_melt <- melt(data_box_plot, as.is = TRUE)

ggplot(data_box_plot_melt, aes(x=X2, y=value, fill=X2)) +
  geom_boxplot()


