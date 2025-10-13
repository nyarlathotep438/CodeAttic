library(readr)
library(dplyr)
library(stringr)

# Data
Mel_cluster_top15 <- read_csv("./marker_gene/Mel/cluster_genes_res0.4.csv")
Mel_metastases_cluster_top15 <- read_csv("./marker_gene/Mel_metastases/cluster_genes_res0.3.csv")
MBM_cluster_top15 <- read_csv("./marker_gene/MBM/cluster_genes_res0.7.csv")

# Transform to ACT format####
# Mel(JobID: 20250709184111IYGWH5P0D5UMWB)
Mel_cluster_top15 %>%
  mutate(
    genes = str_split(top_genes, "/", simplify = TRUE)[, 1:15] %>%  # Split and take the first 15 columns
      apply(1, \(x) paste(na.omit(x), collapse = ","))  # Handling possible missing values
  ) %>%
  transmute(output = sprintf("cluster%d:%s", cluster, genes)) %>%
  pull(output) %>%
  cat(sep = "\n")

# Mel_metastases(JobID: 20250709183839QLVKLGAAIPJIIF)
Mel_metastases_cluster_top15 %>%
  mutate(
    genes = str_split(top_genes, "/", simplify = TRUE)[, 1:15] %>%  # Split and take the first 15 columns
      apply(1, \(x) paste(na.omit(x), collapse = ","))  # Handling possible missing values
  ) %>%
  transmute(output = sprintf("cluster%d:%s", cluster, genes)) %>%
  pull(output) %>%
  cat(sep = "\n")

# MBM(JobID: 202507091846071EGWCP7N5THTX8)
MBM_cluster_top15 %>%
  mutate(
    genes = str_split(top_genes, "/", simplify = TRUE)[, 1:15] %>%  # Split and take the first 15 columns
      apply(1, \(x) paste(na.omit(x), collapse = ","))  # Handling possible missing values
  ) %>%
  transmute(output = sprintf("cluster%d:%s", cluster, genes)) %>%
  pull(output) %>%
  cat(sep = "\n")

# ACT Result
Mel_ACT <- read.table("./ACT_result/Mel_ACT results.txt", header=TRUE, sep="\t", check.names=FALSE)
Mel_metastases_ACT <- read.table("./ACT_result/Mel_metastases_ACT results.txt",
                                 header=TRUE, sep="\t", check.names=FALSE)
MBM_ACT <- read.table("./ACT_result/MBM_ACT results.txt", header=TRUE, sep="\t", check.names=FALSE)
