suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(FactoMineR))



rm(list = ls())

set.seed(124)

file_geneCounts <- here::here("data", "swab_gene_counts.csv")
file_meta <- here::here("data", "metatable_with_viral_status.csv")
file_gene2Name <- here::here("data/annotation", "gene2name.txt")

####################################################################

sampleInfo <- suppressMessages(readr::read_csv(file = file_meta)) %>% 
  dplyr::rename(sampleId = CZB_ID) %>% 
  dplyr::mutate(
    class = dplyr::case_when(
      viral_status == "SC2" ~ "positive",
      viral_status == "no_virus" ~ "negative",
      viral_status == "other_virus" ~ "negative",
      TRUE ~ NA_character_
    ),
    viral_status = forcats::fct_relevel(
      .f = viral_status, "SC2", "no_virus", "other_virus"
    ),
    class = forcats::fct_relevel(
      .f = class, "positive", "negative"
    )
  )

metadataCols <- setdiff(colnames(sampleInfo), "class")

covidData <- suppressMessages(readr::read_csv(file = file_geneCounts)) %>% 
  dplyr::rename(geneId = "...1") %>% 
  tidyr::pivot_longer(cols = !geneId, names_to = "sampleId", values_to = "count") %>% 
  tidyr::pivot_wider(names_from = geneId, values_from = count) %>% 
  dplyr::left_join(y = sampleInfo, by = "sampleId") %>% 
  dplyr::select(class, !!!metadataCols, everything())

####################################################################

pt_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(face = "bold", size = 14, color = "black"),
    axis.title = element_text(face = "bold", size = 13, color = "black")
  )


## class wise count plot
pt_stats <- ggplot(
  data = sampleInfo,
  mapping = aes(x = viral_status, fill = gender)
) +
  geom_bar(
    stat = "count"
  ) +
  geom_text(
    mapping = aes(label = ..count..), stat = "count",
    position = position_stack(vjust = 0.5),
    color = "black", size = 6, fontface = "bold"
  ) +
  labs(title = "Class distribution of data") +
  scale_fill_manual(
    values = c("F" = "#F06C9B", "M" = "#61A0AF")
  ) +
  scale_x_discrete(
    labels = dplyr::count(sampleInfo, viral_status) %>% 
      dplyr::mutate(label = paste(viral_status, "\n(", n, ")", sep = "")) %>% 
      dplyr::pull(var = label, name = viral_status)
  ) +
  pt_theme

ggplot2::ggsave(
  filename = file.path(outDir, "class_stats.png"),
  plot = pt_stats, width = 10, height = 10
)

####################################################################

## normalize counts
countMat <- dplyr::select(covidData, -class, -!!metadataCols) %>% 
  as.matrix()

rownames(countMat) <- covidData$sampleId

normCountMat <- DESeq2::varianceStabilizingTransformation(t(countMat))

## PCA
## select top 2000 rows/genes with highest variation
rv <- matrixStats::rowVars(normCountMat)
keep <- order(rv, decreasing=TRUE)[seq_len(min(2000, length(rv)))]

## remove low count rows
# keep<- rowSums(normCountMat > 1) >= 2

normCountMatFiltered <- normCountMat[keep, ]

pcaData <- tibble::as_tibble(t(normCountMatFiltered), rownames = "sampleId") %>% 
  dplyr::left_join(y = sampleInfo, by = "sampleId") %>% 
  dplyr::select(class, !!!metadataCols, everything()) %>% 
  as.data.frame()

row.names(pcaData) <- pcaData$sampleId

res.pca <- FactoMineR::PCA(
  X = pcaData, graph = FALSE, scale.unit = TRUE,
  quali.sup = 1:ncol(sampleInfo), ncp = 10
)

eig.val <- factoextra::get_eigenvalue(res.pca)

## Graph of individuals
ind <- factoextra::get_pca_ind(res.pca)
var <- factoextra::get_pca_var(res.pca)

####################################################################
## scree plot: variance by PC
factoextra::fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 50))







