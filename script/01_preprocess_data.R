suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidymodels))
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(recipeselectors))
suppressPackageStartupMessages(library(FSelectorRcpp))
suppressPackageStartupMessages(library(doParallel))


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
## set column roles
colRoles <- rep("predictor", ncol(covidData))
names(colRoles) <- colnames(covidData)
colRoles[metadataCols] <- "ID"
colRoles["class"] <- "outcome"

tibble::enframe(x = colRoles, name = "col", value = "role") %>% 
  readr::write_tsv(
    file = here::here("analysis/02_ML_models", "column_roles.txt")
  )

## prepare splits
splits <- rsample::initial_split(data = covidData, prop = 4/5, strata = viral_status)
splits <- rsample::populate(splits)

dataHoldout <- rsample::testing(splits)
dataModelling <- rsample::training(splits)

folds <- rsample::vfold_cv(
  data = dataModelling, v = 5, strata = viral_status
)

folds <- rsample::populate(folds)

## store splits information
splitsDf <- tibble::tibble(
  id = 1:nrow(covidData),
  initial_split = "analysis"
)

splitsDf$initial_split[splits$out_id] <- "assessment"

readr::write_tsv(
  x = splitsDf,
  file = here::here("analysis/02_ML_models", "initial_split.txt")
)


foldsDf <- purrr::transpose(.l = rsample::rsample2caret(folds)) %>% 
  purrr::map_dfr(
    .f = function(x){
      split <- rep("analysis", times = nrow(dataModelling))
      split[x$indexOut] <- "assessment"
      split
    }
  ) %>% 
  dplyr::mutate(id = 1:n()) %>% 
  dplyr::select(id, everything())

readr::write_tsv(
  x = foldsDf,
  file = here::here("analysis/02_ML_models", "cv_folds_5.txt")
)

