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

source(file = here::here("script", "00_recipe_vst_norm.R"))
source(file = here::here("script", "02_utils.R"))

set.seed(124)

file_geneCounts <- here::here("data", "swab_gene_counts.csv")
file_meta <- here::here("data", "metatable_with_viral_status.csv")
file_gene2Name <- here::here("data/annotation", "gene2name.txt")
file_split1 <- here::here("analysis/02_ML_models", "initial_split.txt")
file_folds <- here::here("analysis/02_ML_models", "cv_folds_5.txt")

cores <- 4
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

splitsDf <- suppressMessages(readr::read_tsv(file_split1))  #initial_split
foldsDf <- suppressMessages(readr::read_tsv(file_folds))  #vfold_cv
####################################################################

## extract cross-validation fold information
splits <- make_splits(
  x = split(as.integer(splitsDf$id), splitsDf$initial_split),
  data = covidData,
  class = "initial_split"
)

dataHoldout <- rsample::testing(splits)
dataModelling <- rsample::training(splits)

cvSplits <- dplyr::select(foldsDf, -id) %>% 
  purrr::map(
    .f = ~ split(as.integer(1:length(.x)), .x)
  ) %>% 
  purrr::map(~make_splits(x = .x, data = dataModelling, class = "vfold_split"))

folds <- manual_rset(splits = unname(cvSplits), ids = names(cvSplits))

####################################################################
## basic recipe for count normalization
rec_vst <- recipes::recipe(x = dataModelling) %>% 
  recipes::update_role(tidyselect::everything(), new_role = "predictor") %>% 
  recipes::update_role(class, new_role = "outcome") %>% 
  recipes::update_role(!!!syms(metadataCols), new_role = "ID") %>% 
  step_deseq2_vst(recipes::all_predictors())

# cooked_rcp <- check_if_recipe_cooks(rcp = rec_vst)

## feature selection recipe using information gain 
infogain_rec <- rec_vst %>% 
  recipeselectors::step_select_infgain(
    recipes::all_predictors(), outcome = "class",  
    type = "infogain", threads = !!cores, top_p = 50
  )

# cooked_rcp <- check_if_recipe_cooks(rcp = infogain_rec)

## feature extraction recipe using PCA
pca_rec <- rec_vst %>% 
  recipes::step_normalize(all_predictors()) %>% 
  recipes::step_pca(
    recipes::all_predictors(), threshold = 0.8
  )

# cooked_rcp <- check_if_recipe_cooks(rcp = pca_rec)

####################################################################

## logistic regression: Lasso
lasso_mod <- parsnip::logistic_reg(
  mixture = 1, penalty = tune()
) %>% 
  parsnip::set_engine("glmnet") %>% 
  parsnip::set_mode("classification")

## grid for optimum hyperparameter search
lasso_grid <- dials::grid_regular(
  dials::penalty(),
  levels = 10
)

####################################################################
#### PCA with Lasso
lasso_pca_workflow <- workflows::workflow() %>% 
  workflows::add_model(lasso_mod) %>% 
  workflows::add_recipe(pca_rec)

# doParallel::registerDoParallel(cores = 4)

lasso_pca_res <- lasso_pca_workflow %>% 
  tune::tune_grid(
    resamples = folds,
    grid = lasso_grid,
    control = control_grid(save_pred = TRUE),
    metrics = metric_set(roc_auc)
  )

# stopImplicitCluster()

# tune::collect_metrics(lasso_pca_res)
# tune::show_best(lasso_pca_res)

## select hyperparameters of the best performing model
lasso_pca_best <- fit_store_best_model(
  tune_res = lasso_pca_res, label = "PCA + Lasso",
  outPrefix = here::here("analysis/02_ML_models/lasso_pca"),
  truth_col = "class", pred_cols = ".pred_positive",
  wf = lasso_pca_workflow, evalSplits = splits
)

####################################################################
#### Information gain with Lasso
lasso_info_workflow <- workflows::workflow() %>% 
  workflows::add_model(lasso_mod) %>% 
  workflows::add_recipe(infogain_rec)

# registerDoParallel(cores=4)  
# cl <- makeCluster(4, type="FORK")  

lasso_info_res <- lasso_info_workflow %>% 
  tune::tune_grid(
    resamples = folds,
    grid = lasso_grid,
    control = control_grid(save_pred = TRUE),
    metrics = metric_set(roc_auc)
  )

# stopCluster(cl)

# tune::collect_metrics(lasso_info_res)
# tune::show_best(lasso_info_res)

## select hyperparameters of the best performing model
lasso_info_best <- fit_store_best_model(
  tune_res = lasso_info_res, label = "Information Gain + Lasso",
  outPrefix = here::here("analysis/02_ML_models/lasso_info"),
  truth_col = "class", pred_cols = ".pred_positive",
  wf = lasso_info_workflow, evalSplits = splits
)


