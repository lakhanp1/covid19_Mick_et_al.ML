


fit_store_best_model <- function(
    tune_res, outPrefix, label, truth_col, pred_cols, wf = NULL, evalSplits = NULL
){
  ## select best
  bestModel <- tune::show_best(tune_res, metric = "roc_auc", n = 1)
  readr::write_tsv(
    x = bestModel, file = paste(outPrefix, ".best_model.tsv", sep = "")
  )
  
  ## store AUC
  aucDf <- tune::collect_predictions(x = tune_res, parameters = bestModel) %>% 
    yardstick::roc_curve(!!truth_col, !!!pred_cols) %>%
    dplyr::mutate(model = label)
  
  readr::write_tsv(x = aucDf, file = paste(outPrefix, ".auc.tsv", sep = ""))
  
  finalFit <- NULL
  
  ## fit to initial_split
  if(!is.null(evalSplits) && !is.null(wf)){
    finalWf <- tune::finalize_workflow(x = wf, parameters = bestModel)
    finalFit <- tune::last_fit(object = finalWf, split = evalSplits)
    
    collect_metrics(finalFit) %>% 
      readr::write_tsv(file = paste(outPrefix, ".final_prediction.tsv", sep = ""))
  }
  
  return(
    list(bestModel = bestModel, auc = aucDf, finalFit = finalFit)
  )
}

####################################################################
## a utliity function to test if a recipe works as expected
check_if_recipe_cooks <- function(rcp){
  preped_rec <- recipes::prep(
    x = rcp,
    training = rsample::analysis(x = folds$splits[[1]])
  )
  
  test_baked <- recipes::bake(
    object = preped_rec,
    new_data = rsample::assessment(x = folds$splits[[1]])
  )
  
  return(test_baked)
}

####################################################################