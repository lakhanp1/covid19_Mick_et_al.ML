suppressPackageStartupMessages(library(tidymodels))
suppressPackageStartupMessages(library(DESeq2))



## Create the function
step_deseq2_vst <- function(
    recipe, 
    ..., 
    role = NA,
    trained = FALSE, 
    dispersionFn = NULL,
    geneIds = NULL,
    options = list(blind = TRUE, fitType = "parametric"),
    skip = FALSE,
    id = recipes::rand_id("DESeq2_vst")
) {
  
  # dispersionFn: this variable will store the dispersion function calculated
  # from the training set in prep() and will be applied to new dataset in bake()
  
  ## The variable selectors are not immediately evaluated by using
  ##  the `quos()` function in `rlang`. `ellipse_check()` captures 
  ##  the values and also checks to make sure that they are not empty.  
  terms <- ellipse_check(...) 
  
  recipes::add_step(
    recipe, 
    step_deseq2_vst_new(
      terms = terms, 
      trained = trained,
      role = role, 
      dispersionFn = dispersionFn,
      geneIds = geneIds,
      options = options,
      skip = skip,
      id = id
    )
  )
  
}

## Initialize a new object
step_deseq2_vst_new <- function(terms, role, trained, dispersionFn, geneIds, options, skip, id) {
  recipes::step(
    subclass = "deseq2_vst", 
    terms = terms,
    role = role,
    trained = trained,
    dispersionFn = dispersionFn,
    geneIds = geneIds,
    options = options,
    skip = skip,
    id = id
  )
}


## Create the prep method
#' @export
prep.step_deseq2_vst <- function(x, training, info = NULL, ...) {
  ## Remember that the prep() function does not apply the step to the data;
  ## it only estimates any required values 
  col_names <- recipes::recipes_eval_select(x$terms, training, info)
  recipes::check_type(training[, col_names])
  
  ## use col_names to ensure the same column order for the data
  ## compute the dispersion function from the training data
  dds_train <- DESeq2::DESeqDataSetFromMatrix(
    countData = as.matrix(t(training[, col_names])),
    colData = DataFrame(condition = rep("A", nrow(training))),
    design = ~1
  )
  dds_train <- DESeq2::estimateSizeFactors(dds_train)
  dds_train <- DESeq2::estimateDispersions(dds_train)
  
  dispersionFn <- DESeq2::dispersionFunction(dds_train)
  
  ## Use the constructor function to return the updated object. 
  ## Note that `trained` is now set to TRUE
  ## dispersionFn: has the dispersion function from training data
  step_deseq2_vst_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    dispersionFn = dispersionFn,
    geneIds = col_names,
    options = x$options,
    skip = x$skip,
    id = x$id
  )
}


## Create the bake method
bake.step_deseq2_vst <- function(object, new_data, ...) {
  # object: updated step function that has been through the corresponding prep() code
  # new_data: tibble of data to be processed. 
  # return: a tibble of the modified version of new_data
  
  ## use col_names to ensure the same column order for the data
  col_names <- object$geneIds
  
  recipes::check_new_data(col_names, object, new_data)
  
  dds_newData <- DESeq2::DESeqDataSetFromMatrix(
    countData = as.matrix(t(new_data[, col_names])),
    colData = DataFrame(condition = rep("A", nrow(new_data))),
    design = ~1
  )
  
  dds_newData <- DESeq2::estimateSizeFactors(dds_newData)
  
  ## apply dispersion function from obtained from the training data
  DESeq2::dispersionFunction(dds_newData) <- object$dispersionFn
  
  ## normalize the data
  vstCall <- expr(
    DESeq2::varianceStabilizingTransformation(
      object = dds_newData,
      !!!object$options
    )
  )
  
  vst <- eval(vstCall)
  
  normCounts <- tibble::as_tibble(t(assay(vst)))
  
  ## replace the count data with normalized data
  new_data[, col_names] <- normCounts[, col_names]
  
  new_data
}



