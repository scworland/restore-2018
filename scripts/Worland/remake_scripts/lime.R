# Skip to line #500
# LIME package functions. For some reason the multioutput model only works
# if the functions are loaded this way rather than via the package. This does 
# not need USGS review because I did not write the functions

# lime function ----
lime <- function(x, model, ...) {
  if (is.character(x) && is.image_file(x)) class(x) <- 'imagefile'
  UseMethod('lime', x)
}


model_permutations <- function(x, y, weights, labels, n_labels, n_features, feature_method) {
  if (all(weights[-1] == 0)) {
    stop('All permutations have no similarity to the original observation. Try setting bin_continuous to TRUE and/or increase kernel_size', call. = FALSE)
  }
  if (!is.null(n_labels)) {
    labels <- names(y)[order(y[1,], decreasing = TRUE)[seq_len(n_labels)]]
  }
  x <- x[, colSums(is.na(x)) == 0 & apply(x, 2, var) != 0, drop = FALSE]
  res <- lapply(labels, function(label) {
    
    features <- select_features(feature_method, x, y[[label]], weights, n_features)
    
    # glmnet does not allow n_features=1
    if (length(features) == 1) {
      x_fit = cbind("(Intercept)" = rep(1, nrow(x)), x[, features, drop = FALSE])
      fit <- glm.fit(x = x_fit, y = y[[label]],  weights = weights, family = gaussian())
      r2 <- fit$deviance / fit$null.deviance
      coefs <- coef(fit)
      intercept <- coefs[1]
      coefs <- coefs[-1]
      model_pred <- fit$fitted.values[1]
    } else {
      shuffle_order <- sample(length(y[[label]])) # glm is sensitive to the order of the examples
      fit <- glmnet(x[shuffle_order, features], y[[label]][shuffle_order], weights = weights[shuffle_order], alpha = 0, lambda = 0.001)
      r2 <- fit$dev.ratio
      coefs <- coef(fit)
      intercept <- coefs[1, 1]
      coefs <- coefs[-1, 1]
      model_pred <- predict(fit, x[1, features, drop = FALSE])[1]
    }
    
    data.frame(
      label = label,
      feature = names(coefs),
      feature_weight = unname(coefs),
      model_r2 = r2,
      model_intercept = intercept,
      model_prediction = model_pred,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, res)
}


feature_selection_method <- function() c("auto", "none", "forward_selection", "highest_weights", "lasso_path", "tree")


select_features <- function(method, x, y, weights, n_features) {
  if (n_features >= ncol(x)) {
    return(seq_len(ncol(x)))
  }
  method <- match.arg(method, feature_selection_method())
  switch(
    method,
    auto = if (n_features <= 6) {
      select_features("forward_selection", x, y, weights, n_features)
    } else {
      select_features("highest_weights", x, y, weights, n_features)
    },
    none = seq_len(nrow(x)),
    forward_selection = select_f_fs(x, y, weights, n_features),
    highest_weights = select_f_hw(x, y, weights, n_features),
    lasso_path = select_f_lp(x, y, weights, n_features),
    tree = select_tree(x, y, weights, n_features),
    stop("Method not implemented", call. = FALSE)
  )
}

#' @importFrom glmnet cv.glmnet
select_f_fs <- function(x, y, weights, n_features) {
  features <- c()
  for (i in seq_len(n_features)) {
    max <- -100000
    best <- 0
    for (j in seq_len(ncol(x))) {
      if (j %in% features) next
      if (length(features) == 0) {
        x_fit = cbind("(Intercept)" = rep(1, nrow(x)), x[, j, drop = FALSE])
        fit <- glm.fit(x = x_fit, y = y,  weights = weights, family = gaussian())
        r2 <- fit$deviance / fit$null.deviance
      } else {
        fit <- glmnet(x[, c(features, j), drop = FALSE], y, weights = weights, alpha = 0, lambda = 0)
        r2 <- fit$dev.ratio
      }
      if (is.finite(r2) && r2 > max) {
        max <- r2
        best <- j
      }
    }
    features <- c(features, best)
  }
  features
}

#' @importFrom glmnet glmnet coef.cv.glmnet
#' @importFrom stats coef
#' @importFrom utils head
select_f_hw <- function(x, y, weights, n_features) {
  shuffle_order <- sample(length(y)) # glm is sensitive to the order of the examples
  fit_model <- glmnet(x[shuffle_order,], y[shuffle_order], weights = weights[shuffle_order], alpha = 0, lambda = 0)
  features <- coef(fit_model)[-1, 1]
  features_order <- order(abs(features), decreasing = TRUE)
  head(features_order, n_features)
}

# Tree model for feature selection
# Based on the latest XGBoost version.
# May require the Drat package because of a bug in old version of xgb.model.dt.tree
# @param x the data as a sparse matrix
# @param y the labels
# @param weights distance of the sample with the original datum
# @param n_features number of features to take
#' @importFrom utils packageVersion
select_tree <- function(x, y, weights, n_features) {
  xgb_version <- packageVersion("xgboost")
  if (xgb_version < "0.6.4.6") stop("You need to install latest xgboost (version >= \"0.6.4.6\") from Xgboost Drat repository to use tree mode for feature selection.\nMore info on http://xgboost.readthedocs.io/en/latest/R-package/xgboostPresentation.html")
  number_trees <- max(trunc(log2(n_features)), 2)
  if (log2(n_features) != number_trees) message("In \"tree\" mode, number of features should be a power of 2 and at least of 4 (= deepness of the binary tree of 2), setting was set to [", n_features, "], it has been replaced by [", 2^number_trees, "].")
  mat <- xgboost::xgb.DMatrix(x, label = y, weight = weights)
  bst.bow <- xgboost::xgb.train(params = list(max_depth = number_trees, eta = 1, silent = 1, objective = "binary:logistic"), data = mat, nrounds = 1, lambda = 0)
  dt <- xgboost::xgb.model.dt.tree(model = bst.bow)
  selected_words <- head(dt[["Feature"]], n_features)
  which(colnames(mat) %in% selected_words)
}

#' @importFrom glmnet glmnet coef.glmnet
#' @importFrom stats coef
select_f_lp <- function(x, y, weights, n_features) {
  shuffle_order <- sample(length(y)) # glm is sensitive to the order of the examples
  fit <- glmnet(x[shuffle_order,], y[shuffle_order], weights = weights[shuffle_order], alpha = 1, nlambda = 300)
  # In case that no model with correct n_feature size was found
  if (all(fit$df != n_features)) {
    stop(sprintf("No model with %i features found with lasso_path. Try a different method.", n_features))
  }
  has_value <- apply(coef(fit)[-1, ], 2, function(x) x != 0)
  f_count <- apply(has_value, 2, sum)
  row <- which(f_count >= n_features)[1]
  features <- which(has_value[, row])
  if (length(features) > n_features) {
    features <- features[sample(seq_along(features), n_features)]
  }
  features
}
exp_kernel <- function(width) {
  function(x) sqrt(exp(-(x^2) / (width^2)))
}

# lime data frame function ----
lime.data.frame <- function(x, model, bin_continuous = TRUE, n_bins = 4, quantile_bins = TRUE, use_density = TRUE, ...) {
  explainer <- c(as.list(environment()), list(...))
  explainer$x <- NULL
  explainer$feature_type <- setNames(sapply(x, function(f) {
    if (is.integer(f)) {
      if (length(unique(f)) == 1) 'constant' else 'integer'
    } else if (is.numeric(f)) {
      if (length(unique(f)) == 1) 'constant' else 'numeric'
    } else if (is.character(f)) {
      'character'
    } else if (is.factor(f)) {
      'factor'
    } else if (is.logical(f)) {
      'logical'
    } else if (inherits(f, 'Date') || inherits(f, 'POSIXt')) {
      'date_time'
    } else {
      stop('Unknown feature type', call. = FALSE)
    }
  }), names(x))
  if (any(explainer$feature_type == 'constant')) {
    warning('Data contains numeric columns with zero variance', call. = FALSE)
  }
  explainer$bin_cuts <- setNames(lapply(seq_along(x), function(i) {
    if (explainer$feature_type[i] %in% c('numeric', 'integer')) {
      if (quantile_bins) {
        bins <- quantile(x[[i]], seq(0, 1, length.out = n_bins + 1), na.rm = TRUE)
        bins[!duplicated(bins)]
      } else {
        d_range <- range(x[[i]], na.rm = TRUE)
        seq(d_range[1], d_range[2], length.out = n_bins + 1)
      }
    }
  }), names(x))
  explainer$feature_distribution <- setNames(lapply(seq_along(x), function(i) {
    switch(
      explainer$feature_type[i],
      integer = ,
      numeric = if (bin_continuous) {
        table(cut(x[[i]], unique(explainer$bin_cuts[[i]]), labels = FALSE, include.lowest = TRUE))/nrow(x)
      } else if (use_density) {
        density(x[[i]])
      } else {
        c(mean = mean(x[[i]], na.rm = TRUE), sd = sd(x[[i]], na.rm = TRUE))
      },
      character = ,
      logical = ,
      factor = table(x[[i]])/nrow(x),
      NA
    )
  }), names(x))
  structure(explainer, class = c('data_frame_explainer', 'explainer', 'list'))
}

# explain function ----
explain <- function(x, explainer, labels, n_labels = NULL, n_features,
                    n_permutations = 5000, feature_select = 'auto', n, ...) {
  if (is.character(x) && is.image_file(x)) class(x) <- 'imagefile'
  UseMethod('explain', x)
}
model_type.explainer <- function(x) {
  model_type(x$model)
}
output_type <- function(x) {
  switch(
    model_type(x),
    classification = 'prob',
    regression = 'raw',
    stop(model_type(x), ' models are not supported yet', call. = FALSE)
  )
}

# explain data frame function ----
explain.data.frame <- function(x, explainer, labels = NULL, n_labels = NULL,
                               n_features, n_permutations = 5000,
                               feature_select = 'auto', dist_fun = 'gower',
                               kernel_width = NULL, n, ...) {
  assert_that(is.data_frame_explainer(explainer))
  m_type <- model_type(explainer)
  o_type <- output_type(explainer)
  if (m_type == 'regression') {
    if (!is.null(labels) || !is.null(n_labels)) {
      warning('"labels" and "n_labels" arguments are ignored when explaining regression models')
    }
    n_labels <- 1
    labels <- NULL
  }
  assert_that(is.null(labels) + is.null(n_labels) == 1, msg = "You need to choose between labels and n_labels parameters.")
  assert_that(is.count(n_features))
  assert_that(is.count(n_permutations))
  
  if (is.null(kernel_width)) {
    kernel_width <- sqrt(ncol(x)) * 0.75
  }
  kernel <- exp_kernel(kernel_width)
  
  case_perm <- permute_cases(x, n_permutations, explainer$feature_distribution,
                             explainer$bin_continuous, explainer$bin_cuts,
                             explainer$use_density)
  #case_res <- predict_model(explainer$model, case_perm, type = o_type)
  case_res <- data.frame(predict(object=explainer$model, x=as.matrix(case_perm))[n]) %>% setNames("Response")
  case_res <- set_labels(case_res, explainer$model)
  case_ind <- split(seq_len(nrow(case_perm)), rep(seq_len(nrow(x)), each = n_permutations))
  res <- lapply(seq_along(case_ind), function(ind) {
    i <- case_ind[[ind]]
    if (dist_fun == 'gower') {
      sim <- 1 - gower_dist(case_perm[i[1], , drop = FALSE], case_perm[i, , drop = FALSE])
    }
    perms <- numerify(case_perm[i, ], explainer$feature_type, explainer$bin_continuous, explainer$bin_cuts)
    if (dist_fun != 'gower') {
      sim <- kernel(c(0, dist(feature_scale(perms, explainer$feature_distribution, explainer$feature_type, explainer$bin_continuous),
                              method = dist_fun)[seq_len(n_permutations-1)]))
    }
    res <- model_permutations(as.matrix(perms), case_res[i, , drop = FALSE], sim, labels, n_labels, n_features, feature_select)
    res$feature_value <- unlist(case_perm[i[1], res$feature])
    res$feature_desc <- describe_feature(res$feature, case_perm[i[1], ], explainer$feature_type, explainer$bin_continuous, explainer$bin_cuts)
    guess <- which.max(abs(case_res[i[1], ]))
    res$case <- rownames(x)[ind]
    res$label_prob <- unname(as.matrix(case_res[i[1], ]))[match(res$label, colnames(case_res))]
    res$data <- list(as.list(case_perm[i[1], ]))
    res$prediction <- list(as.list(case_res[i[1], ]))
    res$model_type <- m_type
    res
  })
  res <- do.call(rbind, res)
  res <- res[, c('model_type', 'case', 'label', 'label_prob', 'model_r2', 'model_intercept', 'model_prediction', 'feature', 'feature_value', 'feature_weight', 'feature_desc', 'data', 'prediction')]
  if (m_type == 'regression') {
    res$label <- NULL
    res$label_prob <- NULL
    res$prediction <- unlist(res$prediction)
  }
  res
}
is.data_frame_explainer <- function(x) inherits(x, 'data_frame_explainer')
#' @importFrom stats setNames
numerify <- function(x, type, bin_continuous, bin_cuts) {
  setNames(as.data.frame(lapply(seq_along(x), function(i) {
    if (type[i] %in% c('character', 'factor', 'logical')) {
      as.numeric(x[[i]] == x[[i]][1])
    } else if (type[i] == 'date_time' || type[i] == 'constant') {
      rep(0, nrow(x))
    } else {
      if (bin_continuous) {
        cuts <- bin_cuts[[i]]
        cuts[1] <- -Inf
        cuts[length(cuts)] <- Inf
        xi <- cut(x[[i]], unique(cuts), include.lowest = T)
        as.numeric(xi == xi[1])
      } else {
        x[[i]]
      }
    }
  }), stringsAsFactors = FALSE), names(x))
}
#' @importFrom stats setNames
feature_scale <- function(x, distribution, type, bin_continuous) {
  setNames(as.data.frame(lapply(seq_along(x), function(i) {
    if (type[i] == 'numeric' && !bin_continuous) {
      scale(x[, i], distribution[[i]]['mean'], distribution[[i]]['sd'])
    } else {
      x[, i]
    }
  }), stringsAsFactors = FALSE), names(x))
}
describe_feature <- function(feature, case, type, bin_continuous, bin_cuts) {
  sapply(feature, function(f) {
    if (type[[f]] == 'logical') {
      paste0(f, ' is ', tolower(as.character(case[[f]])))
    } else if (type[[f]] %in% c('character', 'factor')) {
      paste0(f, ' = ', as.character(case[[f]]))
    } else if (bin_continuous) {
      cuts <- bin_cuts[[f]]
      cuts[1] <- -Inf
      cuts[length(cuts)] <- Inf
      bin <- cut(case[[f]], unique(cuts), labels = FALSE, include.lowest = TRUE)
      cuts <- trimws(format(cuts, digits = 3))
      if (bin == 1) {
        paste0(f, ' <= ', cuts[bin + 1])
      } else if (bin == length(cuts) - 1) {
        paste0(cuts[bin], ' < ', f)
      } else {
        paste0(cuts[bin], ' < ', f, ' <= ', cuts[bin + 1])
      }
    } else {
      f
    }
  })
}

# permute cases ----
permute_cases <- function(cases, n_permutations, ...) {
  UseMethod('permute_cases')
}
#' @importFrom stats rnorm runif
permute_cases.data.frame <- function(cases, n_permutations, feature_distribution, bin_continuous, bin_cuts, use_density) {
  nrows <- nrow(cases) * n_permutations
  perm <- as.data.frame(lapply(seq_along(cases), function(i) {
    perms <- if (is.numeric(cases[[i]]) && is.na(feature_distribution[[i]])) {
      rep(cases[[i]], each = n_permutations)
    } else if (is.numeric(cases[[i]]) && bin_continuous) {
      bin <- sample(seq_along(feature_distribution[[i]]), nrows, TRUE, as.numeric(feature_distribution[[i]]))
      diff(bin_cuts[[i]])[bin] * runif(nrows) + bin_cuts[[i]][bin]
    } else if (is.numeric(cases[[i]]) && use_density) {
      width <- diff(feature_distribution[[i]]$x[c(1, 2)]) / 2
      sample(feature_distribution[[i]]$x, nrows, TRUE, feature_distribution[[i]]$y) + runif(nrows, -width, width)
    } else if (is.numeric(cases[[i]])) {
      rnorm(nrows) * feature_distribution[[i]]['sd'] + feature_distribution[[i]]['mean']
    } else if (is.character(cases[[i]])) {
      sample(names(feature_distribution[[i]]), nrows, TRUE, as.numeric(feature_distribution[[i]]))
    } else if (is.logical(cases[[i]])) {
      sample(as.logical(names(feature_distribution[[i]])), nrows, TRUE, as.numeric(feature_distribution[[i]]))
    } else if (is.factor(cases[[i]])) {
      x <- sample(names(feature_distribution[[i]]), nrows, TRUE, as.numeric(feature_distribution[[i]]))
      factor(x, levels = names(feature_distribution[[i]]))
    } else if (inherits(cases[[i]], 'Date') || inherits(cases[[i]], 'POSIXt')) {
      rep(cases[[i]], each = n_permutations)
    }
    if (is.integer(cases[[i]])) {
      as.integer(round(perms))
    } else {
      perms
    }
  }), stringsAsFactors = FALSE)
  names(perm) <- names(cases)
  perm[seq.int(1, by = n_permutations, length.out = nrow(cases)), ] <- cases
  perm
}

# other functions ----
as_classifier <- function(x, labels = NULL) {
  class(x) <- c('lime_classifier', class(x))
  attr(x, 'lime_labels') <- labels
  x
}
#' @rdname as_classifier
#' @export
as_regressor <- function(x) {
  class(x) <- 'lime_regressor'
  x
}
set_labels <- function(res, model) {
  labels <- attr(model, 'lime_labels')
  if (model_type(model) == 'classification' && !is.null(labels)) {
    if (length(labels) != ncol(res)) {
      warning('Ignoring provided class labels as length differs from model output')
    } else {
      names(res) <- labels
    }
  }
  res
}

model_type <- function(x, ...) {
  UseMethod('model_type')
}

model_type.default <- function(x, ...) {
  stop('The class of model must have a model_type method. See ?model_type to get an overview of models supported out of the box', call. = FALSE)
}

model_type.lime_classifier <- function(x, ...) 'classification'

model_type.lime_regressor <- function(x, ...) 'regression'

model_type.train <- function(x, ...) {
  tolower(x$modelType)
}

model_type.WrappedModel <- function(x, ...) {
  switch(
    x$learner$type,
    classif = 'classification',
    regr = 'regression',
    surv = 'survival',
    cluster = 'clustering',
    multilabel = 'multilabel'
  )
}

model_type.xgb.Booster <- function(x, ...) {
  obj <- x$params$objective
  type <- strsplit(obj, ':')[[1]][1]
  switch(
    type,
    reg = 'regression',
    binary = 'classification',
    multi = 'classification',
    stop('Unsupported model type', call. = FALSE)
  )
}

model_type.lda <- function(x, ...) 'classification'

model_type.keras.engine.training.Model <- function(x, ...) {
  if (!requireNamespace('keras', quietly = TRUE)) {
    stop('The keras package is required for predicting keras models')
  }
  if (keras::get_layer(x, index = -1)$activation$func_name %in% c('relu','linear')) {
    'regression'
  } else {
    'classification'
  }
}

model_type.H2OModel <- function(x, ...) {
  h2o_model_class <- class(x)[[1]]
  if (h2o_model_class %in% c("H2OBinomialModel", "H2OMultinomialModel")) {
    return('classification')
  } else if (h2o_model_class == "H2ORegressionModel") {
    return('regression')
  } else {
    stop('This h2o model is not currently supported.')
  }
}

predict_model <- function(x, newdata, type, ...) {
  UseMethod('predict_model')
}
#' @export
predict_model.default <- function(x, newdata, type, ...) {
  p <- predict(x, newdata = newdata, type = type, ...)
  if (type == 'raw') p <- data.frame(Response = p, stringsAsFactors = FALSE)
  as.data.frame(p)
}
#' @export
predict_model.WrappedModel <- function(x, newdata, type, ...) {
  if (!requireNamespace('mlr', quietly = TRUE)) {
    stop('mlr must be available when working with WrappedModel models')
  }
  p <- predict(x, newdata = newdata, ...)
  type2 <- switch(
    type,
    raw = data.frame(Response = mlr::getPredictionResponse(p), stringsAsFactors = FALSE),
    prob = mlr::getPredictionProbabilities(p, p$task.desc$class.levels),
    stop('Type must be either "raw" or "prob"', call. = FALSE)
  )
}

# S Worland function

# set keras models as 'regression' for lime
# model_type.keras.engine.training.Model <- function(x, ...) {
#   "regression"
# }

# predict_model.keras.engine.training.Model <- function(x, newdata, type='raw', n, ...) {
#   if (!requireNamespace('keras', quietly = TRUE)) {
#     stop('The keras package is required for predicting keras models')
#   }
#   
#   res <- predict(object=x, x=as.matrix(newdata)) %>%
#     data.frame()
#   
#   if (type == 'raw') {
#     return(data.frame(Response = res[, n]))
#   } else {
#     if (ncol(res) == 1) {
#       res <- cbind(1 - res, res)
#     }
#     colnames(res) <- as.character(seq_len(ncol(res)))
#     as.data.frame(res, check.names = FALSE)
#   }
# }


# apply lime to multioutput model
sw_lime <- function(final_fdc_model, gage_all, huc12_covariates) {
  
  d <- gage_all
  
  #model <- final_fdc_model
  
  f27 <- c("f0.02","f0.05","f0.1","f0.2","f0.5","f01","f02","f05",
           "f10","f20","f25","f30","f40","f50","f60","f70","f75",
           "f80","f90","f95","f98","f99","f99.5","f99.8","f99.9",
           "f99.95","f99.98")
  
  # covariates
  X <- select(d, ppt_mean, temp_mean, hay_pasture, nid_storage) %>%
    mutate_all(funs(as.numeric(as.factor(.)))) %>%
    mutate_all(funs(as.vector(scale(.)))) %>%
    select_if(~!any(is.na(.)))  %>%
    as.matrix()
  
  # X <- select(d,major:flood_storage, -basin_area) %>%
  #   mutate_all(funs(as.numeric(as.factor(.)))) %>%
  #   mutate_all(funs(as.vector(scale(.)))) %>%
  #   select_if(~!any(is.na(.)))  %>%
  #   as.matrix()
  
  # quantiles
  Y <- select(d,f27) %>%
    mutate_all(funs(.+0.001)) %>%
    mutate_all(funs(./area)) %>%
    mutate_all(funs(log10)) 
  
  # train NN model 
  input <- layer_input(shape=dim(X)[2],name="basinchars")

  base_model <- input  %>%
    layer_dense(units = 10,activation="relu") %>%
    layer_dropout(rate=0.5) %>%
    layer_dense(units = 10,activation="relu") %>%
    layer_dropout(rate=0.5)

  for(i in 1:dim(Y)[2]){
    y <- colnames(Y)[i]
    outstring <- paste0(
      sprintf("%s <- base_model %%>%%", y),
      sprintf(" layer_dense(units = 1, activation='linear', name='%s')",y)
    )
    eval(parse(text=outstring))
  }

  # not very sensitive, just make it small
  loss_weights <- c(rep(1,10),rep(1e-2,17))

  Ylist <- paste0("list(",paste(colnames(Y),sep="",collapse=","),")")
  model <- keras_model(input,eval(parse(text=Ylist))) %>%
    compile(optimizer = optimizer_rmsprop(lr = 0.001),
            loss_weights = loss_weights,
            loss="mse",
            metrics="mae")

  model_fit <- model %>%
    fit(x=X,
        y=Y,
        epochs=150,
        batch_size = 450,
        validation_split=0.1,
        verbose=0)
  
  # extract biases in last layer
  model_weights <- keras::get_weights(model)
  
  bias <- NULL
  for(i in 1:length(model_weights)){
    weight <- model_weights[[i]]
    
    if(dim(weight)==1){
      bias <- rbind(bias,weight)
    }
  }
  
  # make tibble for LIME
  X_explainer <- as.tibble(X)
  
  # prepare covariates for model predictions
  all_X <- huc12_covariates %>% 
    distinct(comid,decade,.keep_all=T) %>% # drop duplicate comids
    select(major:flood_storage, -basin_area) %>%
    filter(decade==2000) %>%
    select(ppt_mean, temp_mean, hay_pasture, nid_storage) %>%
    mutate_all(funs(as.numeric(as.factor(.)))) %>%
    mutate_all(funs(as.vector(scale(.)))) %>%
    select_if(~!any(is.na(.))) %>%
    as.tibble()
  
  explanation_all <- NULL
  for (j in 1:length(f27)) {
  #for (j in c(1,9,14,19,27)) {
    
    print(paste0("Calculating LIME for quantile number ",j," out of ", length(f27)))
    
    # predict_model.keras.engine.training.Model <- function(x, newdata, type='raw', n=j, ...) {
    #   if (!requireNamespace('keras', quietly = TRUE)) {
    #     stop('The keras package is required for predicting keras models')
    #   }
    #   
    #   res <- predict(object=x, x=as.matrix(newdata)) %>%
    #     data.frame()
    #   
    #   if (type == 'raw') {
    #     return(data.frame(Response = res[, n]))
    #   } else {
    #     if (ncol(res) == 1) {
    #       res <- cbind(1 - res, res)
    #     }
    #     colnames(res) <- as.character(seq_len(ncol(res)))
    #     as.data.frame(res, check.names = FALSE)
    #   }
    # }

    # run lime() on training set
    explainer <- lime::lime(x = X_explainer, 
                            model = model, 
                            bin_continuous = FALSE)
    
    # run explain() on huc12 pour point data
    explanation <- lime::explain(x = all_X, 
                                 explainer = explainer, 
                                 n_features = 5,
                                 kernel_width = 0.2,
                                 n_permutations = 500,
                                 feature_select = "highest_weights",
                                 n = j) %>%
      mutate(quantile = f27[j])
    
    explanation_all <- rbind(explanation_all, explanation)
  }
  
  
  result <- list(explanation_all=explanation_all,
                  bias=bias)
  
  return(result)
  
}

# gage_all <- remake::fetch('gage_all')

