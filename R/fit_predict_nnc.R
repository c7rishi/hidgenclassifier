# # source('fit_nn_helper_fns.R')
#
#' Generate an object from the class "nn"
#'
#' See description of output of fit_nn() for details
#'
new_nn = function(map_df, model, model_raw, ind_val, tuning_results, preproc) {

  stopifnot(
    # isS4(X) | is.matrix(X) | is.vector(X) | is(X, "Matrix"),
    # is.vector(Y),
    is.list(map_df),
    typeof(model)=="closure",
    typeof(model_raw)=="raw",
    is.list(tuning_results),
    is.list(preproc)
  )
  # return(structure(list(X=X, Y=Y, map_df=map_df, model=model, ind_val=ind_val, tuning_results=tuning_results, preproc=preproc), class="nn"))
  return(
    structure(
      list(
        map_df=map_df,
        model=model,
        model_raw = model_raw,
        ind_val=ind_val,
        tuning_results=tuning_results,
        preproc=preproc
      ),
      class="nn"
    )
  )
}


#' Train a fully-connected multi-class neural network
#'
#'
#' @description
#' This function first splits the data into a training and validation set and tunes hyperparameters using Bayesian optimization (similar to the approach used in Jiao et al. 2020), then uses the best hyperparameters to train a model on the entire dataset.
#'
#' @inheritParams fit_smlc
#' @param val_split Fraction of data to be used as validation set for hyperparameters
#' @param trials Number of trials for hyperparameter tuning
#' @param epochs Number of training epochs
#' @param verbose_mbo Bayesian optimization verbosity mode (logical)
#' @param seed Random seed
#' @param ... Unused
#'
#' @return Object of class "nn", a named list of length 7 with the components of the neural network training process
#' \item{X}{Input matrix}
#' \item{Y}{Response vector}
#' \item{map_df}{Dataframe with columns "original" and "numeric". The "original" column contains the original class names in Y and the "numeric" column contains the numeric representation of the classes used during training}
#' \item{model}{Final Keras model trained on X and Y (see https://keras.rstudio.com/articles/about_keras_models.html for more details)}
#' \item{ind_val}{Vector of indices of X corresponding to validation set used to tune hyperparameters}
#' \item{tuning_results}{Named list with the results from the hyperparameter search (output of mbo() from mlrMBO). The list elements include "x", a named list with the best hyperparameters found, and "y", the validation accuracy corresponding to the best hyperparameters. See description of MBOSingleObjResult from mlrMBO for more details.}
#' \item{preproc}{Named list with the parameters of the min-max pre-processing transformation applied to X prior to training (output of preProcess() from caret)}
#'
#' @references
#' Jiao W, Atwal G, Polak P, Karlic R, Cuppen E, Danyi A, De Ridder J, van Herpen C, Lolkema MP, Steeghs N, Getz G. A deep learning system accurately classifies primary and metastatic cancers using passenger mutation patterns. Nature communications. 2020 Feb 5;11(1):1-2.
#'
#' @author Zoe Guan. Email: guanZ@mskcc.org
#' @export
fit_nnc = function(X, Y, val_split=1/3, trials=200, epochs=50, verbose_mbo=T, seed=1) {

  ### define mapping from original labels to numeric representation
  map_df = data.frame(name=unique(Y), num=as.numeric(factor(unique(Y))))


  ### apply one-hot encoding to Y
  num_classes = length(unique(Y))
  Y_numeric = as.numeric(factor(Y))
  Y_onehot = to_categorical(Y_numeric-1)

  ### split data into training and validation set for hyperparameter tuning
  train_size = round((1-val_split)*length(Y))
  # set.seed(seed)
  ind_train = sample(1:length(Y), train_size)
  ind_val = setdiff(1:length(Y), ind_train)
  X_train = as.matrix(X[ind_train,])
  Y_train_onehot = Y_onehot[ind_train, ]
  X_val = as.matrix(X[ind_val,])
  Y_val_onehot = Y_onehot[ind_val, ]

  ### tune hyperparameters
  # normalize data using min-max normalization
  # estimate pre-processing transformation from training set
  preproc_tune = caret::preProcess(X_train, method='range')
  # apply transformation to training and validation sets
  X_train_norm = predict(preproc_tune, X_train)
  X_val_norm = predict(preproc_tune, X_val)
  # run hyperparameter search
  tuning_results = tune_model(X_train_norm, Y_train_onehot, X_val_norm, Y_val_onehot, trials=trials, epochs=epochs, verbose_keras=0, verbose_mbo=verbose_mbo, seed=seed)

  ### fit final model using training and validation data
  # normalize data using min-max normalization
  preproc = caret::preProcess(as.matrix(X), method='range')
  X_norm = predict(preproc, as.matrix(X))
  # train model
  input_shape = ncol(X)
  model = create_model(learning_rate = 10^tuning_results$x$learning_rate, weight_decay = 10^tuning_results$x$weight_decay, dropout = 10^tuning_results$x$dropout, num_dense_layers = tuning_results$x$num_dense_layers, num_dense_nodes = tuning_results$x$num_dense_nodes, activation = tuning_results$x$activation, seed=seed, input_shape, num_classes)
  history = fit(model, X_norm, Y_onehot, epochs=epochs, batch_size=32, verbose=0)
  # evaluate(model, X_norm, Y_onehot)

  # return list with X, Y, mapping for Y, final neural network model, indices of validation samples, tuning results, min-max transformation for X
  fit <- new_nn(
    # X=X,
    # Y=Y,
    map_df=map_df,
    model=model,
    model_raw = serialize_model(model),
    ind_val=ind_val,
    tuning_results=tuning_results,
    preproc=preproc
  )

  # return()
  out <- list(
    X = X,
    Y = Y,
    fit = fit
  )

  out
}



#' Get neural network predictions
#'
#' @inheritParams fit_smlc
#' @param fit Fitted neural network hidden genome classifier
#' (output of fit_nnc())
#'
#' @return list with entries:
#'
#' @export
predict_nnc = function(fit,
                       Xnew,
                       Ynew = NULL, ...) {
  # if (!inherits(object, "nn")) {
  #   warning("\"object\" should be of a class inheriting from \"nn\"")
  # }
  # apply transformation

  object <- fit$fit
  Xold_names <- colnames(fit$X)
  newdata <- Xnew_adj <- Xnew %>%
    fill_sparsemat_zero(
      rownames = rownames(.),
      colnames = Xold_names
    )

  newdata_norm = predict(object$preproc, as.matrix(newdata))
  # get predicted probabilities
  pred_prob = object$model_raw %>%
    unserialize_model() %>%
    predict(newdata_norm) %>%
    magrittr::set_colnames(
      as.character(
        object$map_df$name[order(object$map_df$num)]
      )
    ) %>%
    magrittr::set_rownames(
      rownames(Xnew_adj)
    ) %>%
    .[, sort(colnames(.))]


  pred_class <- apply(
    pred_prob,
    1,
    function(x) names(x)[which.max(x)]
  )


  list("predicted" = pred_class,
       "probs_predicted" = pred_prob,
       "observed" = Ynew)

  # if (type=="response") {
  #   return(pred_prob)
  # } else if (type=="class") {
  #   # get class assignments
  #   pred_class = max.col(pred_prob, ties.method="first")
  #   pred_df = data.frame(num=pred_class)
  #   pred_df = left_join(pred_df, object$map_df, by="num")
  #   return(pred_df$name)
  # } else {
  #   stop("\"type\" must be \"response\" or \"class\"")
  # }

}


#' #' Create "predict" method for objects of class "nn"
#' setOldClass("nn")
#' setMethod("predict", signature="nn",
#'           function(object, newdata, type) {
#'             predict_nn(object, newdata, type)
#'           })

