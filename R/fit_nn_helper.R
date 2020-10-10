#
# # functions for neural network hyperparameter tuning, based on the approach used in Jiao et al. 2020 (https://www.nature.com/articles/s41467-019-13825-8; https://github.com/ICGC-TCGA-PanCancer/TumorType-WGS/blob/master/DNN-Model/train_models_tumour_classifier.py)
#
# library(caret)
# library(mlrMBO)
# library(lhs)
# library(dplyr)
#
# set.seed(1)
# library(keras)
# use_session_with_seed(seed=1, disable_gpu=T, disable_parallel_cpu=T)

#' Initialize a fully-connected multi-class neural network
#'
#' @param learning_rate Learning rate for Adam optimizer
#' @param weight_decay L2 regulariation paramter
#' @param dropout Dropout rate
#' @param num_dense_layers Number of hidden layers
#' @param num_dense_nodes Number of nodes per hidden layer (all hidden layers have the same size)
#' @param activation Name of activation function ("relu", "softmax", "softplus", "elu", "tanh", "sigmoid", etc)
#' @param seed Seed for weight initialization and dropout
#' @param input_shape Number of features
#' @param num_classes Number of classes
#'
#' @return A compiled Keras model
#' @author Zoe Guan
#' @examples
#' model <- create_model(learning_rate = 0.001, weight_decay = 0.5, dropout = 0.2, num_dense_layers = 2, num_dense_nodes = 20, activation = "relu", seed = 1, input_shape = 1324, num_classes = 8)
#' summary(model)
create_model <- function(learning_rate = 0.001,
                         weight_decay = 0,
                         dropout = 0,
                         num_dense_layers = 1,
                         num_dense_nodes = 5,
                         activation = "relu",
                         seed = NULL,
                         input_shape,
                         num_classes) {
  model <- keras_model_sequential()

  if (num_dense_layers == 0) {

    # output layer
    model %>%
      layer_dense(
        units = num_classes,
        activation = "softmax",
        kernel_initializer = initializer_glorot_uniform(seed),
        input_shape = input_shape
      )
  } else {

    # hidden layers
    for (i in 1:num_dense_layers) {
      if (i == 1) {
        model %>%
          layer_dense(
            units = num_dense_nodes,
            activation = activation,
            kernel_regularizer = regularizer_l2(weight_decay),
            kernel_initializer = initializer_glorot_uniform(seed),
            input_shape = input_shape
          ) %>%
          # dropout
          layer_dropout(rate = dropout, seed = seed)
      } else {
        model %>%
          layer_dense(
            units = num_dense_nodes,
            activation = activation,
            kernel_regularizer = regularizer_l2(weight_decay),
            kernel_initializer = initializer_glorot_uniform(seed)
          ) %>%
          # dropout
          layer_dropout(rate = dropout, seed = seed)
      }
    }

    # output layer
    model %>%
      layer_dense(
        units = num_classes,
        kernel_initializer = initializer_glorot_uniform(seed),
        activation = "softmax"
      )
  }

  # optimizer
  optimizer <- optimizer_adam(lr = learning_rate)

  # compile model
  model %>%
    compile(
      loss = "categorical_crossentropy",
      optimizer = optimizer,
      metrics = "accuracy"
    )

  # export model as an R object
  # model <- serialize_model(model)

  return(model)
}


#' Fit and validate neural network
#'
#' @param X_train_norm Normalized input matrix for training set
#' @param Y_train_onehot One-hot representation of response for training set
#' @param X_val_norm Normalized input matrix for validation set
#' @param Y_val_onehot One-hot representation of response for validation set
#' @param epochs Number of training epochs
#' @param batch_size Batch size
#' @param verbose Keras verbosity mode (0 = silent, 1 = progress bar, 2 = one line per epoch)
#'
#' @return Validation accuracy
#'
fitness <- function(X_train_norm,
                    Y_train_onehot,
                    X_val_norm,
                    Y_val_onehot,
                    learning_rate,
                    weight_decay,
                    dropout,
                    num_dense_layers,
                    num_dense_nodes,
                    activation,
                    epochs = 50,
                    batch_size = 128,
                    verbose = 0,
                    seed = NULL) {
  input_shape <- ncol(X_train_norm)
  num_classes <- ncol(Y_train_onehot)

  model <- create_model(
    learning_rate,
    weight_decay,
    dropout,
    num_dense_layers,
    num_dense_nodes,
    activation,
    seed,
    input_shape,
    num_classes
  )

  # train model
  history <- fit(
    model,
    X_train_norm,
    Y_train_onehot,
    epochs = epochs,
    batch_size = batch_size,
    validation_data = list(
      X_val_norm,
      Y_val_onehot
    ),
    verbose = verbose
  )

  # get accuracy in validation set
  accuracy <- history$metrics$val_acc[length(history$metrics$val_acc)]

  return(accuracy)
}



#' Tune neural network hyperparameters using Bayesian optimization
#'
#' @param trials Number of trials
#' @param verbose_keras Keras verbosity mode (0 = silent, 1 = progress bar, 2 = one line per epoch)
#' @param verbose_mbo Bayesian optimization verbosity mode (logical)
#'
#' @return See description of mbo() function from mlrMBO
#'
tune_model <- function(X_train_norm,
                       Y_train_onehot,
                       X_val_norm,
                       Y_val_onehot,
                       trials = 200,
                       epochs = 50,
                       verbose_keras = 0,
                       verbose_mbo = T,
                       seed = NULL,
                       batch_size = 128) {

  # hyperparameters to tune
  param_set <- ParamHelpers::makeParamSet(
    # tune learning_rate, weight_decay, dropout on log scale
    ParamHelpers::makeNumericParam(
      "learning_rate",
      lower = -4, upper = -2, trafo = function(x) 10^x
    ),
    ParamHelpers::makeNumericParam(
      "weight_decay",
      lower = -3,
      upper = log10(0.5),
      trafo = function(x) 10^x
    ),
    ParamHelpers::makeNumericParam(
      "dropout",
      lower = -6,
      upper = log10(0.5),
      trafo = function(x) 10^x
    ),
    ParamHelpers::makeIntegerParam(
      "num_dense_layers",
      lower = 0,
      upper = 5
    ),
    ParamHelpers::makeIntegerParam(
      "num_dense_nodes",
      lower = 5,
      upper = 1024
    ),
    ParamHelpers::makeDiscreteParam(
      "activation",
      values = c("relu", "softplus")
    )
  )

  # objective function
  obj_fn <- smoof::makeSingleObjectiveFunction(
    fn = function(x) {
      fitness(
        X_train_norm,
        Y_train_onehot,
        X_val_norm,
        Y_val_onehot,
        learning_rate = x$learning_rate,
        weight_decay = x$weight_decay,
        dropout = x$dropout,
        num_dense_layers = x$num_dense_layers,
        num_dense_nodes = x$num_dense_nodes,
        activation = as.character(x$activation),
        epochs = epochs,
        verbose = verbose_keras,
        seed = seed,
        batch_size = batch_size
      )
    },
    par.set = param_set,
    minimize = F,
    has.simple.signature = F
  )

  control <- mlrMBO::makeMBOControl() %>%
    # aquisition function: expected improvement
    mlrMBO::setMBOControlInfill(
      crit = mlrMBO::makeMBOInfillCritEI()
    ) %>%
    # set number of trials
    mlrMBO::setMBOControlTermination(
      iters = trials
    )

  # surrogate model: Gaussian process
  surrogate <- mlr::makeLearner(
    "regr.gausspr",
    predict.type = "se",
    config = list(show.learner.output = FALSE)
  )

  # initial random points
  initial_points <- ParamHelpers::generateDesign(
    par.set = ParamHelpers::getParamSet(obj_fn),
    fun = lhs::randomLHS
  )

  # run search
  run <- mlrMBO::mbo(
    fun = obj_fn,
    design = initial_points,
    learner = surrogate,
    control = control,
    show.info = verbose_mbo
  )

  return(run)
}
