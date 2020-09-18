#' Cross validation smlc experiments with n_folds
#' INPUT:
#' maf_dt_in: full mutation annotation file (maf)
#' all_covmats: list of covariance matrices
#' (XU combined) for patients
#' with colum names
#' g = gene, v = variants, c = cancer type,
#' and p = patient (tumor) id. Will be used to compute the
#' design matrix of recorded variants in each cross validation fold.
#' XU_g, XU_ms: desgin-meta-design matrix with genes and mutation
#' signatures respectively.
#' @export

do_cv_multiclass_allmodels <- function(
  maf_dt_in,
  all_covmats,
  covmat_append_unsc_vardesign =
    rep(FALSE, length(all_covmats)),
  covmat_append_sc_vardesign =
    rep(FALSE, length(all_covmats)),
  covmat_append_sqrt_TMB =
    rep(FALSE, length(all_covmats)),
  Y,
  sqrt_TMB,
  g_list,
  v_rank_thresh_full = 800,
  n_fold = 5,
  alpha = 1,
  grouped = FALSE,
  mi_normalized = TRUE,
  recorded_var_classifier = FALSE,
  scale_recorded_var_classifer = FALSE,
  ...
)
{

  # browser()
  c_all <- unique(maf_dt_in$c)
  tcga_dt <- data.table::as.data.table(maf_dt_in)


  # find a random partition via stratified random sampling
  folds <- unique(
    tcga_dt[, .(p, c)]
  )[,
    .SD[sample(1:.N)],
    by = c
    ][,
      fold_indic := rep(1:n_fold, length.out = .N),
      by = c
      ][,
        .(p = .(unique(p))),
        by = fold_indic
        ]$p


  allfits <- allpredicts <- vector("list", n_fold)

  # on each fold, first find the training and test sets
  # then fit on the training set, predict on the test set

  for (ii in seq_len(n_fold)) {

    cat("\n***********************************************")
    cat("\nfold =", ii)
    cat("\n***********************************************")

    test_p <- folds[ii] %>% unlist()
    train_p <- folds[-ii] %>% unlist()

    dt_train <- tcga_dt[p %in% train_p]
    dt_test <- tcga_dt[p %in% test_p]


    cat("\n\nscreening v based on NMI ranking..")

    # browser()

    v_mi <- variant_screen_mi(dt_train, v_rank_thresh, normalize_mi = mi_normalized)

    v_train_list_hier <- v_mi$probs_mi[nmi_rank <= v_rank_thresh]$v
    v_train_list_recorded <- v_mi$probs_mi[nmi_rank <= v_rank_thresh_full]$v


    cat("\n\nfinding training design matrices..")

    this_Xv_hier_train_unsc <- extract_design(
      dt_train,
      "v",
      col_select = v_train_list_hier
    )[train_p, ]

    this_Xv_hier_train_sc <- this_Xv_hier_train_unsc * (1/sqrt_TMB[train_p])

    this_Xv_recorded_train_unsc <- extract_design(
      dt_train,
      "v",
      col_select = v_train_list_recorded
    )[train_p, ]


    this_train_covmats <- purrr::pmap(
      list(all_covmats,
           covmat_append_sc_vardesign,
           covmat_append_unsc_vardesign,
           covmat_append_sqrt_TMB),
      function(this_covmat, this_append_sc_vardesign,
               this_append_unsc_vardesign, this_append_sqrt_TMB) {
        out <- this_covmat[train_p, ]

        if (this_append_sc_vardesign) {
          out <- cbind(out, this_Xv_hier_train_sc[train_p, ])
        }

        if (this_append_unsc_vardesign) {
          out <- cbind(out, this_Xv_hier_train_unsc[train_p, ])
        }

        if (this_append_sqrt_TMB) {
          out <- cbind(out, sqrt_TMB = sqrt_TMB[train_p])
        }

        out
      }
    )


    # browser()

    if (recorded_var_classifier) {
      if (scale_recorded_var_classifer) {
        this_train_covmats$recorded_vars <-
          this_Xv_recorded_train_unsc[train_p, ] * (1/sqrt_TMB[train_p])
      } else {
        this_train_covmats$recorded_vars <-
          this_Xv_recorded_train_unsc[train_p, ]
      }

    }


    # fit smlc models on the train set design matrices

    # allfits[[ii]]
    rand_foldid <- get_rand_foldid(Y[train_p], 10)
    this_allfits <- purrr::imap(
      this_train_covmats,
      function(covmat, mod) {
        cat("\n\n")
        cat(
          glue::glue(
            "training model = \'{mod}\' in fold = {ii}..."
          )
        )


        out <- fit_smlc(X = covmat[train_p, ],
                        Y = Y[train_p],
                        alpha = alpha,
                        grouped = grouped,
                        foldid = rand_foldid,
                        # no_penalty_vars = no_penalty_vars,
                        ...)
        out

      }
    )


    # cat("\n\n")
    # cat(
    #   glue::glue(
    #     "training model = \'rf_g\' in fold = {ii}..."
    #   )
    # )
    # allfits[[ii]]$rf_g <- fit_rfc(
    #   X = train_covmats$sml_g,
    #   Y = Y[train_p]
    # )
    #
    #
    #
    # cat("\n\nfinding test design matrix for g..")
    # Xg_test <- extract_design(dt_test, "g", g_list)
    #
    # test_covmats <- list(
    #   sml_v = Xv_test[test_p, ],
    #   sml_g = Xg_test[test_p, ],
    #   smml_v_g_ms = cbind(Xv_test[test_p, ],
    #                       XU_g[test_p, ],
    #                       XU_ms[test_p, ])
    # )

    # browser()


    # test set design matrics

    cat("\n\nfinding test design matrices..")

    test_train_common_v <- v_train_list_recorded %>%
      union(v_train_list_hier) %>%
      intersect(dt_test$v %>% unique())

    this_Xv_hier_test_unsc <- extract_design(
      dt_test,
      "v",
      col_select = test_train_common_v
    )
    this_Xv_hier_test_sc <- this_Xv_hier_test_unsc * (1/sqrt_TMB[test_p])


    this_test_covmats <- purrr::pmap(
      list(all_covmats,
           covmat_append_sc_vardesign,
           covmat_append_unsc_vardesign,
           covmat_append_sqrt_TMB),
      function(this_covmat, this_append_sc_vardesign,
               this_append_unsc_vardesign, this_append_sqrt_TMB) {
        out <- this_covmat[test_p, ]

        if (this_append_sc_vardesign) {
          out <- cbind(out, this_Xv_hier_test_sc[test_p, ])
        }

        if (this_append_unsc_vardesign) {
          out <- cbind(out, this_Xv_hier_test_unsc[test_p, ])
        }

        if (this_append_sqrt_TMB) {
          out <- cbind(out, sqrt_TMB = sqrt_TMB[test_p])
        }

        out
      }
    )

    # this_test_covmats[[1]] <- cbind(
    #   sqrt_TMB = sqrt(sqrt_TMB[test_p]),
    #   Xv_test[test_p, ] * 1/sqrt_TMB[test_p],
    #   this_test_covmats[[1]][test_p, ]
    # )

    if (recorded_var_classifier) {
      if (scale_recorded_var_classifer) {
        this_test_covmats$recorded_vars <-
          this_Xv_hier_test_unsc[test_p, ] * (1/sqrt_TMB[test_p])
      } else {
        this_test_covmats$recorded_vars <- this_Xv_hier_test_unsc[test_p, ]
      }
    }
    # finding mutation burden for testing set tumors
    # tcga_MB_test <- dt_test[
    #   ,
    #   .(MB = length(unique(v))),
    #   by = p
    #   ]
    #
    # MB_p_test <- tcga_MB_test$MB %>%
    #   magrittr::set_names(tcga_MB_test$p)
    #
    # test_raw_covmats <- list(
    #   # raw_v = Xv_test_big[test_p, ],
    #   # raw_g = Xg_test[test_p, ],
    #   raw_v_g_sbs = cbind(Xv_test[test_p, ],
    #                       XU_g[test_p, ],
    #                       XU_ms[test_p, ])
    # )
    #
    # test_raw_covmats$raw_v_g_sbs_window <- cbind(test_raw_covmats$raw_v_g_sbs,
    #                                              XU_window[test_p, ])
    #
    #
    # test_rescaled_covmats <- test_raw_covmats %>%
    #   lapply(
    #     function(this_mat) {
    #       cbind(
    #         # MB = MB_p_test[test_p],
    #         # # sweep(this_mat, 1, MB_p_test[test_p], "/")
    #         # this_mat * (1/MB_p_test[test_p])
    #         sqrt_MB = sqrt(MB_p_test[test_p]),
    #         this_mat * (1/sqrt(MB_p_test[test_p]))
    #       )
    #     }
    #   ) %>%
    #   magrittr::set_names(
    #     names(.) %>%
    #       stringr::str_replace("raw_", "rescaled_MB_")
    #   )

    # test_covmats <- c(test_raw_covmats, test_rescaled_covmats)


    # browser()

    allpredicts[[ii]] <- purrr::pmap(
      list(names(this_test_covmats),
           this_test_covmats[names(this_test_covmats)],
           this_allfits[names(this_test_covmats)]),
      function(mod, covmat, fit) {
        # covmat <- test_covmats[[mod]]
        # fit <- allfits[[ii]][[mod]]
        cat("\n\n")
        cat(
          glue::glue(
            "predicting model = \'{mod}\' in fold = {ii}..."
          )
        )
        predict_smlc(
          Xnew = covmat,
          fit = fit,
          Ynew = Y[test_p]
        )
      }
    ) %>%
      magrittr::set_names(names(this_test_covmats))


    # cat("\n\n")
    # cat(
    #   glue::glue(
    #     "predicting model = \'rf_g\' in fold = {ii}..."
    #   )
    # )
    #
    # allpredicts[[ii]]$rf_g <- predict_rfc(
    #   Xnew = test_covmats$sml_g,
    #   Ynew = Y[test_p],
    #   fit = allfits[[ii]]$rf_g
    # )


    res <- allpredicts[[ii]] %>%
      sapply(function(xx) mean(xx$predicted == xx$observed))


    cat("\n\n")
    cat("\nHard prediction accuracy in this fold:\n\n")
    glue::glue(
      "{names(res)}: {round(res, 3)}"
    ) %>%
      cat(sep = "\n")

    cat("\n***********************************************")

    cat("\n\n")
  }


  models <- names(allpredicts[[1]])

  out <-  lapply(
    models,
    function(x) {
      pred_list <-  lapply(allpredicts, "[[", x)
      predicted <- lapply(pred_list, '[[', "predicted") %>%
        unlist()
      observed <- lapply(pred_list, '[[', "observed") %>%
        unlist()
      probs_predicted <- lapply(pred_list, '[[', "probs_predicted") %>%
        do.call(rbind, .)
      allfits <- NULL #lapply(allfits, "[[", x)


      list(predicted = predicted,
           observed = observed,
           probs_predicted = probs_predicted,
           fit = allfits
      )
    }
  )

  names(out) <- models

  out
}







# do_cv_multiclass <- function(maf_dt_in,
#                              XU_g,
#                              XU_ms,
#                              Y,
#                              g_list,
#                              v_rank_thresh = 100,
#                              v_rank_thresh_full = 800,
#                              n_fold = 5,
#                              alpha = 1,
#                              grouped = FALSE,
#                              ...)
# {
#   c_all <- unique(maf_dt_in$c)
#   tcga_dt <- data.table::as.data.table(maf_dt_in)
#
#
#   # find a random partition via stratified random sampling
#   folds <- unique(
#     tcga_dt[, .(p, c)]
#   )[,
#     .SD[sample(1:.N)],
#     by = c
#     ][,
#       fold_indic := rep(1:n_fold, length.out = .N),
#       by = c
#       ][,
#         .(p = .(unique(p))),
#         by = fold_indic
#         ]$p
#
#
#   allfits <- allpredicts <- vector("list", n_fold)
#
#   # on each fold, first find the training and test sets
#   # then fit on the training set, predict on the test set
#
#   for (ii in seq_len(n_fold)) {
#
#     cat("\n***********************************************")
#     cat("\nfold =", ii)
#     cat("\n***********************************************")
#
#     test_p <- folds[ii] %>% unlist()
#     train_p <- folds[-ii] %>% unlist()
#
#     dt_train <- tcga_dt[p %in% train_p]
#     dt_test <- tcga_dt[p %in% test_p]
#
#
#     cat("\n\nscreening v based on NMI ranking..")
#
#     v_mi <- variant_screen_mi(dt_train, 100)
#
#
#     cat("\n\nfinding training design matrix for v..")
#     # for multilevel model
#     Xv_train <- extract_design(dt_train, "v", v_mi$v_select)
#     # for recoreded variant only model
#     Xv_train_big <- extract_design(
#       dt_train, "v",
#       v_mi$probs_mi[nmi_rank <= v_rank_thresh_full]$v
#     )
#
#
#     cat("\n\nfinding training design matrix for g..")
#     Xg_train <- extract_design(dt_train, "g", g_list)
#
#     # names:
#     # sml_v = sparse multinomial logisitic with recorded variants
#     # sml_g = sparse multinomial logisitic with genes
#     # smml_v_g_ms = sparse multilevel multinomial logisitic with
#     # variants, genes, mutation signatures
#     # rf_g = random forest with genes
#
#     train_covmats <- list(
#       sml_v = Xv_train_big[train_p, ],
#       sml_g = Xg_train[train_p, ],
#       smml_v_g_ms = cbind(Xv_train[train_p, ],
#                           XU_g[train_p, ],
#                           XU_ms[train_p, ])
#     )
#
#
#     allfits[[ii]] <- purrr::imap(
#       train_covmats,
#       function(covmat, mod) {
#         cat("\n\n")
#         cat(
#           glue::glue(
#             "training model = \'{mod}\' in fold = {ii}..."
#           )
#         )
#         fit_smlc(X = covmat,
#                  Y = Y[train_p],
#                  alpha = alpha,
#                  grouped = grouped,
#                  ...)
#
#       }
#     )
#
#
#     cat("\n\n")
#     cat(
#       glue::glue(
#         "training model = \'rf_g\' in fold = {ii}..."
#       )
#     )
#     allfits[[ii]]$rf_g <- fit_rfc(
#       X = train_covmats$sml_g,
#       Y = Y[train_p]
#     )
#
#
#     cat("\n\nfinding test design matrix for v..")
#     Xv_test <- extract_design(dt_test, "v")
#
#     cat("\n\nfinding test design matrix for g..")
#     Xg_test <- extract_design(dt_test, "g", g_list)
#
#     test_covmats <- list(
#       sml_v = Xv_test[test_p, ],
#       sml_g = Xg_test[test_p, ],
#       smml_v_g_ms = cbind(Xv_test[test_p, ],
#                           XU_g[test_p, ],
#                           XU_ms[test_p, ])
#     )
#
#
#
#
#     allpredicts[[ii]] <- lapply(
#       names(test_covmats),
#       function(mod) {
#         covmat <- test_covmats[[mod]]
#         fit <- allfits[[ii]][[mod]]
#         cat("\n\n")
#         cat(
#           glue::glue(
#             "predicting model = \'{mod}\' in fold = {ii}..."
#           )
#         )
#         predict_smlc(
#           Xnew = covmat,
#           fit = fit,
#           Ynew = Y[test_p]
#         )
#       }
#     )
#     names(allpredicts[[ii]]) <- names(test_covmats)
#
#
#     cat("\n\n")
#     cat(
#       glue::glue(
#         "predicting model = \'rf_g\' in fold = {ii}..."
#       )
#     )
#
#     allpredicts[[ii]]$rf_g <- predict_rfc(
#       Xnew = test_covmats$sml_g,
#       Ynew = Y[test_p],
#       fit = allfits[[ii]]$rf_g
#     )
#
#
#     cat("\n***********************************************")
#
#     cat("\n\n")
#   }
#
#
#   models <- names(allpredicts[[1]])
#
#   out <-  lapply(
#     models,
#     function(x) {
#       pred_list <-  lapply(allpredicts, "[[", x)
#       predicted <- lapply(pred_list, '[[', "predicted") %>%
#         unlist()
#       observed <- lapply(pred_list, '[[', "observed") %>%
#         unlist()
#       probs_predicted <- lapply(pred_list, '[[', "probs_predicted") %>%
#         do.call(rbind, .)
#       allfits <- lapply(allfits, "[[", x)
#
#
#       list(predicted = predicted,
#            observed = observed,
#            probs_predicted = probs_predicted,
#            fit = allfits
#       )
#     }
#   )
#
#   names(out) <- models
#
#   out
# }
#
#
#
# do_cv_multiclass_impute <- function(maf_dt_in,
#                                     g_list,
#                                     g_target,
#                                     v_rank_thresh = 250,
#                                     v_rank_thresh_full = 800,
#                                     n_fold = 5,
#                                     gaussian_alpha = 1,
#                                     alpha = 1,
#                                     grouped = FALSE,
#                                     ...)
# {
#   c_all <- unique(maf_dt_in$c)
#   tcga_dt <- data.table::as.data.table(maf_dt_in)
#
#
#   # find a random partition via stratified random sampling
#   folds <- unique(
#     tcga_dt[, .(p, c)]
#   )[,
#     .SD[sample(1:.N)],
#     by = c
#     ][,
#       fold_indic := rep(1:n_fold, length.out = .N),
#       by = c
#       ][,
#         .(p = .(unique(p))),
#         by = fold_indic
#         ]$p
#
#
#   allfits <- allpredicts <- vector("list", n_fold)
#
#   # on each fold, first find the training and test sets
#   # then fit on the training set, predict on the test set
#
#   for (ii in seq_len(n_fold)) {
#
#     cat("\n***********************************************")
#     cat("\nfold =", ii)
#     cat("\n***********************************************")
#
#     test_p <- folds[ii] %>% unlist() %>% sort
#     train_p <- folds[-ii] %>% unlist() %>% sort
#
#     dt_train <- tcga_dt[p %in% train_p]
#     dt_test <- tcga_dt[p %in% test_p]
#
#     cat("\n\nscreening v based on NMI ranking..")
#
#     v_mi_list <- list(
#       train_full = dt_train,
#       train_target = dt_train[g %in% g_target]
#     ) %>%
#       lapply(
#         function(this_dt)
#           variant_screen_mi(this_dt, v_rank_thresh)
#       )
#
#
#     cat("\n\nfinding training design matrix for v..")
#     # for multilevel model
#
#     Xv_train_list <- purrr::pmap(
#       list(
#         list(train_full = dt_train,
#              train_target = dt_train[g %in% g_target]),
#         v_mi_list,
#         # grepl("full", names(v_mi_list))
#         rep(FALSE, length(v_mi_list))
#       ),
#       function(this_dt, this_v_mi, this_logical) {
#         tmp <- extract_design(
#           dat = dt_train,
#           column =  "v",
#           col_select = this_v_mi$v_select,
#           calc_nontarget_p = this_logical,
#           target_g = g_target
#         )
#
#         nontarget_p <- attr(tmp, "p_v")
#
#         absent_p <- setdiff(train_p, rownames(tmp))
#         if (length(absent_p) > 0) {
#           zeromat <- Matrix::Matrix(
#             data = 0,
#             nrow = length(absent_p),
#             ncol = ncol(tmp),
#             dimnames = list(absent_p, colnames(tmp))
#           )
#           mat <- rbind(tmp, zeromat)[train_p, ]
#         } else {
#           mat <- tmp[train_p, ]
#         }
#
#
#         mat
#       }
#     )
#
#     Xv_train_list$impute_sbs <- Xv_train_list$train_full_target_g <- Xv_train_list$train_full
#
#     # nontarget_p <- Xv_train_list$train_full$nontarget_p
#
#
#     cat("\n\nfinding training design-metadesign g..")
#
#     XU_g_train_list <-  list(
#       train_full = dt_train,
#       train_target = dt_train[g %in% g_target]
#     ) %>%
#       purrr::map2(
#         list(g_list, g_target),
#         function(this_dt, this_g_list) {
#           tmp <- extract_design_metadesign_metacat(
#             dat = this_dt,
#             metafeat_col = "g",
#             metafeat_levels =  this_g_list,
#             p_list = train_p,
#             target_g = g_target
#           )
#
#           mat <- Matrix::Matrix(0, nrow = length(train_p),
#                                 ncol = length(this_g_list),
#                                 dimnames = list(train_p, this_g_list),
#                                 sparse = TRUE)
#
#           mat[rownames(tmp), colnames(tmp)] <- tmp[rownames(tmp), colnames(tmp)]
#
#           colnames(mat) <- paste(
#             "g",
#             colnames(mat),
#             sep = "_"
#           )
#
#           mat
#         }
#       )
#
#     XU_g_train_list$train_full_target_g <- XU_g_train_list$impute_sbs <- XU_g_train_list$train_target
#
#
#     SBS_list <- unique(dt_train$SBS) %>% sort()
#
#
#
#     XU_SBS_train_list <-  purrr::map2(
#       list(
#         train_full = dt_train,
#         train_target = dt_train[g %in% g_target]
#       ),
#       list("all", g_target),
#
#       function(this_dt, this_g_subset) {
#         tmp <- extract_design_metadesign_metacat(
#           dat = this_dt,
#           metafeat_col = "SBS",
#           metafeat_levels =  "all",
#           calc_nontarget_p.U = FALSE,
#           target_g = this_g_subset
#         )
#
#         mat <- Matrix::Matrix(0, nrow = length(train_p),
#                               ncol = length(SBS_list),
#                               dimnames = list(train_p, SBS_list),
#                               sparse = TRUE)
#
#         mat[rownames(tmp), colnames(tmp)] <- tmp
#
#         colnames(mat) <- paste(
#           "SBS",
#           colnames(mat),
#           sep = "_"
#         )
#
#         mat
#       }
#     )
#
#
#     XU_SBS_train_list$train_full_target_g <- XU_SBS_train_list$impute_sbs <- XU_SBS_train_list$train_full
#
#
#
#     cat("\n \ntraining linear imputation model for SBS..")
#     # TMB from target, predictor in imputation
#     train_Zmat <- XU_SBS_train_list$train_target
#
#     # TMB from non target, response in imputation
#     train_Ymat <- XU_SBS_train_list$train_full %>%
#       .[rownames(train_Zmat), colnames(train_Zmat)] -
#       train_Zmat
#
#     train_impute <- fit_lm_SBS(Ymat = train_Ymat,
#                                Zmat = train_Zmat)
#
#
#
#     cat("\n\nfinding training cancer labels..")
#     Y_train <- extract_cancer(dt_train)
#
#
#     this_allfits <- purrr::pmap(
#       list(Xv_train_list[c('train_full', 'train_full_target_g', 'train_target')],
#            XU_g_train_list[c('train_full', 'train_full_target_g', 'train_target')],
#            XU_SBS_train_list[c('train_full', 'train_full_target_g', 'train_target')],
#            c('train_full', 'train_full_target_g', 'train_target')),
#       function(Xv, XU_g, XU_SBS, method) {
#         cat("\n\n")
#         cat(
#           glue::glue(
#             "training model = \'{method}\' in fold = {ii}..."
#           )
#         )
#         fit_smlc(X = cbind(Xv[train_p, ],
#                            XU_g[train_p, ],
#                            XU_SBS[train_p, ]),
#                  Y = Y_train[train_p],
#                  alpha = alpha,
#                  grouped = grouped, ...)
#
#
#       }
#     )
#
#     cat("\n\n")
#
#     cat(
#       glue::glue(
#         "trained model \'impute_sbs\' is the same as \'train_full_target_g\'..."
#       )
#     )
#
#     this_allfits$impute_sbs <- this_allfits$train_full_target_g
#
#
#     allfits[[ii]] <- this_allfits
#
#
#
#     #
#
#
#
#     cat("\n\nfinding test design matrix for v..")
#     Xv_test_list <- list(
#       train_full = dt_test,
#       train_target = dt_test[g %in% g_target]
#     ) %>%
#       lapply(
#         function(this_dt) {
#
#           test_v_list <- unique(this_dt$v) %>% unique()
#
#           tmp <- extract_design(
#             dat = this_dt,
#             column = "v",
#             calc_nontarget_p = FALSE
#           )
#
#           mat <- Matrix::Matrix(0, nrow = length(test_p),
#                                 ncol = length(test_v_list),
#                                 dimnames = list(test_p, test_v_list),
#                                 sparse = TRUE)
#
#           mat[rownames(tmp), colnames(tmp)] <- tmp
#
#           mat
#         }
#       )
#
#     Xv_test_list$impute_sbs <- Xv_test_list$train_target
#     Xv_test_list$train_full_target_g <- Xv_test_list$train_full
#
#     cat("\n\nfinding test design metadesign matrix for g..")
#     XU_g_test_list <- list(
#       train_full = dt_test,
#       train_target = dt_test[g %in% g_target]
#     ) %>%
#       purrr::map2(
#         list(g_list, g_target),
#         function(this_dt, this_g_list) {
#           tmp <- extract_design_metadesign_metacat(
#             dat = this_dt,
#             metafeat_col = "g",
#             metafeat_levels =  this_g_list,
#             calc_nontarget_p.U = FALSE,
#             target_g = g_target
#           )
#
#
#           mat <- Matrix::Matrix(0, nrow = length(test_p),
#                                 ncol = length(this_g_list),
#                                 dimnames = list(test_p, this_g_list),
#                                 sparse = TRUE)
#
#           mat[rownames(tmp), colnames(tmp)] <- tmp
#
#           colnames(mat) <- paste(
#             "g",
#             colnames(mat),
#             sep = "_"
#           )
#
#           mat
#         }
#       )
#
#     XU_g_test_list$train_full_target_g <- XU_g_test_list$impute_sbs <- XU_g_test_list$train_target
#
#
#
#     cat("\n\nfinding test design metadesign matrix for SBS..")
#
#     sbs_list_train <- dt_train$SBS %>% unique %>% sort
#     XU_SBS_test_list <- list(
#       train_full = dt_test,
#       train_target = dt_test[g %in% g_target]
#     ) %>%
#       lapply(
#         function(this_dt) {
#           tmp <- extract_design_metadesign_metacat(
#             dat = this_dt,
#             metafeat_col = "SBS",
#             metafeat_levels =  "all",
#             calc_nontarget_p.U = FALSE
#           )
#
#           mat <- Matrix::Matrix(0, nrow = length(test_p),
#                                 ncol = length(sbs_list_train),
#                                 dimnames = list(test_p, sbs_list_train),
#                                 sparse = TRUE)
#
#           mat[rownames(tmp), colnames(tmp)] <- tmp
#
#           colnames(mat) <- paste(
#             "SBS",
#             colnames(mat),
#             sep = "_"
#           )
#
#           mat
#         }
#       )
#
#     XU_SBS_test_list$train_full_target_g <- XU_SBS_test_list$train_full
#
#     cat("\n\nExtract cancer label from test data for validation")
#     Y_test <- extract_cancer(dt_test)
#
#
#     cat("\n\nImputing desgian-meta-design for SBS in test set\n based on trained imputation model...")
#     Zmat_test <- XU_SBS_test_list$train_target[, colnames(train_impute$Xmat)]
#     Ymat_test <- lm_predict_SBS(fit = train_impute,
#                                 Zmat_new = Zmat_test) %>%
#       .[rownames(Zmat_test), colnames(Zmat_test)]
#
#     XU_SBS_test_list$impute_sbs <- Matrix::Matrix(Zmat_test + Ymat_test, sparse = TRUE)
#
#
#     allpredicts[[ii]] <- this_allpredicts <- this_allfits %>%
#       list(
#         .,
#         Xv_test_list[names(.)],
#         XU_g_test_list[names(.)],
#         XU_SBS_test_list[names(.)],
#         names(.)
#       ) %>%
#       purrr::pmap(
#         function(fit, Xv, XU_g, XU_SBS, method) {
#           cat("\n\n")
#           cat(
#             glue::glue(
#               "predicting model = \'{method}\' in fold = {ii}..."
#             )
#           )
#
#           predict_smlc(
#             Xnew = cbind(Xv[test_p, ], XU_g[test_p, ], XU_SBS[test_p, ]),
#             fit = fit,
#             Ynew = Y_test[test_p],
#             impute = FALSE,
#           )
#
#         }
#
#       )
#
#     res <- allpredicts[[ii]] %>%
#       sapply(function(xx) mean(xx$predicted == xx$observed))
#
#     cat("\nHard prediction accuracy in this fold:\n\n")
#     glue::glue(
#       "{names(res)}: {round(res, 3)}"
#     ) %>%
#       cat(sep = "\n")
#
#     #
#     #
#     # cat("\n\n")
#     # cat(
#     #   glue::glue(
#     #     "predicting model = \'rf_g\' in fold = {ii}..."
#     #   )
#     # )
#     #
#     # allpredicts[[ii]]$rf_g <- predict_rfc(
#     #   Xnew = test_covmats$sml_g,
#     #   Ynew = Y[test_p],
#     #   fit = allfits[[ii]]$rf_g
#     # )
#     #
#
#     cat("\n***********************************************")
#
#     cat("\n\n")
#   }
#
#
#   models <- names(allpredicts[[1]])
#
#   out <- sapply(
#     models,
#     function(x) {
#       pred_list <-  lapply(allpredicts, "[[", x)
#       predicted <- lapply(pred_list, '[[', "predicted") %>%
#         unlist()
#       observed <- lapply(pred_list, '[[', "observed") %>%
#         unlist()
#       probs_predicted <- lapply(pred_list, '[[', "probs_predicted") %>%
#         do.call(rbind, .)
#       allfits <- lapply(allfits, "[[", x)
#
#
#       list(predicted = predicted,
#            observed = observed,
#            probs_predicted = probs_predicted,
#            fit = allfits
#       )
#     },
#     simplify = FALSE,
#     USE.NAMES = TRUE
#   )
#
#
#   out
# }
#
#
#
# do_cv_multiclass_scaled_window <- function(maf_dt_in,
#                                            XU_g,
#                                            XU_ms,
#                                            XU_window,
#                                            Y,
#                                            g_list,
#                                            v_rank_thresh = 100,
#                                            v_rank_thresh_full = 800,
#                                            n_fold = 5,
#                                            alpha = 1,
#                                            grouped = FALSE,
#                                            ...)
# {
#   c_all <- unique(maf_dt_in$c)
#   tcga_dt <- data.table::as.data.table(maf_dt_in)
#
#
#   # find a random partition via stratified random sampling
#   folds <- unique(
#     tcga_dt[, .(p, c)]
#   )[,
#     .SD[sample(1:.N)],
#     by = c
#     ][,
#       fold_indic := rep(1:n_fold, length.out = .N),
#       by = c
#       ][,
#         .(p = .(unique(p))),
#         by = fold_indic
#         ]$p
#
#
#   allfits <- allpredicts <- vector("list", n_fold)
#
#   # on each fold, first find the training and test sets
#   # then fit on the training set, predict on the test set
#
#   for (ii in seq_len(n_fold)) {
#
#     cat("\n***********************************************")
#     cat("\nfold =", ii)
#     cat("\n***********************************************")
#
#     test_p <- folds[ii] %>% unlist()
#     train_p <- folds[-ii] %>% unlist()
#
#     dt_train <- tcga_dt[p %in% train_p]
#     dt_test <- tcga_dt[p %in% test_p]
#
#
#     cat("\n\nscreening v based on NMI ranking..")
#
#     v_mi <- variant_screen_mi(dt_train, 100)
#
#     cat("\n\nfinding training design matrices..")
#     # for multilevel model
#     Xv_train <- extract_design(
#       dt_train,
#       "v",
#       col_select = v_mi$v_select
#     )
#
#
#     # browser()
#     # # for recoreded variant only model
#     # Xv_train_big <- extract_design(
#     #   dt_train, "v",
#     #   v_mi$probs_mi[nmi_rank <= v_rank_thresh_full]$v
#     # )
#
#
#
#
#     # cat("\n\nfinding training design matrix for g..")
#     # Xg_train <- extract_design(dt_train, "g", g_list)
#
#     # names:
#     # sml_v = sparse multinomial logisitic with recorded variants
#     # sml_g = sparse multinomial logisitic with genes
#     # smml_v_g_ms = sparse multilevel multinomial logisitic with
#     # variants, genes, mutation signatures
#     # rf_g = random forest with genes
#
#
#     # training set design matrics
#
#     # finding mutation burden for training set tumors
#     tcga_MB_train <- dt_train[
#       ,
#       .(MB = length(unique(v))),
#       by = p
#       ]
#
#     MB_p_train <- tcga_MB_train$MB %>%
#       magrittr::set_names(tcga_MB_train$p)
#
#     train_raw_covmats <- list(
#       # raw_v = Xv_train_big[train_p, ],
#       # raw_g = Xg_train[train_p, ],
#       raw_v_g_sbs = cbind(Xv_train[train_p, ],
#                           XU_g[train_p, ],
#                           XU_ms[train_p, ])
#     )
#
#     train_raw_covmats$raw_v_g_sbs_window <- cbind(train_raw_covmats$raw_v_g_sbs,
#                                                   XU_window[train_p, ])
#
#
#     train_rescaled_covmats <- train_raw_covmats %>%
#       lapply(
#         function(this_mat) {
#           cbind(
#             sqrt_MB = sqrt(MB_p_train[train_p]),
#             # MB = MB_p_train[train_p],
#             # sweep(this_mat, 1, MB_p_train[train_p], "/")
#             this_mat * (1/sqrt(MB_p_train[train_p]))
#           )
#         }
#       ) %>%
#       magrittr::set_names(
#         names(.) %>%
#           stringr::str_replace("raw_", "rescaled_MB_")
#       )
#
#     train_covmats <- c(train_raw_covmats, train_rescaled_covmats)
#
#
#     # fit smlc models on the train set design matrices
#
#     allfits[[ii]] <- purrr::imap(
#       train_covmats,
#       function(covmat, mod) {
#         cat("\n\n")
#         cat(
#           glue::glue(
#             "training model = \'{mod}\' in fold = {ii}..."
#           )
#         )
#
#
#         no_penalty_vars <- "sqrt_MB" %>%
#           .[grepl(., colnames(covmat))]
#
#         # if (any(penalty_factor == 0)) browser()
#
#         out <- fit_smlc(X = covmat,
#                         Y = Y[train_p],
#                         alpha = alpha,
#                         grouped = grouped,
#                         # no_penalty_vars = no_penalty_vars,
#                         ...)
#         out
#
#       }
#     )
#
#
#
#
#     # test set design matrics
#
#     cat("\n\nfinding test design matrices..")
#
#     Xv_test <- extract_design(dt_test, "v")
#
#     # finding mutation burden for testing set tumors
#     tcga_MB_test <- dt_test[
#       ,
#       .(MB = length(unique(v))),
#       by = p
#       ]
#
#     MB_p_test <- tcga_MB_test$MB %>%
#       magrittr::set_names(tcga_MB_test$p)
#
#     test_raw_covmats <- list(
#       # raw_v = Xv_test_big[test_p, ],
#       # raw_g = Xg_test[test_p, ],
#       raw_v_g_sbs = cbind(Xv_test[test_p, ],
#                           XU_g[test_p, ],
#                           XU_ms[test_p, ])
#     )
#
#     test_raw_covmats$raw_v_g_sbs_window <- cbind(test_raw_covmats$raw_v_g_sbs,
#                                                  XU_window[test_p, ])
#
#
#     test_rescaled_covmats <- test_raw_covmats %>%
#       lapply(
#         function(this_mat) {
#           cbind(
#             # MB = MB_p_test[test_p],
#             # # sweep(this_mat, 1, MB_p_test[test_p], "/")
#             # this_mat * (1/MB_p_test[test_p])
#             sqrt_MB = sqrt(MB_p_test[test_p]),
#             this_mat * (1/sqrt(MB_p_test[test_p]))
#           )
#         }
#       ) %>%
#       magrittr::set_names(
#         names(.) %>%
#           stringr::str_replace("raw_", "rescaled_MB_")
#       )
#
#     test_covmats <- c(test_raw_covmats, test_rescaled_covmats)
#
#
#
#     allpredicts[[ii]] <- purrr::pmap(
#       list(names(test_covmats),
#            test_covmats[names(test_covmats)],
#            allfits[[ii]][names(test_covmats)]),
#       function(mod, covmat, fit) {
#         # covmat <- test_covmats[[mod]]
#         # fit <- allfits[[ii]][[mod]]
#         cat("\n\n")
#         cat(
#           glue::glue(
#             "predicting model = \'{mod}\' in fold = {ii}..."
#           )
#         )
#         predict_smlc(
#           Xnew = covmat,
#           fit = fit,
#           Ynew = Y[test_p]
#         )
#       }
#     ) %>%
#       magrittr::set_names(names(test_covmats))
#
#
#     # cat("\n\n")
#     # cat(
#     #   glue::glue(
#     #     "predicting model = \'rf_g\' in fold = {ii}..."
#     #   )
#     # )
#     #
#     # allpredicts[[ii]]$rf_g <- predict_rfc(
#     #   Xnew = test_covmats$sml_g,
#     #   Ynew = Y[test_p],
#     #   fit = allfits[[ii]]$rf_g
#     # )
#
#
#     res <- allpredicts[[ii]] %>%
#       sapply(function(xx) mean(xx$predicted == xx$observed))
#
#     cat("\nHard prediction accuracy in this fold:\n\n")
#     glue::glue(
#       "{names(res)}: {round(res, 3)}"
#     ) %>%
#       cat(sep = "\n")
#
#     cat("\n***********************************************")
#
#     cat("\n\n")
#   }
#
#
#   models <- names(allpredicts[[1]])
#
#   out <-  lapply(
#     models,
#     function(x) {
#       pred_list <-  lapply(allpredicts, "[[", x)
#       predicted <- lapply(pred_list, '[[', "predicted") %>%
#         unlist()
#       observed <- lapply(pred_list, '[[', "observed") %>%
#         unlist()
#       probs_predicted <- lapply(pred_list, '[[', "probs_predicted") %>%
#         do.call(rbind, .)
#       allfits <- lapply(allfits, "[[", x)
#
#
#       list(predicted = predicted,
#            observed = observed,
#            probs_predicted = probs_predicted,
#            fit = allfits
#       )
#     }
#   )
#
#   names(out) <- models
#
#   out
# }

