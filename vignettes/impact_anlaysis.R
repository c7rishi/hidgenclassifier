## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----import-------------------------------------------------------------------
set.seed(42)
library(magrittr)
library(hidgenclassifier)
data("impact")

## ---- view_data---------------------------------------------------------------
impact

## ----defn, include=FALSE------------------------------------------------------
n_gene <- length(unique(impact$Hugo_Symbol))
n_pid <- length(unique(impact$patient_id))

## ----cancer_sites-------------------------------------------------------------
unique(impact$CANCER_SITE)

## ----extract_cancer-----------------------------------------------------------
canc_resp <- extract_cancer_response(
  maf = impact,
  cancer_col = "CANCER_SITE",
  sample_id_col = "patient_id"
)
pid <- names(canc_resp)

## ----train_test_pid-----------------------------------------------------------
set.seed(42)
folds <- data.table::data.table(
  resp = canc_resp
)[,
  foldid := sample(rep(1:5, length.out = .N)),
  by = resp
]$foldid

# 80%-20% stratified separation of training and
# test set tumors
pid_train <- pid[folds != 5]
pid_test <- pid[folds == 5]

## ----screen_variant-----------------------------------------------------------
top_v <- variant_screen_mi(
  maf = impact[patient_id %in% pid_train],
  variant_col = "Variant",
  cancer_col = "CANCER_SITE",
  sample_id_col = "patient_id",
  mi_rank_thresh = 250,
  return_prob_mi = FALSE,
  do_freq_screen = FALSE
)

## ----design_variant-----------------------------------------------------------
X_variant <- extract_design(
  maf = impact,
  variant_col = "Variant",
  sample_id_col = "patient_id",
  variant_subset = top_v
)
dim(X_variant)

## ----mdesign_gene-------------------------------------------------------------
XU_gene <- extract_design_mdesign_mcat(
  maf = impact,
  variant_col = "Variant",
  mfeat_col = "Hugo_Symbol",
  sample_id_col = "patient_id",
  mfeat_subset = NULL
) %>% 
  magrittr::set_colnames(
    paste("Gene_", colnames(.))
  )
dim(XU_gene)

## ----mdesign_sbs96------------------------------------------------------------
XU_sbs96 <- extract_design_mdesign_sbs96(
  maf = impact,
  chromosome_col = "Chromosome",
  start_position_col = "Start_Position",
  end_position_col = "End_Position",
  ref_col = "Reference_Allele",
  alt_col = "Tumor_Seq_Allele2",
  sample_id_col = "patient_id"
) %>% 
  magrittr::set_colnames(
    paste("SBS_", colnames(.))
  )
dim(XU_sbs96)

## ----tmb----------------------------------------------------------------------
tmb <- extract_tmb(
  maf = impact,
  variant_col = "Variant",
  sample_id_col = "patient_id"
)

## ----combine_pred-------------------------------------------------------------
predictor_mat <- cbind(
  X_variant[pid, ],
  XU_gene[pid, ],
  XU_sbs96[pid, ],
  tmb = tmb[pid]
) %>% 
  divide_rows(sqrt(tmb[pid]))

## ----fit, eval = FALSE--------------------------------------------------------
#  fit_impact <- fit_mlogit(
#    X = predictor_mat[pid_train, ],
#    Y = canc_resp[pid_train]
#  )

## ----predict, eval = FALSE----------------------------------------------------
#  pred_impact <- predict_mlogit(
#    fit = fit_impact,
#    Xnew = predictor_mat[pid_test, ]
#  )

## ----odds_ratio, eval=FALSE---------------------------------------------------
#  
#  or <- odds_ratio_mlogit(
#    fit = fit_impact,
#    type = "one-vs-rest",
#    log = TRUE
#  )

