#' Print smlc fit results
#' @export
print_single <- function(out) {

  if (length(unique(na.omit(out$obs))) == 2) {
    a0 = pROC::auc(out$obs, out$pred)
    print(a0)
  }
  if (all(is.numeric(out$pred), is.numeric(out$obs))) {
    a1 = caret::confusionMatrix(data = round(out$pred) %>%
                                  factor(levels = c(1, 0)),
                                reference = out$obs %>%
                                  factor(levels = c(1, 0)))
    print(a1)
  } else {
    unique_obs <- sort(unique(out$obs))
    obs <- out$obs %>% unname() %>% factor(., levels = unique_obs)
    pred <- out$pred %>% unname() %>% factor(., levels = unique_obs)
    a1 = caret::confusionMatrix(data = pred,
                                reference = obs)

    ss = capture.output(print(a1))
    cat(ss[1:23], sep = "\n")
    t(a1$byClass) %>%
      round(3) %>%
      print()

  }
}
