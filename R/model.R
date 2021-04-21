data_for_model <- function(data) {
  mutate(data, primer_short = factor(primer_short),
             sex = factor(sex))
}

fit_model <- function(data_model, file_id) {
  priors = c(prior(normal(0,2), class = "sd"),
             prior(normal(0,2), class = "b"))

  # fit1 <- brm(Cq ~ primer_short * mo(genotype_ord) + (1 || run) + (primer_short || sex + animal_no + litter + animal_plate + animal_replicate), data = d17_4w, prior = priors)
  brm(
    bf(Cq ~ primer_short * genotype + sex*primer_short + (1 || run + animal_no + animal_no:primer_short + litter:primer_short),
       sigma ~ (1 || primer_short + run)
    ),
    data = data_model,
    prior = priors, refresh = 500,
    file = here::here(paste0("local_temp_data/fit_", file_id)), file_refit = "on_change")
}

show_pp_checks <- function(data_model, fit) {
  pred <- posterior_predict(fit, cores = 1)
  pp_theme <- theme(axis.text = element_text(size = 6), strip.text = element_text(size = 9))
  suppressWarnings({
    print(ppc_stat_grouped(data_model$Cq, pred, group = data_model$primer_short, stat = "sd") + pp_theme)
    print(ppc_stat_grouped(data_model$Cq, pred, group = data_model$genotype_ord, stat = "sd") + pp_theme)
    for(g in levels(data_model$genotype_ord)) {
      indices = data_model$genotype_ord == g
      print(
        ppc_stat_grouped(data_model$Cq[indices], pred[,indices], group = data_model$primer_short[indices],stat = "sd") +
          ggtitle(paste0("Genotype = ", g)) + pp_theme)
      print(ppc_stat_grouped(data_model$Cq[indices], pred[,indices], group = data_model$primer_short[indices],stat = "mean") +
          ggtitle(paste0("Genotype = ", g)) + pp_theme)
    }
    print(ppc_stat_grouped(data_model$Cq, pred, group = data_model$sex, stat = "sd") + pp_theme)
    print(ppc_stat_grouped(data_model$Cq, pred, group = data_model$run, stat = "sd") + pp_theme)
    print(ppc_stat_grouped(data_model$Cq, pred, group = data_model$run, stat = "mean")  + pp_theme)
  })
}


predict_comparisons <- function(fit, comparisons_eff) {
  pred_data_numerator_wt <- data.frame(primer_short = comparisons_eff$numerator,
                                       genotype = "wt",
                                       sex = "F"
  )
  pred_data_numerator_d17 <- data.frame(primer_short = comparisons_eff$numerator,
                                        genotype = "d17",
                                        sex = "F"
  )
  pred_data_denominator_wt <- data.frame(primer_short = comparisons_eff$denominator,
                                         genotype = "wt",
                                         sex = "F"
  )
  pred_data_denominator_d17 <- data.frame(primer_short = comparisons_eff$denominator,
                                          genotype = "d17",
                                          sex = "F"
  )

  pred_numerator_wt <- posterior_epred(fit, newdata = pred_data_numerator_wt, re_formula = NA, cores = 1)
  pred_numerator_d17 <- posterior_epred(fit, newdata = pred_data_numerator_d17, re_formula = NA, cores = 1)
  pred_denominator_wt <- posterior_epred(fit, newdata = pred_data_denominator_wt, re_formula = NA, cores = 1)
  pred_denominator_d17 <- posterior_epred(fit, newdata = pred_data_denominator_d17, re_formula = NA, cores = 1)

  pred_numerator <- pred_numerator_wt - pred_numerator_d17
  pred_denominator <- pred_denominator_wt - pred_denominator_d17

  comparisons_pred_numerator <- tidybayes::add_draws(comparisons_eff, pred_numerator, value = "pred_numerator")
  comparisons_pred_denominator <- tidybayes::add_draws(comparisons_eff, pred_denominator, value = "pred_denominator")

  comparisons_pred <- comparisons_pred_numerator %>%
    inner_join(comparisons_pred_denominator, by = setdiff(names(comparisons_pred_numerator), "pred_numerator")) %>%
    mutate(log_ratio = log(eff_numerator) * pred_numerator - log(eff_denominator) * pred_denominator)


  if(nrow(comparisons_pred) != nrow(comparisons_pred_numerator) || nrow(comparisons_pred) != nrow(comparisons_pred_denominator)) {
    stop("Bad join")
  }

  comparisons_pred
}
