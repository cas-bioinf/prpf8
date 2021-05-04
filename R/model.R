data_for_model <- function(data, housekeeping) {
  #data <- dplyr::filter(data, !(primer_short %in% housekeeping))
  data <- dplyr::mutate(data, primer_short = factor(primer_short),
             sex = factor(sex))

  data
}

fit_model <- function(data_model, file_id) {
  priors = c(prior(normal(0,2), class = "sd"),
             prior(normal(0,2), class = "b"))

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
      if(!any(indices)) {
        next
      }
      print(
        ppc_stat_grouped(data_model$Cq[indices], pred[,indices], group = data_model$primer_short[indices],stat = "sd") +
          ggtitle(paste0("Genotype = ", g)) + pp_theme)
      print(ppc_stat_grouped(data_model$Cq[indices], pred[,indices], group = data_model$primer_short[indices],stat = "mean") +
          ggtitle(paste0("Genotype = ", g)) + pp_theme)
    }
    print(ppc_stat_grouped(data_model$Cq, pred, group = data_model$sex, stat = "sd") + pp_theme)
    for(s in levels(data_model$sex)) {
      indices = data_model$sex == s
      print(
        ppc_stat_grouped(data_model$Cq[indices], pred[,indices], group = data_model$primer_short[indices],stat = "sd") +
          ggtitle(paste0("Sex = ", s)) + pp_theme)
      print(ppc_stat_grouped(data_model$Cq[indices], pred[,indices], group = data_model$primer_short[indices],stat = "mean") +
              ggtitle(paste0("Sex = ", s)) + pp_theme)
    }
    print(ppc_stat_grouped(data_model$Cq, pred, group = data_model$sex, stat = "sd") + pp_theme)
    print(ppc_stat_grouped(data_model$Cq, pred, group = data_model$run, stat = "sd") + pp_theme)
    print(ppc_stat_grouped(data_model$Cq, pred, group = data_model$run, stat = "mean")  + pp_theme)
    print(ppc_stat_grouped(data_model$Cq, pred, group = data_model$litter, stat = "sd") + pp_theme)
    print(ppc_stat_grouped(data_model$Cq, pred, group = data_model$litter, stat = "mean")  + pp_theme)
  })
}


predict_comparisons <- function(fit, data, comparisons_eff, genotypes1, genotypes2) {
  res_list <- list()
  if(length(genotypes1) != length(genotypes2)) {
    stop("Length must match")
  }
  all_genotypes <- paste0(genotypes2, "/", genotypes1)
  for(i in 1:length(genotypes1)) {
    res_list[[i]] <- predict_comparisons_genotype(fit, data, comparisons_eff, genotypes1[i], genotypes2[i])
    res_list[[i]]$genotypes = factor(all_genotypes[i], levels = all_genotypes)
  }

  do.call(rbind, res_list)
}

predict_comparisons_genotype <- function(fit, data, comparisons_eff, genotype1, genotype2) {
  primers_present_1 <- unique(dplyr::filter(data, genotype == genotype1)$primer_short)
  primers_present_2 <- unique(dplyr::filter(data, genotype == genotype2)$primer_short)
  primers_both <- intersect(primers_present_1, primers_present_2)

  comparisons_eff <- dplyr::filter(comparisons_eff,
                                   numerator %in% primers_both,
                                   denominator %in% primers_both)

  pred_data_numerator_1 <- data.frame(primer_short = comparisons_eff$numerator,
                                       genotype = genotype1,
                                       sex = "F"
  )
  pred_data_numerator_2 <- data.frame(primer_short = comparisons_eff$numerator,
                                        genotype = genotype2,
                                        sex = "F"
  )
  pred_data_denominator_1 <- data.frame(primer_short = comparisons_eff$denominator,
                                         genotype = genotype1,
                                         sex = "F"
  )
  pred_data_denominator_2 <- data.frame(primer_short = comparisons_eff$denominator,
                                          genotype = genotype2,
                                          sex = "F"
  )

  pred_numerator_1 <- posterior_epred(fit, newdata = pred_data_numerator_1, re_formula = NA, cores = 1)
  pred_numerator_2 <- posterior_epred(fit, newdata = pred_data_numerator_2, re_formula = NA, cores = 1)
  pred_denominator_1 <- posterior_epred(fit, newdata = pred_data_denominator_1, re_formula = NA, cores = 1)
  pred_denominator_2 <- posterior_epred(fit, newdata = pred_data_denominator_2, re_formula = NA, cores = 1)

  pred_numerator <- pred_numerator_1 - pred_numerator_2
  pred_denominator <- pred_denominator_1 - pred_denominator_2

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

compute_ratios_observed <- function(data_model, comparisons) {
  replicate_averaged <- data_model %>%
    group_by(run, animal_no, sex, genotype, primer_short) %>%
    summarise(avg_Cq = mean(Cq), .groups = "drop")


  ratios_observed <- comparisons %>%
    inner_join(replicate_averaged, by = c("numerator" = "primer_short")) %>%
    inner_join(replicate_averaged, by = c("denominator" = "primer_short", "run", "animal_no", "sex", "genotype"), suffix = c("_num", "_denom")) %>%
    mutate(Cq_denom_num = avg_Cq_denom - avg_Cq_num)

  ratios_observed
}

compute_comparisons_obs <- function(data_model, comparisons_eff, genotype_1, genotype_2) {
  averaged <- data_model %>%
    group_by(genotype, primer_short) %>%
    summarise(avg_Cq = mean(Cq), sd_Cq = sd(Cq), sem_Cq = sd_Cq / n(), .groups = "drop")

  averaged_comparisons <- comparisons_eff %>%
    inner_join(averaged, by = c("numerator" = "primer_short")) %>%
    rename(numerator_avg_Cq = avg_Cq, numerator_sem_Cq = sem_Cq) %>%
    inner_join(averaged, by = c("denominator" = "primer_short", "genotype")) %>%
    rename(denominator_avg_Cq = avg_Cq, denominator_sem_Cq = sem_Cq) #%>%

  comparisons_obs <- averaged_comparisons %>%
    inner_join(averaged_comparisons, by = c("label", "comparison_id", "eff_numerator", "eff_denominator", "eff_label"), suffix = c("_2","_1")) %>%
    filter(genotype_2 == !!genotype_2, genotype_1 == !!genotype_1) %>%
    mutate(log_ratio = log(eff_numerator) * (numerator_avg_Cq_1 - numerator_avg_Cq_2) -  log(eff_denominator) * (denominator_avg_Cq_1 - denominator_avg_Cq_2),
           sem_ratio = log(eff_numerator) * (numerator_sem_Cq_1 + numerator_sem_Cq_2) -  log(eff_denominator) * (denominator_sem_Cq_1 + denominator_sem_Cq_2))

  comparisons_obs
}


summarise_comparisons <- function(comparisons_pred) {
  res <- comparisons_pred
  res <- dplyr::group_by(res, label_variant, genotypes, eff_label)
  res <- dplyr::summarise(res, low = quantile(log_ratio, prob = 0.025),
                          high = quantile(log_ratio, prob = 0.975),
                          sign = case_when(low > 0 ~ "+",
                                           high < 0 ~ "-",
                                           TRUE ~ "0"),
                          .groups = "drop_last")

  res <- dplyr::summarise(res,
                          robust = length(unique(sign)) == 1,
                          low = low[eff_label == "Neutral"],
                          high = high[eff_label == "Neutral"],
                          sign = sign[eff_label == "Neutral"]
                          )

  res
}
