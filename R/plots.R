genes_boxplot <- function(data, primers, step_size = 8) {
  for(step in 1:ceiling(length(primers) / step_size)) {
    primers_to_show <- primers[((step - 1) * step_size + 1):(step * step_size)]
    # print(d17_4w %>% filter(primer %in% primers_to_show) %>%
    #   ggplot(aes(x = ct)) + geom_histogram(binwidth = 0.2) + facet_grid(genotype ~ primer, scales = "free_y")
    # )
    print(data %>% filter(primer_short %in% primers_to_show) %>%
            ggplot(aes(y = Cq, x = genotype)) + geom_boxplot(outlier.shape = NA) +
            geom_point(aes(color = sex, shape = sex), position = position_jitter(width = 0.25), alpha = 0.7) +
            facet_wrap( ~ primer_short, scales = "free", nrow = 2)
    )

  }
}

genes_parcoord <- function(data, primers) {
  data %>% filter(primer_short %in% primers) %>%
    group_by(animal_plate, genotype, primer_short, sex) %>%
    summarise(Cq = mean(Cq), .groups = "drop") %>%
    ggplot(aes(y = Cq, x = primer_short, color = sex, group = animal_plate)) +
    geom_line(alpha = 0.7) +
    facet_wrap( ~ genotype, scales = "free")
}

comparison_label_theme <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.3), axis.title.x = element_blank())

scale_color_type <- scale_color_brewer(type = "qual")

plot_detailed_comparisons <- function(comparisons, comparisons_pred, ratios_observed, n_to_show = 4) {

  all_comparisons <- sort(unique(comparisons$label))

  comparison_groups <- list()
  comparison_groups[[1]] <- c()
  for(comp in all_comparisons) {
    variants <- comparisons %>% filter(label == comp) %>% pull(label_variant)
    n_variants <- length(variants)
    last <- length(comparison_groups)
    if(last == 0 || length(comparison_groups[[last]]) + n_variants > n_to_show) {
      comparison_groups[[last + 1]] <- variants
    } else {
      comparison_groups[[last]] <- c(comparison_groups[[last]], variants)
    }
  }



  for(i in 1:length(comparison_groups)) {
    comparisons_to_show <- comparison_groups[[i]]

    nudge_width <- 0.2

    plot_pred <- comparisons_pred %>%
      filter(label_variant %in% comparisons_to_show, eff_label == "Neutral") %>%
      ggplot(aes(x = label_variant, y = exp(log_ratio), color = genotypes, shape = genotypes)) +
      geom_hline(yintercept = 1, color = "blue") +
      expand_limits(y = c(0.5, 2)) +
      stat_pointinterval(position = position_dodge(width = 0.5)) +
      scale_color_type +
      scale_y_log10("ratio of ratios") + comparison_label_theme

    plot_observed <- ratios_observed %>%
      filter(label_variant %in% comparisons_to_show) %>%
      ggplot(aes(x = genotype, y = Cq_denom_num)) + geom_boxplot(outlier.shape = NA) + geom_point(aes(color = sex, shape = sex), position = position_jitter(width = 0.3)) + facet_wrap(~label_variant, scales = "free_y", ncol = 2) +
      scale_y_continuous("Cq (denominator) - Cq (numerator)")

    print((plot_observed | plot_pred) + plot_layout(widths = c(2,1)))
  }
}

plot_all_comparisons <- function(comparisons_pred, expand = c(0.5,2)) {
  comparisons_pred %>%
    filter(eff_label == "Neutral") %>%
    ggplot(aes(x = label_variant, y = exp(log_ratio), color = genotypes, shape = genotypes)) +
    geom_hline(yintercept = 1, color = "blue") +
    stat_pointinterval(position = position_dodge(width = 0.4)) +
    expand_limits(y = expand) +
    scale_y_log10("ratio of ratios") + scale_color_type +
    comparison_label_theme
}


plot_comparisons_sensitivity <- function(comparisons_pred, ncol = 1, expand = c(0.5,2)) {
  comparisons_pred %>%
    ggplot(aes(x = label_variant, y = exp(log_ratio), color = eff_label)) +
    geom_hline(yintercept = 1, color = "blue") +
    stat_pointinterval(position = position_dodge(width = 0.5)) +
    expand_limits(y = expand) +
    scale_y_log10("ratio of ratios") + scale_color_type +
    facet_wrap(~genotypes, ncol = ncol, scales = "free_x") +
    comparison_label_theme
}
