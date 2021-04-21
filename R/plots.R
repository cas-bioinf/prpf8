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


