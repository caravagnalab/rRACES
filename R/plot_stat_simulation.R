

plot_VAFs <- function(x) {
  nchr <- x$chr %>% unique() %>% length()
  nmut <- x %>% nrow

  x %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(.data$VAF, fill = .data$cause),
                            binwidth = 0.01) +
    ggplot2::guides(fill = ggplot2::guide_legend("")) +
    ggplot2::xlim(0, 1) +
    my_theme() +
    ggplot2::labs(
      caption = paste(nchr, "chromosomes", nmut, "mutations"),
      y = "Counts",
      title = "Variant Allele Frequencies"
    ) +
    ggplot2::facet_wrap(~.data$type)
  # ggplot2::scale_fill_manual(values = c(`SNV` = "steelblue",
  #                                       `indel` = "indianred3"))
}

plot_depth = function(x) {
  nchr <- x$chr %>% unique() %>% length()
  nmut <- x %>% nrow

  x %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(.data$coverage,
                                         fill = .data$type), bins = 30) +
    my_theme() +
    ggplot2::labs(
      caption = paste(nchr, "chromosomes", nmut, "mutations"),
      y = "Counts",
      title = "Sequencing coverage"
    ) +
    ggplot2::scale_fill_manual(values = c(`SNV` = "steelblue", 
                                          `indel` = "indianred3"))
}

plot_signatures <- function(x) {
  x_signatures <- x %>%
    dplyr::filter(!.data$is_driver) %>%
    dplyr::group_by(.data$context, .data$ref, .data$alt, .data$cause) %>%
    dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(
      substitution = paste0(.data$ref, ">", .data$alt),
      new_substitution = case_when(
        (.data$ref == "G" & .data$alt == "T") ~ "C>A",
        (.data$ref == "G" & .data$alt == "C") ~ "C>G",
        (.data$ref == "G" & .data$alt == "A") ~ "C>T",
        (.data$ref == "A" & .data$alt == "T") ~ "T>A",
        (.data$ref == "A" & .data$alt == "G") ~ "T>C",
        (.data$ref == "A" & .data$alt == "C") ~ "T>G"
      ),
      substitution = ifelse(
        is.na(.data$new_substitution),
        .data$substitution,
        .data$new_substitution
      )
    )

  substr(x_signatures$context, 2, 2) <- "_"

  nchr <- x$chr %>% unique() %>% length()
  nmut <- x %>% nrow

  x_signatures %>%
    ggplot2::ggplot() +
    ggplot2::geom_bar(ggplot2::aes(x = .data$context, y = .data$n,
                                   fill = .data$cause), stat = "identity") +
    my_theme() +
    ggplot2::labs(
      caption = paste(nchr, "chromosomes", nmut, "mutations"),
      y = "Counts",
      title = "Signatures"
    ) +
    ggplot2::facet_grid(.data$cause~.data$substitution, scales = "free_y") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
}