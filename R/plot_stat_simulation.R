

plot_VAFs = function(x)
{
  nchr = x$chr %>% unique() %>% length()
  nmut = x %>% nrow
  
  x %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(VAF, fill = cause), binwidth = 0.01) +
    ggplot2::guides(fill = ggplot2::guide_legend("")) +
    ggplot2::xlim(0,1) +
    my_theme() +
    ggplot2::labs(
      caption = paste(nchr, "chromosomes", nmut, "mutations"),
      y = "Counts",
      title = 'Variant Allele Frequencies'
    ) +
    ggplot2::facet_wrap(~type)
  # ggplot2::scale_fill_manual(values = c(`SNV` = "steelblue", `indel` = 'indianred3'))
}

plot_depth = function(x)
{
  nchr = x$chr %>% unique() %>% length()
  nmut = x %>% nrow
  
  x %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(coverage,fill = type), bins = 30) +
    my_theme() +
    ggplot2::labs(
      caption = paste(nchr, "chromosomes", nmut, "mutations"),
      y = "Counts",
      title = 'Sequencing coverage'
    ) +
    ggplot2::scale_fill_manual(values = c(`SNV` = "steelblue", `indel` = 'indianred3'))
}

plot_signatures = function(x)
{
  x_signatures = x %>%
    dplyr::filter(!is_driver) %>%
    dplyr::group_by(context, ref, alt, cause) %>%
    dplyr::summarise(n = dplyr::n(), .groups = 'drop') %>%
    dplyr::mutate(
      substitution = paste0(ref, '>', alt),
      new_substitution = case_when(
        (ref == "G" & alt == 'T') ~ "C>A",
        (ref == "G" & alt == 'C') ~ "C>G",
        (ref == "G" & alt == 'A') ~ "C>T",
        (ref == "A" & alt == 'T') ~ "T>A",
        (ref == "A" & alt == 'G') ~ "T>C",
        (ref == "A" & alt == 'C') ~ "T>G"
      ),
      substitution = ifelse(
        is.na(new_substitution),
        substitution,
        new_substitution
      )
    )
  
  substr(x_signatures$context, 2, 2) = '_'
  
  nchr = x$chr %>% unique() %>% length()
  nmut = x %>% nrow
  
  x_signatures %>%
    ggplot2::ggplot() +
    ggplot2::geom_bar(ggplot2::aes(x = context, y = n, fill = cause), stat = 'identity') +
    my_theme() +
    ggplot2::labs(
      caption = paste(nchr, "chromosomes", nmut, "mutations"),
      y = "Counts",
      title = 'Signatures'
    ) +
    ggplot2::facet_grid(cause~substitution, scales = 'free_y') +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
}