my_theme = function(){
  ggplot2::theme_linedraw(base_size = 10) +
    ggplot2::theme(
      legend.position = 'bottom'
    )
}

plot_VAFs = function(x)
{
  nchr = x$chr %>% unique() %>% length()
  nmut = x %>% nrow

   x %>%
    ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(VAF, fill = type), binwidth = 0.01) +
    ggplot2::guides(fill = ggplot2::guide_legend("")) +
    ggplot2::xlim(0,1) +
    my_theme() +
    ggplot2::labs(
      caption = paste(nchr, "chromosomes", nmut, "mutations"),
      y = "Counts",
      title = 'Variant Allele Frequencies'
    ) +
     ggplot2::scale_fill_manual(values = c(`SNV` = "steelblue", `indel` = 'indianred3'))
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
