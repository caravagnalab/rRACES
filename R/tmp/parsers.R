parse_VAF_file = function(f = 'nobuild/chr_1_VAF.csv')
{
  x = readr::read_tsv(f, col_types = readr::cols())

  if(!any(grepl('chr', x$chr))) x$chr = paste0('chr', x$chr)

  x %>% sanitiser()
}

sanitiser = function(x)
{
  x %>%
    dplyr::mutate(
      status =
        case_when(
          grepl("SBS", cause) ~ "SBS",
          grepl("DBS", cause) ~ "DBS",
          grepl("IS", cause) ~ "IS",
          TRUE ~ "driver"
        ),
      is_driver = status == 'driver'
    )
}
