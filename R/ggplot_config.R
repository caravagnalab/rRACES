my_theme = function(){
  ggplot2::theme_linedraw(base_size = 10) +
    ggplot2::theme(
      legend.position = 'bottom'
    )
}

get_species_colors = function(x){
  nm = sort(x)
  
  nmc = RColorBrewer::brewer.pal(x %>% length(), "Set1")
  
  if(length(nm) < length(nmc)) 
    nmc = nmc[1:length(nm)]
  
  names(nmc) = nm
  
  return(nmc)
}
