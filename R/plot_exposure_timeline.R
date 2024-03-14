## This file is part of the rRACES (https://github.com/caravagnalab/rRACES/).
## Copyright (C) 2023 - Giulio Caravagna <gcaravagna@units.it>
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <https://www.gnu.org/licenses/>.

preprocess <- function(x) {
  time_points <- unique(x$time)
  a <- apply(x, 1, function(x) time_points[which(time_points == as.numeric(x[1])) + 1] )
  a[is.na(a)] <- max(time_points) #+ mean(diff(time_points))
  x$endtime <- a
  x <- x[, c("time", "signature", "exposure", "type", "endtime")]
  return(x)
}

#' Plot a signatures exposures
#'
#' @description
#' Plots a signatures exposures in different time stamps
#'
#' @param x exposure dataframe
#'
#' @return An editable ggplot plot.
#' @export
plot_exposure_timeline <- function(x) {
  
  signames <- x %>%dplyr::pull(signature) %>% unique
  all_colors <- rRACES:::get_signatures_colors()
  colors <- all_colors[signames]
  
  x_prep <- preprocess(x)
  
  p <- ggplot2::ggplot(data = x_prep, ggplot2::aes(x = time, y = exposure, group = signature, color = signature)) + 
    ggplot2::geom_segment(ggplot2::aes(xend = endtime, yend = exposure), linewidth = 0.5, position = ggplot2::position_jitter(width = 0.04, height = 0.04, seed = 10)) + 
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.04, height = 0.04, seed = 10)) + 
    ggplot2::ylim(0, 1) + 
    ggplot2::labs(x = "Time Point", y = "Signature Exposure", title = "Mutational Signature Exposure Timeline") + 
    #scale_fill_manual(values=colors) + 
    ggplot2::scale_colour_manual(name = "sigs",values = colors) + 
    ggplot2::theme_minimal()  # Apply a minimal theme
  
  #ggplot(data = x, aes(x = time, y = exposure, group = signature, color = signature)) +
  #  geom_step(direction = "hv", linewidth = 1, position = position_jitter(width = 0.05, height = 0.05, seed = 10)) +  # Connect points with horizontal lines
  #  #geom_point(size = 3, position = position_jitter(width = 0.00, height = 0.004, seed = 10)) +  # Add points
  #  ylim(0, 1) + 
  #  labs(x = "Time Point", y = "Signature Exposure", title = "Mutational Signature Exposure Timeline") +
  #  theme_minimal()  # Apply a minimal theme
  
  return(p)
}

