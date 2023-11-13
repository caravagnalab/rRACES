devtools::load_all()

# Create a new simulation and set the simulation data directory
# to be "output_dir". A 1000x1000-cells tissue is automatically
# added to the simulation
x <- create_simulation(
  name = "My Test rRACES",
  output_dir = "some_folder",
  tissue = "A fake liver",
  tissue_size = c(2000, 2000)
)

# Add two genotype, "A" and "B", having epigenetic states
x <- add_genotype(
  x,
  name = "A",
  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
  growth_rates = c("+" = 0.2, "-" = 0.08),
  death_rates = c("+" = 0.1, "-" = 0.01)
)

x <- add_genotype(
  x,
  name = "B",
  epigenetic_rates = c("+-" = 0.01, "-+" = 0.01),
  growth_rates = c("+" = 0.3, "-" = 0.3),
  death_rates = c("+" = 0.02, "-" = 0.02)
)

x <- add_genotype(
  x,
  name = "C",
  epigenetic_rates = c("+-" = 0.001, "-+" = 0.01),
  growth_rates = c("+" = 0.2, "-" = 0.2),
  death_rates = c("+" = 0.0, "-" = 0.0)
)

x <- add_genotype(
  x,
  name = "D",
  epigenetic_rates = c("+-" = 0.001, "-+" = 0.01),
  growth_rates = c("+" = 0.2, "-" = 0.2),
  death_rates = c("+" = 0.0, "-" = 0.0)
)


# Program a timed transition from species "A" to
# species "B" at time 60
x <- add_genotype_evolution(x, "A", "B", 20)
x <- add_genotype_evolution(x, "B", "C", 35)
x <- add_genotype_evolution(x, "B", "D", 45)

# Add one cell to the simulated tissue
x <- set_initial_cell(
  x,
  genotype = "A",
  epistate = "+",
  position = c(1000, 1000)
)

time_series <- lapply(seq(20, 80, by = 1),
                      function(t) {
                        run(x, what = list(time = t))
                        x$simulation$get_counts() %>%
                          mutate(
                            time = x$simulation$get_clock(),
                            species = paste0(genotype, epistate)
                          )
                      }) %>%
  Reduce(f = bind_rows) %>%
  as_tibble()

plot_state(x)
plot_tissue(x)

plot_tissue(x) + facet_grid(genotype ~ epistate)

time_series %>%
  # dplyr::filter(species == 'B+') %>%
  ggplot() +
  my_theme() +
  geom_point(aes(x = time, y = counts, color = genotype, alpha = epistate)) +
  geom_line(aes(x = time, y = counts, color = genotype, alpha = epistate)) +
  ggplot2::scale_color_manual(values =
                                get_species_colors(time_series$genotype %>%
                                                     unique)) +
  ggplot2::scale_alpha_manual(values = c(`+` = 1, `-` = .5))  +
  ggplot2::theme(legend.position = "bottom")


library(ggmuller)

pop_df <- time_series %>%
  rename(Generation = .data$time, Identity = .data$species,
         Population = .data$counts) %>%
  select(Generation, Identity, Population)

pop_df <- pop_df %>% bind_rows(data.frame(Generation = 0,
                                          Identity = "A+", Population = 1))
pop_df <- pop_df %>% bind_rows(data.frame(Generation = 0,
                                          Identity = "WT", Population = 1))

edges <- data.frame(
  Parent = c("A+", "A+", "B+"),
  Identity = c("A-", "B+", "B-")
)

edges <- data.frame(
  Parent = c("WT", "A+", "A+", "B+", "B+", "C+", "B+", "D+"),
  Identity = c("A+", "A-", "B+", "B-", "C+", "C-", "D+", "D-")
)
muller_df <- get_Muller_df(edges, pop_df)

sp_colors <- get_species_colors(c("A", "B", "C", "D"))

sp_colors_p <- sp_colors
sp_colors_m <- alpha(sp_colors, 0.5)

names(sp_colors_p) <- paste0(names(sp_colors_p), "+")
names(sp_colors_m) <- paste0(names(sp_colors_m), "-")

Muller_plot(muller_df,
            add_legend = TRUE,
            xlab = "Time",
            ylab = "Proportion",
            palette = c(`WT` = "gainsboro", sp_colors_p, sp_colors_m)) +
  my_theme()


muller_df_pop <- add_empty_pop(Muller_df)

# list of legend entries, omitting NA
id_list <- sort(unique(Muller_df_pop$Identity))

ggplot(Muller_dmuller_df_popf_pop,
       aes_string(x = "Generation",
                  y = "Population",
                  group = "Group_id",
                  fill = "Identity",
                  colour = "Identity")) +
  geom_area() +
  theme(legend.position = "right") +
  guides(linetype = FALSE, color = FALSE) + 
  scale_fill_manual(name = "Identity", values = c(`WT` = "gainsboro",
                                                  sp_colors_p, sp_colors_m),
                    breaks = id_list) +
  scale_color_manual(values = c(`WT` = "gainsboro", sp_colors_p, sp_colors_m)) +
  my_theme()


