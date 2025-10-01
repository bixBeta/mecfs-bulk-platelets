# extract coldata and ratio-counts ----

dds_setup = lapply(ratios.by.enid.list, function(x) {
  map(c("target", "aux_filtered_ratios_matrix"), ~ pluck(x, .x)) |>
    set_names(c("target", "aux_filtered_ratios_matrix"))
})

# set up dds ----

pca_list = lapply(dds_setup, function(x) {
  getPCA(dds_setup_list = x)
})


ggPcas = lapply(pca_list, function(x) {
  pluck(x, "plot")
})

cowplot::plot_grid(plotlist = ggPcas, nrow = 1)
