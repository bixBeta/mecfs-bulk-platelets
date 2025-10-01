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


png(
  filename = "figures/ratios_pca_dropped_enids.png",
  width = 2920,
  height = 2080,
  res = 150
)
cowplot::plot_grid(
  plotlist = ggPcas,
  nrow = 2,
  labels = c(names(ggPcas)),
  label_size = 9,
  hjust = -0.1,
  vjust = 70,
  rel_widths = 1080,
  rel_heights = 1080
)
dev.off()
