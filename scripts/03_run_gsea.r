key <- readRDS(
  "~/Desktop/4052_plusMore/Platelet-4052-RNA-seq-Paired-Analysis/data/key.RDS"
)
rank_list = lapply(ratios.by.enid.list, function(x) {
  annotateRanks(ratio_ = x)
})


database_list = list()

database = c("H", "C2", "C5")

# one-time-run ----
# for (i in database) {
#   database_list[[i]] <- msigdbr(species = "Homo sapiens", category = i) %>%
#     dplyr::select(gs_name, gene_symbol, gs_cat, gs_subcat, gs_description)
# }
#
# c2 = database_list$C2
# c2.noncgp = c2 %>% filter(., gs_subcat != "CGP")
# database_list$C2 <- c2.noncgp
#
# saveRDS(database_list, "data/msigdb.RDS")

database_list <- readRDS("data/msigdb.RDS")

gse = list()

for (db in database) {
  message(paste0("Processing ----", db))
  gse[[db]] <- lapply(rank_list, FUN = function(x) {
    GSEA(
      x,
      TERM2GENE = database_list[[db]],
      pAdjustMethod = "BH",
      pvalueCutoff = 1
    )
  })
}

saveRDS(gse, "data/gsea_results_droppedEnids.RDS")

test = gse$C5$D2.PRE_vs_D1.PRE

enrichplot::gseaplot2(
  test,
  geneSetID = 5,
  title = test$Description[5]
)
GSEAPlot(
  test,
  gene_sets = attr(test, "geneSets")["GOBP_CYTOPLASMIC_TRANSLATION"],
  gene_ranks = test@geneList
)

test = gseaRes(gse_ = gse)
