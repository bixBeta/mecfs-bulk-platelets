source("scripts/packages.r")
source("scripts/functions.r")
load("../../Platelet-4052-RNA-seq-Paired-Analysis/data/4052.dds_SARtools.RData")

dds = out.DESeq2$dds

raw.counts = counts(dds, normalized = F)
target = target |> filter(enid != "E759" & enid != "E576")
target
raw.counts = raw.counts |> data.frame() |> select(all_of(target$label))
table(target$label == colnames(raw.counts))

dds = DESeqDataSetFromMatrix(
  countData = raw.counts,
  colData = target,
  design = ~group
)
dds = DESeq(dds)
vst = varianceStabilizingTransformation(dds, blind = T)
plotPCA(vst, intgroup = "group")

xl = extract.filter.normCounts(dds = dds)
saveRDS(xl, "data/xl_droppedEnids.RDS")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

xl <- readRDS("data/xl_droppedEnids.RDS")

combn.df = read.csv("data/combn.csv", header = T)


ratios.by.enid.list = list()

for (i in 1:nrow(combn.df)) {
  ratios.by.enid.list[[i]] <- extractRatios(
    xl = xl,
    num_ = combn.df$num[i],
    denom_ = combn.df$denom[i]
  )

  names(ratios.by.enid.list)[[i]] <- paste0(
    combn.df$num[i],
    "_vs_",
    combn.df$denom[i]
  )
}

saveRDS(ratios.by.enid.list, "data/ratios_by_enid_list_DroppedENID.RDS")
