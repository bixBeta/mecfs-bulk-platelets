source("scripts/packages.r")
source("scripts/functions.r")
load("../Platelet-4052-RNA-seq-Paired-Analysis/data/4052.dds_SARtools.RData")

dds = out.DESeq2$dds

raw.counts = counts(dds, normalized = F)
norm.counts = counts(dds, normalized = T)

xl = extract.filter.normCounts(dds = dds)
saveRDS(xl, "data/xl.RDS")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

xl <- readRDS("data/xl.RDS")
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

saveRDS(ratios.by.enid.list, "data/ratios_by_enid_list.RDS")
