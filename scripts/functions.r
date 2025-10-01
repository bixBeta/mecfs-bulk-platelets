#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####
# Quantile Filtering on norm counts + density plots ----
####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

extract.filter.normCounts = function(dds) {
  # extract coldata
  meta.data = colData(dds) |> data.frame()

  # add enid.id column for grepping later

  # meta.data$orig.ident = paste0("X", meta.data$label)
  # meta.data$orig.ident = gsub(pattern = "-", replacement = ".", x = meta.data$orig.ident)
  colnames(meta.data)[which(colnames(meta.data) == 'label')] <- 'orig.ident'
  meta.data$enid.id = paste0(meta.data$orig.ident, "_", meta.data$enid)
  # extract normalized counts
  norm.counts = counts(dds, normalized = T) |> data.frame()

  # add row medians to norm.counts
  norm.counts$row.medians = apply(X = norm.counts, MARGIN = 1, FUN = median)

  # get the quantile profile of the median
  q.tile = quantile(norm.counts$row.medians)

  # bins.quantiles(norm.counts$row.medians, target.bins = 4, max.breaks = 15)$binct
  q.tile[4] <- 50
  q75.filtered.norm.counts = norm.counts[norm.counts$row.medians > q.tile[4], ]

  # plot norm counts pre and post filter

  counts.df.no.filter = norm.counts |> select(-row.medians)
  counts.df.stacked.no.filter = stack(counts.df.no.filter)
  gg.df.no.filter = left_join(
    counts.df.stacked.no.filter,
    meta.data,
    by = c("ind" = "orig.ident")
  )

  p.pre.filter <- ggplot(gg.df.no.filter, aes(x = .data$values + 1)) +
    stat_density(
      aes(group = .data$ind, color = .data$day),
      position = "identity",
      geom = "line",
      show.legend = TRUE
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(~ 10^.x))
    ) +
    labs(color = "") +
    xlab(paste0("normCounts_", deparse(substitute(dds)))) +
    ylab("Density") +
    ggtitle("Density of counts distribution") +
    theme_gray() +
    facet_wrap("enid")

  counts.df.post.filter = q75.filtered.norm.counts |> select(-row.medians)
  counts.df.stacked.post.filter = stack(counts.df.post.filter)
  gg.df.post.filter = left_join(
    counts.df.stacked.post.filter,
    meta.data,
    by = c("ind" = "orig.ident")
  )

  p.post.filter <- ggplot(gg.df.post.filter, aes(x = .data$values + 1)) +
    stat_density(
      aes(group = .data$ind, color = .data$day),
      position = "identity",
      geom = "line",
      show.legend = TRUE
    ) +
    scale_x_continuous(
      trans = log10_trans(),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(~ 10^.x))
    ) +
    labs(color = "") +
    xlab(paste0(
      "normCounts_",
      deparse(substitute(dds)),
      "__",
      "Row_median >",
      round(q.tile[4], 2)
    )) +
    ylab("Density") +
    ggtitle("Density of counts distribution") +
    theme_gray() +
    facet_wrap("ind")

  #return(p.pre.filter)
  return(list(
    colData = meta.data,
    norm.counts.No.filter = norm.counts,
    Quantile = q.tile,
    q75.filtered.norm.counts = q75.filtered.norm.counts,
    ggplots = list(
      ggplot.pre.filter = p.pre.filter,
      ggplot.post.filter = p.post.filter
    )
  ))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####
# Get Ratios by ENID ----
####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

extractRatios = function(xl, num_, denom_) {
  normCounts_ = as.matrix(
    xl$q75.filtered.norm.counts |> dplyr::select(-row.medians)
  )
  meta_ = xl$colData

  stack_normCounts_ = as_tibble(stack(normCounts_))
  genes_ = levels(stack_normCounts_$row)
  enid.ids_ = levels(stack_normCounts_$col)

  stack_normCounts_[stack_normCounts_ == 0] <- NA

  levels.enid = list()
  for (i in 1:length(enid.ids_)) {
    levels.enid[[i]] <- stack_normCounts_ |> filter(grepl(enid.ids_[i], col))
    names(levels.enid)[[i]] <- enid.ids_[i]
  }

  enids = meta_ |>
    data.frame() |>
    pull(enid) |>
    unique()

  temp_ = list()
  for (i in 1:length(enids)) {
    t = levels.enid[grep(enids[i], names(levels.enid))]
    temp_[[i]] <- do.call(cbind, t)
    names(temp_)[[i]] <- enids[i]
  }

  dropNA = lapply(temp_, function(x) {
    x = x |> select(1, matches(".value"))
    x = x |> column_to_rownames(colnames(x)[1])
    x = x |> drop_na()
    x = x |> rownames_to_column("gene")
  })

  dropNA.joined = dropNA |> reduce(full_join, by = "gene")
  dropNA.joined = dropNA.joined |> column_to_rownames("gene")

  # # test.joined$gm = apply(X = test.joined, MARGIN = 1, FUN = geometric.mean)
  message("Flooring")

  # flooring
  dropNA.joined[dropNA.joined >= 0 & dropNA.joined < 1] <- 1

  denom.counts = dropNA.joined |>
    rownames_to_column("gene") |>
    select(matches(denom_), "gene") |>
    column_to_rownames("gene")
  num.counts = dropNA.joined |>
    rownames_to_column("gene") |>
    select(matches(num_), "gene") |>
    column_to_rownames("gene")

  ratios = list()

  for (i in 1:length(enids)) {
    ratios[[i]] <- as.data.frame(
      num.counts[, grep(enids[i], colnames(num.counts))] /
        denom.counts[, grep(enids[i], colnames(denom.counts))]
    )

    names(ratios)[[i]] <- enids[i]

    rownames(ratios[[i]]) <- rownames(num.counts)

    # ifelse(nrow(ratios[[i]]) > 1, colnames(ratios[[i]]) <- enids[i], NA)

    colnames(ratios[[i]]) <- enids[i]
  }

  ratios.matrix = as.matrix(do.call(cbind, ratios))

  # separate numerator and denominator enids

  cfs_enids = meta_ |>
    filter(group == "CFS") |>
    filter(enid %in% enids) |>
    pull(enid) |>
    unique()
  ctrl_enids = meta_ |>
    filter(group == "control") |>
    filter(enid %in% enids) |>
    pull(enid) |>
    unique()

  ratios.matrix.aux = ratios.matrix |>
    data.frame() |>
    rownames_to_column("gene") |>
    rowwise() |>
    mutate(
      n_cfs_samples_with_non_zero_value = sum(
        c_across(cols = matches(cfs_enids)) > 0,
        na.rm = T
      ),

      n_ctrl_samples_with_non_zero_value = sum(
        c_across(cols = matches(ctrl_enids)) > 0,
        na.rm = T
      )
    ) |>
    ungroup() |>
    filter(
      n_cfs_samples_with_non_zero_value > 2 &
        n_ctrl_samples_with_non_zero_value > 2
    ) |>
    column_to_rownames("gene")

  cfs_ratios = ratios.matrix.aux |> data.frame() |> select(matches(cfs_enids))
  cfs_ratios$geomMean.cfs = apply(
    X = cfs_ratios,
    MARGIN = 1,
    FUN = geometric.mean
  )

  ctrl_ratios = ratios.matrix.aux |> data.frame() |> select(matches(ctrl_enids))
  ctrl_ratios$geomMean.ctrl = apply(
    X = ctrl_ratios,
    MARGIN = 1,
    FUN = geometric.mean
  )

  log2FC_geomMean = log2(cfs_ratios$geomMean.cfs) -
    log2(ctrl_ratios$geomMean.ctrl)
  names(log2FC_geomMean) <- rownames(cfs_ratios)

  log2FC_geomMean = sort(log2FC_geomMean, decreasing = T)

  message("DONEEEE x) ")

  return(list(
    ratios_matrix = ratios.matrix,
    aux_filtered_ratios_matrix = ratios.matrix.aux,
    num.counts = num.counts,
    denom.counts = denom.counts,
    target = meta_,
    cfs_ratios = cfs_ratios,
    ctrl_ratios = ctrl_ratios,
    log2FC_geomMean = log2FC_geomMean
  ))
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####
# dds setup for ratios ----
####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

getPCA = function(dds_setup_list) {
  target = dds_setup_list[["target"]]
  target = target |> select(enid, group) |> unique()
  rownames(target) <- target$enid

  counts = dds_setup_list[["aux_filtered_ratios_matrix"]]
  counts = counts |> select(all_of(target$enid))
  counts = drop_na(counts)

  rv <- rowVars(as.matrix(counts))
  select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
  pca = prcomp(t(as.matrix(counts)[select, ]))
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  pVar.df <- as.data.frame(percentVar)
  pVar.df$x = as.factor(paste0("PC", rownames(pVar.df)))

  pVar.df = pVar.df[, order(names(pVar.df))]
  pVar.df$percentVar = pVar.df$percentVar * 100
  pVar.df$percentVar = round(pVar.df$percentVar, digits = 2)

  d <- data.frame(pca$x, label = rownames(pca$x))
  d2 <- left_join(d, target, by = c("label" = "enid"))

  suppressPackageStartupMessages(library(ggrepel))
  suppressPackageStartupMessages(library(ggplot2))

  pc1 = ggplot(d2, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 5, aes(shape = group)) +
    geom_label_repel(
      aes(label = label),
      box.padding = 0.8,
      point.padding = 0.5,
      segment.color = 'grey55',
      show.legend = F
    ) +
    xlab(paste0(pVar.df$x[1], "  ", pVar.df$percentVar[1], "%")) +
    ylab(paste0(pVar.df$x[2], "  ", pVar.df$percentVar[2], "%")) +
    theme_light() +
    theme(legend.position = "bottom") +
    scale_color_manual(values = colors2)
  # geom_mark_hull(aes(fill = batch),
  #                concavity = 20, alpha = 0.08
  #)

  # png(paste0(pin, "_PC1_PC2.png"), width = 1200, height = 1200, res = 150)
  # pc1
  # dev.off()

  return(list(
    prcomp.out = pca,
    Variance.df = pVar.df,
    colData = target,
    PCA.df = d2,
    plot = pc1
  ))
}

#
colors2 <- c(
  "red",
  "#1f78b4",
  "#1b9e77",
  "purple3",
  "khaki4",
  "#E9A3C9",
  "#A1D76A",
  "#EF8A62",
  "grey"
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####
# getRanklist - annotated  ----
####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

annotateRanks = function(ratio_) {
  l = ratio_[["log2FC_geomMean"]]
  l = l |> data.frame() |> rownames_to_column("id")
  l = left_join(l, key, by = c("id" = "gene"))
  l = l |> distinct(symbol, .keep_all = T)
  l = l |> drop_na()
  l2 = l |> pull(l)
  names(l2) <- l$symbol
  return(l2)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####
# extractGSEA ----
####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

gseaRes = function(gse_) {
  ul = unlist(gse_)
  ul_res = lapply(ul, function(x) {
    pluck(x, "result")
  })

  ul_res2 = do.call(rbind, ul_res)

  return(ul_res2)
}
