test = function(xl, num_, denom_) {
  # get median filtered counts

  normCounts_ = as.matrix(
    xl$q75.filtered.norm.counts |> dplyr::select(-row.medians)
  )

  # get target meta

  meta_ = xl$colData |> filter(group2 == num_ | group2 == denom_)

  # table normalize to 1-1

  stack_normCounts_ = as_tibble(stack(normCounts_))
  genes_ = levels(stack_normCounts_$row)
  enid.ids_ = levels(stack_normCounts_$col)

  # swap zero counts to NA

  stack_normCounts_[stack_normCounts_ == 0] <- NA

  # get 1-1 per enid level

  levels.enid = list()
  for (i in 1:length(enid.ids_)) {
    levels.enid[[i]] <- stack_normCounts_ |> filter(grepl(enid.ids_[i], col))
    names(levels.enid)[[i]] <- enid.ids_[i]
  }

  # get enids of interest

  enids = meta_ |>
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

  message("Flooring")

  # flooring (replace counts between 0 and 1 to 1 )

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

    # names(ratios)[[i]] <- enids[i]
    # rownames(ratios[[i]]) <- rownames(num.counts)
    # colnames(ratios[[i]]) <- enids[i]
  }

  # Get Matrix of the num/denom ratios

  ratios.matrix = as.matrix(do.call(cbind, ratios))

  case.enids = meta_ |> filter(group == "CFS") |> pull(enid) |> unique()
  case.enids = case.enids[case.enids %in% enids]

  ctrl.enids = meta_ |> filter(group == "control") |> pull(enid) |> unique()
  ctrl.enids = ctrl.enids[ctrl.enids %in% enids]

  ratios.matrix.aux = ratios.matrix |>
    data.frame() |>
    rownames_to_column("gene") |>
    rowwise() |>
    mutate(
      #First condition, at least 2 exp columns > 0
      aux_condition_g2 = sum(
        c_across(cols = matches(case.enids)) > 0,
        na.rm = T
      ),
      #First condition, at least 2 stat columns > 0
      aux_condition_g1 = sum(
        c_across(cols = matches(ctrl.enids)) > 0,
        na.rm = T
      )
      # #Second condition, at least 1 stat column > 2
      # aux_condition3 = sum(c_across(cols = contains("D1")) > 2),
      # #Second condition, at least 1 exp column > 2
      # aux_condition4 = sum(c_across(cols = contains("D2")) > 2)
    ) |>
    ungroup() |>
    filter(
      aux_condition_g2 > 2 & aux_condition_g1 > 2
    ) |>
    column_to_rownames("gene")
}


getAllRatios = function(xl, num_, denom_) {
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
