library(biomaRt)

# List available Ensembl marts
listEnsembl()

# Connect to the main Ensembl genes mart
ensembl <- useEnsembl(biomart = "genes")

# List available datasets within the mart (e.g., human, mouse)
listDatasets(ensembl)

# Choose the human gene dataset
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
# A vector of example Ensembl Gene IDs
ensembl_ids <- rownames(xl$q75.filtered.norm.counts)

# Query the database
gene_annotation <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = ensembl_ids,
  mart = ensembl
)

# View the result
print(gene_annotation)

dim(xl$q75.filtered.norm.counts)

filtered_counts = xl$q75.filtered.norm.counts |> rownames_to_column("id")
filtered_counts = left_join(
  filtered_counts,
  gene_annotation,
  by = c("id" = "ensembl_gene_id")
)
filtered_counts.annotation = filtered_counts |>
  dplyr::select(id, external_gene_name)

saveRDS(filtered_counts.annotation, "data/ensembl_to_symbols.RDS")
