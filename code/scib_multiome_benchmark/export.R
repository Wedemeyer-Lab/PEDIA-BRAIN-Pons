# Export RNA/ATAC counts, metadata, reductions, and WNN graph


library(Seurat)
library(Signac)
library(Matrix)

# Configuration


outdir <- "scib_multiome_benchmark"

dir.create(outdir, showWarnings = FALSE)
dir.create(file.path(outdir, "data"), showWarnings = FALSE)
dir.create(file.path(outdir, "embeddings"), showWarnings = FALSE)


sobj <- readRDS("output/objects/atlas_new_celltype.RDS")

required_meta <- c(
  "sample_id",
  "Cell_Type"
)

missing_meta <- setdiff(required_meta, colnames(sobj[[]]))

if (length(missing_meta) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_meta, collapse = ", ")
  )
}

required_reductions <- c(
  "pca",
  "rna_harmony",
  "lsi",
  "atac_harmony"
)

missing_reductions <- setdiff(
  required_reductions,
  Reductions(sobj)
)

if (length(missing_reductions) > 0) {
  stop(
    "Missing reductions: ",
    paste(missing_reductions, collapse = ", ")
  )
}

required_graphs <- c("wknn", "wsnn")

missing_graphs <- setdiff(
  required_graphs,
  names(sobj@graphs)
)

if (length(missing_graphs) > 0) {
  stop(
    "Missing WNN graphs: ",
    paste(missing_graphs, collapse = ", ")
  )
}


# Metadata


meta <- sobj@meta.data

meta$sample_id <- as.character(meta$sample_id)
meta$Cell_Type <- as.character(meta$Cell_Type)

if ("orig.ident" %in% colnames(meta)) {
  meta$orig.ident <- as.character(meta$orig.ident)
}

if ("donor" %in% colnames(meta)) {
  meta$donor <- as.character(meta$donor)
}

meta$cell_id <- rownames(meta)

write.csv(
  meta,
  file = file.path(outdir, "data", "cell_metadata.csv"),
  row.names = FALSE,
  quote = TRUE
)


# RNA counts


message("Extracting RNA counts...")

rna_counts <- GetAssayData(
  sobj,
  assay = "RNA",
  layer = "counts"
)

writeMM(
  rna_counts,
  file = file.path(outdir, "data", "rna_counts.mtx")
)

write.csv(
  data.frame(
    gene = rownames(rna_counts)
  ),
  file = file.path(outdir, "data", "rna_features.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(
    cell_id = colnames(rna_counts)
  ),
  file = file.path(outdir, "data", "rna_cells.csv"),
  row.names = FALSE,
  quote = FALSE
)


# ATAC counts


message("Extracting ATAC counts...")

atac_counts <- GetAssayData(
  sobj,
  assay = "ATAC",
  layer = "counts"
)

writeMM(
  atac_counts,
  file = file.path(outdir, "data", "atac_counts.mtx")
)

write.csv(
  data.frame(
    peak = rownames(atac_counts)
  ),
  file = file.path(outdir, "data", "atac_features.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(
    cell_id = colnames(atac_counts)
  ),
  file = file.path(outdir, "data", "atac_cells.csv"),
  row.names = FALSE,
  quote = FALSE
)


# Helper for exporting reductions


export_reduction <- function(object, reduction, filename) {
  
  message("Exporting reduction: ", reduction)
  
  emb <- Embeddings(
    object,
    reduction = reduction
  )
  
  df <- data.frame(
    cell_id = rownames(emb),
    emb,
    check.names = FALSE
  )
  
  write.csv(
    df,
    file = file.path(
      outdir,
      "embeddings",
      filename
    ),
    row.names = FALSE,
    quote = FALSE
  )
  
  invisible(emb)
}

# RNA representations


rna_pca <- export_reduction(
  sobj,
  "pca",
  "rna_pca.csv"
)

rna_harmony <- export_reduction(
  sobj,
  "rna_harmony",
  "rna_harmony.csv"
)

# ATAC representations


atac_lsi <- export_reduction(
  sobj,
  "lsi",
  "atac_lsi.csv"
)

atac_harmony <- export_reduction(
  sobj,
  "atac_harmony",
  "atac_harmony.csv"
)

# WNN graphs



wknn <- sobj@graphs$wknn
wsnn <- sobj@graphs$wsnn

# Make sure graph dimensions correspond to cells
if (nrow(wknn) != ncol(sobj) ||
    ncol(wknn) != ncol(sobj)) {
  
  stop("wknn dimensions do not match number of cells.")
}

if (nrow(wsnn) != ncol(sobj) ||
    ncol(wsnn) != ncol(sobj)) {
  
  stop("wsnn dimensions do not match number of cells.")
}

writeMM(
  wknn,
  file = file.path(outdir, "data", "wknn.mtx")
)

writeMM(
  wsnn,
  file = file.path(outdir, "data", "wsnn.mtx")
)

write.csv(
  data.frame(
    cell_id = colnames(sobj)
  ),
  file = file.path(outdir, "data", "wnn_cells.csv"),
  row.names = FALSE,
  quote = FALSE
)


# Cell-order checks


message("Checking cell ordering...")

stopifnot(
  identical(
    rownames(meta),
    colnames(rna_counts)
  )
)

stopifnot(
  identical(
    rownames(meta),
    colnames(atac_counts)
  )
)

stopifnot(
  identical(
    rownames(rna_pca),
    rownames(rna_harmony)
  )
)

stopifnot(
  identical(
    rownames(rna_pca),
    rownames(atac_lsi)
  )
)

stopifnot(
  identical(
    rownames(rna_pca),
    rownames(atac_harmony)
  )
)

