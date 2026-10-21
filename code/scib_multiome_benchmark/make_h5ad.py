
# Create AnnData files from the files exported by export.R.

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread
from scipy.sparse.linalg import eigsh

import anndata as ad
import scanpy as sc


# Paths


ROOT = Path("scib_multiome_benchmark")
DATA = ROOT / "data"
EMB = ROOT / "embeddings"



def read_mtx(filename):
    print(f"Reading {filename}")
    return mmread(filename).tocsr()


def read_embedding(filename):
    df = pd.read_csv(filename)

    cells = df.iloc[:, 0].astype(str).values
    values = df.iloc[:, 1:].to_numpy(dtype=np.float32)

    return cells, values


def align_embedding(adata, filename):
    cells, values = read_embedding(filename)

    lookup = pd.Index(cells)

    missing = adata.obs_names.difference(lookup)

    if len(missing) > 0:
        raise ValueError(
            f"{len(missing)} cells in AnnData are missing "
            f"from {filename}"
        )

    indices = lookup.get_indexer(adata.obs_names)

    if np.any(indices < 0):
        raise ValueError(
            f"Could not align all cells for {filename}"
        )

    return values[indices]


# Metadata

print("Reading metadata...")

meta = pd.read_csv(
    DATA / "cell_metadata.csv"
)

meta["cell_id"] = meta["cell_id"].astype(str)

meta = meta.set_index("cell_id")

for column in meta.columns:
    meta[column] = meta[column].astype(str)


# RNA




rna_counts = read_mtx(
    DATA / "rna_counts.mtx"
)

rna_features = pd.read_csv(
    DATA / "rna_features.csv"
)["gene"].astype(str).values

rna_cells = pd.read_csv(
    DATA / "rna_cells.csv"
)["cell_id"].astype(str).values

if rna_counts.shape != (
    len(rna_features),
    len(rna_cells),
):
    raise ValueError(
        "RNA matrix dimensions do not match "
        "feature/cell files."
    )

rna = ad.AnnData(
    X=rna_counts.T.tocsr()
)

rna.obs_names = rna_cells
rna.var_names = rna_features

rna.obs = meta.loc[rna.obs_names].copy()

rna.var_names_make_unique()
rna.obs_names_make_unique()

rna.layers["counts"] = rna.X.copy()

# Normalize RNA for scIB's non-integrated X
print("Normalizing RNA...")

sc.pp.normalize_total(
    rna,
    target_sum=1e4
)

sc.pp.log1p(rna)

# Identify HVGs.

sc.pp.highly_variable_genes(
    rna,
    n_top_genes=3000,
    flavor="seurat"
)


# RNA embeddings


rna.obsm["X_pca"] = align_embedding(
    rna,
    EMB / "rna_pca.csv"
)

rna.obsm["X_rna_harmony"] = align_embedding(
    rna,
    EMB / "rna_harmony.csv"
)


# Save RNA


rna.write_h5ad(
    DATA / "RNA.h5ad",
    compression="gzip"
)



# ATAC


atac_counts = read_mtx(
    DATA / "atac_counts.mtx"
)

atac_features = pd.read_csv(
    DATA / "atac_features.csv"
)["peak"].astype(str).values

atac_cells = pd.read_csv(
    DATA / "atac_cells.csv"
)["cell_id"].astype(str).values

if atac_counts.shape != (
    len(atac_features),
    len(atac_cells),
):
    raise ValueError(
        "ATAC matrix dimensions do not match "
        "feature/cell files."
    )

atac = ad.AnnData(
    X=atac_counts.T.tocsr()
)

atac.obs_names = atac_cells
atac.var_names = atac_features

atac.obs = meta.loc[atac.obs_names].copy()

atac.var_names_make_unique()
atac.obs_names_make_unique()

atac.layers["counts"] = atac.X.copy()

# ATAC embeddings


atac.obsm["X_lsi"] = align_embedding(
    atac,
    EMB / "atac_lsi.csv"
)

atac.obsm["X_atac_harmony"] = align_embedding(
    atac,
    EMB / "atac_harmony.csv"
)

# Save ATAC


atac.write_h5ad(
    DATA / "ATAC.h5ad",
    compression="gzip"
)



# WNN graph


print("Reading WNN graph...")

wsnn = read_mtx(
    DATA / "wsnn.mtx"
)

wnn_cells = pd.read_csv(
    DATA / "wnn_cells.csv"
)["cell_id"].astype(str).values

if wsnn.shape != (
    len(wnn_cells),
    len(wnn_cells),
):
    raise ValueError(
        "WNN graph dimensions do not match cell list."
    )

# Align graph to RNA cell order
wnn_lookup = pd.Index(wnn_cells)
indices = wnn_lookup.get_indexer(rna.obs_names)

if np.any(indices < 0):
    raise ValueError(
        "Could not align WNN graph cells with RNA cells."
    )

wsnn = wsnn[indices][:, indices]



# Symmetrize graph


print("Symmetrizing WNN graph...")

wsnn = (wsnn + wsnn.T) * 0.5

wsnn = wsnn.tocsr()

# Remove self edges
wsnn.setdiag(0)
wsnn.eliminate_zeros()



# Spectral embedding


print("Computing spectral WNN embedding...")

degree = np.asarray(
    wsnn.sum(axis=1)
).ravel()

degree[degree == 0] = 1.0

inv_sqrt_degree = 1.0 / np.sqrt(degree)

D_inv_sqrt = sp.diags(
    inv_sqrt_degree
)

normalized_adjacency = (
    D_inv_sqrt @ wsnn @ D_inv_sqrt
)

# Number of dimensions for WNN benchmark
n_components = 30

if normalized_adjacency.shape[0] <= n_components + 1:
    raise ValueError(
        "Not enough cells for a 30-dimensional "
        "WNN spectral embedding."
    )

eigenvalues, eigenvectors = eigsh(
    normalized_adjacency,
    k=n_components + 1,
    which="LA"
)

# Sort eigenvectors from largest eigenvalue
order = np.argsort(
    eigenvalues
)[::-1]

eigenvectors = eigenvectors[:, order]

# Remove the trivial first eigenvector
wnn_embedding = eigenvectors[:, 1:n_components + 1]

wnn_embedding = np.asarray(
    wnn_embedding,
    dtype=np.float32
)

print(
    "WNN embedding shape:",
    wnn_embedding.shape
)



# Create benchmark AnnData


benchmark = rna.copy()

benchmark.obsm["X_rna_pca"] = rna.obsm["X_pca"]
benchmark.obsm["X_rna_harmony"] = rna.obsm["X_rna_harmony"]

benchmark.obsm["X_atac_lsi"] = (
    atac.obsm["X_lsi"]
)

benchmark.obsm["X_atac_harmony"] = (
    atac.obsm["X_atac_harmony"]
)

benchmark.obsm["X_wnn"] = wnn_embedding

# Store graph too, for reference
benchmark.obsp["wnn_connectivities"] = wsnn

benchmark.uns["benchmark_info"] = {
    "batch_key": "sample_id",
    "label_key": "Cell_Type",
    "wnn_embedding": (
        "30-dimensional spectral embedding "
        "derived from the Seurat wsnn graph"
    )
}

benchmark.write_h5ad(
    DATA / "benchmark.h5ad",
    compression="gzip"
)

print("Benchmark embeddings:")
print(
    list(benchmark.obsm.keys())
)
