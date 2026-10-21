
# Run scIB-metrics on RNA, ATAC, and WNN representations.

from pathlib import Path

import pandas as pd
import scanpy as sc
import numpy as np

from scib_metrics.benchmark import Benchmarker

ROOT = Path("scib_multiome_benchmark")
DATA = ROOT / "data"
RESULTS = ROOT / "results"

RESULTS.mkdir(
    parents=True,
    exist_ok=True
)



#Loading benchmark.h5ad

adata = sc.read_h5ad(
    DATA / "benchmark.h5ad"
)

print("")
print("Cells:", adata.n_obs)
print("Genes:", adata.n_vars)

print("")
print("Metadata:")
print(
    adata.obs[
        ["orig.ident", "Cell_Type"]
    ].head()
)
## all
adata_subset = adata.copy()


for key in adata_subset.obsm.keys():
    print(
        f"  {key}: "
        f"{adata_subset.obsm[key].shape}"
    )



# Configuration


BATCH_KEY = "orig.ident"
LABEL_KEY = "Cell_Type"

EMBEDDINGS = [
    "X_rna_pca",
    "X_rna_harmony",
    "X_atac_lsi",
    "X_atac_harmony",
    "X_wnn"
]


# Verify metadata


if BATCH_KEY not in adata_subset.obs:
    raise ValueError(
        f"Missing batch key: {BATCH_KEY}"
    )

if LABEL_KEY not in adata_subset.obs:
    raise ValueError(
        f"Missing label key: {LABEL_KEY}"
    )


# Verify embeddings


for key in EMBEDDINGS:

    if key not in adata_subset.obsm:
        raise ValueError(
            f"Missing embedding: {key}"
        )

    if adata_subset.obsm[key].shape[0] != adata_subset.n_obs:
        raise ValueError(
            f"{key} has incorrect number of cells."
        )



# Print batch/label structure


print("")
print("Batch counts:")
print(
    adata_subset.obs[BATCH_KEY]
    .value_counts()
)

print("")
print("Cell-type counts:")
print(
    adata_subset.obs[LABEL_KEY]
    .value_counts()
)


# Benchmark



bm = Benchmarker(
    adata_subset,
    batch_key=BATCH_KEY,
    label_key=LABEL_KEY,
    embedding_obsm_keys=EMBEDDINGS,
    n_jobs=-1,
    progress_bar=True,
)


bm.prepare()

print("Running benchmark...")

bm.benchmark()



# Results

results = bm.get_results(
    min_max_scale=False,
    clean_names=True
)


print(results)


# ------------------------------------------------------------
# Save results
# ------------------------------------------------------------

results.to_csv(
    RESULTS / "scib_results_orig.ident.csv"
)

# Also save transposed version
results.T.to_csv(
    RESULTS / "scib_results_transposed_orig.ident.csv"
)

## Fetal
adata_subset = adata[adata.obs["postnatal"] == "Fetal"].copy()



for key in adata_subset.obsm.keys():
    print(
        f"  {key}: "
        f"{adata_subset.obsm[key].shape}"
    )



# Configuration

BATCH_KEY = "orig.ident"
LABEL_KEY = "Cell_Type"

EMBEDDINGS = [
    "X_rna_pca",
    "X_rna_harmony",
    "X_atac_lsi",
    "X_atac_harmony",
    "X_wnn"
]



# Verify metadata

if BATCH_KEY not in adata_subset.obs:
    raise ValueError(
        f"Missing batch key: {BATCH_KEY}"
    )

if LABEL_KEY not in adata_subset.obs:
    raise ValueError(
        f"Missing label key: {LABEL_KEY}"
    )



# Verify embeddings


for key in EMBEDDINGS:

    if key not in adata_subset.obsm:
        raise ValueError(
            f"Missing embedding: {key}"
        )

    if adata_subset.obsm[key].shape[0] != adata_subset.n_obs:
        raise ValueError(
            f"{key} has incorrect number of cells."
        )


# Benchmark


bm = Benchmarker(
    adata_subset,
    batch_key=BATCH_KEY,
    label_key=LABEL_KEY,
    embedding_obsm_keys=EMBEDDINGS,
    n_jobs=-1,
    progress_bar=True,
)


bm.prepare()

print("Running benchmark...")

bm.benchmark()



# Results


results = bm.get_results(
    min_max_scale=False,
    clean_names=True
)


print(results)


# Save results

results.to_csv(
    RESULTS / "fetal_scib_results_orig.ident.csv"
)

# Also save transposed version
results.T.to_csv(
    RESULTS / "fetal_scib_results_transposed_orig.ident.csv"
)

## postnatal
adata_subset = adata[adata.obs["postnatal"] == "Postnatal"].copy()


for key in adata_subset.obsm.keys():
    print(
        f"  {key}: "
        f"{adata_subset.obsm[key].shape}"
    )


# Configuration



EMBEDDINGS = [
    "X_rna_pca",
    "X_rna_harmony",
    "X_atac_lsi",
    "X_atac_harmony",
    "X_wnn"
]


# Verify metadata


if BATCH_KEY not in adata_subset.obs:
    raise ValueError(
        f"Missing batch key: {BATCH_KEY}"
    )

if LABEL_KEY not in adata_subset.obs:
    raise ValueError(
        f"Missing label key: {LABEL_KEY}"
    )

# Verify embeddings
for key in EMBEDDINGS:

    if key not in adata_subset.obsm:
        raise ValueError(
            f"Missing embedding: {key}"
        )

    if adata_subset.obsm[key].shape[0] != adata_subset.n_obs:
        raise ValueError(
            f"{key} has incorrect number of cells."
        )

# Benchmark



bm = Benchmarker(
    adata_subset,
    batch_key=BATCH_KEY,
    label_key=LABEL_KEY,
    embedding_obsm_keys=EMBEDDINGS,
    n_jobs=-1,
    progress_bar=True,
)

bm.prepare()

print("Running benchmark...")

bm.benchmark()

# Results


results = bm.get_results(
    min_max_scale=False,
    clean_names=True
)

print(results)


# Save results


results.to_csv(
    RESULTS / "postnatal_scib_results_orig.ident.csv"
)

# Also save transposed version
results.T.to_csv(
    RESULTS / "postnatal_scib_results_transposed_orig.ident.csv"
)
