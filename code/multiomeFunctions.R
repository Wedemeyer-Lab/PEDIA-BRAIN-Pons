################################################################################
## Multiome Functions
## Desc: Functions used to create the multiome data in the pons atlas.
## Last Updated: 9/18/24
################################################################################

# move and rename files function
move.and.rename <- function(from, to) {
  todir <- dirname(to)
  if (!isTRUE(file.info(todir)$isdir)) dir.create(todir, recursive=TRUE)
  file.rename(from = from,  to = to)
}
  


# Create seurat objects function
create_sobj <- function(name, input.path, out.base, annotations, input.format = "cellranger", project = "Seurat_Project") {
  
  # create paths
  if (!endsWith(out.base, "/")) {
    out.base <- paste0(out.base, "/")
  }
  # create paths
  if (!endsWith(input.path, "/")) {
    input.path <- paste0(input.path, "/")
  }
  
  # define paths
  
  if (input.format == "cellranger") {
    dir <- paste0(input.path, name, "/outs")
  } else {
    dir <- paste0(input.path, name)
  }
  
  output.path <- paste0(out.base, "individual/", name)
  frag.path <- paste(dir, "atac_fragments.tsv.gz", sep = "/")
  mat.path <- paste(dir, "filtered_feature_bc_matrix", sep = "/")
  metadata.atac.path <- paste(dir, "per_barcode_metrics.csv", sep = "/")
  
  # read data
  data <- Read10X(mat.path)
  atac_counts <- data$Peaks
  rna_counts <- data$`Gene Expression`
  
  
  # Get fragment file
  frag.file <- CreateFragmentObject(
    frag.path,
    cells = NULL,
    validate.fragments = T, 
    verbose = TRUE
  )
  
  # use only peaks mapping to standard chromosomes
  ## limit to barcodes present in the out of transcriptome step
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use),]
  
  # Load ATAC-seq per barcode (cell) metadata
  metadata.atac <- read.csv(
    file = metadata.atac.path,
    header = TRUE,
    row.names = 1
  )
  
  # Create chromatin assay
  atac_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = 'hg38',
    fragments = frag.file,
    min.cells = 10,
    annotation = annotations,
    meta.data = metadata.atac
  )
  
  # Create Seurat object
  sobj <- CreateSeuratObject(
    counts = atac_assay,
    assay = 'ATAC',
    project = project
  )
  
  Annotation(sobj) <- annotations
  
  # add assay to seurat object
  rna_assay <- CreateAssayObject(rna_counts)
  sobj[["RNA"]] <- rna_assay
  sobj@meta.data$sample_id <- name
  
  # Save RDS of seurat object
  output <- paste(output.path, ".rds", sep = "")
  saveRDS(object= sobj, file = output)
  
  message(paste0(name, " RDS saved."))
  return(sobj)
} # end create_sobj

#########################################################################################################
# cellCycleScoring

cellCycleScoring <- function(sobj) {
  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  DefaultAssay(sobj) <- "RNA"
  sobj <- NormalizeData(sobj)
  sobj <- CellCycleScoring(sobj, s.features = s.genes, g2m.features = g2m.genes)
  sobj$CC.Difference <- sobj$S.Score - sobj$G2M.Score
  return(sobj)
} # end cellCycleScoring

#########################################################################################################
# sobj_qc

plotQC <-function(sobj, frag.base.path = NULL, frag.full.path = NULL, out ="./", 
                  sample=NULL, sample.column = NULL, find.markers= TRUE, annotations = NULL,
                  run.macs2 = FALSE, macs2.path = NULL, group.by.macs2 = "seurat_clusters", merged = FALSE,
                  input.format = "cellranger") {
  
  # need to have column "sample_id"
  # Takes in list of samples and generate gene expression quality control plots and metrics for each sample.
  if (!merged) {
    if (is.null(sample) & !is.null(sample.column)) {
      sample <- sobj@meta.data[,sample.column][1]
    } else if (is.null(sample) & is.null(sample.column)) {
      stop("Need sample or sample.column arguments.")
    }
    
    if (run.macs2 & is.null(macs2.path)) {
      stop("When running MACS2, the 'macs2.path' argument is required.")
    }
    
    # create paths
    if (!endsWith(out, "/")) {
      out <- paste0(out, "/")
    }
    
    out.path <- paste0(out, sample)
    dir.create(out.path)
    if (!is.null(frag.base.path) & !is.null(frag.full.path)) {
      stop("Need frag.base.path or frag.full.path arguments.")
    } else if (!(is.null(frag.base.path) | is.null(frag.full.path))) {
      stop("Need only frag.base.path or frag.full.path arguments.")
    } else if (!is.null(frag.base.path)) {
      if (!endsWith(frag.base.path, "/")) {
        frag.base.path <- paste0(frag.base.path, "/")
      }
      if (input.format == "cellranger") {
        frag.path <- paste0(frag.base.path, sample, "/outs/atac_fragments.tsv.gz")
      } else {
        frag.path <- paste0(frag.base.path, sample, "/atac_fragments.tsv.gz")
      }
    } else if (!is.null(frag.full.path)) {
      frag.path <- frag.full.path
    }
    
    message(paste0("Running ", sample, "."))
    
    # QC for RNA
    
    DefaultAssay(sobj) <- "RNA"
    
    # store mitochondrial percentage in object meta data
    sobj <- PercentageFeatureSet(sobj, pattern = "MT-", col.name = "percent.mt")
    
    
    # QC for snATAC-seq
    
    DefaultAssay(sobj) <- "ATAC"
    
    sobj <- NucleosomeSignal(sobj)
    sobj <- TSSEnrichment(sobj, fast = FALSE)
    
    # get blacklist ratio
    sobj$blacklist_ratio <- FractionCountsInRegion(
      object = sobj, 
      assay = 'ATAC',
      regions = blacklist_hg38
    )
    
    # Calcualte the fragments in peaks
    
    total_fragments <- CountFragments(frag.path)
    rownames(total_fragments) <- total_fragments$CB
    sobj$atac_fragments <- total_fragments[colnames(sobj), "frequency_count"]
    
    sobj <- FRiP(
      object = sobj,
      assay = 'ATAC',
      total.fragments = 'atac_fragments'
    )
    
    VlnPlot(
      object = sobj,
      features = c('nFeature_ATAC', 'nCount_ATAC', 'TSS.enrichment','nucleosome_signal','blacklist_ratio','FRiP'),
      pt.size = 0.1,
      ncol = 3
    )
    ggsave(paste0(out.path,"/ATACviolinplot.pdf"), height = 8, width = 10)
    
    
    # Create density scatter plot
    DefaultAssay(sobj) <- "ATAC"
    DensityScatter(sobj, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)+ ylim(c(0,15))
    ggsave(paste0(out.path,"/DensityScatter.pdf"), height = 8, width = 10)
    
    # Create TSS plot
    sobj$high.tss <- ifelse(sobj$TSS.enrichment > 2.5, 'High', 'Low')
    TSSPlot(sobj, group.by = 'high.tss') + NoLegend()
    ggsave(paste0(out.path,"/TSSplot.pdf"), height = 8, width = 10)
    
    # Create Fragment Histogram
    sobj$nucleosome_group <- ifelse(sobj$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
    FragmentHistogram(object = sobj, group.by = 'nucleosome_group')
    ggsave(paste0(out.path,"/FragmentHistogram.pdf"), height = 8, width = 10)
    
    # Create Violin plot
    metrics <- c("TSS.enrichment","nucleosome_signal","blacklist_ratio","FRiP",
                 "nFeature_RNA","nCount_RNA", "nFeature_ATAC","nCount_ATAC", "percent.mt")
    VlnPlot(object = sobj,
            features = metrics,
            pt.size = 0.1,
            ncol = 3,
    )
    ggsave(paste0(out.path,"/violinPlot.pdf"), units = "in",height = 10, width = 12)
    
    ##############################################################################
    # pre-processing and dimensional reduction on both assays independently, 
    # using standard approaches for RNA and ATAC-seq data.
    
    # RNA analysis
    DefaultAssay(sobj) <- "RNA"
    sobj <- sobj %>%
      SCTransform(vst.flavor = 'v2') %>% 
      RunPCA(reduction.name = "pca", reduction.key = "PCArna_") %>%
      RunUMAP(reduction="pca", reduction.name = "umap.rna", reduction.key = "UMAPrna_", dims = 1:30) %>%
      FindNeighbors(dims = 1:30)
    
    # test cluster resolution
    sobj <- FindClusters(sobj, resolution=0.5)
    DimPlot(sobj, label = TRUE)
    ggsave((paste0(out.path,"/UMAP_RNA.pdf")))
    
    # Find all markers
    
    if (find.markers) {
      DefaultAssay(sobj) <- "RNA"
      markers.rna <- FindAllMarkers(sobj, only.pos = TRUE)
      write.csv(markers.rna,paste0(out.path, "/markers.rna.csv"))
    }
    
    # Find variable features
    DefaultAssay(sobj) <- "RNA"
    sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
    # plot variable features 
    LabelPoints(VariableFeaturePlot(sobj), points = head(VariableFeatures(sobj), 10), repel = TRUE)
    ggsave((paste0(out.path,"/variableFeaturesRNA.pdf")),height = 8.5, width = 11)
    # Feature plot with QC metrics
    
    sobj_copy <- sobj
    sobj_copy$nucleosome_signal[is.infinite(sobj_copy$nucleosome_signal)] <- 0
    FeaturePlot(sobj_copy, 
                features = metrics,
                pt.size = 0.4,
                label = TRUE
    )
    ggsave((paste0(out.path,"/featurePlot.pdf")), height = 8.5, width = 11)
    rm(sobj_copy)
    
    #top10 genes
    if(find.markers) {
      DefaultAssay(sobj) <- "SCT"
      markers.sobj <- markers.rna
      markers.sobj %>%
        group_by(cluster) %>%
        dplyr::filter(avg_log2FC > 1) %>%
        slice_head(n = 10) %>%
        ungroup() -> top10
      DoHeatmap(sobj, features = top10$gene) + NoLegend()
      ggsave((paste0(out.path,"/topHeatmap.pdf")),width = 28,height = 14)
      
      
      markers.rna %>%
        group_by(cluster) %>%
        slice_head(n = 5) %>%
        ungroup() -> top
      rna.list <- lapply(levels(sobj), function(cluster) {
        cluster.markers <- markers.rna[which(markers.rna$cluster==cluster),"gene"]
      })
      names(rna.list) <- paste0("cluster", 0:(length(unique(sobj$seurat_clusters)) - 1))
      saveRDS(rna.list,paste0(out.path, "/rna.list.rds"))
    }
    
    # ATAC analysis
    # We exclude the first dimension as this is typically correlated with sequencing depth
    DefaultAssay(sobj) <- "ATAC"
    
    sobj <- RunTFIDF(sobj)
    sobj <- FindTopFeatures(sobj, min.cutoff = 'q0')
    sobj <- RunSVD(sobj)
    sobj <- RunUMAP(sobj,
                    reduction = 'lsi',
                    dims = 2:50,
                    reduction.name = "umap.atac", 
                    reduction.key = "UMAPatac_")
    
    
    if (find.markers) {
      DefaultAssay(sobj) <- "ATAC"
      markers.atac <- FindAllMarkers(object = sobj, only.pos = TRUE)
      
      markers.atac %>%
        group_by(cluster) %>%
        slice_head(n = 5) %>%
        ungroup() -> topATAC
      write.csv(markers.atac,paste0(out.path, "/markers.atac.csv"))
      write.csv(topATAC, paste0(out.path, "/topMarkers.atac.csv"))
    }
    
    
    # Find variable features
    DefaultAssay(sobj) <- "ATAC"
    sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
    # plot variable features
    LabelPoints(VariableFeaturePlot(sobj), points = head(VariableFeatures(sobj), 10), repel = TRUE)
    ggsave(paste0(out.path,"/variableFeaturesATAC.pdf"),height = 8.5, width = 11)
    
    # plot both umaps
    umap.rna.plot <- DimPlot(sobj, reduction="umap.rna", label=T) + NoLegend() 
    umap.atac.plot <- DimPlot(sobj, reduction="umap.atac", label=T) + NoLegend() 
    ggsave(paste0(out.path, "/umapRNA_ATAC.pdf"), 
           plot = umap.rna.plot + umap.atac.plot,
           units = "in",
           height = 8.5,
           width = 11)
    message(paste0("QC plots saved to ", out.path))
    if (run.macs2) {
      runMACS2(sobj, macs2.path=macs2.path, frag.path = frag.path,out.path=out.path,
               group.by=group.by.macs2, sample=sample, find.markers= find.markers, annotations=annotations)
    }
  } else {
    # create paths
    if (!endsWith(out, "/")) {
      out <- paste0(out, "/")
    }
    
    out.path <- out
    dir.create(out.path)
    # Perform log-normalization and feature selection, as well as SCT normalization on global object
    DefaultAssay(sobj) <- "RNA"
    sobj <- sobj %>%
      NormalizeData() %>%
      FindVariableFeatures(selection.method = "vst", nfeatures=3000) %>%
      ScaleData() %>%
      SCTransform(vars.to.regress = NULL)
    
    # Calculate PCs using variable features determined by SCTransform (3000 by default)
    sobj <- RunPCA(sobj, assay = "SCT", npcs = 50)

    sobj <- sobj  %>%
      RunUMAP(assay = "SCT", dims = 1:50) %>%
      FindNeighbors() %>%
      FindClusters(resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2), verbose = TRUE)
    # QC for snATAC-seq
    
    DefaultAssay(sobj) <- "ATAC"
    
    sobj <- NucleosomeSignal(sobj)
    sobj <- TSSEnrichment(sobj, fast = FALSE)
    
    # get blacklist ratio
    sobj$blacklist_ratio <- FractionCountsInRegion(
      object = sobj,
      assay = 'ATAC',
      regions = blacklist_hg38
    )
    
    # store mitochondrial percentage in object meta data
    DefaultAssay(sobj) <- "RNA"
    sobj <- PercentageFeatureSet(sobj, pattern = "MT-", col.name = "percent.mt")
    
    
    DefaultAssay(sobj) <- "RNA"
    
    
    # Feature plot with QC metrics
    metrics <- c("TSS.enrichment","nucleosome_signal","blacklist_ratio",
                 "nFeature_RNA","nCount_RNA","percent.mt", "CC.Difference", "S.Score", "G2M.Score")
    
    FeaturePlot(sobj, 
                features = metrics,
                pt.size = 0.4,
                label = TRUE
    )
    ggsave(paste0(out.path,"/featurePlot.pdf"),
           units = "in",
           height = 8.5,
           width = 11)
    DimPlot(sobj, group.by ="SCT_snn_res.0.6",label = T) + ggtitle("res0.6 umap")
    ggsave(paste0(out.path,"/umap_merge.pdf"),
           units = "in",
           height = 8.5,
           width = 11)
    DimPlot(sobj, group.by ="SCT_snn_res.0.6",label = T) + 
      ggtitle("res0.6 UMAP and Feature plot of QC metrics") + 
      theme(legend.position = "none")+ 
      FeaturePlot(sobj, 
                  features = metrics,
                  pt.size = 0.4,
                  label = F)
    ggsave(paste0(out.path,"/metrics_features_merge.pdf"), 
           units = "in",
           height = 12,
           width = 20)
    
    # Create density scatter plot
    DefaultAssay(sobj) <- "ATAC"
    DensityScatter(sobj, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)+ ylim(c(0,15))
    ggsave(paste0(out.path,"/density_scatterplot_merge.pdf"), 
           units = "in",
           height = 8.5,
           width = 11
    )
    
    # Create TSS plot
    sobj$high.tss <- ifelse(sobj$TSS.enrichment > 2.5, 'High', 'Low')
    TSSPlot(sobj, group.by = 'high.tss') + NoLegend()
    ggsave(paste0(out.path,"/TSS_plot_merge.pdf"), 
           units = "in",
           height = 8.5,
           width = 11)
    
    # Create Fragment Histogram
    sobj$nucleosome_group <- ifelse(sobj$nucleosome_signal > 2, 'NS > 2', 'NS < 2')
    FragmentHistogram(object = sobj, group.by = 'nucleosome_group')
    ggsave(paste0(out.path,"/Fragment_Histogram_merge.pdf"), 
           units = "in",
           height = 8.5,
           width = 11)
    
    # Create Violin plot
    VlnPlot(
      object = sobj,
      features = c("TSS.enrichment","nucleosome_signal","blacklist_ratio",
                   "nFeature_RNA","nCount_RNA","percent.mt", "CC.Difference", "S.Score", "G2M.Score"),
      pt.size = 0.1,
      ncol = 3,
    )
    ggsave(paste0(out.path,"/violin_plot_merge.pdf"),
           units = "in",
           height = 10,
           width = 12
    )
    
    DefaultAssay(sobj) <- "RNA"
    sobj@meta.data$sample_trunc <- sub("_MW_.*", "", sobj@meta.data$sample)
    ggsave(paste0(out.path,"/violin_plot0.6_merge.pdf"), 
           plot = VlnPlot(
             object = sobj,
             features = c("TSS.enrichment","nucleosome_signal","blacklist_ratio",
                          "nFeature_ATAC","nCount_ATAC","percent.mt", "CC.Difference", "S.Score", "G2M.Score"),
             group.by ="sample_trunc",
             pt.size = 0.1,
             ncol = 3
           ),
           units = "in",
           height = 10,
           width = 12
    )
    sobj@meta.data$sample_trunc <- NULL
    
  }
  
  return(sobj)
}

################################################################################
# runMACS22

runMACS2<- function(sobj, macs2.path, frag.base.path = NULL, frag.full.path = NULL, out=NULL, out.path=NULL,
                    group.by = "seurat_clusters", sample=NULL, 
                    sample.column = NULL,find.markers=TRUE, annotations=NULL) {
  # call peaks using MACS2 grouping by clusters
  # perform this step after clustering to peak calling by cluster
  
  # call peaks using MACS2 grouping by clusters
  # perform this step after clustering to peak calling by cluster
  
  if (is.null(out) & is.null(out.path)) {
    stop("Need 'out' or 'out.path' arguments.")
  } else if (!is.null(out)) {
    if (is.null(sample) & !is.null(sample.column)) {
      sample <- sobj@meta.data[,sample.column][1]
    } else if (is.null(sample) & is.null(sample.column)) {
      stop("When specifying 'out' argument, need 'sample' or 'sample.column' arguments.")
    }
    if (!endsWith(out, "/")) {
      out <- paste0(out, "/")
    }
    out.path <- paste0(out, sample)
    dir.create(out.path)
  }
  
  out.path <- paste0(out, sample)
  dir.create(out.path)
  if (!is.null(frag.base.path) & !is.null(frag.full.path)) {
    stop("Need frag.base.path or frag.full.path arguments.")
  } else if (!(is.null(frag.base.path) | is.null(frag.full.path))) {
    stop("Need only frag.base.path or frag.full.path arguments.")
  } else if (!is.null(frag.base.path)) {
    if (!endsWith(frag.base.path, "/")) {
      frag.base.path <- paste0(frag.base.path, "/")
    }
    frag.path <- paste0(frag.base.path, sample, "/outs/atac_fragments.tsv.gz")
  } else if (!is.null(frag.full.path)) {
    frag.path <- frag.full.path
  }
  
  if(is.null(annotations)) {
    # get annotation 
    annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
    seqlevelsStyle(annotations) <- 'UCSC'
    genome(annotations) <- "hg38"
  }
  
  #call peaks
  message(paste0("Running MACS2 on ", sample, "."))
  DefaultAssay(sobj) = "ATAC"
  peaksMACS2 <- CallPeaks(sobj, 
                          group.by = group.by,
                          macs2.path=macs2.path)
  
  # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  peaksMACS2 <- keepStandardChromosomes(peaksMACS2, pruning.mode = "coarse")
  peaksMACS2 <- subsetByOverlaps(x = peaksMACS2, ranges = blacklist_hg38_unified, invert = TRUE)
  
  # quantify counts in each peak
  
  frag.file <- CreateFragmentObject(
    frag.path,
    cells = NULL,
    validate.fragments = T, 
    verbose = TRUE
  )
  
  
  macs2_counts <- FeatureMatrix(
    fragments = frag.file,
    features = peaksMACS2,
    cells = colnames(sobj)
  )
  
  
  # create a new assay using the MACS2 peak set and add it to the Seurat object
  sobj[["peaksMACS2"]] <- CreateChromatinAssay(
    counts = macs2_counts,
    fragments = frag.file,
    annotation = annotations
  )
  
  if (find.markers) {
    # Find peaks that are significantly different between clusters
    DefaultAssay(sobj) <- "peaksMACS2"
  markers.macs2 <- FindAllMarkers(object = sobj, only.pos = T, min.pct = 0.25)
  head(markers.macs2)
  
  markers.macs2 %>%
    group_by(cluster) %>%
    slice_head(n = 5) %>%
    ungroup() -> topMACS2
  write.csv(markers.macs2, paste0(out.path, "/markers.macs2.csv"))
  write.csv(topMACS2, paste0(out.path, "/topMarkers.macs2.csv"))
  }
  
  # Find variable features
  DefaultAssay(sobj) <- "peaksMACS2"
  sobj <- FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
  # plot variable features 
  LabelPoints(VariableFeaturePlot(sobj), points = head(VariableFeatures(sobj), 10), repel = TRUE)
  ggsave(paste0(out.path, "/macs2VariableFeatures.pdf"), 
         units = "in",
         height = 8.5,
         width = 11)  
  message("Plots saved to ", out.path,".")
  return(sobj)
}



