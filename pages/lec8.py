import streamlit as st

col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("Single Cell Analysis", divider=True)
    st.markdown("#### Part 1")

    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Set up")
        
        st.code("""
        library(data.table)
        library(tidyverse)
        library(Seurat)
        library(patchwork)
        library(reactable)
        # set your working directory
        setwd("~/Documents/HGEN_663/extra/lec8")
        """, language="r")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### MacCarroll retina")
        
        st.markdown("""
        Here we will analyze single-cell RNA-seq data from mouse retina, obtained at postnatal 
        day 14. Fresh tissue was dissected, dissociated, and the cell suspension was immediately
        subjected to the Drop-seq protocol for cell capture and sequencing.
        
        We will start from the Drop-seq output (a matrix of read counts per gene, per cell) 
        of the MacCarrol sample, and perform basic quality control, dimensionality reduction 
        and clustering. We will then analyze the cell populations obtained, infer their cell-type 
        identity and obtain their gene signatures.
        """)

        st.markdown("#### Import")
        st.markdown("Load gene count estimates per cell, and build the basic Seurat object we’ll be using for analysis.")
        st.code("""
        set.seed(42)
        DGEmatrix <- fread('Dropseq_p14_retina_Mccarroll.txt', sep = "\t") %>% 
        data.frame(row.names = 1)

        # create our base Seurat object, discarding cells that have less than 100 genes
        # and discarding genes that are expressed in less than 3 cells. 
        Seurat_object <- CreateSeuratObject(DGEmatrix, assay = "RNA",
                                            min.cells = 3, min.features = 100,
                                            project = "Mccarroll_retina")

        # number of genes and cells
        dim(DGEmatrix)
        """, language="r")
        st.image("images/lec8.import.png")

        st.code("""
        # check DGE matrix
        DGEmatrix[1:3, 1:3]
        """, language="r")
        st.image("images/lec8.matrix.png")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### QC")
        
        st.markdown("We will compute some QC metrics that reflect the overall quality of the sample: coverage (# of genes, # of UMIs) and mitochondrial read proportions – which is often an indicator of stress/cell damage")
        st.code("""
        # get a list of mitochondrial genes: all genes starting with 'MT'
        mito.genes <- grep("^MT-", rownames(x = Seurat_object@assays$RNA), value = TRUE)

        # compute proportions
        percent.mito <- Matrix::colSums(Seurat_object@assays$RNA[mito.genes,]) /
        Matrix::colSums(Seurat_object@assays$RNA)

        # add the information back to the seurat object
        Seurat_object <- AddMetaData(object = Seurat_object, metadata = percent.mito,
                                    col.name = "percent.mito")

        # get a list of crystal genes (all genes starting by 'CRY') to remove crystalin contamination
        crystal.genes <- grep("^CRY[AB]", rownames(x = Seurat_object@assays$RNA), value = TRUE)

        # compute proportions
        percent.crystal <- Matrix::colSums(Seurat_object@assays$RNA[crystal.genes, ]) / 
        Matrix::colSums(Seurat_object@assays$RNA)

        # add the information back to the seurat object
        Seurat_object <- AddMetaData(object = Seurat_object, metadata = percent.crystal,
                                    col.name = "percent.crystal")

        # display distribution of some metrics:
        # # of genes, # of UMIs, and mitochondrial/crystal proportion
        fts <- c("nFeature_RNA", "nCount_RNA", "percent.mito", "percent.crystal")
        VlnPlot(object = Seurat_object, ncol = 4, pt.size = 0, features = fts)
        """, language="r")
        st.image("images/lec8.QC.png")

        st.code("""
        # Get some summary stats for the sample:
        summary_stats_before_filtering <- tibble(
        total_cells  = nrow(Seurat_object@meta.data),
        mean_n_genes = mean(Seurat_object@meta.data$nFeature_RNA),
        sd_n_genes   = sd(Seurat_object@meta.data$nFeature_RNA),
        max_n_genes  = max(Seurat_object@meta.data$nFeature_RNA),
        min_n_genes  = min(Seurat_object@meta.data$nFeature_RNA),
        mean_UMI     = mean(Seurat_object@meta.data$nCount_RNA),
        sd_UMI       = sd(Seurat_object@meta.data$nCount_RNA),
        max_UMI      = max(Seurat_object@meta.data$nCount_RNA),
        min_UMI      = min(Seurat_object@meta.data$nCount_RNA)
        ) %>% mutate_all(function(x) round(x, 2))

        summary_stats_before_filtering %>% reactable::reactable()
        """, language="r")

        st.image("images/lec8.summary.png")
                
        st.markdown("")
        st.code("""
        print(
            paste0("Number of filtered cells: ",
                summary_stats_before_filtering$total_cells -
                summary_stats_after_filtering$total_cells
            )
        )
        """, language="r")
        st.image("images/lec8.filt.cells.png")

        st.markdown("Not much was removed, only a few outlier cells, because this was a quite good quality data to begin with. In other datasets, the before/after QC metrics will look much more different.")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Normalization")
        
        st.markdown("Normalize the dataset using [SCTransform](https://satijalab.org/seurat/articles/sctransform_vignette.html)")
        st.code("""
        # before proceeding, we note markers for cell cycle
        Seurat_object <- CellCycleScoring(Seurat_object, set.ident = FALSE,
                                        s.features = cc.genes$s.genes,
                                        g2m.features = cc.genes$g2m.genes)

        # we then go on to apply a transforms the data and regresses out unwanted variation
        vs <- c("nFeature_RNA", "percent.mito", "S.Score", "G2M.Score")
        Seurat_object <- SCTransform(Seurat_object, verbose = T,
                                    vars.to.regress = vs)

        FeatureScatter(Seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mito") +
        FeatureScatter(Seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") &
        theme(legend.position = 'none')
        """, language="r")
        st.image("images/lec8.normalization.png")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Dimension reduction")
        
        st.markdown("##### PCA")
        st.markdown("We will perform a first dimension reduction by using only the most variable genes. Then, we’ll transform the data using principal component analysis (PCA), and use only the first few components for downstream analyses.")
        st.code("""
        # a bit of wrangling to prepare for VariableFeaturePlot
        Seurat_object[['SCT']]@meta.features <- SCTResults(Seurat_object[['SCT']], slot = 'feature.attributes')[, c('gmean', 'variance', 'residual_variance')]
        Seurat_object[['SCT']]@meta.features$variable <- F
        Seurat_object[['SCT']]@meta.features[VariableFeatures(Seurat_object[['SCT']]), 'variable'] <- F
        colnames(Seurat_object[["SCT"]]@meta.features) <- paste0("sct.", colnames(Seurat_object[["SCT"]]@meta.features))

        # label top 10 most variable features
        VariableFeaturePlot(Seurat_object, selection.method = 'sct', assay = 'SCT') %>%
        LabelPoints(points = head(VariableFeatures(Seurat_object), 10), repel = T) &
        theme(legend.position = 'none')
        """, language="r")
        st.image("images/lec8.sc.PCA.png")

        st.code("""
        # perform PCA
        Seurat_object <- RunPCA(Seurat_object, pcs.compute = 100, do.print = F)

        # display genes correlated with PCs 1 & 2
        VizDimLoadings(Seurat_object, dims = 1:2, reduction = "pca")
        """, language="r")

        st.code("""
        # alternative representation of genes highly correlated with PC1
        DimHeatmap(Seurat_object, dims = 1:15, cells = 500, balanced = T)
        """, language="r")
        st.image("images/lec8.Dim.Heatmap.png")
        
        st.code("""
        # inspect the standard deviation of each PC
        ElbowPlot(Seurat_object)
        """, language="r")
        st.image("images/lec8.elbowPlot.png")

        st.code("""
        # assess effect of cell cycle
        DimPlot(Seurat_object, reduction = "pca", group.by = "Phase")
        """, language="r")
        st.image("images/lec8.DimPlot.png")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Non-linear approaches")
        
        st.code("""
        Seurat_object <- RunUMAP(Seurat_object, dims = 1:20)
        Seurat_object <- RunTSNE(Seurat_object, dims = 1:20)
        """, language="r")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Clustering")
        
        st.code("""
        # cluster the cells based on UMAP reduction
        Seurat_object <- FindNeighbors(Seurat_object, 
                                    dims = 1:2, 
                                    reduction = "umap")
        # the resolution can be adjusted to tweak clustering accuracy
        Seurat_object <- FindClusters(Seurat_object, 
                                    resolution = 0.5, 
                                    reduction = "umap")
        """, language="r")
        st.image("images/lec8.clust.png")

        st.markdown("UMAP")
        st.code("""
        DimPlot(Seurat_object, reduction = "umap", label = T)
        """, language="r")
        st.image("images/lec8.UMAP.png")

        st.markdown("tSNE")
        st.code("""
        DimPlot(Seurat_object, reduction = "tsne", label = T)
        """, language="r")
        st.image("images/lec8.tSNE.png")

        st.markdown("Metrics")
        st.code("""
        VlnPlot(object = Seurat_object, pt.size = 0.1, ncol = 1, features = fts)
        """, language="r")
        st.image("images/lec8.metrics.png")
    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Identifying cell types")
        
        st.markdown("We now have some clusters of cells (cell populations), but we still don’t know what these cells are. One can apply differential expression analysis to find marker genes higher expressed in every given cluster as compared to all remaining cells, and then infer the cell type based on current knowledge on those genes.")
        st.code("""
        Seurat_object.markers <- FindAllMarkers(Seurat_object, only.pos = TRUE,
                                        min.pct = 0.25, logfc.threshold = 0.25)
        head(Seurat_object.markers) %>% reactable()
        """, language="r")
        st.image("images/lec8.table.png")

        st.markdown("We can take a look at the top 10 markers in each cluster")
        st.code("""
        top10 <- Seurat_object.markers %>% 
        group_by(cluster) %>% 
        slice_max(avg_log2FC, n = 10) %>%
        ungroup()
        DotPlot(
        Seurat_object,
        assay = NULL,
        unique(top10$gene),
        cols = c("blue", "red"),
        col.min = -2.5,
        col.max = 2.5,
        dot.min = 0,
        dot.scale = 6,
        group.by = NULL,
        split.by = NULL,
        scale.by = "radius",
        scale.min = NA,
        scale.max = NA
        ) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
        """, language="r")
        st.image("images/lec8.top10.png")

        st.markdown("Alternatively, one could directly use known marker genes or published expression profiles of purified / sorted cells.")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Known marker genes")
        
        st.code("""
        markers.retina.dotplot <- toupper(c(
        "Rho",  #Rods
        "Opn1mw", #Cones
        "Trpm1", #Bipolar cells
        "Scgn", #Bipolar cells
        "Celf4", #Amacrine cells
        "Rgs5", #Pericytes
        "Cldn5", #Endothelial cells
        "Lyz2", #Immune cells
        "Glul", # Muller glia
        "Gfap", #Astrocytes
        "Optc", #Lens cells
        "Lhx1"  #Horizontal cells
        ))
        """, language="r")

        st.markdown("Example")
        st.code("""
        DefaultAssay(Seurat_object) <- 'RNA'
        FeaturePlot(Seurat_object, features = 'TRPM1', reduction = 'umap')
        """, language="r")
        st.image("images/lec8.known.markers.png")

        st.markdown("Dot plots")
        st.code("""
        DefaultAssay(Seurat_object) <- 'SCT'
        DotPlot(
        Seurat_object,
        assay = NULL,
        markers.retina.dotplot,
        cols = c("blue", "red"),
        col.min = -2.5,
        col.max = 2.5,
        dot.min = 0,
        dot.scale = 6,
        group.by = NULL,
        split.by = NULL,
        scale.by = "radius",
        scale.min = NA,
        scale.max = NA
        ) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
        """, language="r")
        st.image("images/lec8.markers.dotplots.png")

        st.markdown("Violin plots")
        st.code("""
        DefaultAssay(Seurat_object) <- 'SCT'
        VlnPlot(object = Seurat_object, features = markers.retina.dotplot, pt.size = 0.1, stack=TRUE)
        """, language='r')
        st.image("images/lec8.violin.plots.png")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Annotation")
        
        st.markdown("Activity 1: Assign as many cluster labels as you can based on the provided information")
        st.code("""
        DefaultAssay(Seurat_object) <- 'SCT'
        Seurat_object <- SetIdent(Seurat_object, cells = NULL, value="seurat_clusters")
        new.cluster.ids <- c("",     # cluster 0
                            "",     # cluster 1
                            "",     # cluster 2
                            "",     # cluster 3
                            "",     # cluster 4
                            "",     # cluster 5
                            "",     # cluster 6
                            "",     # cluster 7
                            "",     # cluster 8
                            "",     # cluster 9
                            "",     # cluster 10
                            "",     # cluster 11
                            "",     # cluster 12
                            "",     # cluster 13
                            "",     # cluster 14
                            "",     # cluster 15
                            "",     # cluster 16
                            "",     # cluster 17
                            "",     # cluster 18
                            "",     # cluster 19
                            ""      # cluster 20
        )
        names(new.cluster.ids) <- levels(Seurat_object)
        Seurat_object <- RenameIdents(Seurat_object, new.cluster.ids)

        Seurat_object[["Cell_Type"]] <- Idents(object = Seurat_object)

        DimPlot(Seurat_object, reduction = "umap", label = T)
        """, language="r")
    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Joyal retina")
        
        st.markdown("Now that you know how to process single-cell samples, perform QC, identify cell populations and assign identities, you are ready to analyze a dataset on your own! We have another mouse retinal sample from the Joyal Lab.")
        st.code("""
        set.seed(42)
        DGEmatrix <- fread('Dropseq_p14_retina_Mccarroll.txt', sep = "\t") %>% 
        data.frame(row.names = 1)
        """, language="r")
    st.divider()
    #############################################################################################
    st.markdown("#### Part 2")
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Set up")
        
        st.code("""
        library(data.table)
        library(tidyverse)
        library(Seurat)
        library(patchwork)
        library(reactable)
        # set your working directory
        setwd("~/Documents/HGEN_663/extra/lec8")
        """, language="r")
    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Integrating datasets")
        st.markdown("We’ll start with pre-computed Seurat objects.")

        st.markdown("#### Merging")

        st.code("""
        # import MacCarroll dataset
        Seurat_object_MacCarroll <- readRDS("Seurat_object_MacCarroll.rds")
        Seurat_object_MacCarroll[["Dataset"]] <- "MacCarroll"

        # import Joyal dataset
        Seurat_object_Joyal <- readRDS("Seurat_object_Joyal.rds")
        Seurat_object_Joyal[["Dataset"]] <- "Joyal"

        # merge the 2 dataset
        Seurat_object <- merge(Seurat_object_MacCarroll,
                            y = Seurat_object_Joyal,
                            add.cell.ids = NULL,
                            project = "Lab_integration")
        # before proceeding, we note markers for cell cycle
        Seurat_object <- CellCycleScoring(Seurat_object, set.ident = FALSE,
                                        s.features = cc.genes$s.genes,
                                        g2m.features = cc.genes$g2m.genes)
        # we then go on to apply a transformation to the data and regresses out unwanted variation
        vs <- c("nFeature_RNA", "percent.mito", "S.Score", "G2M.Score")
        Seurat_object.list <- SplitObject(Seurat_object, split.by = "Dataset") %>%
        lapply(SCTransform, verbose = F, vars.to.regress = vs)
        """, language="r")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Integration")
        
        st.code("""
        # select most consistently variable features
        Seurat_object.features <- SelectIntegrationFeatures(object.list = Seurat_object.list,
                                                            nfeatures = 3000)

        # some data wrangling
        Seurat_object.list <- PrepSCTIntegration(object.list = Seurat_object.list,
                                                anchor.features = Seurat_object.features)

        # identify anchors shared by the datasets
        Seurat_object.anchors <- FindIntegrationAnchors(object.list = Seurat_object.list,
                                                        normalization.method = "SCT", 
                                                        anchor.features = Seurat_object.features)

        # proceed with integration
        Seurat_object <- IntegrateData(anchorset = Seurat_object.anchors,
                                    normalization.method = "SCT")
        """, language="r")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Dimension reduction")
        
        st.code("""
        Seurat_object <- RunPCA(Seurat_object) %>%
                         RunUMAP(dims = 1:20) %>%
                         RunTSNE(dims = 1:20)
        """, language="r")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Clustering")
        
        st.code("""
        Seurat_object <- FindNeighbors(Seurat_object, 
                               dims = 1:2, 
                               reduction = "umap")
        Seurat_object <- FindClusters(Seurat_object, 
                                    resolution = 0.5, 
                                    reduction = "umap")
        """, language="r")

        st.markdown("##### UMAP")
        st.markdown("Together")
        st.code("""
        DimPlot(Seurat_object, reduction = "umap", label=TRUE)
        """, language="r")
        st.image("images/lec8.2.umap.together.png")

        st.markdown("Split")
        st.code("""
        DimPlot(Seurat_object, reduction = "umap", label=TRUE, split.by="Dataset")
        """, language="r")
        st.image("images/lec8.2.umap.split.png")

        st.markdown("##### tSNE")
        st.markdown("Together")
        st.code("""
        DimPlot(Seurat_object, reduction = "tsne", label=TRUE)
        """, language="r")
        st.image("images/lec8.2.tsne.together.png")

        st.markdown("Split")
        st.code("""
        DimPlot(Seurat_object, reduction = "tsne", label=TRUE, split.by="Dataset")
        """, language="r")
        st.image("images/lec8.2.tsne.split.png")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Identifying cell types")
        
        st.markdown("We could again identify markers genes as before")
        st.code("""
        DefaultAssay(Seurat_object) <- "integrated"
        Seurat_object.markers <- FindAllMarkers(Seurat_object, only.pos = TRUE,
                                                min.pct = 0.25, logfc.threshold = 0.25)
        head(Seurat_object.markers) %>% reactable()
        """, language="r")
        st.image('images/lec8.2.cell.types.png')

        st.code("""
        top10 <- Seurat_object.markers %>% 
        group_by(cluster) %>% 
        slice_max(avg_log2FC, n = 10) %>%
        ungroup()
        DotPlot(
        Seurat_object,
        assay = NULL,
        unique(top10$gene),
        cols = c("blue", "red"),
        col.min = -2.5,
        col.max = 2.5,
        dot.min = 0,
        dot.scale = 6,
        group.by = NULL,
        split.by = NULL,
        scale.by = "radius",
        scale.min = NA,
        scale.max = NA
        ) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
        """, language="r")
        st.image("images/lec8.2.dotplot.png")

        st.markdown("Or again adopt known markers")
        st.code("""
        markers.retina.dotplot <- toupper(c(
        "Rho",  #Rods
        "Opn1mw", #Cones
        "Trpm1", #Bipolar cells
        "Scgn", #Bipolar cells
        "Celf4", #Amacrine cells
        "Rgs5", #Pericytes
        "Cldn5", #Endothelial cells
        "Lyz2", #Immune cells
        "Glul", # Muller glia
        "Gfap", #Astrocytes
        "Optc", #Lens cells
        "Lhx1"  #Horizontal cells
        ))
        DotPlot(
        Seurat_object,
        assay = NULL,
        markers.retina.dotplot,
        cols = c("blue", "red"),
        col.min = -2.5,
        col.max = 2.5,
        dot.min = 0,
        dot.scale = 6,
        group.by = NULL,
        split.by = NULL,
        scale.by = "radius",
        scale.min = NA,
        scale.max = NA
        ) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
        """, language="r")
        st.image("images/lec8.2.dotplot2.png")

        st.markdown("Or make use of pre-existing labels annotated independently")
        st.code("""
        DimPlot(Seurat_object, reduction = "umap", label = TRUE, 
        pt.size = 0.5, group.by="Cell_Type") + NoLegend()
        """, language="r")
        st.image("images/lec8.2.celltype.png")

    st.divider()
    with st.container(border=True):

        st.markdown("#### Annotation")
        st.markdown("**Activity 3**: Assign as many cluster labels as you can based on the provided information")
        st.code("""
        Seurat_object <- SetIdent(Seurat_object, cells = NULL, value="seurat_clusters")
        new.cluster.ids <- c("",     # cluster 0
                            "",     # cluster 1
                            "",     # cluster 2
                            "",     # cluster 3
                            "",     # cluster 4
                            "",     # cluster 5
                            "",     # cluster 6
                            "",     # cluster 7
                            "",     # cluster 8
                            "",     # cluster 9
                            "",     # cluster 10
                            "",     # cluster 11
                            "",     # cluster 12
                            "",     # cluster 13
                            "",     # cluster 14
                            "",     # cluster 15
                            "",     # cluster 16
                            "",     # cluster 17
                            "",     # cluster 18
                            "",     # cluster 19
                            "",     # cluster 20
                            "",     # cluster 21
                            "",     # cluster 22
                            "",     # cluster 23
                            "",     # cluster 24
                            "",     # cluster 25
                            "",     # cluster 26
                            "",     # cluster 27
                            ""      # cluster 28
        )
        names(new.cluster.ids) <- levels(Seurat_object)
        Seurat_object <- RenameIdents(Seurat_object, new.cluster.ids)

        Seurat_object[["Cell_Type"]] <- Idents(object = Seurat_object)

        DimPlot(Seurat_object, reduction = "umap", label = T)
        """, language="r")

    st.divider()