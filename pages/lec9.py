import streamlit as st

col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("Epigenetics, DNA methylation and analysis", divider=True)

    #############################################################################################

    with st.container(border=True):
        st.markdown("Set up directory and copy over files")

        st.code("data=/project/60006/hgen_share/lec9", language="bash")

        st.markdown("#### Use dnmtools to identify hypo-methylated regions (HMRs), differentially methylated regions (DMRs) and partially methylated domains (PMDs)")
        st.markdown("""
        Check out the first few lines from GSC.chr3.meth
        """)
        st.code("""
        head ${data}/GSC.chr3.meth
        """, language="bash")

        st.markdown("Calculate differential methylation with `dnmtools`")
        st.code("""
        module load apptainer
        singularity run $data/dnmtools.sif diff -o chr3.methdiff ${data}/EpiLC.chr3.meth ${data}/GSC.chr3.meth
        cp $data/DMR* .
        """, language='bash')

        st.markdown("Identify HMRs with `hmr`. Do the same for EpiLC")
        st.code("""
        singularity run $data/dnmtools.sif hmr -o GSC.chr3.hmr ${data}/GSC.chr3.meth
        """, language='bash')

        st.markdown("Call DMRs usign `dmr`")
        st.code("""
        singularity run $data/dnmtools.sif dmr chr3.methdiff ${data}/EpiLC.chr3.hmr GSC.chr3.hmr ${data}/DMR_EpiLC_lt_GSC.bed ${data}/DMR_GSC_lt_EpiLC.bed
        """, language='bash')

        st.markdown("Determine PMDs with `pmd`. Do the same for EpiLC")
        st.code("""
        singularity run $data/dnmtools.sif pmd -o GSC.chr3.pmd.bed $data/GSC.chr3.meth
        """, language='bash')

        st.markdown("Add a .bed extension to the .hmr files. Download .bw & .bed result files to your local computer")

    st.divider()
    #############################################################################################

    st.markdown("### Let's jump to R!!!")

    with st.container(border=True):
        st.markdown("#### Set up")
        
        st.code("""
        library(data.table)
        library(tidyverse)
        library(SummarizedExperiment)
        library(pheatmap)
        library(knitr)

        # set your working directory
        setwd("~/Documents/HGEN_663/extra/lec9")
        """, language="r")
    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Prepare data")
        
        st.markdown("As opposed to working from raw data, today we will be working with those pre-processed by TCGA")

        st.markdown("#### Downloading from TCGA")
        st.markdown("Rather than retrieving the dataset piece-by-piece manually, we can access it more programmatically. Since it may take a while to download, you’ve been provided with the output of this chunk")
        st.code("""
        library(FEM)
        library(data.table)
        library(tidyverse)
        library(SummarizedExperiment)
        library(pheatmap)
        library(TCGAbiolinks)
        library(biomaRt)

        #downloading from microarray from TCGA
        library(TCGAbiolinks) #retrieval process from TGCA is simplified using this package
        library(FEM)
        library(biomaRt)
        samples <- data.frame(barcodes = c("TCGA-BB-A5HZ","TCGA-CN-4739","TCGA-CN-4727","TCGA-CV-5441",
                                        "TCGA-CV-6441","TCGA-P3-A6T6","TCGA-UF-A71D","TCGA-BA-A4IF",
                                        "TCGA-CR-5247","TCGA-BA-5555","TCGA-BB-4217","TCGA-CQ-5323",
                                        "TCGA-DQ-5624","TCGA-CN-4722","TCGA-P3-A5Q5","TCGA-CV-5430",
                                        "TCGA-BA-A4IG","TCGA-HD-7229","TCGA-CR-5243","TCGA-BA-4075"),
                            names = c("NSD1m_1","NSD1m_2","NSD1m_3","NSD1m_4","NSD1m_5","NSD1m_6",
                                        "NSD1m_7","NSD1m_8","NSD1m_9","NSD1m_10","noNSD1m_1",
                                        "noNSD1m_2","noNSD1m_3","noNSD1m_4","noNSD1m_5","noNSD1m_6",
                                        "noNSD1m_7","noNSD1m_8","noNSD1m_9","noNSD1m_10"),
                            stringsAsFactors = F) 

        query.meth <- GDCquery(project = "TCGA-HNSC",
                       data.category = "DNA Methylation",
                       platform = "Illumina Human Methylation 450",
                       data.type = "Methylation Beta Value",
                       barcode = samples$barcodes)

        GDCdownload(query.meth)
        data.meth <- GDCprepare(query = query.meth, save = TRUE,
                                save.filename = "meth.RData")

        query.exp <- GDCquery(project = "TCGA-HNSC",
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification",
                            workflow.type = "STAR - Counts",
                            barcode = samples$barcodes)
        GDCdownload(query.exp)
        data.exp <- GDCprepare(query = query.exp, save = TRUE,
                        save.filename = "exp.RData")

        # Additional annotation
        data("probe450kfemanno")
        probes <- as.data.frame(simplify2array(probe450kfemanno))

        # Ensembl Gene ID <-> Probe ID
        ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
        tx <- unique(sub('\\..*$','', rowData(data.meth)$gene))

        map <- getBM(attributes = c('ensembl_transcript_id','ensembl_gene_id'),
                    filters = 'ensembl_transcript_id', values = tx,
                    mart = ensembl)
        
        dict <- data.frame(probe = rownames(rowData(data.meth)),
                   transcript = sub('\\..*$','',rowData(data.meth)$gene))
        dict$gene <- plyr::mapvalues(dict$transcript, map$ensembl_transcript_id, map$ensembl_gene_id)
        """, language="r")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Import")
        
        st.markdown("We can now load in the prepared data and make some further adjustments")
        st.code("""
        load('lec9.RData')

        # Discard normal tissue samples
        data.meth <- subset(data.meth, select = data.meth$shortLetterCode != "NT")

        # Replace barcodes with readable names
        colnames(data.meth) <- colnames(data.exp) <- samples$names

        # Discard probes with all NA values
        mat.meth <- assay(data.meth) %>%
        subset(rowSums(is.na(.)) != ncol(.))
        """, language="r")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("### Distribution")
        st.markdown("#### Genome-wide")

        st.markdown("Individual samples")

        st.code("""
        #Genome-wide: One of the first things to check would be the distribution of methylation levels across the genome
        # Reshape for plotting
        data.plt <- reshape2::melt(mat.meth, varnames = c("probe","sample")) %>%
        mutate(NSD1 = ifelse(grepl("noNSD1",sample), "+", "-"))

        # Plot individual samples
        to_plot <- data.plt[data.plt$sample %in% c("NSD1m_6","noNSD1m_6"),]
        to_plot <- na.omit(to_plot)

        # Run this line in the console
        ggplot(to_plot, aes(x = value)) + geom_histogram() + labs (x = "Beta value", y = "Count") + facet_grid(. ~ sample)

        # Plot all samples
        ggplot(data.plt, aes(x = sample, y = value)) +
        geom_violin(aes(color = NSD1), show.legend = F) + 
        stat_summary(aes(color = NSD1), fun.y = median, geom = "point", show.legend = F) +
        labs(x = "Sample", y = "Beta value") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        """, language="r")
        st.image("images/lec9.all.samples.png")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Specific regions")
        
        st.code("""
        # Split probes based on feature type
        ft <- split(as.character(rowData(data.meth)$Composite.Element.REF),
            f = rowData(data.meth)$Feature_Type)

        # Split probes based on gene group
        gg <- split(as.character(probes$probeID), f = probes$GeneGroup %>%
                    as.numeric() %>%
                    replace(., is.na(.), 7))
        """, language="r")

        st.markdown("Individual samples")
        st.code("""
        # Pick a particular region
        data.plt.reg <- data.plt[data.plt$probe %in% gg$`5`,]#5 is gene body
        # Individual samples
        ggplot(data.plt.reg[data.plt.reg$sample %in% c("NSD1m_2","noNSD1m_6"),], aes(x = value)) + 
        geom_histogram() + labs (x = "Beta value", y = "Count") + facet_grid(. ~ sample)

        #select another region - subset of probes exhibiting unique pattern
        data.plt.reg <- data.plt[data.plt$probe %in% gg$`6`,]
        ggplot(data.plt.reg[data.plt.reg$sample %in% c("NSD1m_2","noNSD1m_6"),], aes(x = value)) + 
        geom_histogram() + labs (x = "Beta value", y = "Count") + facet_grid(. ~ sample)
        """, language="r")
        st.image("images/lec9.specific.regions.ind.samples.png")

        st.markdown("All samples")
        st.code("""
        # All samples
        ggplot(data.plt.reg, aes(x = sample, y = value)) +
        geom_violin(aes(color = NSD1), show.legend = F) +
        stat_summary(aes(color = NSD1), fun.y = median, geom="point", show.legend = F) +
        labs(x = "Sample", y = "Beta value") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
        """, language="r")
        st.image("images/lec9.specific.regions.all.samples.png")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Clustering")
        
        st.markdown("Given the beta value matrix, we can adopt similar approaches as we’ve used for RNA-seq expression matrices")


        st.markdown("#### Methylation-based")
        st.markdown("For example, both hierarchical clustering and PCA can be directly applied")
        st.code("""
        # HCA on the n most variable probes
        sds <- apply(mat.meth, 1, sd)
        keep <- order(sds, decreasing = T)[1:1000]
        anno.col <- data.frame(NSD1=rep(c("-","+"),each=10))
        rownames(anno.col) <- colnames(mat.meth)
        pheatmap(
            as.matrix(mat.meth[keep,]),
            treeheight_row = 0,
            show_rownames = F,
            annotation_col = anno.col,
            annotation_names_col = F,
            clustering_distance_cols = "correlation",
            clustering_distance_rows = "correlation"
        )
        """, language="r")
        st.image("images/lec9.heatmap.png")

        st.code("""
        # PCA
        res.pca <- prcomp(na.omit(mat.meth))
        comps <- summary(res.pca)$importance[2,1:2] * 100
        res.pca$rotation %>%
        data.frame %>%
        mutate(samples = rownames(.)) %>%
        mutate(condition = ifelse(grepl("noNSD1",samples),"+","-")) %>%
        ggplot(aes(x = PC1, y = PC2, color = condition)) +
            geom_point(size = 5, show.legend = F) +
            geom_text(aes(label = samples), show.legend = F) +
            xlab(sprintf("PC1 (%d %%)",round(comps[2]))) +
            ylab(sprintf("PC2 (%d %%)",round(comps[1])))
        """, language="r")
        st.image("images/lec9.PCA.png")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Interrelation")
        
        st.markdown("Now we proceed and investigate the relationship between expression and methylation")
        st.code("""
        # expression-based
        mat.exp <- log2(assay(data.exp) + 1)

        # Only focus on genes and probes with valid mapping
        mat.meth.ok <- mat.meth[rownames(mat.meth) %in% dict$probe,]
        mat.exp.ok <- mat.exp[rownames(mat.exp) %in% dict$gene,]

        # Look at a particular type of CpG sites
        reg <- gg$`2`
        reg.name <- "TSS"

        # Merge methylation with expression data
        dat.reg <- reshape2::melt(mat.meth.ok[rownames(mat.meth.ok) %in% reg,],
                        varnames = c("probe","samples"), value.name = "beta") %>%
                left_join(., dict, by = "probe") %>%
                left_join(., reshape2::melt(mat.exp.ok,varnames=c("gene","samples")), by = c("gene","samples"))

        # Get mean methylation level of all probes associated with each gene
        samp <- "NSD1m_6"
        dat.plt <- dat.reg[dat.reg$samples == samp,] %>%
        group_by(gene) %>%
        summarise(beta = mean(beta), expr = mean(value))
        smoothScatter(x = dat.plt$expr, y = dat.plt$beta, xlab = "Expression",
                    ylab = "Methylation", main = paste0(reg.name," of ",samp))
        """, language="r")
        st.image("images/lec9.interrelation.png")

        st.code("""
        # Aggregate samples within conditions by taking the mean
        # Aggregate samples within conditions by taking the mean
        dat.wt <- dat.reg[grepl("noNSD1",dat.reg$samples),] %>%
            group_by(gene) %>%
            summarise(beta = mean(beta), expr = mean(value))
        dat.mt <- dat.reg[!grepl("noNSD1",dat.reg$samples),] %>%
            group_by(gene) %>%
            summarise(beta = mean(beta), expr = mean(value))
        dat.all <- dat.reg %>%
            group_by(gene) %>%
            summarise(beta = mean(beta), expr = mean(value))

        par(mfrow=c(1,3),mar=c(5,5,5,1))
        smoothScatter(x = dat.wt$expr, y = dat.wt$beta, xlab = "Expression", ylab = "Methylation", main = "NSD1+")
        smoothScatter(x = dat.mt$expr, y = dat.mt$beta, xlab = "Expression", ylab = "Methylation", main = "NSD1-")
        title(reg.name, line = -1, outer = TRUE)
        """, language="r")
        st.image("images/lec9.smooth.scatter.png")

    st.divider()
    #############################################################################################

    st.markdown("### Your turn!!!")
    with st.container(border=True):
        st.image("images/lec9.yourturn.png")
    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Distribution")
        
        st.markdown("For gene body, generate the distribution of methylation across NSD1 mutant and control samples.")

        st.markdown("#### Plot Individual samples")
        st.markdown("See above for the parameters used for ggplot()")
        st.code("""
        ### look at gene body ###
        # Look at a particular type of CpG sites
        reg <- gg$`???`
        reg.name <- "gene body"

        # Merge methylation with expression data
        dat.reg <-  %>%
        left_join() %>%
        left_join()

        # Get mean methylation level of all probes associated with each gene
        samp <- "NSD1m_6"
        dat.plt <- dat.reg[] %>%
        group_by(gene) %>%
        summarise(beta = , expr =)

        smoothScatter()
        """, language="r")

        st.markdown("What are your conclusions?")

        st.markdown("#### Plot All samples")

        st.code("""
        # Aggregate samples within conditions by taking the mean
        dat.wt
        dat.mt
        dat.all
        
        # Complete the lines below:
        smoothScatter( )
        smoothScatter( )
        """, language="r")
        
    st.divider()
