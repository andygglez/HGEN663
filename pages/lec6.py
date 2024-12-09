import streamlit as st

col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("Transcriptome Analysis", divider=True)

    with st.container(border=True):
            st.markdown("#### Set up")
            
            st.code("""
            library(DESeq2)
            library(apeglm)
            library(gprofiler2)
            library(fgsea)
            library(msigdbr)
            library(org.Hs.eg.db)
            library(limma)
            library(patchwork)
            library(ggrepel)
            library(reactable)
            library(shiny)
            library(tippy)
            library(ggdist)
            library(dbplyr)
            library(tidyverse)

            # set your working directory
            setwd("~/Documents/HGEN_663/extra/lec6")
            """, language="r")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("### Part 1")
            
            st.markdown("#### Import")
            st.code("""
            # load data
            load('lec6_1.rda')

            # create DESeq object
            dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~condition)

            # what samples are we working with?
            coldata
            """, language="r")
            st.image("images/lec6.import.png")

            st.markdown("#### Run analysis")
            st.markdown("How does normalization show up on the MA plot?")
            st.code("""
            dds <- DESeq(dds)

            par(mfrow = c(1, 2))
            limma::plotMA(log2(counts(dds) + 1), array = 8, main = "raw")
            abline(h = 0, col = "grey")
            limma::plotMA(log2(counts(dds, normalized = T)  + 1), array = 8, main = "normalized")
            abline(h = 0, col = "grey")
            """, language="r")
            st.image("images/lec6.run.analysis.png")


    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Get results")
            
            st.code("""
            res <- results(dds)
            head(res)
            """, language="r")
            st.image("images/lec6.get.results.png")


            st.markdown("#### Shrink fold change estimate")
            st.code("""
            resLFC <- lfcShrink(dds, coef = 2, type = "apeglm", res = res)
            """, language="r")

            st.markdown("#### Compare effect of shrinkage")

            st.markdown("#### MA plot")
            st.code("""
            par(mfrow = c(1, 2))
            DESeq2::plotMA(res, ylim = c(-5, 5), maiin = 'Original')
            DESeq2::plotMA(resLFC, ylim = c(-5, 5), main = 'Shrunken')
            """, language="r")
            st.image("images/lec6.MA.plot.png")

            st.markdown("#### Volcano plot")
            st.code("""
            volc <- function(r, x, y, ttl) {
            d <- as.data.frame(r) %>%
                dplyr::rename(x = !!x, y = !!y) %>%
                mutate(kind = case_when(abs(x) > 2 & y < .01 ~ 'DE',
                                        abs(x) > 2 ~ 'big',
                                        y < .01 ~ 'sig',
                                        T ~ 'NS'),
                    y = -log10(y))

            ct <- na.omit(d) %>% 
                dplyr::filter(kind == 'DE') %>%
                mutate(up = x > 0) %>%
                dplyr::count(up) %>%
                mutate(x = ifelse(up, Inf, -Inf),
                    y = Inf,
                    h = as.numeric(up))

            ggplot(d, aes(x, y, color = kind)) +
                geom_vline(xintercept = c(-2, 2), linetype = 'dashed') +
                geom_hline(yintercept = 2, linetype = 'dashed') +
                geom_point(alpha = .5) +
                geom_label(aes(x = x, y = y, label = n, hjust = h),
                        vjust = 1, data = ct, inherit.aes = F) +
                scale_y_continuous() +
                scale_color_manual(values = c('forestgreen', 'red2', 'grey30',  'royalblue')) +
                labs(x = x, y = y, title = ttl)
            }

            volc(res, 'log2FoldChange','padj', 'Original') +
            theme(legend.position = 'none') +
            volc(resLFC, 'log2FoldChange', 'padj', 'Shrunken')
            """, language="r")
            st.image("images/lec6.volcano.plot.png")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Grouping")
            st.markdown("Before delving further, let’s make sure the samples we’re comparing for DE make sense with a simple PCA")
            
            st.markdown("#### PCA")
            st.code("""
            rld <- rlog(dds, blind = F)
            pd <- plotPCA(rld, intgroup = "condition", ntop = 500, returnData = T)

            data.frame(extra) %>%
            rownames_to_column('name') %>% 
            merge(pd) %>%
            ggplot(aes(x = PC1, y = PC2)) +
            geom_point(aes(color = condition)) +
            geom_text(aes(label = name)) +
            labs(x = sprintf('PC1: %.1f%% variance', 100 * attr(pd, 'percentVar')[1]),
                y = sprintf('PC2: %.1f%% variance', 100 * attr(pd, 'percentVar')[2])) +
            coord_fixed()
            """, language="r")
            st.image("images/lec6.MA.plot.png")

            st.markdown("#### Heatmap")
            st.code("""
            assay(rld) %>%
                apply(1, sd) %>%
                order(decreasing = TRUE) %>%
                .[1:100] %>%
                assay(rld)[.,] %>%
                heatmap()
            """, language="r")
            st.image("images/lec6.heatmap.png")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Design")
            
            st.code("""
            colData(dds) <- cbind(coldata, extra)
            colData(dds)$condition <- factor(colData(dds)$condition)
            design(dds) <- ~ sex + condition
            dds <- DESeq(dds)
            res.s <- results(dds)

            design(dds) <- ~ sex + age + condition
            dds <- DESeq(dds)
            res.sa <- results(dds)
            """, language="r")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("### Part 2")
            
            st.markdown("#### Import")
            st.code("""
            # clear workspace
            rm(list = setdiff(ls(), 'volc'))

            # load data
            load('lec6_2.rda')

            # create deseq dataset object
            dds <- DESeqDataSetFromMatrix(countData = round(cts),
                                        colData = coldata,
                                        design = ~ condition)

            # run deseq2
            dds <- DESeq(dds)

            # output all results
            res <- results(dds)
            """, language="r")

            st.markdown("#### DE")
            st.code("""
            t10 <- data.frame(res) %>%
                rownames_to_column("gene") %>%
                slice_min(padj, n = 10, with_ties = F) %>%
                arrange(padj) %>%
                pull(gene)

            counts(dds, normalized = TRUE) %>%
                data.frame() %>%
                rownames_to_column('gene') %>%
                dplyr::filter(gene %in% t10) %>%
                pivot_longer(-gene, names_to = 'samp', values_to = 'ct') %>%
                merge(rownames_to_column(data.frame(coldata), 'samp')) %>%
                mutate(gene = factor(gene, t10)) %>%
                ggplot(aes(x = condition, y = ct, color = condition)) +
                geom_point() +
                scale_y_log10() +
                facet_wrap(.~gene)
            """, language="r")
            st.image("images/lec6.DE.png")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Pathway Overrepresentation Analysis")
            
            st.markdown("")
            st.code("""
            rd <- data.frame(res) %>%
                na.omit() %>%
                rownames_to_column('gene') %>%
                mutate(gene = sub('\\..*', '', gene)) 

            g <- gost(rd$gene[rd$padj < .05 & abs(rd$log2FoldChange) > .5],
                    organism = 'hsapiens',
                    custom_bg = rd$gene)

            g$result %>%
                dplyr::select(source, term_name, term_size, query_size,
                                intersection_size, p_value) %>% 
                arrange(source, p_value) %>%
                mutate(p_value = formatC(p_value, digits = 3, format = 'g')) %>%
                reactable()
            """, language="r")
            st.image("images/lec6.pathway.ORA.png")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Gene Set Enrichment Analysis (GSEA)")
            
            st.markdown("#### Run analysis")
            st.code("""
            # mapping Ensembl IDs to gene symbols
            e2s <- AnnotationDbi::select(org.Hs.eg.db,
                                        key = rd$gene,
                                        columns = "SYMBOL",
                                        keytype = "ENSEMBL") %>%
                    na.omit() %>%
                    deframe()

            # use DE test statistic to rank genes
            s <- rd %>% 
                mutate(symb = e2s[gene]) %>%
                na.omit() %>% 
                group_by(symb) %>%
                summarise(stat = mean(stat)) %>%
                deframe()

            # pathways
            p <- msigdbr(species = "Homo sapiens", category = "H") %>%
                mutate(gs_name = sub('^HALLMARK_', '', gs_name)) %>%
                {split(.$gene_symbol, .$gs_name)}

            # GSEA
            r <- fgsea(pathways = p, stats = s, eps = 0.0)

            r %>%
                arrange(NES) %>% 
                mutate(pathway = fct_inorder(pathway)) %>% 
                ggplot(aes(x = pathway, y = NES)) + 
                geom_col(aes(fill = -log10(padj))) +
                scale_fill_viridis_c() +
                coord_flip() +
                labs(x = "Pathway", y = "Normalized enrichment score")
            """, language="r")

            st.markdown("#### Top hit")
            st.code("""
            plotEnrichment(p[[arrange(r, padj)$pathway[1]]], s) +
                ggtitle(arrange(r, padj)$pathway[1])
            """, language="r")
            st.image("images/lec6.GSEA.plot.png")

            st.markdown("#### Table")
            st.code("""
            arrange(r, padj) %>%
                dplyr::select(pathway, padj, ES, NES, size, leadingEdge) %>%
                mutate_at(c('pathway', 'padj', 'ES', 'NES'), function(x) {
                    formatC(x, digits = 3, format = 'g')
                }) %>%
                reactable(
                    searchable = T,
                    highlight = T,
                    wrap = F,
                    resizable = T,
                    striped = T,
                    paginationType = "jump",
                    showPageSizeOptions = T,
                    defaultPageSize = 10,
                    columns = list(
                    leadingEdge = colDef(
                        html = T,
                        cell =  function(value, index, name) {
                        value <- paste(value, collapse = ', ')
                        div(
                            style = "cursor: info;
                                    white-space: nowrap;
                                    overflow: hidden;
                                    text-overflow: ellipsis;",
                            tippy(text = value, tooltip = value)
                        )
                        }
                    )
                    )
                )
            """, language="r")
            st.image("images/lec6.enrichment.table.png")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Your turn")
            
            st.markdown("Run DESeq on the dataset from previous class. Run comparisons for both knockout conditions (ie. PA versus NSD1KO, and PA versus NSD12DKO)")
            st.markdown("#### From previous class")
            st.code("""
            fc <- read.table("~/Documents/HGEN_663/extra/lec5/ex2_featureCounts.txt",header=1)

            df <- fc %>% 
            dplyr::select(-c(2:6)) %>% 
            tibble::column_to_rownames("Geneid")
            mat <- df %>% as.matrix()
            metadata <- data.frame(kind=colnames(mat))
            metadata$condition <- "PA"
            metadata[3:4,"condition"] <- "NSD1KO"
            metadata[5:6,"condition"] <- "NSD12DKO"
            metadata <- tibble::column_to_rownames(metadata,"kind")
            dds <- DESeqDataSetFromMatrix(countData = mat,
                                        colData = metadata,
                                        design = ~condition)
            """, language="r")

            st.markdown("#### Run DESeq")
            st.code("""
            # run deseq2
            dds <- DESeq(dds)

            # check which comparisons can be made
            resultsNames(dds)

            # results output
            res <- dds %>%
            results(contrast = c('condition', 'NSD1KO', 'PA'))
            """, language="r")

            st.markdown("#### Run dispersion analysis. Does the DESeq model have a good fit to the data?")
            st.code("""
            plotDispEsts(dds)
            """, language="r")

            st.markdown("#### Volcano plot")
            st.code("""
            # use LFC shrinkage
            resLFC <- lfcShrink(dds,contrast=c('condition',test_samp,control_samp),type="ashr")
            volc(resLFC, 'log2FoldChange', 'padj', 'Shrunken')
            """, language="r")

            st.markdown("#### Pathway analysis")
            st.markdown("Do pathway analysis for both NSD1KO and NSD12DKO. Do the overrepresented pathways differ?")
            st.code("""
            # dataframe 
            rd <- res %>%
                data.frame() %>%
                na.omit() %>%
                rownames_to_column('gene') %>%
                mutate(gene = sub('\\..*', '', gene))
                # pathway analysis

            g <- gost()
            """, language="r")

            st.markdown("#### GSEA")
            st.markdown("Run GSEA and look at the top hit for each condition compared to parental.")

    st.divider()