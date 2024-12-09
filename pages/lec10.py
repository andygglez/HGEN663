import streamlit as st

col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("Epigenetics, chromatin remodeling, histone modifications", divider=True)

    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Copy the files for today's class from `/home/hgen_share/lec10`")
            
            st.markdown("Set up environment and copy over files")
            st.code("""
            module load StdEnv/2020 samtools/1.11 r/4.0.2
            export PATH="/home/hgen_share/Anaconda/bin:$PATH"
            export PATH="/home/hgen_share/Anaconda/envs/chip_seq_V2/bin:$PATH"
            """, language="bash")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Assess the profile of KO_H3K27ac.chr1.bam near promoters")
            
            st.markdown("Produce coverage track with `bamCoverage`")
            st.code("""
            bamCoverage -b KO_H3K27ac.chr1.bam -o KO_H3K27ac.chr1.bw
            """, language="bash")

            st.markdown("Aggregate signal around promoter regions using `computeMatrix`")
            st.code("""
            computeMatrix reference-point -R promoter.chr1.bed -S KO_H3K27ac.chr1.bw -o KO_H3K27ac.promoter.chr1.mat.gz --referencePoint center -bs 100 -a 10000 -b 10000
            """, language="bash")

            st.markdown("Visualize output through `plotHeatmap`")
            st.code("""
            plotHeatmap -m KO_H3K27ac.promoter.chr1.mat.gz \\
                        -o KO_H3K27ac.promoter.chr1.pdf \\
                        --legendLocation none --yAxisLabel "" \\
                        --xAxisLabel "" --colorMap magma \\
                        --regionsLabel "" --refPointLabel "Promoter"
            """, language="bash")

            st.markdown("Download `.bw` and `.pdf` result files to your local computer")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Call peaks in KO_H3K27ac.chr1.bam with MACS")
            
            st.markdown("Identify significantly enriched regions using `macs3`")
            st.code("""
            macs3 callpeak -t KO_H3K27ac.chr1.bam -c KO_Input.chr1.bam -f BAM -g 2.4e8 -n KO_H3K27ac.chr1 -q 0.01
            """, language="bash")

            st.markdown("Download `.narrowPeak` file to your local computer")
    st.divider()
    #############################################################################################
    st.markdown("### Let's jump to R!!!")

    with st.container(border=True):
            st.markdown("#### Set up")
            
            st.code("""
            library(data.table)
            library(tidyverse)
            library(lattice)
            library(gridExtra)
            library(rtracklayer)
            library(DiffBind)
            library(idr2d)
            library(patchwork)
            library(reactable)

            # set your working directory
            setwd("~/Documents/HGEN_663/extra/lec10")
            """, language="r")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### IDR")

        st.markdown("We read in MACS peak calls and select coordinate columns as well as the p-value column, then run IDR")
        st.code("""
        idr <- c('KO','KO2') %>%
        paste0('_H3K27ac_peaks.narrowPeak') %>%
        lapply(fread, select = c(1:3,8)) %>%
        {estimate_idr1d(.[[1]], .[[2]], value_transformation = 'identity')} %>%
        .$rep1_df

        ggplot(idr, aes(x = rank, y = rep_rank, color = idr)) +
        geom_point(size = .1) + scale_color_gradientn(colors = rainbow(10)) +
        facet_grid(.~'rank') + theme(legend.position = 'none') +
        ggplot(idr, aes(x = value, y = rep_value, color = idr)) +
        geom_point(size = .1) + scale_color_gradientn(colors = rainbow(10)) +
        facet_grid(.~'value')
        """, language="r")
        st.image("images/lec10.idr.png")

    st.divider()

    ##########################################################################################

    with st.container(border=True):
        st.markdown("#### Count matrix")

        st.markdown("We can take the called peaks and merge them to obtain a reference set of genomic intervals – into which we can count reads into just as we’ve done with genes. This can be done easily with DiffBind using a sample sheet containing the locations and BAM / peak files as well as some metadata.")
        st.code("""
        # You don't need to run this chunk in class
        dat <- dba(sampleSheet = "samples.csv")
        dat <- dba.count(dat)
        """, language="r")

        st.markdown("But since the counting process takes a bit of time, we’ll start with the pre-computed object (lec10.RData). As with transcriptomic and methylomic matrices, we can again perform routine analyses such as PCA and hierchical clustering")
        st.code("""
        load('lec10.RData')
        """, language="r")

    st.divider()

    ##########################################################################################

    with st.container(border=True):
        st.markdown("#### PCA")

        st.code("""
        dba.plotPCA(dat, label = DBA_CONDITION)
        """, language="r")
        st.image("images/lec10.PCA.png")

    st.divider()
        
    ##########################################################################################

    with st.container(border=True):
        st.markdown("#### Correlation")

        st.code("""
        dba.plotHeatmap(dat, correlations = T)
        """, language="r")
        st.image("images/lec10.corr.png")
        
    st.divider()
        
    ##########################################################################################

    with st.container(border=True):
        st.markdown("#### Hierarchical clustering")

        st.markdown("")
        st.code("""
        """, language="r")
        st.image("images/lec10.hier.clus.png")

    st.divider()
        
    ##########################################################################################

    with st.container(border=True):
        st.markdown("#### Genome-wide summaries")

        st.markdown("When we have many regions of interest across the genome, we can assess the distribution of signals surrounding such areas and evaluate the overall trend. One tool for performing this is deepTools, and we can examine the results of their computeMatrix module.")
        st.code("""
        # You don't need to run this chunk in class
        library(parallel)
        mats <- mclapply(list.files('mats', full.names = T), function(x)
        colMeans(fread(x, skip = 1, header = F, drop = 1:6),
                na.rm = T), mc.cores = 12) %>%
        setNames(list.files('mats'))
        mat.ex <- fread('mats/KO_H3K27ac.rx.promoter.mat.gz', skip = 1, header = F)
        """, language="r")

    st.divider()
        
    ##########################################################################################

    with st.container(border=True):
        st.markdown("#### Heatmap")

        st.markdown("We’ll first take the KO sample’s H3K27ac enrichment surrounding promoters as an example")
        st.code("""
        # use only the matrix part, throw away the coordinates info
        dat.ex <- mat.ex[,-(1:6)] %>%
        arrange(rowMeans(.)) %>%
        as.matrix

        # plot
        levelplot(t(dat.ex), useRaster = T, xlab = "Position", ylab = "Interval",
                scales = list(draw = F), col.regions = clrs,
                at = seq(-2, 3, length.out = 100), aspect = "fill")
        """, language="r")
        st.image("images/lec10.heatmap.png")

    st.divider()
        
    ##########################################################################################

    with st.container(border=True):
        st.markdown("#### Including expression")

        st.markdown("Given the link between promoter acetylation and the cognate gene’s expression, we further incorporate RNA-seq data. Again, for the sake of time, it’s been imported and provided in lec10.RData")
        st.code("""
        # You don’t need to run this chunk in class
        exps <- fread('BT245-KO-C4.UCSC.featureCounts.txt.counts')
        library(biomaRt)
        ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                        host = "grch37.ensembl.org",
                        path = "/biomart/martservice",
                        dataset = "hsapiens_gene_ensembl")
        res <- getBM(attributes = c("external_gene_name", "chromosome_name",
                                    "start_position", "end_position"),
                    filters = "external_gene_name", values = exps$Geneid,
                    mart = ensembl)
        dict <- left_join(res[!duplicated(res$external_gene_name) &
                                res$external_gene_name %in% exps$Geneid,],
                        exps[exps$Geneid %in% res$external_gene_name],
                        by = c("external_gene_name" = "Geneid")) %>%
        filter(chromosome_name %in% c(1:22, "X", "Y"))
        """, language="r")

        st.markdown("With the expresion data in hand, we now combine the two modalities")
        st.code("""
        # get promoters intervals from matrix
        pmts <- makeGRangesFromDataFrame(mat.ex,
                                        seqnames.field = 'V1',
                                        start.field = 'V2',
                                        end.field = 'V3')

        # get gene intervals overlapping promoters
        dict$chr <- paste0("chr", dict$chromosome_name)
        genes <- makeGRangesFromDataFrame(dict,
                                        start.field = "start_position",
                                        end.field = "end_position",
                                        seqnames.field = "chr",
                                        keep.extra.columns = T)
        pmts.ol <- pmts[overlapsAny(pmts, genes)]
        genes.ol <- genes[overlapsAny(genes, pmts)]

        # partition genes into quintiles by expression
        ranks <- rank(genes.ol$counts/abs((end(genes.ol) - start(genes.ol))),
                    ties.method = "first")
        genes.ol$cut <- as.numeric(cut(ranks, quantile(ranks, probs = seq(0, 1, 0.2)),
                            include.lowest = TRUE))
        pmts.cut <- data.frame(V4 = sprintf("%s:%d-%d", seqnames(pmts.ol),
                                            start(pmts.ol),end(pmts.ol)),
                            cut = genes.ol$cut[findOverlaps(pmts.ol, genes.ol,
                                                            select = "first")])

        # integrate
        dat.ex <- left_join(pmts.cut, mat.ex[,-c(1:3, 5, 6)]) %>%
        select(-V4) %>%
        arrange(cut, rowMeans(select(., -cut))) %>%
        split(., .$cut)
        """, language="r")
        
    st.divider()
        
    ##########################################################################################

    with st.container(border=True):
        st.markdown("#### Heatmaps")

        st.markdown("We can plot the results as heatmaps again")
        st.code("""
        # plot heatmap for each quintile
        lapply(1:5, function(i) {
        levelplot(t(dat.ex[[i]][,2:601]), useRaster = T,
                                xlab = NULL, ylab = NULL, scales = list(draw = F),
                                col.regions = clrs, aspect = "fill",
                                at = seq(-2, 3, length.out = 100),
                                colorkey = F, margin = F)
        }) %>% grid.arrange(grobs = ., ncol = 1)
        """, language="r")
        st.image("images/lec10.grid.heatmaps.png")

    ##########################################################################################

    with st.container(border=True):
        st.markdown("#### Aggregate plots")

        st.markdown("Or alternatively take the column-wise average across said heatmap as a summary")
        st.image("images/lec10.agg.plots.png")

    st.divider()

    ##############################################################################
    with st.container(border=True):
        st.markdown("#### Correlation")

        st.markdown("Or more directly we can just assess the relationship per-promoter")
        st.code("""
        ## Joining with `by = join_by(V4)`
        """, language="r")
        st.image("images/lec10.correlation.png")

    st.divider()
    ##############################################################################
    with st.container(border=True):
        st.markdown("#### Normalization")

        st.markdown("We’ve been only looking at one sample thus far, and to contrast difference samples we’ll need to make sure the scales are meaningfully comparable")

        st.markdown("#### Aggregate plot")
        st.markdown("Taking our previous promoter aggregate plots, different methods would point to rather divergent conclusions")
        st.code("""
        # define groups
        ann <- tibble(f = names(mats)) %>%
        separate(f, c('samps', 'marks', 'norms', 'regs', NA, NA), '[_\\.]', F) 

        # what to plot
        mark <- 'H3K27ac'
        reg <- 'promoter'

        # label groups
        dats <- melt(mats[ann$f[ann$marks == mark & ann$regs == reg]]) %>%
        rownames_to_column("row") %>%
        mutate(row = (as.numeric(row) - 1) %% (n()/6) + 1) %>%
        merge(dplyr::rename(ann, L1 = f)) %>%
        mutate(norms = factor(norms, levels = c('raw', 'inp', 'rx'))) %>%
        select(-L1)

        # label axis
        len <- nrow(dats) / 6 / 3
        if (len * 10 < 1000) {
        lab.l <- sprintf("-%db", len)
        lab.r <- sprintf("+%db", len)
        } else {
        lab.l <- sprintf("-%dkb", len / 100)
        lab.r <- sprintf("+%dkb", len / 100)
        }

        # plot
        ggplot(dats, aes(x = row, y = value, color = samps)) +
        annotate("rect", xmin = len, xmax = len * 2, ymin = -Inf,
                ymax = Inf, alpha = 0.2, fill = "grey") +
        geom_line() +
        labs(y = "Enrichment", title = sprintf("%s @ %s", mark, reg)) +
        scale_x_continuous(breaks = c(1, len * 1:3),
                            labels = c(lab.l, "Start", "End", lab.r),
                            expand = c(0, 0), name = "Position") +
        facet_wrap(. ~ norms, scales = "free") + theme_bw()
        """, language="r")
        st.image("images/lec10.plots.png")

    st.divider()
    ##############################################################################
    with st.container(border=True):
        st.markdown("#### Count in windows")

        st.markdown("As opposed to limiting ourselves to some set of pre-defined regions, we could also just take the counts in uniformly divided windows across the genome, which is stored cts object (a list of bedGraph files imported like below)")
        st.code("""
        # You don’t need to run this chunk in class
        cts <- lapply(list.files('10kb', full.names = T), import.bedGraph) %>%
        setNames(sub('.10kb.bed', '', list.files('10kb')))
        """, language="r")

        st.markdown("We now apply various normalization factors")
        st.code("""
        # read in scaling factors
        facs <- fread('~/Documents/HGEN_663/extra/lec10/rx.csv')
        """, language="r")

        st.code("""
        facs <- facs %>%
        mutate(sample = paste(facs$condition, facs$mark, sep = "_")) %>%
        rowwise() %>%
        mutate(ref = min(hs_chip, hs_input)) %>%
        ungroup() %>%
        mutate(depthfac1 = ref / hs_chip,
                depthfac2 = ref / hs_input,
                rx = (hs_chip / hs_input) / (dm_chip / dm_input),
                scalefac1 = rx * depthfac1)
        facs %>% reactable()
        """, language="r")
        st.image("images/lec10.reactable.png")

        st.code("""
        # normalize
        for (samp in facs$sample) {

        # get ratios
        rx <- facs$scalefac1[facs$sample == samp]
        r1 <- facs$depthfac1[facs$sample == samp]
        r2 <- facs$depthfac2[facs$sample == samp]
        
        # scale by factor
        num.input <- r1 * cts[[samp]]$score + 1
        num.rx <- rx * cts[[samp]]$score + 1
        denom <- r2 * cts[[sub('_.*', '_input', samp)]]$score + 1
        cts[[samp]]$rx <- log2(num.rx / denom)
        cts[[samp]]$input <- log2(num.input / denom)
        }
        """, language="r")

        st.markdown("And take H3K27me3 as an example")
        st.code("""
        # compare counts
        mark <- 'H3K27me3'

        samp1 <- paste('KO', mark, sep = '_')
        samp2 <- paste('K27M', mark, sep = '_')
        smoothScatter(x = log2(cts[[samp1]]$score + 1),
                    y = log2(cts[[samp2]]$score + 1),
                    xlab = samp1, ylab = samp2, main = "Raw")
        abline(0,1)
        """, language="r")
        st.image("images/lec10.K27me3.png")

        st.markdown("Or as MA plots")
        st.code("""
        ma <- list()
        ma[["raw"]] <- data.frame(m = log(cts[[samp1]]$score) - log(cts[[samp2]]$score),
                                a = 0.5 * (log(cts[[samp1]]$score) + log(cts[[samp2]]$score)))
        ma[["input"]] <- data.frame(m = cts[[samp1]]$input - cts[[samp2]]$input,
                                    a = 0.5 * (cts[[samp1]]$input + cts[[samp2]]$input))
        ma[["rx"]] <- data.frame(m = cts[[samp1]]$rx - cts[[samp2]]$rx,
                                    a = 0.5 * (cts[[samp1]]$rx + cts[[samp2]]$rx))

        par(mfrow = c(1,3))
        for(nm in names(ma)) {
        smoothScatter(x = ma[[nm]]$a, y = ma[[nm]]$m, ylab = "M", xlab = "A",
                        main = sprintf("%s - %s (%s)", samp1, samp2, nm))
        abline(h = 0)
        }
        """, language="r")
        st.image("images/lec10.MAplot.png")

    st.divider()
    ##############################################################################
    with st.container(border=True):
        st.markdown("#### Your turn")

        st.markdown("Do other regions or marks exhibit different trends?")
        st.markdown("""
        1) Plot H3K27ac on enhancers.
        2) Plot H3K36me2 on genes, promoters and enhancers.
        3) Plot H3K27me3 on genes, promoters and enhancers.
        """)
        st.code("""
        # define groups
        ann <- tibble(f = names(mats)) %>%
        separate(f, c('samps', 'marks', 'norms', 'regs', NA, NA), '[_\\.]', F) 

        # check what other regulatory regions are avaiable
        ann$regs

        # check what other marks you can plot
        ann$marks

        # what to plot
        mark <- ''
        reg <- ''

        # label groups
        dats <- 

        # label axis
        len <- nrow(dats) / 6 / 3
        if (len * 10 < 1000) {
        lab.l <- sprintf("-%db", len)
        lab.r <- sprintf("+%db", len)
        } else {
        lab.l <- sprintf("-%dkb", len / 100)
        lab.r <- sprintf("+%dkb", len / 100)
        }

        # plot
        ggplot()
        """, language="r")

        st.markdown("Do the other marks also seem adequately normalized with any of the approaches?")
        st.markdown("Other marks: H3K36me2, H3K27ac")
        st.code("""
        # compare counts
        mark <- ''

        samp1 <- paste('KO', mark, sep = '_')
        samp2 <- paste('K27M', mark, sep = '_')
        smoothScatter(x = log2(cts[[samp1]]$score + 1),
                    y = log2(cts[[samp2]]$score + 1),
                    xlab = samp1, ylab = samp2, main = "Raw")
        abline(0,1)
        #
        ma <- list()
        ma[["raw"]] <- data.frame(m = log(cts[[samp1]]$score) - log(cts[[samp2]]$score),
                                a = 0.5 * (log(cts[[samp1]]$score) + log(cts[[samp2]]$score)))
        ma[["input"]] <- data.frame(m = cts[[samp1]]$input - cts[[samp2]]$input,
                                    a = 0.5 * (cts[[samp1]]$input + cts[[samp2]]$input))
        ma[["rx"]] <- data.frame(m = cts[[samp1]]$rx - cts[[samp2]]$rx,
                                    a = 0.5 * (cts[[samp1]]$rx + cts[[samp2]]$rx))

        par(mfrow = c(1,3))
        for(nm in names(ma)) {
        smoothScatter(x = ma[[nm]]$a, y = ma[[nm]]$m, ylab = "M", xlab = "A",
                        main = sprintf("%s - %s (%s)", samp1, samp2, nm))
        abline(h = 0)
        }
        """, language="r")
    st.divider()
