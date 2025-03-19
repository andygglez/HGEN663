import streamlit as st

col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("Epigenetics, chromatin remodeling, histone modifications", divider=True)

    #############################################################################################

    
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Look at the files for the current lecture:")
            
            st.code("""
            data=/project/60006/hgen_share/lec10
            ll $data
            """, language="bash")

            st.markdown("Create a coverage track for the ${data}/KO_H3K27ac.chr1.bam file")

            st.code("""
            module load apptainer
            singularity exec ${data}/deeptools.sif bamCoverage -b ${data}/KO_H3K27ac.chr1.bam -o KO_H3K27ac.chr1.bw
            """)

            st.markdown("Aggregate signal around promoter regions using `computeMatrix`")
            
            st.code("""
            singularity exec ${data}/deeptools.sif \\
                    computeMatrix reference-point \\
                    -R ${data}/promoter.chr1.bed \\
                    -S KO_H3K27ac.chr1.bw \\
                    -o KO_H3K27ac.promoter.chr1.mat.gz \\
                    --referencePoint center \\
                    -bs 100 -a 10000 -b 10000
            """)


            st.markdown("Visualize output through `plotHeatmap`")
            st.code("""
            singularity exec ${data}/deeptools.sif \\
                    plotHeatmap -m KO_H3K27ac.promoter.chr1.mat.gz \\
                    -o KO_H3K27ac.promoter.chr1.pdf \\
                    --legendLocation none --yAxisLabel "" \\
                    --xAxisLabel "" --colorMap magma \\
                    --regionsLabel "" --refPointLabel "Promoter"
            """, language="bash")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Call peaks in KO_H3K27ac.chr1.bam with MACS")
            
            st.markdown("Identify significantly enriched regions using `macs3`")
            st.code("""
            macs3 callpeak -t ${data}/KO_H3K27ac.chr1.bam -c ${data}/KO_Input.chr1.bam -f BAM -g 2.4e8 -n KO_H3K27ac.chr1 -q 0.01
            """, language="bash")

    st.divider()
    #############################################################################################
    st.markdown("### Let's jump to R!!!")

    with st.container(border=True):
            st.markdown("#### Set up")
            
            st.code("""
            library(tidyverse)
            library(data.table)
            library(lattice)
            library(gridExtra)
            library(rtracklayer)
            library(DiffBind)
            library(idr2d)
            library(patchwork)
            library(reshape2)
            library(GenomicRanges)
            library(dplyr)

            # set your working directory
            setwd("")
            """, language="r")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### IDR")

        st.markdown("We read in MACS peak calls and select coordinate columns as well as the p-value column, then run IDR")
        st.code("""
        #IDR
        #We read in MACS peak calls and select coordinate columns as well as the p-value column, then run IDR
        idr <- c('KO','KO2') %>%
            paste0('_H3K27ac.chr1_peaks.narrowPeak') %>%
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
        # Read prepared data
        load("lec10.RData")

        #PCA
        dba.plotPCA(dat, label = DBA_CONDITION)
        """)

        st.image("images/lec10.PCA.png")

        st.code("""
        #HCA
        dba.plotHeatmap(dat, correlations = F)
        """, language="r")

        st.image("images/lec10.corr.png")

    st.divider()

    ##########################################################################################

    with st.container(border=True):
        st.markdown("#### Heatmap")

        st.markdown("We’ll first take the KO sample’s H3K27ac enrichment surrounding promoters as an example")
        st.code("""
        #We’ll first take the KO sample’s H3K27ac enrichment surrounding promoters as an example
        # use only the matrix part, throw away the coordinates info
        dat.ex <- mat.ex[,-(1:6)] %>%
            arrange(rowMeans(.)) %>%
            as.matrix
        """)

    st.divider()
        
    ##########################################################################################

    with st.container(border=True):
        st.markdown("#### Including expression")

        st.markdown("Given the link between promoter acetylation and the cognate gene’s expression, we further incorporate RNA-seq data. Again, for the sake of time, it’s been imported and provided in lec10.RData")
        
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

        # partition genes into quantiles by expression
        # Divide the expression read counts by length of the gene
        ranks <- rank(genes.ol$counts/abs((end(genes.ol) - start(genes.ol))),
                    ties.method = "first")

        # Assign genes to each rank based on the expression
        genes.ol$cut <- as.numeric(cut(ranks, quantile(ranks, probs = seq(0, 1, 0.2)),
                            include.lowest = TRUE))
        # Do the same for promoters
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

        st.markdown("#### What are your conclusions from these plots ????")

    ##########################################################################################

    with st.container(border=True):
        st.markdown("#### Aggregate plots")

        st.markdown("Or alternatively take the column-wise average across said heatmap as a summary")

        st.code("""
        lapply(dat.ex, function(x) colMeans(x[,2:601])) %>%
            melt() %>%
            rownames_to_column("row") %>%
            mutate(row = (as.numeric(row) - 1) %% (n()/5) + 1) %>%
            ggplot(aes(x = row, y = value, color = L1)) +
            annotate("rect", xmin = 200, xmax = 400, ymin = -Inf,
                    ymax = Inf, alpha = 0.2, fill = "grey") +
            geom_line() +
            labs(y = "Enrichment", title = "H3K27ac in promoters",
                color = "Expression level") +
            scale_x_continuous(breaks = c(1, 200 * 1:3),
                                labels = c("-2kb", "Start", "End", "+2kb"),
                                expand = c(0, 0), name = "Position") +
            theme_bw() + theme(legend.position = "bottom")
        """)

        st.image("images/lec10.agg.plots.png")

        st.markdown("#### What are your conclusions from these plots ????")

    st.divider()

    ##############################################################################
    with st.container(border=True):
        st.markdown("#### Correlation")

        st.markdown("Or more directly we can just assess the relationship per-promoter")
        st.code("""
        ## Joining with `by = join_by(V4)`
        ols <- findOverlaps(pmts.ol, genes.ol, select = "first")
        genes.ol$len <- abs(end(genes.ol) - start(genes.ol)) / 1000
        fac <- sum(genes.ol$counts / genes.ol$len) / 1e6
        pmts.cut$tpm <- (genes.ol$counts[ols] / genes.ol$len[ols]) / fac
        dat.ex <- left_join(pmts.cut, mat.ex[,-c(1:3, 5, 6)]) %>%
            select(-c(V4, cut)) %>%
            mutate(rm = rowMeans(select(., 202:401))) %>%
            select(tpm, rm)
        smoothScatter(x = log2(dat.ex$tpm + 1), y = dat.ex$rm,
                    xlab = "Expression", ylab = "H3K27ac")
        """, language="r")
        st.image("images/lec10.correlation.png")

    st.divider()

    ##############################################################################
    with st.container(border=True):
        st.markdown("#### Your Turn!!!")

        st.markdown("""
        Look at the folder :blue[${data}/activity] and download the files. The students will pair up in teams of two
        and each team will analyze the distribution of a mark on lowly and highly expressed genes using `deepTools` commands
        computeMatrix and plotHeatmap.

        Draw conclusions and debate!
        """)
        
        st.markdown("""
        Some technical details:

        - Create just one image to analyze both conditions and both groups of genes (lowly and highly expressed genes)
        - Use a bin size of 250 bp
        - Specify a region upstream the gene of 4500 bp
        - Specify a region downstream the gene of 4500 bp
        - Scale the genic regions to 6000 bp using the `scale-regions` mode
        - Label the samples you use
        - Provide a title for the plot
        - Label the regions you are using (low and high expression gene groups)

        Add the following parameters to `computeMatrix` command:

            --missingDataAsZero --skipZeros -p "max/2"

        Add the following parameters to `plotHeatmap` command:

            --dpi 600 --perGroup --colorMap 'Reds'
        """)

        st.markdown("The plots should look like these: (This is a random example from the web)")
        st.image("images/lec10.example.png")

    st.divider()