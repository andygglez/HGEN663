import streamlit as st

col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("Analysis of epigenomic data", divider=True)

    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Call peaks, run IDR and then pool samples")
            
            st.markdown("Look at the files in the data folder")
            st.code("""
            data=/project/60006/hgen_share/lec11
            ll ${data}
            """, language="bash")

            st.markdown("Call peaks using `MACS3`. Call peaks for the other samples: mut_rep2.bam,rescue_rep1.bam,rescue_rep2.bam")
            st.code("""
            macs3 callpeak -g hs -f BAMPE --extsize 200 \\
                        --shift -100 --nomodel -B -p 1e-3 \\
                        --keep-dup all -t ${data}/mut_rep1.bam \\
                        -n mut_rep1 --outdir .
        
            macs3 callpeak -g hs -f BAMPE --extsize 200 \\
                        --shift -100 --nomodel -B -p 1e-3 \\
                        --keep-dup all -t ${data}/mut_rep2.bam \\
                        -n mut_rep2 --outdir .
            """, language="bash")

            st.markdown("Run IDR for both conditions")
            st.code("""
            idr --samples mut_rep1_peaks.narrowPeak mut_rep2_peaks.narrowPeak \\
                --input-file-type narrowPeak \\
                --rank p.value \\
                --output-file mut-idr.peaks.txt \\
                --output-file-type bed \\
                --idr-threshold 0.05 \\
                --plot \\
                --log-output-file mut.idr.log
            """, language="bash")

            st.markdown("Pool together peaks. Make a .saf file to use with `featureCounts`")
            st.code("""
            module load bedtools

            cat *-idr.peaks.txt | sort -k1,1 -k2,2n | bedtools merge | tee merged.bed | wc -l

            awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' merged.bed > merged.saf
            """, language="bash")

            st.markdown("Run featureCounts")
            st.code("""
            module load subread

            featureCounts -a merged.saf -F SAF -o xie2022.fc.txt \\
                        -p ${data}/mut_rep1.bam ${data}/mut_rep2.bam \\
                        ${data}/rescue_rep1.bam ${data}/rescue_rep2.bam
            """, language="bash")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Assess profile of ATAC-seq signal on different types of reference peaks")
            
            st.markdown("Aggregate signal centered on ATAC-seq peaks using computeMatrix")
            st.code("""
            module load apptainer
            singularity exec ${data}/deeptools.sif \\
                        computeMatrix reference-point -R $data/xie2022.open_peaks.bed \\
                        -S ${data}/mut_rep1.cpm.bw ${data}/mut_rep2.cpm.bw ${data}/rescue_rep1.cpm.bw ${data}/rescue_rep2.cpm.bw \\
                        --referencePoint center \\
                        --binSize 50 -a 2500 -b 2500 \\
                        --missingDataAsZero --skipZeros \\
                        -o ATAC.mat.gz
            

            singularity exec ${data}/deeptools.sif \\
                        plotHeatmap -m ATAC.mat.gz -o heatmap.png \\
                        --colorMap Reds --verbose -x "ATAC" -y "" \\
                        --startLabel "" --endLabel "" --samplesLabel mut_rep1 mut_rep2 rescue_rep1 rescue_rep2

            """, language="bash")

    st.divider()
    #############################################################################################

    st.markdown("#### Let's jump to R!!!")

    with st.container(border=True):
        st.markdown("#### Set up")
        
        st.code("""
        library(data.table)
        library(tidyverse)
        library(lattice)
        library(gridExtra)
        library(rtracklayer)
        library(idr2d)
        library(patchwork)
        library(reactable)
        library(ggplot2)
        library(DESeq2)
        library(TxDb.Hsapiens.UCSC.hg19.knownGene)
        library(org.Hs.eg.db)
        library(plotly)
        library(ChIPseeker)
        library(annotate)
        library(gprofiler2)
        library(clusterProfiler)
        options(ChIPseeker.ignore_1st_exon = TRUE)
        options(ChIPseeker.ignore_1st_intron = TRUE)
        options(ChIPseeker.ignore_downstream = TRUE)
        options(ChIPseeker.ignore_promoter_subcategory = TRUE)

        # set your working directory
        setwd("")
            """, language="r")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Import")
        
        st.code("""
        # Read in and process the feature counts from /home/hgen_share/lec11 as previously shown
        mat <- read.table()
        
        # Do the necessary transformations to the columns
        """, language="r")

        st.code("""
        load("res.RData")
        load("vst.RData")
        load("dds.RData")
        """)
    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### DESeq2")
        
        st.code("""
        mdat <- data.frame(cond=c("mut","mut","rescue","rescue"))

        dds <- DESeqDataSetFromMatrix(countData = ?,
                                colData = ?,
                                design = ?)
        

        dds <- DESeq(dds)
        # get results dataframe
        res <- results(dds,contrast = c("cond","mut","rescue")) %>%
                data.frame() %>%
                na.omit() %>%
                rownames_to_column("Geneid")

        """, language="r")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### PCA")
        
        st.markdown("The code to run a PCA is given. However, you will first need to run variance stabilizing transformation as previously shown.")
        st.code("""
        vsd <- vst(dds)
        # subset data - get the top 500 most variable peaks
        rv <- rowVars(assay(vsd))
        select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
        # get pc scores
        pc <- prcomp(t(assay(vsd)[select,]))
        condition <- mdat$cond
        scores <- data.frame(pc$x, condition)
        percentage <- round(pc$sdev / sum(pc$sdev) * 100, 2)
        percentage <- paste( colnames(scores), "(", paste( as.character(percentage), "%", ")", sep="") )
        # plot
        p <- plot_ly(scores,x=scores$PC1,y=scores$PC2,text=rownames(scores),mode="markers",color=factor(condition),marker=list(size=11),colors = "Dark2")
        p <- layout(p,title="",   xaxis = list(title = percentage[1]),
        yaxis = list(title = percentage[2]))
        p
        """, language="r")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Volcano plot")
        
        st.code("""
        volc <- function(r, x, y, ttl) {
        d <- as.data.frame(r) %>%
        dplyr::rename(x = !!x, y = !!y) %>%
        mutate(kind = case_when(abs(x) > 1 & y < .05 ~ '|log2FC|>1,FDR<0.05',
                                abs(x) < 1 & y < .05 ~ '|log2FC|<1,FDR<0.05',
                                y > .05 ~ 'FDR > 0.05'),
                y = -log10(y))
        ct <- na.omit(d) %>% 
        dplyr::filter(kind == '|log2FC|>1,FDR<0.05') %>%
        mutate(up = x > 0) %>%
        dplyr::count(up) %>%
        mutate(x = ifelse(up, Inf, -Inf),
                y = Inf,
                h = as.numeric(up))
        ggplot(d, aes(x, y, color = kind)) +
        geom_vline(xintercept = c(-1, 1), linetype = 'dashed') +
        geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
        geom_point(alpha = .5) +
        geom_label(aes(x = x, y = y, label = n, hjust = h),
                vjust = 1, data = ct, inherit.aes = F) +
        scale_y_continuous() +
        scale_color_manual(values = c('royalblue', 'red2', 'grey30')) +
        labs(x = x, y = "-log10(FDR)", title = ttl) +
        theme_bw()
        }
        volc(res,'log2FoldChange','padj','mutant (H3K36me3 -) versus rescue (H36K36me3 +)')
        """, language="r")
        st.image("images/lec11.volcano.png")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Peaks processing")
            
        st.code("""
        # filter 'res' object for open and closed ATAC-seq peaks as well as separate the genomic intervals for each peak into 'chr', 'start' and 'end'
        open_peaks <- res[res$log2FoldChange > 1 & res$padj < 0.05,] %>%
        tidyr::separate(Geneid,c("chr","start","end"))
        closed_peaks <- res[res$log2FoldChange < -1 & res$padj < 0.05,] %>%
        tidyr::separate(Geneid,c("chr","start","end"))

        # get rid of NAs and non-numeric values
        value <-  as.numeric( as.character(open_peaks$start) ) # get the numbers
        noninteger <- value %% 1 != 0   # see if there's a fractional part
        noninteger <- noninteger | is.na(noninteger)  # get rid of the NA's 
        open_peaks <- open_peaks[!noninteger,]

        # make granges object
        open_peaks_gr <- makeGRangesFromDataFrame(open_peaks)

        # Do the same for closed peaks
        # get rid of NAs and non-numeric values
        value <-  as.numeric( as.character(closed_peaks$start) ) # get the numbers
        noninteger <- value %% 1 != 0   # see if there's a fractional part
        noninteger <- noninteger | is.na(noninteger)  # get rid of the NA's 
        closed_peaks <- closed_peaks[!noninteger,]
        # make granges object
        closed_peaks_gr <- makeGRangesFromDataFrame(closed_peaks)

        # Lastly,we can do the same for the original pooled peaks
        pooled <- res %>%
        tidyr::separate(Geneid,c("chr","start","end"))
        # get rid of NAs and non-numeric values
        value <-  as.numeric( as.character(pooled$start) ) # get the numbers
        noninteger <- value %% 1 != 0   # see if there's a fractional part
        noninteger <- noninteger | is.na(noninteger)  # get rid of the NA's 
        pooled <- pooled[!noninteger,]
        # make granges object
        pooled_peaks_gr <- makeGRangesFromDataFrame(pooled)
        """, language="r")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Distribution")
        
        st.markdown("We can use ChIPSeeker to annotate our peaks. Below, I’ve shown how to annotate the open peaks for you. Annotate the closed and pooled peaks as well as the open peaks. Combine all three types of peaks and compare their distribution.")
        st.code("""
        # load known genes for hg19
        txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

        # annotate the peaks
        ann <- annotatePeak(open_peaks_gr,TxDb = txdb,tssRegion = c(-3000, 3000))
        """, language="r")

        st.code("""
        # plot the peak distribution
        plotAnnoBar(ann)
        """, language="r")
        st.image("images/lec11.plotAnnoBar.png")

        st.code("""
        # plot the percentage per each annotation, here it's shown for open peaks
        plotAnnoPie(ann)
        """, language="r")
        st.image("images/lec11.annoPie.png")

        st.code("""
        # plot distribution of TF-binding loci relative to TSS. 
        # Since we observed almost most of the open peaks to be in distal intergenic regions, 
        # it’s not surprising that a large percentage of sites are far away from the TSS
        plotDistToTSS(ann)
        """)
        st.image("images/lec11.distToTss.png")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Pathway enrichment analysis")
        
        st.markdown("""We can investigate whether genes nearest each peak are enriched 
        for a pathway. Considering this mutation is associated with kidney cancer, 
        some of these enriched pathways makes sense, such as ‘renin secretion’ and 
        ‘calcium signaling pathway’. Perform KEGG enrichment analysis for the closed peaks.""")
        st.code("""
        load("ann.RData")

        # you will need a background list, which will be all the genes nearest the peaks from your pooled peak list
        pooled <- ann$pooled_peaks %>% as.data.frame()
        bg <- pooled$geneId %>%
        as.character() %>%
        unique()
        # genes for open peaks
        open_peak_df <- ann$open_peaks %>%
        as.data.frame() 
        entrezids <- open_peak_df$geneId %>%
        as.character() %>%
        unique()
        #
        ekegg <- enrichKEGG(gene = entrezids,
                        organism = 'hsa',
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "fdr",
                        universe = bg)
        dotplot(ekegg)
        """, language="r")
        st.image("images/lec11.dotplot.png")

    st.divider()
    #############################################################################################