import streamlit as st

col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("Analysis of epigenomic data", divider=True)

    #############################################################################################

    with st.container(border=True):
            st.markdown("#### For today's class, only copy over files that are not bam files `from/home/hgen_share/lec11`")
            
            st.markdown("Set up environment and copy over files")
            st.code("""
            export PATH="/home/hgen_share/Anaconda/bin:$PATH"
            export PATH="/home/hgen_share/Anaconda/envs/chip_seq_V2/bin:$PATH"
            """, language="bash")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Call peaks, run IDR and then pool samples")
            
            st.markdown("Call peaks using `MACS3`. Call peaks for the other samples: mut_rep2.bam,rescue_rep1.bam,rescue_rep2.bam")
            st.code("""
            macs3 callpeak -g hs -f BAMPE --extsize 200 \\
                        --shift -100 --nomodel -B -p 1e-3 \\
                        --keep-dup all -t /home/hgen_share/lec11/mut_rep1.bam \\
                        -n mut_rep1 --outdir .
            """, language="bash")

            st.markdown("Run IDR for both conditions")
            st.code("""
            idr --samples mut_rep1_peaks.narrowPeak mut_rep2_peaks.narrowPeak \\
                --input-file-type narrowPeak \\
                --rank p.value \\
                --output-file mut-idr \\
                --output-file-type bed \\
                --idr-threshold 0.05 \\
                --plot \\
                --log-output-file mut.idr.log
            """, language="bash")

            st.markdown("Pool together peaks. Make a .saf file to use with `featureCounts`")
            st.code("""
            module load StdEnv/2020 bedtools/2.30.0

            cat *-idr | sort -k1,1 -k2,2n | bedtools merge | tee merged.bed | wc -l
            bedtools intersect -a merged.bed -b blacklist.hg19.bed -v > merged.filtered.bed

            awk 'OFS="\t" {print $1"."$2"."$3, $1, $2, $3, "."}' merged.filtered.bed > merged.saf
            """, language="bash")

            st.markdown("Run featureCounts")
            st.code("""
            module load StdEnv/2020 gcc/9.3.0 subread/2.0.1
            featureCounts -a merged.saf -F SAF -o xie2022.fc.txt \\
                        -p /home/hgen_share/lec11/mut_rep1.bam /home/hgen_share/lec11/mut_rep2.bam \\
                        /home/hgen_share/lec11/rescue_rep1.bam /home/hgen_share/lec11/rescue_rep2.bam
            """)

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Assess profile of ATAC-seq signal on different types of reference peaks")
            
            st.markdown("Aggregate signal centered on ATAC-seq peaks using computeMatrix")
            st.code("""
            computeMatrix reference-point -R open_peaks.bed \\
                -S mut_rep1.cpm.bw mut_rep2.cpm.bw rescue_rep1.cpm.bw rescue_rep2.cpm.bw \\
                -bl blacklist.hg19.bed \\
                --referencePoint center \\
                --binSize 50 -a 2500 -b 2500 \\
                --missingDataAsZero --skipZeros \\
                --samplesLabel mut_1 mut_2 rescue_1 rescue_2 \\
                -o ATAC.mat.gz
            """, language="bash")

            st.markdown("Due to memory/speed limitations, you can cancel it with `Ctrl + C` Instead, the compressed matrix from this analysis has been provided for you in `/home/hgen_share/lec11`")
            st.markdown("""
            Choosing the appropriate reference regions influences how the signal behaves. 
            Plot the output based on different reference regions: `union_peaks.mat.gz`, `open_peaks.mat.gz`, 
            `closed_peaks.mat.gz`, `mut_peaks.mat.gz`, `rescue_peaks.mat.gz`. 
            Visualize the output through `plotHeatmap` and/or `plotProfile`. 
            How do the plots differ using the different reference regions?
            """)
            st.code("""
            plotHeatmap -m ATAC.mat.gz -o heatmap.png --dpi 600 \\
                        --colorMap Reds --verbose -x "ATAC" -y "" \\
                        --startLabel "" --endLabel "" --regionsLabel ""

            plotProfile -m ATAC.mat.gz -out aggregate_plot.png \\
                        --perGroup --refPointLabel "ATAC" -T "ATACseq" \\
                        --dpi 600 --regionsLabel ""
            """, language="bash")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Run HOMER on the open and closed ATAC-seq peaks")
            
            st.markdown("Perform motif enrichment analysis using findMotifsGenome.pl")
            st.code("""
            findMotifsGenome.pl open_peaks.bed hg19 HOMER.open_peaks -size given -mask -bg union_peaks.bed
            """, language="bash")

            st.markdown("You can cancel it with `Ctrl + C`")

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
            setwd("~/Documents/HGEN_663/extra/lec11")
            """, language="r")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Import")
            
            st.code("""
            # Read in and process the feature counts from /home/hgen_share/lec11 as previously shown
            mat <- read.table() %>%
            column_to_rownames()
            """, language="r")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### DESeq2")
            
            st.code("""
            mdat <- data.frame()
            mdat$cond <- c("mut","mut","rescue","rescue") 
            mdat <- mdat %>%
            column_to_rownames("kind")
            dds <- DESeqDataSetFromMatrix(countData = ?,
                                        colData = ?,
                                        design = ?)
            dds <- DESeq(?)
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
            vst <- ?
            # subset data - get the top 500 most variable peaks
            rv <- rowVars(assay(vst))
            select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
            # get pc scores
            pc <- prcomp(t(assay(vst)[select,]))
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

    with st.container(border=True):
            st.markdown("#### Extra")
            st.markdown("#### Split ATAC-seq peaks according to their annotation: promoter, intron or intergenic")
            
            st.code("""
            promoter <- ann$open_peaks %>% as.data.frame() %>%
            dplyr::filter(str_detect(annotation, "Promoter")) %>%
            dplyr::select(1:3)
            intergenic <- ann$open_peaks %>% as.data.frame() %>%
            dplyr::filter(str_detect(annotation, "Distal Intergenic")) %>%
            dplyr::select(1:3)
            intron <- ann$open_peaks %>% as.data.frame() %>%
            dplyr::filter(str_detect(annotation, "Intron")) %>%
            dplyr::select(1:3)
            export.bed(object = promoter,con="open_peaks.promoter.bed")
            export.bed(object = intergenic,con="open_peaks.intergenic.bed")
            export.bed(object = intron,con="open_peaks.intron.bed")
            """, language="r")

            st.markdown("#### Convert peaks into sequence strings. Essentially, going from .bed files to .fasta files")
            st.code("""
            # load required libraries
            library(BSgenome.Hsapiens.UCSC.hg19)
            library(seqinr)

            # load genome
            Hsapiens <- BSgenome.Hsapiens.UCSC.hg19

            # extract sequence from 'granges' object
            seq <- getSeq(Hsapiens,open_peaks_gr)

            # you need to create an 'id' for each sequence, here we are naming them peak1, peak2, etc.
            names <- paste0("peak",1:8228)

            # export fasta file
            write.fasta(seq,names = names,file.out = "open_peaks.fasta")
            """, language="r")

    st.divider()
    #############################################################################################
