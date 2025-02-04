import streamlit as st

col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("Transcription and RNA Sequencing", divider=True)

    with st.container(border=True):
        st.markdown("#### Take a look at the GENCODE GTF file")
        
        st.markdown("Check out the first few lines from the GTF file")
        st.code("""
        data=/project/def-sponsor00/hgen_share/lec5

        head ${data}/gencode.v36.annotation.chr1.gtf | column -t -s $'\t' | less -S
        """, language="bash")
    st.divider()
    #############################################################################################
    with st.container(border=True):
        st.markdown("#### Align `ex.fq` with STAR")
        
        st.markdown("Set up variables")
        st.code("""
        module load gcc/12.3  openmpi/4.1.5 salmon/1.10.2 star/2.7.11b subread/2.0.6

        BASE=/home/hgen_share
        PICARD_JAR=$data/../utils/picard.jar
        """, language="bash")

#       STAR --genomeDir $BASE/star/chr1 \\
        st.markdown("Align with [STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)")

        st.markdown("Build an index. :red[Don't run steps that involve the STAR aligner]")
        st.code("""
        mkdir -p ${data}/STAR_index

        STAR --runThreadN 2 \\
             --runMode genomeGenerate \\
             --genomeDir ${data}/STAR_index \\
             --genomeFastaFiles ${data}/chr1.fa \\
             --sjdbGTFfile$ {data}/gencode.v36.annotation.chr1.gtf &> log.indexing.txt
        """, language="bash")

        st.code("""
        mkdir -p /project/def-sponsor00/hgen_share/lec5/results

        STAR --genomeDir ${data}/STAR_index \\
             --runThreadN 4 \\
             --readFilesIn ${data}/ex.fq \\
             --outFileNamePrefix ${data}/results/ex_ \\
             --outSAMtype BAM SortedByCoordinate \\
             --outSAMunmapped Within \\
             --outSAMattributes Standard \\
             --outFilterScoreMinOverLread 0.3 \\
             --outFilterMatchNminOverLread 0.3 \\
             --twopassMode Basic
        """, language="bash")

        st.markdown("Take a look at the summary statistics with `less`")
        st.code("""
        less ${data}/results/ex_Log.final.out
        """, language="bash")
    st.divider()
    #############################################################################################
    with st.container(border=True):
        st.markdown("#### Quantify expression levels using featureCounts and Salmon")
        
        st.markdown("Run [featureCounts](https://subread.sourceforge.net/SubreadUsersGuide.pdf)")
        st.code("""
        featureCounts -T 4 -s 2 \\
                -a ${data}/gencode.v36.annotation.chr1.gtf \\
                -o ex_featurecounts.txt \\
                ${data}/results/ex_Aligned.sortedByCoord.out.bam
        """, language="bash")

        st.markdown("Run [Salmon](https://salmon.readthedocs.io/en/latest/)")
        st.code("""
        module load salmon/1.10.2
        
        salmon quant -i ${data}/transcriptome_index \\
                -l A \\
                -r ${data}/ex.fq --validateMappings \\
                -o salmon_quant
        """, language="bash")

        st.markdown("Download the output from the gene expression counters")

    st.divider()

    #############################################################################################

    with st.container(border=True):
        
        st.markdown("#### Set up")
        st.code("""
        library(ggplot2)
        library(dplyr)
        library(tximport)
        library(DESeq2)
        library(vsn)
        library(gplots)
        library(hexbin)
        library(eulerr)
        library(GenomicFeatures)
        library(rtracklayer)
        library(patchwork)
        library(tibble)
        """, language="r")

    st.divider()
    #############################################################################################
    with st.container(border=True):
        st.markdown("#### Import data")
        
        st.code("""
        if (!file.exists('gencode.v36.annotation.gff3.gz')) {
        download.file('ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.annotation.gff3.gz', 
                      'gencode.v36.annotation.gff3.gz')
        }

        g <- import.gff3('gencode.v36.annotation.gff3.gz')

        tx2gene <- g[g$type == 'transcript'] %>%
                {tibble(txid = .$transcript_id, gene = .$gene_id)} %>%
                mutate_all(function(x) sub('\\\..*', '', x))
        """, language="r")

        st.code("""
        head(g)
        """, language="r")

        st.image("images/lec5.granges.png")
        st.markdown("## ...")
        st.code("""
        head(tx2gene)
        """, language="r")
        st.image("images/lec5.tibble.png")

        st.divider()
        
        st.markdown("")
        st.code("""
        # read in featureCounts result
        fc <- read.table("ex_featurecounts.txt", header = 1)

        # rename counts column for simplicity
        names(fc)[7] <- 'count'

        # take a look at the table
        head(fc)
        """, language="r")

        st.code("""
        ggplot(fc, aes(x = count)) +
            geom_histogram(bins = 200) +
            scale_y_log10()
        """, language="r")
        st.image("images/lec5.fC.distribution.png")

    st.divider()
    #############################################################################################
    # Sys.setenv(PATH = paste("/Users/padilr1/opt/anaconda3/envs/r_env/bin", 
    #            Sys.getenv("PATH"), sep=":"))
    with st.container(border=True):
        st.markdown("#### Salmon")
        
        st.markdown("Load")
        st.code("""
        # read in raw data
        sf <- read.table('quant.sf', header = 1)

        # import salmon results
        sm <- tximport("quant.sf",
                type = "salmon",
                tx2gene = tx2gene,
                ignoreTxVersion=T)
        """, language="r")

        st.markdown("Per-transcript")
        st.code("""
        head(sf)
        """, language="r")
        st.image("images/lec5.head.sf.png")

        st.markdown("Per-gene")
        st.code("""
        str(sm)
        """, language="r")
        st.image("images/lec5.str.sm.png")
    st.divider()
    #############################################################################################
    with st.container(border=True):
        st.markdown("#### Comparison")
        
        st.markdown("Merge")
        st.code("""
        # merge by ensembl gene id
        both <- fc %>% mutate(Geneid = sub('\\\..*', '', Geneid)) %>% 
                merge(sm$counts, by.x = 'Geneid', by.y = 'row.names') %>%
                dplyr::rename(sm = V1, fc = count)
        """, language="r")

        st.markdown("Unexpressed")
        st.code("""
        fit <- lapply(both[c('sm', 'fc')], `==`, 0) %>%
                bind_cols() %>%
                euler()
                plot(fit, quantities = T)
        """, language="r")
        st.image("images/lec5.comp.fC.sm.png")

        st.markdown("Expressed")
        st.code("""
        both %>% 
            dplyr::filter(fc != 0 & sm != 0) %>% 
            ggplot(aes(x = sm, y = fc)) +
            geom_point() +
            scale_x_log10() + 
            scale_y_log10() + 
            labs(x = 'Salmon', y = 'featureCounts')
        """, language="r")
        st.image("images/lec5.comp.fC.sm2.png")

        # st.markdown("Create DESeq object")
        # st.code("""
        # mat <- as.matrix(both[,c('fc','sm')])

        # cdat <- data.frame(kind = colnames(mat)) %>% 
        #                 `rownames<-`(colnames(mat))

        # dds <- DESeqDataSetFromMatrix(countData = round(mat),
        #                               colData = cdat,
        #                               design = ~kind)
        # dds
        # """, language="r")
        # st.image("images/lec5.deseq.obj.png")
    st.divider()
    #############################################################################################
    with st.container(border=True):
        st.markdown("#### Full dataset")
        
        st.markdown("Load")
        st.code("""
        # load full dataset
        load('lec5.rda')

        # take a look at "data"
        head(data)
        """, language="r")
        st.image("images/lec5.full.dataset.png")

        st.markdown("""We can compare different replicates and conditions of RNA-Seq experiments""")
        st.code("""
        keep <- rowSums(counts(dds)) >= 10
        dds.filt <- dds[keep,]
        """, language="r")

        # st.image("images/lec5.var.stab.png")

        st.markdown("#### Dimension reduction")
        st.markdown("Perform a variance stabilizing transformation")
        st.code("""
        vsd <- vst(dds.filt)
        plotPCA(vsd)
        """, language="r")
        st.image("images/lec5.PCA.png")
    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Hierarchical clustering")
                
        st.markdown("There are too many genes to perform hierarchical clustering quickly, so a subset of data will be looked at instead.")
        st.code("""
        # how can you subset the data?
        select <- assay(vsd) %>%
        apply(1, sd) %>%
        order(decreasing = TRUE) %>%
        .[1:100]
        heatmap(assay(vsd)[select,], margins = c(10, 5))
        """, language="r")
        st.image("images/lec5.heatmap.png")

        st.markdown("#### Spot-check")
        st.markdown("Pick a gene, does the expression level seem different between conditions?")
        st.code("""
        plotCounts(dds.filt, gene = select[1], intgroup = "condition")
        """, language="r")
        st.image("images/lec5.ind.gene.png")

    st.divider()
    #############################################################################################

    # with st.container(border=True):
    #     st.markdown("#### Your turn")
        
    #     st.markdown("Load the featureCount dataset ex2_featureCounts.txt.")
    #     st.markdown("**Load**")
    #     st.code("""
    #     fc <- read.table("~/Documents/HGEN_663/extra/lec5/ex2_featureCounts.txt",header=1)
    #     """, language="r")

    #     st.markdown("**Create DESeq object**")
    #     st.markdown("""Remove unnecessary columns. An example of how to generate
    #     a metadata object is indicated below. Within the metadata dataframe, 
    #     create another column called condition, grouping together the parental (PA),
    #     NSD1 knockout (NSD1KO) and NSD1/2 double knockouts (NSD12DKO) samples. 
    #     Use ~condition in the design parameter when creating your DESeq object.""")
    #     st.code("""
    #     df <- fc %>% 
    #     dplyr::select(-c(2:6)) %>% 
    #             tibble::column_to_rownames("Geneid")

    #     mat <- df %>% as.matrix()

    #     metadata <- data.frame(kind=colnames(mat))

    #     metadata$condition <- "PA"

    #     metadata[3:4,"condition"] <- "NSD1KO"
    #     metadata[5:6,"condition"] <- "NSD12DKO"

    #     metadata <- tibble::column_to_rownames(metadata,"kind")

    #     dds <- DESeqDataSetFromMatrix(countData = mat,
    #                                   colData = metadata,
    #                                   design = ~condition)
    #     """, language="r")

    #     st.markdown("**Apply vst transformation**")
    #     st.code("""
    #     vst <- varianceStabilizingTransformation(dds)
    #     """, language="r")

    #     st.markdown("**PCA**")
    #     st.markdown("Plot only the top 500 most variable genes.")
    #     st.code("""
    #     # subset data
    #     select <- assay(vst) %>%
    #             apply(1, sd) %>%
    #             order(decreasing = TRUE) %>%
    #             .[1:500]
    #     vst_filt <- vst[select,]
        
    #     plotPCA(vst_filt)
    #     """, language="r")
    #     st.image("images/lec5.PCA2.png")

    #     st.markdown("**Hierarchical clustering**")
    #     st.markdown("Do the same as above for hierarchical clustering.")
    #     st.code("""
    #     heatmap(assay(vst_filt), margins = c(10, 5))
    #     """, language="r")
    #     st.image("images/lec5.heatmap2.png")

    #     st.markdown("**Spotcheck**")
    #     st.markdown("""Google the ensembl id for NSD1 (in humans). Plot the expression level.
    #     Does the expression level seem different between the knockout and parental samples?""")
    #     st.code("""
    #     # NSD1 = ENSG00000165671
    #     # NSD2 = ENSG00000109685
    #     plotCounts(dds,gene = "ENSG00000165671",intgroup = "condition")
    #     """, language="r")
    #     st.image("images/lec5.ind.gene2.png")
