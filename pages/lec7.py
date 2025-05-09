import streamlit as st

col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("More Advanced RNA-seq analysis", divider=True)

    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Copy the files for today's class from `data`")
        st.code("data=/project/60006/hgen_share/lec7", language="bash")
        # st.markdown("Set up directory and copy over files")
            
    st.divider()
    #############################################################################################

    # with st.container(border=True):
    #     st.markdown("#### Take a look at the outputs of STAR-FUSION")
        
    #     st.markdown("Check out the first few lines from the abridged results file")
    #     st.code("""
    #     column -t ${data}/fusion/kd3.fusion_predictions.abridged.tsv | less -S
    #     """, language="bash")
    # st.divider()
    #############################################################################################
    with st.container(border=True):
        st.markdown("#### Assemble transcripts using StringTie")
        
        st.markdown("Set up variables")
        st.code("""
        module load StdEnv/2023 stringtie/2.2.3

        awk '$1=="chr6"' ${data}/gencode.gtf > gencode.chr6.gtf
        """, language="bash")
        st.markdown("Apply `stringtie` without and without reference annotations")
        st.code("""
        stringtie --rf -o denovo/ct1.gtf ${data}/ct1.bam
        stringtie --rf -G gencode.chr6.gtf -o guided/ct1.gtf ${data}/ct1.bam
        """, language="bash")
        st.markdown("Combine sample-specific outputs")
        st.code("""
        stringtie --merge -o denovo.gtf denovo/*.gtf
        stringtie --merge -o guided.gtf guided/*.gtf
        """, language="bash")
        
        st.markdown("Compare with reference annotations using `gffcompare`")
        st.code("""
        module load gffcompare/0.12.6
        gffcompare -r gencode.chr6.gtf -o denovo denovo.gtf
        gffcompare -r gencode.chr6.gtf -o guided guided.gtf
        """, language="bash")
        st.markdown("Check out the summary")
        st.code("""
        less denovo.stats
        less guided.stats
        """, language="bash")
        st.markdown("How many novel transcripts did we find?")
        st.code("""
        awk '$3=="transcript"' guided.annotated.gtf | grep 'class_code "u"' | wc -l 
        """, language="bash")
    st.divider()
    #############################################################################################
    with st.container(border=True):
        st.markdown("#### Identify differential splicing events through LeafCutter")
        
        st.markdown("Set up environment")

        st.markdown("Parse out junctions with `regtools`")
        st.code("""
        regtools=/project/60006/hgen_share/utils/regtools/bin/regtools

        for b in ${data}/bam/*.bam; do $regtools junctions extract -a 8 -m 50 -M 500000 -s RF $b -o $(basename $b .bam).junc; done
        """, language="bash")

        st.markdown("Note sample information for downstream use")
        # st.code("""
        # ls -A1 junc/*junc > juncfiles
        # sed 's/.junc//' juncfiles \\
        # | awk -v OFS='\t' '{print $1, substr($1,0,2) == "ct" ? "control" : "knockdown"}' \\
        # > groupfiles
        # """, language="bash")
        # st.code("cp $data/juncfiles /project/60006/$USER", language="bash")

        st.markdown("### Download files to your computer")
        st.code("scp <user>@hgen633.calculquebec.cloud:/project/60006/hgen_share/lec7/download/* .", language="bash" )

        st.markdown("Run [leafCutter](https://davidaknowles.github.io/leafcutter/articles/Usage.html) in your system")
        st.code("""
        python leafcutter_cluster_regtools.py -j juncfiles -m 50 -l 500000
        """, language="bash")
        
        # st.markdown("Run next lines in your local system")
        # st.markdown("First download the data")


        st.code("""
        Rscript leafcutter_ds.R \\
                --num_threads 1 \\
                --exon_file exons.txt.gz \\
                --min_samples_per_intron 3 \\
                leafcutter_perind_numers.counts.gz \\
                groupfiles
        """, language="bash")

        st.markdown("Pack up results for visualization")
        st.code("""
        Rscript prepare_results.R \\
                --meta_data_file groupfiles \\
                --code KD \\
                leafcutter_perind_numers.counts.gz \\
                leafcutter_ds_cluster_significance.txt \\
                leafcutter_ds_effect_sizes.txt \\
                hg38 \\
                -o viz.rda
        """, language="bash")
        st.markdown("Download `viz.rda` to your local computer")
    st.divider()
    #############################################################################################
    with st.container(border=True):
        st.markdown("#### Alternatively, we can run rMATS :red[Don't run this section!]")
        
        st.markdown("Set up directory and prepare files for input into [rMATS](https://github.com/Xinglab/rmats-turbo/blob/v4.1.2/README.md)")
        st.code("""
        mkdir output_rmats
        mkdir tmp_rmats

        echo "${data}/bam/ct1.bam,${data}/bam/ct2.bam,${data}/bam/ct3.bam" > control_list.txt
        echo "${data}/bam/kd1.bam,${data}/bam/kd2.bam,${data}/bam/kd3.bam" > kd_list.txt
        """, language="bash")


        st.markdown("Run rMATS")
        st.code("""
        rmats.py --b1 control_list.txt --b2 kd_list.txt \\
                --gtf ${data}/ref-transcripts.gtf \\
                --od output_rmats -t paired --libType fr-firststrand \\
                --readLength 150 --nthread 1 --tmp tmp_rmats --variable-read-length
        """, language="bash")

        st.divider()
        #############################################################################################

    with st.container(border=True):
        st.markdown("#### Let's jump to R!!!")
        st.divider()
        st.markdown("#### Part 1")
        
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
        library(leafviz)
        # set your working directory
        setwd("~/Documents/HGEN_663/extra/lec7")
        """, language="r")
        st.markdown("#### Install package")
        st.code("""
        # install leafviz
        install.packages("remotes")
        remotes::install_github("jackhump/leafviz")
        library(leafviz)
        """, language="r")
        st.markdown("#### Run Leafviz")
        st.markdown("You can find more details on Leafviz [here](https://github.com/jackhump/leafviz)")
        st.code("""
        leafviz("path/to/viz.rda")
        """, language="r")
    st.divider()
    #############################################################################################
    with st.container(border=True):
        st.markdown("##### Part 2")
        
        st.markdown("#### Import results from rMATS")

        st.code(
        """
# Before starting, make sure that following packages are loaded
library("clusterProfiler")
library("AnnotationDbi")
library("org.Hs.eg.db")

# Step 1: Change the filepath to where the downloaded "splicing_analysis_demo" file is stored on your device
filePath <- "path/to/folder"
        
# Step 2: Set the cutoff criteria for filtering the RMATS results
FDR_cutoff <- 0.05
inclusion_diff_cutoff <- 0.1

# Step 3: Load the RMATS output file titled "RI.MATS.JC.txt"
event_file <- read.delim() ### Complete this line

head(event_file) # see first few lines of the event file

# Step 4: Filter for events with an FDR below the set cutoff and an inclusion level difference about cutoff
event_file_filtered <- event_file[] ### Complete this line

# Step 5: Classify if each event is more included (ie. occurs more) in the control or case condition
event_file_filtered_control <- event_file_filtered[event_file_filtered$IncLevelDifference > inclusion_diff_cutoff,]
event_file_filtered_case <- event_file_filtered[event_file_filtered$IncLevelDifference < -inclusion_diff_cutoff,]


cat("Number of events prior to filtering: ", nrow(event_file),
"\\nNumber of events after filtering on FDR: ", nrow(event_file_filtered),
"\\nNumber of events with more intron retention in control: ", nrow(event_file_filtered_control),
"\\nNumber of events in more intron retention in case: ", nrow(event_file_filtered_case))

# Step 6: Use all events as a background in over-representation analysis
background <- unique(event_file$geneSymbol)
background <- mapIds(org.Hs.eg.db, keys = background, column = "ENTREZID", keytype = "SYMBOL") # converts gene symbols to Entrez ID, required for enrichGO function
background <- background[!is.na(background)] # removes na values (no mapped Entrez ID)

# Step 7: Convert our subset of significant events into a gene list
gene_list <- unique(event_file_filtered$geneSymbol)
gene_list <- mapIds(org.Hs.eg.db, keys = gene_list, column = "ENTREZID", keytype = "SYMBOL") # converts gene symbols to Entrez ID, required for enrichGO function
gene_list <- gene_list[!is.na(names(gene_list))] # removes na values (no mapped Entrez ID)

# Step 8: Run GO results with our background 
GO_results_with_background <- enrichGO(
        
)

# Step 9: Run GO results without our background
GO_results_NO_background <- enrichGO(
        
)

# Try changing the value of "showCategory", what happens?
barplot(GO_results_with_background, showCategory = 10, title = "Results Using Gene Background") + 
    barplot(GO_results_NO_background, showCategory = 10, title = "Results Without Gene Background") # shows the results of enrichGO as barplots

""", language="r")

        # st.markdown("Import the rest of the results from rMATS (ie. those ending with ‘.JC.’, which are junction counts).")
        # st.code("""
        # se <- read.table("path/to/rmats_output_SNRPB/SE.MATS.JC.txt",sep = '\t', header = TRUE, stringsAsFactors = FALSE)
        # ri <- read.table("path/to/rmats_output_SNRPB/RI.MATS.JC.txt",sep = '\t', header = TRUE, stringsAsFactors = FALSE)
        # """, language="r")
        # st.markdown("#### Differential Splicing Events")
        # st.markdown("Which type are the most abundant differentially spliced events (DSEs)? Is there a tendancy for the knockdown samples towards one of the DSEs? (for example exon skipping?)")
        # st.code("""
        # se_dse <- se %>%
        #     filter(FDR < 0.05)
        #     nrow(se_dse[se_dse$IncLevelDifference < 0,]) #differential skipping in control
        #     nrow(se_dse[se_dse$IncLevelDifference > 0,]) #differential skipping in knockdown
        
        # ri_dse <- ri %>%
        #     filter(FDR < 0.05)
        #     nrow(ri_dse[ri_dse$IncLevelDifference > 0,]) #differential retained introns in control
        #     nrow(ri_dse[ri_dse$IncLevelDifference < 0,]) #differential retained introns in knockdown
        # """, language="r")

        # st.markdown("Plot your results!")
        # st.code("""
        # df <- data.frame(SE = c(106,354),
        #         RI = c(26,317),
        #         samp = c("Control","Knockdown")) %>%
        #         pivot_longer(cols = c("SE","RI"),names_to = "DSE",values_to = "count")
        
        # df$samp <- fct_rev(df$samp)

        # df %>%
        # ggplot(aes(x=DSE,y=count,fill=samp)) +
        # geom_bar(position="stack",stat="identity") +
        # coord_flip() +
        # scale_fill_manual(values=c("purple", "red")) +
        # theme_minimal()
        # """, language="r")
                
