import streamlit as st

col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("More Advanced RNA-seq analysis", divider=True)



    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Copy the files for today's class from `/home/hgen_share/lec7`")
        st.markdown("Set up directory and copy over files")
            
    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Take a look at the outputs of STAR-FUSION")
            
        st.markdown("Check out the first few lines from the abridged results file")
        st.code("""
        column -t fusion/kd3.fusion_predictions.abridged.tsv | less -S
        """, language="bash")
    st.divider()
    #############################################################################################
    with st.container(border=True):
        st.markdown("#### Assemble transcripts using StringTie")
        
        st.markdown("Set up variables")
        st.code("""
        module load StdEnv/2020 stringtie/2.1.5
        awk '$1=="chr6"' gencode.gtf > gencode.chr6.gtf
        """, language="bash")
        st.markdown("Apply `stringtie` without and without reference annotations")
        st.code("""
        stringtie --rf -o guided/ct1.gtf ct1.bam
        stringtie --rf -G gencode.chr6.gtf -o denovo/ct1.gtf ct1.bam
        """, language="bash")
        st.markdown("Combine sample-specific outputs")
        st.code("""
        stringtie --merge -o denovo.gtf denovo/*.gtf
        stringtie --merge -o guided.gtf guided/*.gtf
        """, language="bash")
        
        st.markdown("Compare with reference annotations using `gffcompare`")
        st.code("""
        module load nixpkgs/16.09 gffcompare/0.11.6
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
        awk '$3=="transcript"' guided.annotated.gtf| grep 'class_code "u"' | wc -l 
        """, language="bash")
    st.divider()
    #############################################################################################
    with st.container(border=True):
        st.markdown("#### Identify differential splicing events through LeafCutter")
        
        st.markdown("Set up environment")
        st.code("""
        export PATH="/home/hgen_share/Anaconda/bin:$PATH"
        export PATH="/home/hgen_share/Anaconda/envs/r_env/bin:$PATH"
        """, language="bash")
        st.markdown("Parse out junctions with `regtools`")
        st.code("""
        for b in *.bam; do regtools junctions extract -a 8 -m 50 -M 500000 -s RF $b -o $b.junc; done
        """, language="bash")
        st.markdown("Note sample information for downstream use")
        st.code("""
        ls -A1 *junc > juncfiles
        sed 's/.junc//' juncfiles \\
        | awk -v OFS='\t' '{print $1, substr($1,0,2) == "ct" ? "control" : "knockdown"}' \\
        > groupfiles
        """, language="bash")
        st.markdown("Run [leafCutter](https://davidaknowles.github.io/leafcutter/articles/Usage.html)")
        st.code("""
        export PATH="/home/hgen_share/Anaconda/bin:$PATH"
        export PATH="/home/hgen_share/Anaconda/envs/r_env_V3/bin:$PATH"
        python /home/hgen_share/utils/leafcutter/clustering/leafcutter_cluster_regtools.py -j juncfiles -m 50 -l 500000
        Rscript /home/hgen_share/utils/leafcutter/scripts/leafcutter_ds.R \\
                --num_threads 1 \\
                --exon_file exons.txt.gz \\
                --min_samples_per_intron 3 \\
                leafcutter_perind_numers.counts.gz \\
                groupfiles
        """, language="bash")
        st.markdown("Pack up results for visualization")
        st.code("""
        Rscript /home/hgen_share/utils/leafcutter/leafviz/prepare_results.R \\
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
        st.markdown("#### Alternatively, we can run rMATS")
        
        st.markdown("Set up directory and prepare files for input into [rMATS](https://github.com/Xinglab/rmats-turbo/blob/v4.1.2/README.md)")
        st.code("""
        export PATH="/home/hgen_share/Anaconda/bin:$PATH"
        export PATH="/home/hgen_share/Anaconda/envs/rmats/bin:$PATH"
        mkdir output_rmats
        mkdir tmp_rmats
        echo "/home/hgen_share/lec7/ct1.bam,/home/hgen_share/lec7/ct2.bam,/home/hgen_share/lec7/ct3.bam" > control_list.txt
        echo "/home/hgen_share/lec7/kd1.bam,/home/hgen_share/lec7/kd2.bam,/home/hgen_share/lec7/kd3.bam" > kd_list.txt
        """, language="bash")
        st.markdown("Run rMATS")
        st.code("""
        python /home/hgen_share/Anaconda/envs/rmats/rMATS/rmats.py \\
                --b1 control_list.txt --b2 kd_list.txt \\
                --gtf /home/hgen_share/ref/hg38/ref-transcripts.gtf \\
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
        st.markdown("Import the rest of the results from rMATS (ie. those ending with ‘.JC.’, which are junction counts).")
        st.code("""
        se <- read.table("output_rmats/SE.MATS.JC.txt",sep = '\t', 
                        header = TRUE, 
                        stringsAsFactors = FALSE)
        """, language="r")
        st.markdown("#### Differential Splicing Events")
        st.markdown("Which type are the most abundant differentially spliced events (DSEs)? Is there a tendancy for the knockdown samples towards one of the DSEs? (for example exon skipping?)")
        st.code("""
        # Hint: for the differentially spliced events, if they increase in the
        # control samples, then they must decrease in the knockdown samples
        se_dse <- se %>%
        filter(FDR < 0.05)
        """, language="r")
        st.markdown("#### Pathway analysis of genes with DSEs")
        st.markdown("Do pathway analyses for genes for each of the different DSE. Are any pathways over-represented?")
        st.code("""
        # geneIDs and geneSymbols are included in the rMATS output. Remember to remove duplicated genes!
        #
        g <- gost(query = ,
                organism = 'hsapiens',
                custom_bg = ,
                user_threshold = 0.05,
                significant = FALSE,
                correction_method = "fdr")
        gostplot(g)
        """, language="r")
                
    with st.container(border=True):
        st.markdown("#### Extra")
        
        st.markdown("While the dataset above is a toy example, you can take a look at the rMATS output of a real dataset: `/home/hgen_share/lec7_extra/rmats_output_SNRPB`. If you have time, do the same type of analyses as before with this dataset of a SNRPB heterozygous mutant. As previously, sample 1 in the rMATS output are the control samples in triplicates whereas sample 2 are the heterozygous SNRPB mutants, in triplicates as well.")

