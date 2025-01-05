import streamlit as st


col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("Introduction to Genomics", divider=True)

    #############################################################################################

    with st.container(border=True):
        st.markdown("Log in using `ssh`")

        st.code("ssh <username>@hgen633.calculquebec.cloud", language="bash")
    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Set up a work directory at /project/def-sponsor00/$USER")

        st.markdown("Navigate to /project/def-sponsor00/$USER")
        st.code("cd /project/def-sponsor00/$USER", language="bash")

        st.markdown("Create a directory with `mkdir`")
        st.code("mkdir -p lec1", language="bash")
        st.code("cd lec1", language="bash")

    st.divider()
    #############################################################################################

    # with st.container(border=True):
    #     st.markdown("#### Look at the files for today's class from /home/hgen_share/lec1")
    #     st.markdown("Go to the created directory with `cd`")

    #     st.code("ll ~/projects/lec1", language="bash")

    #     # st.markdown("Copy over files with `cp`")
    #     # st.code("cp -R /project/def-sponsor00/hgen_share/lec1/* ./", language="bash")
    # st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Take a look at the files")
        
        st.markdown("Show the list of files with `ls`")
        st.code("ls /project/def-sponsor00/hgen_share/lec1/", language="bash")

        st.markdown("Assign the folder to a variable")
        st.code("data=/project/def-sponsor00/hgen_share/lec1", language="bash")

        st.markdown("Inspect the first few lines of ex* files with `head`")
        st.code("""
        head ${data}/ex1.fa
        head ${data}/ex2.fq
        """, language="bash")
    st.divider()

    #############################################################################################
    with st.container(border=True):
        st.markdown("#### Convert `ex2.fq` to the FASTA format")
        
        st.markdown("Inspect the seqtk module with `module spider`")
        st.code("module spider seqtk/1.4", language="bash")

        st.markdown("Load the seqtk module with `module load`")
        st.code("module load seqtk/1.4", language="bash")

        st.markdown("Inspect loaded modules with `module list`")
        st.code("module list", language="bash")

        st.markdown("Convert FASTQ to FASTA with `seqtk seq`")
        st.code("seqtk seq -a ${data}/ex2.fq > ex2.fa", language="bash")

        st.markdown("Inspect the first few lines of ex2.fa with `head`")
        st.code("head ex2.fa", language="bash")
    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Run FastQC on `ex2.fq` and compare with [another example](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html)")
        st.code("module load fastqc/0.12.1")

        st.markdown("Perform quality control analysis on ex2.fq with fastqc")
        st.code("cp ${data}/ex2.fq .", language="bash")
        st.code("fastqc ex2.fq", language="bash")
        
        st.markdown("Download results to your local computer using `scp`")
        st.code("scp <username>@hgen633.calculquebec.cloud:/directory/fastqc_results.html local_directory", language="bash")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Align the first 10 and then 1000 sequences of `ex2.fa` using BLAT")

        st.markdown("Load the seqkit module")
        st.code("module load seqkit/2.5.1")

        st.markdown("Subset the first 10 and then 1000 sequences of `ex2.fa` using `seqkit`")
        st.code("seqkit head -n 10 ex2.fa > sub_ex2.fa", language="bash")

        st.markdown("Time the alignment against hg19 with blat using a command of the form")
        st.code("""
        module load blat/3.7
        time blat ${data}/hg19_ref/hg19.fa sub_ex2.fa OUTPUT.psl
        """, language="bash")

        st.markdown("More information on the .psl format found [here](https://genomebrowser.wustl.edu/goldenPath/help/blatSpec.html)")
    st.divider()

    #############################################################################################
    
    st.markdown("#####################################################################################################")

    with st.container(border=True):
        st.markdown("#### Align `ex2.fq` using BWA")
        st.markdown("Align against the same reference with `bwa mem`")
        st.code("""
        module load bwa/0.7.17
        bwa mem ${data}/hg19_ref/hg19.fa ${data}/ex2.fq > ex2.sam
        """, language="bash")

        st.markdown("How many reads are there in `ex2.fq`? How much faster is BWA MEM compared to BLAT?")
        st.markdown("Take a look at the `*.sam` file with `less`")
        st.code("less ex2.sam", language="bash")

        st.markdown("Calculate summary statistics with `samtools flagstat`")
        st.code("""
        module load samtools/1.16.1
        samtools flagstat ex2.sam
        """, language="bash")

        st.markdown("More information on interpreting flagstat statistics found [here](https://www.htslib.org/doc/samtools-flagstat.html)")
    st.divider()
    #############################################################################################
    
    with st.container(border=True):
        st.markdown("#### Convert `ex2.sam` to the BAM format and index it")
        st.markdown("Sort and compress in one go with `samtools sort`")
        st.code("samtools sort ex2.sam > ex2.sorted.bam ", language="bash")

        st.markdown("Take a look at the `*.bam` file with `head`")
        st.code("head ex2.sorted.bam", language="bash")

        st.markdown("Compare the size difference after compression with `ls`")
        st.code("ls -lh ex2.*am", language="bash")

        st.markdown("Generate an index file with samtools index")
        st.code("samtools index ex2.sorted.bam", language="bash")
    st.divider()