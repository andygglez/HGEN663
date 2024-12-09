import streamlit as st


col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("Introduction to Genomics", divider=True)

    #############################################################################################

    with st.container(border=True):
        st.markdown("Log in using `ssh`")

        st.code("ssh user@workshop2021a.vhost37.genap.ca", language="bash")
    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Set up a work directory at ~/projects/lec1")
        st.markdown("Create a directory with `mkdir`")

        st.code("mkdir -p ~/projects/lec1", language="bash")
    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Copy the files for today's class from /home/hgen_share/lec1")
        st.markdown("Go to the created directory with `cd`")

        st.code("cd ~/projects/lec1", language="bash")

        st.markdown("Copy over files with `cp`")
        st.code("cp -R /home/hgen_share/lec1/* ./", language="bash")
    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Take a look at the files you copied")
        st.markdown("Show the list of files with `ls`")

        st.code("ls -lh", language="bash")

        st.markdown("Inspect the first few lines of ex* files with `head`")
        st.code("""
        head ex1.fa
        head ex2.fq
        """, language="bash")
    st.divider()

    #############################################################################################
    with st.container(border=True):
        st.markdown("#### Convert `ex2.fq` to the FASTA format")
        
        st.markdown("Inspect the seqtk module with `module spider`")
        st.code("module spider seqtk/1.3", language="bash")

        st.markdown("Load the seqtk module with `module load`")
        st.code("module spider seqtk/1.3", language="bash")

        st.markdown("Inspect loaded modules with `module list`")
        st.code("module list", language="bash")

        st.markdown("Convert FASTQ to FASTA with `seqtk seq`")
        st.code("seqtk seq -a ex2.fq > ex2.fa", language="bash")

        st.markdown("Inspect the first few lines of ex2.fa with `head`")
        st.code("head ex2.fa", language="bash")
    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Set up a conda environment")
        st.markdown("Export conda environment")

        st.code("""
        export PATH="/home/hgen_share/Anaconda/bin:$PATH" 
        export PATH="/home/hgen_share/Anaconda/envs/lec1/bin:$PATH"
        """, language="bash")

        st.markdown("Check environment is in your path")
        st.code("echo $PATH", language="bash")

        st.markdown("More information on how to set up Anaconda environment [here](https://docs.conda.io/projects/conda/en/stable/index.html)")
    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Run FastQC on `ex2.fq` and compare with [another example](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html)")
        st.markdown("Export conda environment")

        st.markdown("Perform quality control analysis on ex2.fq with fastqc")
        st.code("fastqc ex2.fq", language="bash")
        
        st.markdown("Download results to your local computer using `scp`")
        st.code("scp user@workshop2021a.vhost37.genap.ca:/directory/fastqc_results.html local_directory", language="bash")

    st.divider()
    #############################################################################################

    with st.container(border=True):
        st.markdown("#### Align the first 10 and then 1000 sequences of `ex2.fa` using BLAT")
        st.markdown("Inspect the reference assembly")
        st.code("head /home/hgen_share/bwa/hg19.fa", language="bash")

        st.markdown("Assign the reference file to a variable for simplicity")
        st.code("ref=/home/hgen_share/bwa", language="bash")

        st.markdown("Subset the first 10 and then 1000 sequences of `ex2.fa` using `seqkit`")
        st.code("seqkit head -n 10 ex2.fa > sub_ex2.fa", language="bash")

        st.markdown("Time the alignment against hg19 with blat using a command of the form")
        st.code("""
        module load blat/3.5
        time blat $ref/hg19.fa sub_ex2.fa OUTPUT.psl
        """, language="bash")

        st.markdown("More information on the .psl format found [here](https://genomebrowser.wustl.edu/goldenPath/help/blatSpec.html)")
    st.divider()

    #############################################################################################
    
    with st.container(border=True):
        st.markdown("#### Align `ex2.fq` using BWA")
        st.markdown("Align against the same reference with `bwa mem`")
        st.code("""
        module load bwa/0.7.17
        bwa mem $ref/hg19.fa ex2.fq > ex2.sam
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