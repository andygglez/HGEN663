import streamlit as st

col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("Long Read and Single Molecule Sequencing", divider=True)

    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Copy files for part 1 from `/home/hgen_share/lec4_pt1`")
            
            st.markdown("Set up directory and copy over files")
            st.code("""
            module load StdEnv/2020 gcc/9.3.0 openmpi/4.0.3 hdf5/1.10.6 minimap2/2.17 fastqc/0.11.9 samtools/1.11 mugqic/qualimap/2.2.1 mafft/7.471 emboss/6.6.0 fasttree/2.1.11 mugqic/pycoQC/2.5.2
            export PATH="/home/hgen_share/Anaconda/bin:$PATH"
            export PATH="/home/hgen_share/Anaconda/envs/longread_env/bin:$PATH"
            export PATH=/home/hgen_share/ont-guppy-cpu/bin:$PATH
            """, language="bash")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Take a look at ex.fast5")
            
            st.markdown("View the file with `h5ls`")
            st.code("""
            h5ls ex.fast5 | less
            """, language="bash")

            st.markdown("Show structure of one read using `h5ls`")
            st.code("""
            h5ls -r ex.fast5/read_00017b75-18fc-40ae-8a21-a5b3e49dc753
            """, language="bash")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Basecall using Guppy")
            
            st.markdown("This will take a while with just 1 CPU thread, and since the output is provided already, you can terminate it prematurely with `Ctrl-C`")
            st.code("""
            conf="dna_r9.4.1_450bps_hac"

            guppy_basecaller \\
                --compress_fastq \\
                -i ex.fast5 \\
                -s basecall/ \\
                --cpu_threads_per_caller 1 \\
                --num_callers 1 \\
                -c dna_r9.4.1_450bps_hac.cfg
            """, language="bash")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### (Extra) Assess quality with pycoQC & FastQC")
            
            st.markdown("Run `pycoQC` with the provided summary")
            st.code("""
            pycoQC -f sequencing_summary.txt -o pyco.html
            """, language="bash")

            st.markdown("Run `FastQC` with the provided sequences")
            st.code("""
            fastqc ex.fq
            """, language="bash")

            st.markdown("Download `pycoQC` html report")
            st.markdown("Compare with another report [here](https://a-slide.github.io/pycoQC/pycoQC/results/Guppy-2.3_basecall-1D_alignment-DNA.html)")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### (Extra) Rectify quality issues by trimming")
            
            st.markdown("Trim adapters with `Porechop`")
            st.code("""
            porechop -i ex.fq -t 1 -v 2 -o chopped.fq
            """, language="bash")

            st.markdown("Is it better now? Run fastQC on trimmed `fastq`")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### (Extra) Assess error profile after alignment")
            
            st.markdown("Align trimmed reads using `minimap2`")
            st.code("""
            minimap2 -t 1 -x map-ont -a wuhCor1.fa chopped.fq | samtools sort -o ex.bam
            samtools index ex.bam
            """, language="bash")

            st.markdown("Summarize alignment results with `Qualimap`")
            st.code("""
            qualimap bamqc -bam ex.bam --java-mem-size=1G -nw 5000 -nt 1 -c -outdir qualimap/
            """, language="bash")

            st.markdown("Download results to your local computer")
            st.markdown("Compare with the [report](https://rawgit.com/kokonech/kokonech.github.io/master/qualimap/HG00096.chr20_bamqc/qualimapReport.html) from a typical Illumina run")
    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Hybrid assembly with Unicycler")
            
            st.markdown("Set up variables and copy over files for part 2")
            st.code("""
            export PATH=/home/hgen_share/utils/SPAdes-3.15.5-Linux/bin:$PATH
            export PATH=/home/hgen_share/utils/ncbi-blast-2.13.0+/bin:$PATH
            export PATH=/home/hgen_share/utils/racon/build/bin:$PATH
            export PATH="/home/hgen_share/Anaconda/bin:$PATH"
            export PATH="/home/hgen_share/Anaconda/envs/longread_env/bin:$PATH"
            cp /home/hgen_share/lec4_pt2/* ./
            """, language="bash")

            st.markdown("""Run `Unicycler` using boths short- and long-reads. 
            For hybrid assemblies, the direct long-read bridging step of the
            pipeline can take a while to complete. Thus, the output is provided
            for you. You can terminate it prematurely with `Ctrl-C`""")
            st.code("""
            unicycler -1 Illumina_F.fq -2 Illumina_R.fq \\
                    -l ONT.fq -t 20 --min_fasta_length 10000 \\
                    --linear_seqs 1 -o unicycler
            """, language="bash")

            st.markdown("The output can be found in `/home/hgen_share/lec4_pt2/unicycler`")
    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Evaluate the assembled genome's quality")
            
            st.markdown("Use `QUAST` to create a summary report")
            st.code("""
            module load mugqic/Quast/5.0.2
            quast.py -t 1 -o quast -R wuhCor1.fa -g wuhCor1.gff3 unicycler/assembly.fasta
            """, language="bash")

            st.markdown("Compare with a [report](http://cab.cc.spbu.ru/quast/paper/h.sapiens_chr14/) on eukaryotic assemblies")

            st.markdown("[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi) the top 2 contigs (`contig1.fa` and `contig2.fa`) found in the `unicycler` output directory")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Place newly assembled viral genomes in a global context")
            
            st.markdown("Set up variables for nextstrain environment")
            st.code("""
            export PATH="/home/hgen_share/Anaconda/bin:$PATH"
            export PATH="/home/hgen_share/Anaconda/envs/nextstrain_V2/bin:$PATH"
            """, language="bash")

            st.markdown("Copy the tutorial dataset for SARS-Cov2 provided by nextstrain")
            st.code("""
            cp -r /home/hgen_share/ncov ./
            """, language="bash")

            st.markdown("""Take a look at the `/data` sub-directory. Inspect the data 
            files with `less` to get an impression of what data needs to be provided 
            for nextstrain. Unzip the fasta file using `gunzip`""")
            st.markdown("Go back to the main `/ncov` directory and run the basic nextstrain workflow")

            st.code("""
            snakemake --cores 14 --profile ./my_profiles/getting_started
            """, language="bash")

            st.markdown("""On your own computer, not the cluster, establish a tunnel 
            to server. Here you should use your own username.Also, replace the hashtag
            with any two numbers to land on a port that is not busy""")
            st.code("""
            ssh -L 40##:localhost:40## studXX@workshop2021a.vhost37.genap.ca
            """, language="bash")

            st.markdown("""Start the auspice server to view the results using your credentials.
            If the PORT is busy, change to another port by exporting the port number as indicated below""")
            st.code("""
            export PORT=40##
            auspice view --datasetDir /home/studXX/projects/lec4/ncov/auspice
            """, language="bash")

            st.markdown("You can visit the results page at [](http://localhost:40##) in your browser. If you're done inspecting the results, terminate the auspice server with `Ctrl+C`")
            st.markdown("""Include our own recently assembled data. Use `consensus_sequences.fasta`
            and `metadata.tsv` from `/home/hgen_share/lec4_pt2`. Append the fasta file to the example
            sequences (you can also inspect the fasta file with `less`""")

            st.code("""
            cat consensus_sequences.fasta >> data/example_sequences.fasta
            """, language="bash")

            st.markdown("Do the same for the metadata file. Append it to `example_metadata.tsv`")
            st.markdown("Run the same steps as above. Where does our data reside in the phylogeny?")