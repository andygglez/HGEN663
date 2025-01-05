import streamlit as st

col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("Whole Genome Sequencing and analysis", divider=True)



#############################################################################################

    with st.container(border=True):
            st.markdown("#### Take a look at `HG002.g.vcf` from last week")
            
            st.markdown("View the file with `less`")
            st.markdown("Create a new directory and move the `vcf` file and its index to the new directory using `mv`")

    st.divider()
    #############################################################################################
    with st.container(border=True):
            st.markdown("#### Filter variants using GATK recommanded cut-offs")
            
            st.markdown("Set up variables")
            st.code("""
            REF=/home/hgen_share/bwa
            VCF=/home/hgen_share/vcf/b37
            PICARD_JAR=/home/hgen_share/utils/picard.jar
            GATK_JAR=/home/hgen_share/utils/GenomeAnalysisTK-3.8.1.jar
            export PATH="/home/hgen_share/utils/gatk-4.1.9.0:$PATH"
            module load StdEnv/2020 samtools/1.16
            module load nixpkgs/16.09 java/1.8.0_192
            """, language="bash")

            st.markdown("Replace the hashtag (#) with the recommended cut-offs found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)")
            st.code("""
            gatk --java-options "-Xmx4G" \\
                VariantFiltration \\
                -V HG002.g.vcf \\
                -filter "QD < #" --filter-name "QD" \\
                -filter "FS > #" --filter-name "FS" \\
                -filter "SOR > #" --filter-name "SOR" \\
                -filter "MQ < #" --filter-name "MQ" \\
                -filter "MQRankSum < #" --filter-name "MQRankSum" \\
                -filter "ReadPosRankSum < #" --filter-name "ReadPosRankSum" \\
                -O HG002.filt.g.vcf 2> >(grep -v undefined)
            """, language="bash")

            st.markdown("How many variants appear ok?")
            st.code("""
            grep -v '^#' HG002.filt.g.vcf | cut -f7 | grep 'PASS' | wc -l
            """, language="bash")

            st.markdown("How many total variants are called? Hint: use `grep`")
            st.markdown("Set up variables for [DeepVariant](https://github.com/google/deepvariant)")
            st.code("""
            module load StdEnv/2020 apptainer/1.1.3
            mkdir -p "${PWD}/data"
            mkdir -p "${PWD}/output"
            INPUT_DIR="${PWD}/data"
            OUTPUT_DIR="${PWD}/output"
            TMPDIR=.
            cp /home/hgen_share/lec3/* ${INPUT_DIR}
            """, language="bash")

            st.markdown("Run [DeepVariant](https://github.com/google/deepvariant) on the same sample from last week. **Note: this may take a few minutes**")
            st.code("""
            singularity run -B $PWD,/usr/lib/locale/ \\
                /home/hgen_share/deepvariant/deepvariant_1.4.0.sif \\
                /opt/deepvariant/bin/run_deepvariant \\
                --model_type=WGS \\
                --ref="${INPUT_DIR}"/b37.chr20.fa \\
                --reads="${INPUT_DIR}"/HG002.sorted.bam \\
                --regions "20:43200000-43300000" \\
                --output_vcf="${OUTPUT_DIR}"/HG002.dv.vcf.gz \\
                --output_gvcf="${OUTPUT_DIR}"/HG002.dv.g.vcf.gz \\
                --num_shards=1
            """, language="bash")

            st.markdown("""Unzip the `g.vcf` found in the output directory using `gunzip` then
            filter it using similar parameters as above. 
            Compare the number of variants that pass. **Note: you may need to set-up your variables 
            for GATK again**""")

    st.divider()
    #############################################################################################
    with st.container(border=True):
            st.markdown("#### (Extra) Compare hard filtering with `VQSR`")
            
            st.markdown("Copy whole genome variant calls")
            st.code("""
            cp /home/hgen_share/lec3_part2 ./
            """, language="bash")

            st.markdown("Since SNPs and indels need to be trained separately, start with SNPs")
            st.code("""
            gatk --java-options "-Xmx4G" \\
                VariantRecalibrator \\
                -R ${REF}/b37.fa \\
                -V raw.vcf \\
                --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $VCF/b37/hapmap.vcf \\
                --resource:omni,known=false,training=true,truth=false,prior=12.0 $VCF/b37/omni.vcf \\
                --resource:1000G,known=false,training=true,truth=false,prior=10.0 $VCF/b37/1KG.snps.vcf \\
                --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $VCF/b37/dbSNP.vcf \\
                -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \\
                -mode SNP \\
                -O recal_SNP.recal \\
                --tranches-file recal_SNP.tranches \\
                --rscript-file recal_SNP.plots.R
            """, language="bash")



            st.markdown("Apply VQSR to SNPs")

            st.code("""
            gatk --java-options "-Xmx4G" \\
                ApplyVQSR \\
                -R ${REF}/b37.fa \\
                -V raw.vcf \\
                -O raw.snp.vcf \\
                --truth-sensitivity-filter-level 99.5 \\
                --tranches-file recal_SNP.tranches \\
                --recal-file recal_SNP.recal \\
                -mode SNP
            """, language="bash")


            st.markdown("Now repeat for indels")
            st.code("""
            gatk --java-options "-Xmx4G" \\
                VariantRecalibrator \\
                -R ${REF}/b37.fa \\
                -V raw.snp.vcf \\
                --resource:mills,known=false,training=true,truth=true,prior=12.0 $VCF/b37/Mills.indels.vcf \\
                --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $VCF/b37/dbSNP.vcf \\
                --max-gaussians 4 \\
                -an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum \\
                -mode INDEL \\
                -O recal_INDEL.recal \\
                --tranches-file recal_INDEL.tranches \\
                --rscript-file recal_INDEL.plots.R
    
            gatk --java-options "-Xmx4G" \\
                ApplyVQSR \\
                -R ${REF}/b37.fa \\
                -V raw.snp.vcf \\
                -O recal.vcf \\
                --truth-sensitivity-filter-level 99.0 \\
                --tranches-file recal_INDEL.tranches \\
                --recal-file recal_INDEL.recal \\
                -mode INDEL
            """, language="bash")

            st.markdown("How many variants in the same region passed?")
            st.code("""
            gatk SelectVariants \\
                -R ${REF}/b37.fa \\
                -V recal.vcf \\
                -O recal.sub.vcf \\
                -L 20:43200000-43300000
                
            grep -v '^#' recal.sub.vcf | cut -f7 | grep 'PASS' | wc -l
            """, language="bash")

    st.divider()
    #############################################################################################
    with st.container(border=True):
            st.markdown("#### (Extra) Generate a VCF file with only variants that passed")
            
            st.markdown("Set up `bcftools` to filter PASS variants")
            st.code("""
            module load StdEnv/2020 gcc/9.3.0 bcftools/1.16
            """, language="bash")

            st.markdown("Run `bcftools`. The input file can be the VCF file produced from either variant caller")
            st.code("""
            bcftools view -f 'PASS,.' -O v -o output_filename input_file 
            """, language="bash")

    st.divider()
    #############################################################################################
    with st.container(border=True):
            st.markdown("#### (Extra) Use `VariantsToTableto` get a tsv file of variants")
            
            st.markdown("""`VariantsToTable`, by default, only extracts variants that passed 
            filtering. Use the `--show-filtered` parameter to show all variants.""")
            st.code("""
            REF=/home/hgen_share/bwa
            gatk VariantsToTable \\
                -R $REF/b37.fa \\
                -V input_vcf \\
                -F CHROM -F POS -F FILTER -F TYPE -GF AD -GF DP \\
                --show-filtered \\
                -O output_filename.tsv
            """, language="bash")

            st.markdown("You may need to reload `modules` for GATK again")
    st.divider()
    #############################################################################################
    with st.container(border=True):
            st.markdown("#### Add functional annotations with SnpEff")
            
            st.markdown("Copy whole genome variant calls")
            st.code("""
            cp /home/hgen_share/lec3_part2/* ./
            """, language="bash")

            st.markdown("Use `SnpEff` with mostly default parameters")
            st.code("""
            UTILS=/home/partage/utils
            SNPEFF_JAR=/home/hgen_share/utils/snpEff/snpEff.jar

            java -Xmx4G -jar $SNPEFF_JAR \\
                -c /home/hgen_share/utils/snpEff/snpEff.config \\
                -i vcf -o vcf -v hg19 \\
                ex.vcf > ex.snpeff.vcf
            """, language="bash")

            st.markdown("Download .html results file to your local computer using scp")
    st.divider()
    #############################################################################################
    with st.container(border=True):
            st.markdown("#### Annotate with existing resources using GEMINI")
            
            st.markdown("Add the directory containing `gemini`")
            st.code("""
            export PATH=/home/hgen_share/utils/gemini/anaconda/bin:$PATH
            """, language="bash")

            st.markdown("Generate annotated database with gemini")
            st.code("""
            module load mugqic/gemini/0.20.1
            gemini load --cores 4 --tempdir $PWD -t snpEff -v ex.snpeff.vcf ex.snpeff.db
            """, language="bash")

            st.markdown("Dump database to plain-text csv using sqlite")
            st.code("""
            db=ex.snpeff.db
            t=($(sqlite3 $db ".tables"))
            for i in "${t[@]}"; do
                sqlite3 $db<<- EXIT_HERE
                .mode csv
                .headers on
                .output $i.csv
                SELECT * FROM $i;
                .exit
                EXIT_HERE
                echo "$i.csv generated"
            done
            """, language="bash")


            st.markdown("Download `csv` file")
            st.markdown("For more details regarding columns: [GEMINI database](https://gemini.readthedocs.io/en/latest/content/database_schema.html)")

    st.divider()
    #############################################################################################

    st.markdown("### Let's move to R!")

    with st.container(border=True):
            st.markdown("#### Set up")
            
            st.code("""
            using <- function(...) {
            libs <- unlist(list(...))
            need <- libs[!unlist(lapply(libs, require, character.only = T))]
            if(length(need) > 0){ 
                install.packages(need, repos = "https://cloud.r-project.org")
                need <- need[!unlist(lapply(need, require, character.only = T))]
                if (length(need) > 0) {
                if (!requireNamespace("BiocManager", quietly = T))
                    install.packages("BiocManager")
                BiocManager::install(need)
                }
                }
            }
            using("ggplot2", "dplyr", "ggseqlogo", "GenomicRanges",
            "BSgenome.Hsapiens.UCSC.hg19", "biomaRt")
            """, language="r")

            st.markdown("#### Gene annotations")
            st.markdown("""As an R implementation of the BioMart project, biomaRt 
            facilitates access to a variety of resources ranging from Ensembl Genomes
            to the COSMIC cancer database. The chunk below retrieves the Ensembl gene set,
            subsets autosomal genes with valid HGNC symbols, converts the chromosome names
            to UCSC scheme, and finally derives a GRanges object from it.""")
            st.code("""
            # the useMart function connects to a specified
            #  BioMart database and dataset within this database.
            mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
                            host = "grch37.ensembl.org",
                            path = "/biomart/martservice", 
                            dataset = "hsapiens_gene_ensembl")

            # seqinfo extracts the sequence information stored 
            # in a DNAStringSet object
            si <- Seqinfo(genome = "hg19")[paste0("chr", 1:22)]

            # Given a set of filters and corresponding values, 
            # getBM retrieves the user specified attributes from 
            # the BioMart database
            genes <- getBM(attributes = c("hgnc_symbol", "chromosome_name", 
                                          "start_position","end_position", 
                                          "strand", "ensembl_gene_id"),
                           filters = "chromosome_name", 
                           values = as.character(1:22), 
                           mart = mart) %>%
            arrange(chromosome_name, start_position) %>%
            filter(hgnc_symbol != "") %>% # remove any element without a gene symbol
            mutate(strand = ifelse(strand == 1, '+', '-'),
                   chromosome_name = paste0('chr', chromosome_name)) %>%
            makeGRangesFromDataFrame(keep.extra.columns = T,
                                    start.field = "start_position",
                                    end.field = "end_position",
                                    seqinfo = si) # make granges object
            # The GRanges class is a container for the genomic locations
            #  and their associated annotations
            """, language="r")
            st.image("images/lec3.granges.png")

            st.markdown("#### Sequence analysis")
            st.code("""
            set.seed(42)
            genes.sub <- genes[strand(genes) == '+'] %>%
            sample(1000) # sample 1000 genes from positive strand

            # get exons from bioMart database
            exons.sub <- getBM(attributes = c("ensembl_exon_id", 
                                              "exon_chrom_start", "exon_chrom_end",
                                              "chromosome_name", "strand"), 
                               filters = "ensembl_gene_id",
                               values = genes.sub$ensembl_gene_id, 
                               mart = mart) %>%
            arrange(chromosome_name, exon_chrom_start) %>%
            mutate(strand = ifelse(strand == 1, '+', '-'),
                   chromosome_name = paste0('chr', chromosome_name)) %>%
            makeGRangesFromDataFrame(keep.extra.columns = T,
                                    start.field = "exon_chrom_start",
                                    end.field = "exon_chrom_end",
                                    seqinfo = si)

            # resize exons
            resize(exons.sub, fix = 'start', width = 10) %>% # resize returns an object of the same type and length as x containing intervals that have been resized to width width based on the strand(x) values
            resize(fix = 'end', width = 20) %>%
            getSeq(BSgenome.Hsapiens.UCSC.hg19, .) %>%
            consensusMatrix() %>%
            .[c('A', 'C', 'G', 'T'),] %>%
            ggseqlogo() # ggseqlogo is a shortcut for generating sequence logos
            """, language="r")
            st.image("images/lec3.motif.png")

            st.markdown("#### Additional exercise")
            st.markdown("Plot the sequence of splice donor sites using the same **exons.sub** object")
    st.divider()