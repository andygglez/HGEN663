import streamlit as st

col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("Whole Genome Sequencing and analysis", divider=True)

#############################################################################################

    with st.container(border=True):
            st.markdown("#### Take a look at `HG002.g.vcf` from last week")
            
            st.markdown("View the file with `less`")
            st.markdown("Create a new directory and move the `vcf` file and its index to the new directory using `mv`")
            st.code("""
            mkdir -p /project/def-sponsor00/$USER/lec3
            cd /project/def-sponsor00/$USER/lec3
            data=/project/def-sponsor00/hgen_share/lec3
            cp ${data}/HG002.g.vcf .
            """)

    st.divider()
    #############################################################################################
    with st.container(border=True):
            st.markdown("#### Filter variants using GATK recommanded cut-offs")
            
            st.markdown("Set up variables")
            st.code("""
            REF=/project/def-sponsor00/hgen_share/lec3/hg19
            dbSNP=/project/def-sponsor00/hgen_share/lec3/dbSNP.vcf
            Mills=/project/def-sponsor00/hgen_share/lec3/Mills.indels.vcf
            VCF=${data}/HG002.g.vcf
            
            GATK_JAR=/project/def-sponsor00/hgen_share/utils/GenomeAnalysisTK-3.8.1.jar
            PICARD_JAR=/project/def-sponsor00/hgen_share/utils/picard.jar

            module load StdEnv/2023 gatk/4.4.0.0

            """, language="bash")

            st.markdown("The recommended cut-offs are found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)")

            st.code("""
            gatk --java-options "-Xmx4G" \\
                VariantFiltration \\
                -V ${VCF} \\
                -filter "QD < 2.0" --filter-name "QD" \\
                -filter "FS > 60.0" --filter-name "FS" \\
                -filter "SOR > 3.0" --filter-name "SOR" \\
                -filter "MQ < 40.0" --filter-name "MQ" \\
                -filter "MQRankSum < -12.5" --filter-name "MQRankSum" \\
                -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum" \\
                -O HG002.filt.g.vcf &> log.txt
            """, language="bash")

            st.markdown("How many variants appear ok?")
            st.code("""
            grep -v '^#' HG002.filt.g.vcf | cut -f7 | grep 'PASS' | wc -l
            """, language="bash")

            st.markdown("Run [DeepVariant](https://github.com/google/deepvariant) on the same sample from last week. **Note: this may take a few minutes**")
            st.code("""
            module load apptainer/1.3.5

            singularity run -B $PWD,/usr/lib/locale/ \\
                        /project/def-sponsor00/hgen_share/lec3/deepvariant_1.4.0.sif \\
                        /opt/deepvariant/bin/run_deepvariant \\
                        --model_type=WGS \\
                        --ref=$data/hg19/hg19.chr20.fa \\
                        --reads=${data}/HG002.sorted.bam \\
                        --regions "chr20:43200000-43300000" \\
                        --output_vcf=/project/def-sponsor00/$USER/lec3/HG002.dv.vcf.gz \\
                        --output_gvcf=/project/def-sponsor00/$USER/lec3/HG002.dv.g.vcf.gz \\
                        --num_shards=1
            """, language="bash")

            st.markdown("""Unzip the `g.vcf` found in the output directory using `gunzip` then
            filter it using similar parameters as above.
            Compare the number of variants that pass. **Note: you may need to set-up your variables 
            for GATK again**""")

            st.code("gunzip HG002.dv.g.vcf.gz", language="bash")

            st.code("""
            gatk --java-options "-Xmx4G" \\
                VariantFiltration \\
                -V HG002.dv.g.vcf \\
                -filter "QD < 2.0" --filter-name "QD" \\
                -filter "FS > 60.0" --filter-name "FS" \\
                -filter "SOR > 3.0" --filter-name "SOR" \\
                -filter "MQ < 40.0" --filter-name "MQ" \\
                -filter "MQRankSum < -12.5" --filter-name "MQRankSum" \\
                -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum" \\
                -O HG002.dv.filtered.g.vcf &> log.txt

            grep -v '^#' HG002.dv.filtered.g.vcf | cut -f7 | grep 'PASS' | wc -l
            """, language="bash")

    st.divider()
    #############################################################################################
    with st.container(border=True):
            st.markdown("#### Compare hard filtering with `VQSR`")
            
            st.markdown("Copy whole genome variant calls")
            st.code("""
            cp ${data}/raw.vcf .
            """, language="bash")

            st.markdown("Since SNPs and indels need to be trained separately, start with SNPs")
            st.code("""
            b37=/project/def-sponsor00/hgen_share/lec3/b37

            dbSNP=${b37}/dbSNP.vcf
            hapmap=${b37}/hapmap.vcf
            omni=${b37}/omni.vcf
            Mills=${b37}/Mills.indels.vcf
            Snps=${b37}/1KG.snps.vcf
            
            ### Note: this might take some time to run
            gatk --java-options "-Xmx4G" \\
                VariantRecalibrator \\
                -R ${b37}/b37.fa \\
                -V ${data}/raw.vcf \\
                --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \\
                --resource:omni,known=false,training=true,truth=false,prior=12.0 ${omni} \\
                --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${Snps} \\
                --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbSNP} \\
                -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum \\
                -mode SNP \\
                -O recal_SNP.recal \\
                --tranches-file ${data}/recal_SNP.tranches
            """, language="bash")

            st.markdown("Apply VQSR to SNPs")

            st.code("""
            gatk --java-options "-Xmx4G" \\
                ApplyVQSR \\
                -R ${b37}/b37.fa \\
                -V raw.vcf \\
                -O raw.snp.vcf \\
                --truth-sensitivity-filter-level 99.5 \\
                --tranches-file  ${data}/recal_SNP.tranches \\
                --recal-file recal_SNP.recal \\
                -mode SNP
            """, language="bash")

            st.markdown("Now repeat for indels")
            st.code("""
            ### Note: this might take some time to run
            gatk --java-options "-Xmx4G" \\
                VariantRecalibrator \\
                -R ${b37}/b37.fa \\
                -V raw.snp.vcf \\
                --resource:mills,known=false,training=true,truth=true,prior=12.0 ${Mills} \\
                --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbSNP} \\
                --max-gaussians 4 \\
                -an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum \\
                -mode INDEL \\
                -O recal_INDEL.recal \\
                --tranches-file ${data}/recal_INDEL.tranches
    
            gatk --java-options "-Xmx4G" \\
                ApplyVQSR \\
                -R ${b37}/b37.fa \\
                -V raw.snp.vcf \\
                -O recal.vcf \\
                --truth-sensitivity-filter-level 99.0 \\
                --tranches-file ${data}/recal_INDEL.tranches \\
                --recal-file recal_INDEL.recal \\
                -mode INDEL
            """, language="bash")

            st.markdown("How many variants in the same region passed?")
            st.code("""
            gatk SelectVariants \\
                -R ${b37}/b37.fa \\
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
            module load StdEnv/2023 gcc/12.3 bcftools/1.19
            """, language="bash")

            st.markdown("Run `bcftools`. The input file can be the VCF file produced from either variant caller")
            st.code("""
            bcftools view -f 'PASS,.' -O v -o recal.passed.vars.vcf recal.vcf 
            """, language="bash")

    st.divider()
    #############################################################################################
    with st.container(border=True):
            st.markdown("#### (Extra) Use `VariantsToTableto` get a tsv file of variants")
            
            st.markdown("""`VariantsToTable`, by default, only extracts variants that passed 
            filtering. Use the `--show-filtered` parameter to show all variants.""")
            st.code("""
            gatk VariantsToTable \\
                -R ${b37}/b37.fa \\
                -V recal.sub.vcf \\
                -F CHROM -F POS -F FILTER -F TYPE -GF AD -GF DP \\
                --show-filtered \\
                -O variants.table.tsv
            """, language="bash")

            st.markdown("You may need to reload `modules` for GATK again")
    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Add functional annotations with SnpEff")
            
            st.markdown("Use `SnpEff` with mostly default parameters")
            st.code("""
            UTILS=/project/def-sponsor00/hgen_share/utils/
            SNPEFF_JAR=${UTILS}/snpEff.jar

            java -Xmx4G -jar $SNPEFF_JAR \\
                -c ${UTILS}/snpEff.config \\
                -i vcf -o vcf -v hg19 \\
                ${data}/ex.vcf > ex.snpeff.vcf
            """, language="bash")

            st.markdown("Download .html results file to your local computer using scp")
    st.divider()
