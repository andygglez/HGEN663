import streamlit as st

col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("Whole Genome Sequencing and analysis", divider=True)

#############################################################################################

    with st.container(border=True):
            st.markdown("#### Take a look at `HG002.g.vcf` from last week")
            
            st.markdown("We will convert it to vcf first")
            st.code("""
            data=/project/def-sponsor00/hgen_share/lec3

            mkdir -p /project/def-sponsor00/$USER/lec3
            cd /project/def-sponsor00/$USER/lec3
            cp ${data}/HG002.sorted.dup.recal.ba* .
            cp ${data}/HG002.g.vcf .
            """, language="bash")

            st.code("""
            module load StdEnv/2023 gatk/4.4.0.0

            gatk GenotypeGVCFs \\
                --create-output-variant-index \\
                -R ${data}/hg19/hg19.chr20.fa \\
                -V HG002.g.vcf \\
                -O HG002.HaploCaller.vcf

            rm HG002.g.vcf
            """, language="bash")

            st.markdown("Examine the vcf file with `less`")
            st.code("less HG002.HaploCaller.vcf", language="bash")

    st.divider()
    #############################################################################################
    with st.container(border=True):
            st.markdown("#### We will call variants in a different way: by using DeepVariant")
            

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

            st.markdown("Now download the vcf files and look at them in IGV, together with the bam file")
            st.markdown("What are your conclusions? Go to the region chr20:43,277,043-43,300,080")
    
    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Place the variants in a table")
            
            st.code("""
            gatk VariantsToTable \\
                -R $data/hg19/hg19.chr20.fa \\
                -V /project/def-sponsor00/$USER/lec3/HG002.dv.vcf.gz \\
                -F CHROM -F POS -F FILTER -F TYPE -GF AD -GF DP \\
                --show-filtered \\
                -O variants.table.tsv
            """, language="bash")

    st.divider()

    with st.container(border=True):
            st.markdown("#### Add functional annotations with SnpEff")
            
            st.markdown("Use `SnpEff` with mostly default parameters")
            st.code("""
            UTILS=/project/def-sponsor00/hgen_share/utils/
            SNPEFF_JAR=${UTILS}/snpEff.jar

            java -Xmx4G -jar $SNPEFF_JAR \\
                -c ${UTILS}/snpEff.config \\
                -i vcf -o vcf -v hg19 \\
                HG002.dv.vcf > HG002.dv.snpeff.vcf
            
            mkdir -p dv_snpeff
            mv snpEff* dv_snpeff
            """, language="bash")

            st.divider()

            st.markdown("Let's repeat the analysis with another vcf file")

            st.code("""
            java -Xmx4G -jar $SNPEFF_JAR \\
                -c ${UTILS}/snpEff.config \\
                -i vcf -o vcf -v hg19 \\
                ${data}/ex.vcf > ex.snpeff.vcf
            """, language="bash")

            st.markdown("Download .html results file to your local computer using scp")

    st.divider()
 