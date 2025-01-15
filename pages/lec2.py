import streamlit as st

col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("Human Genetic Variation", divider=True)

    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Look at the files for today's class from `/project/def-sponsor00/hgen_share/lec2`")
            st.code("ls /project/def-sponsor00/hgen_share/lec2")

            st.markdown("Store the location in a variable")
            st.code("data=/project/def-sponsor00/hgen_share/lec2")

            st.markdown("Create a directory in your folder")
            st.code("""
            mkdir -p /project/def-sponsor00/$USER/lec2
            cd /project/def-sponsor00/$USER/lec2
            cp ${data}/HG002.sorted.bam .
            """)

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Visualize `HG002.sorted.bam` in IGV")
            st.markdown("Download the `bam` file to your local computer as previously shown using `scp`")
            st.markdown("Check what assembly it was aligned to")

            st.code("""
            module load StdEnv/2023
            module load samtools/1.20
            
            samtools view -H ${data}/HG002.sorted.bam | tail -n 4""", language="bash")

            st.markdown("Open it in IGV with the appropriate reference around chr20:43200000-43300000") #####################################3
    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Re-align `HG002.sorted.bam` around susceptible regions with GATK")
            st.markdown("Assign paths to variables for simplicity")
            st.code("""
            GATK_JAR=/project/def-sponsor00/hgen_share/utils/GenomeAnalysisTK-3.8.1.jar
            PICARD_JAR=/project/def-sponsor00/hgen_share/utils/picard.jar

            data=/project/def-sponsor00/hgen_share/lec2
            REF=$data/hg19/hg19.chr20.fa

            module load StdEnv/2020
            module load java/1.8.0_192
            """, language="bash")

            st.markdown("Find regions with indels and dense in SNPs with `RealignerTargetCreator`")
            st.code("""
            java -Xmx4G -jar $GATK_JAR \\
                    -T RealignerTargetCreator \\
                    -R $REF \\
                    -L chr20:43200000-43300000 \\
                    -o realign.intervals \\
                    -I $data/HG002.sorted.bam
            """, language="bash")

            st.markdown("Re-align around the identified regions with `IndelRealigner`")
            st.code("""
            java -Xmx4G -jar $GATK_JAR \\
                    -T IndelRealigner \\
                    -R $REF \\
                    -targetIntervals realign.intervals \\
                    -o HG002.realigned.sorted.bam \\
                    -I $data/HG002.sorted.bam
            """, language="bash")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Visualize re-aligned reads in IGV")
            
            st.markdown("""Retrieve re-aligned reads based with "OC" tag ("original CIGAR", implying that a new CIGAR string was generated)""")
            st.code("""
            module load StdEnv/2023
            module load samtools/1.20

            samtools view HG002.realigned.sorted.bam | grep 'OC' | cut -f1 | sort > oc.txt
            """, language="bash")

            st.markdown("Sort reads by name using `picard SortSam` to allow efficient searching and subset reads with name matching those found in `oc.txt` with `picard FilterSamReads`")

            st.code("""
            module load StdEnv/2020
            module load java/1.8.0_192
            """)

            st.code("""
            java -Xmx4G -jar $PICARD_JAR SortSam \\
                    I=HG002.realigned.sorted.bam \\
                    O=HG002.realigned.qs.bam \\
                    SORT_ORDER=queryname

            java -Xmx4G -jar $PICARD_JAR FilterSamReads \\
                    I=HG002.realigned.qs.bam \\
                    O=HG002.new.bam \\
                    FILTER=includeReadList \\
                    RLF=oc.txt \\
                    SORT_ORDER=coordinate \\
                    CREATE_INDEX=true

            java -Xmx4G -jar $PICARD_JAR SortSam \\
                    I=$data/HG002.sorted.bam \\
                    O=HG002.qs.bam \\
                    SORT_ORDER=queryname

            java -Xmx4G -jar $PICARD_JAR FilterSamReads \\
                    I=HG002.qs.bam \\
                    O=HG002.old.bam \\
                    FILTER=includeReadList \\
                    RLF=oc.txt \\
                    SORT_ORDER=coordinate \\
                    CREATE_INDEX=true
            """, language="bash")

            st.markdown("Download the new and old `bam` files. What else must you download with your `bam` file for visualization?")
            st.markdown("Go to chr20:43,257,014-43,257,372 to see an obvious example.")
    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Mark PCR duplicates with Picard")
            
            st.markdown("Mark PCR duplicate reads with `MarkDuplicates` ")
            st.code("""
            java -Xmx4G -jar ${PICARD_JAR} MarkDuplicates \\
                    REMOVE_DUPLICATES=false \\
                    CREATE_INDEX=true \\
                    I=HG002.realigned.sorted.bam \\
                    O=HG002.sorted.dup.bam \\
                    METRICS_FILE=HG002.sorted.dup.metrics
            """, language="bash")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Apply base quality score recalibration with GATK")
            
            st.markdown("Establish a set of covariates for correction based on context and known SNPs with `BaseRecalibrator`")

            st.markdown("Use GATK4 in this case")

            st.code("""
            module load StdEnv/2023
            module load gatk/4.4.0.0
            """)

            st.code("""
            gatk --java-options "-Xmx4G" BaseRecalibrator \\
                -R $REF \\
                -known-sites $data/dbSNP.vcf \\
                -known-sites $data/Mills.indels.vcf \\
                -L chr20:43200000-43300000 \\
                -O HG002.sorted.dup.recalibration_report.grp \\
                -I HG002.sorted.dup.bam
            """, language="bash")

            st.markdown("Recalibrate using the metrics calculated above with `ApplyBQSR`")
            st.code("""
            gatk --java-options "-Xmx4G" ApplyBQSR \\
                -R $REF \\
                --bqsr-recal-file HG002.sorted.dup.recalibration_report.grp \\
                -O HG002.sorted.dup.recal.bam \\
                -I HG002.sorted.dup.bam
            """, language="bash")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Call variants with GATK and visualize results in IGV")
            st.markdown("Perform variant calling with `HaplotypeCaller`")

            st.code("""
            gatk --java-options "-Xmx4G" \\
                HaplotypeCaller \\
                -R $REF \\
                -I HG002.sorted.dup.recal.bam \\
                -ERC GVCF \\
                -O HG002.g.vcf \\
                -L chr20:43200000-43300000 \\
                --native-pair-hmm-threads 1
            """, language="bash")

            st.markdown("Download the `vcf` and recalibrated `bam` file to your local computer")
            st.markdown("Guide to viewing `vcf` files in IGV found [here](https://igv.org/doc/desktop/)")
            st.markdown("More information, including parameter details, for GATK HaplotypeCaller found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)")

    st.divider()