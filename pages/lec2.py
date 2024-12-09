import streamlit as st

col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("Human Genetic Variation", divider=True)

    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Copy the files for today's class from `/home/hgen_share/lec2`")
            st.markdown("Set up directory and copy over files")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Compare `HG002.sorted.bam` with what you generated last week")
            st.markdown("**Remember to load modules**")
            st.markdown("Index the `bam` file")
            st.markdown("Run `samtools flagstat` on both files")
            st.markdown("What are some of the similarities and differences?")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Visualize `HG002.sorted.bam` in IGV")
            st.markdown("Download the `bam` file to your local computer as previously shown using `scp`")
            st.markdown("Check what assembly it was aligned to")

            st.code("samtools view -H HG002.sorted.bam | tail -n 3", language="bash")

            st.markdown("Open it in IGV with the appropriate reference around 20:43200000-43300000")
    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Re-align `HG002.sorted.bam` around susceptible regions with GATK")
            st.markdown("Assign paths to variables for simplicity")
            st.code("""
            REF=/home/hgen_share/bwa
            VCF=/home/hgen_share/vcf/b37
            PICARD_JAR=/home/hgen_share/utils/picard.jar
            GATK_JAR=/home/hgen_share/utils/GenomeAnalysisTK-3.8.1.jar
            export PATH="/home/hgen_share/utils/gatk-4.1.9.0:$PATH"
            module load StdEnv/2020 samtools/1.16
            module load nixpkgs/16.09 java/1.8.0_192
            """, language="bash")

            st.markdown("Find regions with indels and dense in SNPs with `RealignerTargetCreator`")
            st.code("""
            java -Xmx4G -jar $GATK_JAR \\
                -T RealignerTargetCreator \\
                -R $REF/b37.fa \\
                -L 20:43200000-43300000 \\
                --known $VCF/dbSNP.vcf \\
                --known $VCF/Mills.indels.vcf \\
                -o realign.intervals \\
                -I HG002.sorted.bam
            """, language="bash")

            st.markdown("Re-align around the identified regions with `IndelRealigner`")
            st.code("""
            java -Xmx4G -jar $GATK_JAR \\
                -T IndelRealigner \\
                -R $REF/b37.fa \\
                -known $VCF/dbSNP.vcf \\
                -known $VCF/Mills.indels.vcf \\
                -targetIntervals realign.intervals \\
                -o HG002.realigned.sorted.bam \\
                -I HG002.sorted.bam
            """, language="bash")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Visualize re-aligned reads in IGV")
            
            st.markdown("""Retrieve re-aligned reads based with "OC" tag ("original CIGAR", implying that a new CIGAR string was generated)""")
            st.code("samtools view HG002.realigned.sorted.bam | grep 'OC' | cut -f1 | sort > oc.txt", language="bash")

            st.markdown("Sort reads by name using `picard SortSam` to allow efficient searching")
            st.code("""
            java -Xmx4G -jar $PICARD_JAR SortSam \\
                -I HG002.realigned.sorted.bam \\
                -O HG002.realigned.qs.bam \\
                -SO queryname
            
            java -Xmx2g -jar $PICARD_JAR SortSam \\
                -I HG002.sorted.bam \\
                -O HG002.qs.bam \\
                -SO queryname
            """, language="bash")

            st.markdown("Subset reads with name matching those found in `oc.txt` with `picard FilterSamReads`")
            st.code("""
            java -Xmx4G -jar $PICARD_JAR FilterSamReads \\
                -I HG002.realigned.qs.bam \\
                -O HG002.new.bam \\
                --FILTER includeReadList \\
                -RLF oc.txt \\
                -SO coordinate \\
                -CREATE_INDEX true
            
            java -Xmx2g -jar $PICARD_JAR FilterSamReads \\
                -I HG002.qs.bam \\
                -O HG002.old.bam \\
                --FILTER includeReadList \\
                -RLF oc.txt \\
                -SO coordinate \\
                -CREATE_INDEX true
            """, language="bash")


            st.markdown("Download the new and old `bam` files. What else must you download with your `bam` file for visualization?")
            st.markdown("Go to 20:43,257,030-43,257,326 to see an obvious example.")
    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Mark PCR duplicates with Picard")
            
            st.markdown("Mark PCR duplicate reads with `MarkDuplicates` ")
            st.code("""
            java -Xmx4G -jar ${PICARD_JAR} MarkDuplicates \\
                -REMOVE_DUPLICATES false \\
                -CREATE_INDEX true \\
                -I HG002.realigned.sorted.bam \\
                -O HG002.sorted.dup.bam \\
                -M HG002.sorted.dup.metrics
            """, language="bash")

            st.markdown("Check the results")
            st.code("""
            grep 'METRICS CLASS' -A 2 HG002.sorted.dup.metrics | column -t | less -S
            """, language="bash")

    st.divider()
    #############################################################################################

    with st.container(border=True):
            st.markdown("#### Apply base quality score recalibration with GATK")
            
            st.markdown("Establish a set of covariates for correction based on context and known SNPs with `BaseRecalibrator`")
            st.code("""
            gatk --java-options "-Xmx4G" \
                BaseRecalibrator \\
                -R $REF/b37.fa \\
                -known-sites $VCF/dbSNP.vcf \\
                -known-sites $VCF/Mills.indels.vcf \\
                -L 20:43200000-43300000 \\
                -O HG002.sorted.dup.recalibration_report.grp \\
                -I HG002.sorted.dup.bam
            """, language="bash")

            st.markdown("Recalibrate using the metrics calculated above with `PrintReads`")
            st.code("""
            gatk --java-options "-Xmx4G" \\
                ApplyBQSR \\
                -R $REF/b37.fa \\
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
                -R $REF/b37.fa \\
                -I HG002.sorted.dup.recal.bam \\
                -ERC GVCF \\
                -O HG002.g.vcf \\
                -L 20:43200000-43300000 \\
                --native-pair-hmm-threads 1
            """, language="bash")

            st.markdown("Download the `vcf` and recalibrated `bam` file to your local computer")
            st.markdown("Guide to viewing `vcf` files in IGV found [here](https://igv.org/doc/desktop/)")
            st.markdown("More information, including parameter details, for GATK HaplotypeCaller found [here](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller)")

    st.divider()