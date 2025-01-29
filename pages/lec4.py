import streamlit as st

col1, col2, col3 = st.columns([2,8,2])

with col2:
    st.header("Long Read and Single Molecule Sequencing", divider=True)

    #############################################################################################

    with st.container(border=True):
        st.markdown("#### In this first part of the class we give a general description of the tools that were used to call variants with the PacBio data")
        st.markdown("NOTE: :red[You won't be able to run this code!!!]")

        st.markdown("##### Alignment")
        st.code("""
        pbmm2 align -j 10 ${REF_FILE} \\
                ${BASEDIR}/RS_PacBio_WGS_202X/controls/fastq/R0040235.fastq R0040235.sorted.bam \\
                --sort --rg '@RG\tID:idR0040235\tSM:myR0040235'
        """, language="bash")
        
        st.markdown("##### Call variants")
        st.code("""
        singularity exec --bind /usr/lib/locale/:/usr/lib/locale/ \\
                --bind /lb/scratch/ebareke/tmp/:/lb/scratch/ebareke/tmp/ \\
                --bind ${BASEDIR}:/data/ \\
                --workdir /data/RS_PacBio_WGS_202X/tmp/R0040235 \\
                --writable-tmpfs \\
                ${BASEDIR}/RS_PacBio_WGS_202X/images/deepvariant_1.8.0.sif /opt/deepvariant/bin/run_deepvariant \\
                --model_type=PACBIO \\
                --ref /data/RS_PacBio_WGS_202X/ref/hg38.fasta \\
                --reads /data/RS_PacBio_WGS_202X/reads/ALIGNED/R0040235.sorted.bam \\
                --vcf_stats_report=true \\
                --output_vcf=/data/RS_PacBio_WGS_202X/final/hg38/SNV/${SAMPLE}.vcf.gz \\
                --output_gvcf=/data/RS_PacBio_WGS_202X/final/hg38/SNV/${SAMPLE}.g.vcf.gz \\
                --intermediate_results_dir=/data/RS_PacBio_WGS_202X/intermediate_results_dir/${SAMPLE} \\
                --num_shards=12
        """, language="bash")
           
        st.markdown("##### Call CNVs")
        st.code("""
        hificnv --bam ${INPUT_DIR}/${SAMPLE}.sorted.bam \\
                --maf ${MAF_FILE} \\
                --ref ${REF_FILE} \\
                --exclude ${EXCLUDE_REGION} \\
                --expected-cn ${EXPECTED_CN} \\
                --threads 4 \\
                --output-prefix ${SAMPLE}
        """, language="bash")
        
        st.markdown("##### Call SVs")
        st.code("""
        pbsv discover ${BASEDIR}/RS_PacBio_WGS_202X/reads/ALIGNED/R0040235.sorted.bam R0040235.svsig.gz
        tabix -c '#' -s 3 -b 4 -e 4 R0040235.svsig.gz
        pbsv call -j 10 ${REF_FILE} R0040235.svsig.gz R0040235.SV.vcf
        """, language="bash")

        st.markdown("##### Genotype repeats")
        st.code("""
        trgt genotype --genome ${REF_FILE} \\
                --repeats ${REPEAT_DEF} \\
                --reads ${INPUT_DIR}/${SAMPLE}.sorted.bam \\
                --output-prefix ${SAMPLE}
        """, language="bash")


    if "authenticated" not in st.session_state:
        st.session_state.authenticated = False
        
    if not st.session_state.authenticated:
        st.markdown("#### Password Protected Section")
        password = st.text_input("Enter the password: ", type="password")
        PASSWORD="bcftools"

        if st.button("Submit"):
            if password == PASSWORD:
                st.session_state.authenticated = True
                st.success("Access granted!")
            else:
                st.error("Invalid password. Please try again.")

    else:

        st.divider()
        with st.container(border=True):
            st.markdown("#### Compare Illumina vs PacBio SNVs: ")

            st.markdown("Load modules")
            st.code("module load bcftools/1.19")

            st.code("""
            ill_snv=/project/def-sponsor00/hgen_share/lec4/HG002.Illumina.dv.vcf.gz
            pacbio_snv=/project/def-sponsor00/hgen_share/lec4/HG002_PacBio_GRCh38.deepvariant.phased.chr22.vcf.gz

            # Only present in the first file
            user=/project/def-sponsor00/$USER/lec4

            bcftools isec -C ${ill_snv} ${pacbio_snv} -o ${user}/snv.illumina.uniq.txt
            bcftools isec -C ${pacbio_snv} ${ill_snv} -o ${user}/snv.pacbio.uniq.txt

            # Present in both
            bcftools isec -n=2 ${pacbio_snv} ${ill_snv} -o ${user}/snv.common.txt
            """, language="bash")

            st.markdown("#### Now compare Illumina vs PacBio SVs: ")

            st.code("""
            ill_sv=/project/def-sponsor00/hgen_share/lec4/HG002_Illumina_GRCh38.manta.chr22.vcf.gz
            pacbio_sv=/project/def-sponsor00/hgen_share/lec4/HG002_PacBio_GRCh38.pbsv.phased.chr22.vcf.gz

            # Only present in the first file
            bcftools isec -C ${ill_sv} ${pacbio_sv} -o ${user}/sv.illumina.uniq.txt
            bcftools isec -C ${pacbio_sv} ${ill_sv} -o ${user}/sv.pacbio.uniq.txt

            # Present in both
            bcftools isec -n=2 ${ill_sv} ${pacbio_sv} -o ${user}/sv.common.txt
            """, language="bash")

            st.markdown("#### Now inspect the files with unique and common sites:")

            st.code("less ${user}/snv.illumina.uniq.txt", language="bash")
            st.markdown("What can you see?")
            st.markdown("Look at the following regions and describe what you see:")

            st.markdown("""
            :blue[#1]- chr22:36,985,075-36,985,194 \\
            :blue[#2]- chr22:36,742,174-36,757,487 \\
            :blue[#3]- chr22:36,447,091-36,447,130 \\
            :blue[#4]- chr22:36,433,926-36,434,508 \\
            :blue[#5]- chr22:36,219,565-36,219,779 \\
            :blue[#6]- chr22:10,699,393-10,699,432
            """)