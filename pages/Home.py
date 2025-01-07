
import py3Dmol
from stmol import showmol


st.image("images/home-genomics.jpg", use_column_width=True)

st.markdown("# Welcome to HGEN663")

st.html("""<p class='justified-text'>
This website is designed to support your journey through HGEN663 by providing access to well-organized, practical coding resources. Whether you're a beginner or looking to sharpen your programming skills in biological data analysis, you'll find materials tailored to our class topics.
</p>""")

# st.markdown("##### ")
st.markdown("# Instructor")
st.markdown("##### ")

col1, col2, col3, col4, col5 = st.columns([1,5,1,5,1])

with col2:
    st.image("images/Jacek.jpg", caption="Professor Jacek Majewski, PhD")

with col4:
    st.html("""<p class='justified-text'>
    Professor Jacek Majewski received a PhD in Evolutionary Biology from 
    Wesleyan University (1999), a Master's in Electrical Engineering
    from Stanford University (1991) and a Bachelor in Physics, also from Stanford 
    (1990). He is Associate Professor, Department of Human Genetics, 
    Faculty of Medicine at McGill University and holds the Canada Research 
    Chair in Statistical Genetics. Prof. Majewski’s research is based on genomics 
    and bioinformatics analysis of high throughput data. The recent revolution 
    in massively parallel DNA sequencing has opened new venues 
    into numerous biomedical problems.</p>
    """)

st.markdown("# What do we do?")
st.markdown("##### ")
col1, col2, col3, col4, col5 = st.columns([1,5,1,5,1])
with col2:
    ######################
    cartoon_style = {"cartoon": {
                        "color": "spectrum",
                        "thickness": 0.4,
                        "spin_on": True
                    }}
    sphere_style = {'sphere':{}}

    if 'view' not in st.session_state:
        st.session_state.view = py3Dmol.view(query="pdb:3av2", width=480, height=480)
        st.session_state.view.setStyle(cartoon_style)
        st.session_state.view.spin(True)

    if 'current_style' not in st.session_state:
        st.session_state.current_style = 'cartoon'

    def change_mol_style():
        if st.session_state.current_style == 'cartoon':
            st.session_state.view.setStyle(sphere_style)
            st.session_state.current_style = 'sphere'
        else:
            st.session_state.view.setStyle(cartoon_style)
            st.session_state.current_style = 'cartoon'

    st.html("""<p class='justified-text'>
        At <a href="https://majewskilab.github.io/">Majewski's Lab</a> we aim to elucidate the role of
        lysine 36 methyltransferases in diseases and developmental syndromes. Using 
        Next-Generation Sequencing (NGS) technologies, we investigate how alterations 
        in these enzymes impact the transcriptome, epigenome, and chromatin conformation.</p>
    """)

    st.button(label="Change Style", on_click=change_mol_style)
    
with col4:
    showmol(st.session_state.view, height=450, width=500)

st.markdown("""
## What You’ll Find Here:
""")
st.markdown("Scripts and Code Examples:")
st.html("""<p class='justified-text'>
Ready-to-use scripts in class, written in Python, R, and bash, 
focused on analyzing genomic data, and exploring bioinformatics tools.</p>
""")
st.markdown("#### Hands-On Tutorials:")
st.html("""<p class='justified-text'>
Step-by-step guides that walk you through key computational techniques used in modern biology.
</p>""")

st.markdown("#### Project Resources: ")
st.html("""<p class='justified-text'>
Our website will provide you with toy datasets for class 
projects to help you hit the ground running.</p>""")

st.markdown("#### Reference Material:")
st.html("""<p class='justified-text'>
Links to relevant documentation and further reading to deepen
 your understanding of the methods and tools we use.</p>""")

st.markdown("#### Why Use This Site?")
st.html("""<p class='justified-text'>
Computational biology is a fast-growing field that bridges biology,
 computer science, and data analysis. By engaging with the resources here, 
 you’ll develop the coding skills needed to tackle real-world biological questions, 
 from differential gene expression analysis to evaluating changes in chromatin conformation.
Feel free to explore, experiment, and expand your coding expertise as you dive into the fascinating
 intersection of computation and biology.</p>""")
