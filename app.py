import streamlit as st
from streamlit_option_menu import option_menu
# import py3Dmol
# from stmol import showmol

st.set_page_config(
    page_title="HGEN663",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
    <style>
        .css-79elbk{
            display: none;
        } /**/
        .st-emotion-cache-79elbk{
            display: none;
        }
        .e1fqkh3o4{
            margin-top: 30px;
        }

        .st-emotion-cache-bjn8wh.eczjsme17{
        display: none;        
        }
        
        .justified-text {
        text-align: justify;
        min-width: 450px;
        font-size: 20px; /* Adjust size as needed */
        line-height: 1.6; /* Line height for better readability */
        }
        p {
        text-wrap: wrap;
        font-size: 20px;
        }
        .hr-header {
        padding-top: 0 px;
        }
        .st-emotion-cache-14noe3u {
        font-size: 22px;
        }
        pre code {
        font-size: 18px;
        }
        .st-emotion-cache-1fmfajh{
        border: 1px solid #1c7fd6;
        }
        code {
        color: #1c7fd6;
        }
    </style>
    """, unsafe_allow_html=True)

lectures = ["Lecture "+str(i) for i in range(1,13)]
options = ["HGEN663"]
options.extend(lectures)
options.append("---")
options.append("External Resources")
options.append("Student Selector")
options.append("---")

icons = ['house']
icons.extend([str(i)+"-square" for i in range(1, len(lectures)+1)])

with st.sidebar:
    page_selected = option_menu(None, 
        options=options, icons=['house',"","","","",
                                        "","","",""
                                        "","","","",
                                        "","","box-arrow-up-right", "bullseye"],default_index=0)

if page_selected == "HGEN663":
    exec(open("pages/Home.py").read())

if page_selected == "Lecture 1":
    exec(open('pages/lec1.py').read())

if page_selected == "Lecture 2":
    exec(open('pages/lec2.py').read())

if page_selected == "Lecture 3":
    exec(open('pages/lec3.py').read())

if page_selected == "Lecture 4":
    exec(open('pages/lec4.py').read())

if page_selected == "Lecture 5":
    exec(open('pages/lec5.py').read())

if page_selected == "Lecture 6":
    exec(open('pages/lec6.py').read())

if page_selected == "Lecture 7":
    exec(open('pages/lec7.py').read())

if page_selected == "Lecture 8":
    exec(open('pages/lec8.py').read())

if page_selected == "Lecture 9":
    exec(open('pages/lec9.py').read())

if page_selected == "Lecture 10":
    exec(open('pages/lec10.py').read())

if page_selected == "Lecture 11":
    exec(open('pages/lec11.py').read())

if page_selected == "Lecture 12":
    exec(open('pages/lec12.py').read())

if page_selected == "External Resources":
    exec(open('pages/external.resources.py').read())

if page_selected == "Student Selector":
    exec(open('pages/student.selector.py').read())