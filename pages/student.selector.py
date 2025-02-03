import streamlit as st
import random

# Configure page
# st.set_page_config(page_title="Random Student Selector", layout="centered")

# Custom CSS for styling
st.markdown("""
    <style>
        .fancy-text {
            font-size: 3em;
            font-weight: bold;
            color: #FF4B4B;
            text-align: center;
            margin: 20px;
            padding: 20px;
            background: linear-gradient(45deg, #FFF3B0, #FFD700);
            border-radius: 15px;
            box-shadow: 0 4px 8px rgba(0,0,0,0.2);
            animation: fadeIn 1s ease-in;
        }
        @keyframes fadeIn {
            from { opacity: 0; }
            to { opacity: 1; }
        }
    </style>
""", unsafe_allow_html=True)

# Student list handling
st.title("🎉 Random Student Selector")

# Input for student names
students_input = st.text_area(
    "Enter student names (comma-separated)",
    "Alexander, Fiona, Yanchen, Cedric, Lang, Miranda, Seol, Zhan",
    height=150
)

# Process student names
student_list = [name.strip() for name in students_input.split(',') if name.strip()]

# Button to select random student
if st.button("🎯 Pick Random Student!", use_container_width=True):
    if len(student_list) > 0:
        selected_student = random.choice(student_list)
        st.markdown(
            f'<div class="fancy-text">✨ {selected_student} ✨</div>', 
            unsafe_allow_html=True
        )
    else:
        st.warning("Please enter some student names first!")