# External Modules
import pickle
from pathlib import Path

import streamlit as st
import pandas as pd
from pandasai import SmartDataframe
from pandasai.llm import OpenAI

import streamlit_authenticator as stauth

# Internal Modules
from genesys.visuals import render_protein_file
from genesys.ai import run_conversation
from genesys.DNAToolKit import sequence_type
from genesys.client import upload_content_to_s3, get_s3_url

# Standard Library
import os
from genesys.env import load_dotenv
import logging
from io import StringIO
from time import time
import random

# --- USER AUTHENTICATION

names = ["Charlie", "Eshcol", "Sandeep", "Raj", "Ashish"]
usernames = ["charlie", "eshcol", "sandeep", "raj", "ashish"]

# Load hashed passwords
file_path = Path("genesys/hashed_pw.pkl")
with file_path.open("rb") as file:
    hashed_passwords = pickle.load(file)

authenticator = stauth.Authenticate(
    names, usernames, hashed_passwords, "sales_dashboard", "abcdef", cookie_expiry_days=0)

header = """
<div class="fade-in" style='text-align: center;'>
    <h1>🧬 GeneSys AI 🧬</h1>
</div>
"""

subheader = """
<div class="fade-in" style='text-align: center;'>
    <p><em>Making it as easy as AUG</em></p>
</div>
"""

# Inject custom CSS and content
st.markdown(
    f"""
    <style>
        .fade-in {{
            opacity: 0;
            animation: fade-in 1s ease-in-out forwards;
        }}
        @keyframes fade-in {{
            0% {{
                opacity: 0;
            }}
            100% {{
                opacity: 1;
            }}
        }}
    </style>
    {header}
    {subheader}
    """,
    unsafe_allow_html=True,
)

name, authentication_status, username = authenticator.login("Login", "main")

if authentication_status == False:
    st.error("Username/password is incorrect")

if authentication_status == None:
    st.warning("Please enter your username and password")

if authentication_status:
    authenticator.logout("Logout", "sidebar")
    st.sidebar.title(f"Welcome {name}")
    st.sidebar.write("GeneSys AI is your all-in-one genomics and bioinformatics companion, designed to empower a wide range of users, from bioinformaticians and researchers to clinicians and students.")

    st.sidebar.title(f"Instructions")
    st.sidebar.write(
        """
        1. Upload file.
        2. Ask a question about your data.
        3. Get an answer.
    """
    )
    logging.basicConfig(
        format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)
    load_dotenv()

    def determine_file_type(file):
        if file is not None:
            file_extension = file.name.split('.')[-1].lower()
            if file_extension == "fasta":
                return "FASTA"
            elif file_extension == "csv":
                return "CSV"
            elif file_extension == "pdb":
                return "PDB"
        return None

    file = st.file_uploader("Drop a file here")

    data_type = determine_file_type(file)

    if data_type:
        st.write(f"What would you like to do?")
    else:
        st.write("Upload a file to get started.")
    # Another reason we don't want this to be our implementation is we are going to hopefully have 10 --> 20 file types in the future. Which will be a lot of ifs.

    if data_type == "FASTA":
        if data_type == "FASTA":
            fasta_file = file

            if fasta_file is not None:
                temp_file_path = os.path.join("/tmp", fasta_file.name)
                with open(temp_file_path, "wb") as temp_file:
                    temp_file.write(fasta_file.read())

                st.success(f"File uploaded successfully!")

        if fasta_file is not None:

            # url = get_s3_url(filename=fasta_file.name)
            # upload_content_to_s3(url, fasta_content)

            user_input = st.chat_input("What is up?")

            col1, col2 = st.columns(2)

            if user_input:
                st.write(run_conversation(user_input, temp_file_path))

            sequence_type = sequence_type(temp_file_path)

            with col1:
                msa_button = None
                orf_button = None
                restriction_button = None
                translate_button = None
                gc_button = None
                isoelectric_button = None
                if sequence_type == "DNA":
                    orf_button = st.button("What are the open reading frames?")
                    restriction_button = st.button("Find restriction sites on the first sequence")
                elif sequence_type == "RNA":
                    translate_button = st.button("Translate the given sequence")
                    gc_button = st.button("Calculate the GC content(s)")
                elif sequence_type == "Protein":
                    msa_button = st.button("Perform MSA")
                    isoelectric_button = st.button("Calculate isoelectric points")

            with col2:
                mass_button = None
                transcription_button = None
                gc_button = None
                reverse_button = None
                complement_button = None
                phylogenetic_button = None

                if sequence_type == "DNA":
                    mass_button = st.button("What is the mass of the given sequence?")
                    transcription_button = st.button("Generate the mRNA transcript")
                elif sequence_type == "RNA":
                    reverse_button = st.button("Reverse the given sequence")
                    complement_button = st.button("Generate the complement(s)")
                elif sequence_type == "Protein":
                    mass_button = st.button("What is the mass of the given sequence?")
                    phylogenetic_button = st.button("Generate a phylogenetic tree")
                    

            if msa_button:
                st.write(run_conversation("Perform MSA on the given FASTA file", temp_file_path))
            elif mass_button:
                st.write(run_conversation("What is the mass of the given sequence?", temp_file_path))
            elif orf_button:
                st.write(run_conversation("What are the ORFs for the given file?", temp_file_path))
            elif restriction_button:
                st.write(run_conversation("What are restriction sites on the first sequence?", temp_file_path))
            elif transcription_button:
                st.write(run_conversation("Generate the mRNA transcript", temp_file_path))
            elif translate_button:
                st.write(run_conversation("Translate the given sequence", temp_file_path))
            elif gc_button:
                st.write(run_conversation("Calculate the GC content(s)", temp_file_path))
            elif reverse_button:
                st.write(run_conversation("Reverse the given sequence", temp_file_path))
            elif complement_button:
                st.write(run_conversation("Generate the complement of the given sequence", temp_file_path))
            elif isoelectric_button:
                st.write(run_conversation("What is the isoelectric point of the given sequence?", temp_file_path))
            elif phylogenetic_button:
                st.write(run_conversation("Generate a phylogenetic tree", temp_file_path))
            else:
                st.write("Please upload a FASTA file.")

    elif data_type == "PDB":
        pdb_file = file

        if pdb_file is not None:
            pdb_content = pdb_file.read().decode("utf-8")
            url = get_s3_url(filename=pdb_file.name)
            upload_content_to_s3(url, pdb_content)

            pdb_user_input = st.chat_input("What is up?")

            if pdb_user_input:
                st.write(render_protein_file(pdb_content))

        else:
            st.write("Please upload a PDB file.")

    elif data_type == "CSV":
        csv_file = file

        if csv_file is not None:
            st.write("CSV file uploaded. Displaying DataFrame:")
            df = pd.read_csv(csv_file)

            url = get_s3_url(filename=csv_file.name)
            upload_content_to_s3(url, df.to_csv(StringIO()))

            st.dataframe(df)

            llm = OpenAI(api_token=os.getenv("OPEN_API_KEY"))
            sdf = SmartDataframe(df, config={"llm": llm})

            csv_user_input = st.chat_input("What is up?")

            url = get_s3_url(
                filename=f"ChatSession/query-{str(int(time()))}.txt")
            upload_content_to_s3(url, csv_user_input)

            col1, col2 = st.columns(2)

            with col1:
                shape_button = st.button("What is the shape of the dataframe?")
                columns_button = st.button("What are the columns of the dataframe?")
            with col2:
                describe_button = st.button("Tell me the descriptive statistics")
                head_button = st.button("Show the first 5 rows of the dataframe")

            if shape_button:
                csv_user_input = "What is the shape of the dataframe?"
                st.write(csv_user_input)
            elif columns_button:
                csv_user_input = "What are the columns of the dataframe?"
                st.write(csv_user_input)
            elif describe_button:
                csv_user_input = "Tell me the descriptive statistics of the dataframe"
                st.write(csv_user_input)
            elif head_button:
                csv_user_input = "Show the first 5 rows of the dataframe"
                st.write(csv_user_input)

            if csv_user_input:
                st.write(csv_user_input)

                response = sdf.chat(csv_user_input)
                st.write(response)

                if response == None:
                    st.image("exports/charts/temp_chart.png",
                             caption="Chart Image", use_column_width=True)
            