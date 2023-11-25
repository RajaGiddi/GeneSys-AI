import streamlit as st
from genesys.assistants import bioinformatician_assistant
from genesys.openai import openai_client as client
import time
import tempfile
import os

st.title("Bioinformatics")
st.write("Upload your data to perform bioinformatics analysis")

if "thread" not in st.session_state:
    st.session_state.thread = client.beta.threads.create()

uploaded_file = st.file_uploader("Upload your .fasta file", type=["fasta"])
if uploaded_file is not None:
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as tmp_file:
        tmp_file.write(uploaded_file.getvalue())
        st.session_state.fasta_file_path = tmp_file.name

for message in client.beta.threads.messages.list(st.session_state.thread.id).data:
    with st.chat_message(message.role):
        st.markdown(message.content[0].text.value)

if prompt := st.chat_input("Say something"):
    with st.chat_message("user"):
        st.markdown(prompt)
        client.beta.threads.messages.create(
            st.session_state.thread.id,
            role="user",
            content=prompt + tmp_file.name,
        )
        

    with st.chat_message("assistant"):
        run = client.beta.threads.runs.create(
            thread_id=st.session_state.thread.id,
            assistant_id=bioinformatician_assistant.id,
        )
        while True:
            time.sleep(1)
            run  = client.beta.threads.runs.retrieve(
                run_id=run.id, thread_id=run.thread_id
            )
            
            if run.status == "requires_action" and run.required_action is not None:
                tool_outputs = bioinformatician_assistant.get_tool_outputs(run)
                res = client.beta.threads.runs.submit_tool_outputs(
                    run_id=run.id,
                    thread_id=run.thread_id,
                    tool_outputs=tool_outputs,
                )
            elif run.status == "completed":
                break

        messages = client.beta.threads.messages.list(thread_id=run.thread_id).data
            
        st.markdown(messages[-1].content[0].text.value)