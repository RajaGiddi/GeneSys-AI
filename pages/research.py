import streamlit as st
from genesys.assistants import research_assistant
from genesys.openai import openai_client as client
import time

st.title("Research")

if "thread" not in st.session_state:
    st.session_state.thread = client.beta.threads.create()
    
for message in client.beta.threads.messages.list(st.session_state.thread.id).data[::-1]:
    with st.chat_message(message.role):
        st.markdown(message.content[0].text.value)

if prompt := st.chat_input("Say something"):
    with st.chat_message("user"):
        st.markdown(prompt)
        client.beta.threads.messages.create(
            st.session_state.thread.id,
            role="user",
            content=prompt,
        )
        

    with st.chat_message("assistant"):
        run = client.beta.threads.runs.create(
            thread_id=st.session_state.thread.id,
            assistant_id=research_assistant.id,
        )
        while True:
            time.sleep(1)
            run  = client.beta.threads.runs.retrieve(
                run_id=run.id, thread_id=run.thread_id
            )
            
            if run.status == "requires_action" and run.required_action is not None:
                tool_outputs = research_assistant.get_tool_outputs(run)
                res = client.beta.threads.runs.submit_tool_outputs(
                    run_id=run.id,
                    thread_id=run.thread_id,
                    tool_outputs=tool_outputs,
                )
            elif run.status == "completed":
                break

        messages = client.beta.threads.messages.list(thread_id=run.thread_id).data
            
        st.markdown(messages[-1].content[0].text.value)