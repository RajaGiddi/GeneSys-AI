import streamlit as st
from genesys.assistants import research_assistant
from genesys.openai import openai_client as client
import time

st.title("Research")

if "research_thread" not in st.session_state:
    st.session_state.research_thread = client.beta.threads.create()

for message in client.beta.threads.messages.list(
    st.session_state.research_thread.id, order="asc"
).data:
    with st.chat_message(message.role):
        st.markdown(message.content[0].text.value)

if prompt := st.chat_input("Say something"):
    with st.chat_message("user"):
        st.markdown(prompt)
        client.beta.threads.messages.create(
            st.session_state.research_thread.id,
            role="user",
            content=prompt,
        )
        

    with st.chat_message("assistant"):
        run = client.beta.threads.runs.create(
            thread_id=st.session_state.research_thread.id,
            assistant_id=research_assistant.id,
        )

        with st.spinner('Waiting for a response...'):
            while True:
                time.sleep(1)
                run = client.beta.threads.runs.retrieve(
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

            res = client.beta.threads.messages.list(thread_id=run.thread_id, limit=1)

            if (content := res.data[0].content[0]).type == "text":
                st.markdown(content.text.value)
