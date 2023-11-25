import pytest
from genesys.openai import openai_client as client
from genesys.assistants import bioinformatician_assistant
import time

def test__ask_MSA():

    run = client.beta.threads.create_and_run(
        assistant_id=bioinformatician_assistant.id,
        thread={
            "messages": [
                {"role": "user",
                 "content": f"Can you perform MSA on tests/fixtures/covid_sequences.fasta ?",}
            ]
        }
    )

    while True:
        time.sleep(1)
        run = client.beta.threads.runs.retrieve(
            run_id=run.id, thread_id=run.thread_id)

        if run.status == "requires_action" and run.required_action is not None:
            tool_outputs = bioinformatician_assistant.get_tool_outputs(run)
            res = client.beta.threads.runs.submit_tool_outputs(
                run_id=run.id, thread_id=run.thread_id, tool_outputs=tool_outputs)
        elif run.status == "completed":
            break

    # Retrieve and assert messages
    messages = client.beta.threads.messages.list(thread_id=run.thread_id)
    assert run.status == "completed"
    print(messages.model_dump_json(indent=2))

def test_ask_restriction_sites():
    run = client.beta.threads.create_and_run(
        assistant_id=bioinformatician_assistant.id,
        thread={
            "messages": [
                {"role": "user",
                 "content": f"Can you find the restriction sites on the 2nd sequence tests/fixtures/sequence.fasta ?",}
            ]
        }
    )

    while True:
        time.sleep(1)
        run = client.beta.threads.runs.retrieve(
            run_id=run.id, thread_id=run.thread_id)

        if run.status == "requires_action" and run.required_action is not None:
            tool_outputs = bioinformatician_assistant.get_tool_outputs(run)
            res = client.beta.threads.runs.submit_tool_outputs(
                run_id=run.id, thread_id=run.thread_id, tool_outputs=tool_outputs)
        elif run.status == "completed":
            break

    # Retrieve and assert messages
    messages = client.beta.threads.messages.list(thread_id=run.thread_id)
    assert run.status == "completed"
    print(messages.model_dump_json(indent=2))

def test_ask_isoelectric_points():
    run = client.beta.threads.create_and_run(
        assistant_id=bioinformatician_assistant.id,
        thread={
            "messages": [
                {"role": "user",
                 "content": f"Can you find the isoelectric points on tests/fixtures/rand_gen.fasta ?",}
            ]
        }
    )

    while True:
        time.sleep(1)
        run = client.beta.threads.runs.retrieve(
            run_id=run.id, thread_id=run.thread_id)

        if run.status == "requires_action" and run.required_action is not None:
            tool_outputs = bioinformatician_assistant.get_tool_outputs(run)
            res = client.beta.threads.runs.submit_tool_outputs(
                run_id=run.id, thread_id=run.thread_id, tool_outputs=tool_outputs)
        elif run.status == "completed":
            break

    # Retrieve and assert messages
    messages = client.beta.threads.messages.list(thread_id=run.thread_id)
    assert run.status == "completed"
    print(messages.model_dump_json(indent=2))