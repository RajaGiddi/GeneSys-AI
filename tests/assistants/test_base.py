from genesys.assistants.base import BaseAssistant
from genesys.openai import openai_client as client

def test_assistant_without_functions():
    test_assistant = BaseAssistant(
        name="Test Assistant",
        description="Used for testing.",
        metadata={"test": True}
    )
    also_test_assistant = client.beta.assistants.retrieve(test_assistant.id)
    assert test_assistant.id == also_test_assistant.id

def test_assistant_with_functions():
    from time import time
    from genesys.tools import pubmed

    test_assistant = BaseAssistant(
        name="Another Test Assistant",
        description="An assistant with functions.",
        functions_module=pubmed,
        metadata={"test": True, "time": time()}
    )

    retrieved_assistant = client.beta.assistants.retrieve(test_assistant.id)
    print(retrieved_assistant.model_dump_json(indent=2))
    assert retrieved_assistant.id == test_assistant.id

def test_assistant_with_multiple_mods():
    from time import time
    from genesys.tools import pubmed, sequence

    test_assistant = BaseAssistant(
        name="Another Test Assistant",
        description="An assistant with a lot of functions.",
        functions_modules=[pubmed, sequence],
        metadata={"test": True, "time": time()}
    )

    retrieved_assistant = client.beta.assistants.retrieve(test_assistant.id)
    print(retrieved_assistant.model_dump_json(indent=2))
    assert retrieved_assistant.id == test_assistant.id
