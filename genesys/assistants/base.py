import json
from types import ModuleType
from typing import NotRequired, Unpack, cast

from openai.types.beta.assistant_create_params import AssistantCreateParams
from openai.types.beta.threads.run import Run
from openai.types.beta.assistant import Assistant
from openai.types.beta.assistant_update_params import AssistantUpdateParams

from ..utils import gen_tools_schema
from ..openai import openai_client as client

class _BaseAssistantInitParams(AssistantCreateParams):
    model: NotRequired[str]
    functions_module: NotRequired[ModuleType]
    """A module containing functions for the assistant to call."""

class BaseAssistant:
    model = "gpt-4-1106-preview"

    def __init__(self, **kwargs: Unpack[_BaseAssistantInitParams]):
        kwargs.setdefault("model", self.model)

        tools = kwargs.setdefault("tools", [])

        if mod := kwargs.pop("functions_module", None):
            self.functions = mod.__dict__
            tools.extend(gen_tools_schema(mod))

        self.assistant = client.beta.assistants.create(
            **cast(AssistantCreateParams, kwargs)
        )
        
    def __getattr__(self, name: str):
        return getattr(self.assistant, name)

    def __repr__(self) -> str:
        return f"<BaseAssistant: {self.assistant.id}>"
    
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, BaseAssistant):
            return NotImplemented
        return self.assistant.id == other.assistant.id

    def __hash__(self) -> int:
        return hash(self.assistant.id)

    def delete(self):
        return client.beta.assistants.delete(self.assistant.id)

    def update(self, **kwargs: Unpack[AssistantUpdateParams]) -> Assistant:
        return client.beta.assistants.update(self.assistant.id, **kwargs)

    def get_tool_outputs(self, run: Run):
        assert(
            run.status == "requires_action"
            and run.required_action is not None
            and run.required_action.type == "submit_tool_outputs"
        )
        
        tool_outputs = []

        for action in run.required_action.submit_tool_outputs.tool_calls:
            if (fn_name := action.function.name) in self.functions:
                function_to_call = self.functions[fn_name]
                args = json.loads(action.function.arguments)
                ret = function_to_call(**args)
                tool_outputs.append({
                    "tool_call_id": action.id,
                    "output": json.dumps(str(ret))
                })
            else:
                raise ValueError(f"Unknown tool: {fn_name}")
                
        return tool_outputs
            
