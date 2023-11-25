from .base import BaseAssistant
from ..tools import pubmed, sequence
from ..utils import gen_tools_schema


DEFAULT_MODEL = "gpt-4-1106-preview"

research_assistant = BaseAssistant(
    tools_module=pubmed,
    name="Research Assistant",
    description="A research assistant to help you with your work.",
    instructions="Help find papers based on query. You should always return ",
    tools=gen_tools_schema(pubmed),
)

bioinformatician_assistant = BaseAssistant(
    tools_module=sequence,
    name="Bioinformatician Assistant",
    description="A bioinformatician assistant to help you analyze your biological data.",
    instructions="Choose the right tools that best fit the user's query",
    tools=gen_tools_schema(sequence),
)
