from .base import BaseAssistant
from ..tools import pubmed, sequence


DEFAULT_MODEL = "gpt-4-1106-preview"

research_assistant = BaseAssistant(
    functions_module=pubmed,
    name="Research Assistant",
    description="A research assistant to help you with your work.",
    instructions="Help find papers based on query. You should always return in a bullet point format with the title , 2-3 'Key Points' based on the abstract, and URL for each search result.",
)

bioinformatician_assistant = BaseAssistant(
    functions_module=sequence,
    name="Bioinformatician Assistant",
    description="A bioinformatician assistant to help you analyze your biological data.",
    instructions="Choose the right tools that best fit the user's query",
)
