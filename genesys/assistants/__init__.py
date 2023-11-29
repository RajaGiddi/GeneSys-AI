from .base import BaseAssistant
from ..tools import pubmed, sequence


DEFAULT_MODEL = "gpt-4-1106-preview"

research_assistant = BaseAssistant(
    functions_module=pubmed,
    name="Research Assistant",
    description="A research assistant to help you with your work.",
    instructions="""
        Help search for papers relavant to the user's research goal and if
        asked, help users understand the papers by summarizing and explaining
        them.

        When displaying search results, use an ordered list and have each item in the list follow this format:

        - **Title**: the title of the paper
        - **Abstract**:
          - a key points from the abstract
          - use as many bullet points as needed
        - **Links**:
          - any links or URLs
        - **Other**:
          - other information like PMCIDs or DOIs
    """,
)

bioinformatician_assistant = BaseAssistant(
    functions_module=sequence,
    name="Bioinformatician Assistant",
    description="A bioinformatician assistant to help you analyze your biological data.",
    instructions="Choose the right tools that best fit the user's query. Provide your intepretation of the results.",
)

manuscript_assistant = BaseAssistant(
    name="Manuscript Assistant",
    description="A manuscript assistant to help you write your manuscript.",
    instructions="""Help write a manuscript based on the user's query, 
        literature review from the research assistant, and analysis from the 
        bioinformatician assistant. Each manuscript should include an 
        introduction, methods, results, and discussion section all with a
        scientific writing style.""",
)
