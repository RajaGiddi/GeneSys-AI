from .base import BaseAssistant
from ..tools import pubmed, sequence


DEFAULT_MODEL = "gpt-4-1106-preview"

research_assistant = BaseAssistant(
    functions_module=pubmed,
    name="Research Assistant",
    description="A research assistant to help you with your work.",
    instructions="""Help find papers based on query. When fetching papers, display results in 
        a bullet point format with the title , 2-3 'Key Points' based on the abstract, 
        and the PMID, PMCID, and DOI urls.""",
)

bioinformatician_assistant = BaseAssistant(
    functions_module=sequence,
    name="Bioinformatician Assistant",
    description="A bioinformatician assistant to help you analyze your biological data.",
    instructions="Choose the right tools that best fit the user's query",
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
