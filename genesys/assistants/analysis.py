from genesys.assistants.base import *
import os

# Instantiate the ResearchAssistant object
assistant = ResearchAssistant(api_key=os.environ['OPENAI_API_KEY'])

user_input= input("What do you need help with? ")

assistant_name = "BioinformaticsGenie"

instructions = "Create a Python function to perform Multiple Sequence Alignment (MSA) using Biopython. The function should be well-commented, handle errors gracefully, and be efficient in handling multiple sequences. Ensure the function is easy to understand and use."

thread_instructions = "You need to create a Python function that performs Multiple Sequence Alignment (MSA). The function should accept a list of sequences in FASTA format and return the aligned sequences. Use Biopython's AlignIO and MultipleSeqAlignment modules for the alignment process. Make sure the code is clean, well-commented, and includes error handling. Additionally, you are to output ONLY the generated code, no additional text."

file_path = "genesys/assistants/msa.py"

# Call the create_assistant_and_run method
assistant.create_assistant_and_run(user_input=user_input, assistant_name=assistant_name, instructions=instructions, tools_list=[], thread_instructions=thread_instructions, file_path=file_path)
