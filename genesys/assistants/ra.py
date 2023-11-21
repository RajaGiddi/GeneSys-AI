from openai import OpenAI
import time
import json
from genesys.env import load_dotenv


load_dotenv()

class ResearchAssistant:
    def __init__(self, api_key, model="gpt-4-1106-preview"):
        self.client = OpenAI(api_key=api_key)
        self.model = model

    def create_assistant_and_run(self, user_input, assistant_name, instructions, tools_list, thread_instructions, toolkit, file_path=None):
        assistant = self.client.beta.assistants.create(
            name=assistant_name,
            instructions=instructions,
            tools=tools_list,
            model=self.model
        )
        thread = self.client.beta.threads.create()
        self.client.beta.threads.messages.create(thread_id=thread.id, role="user", content=user_input)
        run = self.client.beta.threads.runs.create(thread_id=thread.id, assistant_id=assistant.id, instructions=thread_instructions)
        self.check_run_status(thread.id, run.id, toolkit)

    def check_run_status(self, thread_id, run_id, toolkit):
        while True:
            time.sleep(5)
            run_status = self.client.beta.threads.runs.retrieve(thread_id=thread_id, run_id=run_id)

            if run_status.status == 'completed':
                messages = self.client.beta.threads.messages.list(thread_id=thread_id)
                for msg in messages.data:
                    role = msg.role
                    content = msg.content[0].text.value
                    print(f"{role.capitalize()}: {content}")
                break
            elif run_status.status == 'requires_action':
                self.process_required_actions(thread_id, run_id, run_status.required_action.submit_tool_outputs.model_dump(), toolkit)
            else:
                print("Waiting for the Assistant to process...")

    def process_required_actions(self, thread_id, run_id, required_actions, toolkit):
        tool_outputs = []

        for action in required_actions["tool_calls"]:
            func_name = action['function']['name']
            arguments = json.loads(action['function']['arguments'])

            try:
                # Dynamically get the function from the specified toolkit
                function_to_call = getattr(toolkit, func_name)
                

                # Call the function with unpacked arguments
                output = function_to_call(**arguments)

                # Convert the output to a string and append it to tool_outputs
                output_str = json.dumps(output)
                tool_outputs.append({
                    "tool_call_id": action['id'],
                    "output": output_str
                })

            except AttributeError:
                raise ValueError(f"Unknown function: {func_name}")
            except json.JSONDecodeError:
                raise ValueError("Error decoding arguments for the function call")

        # Submitting outputs back to the Assistant
        print("Submitting outputs back to the Assistant...")
        self.client.beta.threads.runs.submit_tool_outputs(thread_id=thread_id, run_id=run_id, tool_outputs=tool_outputs)