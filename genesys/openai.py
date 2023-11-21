import os

from openai import OpenAI

from .env import load_dotenv

load_dotenv()

openai_client = OpenAI(
    api_key=os.environ["OPENAI_API_KEY"]
)
