import inspect
import math
from typing import Annotated

from typing_extensions import Doc

from genesys import utils

def f(x: float, y: float) -> float:
    """A simple function without documentation in annotated metadata."""
    return x + y

def greet(name: str = "friend") -> str:
    """Returns a friendly greeting."""
    return f"Hello there {name}!"

def pythagoras(
    a: Annotated[float, Doc("The length of side A.")],
    b: Annotated[float, Doc("The length of side B.")],
) -> Annotated[float, Doc("The length of side C.")]:
    """Calculate the length of the hypotenuse of a right triangle."""
    return math.sqrt(math.pow(a, 2) + math.pow(b, 2))

def test_get_function_description():
    schema = utils.gen_function_schema(pythagoras)
    assert (
        schema.get("description")
        == "Calculate the length of the hypotenuse of a right triangle."
    )


def test_get_function_name():
    schema = utils.gen_function_schema(pythagoras)
    assert schema.get("name") == "pythagoras"

def test_annotated_param():
    sign = inspect.signature(pythagoras)
    assert utils.is_annotated(sign.parameters["a"].annotation)

def test_get_param_metadata():
    schema = utils.gen_function_schema(pythagoras)
    assert schema["parameters"]["properties"]["a"]["type"] == "number"
    assert schema["parameters"]["properties"]["b"]["type"] == "number"
    assert (
        schema["parameters"]["properties"]["a"]["description"]
        == "The length of side A."
    )
    assert (
        schema["parameters"]["properties"]["b"]["description"]
        == "The length of side B."
    )

def test_get_simple_function_schema():
    schema = utils.gen_function_schema(f)
    assert schema["name"] == "f"
    assert (
        schema["description"]
        == "A simple function without documentation in annotated metadata."
    )
    required_params: list = schema["parameters"]["required"]
    assert len(required_params) == 2
    assert "x" in required_params
    assert "y" in required_params


def test_get_simple_param_annotations():
    schema = utils.gen_function_schema(f)
    assert schema["parameters"]["properties"]["x"]["type"] == "number"
    assert schema["parameters"]["properties"]["y"]["type"] == "number"
    assert schema["parameters"]["properties"]["x"].get("description") is None
    assert schema["parameters"]["properties"]["y"].get("description") is None

def test_get_default_param_value():
    schema = utils.gen_function_schema(greet)
    assert schema["parameters"]["properties"]["name"]["type"] == "string"
    assert schema["parameters"]["properties"]["name"]["default"] == "friend"

def test_gen_tools_schema():
    from genesys.tools import pubmed

    schema = utils.gen_tools_schema(pubmed)
    assert isinstance(schema, list)
    
    for tool in schema:
        assert tool["type"] == "function"

def test_gen_sequence_tools_schema():
    import json
    from genesys.tools import sequence

    schema = utils.gen_tools_schema(sequence)
    print(json.dumps(schema, indent=2))
    assert isinstance(schema, list)

    for tool in schema:
        assert tool["type"] == "function"
