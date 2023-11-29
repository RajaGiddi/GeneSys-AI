"""
This testing module is to make sure that the prompts and the function
descriptions defined in `tools.pubmed` results in the behavior that we expect.
"""

from genesys.tools.pubmed import fetch_papers, search_pmc
import json
import time

def test_perform_search():
    time.sleep(1)
    query = "genomics"
    result = fetch_papers(query=query)

    assert isinstance(result, list)
    for item in result:
        assert isinstance(item, dict)
        assert "pubmed_url" in item
        assert "pubmed_id" in item
        assert "title" in item
        assert "abstract" in item

def test_search_pmc():
    result = search_pmc("cancer")
    print(json.dumps(result, indent=2))
