from typing import Annotated
import requests
import xml.etree.ElementTree as ET

from typing_extensions import Doc

base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


def _parse_search_results(search_result: str) -> list[str]:
    """
    Parses the XML response from PubMed's ESearch API and returns a list of PMIDs.
    """
    root = ET.fromstring(search_result)
    return [id_tag.text for id_tag in root.findall('.//IdList/Id')]


def fetch_papers(
    query: Annotated[str, Doc("Search query to submit to PubMed.")],
    max_results: int = 10
):
    """Search PubMed and get the title, abstract, URL, and PMCID for each search result."""

    search_url = f"{base_url}esearch.fcgi?db=pubmed&term={query}&retmode=xml"
    search_response = requests.get(search_url)
    id_list = _parse_search_results(search_response.text)

    ids = ','.join(id_list[:max_results])
    fetch_url = f"{base_url}efetch.fcgi?db=pubmed&id={ids}&retmode=xml"
    fetch_response = requests.get(fetch_url)
    root = ET.fromstring(fetch_response.text)

    papers = []
    for article in root.findall('.//PubmedArticle'):
        pmid = article.find('.//PMID').text
        title = article.find('.//ArticleTitle').text
        abstract = article.find('.//Abstract/AbstractText').text

        # Extracting PMCID
        pmcid = None
        doi = None
        article_ids = article.findall('.//ArticleIdList/ArticleId')
        for aid in article_ids:
            if aid.get('IdType') == 'pmc':
                pmcid = aid.text
            elif aid.get('IdType') == 'doi':
                doi = aid.text

        pmid_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
        pmcid_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmcid}/"
        doi_url = f"https://doi.org/{doi}"

        papers.append({
            "pmid_url": pmid_url,
            "pmcid_url": pmcid_url,
            "doi_url": doi_url,
            "title": title,
            "abstract": abstract,
        })

    return papers
