from elsapy.elsclient import ElsClient
from elsapy.elsprofile import ElsAuthor, ElsAffil
from elsapy.elsdoc import FullDoc, AbsDoc
from elsapy.elssearch import ElsSearch
import json
from Bio import Entrez
import requests
import os
import httpx
from bs4 import BeautifulSoup
from urllib.parse import urljoin


import time


def run(api_key, download_dir_full, download_dir):
    # Provide your email address to PubMed
    Entrez.email = 'acuna.mario@gmail.com'

    # Define the search terms
    search_terms = 'Chronic pain mice brain treatment'
    search_filter = '("2022/01/01"[Date - Create] : "2023/03/03"[Date - Create]) AND ("review"[Publication Type] OR "article"[Publication Type])'

    # Search PubMed
    handle = Entrez.esearch(db='pubmed', term=search_terms, retmax=5, mindate='2022/01/01', maxdate='2023/04/04',
                            datetype='pdat', retmode='xml', sort='date', usehistory='y', idtype='acc')
    results = Entrez.read(handle)

    # Fetch the abstracts, titles, and DOIs of the articles
    id_list = results['IdList']
    batch_size = 5  # 100
    start = 0
    end = batch_size

    while start < len(id_list):
        id_batch = id_list[start:end]
        fetch_handle = Entrez.efetch(db='pubmed', id=id_batch, rettype='abstract', retmode='xml')
        fetch_results = Entrez.read(fetch_handle)
        try:
            for article in fetch_results['PubmedArticle']:
                pmid = article['MedlineCitation']['PMID']
                title = article['MedlineCitation']['Article']['ArticleTitle']
                abstract = article['MedlineCitation']['Article']['Abstract']['AbstractText']
                doi = next(
                    (id for id in article['PubmedData']['ArticleIdList'] if id.attributes.get('IdType') == 'doi'), None)
                article_year = article['MedlineCitation']['Article']['ArticleDate'][0]['Year']
                # Save the abstract to a text file
                # Create the download

                filename = f"{article_year}_{pmid}_abstract.txt"
                # if exist(filename):  # check this syntax ## TODO
                #   continue

                with open(os.path.join(download_dir, filename), 'w', encoding='UTF-8') as f:
                    f.write(f'Title: {title}\n')
                    f.write(f'pmid: {pmid}\n')
                    f.write(f'DOI: {doi}\n')
                    f.write(f'Abstract: {abstract}\n')
                # print(f'pmid: {pmid}\n')

                # Download the PDF file associated with the article
                if doi:

                    # y = download_full('10.1016/j.solmat.2021.111326', api_key)
                    y = download_elsevier(doi, api_key)
                    if y.status_code == 200:
                        json_acceptable_string = y.text
                        d = json.loads(json_acceptable_string)
                        # Print document
                        print(f"{d['full-text-retrieval-response']['coredata']['dc:description']} \n")
                        # Print document
                        print(d['full-text-retrieval-response']['originalText'])

                    if y.status_code == 404:
                        print(f'PMID: {pmid} not found in elsevier')
                        try:
                            # Get the DOI webpage
                            doi_url = f"https://doi.org/{doi}"
                            doi_response = requests.get(doi_url)
                            real_url = doi_response.url
                            # Parse the HTML of the DOI webpage using BeautifulSoup
                            doi_soup = BeautifulSoup(doi_response.text, "html.parser")
                            pdf_links = doi_soup.select("a[href$='pdf']")

                            if pdf_links:
                                pdf_link = pdf_links[0]['href']
                                pdf_url = pdf_link if pdf_link.startswith('http') else urljoin(real_url, pdf_link)
                                #download_manual(article_doi, pdf_url, article_pmid, article_year)
                                print('found fast manual pdf')
                            else:
                                print('maybe try slow pdf')
                                #citation_pdf_url = doi_soup.find("meta", attrs={"name": "citation_pdf_url"})
                                #download_manual

                        except:
                            print('')


                    ## TODO: set the download
                else:
                    print('')

            start += batch_size
            end += batch_size
        except Exception as e:
            print(f"{pmid} Error: {str(e)}")


def download_elsevier(doi, key):
    # Try with url
    HEADERS = {
        'X-ELS-APIKey': key,
        "Accept": 'application/json'
        # 'Accept': 'text / xml'
        # 'Accept': 'text / xml'
        # 'Accept': 'text / html'
    }

    timeout = httpx.Timeout(10.0, connect=60.0)
    client = httpx.Client(timeout=timeout, headers=HEADERS)
    # query="&view=FULL"
    url = f"https://api.elsevier.com/content/search/scopus"
    url = f"https://api.elsevier.com/content/article/doi/{doi}"
    #f"https://api.elsevier.com/content/article/pubmed_id/{pubmed_id}" #28791431
    r = client.get(url,follow_redirects=True)
    # print(response)
    # if response.status_code == 200:
    #    r = requests.get(url, stream=True, headers=HEADERS)
    # else:
    #    r = False
    #    print(f"Did not find response for: {doi} \n Status: {response.status.code}")
    return r

def download_manual():
    print('')

def seach_scopus():

    # replace with your own API key
    API_KEY = api_key

    # base URL for the Scopus Search API
    BASE_URL = "https://api.elsevier.com/content/search/scopus"

    # search query
    query = "all(\"machine learning\")"

    # request parameters
    params = {
        "query": query,
        "apiKey": API_KEY,
        "view": "COMPLETE",
        "count": 25,
        "start": 0,
        "sort": "relevance"
    }

    # send the API request
    response = requests.get(BASE_URL, params=params)

    # check if the request was successful
    if response.status_code == 200:
        # extract the results from the response
        results = response.json()["search-results"]["entry"]
        for result in results:
            # extract the metadata from the result
            title = result["dc:title"]
            authors = [author["authname"] for author in result.get("author", [])]
            publication = result.get("prism:publicationName", "Unknown")
            doi = result.get("prism:doi", "N/A")
            print(f"Title: {title}\nAuthors: {authors}\nPublication: {publication}\nDOI: {doi}\n")
    else:
        print(f"Error: {response.status_code} - {response.text}")


'''
    response = requests.get(url, stream=True, headers=HEADERS, allow_redirects=True)
    if response.status_code == 200:
        r = requests.get(url, stream=True, headers=HEADERS)
    else:
        r = False
        print(f"Did not find response for: {doi} \n Status: {response.status.code}")
    return r
'''

# Get document


if __name__ == '__main__':
    ## Load the api from text
    with open(
            '/Users/marioacuna/Library/Mobile Documents/com~apple~CloudDocs/_todo/ChronicPainLLMs/elsevier_api_key.txt',
            'r') as file:
        api_key = file.read().replace('\n', '')
    target_folder_full = 'download/full_text'
    target_folder_abstract = 'download/abstract'
    if not os.path.exists(target_folder_full):
        os.makedirs(target_folder_full)
    if not os.path.exists(target_folder_abstract):
        os.makedirs(target_folder_abstract)

    run(api_key, target_folder_full, target_folder_abstract)
