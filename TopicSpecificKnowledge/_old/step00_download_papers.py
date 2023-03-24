import os
import requests
import zipfile
from pymed import PubMed
from bs4 import BeautifulSoup

# Set up the search terms
search_terms = "chronic pain neuropathic pain rodents mice"
search_results = []

# Set up the PubMed query
pubmed = PubMed(tool="MyResearch", email="myemail@mydomain.com")

# Perform the search and get the results
results = pubmed.query(search_terms, max_results=5)

# Loop through the search results and save the PDF files
for article in results:
    try:
        article_title = article.title
        article_doi = article.doi
        article_pmid = article.pubmed_id
        article_year = article.publication_date.year

        # it can happen that there are more than 1 doi
        article_doi = article.doi.split('\n')[0]
        # Construct the DOI URL
        doi_url = f"https://doi.org/{article_doi}"

        # Get the DOI webpage
        doi_response = requests.get(doi_url)
        real_url = doi_response.url
        #pdf_url = real_url+'.pdf'
        ## Parse the HTML of the DOI webpage using BeautifulSoup
        #doi_soup = BeautifulSoup(doi_response.text, "html.parser")
        response = requests.get(real_url, allow_redirects=False)
        if response.status_code == 200:
            soup = BeautifulSoup(response.content, "html.parser")
            if len(soup.select("a[href$='.pdf']")) > 0:
                for link in soup.select("a[href$='.pdf']"):
                    url_pdf = 'https://' + real_url.split('/')[2] + link['href']
                    response = requests.get(url_pdf)
            else:
                print('Skipped: ' + article_doi)
                continue

            # Create the download directory if it doesn't exist
            download_dir = "downloads"
            if not os.path.exists(download_dir):
                os.makedirs(download_dir)

            # Save the PDF file to the downloads directory
            filename = f"{article_year}_{article_pmid}.pdf"
            with open(os.path.join(download_dir, filename), "wb") as f:
                f.write(response.content)

            print(article_doi + ' Downloaded')
        else:
            print('Skipped: ' + article_doi)

    except Exception as e:
        print(f"Error: {str(e)}")
