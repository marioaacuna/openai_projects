import os
import requests
from pymed import PubMed
from bs4 import BeautifulSoup
from urllib.parse import urljoin


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
        article_doi = article_doi.split('\n')[0]
        article_pmid = article.pubmed_id.split('\n')[0]
        # Construct the DOI URL
        doi_url = f"https://doi.org/{article_doi}"

        # Get the DOI webpage
        doi_response = requests.get(doi_url)
        real_url = doi_response.url

        # Check if real_url redirects to another url
        while real_url != doi_url:
            doi_url = real_url
            doi_response = requests.get(doi_url)
            real_url = doi_response.url

        # Parse the HTML of the DOI webpage using BeautifulSoup
        doi_soup = BeautifulSoup(doi_response.text, "html.parser")

        # Find the link to the PDF file
        pdf_links = doi_soup.select("a[href$='.pdf']")
        if pdf_links:
            pdf_link = pdf_links[0]['href']
            pdf_url = pdf_link if pdf_link.startswith('http') else urljoin(real_url, pdf_link)

            # Download the PDF file
            response = requests.get(pdf_url, allow_redirects=False)
            if response.status_code == 200 or response.status_code == 303:
                # Create the download directory if it doesn't exist
                download_dir = "downloads"
                if not os.path.exists(download_dir):
                    os.makedirs(download_dir)

                # Save the PDF file to the downloads directory
                filename = f"{article_year}_{article_pmid}.pdf"
                with open(os.path.join(download_dir, filename), "wb") as f:
                    f.write(response.content)

                print(f"{article_doi} Downloaded")
            else:
                print(f"{article_doi} Skipped: Failed to download PDF")
        else:
            print(f"{article_doi} Skipped: No PDF link found")

    except Exception as e:
        print(f"{article_doi} Error: {str(e)}")
