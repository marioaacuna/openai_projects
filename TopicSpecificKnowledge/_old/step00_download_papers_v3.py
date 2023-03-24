import os
import requests
from pymed import PubMed
from bs4 import BeautifulSoup
from urllib.parse import urljoin


def download_pdf(article_doi, pdf_url, article_pmid, article_year):
    try:
        # Download the PDF file
        response = requests.get(pdf_url, allow_redirects=False)
        if response.status_code == 200 or response.status_code == 302 or response.status_code == 303:
            if response.status_code == 303 or response.status_code == 302:
                new_url = response.headers['Location']
                response = requests.get(new_url)
            # Create the download directory if it doesn't exist
            download_dir = "downloads/pdfs"
            if not os.path.exists(download_dir):
                os.makedirs(download_dir)

            # Save the PDF file to the downloads directory
            filename = f"{article_year}_{article_pmid}.pdf"
            with open(os.path.join(download_dir, filename), "wb") as f:
                f.write(response.content)

            print(f"{article_pmid} Downloaded")
        else:
            print(f"{article_pmid} Skipped: Failed to download PDF")
    except Exception as e:
        print(f"{article_pmid} Error: {str(e)}")


# Set up the search terms
search_terms = "chronic pain neuropathic pain rodents mice"
search_results = []

# Set up the PubMed query
pubmed = PubMed(tool="MyResearch", email="acuna.mario@gmail.com")

# Perform the search and get the results
results = pubmed.query(search_terms, max_results=20)

# Loop through the search results and save the PDF files
for article in results:
    try:
        article_title = article.title
        article_doi = article.doi
        article_pmid = article.pubmed_id
        article_year = article.publication_date.year

        if not article_doi:
            print(f"No DOI found for {article_title}")
            continue

        # it can happen that there are more than 1 doi
        article_doi = article_doi.split('\n')[0]
        article_pmid = article.pubmed_id.split('\n')[0]

        # if not article_doi.startswith("10."):
        #    article_doi = f"10.{article_doi}"

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
        #if article_doi == '10.1089/neu.2022.0482':
        #    print('')
        # Find the link to the PDF file
        pdf_links = doi_soup.select("a[href$='.pdf']")


        if pdf_links:
            pdf_link = pdf_links[0]['href']
            pdf_url = pdf_link if pdf_link.startswith('http') else urljoin(real_url, pdf_link)
            download_pdf(article_doi, pdf_url, article_pmid, article_year)
        else:
            # citation_pdf_url = doi_soup.find("meta", attrs={"name": "citation_pdf_url"})
            # pdf_url = citation_pdf_url.get("content")
            pdf_url = None
            pdf_links = doi_soup.find_all("a", href=True)
            done = False
            for link in pdf_links:
                link_class = link.get("class")
                if link_class:
                    #print(link_class)
                    for c in link_class:
                        if "pdf" in c or "download" in c or 'PDF' in c or 'UD_ArticlePDF' in c:
                            pdf_link = link.get("href")
                            pdf_url = pdf_link if pdf_link.startswith('http') else urljoin(real_url, pdf_link)
                            try:
                                download_pdf(article_doi, pdf_url, article_pmid, article_year)
                                done = True
                            except:
                                done = False

                            if done:
                                break
                if done:
                    break
            #if pdf_url:
            #    download_pdf(article_doi, pdf_url, article_pmid, article_year)

            #if not pdf_url:
        abstract = None
        meta_tags = doi_soup.find_all('meta', attrs={'name': 'DC.Description'})
        content = ''
        for meta_tag in meta_tags:
            content = meta_tag.get('content')

        if not content:
            print(f"{article_doi} no content")
            abstract = article.abstract
            if not abstract:
                print(f"{article_doi} Skipped: No pubmed abstract found")
                continue
        else:
            abstract = content
        print(abstract)
        abstract = article.title

        # Create the download
        download_dir = "downloads/abstracts"
        if not os.path.exists(download_dir):
            os.makedirs(download_dir)

    #    Save the abstract to a text file in the downloads directory
        filename = f"{article_year}_{article_pmid}_abstract.txt"
        with open(os.path.join(download_dir, filename), "w", encoding="utf-8") as f:
            f.write(abstract)

    except Exception as e:
        print(f"{article_doi} Error: {str(e)}")
