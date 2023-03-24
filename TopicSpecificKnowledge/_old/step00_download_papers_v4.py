from Bio import Entrez
import requests
import os

# Provide your email address to PubMed
Entrez.email = 'acuna.mario@gmail.com'

# Define the search terms
search_terms = 'Chronic pain mice brain treatment'
search_filter = '("2022/01/01"[Date - Create] : "2023/03/03"[Date - Create]) AND ("review"[Publication Type] OR "article"[Publication Type])'

# Search PubMed
handle = Entrez.esearch(db='pubmed', term=search_terms, retmax=100000, mindate='2010/01/01', maxdate='3000',
                        datetype='pdat', retmode='xml', sort='date', usehistory='y', idtype='acc')
results = Entrez.read(handle)

# Fetch the abstracts, titles, and DOIs of the articles
id_list = results['IdList']
batch_size = 100
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
            doi = next((id for id in article['PubmedData']['ArticleIdList'] if id.attributes.get('IdType') == 'doi'), None)
            article_year = article['MedlineCitation']['DateCompleted']['Year']
            # Save the abstract to a text file
            # Create the download
            download_dir = "downloads/abstracts"
            if not os.path.exists(download_dir):
                os.makedirs(download_dir)

            filename = f"{article_year}_{pmid}_abstract.txt"
            if exist(filename): #check this syntax
                continue
            with open(os.path.join(download_dir, filename), 'w', encoding='UTF-8') as f:
                f.write(f'Title: {title}\n')
                f.write(f'pmid: {pmid}\n')
                f.write(f'DOI: {doi}\n')
                f.write(f'Abstract: {abstract}\n')
            #print(f'pmid: {pmid}\n')


            # Download the PDF file associated with the article
            if False:
                headers = {'User-Agent': 'Mozilla/5.0'}
                url = f'https://doi.org/{doi}'
                response = requests.get(url)
                real_url = response.url
                response = requests.get(real_url, )
                if response.ok:
                    content_type = response.headers.get('Content-Type')
                    if 'pdf' in content_type:
                        with open(f'{pmid}.pdf', 'wb') as f:
                            f.write(response.content)
                    else:
                        # Download the full text from the HTML page
                        from bs4 import BeautifulSoup

                        soup = BeautifulSoup(response.content, 'html.parser')
                        full_text_url = soup.find('meta', {'name': 'citation_fulltext_html_url'})['content']
                        full_text_response = requests.get(full_text_url, headers=headers)
                        if full_text_response.ok:
                            with open(f'{pmid}.html', 'wb') as f:
                                f.write(full_text_response.content)
        start += batch_size
        end += batch_size
    except Exception as e:
        print(f"{pmid} Error: {str(e)}")
