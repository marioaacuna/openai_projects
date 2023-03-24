import pandas as pd
import os
import requests

with open('/Users/marioacuna/Library/Mobile Documents/com~apple~CloudDocs/_todo/ChronicPainLLMs/elsevier_api_key.txt',
          'r') as file:
    api_key = file.read().replace('\n', '')

HEADERS = {
    'X-ELS-APIKEY': api_key,
    'Accept': 'application/json'
}


def get_pdf(doi, file_to_save_to):
    url = f"http://api.elsevier.com/content/article/doi:{doi}?view=FULL"
    with requests.get(url, stream=True, headers=HEADERS) as r:
        if r.status_code == 200:
            for chunk in r.iter_content(chunk_size=1024*1024):
                file_to_save_to.write()





file_path = open(f'downloads/pdfs/{doi}.pdf', 'wb')
file = f'downloads/pdfs/{pdffile}.pdf'
get_pdf(doi, file_path)

doi_list = pd.read_excel('list.xls')
doi_list.columns = ['DOIs']
count = 0
for doi in doi_list['DOIs']:
    doi = doi.replace('DOI:','')
    pdffile = doi.replace('/','_')
    if not os.path.exists(f'path/{pdf}.pdf'):
        with open(f'path/{pdf}.pdf', 'wb') as file:
            get_pdf(doi, file)
        count += 1
        print(f"Dowloaded: {count} of {len(doi_list['DOIs'])} articles")