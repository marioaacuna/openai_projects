import os
import json
import requests

'''
This code will search articles in scopus and download json files containing full text (when found) and abstracts 
related to chronic pain
The number pf search items can be modified
'''


def run(query):
    # create folders if they don't exist
    global total_results
    target_folder_json_abstract = 'download/scopus/JSON/abstracts'
    target_folder_json_full = 'download/scopus/JSON/full'
    # target_folder_pdf = 'download/pdf'
    os.makedirs(target_folder_json_abstract, exist_ok=True)
    os.makedirs(target_folder_json_full, exist_ok=True)

    # replace with your own API key
    with open(
            '/Users/marioacuna/Library/Mobile Documents/com~apple~CloudDocs/_todo/ChronicPainLLMs/elsevier_api_key.txt',
            'r') as file:
        api_key = file.read()

    BASE_URL = "https://api.elsevier.com/content/search/scopus"
    start = 0
    count = 25
    max_results = 100
    params = {
        "query": f"{query} ",
        "apiKey": api_key,
        "view": "COMPLETE",
        "count": count,
        "start": start,
        "sort": "date",
    }

    response = requests.get(BASE_URL, params=params)
    if response.status_code == 200:
        total_results = int(response.json()["search-results"]["opensearch:totalResults"])
        print(f"Total number of results: {total_results}")
    else:
        print(f"Error: {response.status_code} - {response.text}")
        exit()

    while start < total_results and start < max_results:
        params["start"] = start
        response = requests.get(BASE_URL, params=params)
        results = response.json()["search-results"]["entry"]
        for result in results:
            title = result["dc:title"]
            authors = [author["authname"] for author in result.get("author", [])]
            publication = result.get("prism:publicationName", "Unknown")
            doi = result.get("prism:doi", "N/A")
            description = result['dc:description']
            print(doi)

            # download the JSON with the full text
            json_headers = {
                "X-ELS-APIKey": api_key,
                "Accept": "application/json"
            }
            json_url = f"https://api.elsevier.com/content/article/doi/{doi}?view=FULL"
            json_response = requests.get(json_url, headers=json_headers)
            print(json_response.status_code)
            json_data = json_response.json()

            # New version 2
            jsonfile = doi.replace('/', '_')

            # check if full-text is available
            if 'full-text-retrieval-response' in json_data:
                # extract full-text from the json data
                main_text = json_data['full-text-retrieval-response']['originalText']
                jsonfile_full = os.path.join(target_folder_json_full, jsonfile)
                # create a dictionary with the desired fields
                json_dict = {'title': title, 'doi': doi, 'main_text': main_text}

                # convert dictionary to json string
                json_file = json.dumps(json_dict)

                # write it
                with open(f"{jsonfile_full}.json", "w") as f:
                    f.write(json.dumps(json_file))

            # extract abstract from the json data
            jsonfile_abs = os.path.join(target_folder_json_abstract, jsonfile)

            # create a dictionary with the desired fields
            json_dict = {'title': title, 'doi': doi, 'main_text': description}

            # convert dictionary to json string
            json_file = json.dumps(json_dict)

            # write it
            with open(f"{jsonfile_abs}.json", "w") as f:
                f.write(json.dumps(json_file))
        start += count


if __name__ == '__main__':
    q = 'TITLE-ABS-KEY ( "chronic pain" OR "neuropathic pain" AND mice  )  AND  PUBYEAR  >  2010  AND  PUBYEAR  >  ' \
        '2010'
    run(q)
