'''
This script and idea was taken from David Shapiro:
https://github.com/daveshap/ChatGPT_QA_Regenerative_Medicine/blob/main/step02_embed_papers.py
I modified it to read my data
'''

import os
import openai
import json
import pandas as pd
from time import time, sleep


def open_file(filepath):
    with open(filepath, 'r', encoding='utf-8') as infile:
        return infile.read()


def save_file(filepath, content):
    with open(filepath, 'w', encoding='utf-8') as outfile:
        outfile.write(content)


def load_json(filepath):
    with open(filepath, 'r', encoding='utf-8') as infile:
        return json.load(infile)


def save_json(filepath, payload):
    with open(filepath, 'w', encoding='utf-8') as outfile:
        json.dump(payload, outfile, ensure_ascii=False, sort_keys=True, indent=2)


def gpt3_embedding(content, engine='text-embedding-ada-002'):
    # try:
    print('CONTENT TO EMBED:', content)
    content = content.encode(encoding='ASCII', errors='ignore').decode()  # fix any UNICODE errors
    # TODO -  get the timerate limiting fixed
    '''
    import openai
    from tenacity import (
        retry,
        stop_after_attempt,
        wait_random_exponential,
    )  # for exponential backoff
     
    @retry(wait=wait_random_exponential(min=1, max=60), stop=stop_after_attempt(6))
    def completion_with_backoff(**kwargs):
        return openai.Completion.create(**kwargs)
     
    completion_with_backoff(model="text-davinci-003", prompt="Once upon a time,")'''
    response = openai.Embedding.create(input=content, engine=engine)
    vector = response['data'][0]['embedding']  # this is a normal list
    print('VECTOR:', vector)
    return vector


# except:
#    return None


def process_text_files(dir_path, out_path):
    # Create the output directory if it doesn't exist
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    # Loop through all files in the directory
    for filename in os.listdir(dir_path):
        if filename.endswith(".json"):
            # Read in the text file
            with open(dir_path + filename, "r", encoding="utf-8") as f:
                text = f.read()

            # Split the text into pages
            pages = text.split("NEW PAGE")

            # Generate embeddings for each page
            embeddings = [gpt3_embedding(page) for page in pages]

            # Create a dictionary with the filename and the pages and embeddings
            output_dict = {"original_filename": filename,
                           "pages": [{"page_number": i + 1,
                                      "text": page,
                                      "embedding": embedding} for i, (page, embedding) in
                                     enumerate(zip(pages, embeddings))]}

            # Save the dictionary to a JSON file
            # with open(os.path.join(out_path, filename.replace(".txt", ".json")), "w") as f:
            #    json.dump(output_dict, f)
            save_json(os.path.join(out_path, filename.replace(".txt", ".json")), output_dict)


def process_json_files(dir_path, out_path):
    # Create the output directory if it doesn't exist
    if not os.path.exists(out_path):
        os.mkdir(out_path)

    # Loop through all files in the directory
    for filename in os.listdir(dir_path):
        if filename.endswith(".json"):
            # Read in the JSON file
            with open(os.path.join(dir_path, filename), "r") as f:
                json_data = json.load(f)

            # TODO - check if the file exists to not embed it again.
            if os.path.exists(os.path.join(out_path, filename)):
                print(f"{filename} already exists in {out_path}. Skipping...")
                continue

            # Extract the main_text key
            # See first if the input is a string (that happened with the old version of the step01)
            if type(json_data) == str:
                json_data = json.loads(json_data)
                # Load the JSON string into a dictionary

            main_text = json_data["main_text"]

            # TODO: check here the length or alternatively switch to either full or abstract. We need to set a nr of
            #  words and split the text, otherwise it wont work because of the number of tokens

            # Generate embeddings for the main_text
            embedding = gpt3_embedding(main_text)

            # Create a dictionary with the filename and the main_text and embedding
            output_dict = {"original_filename": filename,
                           "main_text": main_text,
                           "embedding": embedding}

            # Save the dictionary to a JSON file
            save_json(os.path.join(out_path, filename), output_dict)


if __name__ == "__main__":
    # Read Key
    T = pd.read_csv('../API_key.csv')
    key = T['OpenAI API key']
    key = key[0]
    openai.api_key = key

    # Define the directories where the text files and output JSON files are located
    dir_path = "download/scopus/JSON/abstracts/"
    out_path = '/Users/marioacuna/Library/Mobile Documents/com~apple~CloudDocs/AI_results/Knowledge/embeddings_json'
    # out_path = "embeddings_json/"

    # Process the json files
    process_json_files(dir_path, out_path)
