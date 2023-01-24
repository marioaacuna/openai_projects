import csv

import openai
import pandas as pd

# Read Key
T = pd.read_csv('/Users/marioacuna/Documents/openai_projects/API_key.csv')
key = T['OpenAI API key']
key = key[0]
openai.api_key = key


def opt3_completion(promt, engine='text-davinci-003', temp=0.2, top_p=1.0, tokens=400, freq_pen=0.0, pres_pen=0.0,
                    stop=['<<END>>']):
    prmt = promt.encode('ASCII', errors='ignore').decode()
    r = openai.Completion.create(
        engine=engine,
        prompt=prmt,
        temperature=temp,
        max_tokens=tokens,
        top_p=top_p,
        frequency_penalty=freq_pen,
        presence_penalty=pres_pen,
        stop=stop)
    text = r['choices'][0]['text'].strip()
    return text


if __name__ == '__main__':
    # Enter metadata
    METADATA = pd.read_csv('/Users/marioacuna/Library/Mobile '
                           'Documents/com~apple~CloudDocs/Project_AcuNez/Studies_December2022'
                           '/Studies_December2022_CheckedMA.csv')

    iid = 0
    responses = []
    for itext in METADATA['summary']:
        iid =iid+ 1
        init_text = "From the following text, infer from where the data were taken (leucocytes?, PBMC?, etc), providing  " \
                    "me with the origin of the samples where the " \
                    "and what kind of patients (AD, PD, MS, etc). Nnothing more, be as concise as possible." \
                    " TEXT: \n <<BLOCK>>"

        text = itext
        #print(text)

        prompt = init_text.replace('<<BLOCK>>', text)
        #prompt = prompt + 'JAX: '
        response = opt3_completion(prompt)
        print(str(iid) + " " + response)
        responses.append(response)
        #conversation.append('JAX: %s' % response)

    with open('output.csv', 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['iteration', 'output'])
        for i, output in enumerate(responses):
            writer.writerow([i, output])