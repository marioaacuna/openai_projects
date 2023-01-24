import openai
import pandas as pd

# Read Key
T = pd.read_csv('../API_key.csv')
key = T['OpenAI API key']
key = key[0]
openai.api_key = key






def opt3_completion(promt, engine='text-davinci-003', temp=0.7, top_p=1.0, tokens=400, freq_pen=0.0, pres_pen=0.0,
                    stop=['<<EMD>>']):
    prompt = promt.encode('ASCII', errors='ignore').decode()
    response = openai.Completion.create(
        engine=engine,
        prompt=prompt,
        temperature=temp,
        max_tokens=tokens,
        top_p=top_p,
        frequency_penalty=freq_pen,
        presence_penalty=pres_pen,
        stop=stop)
    text = response['choices'][0]['text'].strip()
    return text


if __name__ == '__main__':
    prompt = 'Write a list of famous chilean actors:'
    response = opt3_completion(prompt)
    print(response)
