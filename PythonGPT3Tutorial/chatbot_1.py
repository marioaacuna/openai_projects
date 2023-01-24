import openai
import pandas as pd

# Read Key

T = pd.read_csv('../API_key.csv')
key = T['OpenAI API key']
key = key[0]
openai.api_key = key


def opt3_completion(promt, engine='text-davinci-003', temp=0.7, top_p=1.0, tokens=400, freq_pen=0.0, pres_pen=0.0,
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
    init_text = 'The following is a conversation between USER and JAX: \n <<BLOCK>>'
    conversation = list()
    while True:
        # Create a conversation
        user_input = input('USER: ')
        if user_input == 'exit':
            print('Exiting...')
            break
        conversation.append('USER: %s' % user_input)
        text_block = '\n'.join(conversation)
        prompt = init_text.replace('<<BLOCK>>', text_block)
        prompt = prompt + 'JAX: '
        response = opt3_completion(prompt)
        print('JAX: ', response)
        conversation.append('JAX: %s' % response)
