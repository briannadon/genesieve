import re
from sys import argv
#import smart_open

def remove_text_inside_brackets(text, brackets="()[]{}"):
    count = [0] * (len(brackets) // 2) # count open/close brackets
    saved_chars = []
    for character in text:
        for i, b in enumerate(brackets):
            if character == b: # found bracket
                kind, is_close = divmod(i, 2)
                count[kind] += (-1)**is_close # `+1`: open, `-1`: close
                if count[kind] < 0: # unbalanced bracket
                    count[kind] = 0  # keep it
                else:  # found bracket to remove
                    break
        else: # character is not a [balanced] bracket
            if not any(count): # outside brackets
                saved_chars.append(character)
    return ''.join(saved_chars)

def sanitize_text(string):
    regex = re.compile('[^a-zA-Z\.\,\ ]')
    string = remove_text_inside_brackets(string)
    string = regex.sub('',string)
    return string
