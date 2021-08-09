

import nltk
import re
import json
import operator
from nltk.stem import PorterStemmer
queryLength=150


query_dict={}
with open("MyRoc\disk45\Disk45_query_dict301-450.json",'r') as load_f:        
	query_dict = json.load(load_f)

with open("MyRoc\disk45\Disk45_bert_senten_scoN=6.json",'r') as load_f:        
    sente_sco_list = json.load(load_f)

stopwords = [line.strip() for line in open('stopwords\stopwords_418.txt', 'r').readlines()]

def remove_characters_before_tokenization(sentence,
										  keep_apostrophes=False):
	sentence = sentence.strip()
	if keep_apostrophes:                                       
		PATTERN = r"[-|{|}|\\|_|$|&|*|%|@|(|)|~|\n]"
		filtered_sentence = re.sub(PATTERN, r'', sentence)
	else:
		PATTERN = r'[^a-zA-Z0-9 ]'
		filtered_sentence = re.sub(PATTERN, r'', sentence)
	return filtered_sentence


def Stemmer(tokens):
	ps = PorterStemmer()        
	
	filtered_tokens=[ps.stem(token)  for token in tokens]
	return filtered_tokens




termSco_list=[]

for Dict in sente_sco_list:
    queryTerm_dict={}
    query=Dict.get("query")
    queryID=list(query_dict.keys())[list(query_dict.values()).index(query)]
    queryTerm_dict["queryID"] = queryID
    termSco_dict={}
    Term_list = []
    articleList=[]
    articleList=Dict.get("article_list")
    for dict in articleList:
        sentenceList=[]
        sentenceList=dict.get("sentence_list")
        for ddict in sentenceList:
            sentence=''
            L1=[]
            sentence=ddict.get("sentence")
            filtered_sentence = remove_characters_before_tokenization(sentence)
            filtered_sentence = filtered_sentence.lower()      
            L1 = nltk.word_tokenize(filtered_sentence)           
            L2=[]
            for word in L1:
                if word not in stopwords:
                    L2.append(word)
            L3 = Stemmer(L2)
            L4 = [token.encode('utf-8') for token in L3]
            for term in L4:
                if termSco_dict.has_key(term):
                    termSco_dict[term]=float(termSco_dict[term])+float(ddict.get("score"))
                else:
                    termSco_dict[term]=float(ddict.get("score"))

    sorted_Sco_list=sorted(termSco_dict.items(),key=operator.itemgetter(1),reverse=True)
    
    queryTerm_dict["termSco_dict"] = sorted_Sco_list
    termSco_list.append(queryTerm_dict)



'''
print termSco_list
json_str=json.dumps(termSco_list)
with open('termSco_list.json','w') as json_file:
	json_file.write(json_str)
'''
NewTermScoList=[]
for i in range(queryLength+1):
    for DDict in termSco_list:
        if DDict["queryID"].encode('utf-8')==str(i+1):
            NewTermScoList.append(DDict["termSco_dict"])


json_str=json.dumps(NewTermScoList)
with open('Disk45_bert_TermScoListN=6.json','w') as json_file:
	json_file.write(json_str)
