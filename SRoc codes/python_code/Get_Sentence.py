

import os
import nltk
import re
import json
from contractions import  CONTRACTION_MAP


N=6         
queryl=150   


query_content=[]
with open("MyRoc\disk45\Disk45_NO_CRhtmlContent.json",'r') as load_f:   
    query_content = json.load(load_f)

docnum_list=[]
with open('MyRoc\disk45\Disk45_Roc_docnum_first_retrievalN=6_b_0.35.txt', 'r') as f: 
	for line in f.readlines():  
		lineStr = line.strip()  
		docnum_list.append(lineStr)
docnum_list1=[docnum_list[x:x+6] for x in range(0,len(docnum_list),6)] 


with open("MyRoc\disk45\Disk45_query_dict301-450.json",'r') as load_f:        
	query_dict = json.load(load_f)







def expend_contractions(sentence,contraction_mapping):
	contractions_pattern=re.compile('({})'.format('|'.join(contraction_mapping.keys())),flags=re.IGNORECASE|re.DOTALL)
	def expand_match(contraction):
		match=contraction.group(0)
		first_char=match[0]
		expanded_contraction=contraction_mapping.get(match)\
		                      if contraction_mapping.get(match)\
			                  else contraction_mapping.get(match.lower())
		expanded_contraction=first_char+expanded_contraction[1:]
		return expanded_contraction
	expend_sentence= contractions_pattern.sub(expand_match,sentence)
	return expend_sentence


def remove_characters_before_tokenization(sentence,
										  keep_apostrophes=False):
	sentence = sentence.strip()
	if keep_apostrophes:                                       
		PATTERN = r"[-|{|}|\\|_|$|&|*|%|@|(|)|~|\n|\f]"
		filtered_sentence = re.sub(PATTERN, r'', sentence)
	else:
		PATTERN = r'[^a-zA-Z0-9 ]'
		filtered_sentence = re.sub(PATTERN, r'', sentence)
	return filtered_sentence


def right_shangyinhao(sentence):
	shangyinhao = '`'
	shuangyinhao = '``'
	left='['
	right=']'
	if shuangyinhao in sentence:
		sentence = sentence.replace(shuangyinhao, "\"")  
	if shangyinhao in sentence:
		sentence = sentence.replace(shangyinhao, "\'")
	if 	left in sentence:
		sentence =sentence.replace(left,'')
	if 	right in sentence:
		sentence =sentence.replace(right,'')
	return sentence


List=[]
for i in range(queryl):
	qd_dict = {}
	print i
	for each_qd_dict in query_content:
		if each_qd_dict["query"]==query_dict[str(i+1)]:
			qd_dict["query"] = each_qd_dict["query"]
			articlelist = []
			for n in range(N):
				for doc_dict in each_qd_dict['article_list']:
					if doc_dict["articleId"].encode('utf-8')==str(docnum_list1[i][n]):
						doc_dict1 = {}
						doc_dict1["article_id"]=str(docnum_list1[i][n])
						
						content1=''
						content1=doc_dict["content"]
						sentence_list=[]
						if doc_dict["content"].strip()=='':  
							print '空文档编号'
							print doc_dict["articleId"]
						content1=re.sub(' +',' ',content1)

						text_sentence = nltk.sent_tokenize(content1)
						expend_text_sentence = [expend_contractions(sentence, CONTRACTION_MAP) for sentence in
												text_sentence]  
						filtered_sentence = [remove_characters_before_tokenization(sentence, keep_apostrophes=True) for sentence in
											 expend_text_sentence]  
						for j in range(len(filtered_sentence)):  
							if len(filtered_sentence[j]) > 1 and len(filtered_sentence[j]) < 400:
								filtered_sentence[j]=right_shangyinhao(filtered_sentence[j])

								sentence_list.append(filtered_sentence[j])
							if len(filtered_sentence[j]) > 400:    
								filtered_sentence[j] = right_shangyinhao(filtered_sentence[j])

								termList = filtered_sentence[j].split()
								for term in termList:
									if len(term) > 45:
										termList.remove(term)

								sen = ''
								for leg in range(len(termList)):
									s = leg // 20  
									sen = sen + termList[s * 20 + leg % 20] + ' '  
									if (leg+1) % 20==0 or leg==len(termList)-1:
										sentence_list.append(sen)
										
										sen = ''
						doc_dict1["sentence_list"]=sentence_list
						articlelist.append(doc_dict1)
	qd_dict['article_list']=articlelist
	List.append(qd_dict)
	






json_str=json.dumps(List)
with open('Disk_query_senten_listN=6.json','w') as json_file:
	json_file.write(json_str)
















