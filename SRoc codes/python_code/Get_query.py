


import os
import nltk
import re
import string
from pprint import pprint
from contractions import  CONTRACTION_MAP
from nltk.stem import PorterStemmer
import json
title_tag = "<title>"

query_dict ={}
queryCode_Num = 0
queryCode = ""






queryfilePath = 'MyRoc\disk45\Topics.disk45.301-450'
with open(queryfilePath, 'r') as f:  
      isFindtitle = False  
      isStartReadtitle = False  
      for line in f.readlines():  
           lineStr = line.strip()  
           if title_tag in lineStr:  
                isFindtitle = True
                lineStr = line.strip(title_tag)
                queryCode = lineStr.strip()
                
                queryCode_Num = queryCode_Num + 1  
                query_dict[str(queryCode_Num)] = queryCode  

                




json_str1=json.dumps(query_dict)
with open('Disk45_query_dict301-450.json','w') as json_file:
	json_file.write(json_str1)
'''
调用json的方式
with open("query_dict1-50.json",'r') as load_f:
    load_dict1 = json.load(load_f)

'''
print len(query_dict)
