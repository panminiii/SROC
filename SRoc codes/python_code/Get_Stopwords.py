





stopwordsfilePath = 'E:/Text-analytics-with-python-master/Old_Edition_v1/Notebooks/Ch03_Processing_and_Understanding_Text/stopwords_418.txt'
with open(stopwordsfilePath, 'r') as f:  
      stopword = ''
      for line in f.readlines():  
           lineStr = line.strip()  
           if lineStr :
               stopword = stopword+lineStr+', '
           else:
               break
f1 = open('stopwords418.txt','w')
f1.write(stopword)
f1.close()
