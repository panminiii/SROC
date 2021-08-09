

import os
doc_start_tag = '<DOC>'
doc_end_tag = '</DOC>'
doc_tag = "<DOCNO>"
head_start_tag = "<HEAD>"
head_end_tag = "</HEAD>"
start_text_tag = "<TEXT>"
end_text_tag = "</TEXT>"

num_docCode_dict = {}    
docCode_text_dict = {}   
docCode_Num = 0
docCode = ""

def traverse_Directory(directoryPath):
	for root, dirs, files in os.walk(directoryPath):                      
		global docCode_Num											      
		global docCode                                                    
		for file in files:											      
			filePath = os.path.join(root, file)						      
			with open(filePath, 'r') as f:							      
				isFindDoc = False
				isStartReadText = False        						      
				textContent = ""                                          
				for line in f.readlines():                                
					lineStr = line.strip()                               
					if (lineStr==doc_end_tag):							  
						isStartReadText = False
						docCode_text_dict[docCode] = textContent          
						textContent = ""							      

					if isFindDoc and (lineStr==end_text_tag):
						isStartReadText = False	

					if isStartReadText:
						if head_end_tag in line:
							textContent = textContent + lineStr.replace(head_end_tag, '.') + '\n'
							isStartReadText = False
						else:	
							textContent = textContent + line              
					
					if doc_start_tag in lineStr:
						isFindDoc = True

					if doc_end_tag in lineStr:
						isFindDoc = False	

					if isFindDoc:
						if doc_tag in lineStr:                                
							dataArray = lineStr.split()                       
							docCode = dataArray[1]						      
							docCode_Num = docCode_Num + 1                     
							num_docCode_dict[str(docCode_Num)] = docCode      
					
						if head_start_tag in lineStr:						  
							head_text = lineStr.replace(head_start_tag,'')    
							if head_end_tag in lineStr:						  
								head_text = head_text.replace(head_end_tag, '.')   
								textContent = textContent + head_text + '\n'
							else:											   
								textContent = textContent + head_text + '\n'
								isStartReadText = True
						if lineStr==start_text_tag:                           
							isStartReadText = True
					


def searchTextWithDocCode(doc_code):
	
	return docCode_text_dict.get(doc_code)

traverse_Directory("E:\Text-analytics-with-python-master\Old_Edition_v1\Notebooks\Ch03_Processing_and_Understanding_Text\DataSet_AP90")

'''

for key in num_docCode_dict:
	tags = [doc_tag,head_start_tag, head_end_tag, start_text_tag, end_text_tag]
	doc_code = num_docCode_dict[key]
	content = searchTextWithDocCode(doc_code)
	if (len(content) == 0):  
		print("异常！读取到的内容为空，文档编号:" + doc_code)
	for tag in tags:
		if tag in content:
			print("异常!读取到标签了，文档编号:" + doc_code)	
	
'''


text = searchTextWithDocCode("AP900105-0226")
print("-----------------------  查询结果  -----------------------")
print(text)
