# -*- coding: utf-8 -*-

import json
import tensorflow as tf
import os,shutil

class InputExample(object):
  """A single training/test example for simple sequence classification."""

  def __init__(self, guid, text_a, text_b=None, label=None):
    """Constructs a InputExample.

    Args:
      guid: Unique id for the example.
      text_a: string. The untokenized text of the first sequence. For single
        sequence tasks, only this sequence must be specified.
      text_b: (Optional) string. The untokenized text of the second sequence.
        Only must be specified for sequence pair tasks.
      label: (Optional) string. The label of the example. This should be
        specified for train and dev examples, but not for test examples.
    """
    self.guid = guid
    self.text_a = text_a
    self.text_b = text_b
    self.label = label


class DataExample(object):

    def __init__(self): # known special case of object.__init__
        pass

    def gettestexample(self, datadir, json, max_seq_length=128):
        with open(datadir) as jsf:
            jsonFile = json.load(jsf)
            i = 0
            examples = []
            for json in jsonFile: #遍历查询数组
                article_list = json.get("article_list")
                query = json.get('query')
                if len(query) > max_seq_length - 2:
                    query = query[0:(max_seq_length - 2)]
                for article in article_list:
                    i = i + 1
                    j = 0
                    sentence_list = article.get("sentence_list")
                    for sentence in sentence_list:
                        my_sentence = sentence.get("sentence")
                        if len(my_sentence) > max_seq_length - 2:
                            my_sentence = my_sentence[0:(max_seq_length - 2)]
                        # tf.logging.info("query is %s ", query)
                        # tf.logging.info("sentence is %s ", my_sentence)
                        examples.append(InputExample(guid="%s-%s" % (i, j), text_a=my_sentence, text_b=query, label = "contradiction"))
                        j = j + 1
                    # if i % 100 == 0:
                        # tf.logging.info("get a new article id is %s ", article.get("article_id"))

        return examples

    def settestexample(self, datadir, json, destdir, scores=None):
        row = 0

        with open(datadir) as jsf:
            json_file = json.load(jsf)
            tf.logging.info("load jsonfile@@@@@@@@@@@")
            for json2 in json_file:
                # tf.logging.info("json2 in json_file@@@@@@@@@@@")
                article_list = json2.get("article_list")
                i = 0
                for article in article_list:
                    i = i + 1
                    # if i % 100 == 0:
                    #     tf.logging.info("@@output a new article id is %s ", article.get("article_id"))
                    sentence_list = article.get("sentence_list")
                    for sentence in sentence_list:
                        del sentence['sentence']
                        sentence['sentence_sco'] = scores[row]
                        row = row + 1

        with open(destdir, 'w') as jtf:
            tf.logging.info("@@@destdir is %s ", destdir)
            jtf.writelines(json.dumps(json_file, ensure_ascii=False))

        fpath, fname = os.path.split(datadir)    #分离文件名和路径
        out_dir = os.path.join(fpath, "out", fname)
        move_file(datadir, out_dir)


def move_file(srcfile, dstfile):
    if not os.path.isfile(srcfile):
        print ("%s not exist!" % (srcfile))
    else:
        fpath, fname = os.path.split(dstfile)    #分离文件名和路径
        if not os.path.exists(fpath):
            os.makedirs(fpath)                #创建路径
        shutil.move(srcfile, dstfile)          #移动文件


# examples = gettestExample("E:/DATA/nlp/IR/BERT_AP90_sentenceScoListbak.json", json)
if __name__ == "__main__":
    
    examples = DataExample().gettestexample("F:/data128/out/DataSet_DISK45-NO-CR query=0htmlContent.json", json)
    #examples = DataExample().gettestexample("E:/DATA/nlp/IR/BERT_AP90_sentenceScoList_bak.json", json)
    #inputExample.settestExample("E:/DATA/nlp/IR/BERT_AP90_sentenceScoList_bak.json", json, "E:/DATA/nlp/IR/BERT_AP90_sentenceScoListbak_output.json")

    for example in examples:
       # print(example)
        print(example.__getattribute__("text_a"))
        print(example.__getattribute__("text_b"))
