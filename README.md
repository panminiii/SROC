# SRoc:Combination of Sentence-Level and Term-Level Weights

>A Probabilistic Framework for Integrating Sentence-level Semantics via BERT into Pseudo-relevance Feedback<br/>
a probabilistic framework that incorporates sentence-level semantics via Bidirectional Encoder Representations from Transformers (BERT) into PRF.
First, we obtain the importance of terms at the term level. Then, we use BERT to interactively encode the query and sentences in the feedback document to acquire the semantic similarity score of a sentence and the query.
 Next, the semantic scores of different sentences are summed as the term score at the sentence level.
Finally, we balance the term-level and sentence-level weights by adjusting factors and combine the terms with the top-k scores to form a new query for the next-round processing. 

//The two folders are processing code for neural IR part and processing code for traditional model part.

>基于BERT将句子级语义集成到伪关联反馈的概率框架<br/>
一个概率框架，通过从转换器的双向编码器表示(BERT)将句子级语义合并到PRF中。
首先，我们得到项在项级的重要性。然后，我们使用BERT对反馈文档中的查询和句子进行交互编码，获得句子和查询的语义相似度评分。
然后，将不同句子的语义得分归纳为句子级的术语得分。
最后，我们通过调整因子来平衡词汇级和句子级的权重，并将词汇与排名前k的分数结合起来，形成一个新的查询，以进行下一轮的处理。

//两个文件夹分别为神经IR部分的处理代码和传统模型部分的处理代码。


