# SRoc:Combination of Sentence-Level and Term-Level Weights

>A Probabilistic Framework for Integrating Sentence-level Semantics via BERT into Pseudo-relevance Feedback<br/>
a probabilistic framework that incorporates sentence-level semantics via Bidirectional Encoder Representations from Transformers (BERT) into PRF.
First, we obtain the importance of terms at the term level. Then, we use BERT to interactively encode the query and sentences in the feedback document to acquire the semantic similarity score of a sentence and the query.
 Next, the semantic scores of different sentences are summed as the term score at the sentence level.
Finally, we balance the term-level and sentence-level weights by adjusting factors and combine the terms with the top-k scores to form a new query for the next-round processing. 

>The two folders are processing code for neural IR part and processing code for traditional model part.



