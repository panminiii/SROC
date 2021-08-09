package sigir;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.apache.lucene.benchmark.quality.Judge;
import org.apache.lucene.benchmark.quality.QualityQuery;
import org.apache.lucene.benchmark.quality.QualityQueryParser;
import org.apache.lucene.benchmark.quality.trec.TrecJudge;
import org.apache.lucene.benchmark.quality.trec.TrecTopicsReader;
import org.apache.lucene.benchmark.quality.utils.DocNameExtractor;
import org.apache.lucene.benchmark.quality.utils.SubmissionReport;
import org.apache.lucene.index.AtomicReader;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.DocsEnum;
import org.apache.lucene.index.Fields;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.MultiFields;
import org.apache.lucene.index.SlowCompositeReaderWrapper;
import org.apache.lucene.index.Term;
import org.apache.lucene.index.Terms;
import org.apache.lucene.index.TermsEnum;
import org.apache.lucene.search.DocIdSetIterator;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.ScoreDoc;
import org.apache.lucene.search.TopDocs;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.FSDirectory;
import org.apache.lucene.util.BytesRef;

public class Retrieval_BM25_test {                                                         
	 
	public static double k1 = 1.2;                                                   
	public static double b = 0.35;                                                    

	public static String[] topics={"topics  N.51-150","topics  
	public static String[] qrels={"qrels  els.WSJ.151-200"}; 
	public static String[] index={"index_AP90","index_AP8889","index_disk12","index_disk45","index_wt2g","index_wt10g","index_FT","index_FBIS","index_LA","index_SJMN","index_WSJ"}; 
	public static String[] dataName={"AP90_BM25","AP8889_BM25","disk12_BM25","disk45_BM25","wt2g_BM25","wt10g_BM25","FT_BM25","FBIS_BM25","LA_BM25","SJMN_BM25","WSJ_BM25"};
	
	String docNameField = "docno";	                                                  
	String docContentField = "contents";                                              
	
	int numDocs;
	int maxResults = 1000;
	long total_term_freq = 0;                                                         
	double avg_doc_length=0.0;                                                        
	static int name_i=0;
	QualityQuery[] qualityQueries;	
	QualityQueryParser qqParser;	                                                  
	IndexReader reader;	
	IndexSearcher searcher;			
	PrintWriter qualityLog;	
	SubmissionReport submitRep;	
	Judge judge;	
	
	
	
	HashMap<String, Long> term_total_freq_map = new HashMap<String, Long>();		  
	HashMap<String, Integer> term_doc_freq_map = new HashMap<String, Integer>();	  
	HashMap<Integer, Integer> doc_length_map = new HashMap<Integer, Integer>();		  
	HashMap<Integer, Double> doc_avg_tf_map = new HashMap<Integer, Double>();         

	  
	public static double log2(double n) {
		return (Math.log(n)   
	}

	
	  
	public static void main(String[] args) throws Exception{
		  
		
		
		Param_Test();
	}
	
	
	  
	public static void Param_Test() throws Exception 
	{
		for (int i = 0; i <= 0; i++) {
			switch(i)
			{                       
				case 0:b=0.55;break;
				case 1:b=0.0 ;break;
				case 2:b=0.35;break;
				case 3:b=0.35;break;
				case 4:b=0.25;break;
				case 5:b=0.2;break;
				case 6:b=0.3;break;
				case 7:b=0.05;break;
				case 8:b=0.3;break;
				case 9:b=0.55;break;
				case 10:b=0.3;break;
			}   
			    name_i=i;
				String resultFile="result  
				String reportFile="result  
				System.out.println(resultFile);
				new Retrieval_BM25_test().RunBM25_test(topics[i], qrels[i], index[i], resultFile, reportFile);
		}
	}

	
	  
	public void RunBM25_test(String topicFile, String qrelsFile, String indexFile, String resultFile, String reportFile) throws Exception {
		  
		Directory directory = FSDirectory.open(new File(indexFile));                  
		reader = DirectoryReader.open(directory);                                     
		searcher = new IndexSearcher(reader);                                         
		docNameField = "docno";                                                       
		docContentField = "contents";                                                 
		qualityLog = new PrintWriter(new File(resultFile), "UTF-8");                  
		TrecTopicsReader qReader = new TrecTopicsReader();	                          
		qualityQueries = qReader.readQueries(new BufferedReader(new FileReader(new File(topicFile))));	  
	                                                                                  
		judge = new TrecJudge(new BufferedReader(new FileReader(new File(qrelsFile))));	  
														                              
		judge.validateData(qualityQueries, qualityLog);                               
		                                                                              
		qqParser = new MyQQParser("title", "contents");                               
  
  
		execute();                                                                    
		directory.close();                                                            
		qualityLog.close();                                                           
  
	}	
	
	  
	public void execute() throws Exception {
		termStats();
		docStats();
		search();
	}	
	
	  
	public void termStats() throws Exception {	                                      
		Fields fields = MultiFields.getFields(reader);                                
		Terms terms = fields.terms(docContentField);                                  
		TermsEnum iterator = terms.iterator(null);                                    
		  
		BytesRef byteRef = null;
		while ((byteRef = iterator.next()) != null) {                                 
			String text = new String(byteRef.bytes, byteRef.offset, byteRef.length);  
			term_total_freq_map.put(text, iterator.totalTermFreq());                  
			term_doc_freq_map.put(text, iterator.docFreq());	                      
			
			
			  
			
			
			total_term_freq += iterator.totalTermFreq();                              
		}
	}	

	  
	public void docStats() throws Exception {	
		long total_dl = 0;
		for (int j = 0; j < reader.numDocs(); j++) {                                  
			int docLen = 0;		
			int term_num = 0;	
			Terms terms = reader.getTermVector(j, docContentField);	                  
			if (terms != null && terms.size() > 0) {                                  
				TermsEnum termsEnum = terms.iterator(null);		                      
				while ((termsEnum.next()) != null) {
					int freq = (int) termsEnum.totalTermFreq();                       
					docLen += freq;                                                   
					term_num++;                                                       
				}
			}
			total_dl += docLen;
			doc_length_map.put(j, docLen);                                            
			double avg_tf = (term_num == 0) ? 0 : ((double) docLen)   
			doc_avg_tf_map.put(j, avg_tf);									          
		}
		avg_doc_length = ((double) total_dl)   
	}

	  
	public void search() throws Exception {
  
		AtomicReader atomicReader = SlowCompositeReaderWrapper.wrap(reader);                           
		QualityStats stats[] = new QualityStats[qualityQueries.length];                                
		for (int i = 0; i < qualityQueries.length; i++) {                             
			QualityQuery qq = qualityQueries[i];                                      
			Query query = qqParser.parse(qq);	                                      
			HashSet<Term> term_set = new HashSet<Term>();                             
			query.extractTerms(term_set);	                                          
			ArrayList<Term> term_list = new ArrayList<Term>();                      
			for (Term term : term_set) {                                              
				DocsEnum docsEnum = atomicReader.termDocsEnum(term);	              
				if(docsEnum!=null)
				term_list.add(term);                                                  
			}
			Set<Integer> doc_set = new HashSet<Integer>();	
			int[][] freq_array = new int[term_list.size()][reader.numDocs()];         

			for (int k = 0; k < term_list.size(); k++) {
				Term term = term_list.get(k);
				  
  
				DocsEnum docsEnum = atomicReader.termDocsEnum(term);	              
				while ((docsEnum.nextDoc()) != DocIdSetIterator.NO_MORE_DOCS) {
					int doc_id = docsEnum.docID();                                    
					  
					doc_set.add(doc_id);                                              
					freq_array[k][doc_id] = docsEnum.freq();                          
				}
			}
  
			numDocs = reader.numDocs();                                               
			  
			  
			ScoreDoc[] compact_score_array=QueryExe(doc_set,term_list,freq_array);  
			
			  
			
			
			
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
              
  
  
  
  
            
  
  
  
  
  
  
  
  
  
  
  
  
  
            
  
  
            
            

			  
			CreateResult(compact_score_array,qualityQueries,i,query,stats);
		}
		  
		WriteResult(stats);           
	}
	
	  
	public ScoreDoc[] QueryExe(Set<Integer> doc_set,ArrayList<Term> term_list,int[][] freq_array){
		int kk = 0;	
		ScoreDoc[] compact_score_array=new ScoreDoc[doc_set.size()];
		for (int j = 0; j < numDocs; j++) {                                       
			if (doc_set.contains(j)) {                                            
				double total_score = 0.0f;                                        
				int doc_length = doc_length_map.get(j);                           
				double K = k1 * ((1 - b) + b * doc_length   
				for (int k = 0; k < term_list.size(); k++) {                      
					int df = term_doc_freq_map.get(term_list.get(k).text());      
					double tf = freq_array[k][j];	                              
					double TF = (k1 + 1) * tf   
					double IDF = log2((numDocs - df + 0.5)   
					total_score += TF * IDF;
					  
				}   
			compact_score_array[kk++] = new ScoreDoc(j, (float) total_score);     
			}
		}
		Arrays.sort(compact_score_array, new ByWeightComparator());               
		return compact_score_array;
	}
	
	  
	public void CreateResult(ScoreDoc[] compact_score_array1,QualityQuery[] qualityQueries,int i,Query query,QualityStats[] stats) throws IOException{
		int max_result = Math.min(maxResults, compact_score_array1.length);
		ScoreDoc[] score_doc = new ScoreDoc[max_result];
		System.arraycopy(compact_score_array1, 0, score_doc, 0, max_result);	      
		TopDocs td = new TopDocs(max_result, score_doc, (float) score_doc[0].score);  
		stats[i] = analyzeQueryResults(qualityQueries[i], query, td, judge, qualityLog, 1);  
	}
	
	  
	public void WriteResult(QualityStats[] stats){
		QualityStats avg = QualityStats.average(stats);		                  
		avg.log("SUMMARY", 2, qualityLog, "  ");
	}

	  
	private QualityStats analyzeQueryResults(QualityQuery qq, Query q, TopDocs td, Judge judge, PrintWriter logger, long searchTime) throws IOException {
		QualityStats stts = new QualityStats(judge.maxRecall(qq), searchTime);
		long t1 = System.currentTimeMillis();	
		ScoreDoc[] scoreDocs = td.scoreDocs;
		DocNameExtractor xt = new DocNameExtractor(docNameField);
		for (int i = 0; i < scoreDocs.length; i++) {
			String docName = xt.docName(searcher, scoreDocs[i].doc);	
			long docNameExtractTime = System.currentTimeMillis() - t1;
			t1 = System.currentTimeMillis();
			boolean isRelevant = judge.isRelevant(docName, qq);		
			stts.addResult(i + 1, isRelevant, docNameExtractTime);
		}
		if (logger != null) {
			logger.println(qq.getQueryID() + "  -  " + q);
			stts.log(qq.getQueryID() + " Stats:", 1, logger, "  ");
		}
		return stts;
	}
}
