package sigir;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
  
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map.Entry;
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

public class Retrieval_RM3 {                                                       

	public static double alpha=0.55;                                                   
	public static int N = 10;                                                         
	public static int N1 = 30;                                                        
	public static int mu = 500;                                                      

  
  
  
	public static String result="";                                  
	public static String report="";                              
	
	public static String[] topics={"topics  N.51-150","topics  
	public static String[] qrels={"qrels  els.WSJ.151-200"}; 
	public static String[] index={"index_AP90","index_AP8889","index_disk12","index_disk45","index_wt2g","index_wt10g","index_FT","index_FBIS","index_LA","index_SJMN","index_WSJ"}; 
	public static String[] dataName={"AP90_RM3","AP8889_RM3","disk12_RM3","disk45_RM3","wt2g_RM3","wt10g_RM3","FT_RM3","FBIS_RM3","LA_RM3","SJMN_RM3","WSJ_RM3"};
	
	String docNameField = "docno";	                                                  
	String docContentField = "contents";                                              
	
	int numDocs;
	int maxResults = 1000;
	long total_term_freq = 0;                                                         
	double avg_doc_length=0.0;                                                        
	
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

	  
	public static void main(String[] args) throws Exception{
		Param_Test();
	}
	
	
	  
	public static void Param_Test() throws Exception {
		for (int i = 6; i <= 10; i++){
			for (N1 = 20; N1 <= 50; N1 += 10) {
				if(N1==40) continue;
				switch(i)
				{
					case 0:mu=600;if(N1==10) {alpha=0.7;}if(N1==20) {alpha=0.7;}if(N1==30) {alpha=0.8;}if(N1==50) {alpha=0.8;}break;
					case 1:mu=750;if(N1==10) {alpha=0.5;}if(N1==20) {alpha=0.7;}if(N1==30) {alpha=0.6;}if(N1==50) {alpha=0.7;}break;
					case 2:mu=850;if(N1==10) {alpha=0.5;}if(N1==20) {alpha=0.7;}if(N1==30) {alpha=0.7;}if(N1==50) {alpha=0.7;}break;
					case 3:mu=300;if(N1==10) {alpha=0.6;}if(N1==20) {alpha=0.6;}if(N1==30) {alpha=0.7;}if(N1==50) {alpha=0.7;}break;
					case 4:mu=1400;if(N1==10) {alpha=0.3;}if(N1==20) {alpha=0.4;}if(N1==30) {alpha=0.5;}if(N1==50) {alpha=0.5;}break;
					case 5:mu=1600;if(N1==10) {alpha=0.4;}if(N1==20) {alpha=0.3;}if(N1==30) {alpha=0.4;}if(N1==50) {alpha=0.4;}break;
					case 6:mu=1450;if(N1==10) {alpha=0.8;}if(N1==20) {alpha=0.3;}if(N1==30) {alpha=0.4;}if(N1==50) {alpha=0.4;}break;
					case 7:mu=350;if(N1==10) {alpha=0.5;}if(N1==20) {alpha=0.3;}if(N1==30) {alpha=0.4;}if(N1==50) {alpha=0.4;}break;
					case 8:mu=850;if(N1==10) {alpha=0.6;}if(N1==20) {alpha=0.3;}if(N1==30) {alpha=0.4;}if(N1==50) {alpha=0.4;}break;
					case 9:mu=350;if(N1==10) {alpha=0.8;}if(N1==20) {alpha=0.3;}if(N1==30) {alpha=0.4;}if(N1==50) {alpha=0.4;}break;
					case 10:mu=1200;if(N1==10) {alpha=0.7;}if(N1==20) {alpha=0.3;}if(N1==30) {alpha=0.4;}if(N1==50) {alpha=0.4;}break;
					
				}
				for(alpha=0.5;alpha<0.9;alpha=Math.round((alpha+0.1)*10.0)  
				{  
					Param_Name(N, N1, alpha, i);                                  
					new Retrieval_RM3().RunLM_PFB(topics[i], qrels[i], index[i], result, report);
				}
			}
		}
	}


	  
	public static void Param_Name(int N, int N1, double alpha,int i) {
		StringBuffer Name = new StringBuffer("result  
		Name.append(dataName[i]).append("_mu_").append(String.format("%d", mu)).append("_N_").append(String.format("%2d", N));
		Name.append("_N1_").append(String.format("%2d", N1)).append("_alpha_").append(String.format("%.2f", alpha));
		StringBuffer Name1 = new StringBuffer(Name);
		result = String.valueOf(Name.append("_c.txt"));
		report = String.valueOf(Name1.append("_report.txt"));
		System.out.println(result);                                                  
	}
	
	  
	public void RunLM_PFB(String topicFile, String qrelsFile, String indexFile, String resultFile, String reportFile) throws Exception {
		System.out.println("RunLM_PFB...");
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
		total_term_freq=0;
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
				if(docsEnum!=null)
				while ((docsEnum.nextDoc()) != DocIdSetIterator.NO_MORE_DOCS) {
					int doc_id = docsEnum.docID();                                    
					  
					doc_set.add(doc_id);                                              
					freq_array[k][doc_id] = docsEnum.freq();                          
				}
			}
  
			numDocs = reader.numDocs();                                               
			
			HashMap<String,Double> q_Doc=new HashMap<String,Double>();
			for(Term t:term_list)              
			{
				q_Doc.put(t.text(), (1.0)  
			}
			
			  
			ScoreDoc[] compact_score_array=QueryExe(doc_set,term_list,freq_array,q_Doc);  
			  
			int min = Math.min(N, compact_score_array.length);

			
			double doc_score_Sum=0.0;                                                  
			for(int n=0;n<min;n++)
			{
				  
				doc_score_Sum+=Math.pow(2, compact_score_array[n].score);
				  
			}
			  
			  
			ArrayList<HashMap<String,Double>> DocN_list = new ArrayList<HashMap<String,Double>>(); 
			for(int n=0;n<min;n++)
			{
				HashMap<String,Double> Vecter_DocN = new HashMap<String,Double>();
				int docNum=compact_score_array[n].doc;                                
				double doc_score=Math.pow(2, compact_score_array[n].score)  
				  
				int dl = doc_length_map.get(docNum);                                  
				Terms terms = reader.getTermVector(docNum, docContentField);	      
				TermsEnum iterator = terms.iterator(null);                                  
				BytesRef byteRef = null;
				while ((byteRef = iterator.next()) != null) {  
				    String term = new String(byteRef.bytes, byteRef.offset, byteRef.length);
					double tf=iterator.totalTermFreq(); 
					Long ctf= term_total_freq_map.get(term);
					double term_score=1.0*tf  
					Vecter_DocN.put(term, term_score*doc_score);                     
					
					  
				}
				DocN_list.add(Vecter_DocN);                                           
			}
			  
			
			  
			HashMap<String,Double> sum_DocN=new HashMap<String,Double>();
			for(int m=0;m<DocN_list.size();m++)
			{
				sum_DocN=addVecter(sum_DocN, DocN_list.get(m));                       
			}			
			
			  
			HashMap<String,Double> pre_DocN=new HashMap<String,Double>();
			pre_DocN=sortVecter(sum_DocN);                                            
			
			pre_DocN=addVecter_alpha(q_Doc,pre_DocN);                                 
			  
			
			  
			  
			ScoreDoc[] compact_score_array1 = PFB_QueryExe(atomicReader , pre_DocN);

			  
			CreateResult(compact_score_array1,qualityQueries,i,query,stats);
		}
		WriteResult(stats);                                                           
	}
	
	  
	public HashMap<String,Double> addVecter_alpha(HashMap<String,Double> A, HashMap<String,Double> B) throws Exception {
		HashMap<String,Double> C=new HashMap<String,Double>();
		for(Entry<String,Double> a : A.entrySet())
		{
			if(B.containsKey(a.getKey()))
			{				
				C.put(a.getKey(), ((1-alpha)*a.getValue()+alpha*B.get(a.getKey())));
			}
			else C.put(a.getKey(), (1-alpha)*a.getValue());			
		}
		for(Entry<String,Double> b : B.entrySet())
		{
			if(!C.containsKey(b.getKey()))
			{				
				C.put(b.getKey(), alpha*b.getValue());
			}					
		}
		return C;
	}
	
	  
	public ScoreDoc[] QueryExe(Set<Integer> doc_set,ArrayList<Term> term_list,int[][] freq_array,HashMap<String,Double> q_Doc){
	  
		int kk = 0;	
		ScoreDoc[] compact_score_array=new ScoreDoc[doc_set.size()];
		for (int j = 0; j < numDocs; j++) {                                           
			if (doc_set.contains(j)) {                                                
				double total_score = 0.0;                                             
				
				int dl = doc_length_map.get(j);                                       
				for (int k = 0; k < term_list.size(); k++) {                          
					  
					  
					Long ctf= term_total_freq_map.get(term_list.get(k).text());       
					if(ctf==null) ctf=(long) 0;
					double tf = freq_array[k][j];	                                  
					  
					total_score+=log2((1.0*tf  
					                                                                  
				} 
				  
			compact_score_array[kk++] = new ScoreDoc(j, (float) total_score);         
			}
		}
		Arrays.sort(compact_score_array, new ByWeightComparator());                   
		return compact_score_array;
	}
	                                                                                                                                                                                        
	  
	  
	public ScoreDoc[] PFB_QueryExe(AtomicReader atomicReader ,HashMap<String,Double>  pre_DocN) throws IOException{
		int k=0;
		int[][] freq_array = new int[pre_DocN.size()][reader.numDocs()];
		Set<Integer> doc_set = new HashSet<Integer>();
		ArrayList<String> query_list = new ArrayList<String>(); 
		for(Entry<String,Double> pre : pre_DocN.entrySet()){
			Term term = new Term(docContentField,pre.getKey());
			query_list.add(pre.getKey());
			  
			DocsEnum docsEnum = atomicReader.termDocsEnum(term);	                  
			if(docsEnum!=null)
			while ((docsEnum.nextDoc()) != DocIdSetIterator.NO_MORE_DOCS) {
				int doc_id = docsEnum.docID();                                        
				doc_set.add(doc_id);                                                  
				freq_array[k][doc_id] = docsEnum.freq();                              
			}
			k++;
		}
		
		int kk = 0;
		ScoreDoc[] compact_score_array=new ScoreDoc[doc_set.size()];
		for (int j = 0; j < numDocs; j++) {                                           
			if (doc_set.contains(j)) {                                                
				double total_score = 0.0;                                             
				int dl = doc_length_map.get(j);
				for (int i = 0; i < query_list.size(); i++) {                         
					double qtf=pre_DocN.get(query_list.get(i));
					Long ctf= term_total_freq_map.get(query_list.get(i));
					if(ctf==null) ctf=(long) 0;
					double tf = freq_array[i][j];                                     
					  
					total_score +=log2(1.0*tf  
					  
				}  
				  
			compact_score_array[kk++] = new ScoreDoc(j, (float) (total_score));       
			}
		}
		Arrays.sort(compact_score_array, new ByWeightComparator());                   
		return compact_score_array;
	}
	
	  
	public HashMap<String, Double> sortVecter(HashMap<String,Double> A)  {
		List<HashMap.Entry<String,Double>> list = new ArrayList<HashMap.Entry<String,Double>>(A.entrySet());
		Collections.sort(list,new Comparator<HashMap.Entry<String,Double>>() {
              
            public int compare(Entry<String, Double> o1,
                    Entry<String, Double> o2) {
                 return (o1.getValue().compareTo(o2.getValue()))*-1;                  
            }
		});
		
		HashMap<String,Double> C=new HashMap<String,Double>();
		double sum=0.0;
		for(int i=0;i<N1;i++)
		{
			sum+= list.get(i).getValue();	                                          
		}
		for(int i=0;i<N1;i++)
		{
			if(sum==0)  
			{
				C.put(list.get(i).getKey(), 0.0);	            	
			}
			else 
			{
				C.put(list.get(i).getKey(), list.get(i).getValue()  
				
				  
			}
			
			
			  
			
		}
		
		return C;
	}
	
	
	
	  
	public HashMap<String,Double> addVecter(HashMap<String,Double> A, HashMap<String,Double> B) throws Exception {
		HashMap<String,Double> C=new HashMap<String,Double>();
		for(Entry<String,Double> a : A.entrySet())
		{
			if(B.containsKey(a.getKey()))
			{				
				C.put(a.getKey(), (a.getValue()+B.get(a.getKey())));
			}
			else C.put(a.getKey(), a.getValue());			
		}
		for(Entry<String,Double> b : B.entrySet())
		{
			if(!C.containsKey(b.getKey()))
			{				
				C.put(b.getKey(), b.getValue());
			}					
		}
		return C;
	}
	
	  
	public void CreateResult(ScoreDoc[] compact_score_array1,QualityQuery[] qualityQueries,int i,Query query,QualityStats[] stats) throws IOException{
		int max_result = Math.min(maxResults, compact_score_array1.length);
		ScoreDoc[] score_doc = new ScoreDoc[max_result];
		System.arraycopy(compact_score_array1, 0, score_doc, 0, max_result);	      
		TopDocs td = new TopDocs(max_result, score_doc, (float) score_doc[0].score);  
		stats[i] = analyzeQueryResults(qualityQueries[i], query, td, judge, qualityLog, 1);  
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
	
	  
	public void WriteResult(QualityStats[] stats){
		QualityStats avg = QualityStats.average(stats);		                  
		avg.log("SUMMARY", 2, qualityLog, "  ");
	}
	
	  
	public void printVecter(HashMap<String,Double> A){
		for(Entry<String,Double> a : A.entrySet())
		{
			System.err.println("name=:"+a.getKey()+" value=:"+a.getValue());
		}
	}
	
	  
	public static double log2(double n) {
		return (Math.log(n)   
	}
	
}


