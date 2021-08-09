package sigir;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
//import java.text.Format;
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

public class Retrieval_DLM {                                                     

	public static double alpha=0.55;                                                
	public static int N = 10;                                                       
	public static int N1 = 30;                                                      
	public static int mu = 2000;                                                    

//	public static String[] topics={"topics/topics.AP90.51-100","topics/topics.AP8889.51-100","topics/topics.disk12.51-200","topics/topics.disk45.301-450","topics/topics.wt2g","topics/topics.wt10g","topics/topics.FT.301-400","topics/topics.FBIS.351-450","topics/topics.LA.301-400","topics/topics.SJMN.51-150","topics/topics.WSJ.151-200"}; 
//	public static String[] qrels={"qrels/qrels.AP90.51-100","qrels/qrels.AP8889.51-100","qrels/qrels.disk12.51-200","qrels/qrels.disk45.301-450","qrels/qrels.wt2g","qrels/qrels.wt10g","qrels/qrels.FT.301-400","qrels/qrels.FBIS.351-450","qrels/qrels.LA.301-400","qrels/qrels.SJMN.51-150","qrels/qrels.WSJ.151-200"}; 
//	public static String[] index={"index_AP90","index_AP8889","index_disk12","index_disk45","index_WT2G","index_WT10G","index_FT","index_FBIS","index_LA","index_SJMN","index_WSJ"}; 
//	public static String[] dataName={"AP90_DLM","AP8889_DLM","disk12_DLM","disk45_DLM","WT2G_DLM","WT10G_DLM","FT_DLM","FBIS_DLM","LA_DLM","SJMN_DLM","WSJ_DLM"};
	
	public static String[] topics={"topics/topics.AP89.1-50","topics/topics.AP8889.101-150","topics/topics.AP8889.151-200","topics/topics.disk45.301-450","topics/topics.wt2g","topics/topics.wt10g","topics/topics.FT.301-400","topics/topics.FBIS.351-450","topics/topics.LA.301-400","topics/topics.SJMN.51-150","topics/topics.WSJ.151-200"}; 
	public static String[] qrels={"qrels/qrels.AP89.1-50","qrels/qrels.AP8889.101-150","qrels/qrels.AP8889.151-200","qrels/qrels.disk45.301-450","qrels/qrels.wt2g","qrels/qrels.wt10g","qrels/qrels.FT.301-400","qrels/qrels.FBIS.351-450","qrels/qrels.LA.301-400","qrels/qrels.SJMN.51-150","qrels/qrels.WSJ.151-200"}; 
	public static String[] index={"index_AP89","index_AP8889","index_AP8889","index_disk45","index_WT2G","index_WT10G","index_FT","index_FBIS","index_LA","index_SJMN","index_WSJ"}; 
	public static String[] dataName={"1AP89_DLM","2AP8889_DLM","3AP8889_DLMM","disk45_DLM","WT2G_DLM","WT10G_DLM","FT_DLM","FBIS_DLM","LA_DLM","SJMN_DLM","WSJ_DLM"};
	
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
		for (int i = 4; i <= 4; i++){
			for (N1 = 10; N1 <= 10; N1 += 10) {
				if(N1==40) continue;

				for(alpha=0.0;alpha<=0.0;alpha=Math.round((alpha+0.1)*10.0)/10.0)
				{
					Param_Name(N, N1, alpha, i);                                
					new Retrieval_DLM().RunLM_PFB(topics[i], qrels[i], index[i], result, report);
				}
			}
		}
	}


	
	public static void Param_Name(int N, int N1, double alpha,int i) {
		StringBuffer Name = new StringBuffer("result/DLM/s/");
		Name.append(dataName[i]).append("_mu_").append(String.format("%d", mu)).append("_N_").append(String.format("%2d", N));
		Name.append("_N1_").append(String.format("%2d", N1)).append("_alpha_").append(String.format("%.2f", alpha));
		StringBuffer Name1 = new StringBuffer(Name);
		result = String.valueOf(Name.append("_p1.txt"));
		report = String.valueOf(Name1.append("_report.txt"));
	}
	
	
	public void RunLM_PFB(String topicFile, String qrelsFile, String indexFile, String resultFile, String reportFile) throws Exception {
		System.out.println("RunLM_PFB...");
		Directory directory = FSDirectory.open(new File(indexFile));                
		reader = DirectoryReader.open(directory);                                   
		searcher = new IndexSearcher(reader);                                       
		docNameField = "docno";                                                     
		docContentField = "contents";                                              contents=title+text
		qualityLog = new PrintWriter(new File(resultFile), "UTF-8");                
		TrecTopicsReader qReader = new TrecTopicsReader();	                        
		qualityQueries = qReader.readQueries(new BufferedReader(new FileReader(new File(topicFile))));	
	                                                                               
		judge = new TrecJudge(new BufferedReader(new FileReader(new File(qrelsFile))));	
														                           
		judge.validateData(qualityQueries, qualityLog);                             
		                                                                           
		qqParser = new MyQQParser("title", "contents");                             
                                                             
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
			double avg_tf = (term_num == 0) ? 0 : ((double) docLen) / term_num;     
			doc_avg_tf_map.put(j, avg_tf);									       
		}
		avg_doc_length = ((double) total_dl) / reader.numDocs();                    
	}

	
	public void search() throws Exception {
		AtomicReader atomicReader = SlowCompositeReaderWrapper.wrap(reader);        
		QualityStats stats[] = new QualityStats[qualityQueries.length];             
		for (int i = 0; i < qualityQueries.length; i++) {                           
			QualityQuery qq = qualityQueries[i];                                    
			Query query = qqParser.parse(qq);	                                    
			HashSet<Term> term_set = new HashSet<Term>();                            
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
			HashMap<String,Double> q_Doc=new HashMap<String,Double>();
			for(Term t:term_list)              
			{
				q_Doc.put(t.text(), (1.0)/term_list.size());                          
			}
			
			
			ScoreDoc[] compact_score_array=QueryExe(doc_set,term_list,freq_array,q_Doc);  
		
			
			CreateResult(compact_score_array,qualityQueries,i,query,stats);
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
	//public ScoreDoc[] QueryExe(Set<Integer> doc_set,ArrayList<Term> term_list,int[][] freq_array,HashMap<String,Double> q_Doc){
		int kk = 0;	
		ScoreDoc[] compact_score_array=new ScoreDoc[doc_set.size()];
		for (int j = 0; j < numDocs; j++) {                                         
			if (doc_set.contains(j)) {                                              
				double total_score = 0.0;                                           
				
				int dl = doc_length_map.get(j);                                     
				for (int k = 0; k < term_list.size(); k++) {                        
					//double qtf=q_Doc.get(term_list.get(k).text());
					//System.out.println("qtf="+qtf);					
					Long ctf= term_total_freq_map.get(term_list.get(k).text());     
					double tf = freq_array[k][j];	                                		
					//System.out.println(tf+"  "+dl+"  "+ctf+"  "+total_term_freq);
					total_score+=log2((1.0*tf/(dl+mu)+1.0*mu/(dl+mu)*ctf/total_term_freq)); 
					                                                                
				} 
				//System.err.println(String.format("%.20f",total_score));
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
			//processedQueriesWriter.write(termText+" ");                           
			DocsEnum docsEnum = atomicReader.termDocsEnum(term);	                
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
					double tf = freq_array[i][j];  				
					
					total_score +=log2(1.0*tf/(dl+mu)+1.0*mu/(dl+mu)*ctf/total_term_freq)*qtf;     
					
				}  
				
			compact_score_array[kk++] = new ScoreDoc(j, (float) (total_score));     
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
				C.put(list.get(i).getKey(), list.get(i).getValue()/sum);	            
				
				
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
		return (Math.log(n) / Math.log(2));
	}
	
}


