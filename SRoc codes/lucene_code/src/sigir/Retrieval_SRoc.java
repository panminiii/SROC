package sigir;
import org.apache.commons.io.FileUtils; 
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
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
import org.json.JSONArray;
public class SRoc {
	  
	public static double k1 = 1.2;                                                    
	public static double b = 0.55;                                                    
	public static double k3 = 8.0;
	public static int N = 10;          
	  
	  
	public static int N2 = 100;
	public static int N1 = 30;                                                        
	public static double beta =0.2;                                                   
	public static double alpha=0.5;                                                   
	
	public static String[] topics ={"topics  
    public static String[] qrels={"qrels  
    public static String[] index={"index_AP90","index_AP8889","index_WT2G","index_disk45"};
    public static String[] dataName={"AP90","AP8889","WT2G","DISK45"};
     
     
  
    public static String result="";
	public static String report="";
	public static String jsonName="";
	
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
	
	  
	public static void Param_Test() throws Exception
	{
		for (int i=0;i<=0;i++) {
			  
			for (N1=10;N1<=10;N1+=10) {
				if(N1==40) continue;
				  
				{                       
					case 0:b=0.55;alpha=0.5;break;
					case 1:b=0.0 ;alpha=0.3;break;
					case 2:b=0.25;alpha=0.1;break;
					case 3:b=0.35;alpha=0.3;break;
					case 4:b=0.2;alpha=0.4;break;
					case 5:b=0.4;alpha=0.25;break;
				
				} *  
			for(N=10;N<=10;N=N+1){
				for(alpha=0.1;alpha<=0.5;alpha=round(alpha+0.1)){
				   for(beta=0.1;beta<=0.1;beta=round(beta+0.1)){
					   jsonName="json  
					   String result="result1  
						String report="report1  
						System.out.println(result+"  "+jsonName);    
						new SRoc().RunBM25_PFB(topics[i],qrels[i],index[i], result, report);
						}
				}
					}
				}
			}
		}
	}

  
  
  
  
  
  
  
  
  
  
  
	  
	public void RunBM25_PFB(String topicFile, String qrelsFile, String indexFile, String resultFile, String reportFile) throws Exception {
		  
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
  
		System.out.println("ok!");
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

		for (int i = 0; i <qualityQueries.length; i++) {                     
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
			
			  
			  
			for(int qi=0;qi<N;qi++)
			{

				docnum_array[i][qi]=compact_score_array[qi].doc;
				  
			   printTextToFile(docnum_array[i][qi]);     
			}
			
			*  
			  
			int min = Math.min(N, compact_score_array.length);        
					
			  
			ArrayList<HashMap<String,Double>> DocN_list = new ArrayList<HashMap<String,Double>>(); 
			for(int n=0;n<min;n++)
			{
				HashMap<String,Double> Vecter_DocN = new HashMap<String,Double>();
				int docNum=compact_score_array[n].doc;                                  
				double dl= doc_length_map.get(docNum); 
				Terms terms = reader.getTermVector(docNum, docContentField);	        
				TermsEnum iterator = terms.iterator(null);                                  
				BytesRef byteRef = null;
				while ((byteRef = iterator.next()) != null) {  
				    String term = new String(byteRef.bytes, byteRef.offset, byteRef.length);
					double tf=iterator.totalTermFreq();                                 
					tf=tf  
					double df=term_doc_freq_map.get(term);                              
					double IDF = log2((numDocs - df + 0.5)   
					Vecter_DocN.put(term, tf*IDF  
				}
				DocN_list.add(Vecter_DocN);                                             
			}
			  
			HashMap<String,Double> sum_DocN=new HashMap<String,Double>();
			for(int m=0;m<DocN_list.size();m++)
			{
				sum_DocN=addVecter(sum_DocN, DocN_list.get(m),1.0,1.0);                         
			}			
						
			  
			HashMap<String,Double> pre_DocN2=new HashMap<String,Double>();
			pre_DocN2=sortVecter(sum_DocN);                                                     
			  
			  
			  
			  
			File file=new File(jsonName);   
			String content= FileUtils.readFileToString(file,"UTF-8"); 
			JSONArray jsonArray=new JSONArray(content);        			
			JSONArray termsco=jsonArray.getJSONArray(i);
			HashMap<String,Double> pre_DocN1 = new HashMap<String,Double>();
			
			
			for(int t=0;t<termsco.length();t++){
				
				JSONArray array=termsco.getJSONArray(t);
				
				Term term = new Term(docContentField,array.getString(0));     
				if(atomicReader.termDocsEnum(term)!=null)                   
				{				
					pre_DocN1.put(array.getString(0),array.getDouble(1));
				}
				else
				{
					  
				}
				
				  
			}
			  
			  

			pre_DocN2=guiYiHua1(pre_DocN2); 
			pre_DocN2=sortVecter(pre_DocN2);    
			HashMap<String,Double> pre_DocN = new HashMap<String,Double>();
			pre_DocN1=guiYiHua1(pre_DocN1);      
			pre_DocN=searchVecter(pre_DocN2, pre_DocN1, 1.0-beta, beta);  
			
			pre_DocN=sortVecter1(pre_DocN);    
			
			HashMap<String,Double> q_Doc=new HashMap<String,Double>();
			for(Term t:term_list)              
			{
				q_Doc.put(t.text(), 1.0);                                               
			}
			 
			pre_DocN=addVecter(q_Doc,pre_DocN,1-alpha,alpha);
			pre_DocN=guiYiHua1(pre_DocN);
			
			
			
			  
			ScoreDoc[] compact_score_array2 = PFB_QueryExe(atomicReader , pre_DocN);    
			
			CreateResult(compact_score_array2,qualityQueries,i,query,stats);            
			  
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
	 
	public HashMap<String, Double> sortVecter1(HashMap<String,Double> A)  {
		List<HashMap.Entry<String,Double>> list = new ArrayList<HashMap.Entry<String,Double>>(A.entrySet());
		Collections.sort(list,new Comparator<HashMap.Entry<String,Double>>() {
            public int compare(Entry<String, Double> o1,
                    Entry<String, Double> o2) {
                return (o1.getValue().compareTo(o2.getValue()))*-1;   
            }
		});
		int min=Math.min(N1, list.size());
		HashMap<String,Double> C=new HashMap<String,Double>();
    	double dmax = list.get(0).getValue();
    	  
		for(int i=0;i<min;i++)    
		{
			  
			if (Double.isNaN(dmax)){
				C.put(list.get(i).getKey(), 0.0);	
			}else{
				  
				C.put(list.get(i).getKey(), list.get(i).getValue()  
			}
		} 
		return C;
	}
	
	
	public HashMap<String, Double> guiYiHua1(HashMap<String,Double> A)  {
		List<HashMap.Entry<String,Double>> list = new ArrayList<HashMap.Entry<String,Double>>(A.entrySet());
		Collections.sort(list,new Comparator<HashMap.Entry<String,Double>>() {
            public int compare(Entry<String, Double> o1,
                    Entry<String, Double> o2) {
                return (o1.getValue().compareTo(o2.getValue()))*-1;   
            }
		});
		HashMap<String,Double> C=new HashMap<String,Double>();
		double dmax = list.get(0).getValue();
		for(int i=0;i<A.size();i++)    
		{
			  
			if (Double.isNaN(dmax)){
				C.put(list.get(i).getKey(), 0.0);	
			}else{
				  
				C.put(list.get(i).getKey(), list.get(i).getValue()  
			}
		} 
		return C;
	}
	
	public HashMap<String, Double> guiYiHua2(HashMap<String,Double> A)  {
		List<HashMap.Entry<String,Double>> list = new ArrayList<HashMap.Entry<String,Double>>(A.entrySet());
		Collections.sort(list,new Comparator<HashMap.Entry<String,Double>>() {
            public int compare(Entry<String, Double> o1,
                    Entry<String, Double> o2) {
                return (o1.getValue().compareTo(o2.getValue()))*-1;   
            }
		});
		int min=list.size();
		HashMap<String,Double> C=new HashMap<String,Double>();
		double dmax = list.get(0).getValue();
	  	double dmin = list.get(min-1).getValue();
		for(int i=0;i<A.size();i++)    
		{
			  
			if (Double.isNaN(dmax)){
				C.put(list.get(i).getKey(), 0.0);	
			}else{
				C.put(list.get(i).getKey(), (list.get(i).getValue()-dmin)  
				  
			}
		} 
		return C;
	}

	
	public HashMap<String, Double> guiYiHua3(HashMap<String,Double> A)  {
		List<HashMap.Entry<String,Double>> list = new ArrayList<HashMap.Entry<String,Double>>(A.entrySet());
		Collections.sort(list,new Comparator<HashMap.Entry<String,Double>>() {
            public int compare(Entry<String, Double> o1,
                    Entry<String, Double> o2) {
                return (o1.getValue().compareTo(o2.getValue()))*-1;   
            }
		});
		int min=list.size();
		HashMap<String,Double> C=new HashMap<String,Double>();
		double sum=0.0;
		for(int i=0;i<min;i++) 
		{
			sum+=list.get(i).getValue()*list.get(i).getValue();
		}
		sum=Math.sqrt(sum);
		for(int i=0;i<min;i++)    
		{
			C.put(list.get(i).getKey(), list.get(i).getValue()  
		}    
		return C;
	}
	
	
	  
	public ScoreDoc[] PFB_QueryExe(AtomicReader atomicReader ,HashMap<String,Double>  pre_DocN) throws IOException{
		int k=0;
		int[][] freq_array = new int[pre_DocN.size()][reader.numDocs()];
		Set<Integer> doc_set = new HashSet<Integer>();
		ArrayList<String> query_list = new ArrayList<String>(); 
		for(Entry<String,Double> pre : pre_DocN.entrySet()){
			Term term = new Term(docContentField,pre.getKey());
			query_list.add(pre.getKey());
			DocsEnum	docsEnum = atomicReader.termDocsEnum(term);	                  
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
				double total_score = 0.0f;                                            
				int doc_length = doc_length_map.get(j);                               
				double K = k1 * ((1 - b) + b * doc_length   
				for (int i = 0; i < query_list.size(); i++) {                         
					double qtf=pre_DocN.get(query_list.get(i));
					int df = term_doc_freq_map.get(query_list.get(i));                
					double tf = freq_array[i][j];                                     
					double TF = (k1 + 1) * tf   
					double IDF = log2((numDocs - df + 0.5)   
					total_score += TF * IDF*qtf;                                      
				} 
				  
			compact_score_array[kk++] = new ScoreDoc(j, (float) total_score);         
			}
		}
		Arrays.sort(compact_score_array, new ByWeightComparator());                   
		return compact_score_array;
	}
	
	  
	
	  
	public void CreateResult(ScoreDoc[] compact_score_array1, QualityQuery[] qualityQueries,int i,Query query,QualityStats[] stats) throws IOException{
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
	
	
	  
	public HashMap<String, Double> sortVecter(HashMap<String,Double> A)  {
		List<HashMap.Entry<String,Double>> list = new ArrayList<HashMap.Entry<String,Double>>(A.entrySet());
		Collections.sort(list,new Comparator<HashMap.Entry<String,Double>>() {
              
            public int compare(Entry<String, Double> o1,
                    Entry<String, Double> o2) {
                return (o1.getValue().compareTo(o2.getValue()))*-1;
            }
		});
		HashMap<String,Double> C=new HashMap<String,Double>();
		int min=Math.min(N2, list.size());
		  
		for(int i=0;i<min;i++)
		{
			 C.put(list.get(i).getKey(), list.get(i).getValue());			
			   
		}
		
		return C;
	}
	
	  
	public HashMap<String,Double> addVecter(HashMap<String,Double> A, HashMap<String,Double> B,double alpha,double beta) throws Exception {
		HashMap<String,Double> C=new HashMap<String,Double>();
		for(Entry<String,Double> a : A.entrySet())
		{
			if(B.containsKey(a.getKey()))
			{	
				if(a.getValue()*alpha+B.get(a.getKey())*beta!=0)
				C.put(a.getKey(), (a.getValue()*alpha+B.get(a.getKey())*beta));
			}
			else 
				if(a.getValue()*alpha!=0)
				C.put(a.getKey(), a.getValue()*alpha);			
		}
		for(Entry<String,Double> b : B.entrySet())
		{
			if(!C.containsKey(b.getKey()))
			{		
				if(b.getValue()*beta!=0)
				C.put(b.getKey(), b.getValue()*beta);
			}					
		}
		return C;
	}
	
	public HashMap<String,Double> searchVecter(HashMap<String,Double> A, HashMap<String,Double> B,double al,double be) throws Exception {
		HashMap<String,Double> C=new HashMap<String,Double>();
		for(Entry<String,Double> a : A.entrySet())
		{
			if(B.containsKey(a.getKey()))
			{	
				if(a.getValue()*al+B.get(a.getKey())*be!=0)
				C.put(a.getKey(), (a.getValue()*al+B.get(a.getKey())*be));
			}
			else 
				if(a.getValue()*al!=0)
				C.put(a.getKey(), a.getValue()*al);			
		}

		return C;
	}
	
	
	
	
	
	public static double round(double n) {
		return (Math.round(n*100)  
	}
	
    
	  
	public static double log2(double n) {
		return (Math.log(n)   
	}
	
	  
	public static void printTextToFile(int docnum){
        try{
            File file = new File("Roc_docnum_first_retrieval  
			PrintStream printStream = new PrintStream(new FileOutputStream(file,true),true);
			printStream.println(docnum);
        }catch(Exception ex){
            ex.printStackTrace();
        }
    }
    
}
