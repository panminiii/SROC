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

public class MyRocchio {                                                      

	public  double k1 = 1.2;                                                
	public  double b = 0.35;                                               
	public  double k3 = 8.0;
	public  int N = 10;                                                       
	public  double beta =0.25;                                                  
	public  double aphla=1.0;                                                     
	public  int N1 = 30;                                                         
	public  double K = 0.6;                                                    
	public  double L = 1;                                                      


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

	  
	public static double log2(double n) {
		return (Math.log(n)   
	}

	public MyRocchio() {
	}

	public MyRocchio(double b) {
		this.b = b;
	}
	
	public MyRocchio(double b, int N, int N1, double beta) {
		this.b = b;
		this.N = N;
		this.N1 = N1;
		this.beta = beta;
	}
	  
	public static void main(String[] args) throws Exception{
		  
		String usage =
				"Usage:\tjava BatchSearch [-index indexFile] [-topics topicsFile] [-qrels qrelsFile] [-data dataName] [-b b_value]"
				+ "[-docN N_value][-termN N1_value][-beta beta_value]";
		if (args.length > 0 && ("-h".equals(args[0]) || "-help".equals(args[0]))) {
			System.out.println(usage);
			System.exit(0);
		}
		
		String topicsFile = "topics  
		String qrelsFile = "qrels  
		String indexFile = "index_wt2g";
		String dataName = "wt2g";
		double b_value=0.25;
		int N_value=10;
		int N1_value=30;
		double beta_value=0.9;
		
		for(int i = 0;i < args.length;i++) {
			if ("-b".equals(args[i])) {
				b_value =  Double.parseDouble(args[i+1]);
				i++;
			}else if ("-beta".equals(args[i])) {
				beta_value =  Double.parseDouble(args[i+1]);
				i++;
			}else if ("-docN".equals(args[i])) {
				N_value =  Integer.parseInt(args[i+1]);
				i++;
			}else if ("-termN".equals(args[i])) {
				N1_value =  Integer.parseInt(args[i+1]);
				i++;
			} 		
			else if ("-index".equals(args[i])) {
				indexFile = args[i+1];
				i++;
			} else if ("-topics".equals(args[i])) {
				topicsFile = args[i+1];
				i++;
			} else if ("-qrels".equals(args[i])) {
				qrelsFile = args[i+1];
				i++;
			} else if ("-data".equals(args[i])) {
				dataName = args[i+1];
				i++;
			}
		}
		String resultFile = "retrievalresult  
		String reportFile = "retrievalresult  
		new MyRocchio(b_value,  N_value,  N1_value,  beta_value).RunRocchio(topicsFile, qrelsFile, indexFile,  resultFile,  reportFile);
	}
	  
	public void RunRocchio(String topicFile, String qrelsFile, String indexFile, String resultFile, String reportFile) throws Exception {
		System.out.println("RunRocchio...");
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
		PrintWriter pw=new PrintWriter(new File(reportFile), "UTF-8");                
		submitRep = new SubmissionReport(pw, "TEST");                                 
		execute();                                                                    
		directory.close();                                                            
		qualityLog.close();                                                           
		pw.close();                                                                   
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
				if(docsEnum!=null)
				while ((docsEnum.nextDoc()) != DocIdSetIterator.NO_MORE_DOCS) {
					int doc_id = docsEnum.docID();                                    
					  
					doc_set.add(doc_id);                                              
					freq_array[k][doc_id] = docsEnum.freq();                          
				}
			}
			  
			numDocs = reader.numDocs();                                               

			 
			  
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
					tf = tf  
					double df=term_doc_freq_map.get(term);                              
					double IDF = log2((numDocs - df + 0.5)   
					Vecter_DocN.put(term, tf*IDF*beta  
					  
					  
					  

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
			HashMap<String,Double> q_Doc=new HashMap<String,Double>();
			for(Term t:term_list)              
			{
				q_Doc.put(t.text(), 1.0*aphla);                                           
				  
			}

			  
			  
			  
			pre_DocN=addVecter(pre_DocN,q_Doc);                                     
			  
			  
			  
			ScoreDoc[] compact_score_array1 = PFB_QueryExe(atomicReader , pre_DocN);

			  
			CreateResult(compact_score_array1,qualityQueries,i,query,stats);
		}
		  
		WriteResult(stats);           
	}
	
	

	  
	public void printVecter(HashMap<String,Double> A){
		for(Entry<String,Double> a : A.entrySet())
		{
			System.err.println("name=:"+a.getKey()+" value=:"+a.getValue());
		}
	}

	  
	 * @return **  
		List<HashMap.Entry<String,Double>> list = new ArrayList<HashMap.Entry<String,Double>>(A.entrySet());
		Collections.sort(list,new Comparator<HashMap.Entry<String,Double>>() {
			  
			public int compare(Entry<String, Double> o1,
					Entry<String, Double> o2) {
				return (o1.getValue().compareTo(o2.getValue()))*-1;
			}
		});
  
        int min = Math.min(N1, list.size());
		HashMap<String,Double> C=new HashMap<String,Double>();
		for(int i=0;i<min;i++)
		{
			C.put(list.get(i).getKey(), list.get(i).getValue());			
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

	  
	public ScoreDoc[] QueryExe(Set<Integer> doc_set,ArrayList<Term> term_list,int[][] freq_array){
		int kk = 0;	
		ScoreDoc[] compact_score_array=new ScoreDoc[doc_set.size()];
		for (int j = 0; j < numDocs; j++) {                                       
			if (doc_set.contains(j)) {                                            
				double total_score = 0.0f;                                        
				int doc_length = doc_length_map.get(j);                           
				double K = k1 * ((1 - b) + b * doc_length   
				for (int k = 0; k < term_list.size(); k++) {                      
					
					Integer df = term_doc_freq_map.get(term_list.get(k).text());      
					if(df!=null&&df!=0){
						double tf = freq_array[k][j];	                              
					    double TF = (k1 + 1) * tf   
					    double IDF = log2((numDocs - df + 0.5)   
					    total_score += TF * IDF;
					      
				    }
				}   
				compact_score_array[kk++] = new ScoreDoc(j, (float) total_score);     
			}
		}
		Arrays.sort(compact_score_array, new ByWeightComparator());               
		return compact_score_array;
	}

	  
	public double[] ComputeQueryWeight(ArrayList<Term> term_list,HashMap<Integer, Double> Doc2DocVecterd,int[][] freq_array){
		double[] queryWeight = new double[term_list.size()];
		for (int k = 0; k < term_list.size(); k++) {
			queryWeight[k]=0.0;
			for (Entry<Integer, Double> docSim : Doc2DocVecterd.entrySet()) {
				int j=docSim.getKey();
				int doc_length = doc_length_map.get(j);
				double tf = freq_array[k][j];                                     
				int df = term_doc_freq_map.get(term_list.get(k).text());          
				double K = k1 * ((1 - b) + b * doc_length   
				double TF = (k1 + 1) * tf   
				double IDF = log2((numDocs - df + 0.5)   
				queryWeight[k] += TF * IDF*docSim.getValue();                     
			}
			  
			  
		}
		return queryWeight;
	}

	  
	  
	 * @param pre_Docn 
	 * @throws IOException **  

		int[][] freq_array = new int[pre_DocN.size()][reader.numDocs()];
		int k=0;
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
				double total_score = 0.0f;                                        
				int doc_length = doc_length_map.get(j);                           
				double K = k1 * ((1 - b) + b * doc_length   
				for (int i = 0; i < query_list.size(); i++) {                      
					double qtf=pre_DocN.get(query_list.get(i));
					Integer df = term_doc_freq_map.get(query_list.get(i));            
					if(df!=null&&df!=0){
						double tf = freq_array[i][j];                                 
					    double TF = (k1 + 1) * tf   
					    double IDF = log2((numDocs - df + 0.5)   

					      
					    total_score += TF * IDF*qtf;
					  
				    }
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
		submitRep.report(qualityQueries[i], td, docNameField, searcher);
		submitRep.flush();
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
