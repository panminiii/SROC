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
import java.util.Set;
import java.util.Map.Entry;

import org.apache.commons.io.FileUtils;
import org.apache.lucene.benchmark.quality.Judge;
import org.apache.lucene.benchmark.quality.QualityQuery;
import org.apache.lucene.benchmark.quality.QualityQueryParser;
  
import org.apache.lucene.benchmark.quality.trec.TrecJudge;
import org.apache.lucene.benchmark.quality.trec.TrecTopicsReader;
import org.apache.lucene.benchmark.quality.utils.DocNameExtractor;
import org.apache.lucene.index.AtomicReader;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.DocsAndPositionsEnum;
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

import sigir.kernels.kernel;

public class Retrieval_KRoc {
	
	private static final String docNameField = "docno";
	private static final String docContentField = "contents";
	QualityQuery[] qualityQueries;
	QualityQueryParser qqParser;	  
	IndexReader reader;
	IndexSearcher searcher;
	Judge judge;
	kernel myKernel;
	  
	PrintWriter qualityLogCRoc;
	
	  
  
  
  
  
  
	public static String[] topics={"topics  N.51-150","topics  
	public static String[] qrels={"qrels  els.WSJ.151-200","qrels  
	public static String[] index={"index_AP90","index_AP8889","index_disk12","index_disk45","index_wt2g","index_wt10g","index_FT","index_FBIS","index_LA","index_SJMN","index_WSJ","index_GOV2_stem"}; 
	public static String[] dataName={"3AP90_KRoc","4AP8889_KRoc","disk12_KRoc","disk45_KRoc","wt2g_KRoc","wt10g_KRoc","FT_KRoc","FBIS_KRoc","LA_KRoc","SJMN_KRoc","WSJ_KRoc","GOV2_KRoc"};
	
	int numDocs;
	int maxResults = 1000;
	double avg_doc_length;
	static int N=10;
	  
	static int N1=30;
	  
	static double alpha=1.0;
	static double beta=0.1;
	static double deta=0.1;
	static double k1 = 1.2;
	static double b = 0.25;
	static double k3 = 8.0;
	static double sigma=50;
	
	public Retrieval_KRoc() {
	}

	public static void main(String[] args) throws Throwable {
		String kernelName="gaussKernel";         
		kernel newKernel=(kernel)Class.forName("sigir.kernels."+kernelName).newInstance();
		for (int i=1;i<=1;i+=1) {
			for (N1=10;N1<=50;N1+=10) {					
						b=0.55 ;
						  
						  
						sigma=25;					
					for(alpha=0.4;alpha<=0.4;alpha=round(alpha+0.1)){
						for(beta=0.4;beta<=0.4;beta=round(beta+0.1)){
							for(deta=0.0;deta<=1.0;deta=round(deta+0.1)){
								b=0.0;k1=1.2;N=8;
							String resultFile="test1  
							String reportFile="report  
							  
							System.out.println(resultFile);    
							new Retrieval_KRoc().RunCRTER2(topics[i],qrels[i],index[i],newKernel, sigma,resultFile, reportFile);
						}
					}
				}	
			}
		}
	}

	public void RunCRTER2(String topicFile, String qrelsFile, String indexFile,kernel myKernel, double sigma,String resultFile, String reportFile) throws Exception {
		this.myKernel=myKernel;
		myKernel.setParameter(sigma);
		Directory directory = FSDirectory.open(new File(indexFile));
		reader = DirectoryReader.open(directory);
		searcher = new IndexSearcher(reader);
		qualityLogCRoc = new PrintWriter(new File(resultFile), "UTF-8");
		TrecTopicsReader qReader = new TrecTopicsReader();	  
		qualityQueries = qReader.readQueries(new BufferedReader(new FileReader(new File(topicFile))));	  
		  
		judge = new TrecJudge(new BufferedReader(new FileReader(new File(qrelsFile))));	  
		judge.validateData(qualityQueries, qualityLogCRoc);   
		qqParser = new MyQQParser("title", "contents");   
		  
		execute();      
		directory.close();
		qualityLogCRoc.close();
	}

	HashMap<String, Long> term_total_freq_map = new HashMap<String, Long>();		  
	HashMap<String, Integer> term_doc_freq_map = new HashMap<String, Integer>();	  
	HashMap<Integer, Integer> doc_length_map = new HashMap<Integer, Integer>();		  
	HashMap<Integer, Double> doc_avg_tf_map = new HashMap<Integer, Double>();		  
	HashMap<String, HashMap<Integer, ArrayList<Integer>>> within_query_freq_map = new HashMap<String, HashMap<Integer, ArrayList<Integer>>>();   
	
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
			String term = new String(byteRef.bytes, byteRef.offset, byteRef.length);
			term_total_freq_map.put(term, iterator.totalTermFreq());
			term_doc_freq_map.put(term, iterator.docFreq());
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

	  
	public void CrosstermStats(ArrayList<Term> term_list,Set<Integer> BM25_doc_set) throws Exception {
		AtomicReader atomicReader = SlowCompositeReaderWrapper.wrap(reader);
			for (int k = 0; k < term_list.size(); k++) {
				Term term = term_list.get(k);                                   
				String termText = term.text();
				DocsAndPositionsEnum docsAndPositionsEnum = atomicReader.termPositionsEnum(term);     
				if (docsAndPositionsEnum == null) {
					continue;
				}               
				int doc_id;
				HashMap<Integer, ArrayList<Integer>> Docid_position_map = new HashMap<Integer, ArrayList<Integer>>();       
				while ((doc_id = docsAndPositionsEnum.nextDoc()) != DocIdSetIterator.NO_MORE_DOCS) {
					if(!BM25_doc_set.contains(doc_id)) continue;    
					int freq = docsAndPositionsEnum.freq();         
					int position;
					ArrayList<Integer> query_terms_position = new ArrayList<Integer>();    
					for (int j = 0; j < freq; j++) {
						position = docsAndPositionsEnum.nextPosition();
						query_terms_position.add(position);                            
					}
					Docid_position_map.put(doc_id, query_terms_position);         
				}
				within_query_freq_map.put(termText, Docid_position_map);       
			}
	}

	public void search() throws Exception {
		AtomicReader atomicReader = SlowCompositeReaderWrapper.wrap(reader);
		QualityStats statsCRoc[] = new QualityStats[qualityQueries.length];
		  
		for (int i = 0; i < qualityQueries.length; i++) {
			QualityQuery qq = qualityQueries[i];
			  
			Query query = qqParser.parse(qq);
			  
			HashSet<Term> term_set = new HashSet<Term>();
			query.extractTerms(term_set);
			ArrayList<Term> term_list = new ArrayList<Term>();
			for (Term term : term_set) {
				DocsEnum docsEnum = atomicReader.termDocsEnum(term);	              
				if(docsEnum!=null)
				{
					term_list.add(term);
				  
				}
				
			}
			  
			Set<Integer> doc_set = new HashSet<Integer>();
			int[][] freq_array = new int[term_list.size()][reader.numDocs()];
			for (int k = 0; k < term_list.size(); k++) {
				Term term = term_list.get(k);
				DocsEnum docsEnum = atomicReader.termDocsEnum(term);
				if (docsEnum == null) {
					continue;	
				}
				while (docsEnum.nextDoc() != DocIdSetIterator.NO_MORE_DOCS) {
					int doc_id = docsEnum.docID();
					doc_set.add(doc_id);
					freq_array[k][doc_id] = docsEnum.freq();
				}
			}
			numDocs = reader.numDocs();
			
			  
			ScoreDoc[] compact_score_array_BM25=QueryExe(doc_set,term_list,freq_array);     
			int doc_min=Math.min(N, compact_score_array_BM25.length);          
			Set<Integer> BM25_doc_set = new HashSet<Integer>();                
			for(int doc_i=0;doc_i < doc_min ;doc_i++)
			{
				BM25_doc_set.add(compact_score_array_BM25[doc_i].doc);
			}
			
			  
			ArrayList<Term> term_list1=new ArrayList<Term>();                                
			for(int n=0;n<doc_min;n++)
			{
				int docNum=compact_score_array_BM25[n].doc;                                
				Terms terms = reader.getTermVector(docNum, docContentField);	           
				TermsEnum iterator = terms.iterator(null);                                  
				BytesRef byteRef = null;
				while ((byteRef = iterator.next()) != null) {  
				    String term_text = new String(byteRef.bytes, byteRef.offset, byteRef.length);
				    Term term = new Term(docContentField,term_text );
				    if(!term_list1.contains(term))
				    	term_list1.add(term);                                                 
				}
			}
			
			  
			  
			ArrayList<HashMap<String,Double>> DocN_list1 = new ArrayList<HashMap<String,Double>>(); 
			for(int n=0;n<doc_min;n++)
			{
				HashMap<String,Double> Vecter_DocN1 = new HashMap<String,Double>();
				int j=compact_score_array_BM25[n].doc;                                
				double dl= doc_length_map.get(j); 
				  
				Terms terms = reader.getTermVector(j, docContentField);	      
				TermsEnum iterator = terms.iterator(null);                                  
				BytesRef byteRef = null;
				while ((byteRef = iterator.next()) != null) {  
				    String term = new String(byteRef.bytes, byteRef.offset, byteRef.length);
					double tf=iterator.totalTermFreq();                                 
					tf=tf  
					double df=term_doc_freq_map.get(term);                              
					  
					double IDF = log2((numDocs - df + 0.5)   
					Vecter_DocN1.put(term, tf*IDF  
					  
					
				}
				DocN_list1.add(Vecter_DocN1);                                             
			}
			HashMap<String,Double> sum_DocN1=new HashMap<String,Double>();
			for(int m=0;m<DocN_list1.size();m++)
			{
				sum_DocN1=addVecter(sum_DocN1, DocN_list1.get(m),1.0,1.0);                         
			}			
						
			  
			  
			  
			CrosstermStats(term_list1,BM25_doc_set);                                      
			  
			ArrayList<HashMap<String,Double>> DocN_list = new ArrayList<HashMap<String,Double>>();    
			  
			for (int n = 0; n < doc_min; n++) {                                            
				HashMap<String,Double> Vecter_DocN = new HashMap<String,Double>();
				int j=compact_score_array_BM25[n].doc;  
				
				int dl = doc_length_map.get(j);
				double K = k1 * ((1 - b) + b * dl   
				ArrayList<Term> term_list_tmp=new ArrayList<Term>();
				Terms terms = reader.getTermVector(j, docContentField);	      
				TermsEnum iterator = terms.iterator(null);                                  
				BytesRef byteRef = null;
				while ((byteRef = iterator.next()) != null) {  
					String term_text = new String(byteRef.bytes, byteRef.offset, byteRef.length);
				    Term term = new Term(docContentField,term_text );
				    if(!term_list_tmp.contains(term))
				    	term_list_tmp.add(term);                                          
				}
				
				for (int ti = 0; ti < term_list_tmp.size(); ti++) {                       
					double ptf = 0.0;
					int tftiD = within_query_freq_map.get(term_list_tmp.get(ti).text()).get(j).size();
					
					for (int qj = 0; qj < term_list.size(); qj++) {                       
						double tftiqjD = 0;                                                   
						int tfqjD = freq_array[qj][j];                          
						double qtf = myKernel.intersect(1);                     
						
						if (tftiD == 0 || tfqjD == 0) {
							continue;
						}
						
						for (int tk = 0; tk < tftiD; tk++) {                    
							int positionkti = within_query_freq_map.get(term_list_tmp.get(ti).text()).get(j).get(tk);
							for (int qk = 0; qk < tfqjD; qk++) {
								int positionkqj = within_query_freq_map.get(term_list.get(qj).text()).get(j).get(qk);
								double kerneltermp = myKernel.intersect(Math.abs(positionkti - positionkqj));   
								if (kerneltermp >= Double.MIN_VALUE) {                
									tftiqjD += kerneltermp;
								} 
							}     
						} 
						  
						  
						double TFtiqjD = tftiqjD  
						  
						int qjdf=term_doc_freq_map.get(term_list.get(qj).text());
					    double IDFqj = log2((numDocs - qjdf+ 0.5)   
						double QTF = (k3 + 1) * qtf   
						ptf += TFtiqjD * IDFqj * QTF;                           
						
					}
					if(ptf  
						Vecter_DocN.put(term_list_tmp.get(ti).text(), ptf  
					}
				}
				DocN_list.add(Vecter_DocN);                                                
			}
			HashMap<String,Double> sum_DocN=new HashMap<String,Double>();
			for(int m=0;m<DocN_list.size();m++)
			{
				sum_DocN=addVecter(sum_DocN, DocN_list.get(m),1.0,1.0);                
			}			
			
			
			HashMap<String,Double> q_Doc=new HashMap<String,Double>();
			for(Term t:term_list)              
			{
				q_Doc.put(t.text(), 1.0);                                               
			}
			
			HashMap<String,Double> pre_DocN1=new HashMap<String,Double>();
			HashMap<String,Double> pre_DocN=new HashMap<String,Double>();
			pre_DocN1=sortVecter1(sum_DocN1);  
			pre_DocN=sortVecter1(sum_DocN);                                             
			  
			File file=new File("json  
			  
			String content= FileUtils.readFileToString(file,"UTF-8"); 
			JSONArray jsonArray=new JSONArray(content);        			
			JSONArray termsco=jsonArray.getJSONArray(i);
			HashMap<String,Double> pre_DocN2 = new HashMap<String,Double>();
			for(int t=0;t<termsco.length();t++){
				JSONArray array=termsco.getJSONArray(t);
				Term term = new Term(docContentField,array.getString(0));     
				if(atomicReader.termDocsEnum(term)!=null)                   
				{				
					pre_DocN2.put(array.getString(0),array.getDouble(1));
				}
				else
				{
					  
				}
				  
			}
			pre_DocN2=sortVecter1(pre_DocN2); 
			
			pre_DocN1=addVecter(pre_DocN1,pre_DocN,1.0,beta);                      
			pre_DocN1=addVecter(pre_DocN1,pre_DocN2,1.0,deta); 
			pre_DocN=addVecter(q_Doc,pre_DocN1,1.0,alpha);
			  
			  
			ScoreDoc[] compact_score_arrayKRoc=PFB_QueryExe(atomicReader , pre_DocN);

			  
			Arrays.sort(compact_score_arrayKRoc, new ByWeightComparator());
			int max_resultCRoc = Math.min(maxResults, compact_score_arrayKRoc.length);
			ScoreDoc[] score_docCRoc = new ScoreDoc[max_resultCRoc];
			System.arraycopy(compact_score_arrayKRoc, 0, score_docCRoc, 0, max_resultCRoc);
			TopDocs tdCRoc = new TopDocs(max_resultCRoc, score_docCRoc, (float) score_docCRoc[0].score);
			statsCRoc[i] = analyzeQueryResults(qualityQueries[i], query, tdCRoc, judge, qualityLogCRoc, 1);
			  
			  
			  
		}

		QualityStats avgKRoc = QualityStats.average(statsCRoc);
		avgKRoc.log("SUMMARY", 2, qualityLogCRoc, "  ");
	}

	public void printVecter(HashMap<String, Double> A) {
		for (Entry<String, Double> a : A.entrySet()) {
			System.err.println("name=:" + a.getKey() + " value=:" + a.getValue());
		}
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
	
	  
	public ScoreDoc[] QueryExe(Set<Integer> doc_set,ArrayList<Term> term_list,int[][] freq_array){
		int kk = 0;	
		ScoreDoc[] compact_score_array=new ScoreDoc[doc_set.size()];
		for (int j = 0; j < numDocs; j++) {                                       
			if (doc_set.contains(j)) {                                            
				double total_score = 0.0f;                                        
				int doc_length = doc_length_map.get(j);                           
				double K = k1 * ((1 - b) + b * doc_length   
				for (int k = 0; k < term_list.size(); k++) {
					String termText = term_list.get(k).text();
					double tf = freq_array[k][j];
					Integer dfo = term_doc_freq_map.get(termText);
					double df = (dfo == null) ? 0 : dfo;
					double TF = (k1 + 1) * tf   
					double IDF = log2((numDocs - df+ 0.5)   
					total_score += TF * IDF;
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
	
	  
	public HashMap<String, Double> sortVecter2(HashMap<String,Double> A)  {
		List<HashMap.Entry<String,Double>> list = new ArrayList<HashMap.Entry<String,Double>>(A.entrySet());
		Collections.sort(list,new Comparator<HashMap.Entry<String,Double>>() {
              
            public int compare(Entry<String, Double> o1,
                    Entry<String, Double> o2) {
                return (o1.getValue().compareTo(o2.getValue()))*-1;
            }
		});
		
		HashMap<String,Double> C=new HashMap<String,Double>();
		int min=Math.min(N1, list.size());
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
	
	  
	public HashMap<String, Double> sortVecter3(HashMap<String,Double> A)  {
		return sortVecter2(sortVecter1(A));
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
		
		int min=Math.min(N1, list.size());
		for(int i=0;i<min;i++)
		{
			 C.put(list.get(i).getKey(), list.get(i).getValue());			
		}
		
		return C;
	}
	
	  
	public static double log2(double n) {
		return (Math.log(n)   
	}
	
	  
	public static double round(double n) {
		return (Math.round(n*100)  
	}
	
	  
	public static double norm(double n) {
		return (n  
	}
}

