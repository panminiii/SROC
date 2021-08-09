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

import sigir.kernels.kernel;

public class Retrieval_KRM3 {
	
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
	public static String[] dataName={"AP90_KRM3","AP8889_KRM3","disk12_KRM3","disk45_KRM3","wt2g_KRM3","wt10g_KRM3","FT_KRM3","FBIS_KRM3","LA_KRM3","SJMN_KRM3","WSJ_KRM3","GOV2_KRM3"};

	
	
	int numDocs;
	int maxResults = 1000;
	static int N=10;
	  
	static int N1=20;      
	  
	static double alpha=0.1;
	public static double beta=0.1;
	static double k1 = 1.2;
	static double b = 0.25;
	static double k3 = 8.0;
	static double sigma;
	
	long total_term_freq = 0;
	double avg_doc_length;
	
	static int mu=600;
	
	public Retrieval_KRM3() {
	}

	public static void main(String[] args) throws Throwable {
		String kernelName="gaussKernel";         
		kernel newKernel=(kernel)Class.forName("sigir.kernels."+kernelName).newInstance();
		for (int i=0;i<1;i+=1) {
			for (N1=30;N1<=30;N1+=10) {
				if (N1==40)	continue;
				for (int j=1;j<=1;j+=1) {
					switch(i)
					{
					  
					case 0:mu=400;b=0.0 ;break;  
					case 1:mu=750;b=0.0 ;alpha=0.6;beta=0.2;break;
					case 2:mu=750;b=0.0 ;alpha=0.6;beta=0.2;break;
					  
					case 3:mu=300;b=0.35;alpha=0.6;beta=0.3;break;
					case 4:mu=1400;b=0.25;alpha=0.5;beta=0.5;break;
					case 5:mu=1600;b=0.2;alpha=0.2;beta=1.0;break;
					case 6:mu=1450;b=0.3;alpha=0.9;beta=0.1;break;
					case 7:mu=350;b=0.05 ;alpha=011.6;beta=0.3;break;
					case 8:mu=850;b=0.3;alpha=0.8;beta=0.1;break;
					case 9:mu=350;b=0.55;alpha=0.6;beta=0.4;break;
					case 10:mu=1200;b=0.3;alpha=0.7;beta=0.2;break;
					case 11:mu=600;b=0.4;alpha=0.3;beta=0.7;break;
					}
					switch(j)
					{                       
						case 0:sigma=10;break;
						case 1:sigma=25;break;
						case 2:sigma=50;break;
						case 3:sigma=80;break;
						case 4:sigma=100;break;
						case 5:sigma=200;break;
						case 6:sigma=500;break;
						case 7:sigma=1000;break;
						case 8:sigma=1500;break;
					}    
					
					for(alpha=0.8;alpha<=0.8;alpha=round(alpha+0.1)){
						for(beta=0.0;beta<=1.0;beta=round(beta+0.05)){
							String resultFile="result4  
							String reportFile="report  
							System.out.println(resultFile);    
							new Retrieval_KRM3().RunCRTER3(topics[i],qrels[i],index[i],newKernel, sigma,resultFile, reportFile);
						}
					}	
				}
			}
		}
	}

	public void RunCRTER3(String topicFile, String qrelsFile, String indexFile,kernel myKernel, double sigma,String resultFile, String reportFile) throws Exception {
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
		total_term_freq=0;
		Fields fields = MultiFields.getFields(reader);         
		Terms terms = fields.terms(docContentField);           
		TermsEnum iterator = terms.iterator(null);             
		BytesRef byteRef = null;
		while ((byteRef = iterator.next()) != null) {          
			String term = new String(byteRef.bytes, byteRef.offset, byteRef.length);
			term_total_freq_map.put(term, iterator.totalTermFreq());
			term_doc_freq_map.put(term, iterator.docFreq());
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

	  
	public void CrosstermStats(ArrayList<Term> term_list,Set<Integer> BM25_doc_set) throws Exception {
		AtomicReader atomicReader = SlowCompositeReaderWrapper.wrap(reader);
			within_query_freq_map.clear();
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
				term_list.add(term);
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
			
			  
			ScoreDoc[] compact_score_array_LM=QueryExe(doc_set,term_list,freq_array);     
			int doc_min=Math.min(N, compact_score_array_LM.length);          
			
			  
			double doc_score_Sum=0.0;                                                  
			for(int n=0;n<doc_min;n++)
			{
				doc_score_Sum+=Math.pow(2, compact_score_array_LM[n].score);
				
			}
			  
			ArrayList<HashMap<String,Double>> DocN_list1 = new ArrayList<HashMap<String,Double>>(); 
			
			for(int n=0;n<doc_min;n++)
			{
				HashMap<String,Double> Vecter_DocN = new HashMap<String,Double>();
				int docNum=compact_score_array_LM[n].doc;                                
				double doc_score=Math.pow(2, compact_score_array_LM[n].score)  
				
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
				DocN_list1.add(Vecter_DocN);                                           
			}
			
			HashMap<String,Double> sum_DocN1=new HashMap<String,Double>();
			for(int m=0;m<DocN_list1.size();m++)
			{
				sum_DocN1=addVecter(sum_DocN1, DocN_list1.get(m),1.0,1.0);                         
			}
			
			DocN_list1.clear();
			Set<Integer> LM_doc_set = new HashSet<Integer>();                
			for(int doc_i=0;doc_i < doc_min ;doc_i++)
			{
				LM_doc_set.add(compact_score_array_LM[doc_i].doc);
			}
			
			  
			ArrayList<Term> term_list1=new ArrayList<Term>();                                
			for(int n=0;n<doc_min;n++)
			{
				int docNum=compact_score_array_LM[n].doc;                                
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
			  
									
			  
			  
			  
			CrosstermStats(term_list1,LM_doc_set);                                      
			  
			ArrayList<HashMap<String,Double>> DocN_list = new ArrayList<HashMap<String,Double>>();    
			  
			for (int n = 0; n < doc_min; n++) {                                            
				HashMap<String,Double> Vecter_DocN = new HashMap<String,Double>();
				int j=compact_score_array_LM[n].doc;  
				
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
			DocN_list.clear();
			
			HashMap<String,Double> q_Doc=new HashMap<String,Double>();
			for(Term t:term_list)              
			{
				
				q_Doc.put(t.text(), 1.0  
			}
			
			HashMap<String,Double> pre_DocN1=new HashMap<String,Double>();
			HashMap<String,Double> pre_DocN=new HashMap<String,Double>();
			pre_DocN1=sortVecter5(sum_DocN1); 
			pre_DocN=sortVecter2(sum_DocN);                                             
			
			pre_DocN1=addVecter(pre_DocN1,pre_DocN,(1.0-beta),beta);                      
			pre_DocN1=sortVecter(pre_DocN1);

			HashMap<String,Double> pre_Doc=new HashMap<String,Double>();
			pre_Doc=addVecter(q_Doc,pre_DocN1,(1.0-alpha),alpha);

			  
			ScoreDoc[] compact_score_arrayKRoc=PFB_QueryExe(atomicReader , pre_Doc);
		
		
			  
			int max_resultKRoc = Math.min(maxResults, compact_score_arrayKRoc.length);
			ScoreDoc[] score_docKRoc = new ScoreDoc[max_resultKRoc];
			System.arraycopy(compact_score_arrayKRoc, 0, score_docKRoc, 0, max_resultKRoc);
			
			TopDocs tdCRoc = new TopDocs(max_resultKRoc, score_docKRoc, (float) score_docKRoc[0].score);
			
			statsCRoc[i] = analyzeQueryResults(qualityQueries[i], query, tdCRoc, judge, qualityLogCRoc, 1);
			  
			  
		}

		QualityStats avgKRoc = QualityStats.average(statsCRoc);
		avgKRoc.log("SUMMARY", 2, qualityLogCRoc, "  ");
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
				double total_score = 0.0;                                             
				
				int dl = doc_length_map.get(j);                                       
				for (int k = 0; k < term_list.size(); k++) {                          
					Long ctf= term_total_freq_map.get(term_list.get(k).text());       
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
					total_score +=log2(1.0*tf  
					  
					  
				}  
			compact_score_array[kk++] = new ScoreDoc(j, (float) (total_score));       
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
	
	  
	public HashMap<String, Double> sortVecter5(HashMap<String,Double> A)  {
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
			sum+= list.get(i).getValue();	                                          
		}
		for(int i=0;i<min;i++)
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

