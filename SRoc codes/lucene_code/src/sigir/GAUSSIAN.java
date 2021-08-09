package sigir;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
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

public class GAUSSIAN {                                                      

	private double k1 = 1.2;                                                 
	private double b = 0.35; 
	private double l = 280;
	private double s = 26;
	
	public GAUSSIAN( ) { //************************师兄最新修改************************
	}
	public GAUSSIAN(double b) {
		this.b = b;
	}
	public GAUSSIAN(double l, double s) {
		this.l = l;
		this.s = s;
	}

	public GAUSSIAN(double b, double l, double s) {
		this.b = b;
		this.l = l;
		this.s = s;
	}
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
	//SubmissionReport submitRep;	
	Judge judge;	



	HashMap<String, Long> term_total_freq_map = new HashMap<String, Long>();		
	HashMap<String, Integer> term_doc_freq_map = new HashMap<String, Integer>();	
	HashMap<Integer, Integer> doc_length_map = new HashMap<Integer, Integer>();		
	HashMap<Integer, Double> doc_avg_tf_map = new HashMap<Integer, Double>();       

	/** log2函数**/
	public static double log2(double n) {
		return (Math.log(n) / Math.log(2));
	}



	/** 代码步骤0，直接调用核心算法BM25**/
	public static void main(String[] args) throws Exception{
		String[] topics = {"topics/topics.AP90.51-100","topics/topics.AP8889.51-100","topics/topics.disk12.51-200","topics/topics.disk45.301-450","topics/topics.trec8.401-450","topics/topics.301_450.601_700.robust2004","topics/topics.wt2g","topics/topics.wt10g"}; 
		String[] qrels = {"qrels/qrels.AP90.51-100","qrels/qrels.AP8889.51-100","qrels/qrels.disk12.51-200","qrels/qrels.disk45.301-450","qrels/qrels.trec8.401-450","qrels/qrels.301_450.601_700.robust2004","qrels/qrels.wt2g","qrels/qrels.wt10g"}; 
		String[] index = {"index_AP90","index_AP8889","index_disk12","index_disk45","index_trec8","index_robust2004","index_wt2g","index_wt10g"}; 
		String[] result = {"AP90_GAUSSIAN","AP8889_GAUSSIAN","disk12_GAUSSIAN","disk45_GAUSSIAN","TREC8_GAUSSIAN","ROBUST2004_GAUSSIAN","WT2G_GAUSSIAN","WT10G_GAUSSIAN"};
		for (int i =0;i <=0; i++) {
			for (double b = 0.35; b <=0.35; b += 0.05) {
				for (double l =100; l <=900; l += 100) {
					for (double s =18; s <=18 ; s +=1){
						GAUSSIAN myDSP = new GAUSSIAN();
						myDSP.b = b;
						myDSP.s = s;
						myDSP.l = l;
						
						final String resultFile = String.format("_GAUSSIAN_DSP/"+result[i] + "_k1_%.1f_b_%.2f_l_%.0f_s_%.0f.txt", myDSP.k1, myDSP.b, myDSP.l, myDSP.s );
						final String reportFile = String.format("_GAUSSIAN_DSP/"+result[i] + "_k1_%.1f_b_%.2f_l_%.0f_s_%.0f-report.txt", myDSP.k1, myDSP.b, myDSP.l, myDSP.s );
						
						//final String resultFile = String.format("result/GAUSSIAN_DSP/"+result[i] + "_k1_%.1f_b_%.2f_l_%.0f_s_%.0f.txt", myDSP.k1, myDSP.b, myDSP.l, myDSP.s );
						//final String reportFile = String.format("result/GAUSSIAN_DSP/"+result[i] + "_k1_%.1f_b_%.2f_l_%.0f_s_%.0f-report.txt", myDSP.k1, myDSP.b, myDSP.l, myDSP.s );
						
						myDSP.RunLSPR(topics[i], qrels[i], index[i], resultFile, reportFile );
					}
				}
			}
		}   
	}


	/** 代码步骤1.2，核心算法LSPR(初始化和对索引的处理),各个参数含义 见最上面的static定义**/
	public void RunLSPR(String topicFile, String qrelsFile, String indexFile, String resultFile, String reportFile ) throws Exception {
		File result = new File(resultFile);
		if (result.exists()) {
			System.err.println(String.format("结果文件[%s]已存在，跳过!!!", result.getName()));
			return;
		} else {
			result.getParentFile().mkdirs();
			System.out.println(String.format("--------------------%s--------------------", result.getName().replace(".txt", "")));
		}
		Directory directory = FSDirectory.open(new File(indexFile));                // 打开磁盘已建好的索引
		reader = DirectoryReader.open(directory);                                   // IndexReader
		searcher = new IndexSearcher(reader);                                       // IndexSearcher
		docNameField = "docno";                                                     // 文档名
		docContentField = "contents";                                               // 文档题目和内容，contents=title+text
		qualityLog = new PrintWriter(new File(resultFile), "UTF-8");                // 创建文件写入流AP90-BM25.txt，准备写入结果文件(包括文件丢失信息和每个查询的各项统计指标)，以utf-8的方式
		TrecTopicsReader qReader = new TrecTopicsReader();	                        // #1  定义内容读取器（benchmark自带）
		qualityQueries = qReader.readQueries(new BufferedReader(new FileReader(new File(topicFile))));	// #1 Read TREC topics as QualityQuery[],#1 将标准查询内容（topics.51-100 ）的内容按标签一对对的读取到qualityQueries[]数组中（按标签读取）,(queryID，HashMap<标签,内容>)
		// 数组中存放了 ID和（title，description，narrative）3个键的键值对
		judge = new TrecJudge(new BufferedReader(new FileReader(new File(qrelsFile))));	// judge得到人工标准的Q与D相关的结果qrels.51-100.disk3，（函数处理了第1列（查询序号），第3列（文档名），第4（是否相关）列
		// judge的judgements存储了HashMap<queryID,QRelJudgement>键值对,而QRelJudgement为(queryID，HashMap<docName,docName>)
		judge.validateData(qualityQueries, qualityLog);                             // #3 Verify query and Judge match  验证qualityQueries与judgements的queryID是否对应，
		//如果后者缺少对应的queryID人工标注数据会在qualityLog（AP90-BM25.txt）中写入缺失的文件信息
		qqParser = new MyQQParser("title", "contents");                             // #4 Create parser to translate queries into Lucene queries
		//		PrintWriter pw=new PrintWriter(new File(reportFile), "UTF-8");              // pan add 将下一句的参数单独拿出来为了更好的关闭
		//		submitRep = new SubmissionReport(pw, "TEST");                               // 创建文件写入流AP90-BM25-report.txt，准备写入结果文件(包括每个查询的命中文档的具体得分)，以utf-8的方式
		execute();                                                                  //termStats(); docStats();search(); 执行具体的参数统计和查询search
		directory.close();                                                          //关闭索引文件目录读取流
		qualityLog.close();                                                         //关闭统计文件写入流
		//		pw.close();                                                                 //pan add 关闭具体结果文件写入流(更加稳定)
	}	

	/** 代码步骤1.2.1，依次计算各个参数，然后进行查询**/
	public void execute() throws Exception {
		termStats();           
		docStats();
		search();
	}	

	/** 代码步骤1.2.1.1，以term的角度来进行处理，求出cf,df,tf**/
	public void termStats() throws Exception {	                                    // directory = FSDirectory.open(new File(indexFile)),indexFile="index_AP90_stem"
		Fields fields = MultiFields.getFields(reader);                              // reader(IndexReader) = DirectoryReader.open(directory)
		Terms terms = fields.terms(docContentField);                                // terms 是所有索引文件中"contents"域 的个数 共有78305个，String docContentField = "contents"   
		TermsEnum iterator = terms.iterator(null);                                  // terms.getSumTotalTermFreq(),16322984 是  所有词条（174836个）的词频ctf之和
		//terms.size();
		BytesRef byteRef = null;
		while ((byteRef = iterator.next()) != null) {                               // 循环执行 174836 次  是所有的词条信息，每次取一条（倒排）
			String text = new String(byteRef.bytes, byteRef.offset, byteRef.length);// 每一个词条的信息
			term_total_freq_map.put(text, iterator.totalTermFreq());                // 每一个词条和它的在所有文档中的出现次数cf			
			term_doc_freq_map.put(text, iterator.docFreq());	                    // 每一个词条和出现它文章的总数df	
			total_term_freq += iterator.totalTermFreq();                            // 所有词条和它的在所有文档中的出现次数tf
		}
	}	

	/** 代码步骤1.2.1.2，以doc的角度来进行处理**/
	public void docStats() throws Exception {	
		long total_dl = 0;
		for (int j = 0; j < reader.numDocs(); j++) {                                // "index_AP90_stem" reader.numDocs()索引中文档总数78305
			int docLen = 0;		
			int term_num = 0;	
			Terms terms = reader.getTermVector(j, docContentField);	                // docContentField="content",getTermVectors(docID).terms("content"),得到某一篇文档content中的倒排词典，第1,2,3,4篇各有259，88，240，156个词条
			if (terms != null && terms.size() > 0) {                                // 如果词典不为空
				TermsEnum termsEnum = terms.iterator(null);		                    // 对词典的每个词进行遍历
				while ((termsEnum.next()) != null) {
					int freq = (int) termsEnum.totalTermFreq();                     // 当前词条在当前文档中的出现次数cf
					docLen += freq;                                                 // 该文档的总词数（长度）
					term_num++;                                                     // 该文档词条的个数 等价于terms.size()
				}
			}
			total_dl += docLen;
			doc_length_map.put(j, docLen);                                          // 保存该文档的长度dl
			double avg_tf = (term_num == 0) ? 0 : ((double) docLen) / term_num;     // 求该文档的所有词的平均出现次数atf  
			doc_avg_tf_map.put(j, avg_tf);									        // 保存该文档的所有词的平均出现次数atf
		}
		avg_doc_length = ((double) total_dl) / reader.numDocs();                    // 所有文档的平均长度
	}

	/** 代码步骤1.2.1.3，计算各个查询Query与各个文档Doc的得分（核心步骤！）**/
	public void search() throws Exception {
		numDocs = reader.numDocs();                                             // 获得第文档总数N
		AtomicReader atomicReader = SlowCompositeReaderWrapper.wrap(reader);                         // fields = MultiFields.getFields(reader);liveDocs = MultiFields.getLiveDocs(reader)
		QualityStats stats[] = new QualityStats[qualityQueries.length];                              // stats为存放各个查询统计结果的数组，qualityQueries=》(queryID，HashMap<标签,内容>) 50个
		for (int i = 0; i < qualityQueries.length; i++) {                           // 对每一个查询做出如下处理
			System.err.println(i+" ");
			QualityQuery qq = qualityQueries[i];                                    // 获得上面读取好的 50条查询信息
			Query query = qqParser.parse(qq);	                                    // 构建query查询对象,System.err.println(qq.getQueryID()+qq.getValue("title"));  打印qq数组中一个元素的ID和1个键值对 
			HashSet<Term> term_set = new HashSet<Term>();                           // 构建一个词典set
			query.extractTerms(term_set);	                                        // 将每一个query 进行词条化然后将结果存放在 term_set中，term_set种就是一个查询的几个分词
			ArrayList<Term> term_list = new ArrayList<Term>();                      
			for (Term term : term_set) {
				DocsEnum docsEnum = atomicReader.termDocsEnum(term);
				if(docsEnum!=null){
					term_list.add(term);
				}
			}
			Set<Integer> doc_set = new HashSet<Integer>();	
			int[][] freq_array = new int[term_list.size()][reader.numDocs()];       // 构建一个二维数组有 词条的数目行（第一个查询2行），有78305列，存放每个词条在每个文档中出现的次数
			double maxIDF_query = 0.0;

			for (int k = 0; k < term_list.size(); k++) {
				Term term = term_list.get(k);
				String termText = term_list.get(k).text() ;
				double df = term_doc_freq_map.get(termText) ;

				double IDF = (log2((numDocs-df+0.5) / (df+0.5)));
				if (maxIDF_query<IDF) {
					maxIDF_query = IDF;
				}
				DocsEnum docsEnum = atomicReader.termDocsEnum(term);	            // 某个词条的按文档的倒排索引 例如  word-》1-》3-》4... word代表查询的词条，1,3,4代表出现该词的文档
				while ((docsEnum.nextDoc()) != DocIdSetIterator.NO_MORE_DOCS) {
					int doc_id = docsEnum.docID();                                  // 某个词条的在某个文档中出现
					doc_set.add(doc_id);                                            // 将这个出现该词的文档的编号存放在doc_set中
					freq_array[k][doc_id] = docsEnum.freq();                        // 某个词条的在某个文档中出现的次数tf
				}
			}

			/**进行第一次查询 ，得到根据文档得分排好顺序的文档数组，调用函数1.2.1.3.1**/
			ScoreDoc[] compact_score_array = retrievalByDSP(doc_set,term_list,freq_array, maxIDF_query );  ///////////////////***************////////////////////

			//System.out.println(l);//////////////////////////////**************************************/////////////////////////////////

			/**进行后续的结果统计，调用函数1.2.1.3.5**/ 
			CreateResult(compact_score_array,qualityQueries,i,query,stats);
		}
		/**将总的结果统计写入文件，调用函数1.2.1.3.6**/ 
		WriteResult(stats);           
	}

	/** 代码步骤1.2.1.3.1,进行第一次正常查询**/
	public ScoreDoc[] retrievalByDSP(Set<Integer> doc_set,ArrayList<Term> term_list,int[][] freq_array, double maxIDF_query ){

		int seq_n=term_list.size(),seq_i=0,seq_j=0;
		int seq_number= Fun_number(seq_n);
		double [] idf_i=new double[seq_n];
		double [] SSS=new double[seq_number];
		Double seq_s=0.0;
		for(seq_j=0;seq_j<seq_number;seq_j++)
		{
			seq_s=0.0;
			for(seq_i=0;seq_i<seq_n;seq_i++)
			{
				String termText = term_list.get(seq_i).text() ;
				double df = term_doc_freq_map.get(termText);
				double IDF = log2((numDocs - df + 0.5) / (df + 0.5));
				idf_i[seq_i]= IDF;
				//seq_s += idf_i[seq_i]* Math.exp(-(Math.pow(((seq_j-300*(seq_i+1))/Math.sqrt(2*l)),2)));
				
				seq_s += idf_i[seq_i]* Math.exp((Math.pow(((seq_j-300*(seq_i+1))/l),2)));///////。。。。。。。核函数位置。。。。。
			}
			SSS[seq_j]=seq_s;
		}

		int kk = 0;	

		ScoreDoc[] compact_score_array=new ScoreDoc[doc_set.size()];
		for (int j = 0; j < numDocs; j++) { 
			if (doc_set.contains(j)) {  
				double []TEMP = Arrays.copyOf(SSS, SSS.length);

				double total_score = 0.0; 

				int doc_length = doc_length_map.get(j);

				double K =k1* ((1-b) + b*doc_length / avg_doc_length);

				for (int k = 0; k < term_list.size(); k++) {

					double tf = freq_array[k][j];// BM25
					double TF = (k1 + 1) * tf / (K + tf);// BM25
					int df = term_doc_freq_map.get(term_list.get(k).text());// BM25
					double IDF = log2((numDocs - df + 0.5) / (df + 0.5));// BM25
					
					
					//double tf = freq_array[k][j];//LOG BM25
					//double TF = log2((1+(k1 + 1) * tf / (K + tf)));//LOG BM25
					//int df = term_doc_freq_map.get(term_list.get(k).text());//LOG BM25
					//double IDF = log2(1+log2((numDocs - df + 0.5) / (df + 0.5)));// LOG BM25  NOT VERY GOOD
					
					
					
					//使用方法
					//global_maxTermTF         所有文档所有词的最大TF
					//global_maxTermIDF        所有文档所有词的最大IDF
					//doc_maxTermTF_map.get(j)     文档j的所有词的最大TF
					//doc_maxTermIDF_map.get(j)    文档j的所有词的最大IDF
					//maxIDF_query   某个查询表达内的若干词内的最大IDF
					//avg_TF                   所有词的平均TF
					//avg_IDF                  所有词的平均IDF
					//System.out.println("global_maxTermTF  "+global_maxTermTF+"   global_maxTermIDF  "+global_maxTermIDF);
					//System.out.println("doc_maxTermTF_map.get(j)  "+doc_maxTermTF_map.get(j)+"   doc_maxTermIDF_map.get(j)  "+doc_maxTermIDF_map.get(j));
					
					if (TF*IDF/maxIDF_query>0) {
						int ZL=300*(k+1);
						//int ZL=300*(k+1);
						int ZR=ZL+1;
						
                        //System.out.println(s);//////////////////////////////////////////////////
						
                        //int amplitude = (int) Math.round(s*TF*IDF/maxIDF);//MAIN 0.2765__
						
						//int amplitude = (int) Math.round(s*TF*log2(1+IDF)/log2(1+doc_maxTermIDF_map.get(j)));//MAIN LOG 0.2771__						
						
						//int amplitude = (int) Math.round(s*log2(1+TF)*log2(1+IDF)/log2(1+global_maxTermIDF));//0.2764

						//int amplitude = (int) Math.round(s*log2(1+TF)*log2(1+IDF)/log2(1+maxIDF_query));//IMPORTANT-2:0.2790
						
						//double amplitude = s*log2(1+TF)*log2(1+IDF)/log2(1+maxIDF);

                        //int amplitude = (int) Math.round(s*log2(1+TF)/log2(1+global_maxTermTF)*log2(1+IDF)/log2(1+maxIDF_query));//0.2759
						
						//int amplitude = (int) Math.round(s*log2(1+TF)/log2(1+doc_maxTermTF_map.get(j))*log2(1+IDF)/log2(1+maxIDF_query));//0.2712
						
						//int amplitude = (int) Math.round(s*log2(1+TF)*log2(1+IDF)/log2(1+maxIDF_query));
						
						//System.out.println(s*log2(1+TF)*log2(1+IDF));
						
						int amplitude = (int) Math.round(s*log2(1+TF)*log2(1+IDF));//IMPORTANT-1
						
                        //double amplitude = s*log2(1+TF)*log2(1+IDF);//0.2763 //THE PERFORMANCE IS BAD WHEN CHANGING INT FOR DOUBLE
						
						//int amplitude = (int) Math.round(s  *  log2(1+TF)/log2(1+global_maxTermTF) * log2(1+IDF)/log2(1+global_maxTermIDF));//0.2759 
						
						//int amplitude = (int) Math.round(s  *  log2(1+TF)/log2(1+doc_maxTermTF_map.get(j)) * log2(1+IDF)/log2(1+doc_maxTermIDF_map.get(j)));//NOT VERY GOOD
						
						// amplitude = (int) Math.round(s  *  log2(1+TF)/log2(1+doc_maxTermTF_map.get(j)) * log2(1+IDF)/log2(1+global_maxTermIDF));//0.2722
						
						//System.out.println(s);
						
						//System.out.println(24*TF*IDF/maxIDF);//////////////////////////////////////////////////////////////////////////
					
						
						//int amplitude = (int) Math.round(24*TF*IDF*log2(1+doc_length/avg_doc_length)/maxIDF);
						//int amplitude = (int) Math.round(24*TF*IDF*(avg_doc_length/doc_length)/maxIDF);	
						
						
						/*if (s*log2(1+TF)*log2(1+IDF)>132){   //continue 
							continue;
							}
							else {
								int amplitude = (int) Math.round(s*log2(1+TF)*log2(1+IDF));
								for (int r = 0; r <=amplitude; r++) {
									TEMP[ZL-r] = TEMP[ZL-r] *r/amplitude;
									TEMP[ZR+r] = TEMP[ZR+r] *r/amplitude;
								}*/
						
						
						/*if (s*log2(1+TF)*log2(1+IDF)<126){   //else
						int amplitude = (int) Math.round(s*log2(1+TF)*log2(1+IDF));
						for (int r = 0; r <=amplitude; r++) {
							TEMP[ZL-r] = TEMP[ZL-r] *r/amplitude;
							TEMP[ZR+r] = TEMP[ZR+r] *r/amplitude;
						}
						}
						else {
							int amplitude =125;
							for (int r = 0; r <=amplitude; r++) {
								TEMP[ZL-r] = TEMP[ZL-r] *r/amplitude;
								TEMP[ZR+r] = TEMP[ZR+r] *r/amplitude;
						}*/
						
						//int amplitude = (int) Math.round(s*TF*IDF/maxIDF_query);
						
						
						
						
						for (int r = 0; r <=amplitude; r++) {
							TEMP[ZL-r] = TEMP[ZL-r] *r/amplitude;
							TEMP[ZR+r] = TEMP[ZR+r] *r/amplitude;
						}
					}
				}

				for(seq_j=0;seq_j<seq_number;seq_j++)
				{
					total_score +=TEMP[seq_j];
				}

				compact_score_array[kk++] = new ScoreDoc(j, (float)(-total_score));
			}
		}  	
		Arrays.sort(compact_score_array, new ByWeightComparator());             // 按重写的排序方法进行排序（从大到小的顺序-Float.compare(a, b);），比较两个ScoreDoc的score值（float型）   
		return compact_score_array;	
		
	}

	/** 代码步骤1.2.1.3.5,进行后续的结果统计**/
	public void CreateResult(ScoreDoc[] compact_score_array1,QualityQuery[] qualityQueries,int i,Query query,QualityStats[] stats) throws IOException{
		int max_result = Math.min(maxResults, compact_score_array1.length);
		ScoreDoc[] score_doc = new ScoreDoc[max_result];
		System.arraycopy(compact_score_array1, 0, score_doc, 0, max_result);	    // 硬拷贝要求大小的相关得分，// 下一句：构建TopDocs，它有3要素 （总数、得分数组地址、最大得分）		
		TopDocs td = new TopDocs(max_result, score_doc, (float) score_doc[0].score);  
		stats[i] = analyzeQueryResults(qualityQueries[i], query, td, judge, qualityLog, 1);  
	}

	/** 代码步骤1.2.1.3.6,写入各个统计结果**/
	public void WriteResult(QualityStats[] stats){
		QualityStats avg = QualityStats.average(stats);		                  
		avg.log("SUMMARY", 2, qualityLog, "  ");
	}

	/** 代码步骤1.2.1.3.5.1，统计结果中的各项指标的值**/
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

	public static int Fun_number(int n)
	{
		return (int) Math.pow(2, (int) log2(n*300)+1);//Index-array-bound-exception-just enlarge the number here, not everywhere
}
}