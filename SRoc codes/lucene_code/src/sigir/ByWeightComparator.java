package sigir;

import java.util.Comparator;

import org.apache.lucene.search.ScoreDoc;

class ByWeightComparator implements Comparator<Object> {
	public final int compare(Object pFirst, Object pSecond) {
		float a = ((ScoreDoc) pFirst).score;
		float b = ((ScoreDoc) pSecond).score;
		return -Float.compare(a, b);
	}
}

