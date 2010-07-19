package edu.mit.csail.psrg.georg.HDPTopicModel;

public class GeneHolder implements Comparable {
	String[] names = null;
	double rank = 0.0;
	double[] counts = null;
	
	public GeneHolder(String[] nms,double r,double[] c) {
		int i = 0;
		names = nms;
		rank = r;
		counts = c;
	}
	
	public int compareTo(Object o) {
		GeneHolder holder = (GeneHolder) o;
		if (rank == holder.rank)
			return 0;
		if (rank < holder.rank)
			return 1;
		return -1;
	}
}

