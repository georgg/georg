package edu.mit.csail.psrg.georg.HDP2;

public class OverlapTissue implements Comparable {
	// the percent load of the controlling tissue (the percentage of the tissue's overlapping genes
	// out of the total genes)
	public double tissueUse = 0.0;
	// the percentage of the tissue's overlaping genes out of the topic's genes
	public double topicUse = 0.0;
	// number of times the tissue uses the topic for upregulated genes
	public int numUp = 0;
	// number of times the tissue significantly uses the topic
	public int numSignificant = 0;
	String tissueName;
	int code = 0;
	boolean sortTissueUse = false;
	
	public OverlapTissue(String tn,double tiu,double tou,boolean up,boolean isSig) {
		tissueUse = tiu;
		topicUse = tou;
		if (up) {
			numUp = 1;
		}
		if (isSig) {
			numSignificant = 1;
		}
		tissueName = tn;
		
		code = tissueName.hashCode();
	}
	
	public int hashCode() {
		return code;
	}
	
	public void add(OverlapTissue t2) {
		tissueUse += t2.tissueUse;
		topicUse += t2.topicUse;
		numUp += t2.numUp;
		numSignificant += t2.numSignificant;
	}
	
/*	public String getString(int numSamples) {
		String s = tissueName;
		String t = null;
		
		t = String.format("%.3f",tissueUse/((double) numSamples));
		s = s + " (" + t + "/";
		t = String.format("%.3f",topicUse/((double) numSamples));
		s = s + t + "/";
		t = String.format("%.3f",((double) numSignificant)/((double) numSamples));
		s = s + t + ")";
		
		double u = ((double) numUp)/((double) numSamples);
		t = String.format("%.3f",u);
		if (u >= 0.50) {
			s = s + " UP[" + t + "]";
		} else {
			s = s + " DN[" + t + "]";
		}
		
		return s;
	} */
	
	public int compareTo(Object o2) {
		OverlapTissue ct2 = (OverlapTissue) o2;
		
		if (sortTissueUse) {
			if (tissueUse > ct2.tissueUse)
				return -1;
			
			if (tissueUse == ct2.tissueUse)
				return 0;
	
			return 1;
		} else {
			if (topicUse > ct2.topicUse)
				return -1;
			
			if (topicUse == ct2.topicUse)
				return 0;
	
			return 1;
		}
	}
	
	public boolean equals(Object o) {
		if (((ControllingTissue) o).code == code) {
			return true;
		}
		
		return false;
	}
}
