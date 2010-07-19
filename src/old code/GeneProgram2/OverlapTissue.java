package edu.mit.csail.psrg.georg.GeneProgram2;

public class OverlapTissue implements Comparable {
	// the percent load of the controlling tissue (the percentage of the tissue's overlapping genes
	// out of the total genes)
	public double tissueUse = 0.0;
	// the percentage of the tissue's overlaping genes out of the topic's genes
	public double topicUse = 0.0;
	// number of times the tissue uses the topic for upregulated genes
	public int numUp = 0;
	// number of times the tissue uses a topic with a given modifier value
	public int[][] numModifier = null;
	// number of times the tissue significantly uses the topic
	public int numSignificant = 0;
	String tissueName;
	int code = 0;
	boolean sortTissueUse = false;
	double weight = 0.0;
	boolean passThresh = true;
	
	public OverlapTissue(String tn,double tiu,double tou,boolean up,boolean isSig,int[] modifiers,HDP myHDP) {
		tissueUse = tiu;
		topicUse = tou;
		
		weight = tou;
		
		if (up) {
			numUp = 1;
		}
		if (isSig) {
			numSignificant = 1;
		}
		tissueName = tn;
		
		code = tissueName.hashCode();
		
		if (myHDP.modifierLevels != null) {
			int i = 0;
			
			numModifier = new int[myHDP.modifierLevels.length][];
			for (i=0;i<numModifier.length;i++) {
				numModifier[i] = new int[myHDP.modifierLevels[i]];
			}
			
			for (i=0;i<modifiers.length;i++) {
				numModifier[i][modifiers[i]]++;
			}
		}
	}
	
	public int hashCode() {
		return code;
	}
	
	public void add(OverlapTissue t2) {
		tissueUse += t2.tissueUse;
		topicUse += t2.topicUse;
		numUp += t2.numUp;
		numSignificant += t2.numSignificant;
		
		if (numModifier != null) {
			int i = 0;
			int j = 0;
			for (i=0;i<numModifier.length;i++) {
				for (j=0;j<numModifier[i].length;j++) {
					numModifier[i][j] += t2.numModifier[i][j];
				}
			}
		}
	}
	
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
		if (((OverlapTissue) o).code == code) {
			return true;
		}
		
		return false;
	}
}
