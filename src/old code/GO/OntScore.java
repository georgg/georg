package edu.mit.csail.psrg.georg.GO;

import java.util.HashSet;

public class OntScore implements Comparable {
	public HashSet<OntGene> genes = new HashSet<OntGene>();
	public double fraction = 0.0;
	public double pval = 0.0;
	public boolean acceptPVal = false;
	int rank = 0;
	
	public void addGene(OntGene gene) {
		genes.add(gene);
	}

	public int compareTo(Object o) {
		OntScore score = (OntScore) o;
		
		if (pval == score.pval)
			return 0;
		
		if (pval < score.pval)
			return -1;
		
		return 1;
	}
}
