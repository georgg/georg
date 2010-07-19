package edu.mit.csail.psrg.georg.Annotation;

import java.util.HashSet;

// used to score terms in the ontology
public class OntTermScore implements Comparable {
	public HashSet<OntGene> genes = new HashSet<OntGene>();
	public double fraction = 0.0;
	public double pval = 0.0;
	public boolean acceptPVal = false;
	int rank = 0;
	
	public void addGene(OntGene gene) {
		genes.add(gene);
	}

	public int compareTo(Object o) {
		OntTermScore score = (OntTermScore) o;
		
		if (pval == score.pval)
			return 0;
		
		if (pval < score.pval)
			return -1;
		
		return 1;
	}
}
