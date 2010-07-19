package edu.mit.csail.psrg.georg.GO;

import java.util.ArrayList;
import java.util.HashSet;

public class OntGene {
	public String ID;
	public String Name;
	public String Description;
	public String REFSEQ;
	HashSet<OntTerm> terms = null;
	
	public HashSet<OntTerm> getTerms() {
		return terms;
	}
}
