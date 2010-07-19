package edu.mit.csail.psrg.georg.Annotation;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;

public class OntGene {
	public String ID;
	public String Name;
	public String Description;
	
	HashSet<OntTerm> terms = null;
	
	public HashSet<OntTerm> getTerms() {
		return terms;
	}
	
	public String getTermString(String namespace) {
		String sterms = "";
		OntTerm term = null;
		
		if (terms == null) {
			return null;
		}
		
		Iterator<OntTerm> termIter = terms.iterator();
		
		while(termIter.hasNext()) {
			term = termIter.next();
			if (term.namespace.equals(namespace)) {
				sterms = term.name + " ";
			}
		}
		
		return sterms;
	}
}
