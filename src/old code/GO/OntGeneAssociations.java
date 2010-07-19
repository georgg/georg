package edu.mit.csail.psrg.georg.GO;

import java.util.HashMap;

public class OntGeneAssociations {
	HashMap<String,OntGene> genes = new HashMap<String,OntGene>();
	
	public void addGene(OntGene gene) {
		genes.put(gene.ID,gene);
	}
	
	public boolean checkGene(OntGene gene) {
		return genes.containsKey(gene.ID);
	}
	
	public OntGene getGene(String name) {
		return genes.get(name);
	}
}
