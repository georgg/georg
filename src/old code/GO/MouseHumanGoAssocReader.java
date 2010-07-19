package edu.mit.csail.psrg.georg.GO;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class MouseHumanGoAssocReader extends OntAssociationReader {
	public OntGeneAssociations readFile(OntDAG DAG, HashSet<String> useGenes) throws IOException {
		
		BufferedReader is = new BufferedReader(new FileReader(fName));
		String line = "";
		String s = "";
		int i = 0;
		int f = 0;
		
		OntGene gene = null;
		OntGeneAssociations assoc = new OntGeneAssociations();
		ArrayList<String> GoIDs = null;
		int numGenes = 0;
		boolean useGene = false;
		String geneName;
		String cats[] = null;
		line = is.readLine();
		String cols[] = null;
		while(line != null) {
			cols = line.split("\t");
			
			gene = new OntGene();
			gene.ID = cols[0];
			
			if (cols.length == 4) {
				gene.Name = cols[2];
				gene.Description = cols[3];
			} else {
				if (cols.length > 1) {
					gene.Name = cols[1];
				}
				if (cols.length > 2) {
					gene.Description = cols[2];
				}
			}
			
			// see if the gene is in the filter list
			useGene = true;
			if (useGenes != null) {
				useGene = useGenes.contains(gene.ID);
			}
			
			if (!assoc.checkGene(gene) && useGene) {
				assoc.addGene(gene);
				numGenes++;
			}
			
			if (cols.length == 4) {
				cats = cols[1].split(";");
			
				if (useGene) {
					addGoTerms(gene,cats,DAG);
				}
			}
			
			line = is.readLine();
		}
		
		is.close();
		
		DAG.propagateGeneCounts();
		return assoc;
	}
	
	void addGoTerms(OntGene gene,String[] terms,OntDAG DAG) {
		int i = 0;
		OntTerm parent = null;
		for (i=0;i<terms.length;i++) {
			parent = DAG.getTerm(terms[i]);
			if (parent != null) {
				parent.addOntTerm(gene);
			}
		}
	}

}
