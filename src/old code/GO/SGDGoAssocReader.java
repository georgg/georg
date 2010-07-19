package edu.mit.csail.psrg.georg.GO;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

public class SGDGoAssocReader extends OntAssociationReader {

	public SGDGoAssocReader() {
		setFile("C:\\research_active\\go_data\\gene_association.sgd"); 
	}
	
	public OntGeneAssociations readFile(OntDAG DAG, HashSet<String> useGenes) throws IOException {
		BufferedReader is = new BufferedReader(new FileReader(fName));
		String line = "";
		String s = "";
		line = is.readLine();
		int i = 0;
		
		String[] cols = null;
		String[] cols2 = null;
		
		while(line.indexOf("!") == 0) {
			line = is.readLine();
		}
		
		OntGene gene = null;
		OntGeneAssociations assoc = new OntGeneAssociations();
		OntTerm term = null;
		int numGenes = 0;
		boolean useGene = false;
		String geneID = "";
		String qualifier = "";
		while(line != null) {
			cols = line.split("\t");
			geneID = trimGeneID(cols[10]);
			
			// see if the gene is in the filter list
			useGene = true;
			if (useGenes != null) {
				useGene = useGenes.contains(geneID);
			}
			
			if (useGene) {
				gene = assoc.getGene(geneID);	
				if (gene == null) {
					gene = new OntGene();
					gene.ID = geneID;
					gene.Description = cols[9];
					gene.Name = cols[2];
					assoc.addGene(gene);
					numGenes++;
				}
				
				qualifier = cols[3];
				if (qualifier == null) {
					qualifier = "none";
				}
				
				if (!qualifier.contentEquals("NOT")) {
					term = DAG.getTerm(cols[4]);
					if (term != null) {
						term.addOntTerm(gene);
					}
				}
			}
			
			line = is.readLine();
		}
		
		is.close();
		DAG.propagateGeneCounts();
		
		return assoc;
	}

	public String trimGeneID(String ID) {
		int fi = ID.indexOf("|");
		if (fi > 0) {
			return ID.substring(0,fi);
		}
		
		return ID;
	}
}
