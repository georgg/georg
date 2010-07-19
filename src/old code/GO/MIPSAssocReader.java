package edu.mit.csail.psrg.georg.GO;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;

public class MIPSAssocReader extends OntAssociationReader {
	public MIPSAssocReader() {
		setFile("C:\\research_data\\mips_data\\funcat-2.0_data_20062005");
	}
	
	public OntGeneAssociations readFile(OntDAG DAG, HashSet<String> useGenes) throws IOException {
		OntGeneAssociations assoc = new OntGeneAssociations();
		
		BufferedReader is = new BufferedReader(new FileReader(fName));
		String line = "";
		String geneID = "";
		String catID = "";
		int fi = 0;
		int fi2 = 0;
		OntGene gene = null;
		OntTerm cat = null;
		boolean useGene = true;
		
		line = is.readLine();
		
		while(line != null) {
			fi = line.indexOf("|");
			geneID = line.substring(0,fi);
			geneID = geneID.toUpperCase();
			catID = line.substring(fi+1);
			fi2 = catID.indexOf("|");
			catID = catID.substring(0,fi2);
			
			useGene = true;
			if (useGenes != null) {
				useGene = useGenes.contains(geneID);
			}
			
			if (useGene) {
				gene = assoc.getGene(geneID);
				if (gene == null) {
					gene = new OntGene();
					gene.ID = geneID;
					assoc.addGene(gene);
				}
				cat = DAG.getTerm(catID);
				cat.addOntTerm(gene);
			}
			
			line = is.readLine();
		}
		
		DAG.propagateGeneCounts();
		return assoc;
	}

}
