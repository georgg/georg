package edu.mit.csail.psrg.georg.GO;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.StringTokenizer;

public class AffyGoAssocReader extends OntAssociationReader {
	boolean useREFSEQ = false;
	
	public AffyGoAssocReader() {
		setFile("C:\\research_data\\go_data\\GPL339.txt");
	}
	
	public void setUseREFSEQ(boolean u) {
		useREFSEQ = u;
	}
	
	public OntGeneAssociations readFile(OntDAG DAG,HashSet<String> useGenes) throws IOException {
		BufferedReader is = new BufferedReader(new FileReader(fName));
		String line = "";
		String s = "";
		line = is.readLine();
		int i = 0;
		
		String[] colNames = new String[15];
		String[] cols = null;
		String[] cols2 = null;
		
		while(!line.contentEquals("!platform_table_begin")) {
			line = is.readLine();
		}
		
		line = is.readLine();
		colNames = line.split("\t");
		
		line = is.readLine();
		OntGene gene = null;
		OntGeneAssociations assoc = new OntGeneAssociations();
		ArrayList<String> GoIDs = null;
		int numGenes = 0;
		boolean useGene = false;
		while(!line.contentEquals("!platform_table_end")) {
			cols = line.split("\t");
			cols2 = new String[15];
			for (i=0;i<cols.length;i++) {
				cols2[i] = "";
			}
			for (i=0;i<cols.length;i++) {
				cols2[i] = cols[i];
			}
			cols = cols2;
			
			gene = new OntGene();
			if (!useREFSEQ) {
				gene.ID = cols[0];
			} else {
				if (cols[3] != null) {
					gene.ID = cols[3];
				} else {
					gene.ID = cols[0];
				}
			}
			
			// see if the gene is in the filter list
			useGene = true;
			if (useGenes != null) {
				useGene = useGenes.contains(gene.ID);
			}
			
			if (!assoc.checkGene(gene) && useGene) {		
				gene.Description = cols[7];
				gene.Name = cols[8];
				gene.REFSEQ = cols[3];
			
				GoIDs = parseGoLine(cols[12]);
				addGoTerms(gene,GoIDs,DAG);
				
				GoIDs = parseGoLine(cols[13]);
				addGoTerms(gene,GoIDs,DAG);
				
				GoIDs = parseGoLine(cols[14]);
				addGoTerms(gene,GoIDs,DAG);
							
				assoc.addGene(gene);
				numGenes++;
			}
			
			line = is.readLine();
		}
		
		is.close();
		
		DAG.propagateGeneCounts();
		return assoc;
	}
	
	ArrayList<String> parseGoLine(String line) {
		ArrayList<String> goIDs = new ArrayList<String>();
		String s = "";
		
		if (line == null) {
			return goIDs;
		}
		
		int fi = line.indexOf(" ");
		int fi2 = 0;
		if (fi >= 0) {
			s = line.substring(0,fi);
			goIDs.add(padGOID(s));
			fi = line.indexOf("///",fi);
			while(fi >= 0) {
				fi2 = line.indexOf(" ",fi+4);
				if (fi2 >= 0) {
					s = line.substring(fi+4,fi2);
					goIDs.add(padGOID(s));
					fi = fi2+1;
				}
				fi = line.indexOf("///",fi);
			}
		}
		
		return goIDs;
	}
	
	String padGOID(String ID) {
		String s = ID;
		while(s.length() < 7) {
			s = "0" + s;
		}
		s = "GO:" + s;
		return s;
	}
	
	void addGoTerms(OntGene gene,ArrayList<String> terms,OntDAG DAG) {
		int i = 0;
		OntTerm parent = null;
		for (i=0;i<terms.size();i++) {
			parent = DAG.getTerm(terms.get(i));
			if (parent != null) {
				parent.addOntTerm(gene);
			}
		}
	}
	
	public LinkedHashMap<String,ArrayList<Integer>> uniqueGeneMap(OntGeneAssociations assoc,String[] geneNames) {
		LinkedHashMap<String,ArrayList<Integer>> geneCounts = new LinkedHashMap<String,ArrayList<Integer>>();
		
		int i = 0;
		ArrayList<Integer> geneList = null;
		OntGene gene = null;
		for (i=0;i<geneNames.length;i++) {
			gene = assoc.getGene(geneNames[i]);
			if (!geneCounts.containsKey(gene.REFSEQ)) {
				geneList = new ArrayList<Integer>();
				geneList.add(i);
				geneCounts.put(gene.REFSEQ,geneList);
			} else {
				geneList = geneCounts.get(gene.REFSEQ);
				geneList.add(i);
			}
		}
		
		return geneCounts;
	}
}
