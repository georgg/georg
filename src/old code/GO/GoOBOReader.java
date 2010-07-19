package edu.mit.csail.psrg.georg.GO;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.StringTokenizer;

public class GoOBOReader extends DAGReader {
	public String OBOFileName = "C:\\research_data\\go_data\\gene_ontology.obo";
	
	public OntDAG readFile() throws IOException, DAGException {
		BufferedReader is = new BufferedReader(new FileReader(OBOFileName));
		String line = "";
		OntDAG DAG = new OntDAG();
		line = is.readLine();
		OntTerm term = null;
		String tag = "";
		String content = "";
		int fi = 0;
		int terms = 0;
		
		while(line != null) {
			if (line.equals("[Term]")) {
				if (term != null) {
					DAG.addTerm(term);
				}
				term = new OntTerm();
				terms++;
			} else {
				if (terms > 0) {
					fi = line.indexOf(":");
					if (fi >= 0) {
						tag = line.substring(0,fi);
						content = line.substring(fi+2);
					}
					term.setField(tag,content);
				}
			}
			line = is.readLine();
		}
		
		is.close();
		DAG.associateParents();
		DAG.topologicalSort();
		return DAG;
	}
}
