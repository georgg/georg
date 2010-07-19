package edu.mit.csail.psrg.georg.GO;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class MIPSDAGReader extends DAGReader {
	String schemeName = "C:\\research_data\\mips_data\\funcat-2.0_SCHEME";
	
	public MIPSDAG readFile() throws IOException, DAGException {
		BufferedReader is = new BufferedReader(new FileReader(schemeName));
		String line = "";
		String desc = "";
		String ID = "";
		int fi = 0;
		MIPSDAG DAG = new MIPSDAG();
		OntTerm term = null;
		String pID = "";
		
		line = is.readLine();
		
		while(line != null) {
			fi = line.indexOf(" ");
			ID = line.substring(0,fi);
			desc = line.substring(fi+1);
			pID = getTermParent(ID);
			
			term = new OntTerm();
			term.ID = ID;
			term.name = desc;
			DAG.addTerm(term);
			
			if (pID != null) {
				term.addParentName(pID);
			}
			line = is.readLine();
		}
		
		DAG.associateParents();
		DAG.topologicalSort();
		return DAG;
	}
	
	public String getTermParent(String termID) {
		String[] parents = null;
		int fi = termID.lastIndexOf(".");
		if (fi < 0) {
			return null;
		}
		return termID.substring(0,fi);
	}
}
