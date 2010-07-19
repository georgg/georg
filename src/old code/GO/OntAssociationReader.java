package edu.mit.csail.psrg.georg.GO;

import java.io.IOException;
import java.util.HashSet;

public abstract class OntAssociationReader {
	String fName = "";
	
	public void setFile(String f) {
		fName = f;
	}
	
	public abstract OntGeneAssociations readFile(OntDAG DAG,HashSet<String> useGenes) throws IOException;
}
