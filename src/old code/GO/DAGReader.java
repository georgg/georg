package edu.mit.csail.psrg.georg.GO;

import java.io.IOException;

public abstract class DAGReader {
	public abstract OntDAG readFile() throws IOException, DAGException;
}
