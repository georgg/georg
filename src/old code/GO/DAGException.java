package edu.mit.csail.psrg.georg.GO;

public class DAGException extends Exception {
	public String toString() {
		return "Graph is not a DAG";
	}
}
