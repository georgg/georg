package edu.mit.csail.psrg.georg.Annotation;

import java.io.IOException;

public class Test {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Annotations DAG = new Annotations();
		
		try {
			DAG.readFiles("c:\\research_data\\Jenner\\infection_annotations.txt","C:\\research_data\\go_data\\gene_ontology.obo",null);
			System.out.println("Done");
		} catch(IOException e) {
			System.out.println(e);
		} catch(DAGException e) {
			System.out.println(e);
		}

	}

}
