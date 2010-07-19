package edu.mit.csail.psrg.georg.GO;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import edu.mit.csail.psrg.georg.StatUtil.*;

public class TestGO {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	//	testSGDReader();
	//	testMIPS();
		testMouseHumanReader();
	}
	
	public static void testMIPS() {
		MIPSDAGReader reader = new MIPSDAGReader();
		MIPSAssocReader assocReader = new MIPSAssocReader();
		OntGeneAssociations assoc = null;
		MIPSDAG DAG = null;
		
		try {
			DAG = reader.readFile();
			System.out.println("Loaded MIPS terms");
			assoc = assocReader.readFile(DAG,null);
			System.out.println("Loaded MIPS associations");
		} catch(IOException e) {
			System.out.println(e);
		} catch(DAGException e) {
			System.out.println(e);
		}
	}
	
	public static void testSGDReader() {
		GoOBOReader reader = new GoOBOReader();
		SGDGoAssocReader SGDReader = new SGDGoAssocReader();
		OntDAG DAG = null;
		OntGeneAssociations assoc = null;
		
		try {
			DAG = reader.readFile();
			System.out.println("Loaded GO terms");
			assoc = SGDReader.readFile(DAG,null);
			System.out.println("Loaded gene associations");
		} catch(IOException e) {
			System.out.println(e);
		} catch(DAGException e) {
			System.out.println(e);
		}
	}
	
	public static void testAffyReader() {
		GoOBOReader reader = new GoOBOReader();
		AffyGoAssocReader affyReader = new AffyGoAssocReader();
		OntDAG DAG = null;
		OntGeneAssociations assoc = null;
		
		OntTermScores scores = new OntTermScores();
		ArrayList<String> geneNames = new ArrayList<String>();
		geneNames.add("1415675_at");
		geneNames.add("1415676_a_at");
		geneNames.add("1415677_at");
		geneNames.add("1415678_at");
		geneNames.add("1415679_at");
		
		try {
			DAG = reader.readFile();
			System.out.println("Loaded GO terms");
			assoc = affyReader.readFile(DAG,null);
			System.out.println("Loaded gene associations");
			scores.scoreGeneSet(geneNames,null,assoc,5,200,0.05,OntTerm.AnyCategory);
			System.out.println("Scored genes");
			FileWriter outFile = new FileWriter("C:\\CPPDataFiles\\test_go.txt");
			scores.outputScoresLong(outFile,true);
			outFile.close();
		} catch(IOException e) {
			System.out.println(e);
		} catch(DAGException e) {
			System.out.println(e);
		}

	}
	
	public static void testMouseHumanReader() {
		GoOBOReader reader = new GoOBOReader();
		MouseHumanGoAssocReader assocReader = new MouseHumanGoAssocReader();
		assocReader.setFile("c:\\research_data\\mouse_human\\Hs.genelist.go");
		OntDAG DAG = null;
		OntGeneAssociations assoc = null;
		
		OntTermScores scores = new OntTermScores();
		ArrayList<String> geneNames = new ArrayList<String>();
		geneNames.add("NM_002079");
		geneNames.add("NM_001308");
		
		try {
			DAG = reader.readFile();
			System.out.println("Loaded GO terms");
			assoc = assocReader.readFile(DAG,null);
			System.out.println("Loaded gene associations");
			scores.scoreGeneSet(geneNames,null,assoc,5,200,0.05,OntTerm.AnyCategory);
			System.out.println("Scored genes");
			FileWriter outFile = new FileWriter("C:\\CPPDataFiles\\test_go.txt");
			scores.outputScoresLong(outFile,true);
			outFile.close();
		} catch(IOException e) {
			System.out.println(e);
		} catch(DAGException e) {
			System.out.println(e);
		}
	}

}
