package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.GO.*;

public class OutputTopics {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	//	extractTopics("C:\\CPPDataFiles\\yeast_genedocs_topics.txt","C:\\CPPDataFiles\\yeast_genedocs_topics_summary.txt");
		topicCategories("C:\\CPPDataFiles\\yeast_genedocs_dps.txt","C:\\CPPDataFiles\\yeast_genedocs_dps_mips_short.txt","C:\\CPPDataFiles\\yeast_genedocs_dps_mips.txt");
	}
	
	public static void topicCategories(String fInName,String fOutNameShort,String fOutNameLong) {
//		String fInName = "C:\\CPPDataFiles\\macklis_topics.txt";
		//	String fOutNameShort = "C:\\CPPDataFiles\\macklis_topics_GO_short.txt";
		//	String fOutNameLong = "C:\\CPPDataFiles\\macklis_topics_GO_long.txt";
		//	String fAssocName = "C:\\research_active\\go_data\\GPL339.txt";
		//	String fInName = "C:\\CPPDataFiles\\yeast_topics_test.txt";
		//	String fInName = "C:\\CPPDataFiles\\yeast_genedocs_dps.txt";
		//	String fOutNameShort = "C:\\CPPDataFiles\\yeast_topics_GO_short.txt";
		//	String fOutNameLong = "C:\\CPPDataFiles\\yeast_topics_GO_long.txt";
		//	String fAssocName = "C:\\research_active\\go_data\\GPL339.txt";
			OntGeneAssociations assoc = null;
		//	OntDAG DAG = null;
			MIPSDAG DAG = null;
			
			int minGenes = 5;
			int maxGenes = 500;
			double FDR = 0.05;
			int namespace = OntTerm.AnyCategory;
			
			MicroArrayData myData = new MicroArrayData();
			MIPSDAGReader MIPSReader = new MIPSDAGReader();
			MIPSAssocReader assocReader = new MIPSAssocReader();
//			GoOBOReader goReader = new GoOBOReader();
//			AffyGoAssocReader assocReader = new AffyGoAssocReader();
//			assocReader.setFile(fAssocName);
//			SGDGoAssocReader assocReader = new SGDGoAssocReader();
			try {
				myData.readFile(fInName);
			} catch(IOException e) {
				System.out.println(e);
			}
			
			HashSet<String> useGenes = new HashSet<String>();
			int i = 0;
			int j = 0;
			
			for (i=0;i<myData.numRows;i++) {
				useGenes.add(myData.geneNames[i]);
			}
			try {
			//	DAG = goReader.readFile();
				DAG = MIPSReader.readFile();
				System.out.println("Loaded DAG");
				assoc = assocReader.readFile(DAG,useGenes);
				System.out.println("Loaded gene associations");
			} catch(IOException e) {
				System.out.println(e);
			} catch(DAGException e) {
				System.out.println(e);
			}
			
			ArrayList<String> genes = new ArrayList<String>();
			
			OntClusters clusters = new OntClusters(assoc,minGenes,maxGenes,FDR,namespace);
			for(j=0;j<myData.numCols;j++) {
				genes.clear();
				for (i=0;i<myData.numRows;i++) {
					if (myData.values[i][j] > 0) {
						genes.add(myData.geneNames[i]);
					}
				}
	//			clusters.addCluster(genes);
			}
			clusters.clusterSignif();
			
			try {
				clusters.outputClustersLong(fOutNameLong);
			} catch(IOException e) {
				System.out.println(e);
			}
	}
	
	public static void extractTopics(String fInName,String fOutName) {
		MicroArrayData myData = new MicroArrayData();
		myData.setDiscrete();
		try {
			myData.readFile(fInName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		int j = 0;
		int minStrength = 5;
		int totalStrength = 0;
		double percent = 0.0;
		String s = "";
		
		try {
			FileWriter outFile = new FileWriter(fOutName);
			
			for(j=0;j<myData.numCols;j++) {
				outFile.write(myData.experimentNames[j]);
				outFile.write("\n");
				for (i=0;i<myData.numRows;i++) {
					totalStrength += myData.dvalues[i][j];
				}
				for (i=0;i<myData.numRows;i++) {
					if (myData.dvalues[i][j] >= minStrength) {
						percent = ((double) myData.dvalues[i][j])/((double) totalStrength);
						outFile.write(myData.geneNames[i]);
						outFile.write("\t");
						outFile.write((new Double(percent)).toString());
						outFile.write("\n");
					}
				}
				outFile.write("\n\n");
			}
			
			outFile.close();
		} catch(IOException e) {
			System.out.println(e);
		}
		
	}

}
