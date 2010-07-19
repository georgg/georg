package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.GO.DAGException;
import edu.mit.csail.psrg.georg.GO.MIPSAssocReader;
import edu.mit.csail.psrg.georg.GO.MIPSDAG;
import edu.mit.csail.psrg.georg.GO.MIPSDAGReader;
import edu.mit.csail.psrg.georg.GO.OntClusters;
import edu.mit.csail.psrg.georg.GO.OntGeneAssociations;
import edu.mit.csail.psrg.georg.GO.OntTerm;

public class OutputYeastTopics {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String dirPath = "C:\\research_data\\combine_yeast_bind_expression\\";
		String topicInFile = dirPath + "yeast_combine_topics.txt";
		String topicOutFileL = dirPath + "yeast_combine_mips_long.txt";
		String topicOutFileS = dirPath + "yeast_combine_mips_short.txt";
		String normDPsInFile = dirPath + "yeast_combine_dps.txt";
		String normDPsOutFile = dirPath + "yeast_combine_dps_norm.txt";
		
		String snapFileName = dirPath + "yeast_combine_snap_19999.persist";
		HierDirichletProcess myDP = HierDirichletProcess.restoreFromFile(snapFileName);
		try {
			myDP.outputClustersToFile(dirPath + "yeast_combine");
		} catch(IOException e) {
			System.out.println(e);
		}
		
		MicroArrayData myData = new MicroArrayData();
		try {
			myData.setDiscrete();
			myData.readFile(topicInFile);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int minGenes = 5;
		int maxGenes = 500;
		double FDR = 0.05;
		int namespace = OntTerm.AnyCategory;
			
		MIPSDAGReader MIPSReader = new MIPSDAGReader();
		MIPSAssocReader assocReader = new MIPSAssocReader();
		
		TopicSummarizer topicSumm = new TopicSummarizer(myData,MIPSReader,assocReader);
		topicSumm.maxGenes = minGenes;
		topicSumm.maxGenes = maxGenes;
		topicSumm.FDR = FDR;
		topicSumm.namespace = namespace;
		topicSumm.outputTopics(topicOutFileL,topicOutFileS);
		
		double norm = 0.0;
		int i = 0;
		int j = 0;
		myData = new MicroArrayData();
		try {
			myData.readFile(normDPsInFile);
			
			for (i=0;i<myData.numRows;i++) {
				norm = 0.0;
				for (j=0;j<myData.numCols;j++) {
					norm = norm + myData.values[i][j];
				}
				for (j=0;j<myData.numCols;j++) {
					myData.values[i][j] = myData.values[i][j]/norm;
				}
			}
			
			myData.writeFile(normDPsOutFile);
		} catch(IOException e) {
			System.out.println(e);
		}

	}

}
