package edu.mit.csail.psrg.georg.Graph;

import java.io.FileWriter;
import java.io.IOException;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

public class OutputYeastTFGraph {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		outputYFilesGML();
	}

	public static void outputYFilesGML() {
		String dirPath = "C:\\research_data\\combine_yeast_bind_expression\\tfs_only\\";
		String humanTopicOutFileS = dirPath + "yeast_combine_mips_short.txt";
		String mouseTopicOutFileS = dirPath + "yeast_combine_mips_long.txt";
		String normDPsFile = dirPath + "yeast_combine_dps_norm.txt";
		String graphName = dirPath + "yeast_combine_graph.gml";
		
		MicroArrayData myData = new MicroArrayData();
		
		try {
			myData.readFile(normDPsFile);
		} catch (IOException e) {
			System.out.println(e);
		}
		
		String nodeLabel = "";
		String nodeIDS = "";
		String node2IDS = "";
		int i = 0;
		int j = 0;
		int k = 0;
		int f = 0;
		
		// minimal percentage use of a topic
		double loadThreshold = 0.10;
		
		String[] conditionNames = {"YPD","SM","RAPA","H2O2Lo","Acid","Alpha","BUT","HEAT"};
		int[] conditionMap = new int[myData.numRows];
		
		for (i=0;i<myData.numRows;i++) {
			nodeLabel = myData.geneNames[i];
			j = 0;
			while (j < conditionNames.length) {
				f = nodeLabel.indexOf(conditionNames[j]);
				if (f >= 0) {
					break;
				}
				j++;
			}
			if (j == conditionNames.length) {
				conditionMap[i] = -1;
			} else {
				conditionMap[i] = j;
			}
		}
		
		int[] topicCode = new int[myData.numCols];
		int[] topicMap = new int[myData.numCols];
		int numTopics = 0;
		double[] maxLoad = new double[conditionNames.length];
		double mx = 0.0;
		int maxCondition = 0;
		int numMaxCondition = 0;
		int gid = 0;
		
		for (j=0;j<myData.numCols;j++) {
			for (k=0;k<maxLoad.length;k++) {
				maxLoad[k] = 0.0;
			}
			for (i=0;i<myData.numRows;i++) {
				k = conditionMap[i];
				if (k > -1) {
					if (myData.values[i][j] > maxLoad[k] & myData.values[i][j] >= loadThreshold) {
						maxLoad[k] = myData.values[i][j];
					}
				}
			}
			mx  = 0.0;
			numMaxCondition = 0;
			for (k=0;k<maxLoad.length;k++) {
				if (maxLoad[k] > 0.0 & maxLoad[k] >= loadThreshold) {
					if (maxLoad[k] > mx) {
						maxCondition = k;
						mx = maxLoad[k];
					}
					numMaxCondition++;
				}
			}
			if (numMaxCondition > 1) {
				topicCode[j] = conditionNames.length;
			//	topicCode[j] = maxCondition;
			} else {
				if (numMaxCondition == 1) {
					topicCode[j] = maxCondition;
				}
				if (numMaxCondition == 0) {
					topicCode[j] = -1;
				}
			}
			
			if (topicCode[j] > -1) {
				numTopics++;
				topicMap[j] = numTopics;
			}
		}
		
		try {
			FileWriter outFile = new FileWriter(graphName);
			outFile.write("Creator\t\"yFiles\"\nVersion	2.2\n");
			outFile.write("graph\n[\n");
			outFile.write("directed 1\n");
			outFile.write("hierarchic 1\n");
			
		/*	for (i=0;i<conditionNames.length+1;i++) {
				if (i < conditionNames.length) {
					gid = i + numTopics + myData.numRows;
				} else {
					gid = -1;
				}
				outputTopics(outFile,topicCode,topicMap,i,gid);
			} */
			outputTopics(outFile,topicCode,topicMap,conditionNames.length,-1);
			
			for (i=0;i<conditionNames.length;i++) {
				outputTFs(outFile,myData,loadThreshold,topicMap,topicCode,conditionNames.length,numTopics,conditionNames[i],i + numTopics + myData.numRows);
			}
			
			outFile.write("]\n");
			outFile.close();
			
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void outputTopics(FileWriter outFile,int[] topicCode,int[] topicMap,int selectCode,int gid) throws IOException {
		// output shared topics
		int i = 0;
		String nodeLabel = "";
		for (i=0;i<topicCode.length;i++) {
			if (topicCode[i] == selectCode) {
				outFile.write("node\n[\n");
				outFile.write("id " + (new Integer(topicMap[i])).toString() + "\n");
				nodeLabel = (new Integer(i+1)).toString();
				outFile.write("label " + "\"" + nodeLabel + "\"\n");
				outFile.write("graphics\n[\n");
				outFile.write("type \"circle\"\n");
				outFile.write("]\n");
				if (gid > -1) {
					outFile.write("gid " + (new Integer(gid)).toString() + "\n");
				}
				outFile.write("]\n");
			}
		}
	}
	
	public static void outputTFs(FileWriter outFile,MicroArrayData myData,double loadThreshold,int[] topicMap,int[] topicCode,int maxTopicCode,int numTopics,String conditionName,int gid) throws IOException {
		int i = 0;
		int j = 0;
		int f = 0;
		String nodeLabel = "";
		String nodeIDS = "";
		String node2IDS = "";
		
		outFile.write("node\n[\n");
		outFile.write("id " + (new Integer(gid)).toString() + "\n");
		outFile.write("label " + "\"" + conditionName + "\"\n");
		outFile.write("isGroup 1\n");
		outFile.write("]\n");
		
		for (i=0;i<myData.numRows;i++) {
			f = myData.geneNames[i].indexOf("_" + conditionName);
			if (f >= 0) {
				nodeLabel = myData.geneNames[i].substring(0,f);
				nodeIDS = (new Integer(i+numTopics)).toString();
				outFile.write("node\n[\n");
				outFile.write("id " + nodeIDS + "\n");
				outFile.write("label " + "\"" + nodeLabel + "\"\n");
				outFile.write("gid " + (new Integer(gid)).toString() + "\n");
				outFile.write("]\n");
				for (j=0;j<myData.numCols;j++) {
				//	if (myData.values[i][j] >= loadThreshold) {
					if (myData.values[i][j] >= loadThreshold & topicCode[j] == maxTopicCode) {
						node2IDS = (new Integer(topicMap[j])).toString();
						outFile.write("edge\n[\n");
						outFile.write("source " + nodeIDS + "\n");
						outFile.write("target " + node2IDS + "\n");
						outFile.write("]\n");
					}
				}
			}
		}
	}
}
