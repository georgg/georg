package edu.mit.csail.psrg.georg.Graph;
import java.io.FileWriter;
import java.io.IOException;

import edu.mit.csail.psrg.georg.DataAccess.*;

public class OutputMouseHumanGraph {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	//	outputGraphViz();
		outputYFilesGML();
	}
	
	public static void outputYFilesGML() {
		String dirPath = "C:\\research_data\\mouse_human\\";
		String humanTopicOutFileS = dirPath + "tissue_as_docs_topics_go_human_short.txt";
		String mouseTopicOutFileS = dirPath + "tissue_as_docs_topics_go_mouse_short.txt";
		String normDPsFile = dirPath + "tissue_as_docs_dps_norm.txt";
		String graphName = dirPath + "tissue_as_docs_graph.gml";
		
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
		int f = 0;
		
		// minimal percentage use of a topic by any given tissue
		double loadThreshold = 0.10;
		
		// 0 = discard the topic (no significant loading)
		// 1 = topic in human only
		// 2 = topic in mouse only
		// 3 = topic shared
		int[] topicCode = new int[myData.numCols];
		int[] topicMap = new int[myData.numCols];
		int numTopics = 0;
		double maxLoadHuman = 0.0;
		double maxLoadMouse = 0.0;
		
		for (j=0;j<myData.numCols;j++) {
			maxLoadHuman = 0.0;
			maxLoadMouse = 0.0;
			for (i=0;i<myData.numRows;i++) {
				f = myData.geneNames[i].indexOf("human_");
				if (f >= 0) {
					if (myData.values[i][j] > maxLoadHuman) {
						maxLoadHuman = myData.values[i][j];
					}
				}
				f = myData.geneNames[i].indexOf("mouse_");
				if (f >= 0) {
					if (myData.values[i][j] > maxLoadMouse) {
						maxLoadMouse = myData.values[i][j];
					}
				}
			}
			if (maxLoadHuman >= loadThreshold & maxLoadMouse >= loadThreshold) {
				topicCode[j] = 3;
			} else {
				if (maxLoadHuman >= loadThreshold) {
					topicCode[j] = 1;
				}
				if (maxLoadMouse >= loadThreshold) {
					topicCode[j] = 2;
				}
			}
			if (topicCode[j] > 0) {
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
			
			int humanGID = numTopics + myData.numRows + 1;
			int mouseGID = humanGID + 1;
			
			// output shared topics
			outputTopics(outFile,topicCode,topicMap,3,-1);
			
			// output human specific topics
			outputTopics(outFile,topicCode,topicMap,1,humanGID);
			// output human tissues
			outputTissues(outFile,myData,loadThreshold,topicMap,numTopics,"human",humanGID);
			
			// output mouse specific topics
			outputTopics(outFile,topicCode,topicMap,2,mouseGID);
			// output mouse tissues
			outputTissues(outFile,myData,loadThreshold,topicMap,numTopics,"mouse",mouseGID);
			
			outFile.write("]\n");
			
			outFile.close();
			
		} catch (IOException e) {
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
	
	public static void outputTissues(FileWriter outFile,MicroArrayData myData,double loadThreshold,int[] topicMap,int numTopics,String speciesName,int gid) throws IOException {
		int i = 0;
		int j = 0;
		int f = 0;
		String nodeLabel = "";
		String nodeIDS = "";
		String node2IDS = "";
		
		outFile.write("node\n[\n");
		outFile.write("id " + (new Integer(gid)).toString() + "\n");
		outFile.write("label " + "\"" + speciesName + "\"\n");
		outFile.write("isGroup 1\n");
		outFile.write("]\n");
		
		for (i=0;i<myData.numRows;i++) {
			f = myData.geneNames[i].indexOf(speciesName+"_");
			if (f >= 0) {
				nodeLabel = myData.geneNames[i].substring(f+6);
				nodeIDS = (new Integer(i+numTopics)).toString();
				outFile.write("node\n[\n");
				outFile.write("id " + nodeIDS + "\n");
				outFile.write("label " + "\"" + nodeLabel + "\"\n");
				outFile.write("gid " + (new Integer(gid)).toString() + "\n");
				outFile.write("]\n");
				for (j=0;j<myData.numCols;j++) {
					if (myData.values[i][j] >= loadThreshold) {
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
	
	public static void outputGraphViz() {
		String dirPath = "C:\\research_data\\mouse_human\\";
		String humanTopicOutFileS = dirPath + "tissue_as_docs_topics_go_human_short.txt";
		String mouseTopicOutFileS = dirPath + "tissue_as_docs_topics_go_mouse_short.txt";
		String normDPsFile = dirPath + "tissue_as_docs_dps_norm.txt";
		String graphName = dirPath + "tissue_as_docs_graph.dot";
		
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
		int f = 0;
		
		// minimal percentage use of a topic by any given tissue
		double loadThreshold = 0.05;
		
		// 0 = discard the topic (no significant loading)
		// 1 = topic in human only
		// 2 = topic in mouse only
		// 3 = topic shared
		int[] topicCode = new int[myData.numCols];
		double maxLoadHuman = 0.0;
		double maxLoadMouse = 0.0;
		
		for (j=0;j<myData.numCols;j++) {
			maxLoadHuman = 0.0;
			maxLoadMouse = 0.0;
			for (i=0;i<myData.numRows;i++) {
				f = myData.geneNames[i].indexOf("human_");
				if (f >= 0) {
					if (myData.values[i][j] > maxLoadHuman) {
						maxLoadHuman = myData.values[i][j];
					}
				}
				f = myData.geneNames[i].indexOf("mouse_");
				if (f >= 0) {
					if (myData.values[i][j] > maxLoadMouse) {
						maxLoadMouse = myData.values[i][j];
					}
				}
			}
			if (maxLoadHuman >= loadThreshold & maxLoadMouse >= loadThreshold) {
				topicCode[j] = 3;
			} else {
				if (maxLoadHuman >= loadThreshold) {
					topicCode[j] = 1;
				}
				if (maxLoadMouse >= loadThreshold) {
					topicCode[j] = 2;
				}
			}
		}
		
		try {
			FileWriter outFile = new FileWriter(graphName);
			outFile.write("diGraph G {\n");
			
			// output shared topics
			for (i=0;i<topicCode.length;i++) {
				if (topicCode[i] == 3) {
					nodeLabel = (new Integer(i)).toString();
					nodeIDS = "topic" + (new Integer(i)).toString();
					outFile.write(nodeIDS + " [label=\"" + nodeLabel + "\"];\n");
				}
			}
			
			// output human tissues
			outFile.write("subgraph cluster0 {\n");
			
			// output human only topics
			for (i=0;i<topicCode.length;i++) {
				if (topicCode[i] == 1) {
					nodeLabel = (new Integer(i)).toString();
					nodeIDS = "topic" + (new Integer(i)).toString();
					outFile.write(nodeIDS + " [label=\"" + nodeLabel + "\"];\n");
				}
			}
			
			for (i=0;i<myData.numRows;i++) {
				f = myData.geneNames[i].indexOf("human_");
				if (f >= 0) {
					nodeLabel = myData.geneNames[i].substring(f+6);
					nodeIDS = "human" + (new Integer(i)).toString();
					outFile.write(nodeIDS + " [label=\"" + nodeLabel + "\",shape=rectangle];\n");
					for (j=0;j<myData.numCols;j++) {
						if (myData.values[i][j] >= loadThreshold) {
							node2IDS = "topic" + (new Integer(j)).toString();
							outFile.write(nodeIDS + " -> " + node2IDS + ";\n");
						}
					}
				}
			}
			
			outFile.write("color=black;\n");
			outFile.write("label=\"H. sapiens\";\n");
			outFile.write("}\n");
			
			//	output mouse tissues
			outFile.write("subgraph cluster1 {\n");
			
			// output mouse only topics
			for (i=0;i<topicCode.length;i++) {
				if (topicCode[i] == 2) {
					nodeLabel = (new Integer(i)).toString();
					nodeIDS = "topic" + (new Integer(i)).toString();
					outFile.write(nodeIDS + " [label=\"" + nodeLabel + "\"];\n");
				}
			}
			
			for (i=0;i<myData.numRows;i++) {
				f = myData.geneNames[i].indexOf("mouse_");
				if (f >= 0) {
					nodeLabel = myData.geneNames[i].substring(f+6);
					nodeIDS = "mouse" + (new Integer(i)).toString();
					outFile.write(nodeIDS + " [label=\"" + nodeLabel + "\",shape=rectangle];\n");
					for (j=0;j<myData.numCols;j++) {
						if (myData.values[i][j] >= loadThreshold) {
							node2IDS = "topic" + (new Integer(j)).toString();
							outFile.write(nodeIDS + " -> " + node2IDS + ";\n");
						}
					}
				}
			}
			
			outFile.write("color=black;\n");
			outFile.write("label=\"M. muscularis\";\n");
			outFile.write("}\n");
			
			
			outFile.write("}\n");
			
			outFile.close();
		} catch (IOException e) {
			
		}
			
	}

}
