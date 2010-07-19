package edu.mit.csail.psrg.georg.Graph;

import java.awt.Container;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JFrame;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

public class OutputMouseHumanTabular extends JFrame {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		createTabular();
	}
	
	public OutputMouseHumanTabular(TabularComponent sheet) {
		Container cp = getContentPane();
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		
	//	cp.add(sheet);
		pack();
	}
	
	public static void createTabular() {
		String dirPath = "C:\\research_data\\mouse_human\\";
		String humanTopicOutFileS = dirPath + "tissue_as_docs_topics_go_human_short.txt";
		String mouseTopicOutFileS = dirPath + "tissue_as_docs_topics_go_mouse_short.txt";
		String normDPsFile = dirPath + "tissue_as_docs_dps_norm.txt";
		String tableNameJoint = dirPath + "tissue_as_docs_joint_topics.svg";
		
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
		
	//	outputWholeMatrix(myData,loadThreshold,topicCode,numTopics);
		outputTopicType(tableNameJoint,myData,loadThreshold,topicCode,3);
	}
	
	public static void outputTopicType(String fName,MicroArrayData myData,double loadThreshold,int[] topicCode,int topicType) {
		int numRows = 0;
		String speciesName = "";
		int i = 0;
		int j = 0;
		int k = 0;
		int humanStart = 3;
		int mouseStart = humanStart + 31;
		ArrayList<ArrayList<Integer>> topics = new ArrayList<ArrayList<Integer>>();
		ArrayList<Integer> topic = null;
		
		if (topicType == 3) {
			for (j=0;j<topicCode.length;j++) {
				if (topicCode[j] == 3) {
					topic = new ArrayList<Integer>();
					topics.add(topic);
					for (i=humanStart;i<mouseStart;i++) {
						if (myData.values[i][j] >= loadThreshold | myData.values[i+31][j] >= loadThreshold) {
							topic.add(i);
						}
					}
				}
			}
		}
		
		if (topicType == 2) {
			for (j=0;j<topicCode.length;j++) {
				if (topicCode[j] == 2) {
					topic = new ArrayList<Integer>();
					topics.add(topic);
					for (i=mouseStart;i<myData.numRows;i++) {
						if (myData.values[i][j] >= loadThreshold) {
							topic.add(i);
						}
					}
				}
			}
		}
		
		if (topicType == 1) {
			for (j=0;j<topicCode.length;j++) {
				if (topicCode[j] == 1) {
					topic = new ArrayList<Integer>();
					topics.add(topic);
					for (i=humanStart;i<mouseStart;i++) {
						if (myData.values[i][j] >= loadThreshold) {
							topic.add(i);
						}
					}
				}
			}
		}
		
		for (i=0;i<topics.size();i++) {
			topic = topics.get(i);
			numRows += topic.size();
			numRows++;
		}
		
		TabularComponent sheet = null;
		
		if (topicType == 3) {
			sheet = new TabularComponent(numRows,3,30);
			int[] colWidths = {100,10,10};
			sheet.setColWidths(colWidths);
		} else {
			sheet = new TabularComponent(numRows,2,30);
			int[] colWidths = {100,10};
			sheet.setColWidths(colWidths);
		}
		
		numRows = 0;
		String s = "";
		TabText cell = null;
		int numTopics = 0;
		double v = 0.0;
		for (j=0;j<topicCode.length;j++) {
			if (topicCode[j] == topicType) {
				topic = topics.get(numTopics);
				numTopics++;
				for (k=0;k<topic.size();k++) {
					i = topic.get(k);
					s = myData.geneNames[i].substring(6);
					cell = new TabText(s);
					sheet.setCell(numRows,0,cell);
					
					v = myData.values[i][j];
					cell = new TabText("");
					if (v < loadThreshold) {
						v = 0.0;
					}
					cell.setFillColor((float) (0.0),(float) (0.0),(float) (Math.pow(v,1.0/3.0)));
					cell.setFill(true);
					sheet.setCell(numRows,1,cell);
					
					if (topicType == 3) {
						v = myData.values[i+31][j];
						cell = new TabText("");
						if (v < loadThreshold) {
							v = 0.0;
						}
						cell.setFillColor((float) (0.0),(float) (0.0),(float) (Math.pow(v,1.0/3.0)));
						cell.setFill(true);
						sheet.setCell(numRows,2,cell);
					}
					
					numRows++;
				}
				numRows++;
			}
		}
	//	OutputMouseHumanTabular frame = new OutputMouseHumanTabular(sheet);
	//	frame.setVisible(true);
		try {
			sheet.outputSVG(fName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void outputWholeMatrix(MicroArrayData myData,double loadThreshold,int[] topicCode,int numTopics) {
		int i = 0;
		int j = 0;
		TabularComponent sheet = new TabularComponent(myData.numRows,numTopics+1,30);
		int[] colWidths = new int[numTopics+1];
		colWidths[0] = 150;
		for (i=1;i<numTopics+1;i++) {
			colWidths[i] = 10;
		}
		sheet.setColWidths(colWidths);
		TabText cell = null;
		double v = 0.0;
		for (i=0;i<myData.numRows;i++) {
			cell = new TabText(myData.geneNames[i]);
			sheet.setCell(i,0,cell);
			numTopics = 0;
			for (j=0;j<myData.numCols;j++) {
				v = myData.values[i][j];
				if (topicCode[j] > 0) {
					cell = new TabText("");
					if (v < loadThreshold) {
						v = 0.0;
					}
					cell.setFillColor((float) (0.0),(float) (0.0),(float) (Math.pow(v,1.0/3.0)));
					cell.setFill(true);
					sheet.setCell(i,numTopics+1,cell);
					numTopics++;
				}
			}
		}
		
		OutputMouseHumanTabular frame = new OutputMouseHumanTabular(sheet);
		frame.setVisible(true);
	}

}
