package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

public class OutputGroupsGraph {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	//	outputGraph();
		outputYFilesGraph();
	//	outputMatrix();
	}

	public static void outputMatrix() {
		String dirPath = "C:\\research_data\\mouse_human\\b47mm7\\geo\\";
		String snapFileName = dirPath + "tissue_as_docs_snap_149999.persist";
		String matrixFileName = dirPath + "tissues_prob_matrix.txt";
		
		boolean useProbs = true;
		HierDirichletProcess myHDP = HierDirichletProcess.restoreFromFile(snapFileName);
		int g = 0;
		int ii = 0;
		int jj = 0;
		int ng1 = 0;
		int ng2 = 0;
		
		DPDistanceCollector myDistanceCollector = null;
		NumClustersCollector myNumClustersCollector = null;
		
		if (!useProbs) {
			myDistanceCollector = (DPDistanceCollector) myHDP.getStatCollector(DPDistanceCollector.class);
			myNumClustersCollector = (NumClustersCollector) myHDP.getStatCollector(NumClustersCollector.class);
			for (ii=0;ii<myDistanceCollector.pairProb.length;ii++) {
				for (jj=0;jj<myDistanceCollector.pairProb[ii].length;jj++) {
					myDistanceCollector.pairProb[ii][jj] = myDistanceCollector.pairProb[ii][jj]/((double) myNumClustersCollector.numClustersSamples.size());
				}
			}
		}
		
		DPGroup myGroup = null;
		double v = 0.0;
		int numTissues = 0;
		
		if (useProbs) {
			myGroup = myHDP.DPGroups.get(0);
			myGroup.normPairProbs();
			numTissues = myGroup.pairProb.length + 1;
		} else {
			numTissues = myDistanceCollector.pairProb.length + 1;
		}
		
		double[][] values = new double[numTissues][numTissues];
		
		if (useProbs) {
			ng1 = myGroup.pairProb.length;
		} else {
			ng1 = myDistanceCollector.pairProb.length;
		}
		for (ii=0;ii<ng1;ii++) {
			if (useProbs) {
				ng2 = myGroup.pairProb[ii].length;
			} else {
				ng2 = myDistanceCollector.pairProb[ii].length;
			}
			for (jj=0;jj<ng2;jj++) {
				if (useProbs) {
					v = myGroup.pairProb[ii][jj];
				} else {
					v = 1.0-myDistanceCollector.pairProb[ii][jj];
				}
				values[ii][ii+jj+1] = v;
				values[ii+jj+1][ii] = v;
			}
		}
		
		MicroArrayData matrix = new MicroArrayData();
		matrix.values = values;
		matrix.numRows = values.length;
		matrix.numCols = values[0].length;
		
		try {
			matrix.writeFile(matrixFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void outputGraph() {
		String dirPath = "C:\\research_data\\mouse_human\\b47mm7\\geo\\";
		String snapFileName = dirPath + "tissue_as_docs_snap_149999.persist";
		String graphFileName = dirPath + "groups_biolayout.net";
		
		HierDirichletProcess myHDP = HierDirichletProcess.restoreFromFile(snapFileName);
		int g = 0;
		int ii = 0;
		int jj = 0;
		int ng1 = 0;
		int ng2 = 0;
		DPGroup myGroup = null;
		double thresh = 0.50;
		String name1;
		String name2;
		double v = 0.0;
		int species = 0;
		
		DPDistanceCollector myDistanceCollector = (DPDistanceCollector) myHDP.getStatCollector(DPDistanceCollector.class);
		NumClustersCollector myNumClustersCollector = (NumClustersCollector) myHDP.getStatCollector(NumClustersCollector.class);
		
		for (ii=0;ii<myDistanceCollector.pairProb.length;ii++) {
			for (jj=0;jj<myDistanceCollector.pairProb[ii].length;jj++) {
				myDistanceCollector.pairProb[ii][jj] = myDistanceCollector.pairProb[ii][jj]/((double) myNumClustersCollector.numClustersSamples.size());
			}
		}
		
		try {
			FileWriter file = new FileWriter(graphFileName);
			
			for (g=0;g<myHDP.DPGroups.size();g++) {
				myGroup = myHDP.DPGroups.get(g);
				myGroup.normPairProbs();
				
				file.write("*Vertices ");
				file.write((new Integer(myGroup.groupableDPs.size())).toString());
				file.write("\n");
				
				for (ii=0;ii<myGroup.groupableDPs.size();ii++) {
					file.write((new Integer(ii+1)).toString());
					name1 = stripString(myGroup.groupableDPs.get(ii).label);
					file.write(" \"" + name1 + "\" 0.0 0.0 0.0 ic ");
					species = speciesCode(myGroup.groupableDPs.get(ii).label);
					if (species == 1) {
						file.write("Red");
					}
					if (species == 2) {
						file.write("Green");
					}
					file.write("\n");
				}
				
				file.write("*Edges\n");
				
			//	ng1 = myGroup.pairProb.length;
				ng1 = myDistanceCollector.pairProb.length;
				for (ii=0;ii<ng1;ii++) {
					ng2 = myDistanceCollector.pairProb[ii].length;
					for (jj=0;jj<ng2;jj++) {
					//	v = 1.0-myDistanceCollector.pairProb[ii][jj];
						v = myDistanceCollector.pairProb[ii][jj];
						if (v >= thresh) {
							file.write((new Integer(ii+1)).toString());
							file.write(" ");
							file.write((new Integer(ii+jj+1+1)).toString());
							file.write(" ");
							file.write((new Double(v)).toString());
							file.write("\n");
						}
					}
				}
			}
			
			file.close();
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static String stripString(String s) {
		int f = 0;
		f = s.indexOf("_");
		return s.substring(f+1);
	}
	
	public static int speciesCode(String s) {
		int f = 0;
		f = s.indexOf("m_");
		if (f == 0)
			return 2;
		
		return 1;
	}
	
	public static void outputYFilesGraph() {
		String dirPath = "C:\\research_data\\mouse_human\\b47mm7\\geo\\";
		String snapFileName = dirPath + "tissue_as_docs_snap_149999.persist";
		String outputFile = dirPath + "tissues_structure.gml";
		String groupsFile = dirPath + "tissue_load.txt";
		
		int ii = 0;
		int jj = 0;
		HierDirichletProcess myHDP = HierDirichletProcess.restoreFromFile(snapFileName);
		DPDistanceCollector myDistanceCollector = (DPDistanceCollector) myHDP.getStatCollector(DPDistanceCollector.class);
		NumClustersCollector myNumClustersCollector = (NumClustersCollector) myHDP.getStatCollector(NumClustersCollector.class);
		int numGroups = 0;
		DPGroup myGroup = myHDP.DPGroups.get(0);
		myGroup.normPairProbs();
		
		for (ii=0;ii<myDistanceCollector.pairProb.length;ii++) {
			for (jj=0;jj<myDistanceCollector.pairProb[ii].length;jj++) {
				myDistanceCollector.pairProb[ii][jj] = myDistanceCollector.pairProb[ii][jj]/((double) myNumClustersCollector.numClustersSamples.size());
			}
		}
		
		try {
			HashMap<String,Integer> groupsIDMap = readGroupsIDs(groupsFile);
			FileWriter file = new FileWriter(outputFile);
			outputYFilesHeader(file);
			numGroups = outputNodes(file,groupsIDMap,myDistanceCollector.DPs);
		//	outputEdges(file,numGroups,myDistanceCollector.pairProb,groupsIDMap,myDistanceCollector.DPs);
			outputEdges(file,numGroups,myGroup.pairProb,groupsIDMap,myDistanceCollector.DPs);
			file.write("]\n");
			file.close();
		} catch(IOException e) {
			System.out.println(e);
		}
		
	}
	
	public static HashMap<String,Integer> readGroupsIDs(String fName) throws IOException {
		HashMap<String,Integer> groupIDMap = new HashMap<String,Integer>();
		MicroArrayData dataFile = new MicroArrayData();
		
		dataFile.readFile(fName);
		
		int i = 0;
		for (i=0;i<dataFile.numRows;i++) {
			groupIDMap.put(dataFile.geneNames[i],new Double(dataFile.values[i][0]).intValue());
		}
		
		return groupIDMap;
	}
	
	public static void outputYFilesHeader(FileWriter file) throws IOException {
		file.write("Creator\t\"yFiles\"\nVersion	2.4\n");
		file.write("graph\n[\n");
		file.write("directed 0\n");
		file.write("hierarchic 1\n");
	}
	
	public static int outputNodes(FileWriter file,HashMap<String,Integer> groupIDMap,DirichletProcess[] DPs) throws IOException {
		int n = 0;
		int i = 0;
		int speciesCode = 0;
		Iterator<Integer> iter = groupIDMap.values().iterator();
		while(iter.hasNext()) {
			i = iter.next();
			if (i > n) {
				n=i;
			}
		}
		n++;
		
		for (i=0;i<n;i++) {
			file.write("node\n[\n");
			file.write("id " + (new Integer(i+1)).toString() + "\n");
			file.write("label " + "\"" + "group#" + ((new Integer(i+1))).toString() + "\"\n");
			file.write("isGroup 1\n");
			file.write("]\n");
		}
		
		String label = "";
		String color = "";
		int group = 0;
		for (i=0;i<DPs.length;i++) {
			label = DPs[i].label;
			speciesCode = speciesCode(label);
			group = groupIDMap.get(label);
			file.write("node\n[\n");
			file.write("id " + (new Integer(i+n+1)).toString() + "\n");
			file.write("label " + "\"" + stripString(label) + "\"\n");
			file.write("gid " + (new Integer(group+1)).toString() + "\n");
			file.write("graphics [\n");
			file.write("type \"circle\"\n");
			
			if (speciesCode == 1) {
				color = "#0000FF";
			}
			if (speciesCode == 2) {
				color = "#FF0000";
			}
			file.write("fill \"" + color + "\"\n");
			file.write("h 10.0\n");
			file.write("w 10.0\n");
			
			file.write("]\n");
			file.write("]\n");
		}
		
		return n;
	}
	
	public static void outputEdges(FileWriter file,int numGroups,double[][] pairProb, HashMap<String,Integer> groupIDMap,DirichletProcess[] DPs) throws IOException {
		int ii = 0;
		int jj = 0;
		double thresh = 0.50;
		
		int ng1 = pairProb.length;
		int ng2 = 0;
		double v = 0.0;
		int source = 0;
		int destination = 0;
		int sourceGroup = 0;
		int destinationGroup = 0;
		
		double maxEdgeWidth = 3.0;
		int edgeWidth = 0;
		
		for (ii=0;ii<ng1;ii++) {
			ng2 = pairProb[ii].length;
			source = ii;
			sourceGroup = groupIDMap.get(DPs[source].label);
			for (jj=0;jj<ng2;jj++) {
				v = 1.0-pairProb[ii][jj];
				destination = ii + jj + 1;
				destinationGroup = groupIDMap.get(DPs[destination].label);
				if (v >= thresh & (sourceGroup != destinationGroup)) {
					v = (v-thresh)/thresh;
					edgeWidth = (int) Math.round(maxEdgeWidth * v) + 1;
					file.write("edge\n[\n");
					file.write("source " + (new Integer(ii+1+numGroups)).toString() + "\n");
					file.write("target " + (new Integer(ii+jj+1+1+numGroups)).toString() + "\n");
					file.write("graphics [\n");
					file.write("width " + (new Integer(edgeWidth)).toString() + "\n");
					file.write("]\n");
					file.write("]\n");
				}
			}
		}
	}
}
