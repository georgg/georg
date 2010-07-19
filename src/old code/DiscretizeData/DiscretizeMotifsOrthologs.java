package edu.mit.csail.psrg.georg.DiscretizeData;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

public class DiscretizeMotifsOrthologs {
	
	public static void main(String[] args) {
	//	processData();
	//	discretizeFiles();
	//	mergeData();
		mergeMotifs();
	}
	
	public static void processData() {
		String humanExpressionInFileName = "C:\\research_data\\mouse_human\\homo_b47_data\\GNF1H.gcrma.snco.txt";
		String mouseExpressionInFileName = "C:\\research_data\\mouse_human\\homo_b47_data\\GNF1M.gcrma.snco.txt";
		String humanMotifInFileName = "C:\\research_data\\mouse_human\\motifs\\human_2000bp_6ord_ratios.txt.trimmed";
		String mouseMotifInFileName = "C:\\research_data\\mouse_human\\motifs\\mouse_2000bp_6ord_ratios.txt.trimmed";
		
		String humanMotifOutFileName = "C:\\research_data\\mouse_human\\motifs_orthologs\\human_2000bp_6ord_discretized.txt";
		String mouseMotifOutFileName = "C:\\research_data\\mouse_human\\motifs_orthologs\\mouse_2000bp_6ord_discretized.txt";
		
		String[] humanENames = null;
		String[] mouseENames = null;
		String[] humanMNames = null;
		String[] mouseMNames = null;
		
		MicroArrayData humanExpression = new MicroArrayData();
		MicroArrayData mouseExpression = new MicroArrayData();
		MicroArrayData humanMotifs = new MicroArrayData();
		MicroArrayData mouseMotifs = new MicroArrayData();
		try {
			humanExpression.readFile(humanExpressionInFileName);
			humanENames = humanExpression.geneNames;
			humanExpression = null;
			mouseExpression.readFile(mouseExpressionInFileName);
			mouseENames = mouseExpression.geneNames;
			mouseExpression = null;
			humanMotifs.readFileBlind(humanMotifInFileName);
			humanMNames = humanMotifs.geneNames;
			mouseMotifs.readFileBlind(mouseMotifInFileName);
			mouseMNames = mouseMotifs.geneNames;
		} catch(IOException e) {
			System.out.println(e);
		}
		
		ArrayList<Integer> humanGenes = new ArrayList<Integer>();
		ArrayList<Integer> mouseGenes = new ArrayList<Integer>();
		
		int i = 0;
		int j = 0;
		HashMap<String,String> orthologMap = new HashMap<String,String>();
		for (i=0;i<humanENames.length;i++) {
			orthologMap.put(humanENames[i],mouseENames[i]);
		}
		
		HashMap<String,Integer> mouseGeneMap = new HashMap<String,Integer>();
		for (i=0;i<mouseMNames.length;i++) {
			mouseGeneMap.put(mouseMNames[i],i);
		}
		String orthoName;
		for (i=0;i<humanMNames.length;i++) {
			orthoName = orthologMap.get(humanMNames[i]);
			if (mouseGeneMap.containsKey(orthoName)) {
				humanGenes.add(i);
				mouseGenes.add(mouseGeneMap.get(orthoName));
			}
		}
		
		String[] humanMNames2 = new String[humanGenes.size()];
		double[][] humanValues = new double[humanGenes.size()][humanMotifs.numCols];
		
		String[] mouseMNames2 = new String[humanGenes.size()];
		double[][] mouseValues = new double[humanGenes.size()][mouseMotifs.numCols];
		
		int m_i;
		int h_i;
		for (i=0;i<humanGenes.size();i++) {
			h_i = humanGenes.get(i);
			m_i = mouseGenes.get(i);
			humanMNames2[i] = humanMNames[h_i];
			mouseMNames2[i] = mouseMNames[m_i];
			System.arraycopy(humanMotifs.values[h_i],0,humanValues[i],0,humanMotifs.numCols);
			System.arraycopy(mouseMotifs.values[m_i],0,mouseValues[i],0,mouseMotifs.numCols);
		}
		
		DiscretizeMouseHumanMotifs.thresholdMotifs(humanValues);
		DiscretizeMouseHumanMotifs.thresholdMotifs(mouseValues);
		
		humanMotifs.values = humanValues;
		humanMotifs.geneNames = humanMNames2;
		humanMotifs.numRows = humanValues.length;
		mouseMotifs.values = mouseValues;
		mouseMotifs.geneNames = mouseMNames2;
		mouseMotifs.numRows = mouseValues.length;
		
		try {
			humanMotifs.writeFile(humanMotifOutFileName);
			mouseMotifs.writeFile(mouseMotifOutFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void discretizeFiles() {
		String humanMotifFileName = "C:\\research_data\\mouse_human\\motifs_orthologs\\human_2000bp_6ord_discretized.txt";
		String humanMIFileName = "C:\\research_data\\mouse_human\\motifs_orthologs\\human_MI.txt";
		String humanLevelsFileName = "C:\\research_data\\mouse_human\\motifs_orthologs\\human_levels.txt";
		String mouseMIFileName = "C:\\research_data\\mouse_human\\motifs_orthologs\\mouse_MI.txt";
		String mouseLevelsFileName = "C:\\research_data\\mouse_human\\motifs_orthologs\\mouse_levels.txt";
		String mouseMotifFileName = "C:\\research_data\\mouse_human\\motifs_orthologs\\mouse_2000bp_6ord_discretized.txt";
		
		discretizeData(humanMotifFileName,humanMotifFileName,humanMIFileName,humanLevelsFileName);
		discretizeData(mouseMotifFileName,mouseMotifFileName,mouseMIFileName,mouseLevelsFileName);
	}
	
	public static void discretizeData(String inFile,String outFile,String MIFileName,String levelsFileName) {
		MicroArrayData file = new MicroArrayData();
		
		try {
			file.readFile(inFile);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		DiscretizeExpression discretizer = new DiscretizeExpression(file.values,10);
		
		try {
			discretizer.discretizeDownTree(MIFileName);
			discretizer.mergeDownLevels(2,levelsFileName);
			discretizer.transposeDExpression();
			file.setDiscrete();
			file.dvalues = discretizer.dExpression;
			file.writeFile(outFile);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void mergeData() {
		String humanMotifFileName = "C:\\research_data\\mouse_human\\motifs_orthologs\\human_2000bp_6ord_discretized.txt";
		String mouseMotifFileName = "C:\\research_data\\mouse_human\\motifs_orthologs\\mouse_2000bp_6ord_discretized.txt";
		String combinedMotifFileName = "C:\\research_data\\mouse_human\\motifs_orthologs\\combined_2000bp_6ord_discretized.txt";
		String motifMapFileName = "C:\\research_data\\mouse_human\\motifs\\motif_index.txt";
		
		MicroArrayData humanMotifs = new MicroArrayData();
		humanMotifs.setDiscrete();
		MicroArrayData mouseMotifs = new MicroArrayData();
		mouseMotifs.setDiscrete();
		
		HashMap<String,String> motifMap = null;
		
		try {
			humanMotifs.readFile(humanMotifFileName);
			mouseMotifs.readFile(mouseMotifFileName);
			motifMap = DiscretizeMouseHumanMotifs.loadMotifMap(motifMapFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int[][] combinedValues = new int[humanMotifs.numRows*2][humanMotifs.numCols];
		String[] combinedNames = new String[humanMotifs.numRows*2];
		int i = 0;
		
		for (i=0;i<humanMotifs.numRows;i++) {
			System.arraycopy(humanMotifs.dvalues[i],0,combinedValues[i],0,humanMotifs.numCols);
			System.arraycopy(mouseMotifs.dvalues[i],0,combinedValues[i+humanMotifs.numRows],0,humanMotifs.numCols);
			combinedNames[i] = humanMotifs.geneNames[i];
			combinedNames[i+humanMotifs.numRows] = mouseMotifs.geneNames[i];
		}
		
		humanMotifs.dvalues = combinedValues;
		humanMotifs.geneNames = combinedNames;
		humanMotifs.numRows = humanMotifs.geneNames.length;
		
		for (i=0;i<humanMotifs.numCols;i++) {
			humanMotifs.experimentNames[i] = motifMap.get(humanMotifs.experimentNames[i]);
		}
		
		try {
			humanMotifs.writeFile(combinedMotifFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void mergeMotifs() {
		String combinedMotifFileName = "C:\\research_data\\mouse_human\\motifs_orthologs\\combined_2000bp_6ord_discretized.txt";
		double presentThresh = 0.50;
		
		MicroArrayData motifData = new MicroArrayData();
		motifData.setDiscrete();
		try {
			motifData.readFile(combinedMotifFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		HashMap<String,ArrayList<Integer>> motifMap = new HashMap<String,ArrayList<Integer>>();
		ArrayList<Integer> motifIDXs = null;
		
		int i = 0;
		int j = 0;
		int k = 0;
		int f = 0;
		String s;
		
		for (i=0;i<motifData.experimentNames.length;i++) {
			s = motifData.experimentNames[i];
			f = s.indexOf("_");
			if (f == -1) {
				s = s.substring(2);
			} else {
				s = s.substring(2,f);
			}
			if (!motifMap.containsKey(s)) {
				motifIDXs = new ArrayList<Integer>();
				motifMap.put(s,motifIDXs);
			}
			motifIDXs = motifMap.get(s);
			motifIDXs.add(i);
		}
		
		String[] uMotifs = new String[motifMap.size()];
		uMotifs = motifMap.keySet().toArray(uMotifs);
		int maxValue = 0;
		int[][] values = new int[motifData.numRows][motifMap.size()];
		for (i=0;i<values.length;i++) {
			for (j=0;j<uMotifs.length;j++) {
				motifIDXs = motifMap.get(uMotifs[j]);
				maxValue = 0;
				for (k=0;k<motifIDXs.size();k++) {
					if (motifData.dvalues[i][motifIDXs.get(k)] > maxValue) {
						maxValue = motifData.dvalues[i][motifIDXs.get(k)];
					}
				}
				values[i][j] = maxValue;
			}
		}
		
		ArrayList<Integer> keepMotifs = new ArrayList<Integer>();
		int numPresent = 0;
		int maxNumPresent = (int) Math.round(((double) values.length)*presentThresh);
		
		for (j=0;j<values[0].length;j++) {
			numPresent = 0;
			for (i=0;i<values.length;i++) {
				if (values[i][j] > 0) {
					numPresent++;
				}
			}
			if (numPresent < maxNumPresent) {
				keepMotifs.add(j);
			}
		}
		int[][] values2 = new int[values.length][keepMotifs.size()];
		String[] uMotifs2 = new String[keepMotifs.size()];
		for (k=0;k<keepMotifs.size();k++) {
			j = keepMotifs.get(k);
			uMotifs2[k] = uMotifs[j];
			for (i=0;i<values.length;i++) {
				values2[i][k] = values[i][j];
			}
		}
		values = values2;
		uMotifs = uMotifs2;
		
		motifData.dvalues = values;
		motifData.numCols = uMotifs.length;
		motifData.experimentNames = uMotifs;
		
		try {
			motifData.writeFile(combinedMotifFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
}
