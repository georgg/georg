package edu.mit.csail.psrg.georg.DiscretizeData;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

public class DiscretizeMouseHuman2 {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	//	mergeExpressionData();
	//	thresholdMotifs();
		thresholdExpression();
		discretizeData();
	//	mergeSameMotifs();
	}
	
	public static void discretizeData() {
	//	String outputDir = "C:\\research_data\\mouse_human\\b47mm7\\";
	//	String discreteEName = "C:\\research_data\\mouse_human\\b47mm7\\expression_discrete.txt";
		String outputDir = "C:\\research_data\\mouse_human\\simulated\\";
		String discreteEName = "C:\\research_data\\mouse_human\\simulated\\expression_discrete.txt";
		
	//	String discreteMotifsName = "C:\\research_data\\mouse_human\\motifs\\motifs_discrete.txt";
		
		String combinedName = outputDir + "combined_discretized.txt";
		String levelsFileName = outputDir + "discretized_levels.txt";
		String MIFileName = outputDir + "discretized_MI.txt";
		
		MicroArrayData expression = new MicroArrayData();
		MicroArrayData motifs = new MicroArrayData();
		
		try {
			expression.readFile(discreteEName);
		//	motifs.readFile(discreteMotifsName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
	//	MicroArrayData mergedData = MicroArrayData.mergeCols(expression,motifs);
		MicroArrayData mergedData = expression;
	//	expression = null;
		motifs = null;
		
		DiscretizeExpression discretizer = new DiscretizeExpression(mergedData.values,10);
		
		try {
			discretizer.discretizeDownTree(MIFileName);
			discretizer.mergeDownLevels(2,levelsFileName);
			discretizer.transposeDExpression();
			mergedData.setDiscrete();
			mergedData.dvalues = discretizer.dExpression;
			mergedData.writeFile(combinedName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	// combine replicate PMA values (assumes two replicates, adjacent columns)
	public static MicroArrayData combinePMAReplicates(String fName) {
		MicroArrayData PMA = new MicroArrayData();
		PMA.setDiscrete();
		
		try {
			PMA.readFileBlind(fName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int[][] PMAc = new int[PMA.numRows][PMA.numCols/2];
		String[] experimentNames = new String[PMA.numCols/2];
		
		int i = 0;
		int j = 0;
		int k = 0;
		int v = 0;
		k = 0;
		j = 0;
		while (j<PMA.numCols-1) {
			for (i=0;i<PMAc.length;i++) {
				v = 0;
				if (PMA.dvalues[i][j] > 0 | PMA.dvalues[i][j+1] > 0) {
					v = 1;
				}
				PMAc[i][k] = v;
			}
			experimentNames[k] = PMA.experimentNames[j];
			k++;
			j = j + 2;
		}
		
		PMA.experimentNames = experimentNames;
		PMA.dvalues = PMAc;
		PMA.numCols = PMA.numCols/2;
		
		return PMA;
	}
	
	public static void mergeExpressionData() {
		String dirName = "C:\\research_data\\mouse_human\\b47mm7\\";
		String humanEName = "GNF1H.gcRMA.mean.txt";
		String humanPMAName = "GNF1H.pma.txt";
		String mouseEName = "GNF1M.gcRMA.mean.txt";
		String mousePMAName = "GNF1M.pma.txt";
		
		String mergedEName = "GNF1_merged_expression.txt";
		String mergedPMAName = "GNF1_merged_PMA.txt";
		
		MicroArrayData humanPMA = combinePMAReplicates(dirName + humanPMAName);
		MicroArrayData mousePMA = combinePMAReplicates(dirName + mousePMAName);
		int i = 0;
		
		for (i=0;i<humanPMA.numCols;i++) {
			humanPMA.experimentNames[i] = "h_" + humanPMA.experimentNames[i];
		}
		
		for (i=0;i<mousePMA.numCols;i++) {
			mousePMA.experimentNames[i] = "m_" + mousePMA.experimentNames[i];
		}
		
		MicroArrayData combinedPMA = MicroArrayData.mergeCols(humanPMA,mousePMA);
		
		for (i=0;i<combinedPMA.numRows;i++) {
			combinedPMA.geneNames[i] = humanPMA.geneNames[i] + "|" + mousePMA.geneNames[i];
		}
		
		try {
			combinedPMA.writeFile(dirName + mergedPMAName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		mergeMotifData(humanPMA.geneNames,mousePMA.geneNames);
		
		humanPMA = null;
		mousePMA = null;
		combinedPMA = null;
		
		MicroArrayData humanE = new MicroArrayData();
		MicroArrayData mouseE = new MicroArrayData();
				
		try {
			humanE.readFileBlind(dirName + humanEName);
			mouseE.readFileBlind(dirName + mouseEName);
			
			humanE.experimentCode = new int[humanE.numCols];
			for (i=0;i<humanE.numCols;i++) {
				humanE.experimentNames[i] = "h_" + humanE.experimentNames[i];
				humanE.experimentCode[i] = 1;
			}
			
			mouseE.experimentCode = new int[mouseE.numCols];
			for (i=0;i<mouseE.numCols;i++) {
				mouseE.experimentNames[i] = "m_" + mouseE.experimentNames[i];
				mouseE.experimentCode[i] = 2;
			}
			
			MicroArrayData combinedE = MicroArrayData.mergeCols(humanE,mouseE);
			for (i=0;i<combinedE.numRows;i++) {
				combinedE.geneNames[i] = humanE.geneNames[i] + "|" + mouseE.geneNames[i];
			}
			
			humanE = null;
			mouseE = null;
			combinedE.writeFile(dirName+mergedEName);
			combinedE = null;
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void mergeMotifData(String[] humanGeneNames,String[] mouseGeneNames) {
		String dirName = "C:\\research_data\\mouse_human\\motifs\\";
		String motifIndexName = dirName + "motif_index.txt";
		String humanMotifsName = dirName + "human_cons_2000up_1000down_3ord_ratios.txt";
		String mouseMotifsName = dirName + "mouse_cons_2000up_1000down_3ord_ratios.txt";
		String mergeMotifsName = dirName + "motifs_cons_merged.txt";
		
		MicroArrayData humanMotifs = loadMotifFile(humanMotifsName,humanGeneNames);
		MicroArrayData mouseMotifs = loadMotifFile(mouseMotifsName,mouseGeneNames);
		
		int i = 0;
		HashMap<String,String> motifMap = null;
		
		try {
			motifMap = loadMotifMap(motifIndexName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		humanMotifs.experimentCode = new int[humanMotifs.numCols];
		for (i=0;i<humanMotifs.numCols;i++) {
			humanMotifs.experimentNames[i] = "h_mot_" + motifMap.get(humanMotifs.experimentNames[i]);
			humanMotifs.experimentCode[i] = 3;
		}
		mouseMotifs.experimentCode = new int[mouseMotifs.numCols];
		for (i=0;i<mouseMotifs.numCols;i++) {
			mouseMotifs.experimentNames[i] = "m_mot_" + motifMap.get(mouseMotifs.experimentNames[i]);
			mouseMotifs.experimentCode[i] = 4;
		}
		
		MicroArrayData mergeMotifs = MicroArrayData.mergeCols(humanMotifs,mouseMotifs);
		for (i=0;i<mergeMotifs.numRows;i++) {
			mergeMotifs.geneNames[i] = humanGeneNames[i] + "|" + mouseGeneNames[i];
		}
		
		try {
			mergeMotifs.writeFile(mergeMotifsName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static MicroArrayData loadMotifFile(String fName,String[] geneNames) {
		MicroArrayData motifs = new MicroArrayData();
		
		try {
			motifs.readFileBlind(fName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		double[][] values = new double[geneNames.length][motifs.numCols];
		HashMap<String,Integer> geneMap = new HashMap<String,Integer>();
		int i = 0;
		int ii = 0;
		int j = 0;
		
		for (i=0;i<motifs.geneNames.length;i++) {
			geneMap.put(motifs.geneNames[i],i);
		}
		
		for (i=0;i<geneNames.length;i++) {
			if (geneMap.containsKey(geneNames[i])) {
				ii = geneMap.get(geneNames[i]);
				System.arraycopy(motifs.values[ii],0,values[i],0,motifs.numCols);
			} else {
				for (j=0;j<motifs.numCols;j++) {
					values[i][j] = 0.0/0.0;
				}
			}
		}
		
		motifs.geneNames = geneNames;
		motifs.values = values;
		motifs.numRows = motifs.geneNames.length;
		
		return motifs;
	}
	
	public static void thresholdExpression() {
	/*	String dirName = "C:\\research_data\\mouse_human\\b47mm7\\";
		String mergedEName = dirName + "GNF1_merged_expression.txt";
		String mergedPMAName = dirName + "GNF1_merged_PMA.txt";
		String discreteEName = dirName + "expression_discrete.txt"; */
		
		String dirName = "C:\\research_data\\mouse_human\\simulated\\";
		String mergedEName = dirName + "overlap1.txt";
		String discreteEName = dirName + "expression_discrete.txt";
		
	//	double expressionThresh = 1.5;
		double expressionThresh = 7.24;
		
	//	MicroArrayData PMA = new MicroArrayData();
	//	PMA.setDiscrete();
		
		MicroArrayData expression = new MicroArrayData();
		
		try {
		//	PMA.readFile(mergedPMAName);
			expression.readFile(mergedEName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		int j = 0;
		
		for (i=0;i<expression.numRows;i++) {
			for (j=0;j<expression.numCols;j++) {
		//		if (expression.values[i][j] < expressionThresh | PMA.dvalues[i][j] == 0) {
				if (expression.values[i][j] < expressionThresh) {
					expression.values[i][j] = 0.0/0.0;
				}
			}
		}
		
		try {
			expression.writeFile(discreteEName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void thresholdMotifs() {
		String dirName = "C:\\research_data\\mouse_human\\motifs\\";
		String mergedMotifsName = dirName + "motifs_cons_merged.txt";
		String discreteMotifsName = dirName + "motifs_discrete.txt";
		
		double percentMax = 0.60;
		double presentThresh = 0.25;
		
		MicroArrayData motifs = new MicroArrayData();
		try {
			motifs.readFile(mergedMotifsName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int[] numHits = new int[motifs.numCols];
		Vector<Double> v = new Vector<Double>();
		int maxNumPresent = (int) Math.round(((double) motifs.numRows)*presentThresh);
		int i = 0;
		int j = 0;
		double maxV = 0.0;
		
		for (j=0;j<motifs.numCols;j++) {
			maxV = -10e8;
			for (i=0;i<motifs.numRows;i++) {
				if (motifs.values[i][j] > maxV) {
					maxV = motifs.values[i][j];
				}
			}
			maxV = percentMax * maxV;
			for (i=0;i<motifs.numRows;i++) {
				if (motifs.values[i][j] >= maxV) {
					numHits[j]++;
				} else {
					motifs.values[i][j] = 0.0/0.0;
				}
				
			}
		}
		
		HashSet<Integer> delMotifs = new HashSet<Integer>();
		
		for (j=0;j<motifs.numCols;j++) {
			if (numHits[j] > maxNumPresent | numHits[j] == 0) {
				delMotifs.add(j);
			}
		}
		motifs.deleteCols(delMotifs);
		
		try {
			motifs.writeFile(discreteMotifsName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static HashMap<String,String> loadMotifMap(String motifMapName) throws IOException {
		HashMap<String,String> motifMap = new HashMap<String,String>();
		String[] lineSplit = null;
		BufferedReader is = new BufferedReader(new FileReader(motifMapName));
		String line = "";		
		line = is.readLine();
		while(line != null) {
			lineSplit = line.split("\t");
			motifMap.put(lineSplit[0],lineSplit[1]);
			line = is.readLine();
		}
		is.close();
		return motifMap;
	}
	
	public static void mergeSameMotifs() {
		String outputDir = "C:\\research_data\\mouse_human\\b47mm7\\";
		String combinedName = outputDir + "combined_discretized.txt";
		String combinedNameFinal = outputDir + "GNF1_discretized_data.txt";
		
		MicroArrayData combinedData = new MicroArrayData();
		combinedData.setDiscrete();
		try {
			combinedData.readFile(combinedName);
		} catch(IOException e) { System.out.println(e); }
		
		ArrayList<Integer> humanMotifs = new ArrayList<Integer>();
		ArrayList<Integer> mouseMotifs = new ArrayList<Integer>();
		HashSet<Integer> delCols = new HashSet<Integer>();
		
		int f = 0;
		int i = 0;
		
		for (i=0;i<combinedData.numCols;i++) {
			f = combinedData.experimentNames[i].indexOf("h_mot_");
			if (f > -1) {
				humanMotifs.add(i);
				delCols.add(i);
			}
			f = combinedData.experimentNames[i].indexOf("m_mot_");
			if (f > -1) {
				mouseMotifs.add(i);
				delCols.add(i);
			}
		}
		
		MicroArrayData humanMotifData = mergeSameMotifsUtil(combinedData,humanMotifs);
		for (i=0;i<humanMotifData.numCols;i++) {
			humanMotifData.experimentCode[i] = 3;
			humanMotifData.experimentNames[i] = "h_mot_" + humanMotifData.experimentNames[i];
		}
		
		MicroArrayData mouseMotifData = mergeSameMotifsUtil(combinedData,mouseMotifs);
		for (i=0;i<mouseMotifData.numCols;i++) {
			mouseMotifData.experimentCode[i] = 4;
			mouseMotifData.experimentNames[i] = "m_mot_" + mouseMotifData.experimentNames[i];
		}
		
		combinedData.deleteCols(delCols);
		combinedData = MicroArrayData.mergeCols(combinedData,humanMotifData);
		combinedData = MicroArrayData.mergeCols(combinedData,mouseMotifData);
		
		try {
			combinedData.writeFile(combinedNameFinal);
		} catch(IOException e) { System.out.println(e); }
	}
	
	public static String getMotifName(String motifName) {
		int f = 0;
		f = motifName.indexOf("V$");
		return motifName.substring(f);
	}
	
	public static MicroArrayData mergeSameMotifsUtil(MicroArrayData data,ArrayList<Integer> colList) {
		HashMap<String,ArrayList<Integer>> motifMap = new HashMap<String,ArrayList<Integer>>();
		ArrayList<Integer> motifIDXs = null;
		
		int i = 0;
		int ii = 0;
		int j = 0;
		int k = 0;
		int f = 0;
		String s;
		
		for (ii=0;ii<colList.size();ii++) {
			i = colList.get(ii);
			s = getMotifName(data.experimentNames[i]);
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
		int[][] values = new int[data.numRows][motifMap.size()];
		for (i=0;i<values.length;i++) {
			for (j=0;j<uMotifs.length;j++) {
				motifIDXs = motifMap.get(uMotifs[j]);
				maxValue = 0;
				for (k=0;k<motifIDXs.size();k++) {
					if (data.dvalues[i][motifIDXs.get(k)] > maxValue) {
						maxValue = data.dvalues[i][motifIDXs.get(k)];
					}
				}
				values[i][j] = maxValue;
			}
		}
		
		MicroArrayData newData = new MicroArrayData();
		
		newData.dvalues = values;
		newData.discreteData = true;
		newData.numCols = uMotifs.length;
		newData.experimentNames = uMotifs;
		newData.geneNames = data.geneNames;
		newData.numRows = data.numRows;
		newData.experimentCode = new int[newData.numCols];
		
		return newData;
	}
}
