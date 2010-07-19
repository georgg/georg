package edu.mit.csail.psrg.georg.DiscretizeData;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.StatUtil.VectorUtil;

public class DiscretizeMouseHumanMotifs {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	//	reformatMotifFile();
		processData();
		discretizeData();
		mergeMotifs();
	}
	
	public static void discretizeData() {
		String humanMotifOutFileName = "C:\\research_data\\mouse_human\\motifs\\human_cons_discretized.txt";
		String mouseMotifOutFileName = "C:\\research_data\\mouse_human\\motifs\\mouse_cons_discretized.txt";
		String combinedMotifOutFileName = "C:\\research_data\\mouse_human\\motifs\\combined_cons_discretized.txt";
		String levelsFileName = "C:\\research_data\\mouse_human\\motifs\\discretized_levels.txt";
		String MIFileName = "C:\\research_data\\mouse_human\\motifs\\discretized_MI.txt";
		
		MicroArrayData humanMotifData = new MicroArrayData();
		MicroArrayData mouseMotifData = new MicroArrayData();
		
		try {
			humanMotifData.readFileBlind(humanMotifOutFileName);
			mouseMotifData.readFileBlind(mouseMotifOutFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		double[][] values2 = new double[humanMotifData.numRows + mouseMotifData.numRows][humanMotifData.numCols];
		String[] geneNames = new String[humanMotifData.numRows + mouseMotifData.numRows];
		int i = 0;
		int j = 0;
		for (i=0;i<humanMotifData.numRows;i++) {
			values2[i] = humanMotifData.values[i];
			geneNames[i] = humanMotifData.geneNames[i];
		}
		for (i=0;i<mouseMotifData.numRows;i++) {
			values2[i+humanMotifData.numRows] = mouseMotifData.values[i];
			geneNames[i+humanMotifData.numRows] = mouseMotifData.geneNames[i];
		}
		
		DiscretizeExpression discretizer = new DiscretizeExpression(values2,10);
		MicroArrayData mergedData = new MicroArrayData();
		mergedData.geneNames = geneNames;
		mergedData.experimentNames = humanMotifData.experimentNames;
		mergedData.numCols = mergedData.experimentNames.length;
		mergedData.numRows = geneNames.length;
		
	/*	int[][] discreteData = new int[values2.length][values2[0].length];
		for (i=0;i<values2.length;i++) {
			for (j=0;j<values2[0].length;j++) {
				if (!Double.isNaN(values2[i][j])) {
					discreteData[i][j] = 1;
				}
			}
		} */
		
		try {
			discretizer.discretizeDownTree(MIFileName);
			discretizer.mergeDownLevels(2,levelsFileName);
			discretizer.transposeDExpression();
			mergedData.setDiscrete();
			mergedData.dvalues = discretizer.dExpression;
		//	mergedData.dvalues = discreteData;
			mergedData.writeFile(combinedMotifOutFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void reformatMotifFile() {
		String humanMotifInFileName = "C:\\research_data\\mouse_human\\motifs\\human_cons_2000up_1000down_3ord_ratios.txt";
		String mouseMotifInFileName = "C:\\research_data\\mouse_human\\motifs\\mouse_cons_2000up_1000down_3ord_ratios.txt";
		String humanMotifOutFileName = "C:\\research_data\\mouse_human\\motifs\\human_cons_2000up_1000down_3ord_ratios_rename.txt";
		String mouseMotifOutFileName = "C:\\research_data\\mouse_human\\motifs\\mouse_cons_2000up_1000down_3ord_ratios_rename.txt";
		String motifMapFileName = "C:\\research_data\\mouse_human\\motifs\\motif_index.txt";	
		
		HashMap<String,String> motifMap = null;
		MicroArrayData humanMotifData = new MicroArrayData();
		MicroArrayData mouseMotifData = new MicroArrayData();
		try {
			motifMap = loadMotifMap(motifMapFileName);
			humanMotifData.readFileBlind(humanMotifInFileName);
			mouseMotifData.readFileBlind(mouseMotifInFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		
		for (i=0;i<humanMotifData.numCols;i++) {
		//	humanMotifData.experimentNames[i] = motifMap.get(humanMotifData.experimentNames[i]);
		}
		
		for (i=0;i<mouseMotifData.numCols;i++) {
		//	mouseMotifData.experimentNames[i] = motifMap.get(mouseMotifData.experimentNames[i]);
		}
		
		try {
			humanMotifData.writeFile(humanMotifOutFileName);
			mouseMotifData.writeFile(mouseMotifOutFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void processData() {
		String topGenesFileName = "c:\\research_data\\mouse_human\\b47mm7\\geo\\top_genes_50.txt";
		String humanMotifInFileName = "C:\\research_data\\mouse_human\\motifs\\human_cons_2000up_1000down_3ord_ratios.txt";
		String mouseMotifInFileName = "C:\\research_data\\mouse_human\\motifs\\mouse_cons_2000up_1000down_3ord_ratios.txt";
		String humanMotifOutFileName = "C:\\research_data\\mouse_human\\motifs\\human_cons_discretized.txt";
		String mouseMotifOutFileName = "C:\\research_data\\mouse_human\\motifs\\mouse_cons_discretized.txt";
		String motifMapFileName = "C:\\research_data\\mouse_human\\motifs\\motif_index.txt";	
		
		String[] geneNames = null;
		HashMap<String,String> motifMap = null;
		MicroArrayData humanMotifData = new MicroArrayData();
		MicroArrayData mouseMotifData = new MicroArrayData();
		try {
			geneNames = getGenes(topGenesFileName);
			motifMap = loadMotifMap(motifMapFileName);
			humanMotifData.readFileBlind(humanMotifInFileName);
			mouseMotifData.readFileBlind(mouseMotifInFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		int j = 0;
		
		thresholdMotifs(humanMotifData.values);
		thresholdMotifs(mouseMotifData.values);
		
		HashMap<String,Integer> humanGeneMap = new HashMap<String,Integer>();
		HashMap<String,Integer> mouseGeneMap = new HashMap<String,Integer>();
		
		for (i=0;i<humanMotifData.numRows;i++) {
			humanGeneMap.put(humanMotifData.geneNames[i],i);
		}
		for (i=0;i<mouseMotifData.numRows;i++) {
			mouseGeneMap.put(mouseMotifData.geneNames[i],i);
		}
		
		for (i=0;i<humanMotifData.numCols;i++) {
			humanMotifData.experimentNames[i] = motifMap.get(humanMotifData.experimentNames[i]);
		}
		
		for (i=0;i<mouseMotifData.numCols;i++) {
			mouseMotifData.experimentNames[i] = motifMap.get(mouseMotifData.experimentNames[i]);
		}
		
		try {
			outputMotifFile(humanMotifOutFileName,humanMotifData.values,geneNames,humanGeneMap,humanMotifData.experimentNames);
			outputMotifFile(mouseMotifOutFileName,mouseMotifData.values,geneNames,mouseGeneMap,mouseMotifData.experimentNames);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void outputMotifFile(String fOutName,double[][] values,String[] geneNames,HashMap<String,Integer> geneMap,String[] experimentNames) throws IOException {
		FileWriter file = new FileWriter(fOutName);
		int i = 0;
		int j = 0;
		
		file.write("genes");
		for (i=0;i<experimentNames.length;i++) {
			file.write("\t");
			file.write(experimentNames[i]);
		}
		file.write("\n");
		
		int gg = 0;
		
		for (i=0;i<geneNames.length;i++) {
			if (geneMap.containsKey(geneNames[i])) {
				gg = geneMap.get(geneNames[i]);
				file.write(geneNames[i]);
				for (j=0;j<values[0].length;j++) {
					file.write("\t");
					file.write((new Double(values[gg][j])).toString());
				}
				file.write("\n");
			}
		}
		
		file.close();
	}
	
	public static String[] getGenes(String topGenesFileName) throws IOException {
		HashSet<String> genes = new HashSet<String>();
		int i = 0;
		String[] lineSplit = null;
		BufferedReader is = new BufferedReader(new FileReader(topGenesFileName));
		String line = "";		
		line = is.readLine();
		while(line != null) {
			lineSplit = line.split("\t");
			
			if (lineSplit.length > 1) {
				for (i=1;i<lineSplit.length;i++) {
					genes.add(lineSplit[i]);
				}
			}
			line = is.readLine();
		}
		is.close();
		
		String[] geneNames = new String[genes.size()];
		geneNames = genes.toArray(geneNames);
		
		return geneNames;
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
	
	public static void thresholdMotifs(double[][] values) {
		Vector<Double> mv = new Vector<Double>();
		double[] maxV = new double[values[0].length];
		double[] rankV = new double[values[0].length];
		int i = 0;
		int j = 0;
		
		double percentile = 0.75;
		double percentMax = 0.60;
	//	double percentMax = 0.70;
		for (j=0;j<values[0].length;j++) {
			mv.clear();
			for (i=0;i<values.length;i++) {
				mv.add(values[i][j]);
			}
			maxV[j] = percentMax*VectorUtil.max(mv);
			rankV[j] = VectorUtil.percentile(mv,percentile);
		}
		
		for (i=0;i<values.length;i++) {
			for (j=0;j<values[0].length;j++) {
		//		if (values[i][j] < maxV[j] | values[i][j] < rankV[j]) {
				if (values[i][j] < maxV[j]) {
					values[i][j] = 0.0/0.0;
				}
			}
		}
	}
	
	public static void mergeMotifs() {
		String combinedMotifOutFileName = "C:\\research_data\\mouse_human\\motifs\\combined_cons_discretized.txt";
		String mergedOutFileName = "C:\\research_data\\mouse_human\\motifs\\combined_cons_merged.txt";
		double presentThresh = 0.25;
		
		MicroArrayData motifData = new MicroArrayData();
		motifData.setDiscrete();
		try {
			motifData.readFile(combinedMotifOutFileName);
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
			motifData.writeFile(mergedOutFileName);
		} catch(IOException e) {
			
		}
	}
}
