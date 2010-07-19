package edu.mit.csail.psrg.georg.DiscretizeData;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

public class DiscretizeMouseHuman {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	//	processData();
	//	resortMotifData();
		processSingleSpecies();
	}
	
	public static void processSingleSpecies() {
		int speciesCode = 2;
		String speciesName = "";
		String dirName = "C:\\research_data\\mouse_human\\b47mm7\\geo\\mouse_only_all\\";
		String EName = dirName + "GNF1M.gcRMA.geo.ugene.txt";
		String PMAName = dirName + "GNF1M.pma.ugene.txt";
		String combinedName = dirName + "combined_discretized.txt";
		String levelsFileName = dirName + "discretized_levels.txt";
		String MIFileName = dirName + "discretized_MI.txt";
		double expressionThresh = 2.0;
	
		MicroArrayData PMA = combinePMAReplicates(PMAName);
		int i = 0;
		int j = 0;
		
		if (speciesCode == 1) {
			speciesName = "h_";
		} else {
			speciesName = "m_";
		}
		
		for (i=0;i<PMA.numCols;i++) {
			PMA.experimentNames[i] = speciesName + PMA.experimentNames[i];
		}
		
		MicroArrayData expression = new MicroArrayData();
				
		try {
			expression.readFileBlind(EName);
			
			expression.experimentCode = new int[expression.numCols];
			for (i=0;i<expression.numCols;i++) {
				expression.experimentNames[i] = speciesName + expression.experimentNames[i];
				expression.experimentCode[i] = speciesCode;
			}
			
			for (i=0;i<expression.numRows;i++) {
				if (speciesCode == 1) {
					expression.geneNames[i] = expression.geneNames[i] + "|" + "NM_none";
				} else {
					expression.geneNames[i] = "NM_none" + "|" + expression.geneNames[i];
				}
			}
		} catch(IOException e) {
			System.out.println(e);
		}
		
		for (i=0;i<expression.numRows;i++) {
			for (j=0;j<expression.numCols;j++) {
				if (expression.values[i][j] < expressionThresh | PMA.dvalues[i][j] == 0) {
					expression.values[i][j] = 0.0/0.0;
				}
			}
		}
		
		DiscretizeExpression discretizer = new DiscretizeExpression(expression.values,10);
		
		try {
			discretizer.discretizeDownTree(MIFileName);
			discretizer.mergeDownLevels(3,levelsFileName);
			discretizer.transposeDExpression();
			expression.setDiscrete();
			expression.dvalues = discretizer.dExpression;
			expression.writeFile(combinedName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
//	 combine replicate PMA values (assumes two replicates, adjacent columns)
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
	
	public static void processData() {
	//	String dirName = "C:\\research_data\\mouse_human\\";
		String dirName = "C:\\research_data\\mouse_human\\homo_b47_data\\";
	//	String GNF1H_fName = "GNF1H_genes_chopped";
	//	String GNF1M_fName = "GNF1M_genes_chopped";
		String GNF1H_fName = "GNF1H";
		String GNF1M_fName = "GNF1M";
		String mergedValues_fName = "GNF1_merged_expression.txt";
		String mergedPMA_fName = "GNF1_merged_PMA.txt";
		String discretizedMI_fName = "GNF1_discretized_MI.txt";
		String discretizedLevels_fName = "GNF1_discretized_levels.txt";
		String discretizedData_fName = "GNF1_discretized_data.txt";
		
		boolean useMotifs = false;
		MicroArrayData humanMotifs = null;
		MicroArrayData mouseMotifs = null;
		
		MicroArrayData GNF1H_value = new MicroArrayData();
		MicroArrayData GNF1H_PMA = new MicroArrayData();
		GNF1H_PMA.setDiscrete();
		MicroArrayData GNF1M_value = new MicroArrayData();
		MicroArrayData GNF1M_PMA = new MicroArrayData();
		GNF1M_PMA.setDiscrete();
		
		try {
		/*	GNF1H_value.readFile(dirName+GNF1H_fName+".gcrma.txt");
			GNF1H_PMA.readFile(dirName+GNF1H_fName+".pma");
			GNF1M_value.readFile(dirName+GNF1M_fName+".gcrma.txt");
			GNF1M_PMA.readFile(dirName+GNF1M_fName+".pma"); */
		/*	GNF1H_value.readFile(dirName+GNF1H_fName+".gcrma.txt.homob44");
			GNF1H_PMA.readFile(dirName+GNF1H_fName+".pma.txt.homob44");
			GNF1M_value.readFile(dirName+GNF1M_fName+".gcrma.txt.homob44");
			GNF1M_PMA.readFile(dirName+GNF1M_fName+".pma.txt.homob44"); */
			GNF1H_value.readFile(dirName+GNF1H_fName+".gcrma.snco.txt");
			GNF1H_PMA.readFile(dirName+GNF1H_fName+".pma.sn.txt");
			GNF1M_value.readFile(dirName+GNF1M_fName+".gcrma.snco.txt");
			GNF1M_PMA.readFile(dirName+GNF1M_fName+".pma.sn.txt");
		} catch(IOException e) {
			System.out.println(e);
		}
		
		if (useMotifs) {
			humanMotifs = new MicroArrayData();
			mouseMotifs = new MicroArrayData();
			loadMotifs(humanMotifs,GNF1H_value.geneNames,mouseMotifs,GNF1M_value.geneNames);
		}
		
		// combine replicates in PMA values
		int[][] H_PMA = new int[GNF1H_PMA.numRows][GNF1H_PMA.numCols/2];
		int[][] M_PMA = new int[GNF1M_PMA.numRows][GNF1M_PMA.numCols/2];
		
		int i = 0;
		int j = 0;
		int k = 0;
		int v = 0;
		k = 0;
		j = 0;
		while (j<GNF1H_PMA.numCols-1) {
			for (i=0;i<H_PMA.length;i++) {
				v = 0;
			//	if (GNF1H_PMA.dvalues[i][j] > 0 & GNF1H_PMA.dvalues[i][j+1] > 0) {
				if (GNF1H_PMA.dvalues[i][j] > 0 | GNF1H_PMA.dvalues[i][j+1] > 0) {
					v = 1;
				}
				H_PMA[i][k] = v;
			}
			k++;
			j = j + 2;
		}
		
		j = 0;
		k = 0;
		while (j<GNF1M_PMA.numCols-1) {
			for (i=0;i<M_PMA.length;i++) {
				v = 0;
			//	if (GNF1M_PMA.dvalues[i][j] > 0 & GNF1M_PMA.dvalues[i][j+1] > 0) {
				if (GNF1M_PMA.dvalues[i][j] > 0 | GNF1M_PMA.dvalues[i][j+1] > 0) {
					v = 1;
				}
				M_PMA[i][k] = v;
			}
			k++;
			j = j + 2;
		}
		
		int numMatched = 31;
		
	/*	GNF1H_value.numCols = numMatched;
		GNF1M_value.numCols = numMatched; */
		
		int[][] matchPairs = new int[numMatched][2];
		for(i=0;i<numMatched;i++) {
			matchPairs[i][0] = i;
			matchPairs[i][1] = i + GNF1H_value.numCols;
		}
		
	//	DiscretizeAffyData H_discretizer = new DiscretizeAffyData(GNF1H_value.values,H_PMA,0);
	//	H_discretizer.adjustIntensities();
		transformRelativeAbundance(GNF1H_value.values,H_PMA,GNF1H_value.numCols);
		
	//	DiscretizeAffyData M_discretizer = new DiscretizeAffyData(GNF1M_value.values,M_PMA,0);
	//	M_discretizer.adjustIntensities();
		transformRelativeAbundance(GNF1M_value.values,M_PMA,GNF1M_value.numCols);
		
		double[][] mergedExpression = new double[GNF1H_value.numRows][GNF1H_value.numCols+GNF1M_value.numCols];
		int[][] mergedPMA = new int[GNF1H_value.numRows][GNF1H_value.numCols+GNF1M_value.numCols];
		int[] species = null;
		if (!useMotifs) {
			species = new int[GNF1H_value.numCols+GNF1M_value.numCols];
		} else {
			species = new int[GNF1H_value.numCols+GNF1M_value.numCols + humanMotifs.numCols + mouseMotifs.numCols];
		}
		
		for (i=0;i<GNF1H_value.numRows;i++) {
			System.arraycopy(GNF1H_value.values[i],0,mergedExpression[i],0,GNF1H_value.numCols);
			System.arraycopy(H_PMA[i],0,mergedPMA[i],0,GNF1H_value.numCols);	
		}
		for (i=0;i<GNF1M_value.numRows;i++) {
			System.arraycopy(GNF1M_value.values[i],0,mergedExpression[i],GNF1H_value.numCols,GNF1M_value.numCols);
			System.arraycopy(M_PMA[i],0,mergedPMA[i],GNF1H_value.numCols,GNF1M_value.numCols);
		}
		
		// concatenate experiment and gene names
		for (i=0;i<GNF1H_value.numCols;i++) {
			GNF1H_value.experimentNames[i] = "h_" + GNF1H_value.experimentNames[i];
			species[i] = 1;
		}
		for (i=0;i<GNF1M_value.numCols;i++) {
			GNF1M_value.experimentNames[i] = "m_" + GNF1M_value.experimentNames[i];
			species[i+GNF1H_value.numCols] = 2;
		}
		
		String[] mergedExperimentNames = null;
		if (!useMotifs) {
			mergedExperimentNames = new String[GNF1H_value.numCols+GNF1M_value.numCols];
		} else {
			mergedExperimentNames = new String[GNF1H_value.numCols+GNF1M_value.numCols + humanMotifs.numCols + mouseMotifs.numCols];
		}
		System.arraycopy(GNF1H_value.experimentNames,0,mergedExperimentNames,0,GNF1H_value.numCols);
		System.arraycopy(GNF1M_value.experimentNames,0,mergedExperimentNames,GNF1H_value.numCols,GNF1M_value.numCols);
		if (useMotifs) {
			System.arraycopy(humanMotifs.experimentNames,0,mergedExperimentNames,GNF1H_value.numCols+GNF1M_value.numCols,humanMotifs.numCols);
			System.arraycopy(mouseMotifs.experimentNames,0,mergedExperimentNames,GNF1H_value.numCols+GNF1M_value.numCols+humanMotifs.numCols,mouseMotifs.numCols);
			for (i=0;i<humanMotifs.numCols;i++) {
				species[i + GNF1H_value.numCols+GNF1M_value.numCols] = 3;
			}
			for (i=0;i<mouseMotifs.numCols;i++) {
				species[i + GNF1H_value.numCols+GNF1M_value.numCols + humanMotifs.numCols] = 4;
			}
		}
		
		String[] mergedGeneNames = new String[GNF1H_value.numRows];
		for (i=0;i<GNF1M_value.numRows;i++) {
			mergedGeneNames[i] = GNF1H_value.geneNames[i] + "|" + GNF1M_value.geneNames[i];
		}
		
		int[] filterList = new int[GNF1M_value.numRows];
		int numFiltered = 0;
		double maxExpressedPercent = 1.25;
		int maxExpressed_H = (int) Math.floor(maxExpressedPercent*((double) GNF1H_value.numCols));
		int maxExpressed_M = (int) Math.floor(maxExpressedPercent*((double) GNF1M_value.numCols));
		int numExpressed_H = 0;
		int numExpressed_M = 0;
		
		for (i=0;i<GNF1H_value.numRows;i++) {
			numExpressed_H = 0;
			for (j=0;j<GNF1H_value.numCols;j++) {
				if (GNF1H_PMA.dvalues[i][j] > 0) {
			//	if (!Double.isNaN(GNF1H_value.values[i][j])) {
			//	if (GNF1H_value.values[i][j] > 1.5) {
					numExpressed_H++;
				}
			}
			numExpressed_M = 0;
			for (j=0;j<GNF1M_value.numCols;j++) {
				if (GNF1M_PMA.dvalues[i][j] > 0) {
			//	if (GNF1H_value.values[i][j] > 1.5) {
			//	if (!Double.isNaN(GNF1M_value.values[i][j])) {
					numExpressed_M++;
				}
			}
			
			if (numExpressed_M >= maxExpressed_M | numExpressed_H >= maxExpressed_H) {
				filterList[i] = 1;
				numFiltered++;
			}
		}
		
		GNF1H_value = null;
		GNF1H_PMA = null;
		GNF1M_value = null;
		GNF1M_PMA = null;
		
		double[][] mergedExpression2 = null;
		if (numFiltered > 0) {
			k = 0;
			int[][] mergedPMA2 = new int[mergedPMA.length-numFiltered][mergedPMA[0].length];
			if (!useMotifs) {
				mergedExpression2 = new double[mergedExpression.length-numFiltered][mergedExpression[0].length];
			} else {
				mergedExpression2 = new double[mergedExpression.length-numFiltered][mergedExpression[0].length + humanMotifs.numCols + mouseMotifs.numCols];
			}
			String[] mergedGeneNames2 = new String[mergedGeneNames.length-numFiltered];
			for (i=0;i<filterList.length;i++) {
				if (filterList[i] == 0) {
					mergedPMA2[k] = mergedPMA[i];
					if (!useMotifs) {
						mergedExpression2[k] = mergedExpression[i];
					} else {
						System.arraycopy(mergedExpression[i],0,mergedExpression2[k],0,mergedExpression[i].length);
						System.arraycopy(humanMotifs.values[i],0,mergedExpression2[k],mergedExpression[i].length,humanMotifs.numCols);
						System.arraycopy(mouseMotifs.values[i],0,mergedExpression2[k],humanMotifs.numCols+mergedExpression[i].length,mouseMotifs.numCols);
					}
					mergedGeneNames2[k] = mergedGeneNames[i];
					k++;
				}
			}
			mergedPMA = mergedPMA2;
			mergedExpression = mergedExpression2;
			mergedGeneNames = mergedGeneNames2;
		} else {
			if (useMotifs) {
				mergedExpression2 = new double[mergedExpression.length][mergedExpression[0].length + humanMotifs.numCols + mouseMotifs.numCols];
				for (i=0;i<mergedExpression.length;i++) {
					System.arraycopy(mergedExpression[i],0,mergedExpression2[i],0,mergedExpression[i].length);
					System.arraycopy(humanMotifs.values[i],0,mergedExpression2[i],mergedExpression[i].length,humanMotifs.numCols);
					System.arraycopy(mouseMotifs.values[i],0,mergedExpression2[i],humanMotifs.numCols+mergedExpression[i].length,mouseMotifs.numCols);
				}
				mergedExpression = mergedExpression2;
			}
		}
		
		DiscretizeAffyPairedData discretizer = new DiscretizeAffyPairedData(mergedExpression,mergedPMA,matchPairs,20);
		MicroArrayData mergedData = new MicroArrayData();
		mergedData.values = discretizer.expression;
		mergedData.geneNames = mergedGeneNames;
		mergedData.experimentNames = mergedExperimentNames;
		mergedData.experimentCode = species;
		mergedData.numCols = discretizer.numExperiments;
		mergedData.numRows = discretizer.numGenes;
		
		try {
			mergedData.writeFile(dirName+mergedValues_fName);
			discretizer.discretizeDownTree(dirName+discretizedMI_fName);
			discretizer.mergeDownLevels(3,dirName+discretizedLevels_fName);
			discretizer.transposeDExpression();
			mergedData.setDiscrete();
			mergedData.values = null;
			mergedData.dvalues = discretizer.dExpression;
			mergedData.writeFile(dirName+discretizedData_fName);
	
			mergedData.dvalues = mergedPMA;
			mergedData.numCols = mergedPMA[0].length;
			mergedData.writeFile(dirName+mergedPMA_fName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
	/*	MicroArrayData mergedData = new MicroArrayData();
		mergedData.values = mergedExpression;
		mergedData.geneNames = mergedGeneNames;
		mergedData.experimentNames = mergedExperimentNames;
		mergedData.experimentCode = species;
		mergedData.numCols = mergedExpression[0].length;
		mergedData.numRows = mergedExpression.length;
		try {
			mergedData.writeFile(dirName+mergedValues_fName);
			mergedData.setDiscrete();
			mergedData.values = null;
		//	mergedData.dvalues = simpleDiscretization(mergedExpression,mergedPMA);
			mergedData.dvalues = simpleDiscretization2(mergedExpression,species);
			mergedData.writeFile(dirName+discretizedData_fName);
	
			mergedData.dvalues = mergedPMA;
			mergedData.numCols = mergedPMA[0].length;
			mergedData.writeFile(dirName+mergedPMA_fName);
		} catch(IOException e) {
			System.out.println(e);
		} */
		
	}
	
	public static void transformRelativeAbundance(double[][] expression,int[][] PMA,int numPairs) {
		double norm = 0.0;
		int i = 0;
		int j = 0;
		
		for (i=0;i<expression.length;i++) {
			norm = 0.0;
			for (j=0;j<numPairs;j++) {
				norm = norm + expression[i][j];
			}
			for (j=0;j<expression[0].length;j++) {
			//	expression[i][j] = expression[i][j]/norm;
			//	expression[i][j] = expression[i][j]*((double) numPairs);
				expression[i][j] = expression[i][j]*31.0;
				//	expression[i][j] = Math.log(expression[i][j]);
			//	if (PMA[i][j] == 0 | expression[i][j] <= 1.0/((double)numPairs)) {
				if (PMA[i][j] == 0 | expression[i][j] <= 1.5) {
			//	if (PMA[i][j] == 0 | expression[i][j] <= 2.0) {
					expression[i][j] = 0.0/0.0;
				} 
			} 
		}
	}
	
	public static int[][] simpleDiscretization2(double[][] expression,int[] experimentCode) {
		int[][] DExpression = new int[expression.length][expression[0].length];
		int i = 0;
		int j = 0;
		double v = 0.0;
		double minFoldChange = 1.5;
		
		for (i=0;i<DExpression.length;i++) {
			for (j=0;j<DExpression[0].length;j++) {
				v = expression[i][j];
				if (!Double.isNaN(v) & v > minFoldChange & experimentCode[j] <= 2) {
					v = Math.log(v/minFoldChange);
					DExpression[i][j] = ((int) Math.round(v)) + 1;
				} else {
					if (experimentCode[j] > 2 & !Double.isNaN(v) & v > 12.0) {
						DExpression[i][j] = 2;
					}
				}
			}
		}
		
		return DExpression;
	}
	
	public static int[][] simpleDiscretization(double[][] expression,int[][] PMA) {
		int[][] DExpression = new int[expression.length][expression[0].length];
		int i = 0;
		int j = 0;
		double v = 0.0;
		double minFoldChange = 1.5;
		
		for (i=0;i<DExpression.length;i++) {
			for (j=0;j<DExpression[0].length;j++) {
				v = expression[i][j];
				if (PMA[i][j] == 1 & v > minFoldChange) {
				//	v = Math.log(v/minFoldChange)/Math.log(minFoldChange);
					v = Math.log(v/minFoldChange);
					DExpression[i][j] = ((int) Math.round(v)) + 1;
				}
			/*	if (PMA[i][j] == 1 & v < 1/minFoldChange) {
					v = Math.log(v/minFoldChange)/Math.log(minFoldChange);
					DExpression[i][j] = ((int) Math.round(v)) + 1;
				} */
			}
		}
		
		return DExpression;
	}
	
	public static void loadMotifs(MicroArrayData humanMotifs,String[] humanGeneNames,MicroArrayData mouseMotifs,String[] mouseGeneNames) {
		String dirName = "C:\\research_data\\mouse_human\\homo_b47_data\\";
		String humanMotifDName = "human_5kb.out";
		String humanMotifVName = "human_5kb_raw.out";
		String mouseMotifDName = "mouse_5kb.out";
		String mouseMotifVName = "mouse_5kb_raw.out";
		
		MicroArrayData humanMotifsD = new MicroArrayData();
		humanMotifsD.setDiscrete();
		
		MicroArrayData mouseMotifsD = new MicroArrayData();
		mouseMotifsD.setDiscrete();
		
		try {
			humanMotifs.readFile(dirName + humanMotifVName);
			humanMotifsD.readFile(dirName + humanMotifDName);
			mouseMotifs.readFile(dirName + mouseMotifVName);
			mouseMotifsD.readFile(dirName + mouseMotifDName);
		} catch (IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		int j = 0;
		
		HashMap<String,Integer> humanGeneMap = new HashMap<String,Integer>();
		for (i=0;i<humanGeneNames.length;i++) {
			humanGeneMap.put(humanGeneNames[i],i);
		}
		
		HashMap<String,Integer> mouseGeneMap = new HashMap<String,Integer>();
		for (i=0;i<mouseGeneNames.length;i++) {
			mouseGeneMap.put(mouseGeneNames[i],i);
		}
		
		for (j=0;j<humanMotifs.numCols;j++) {
			humanMotifs.experimentNames[j] = "h_mot_" + humanMotifs.experimentNames[j];
			for (i=0;i<humanMotifs.numRows;i++) {
				if (humanMotifsD.dvalues[i][j] == 0) {
					humanMotifs.values[i][j] = 0.0/0.0;
				}
			}
		}
		
		for (j=0;j<mouseMotifs.numCols;j++) {
			mouseMotifs.experimentNames[j] = "m_mot_" + mouseMotifs.experimentNames[j];
			for (i=0;i<mouseMotifs.numRows;i++) {
				if (mouseMotifsD.dvalues[i][j] == 0) {
					mouseMotifs.values[i][j] = 0.0/0.0;
				}
			}
		}
		
		double maxPercent = 0.25;
		int maxPresent = (int) Math.round(maxPercent*humanMotifs.numRows);
		int minPresent = 99;
		ArrayList<Integer> retainHuman = new ArrayList<Integer>();
		ArrayList<Integer> retainMouse = new ArrayList<Integer>();
		
		int numP_human = 0;
		int numP_mouse = 0;
		int numHuman = 0;
		int numMouse = 0;
		for (j=0;j<humanMotifs.numCols;j++) {
			numP_human = 0;
			for (i=0;i<humanMotifs.numRows;i++) {
				if (humanMotifsD.dvalues[i][j] == 1)
					numP_human++;
			}
			if (numP_human <= maxPresent & numP_human > minPresent) {
				retainHuman.add(j);
				numHuman++;
			}
		}
		
		for (j=0;j<mouseMotifs.numCols;j++) {
			numP_mouse = 0;
			for (i=0;i<mouseMotifs.numRows;i++) {
				if (mouseMotifsD.dvalues[i][j] == 1)
					numP_mouse++;
			}
			if (numP_mouse <= maxPresent & numP_mouse > minPresent) {
				retainMouse.add(j);
				numMouse++;
			}
		}
		
		double[][] humanV = new double[humanMotifs.numRows][numHuman];
		String[] humanE = new String[numHuman];
		double[][] mouseV = new double[humanMotifs.numRows][numMouse];
		String[] mouseE = new String[numMouse];
		
		int k = 0;
		int gg = 0;
		
		for (k=0;k<retainHuman.size();k++) {
			j = retainHuman.get(k);
			humanE[k] = humanMotifs.experimentNames[j];
			for (i=0;i<humanMotifs.numRows;i++) {
				gg = humanGeneMap.get(humanMotifs.geneNames[i]);
				humanV[gg][k] = humanMotifs.values[i][j];
			}
		}
		humanMotifs.values = humanV;
		humanMotifs.experimentNames = humanE;
		humanMotifs.numCols = humanV[0].length;
		
		for (k=0;k<retainMouse.size();k++) {
			j = retainMouse.get(k);
			mouseE[k] = mouseMotifs.experimentNames[j];
			for (i=0;i<mouseMotifs.numRows;i++) {
				gg = mouseGeneMap.get(mouseMotifs.geneNames[i]);
				mouseV[gg][k] = mouseMotifs.values[i][j];
			}
		}
		mouseMotifs.values = mouseV;
		mouseMotifs.experimentNames = mouseE;
		mouseMotifs.numCols = mouseV[0].length;
	}

	static void resortMotifData() {
		String dirName = "C:\\research_data\\mouse_human\\homo_b47_data\\";
		String humanMotifDName = "human_5kb.out";
		String mouseMotifDName = "mouse_5kb.out";
		String humanMotifDNameOut = "human_5kb_conformed.out";
		String mouseMotifDNameOut = "mouse_5kb_conformed.out";
		String humanExpressionName = "GNF1H.gcrma.snco.txt";
		String mouseExpressionName = "GNF1M.gcrma.snco.txt";
		
		MicroArrayData humanMotifsD = new MicroArrayData();
		humanMotifsD.setDiscrete();
		
		MicroArrayData mouseMotifsD = new MicroArrayData();
		mouseMotifsD.setDiscrete();
		
		MicroArrayData humanExpression = new MicroArrayData();
		MicroArrayData mouseExpression = new MicroArrayData();
		String[] humanGeneNames = null;
		String[] mouseGeneNames = null;
		
		try {
			humanExpression.readFile(dirName + humanExpressionName);
			humanGeneNames = humanExpression.geneNames;
			humanExpression = null;
			
			mouseExpression.readFile(dirName + mouseExpressionName);
			mouseGeneNames = mouseExpression.geneNames;
			mouseExpression = null;
			
			humanMotifsD.readFile(dirName + humanMotifDName);
			mouseMotifsD.readFile(dirName + mouseMotifDName);
		} catch (IOException e) {
			System.out.println(e);
		}
		
		HashMap<String,String> motifNameMap = readMotifNames();
		
		LinkedHashSet<String> unionExperimentNames = new LinkedHashSet<String>();
		int i = 0;
		int j = 0;
		
		for (j=0;j<mouseMotifsD.numCols;j++) {
			unionExperimentNames.add(mouseMotifsD.experimentNames[j]);
		}
		for (j=0;j<humanMotifsD.numCols;j++) {
			unionExperimentNames.add(humanMotifsD.experimentNames[j]);
		}
		
		HashMap<String,Integer> experimentOrder = new HashMap<String,Integer>();
		Iterator<String> iter = unionExperimentNames.iterator();
		j = 0;
		while(iter.hasNext()) {
			experimentOrder.put(iter.next(),j);
			j++;
		}
		
		HashMap<String,Integer> humanGeneMap = new HashMap<String,Integer>();
		HashMap<String,Integer> mouseGeneMap = new HashMap<String,Integer>();
		
		for (i=0;i<mouseGeneNames.length;i++) {
			mouseGeneMap.put(mouseGeneNames[i],i);
		}
		
		for (i=0;i<humanGeneNames.length;i++) {
			humanGeneMap.put(humanGeneNames[i],i);
		}
		
		int[][] humanMotifs = new int[humanGeneNames.length][experimentOrder.size()];
		int[][] mouseMotifs = new int[mouseGeneNames.length][experimentOrder.size()];
		
		int ii = 0;
		int jj = 0;
		for (ii=0;ii<mouseGeneNames.length;ii++) {
			i = mouseGeneMap.get(mouseMotifsD.geneNames[ii]);
			for (jj=0;jj<mouseMotifsD.numCols;jj++) {
				if (experimentOrder.containsKey(mouseMotifsD.experimentNames[jj])) {
					j = experimentOrder.get(mouseMotifsD.experimentNames[jj]);
					mouseMotifs[i][j] = mouseMotifsD.dvalues[ii][jj];
				}
			}
		}
		
		String[] experimentNames = new String[experimentOrder.size()];
		experimentNames = experimentOrder.keySet().toArray(experimentNames);
		for (j=0;j<experimentNames.length;j++) {
			if (motifNameMap.containsKey(experimentNames[j])) {
			//	experimentNames[j] = motifNameMap.get(experimentNames[j]);
			}
		}
		
		mouseMotifsD.experimentNames = experimentNames;
		mouseMotifsD.geneNames = mouseGeneNames;
		mouseMotifsD.dvalues = mouseMotifs;
		mouseMotifsD.numCols = experimentNames.length;
		
		for (ii=0;ii<humanGeneNames.length;ii++) {
			i = humanGeneMap.get(humanMotifsD.geneNames[ii]);
			for (jj=0;jj<humanMotifsD.numCols;jj++) {
				if (experimentOrder.containsKey(humanMotifsD.experimentNames[jj])) {
					j = experimentOrder.get(humanMotifsD.experimentNames[jj]);
					humanMotifs[i][j] = humanMotifsD.dvalues[ii][jj];
				}
			}
		}
		
		humanMotifsD.experimentNames = experimentNames;
		humanMotifsD.geneNames = humanGeneNames;
		humanMotifsD.dvalues = humanMotifs;
		humanMotifsD.numCols = experimentNames.length;
		
		try {
			humanMotifsD.writeFile(dirName + humanMotifDNameOut);
			mouseMotifsD.writeFile(dirName + mouseMotifDNameOut);
		} catch (IOException e) {
			System.out.println(e);
		}
	}
	
	public static HashMap<String,String> readMotifNames() {
		String dirName = "C:\\research_data\\mouse_human\\homo_b47_data\\";
		String motifNamesFile = "tf_ids_and_names.txt";
		String line = null;
		String[] info = null;
		String tfName = null;
		int f = 0;
		HashMap<String,String> motifNames = new HashMap<String,String>();
		
		try {
			BufferedReader is = new BufferedReader(new FileReader(dirName + motifNamesFile));
			line = is.readLine();
			while (line != null) {
				info = line.split("\t");
				f = info[1].indexOf(")");
				if (f > 0) {
					if (f+2 < info[1].length() - 1) {
						tfName = info[1].substring(f+2);
					} else {
						tfName = info[1].substring(f+1);
					}
				} else {
					tfName = info[1];
				}
				motifNames.put(info[0],tfName);
				line = is.readLine();
			}
		} catch (IOException e) {
			System.out.println(e);
		}
		return motifNames;
	}
}
