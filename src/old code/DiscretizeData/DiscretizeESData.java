package edu.mit.csail.psrg.georg.DiscretizeData;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

public class DiscretizeESData {
	public static void main(String[] args) {
		discretizeData();
		mergeReplicates();
	}
	
	public static void discretizeData() {
		String dirPath = "c:\\research_data\\mouse_human\\ESNovartis\\";
		String humanExpressionName = dirPath + "ESNovartis.geo.txt";
		String humanPMAName = dirPath + "ESNovartis.pma.txt";
		String mouseExpressionName = dirPath + "GNF1M.geo.txt";
		String mousePMAName = dirPath + "GNF1M.pma.txt";
		
		String combinedName = dirPath + "ESNovartis_reps_combined_discretized.txt";
		String levelsFileName = dirPath + "discretized_levels.txt";
		String MIFileName = dirPath + "discretized_MI.txt";
		
		MicroArrayData humanExpression = new MicroArrayData();
		MicroArrayData mouseExpression = new MicroArrayData();
		
		MicroArrayData humanPMA = new MicroArrayData();
		humanPMA.setDiscrete();
		MicroArrayData mousePMA = new MicroArrayData();
		mousePMA.setDiscrete();
		
		try {
			humanExpression.readFileBlindContinuous(humanExpressionName);
			mouseExpression.readFileBlindContinuous(mouseExpressionName);
			humanPMA.readFileBlindDiscrete(humanPMAName);
			mousePMA.readFileBlindDiscrete(mousePMAName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		int j = 0;
		
		humanExpression.experimentCode = new int[humanExpression.numCols];
		for (j=0;j<humanExpression.numCols;j++) {
			humanExpression.experimentNames[j] = "h_" + humanExpression.experimentNames[j];
			humanExpression.experimentCode[j] = 1;
		}
		for (i=0;i<humanExpression.values.length;i++) {
			for (j=0;j<humanExpression.values[0].length;j++) {
				if (humanPMA.dvalues[i][j] == 0 | humanExpression.values[i][j] < 2.0) {
					humanExpression.values[i][j] = 0.0/0.0;
				}
			}
		}
		
		mouseExpression.experimentCode = new int[mouseExpression.numCols];
		for (j=0;j<mouseExpression.numCols;j++) {
			mouseExpression.experimentNames[j] = "m_" + mouseExpression.experimentNames[j];
			mouseExpression.experimentCode[j] = 2;
		}
		for (i=0;i<mouseExpression.values.length;i++) {
			for (j=0;j<mouseExpression.values[0].length;j++) {
				if (mousePMA.dvalues[i][j] == 0 | mouseExpression.values[i][j] < 2.0) {
					mouseExpression.values[i][j] = 0.0/0.0;
				}
			}
		}
		
		String[] combinedGenes = new String[humanExpression.numRows];
		for (i=0;i<humanExpression.numRows;i++) {
			combinedGenes[i] = humanExpression.geneNames[i] + "|" + mouseExpression.geneNames[i];
		}
		humanExpression.geneNames = combinedGenes;
		mouseExpression.geneNames = combinedGenes;
		MicroArrayData mergedData = MicroArrayData.mergeCols(humanExpression,mouseExpression);
		
		humanExpression = null;
		mouseExpression = null;
		
		DiscretizeExpression discretizer = new DiscretizeExpression(mergedData.values,10);
		
		try {
			discretizer.discretizeDownTree(MIFileName);
			discretizer.mergeDownLevels(3,levelsFileName);
			discretizer.transposeDExpression();
			mergedData.setDiscrete();
			mergedData.dvalues = discretizer.dExpression;
			mergedData.writeFile(combinedName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void mergeReplicates() {
		String dirPath = "c:\\research_data\\mouse_human\\ESNovartis\\";
		String combinedName = dirPath + "ESNovartis_reps_combined_discretized.txt";
		String mergedName = dirPath + "ESNovartis_combined_discretized.txt";
		
		MicroArrayData combinedData = new MicroArrayData();
		combinedData.setDiscrete();
		
		try {
			combinedData.readFile(combinedName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		LinkedHashMap<String,ArrayList<Integer>> tissueMap = new LinkedHashMap<String,ArrayList<Integer>>();
		ArrayList<Integer> tissueList = null;
		
		int i = 0;
		int j = 0;
		int jj = 0;
		int k = 0;
		for (i=0;i<combinedData.numCols;i++) {
			tissueList = tissueMap.get(combinedData.experimentNames[i]);
			if (tissueList == null) {
				tissueList = new ArrayList<Integer>();
				tissueMap.put(combinedData.experimentNames[i],tissueList);
			}
			tissueList.add(i);
		}
		
		int[][] values = new int[combinedData.numRows][tissueMap.size()];
		int[] experimentCode = new int[tissueMap.size()];
		String[] tissueNames = new String[tissueMap.size()];
		Iterator<String> nameIter = tissueMap.keySet().iterator();
		
		int maxValue = 0;
		j = 0;
		while(nameIter.hasNext()) {
			tissueNames[j] = nameIter.next();
			tissueList = tissueMap.get(tissueNames[j]);
			for (i=0;i<combinedData.numRows;i++) {
				maxValue = 0;
				for (k=0;k<tissueList.size();k++) {
					jj = tissueList.get(k);
					if (combinedData.dvalues[i][jj] > maxValue) {
						maxValue = combinedData.dvalues[i][jj];
					}
					experimentCode[j] = combinedData.experimentCode[jj];
				}
				values[i][j] = maxValue;
			}
			j++;
		}
		
		combinedData.dvalues = values;
		combinedData.numCols = tissueNames.length;
		combinedData.experimentCode = experimentCode;
		combinedData.experimentNames = tissueNames;
		
		try {
			combinedData.writeFile(mergedName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
}
