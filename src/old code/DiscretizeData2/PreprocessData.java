package edu.mit.csail.psrg.georg.DiscretizeData2;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;
import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

// Do simple filtering on data: remove genes w/ insufficient
// fold-change, etc.
// also, optionally log transform the data and
// thresh-hold values
public class PreprocessData {
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		doPreprocess();
	}
	
	public static void doPreprocess() {
		String outputDir = "C:\\research_data\\zon\\";
		
		String inFName = outputDir + "zfish_gcrma.txt";
		String outFName = outputDir + "zfish_gcrma_filtered.txt";
		String outputFile = outputDir + "zfish_gcrma_discrete.txt";
		String selectFName = null;
		
		// threshold for considering a gene up or down
		double threshold = 0.5850;
	//	double logBase = 2.0;
		double logBase = 0.0;
		filterFile(inFName,outFName,selectFName,threshold,logBase);
		
		doDiscretization(outputDir,outFName,outputFile);
	}
	
	public static void doDiscretization(String outputDir,String inputFile,String outputFile) {
		String combinedName = outputDir + "combined_discretized.txt";
		String levelsFileName = outputDir + "discretized_levels";
		String MIFileName = outputDir + "discretized_MI.txt";
		
		int numLevels_init = 10;
		int numLevels_final = 3;
		
		MicroArrayData expression = new MicroArrayData();
		
		try {
			expression.readFile(inputFile);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		MIDiscretize discretizer = new MIDiscretize(expression.values,numLevels_init);
		
		try {
			discretizer.discretizeDownTree(MIFileName);
			discretizer.mergeDownLevels(numLevels_final,levelsFileName);
			discretizer.transposeDExpression();
			expression.setDiscrete();
			expression.dvalues = discretizer.dExpression;
			expression.writeFile(combinedName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void filterFile(String inFName,String outFName,String selectFName,double threshold,double logBase) {
		int minChangeNum = 1;
		
		MicroArrayData data = new MicroArrayData();
		
		try {
			data.readFileBlindContinuous(inFName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		int j = 0;
		
		if (logBase > 0.0) {
			for (i=0;i<data.numRows;i++) {
				for (j=0;j<data.numCols;j++) {
					data.values[i][j] = Math.log(data.values[i][j])/Math.log(logBase);
				}
			}
		}
		
		int numChange = 0;
		HashSet<Integer> drows = new HashSet<Integer>();
		
		for (i=0;i<data.numRows;i++) {
			numChange = 0;
			for (j=0;j<data.numCols;j++) {
				if (Math.abs(data.values[i][j]) >= threshold) {
					numChange++;
				}
			}
			if (numChange < minChangeNum) {
				drows.add(i);
			}
		}
		
		ClusterReader selectReader = new ClusterReader();
		
		if (selectFName != null) {
			try {
				selectReader.readFile(selectFName);
			} catch(IOException e) {
				System.out.println(e);
			}
			
		
			HashSet<String> selectGenes = new HashSet<String>();
			ArrayList<String> genes = selectReader.clusters.get(0);
			for (i=0;i<genes.size();i++) {
				selectGenes.add(genes.get(i));
			}
		
			for (i=0;i<data.numRows;i++) {
				if (!selectGenes.contains(data.geneNames[i])) {
					drows.add(i);
				}
			}
		}
	
		data.deleteRows(drows);
		
		for (i=0;i<data.numRows;i++) {
			for (j=0;j<data.numCols;j++) {
				if (Math.abs(data.values[i][j]) < threshold) {
					data.values[i][j] = 0.0/0.0;
				}
			}
		}
		
		try {
			data.writeFile(outFName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}

}
