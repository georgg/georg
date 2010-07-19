package edu.mit.csail.psrg.georg.GeneProgram2;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;
import edu.mit.csail.psrg.georg.DataAccess.CommandLineParser;
import edu.mit.csail.psrg.georg.DataAccess.CommandLineParserException;
import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

// Do simple filtering on data: remove genes w/ insufficient
// fold-change, etc.
// also, optionally log transform the data and
// thresh-hold values
public class PreprocessData {
	
	/**
	 * @param args
	 */
	public static CommandLineParser constructParser() {
		CommandLineParser parser = new CommandLineParser();
		try {
			parser.addSwitch("--controlFile","String","file with parameter settings");
			parser.addSwitch("--inputFileName","String","file name for input data file");
			parser.addSwitch("--baseOutputFileName","String","base name for output files");
			parser.addSwitch("--selectGenesFileName","String","file specifying genes to use");
			parser.addSwitch("--threshold","Double","threshold (absolute value) to consider gene differentially expressed");
			parser.addSwitch("--logBase","Double","take the log of the data to using this base (0 = don't take the log)");
			parser.addSwitch("--noUpDown","NoArgument","discretize data using induction only");
			parser.addSwitch("--simpleUpDown","NoArgument","discretize data using symmetric up/down discretization");
			parser.addSwitch("--numLevels_start","Integer","number of discretization levels to start with");
			parser.addSwitch("--numLevels_end","Integer","number of discretization levels to end with");
			parser.addSwitch("--minChange","Integer","minimum number of tissues a gene must show change in to be included");
			parser.addSwitch("--baseModifierFileName","String","file name for input modifier files");
			parser.addSwitch("--numModifiers","Integer","number of modifiers");
		} catch(CommandLineParserException e) {
			System.out.println(e);
		}
		return parser;
	}
	
	public static void main(String[] args) {
		doPreprocess(args);
	}
	
	public static void doPreprocess(String[] args) {
		String controlFile = null;
		String inputFileName = null;
		String baseOutputFileName = null;
		String selectGenesFileName = null;
		// fold change of 1.5 in log2 space
		double threshold = 0.5850;
		double logBase = 0.0;
		int numLevels_start = 10;
		int numLevels_end = 3;
		int minChange = 1;
		String baseModifierFileName = null;
		int numModifiers = 0;
		boolean noUpDown = false;
		boolean simpleUpDown = false;
		
		try {
			CommandLineParser parser = constructParser();
			if (args.length > 0) {
				if (args[0].equals("--help")) {
					parser.displayDescriptions();
					return;
				}
			}
			HashMap<String,Object> parseMap = parser.parseCommandLine(args);
			
			if (parseMap.containsKey("--controlFile")) {
				controlFile = (String) parseMap.get("--controlFile");
			}
			
			if (controlFile != null) {
				parseMap = parser.parseCommandFile(controlFile);
				HashMap<String,Object> parseMap2 = parser.parseCommandLine(args);
				parseMap.putAll(parseMap2);
			}
			
			if (parseMap.containsKey("--inputFileName")) {
				inputFileName = (String) parseMap.get("--inputFileName");
			}
			
			if (parseMap.containsKey("--baseOutputFileName")) {
				baseOutputFileName = (String) parseMap.get("--baseOutputFileName");
			}
			
			if (parseMap.containsKey("--selectGenesFileName")) {
				selectGenesFileName = (String) parseMap.get("--selectGenesFileName");
			}
			
			if (parseMap.containsKey("--threshold")) {
				threshold = (Double) parseMap.get("--threshold");
			}
			
			if (parseMap.containsKey("--noUpDown")) {
				noUpDown = true;
			}
			
			if (parseMap.containsKey("--simpleUpDown")) {
				simpleUpDown = true;
			}
			
			if (parseMap.containsKey("--logBase")) {
				logBase = (Double) parseMap.get("--logBase");
			}
			
			if (parseMap.containsKey("--numLevels_start")) {
				numLevels_start = (Integer) parseMap.get("--numLevels_start");
			}
			
			if (parseMap.containsKey("--numLevels_end")) {
				numLevels_end = (Integer) parseMap.get("--numLevels_end");
			}
			
			if (parseMap.containsKey("--minChange")) {
				minChange = (Integer) parseMap.get("--minChange");
			}
			
			if (parseMap.containsKey("--baseModifierFileName")) {
				baseModifierFileName = (String) parseMap.get("--baseModifierFileName");
			}
			
			if (parseMap.containsKey("--numModifiers")) {
				numModifiers = (Integer) parseMap.get("--numModifiers");
			}
			
			parser.outputMap(parseMap);
			
		} catch(CommandLineParserException e) {
			System.out.println(e);
			return;
		} catch(IOException e) {
			System.out.println(e);
			return;
		}
		
		String outFilteredName = baseOutputFileName + "_filtered.txt";
		
		filterFile(inputFileName,outFilteredName,selectGenesFileName,threshold,logBase,minChange,baseModifierFileName,numModifiers,noUpDown,simpleUpDown);
		
		doDiscretization(baseOutputFileName,outFilteredName,numLevels_start,numLevels_end,noUpDown,simpleUpDown);
	}
	
	public static void doDiscretization(String baseOutputFileName,String inputFile,int numLevels_init,int numLevels_final,boolean noUpDown,boolean simpleUpDown) {
		String combinedName = baseOutputFileName + "_discretized.txt";
		String levelsFileName = baseOutputFileName + "_discretized_levels";
		String MIFileName = baseOutputFileName + "_discretized_MI.txt";
		
		MicroArrayData expression = new MicroArrayData();
		
		try {
			expression.readFile(inputFile);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		MIDiscretizeUpDown discretizerUD = null;
		MIDiscretize discretizer = null;
		
		if (noUpDown | simpleUpDown) {
			discretizer = new MIDiscretize(expression.values,numLevels_init);
		} else {
			discretizerUD = new MIDiscretizeUpDown(expression.values,numLevels_init);
		}
		
		try {
			if (noUpDown | simpleUpDown) {
				discretizer.discretizeDownTree(MIFileName);
				discretizer.mergeDownLevels(numLevels_final,levelsFileName);
				discretizer.transposeDExpression();
				expression.setDiscrete();
				expression.dvalues = discretizer.dExpression;
			} else {
				discretizerUD.discretizeDownTree(MIFileName);
				discretizerUD.mergeDownLevels(numLevels_final,levelsFileName);
				discretizerUD.transposeDExpression();
				expression.setDiscrete();
				expression.dvalues = discretizerUD.dExpression;
			}
			expression.writeFile(combinedName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void filterFile(String inFName,String outFName,String selectFName,double threshold,double logBase,int minChangeNum,String baseModifierFileName,int numModifiers,boolean noUpDown,boolean simpleUpDown) {
		
		MicroArrayData data = new MicroArrayData();
		
		try {
			data.readFileBlindContinuous(inFName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		int j = 0;
		
		MicroArrayData simpleMod = null;
		if (simpleUpDown) {
			simpleMod = new MicroArrayData();
			simpleMod.numCols = data.numCols;
			simpleMod.numRows = data.numRows;
			simpleMod.discreteData = true;
			simpleMod.dvalues = new int[simpleMod.numRows][simpleMod.numCols];
			simpleMod.geneNames = new String[data.geneNames.length];
			for (i=0;i<data.geneNames.length;i++) {
				simpleMod.geneNames[i] = data.geneNames[i];
			}
			
			simpleMod.experimentNames = data.experimentNames;
			
			if (logBase > 0.0) {
				for (i=0;i<data.numRows;i++) {
					for (j=0;j<data.numCols;j++) {
						data.values[i][j] = Math.log(data.values[i][j])/Math.log(logBase);
					}
				}
			}
			
			for (i=0;i<data.numRows;i++) {
				for (j=0;j<data.numCols;j++) {
					if (data.values[i][j] >= 0.0) {
						simpleMod.dvalues[i][j] = 1;
					}
					data.values[i][j] = Math.abs(data.values[i][j]);
				}
			}
		} else {
			if (noUpDown) {
				for (i=0;i<data.numRows;i++) {
					for (j=0;j<data.numCols;j++) {
						data.values[i][j] = Math.abs(data.values[i][j]);
					}
				}
			}
		
			if (logBase > 0.0) {
				for (i=0;i<data.numRows;i++) {
					for (j=0;j<data.numCols;j++) {
						data.values[i][j] = Math.log(data.values[i][j])/Math.log(logBase);
					}
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
		
		if (simpleUpDown) {
			String mfn = baseModifierFileName + "1_filtered.txt";
			simpleMod.deleteRows(drows);
			try {
				simpleMod.writeFile(mfn);
			} catch(IOException e) {
				System.out.println(e);
			}
		}
		
		if (numModifiers > 0) {
			MicroArrayData modifyFile = null;
			String mfn = "";
			
			for (i=0;i<numModifiers;i++) {
				modifyFile = new MicroArrayData();
				mfn = baseModifierFileName + (new Integer(i+1)).toString()+".txt";
				try {
					modifyFile.readFileBlindDiscrete(mfn);
					modifyFile.deleteRows(drows);
					mfn = baseModifierFileName + (new Integer(i+1)).toString()+"_filtered.txt";
					modifyFile.writeFile(mfn);
				} catch(IOException e) {
					System.out.println(e);
				}
			}
		}
	}

}
