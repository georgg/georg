package edu.mit.csail.psrg.georg.GeneProgram2;

import java.io.IOException;
import java.util.HashMap;

import edu.mit.csail.psrg.georg.DataAccess.CommandLineParser;
import edu.mit.csail.psrg.georg.DataAccess.CommandLineParserException;

public class OutputREPs {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		doOutput(args);
	}
	
	public static CommandLineParser constructParser() {
		CommandLineParser parser = new CommandLineParser();
		try {
			parser.addSwitch("--controlFile","String","file with parameter settings");
			parser.addSwitch("--basePersistFileName","String","base file name for persist files to load");
			parser.addSwitch("--baseOutputFileName","String","base file name for output files");
			parser.addSwitch("--sampleStart","Integer","starting sample number");
			parser.addSwitch("--sampleIncrement","Integer","sample increment");
			parser.addSwitch("--numSamples","Integer","total number of samples to load");
			parser.addSwitch("--modifierFileName","String","file name for modifier names");
			parser.addSwitch("--minGenesInREP","Integer","minimum number of genes in a REP");
			parser.addSwitch("--tissueUseGenesThreshold","Double","percentage of the REP's genes that must be used by the tissue in a sample to be deemed significant");
			parser.addSwitch("--minSimilarity","Double","minimum similarity to merge REPs");
			parser.addSwitch("--minGeneOccurPercent","Double","minimum percentage occurance for a gene in a REP");
			parser.addSwitch("--minREPOccurPercent","Double","the minimum percentage of samples the REP must occur in");
			parser.addSwitch("--minREPPercentSignificance","Double","the minimum percentage of samples the REP must be significant in");
			parser.addSwitch("--minGroupMerge","Double","the minimum similarity to merge tissues into the same consensus group");
			parser.addSwitch("--GOAssocFileName","String","file name for gene to GO categories associations");
			parser.addSwitch("--GOHierarchyFileName","String","file name for the GO DAG (OBO file)");
			parser.addSwitch("--FDR","Double","false discovery rate to use for determining significance of REP GO category enrichment");
			parser.addSwitch("--minGOGenes","Integer","minimum number of genes for a GO category to be tested for enrichment");
			parser.addSwitch("--maxGOGenes","Integer","maximum number of genes for a GO category to be tested for enrichment");
		} catch(CommandLineParserException e) {
			System.out.println(e);
		}
		return parser;
	}
	
	public static void doOutput(String[] args) {
		REPBuilder builder = new REPBuilder();
		
		String basePersistFileName = null;
		String baseOutputFileName = null;
		int sampleStart = 0;
		int sampleIncrement = 0;
		int numSamples = 0;
		String modifierFileName = null;
		String controlFile = null;
		
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
			
			if (parseMap.containsKey("--basePersistFileName")) {
				basePersistFileName = (String) parseMap.get("--basePersistFileName");
			}
			
			if (parseMap.containsKey("--baseOutputFileName")) {
				baseOutputFileName = (String) parseMap.get("--baseOutputFileName");
			}
			
			if (parseMap.containsKey("--sampleStart")) {
				sampleStart = (Integer) parseMap.get("--sampleStart");
			}
			
			if (parseMap.containsKey("--sampleIncrement")) {
				sampleIncrement = (Integer) parseMap.get("--sampleIncrement");
			}
			
			if (parseMap.containsKey("--numSamples")) {
				numSamples = (Integer) parseMap.get("--numSamples");
			}
			
			if (parseMap.containsKey("--modifierFileName")) {
				modifierFileName = (String) parseMap.get("--modifierFileName");
			}
			
			if (parseMap.containsKey("--minGenesInREP")) {
				builder.minGenesInREP = (Integer) parseMap.get("--minGenesInREP");
			}
			
			if (parseMap.containsKey("--tissueUseGenesThreshold")) {
				builder.tissueUseGenesThreshold = (Double) parseMap.get("--tissueUseGenesThreshold");
			}
			
			if (parseMap.containsKey("--minSimilarity")) {
				builder.minSimilarity = (Double) parseMap.get("--minSimilarity");
			}
			
			if (parseMap.containsKey("--minGeneOccurPercent")) {
				builder.minGeneOccurPercent = (Double) parseMap.get("--minGeneOccurPercent");
			}
			
			if (parseMap.containsKey("--minREPOccurPercent")) {
				builder.minREPOccurPercent = (Double) parseMap.get("--minREPOccurPercent");
			}
			
			if (parseMap.containsKey("--minREPPercentSignificance")) {
				builder.minREPPercentSignificance = (Double) parseMap.get("--minREPPercentSignificance");
			}
			
			if (parseMap.containsKey("--minGroupMerge")) {
				builder.minGroupMerge = (Double) parseMap.get("--minGroupMerge");
			}
			
			if (parseMap.containsKey("--GOAssocFileName")) {
				builder.GOAssocFileName = (String) parseMap.get("--GOAssocFileName");
			}
			
			if (parseMap.containsKey("--GOHierarchyFileName")) {
				builder.GOHierarchyFileName = (String) parseMap.get("--GOHierarchyFileName");
			}
	
			if (parseMap.containsKey("--FDR")) {
				builder.FDR = (Double) parseMap.get("--FDR");
			}
			
			if (parseMap.containsKey("--minGOGenes")) {
				builder.minGOGenes = (Integer) parseMap.get("--minGOGenes");
			}
			
			if (parseMap.containsKey("--maxGOGenes")) {
				builder.maxGOGenes = (Integer) parseMap.get("--maxGOGenes");
			}
			
			parser.outputMap(parseMap);
			
		} catch(CommandLineParserException e) {
			System.out.println(e);
			return;
		} catch(IOException e) {
			System.out.println(e);
			return;
		}
		
		builder.processSamples(basePersistFileName,baseOutputFileName,modifierFileName,sampleStart,sampleIncrement,numSamples);
	}

}
