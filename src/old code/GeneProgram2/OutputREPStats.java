package edu.mit.csail.psrg.georg.GeneProgram2;

import java.io.IOException;
import java.util.HashMap;

import edu.mit.csail.psrg.georg.DataAccess.CommandLineParser;
import edu.mit.csail.psrg.georg.DataAccess.CommandLineParserException;

public class OutputREPStats {

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
			parser.addSwitch("--baseFileName","String","base file name for files to load and output");
			parser.addSwitch("--minGenesInREP","Integer","minimum number of genes in a REP");
			parser.addSwitch("--nameSpace","String","the type of annotation category to test");
			parser.addSwitch("--tissueUseGenesThreshold","Double","percentage of the REP's genes that must be used by the tissue in a sample to be deemed significant");
			parser.addSwitch("--minGeneOccurPercent","Double","minimum percentage occurance for a gene in a REP");
			parser.addSwitch("--minREPOccurPercent","Double","the minimum percentage of samples the REP must occur in");
			parser.addSwitch("--GeneAssocFileName","String","file name for gene to associations");
			parser.addSwitch("--GOHierarchyFileName","String","file name for the GO DAG (OBO file)");
			parser.addSwitch("--FDR","Double","false discovery rate to use for determining significance of REP GO category enrichment");
			parser.addSwitch("--minGOGenes","Integer","minimum number of genes for a GO category to be tested for enrichment");
			parser.addSwitch("--maxGOGenes","Integer","maximum number of genes for a GO category to be tested for enrichment");
			parser.addSwitch("--outputAnnotations","NoArgument","output annotations for each gene");
		} catch(CommandLineParserException e) {
			System.out.println(e);
		}
		return parser;
	}
	
	public static void doOutput(String[] args) {
		REPStatBuilder builder = new REPStatBuilder();
		
		String baseFileName = null;
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
			
			if (parseMap.containsKey("--baseFileName")) {
				baseFileName = (String) parseMap.get("--baseFileName");
			}
			
			if (parseMap.containsKey("--nameSpace")) {
				builder.namespace = (String) parseMap.get("--nameSpace");
			}
			
			if (parseMap.containsKey("--minGenesInREP")) {
				builder.minGenesInREP = (Integer) parseMap.get("--minGenesInREP");
			}
			
			if (parseMap.containsKey("--tissueUseGenesThreshold")) {
				builder.tissueUseGenesThreshold = (Double) parseMap.get("--tissueUseGenesThreshold");
			}
			
			if (parseMap.containsKey("--minGeneOccurPercent")) {
				builder.minGeneOccurPercent = (Double) parseMap.get("--minGeneOccurPercent");
			}
			
			if (parseMap.containsKey("--minREPOccurPercent")) {
				builder.minREPOccurPercent = (Double) parseMap.get("--minREPOccurPercent");
			}
			
			if (parseMap.containsKey("--GeneAssocFileName")) {
				builder.GeneAssocFileName = (String) parseMap.get("--GeneAssocFileName");
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
			
			if (parseMap.containsKey("--outputAnnotations")) {
				builder.outputAnnotations = true;
			}
			
			parser.outputMap(parseMap);
			
		} catch(CommandLineParserException e) {
			System.out.println(e);
			return;
		} catch(IOException e) {
			System.out.println(e);
			return;
		}
		
		try {
			builder.outputREPStats(baseFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}

}
