package edu.mit.csail.psrg.georg.Annotation;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;
import edu.mit.csail.psrg.georg.DataAccess.CommandLineParser;
import edu.mit.csail.psrg.georg.DataAccess.CommandLineParserException;

public class ScoreGeneSets {
	//	 parameters for GO categories
	public int minGOGenes = 5;
	public int maxGOGenes = 200;
	public double FDR = 0.05;
	public String namespace = "biological_process";
	public String GeneAssocFileName = "C:\\research_data\\mouse_human\\b47mm7\\Hs.genelist.go.plus";
	public String GOHierarchyFileName = "C:\\research_data\\go_data\\gene_ontology_plus.obo";
	public String inputFile;
	public String outputFile;
	
	public static CommandLineParser constructParser() {
		CommandLineParser parser = new CommandLineParser();
		try {
			parser.addSwitch("--controlFile","String","file with parameter settings");
			parser.addSwitch("--inputFile","String","file with gene sets");
			parser.addSwitch("--outputFile","String","file to output significant gdfmx");
			parser.addSwitch("--nameSpace","String","the type of annotation category to test");
			parser.addSwitch("--GeneAssocFileName","String","file name for gene to associations");
			parser.addSwitch("--GOHierarchyFileName","String","file name for the GO DAG (OBO file)");
			parser.addSwitch("--FDR","Double","false discovery rate to use for determining significance of GO category enrichment");
			parser.addSwitch("--minGOGenes","Integer","minimum number of genes for a GO category to be tested for enrichment");
			parser.addSwitch("--maxGOGenes","Integer","maximum number of genes for a GO category to be tested for enrichment");
		} catch(CommandLineParserException e) {
			System.out.println(e);
		}
		return parser;
	}
	
	public static void main(String[] args) {
		
		ScoreGeneSets scorer = new ScoreGeneSets();
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
			
			if (parseMap.containsKey("--inputFile")) {
				scorer.inputFile = (String) parseMap.get("--inputFile");
			}
			
			if (parseMap.containsKey("--outputFile")) {
				scorer.outputFile = (String) parseMap.get("--outputFile");
			}
			
			if (parseMap.containsKey("--nameSpace")) {
				scorer.namespace = (String) parseMap.get("--nameSpace");
			}
			
			
			if (parseMap.containsKey("--GeneAssocFileName")) {
				scorer.GeneAssocFileName = (String) parseMap.get("--GeneAssocFileName");
			}
			
			if (parseMap.containsKey("--GOHierarchyFileName")) {
				scorer.GOHierarchyFileName = (String) parseMap.get("--GOHierarchyFileName");
			}
	
			if (parseMap.containsKey("--FDR")) {
				scorer.FDR = (Double) parseMap.get("--FDR");
			}
			
			if (parseMap.containsKey("--minGOGenes")) {
				scorer.minGOGenes = (Integer) parseMap.get("--minGOGenes");
			}
			
			if (parseMap.containsKey("--maxGOGenes")) {
				scorer.maxGOGenes = (Integer) parseMap.get("--maxGOGenes");
			}
			
			parser.outputMap(parseMap);
			
		} catch(CommandLineParserException e) {
			System.out.println(e);
			return;
		} catch(IOException e) {
			System.out.println(e);
			return;
		}
		
		scorer.scoreCategories();
	}
	
	public void scoreCategories() {
		Annotations myAnnotations = new Annotations();
		
		ClusterReader genesFile = new ClusterReader();
		try {
			genesFile.readFile(inputFile);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		HashSet<String> useGenes = new HashSet<String>();
		int i = 0;
		int j = 0;
		ArrayList<String> genes = null;
		for (i=0;i<genesFile.clusters.size();i++) {
			genes = genesFile.clusters.get(i);
			for (j=0;j<genes.size();j++) {
				useGenes.add(genes.get(j));
			}
		}
		
		try {
			myAnnotations.readFiles(GeneAssocFileName,GOHierarchyFileName,useGenes);
			System.out.println("Loaded annotations");
		} catch(IOException e) {
			System.out.println(e);
		} catch(DAGException e) {
			System.out.println(e);
		}
		
		GeneSetCollection geneSets = new GeneSetCollection(myAnnotations,minGOGenes,maxGOGenes,FDR,namespace);
		OntGene oGene = null;
		
		for (i=0;i<genesFile.clusters.size();i++) {
			genes = genesFile.clusters.get(i);
			geneSets.addGeneSet(genes);
		}
		
		geneSets.geneSetSignif();
		
		try {
			FileWriter file = new FileWriter(outputFile);
		
			for (i=0;i<geneSets.geneSets.size();i++) {
				geneSets.geneSets.get(i).outputScoresOneLine(file);
			}
			file.close();
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
}
