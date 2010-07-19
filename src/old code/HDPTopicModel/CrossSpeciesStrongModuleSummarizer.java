package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

import edu.mit.csail.psrg.georg.GO.DAGException;
import edu.mit.csail.psrg.georg.GO.DAGReader;
import edu.mit.csail.psrg.georg.GO.GoOBOReader;
import edu.mit.csail.psrg.georg.GO.MouseHumanGoAssocReader;
import edu.mit.csail.psrg.georg.GO.OntAssociationReader;
import edu.mit.csail.psrg.georg.GO.OntClusters;
import edu.mit.csail.psrg.georg.GO.OntDAG;
import edu.mit.csail.psrg.georg.GO.OntGene;
import edu.mit.csail.psrg.georg.GO.OntGeneAssociations;
import edu.mit.csail.psrg.georg.GO.OntTerm;
import edu.mit.csail.psrg.georg.HierarchicalCluster.HierarchicalCluster;
import edu.mit.csail.psrg.georg.StatUtil.VectorUtil;

public class CrossSpeciesStrongModuleSummarizer {
	int minGenes = 5;
	int maxGenes = 200;
	// minimum count of gene in topic
	int minOccurInTopic = 2;
	double FDR = 0.05;
	// minimum percentage occurance for a gene in a module
	double minGeneOccurPercent = 0.05;
	double minModuleOccurPercent = 0.05;
	double minGeneOverlap = 0.50;
	double finalFilterModuleOccurPercent = 0.50;
	int namespace = OntTerm.BiologicalProcess;
	HierDirichletProcess myHDP = null;
	int numSpecies = 2;
	
	private class geneHolder implements Comparable {
		String[] names = new String[numSpecies];
		double rank = 0.0;
		double[] counts = new double[numSpecies*2];
		
		public geneHolder(String[] nms,double r,double[] c) {
			int i = 0;
			names = nms;
			rank = r;
			counts = c;
		}
		
		public int compareTo(Object o) {
			geneHolder holder = (geneHolder) o;
			if (rank == holder.rank)
				return 0;
			if (rank < holder.rank)
				return 1;
			return -1;
		}
	}
	
	private class geneList {
		ArrayList<geneHolder> genes = new ArrayList<geneHolder>();
		
		void addGene(String[] names,double rank,double[] counts) {
			genes.add(new geneHolder(names,rank,counts));
		}
		
		int size() {
			return genes.size();
		}
		
		double[] getRanks() {
			double[] ranks = new double[genes.size()];
			int i = 0;
			for (i=0;i<genes.size();i++) {
				ranks[i] = genes.get(i).rank;
			}
			return ranks;
		}
		
		void sort() {
			ArrayList<geneHolder> genes2 = new ArrayList<geneHolder>();
			geneHolder[] garray = new geneHolder[genes.size()];
			garray = genes.toArray(garray);
			List<geneHolder> glist = Arrays.asList(garray);
			Collections.sort(glist);
			
			double[] ranks = new double[genes.size()];
			int i = 0;
			Iterator<geneHolder> iter = glist.iterator();
			geneHolder gene = null;
			while(iter.hasNext()) {
				gene = iter.next();
				ranks[i] = gene.rank;
				genes2.add(gene);
				i++;
			}
			
			genes = genes2;
		}
	}
	
	public CrossSpeciesStrongModuleSummarizer(HierDirichletProcess hdp) {
		myHDP = hdp;
	}
	
//	 create strong modules from existing modules
	ArrayList<StrongModule> createStrongModules() {
		int i = 0;
  		int j = 0;
  		int k = 0;
  		double p = 0;
  		DirichletProcess dp = null;
  		Document[] docs = null;
  		Document doc = null;
  		StrongModule[] modules = new StrongModule[myHDP.numClusters];
  		
  		StrongModuleCollector myStrongModuleCollector = (StrongModuleCollector) myHDP.getStatCollector(StrongModuleCollector.class);
  		
  		for (i=0;i<myHDP.numClusters;i++) {
  			modules[i] = new CrossSpeciesStrongModule(myHDP.totalGenes,myHDP.posCounts[i]);
  		}
  		
  		for (j=0;j<myHDP.DP.length;j++) {
  			dp = myHDP.DP[j];
  			if (dp.state != DirichletProcess.HELDOUT) {
				docs = dp.documents;
				if (docs != null) {
					for (k=0;k<docs.length;k++) {
						doc = docs[k];
						doc.allocateGenesToClusterMap(myHDP.tempClusterMap,myHDP.tempClusterArray,myHDP.numClusters);
						for (i=0;i<myHDP.numClusters;i++) {
							p = modules[i].controlPVal(myHDP.tempClusterMap[i],doc,myHDP.totalGenes,minGenes);
	  						if (p <= myStrongModuleCollector.documentControlPvalThresh) {
	  							modules[i].addControllingDocument(doc,myHDP.tempClusterMap[i],myHDP.tempClusterArray[i]);
	  						}
						}
					}
				}
  			}
		}
  		
  		ArrayList<StrongModule> moduleArray = new ArrayList<StrongModule>();
  		for (i=0;i<myHDP.numClusters;i++) {
  			if (modules[i].controllingDocuments.size() >= 1 & modules[i].totalGenesInModule >= myStrongModuleCollector.minGenesInModule) {
  				moduleArray.add(modules[i]);
  			}
  		}
  		return moduleArray;
	}
	
	// merge together modules w/ overlapping genes
	ArrayList<StrongModule> collapseModules() {
		int i = 0;
		int j = 0;
		StrongModuleCollector myStrongModuleCollector = (StrongModuleCollector) myHDP.getStatCollector(StrongModuleCollector.class);
		int nc = myStrongModuleCollector.strongModuleMap.size();
		ArrayList<StrongModule> modules = new ArrayList<StrongModule>();
		ArrayList<HashSet<Integer>> geneSets = new ArrayList<HashSet<Integer>>();
		HashSet<Integer> genes = null;
		Iterator<StrongModule> iter = myStrongModuleCollector.strongModuleMap.values().iterator();
		CrossSpeciesStrongModule module = null;
		double v1 = 0.0;
		double v2 = 0.0;
		double norm = 0.0;
		int geneID = 0;
		Integer[] counts = null;
		Iterator<Integer> gIter = null;
		
		while(iter.hasNext()) {	
			module = (CrossSpeciesStrongModule) iter.next();
			v2 = ((double) module.numOccur)/((double) myStrongModuleCollector.numSamples);
			
			if (v2 >= minModuleOccurPercent) {
				genes = new HashSet<Integer>();
		//		norm = (double) module.numOccur*module.controllingDocuments.size();
				
				gIter = module.geneMap.keySet().iterator();
				while (gIter.hasNext()) {
					geneID = gIter.next();
					counts = module.geneMap.get(geneID);
					v1 = ((double) (counts[0] + counts[2]))/((double) (module.numOccurInSpecies[0] + module.numOccurInSpecies[1]));
				//	v1 = ((double) (counts[1] + counts[3]))/((double) (module.numOccurInSpecies[0] + module.numOccurInSpecies[1]));
					if (v1 >= minGeneOccurPercent) {
				//	if (v1 >= 1.0) {
						genes.add(geneID);
					}
				}
				if (genes.size() >= minGenes) {
					geneSets.add(genes);
					modules.add(module);
				}
			}
		}
		
		nc = modules.size();
		double[][] distances = new double[nc-1][];
		int numUnion = 0;
		int numIntersect = 0;
		genes = new HashSet<Integer>();
		
		for (i=0;i<nc-1;i++) {
			distances[i] = new double[nc-i-1];
			for (j=(i+1);j<nc;j++) {
				genes.clear();
				genes.addAll(geneSets.get(i));
				genes.addAll(geneSets.get(j));
				numUnion = genes.size();
				genes.clear();
				genes.addAll(geneSets.get(i));
				genes.retainAll(geneSets.get(j));
				numIntersect = genes.size();
				distances[i][j-i-1] = 1.0 - ((double) numIntersect)/((double) numUnion);
			}
		}
		
		HierarchicalCluster hierCluster = new HierarchicalCluster();
		hierCluster.cluster(distances,1.0-minGeneOverlap);
		
		ArrayList<StrongModule> modules2 = new ArrayList<StrongModule>();
		CrossSpeciesStrongModule module2 = null;
		
		int moduleNum = 0;
		for (i=0;i<hierCluster.clusters.size();i++) {
			moduleNum = hierCluster.clusters.get(i).get(0);
			module = (CrossSpeciesStrongModule) modules.get(moduleNum);
			modules2.add(module);
			if (hierCluster.clusters.get(i).size() > 1) {
				for (j=1;j<hierCluster.clusters.get(i).size();j++) {
					moduleNum = hierCluster.clusters.get(i).get(j);
					module2 = (CrossSpeciesStrongModule) modules.get(moduleNum);
					module.controllingDocuments.retainAll(module2.controllingDocuments);
				}
				for (j=1;j<hierCluster.clusters.get(i).size();j++) {
					moduleNum = hierCluster.clusters.get(i).get(j);
					module2 = (CrossSpeciesStrongModule) modules.get(moduleNum);
					module.merge(module2);
				}
			}
		}
		
		ArrayList<StrongModule> modules3 = new ArrayList<StrongModule>();
		// perform a final filtering on merged modules
		for (i=0;i<modules2.size();i++) {
			module = (CrossSpeciesStrongModule) modules2.get(i);
			v2 = ((double) module.numOccur)/((double) myStrongModuleCollector.numSamples);
			if (v2 >= finalFilterModuleOccurPercent) {
				modules3.add(module);
			}
		}
		
		return modules3;
	}
	
	public void outputTopics(String fOutNameLong,String fOutNameShort,String fOutTabularName,String fOutNameGeneCross,String fOutNameTopGenes,String tissueMatrixLoadFile,String topicGOCountFile,boolean doCollapseModules) {
		int j = 0;
		int i = 0;
		int k = 0;
		ArrayList<String> genes = new ArrayList<String>();
		double norm = 0.0;
		
		StrongModuleCollector myStrongModuleCollector = (StrongModuleCollector) myHDP.getStatCollector(StrongModuleCollector.class);
		
		ArrayList<StrongModule> modules = null;
		
		if (doCollapseModules) {
			modules = collapseModules();
		} else {
			modules = createStrongModules();
		}
		
		int[][] geneCross = new int[myHDP.geneNames.length][modules.size()];
		
		String[] humanGeneNames = new String[myHDP.geneNames.length];
		String[] mouseGeneNames = new String[myHDP.geneNames.length];
		String s1;
		String s2;
		int f = 0;
		for (i=0;i<myHDP.geneNames.length;i++) {
			f = myHDP.geneNames[i].indexOf("|");
			s1 = myHDP.geneNames[i].substring(0,f);
			s2 = myHDP.geneNames[i].substring(f+1);
			humanGeneNames[i] = s1;
			mouseGeneNames[i] = s2;
		}
		String refDirPath = "C:\\research_data\\mouse_human\\b47mm7\\";
		String humanAssocFile = refDirPath + "Hs.genelist.go.plus";
		String mouseAssocFile = refDirPath + "Mm.genelist.go.plus";
		
		GoOBOReader goReader = new GoOBOReader();
		MouseHumanGoAssocReader humanAssocReader = new MouseHumanGoAssocReader();
		humanAssocReader.setFile(humanAssocFile);
		OntDAG humanDAG = null;
		OntGeneAssociations humanAssoc = null;
		
		MouseHumanGoAssocReader mouseAssocReader = new MouseHumanGoAssocReader();
		mouseAssocReader.setFile(mouseAssocFile);
		OntDAG mouseDAG = null;
		OntGeneAssociations mouseAssoc = null;
		
		HashSet<String> humanUseGenes = new HashSet<String>();
		HashSet<String> mouseUseGenes = new HashSet<String>();
		for (i=0;i<humanGeneNames.length;i++) {
			humanUseGenes.add(humanGeneNames[i]);
			mouseUseGenes.add(mouseGeneNames[i]);
		}
		
		try {
			humanDAG = goReader.readFile();
			System.out.println("Loaded human DAG");
			humanAssoc = humanAssocReader.readFile(humanDAG,humanUseGenes);
			System.out.println("Loaded human associations");
			
			mouseDAG = goReader.readFile();
			System.out.println("Loaded mouse DAG");
			mouseAssoc = mouseAssocReader.readFile(mouseDAG,mouseUseGenes);
			System.out.println("Loaded mouse associations");
		} catch(IOException e) {
			System.out.println(e);
		} catch(DAGException e) {
			System.out.println(e);
		}
		
		OntClusters humanClusters = new OntClusters(humanAssoc,minGenes,maxGenes,FDR,namespace);
		OntClusters mouseClusters = new OntClusters(mouseAssoc,minGenes,maxGenes,FDR,namespace);
		
		try {
			OutputStrongModulesTabular tab = new OutputStrongModulesTabular(myHDP,modules,fOutTabularName,tissueMatrixLoadFile);
		} catch(Exception e) {
			System.out.println(e);
		}
		Iterator<StrongModule> iter = modules.iterator();
		geneList[] geneLists = new geneList[modules.size()];
		
		Iterator<Document> docIter = null; 
		CrossSpeciesStrongModule module = null;
		Iterator<Integer> gIter = null;
		Integer[] counts = null;
		int geneID = 0;
		
		double v1 = 0.0;
		double v2 = 0.0;
		double p_H = 0.0;
		double p_M = 0.0;
		int code = 0;
		String description = "";
		Document doc = null;
		int gene = 0;
		String[] geneNames = null;
		double[] normCounts = null;
		ArrayList<String> humanGeneSetNames = null;
		ArrayList<String> mouseGeneSetNames = null;
		while(iter.hasNext()) {
			genes.clear();
			geneLists[k] = new geneList();
			module = (CrossSpeciesStrongModule) iter.next();
			gIter = module.geneMap.keySet().iterator();
			while(gIter.hasNext()) {
				geneID = gIter.next();
				counts = module.geneMap.get(geneID);
				v1 = ((double) (counts[0] + counts[2]))/((double) (module.numOccurInSpecies[0] + module.numOccurInSpecies[1]));
				v2 = ((double) (counts[1] + counts[3]))/((double) (module.numOccurInSpecies[0] + module.numOccurInSpecies[1]));
				p_H = ((double) counts[0])/((double) module.numOccurInSpecies[0]);
				p_M = ((double) counts[2])/((double) module.numOccurInSpecies[1]);
				
			//	if (p_H >= 0.50 | p_M >= 0.50) {
				if (v1 >= minGeneOccurPercent) {
					geneNames = new String[2];
					normCounts = new double[module.numSpecies*2];
					for(i=0;i<module.numSpecies;i++) {
						norm = (double) module.numOccurInSpecies[i];
						normCounts[i*numSpecies] = ((double) counts[i*numSpecies])/norm;
						normCounts[i*numSpecies+1] = ((double) counts[i*numSpecies+1])/norm;
					}
					gene = myHDP.reverseGeneMap[geneID];
					geneNames[0] = humanGeneNames[gene];
					geneNames[1] = mouseGeneNames[gene];
					geneLists[k].addGene(geneNames,v2,normCounts);
					code = 0;
					if (counts[0] > 0 & counts[2] == 0) {
						code = 1;
					}
					if (counts[0] == 0 & counts[2] > 0) {
						code = 2;
					}
					if (counts[0] > 0 & counts[2] > 0) {
						code = 3;
					}
					geneCross[gene][k] = code;
				}
			}
			description = "";
			docIter = module.controllingDocuments.iterator();
			while(docIter.hasNext()) {
				doc = docIter.next();
			//	v1 = module.documentLoad.get(doc)/((double) module.numOccur);
				v1 = module.documentLoad.get(doc);
				description = description + doc.annotation + " (" + (new Double(v1)).toString() + ") ";
				description = description + " | ";
			}
			geneLists[k].sort();
			v1 = ((double) module.numOccur)/((double) myStrongModuleCollector.numSamples);
			if (v1 >= minModuleOccurPercent & geneLists[k].size() >= minGenes) {
				humanGeneSetNames = new ArrayList<String>();
				mouseGeneSetNames = new ArrayList<String>();
				for (i=0;i<geneLists[k].genes.size();i++) {
					humanGeneSetNames.add(geneLists[k].genes.get(i).names[0]);
					mouseGeneSetNames.add(geneLists[k].genes.get(i).names[1]);
				}
				humanClusters.addCluster(humanGeneSetNames,geneLists[k].getRanks(),description,v1);
				mouseClusters.addCluster(mouseGeneSetNames,geneLists[k].getRanks(),description,v1);
			}
			k++;
		}
		
		humanClusters.clusterSignif();
		mouseClusters.clusterSignif();
		
		int[] humanSignif = humanClusters.numSignificantTerms();
		int[] mouseSignif = mouseClusters.numSignificantTerms();
		try {
			FileWriter gf = new FileWriter(topicGOCountFile);
			for (i=0;i<humanSignif.length;i++) {
				if (humanSignif[i] > 0 | mouseSignif[i] > 0) {
					humanSignif[i] = 1;
				}
				gf.write((new Integer(humanSignif[i])).toString());
				gf.write("\n");
			}
			gf.close();
		} catch(IOException e) {
			System.out.println(e);
		}
		
		try {
			outputClustersLong(fOutNameLong,geneLists,mouseClusters,humanClusters,humanAssoc,mouseAssoc);
			outputGeneCross(fOutNameGeneCross,humanGeneNames,mouseGeneNames,geneCross);
			outputTopGenes(fOutNameTopGenes,geneLists);
		//	clusters.outputClustersShort(fOutNameShort);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public void outputGeneCross(String fOutNameGeneCross,String[] humanGeneNames,String[] mouseGeneNames,int[][] geneCross) throws IOException {
		int i = 0;
		int j = 0;
		FileWriter outFile = new FileWriter(fOutNameGeneCross);
		
		for (i=0;i<humanGeneNames.length;i++) {
			outFile.write(humanGeneNames[i]);
			outFile.write("\t");
			outFile.write(mouseGeneNames[i]);
			for (j=0;j<geneCross[0].length;j++) {
				outFile.write("\t");
				outFile.write((new Integer(geneCross[i][j])).toString());
			}
			outFile.write("\n");
		}
		outFile.close();
	}
	
	public void outputClustersLong(String fName,geneList[] geneLists,OntClusters mouseClusters,OntClusters humanClusters,OntGeneAssociations humanAssoc,OntGeneAssociations mouseAssoc) throws IOException {
		FileWriter outFile = new FileWriter(fName);
		int i = 0;
		for (i=0;i<geneLists.length;i++) {
			outFile.write("Cluster#");
			outFile.write((new Integer(i+1)).toString());
			outFile.write("\t");
			outFile.write((new Integer(mouseClusters.getClusters().get(i).numGenesInSet)).toString());
			outFile.write("\t");
			outFile.write((new Double(mouseClusters.getClusters().get(i).probability)).toString());
			outFile.write("\n");
			if (mouseClusters.getClusters().get(i).description != null) {
				outFile.write(mouseClusters.getClusters().get(i).description);
				outFile.write("\n");
			}
			
			outputGeneInfo(outFile,geneLists[i],humanAssoc,mouseAssoc);
			
			humanClusters.getClusters().get(i).outputScoresLong(outFile,false);
			mouseClusters.getClusters().get(i).outputScoresLong(outFile,false);
		}
		outFile.close();
	}
	
	public void outputGeneInfo(FileWriter outFile,geneList genes,OntGeneAssociations humanAssoc,OntGeneAssociations mouseAssoc) throws IOException {
		int i = 0;
		int j = 0;
		geneHolder gene = null;
		OntGene humanOGene = null;
		OntGene mouseOGene = null;
		
		for (i=0;i<genes.size();i++) {
			gene = genes.genes.get(i);
			outFile.write(gene.names[0]);
			outFile.write("\t");
			outFile.write(gene.names[1]);
			outFile.write("\t");
			outFile.write((new Double(gene.rank)).toString());
			outFile.write("\t");
			for (j=0;j<gene.counts.length;j++) {
				outFile.write((new Double(gene.counts[j])).toString());
				outFile.write("\t");
			}
			humanOGene = humanAssoc.getGene(gene.names[0]);
			mouseOGene = mouseAssoc.getGene(gene.names[1]);
			outputGeneDescription(outFile,humanOGene);
			outFile.write("\t");
			outputGeneDescription(outFile,mouseOGene);
			outFile.write("\n");
		}
	}
	
	public void outputGeneDescription(FileWriter outFile,OntGene gene) throws IOException {
		if (gene == null) {
			outFile.write("none\tnone");
			return;
		}
		if (gene.Name != null) {
			outFile.write(gene.Name);
		} else {
			outFile.write("none");
		}
		outFile.write("\t");
		if (gene.Description != null) {
			outFile.write(gene.Description);
		} else {
			outFile.write("none");
		}
	}
	
	public void outputTopGenes(String fName,geneList[] geneLists) throws IOException {
		int i = 0;
		int j = 0;
		FileWriter outFile = new FileWriter(fName);
		geneHolder gene = null;
		double countThresh = 0.50;
		int maxGenesInCluster = 50;
		Vector<Double> geneRank = new Vector<Double>();
		Vector<Integer> v2 = new Vector<Integer>();
		int[] order = null;
		
		for (i=0;i<geneLists.length;i++) {
			outFile.write((new Integer(i+1)).toString());
			geneRank.clear();
			v2.clear();
		/*	for (j=0;j<geneLists[i].size();j++) {
				gene = geneLists[i].genes.get(j);
				if (gene.counts[0] >= countThresh) {
					geneRank.add(gene.counts[0]);
					v2.add(j);
				}
			}
			if (geneRank.size() >= maxGenesInCluster) {
				order = VectorUtil.sortOrder(geneRank);
				for (j=order.length-1;j>=(order.length-maxGenesInCluster-1);j--) {
					gene = geneLists[i].genes.get(v2.get(order[j]));
					outFile.write("\t");
					outFile.write(gene.names[0]);
				}
			} */
			
			for (j=0;j<geneLists[i].size();j++) {
				gene = geneLists[i].genes.get(j);
				if (gene.counts[0] >= countThresh) {
					outFile.write("\t");
					outFile.write(gene.names[0]);
				}
				if (gene.counts[2] >= countThresh) {
					outFile.write("\t");
					outFile.write(gene.names[1]);
				}
			}
			
			outFile.write("\n");
		}
		
		outFile.close();
	}
}
