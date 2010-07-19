package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.DiscretizeData.DiscretizeMouseHuman;
import edu.mit.csail.psrg.georg.GO.DAGException;
import edu.mit.csail.psrg.georg.GO.DAGReader;
import edu.mit.csail.psrg.georg.GO.OntAssociationReader;
import edu.mit.csail.psrg.georg.GO.OntClusters;
import edu.mit.csail.psrg.georg.GO.OntDAG;
import edu.mit.csail.psrg.georg.GO.OntGeneAssociations;
import edu.mit.csail.psrg.georg.GO.OntTerm;
import edu.mit.csail.psrg.georg.GO.OntTermScores;
import edu.mit.csail.psrg.georg.HierarchicalCluster.HierarchicalCluster;

public class StrongModuleSummarizer {
	OntGeneAssociations assoc = null;
	OntDAG DAG = null;
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
//	int namespace = OntTerm.AnyCategory;
	int namespace = OntTerm.BiologicalProcess;
	HashSet<String> useGenes = new HashSet<String>();
	HierDirichletProcess myHDP = null;
	
	private class Edge {
		int source = 0;
		int destination = 0;
		double strength = 0.0;
		Edge(int s,int d) {
			source = s;
			destination = d;
		}
	}
	
	private class UndirectedEdge extends Edge {
		UndirectedEdge(int s,int d) {
			super(s,d);
		}
		
		public String toString() {
			String s = "";
			if (source > destination) {
				s = (new Integer(destination)).toString() + "," + (new Integer(source)).toString();
			} else {
				s = (new Integer(source)).toString() + "," + (new Integer(destination)).toString();
			}
			return s;
		}
	}
	
	private class geneHolder implements Comparable {
		String name = "";
		double rank = 0.0;
		
		public geneHolder(String n,double r) {
			name = n;
			rank = r;
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
		
		void addGene(String name,double rank) {
			genes.add(new geneHolder(name,rank));
		}
		
		int size() {
			return genes.size();
		}
		
		double[] sort(ArrayList<String> geneNames) {
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
				geneNames.add(gene.name);
				i++;
			}
			
			return ranks;
		}
	}
	
	public StrongModuleSummarizer(HierDirichletProcess init_HDP,DAGReader init_DAGReader, OntAssociationReader init_ontAssociationReader) {
		processGenes(init_HDP);
		readDAG(init_DAGReader);
		readGeneAssocations(init_ontAssociationReader);
	}
	
	public void setMinOccurInTopic(int m) {
		minOccurInTopic = m;
	}
	
	void processGenes(HierDirichletProcess init_HDP) {
		myHDP = init_HDP;
		int i = 0;
		for (i=0;i<myHDP.geneNames.length;i++) {
			useGenes.add(myHDP.geneNames[i]);
		}
	}
	
	void readDAG(DAGReader reader) {
		try {
				DAG = reader.readFile();
				System.out.println("Loaded DAG");
			} catch(IOException e) {
				System.out.println(e);
			} catch(DAGException e) {
				System.out.println(e);
			}
	}
	
	void readGeneAssocations(OntAssociationReader reader) {
		try {
			assoc = reader.readFile(DAG,useGenes);
			System.out.println("Loaded associations");
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public void outputTopics(String fOutNameLong,String fOutNameShort,String fOutGMLName, String fOutMotifsName, boolean doCollapseModules, boolean human) {
		int j = 0;
		int i = 0;
		ArrayList<String> genes = new ArrayList<String>();
		OntClusters clusters = new OntClusters(assoc,minGenes,maxGenes,FDR,namespace);
		geneList glist = null;
		double[] ranks = null;
		double norm = 0.0;
		
		ArrayList<StrongModule> modules = null;
		
		if (doCollapseModules) {
			modules = collapseModules();
		} else {
			modules = createStrongModules();
		}
		
		//	outputYFilesGML(fOutGMLName,modules);
		OutputStrongModulesTabular tab = new OutputStrongModulesTabular(myHDP,modules,fOutGMLName);
		
	//	Iterator<StrongModule> iter = myHDP.strongModuleMap.values().iterator();
		Iterator<StrongModule> iter = modules.iterator();
		
		Iterator<Document> docIter = null; 
		SingleStrongModule module = null;
		double v1 = 0.0;
		double v2 = 0.0;
		String description = "";
		Document doc = null;
		int gene = 0;
		while(iter.hasNext()) {
			genes.clear();
			glist = new geneList();
			module = (SingleStrongModule) iter.next();
			norm = (double) module.numOccur*module.controllingDocuments.size();
			for (i=0;i<module.geneCounts.length;i++) {
				v1 = ((double) module.geneCounts[i])/norm;
				v2 = ((double) module.geneIntensity[i])/norm;
				if (v1 >= minGeneOccurPercent) {
					gene = myHDP.reverseGeneMap[i];
					glist.addGene(myHDP.geneNames[gene],v2);
				}
			}
			description = "";
			docIter = module.controllingDocuments.iterator();
			while(docIter.hasNext()) {
				doc = docIter.next();
				v1 = module.documentLoad.get(doc)/((double) module.numOccur);
				description = description + doc.annotation + " (" + (new Double(v1)).toString() + ") ";
				description = description + " | ";
			}
			ranks = glist.sort(genes);
			v1 = ((double) module.numOccur)/((double) myHDP.numClustersSamples.size());
			if (v1 >= minModuleOccurPercent & glist.size() >= minGenes) {
				clusters.addCluster(genes,ranks,description,v1);
			}
		}
		
		clusters.clusterSignif();
		
		try {
			clusters.outputClustersLong(fOutNameLong);
			clusters.outputClustersShort(fOutNameShort);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		outputMotifs(human,fOutMotifsName,clusters);
	}
	
	public void outputMotifs(boolean human,String foutName,OntClusters clusters) {
		String dirName = "C:\\research_data\\mouse_human\\homo_b47_data\\";
		String humanMotifDName = "human_5kb_conformed.out";
		String mouseMotifDName = "mouse_5kb_conformed.out";
		
		MicroArrayData motifsD = new MicroArrayData();
		
		try {
			if (human) {
				motifsD.readFile(dirName + humanMotifDName);
			} else {
				motifsD.readFile(dirName + mouseMotifDName);
			}
		} catch (IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		int j = 0;
		int k = 0;
		int gg = 0;
		
		HashMap<String,Integer> geneMap = new HashMap<String,Integer>();
		for (i=0;i<motifsD.numRows;i++) {
			geneMap.put(motifsD.geneNames[i],i);
		}
		
		double[][] values = new double[clusters.getClusters().size()][motifsD.numCols];
		double[] mu = new double[motifsD.numCols];
		double[] geneRanks = null;
		OntTermScores score = null;
		double totalV = 0.0;
		
		for (i=0;i<motifsD.numRows;i++) {
			for (j=0;j<motifsD.numCols;j++) {
				mu[j] += motifsD.values[i][j];
			}
		}
		for (j=0;j<motifsD.numCols;j++) {
			mu[j] = mu[j]/((double) motifsD.numRows);
		}
		
		for (i=0;i<values.length;i++) {
			score = clusters.getClusters().get(i);
			geneRanks = score.getGeneRanks();
			totalV = 0.0;
			for (k=0;k<score.getGenes().size();k++) {
				gg = geneMap.get(score.getGenes().get(k).ID);
				if (geneRanks[k] >= 1.5) {
					for (j=0;j<values[0].length;j++) {
					//	values[i][j] += motifsD.values[gg][j]*geneRanks[k];
						values[i][j] += motifsD.values[gg][j];
					}
				//	totalV = totalV + geneRanks[k];
					totalV++;
				}
			}
			for (j=0;j<values[0].length;j++) {
				values[i][j] = values[i][j]/totalV;
				values[i][j] = values[i][j]/mu[j];
			}
		}
		
		HashMap<String,String> motifNameMap = DiscretizeMouseHuman.readMotifNames();
		
		motifsD.values = values;
		motifsD.numCols = values[0].length;
		motifsD.numRows = values.length;
		String[] clusterNames = new String[motifsD.numRows];
		for (i=0;i<clusterNames.length;i++) {
			clusterNames[i] = clusters.getClusters().get(i).getDescription();
		}
		motifsD.geneNames = clusterNames;
		String[] experimentNames = new String[motifsD.experimentNames.length];
		for (j=0;j<experimentNames.length;j++) {
			if (motifNameMap.containsKey(motifsD.experimentNames[j])) {
				motifsD.experimentNames[j] = motifNameMap.get(motifsD.experimentNames[j]);
			}
		}
		
		try {
			motifsD.writeFile(foutName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public void outputYFilesGML(String foutName,ArrayList<SingleStrongModule> modules) {
		FileWriter file = null;
		HashMap<DirichletProcess,Integer> DPMap = null;
		
		try {
			file = new FileWriter(foutName);
			outputYFilesHeader(file);
			DPMap = outputHDPNodes(file);
			outputTopicNodes(file,modules,DPMap);
			outputEdges(file,modules,DPMap);
			
			file.write("]\n");
			file.close();
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public void outputYFilesHeader(FileWriter file) throws IOException {
		file.write("Creator\t\"yFiles\"\nVersion	2.2\n");
		file.write("graph\n[\n");
		file.write("directed 1\n");
		file.write("hierarchic 1\n");
	}
	
	HashMap<DirichletProcess,Integer> outputHDPNodes(FileWriter file) throws IOException {
		int i = 0;
		DirichletProcess DP = null;
		HashMap<DirichletProcess,Integer> DPMap = new HashMap<DirichletProcess,Integer>();
		
		for (i=0;i<myHDP.DP.length;i++) {
			DP = myHDP.DP[i];
			DPMap.put(DP,i+1);
			// DPs w/ no documents serve to group others in the hierarchy
			if (DP.documents == null) {
				file.write("node\n[\n");
				file.write("id " + (new Integer(i+1)).toString() + "\n");
				file.write("label " + "\"" + processHDPLabel(DP.getLabel()) + "\"\n");
				file.write("isGroup 1\n");
				
				if (DP.parent != null) {
					file.write("gid " + (DPMap.get(DP.parent)).toString() + "\n");
				}
				file.write("]\n");
			} else {
				file.write("node\n[\n");
				file.write("id " + (new Integer(i+1)).toString() + "\n");
				file.write("label " + "\"" + processHDPLabel(DP.getLabel()) + "\"\n");
				if (DP.parent != null) {
					file.write("gid " + (DPMap.get(DP.parent)).toString() + "\n");
				}
				file.write("]\n");
			}
		}
		
		return DPMap;
	}
	
	String processHDPLabel(String label1) {
		String label2 = label1;
		int f = label1.indexOf("human_");
		if (f == 0) {
			label2 = label1.substring(f+6);
		} else {
			f = label1.indexOf("mouse_");
			if (f == 0) {
				label2 = label1.substring(f+6);
			} else {
				label2 = label1;
			}
		}
		return label2;
	}
	
	void outputTopicNodes(FileWriter file,ArrayList<SingleStrongModule> modules,HashMap<DirichletProcess,Integer> DPMap) throws IOException {
		int i = 0;
		int gid = 0;
		String nodeLabel = "";
		
		for (i=0;i<modules.size();i++) {
			gid = closestRelative(modules.get(i),DPMap);
			file.write("node\n[\n");
			file.write("id " + (new Integer(i+1+DPMap.size())).toString() + "\n");
			nodeLabel = (new Integer(i+1)).toString();
			file.write("label " + "\"" + nodeLabel + "\"\n");
			file.write("graphics\n[\n");
			file.write("type \"circle\"\n");
			file.write("]\n");
			if (gid > -1) {
				file.write("gid " + (new Integer(gid)).toString() + "\n");
			}
			file.write("]\n");
		}
	}
	
	void outputEdges(FileWriter file,ArrayList<SingleStrongModule> modules,HashMap<DirichletProcess,Integer> DPMap) throws IOException {
		int i = 0;
		int j = 0;
		int k = 0;
		int source = 0;
		int destination = 0;
		int relative = 0;
		
		Document doc = null;
		Iterator<Document> docIter = null;
		SingleStrongModule module = null;
		HashSet<Integer> redundantSource = new HashSet<Integer>();
		Integer[] sources = null;
		HashMap<String,UndirectedEdge> uniqueEdges = new HashMap<String,UndirectedEdge>();
		UndirectedEdge myEdge = null;
		UndirectedEdge myEdge2 = null;
		boolean withinGroup = false;
		String ec = "";
		double maxLoad = 0.0;
		double v = 0.0;
		
		for (i=0;i<modules.size();i++) {
			module = modules.get(i);
			destination = i + 1 + DPMap.size();
			relative = closestRelative(module,DPMap);
			docIter = module.controllingDocuments.iterator();
			maxLoad = 0.0;
			while(docIter.hasNext()) {
				doc = docIter.next();
				v = module.documentLoad.get(doc)/module.numOccur;
				if (v > maxLoad) {
					maxLoad = v;
				}
			}
			docIter = module.controllingDocuments.iterator();
			redundantSource.clear();
			withinGroup = false;
			while(docIter.hasNext()) {
				doc = docIter.next();
			//	source = DPMap.get(doc.parent);
				source = DPMap.get(doc.parent.parent);
				if (relative == source) {
					withinGroup = true;
					source = DPMap.get(doc.parent);
				}
				
			//	if (!redundantSource.contains(source) & withinGroup) {
				if (!redundantSource.contains(source)) {
					if (withinGroup) {
						ec = edgeColor(module.documentLoad.get(doc)/module.numOccur);
					} else {
						ec = edgeColor(maxLoad);
					}
					file.write("edge\n[\n");
					file.write("source " + source + "\n");
					file.write("target " + destination + "\n");
					file.write("graphics [\n");
					file.write("fill \"" + ec + "\"\n");
					file.write("]\n");
					file.write("]\n");
				}
				redundantSource.add(source);
			}
			if (!withinGroup) {
				sources = new Integer[redundantSource.size()];
				sources = redundantSource.toArray(sources);
				if (sources.length > 1) {
					for (j=0;j<sources.length-1;j++) {
						source = sources[j];
						for (k=j+1;k<sources.length;k++) {
							destination = sources[k];
							myEdge = new UndirectedEdge(source,destination);
							myEdge.strength = maxLoad;
							if (uniqueEdges.containsKey(myEdge.toString())) {
								myEdge2 = uniqueEdges.get(myEdge.toString());
								if (myEdge.strength > myEdge2.strength) {
									myEdge2.strength = myEdge.strength;
								}
							} else {
								uniqueEdges.put(myEdge.toString(),myEdge);
							}
						}
					}
				}
			}
		}
		
	/*	Iterator<UndirectedEdge> edgeIter = uniqueEdges.values().iterator();
		while(edgeIter.hasNext()) {
			myEdge = edgeIter.next();
			ec = edgeColor(myEdge.strength);
			file.write("edge\n[\n");
			file.write("source " + myEdge.source + "\n");
			file.write("target " + myEdge.destination + "\n");
			file.write("graphics [\n");
			file.write("fill \"" + ec + "\"\n");
			file.write("]\n");
			file.write("]\n");
		} */
	}
	
	String toHexadecimal(int v) {
		return Integer.toHexString(v).toUpperCase();
	}
	
	String edgeColor(double v) {
		String c = "";
		int cv = (int) Math.round((1-v)*(255.0/1.5));
		if (cv > 255) {
			cv = 255;
		}
		if (cv < 0) {
			cv = 0;
		}
		String c2 = toHexadecimal(cv);
		c = "#" + c2 + c2 + c2;
		return c;
	}
	
	int closestRelative(SingleStrongModule module,HashMap<DirichletProcess,Integer> DPMap) {
		int relative = -1;
		
		Iterator<Document> iter = module.controllingDocuments.iterator();
		Document doc = null;
		HashSet<Integer> relatives = new HashSet<Integer>();
		HashSet<Integer> myAncestors = null;
		
		while(iter.hasNext()) {
			doc = iter.next();
			myAncestors = ancestors(doc,DPMap);
			if (relatives.isEmpty()) {
				relatives.addAll(myAncestors);
			} else {
				relatives.retainAll(myAncestors);
			}
		}
		
		Iterator<Integer> relativeIter = relatives.iterator();
		int curRelative = 0;
		while (relativeIter.hasNext()) {
			curRelative = relativeIter.next();
			if (curRelative > relative) {
				relative = curRelative;
			}
		}
		
		return relative;
	}
	
	HashSet<Integer> ancestors(Document doc,HashMap<DirichletProcess,Integer> DPMap) {
		HashSet<Integer> myAncestors = new HashSet<Integer>();
		
		DirichletProcess parent = null;
		parent = doc.parent.parent;
		
		while(parent != null) {
			myAncestors.add(DPMap.get(parent));
			parent = parent.parent;
		}
		
		return myAncestors;
	}
	
	// create strong modules from existing modules
	ArrayList<StrongModule> createStrongModules() {
		int i = 0;
  		int j = 0;
  		int k = 0;
  		double p = 0;
  		DirichletProcess dp = null;
  		Document[] docs = null;
  		Document doc = null;
  		SingleStrongModule[] modules = new SingleStrongModule[myHDP.numClusters];
  		
  		myHDP.numClustersSamples.clear();
  		myHDP.numClustersSamples.add(1);
  		
  		for (i=0;i<myHDP.numClusters;i++) {
  			modules[i] = new SingleStrongModule(myHDP.totalGenes,myHDP.posCounts[i]);
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
	  						if (p <= myHDP.documentControlPvalThresh) {
						//	if (p < 0.1) {
	  							modules[i].addControllingDocument(doc,myHDP.tempClusterMap[i],myHDP.tempClusterArray[i]);
	  						}
						}
					}
				}
  			}
		}
  		
  		ArrayList<StrongModule> moduleArray = new ArrayList<StrongModule>();
  		for (i=0;i<myHDP.numClusters;i++) {
  			if (modules[i].controllingDocuments.size() >= 1 & modules[i].totalGenesInModule >= myHDP.minGenesInModule) {
  				moduleArray.add(modules[i]);
  			}
  		}
  		return moduleArray;
	}
	
	// merge together modules w/ overlapping genes
	ArrayList<StrongModule> collapseModules() {
		int i = 0;
		int j = 0;
		int nc = myHDP.strongModuleMap.size();
		ArrayList<StrongModule> modules = new ArrayList<StrongModule>();
		ArrayList<HashSet<Integer>> geneSets = new ArrayList<HashSet<Integer>>();
		HashSet<Integer> genes = null;
		Iterator<StrongModule> iter = myHDP.strongModuleMap.values().iterator();
		SingleStrongModule module = null;
		double v1 = 0.0;
		double v2 = 0.0;
		double norm = 0.0;
		
		while(iter.hasNext()) {	
			module = (SingleStrongModule) iter.next();
			v2 = ((double) module.numOccur)/((double) myHDP.numClustersSamples.size());
			
			if (v2 >= minModuleOccurPercent) {
				genes = new HashSet<Integer>();
				norm = (double) module.numOccur*module.controllingDocuments.size();
				
				for (j=0;j<module.geneCounts.length;j++) {
					v1 = ((double) module.geneCounts[j])/norm;
					if (v1 >= minGeneOccurPercent) {
						genes.add(j);
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
		
		ArrayList<SingleStrongModule> modules2 = new ArrayList<SingleStrongModule>();
		SingleStrongModule module2 = null;
		
		int moduleNum = 0;
		for (i=0;i<hierCluster.clusters.size();i++) {
			moduleNum = hierCluster.clusters.get(i).get(0);
			module = (SingleStrongModule) modules.get(moduleNum);
			modules2.add(module);
			if (hierCluster.clusters.get(i).size() > 1) {
				for (j=1;j<hierCluster.clusters.get(i).size();j++) {
					moduleNum = hierCluster.clusters.get(i).get(j);
					module2 = (SingleStrongModule) modules.get(moduleNum);
					module.controllingDocuments.retainAll(module2.controllingDocuments);
				}
				for (j=1;j<hierCluster.clusters.get(i).size();j++) {
					moduleNum = hierCluster.clusters.get(i).get(j);
					module2 = (SingleStrongModule) modules.get(moduleNum);
					module.merge(module2);
				}
			}
		}
		
		ArrayList<StrongModule> modules3 = new ArrayList<StrongModule>();
		// perform a final filtering on merged modules
		for (i=0;i<modules2.size();i++) {
			module = modules2.get(i);
			v2 = ((double) module.numOccur)/((double) myHDP.numClustersSamples.size());
			if (v2 >= finalFilterModuleOccurPercent) {
				modules3.add(module);
			}
		}
		
		return modules3;
	}
	
}
