package edu.mit.csail.psrg.georg.HDP2;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Vector;

import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;
import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.GO.DAGException;
import edu.mit.csail.psrg.georg.GO.GoOBOReader;
import edu.mit.csail.psrg.georg.GO.MouseHumanGoAssocReader;
import edu.mit.csail.psrg.georg.GO.OntClusters;
import edu.mit.csail.psrg.georg.GO.OntDAG;
import edu.mit.csail.psrg.georg.GO.OntGeneAssociations;
import edu.mit.csail.psrg.georg.GO.OntTerm;
import edu.mit.csail.psrg.georg.HierarchicalCluster.HierarchicalCluster;
import edu.mit.csail.psrg.georg.StatUtil.StatUtil;

public class SummarizeInferredModel {
	public int minGenes = 5;
	public int maxGenes = 200;
	// minimum count of gene in topic
	public int minOccurInTopic = 2;
	public double FDR = 0.05;
	// minimum percentage occurance for a gene in a topic
//	public double minGeneOccurPercent = 0.05;
	public double minGeneOccurPercent = 1.0;
	public double minModuleOccurPercent = 0.05;
	public double minGeneOverlap = 0.50;
	public double finalFilterTopicOccurPercent = 0.50;
	public int namespace = OntTerm.BiologicalProcess;
	public HashSet<String> useGenes = new HashSet<String>();
	
	public String refDirPath = "C:\\research_data\\mouse_human\\b47mm7\\";
	public String humanAssocFile = refDirPath + "Hs.genelist.go.plus";
	
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
	
	public void outputFiles(String dirPath,String snapFileName,String classifyFileName) {
			String REPsAnnotatedFile = dirPath + "REPs_GO_annotated.txt";
			String geneCrossOutFile = dirPath + "gene_cross.txt";
			String tissueMatrixLoadFile = dirPath + "tissue_load.txt";
			String topicGOCountFile = dirPath + "topic_GO_count.txt";
			String topicCountFile = dirPath + "numtopics.txt";
			String groupFile = dirPath + "groups.txt";
			String classifyPValuesFileName = dirPath + "classify_pvals.txt";
			
			HDP myHDP = HDP.restoreFromFile(dirPath + snapFileName);
			
			int i = 0;
			int j = 0;
			
			String[] humanGeneNames = new String[myHDP.geneNames.length];
			String[] mouseGeneNames = new String[myHDP.geneNames.length];
			String s1 = "";
			String s2 = "";
			int f = 0;
			
			for (i=0;i<myHDP.geneNames.length;i++) {
		/*		f = myHDP.geneNames[i].indexOf("|");
				s1 = myHDP.geneNames[i].substring(0,f);
				s2 = myHDP.geneNames[i].substring(f+1);
				humanGeneNames[i] = s1;
				mouseGeneNames[i] = s2; */
				humanGeneNames[i] = myHDP.geneNames[i];
			}
			
			if (myHDP.DPGroups.size() > 0) {
				if (myHDP.DPGroups.size() > 0) {
					reconstructGroups(myHDP,groupFile);
				}
			}
			
			outputREPs(myHDP,humanGeneNames,REPsAnnotatedFile,classifyFileName,classifyPValuesFileName);
			
			outputTissueLoadMatrixDirect(myHDP,tissueMatrixLoadFile);

			try {
				myHDP.outputNumTopicsSamples(topicCountFile);
			} catch(IOException e) {
				System.out.println(e);
			}
		}
		
		public static void reconstructGroups(HDP myDP,String groupFName) {
			if (myDP.DPGroups.size() == 0)
				return;
			
			int i = 0;
			int j = 0;
			for (i=0;i<myDP.DPGroups.size();i++) {
				if (myDP.DPGroups.get(i).numGoodIter > 0) {
					myDP.DPGroups.get(i).normPairProbs();
					myDP.DPGroups.get(i).buildConcensusModel();
				}
			}
			myDP.buildDPList();
			
			ParentDP pdp = null;
			TissueDP dp = null;
			TissueGroup group = myDP.DPGroups.get(0);
			
			try {
				FileWriter file = new FileWriter(groupFName);
				
				int dpID = 0;
				for (i=0;i<group.controllingProcess.numChildren();i++) {
					pdp = (ParentDP) group.controllingProcess.children.get(i);
					for (j=0;j<pdp.numChildren();j++) {
						dp = (TissueDP) pdp.children.get(j);
						file.write(dp.label);
						if (j < pdp.numChildren() - 1) {
							file.write("\t");
						}
					}
					file.write("\n");
				}
				file.close();
			} catch(IOException e) {
				System.out.println(e);
			}
			
			System.out.println("reconstructed groups");
		}
		
		// create REPs from single sample
		public ArrayList<REP> createREPsFromSample(HDP myHDP) {
			int i = 0;
	  		int j = 0;
	  		double p = 0;
	  		DirichletProcess dp = null;
	  		TissueDP tissue = null;
	  		REP[] reps = new REP[myHDP.numTopics];
	  		
	  		REPCollector myREPCollector = (REPCollector) myHDP.getStatCollector(REPCollector.class);
	  		if (myREPCollector == null) {
	  			myREPCollector = new REPCollector(myHDP);
	  		}
	  		
	  		for (i=0;i<myHDP.numTopics;i++) {
	  			reps[i] = new GenericREP(myHDP.totalGenes,myHDP.topics.get(i).posCounts);
	  		}
	  		
	  		for (j=0;j<myHDP.DP.length;j++) {
	  			dp = myHDP.DP[j];
	  			if (dp.state != DirichletProcess.HELDOUT & dp.getClass() == TissueDP.class) {
	  				tissue = (TissueDP) dp;
					tissue.allocateGenesToTopicMap();
					for (i=0;i<myHDP.numTopics;i++) {
						p = reps[i].controlPVal(myHDP.topics.get(i).topicMap,tissue,myHDP.totalGenes,minGenes);
		  				if (p <= myREPCollector.tissueControlPvalThresh) {
		  					reps[i].addControllingTissue(tissue,myHDP.topics.get(i).topicMap,myHDP.topics.get(i).topicArray,tissue.upDownTopicChoice[i]);
		  				}
					}
	  			}
			}
	  		
	  		ArrayList<REP> REPArray = new ArrayList<REP>();
	  		for (i=0;i<myHDP.numTopics;i++) {
	  			if (reps[i].controllingTissues.size() >= 1 & reps[i].totalGenesInREP >= myREPCollector.minGenesInREP) {
	  				REPArray.add(reps[i]);
	  			}
	  		}
	  		return REPArray;
		}
		
		public LinkedHashMap<String,String> loadClassifications(String classifyFileName) {
			LinkedHashMap<String,String> classifications = new LinkedHashMap<String,String>();
			ClusterReader reader = new ClusterReader();
			
			try {
				reader.readFile(classifyFileName);
			} catch(IOException e) {
				System.out.println(e);
			}
			
			int i = 0;
			String type = null;
			for (i=0;i<reader.clusters.size();i++) {
				type = reader.clusters.get(i).get(4);
				if (type.equals("D")) {
					classifications.put(reader.clusters.get(i).get(0),reader.clusters.get(i).get(1));
				}
			}
			
			return classifications;
		}
		
		public void outputREPs(HDP myHDP,String[] humanGeneNames,String goOutName,String classifyFileName,String classifyPValuesFileName) {
			int numSamples = 0;
			
			int j = 0;
			int i = 0;
			int k = 0;
			
			LinkedHashMap<String,String> classifications = null;
			LinkedHashMap<String,Integer> classCountMap = null;
			LinkedHashMap<String,Integer> classMap = null;
			String[] classNames = null;
			String tempS = null;
			Iterator<String> siter = null;
			int[] classCountUP = null;
			int[] classCountDN = null;
			double[][] classUPPValues = null;
			double[][] classDNPValues = null;
			int classID = 0;
			int pvx = 0;
			int pvN = 0;
			int pvs = 0;
			int pvn = 0;
			int totalUP = 0;
			int totalDN = 0;
			if (classifyFileName != null) {
				classifications = loadClassifications(classifyFileName);
				classCountMap = new LinkedHashMap<String,Integer>();
				classMap = new LinkedHashMap<String,Integer>();
				siter = classifications.values().iterator();
				while(siter.hasNext()) {
					tempS = siter.next();
					if (!classCountMap.containsKey(tempS)) {
						classCountMap.put(tempS,0);
					}
					i = classCountMap.get(tempS);
					classCountMap.put(tempS,i+1);
				}
				classCountUP = new int[classCountMap.size()];
				classCountDN = new int[classCountMap.size()];
				classNames = new String[classCountMap.size()];
				
				siter = classCountMap.keySet().iterator();
				i = 0;
				while(siter.hasNext()) {
					tempS = siter.next();
					classMap.put(tempS,i);
					classNames[i] = tempS;
					i++;
				}
			}
			
			REPCollector ms = null;
			if (myHDP.statCollectors.containsKey(REPCollector.class)) {
				ms = (REPCollector) myHDP.getStatCollector(REPCollector.class);
				if (ms.numSamples > 0) {
					numSamples = ms.numSamples;
				} else {
					ms = null;
					numSamples = 1;
				}
			}
			
			ArrayList<String> genes = new ArrayList<String>();
			
			GoOBOReader goReader = new GoOBOReader();
			MouseHumanGoAssocReader humanAssocReader = new MouseHumanGoAssocReader();
			humanAssocReader.setFile(humanAssocFile);
			OntDAG humanDAG = null;
			OntGeneAssociations humanAssoc = null;
			
			HashSet<String> humanUseGenes = new HashSet<String>();
			for (i=0;i<humanGeneNames.length;i++) {
				humanUseGenes.add(humanGeneNames[i]);
			}
			
			try {
				humanDAG = goReader.readFile();
				System.out.println("Loaded human GO DAG");
				humanAssoc = humanAssocReader.readFile(humanDAG,humanUseGenes);
				System.out.println("Loaded human GO associations");
			} catch(IOException e) {
				System.out.println(e);
			} catch(DAGException e) {
				System.out.println(e);
			}
			
			OntClusters clusters = new OntClusters(humanAssoc,minGenes,maxGenes,FDR,namespace);
			
			geneList glist = null;
			double[] ranks = null;
			double norm = 0.0;
			
			ArrayList<REP> reps = null;
			
			if (ms != null) {
				reps = collapseREPs(ms);
			} else {
				reps = createREPsFromSample(myHDP);
			}
			
			if (classifyFileName != null) {
				classUPPValues = new double[reps.size()][classCountMap.size()];
				classDNPValues = new double[reps.size()][classCountMap.size()];
				for (i=0;i<classUPPValues.length;i++) {
					for (j=0;j<classUPPValues[0].length;j++) {
						classUPPValues[i][j] = 1.0;
						classDNPValues[i][j] = 1.0;
					}
				}
			}
			
			Iterator<REP> iter = reps.iterator();
			
			Vector<ControllingTissue> tissues = null;
			GenericREP rep = null;
			double v1 = 0.0;
			double v2 = 0.0;
			String description = "";
			ControllingTissue tissue = null;
			int gene = 0;
			Iterator<Integer> intIter = null;
			int geneIDX = 0;
			while(iter.hasNext()) {
				genes.clear();
				glist = new geneList();
				rep = (GenericREP) iter.next();
				norm = (double) rep.numOccur*rep.controllingTissues.size();
				intIter = rep.geneCounts.keySet().iterator();
				while(intIter.hasNext()) {
					geneIDX = intIter.next();
					v1 = ((double) rep.geneCounts.get(geneIDX))/norm;
					v2 = ((double) rep.geneIntensity.get(geneIDX))/norm;
				//	if (v1 >= minGeneOccurPercent) {
					if (v2 >= minGeneOccurPercent) {
						gene = myHDP.reverseGeneMap[geneIDX];
						glist.addGene(humanGeneNames[gene],v2);
					}
				}
				
				description = "";
				tissues = rep.sortTissues();
				
				if (classifyFileName != null) {
					for (i=0;i<classCountUP.length;i++) {
						classCountUP[i] = 0;
						classCountDN[i] = 0;
					}
					
					totalUP = 0;
					totalDN = 0;
					for (i=0;i<tissues.size();i++) {
						tempS = classifications.get(tissues.get(i).dp.label);
						classID = classMap.get(tempS);
						if (tissues.get(i).UDuse) {
							totalUP++;
						} else {
							totalDN++;
						}
						if (tissues.get(i).UDuse) {
							classCountUP[classID]++;
						} else {
							classCountDN[classID]++;
						}
					}
					
					for (i=0;i<classCountUP.length;i++) {
						pvx = classCountUP[i];
						pvN = classifications.size();
						pvs = classCountMap.get(classNames[i]);
						pvn = totalUP;
						// x = # observed successes in the sample
						// N = population size
						// s = # of successes in population (e.g., size of positive set in the population)
						// n = sample size	
						if (pvx >= 4) {
							classUPPValues[k][i] = 1.0-StatUtil.hyperGeometricCDF(pvx-1,pvN,pvs,pvn);
						}
						pvx = classCountDN[i];
						pvn = totalDN;
						if (pvx >= 4) {
							classDNPValues[k][i] = 1.0-StatUtil.hyperGeometricCDF(pvx-1,pvN,pvs,pvn);
						}
						
					}
					k++;
				}
				
				for(i=0;i<tissues.size();i++) {
					tissue = tissues.get(i);
					v1 = tissue.load/((double) rep.numOccur);
					if (!tissue.UDuse) {
						v1 = -v1;
					}
					description = description + tissue.sID + " (" + (new Double(v1)).toString() + ") ";
					description = description + " | ";
				}
				ranks = glist.sort(genes);
				v1 = ((double) rep.numOccur)/((double) numSamples);
				if (v1 >= minModuleOccurPercent & glist.size() >= minGenes) {
					clusters.addCluster(genes,ranks,description,v1);
				}
			}
			
			clusters.clusterSignif();
			
			try {
				clusters.outputClustersLong(goOutName);
				if (classifyFileName != null) {
					outputClassifyPValues(classifyPValuesFileName,classUPPValues,classDNPValues,classNames);
				}
			} catch(IOException e) {
				System.out.println(e);
			}
		}
		
		public static void outputClassifyPValues(String classifyPValuesFileName,double[][] UPpvals,double[][] DNpvals,String[] classNames) throws IOException {
			FileWriter file = new FileWriter(classifyPValuesFileName);
			int i = 0;
			int j = 0;
			
			file.write("GP\t");
			
			for (i=0;i<classNames.length;i++) {
				file.write(classNames[i]);
				file.write("_UP\t");
				file.write(classNames[i]);
				file.write("_DN");
				if (i < classNames.length - 1) {
					file.write("\t");
				} else {
					file.write("\n");
				}
			}
			
			for (i=0;i<UPpvals.length;i++) {
				file.write((new Integer(i+1)).toString());
				file.write("\t");
				for (j=0;j<UPpvals[0].length;j++) {
					file.write((new Double(UPpvals[i][j])).toString());
					file.write("\t");
					file.write((new Double(DNpvals[i][j])).toString());
					if (j < UPpvals[0].length - 1) {
						file.write("\t");
					} else {
						file.write("\n");
					}
				}
			}
			
			file.close();
		}
		
		public static void outputTissueLoadMatrixDirect(HDP myHDP,String tissueMatrixLoadFile) {
			double[][] values = new double[myHDP.DP.length][myHDP.numTopics];
			String[] DPNames = new String[myHDP.DP.length];
			String[] GPNames = new String[myHDP.numTopics];
			
			int dd = 0;
			int tt = 0;
			
			for (tt=0;tt<myHDP.numTopics;tt++) {
				GPNames[tt] = "GP#" + (tt+1);
			}
			
			for (dd=0;dd<myHDP.DP.length;dd++) {
				DPNames[dd] = myHDP.DP[dd].getLabel();
				for (tt=0;tt<myHDP.numTopics;tt++) {
					values[dd][tt] = myHDP.DP[dd].beta[tt];
					if (myHDP.DP[dd].getClass() == TissueDP.class) {
						if (!((TissueDP) myHDP.DP[dd]).upDownTopicChoice[tt]) {
							values[dd][tt] = -values[dd][tt];
						}
					}
				}
			}
			
			MicroArrayData matrix = new MicroArrayData();
			matrix.setContinuous();
			matrix.values = values;
			matrix.numRows = myHDP.DP.length;
			matrix.numCols = myHDP.numTopics;
			matrix.geneNames = DPNames;
			matrix.experimentNames = GPNames;
			
			try {
				matrix.writeFile(tissueMatrixLoadFile);
			} catch(IOException e) {
				System.out.println(e);
			}
		}
		
//		 merge together modules w/ overlapping genes
		ArrayList<REP> collapseREPs(REPCollector myCollector) {
			int i = 0;
			int j = 0;
			
			int nc = myCollector.REPMap.size();
			
			ArrayList<REP> reps = new ArrayList<REP>();
			ArrayList<HashSet<Integer>> geneSets = new ArrayList<HashSet<Integer>>();
			HashSet<Integer> genes = null;
			Iterator<REP> iter = myCollector.REPMap.values().iterator(); 
			GenericREP rep = null;
			double v1 = 0.0;
			double v2 = 0.0;
			double norm = 0.0;
			int geneID = 0;
			int count = 0;
			Iterator<Integer> gIter = null;
			
			while(iter.hasNext()) {	
				rep = (GenericREP) iter.next();
				v2 = ((double) rep.numOccur)/((double) myCollector.numSamples);
				
			//	if (v2 >= minModuleOccurPercent) {
					genes = new HashSet<Integer>();
					gIter = rep.geneCounts.keySet().iterator();
					norm = (double) rep.numOccur*rep.controllingTissues.size();
					while (gIter.hasNext()) {
						geneID = gIter.next();
						count = rep.geneCounts.get(geneID);
						v1 = ((double) count)/norm;
						if (v1 >= minGeneOccurPercent) {
							genes.add(geneID);
						}
					}
					if (genes.size() >= minGenes) {
						geneSets.add(genes);
						reps.add(rep);
					}
			//	}
			}
			
			nc = reps.size();
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
			
			ArrayList<REP> reps2 = new ArrayList<REP>();
			GenericREP rep2 = null;
			
			int repNum = 0;
			for (i=0;i<hierCluster.clusters.size();i++) {
				repNum = hierCluster.clusters.get(i).get(0);
				rep = (GenericREP) reps.get(repNum);
				reps2.add(rep);
				if (hierCluster.clusters.get(i).size() > 1) {
					for (j=1;j<hierCluster.clusters.get(i).size();j++) {
						repNum = hierCluster.clusters.get(i).get(j);
						rep2 = (GenericREP) reps.get(repNum);
						rep.intersectControllingTissues(rep2);
					}
					for (j=1;j<hierCluster.clusters.get(i).size();j++) {
						repNum = hierCluster.clusters.get(i).get(j);
						rep2 = (GenericREP) reps.get(repNum);
						rep.merge(rep2);
					}
				}
			}
			
			ArrayList<REP> reps3 = new ArrayList<REP>();
			// perform a final filtering on merged modules
			for (i=0;i<reps2.size();i++) {
				rep = (GenericREP) reps.get(i);
				v2 = ((double) rep.numOccur)/((double) myCollector.numSamples);
				if (v2 >= finalFilterTopicOccurPercent) {
					reps3.add(rep);
				}
			}
			
			return reps3;
		}
}
