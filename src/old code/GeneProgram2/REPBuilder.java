package edu.mit.csail.psrg.georg.GeneProgram2;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Vector;

import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;
import edu.mit.csail.psrg.georg.GO.DAGException;
import edu.mit.csail.psrg.georg.GO.GoOBOReader;
import edu.mit.csail.psrg.georg.GO.MouseHumanGoAssocReader;
import edu.mit.csail.psrg.georg.GO.OntClusters;
import edu.mit.csail.psrg.georg.GO.OntDAG;
import edu.mit.csail.psrg.georg.GO.OntGeneAssociations;
import edu.mit.csail.psrg.georg.GO.OntTerm;

// this class will take a series of .persist files,
// extract EPs, and then merge EPs from subsequent samples
// to create a coherent set of REPs
public class REPBuilder {
	public ArrayList<REPOverlap> overlaps = new ArrayList<REPOverlap>();
	public int minGenesInREP = 5;
	// the percentage of the REP's genes that must be used by the tissue
	// in a sample to be deemed significant
	public double tissueUseGenesThreshold = 0.05;
	// minimum similarity to merge REPs
	public double minSimilarity = 0.50;
	// minimum percentage occurance for a gene in a REP
	public double minGeneOccurPercent = 0.05;
	// the minimum number of samples the REP must occur in
	public double minREPOccurPercent = 0.50;
	// the minimum percentage of samples the REP must be significant in
	public double minREPPercentSignificance = 0.50;
	
	// the minimum similarity to merge tissues into the same consensus group
	public double minGroupMerge = 0.90;
	
	// parameters for GO categories
	public int minGOGenes = 5;
	public int maxGOGenes = 200;
	public double FDR = 0.05;
	public int namespace = OntTerm.BiologicalProcess;
	public HashSet<String> useGenes = new HashSet<String>();
	public String GOAssocFileName = "C:\\research_data\\mouse_human\\b47mm7\\Hs.genelist.go.plus";
	public String GOHierarchyFileName = "C:\\research_data\\go_data\\gene_ontology.obo";
	
	private class geneHolder implements Comparable {
		String name = "";
		double rank = 0.0;
		double occur = 0.0;
		
		public geneHolder(String n,double r,double o) {
			name = n;
			rank = r;
			occur = o;
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
		
		void addGene(String name,double rank,double occur) {
			genes.add(new geneHolder(name,rank,occur));
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
	
	public void processSamples(String baseFile,String outfName,String modifierFileName, int sampleStart,int sampleInterval,int numSamples) {
		int i = 0;
		String fName = null;
		int sampleID = 0;
		HDP myHDP = null;
		
		fName = baseFile + ((Integer) sampleStart).toString();
		fName = fName + ".persist";
		myHDP = HDP.restoreFromFile(fName);
		LinkedHashMap<String,Integer> tissueMap = new LinkedHashMap<String,Integer>();
		int tissueNum = 0;
		for (i=0;i<myHDP.DP.length;i++) {
			if (myHDP.DP[i].getClass() == TissueDP.class) {
				tissueMap.put(myHDP.DP[i].label,tissueNum);
				tissueNum++;
			}
		} 
		
		for (i=0;i<numSamples;i++) {
			sampleID = sampleStart + i*sampleInterval;
			fName = baseFile + ((Integer) sampleID).toString();
			fName = fName + ".persist";
			myHDP = HDP.restoreFromFile(fName);
			System.out.println("Processing sample " + ((Integer) sampleID).toString());
			addSample(myHDP,tissueMap);
		}
		
		try {
			outputOverlapREPs(outfName,myHDP,numSamples,tissueMap,modifierFileName);
			// use last HDP to output number of expression programs for samples
			myHDP.outputNumExpressionProgramSamples(outfName+"_num_EPs.txt");
			// output consensus tissue groups if groups have been collected
			if (myHDP.DPGroups.size() > 0) {
				if (myHDP.DPGroups.get(0).numGoodIter > 0) {
					reconstructGroups(myHDP,outfName+"_groups.txt");
				}
			}
		} catch(IOException e) {
			System.out.println(e);
		}
		
		System.out.println("Done");
	}
	
	public void addSample(HDP myHDP,LinkedHashMap<String,Integer> tissueMap) {
		int i = 0;
  		int j = 0;
  		int k = 0;
  		double p = 0;
  		DirichletProcess dp = null;
  		REPOverlap[] tempREPs = new REPOverlap[myHDP.numExpressionPrograms];
  		int[] modifiers = null;
  		boolean UDChoice = true;
  		
  		if (myHDP.modifierLevels != null) {
  			modifiers = new int[myHDP.modifierLevels.length];
  		}
  		
  		for (i=0;i<myHDP.numExpressionPrograms;i++) {
  			tempREPs[i] = new REPOverlap(myHDP.totalGenes);
  			tempREPs[i].addGeneCounts(myHDP.expressionPrograms.get(i).posCounts);
  		}
  		
  		for (j=0;j<myHDP.DP.length;j++) {
  			dp = myHDP.DP[j];
  			if (dp.state != DirichletProcess.HELDOUT & dp.getClass() == TissueDP.class) {
				((TissueDP) dp).allocateGenesToExpressionProgramMap();
				for (i=0;i<myHDP.numExpressionPrograms;i++) {
					if (modifiers != null) {
						for (k=0;k<modifiers.length;k++) {
							modifiers[k] = ((TissueDP) dp).modifierExpressionProgramChoice[k][i];
						}
					}
					if (myHDP.useUpDown) {
						UDChoice = ((TissueDP) dp).upDownExpressionProgramChoice[i];
					}
					tempREPs[i].addTissue(myHDP,UDChoice,modifiers,myHDP.expressionPrograms.get(i).expressionProgramMap,myHDP.expressionPrograms.get(i).expressionProgramArray,(TissueDP) dp,minGenesInREP,tissueUseGenesThreshold);
				}
			}
		}
  		
  		for (i=0;i<myHDP.numExpressionPrograms;i++) {
  			tempREPs[i].normCounts();
  		}
  		
  		double[] bestSimilarity = new double[tempREPs.length];
  		int[] bestSimilarityIDX = new int[tempREPs.length];
  		double similarity = 0.0;
  		
  		for (i=0;i<tempREPs.length;i++) {
  			bestSimilarity[i] = 0;
  			bestSimilarityIDX[i] = -1;
  			for (j=0;j<overlaps.size();j++) {
  				similarity = overlaps.get(j).similarity(tempREPs[i]);
  				if (similarity > bestSimilarity[i] & similarity >= minSimilarity) {
  					bestSimilarity[i] = similarity;
  					bestSimilarityIDX[i] = j;
  				}
  			}
  		}
  		
  		for (i=0;i<tempREPs.length;i++) {
  			if (bestSimilarityIDX[i] >= 0) {
  				overlaps.get(bestSimilarityIDX[i]).add(tempREPs[i]);
  			} else {
  				overlaps.add(tempREPs[i]);
  			}
  		}
	}
	
	public void outputOverlapREPs(String fOutName,HDP myHDP,int numSamples,LinkedHashMap<String,Integer> tissueMap,String modifierFileName) throws IOException {
		ArrayList<String> genes = new ArrayList<String>();
		int i = 0;
		int j = 0;
		int k = 0;
		GoOBOReader goReader = new GoOBOReader();
		goReader.OBOFileName = GOHierarchyFileName;
		MouseHumanGoAssocReader humanAssocReader = new MouseHumanGoAssocReader();
		humanAssocReader.setFile(GOAssocFileName);
		OntDAG GODAG = null;
		OntGeneAssociations GOAssoc = null;
		
		String[] humanGeneNames = myHDP.geneNames;
		
		HashSet<String> humanUseGenes = new HashSet<String>();
		for (i=0;i<humanGeneNames.length;i++) {
			humanUseGenes.add(humanGeneNames[i]);
		}
		
		try {
			GODAG = goReader.readFile();
			System.out.println("Loaded GO DAG");
			GOAssoc = humanAssocReader.readFile(GODAG,humanUseGenes);
			System.out.println("Loaded GO associations");
		} catch(IOException e) {
			System.out.println(e);
		} catch(DAGException e) {
			System.out.println(e);
		}
		
		OntClusters clusters = new OntClusters(GOAssoc,minGOGenes,maxGOGenes,FDR,namespace);
		
		geneList glist = null;
		double[] ranks = null;
		double norm = 0.0;
		REPOverlap rep = null;
		double v1 = 0.0;
		double v2 = 0.0;
		Vector<OverlapTissue> tissues = new Vector<OverlapTissue>();
		Iterator<OverlapTissue> tissueIter = null;
		OverlapTissue tissue = null;
		
		FileWriter file = new FileWriter(fOutName+"_REPs.txt");
		ArrayList<REPOverlap> useREPs = new ArrayList<REPOverlap>();
		int gene = 0;
		
		int numGoodREPs = 0;
		for (i=0;i<overlaps.size();i++) {
			rep = overlaps.get(i);
			v1 = ((double) rep.numOccur)/((double) numSamples);
			if (v1 >= minREPOccurPercent) {
				numGoodREPs++;
			}
		}
		
		for (i=0;i<overlaps.size();i++) {
			rep = overlaps.get(i);
			v1 = ((double) rep.numOccur)/((double) numSamples);
			if (v1 >= minREPOccurPercent) {
				tissues.clear();
				tissueIter = rep.tissues.values().iterator();
				while (tissueIter.hasNext()) {
					tissue = tissueIter.next();
					v2 = ((double) tissue.numSignificant)/((double) rep.numOccur);
					if (v2 >= minREPPercentSignificance) {
						tissues.add(tissue);
					}
				}
				if (tissues.size() >= 1) {
					genes.clear();
					glist = new geneList();
					for (j=0;j<rep.geneOccurNormalized.length;j++) {
						v1 = rep.geneOccurNormalized[j]/((double) rep.numOccur);
						v2 = rep.geneIntensityNormalized[j]/((double) rep.numOccur);
						if (v1 >= minGeneOccurPercent) {
							gene = myHDP.reverseGeneMap[j];
							glist.addGene(humanGeneNames[gene],v2,v1);
						}
					}
					
					if (glist.size() >= minGenesInREP) {
						useREPs.add(rep);
						ranks = glist.sort(genes);
						clusters.addCluster(genes,ranks,null,1.0);
					}
				}
			}
		}
		
		clusters.clusterSignif();
		
		ArrayList<ArrayList<String> > modifierNames = null;
		int[][] modifiers = null;
		int mi = 0;
		int mj = 0;
		if (modifierFileName != null) {
			ClusterReader modifierReader = new ClusterReader();
			modifierReader.readFile(modifierFileName);
			modifierNames = modifierReader.clusters;
		}
		
		for (i=0;i<useREPs.size();i++) {
			rep = useREPs.get(i);
			v1 = ((double) rep.numOccur)/((double) numSamples);
			
			tissues.clear();
			tissueIter = rep.tissues.values().iterator();
			while (tissueIter.hasNext()) {
				tissue = tissueIter.next();
				v2 = ((double) tissue.numSignificant)/((double) rep.numOccur);
				if (v2 >= minREPPercentSignificance) {
					tissues.add(tissue);
				}
			}
			Collections.sort(tissues);
			file.write("REP#" + ((new Integer(i+1))).toString());
			file.write("\t" + String.format("%.3f",v1));
			file.write("\n");
			file.write("tissues");
			for (j=0;j<tissues.size();j++) {
				file.write("\t");
				file.write(tissues.get(j).tissueName);
			}
			file.write("\n");
			
			if (myHDP.useUpDown) {
				file.write("UD");
				for (j=0;j<tissues.size();j++) {
					file.write("\t");
					v1 = (double) tissues.get(j).numUp;
					v1 = v1 / ((double) rep.numOccur);
					if (v1 < 0.50) {
						v1 = v1 - 1.0;
					}
					file.write(String.format("%.3f",v1));
				}
				file.write("\n");
			}
			
			file.write("topicUse");
			for (j=0;j<tissues.size();j++) {
				file.write("\t");
				v1 = tissues.get(j).topicUse;
				v1 = v1 / ((double) rep.numOccur);
				file.write(String.format("%.3f",v1));
			}
			file.write("\n");
			file.write("tissueUse");
			for (j=0;j<tissues.size();j++) {
				file.write("\t");
				v1 = tissues.get(j).tissueUse;
				v1 = v1 / ((double) rep.numOccur);
				file.write(String.format("%.3f",v1));
			}
			file.write("\n");
			file.write("signif");
			for (j=0;j<tissues.size();j++) {
				file.write("\t");
				v1 = (double) tissues.get(j).numSignificant;
				v1 = v1 / ((double) rep.numOccur);
				file.write(String.format("%.3f",v1));
			}
			file.write("\n");
			
			if (modifierFileName != null) {
				for (mi=0;mi<modifierNames.size();mi++) {
					for (mj=0;mj<modifierNames.get(mi).size();mj++) {
						file.write(modifierNames.get(mi).get(mj));
						for (j=0;j<tissues.size();j++) {
							v1 = (double) tissues.get(j).numModifier[mi][mj];
							v1 = v1 / ((double) rep.numOccur);
							file.write("\t");
							file.write(String.format("%.3f",v1));
						}
						file.write("\n");
					}
				}
				file.write("\n");
			}
			
			clusters.clusters.get(i).outputScoresLong(file,true);
		}	
		file.close();
	}
	
	public void reconstructGroups(HDP myDP,String groupFName) {
		if (myDP.DPGroups.size() == 0)
			return;
		
		int i = 0;
		int j = 0;
		for (i=0;i<myDP.DPGroups.size();i++) {
			if (myDP.DPGroups.get(i).numGoodIter > 0) {
				myDP.DPGroups.get(i).normPairProbs();
				myDP.DPGroups.get(i).buildConcensusModel(minGroupMerge);
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
}
