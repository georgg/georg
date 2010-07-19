package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;

import edu.mit.csail.psrg.georg.HierarchicalCluster.HierarchicalCluster;
import edu.mit.csail.psrg.georg.StatUtil.StatUtil;

public class StrongModuleTransposeCollector extends HDPStatCollector {
	
	// a list of processes to track strong modules for
	ArrayList<DirichletProcess> clusterDPList = null;
	// thresh-hold for p-value to determine that a motif controls a module
	double motifControlPvalThresh = 0.001;
	// minimum number of genes that use the module
	int minGenesInModule = 5;
	// parameters below are used for merged modules
	// minimum percentage occurance for a gene in a module
	double minClusterLoad = 0.01;
	double minModuleOccurPercent = 0.05;
	double minClusterOverlap = 0.75;
	double finalFilterModuleOccurPercent = 0.50;
	
	// map the gene ID number to the gene DP
	HashMap<Integer,DirichletProcess> geneDPMap = new HashMap<Integer,DirichletProcess>();
	
	// total # of genes that each motif occurs in (length = # of motifs)
	int[] totalMotifGenes = null;
	// total # of events for each motif (length = # of motifs)
	int[] totalMotifEvents = null;
	
	// total number of genes in each module (e.g., that have 1 or more motif event using the module) length = # modules
	ArrayList<HashSet<Integer>> totalGenesModule = new ArrayList<HashSet<Integer>>();
	// # of modules to expand the totalGenesModule array by
	int expandTotalGenesModule = 10;
	// stores gene counts for each module (by motif and cluster)
	ArrayList<HashMap<Integer,HashMap<Integer,HashMap<Integer,Integer>>>> genesInModuleByMotif = new ArrayList<HashMap<Integer,HashMap<Integer,HashMap<Integer,Integer>>>>();
	// temporary vector for computing the cluster loading or avg. # of motif occurances per gene across the cluster that load on the topic (length = # of clusters; those in DPList)
	int[] clusterTotals = null;
	int totalNumGenes = 0;
	int totalNumMotifs = 0;
	
	// number of samples (e.g., number of times stats have been updated)
	int numSamples = 0;
	
	// maps sets of motifs to modules (motifs are indexed by integers)
	LinkedHashMap<HashSet<Integer>,StrongModuleTranspose> strongModuleMap = new LinkedHashMap<HashSet<Integer>,StrongModuleTranspose>();
	// maximum number of strong modules that will be maintained
	int maxNumStrongModules = 2500;
	// interval to "clean" the number of strong modules
	int cleanStrongModulesInterval = 500;
	int cleanStrongModulesCounter = 0;

	public StrongModuleTransposeCollector(HierDirichletProcess h) {
		super(h);
	}
	
	public StrongModuleTransposeCollector(HierDirichletProcess h,ArrayList<DirichletProcess> dpL) {
		super(h);
		clusterDPList = dpL;
		totalNumMotifs = myHDP.totalGenes;
		
		clusterTotals = new int[clusterDPList.size()];
		totalMotifGenes = new int[totalNumMotifs];
		totalMotifEvents = new int[totalNumMotifs];
		
		int i = 0;
		int j = 0;
		int gg = 0;
		Document doc = null;
		DirichletProcess parentDP = null;
		DirichletProcess dp = null;
		
		for (i=0;i<clusterDPList.size();i++) {
			parentDP = clusterDPList.get(i);
			for (j=0;j<parentDP.children.size();j++) {
				dp = parentDP.children.get(j);
				geneDPMap.put(dp.documents[0].docID,dp);
				totalNumGenes++;
			}
		}
		
		// the array is indexed motif X genes
		int[][] motifGenes = new int[myHDP.totalGenes][totalNumGenes];
		int docCode = 0;
		
		for (j=0;j<myHDP.totalGenes;j++) {
			motifGenes[j] = new int[totalNumGenes];
		}
		
		for (i=0;i<clusterDPList.size();i++) {
			parentDP = clusterDPList.get(i);
			for (j=0;j<parentDP.children.size();j++) {
				dp = parentDP.children.get(j);
				doc = dp.documents[0];
				docCode = doc.code;
				for (gg=0;gg<doc.genes.length;gg++) {
					motifGenes[doc.genes[gg]][docCode] = 1;
					totalMotifEvents[doc.genes[gg]]++;
				}
			}
		}
		
		for (j=0;j<motifGenes.length;j++) {
			for (i=0;i<motifGenes[j].length;i++) {
				totalMotifGenes[j] += motifGenes[j][i];
			}
		}
	}
	
	public void allocateTotalGenesModule() {
		if (totalGenesModule == null) {
			totalGenesModule = new ArrayList<HashSet<Integer>>();
		}
		
		if (genesInModuleByMotif == null) {
			genesInModuleByMotif = new ArrayList<HashMap<Integer,HashMap<Integer,HashMap<Integer,Integer>>>>();
		}
		
		while(totalGenesModule.size() <= myHDP.numClusters + expandTotalGenesModule) {
			totalGenesModule.add(new HashSet<Integer>());
		}
		
		while(genesInModuleByMotif.size() <= myHDP.numClusters + expandTotalGenesModule) {
			genesInModuleByMotif.add(new HashMap<Integer,HashMap<Integer,HashMap<Integer,Integer>>>());
		}
	}

	public void updateStats() {
		int mm = 0;
		int i = 0;
  		int j = 0;
  		int gg = 0;
  		int cc = 0;
  		
  		allocateTotalGenesModule();
  		numSamples++;
  		
  		DirichletProcess clusterDP = null;
  		DirichletProcess dp = null;
  		Document doc = null;
  		StrongModuleTranspose[] modules = new StrongModuleTranspose[myHDP.numClusters];
  		
  		for (i=0;i<myHDP.numClusters;i++) {
  			modules[i] = new StrongModuleTranspose();
  		}
  		
  		int docID = 0;
  		int x = 0;
  		int n = 0;
  		int s = 0;
  		double p = 0.0;
  		
  		HashSet<Integer> usedGenes = new HashSet<Integer>();
  		HashMap<Integer,HashMap<Integer,HashMap<Integer,Integer>>> topicMotifMap = null;
  		HashMap<Integer,HashMap<Integer,Integer>> motifToClusters = null;
  		HashMap<Integer,Integer> clusterToGenes = null;
  		Integer oldCount = null;
  		
  		for (cc=0;cc<myHDP.numClusters;cc++) {
  			topicMotifMap = genesInModuleByMotif.get(cc);
  			topicMotifMap.clear();
  			totalGenesModule.get(cc).clear();
  		}
  		
  		for (i=0;i<clusterDPList.size();i++) {
	  		clusterDP = clusterDPList.get(i);
	  		for (j=0;j<clusterDP.children.size();j++) {
	  			dp = clusterDP.children.get(j);
	  			doc = dp.documents[0];
	  			docID = doc.code;
	  			for (gg=0;gg<doc.genes.length;gg++) {
	  				cc = doc.geneClusterAssigns[gg];
	  				mm = doc.genes[gg];
	  				topicMotifMap = genesInModuleByMotif.get(cc);
	  				motifToClusters = topicMotifMap.get(mm);
	  				if (motifToClusters == null) {
	  					motifToClusters = new HashMap<Integer,HashMap<Integer,Integer>>();
	  					topicMotifMap.put(mm,motifToClusters);
	  				}
	  				clusterToGenes = motifToClusters.get(i);
	  				if (clusterToGenes == null) {
	  					clusterToGenes = new HashMap<Integer,Integer>();
	  					motifToClusters.put(i,clusterToGenes);
	  				}
	  				
	  				oldCount = clusterToGenes.get(docID);
	  				if (oldCount == null) {
	  					clusterToGenes.put(docID,1);
	  				} else {
	  					clusterToGenes.put(docID,oldCount+1);
	  				}
	  				totalGenesModule.get(cc).add(docID);
	  			}
	  		}
		}
  		
  		Iterator<Integer> motifIter = null;
  		Iterator<Integer> clusterIDIter = null;
  		int clusterID = 0;
  		for (cc=0;cc<myHDP.numClusters;cc++) {
  			n = totalGenesModule.get(cc).size();
  			if (n > minGenesInModule) {
  				topicMotifMap = genesInModuleByMotif.get(cc);
  				motifIter = topicMotifMap.keySet().iterator();
  				while (motifIter.hasNext()) {
  					mm = motifIter.next();
  					motifToClusters = topicMotifMap.get(mm);
  					clusterIDIter = motifToClusters.keySet().iterator();
  					x = 0;
  					while(clusterIDIter.hasNext()) {
  						clusterID = clusterIDIter.next();
  						clusterToGenes = motifToClusters.get(clusterID);
  						x += clusterToGenes.size();
  					}
  					if (x > 1) {
		  				s = totalMotifGenes[mm];
		  				try {
		  					p = 1.0-StatUtil.hyperGeometricCDF(x-1,totalNumGenes,s,n);
		  				} catch(ArithmeticException e) {
		  					p = 1.0;
		  				}
		  				if (p <= motifControlPvalThresh) {
		  					modules[cc].addControllingMotif(mm,motifToClusters);
		  				}
		  			}
  				}
  			}
  		}
  		
  		for (i=0;i<myHDP.numClusters;i++) {
  			if (modules[i].controllingMotifs.size() >= 1) {
  				modules[i].norm();
  				modules[i].addToModuleMap(strongModuleMap);
  			}
  		}
  		
  		cleanStrongModulesCounter++;
  		if (cleanStrongModulesCounter >= cleanStrongModulesInterval) {
  			cleanStrongModulesCounter = 0;
  			if (strongModuleMap.size() > maxNumStrongModules) {
  				StrongModuleTranspose module = null;
  				StrongModuleTranspose[] marray = new StrongModuleTranspose[strongModuleMap.size()];
				marray = strongModuleMap.values().toArray(marray);
				List<StrongModuleTranspose> mlist = Arrays.asList(marray);
				Collections.sort(mlist);
  				int ss = strongModuleMap.size();
  				Iterator<StrongModuleTranspose> iter = mlist.iterator();
  				while(ss > maxNumStrongModules & iter.hasNext()) {
  					module = iter.next();
  					module = strongModuleMap.remove(module.controllingMotifs);
  					ss--;
  				}
  			}
  		}
	}
	
	// matrix of loadings of motifs on strong topics
	public double[][] generateMotifLoadMatrix(ArrayList<StrongModuleTranspose> modules) {
		double[][] mlm = new double[totalNumMotifs][modules.size()];
		int motifID = 0;
		StrongModuleTranspose module = null;
		Iterator<StrongModuleTranspose> moduleIter = modules.iterator();
		Iterator<Integer> motifIter = null;
		double v = 0.0;
		
		int i = 0;
		int j = 0;
		while (moduleIter.hasNext()) {
			module = moduleIter.next();
			motifIter = module.controllingMotifs.iterator();
			while (motifIter.hasNext()) {
				motifID = motifIter.next();
				v = module.motifLoad.get(motifID);
				mlm[motifID][i] = v/((double) module.numOccur);
			}
			i++;
		}
		
		// normalize the columns
		double norm = 0.0;
		for (j=0;j<mlm[0].length;j++) {
			norm = 0.0;
			for (i=0;i<mlm.length;i++) {
				norm += mlm[i][j];
			}
			if (norm > 0.0) {
				for (i=0;i<mlm.length;i++) {
					mlm[i][j] = mlm[i][j]/norm;
				}
			}
		}
		return mlm;
	}
	
	public void computeFinalClusterTotals() {
		int i = 0;
		int j = 0;
		DirichletProcess clusterDP = null;
		DirichletProcess dp = null;
		for (i=0;i<clusterDPList.size();i++) {
			clusterDP = clusterDPList.get(i);
			clusterTotals[i] = 0;
			for (j=0;j<clusterDP.children.size();j++) {
				dp = clusterDP.children.get(j);
				clusterTotals[i] += dp.documents[0].genes.length;
			}
		}
	}
	
	// matrix of loadings of clusters on strong topics
	public double[][] generateClusterLoadMatrix(ArrayList<StrongModuleTranspose> modules) {
		double[][] clm = new double[clusterDPList.size()][modules.size()];
		int clusterID = 0;
		StrongModuleTranspose module = null;
		Iterator<StrongModuleTranspose> moduleIter = modules.iterator();
		Iterator<Integer> clusterIter = null;
		double v = 0.0;
		
		int i = 0;
		int j = 0;
		int gg = 0;
		
		computeFinalClusterTotals();
		HashMap<Integer,Integer> clusterCounts = null;
		
		i = 0;
		while(moduleIter.hasNext()) {
			module = moduleIter.next();
			clusterCounts = module.addUpClusterCounts();
			clusterIter = clusterCounts.keySet().iterator();
			while(clusterIter.hasNext()) {
				clusterID = clusterIter.next();
				v = ((double) clusterCounts.get(clusterID))/(((double) clusterTotals[clusterID])*((double) module.numOccur));
				clm[clusterID][i] = v;
			}
			i++;
		}
		
		int[] tgc = new int[clusterDPList.size()];
		double[][] signif = new double[clusterDPList.size()][modules.size()];
		int x = 0;
		int N = 0;
		int n = 0;
		int s = 0;
		
		HashSet<Integer> allGenes = new HashSet<Integer>();
		for (i=0;i<modules.size();i++) {
			module = modules.get(i);
			allGenes.addAll(module.allGenes());
		}
		N = allGenes.size();
		
		int docID = 0;
		for (i=0;i<clusterDPList.size();i++) {
			tgc[i] = 0;
			for (j=0;j<clusterDPList.get(i).children.size();j++) {
				docID = clusterDPList.get(i).children.get(j).code;
				if (allGenes.contains(docID)) {
					tgc[i]++;
				}
			}
		}
		
		for (i=0;i<modules.size();i++) {
			module = modules.get(i);
			n = module.allGenes().size();
			clusterIter = module.clusterTotal.keySet().iterator();
			while(clusterIter.hasNext()) {
				clusterID = clusterIter.next();
				x = module.clusterTotal.get(clusterID).size();
				s = tgc[clusterID];
				v = 1.0 - StatUtil.hyperGeometricCDF(x-1,N,s,n);
			//	if (v > 0.05/((double) modules.size())) {
				if (v > 0.05 | x < 5) {
					clm[clusterID][i] = 0.0;
				}
			}
		}
		
		// normalize the rows
		double norm = 0;
		for (i=0;i<clm.length;i++) {
			norm = 0;
			for (j=0;j<clm[0].length;j++) {
				if (!Double.isNaN(clm[i][j])) {
					norm += clm[i][j];
				}
			}
			if (norm > 0.0) {
				for (j=0;j<clm[0].length;j++) {
					clm[i][j] = clm[i][j]/norm;
				}
			}
		}
		
		return clm;
	}
	
	// merge together modules w/ overlapping clusters
	ArrayList<StrongModuleTranspose> collapseModules() {
		int i = 0;
		int j = 0;
		
		int nc = strongModuleMap.size();
		ArrayList<StrongModuleTranspose> modules = new ArrayList<StrongModuleTranspose>();
		Iterator<StrongModuleTranspose> iter = strongModuleMap.values().iterator();
		StrongModuleTranspose module = null;
		double v2 = 0.0;
		
		computeFinalClusterTotals();
		while(iter.hasNext()) {	
			module = (StrongModuleTranspose) iter.next();
			v2 = ((double) module.numOccur)/((double) numSamples);
			if (v2 >= minModuleOccurPercent) {
				modules.add(module);
			}
		}
		
		nc = modules.size();
		double[][] distances = new double[nc-1][];
		
		for (i=0;i<nc-1;i++) {
			distances[i] = new double[nc-i-1];
			for (j=(i+1);j<nc;j++) {
				distances[i][j-i-1] = 1.0 - modules.get(i).geneOverlap(modules.get(j));
			}
		}
		
		HierarchicalCluster hierCluster = new HierarchicalCluster();
		hierCluster.cluster(distances,1.0-minClusterOverlap);
		
		ArrayList<StrongModuleTranspose> modules2 = new ArrayList<StrongModuleTranspose>();
		StrongModuleTranspose module2 = null;
		
		int moduleNum = 0;
		for (i=0;i<hierCluster.clusters.size();i++) {
			moduleNum = hierCluster.clusters.get(i).get(0);
			module = (StrongModuleTranspose) modules.get(moduleNum);
			modules2.add(module);
			if (hierCluster.clusters.get(i).size() > 1) {
			/*	for (j=1;j<hierCluster.clusters.get(i).size();j++) {
					moduleNum = hierCluster.clusters.get(i).get(j);
					module2 = (StrongModuleTranspose) modules.get(moduleNum);
					module.controllingMotifs.retainAll(module2.controllingMotifs);
				} */
				for (j=1;j<hierCluster.clusters.get(i).size();j++) {
					moduleNum = hierCluster.clusters.get(i).get(j);
					module2 = (StrongModuleTranspose) modules.get(moduleNum);
					module.merge(module2);
				}
			}
		}
		
		ArrayList<StrongModuleTranspose> modules3 = new ArrayList<StrongModuleTranspose>();
		// perform a final filtering on merged modules
		for (i=0;i<modules2.size();i++) {
			module = (StrongModuleTranspose) modules2.get(i);
			v2 = ((double) module.numOccur)/((double) numSamples);
			if (v2 >= finalFilterModuleOccurPercent) {
				modules3.add(module);
			}
		}
		
		return modules3;
	}
}
