package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.Serializable;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

public class StrongModuleTranspose implements Serializable, Comparable {
	// motifs that "control" the module are those w/ a significant overlap with
	// the represented genes in the module
	HashSet<Integer> controllingMotifs = new HashSet<Integer>();
	HashMap<Integer,Double> motifLoad = new HashMap<Integer,Double>();
	// gene counts indexed by cluster
	HashMap<Integer,HashMap<Integer,Integer>> clusterTotal = new HashMap<Integer,HashMap<Integer,Integer>>();
	
	// total number of genes in the module
	int numOccur = 1;
	
	// ml = motif load
	public void addControllingMotif(int motifIDX,HashMap<Integer,HashMap<Integer,Integer>> clusterCounts) {
		boolean present = false;
	
		Iterator<Integer> clusterIter = clusterCounts.keySet().iterator();
		Iterator<Integer> geneIter = null;
		HashMap<Integer,Integer> geneToCountsAdd = null;
		HashMap<Integer,Integer> geneToCounts = null;
		int clusterID = 0;
		int geneID = 0;
		Integer count = null;
		int countAdd = 0;
		int motifEvents = 0;
		
		while (clusterIter.hasNext()) {
			clusterID = clusterIter.next();
			geneToCountsAdd = clusterCounts.get(clusterID);
			geneToCounts = clusterTotal.get(clusterID);
			if (geneToCounts == null) {
				geneToCounts = new HashMap<Integer,Integer>();
				clusterTotal.put(clusterID,geneToCounts);
			}
			geneIter = geneToCountsAdd.keySet().iterator();
			while(geneIter.hasNext()) {
				geneID = geneIter.next();
				countAdd = geneToCountsAdd.get(geneID);
				count = geneToCounts.get(geneID);
				if (count == null) {
					count = new Integer(0);
				}
				geneToCounts.put(geneID,count+countAdd);
				motifEvents += countAdd;
			}
		}
			
		double ml = ((double) motifEvents);
		present = controllingMotifs.contains(motifIDX);
		if (!present) {
			controllingMotifs.add(motifIDX);
			motifLoad.put(motifIDX,ml);
		} else {
			motifLoad.put(motifIDX,ml+motifLoad.get(motifIDX));
		}
	}
	
	public void norm() {
		Iterator<Integer> motifIter = motifLoad.keySet().iterator();
		double norm = 0.0;
		int motifID = 0;
		while(motifIter.hasNext()) {
			motifID = motifIter.next();
			norm += motifLoad.get(motifID);
		}
		motifIter = motifLoad.keySet().iterator();
		HashMap<Integer,Double> motifLoad2 = new HashMap<Integer,Double>();
		while (motifIter.hasNext()) {
			motifID = motifIter.next();
			motifLoad2.put(motifID,motifLoad.get(motifID)/norm);
		}
		motifLoad = motifLoad2;
	}
	
	public void merge(StrongModuleTranspose m2) {
		numOccur = numOccur + m2.numOccur;
		
		Iterator<Integer> clusterIter = m2.clusterTotal.keySet().iterator();
		Iterator<Integer> geneIter = null;
		HashMap<Integer,Integer> geneToCountsAdd = null;
		HashMap<Integer,Integer> geneToCounts = null;
		int countAdd = 0;
		Integer count = null;
		int clusterID = 0;
		int geneID = 0;
		int motifID = 0;
		
		while(clusterIter.hasNext()) {
			clusterID = clusterIter.next();
			geneToCountsAdd = m2.clusterTotal.get(clusterID);
			geneToCounts = clusterTotal.get(clusterID);
			if (geneToCounts == null) {
				geneToCounts = new HashMap<Integer,Integer>();
				clusterTotal.put(clusterID,geneToCounts);
			}
			
			geneIter = geneToCountsAdd.keySet().iterator();
			while(geneIter.hasNext()) {
				geneID = geneIter.next();
				countAdd = geneToCountsAdd.get(geneID);
				count = geneToCounts.get(geneID);
				if (count == null) {
					count = new Integer(0);
				}
				geneToCounts.put(geneID,count+countAdd);
			}
		}
		
	//	Iterator<Integer> motifIter = motifLoad.keySet().iterator();
		Iterator<Integer> motifIter = m2.motifLoad.keySet().iterator();
		
		while(motifIter.hasNext()) {
			motifID = motifIter.next();
			if (!motifLoad.containsKey(motifID)) {
				motifLoad.put(motifID,0.0);
			}
			motifLoad.put(motifID,motifLoad.get(motifID) + m2.motifLoad.get(motifID));
		/*	if (m2.motifLoad.containsKey(motifID) & motifLoad.containsKey(motifID)) {
				motifLoad.put(motifID,motifLoad.get(motifID) + m2.motifLoad.get(motifID));
			} */
		}
	}
	
	public int compareTo(Object o) {
		StrongModuleTranspose module2 = (StrongModuleTranspose) o;
		
		if (numOccur == module2.numOccur) {
			return 0;
		}
		
		if (numOccur > module2.numOccur)
			return 1;
		
		return -1;
	}
	
	public void addToModuleMap(HashMap<HashSet<Integer>,StrongModuleTranspose> strongModuleMap) {
		boolean alreadyThere = false;
		alreadyThere = strongModuleMap.containsKey(controllingMotifs);
		if (!alreadyThere) {
			strongModuleMap.put(controllingMotifs,this);
		} else {
			StrongModuleTranspose module = strongModuleMap.get(controllingMotifs);
			module.merge(this);
		}
	}
	
	public HashSet<Integer> allGenes() {
		HashSet<Integer> genes = new HashSet<Integer>();
		HashMap<Integer,Integer> genesToCountsMap = null;
		Iterator<HashMap<Integer,Integer>> clusterMapIter = clusterTotal.values().iterator();
		while(clusterMapIter.hasNext()) {
			genesToCountsMap = clusterMapIter.next();
			genes.addAll(genesToCountsMap.keySet());
		}
		return genes;
	}
	
	public double geneOverlap(StrongModuleTranspose m2) {
		int unionTotal = 0;
		int intersectTotal = 0;
		
		HashSet<Integer> g1 = allGenes();
		HashSet<Integer> g2 = m2.allGenes();
		
		HashSet<Integer> g3 = new HashSet<Integer>();
		g3.addAll(g1);
		g3.addAll(g2);
		unionTotal = g3.size();
		
		g1.retainAll(g2);
		intersectTotal = g1.size();
		
		return ((double) intersectTotal)/((double) unionTotal);
	}
	
	public HashMap<Integer,Integer> addUpClusterCounts() {
		HashMap<Integer,Integer> map = new HashMap<Integer,Integer>();
		int clusterID = 0;
		int clusterCount = 0;
		HashMap<Integer,Integer> geneToCounts = null;
		Iterator<Integer> clusterIter = clusterTotal.keySet().iterator();
		Iterator<Integer> countIter = null;
		
		while(clusterIter.hasNext()) {
			clusterCount = 0;
			clusterID = clusterIter.next();
			geneToCounts = clusterTotal.get(clusterID);
			countIter = geneToCounts.values().iterator();
			while (countIter.hasNext()) {
				clusterCount += countIter.next();
			}
			map.put(clusterID,clusterCount);
		}
		
		return map;
	}
}
