package edu.mit.csail.psrg.georg.GeneProgram2;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;

import edu.mit.csail.psrg.georg.StatUtil.StatUtil;

public class REPOverlap implements Comparable {
	// total counts for each gene
	public int[] geneCounts = null;
	// normalized by the number of significant tissues
	public double[] geneOccurNormalized = null;
	public double[] geneIntensityNormalized = null;
	// the multinomial distribution
	public double[] dist = null;
	public LinkedHashMap<String,OverlapTissue> tissues = new LinkedHashMap<String,OverlapTissue>();
	public LinkedHashMap<ParentDP,ArrayList<OverlapTissue>> groupMap = new LinkedHashMap<ParentDP,ArrayList<OverlapTissue>>();
	
	int numOccur = 1;
	int numUse = 0;
	double numUseWeight = 0.0;
	int totalGenesInREP = 0;
	int totalCountsInREP = 0;
	
	double generality = 0.0;
	
	public REPOverlap(int numGenes) {
		geneCounts = new int[numGenes];
		geneOccurNormalized = new double[numGenes];
		geneIntensityNormalized = new double[numGenes];
		dist = new double[numGenes];
	}
	
	public void addGeneCounts(int[] counts) {
		System.arraycopy(counts,0,geneCounts,0,counts.length);
		updateDist();
	}
	
	// call updateDist before calling this method
	public void addTissue(HDP myHDP,boolean UD,int[] modifiers,HashSet<Integer> tissueOverlapGenes,ArrayList<Integer> tissueGenes,TissueDP dp,int minGenesInREP,double tissueGenesThresh) {
		
		String tissueName = dp.label;
		double tiu = ((double) tissueGenes.size())/((double) dp.genes.length);
	//	double tou = ((double) tissueOverlapGenes.size())/((double) totalGenesInREP);
		double tou = 0.0;
		boolean sig = false;
		
		double p = 0.0;
		
		Iterator<Integer> iter = null;
		int geneID = 0;
		
		if (tissueGenes.size() > 0) {
			iter = tissueGenes.iterator();
			while(iter.hasNext()) {
				geneID = iter.next();
				tou += dist[geneID];
			}
		}
		
		int x = tissueOverlapGenes.size();
		
		if (x < minGenesInREP) {
			p = -1000.0;
		} else {
		//	p = tiu;
			p = tou;
		}
		
		if (p >= tissueGenesThresh) {
			numUse++;
			numUseWeight += tou;
			sig = true;
			iter = tissueOverlapGenes.iterator();
			while(iter.hasNext()) {
				geneID = iter.next();
			//	geneOccurNormalized[geneID]++;
				geneOccurNormalized[geneID] += tou;
			}
			iter = tissueGenes.iterator();
			while(iter.hasNext()) {
				geneID = iter.next();
			//	geneIntensityNormalized[geneID]++;
				geneIntensityNormalized[geneID] += tou;
			}
		}
		
		OverlapTissue ot = new OverlapTissue(tissueName,tiu,tou,UD,sig,modifiers,myHDP);
		tissues.put(tissueName,ot);
		
		if (myHDP.DPGroups.size() > 0) {
			ArrayList<OverlapTissue> tArray = groupMap.get(dp.parent);
			if (tArray == null) {
				tArray = new ArrayList<OverlapTissue>();
				groupMap.put(dp.parent,tArray);
			}
			tArray.add(ot);
		}
		
		if (p >= tissueGenesThresh) {
			ot.passThresh = true;
		} else {
			ot.passThresh = false;
		}
	}
	
	public void normCounts() {
		int i = 0;
		for (i=0;i<geneCounts.length;i++) {
		//	geneOccurNormalized[i] = geneOccurNormalized[i]/((double) numUse);
			geneOccurNormalized[i] = geneOccurNormalized[i]/numUseWeight;
		//	geneIntensityNormalized[i] = geneIntensityNormalized[i]/((double) numUse);
			geneIntensityNormalized[i] = geneIntensityNormalized[i]/numUseWeight;
		}
		
		computeGenerality();
	}
	
	public void updateDist() {
		int i = 0;
		totalGenesInREP = 0;
		totalCountsInREP = 0;
		for (i=0;i<geneCounts.length;i++) {
			totalCountsInREP += geneCounts[i];
			if (geneCounts[i] > 0) {
				totalGenesInREP++;
			}
		}
		for (i=0;i<geneCounts.length;i++) {
			dist[i] = geneCounts[i]/((double) totalCountsInREP);
		}
	}
	

	// calculates the Hellinger similarity between the two probability
	// distributions; other distances can be used
	public double similarity(REPOverlap rep2) {
		double p = 0;
		int i = 0;
		for (i=0;i<dist.length;i++) {
			p += Math.sqrt(dist[i]*rep2.dist[i]);
		}
		return p;
	}
	
	public void add(REPOverlap rep2) {
		int i = 0;
		Iterator<String> iter = tissues.keySet().iterator();
		String tissueName = null;
		OverlapTissue tissue1 = null;
		OverlapTissue tissue2 = null;
		
		while(iter.hasNext()) {
			tissueName = iter.next();
			tissue1 = tissues.get(tissueName);
			tissue2 = rep2.tissues.get(tissueName);
			tissue1.add(tissue2);
		}
		
		for (i=0;i<geneCounts.length;i++) {
			geneCounts[i] += rep2.geneCounts[i];
			geneOccurNormalized[i] += rep2.geneOccurNormalized[i];
			geneIntensityNormalized[i] += rep2.geneIntensityNormalized[i];
		}
		
		generality += rep2.generality;
		
		updateDist();
		numOccur++;
	}
	
	public void computeGenerality() {
		if (groupMap.keySet().isEmpty()) {
			return;
		}
		
		double[] vals = new double[groupMap.keySet().size()];
		double norm = 0.0;
		
		ArrayList<OverlapTissue> tArray = null;
		Iterator<ArrayList<OverlapTissue>> iter = groupMap.values().iterator();
		int i = 0;
		
		int grpNum = 0;
		while(iter.hasNext()) {
			tArray = iter.next();
			for (i=0;i<tArray.size();i++) {
				if (tArray.get(i).passThresh) {
					norm += tArray.get(i).topicUse;
					vals[grpNum] += tArray.get(i).topicUse;
				}
			}
			grpNum++;
		}
		
		generality = 0.0;
		
		if (vals.length > 1) {
			for (i=0;i<vals.length;i++) {
				vals[i] = vals[i]/norm;
				if (vals[i] > 0.0) {
					generality += -(vals[i]*Math.log(vals[i]));
				}
			}
		}
		
		if (vals.length > 0) {
			generality = generality/Math.log((double) vals.length);
		}
		
		if ((new Double(generality)).isNaN()) {
			generality = 0.0;
		}
		
		groupMap.clear();
	}

	public int compareTo(Object o) {
		REPOverlap r2 = (REPOverlap) o;
		
		if (generality > r2.generality)
			return 1;
		
		if (generality < r2.generality)
			return -1;
		
		return 0;
	}
}
