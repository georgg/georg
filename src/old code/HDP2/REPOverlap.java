package edu.mit.csail.psrg.georg.HDP2;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;

import edu.mit.csail.psrg.georg.StatUtil.StatUtil;

public class REPOverlap {
	// total counts for each gene
	public int[] geneCounts = null;
	// normalized by the number of significant tissues
	public double[] geneOccurNormalized = null;
	public double[] geneIntensityNormalized = null;
	// the multinomial distribution
	public double[] dist = null;
	public LinkedHashMap<String,OverlapTissue> tissues = new LinkedHashMap<String,OverlapTissue>();
	
	int numOccur = 1;
	int numUse = 0;
	int totalGenesInREP = 0;
	int totalCountsInREP = 0;
	
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
	public void addTissue(boolean UD,HashSet<Integer> tissueOverlapGenes,ArrayList<Integer> tissueGenes,TissueDP dp,int minGenesInREP,double tissueGenesThresh) {
		
		String tissueName = dp.label;
		double tiu = ((double) tissueGenes.size())/((double) dp.genes.length);
		double tou = ((double) tissueOverlapGenes.size())/((double) totalGenesInREP);
		boolean sig = false;
		
		double p = 0.0;
		
		int x = tissueOverlapGenes.size();
		
		if (x < minGenesInREP) {
			p = -1000.0;
		} else {
			p = tiu;
		}
				
		Iterator<Integer> iter = null;
		int geneID = 0;
		if (p >= tissueGenesThresh) {
			numUse++;
			sig = true;
			iter = tissueOverlapGenes.iterator();
			while(iter.hasNext()) {
				geneID = iter.next();
				geneOccurNormalized[geneID]++;
			}
			iter = tissueGenes.iterator();
			while(iter.hasNext()) {
				geneID = iter.next();
				geneIntensityNormalized[geneID]++;
			}
		}
		
		OverlapTissue ot = new OverlapTissue(tissueName,tiu,tou,UD,sig);
		tissues.put(tissueName,ot);
	}
	
	public void normCounts() {
		int i = 0;
		for (i=0;i<geneCounts.length;i++) {
			geneOccurNormalized[i] = geneOccurNormalized[i]/((double) numUse);
			geneIntensityNormalized[i] = geneIntensityNormalized[i]/((double) numUse);
		}
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
		updateDist();
		numOccur++;
	}
}
