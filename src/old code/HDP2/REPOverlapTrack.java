package edu.mit.csail.psrg.georg.HDP2;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;

public class REPOverlapTrack {
	public int[][] rMatrix = null;
	int numOccur = 1;
	
	public REPOverlapTrack(int numGenes,int numTissues) {
		rMatrix = new int[numGenes][numTissues];
	}
	
	public void addTissue(HDP myHDP, TissueDP dp, LinkedHashMap<String,Integer> tissueMap,ArrayList<Integer> tissueGenes) {
		Iterator<Integer> iter = null;
		int geneID = 0;
		int tissueID = tissueMap.get(dp.label);
		
		iter = tissueGenes.iterator();
		while(iter.hasNext()) {
			geneID = iter.next();
			geneID = myHDP.reverseGeneMap[geneID];
			rMatrix[geneID][tissueID]++;
		}
		
	}
	
	public void add(REPOverlapTrack rep) {
		int i = 0;
		int j = 0;
		
		for (i=0;i<rMatrix.length;i++) {
			for (j=0;j<rMatrix[0].length;j++) {
				rMatrix[i][j] += rep.rMatrix[i][j];
			}
		}
		
		numOccur++;
	}
}
