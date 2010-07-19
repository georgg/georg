package edu.mit.csail.psrg.georg.HDP2;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

import edu.mit.csail.psrg.georg.StatUtil.StatUtil;

public class GenericREP extends REP {
	
	/**
	 * 
	 */
	private static final long serialVersionUID = -3179519025868030859L;
	// number of times gene occurs in at least one tissue (disregarding intensity)
	HashMap<Integer,Integer> geneCounts = null;
	HashMap<Integer,Integer> geneIntensity = null;
	
	public GenericREP(int numGenes,int[] posCounts) {
		super(numGenes,posCounts);
		geneIntensity = new HashMap<Integer,Integer>();
		geneCounts = new HashMap<Integer,Integer>();
	}
	
	public void addControllingTissue(TissueDP dp, HashSet<Integer> tissueOverlapGenes, ArrayList<Integer> tissueGeneRepeats, boolean UDUse) {
		double v = ((double) tissueOverlapGenes.size())/((double) dp.numUniqueGenes);
		ControllingTissue tissue = new ControllingTissue(dp,v,UDUse);
		controllingTissues.put(tissue,tissue);
		
		Object[] g2 = tissueOverlapGenes.toArray();
		int i = 0;
		int idx = 0;
		Integer vc = null;
		for (i=0;i<g2.length;i++) {
			idx = ((Integer) g2[i]).intValue();
			vc = geneCounts.get(idx);
			if (vc == null) {
				geneCounts.put(idx,1);
			} else {
				geneCounts.put(idx,vc+1);
			}
		}
		
		g2 = tissueGeneRepeats.toArray();
		for (i=0;i<g2.length;i++) {
			idx = ((Integer) g2[i]).intValue();
			vc = geneIntensity.get(idx);
			if (vc == null) {
				geneIntensity.put(idx,1);
			} else {
				geneIntensity.put(idx,vc+1);
			}
		}
	}
	
	public void merge(REP rep2) {
		Iterator<Integer> intIter1 = ((GenericREP) rep2).geneCounts.keySet().iterator();
		Iterator<Integer> intIter2 = ((GenericREP) rep2).geneCounts.values().iterator();
		int idx2 = 0;
		int gc2 = 0;
		Integer gc1 = 0;
		
		while(intIter1.hasNext()) {
			idx2 = intIter1.next();
			gc2 = intIter2.next();
			gc1 = geneCounts.get(idx2);
			if (gc1 == null) {
				geneCounts.put(idx2,gc2);
			} else {
				geneCounts.put(idx2,gc1+gc2);
			}
		}
		
		intIter1 = ((GenericREP) rep2).geneIntensity.keySet().iterator();
		intIter2 = ((GenericREP) rep2).geneIntensity.values().iterator();
		
		while(intIter1.hasNext()) {
			idx2 = intIter1.next();
			gc2 = intIter2.next();
			gc1 = geneIntensity.get(idx2);
			if (gc1 == null) {
				geneIntensity.put(idx2,gc2);
			} else {
				geneIntensity.put(idx2,gc1+gc2);
			}
		}
		
		
		Iterator<ControllingTissue> iter = rep2.controllingTissues.keySet().iterator();
		ControllingTissue tissue2 = null;
		ControllingTissue tissue = null;
		
		while(iter.hasNext()) {
			tissue2 = iter.next();
			tissue = controllingTissues.get(tissue2);
			if (tissue != null) {
				tissue.merge(tissue2);
			}
		}
		
		numOccur = numOccur + rep2.numOccur;
	}
}
