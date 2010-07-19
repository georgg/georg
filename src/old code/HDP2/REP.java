package edu.mit.csail.psrg.georg.HDP2;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

import edu.mit.csail.psrg.georg.StatUtil.StatUtil;

public abstract class REP implements Serializable,Comparable {
	private static final long serialVersionUID = 2282033931753592329L;
	// tissues that "control" the REP are those w/ a significant overlap with
	// the represented genes in the REP
	public HashMap<ControllingTissue,ControllingTissue> controllingTissues = new HashMap<ControllingTissue,ControllingTissue>();
	
	// number of times the REP occurred
	public int numOccur = 1;
	
	// total number of genes in the REP
	public int totalGenesInREP = 0;
	
	public REP(int numGenes,int[] posCounts) {
		int i = 0;
		for (i=0;i<posCounts.length;i++) {
			if (posCounts[i] > 0) {
				totalGenesInREP++;
			}
		}
	}
	
	public abstract void addControllingTissue(TissueDP dp, HashSet<Integer> tissueOverlapGenes, ArrayList<Integer> tissueGeneRepeats, boolean UDUse);
	
	public int compareTo(Object o) {
		REP module2 = (REP) o;
		
		if (numOccur == module2.numOccur) {
			return 0;
		}
		
		if (numOccur > module2.numOccur)
			return 1;
		
		return -1;
	}
	
	public double controlPVal(HashSet<Integer> tissueOverlapGenes,TissueDP dp,int totalGenes,int minGenesInREP) {
		double p = 0.0;
		
		int x = tissueOverlapGenes.size();
		
		if (x < minGenesInREP) {
			return 1.0;
		}
		
		int s = totalGenesInREP;
		
		int N = totalGenes;
		int n = dp.numUniqueGenes;
		
//		p = 1.0-StatUtil.hyperGeometricCDF(x-1,N,s,n);
		p = 1.0-StatUtil.hyperGeometricCDF(x-1,N,n,s);
		
		return p;
	}
	
	public abstract void merge(REP module2);
	
	public void addToREPMap(HashMap<HashMap<ControllingTissue,ControllingTissue>,REP> REPMap, int totalNumGenes) {
		boolean alreadyThere = false;
		alreadyThere = REPMap.containsKey(controllingTissues);
		if (!alreadyThere) {
			REPMap.put(controllingTissues,this);
		} else {
			REP rep = REPMap.get(controllingTissues);
			rep.merge(this);
		}
	}
	
	public Vector<ControllingTissue> sortTissues() {
		Vector<ControllingTissue> cts = new Vector<ControllingTissue>();
		Iterator<ControllingTissue> iter = controllingTissues.keySet().iterator();
		while(iter.hasNext()) {
			cts.add(iter.next());
		}
		Collections.sort(cts);
		return cts;
	}
	
	public void intersectControllingTissues(REP rep2) {
		HashMap<ControllingTissue,ControllingTissue> cts = new HashMap<ControllingTissue,ControllingTissue>();
		Iterator<ControllingTissue> iter = rep2.controllingTissues.keySet().iterator();
		ControllingTissue ct = null;
		while(iter.hasNext()) {
			ct = iter.next();
			ct = controllingTissues.get(ct);
			if (ct != null) {
				cts.put(ct,ct);
			}
		}
		controllingTissues = cts;
	}
}
