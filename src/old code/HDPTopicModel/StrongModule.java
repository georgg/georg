package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import edu.mit.csail.psrg.georg.StatUtil.StatUtil;

public abstract class StrongModule implements Serializable,Comparable {
	// documents that "control" the module are those w/ a significant overlap with
	// the represented genes in the module
	HashSet<Document> controllingDocuments = new HashSet<Document>();
	HashMap<Document,Double> documentLoad = new HashMap<Document,Double>();
	
	// number of times the strong module occurred
	int numOccur = 1;
	
	// total number of genes in the module
	int totalGenesInModule = 0;
	
	public StrongModule(int numGenes,int[] posCounts) {
		int i = 0;
		for (i=0;i<posCounts.length;i++) {
			if (posCounts[i] > 0) {
				totalGenesInModule++;
			}
		}
	}
	
	public abstract void addControllingDocument(Document doc, HashSet documentOverlapGenes, ArrayList documentGeneRepeats);
	
	public int compareTo(Object o) {
		StrongModule module2 = (StrongModule) o;
		
		if (numOccur == module2.numOccur) {
			return 0;
		}
		
		if (numOccur > module2.numOccur)
			return 1;
		
		return -1;
	}
	
	public double controlPVal(HashSet documentOverlapGenes,Document doc,int totalGenes,int minGenesInModule) {
		double p = 0.0;
		
		int x = documentOverlapGenes.size();
		
		if (x < minGenesInModule) {
			return 1.0;
		}
		
		int s = totalGenesInModule;
		
		int N = totalGenes;
		int n = doc.numUniqueGenes;
		
//		p = 1.0-StatUtil.hyperGeometricCDF(x-1,N,s,n);
		p = 1.0-StatUtil.hyperGeometricCDF(x-1,N,n,s);
		
		return p;
	}
	
	public abstract void merge(StrongModule module2);
	
	public void addToModuleMap(HashMap<HashSet<Document>,StrongModule> strongModuleMap, int totalNumGenes) {
		boolean alreadyThere = false;
		alreadyThere = strongModuleMap.containsKey(controllingDocuments);
		if (!alreadyThere) {
			strongModuleMap.put(controllingDocuments,this);
		} else {
			StrongModule module = strongModuleMap.get(controllingDocuments);
			module.merge(this);
		}
	}
}
