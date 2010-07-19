package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.TreeSet;

import edu.mit.csail.psrg.georg.StatUtil.StatUtil;

public class SingleStrongModule extends StrongModule {
	// documents that "control" the module are those w/ a significant overlap with
	// the represented genes in the module
	HashSet<Document> controllingDocuments = new HashSet<Document>(10);
	HashMap<Document,Double> documentLoad = new HashMap<Document,Double>(10);
	
	// number of times gene occurs in at least one document (disregarding intensity)
	int[] geneCounts = null;
	int[] geneIntensity = null;
	
	// total number of genes in the module
	int totalGenesInModule = 0;
	int numOccur = 1;
	
	public SingleStrongModule(int numGenes,int[] posCounts) {
		super(numGenes,posCounts);
		geneIntensity = new int[numGenes];
		geneCounts = new int[numGenes];
	}
	
	public void addControllingDocument(Document doc, HashSet documentOverlapGenes, ArrayList documentGeneRepeats) {
		double v = ((double) documentOverlapGenes.size())/((double) doc.numUniqueGenes);
		controllingDocuments.add(doc);
		documentLoad.put(doc,v);
		Object[] g2 = documentOverlapGenes.toArray();
		int i = 0;
		int idx = 0;
		for (i=0;i<g2.length;i++) {
			idx = ((Integer) g2[i]).intValue();
			geneCounts[idx]++;
		}
		
		g2 = documentGeneRepeats.toArray();
		for (i=0;i<g2.length;i++) {
			idx = ((Integer) g2[i]).intValue();
			geneIntensity[idx]++;
		}
	}
	
	public void merge(StrongModule module2) {
		int i = 0;
		for (i=0;i<geneIntensity.length;i++) {
			geneIntensity[i] += ((SingleStrongModule) module2).geneIntensity[i];
			geneCounts[i] += ((SingleStrongModule) module2).geneCounts[i];
		}
		
		Iterator<Document> iter = controllingDocuments.iterator();
		Document doc = null;
		double v = 0.0;
		
		while(iter.hasNext()) {
			doc = iter.next();
			v = documentLoad.get(doc) + module2.documentLoad.get(doc);
			documentLoad.put(doc,v);
		}
		
		numOccur += module2.numOccur;
	}
}
