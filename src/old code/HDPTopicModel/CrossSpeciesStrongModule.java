package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

public class CrossSpeciesStrongModule extends StrongModule {
	int numSpecies = 2;
	HashMap<Integer,Integer[]> geneMap = new HashMap<Integer,Integer[]>();
	int[] numOccurInSpecies = new int[numSpecies];
	
	public CrossSpeciesStrongModule(int numGenes, int[] posCounts) {
		super(numGenes, posCounts);
	}
	
	public void addControllingDocument(Document doc, HashSet documentOverlapGenes, ArrayList documentGeneRepeats) {
		double v = ((double) documentOverlapGenes.size())/((double) doc.numUniqueGenes);
		controllingDocuments.add(doc);
		documentLoad.put(doc,v);
		Object[] g2 = documentOverlapGenes.toArray();
		int i = 0;
		int j = 0;
		int idx = 0;
		Integer[] counts = null;
		int species = speciesCode(doc);
		
		for (i=0;i<g2.length;i++) {
			idx = ((Integer) g2[i]).intValue();
			counts = geneMap.get(idx);
			if (counts == null) {
				counts = new Integer[2*numSpecies];
				for (j=0;j<counts.length;j++) {
					counts[j] = new Integer(0);
				}
				geneMap.put(idx,counts);
			}
			counts[species*2]++;
		}
		
		g2 = documentGeneRepeats.toArray();
		for (i=0;i<g2.length;i++) {
			idx = ((Integer) g2[i]).intValue();
			counts = geneMap.get(idx);
			if (counts == null) {
				counts = new Integer[2*numSpecies];
				for (j=0;j<counts.length;j++) {
					counts[j] = new Integer(0);
				}
				geneMap.put(idx,counts);
			}
			counts[species*2+1]++;
		}
		numOccurInSpecies[species]++;
	}
	
	public int speciesCode(Document doc) {
		int f = 0;
		f = doc.annotation.indexOf("h_");
		if (f == 0)
			return 0;
		
		return 1;
	}
	
	public void merge(StrongModule module2) {
		int geneID = 0;
		Integer[] counts2 = null;
		Integer[] counts1 = null;
		Iterator<Integer> iter = geneMap.keySet().iterator();
		int i = 0;
		int j = 0;
		
		while(iter.hasNext()) {
			geneID = iter.next();
			counts2 = ((CrossSpeciesStrongModule) module2).geneMap.get(geneID);
			if (counts2 != null) {
				counts1 = geneMap.get(geneID);
				if (counts1 == null) {
					counts1 = new Integer[2*numSpecies];
					for (j=0;j<counts1.length;j++) {
						counts1[j] = new Integer(0);
					}
					geneMap.put(geneID,counts1);
				}
				for (i=0;i<counts2.length;i++) {
					counts1[i] += counts2[i];
				}
			}
		}
		
		Iterator<Document> iter2 = controllingDocuments.iterator();
		Document doc = null;
		double v = 0.0;
		
		while(iter.hasNext()) {
			doc = iter2.next();
			documentLoad.put(doc,new Double(documentLoad.get(doc) + module2.documentLoad.get(doc)));
		}
		
		numOccur += module2.numOccur;
		for (i=0;i<numSpecies;i++) {
			numOccurInSpecies[i] += ((CrossSpeciesStrongModule) module2).numOccurInSpecies[i];
		}
	}
}
