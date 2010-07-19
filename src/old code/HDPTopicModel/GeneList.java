package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

public class GeneList {
	ArrayList<GeneHolder> genes = new ArrayList<GeneHolder>();
	
	void addGene(String[] names,double rank,double[] counts) {
		genes.add(new GeneHolder(names,rank,counts));
	}
	
	int size() {
		return genes.size();
	}
	
	double[] getRanks() {
		double[] ranks = new double[genes.size()];
		int i = 0;
		for (i=0;i<genes.size();i++) {
			ranks[i] = genes.get(i).rank;
		}
		return ranks;
	}
	
	void sort() {
		ArrayList<GeneHolder> genes2 = new ArrayList<GeneHolder>();
		GeneHolder[] garray = new GeneHolder[genes.size()];
		garray = genes.toArray(garray);
		List<GeneHolder> glist = Arrays.asList(garray);
		Collections.sort(glist);
		
		double[] ranks = new double[genes.size()];
		int i = 0;
		Iterator<GeneHolder> iter = glist.iterator();
		GeneHolder gene = null;
		while(iter.hasNext()) {
			gene = iter.next();
			ranks[i] = gene.rank;
			genes2.add(gene);
			i++;
		}
		
		genes = genes2;
	}
}
