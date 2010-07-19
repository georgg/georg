package edu.mit.csail.psrg.georg.HDP2;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;

public class Topic implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = -2156340130982474428L;
	// how many positive counts (for each gene) are associated with each topic
	public int[] posCounts = null;
	public int totalCount = 0;
	// temporary map used for evaluating tissue loading on topics
	public HashSet<Integer> topicMap = new HashSet<Integer>();
	// temporary array for evaluating tissue loading on topics
	public ArrayList<Integer> topicArray = new ArrayList<Integer>();
	
	public Topic(int numGenes) {
		posCounts = new int[numGenes];
	}
	
	public void zeroCounts() {
		int i = 0;
		for (i=0;i<posCounts.length;i++) {
			posCounts[i] = 0;
		}
		totalCount = 0;
	}
	
	public void removeGeneCount(int gene) {
		posCounts[gene]--;
		totalCount--;
	}
	
	public void addGeneCount(int gene) {
		posCounts[gene]++;
		totalCount++;
	}
}
