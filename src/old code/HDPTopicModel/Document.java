package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.util.ArrayList;
import java.util.HashSet;
import java.io.Serializable;

public class Document implements Serializable {
	int[] genes;
	int[] geneClusterAssigns;
	int[] clusterNumData;
	String annotation;
	DirichletProcess parent = null;
	int docID = 0;
	int numUniqueGenes = 0;
	int code = 0;
	
	Document(String s, ArrayList<Integer> glist) {
		annotation = s;
		
		int i = 0;
		
		genes = new int[glist.size()];
		geneClusterAssigns = new int[glist.size()];
		
		for (i=0;i<glist.size();i++) {
			genes[i] = glist.get(i);
			geneClusterAssigns[i] = 0;
		}
		
		HashSet<Integer> uniqueGenes = new HashSet<Integer>();
		for (i=0;i<genes.length;i++) {
			uniqueGenes.add(genes[i]);
		}
		numUniqueGenes = uniqueGenes.size();
	}
	
	public int getCode() { return code; }
	public void setCode(int c) { code = c; }
	
	void deleteCluster(int c, int numClusters) {
		int i = 0;
		
		for (i=0;i<geneClusterAssigns.length;i++) {
			if (geneClusterAssigns[i] > c) {
				geneClusterAssigns[i] = geneClusterAssigns[i] - 1;
			}
		}
		System.arraycopy(clusterNumData,c+1,clusterNumData,c,numClusters-c);
	}
	
	void growClusters(int growClustersAlloc) {
  		int currSize = clusterNumData.length;
  		
  		int[] clusterNumData2 = new int[currSize + growClustersAlloc];
  		System.arraycopy(clusterNumData,0,clusterNumData2,0,currSize);
  		clusterNumData = clusterNumData2;
  	}
	
	void allocateGenesToClusterMap(HashSet[] clusterMap,ArrayList[] clusterArray,int numClusters) {
		int i = 0;
		int gene = 0;
		
		for (i=0;i<numClusters;i++) {
			clusterMap[i].clear();
			clusterArray[i].clear();
		}
		
		for (i=0;i<genes.length;i++) {
			gene = genes[i];
			clusterMap[geneClusterAssigns[i]].add(new Integer(gene));
			clusterArray[geneClusterAssigns[i]].add(new Integer(gene));
		}
	}
}
