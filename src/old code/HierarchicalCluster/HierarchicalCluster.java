package edu.mit.csail.psrg.georg.HierarchicalCluster;

import java.util.ArrayList;

public class HierarchicalCluster {
	public ArrayList< ArrayList<Integer> > clusters;
	
	public void cluster(double[][] distances, double cutDist) {
		int i = 0;
		int j = 0;
		clusters = new ArrayList< ArrayList<Integer> >();
		ArrayList<Integer> arl = null;
		
		// assign each data item to its own cluster
		for (i=0;i<distances.length+1;i++) {
			arl = new ArrayList<Integer>();
			arl.add(i);
			clusters.add(arl);
		}
		
		int merge1 = 0;
		int merge2 = 0;
		int nc = 0;
		double d = 0.0;
		double minD = 100000;
		
		System.out.println("Hierarchical clustering...");
		
		while(merge1 >= 0 & clusters.size() > 1) {
			merge1 = -1;
			merge2 = -1;
			minD = 10000;
			nc = clusters.size();
			System.out.println("Num clusters=" + nc);
			for (i=0;i<nc-1;i++) {
				for (j=(i+1);j<nc;j++) {
					d = totalLinkage(distances,i,j);
					if (d < minD & d < cutDist) {
						merge1 = i;
						merge2 = j;
						minD = d;
					}
				}
			}
			if (merge1 >= 0) {
				merge(merge1,merge2);
			}
		}
	}
	
//	 compute total linkage for the two clusters (max distance between any pair)
	double totalLinkage(double[][] distances, int cluster1, int cluster2) {
		double maxD = 0.0;
		double d = 0.0;
		int i = 0;
		int j = 0;
		int ID1 = 0;
		int ID2 = 0;
		int n1 = clusters.get(cluster1).size();
		int n2 = clusters.get(cluster2).size();
		
		for (i=0;i<n1;i++) {
			ID1 = clusters.get(cluster1).get(i);
			for (j=0;j<n2;j++) {
				ID2 = clusters.get(cluster2).get(j);
				if (ID1 > ID2) {
					d = distances[ID2][ID1-ID2-1];
				} else {
					d = distances[ID1][ID2-ID1-1];
				}
				if (d > maxD) {
					maxD = d;
				}
			}
		}
		return maxD;
	}
	
	void merge(int merge1,int merge2) {
		int gg = 0;
		int ng = 0;
		
		ng = clusters.get(merge2).size();
		
		for (gg=0;gg<ng;gg++) {
			clusters.get(merge1).add(clusters.get(merge2).get(gg));
		}
		
		clusters.remove(merge2);
	}
	
	
}
