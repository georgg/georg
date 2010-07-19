package edu.mit.edu.csail.psrg.georg.AgglomerativeCluster;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class AgglomerativeCluster {
	public int[][] data = null;
	public int maxValue = 0;
	public ArrayList<ACluster> clusters = new ArrayList<ACluster>();
	public int numCols = 0;
	String[] labels = null;
	
	public AgglomerativeCluster(int[][] d, int n, String[] l) {
		data = d;
		numCols = n;
		labels = l;
		init();
	}
	
	void init() {
		int i = 0;
		int j = 0;
		
		for (i=0;i<data.length;i++) {
			clusters.add(new ACluster(this,i));
			for (j=0;j<data[i].length;j++) {
				if (data[i][j] > maxValue) {
					maxValue = data[i][j];
				}
			}
		}
		maxValue++;
	}
	
	public void iterate() {
		int i = 0;
		int j = 0;
		ACluster cluster = null;
		double oldLike = 0.0;
		double newLike = 0.0;
		double sampleSize = (double) data.length;
	//	double numParams = (double) (clusters.size()*(maxValue-1)*numCols);
		double numParams = (double) (clusters.size()*(numCols-1));
		
		for (i=0;i<clusters.size();i++) {
			cluster = clusters.get(i);
			cluster.updateLikelihood();
			newLike += cluster.lik;
		}
		
		double newPenalty = -numParams*Math.log(sampleSize)/2.0;
		double oldPenalty = 0.0;
		
		boolean continueMerge = true;
		int merge1 = 0;
		int merge2 = 0;
		double bestLike = 0.0;
		double l = 0.0;
		ACluster cluster1 = null;
		ACluster cluster2 = null;
		ACluster mergeCluster = null;
		ArrayList<ACluster> newClusters = new ArrayList<ACluster>();
		
		while(continueMerge & clusters.size() > 1) {
			oldLike = newLike;
			oldPenalty = newPenalty;
			
		//	numParams = (double) ((clusters.size()-1)*(maxValue-1)*numCols);
			numParams = (double) ((clusters.size()-1)*(numCols-1));
			newPenalty = -numParams*Math.log(sampleSize)/2.0;
			bestLike = 0.0/0.0;
			merge1 = -1;
			merge2 = -1;
			for (i=0;i<clusters.size()-1;i++) {
				for (j=(i+1);j<clusters.size();j++) {
					cluster1 = clusters.get(i);
					cluster2 = clusters.get(j);
					mergeCluster = ACluster.merge(cluster1,cluster2);
					l = mergeCluster.lik - cluster1.lik - cluster2.lik;
					if (l > bestLike | Double.isNaN(bestLike)) {
						bestLike = l;
						merge1 = i;
						merge2 = j;
					}
				}
			}
			
			newLike = oldLike + bestLike;
			
			if (newLike + newPenalty > oldLike + oldPenalty) {
				newClusters.clear();
				System.out.println("Merge " + merge1 + " " + merge2);
				System.out.println(newLike + newPenalty);
				cluster1 = clusters.get(merge1);
				cluster2 = clusters.get(merge2);
				System.out.println(getClusterName(cluster1) + " <<>> " + getClusterName(cluster2));
				mergeCluster = ACluster.merge(cluster1,cluster2);
				newClusters.add(mergeCluster);
				for (i=0;i<clusters.size();i++) {
					if ((i == merge1) | (i == merge2)) {
					} else {
						newClusters.add(clusters.get(i));
					}
				}
				clusters.clear();
				clusters.addAll(newClusters);
			} else {
				continueMerge = false;
			}
			System.out.println(clusters.size());
		}
	}
	
	public String getClusterName(ACluster cluster) {
		int idx = 0;
		int i = 0;
		String s = "";
		for (idx=0;idx<cluster.dataItems.size();idx++) {
			i = cluster.dataItems.get(idx);
			s = s + labels[i] + " | ";
		}
		return s;
	}
	
	public void outputClusters(String fname) throws IOException {
		FileWriter file = new FileWriter(fname);
		int i = 0;
		int j = 0;
		int k = 0;
		ACluster cluster = null;
		
		for (i=0;i<clusters.size();i++) {
			cluster = clusters.get(i);
			for (j=0;j<cluster.dataItems.size();j++) {
				k = cluster.dataItems.get(j);
				file.write(labels[k]);
				if (j < cluster.dataItems.size() - 1) {
					file.write("\t");
				}
			}
			file.write("\n");
		}
		file.close();
	}
}
