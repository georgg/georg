package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.util.ArrayList;

public class DPDistanceCollector extends HDPStatCollector {
	DirichletProcess[] DPs = null;
	double[][] pairProb = null;
	int numSamples = 0;

	public DPDistanceCollector(HierDirichletProcess h) {
		super(h);
	}
	
	public DPDistanceCollector(HierDirichletProcess h,ArrayList<DirichletProcess> d) {
		super(h);
		int i = 0;
		
		DPs = new DirichletProcess[d.size()];
		DPs = d.toArray(DPs);
		
		pairProb = new double[DPs.length-1][];
		
		for (i=0;i<DPs.length-1;i++) {
			pairProb[i] = new double[DPs.length-i-1];
		}
	}

	public void updateStats() {
		int i = 0;
		int j = 0;
		for (i=0;i<DPs.length-1;i++) {
			for (j=i+1;j<DPs.length;j++) {
				pairProb[i][j-i-1] += hellingerDistance(DPs[i],DPs[j]);
			}
		}
		numSamples++;
	}
	
	public double hellingerDistance(DirichletProcess DP1,DirichletProcess DP2) {
		int cc = 0;
		int numTopics = myHDP.numClusters;
		int norm1 = 0;
		int norm2 = 0;
		double dist = 0.0;
		double norm = 0.0;
		
		dist = 0.0;
		for (cc=0;cc<numTopics;cc++) {
			norm1 += DP1.clusterNumData[cc];
			norm2 += DP2.clusterNumData[cc];
			// compute the Hellinger distance between the distributions
			dist += Math.sqrt((double) (DP1.clusterNumData[cc]*DP2.clusterNumData[cc]));
		}
		norm = Math.sqrt((double) (norm1*norm2));
		dist = dist/norm;
		
		return dist;
	}

}
