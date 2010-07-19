package edu.mit.csail.psrg.georg.BindDP;

import edu.mit.csail.psrg.georg.SampleUtil.LogDensityFn;
import cern.jet.stat.Gamma;

public class BindCountPriorPosteriorLogDensityFn implements LogDensityFn {
	int[] numGenesInCluster = null;
	int[][] numBound = null;
	double[][] lambda = null;
	int numClusters = 0;
	int numExperiments = 0;
	
	public void setParams(int[] ngc, int[][] nb, double[][] lb, int nc) {
		numGenesInCluster = ngc;
		lambda = lb;
		numBound = nb;
		numClusters = nc;
		numExperiments = lambda[0].length;
	}
	
	public double eval(double[] x) {
		double alpha = x[0];
		double beta = x[1];
		double h = 0.0;
		double n1 = 0.0;
		double n2 = 0.0;
		
		if (alpha <= 0.0 || beta <= 0.0) {
			return -1.0/0.0;
		}
		
		h = -(5.0/2.0)*Math.log(alpha+beta);
		h = h + ((double) numExperiments*numClusters)*(Gamma.logGamma(alpha+beta) - Gamma.logGamma(alpha) - Gamma.logGamma(beta));
		
		int cc = 0;
		int j = 0;
		for (cc=0;cc<numClusters;cc++) {
			for (j=0;j<numExperiments;j++) {
				n1 = (alpha-1.0 + (double) numBound[cc][j]);
				h = h + n1*Math.log(lambda[cc][j]);
				n2 = (beta-1.0 + (double) (numGenesInCluster[cc] - numBound[cc][j]));
				h = h + n2*Math.log(1-lambda[cc][j]);
			//	h = h + (Gamma.logGamma(n1+n2) - Gamma.logGamma(n1) - Gamma.logGamma(n2));
			}
		} 
		
		return h;
	}

}
