package com.newmacondo.georg.DynamicDUO;

import java.io.Serializable;

import com.newmacondo.georg.StatUtil.StatUtil;

// each cluster represents the core dynamics, which can then be mapped
// to states via permutation operators or time-shift operators
public class DynamicsCluster implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = -4197338717134677716L;
	public DP myDP = null;
	public int[][] transitionCounts;
	public int[] initCounts;
	public double[][] transitionProbs;
	public double[] initProbs;
	public double[] alpha;
	int numDataAssigned = 0;
	// total data in each state
	int[] totalState;
	
	public DynamicsCluster(DP dp) {
		myDP = dp;
		transitionCounts = new int[myDP.numStates][myDP.numStates];
		initCounts = new int[myDP.numStates];
		initProbs = new double[myDP.numStates];
		transitionProbs = new double[myDP.numStates][myDP.numStates];
		totalState = new int[myDP.numStates];
		alpha = new double[myDP.numStates];
	}
	
	public void sampleTransitionProbsPosterior() {
		int i = 0;
		int j = 0;
		
		for (i=0;i<myDP.numStates;i++) {
			alpha[i] = myDP.transitionPriorAlpha + (double) initCounts[i];
		}
		StatUtil.dirichlet_rnd(initProbs, alpha, myDP.numStates);
		
		for (i=0;i<myDP.numStates;i++) {
			for (j=0;j<myDP.numStates;j++) {
				alpha[j] = myDP.transitionPriorAlpha + (double) transitionCounts[i][j];
			}
			StatUtil.dirichlet_rnd(transitionProbs[i], alpha, myDP.numStates);
		}
	}
	
	public void addOTU(int OTUNum,int locationNum) {
		numDataAssigned++;
		int k = 0;
		int state1 = myDP.states[OTUNum][locationNum][0];
		int state2 = 0;
		initCounts[state1]++;
		for (k=1;k<myDP.numTimepoints;k++) {
			state2 = myDP.states[OTUNum][locationNum][k];
			totalState[state2]++;
			transitionCounts[state1][state2]++;
			state1 = state2;
		}
	}

	public void removeOTU(int OTUNum,int locationNum) {
		numDataAssigned--;
		int k = 0;
		int state1 = myDP.states[OTUNum][locationNum][0];
		int state2 = 0;
		initCounts[state1]--;
		for (k=1;k<myDP.numTimepoints;k++) {
			state2 = myDP.states[OTUNum][locationNum][k];
			totalState[state2]--;
			transitionCounts[state1][state2]--;
			state1 = state2;
		}
	}
	
	// compute the likelihood of the assigned data states (integrating out
	// the transition parameters)
	public double logStateLike() {
		double ll = 0.0;
		int s1 = 0;
		int s2 = 0;
		
		// calculate for initial state
		ll = -cern.jet.stat.Gamma.logGamma(((double) (numDataAssigned*myDP.numTimepoints)) + myDP.transitionPriorAlpha*((double) myDP.numStates));
		ll = ll - ((double) myDP.numStates)*Math.log((double)(numDataAssigned*myDP.numTimepoints));
		for (s1=0;s1<myDP.numStates;s1++) {
			ll += cern.jet.stat.Gamma.logGamma(((double) initCounts[s1])+myDP.transitionPriorAlpha);
			ll += Math.log(((double) initCounts[s1])+myDP.transitionPriorAlpha);
		}
		
		// now calculate for all states
		for (s1=0;s1<myDP.numStates;s1++) {
			ll = ll - cern.jet.stat.Gamma.logGamma(((double) totalState[s1]) + myDP.transitionPriorAlpha*((double) myDP.numStates));
			ll = ll - ((double) myDP.numStates)*Math.log((double)totalState[s1]); 
			for (s2=0;s2<myDP.numStates;s2++) {
				ll += cern.jet.stat.Gamma.logGamma(((double) transitionCounts[s1][s2])+myDP.transitionPriorAlpha);
				ll += Math.log(((double) transitionCounts[s2][s1])+myDP.transitionPriorAlpha);
			}
		}
		
		return ll;
	}
}