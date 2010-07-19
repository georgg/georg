package com.newmacondo.georg.DynamicDUO;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;

import cern.jet.random.Beta;
import cern.jet.random.Gamma;

import com.newmacondo.georg.StatUtil.PermutationGenerator;
import com.newmacondo.georg.StatUtil.StatUtil;

import umontreal.iro.lecuyer.probdist.NegativeBinomialDist;

public class DP implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = 5293768161600539738L;
	int numOTUs;
	int numLocations;
	int numTimepoints;
	int numReplicates;
	
	public String[] OTUNames;
	
	// dim 1 indexes OTUs, dim 2 location, dim 3 indexes time-points, dim 4 indexes replicates at each time-point
	public int[][][][] OTUCounts;
	
	// precomputed log likelihood for each observed value at each state
	public double[][][][] observedLoglike;
	
	// hidden states for each OTU
	public int[][][] states;
	
	// collects observations for each state
	public ArrayList<ArrayList<Integer>> collectedObservations = new ArrayList<ArrayList<Integer>>();
	
	// cluster assignments for each OTU
	public DynamicsCluster[][] OTUAssignment;
	
	// number of states for the HMM
	public int numStates = 4;
	
	// parameters for the negative binomial
	public double[] negBinomR;
	public double[] negBinomP;
	
	// hyper-parameters for the negative binomial value p (assumes a beta prior)
	public double negBinomHyper1;
	public double negBinomHyper2;
	
	// prior for transition probabilities (assumes symmetric Dirichlet prior)
	public double transitionPriorAlpha = 0.1;
	
	ArrayList<DynamicsCluster> clusters = new ArrayList<DynamicsCluster>();
	
	// the concentration parameter for the DP
	double alpha = 1.0;
	// parameters for the gamma prior on alpha
	double alpha_a = 10e-8;
	double alpha_b = 10e-8;
	// numbers of iterations for resampling alpha
	int numAlphaIters = 15;
	
	// precomputed state permutations
	int[][] permutations = null;
	
	public void init() {
		int i = 0;
		int j = 0;
		int s = 0;
		states = new int[numOTUs][numLocations][numTimepoints];
		OTUAssignment = new DynamicsCluster[numOTUs][numLocations];
		negBinomR = new double[numStates];
		negBinomP = new double[numStates];
		observedLoglike = new double[numOTUs][numLocations][numTimepoints][numStates];
		populatePermutations();
		
		DynamicsCluster cluster = new DynamicsCluster(this);
		
		for (i=0;i<numOTUs;i++) {
			for (j=0;j<numLocations;j++) {
				OTUAssignment[i][j] = cluster;
			}
		}
		
		for (s=0;s<numStates;s++) {
			collectedObservations.add(new ArrayList<Integer>());
		}
		
		clusters.add(cluster);
		cluster.numDataAssigned = numOTUs*numLocations;
		resampleTransitionProbs();
		initEmissionProbs();
		recomputeObservedLoglike();
		initStates();
		resampleTransitionProbs();
	}
	
	public void initStates() {
		int i = 0;
		int j = 0;
		int k = 0;
		int l = 0;
		int s = 0;
		int prevState = 0;
		DynamicsCluster cluster = null;
		ArrayList<Integer> stateObservations = null;
		double[] probs = new double[numStates];
		
		for (i=0;i<numOTUs;i++) {
			for (j=0;j<numLocations;j++) {
				cluster = OTUAssignment[i][j];
				for (s=0;s<numStates;s++) {
					probs[s] = Math.exp(observedLoglike[i][j][0][s]);
				}
				states[i][j][0] = StatUtil.multinomial_rnd(probs, numStates);
				cluster.initCounts[states[i][j][0]]++;
				stateObservations = collectedObservations.get(states[i][j][0]);
				for (l=0;l<numReplicates;l++) {
					stateObservations.add(OTUCounts[i][j][0][l]);
				}
				for (k=1;k<numTimepoints;k++) {
					prevState = states[i][j][k-1];
					for (s=0;s<numStates;s++) {
						probs[s] = Math.exp(observedLoglike[i][j][k][s]);
					}
					states[i][j][k] = StatUtil.multinomial_rnd(probs, numStates);
					cluster.totalState[states[i][j][k]]++;
					cluster.transitionCounts[prevState][states[i][j][k]]++;
					stateObservations = collectedObservations.get(states[i][j][k]);
					for (l=0;l<numReplicates;l++) {
						stateObservations.add(OTUCounts[i][j][k][l]);
					}
				}
			}
		}
	}
	
	
	// divide data into quantiles
	public void initEmissionProbs() {
		ArrayList<Integer> allData = new ArrayList<Integer>();
		int i = 0;
		int j = 0;
		int k = 0;
		int l = 0;
		for (i=0;i<numOTUs;i++) {
			for (j=0;j<numLocations;j++) {
				for (k=0;k<numTimepoints;k++) {
					for (l=0;l<numReplicates;l++) {
						allData.add(OTUCounts[i][j][k][l]);
					}
				}
			}
		}
		
		Collections.sort(allData);
		
		int inc = (int) Math.floor(((double) allData.size())/4.0);
		i = 0;
		j = 0;
		k = 0;
		double p = 0.0;
		double r = 0.0;
		int[] v = new int[inc];
		int n = Math.min(inc*4,allData.size());
		double[] params;
		while(i<n) {
			v[j] = allData.get(i);
			i++;
			j++;
			if (j==inc) {
				j = 0;
				params = NegativeBinomialDist.getMLE(v, inc);
				negBinomR[k] = params[0];
				negBinomP[k] = params[1];
				p = negBinomP[k];
				r = negBinomR[k];
				System.out.print("mu");
				System.out.print(k);
				System.out.print("=");
				System.out.print(r*(1.0-p)/p);
				System.out.print(" r");
				System.out.print(k);
				System.out.print("=");
				System.out.print(r);
				System.out.print("\n");
				k++;
			}
		}
	}
	
	public void populatePermutations() {
		int[] indices;
		int numPerm = PermutationGenerator.getSmallFactorial(numStates);
		int i = 0;
		int pn = 0;
		permutations = new int[numPerm][numStates];
		PermutationGenerator x = new PermutationGenerator(numStates);
		while(x.hasMore()) {
			indices = x.getNext();
			for (i=0;i<numStates;i++) {
				permutations[pn][i] = indices[i];
			}
			pn++;
		}
	}
		
	public void iterate(int numIters) {
		int i = 0;
		int nc = 0;
		int cc = 0;
		DynamicsCluster cluster = null;
		int outputFreq = 100;
		int numSinceOutput = 0;
		double r = 0.0;
		double p = 0.0;
		
		for (i=0;i<numIters;i++) {
			sampleOTUAssigns();
			
			nc = clusters.size();
  			for (cc=nc-1;cc>=0;cc--) {
  				cluster = clusters.get(cc);
  				if (cluster.numDataAssigned == 0) {
  					clusters.remove(cc);
  				}
  			}
			
  			resampleStates();
  			resampleTransitionProbs();
  			resampleEmissionProbs();
  			recomputeObservedLoglike();
  			resampleAlpha();
  			
  			if (numSinceOutput == outputFreq) {
  				System.out.print("iter=");
  				System.out.print(i);
  				System.out.print(" # clusters=");
  				System.out.print(clusters.size());
  				for (cc=0;cc<numStates;cc++) {
  					p = negBinomP[cc];
  					r = negBinomR[cc];
  					System.out.print(" mu");
  					System.out.print(cc);
  					System.out.print("=");
  					System.out.print(r*(1.0-p)/p);
  				}
  				System.out.print("\n");
  				numSinceOutput = 0;
  			}
  			
  			numSinceOutput++;
		}
	
	}
	
	void sampleOTUAssigns() {
		int cc = 0;
		int i = 0;
		int j = 0;
	//	int[][] transitionCounts = new int[numStates][numStates];
	//	int totalState[] = new int[numStates];
		
		double denom = (((double) (numOTUs*numLocations)) - 1.0 + alpha);
		double tmp = 0.0;
		double ll = 0.0;
		
		DynamicsCluster newCluster = new DynamicsCluster(this);
		clusters.add(newCluster);
		newCluster.sampleTransitionProbsPosterior();
		DynamicsCluster cluster = null;
		ArrayList<Double> condLike = new ArrayList<Double>();
		for (cc=0;cc<clusters.size();cc++) {
			condLike.add(0.0);
		}
		
		for (i=0;i<numOTUs;i++) {
			for (j=0;j<numLocations;j++) {
				cluster = OTUAssignment[i][j];
				cluster.removeOTU(i, j);
				for (cc=0;cc<clusters.size();cc++) {
					cluster = clusters.get(cc);
					
				/*	if (cc == clusters.size()-1) {
						ll = logDataLikeEmptyCluster(i, j, transitionCounts,totalState);
					} else {
						ll = logDataLikeCluster(i, j, cluster);
					} */
					ll = logDataLikeCluster(i, j, cluster);
					
					if (cluster.numDataAssigned>=1) {
						tmp = ((double) cluster.numDataAssigned);
					} else {
						tmp = alpha;
					}
					ll = Math.exp(ll)*tmp/denom;
					if (Double.isNaN(ll)) {
						ll = 0.0;
					}
					condLike.set(cc, ll);
				}
			
				// now choose the new cluster
				cc = StatUtil.multinomial_rnd(condLike);
				cluster = clusters.get(cc);
				OTUAssignment[i][j] = cluster;
				cluster.addOTU(i,j);
				
				if (cc==clusters.size() - 1) {
					newCluster = new DynamicsCluster(this);
					newCluster.sampleTransitionProbsPosterior();
					clusters.add(newCluster);
					condLike.add(0.0);
				} else {
					clusters.get(clusters.size()-1).sampleTransitionProbsPosterior();
				}
			}
		}
	}
		
	void resampleStates() {
		int i = 0;
		int j = 0;
		int k = 0;
		int l = 0;
		int s = 0;
		double[] condProb = new double[numStates];
		DynamicsCluster cluster = null;
		int prevState = 0;
		ArrayList<Integer> stateObservations = null;
		// clear out the observations, as they'll be reassigned
		for (s=0;s<numStates;s++) {
			collectedObservations.get(s).clear();
		}
		
		for (i=0;i<numOTUs;i++) {
			for (j=0;j<numLocations;j++) {
				cluster = OTUAssignment[i][j];
				for (s=0;s<numStates;s++) {
					condProb[s] = cluster.initProbs[s] * Math.exp(observedLoglike[i][j][0][s]);
				}
				cluster.initCounts[states[i][j][0]]--;
				states[i][j][0] = StatUtil.multinomial_rnd(condProb, numStates);
				cluster.initCounts[states[i][j][0]]++;
				stateObservations = collectedObservations.get(states[i][j][0]);
				for (l=0;l<numReplicates;l++) {
					stateObservations.add(OTUCounts[i][j][0][l]);
				}
				for (k=1;k<numTimepoints;k++) {
					prevState = states[i][j][k-1];
					for (s=0;s<numStates;s++) {
						condProb[s] = cluster.transitionProbs[prevState][s] * Math.exp(observedLoglike[i][j][k][s]);
					}
					cluster.transitionCounts[prevState][states[i][j][k]]--;
					cluster.totalState[states[i][j][k]]--;
					states[i][j][k] = StatUtil.multinomial_rnd(condProb, numStates);
					cluster.transitionCounts[prevState][states[i][j][k]]++;
					cluster.totalState[states[i][j][k]]++;
					stateObservations = collectedObservations.get(states[i][j][k]);
					for (l=0;l<numReplicates;l++) {
						stateObservations.add(OTUCounts[i][j][k][l]);
					}
				}
			}
		}
	}
	
	public void resampleTransitionProbs() {
		int cc = 0;
		DynamicsCluster cluster = null;
		for (cc=0;cc<clusters.size();cc++) {
			cluster = clusters.get(cc);
			cluster.sampleTransitionProbsPosterior();
		}
	}
	
	// this needs to be updated to do sampling (will need to do Metropolis steps
	// for the nonconjugate gamma prior on r
	// at the moment, it just computes MLE estimates
	public void resampleEmissionProbs() {
		int s = 0;
		int i = 0;
		int[] observed = null;
		double[] params = null;
		ArrayList<Integer> stateObservations = null;
		for (s=0;s<numStates;s++) {
			stateObservations = collectedObservations.get(s);
			if (!stateObservations.isEmpty()) {
				observed = new int[stateObservations.size()];
				for (i=0;i<observed.length;i++) {
					observed[i] = stateObservations.get(i);
				}
				params = NegativeBinomialDist.getMLE(observed, observed.length);
			}
			else {
				params = new double[2];
				params[0] = 5;
				params[1] = StatUtil.beta_rnd(5+params[0], 5);
			}
			negBinomR[s] = params[0];
			negBinomP[s] = params[1];
		}
	}
	
	// sample concentration parameter using aux. variable sampling scheme
	// from Escobar and West (1995)
	void resampleAlpha() {
		int i = 0;
		double eta = 0.0;
		int zz = 0;
		double[] w = new double[2];
		int numClusters = clusters.size();
		int numData = numOTUs * numLocations;
		if (numClusters <= 1) {
			return;
		}
		
		for (i=0;i<numAlphaIters;i++) {
			eta = Beta.staticNextDouble(alpha+1.0,numData); 
			if (eta > 0) {
				w[0] = (alpha_a + (double) (numClusters - 1))/(alpha_b - Math.log(eta));
				w[1] = (double) numData;
				zz = StatUtil.multinomial_rnd(w,2);
				w[0] = alpha_b - Math.log(eta);
				if (w[0] >= 0.0) {
					if (zz == 0) {
						alpha = Gamma.staticNextDouble(alpha_a+((double) numClusters),w[0]);
					} else {
						alpha = Gamma.staticNextDouble(alpha_a+((double) numClusters-1),w[0]);
					}
				}
			}
		}
	}
	
	void recomputeObservedLoglike() {
		int i = 0;
		int j = 0;
		int k = 0;
		int l = 0;
		int s = 0;
		for (i=0;i<numOTUs;i++) {
			for (j=0;j<numLocations;j++) {
				for (k=0;k<numTimepoints;k++) {
					for (s=0;s<numStates;s++) {
						observedLoglike[i][j][k][s] = 0.0;
						for (l=0;l<numReplicates;l++) {
							observedLoglike[i][j][k][s] += Math.log(NegativeBinomialDist.prob(negBinomR[s], negBinomP[s], OTUCounts[i][j][k][l]));
						}
					}
				}
			}
		}
	}
	
	// compute log likelihood for data, given state sequence and putative cluster
	public double logDataLikeCluster(int OTUNum, int OTULocation, DynamicsCluster cluster) {
		double ll = 0.0;
		int k = 0;
		int state1 = states[OTUNum][OTULocation][0];
		int state2 = 0;
		double prob = 0.0;
		
		prob = cluster.initProbs[state1];
		ll += Math.log(prob);
		
		for (k=1;k<numTimepoints;k++) {
			state2 = states[OTUNum][OTULocation][k];
			prob = cluster.transitionProbs[state1][state2];
			ll += Math.log(prob);
			state1 = state2;
		}
		
		return ll;
	}
	
	// compute log likelihood for data, given state sequence,
	// integrating out transition probabilities
	public double logDataLikeEmptyCluster(int OTUNum, int OTULocation, int[][] transitionCounts,int[] totalState) {
		double ll = 0.0;
		int s1 = 0;
		int s2 = 0;
		int k = 0;
		int initCounts = 0;
		
		// calculate for initial state
		ll = -cern.jet.stat.Gamma.logGamma(((double) numTimepoints) + transitionPriorAlpha*((double) numStates));
		ll = ll - ((double) numStates)*Math.log((double)numTimepoints + transitionPriorAlpha);
		for (s1=0;s1<numStates;s1++) {
			if (s1 == states[OTUNum][OTULocation][0]) {
				initCounts = 1;
			} else {
				initCounts = 0;
			}
			ll += cern.jet.stat.Gamma.logGamma(((double) initCounts)+transitionPriorAlpha);
			ll += Math.log(((double) initCounts)+transitionPriorAlpha);
		}
		
		for (s1=0;s1<numStates;s1++) {
			totalState[s1] = 0;
			for (s2=0;s2<numStates;s2++) {
				transitionCounts[s1][s2] = 0;
			}
		}
		for (k=0;k<numTimepoints;k++) {
			totalState[states[OTUNum][OTULocation][k]]++;
			if (k>0) {
				transitionCounts[states[OTUNum][OTULocation][k-1]][states[OTUNum][OTULocation][k]]++;
			}
		}
		
		// now calculate for all states
		for (s1=0;s1<numStates;s1++) {
			ll = ll - cern.jet.stat.Gamma.logGamma(((double) totalState[s1]) + transitionPriorAlpha*((double) numStates));
			ll = ll - ((double) numStates)*Math.log((double)totalState[s1] + transitionPriorAlpha); 
			for (s2=0;s2<numStates;s2++) {
				ll += cern.jet.stat.Gamma.logGamma(((double) transitionCounts[s1][s2])+transitionPriorAlpha);
				ll += Math.log(((double) transitionCounts[s2][s1])+transitionPriorAlpha);
			}
		}
		
		return ll;
	}
	
	public void writeClusters(String fname) throws IOException {
		int s1 = 0;
		int s2 = 0;
		int cc = 0;
		double ml = 0.0;
		int bestState = 0;
		int prevState = 0;
		int k = 0;
		FileWriter outFile = new FileWriter(fname);
		String s;
		double p = 0.0;
		double r = 0.0;
		DynamicsCluster cluster = null;
		
		for (s1=0;s1<numStates;s1++) {
			p = negBinomP[s1];
			r = negBinomR[s1];
			s = "mu" + Integer.toString(s1+1) + " = " + Double.toString(r*(1.0-p)/p);
			s = s + "     " + "r" + Integer.toString(s1+1) + " = " + Double.toString(r);
			outFile.write(s+"\n");
		}
		
		outFile.write("\n");
		
		for (cc=0;cc<clusters.size();cc++) {
			cluster = clusters.get(cc);
			s = "DynamicsCluster#" + Integer.toString(cc+1);
			s = s + "\t" + "# OTUs = " + Integer.toString(cluster.numDataAssigned);
			outFile.write(s+"\n");
			outFile.write("MLS");
			ml = 0.0;
			for (s1=0;s1<numStates;s1++) {
				if (cluster.initProbs[s1] > ml) {
					ml = cluster.initProbs[s1];
					bestState = s1;
				}
			}
			prevState = bestState;
			outFile.write("\t"+Integer.toString(bestState+1));
			for (k=1;k<numTimepoints;k++) {
				ml = 0.0;
				for (s1=0;s1<numStates;s1++) {
					if (cluster.transitionProbs[prevState][s1] > ml) {
						ml = cluster.transitionProbs[prevState][s1];
						bestState = s1;
					}
				}
				outFile.write("\t" + Integer.toString(bestState+1));
				prevState = bestState;
			}
			outFile.write("\n");
			outFile.write("TPM\n");
			for (s2=0;s2<numStates;s2++) {
				outFile.write("\t" + Double.toString(cluster.initProbs[s2]));
				for (s1=0;s1<numStates;s1++) {
					outFile.write("\t" + Double.toString(cluster.transitionProbs[s1][s2]));
				}
				outFile.write("\n");
			}
			outFile.write("\n\n");
		}
		outFile.close();
	}
}
