package edu.mit.csail.psrg.georg.BindDP;

import java.util.ArrayList;
import cern.jet.random.Gamma;
import cern.jet.random.Uniform;
import cern.jet.random.Exponential;
import cern.jet.random.Beta;
import edu.mit.csail.psrg.georg.SampleUtil.SliceSampler;
import edu.mit.csail.psrg.georg.StatUtil.StatUtil;

public class BindDP {
	int numGenes = 0;
	int numExperiments = 0;
	int numClusters = 0;
	// allocated size of the cluster array
	int allocClusters = 5;
	// amount to grow cluster allocation by
	int growClusterAlloc = 10;
	
	// the largest possible p-value to consider as a binding event
	double possibleBindThreshold = 0.10;
	
	// the hidden binding variable for each gene
	boolean[][] inferredBinding = null;
	// for storing the mean of the sampled inferred binding values
	double[][] meanInferredBinding = null;
	// p-values for each gene
	double[][] pvals = null;
	// indices of (possible) binding events for each gene
	int[][] geneBindIdx = null;
	// array of genes with NaN values (there will be duplicates)
	int[] NaNGenes = null;
	// array of experiments with NaN values (corresponding to NaNGenes
	int[] NaNExperiments = null;
	
	// the concentration parameter for clustering
	double alpha = 10.0;
	// parameters for the gamma prior on alpha
	double alpha_a = 0.5;	// 10e-10;
	double alpha_b = 0.5;	// 10e-10;
	
	// number of genes assigned to each cluster
	int[] numGenesInCluster = null;
	// number of bound genes in each cluster and experiment
	int[][] numBound = null;
	// lambda parameters for cluster Binomial distributions
	double[][] lambda = null;
	// cluster assignment for each gene
	int[] clusterAssign = null;
	// internal variable used for computing conditional likelihood
	double[] condLike = null;
	// internal array for storing temporary cluster statistics
	double[] clusterTempStat = null;
	
	// parameters for beta prior on binding (shared across all clusters)
//	double priorNoBind = 95;
//	double priorBind = 5;
	double priorNoBind = 0.99;
	double priorBind = 0.01;
	// parameters for exponential densities (for generating p-values); the
	// vector is of dimension = # binding experiments
	double[] beta = null;
	// internval vectors used for sampling beta
	double[] Nj = null;
	double[] sbj = null;
	// vector to hold mean of beta posterior
	double[] meanBeta = null;
	// parameters for gamma prior on exponential densities
	// alpha_beta is like an 'effective sample size' for binding events
	double alpha_beta = 0.0;
	double w_beta = 0.0;
	
	// parameter for the prior on alpha_beta and w_beta
	double K = 0.0;
	
	// used to store the probability of two genes co-occurring in a cluster
	double[][] pairProb;
	
	// bindThreshold is an initial guess at the p-value threshold for
	// binding; nc is the initial number of clusters
	BindDP(double[][] p, double bindThreshold, int nc) {
		pvals = p;
		numGenes = pvals.length;
		numExperiments = pvals[0].length;
		clusterAssign = new int[numGenes];
		beta = new double[numExperiments];
		Nj = new double[numExperiments];
		sbj = new double[numExperiments];
		meanBeta = new double[numExperiments];
		
		compressGenes();
		flagNaNs();
		initParameters(bindThreshold);
		
		numClusters = nc;
		initAllocClusters();
		initClusterAssign();
	}
	
	void flagNaNs() {
		int i = 0;
		int j = 0;
		int n = 0;
		ArrayList<Integer> ng = new ArrayList<Integer>();
		ArrayList<Integer> ne = new ArrayList<Integer>();
		
		for(i=0;i<numGenes;i++) {
			if (pvals[i] != null) {
				n = pvals[i].length;
				for(j=0;j<n;j++) {
					if (Double.isNaN(pvals[i][j])) {
						ng.add(i);
						ne.add(j);
					}
				}
			}
		}
		
		if (!ng.isEmpty()) {
			NaNGenes = new int[ng.size()];
			NaNExperiments = new int[ng.size()];
			for (i=0;i<ng.size();i++) {
				NaNGenes[i] = ng.get(i);
				NaNExperiments[i] = ne.get(i);
			}
		}
	}
	
	void compressGenes() {
		int i = 0;
		int j = 0;
		int n = 0;
		int idx = 0;
		ArrayList<Integer> tempArray = new ArrayList<Integer>();
		geneBindIdx = new int[numGenes][];
		inferredBinding = new boolean[numGenes][];
		meanInferredBinding = new double[numGenes][];
		double[][] pvals2 = new double[numGenes][];
		
		for (i=0;i<numGenes;i++) {
			tempArray.clear();
			for (j=0;j<numExperiments;j++) {
				if (pvals[i][j] < possibleBindThreshold || Double.isNaN(pvals[i][j])) {
					tempArray.add(j);
				}
			}
			n = tempArray.size();
			if (n > 0) {
				geneBindIdx[i] = new int[n];
				pvals2[i] = new double[n];
				inferredBinding[i] = new boolean[n];
				meanInferredBinding[i] = new double[n];
				for (j=0;j<n;j++) {
					idx = tempArray.get(j);
					geneBindIdx[i][j] = idx;
					pvals2[i][j] = pvals[i][idx];
				}
			} else {
				geneBindIdx[i] = null;
				pvals2[i] = null;
				inferredBinding[i] = null;
			}
		}
		pvals = pvals2;
	}
	
	void initParameters(double bindThreshold) {
		int i = 0;
		int j = 0;
		int n = 0;
		w_beta = 0.0;
		double numK = 0.0;
		double p = 0.0;
		
		for (i=0;i<numGenes;i++) {
			if (pvals[i] != null) {
				n = pvals[i].length;
				for (j=0;j<n;j++) {
					p = pvals[i][j];
					if (!Double.isNaN(p) && p <= bindThreshold) {
						inferredBinding[i][j] = true;
						w_beta = w_beta + p;
						numK = numK + 1.0;
					}
					else {
						inferredBinding[i][j] = false;
					}
				}
			}
		}
		
		// this corresponds to the average p-value for a bound gene across all
		// genes and experiments
		w_beta = w_beta/numK;
		K = w_beta;
		// estimate alpha_beta as the average number of binding events per transcription factor
		alpha_beta = numK/((double) numExperiments);
		//alpha_beta = 1.0;
		
		priorBind = numK/((double) numExperiments*numGenes);
		priorNoBind = 1.0 - priorBind;
				
		// now sample beta values from the prior
/*		for (i=0;i<numExperiments;i++) {
			beta[i] = Gamma.staticNextDouble(alpha_beta/2.0,(w_beta*alpha_beta)/2.0);
		} */
		// set beta values to the prior mean, which is 1/K or 1/average p-value for a bound gene
		for (i=0;i<numExperiments;i++) {
			beta[i] = 1.0/K;
		}
		
		// fill-in initial guesses for missing values
		double t = 0.0;
		int gene = 0;
		int expt = 0;
		if (NaNGenes != null) {
			for (i=0;i<NaNGenes.length;i++) {
				gene = NaNGenes[i];
				expt = NaNExperiments[i];
				t = Uniform.staticNextDouble();
				if (t > priorBind) {
					inferredBinding[gene][expt] = false;
				//	pvals[gene][expt] = Uniform.staticNextDouble();
				}
				else {
					inferredBinding[gene][expt] = true;
				//	pvals[gene][expt] = Exponential.staticNextDouble(beta[expt]);
				}
			}
		}
	}
	
	void initClusterAssign() {
		int i = 0;
		int j = 0;
		int cc = 0;
		int n = 0;
		int idx = 0;
		
		for (i=0;i<numGenes;i++) {
			cc = Uniform.staticNextIntFromTo(0,numClusters-1);
			numGenesInCluster[cc]++;
			clusterAssign[i] = cc;
			if (pvals[i] != null) {
				n = pvals[i].length;
				for (j=0;j<n;j++) {
					if (inferredBinding[i][j]) {
						idx = geneBindIdx[i][j];
						numBound[cc][idx]++;
					}
				}
			}
		}
		
		sampleLambda();
	}
	
	void initAllocClusters() {
		int totalAllocate = allocClusters;
		if (allocClusters < numClusters+2) {
			totalAllocate = numClusters+growClusterAlloc;
		}
		
		numGenesInCluster = new int[totalAllocate];
		numBound = new int[totalAllocate][numExperiments];
		lambda = new double[totalAllocate][numExperiments];
		condLike = new double[totalAllocate];
		clusterTempStat = new double[totalAllocate];
				
		allocClusters = totalAllocate;
	}
	
	void addCluster() {
		if (numGenesInCluster.length < numClusters+2) {
			int[] numGenesInCluster2 = new int[allocClusters+growClusterAlloc];
			System.arraycopy(numGenesInCluster, 0, numGenesInCluster2, 0, allocClusters);
			numGenesInCluster = numGenesInCluster2;
			
			int[][] numBound2 = new int[allocClusters+growClusterAlloc][numExperiments];
			System.arraycopy(numBound, 0, numBound2, 0, allocClusters);
			numBound = numBound2;
			
			double[][] lambda2 = new double[allocClusters+growClusterAlloc][numExperiments];
			System.arraycopy(lambda, 0, lambda2, 0, allocClusters);
			lambda = lambda2;
			
			double[] condLike2 = new double[allocClusters+growClusterAlloc];
			System.arraycopy(condLike,0,condLike2,0,allocClusters);
			condLike = condLike2;
			
			double[] clusterTempStat2 = new double[allocClusters+growClusterAlloc];
			System.arraycopy(clusterTempStat,0,clusterTempStat2,0,allocClusters);
			clusterTempStat = clusterTempStat2;
			
			allocClusters = numGenesInCluster.length;
		}
		
		numGenesInCluster[numClusters+1] = 0;
		clusterTempStat[numClusters+1] = 0;
		
		int k = 0;
		for (k=0;k<numExperiments;k++) {
			numBound[numClusters+1][k] = 0;
			lambda[numClusters+1][k] = 0.0;
		}
		
		numClusters++;
	}
	
	void deleteCluster(int c) {
		// the cluster to be deleted is zero-indexed, as is standard in Java
		// this method will copy numClusters+1 entries
		
		System.arraycopy(numGenesInCluster,c+1,numGenesInCluster,c,numClusters-c);
		System.arraycopy(clusterTempStat,c+1,clusterTempStat,c,numClusters-c);
		
		int k = 0;
		
		int[] pc = numBound[c];
		for(k=c;k<numClusters;k++) {
			numBound[k] = numBound[k+1];
		}
		numBound[numClusters] = pc;
		
		double[] ls = lambda[c];
		for(k=c;k<numClusters;k++) {
			lambda[k] = lambda[k+1];
		}
		lambda[numClusters] = ls;
		
		for (k=0;k<clusterAssign.length;k++) {
			if (clusterAssign[k] > c) {
				clusterAssign[k] = clusterAssign[k] - 1;
			}
		}
		
		numClusters--;
	}
	
	void updateClusterTempStat() {
		int cc = 0;
		int j = 0;
		
		for (cc=0;cc<numClusters;cc++) {
			clusterTempStat[cc] = 0.0;
			for(j=0;j<numExperiments;j++) {
				clusterTempStat[cc] = clusterTempStat[cc] + Math.log(1.0-lambda[cc][j]);
			}
		}
	}
	
	void sampleGeneClusterAssigns() {
		int cc = 0;
		int oldCC = 0;
		int i = 0;
		int j = 0;
		int n = 0;
		int idx = 0;
		double denom = ((double) numGenes - 1.0 + alpha);
		double pBind = priorBind/(priorBind+priorNoBind);
		double pBindStat = 0.0;
		double tmp = 0.0;
		
		updateClusterTempStat();
		for (j=0;j<numExperiments;j++) {
			pBindStat = pBindStat + Math.log(1.0-pBind);
		}
		
		for (i=0;i<numGenes;i++) {
			oldCC = clusterAssign[i];
			if (inferredBinding[i] != null) {
				n = inferredBinding[i].length;
			} else { n = 0; }
			for (cc=0;cc<numClusters;cc++) {
				condLike[cc] = clusterTempStat[cc];
				if (n>0) {
					for (j=0;j<n;j++) {
						idx = geneBindIdx[i][j];
						if (inferredBinding[i][j]) {
							condLike[cc] = condLike[cc] + Math.log(lambda[cc][idx]) - Math.log(1.0-lambda[cc][idx]);
						}
					}
				}
				
				// we're considering a cluster that the gene is not already assigned to
				if (cc != oldCC) {
					condLike[cc] = Math.exp(condLike[cc])*((double) numGenesInCluster[cc])/denom;
				} else {
					if (numGenesInCluster[cc] > 1) {
						tmp = ((double) numGenesInCluster[cc] - 1.0);
					} else {
						tmp = alpha;
					}
					condLike[cc] = Math.exp(condLike[cc])*tmp/denom;
				}
			}
			condLike[numClusters] = pBindStat;
			if (n>0) {
				for (j=0;j<n;j++) {
					idx = geneBindIdx[i][j];
					if (inferredBinding[i][j]) {
						condLike[numClusters] = condLike[numClusters] + Math.log(pBind) - Math.log(1.0-pBind);
					}
				}
			}
			condLike[numClusters] = Math.exp(condLike[numClusters])*alpha/denom;
		
			// now choose the new cluster
			cc = StatUtil.multinomial_rnd(condLike,numClusters+1);
		
			if (cc != oldCC) {
				numGenesInCluster[oldCC]--;
				numGenesInCluster[cc]++;
				if (n>0) {
					for (j=0;j<n;j++) {
						idx = geneBindIdx[i][j];
						if (inferredBinding[i][j]) {
							numBound[cc][idx]++;
							numBound[oldCC][idx]--;
						}
					}
				}
				clusterAssign[i] = cc;
			}
			
			if (cc==numClusters) {
				addCluster();
				clusterTempStat[cc] = 0.0;
				for (j=0;j<numExperiments;j++) {
					lambda[cc][j] = Beta.staticNextDouble(priorBind + (double) numBound[cc][j],priorNoBind + (double) (numGenesInCluster[cc] - numBound[cc][j]));
					clusterTempStat[cc] = clusterTempStat[cc] + Math.log(1.0-lambda[cc][j]);
				}
			}
		}
	}
	
	void sampleLambda() {
		int cc = 0;
		int j = 0;
		
		for(cc=0;cc<numClusters;cc++) {
			for(j=0;j<numExperiments;j++) {
				lambda[cc][j] = Beta.staticNextDouble(priorBind + (double) numBound[cc][j],priorNoBind + (double) (numGenesInCluster[cc] - numBound[cc][j]));
			}
		}
	}
	
	void sampleInferredBinding() {
		int i = 0;
		int j = 0;
		int cc = 0;
		int n = 0;
		int idx = 0;
		double tmp = 0.0;
		double tmp2 = 0.0;
		
		for (cc=0;cc<numClusters;cc++) {
			for (j=0;j<numExperiments;j++) {
				numBound[cc][j] = 0;
			}
		}
		
		for (i=0;i<numGenes;i++) {
			cc = clusterAssign[i];
			if (inferredBinding[i] != null) {
				n = inferredBinding[i].length;
				for (j=0;j<n;j++) {
					idx = geneBindIdx[i][j];
					tmp = 0.0;
					tmp2 = 0.0;
					tmp = lambda[cc][idx];
					tmp2 = 1.0 - tmp;
					if (!Double.isNaN(pvals[i][j])) {
						tmp = tmp*beta[idx]*Math.exp(-beta[idx]*pvals[i][j]);
					}
					tmp = tmp/(tmp+tmp2);
					tmp2 = Uniform.staticNextDouble();
					if (tmp2 <= tmp) {
						inferredBinding[i][j] = true;
					} else {
						inferredBinding[i][j] = false;
					}
					if (inferredBinding[i][j]) {
						numBound[cc][idx]++;
					}
				}
			}
		}
		
	//	samplePriorBind();
		estimatePriorBind();
	}
	
/*	void sampleMissingValues() {
		int i = 0;
		int gene = 0;
		int expt = 0;
		if (NaNGenes != null) {
			for (i=0;i<NaNGenes.length;i++) {
				gene = NaNGenes[i];
				expt = NaNExperiments[i];
				if (!inferredBinding[gene][expt]) {
					pvals[gene][expt] = Uniform.staticNextDouble();
				}
				else {
					pvals[gene][expt] = Exponential.staticNextDouble(beta[expt]);
				}
			}
		}
	} */
	
	// sample parameters for exponential binding distributions and priors
	void sampleBeta() {
		int i = 0;
		int j = 0;
		int n = 0;
		int idx = 0;
		double sBeta = 0.0;
		
		for (j=0;j<numExperiments;j++) {
			Nj[j] = 0;
			sbj[j] =0;
		}
		
		for (i=0;i<numGenes;i++) {
			if (inferredBinding[i] != null) {
				n = inferredBinding[i].length;
				for (j=0;j<n;j++) {
					idx = geneBindIdx[i][j];
					if (!Double.isNaN(pvals[i][j]) && inferredBinding[i][j]) {
						Nj[idx] = Nj[idx] + 1.0;
						sbj[idx] = sbj[idx] + pvals[i][j];
					}
				}
			}
		}
		
		for (j=0;j<numExperiments;j++) {
			beta[j] = Gamma.staticNextDouble(Nj[j]+alpha_beta/2.0,alpha_beta*w_beta/2.0 + sbj[j]);
			sBeta = sBeta + beta[j];
		}
		
		double M = (double) numExperiments;
	//	w_beta = Gamma.staticNextDouble(M*alpha_beta/2.0 + 0.5,1.0/(2.0*K) + alpha_beta*sBeta/2.0);
	//	System.out.println("w_beta= " + w_beta);
			
	//	sampleAlpha_beta();
	}
	
	void sampleAlpha_beta() {
	//	if (alpha_beta <= 1.0) {
	//		return;
	//	}
		
		Alpha_BetaPriorPosteriorLogDensityFn myAlpha_BetaPriorPosteriorLogDensityFn = new Alpha_BetaPriorPosteriorLogDensityFn();
		myAlpha_BetaPriorPosteriorLogDensityFn.setParams(w_beta,beta);
		double[] x = new double[1];
		x[0] = alpha_beta;
		double[] wp = new double[1];
		wp[0] = 0.5;
		double[][] S = new double[1][1];
		SliceSampler.sample(1,500,myAlpha_BetaPriorPosteriorLogDensityFn,x,wp,S);
		if (!Double.isInfinite(S[0][0])) {
			alpha_beta = S[0][0];
			System.out.println("alpha_beta=" + S[0][0]);
		}
	}
	
	// set prior binding parameters using method of moments,
	// which will be much faster than sampling it
	void estimatePriorBind() {
		double lambdaMean = 0.0;
		double lambdaVar = 0.0;
		int cc = 0;
		int j = 0;
		
		for (cc=0;cc<numClusters;cc++) {
			for (j=0;j<numExperiments;j++) {
				lambdaMean = lambdaMean + lambda[cc][j];
				lambdaVar = lambdaVar + Math.pow(lambda[cc][j],2.0);
			}
		}
		
		double tot = (double) numClusters*numExperiments;
		lambdaMean = lambdaMean/tot;
		lambdaVar = lambdaVar/tot - Math.pow(lambdaMean,2.0);
		
		// estimate beta distribution parameters
		tot = (lambdaMean*(1.0-lambdaMean)/lambdaVar) - 1.0;
		priorBind = lambdaMean*tot;
		priorNoBind = (1.0-lambdaMean)*tot;
		System.out.println("priorBind=" + priorBind + " " + "priorNoBind=" + priorNoBind);
	}
	
	void samplePriorBind() {
		BindCountPriorPosteriorLogDensityFn myBindCountPriorPosteriorLogDensity = new BindCountPriorPosteriorLogDensityFn();
		myBindCountPriorPosteriorLogDensity.setParams(numGenesInCluster,numBound,lambda,numClusters);
		double[] x = new double[2];
		x[0] = priorBind;
		x[1] = priorNoBind;
		double[] wp = new double[2];
		wp[0] = 0.5;
		wp[1] = 0.5;
		double[][] S = new double[1][2];
		SliceSampler.sample(1,500,myBindCountPriorPosteriorLogDensity,x,wp,S);
		if (!Double.isInfinite(S[0][0]) && !Double.isInfinite(S[0][1])) {
			priorBind = S[0][0];
			priorNoBind = S[0][1];
			System.out.println("priorBind=" + S[0][0] + " " + "priorNoBind=" + S[0][1]);
		}
	}
	
	// sample concentration parameter using aux. variable sampling scheme
	// from Escobar and West (1995)
	void sampleAlpha(int numIters) {
		int i = 0;
		double eta = 0.0;
		double zz = 0.0;
		
		for (i=0;i<numIters;i++) {
			eta = Beta.staticNextDouble(alpha,numGenes);
			if (eta > 0) {
				zz = alpha_b - Math.log(eta);
				if (zz > 0) {
					alpha = Gamma.staticNextDouble(alpha_a+numClusters-1,zz);
				}
			}
		}
	}
	
	void iterate(int numIters,int numBurnin,int numAlphaIters) {
		int i = 0;
		int nc = 0;
		int cc = 0;
		for (i=0;i<numIters;i++) {
			sampleGeneClusterAssigns();
			
			nc = numClusters;
  			for (cc=nc-1;cc>=0;cc--) {
  				if (numGenesInCluster[cc] == 0) {
  					deleteCluster(cc);
  				}
  			}
			
  			sampleLambda();
			sampleInferredBinding();
		//	sampleMissingValues();
			sampleBeta();
			sampleAlpha(numAlphaIters);
			
			if (i >= numBurnin) {
				incMeans();
			}
			
			System.out.println(i + " " + numClusters + " " + alpha);
		}
		
		normMeans(numIters - numBurnin);
	}
	
	void incMeans() {
		int i = 0;
		int j = 0;
		int n = 0;
		
		for (i=0;i<numGenes;i++) {
			if (inferredBinding[i] != null) {
				n = inferredBinding[i].length;
				for (j=0;j<n;j++) {
					if (inferredBinding[i][j]) {
						meanInferredBinding[i][j] = meanInferredBinding[i][j] + 1.0;
					}
				}
			}
		}
		
		for (j=0;j<numExperiments;j++) {
			meanBeta[j] = meanBeta[j] + beta[j];
		}
	}
	
	void normMeans(int ns) {
		double norm = (double) ns;
		int i = 0;
		int j = 0;
		int n = 0;
		
		for (i=0;i<numGenes;i++) {
			if (inferredBinding[i] != null) {
				n = inferredBinding[i].length;
				for (j=0;j<n;j++) {
					meanInferredBinding[i][j] = meanInferredBinding[i][j]/norm;
				}
			}
		}
		
		for (j=0;j<numExperiments;j++) {
			meanBeta[j] = meanBeta[j]/norm;
		}
	}
	
	void expandMeanInferredBinding() {
		int n = 0;
		int idx = 0;
		int i = 0;
		int j = 0;
		double[][] meanInferredBinding2 = new double[numGenes][numExperiments];
		for (i=0;i<numGenes;i++) {
			if (inferredBinding[i] != null) {
				n = inferredBinding[i].length;
				for (j=0;j<n;j++) {
					idx = geneBindIdx[i][j];
					meanInferredBinding2[i][idx] = meanInferredBinding[i][j];
				}
			}
		}
		meanInferredBinding = meanInferredBinding2;
	}
	
	void outputClusters() {
		int cc = 0;
		int j = 0;
		for (cc=0;cc<numClusters;cc++) {
			System.out.println("Cluster #" + cc + " " + numGenesInCluster[cc] + " genes");
			System.out.print("[");
			for (j=0;j<numExperiments;j++) {
				System.out.print(" " + lambda[cc][j]);
			}
			System.out.println(" ]");
		}
		
		System.out.println("Mean beta=[ ");
		for (j=0;j<numExperiments;j++) {
			System.out.print(meanBeta[j]);
			System.out.print(" ");
		}
		System.out.println("]");
	}
}
