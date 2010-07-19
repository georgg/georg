package edu.mit.csail.psrg.georg.IGMM;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Vector;

import cern.jet.random.Gamma;
import cern.jet.random.Uniform;
import cern.jet.random.Beta;
import cern.jet.random.Normal;
import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.HierarchicalCluster.HierarchicalCluster;
import edu.mit.csail.psrg.georg.SampleUtil.*;
import edu.mit.csail.psrg.georg.StatUtil.StatUtil;

public class GaussianDP implements Serializable {
	int numGenes = 0;
	int numExperiments = 0;
	int numClusters = 0;
	// allocated size of the cluster array
	int allocClusters = 5;
	// amount to grow cluster allocation by
	int growClusterAlloc = 10;
	
	// expression values for each gene
	double[][] expression = null;
	// array of genes with NaN values (there will be duplicates)
	int[] NaNGenes = null;
	// array of experiments with NaN values (corresponding to NaNGenes
	int[] NaNExperiments = null;
	
	// the concentration parameter for clustering
	double alpha = 1.0;
	// parameters for the gamma prior on alpha
	double alpha_a = 0.5;	// 10e-10;
	double alpha_b = 1.0;	// 10e-10;
	
	// cluster means (size = clusters X experiments)
	double[][] mu = null;
	// used for sampling mu values
	double[] tempMu = null;
	// mean hyperparameters on mu (size = # experiments)
	double[] lambda = null;
	// varaince hyperparameter on mu
	double r = 0.0;
	// cluster variances (size = # clusters)
	double[] s = null;
	// hyperparameters on s
	double w = 0.0;
	double beta = 0.0;
	// used for sampling beta
//	Vector<Double> x_init = new Vector<Double>();
//	BetaPosteriorSampleFn betaFn = new BetaPosteriorSampleFn();
//	BetaPosteriorLogDensityFn betaFn = new BetaPosteriorLogDensityFn();
	BetaPosteriorSampleFn2 betaFn = new BetaPosteriorSampleFn2();
	AdaptiveRejectionSampleCR ARS = new AdaptiveRejectionSampleCR();
	
	// total means for all expression data (size = # experiments)
	double[] muX = null;
	// total standard deviation for all expression data
	double sX = 0.0;
	
	// number of genes assigned to each cluster
	int[] numGenesInCluster = null;
	// cluster assignment for each gene
	int[] clusterAssign = null;
	// internal variable used for computing conditional likelihood
	double[] condLike = null;
	// used to map clusters to the genes they contain
	HashSet[] clusterMap;
	
	// used to store the probability of two genes co-occurring in a cluster
	double[][] pairProb;
	int numGoodIters = 0;
	boolean capturePairProb = true;
	
	// variables to control snap-shots (saving of the whole DP object
	// every N # of samples)
	int snapShotInterval = 1000;
	String snapShotFileName = null;
	
	String[] geneNames = null;
	
	// array of genes pairs to be tracked
	int[][] pairedGenes = null;
	// array for counting co-clustering of pairs
	int[][] pairCounts = null;
	boolean capturePairCounts = false;
	
	public GaussianDP(MicroArrayData myData,int nc) {
		expression = myData.values;
		numGenes = myData.numRows;
		geneNames = myData.geneNames;
		numExperiments = myData.numCols;
		clusterAssign = new int[numGenes];
		muX = myData.meanCols();
		tempMu = new double[numExperiments];
		lambda = new double[numExperiments];
		numClusters = nc;
		
		flagNaNs();
		initAllocClusters();
		initParameters();
		initClusterAssign();
		resampleMissingValues();
		
		// initialize pairwise distances
		pairProb = new double[numGenes-1][];
		int nn = 0;
		int kk = 0;
		int jj = 0;
		for (jj=0;jj<numGenes-1;jj++) {
			nn = 0;
			for (kk=(jj+1);kk<numGenes;kk++) {
				nn++;
			}
			pairProb[jj] = new double[nn];
		} 
	}
	
	void initParameters() {
		totalS();
		
		// sample the variance hyperparamter for mu
		r = Gamma.staticNextDouble(0.5,Math.pow(sX,2.0)*0.5);
		
		// set hyperparamters for s
	//	beta = ((double) numGenes)/((double) numClusters);
	//	beta = 1.0;
		beta = (double) numExperiments;
		w = Math.pow(sX,2.0);
		// set up initial vector for adaptive rejection sampling
/*		x_init.add(-3.5066);
		x_init.add(-2.9957);
		x_init.add(-2.3026);
		x_init.add(-1.3026);
		x_init.add(-0.3026);
		x_init.add(0.6974);
		x_init.add(1.6974);
		x_init.add(2.6974);
		x_init.add(3.6974);
		x_init.add(4.3175);
		x_init.add(5.0106);
		x_init.add(5.2983);
		x_init.add(5.7038);
		x_init.add(6.2146);
		x_init.add(6.6201);
		x_init.add(6.9078); */
		
		int j = 0;
		for (j=0;j<numExperiments;j++) {
			lambda[j] = Normal.staticNextDouble(muX[j],sX);
		}
		
	}
	
	public void enableSnapShots(int i,String f) {
		snapShotInterval = i;
		snapShotFileName = f;
	}
	
	public void disableSnapShots() {
		snapShotFileName = null;
	}
	
	void totalS() {
		double x = 0.0;
		double tv = 0.0;
		double xt = 0.0;
		int i = 0;
		int j = 0;
		for(i=0;i<numGenes;i++) {
			for(j=0;j<numExperiments;j++) {
				if (!Double.isNaN(expression[i][j])) {
					xt = Math.pow(expression[i][j] - muX[j],2.0);
					x = x + xt;
					tv = tv + 1.0;
				}
			}
		}
		sX = Math.sqrt(x/(tv-1.0));
	}
	
	void initAllocClusters() {
		int totalAllocate = allocClusters;
		if (allocClusters < numClusters+2) {
			totalAllocate = numClusters+growClusterAlloc;
		}
		
		numGenesInCluster = new int[totalAllocate];
		mu = new double[totalAllocate][numExperiments];
		s = new double[totalAllocate];
		condLike = new double[totalAllocate];
		
		clusterMap = new HashSet[totalAllocate];
		int i = 0;
		for (i=0;i<totalAllocate;i++) {
			clusterMap[i] = new HashSet();
		}
				
		allocClusters = totalAllocate;
	}
	
	void initClusterAssign() {
		int i = 0;
		int j = 0;
		int cc = 0;
		double[][] tn = new double[numClusters][numExperiments];
		double[] tns = new double[numClusters];
		
		for (i=0;i<numGenes;i++) {
			cc = Uniform.staticNextIntFromTo(0,numClusters-1);
			numGenesInCluster[cc]++;
			clusterAssign[i] = cc;
			clusterMap[cc].add(new Integer(i));
			for (j=0;j<numExperiments;j++) {
				if (!Double.isNaN(expression[i][j])) {
					tn[cc][j] = tn[cc][j] + 1.0;
					mu[cc][j] = mu[cc][j] + expression[i][j];
				}
			}
		}
		
		for (cc=0;cc<numClusters;cc++) {
			for (j=0;j<numExperiments;j++) {
				mu[cc][j] = mu[cc][j]/tn[cc][j];
			}
		}
		
		for (i=0;i<numGenes;i++) {
			cc = clusterAssign[i];
			for (j=0;j<numExperiments;j++) {
				if (!Double.isNaN(expression[i][j])) {
					s[cc] = s[cc] + Math.pow(expression[i][j]-mu[cc][j],2.0);
					tns[cc] = tns[cc] + 1.0;
				}
			}
		}
		
		for (cc=0;cc<numClusters;cc++) {
			s[cc] = Math.sqrt(s[cc]/tns[cc]);
		}
	}
	
	void flagNaNs() {
		int i = 0;
		int j = 0;
		ArrayList<Integer> ng = new ArrayList<Integer>();
		ArrayList<Integer> ne = new ArrayList<Integer>();
		
		for(i=0;i<numGenes;i++) {
			for(j=0;j<numExperiments;j++) {
			if (Double.isNaN(expression[i][j])) {
				ng.add(i);
				ne.add(j);
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
	
	void addCluster() {
		if (numGenesInCluster.length < numClusters+2) {
			int[] numGenesInCluster2 = new int[allocClusters+growClusterAlloc];
			System.arraycopy(numGenesInCluster, 0, numGenesInCluster2, 0, allocClusters);
			numGenesInCluster = numGenesInCluster2;
			
			double[][] mu2 = new double[allocClusters+growClusterAlloc][numExperiments];
			System.arraycopy(mu, 0, mu2, 0, allocClusters);
			mu = mu2;
			
			double[] s2 = new double[allocClusters+growClusterAlloc];
			System.arraycopy(s, 0, s2, 0, allocClusters);
			s = s2;
			
			double[] condLike2 = new double[allocClusters+growClusterAlloc];
			System.arraycopy(condLike,0,condLike2,0,allocClusters);
			condLike = condLike2;
			
			int jj = 0;
			HashSet[] clusterMap2 = new HashSet[allocClusters+growClusterAlloc];
			System.arraycopy(clusterMap,0,clusterMap2,0,allocClusters);
			for (jj=allocClusters;jj<allocClusters+growClusterAlloc;jj++) {
				clusterMap2[jj] = new HashSet();
			}
			clusterMap = clusterMap2;
			
			allocClusters = numGenesInCluster.length;
		}
		
		numGenesInCluster[numClusters+1] = 0;
		s[numClusters+1] = 0.0;
		clusterMap[numClusters+1].clear();
		
		int k = 0;
		for (k=0;k<numExperiments;k++) {
			mu[numClusters+1][k] = 0.0;
		}
		
		numClusters++;
	}
	
	void deleteCluster(int c) {
		// the cluster to be deleted is zero-indexed, as is standard in Java
		// this method will copy numClusters+1 entries
		
		System.arraycopy(numGenesInCluster,c+1,numGenesInCluster,c,numClusters-c);
		System.arraycopy(s,c+1,s,c,numClusters-c);
		
		int k = 0;
		
		double[] mc = mu[c];
		for(k=c;k<numClusters;k++) {
			mu[k] = mu[k+1];
		}
		mu[numClusters] = mc;
		
		HashSet cm = clusterMap[c];
		for(k=c;k<numClusters;k++) {
			clusterMap[k] = clusterMap[k+1];
		}
		clusterMap[numClusters] = cm;
		
		for (k=0;k<clusterAssign.length;k++) {
			if (clusterAssign[k] > c) {
				clusterAssign[k] = clusterAssign[k] - 1;
			}
		}
		
		numClusters--;
	}
	
	void sampleGeneClusterAssigns() {
		int cc = 0;
		int oldCC = 0;
		int i = 0;
		int j = 0;
		double denom = ((double) numGenes - 1.0 + alpha);
		double tempS = 0.0;
		double tmp = 0.0;
		double minL = 0.0;
		
		for (i=0;i<numGenes;i++) {
			oldCC = clusterAssign[i];
			for (cc=0;cc<numClusters;cc++) {
				condLike[cc] = StatUtil.logNormPDF(mu[cc],s[cc],expression[i]);
			}
			// sample new mu and s
			for (j=0;j<numExperiments;j++) {
				tempMu[j] = Normal.staticNextDouble(lambda[j],1.0/r);
				tempS = Gamma.staticNextDouble(beta*0.5,beta*w*0.5);
				tempS = Math.sqrt(1.0/tempS);
			}
			condLike[numClusters] = StatUtil.logNormPDF(tempMu,tempS,expression[i]);
			
			minL = condLike[0];
			for (cc=1;cc<numClusters;cc++) {
				if (condLike[cc] < minL) {
					minL = condLike[cc];	
				}
			}
			
			for (cc=0;cc<numClusters;cc++) {
			//	condLike[cc] = condLike[cc] - minL;
				// we're considering a cluster that the gene is not already assigned to
				if (cc != oldCC) {
					condLike[cc] = Math.exp(condLike[cc])*((double) numGenesInCluster[cc])/denom;
				} else {
					if (numGenesInCluster[cc] > 1) {
						tmp = ((double) numGenesInCluster[cc]) - 1.0;
					} else {
						tmp = alpha;
					}
					condLike[cc] = Math.exp(condLike[cc])*tmp/denom;
				}
			}
		//	condLike[numClusters] = condLike[numClusters] - minL;
			condLike[numClusters] = Math.exp(condLike[numClusters])*alpha/denom;
		
			// now choose the new cluster
			cc = StatUtil.multinomial_rnd(condLike,numClusters+1);
		
			if (cc != oldCC) {
				numGenesInCluster[oldCC]--;
				numGenesInCluster[cc]++;
				clusterAssign[i] = cc;
				clusterMap[oldCC].remove(new Integer(i));
				clusterMap[cc].add(new Integer(i));
			/*	if (numGenesInCluster[oldCC]>=1) {
					resampleMuSSingle(oldCC);
				} 
				if (numGenesInCluster[cc]>=1) {
					resampleMuSSingle(cc);
				} */
			}
			
			if (cc==numClusters) {
				s[numClusters] = tempS;
				System.arraycopy(tempMu,0,mu[cc],0,numExperiments);
				addCluster();
			}
		}
	}
	
	void resampleMissingValues() {
		// fill-in initial guesses for missing values
		int gene = 0;
		int expt = 0;
		int i = 0;
		int cc = 0;
		if (NaNGenes != null) {
			for (i=0;i<NaNGenes.length;i++) {
				gene = NaNGenes[i];
				expt = NaNExperiments[i];
				cc = clusterAssign[gene];
				expression[gene][expt] = Normal.staticNextDouble(mu[cc][expt],s[cc]);
			}
		}
	}
	
	void resampleMuS() {
		int cc = 0;
		for (cc=0;cc<numClusters;cc++) {
			resampleMuSSingle(cc);
		}
	}
	
	// sample from the posterior for mu and s for cluster cc
	void resampleMuSSingle(int cc) {
		Object[] genes = null;
		int gene = 0;
		int gg = 0;
		int j = 0;
		double nj = 0.0;
		double tempS = 0.0;
		double mujs = 0.0;
		
		genes = clusterMap[cc].toArray();
		nj = (double) numGenesInCluster[cc];
		// sample posterior mu
		if (numGenesInCluster[cc] == 1) {
			gene = (Integer) genes[0];
			System.arraycopy(expression[gene],0,tempMu,0,numExperiments);
		} else {
			for(j=0;j<numExperiments;j++) {
				tempMu[j] = 0.0;
			}
			for (gg=0;gg<genes.length;gg++) {
				gene = (Integer) genes[gg];
				for(j=0;j<numExperiments;j++) {
					tempMu[j] = tempMu[j] + expression[gene][j];
				}
			}
			for(j=0;j<numExperiments;j++) {
				tempMu[j] = tempMu[j]/((double) genes.length);
			}
		}
		mujs = ((Math.pow(s[cc],2.0)/nj)*(1.0/r))/(1.0/r + (Math.pow(s[cc],2.0)/nj));
		mujs = Math.sqrt(mujs);
		for (j=0;j<numExperiments;j++) {
			tempMu[j] = (tempMu[j]/r + (Math.pow(s[cc],2.0)/nj)*lambda[j])/((1.0/r) + (Math.pow(s[cc],2.0)/nj));
			mu[cc][j] = Normal.staticNextDouble(tempMu[j],mujs);
		}
		
		// sample posterior variance
		if (numGenesInCluster[cc] == 1) {
			s[cc] = Gamma.staticNextDouble(0.5*beta,beta*w*0.5);
			s[cc] = Math.sqrt(1.0/s[cc]);
		} else {
			tempS = 0.0;
			for (gg=0;gg<genes.length;gg++) {
				gene = (Integer) genes[gg];
				for(j=0;j<numExperiments;j++) {
					tempS = tempS + Math.pow(expression[gene][j] - mu[cc][j],2.0);
				}
			}
			s[cc] = Gamma.staticNextDouble((((double) numExperiments)*nj+beta)/2.0,0.5*(beta*w + tempS));
			s[cc] = Math.sqrt(1.0/s[cc]);
		}
	}
	
	// resample mu and s hyperparameters (lambda, r, beta, and w)
	void resampleHyperparameters() {
		double lambda_s = 0.0;
		double Q = (double) numClusters;
		int j = 0;
		int cc = 0;
		
		// sample lambda
		lambda_s = (Math.pow(sX,2.0)*(1.0/(r*Q)))/(Math.pow(sX,2.0) + 1.0/(r*Q));
		lambda_s = Math.sqrt(lambda_s);
		for (j=0;j<numExperiments;j++) {
			tempMu[j] = 0.0;
			for (cc=0;cc<numClusters;cc++) {
				tempMu[j] = tempMu[j] + mu[cc][j]/Q;
			}
		}
		for (j=0;j<numExperiments;j++) {
			tempMu[j] = (Math.pow(sX,2.0)*tempMu[j] + muX[j]/(r*Q))/(Math.pow(sX,2.0) + 1.0/(r*Q));
			lambda[j] = Normal.staticNextDouble(tempMu[j],lambda_s);
		}
		
		// sample r
		double lr = 0.0;
		for (cc=0;cc<numClusters;cc++) {
			for (j=0;j<numExperiments;j++) {
				lr = lr + Math.pow(mu[cc][j]-lambda[j],2.0);
			}
		}
		r = Gamma.staticNextDouble((((double) numExperiments)*Q + 1.0)/2.0,(lr + Math.pow(sX,2.0))*0.5);
		
		// sample beta
	//	resampleBeta();
		double beta_next = 0.0;
		betaFn.init(s,w,numClusters);
	//	AdaptiveRejectionSampler sampler = new AdaptiveRejectionSampler(betaFn,x_init,-1.0/0.0,1.0/0.0,1);
	//	AdaptiveRejectionSampler sampler = new AdaptiveRejectionSampler(betaFn,Math.log(beta),0.1,-1.0/0.0,1.0/0.0,1);
	//	beta_next = sampler.samples[0];
		beta_next = ARS.ars(Math.log(beta),0.1,betaFn);
		if (!Double.isInfinite(beta_next) && !Double.isNaN(beta_next)) {
			beta = Math.exp(beta_next);
		}
		
		// sample w
		double ss = 0.0;
		for (j=0;j<numClusters;j++) {
			ss = ss + 1.0/Math.pow(s[j],2.0);
		}
		w = Gamma.staticNextDouble((Q*beta + 1.0)*0.5,(beta*ss +Math.pow(sX,-2.0))*0.5);
	}
	
/*	public void resampleBeta() {
		betaFn.init(s,w,numClusters);
		
		double[] x = new double[1];
		x[0] = beta;
		double[] wp = new double[1];
		wp[0] = 5.0;
		double[][] S = new double[1][1];
		SliceSampler.sample(1,100,betaFn,x,wp,S);
		if (!Double.isInfinite(S[0][0])) {
			beta = S[0][0];
			System.out.println("alpha_beta=" + S[0][0]);
		}
	} */
	
	// sample concentration parameter using aux. variable sampling scheme
	// from Escobar and West (1995)
	void resampleAlpha(int numIters) {
		int i = 0;
		double eta = 0.0;
		int zz = 0;
		double[] w = new double[2];
		if (numClusters <= 1) {
			return;
		}
		
		for (i=0;i<numIters;i++) {
			eta = Beta.staticNextDouble(alpha+1.0,numGenes); 
			if (eta > 0) {
				w[0] = (alpha_a + (double) (numClusters - 1))/(alpha_b - Math.log(eta));
				w[1] = (double) numGenes;
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
	
	
	void iterate(int numIters,int numBurnin,int numAlphaIters) {
		int i = 0;
		int nc = 0;
		int cc = 0;
		String snapShotOutName;
		
		for (i=0;i<numIters;i++) {
			sampleGeneClusterAssigns();
			
			nc = numClusters;
  			for (cc=nc-1;cc>=0;cc--) {
  				if (numGenesInCluster[cc] == 0) {
  					deleteCluster(cc);
  				}
  			}
			
  			resampleMuS();
  			resampleHyperparameters();
  			resampleAlpha(numAlphaIters);
  			resampleMissingValues();
			
			if (i >= numBurnin) {
				numGoodIters++;
				incPairProbs();
			}
			
			if (snapShotFileName != null) {
  				double t = ((double) i+1)/((double) snapShotInterval);
  				if (((double) Math.round(t)) == t) {
  					snapShotOutName = snapShotFileName + "_" + (new Integer(i)).toString() + ".persist";
  					System.out.println("Snap shot at iteration" + i);
  					persistSelf(snapShotOutName);
  				}
  			}
			
			System.out.println(i + " " + numClusters + " alpha=" + alpha + " beta=" + beta + " w=" + w);
		}
	//	normPairProbs(numIters-numBurnin);
	}
	
	public void writeClusterStatsToFile(String fname) throws IOException {
		int cc = 0;
		int j = 0;
		FileWriter outFile = new FileWriter(fname);
		
		for (cc=0;cc<numClusters;cc++) {
			outFile.write(new Double(s[cc]).toString());
			for (j=0;j<mu[cc].length;j++) {
				outFile.write("\t");
				outFile.write(new Double(mu[cc][j]).toString());
			}
			outFile.write("\n");
		}
		
		outFile.close();
	}
	
	public void writeClustersToFile(String fname) throws IOException {
		int cc = 0;
		int gg = 0;
		String s = "";
		FileWriter outFile = new FileWriter(fname);
		
		for (gg=0;gg<numGenes;gg++) {
			outFile.write(geneNames[gg] + "\t");
			outFile.write((new Integer(clusterAssign[gg])).toString());
			outFile.write("\n");
		}
		
		outFile.close();
	}

//used to increment the pair-wise adjacency matrix
	void incPairProbs() {
		int cc = 0;
		int ii = 0;
		int jj = 0;
		int ng = 0;
		int ID1 = 0;
		int ID2 = 0;
		Object[] genes = null;
		
		for (cc=0;cc<numClusters;cc++) {
			genes = clusterMap[cc].toArray();
			ng = genes.length;
			for (ii=0;ii<ng-1;ii++) {
				ID1 = (Integer) genes[ii];
				for (jj=(ii+1);jj<ng;jj++) {
					ID2 = (Integer) genes[jj];
					if (ID2 > ID1) {
						pairProb[ID1][ID2-ID1-1]++;
					} else {
						pairProb[ID2][ID1-ID2-1]++;
					}
				}
			}
		}
	}
	
	void normPairProbs(int goodIters) {
		double norm = (double) goodIters;
		int ii = 0;
		int jj = 0;
		int ng1 = 0;
		int ng2 = 0;
		
		ng1 = pairProb.length;
		
		for (ii=0;ii<ng1;ii++) {
			ng2 = pairProb[ii].length;
			for (jj=0;jj<ng2;jj++) {
				// normalize and convert to a distance
				pairProb[ii][jj] = 1.0 - pairProb[ii][jj]/norm;
			}
		}	
	}
	
//  this will find a concensus model from the pair-wise probabilities
  	void concensusModel() {
  		// perform hierarhical clustering, using 1.0-pair-wise probabilities as
  		// distances, cutting clusters with distances greater than 0.50
  		HierarchicalCluster hier = new HierarchicalCluster();
  		hier.cluster(pairProb,0.75);
  		
  		numClusters = hier.clusters.size();
  		
  		int cc = 0;
  		int gg = 0;
  		int gene = 0;
  		
  		for (cc=0;cc<numClusters;cc++) {
  			for (gg=0;gg<hier.clusters.get(cc).size();gg++) {
  				gene = hier.clusters.get(cc).get(gg);
  				clusterAssign[gene] = cc;
  			}
  		}
  		
  		numGenesInCluster = new int[numClusters];
  		mu = new double[numClusters][numExperiments];
  		s = new double[numClusters];
  		
  		int j = 0;
  		for (cc=0;cc<numClusters;cc++) {
  			s[cc] = 0;
  			numGenesInCluster[cc] = 0;
  			for (j=0;j<mu[cc].length;j++) {
  				mu[cc][j] = 0;
  			}
  		}
  		
  		for (gg=0;gg<numGenes;gg++) {
  			cc = clusterAssign[gg];
  			numGenesInCluster[cc]++;
  			for (j=0;j<mu[cc].length;j++) {
  				mu[cc][j] += expression[gg][j];
  			}
  		}
  		
  		for (cc=0;cc<numClusters;cc++) {
  			for (j=0;j<mu[cc].length;j++) {
  				mu[cc][j] = mu[cc][j]/((double) numGenesInCluster[cc]);
  			}
  		}
  		
  		for (gg=0;gg<numGenes;gg++) {
  			cc = clusterAssign[gg];
  			numGenesInCluster[cc]++;
  			for (j=0;j<mu[cc].length;j++) {
  				s[cc] += Math.pow(expression[gg][j]-mu[cc][j],2.0);
  			}
  		}
  		
  		for (cc=0;cc<numClusters;cc++) {
  			s[cc] = Math.sqrt(s[cc]/((double) numGenesInCluster[cc]*mu[cc].length));
  		}
  		
  	}
  	
  	void persistSelf(String fName) {
  		FileOutputStream fos = null;
  		ObjectOutputStream out = null;
  		try {
  			fos = new FileOutputStream(fName);
  			out = new ObjectOutputStream(fos);
  			out.writeObject(this);
  			out.close();
  			}
  		catch(IOException e) {
  			System.out.println(e);
  		}
  	}
  	
  	public static GaussianDP restoreFromFile(String fName) {
  		GaussianDP myDP = null;
  		FileInputStream fis = null;
  		ObjectInputStream in = null;
  		try {
  			fis = new FileInputStream(fName);
  		    in = new ObjectInputStream(fis);
  		    myDP = (GaussianDP) in.readObject();
  		    in.close();
  		} catch(IOException e) {
  			System.out.println(e);
  		} catch(ClassNotFoundException e) {
  			System.out.println(e);
  		}
  		return myDP;
  	}
}