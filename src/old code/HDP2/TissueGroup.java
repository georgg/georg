package edu.mit.csail.psrg.georg.HDP2;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import cern.jet.random.Beta;
import cern.jet.random.Gamma;
import edu.mit.csail.psrg.georg.HierarchicalCluster.HierarchicalCluster;
import edu.mit.csail.psrg.georg.StatUtil.*;

import java.io.Serializable;

// this class implements DP groupings (e.g., the parent of a DP is uncertain)
// a "controlling" DP is specified, which generates DP's under it
// so, when a new cluster is needed, the "controlling" DP's prior is used to sample
// a new DP
public class TissueGroup implements Serializable {
	/**
	 * 
	 */
	private static final long serialVersionUID = -9032855203386904422L;
	public HDP myModel = null;
	public ParentDP controllingProcess = null;
	// concentration parameters for clustering
	public double alpha = 1.0;
	public double alpha_a = 1.0;
	public double alpha_b = 0.25;
	
	// all the DPs in the group
	public ArrayList<TissueDP> groupableTissues = new ArrayList<TissueDP>();
	
	public int numClusters = 0;
	public int numGroupable = 0;
	public int numGoodIter = 0;
	
	// used for storing pairwise probabilities
	// for hierarchical clustering
	public double[][] pairProb = null;
	
	public TissueGroup(HDP mm, ParentDP c) {
		myModel = mm;
		controllingProcess = c;
		int i = 0;
		int j = 0;
		ParentDP pdp = null;
		TissueDP dp = null;
		
		int dpID = 0;
		for (i=0;i<controllingProcess.numChildren();i++) {
			pdp = (ParentDP) controllingProcess.children.get(i);
			for (j=0;j<pdp.numChildren();j++) {
				dp = (TissueDP) pdp.children.get(j);
				groupableTissues.add(dp);
			}
		}
		numGroupable = groupableTissues.size();
		
		initPairProb();
	}
	
	public void initPairProb() {
		//  initialize pairwise distances
		pairProb = new double[numGroupable-1][];
		int nn = 0;
		int jj = 0;
		int kk = 0;
		for (jj=0;jj<numGroupable-1;jj++) {
			nn = 0;
			for (kk=(jj+1);kk<numGroupable;kk++) {
		  		nn++;
		  	}
		  	pairProb[jj] = new double[nn];
		}
	}
	
	public void sampleDPAssignments(int numTopics) {
		int i = 0;
		int cc = 0;
		int tt = 0;
		numClusters = controllingProcess.numChildren();
		TissueDP dp = null;
		ParentDP parent = null;
		ParentDP oldParent = null;
		double[] condLike = new double[numClusters+2];
		double[] condLikeBeta = new double[myModel.allocTopics];
		double[] newBeta = new double[myModel.allocTopics];
		double[] newTheta = new double[myModel.allocTopics];
		double tmp = 0.0;
		double minCondLike = 0;
		
		for (i=0;i<numGroupable;i++) {
			dp = groupableTissues.get(i);
			if (dp.state == DirichletProcess.ACTIVE) {
				for (cc=0;cc<numClusters;cc++) {
					parent = (ParentDP) controllingProcess.children.get(cc);	
				//	condLike[cc] = dataLikeDP(numTopics,dp.clusterNumData, parent.beta, dp.alpha, dp.upDownTopicChoice, parent.theta);
					condLike[cc] = dataLikeDP(numTopics,parent.beta,parent.theta,dp);
					}
				}
				// sample new beta and theta from the prior
				for (tt=0;tt<numTopics+1;tt++) {
					condLikeBeta[tt] = parent.alpha*controllingProcess.beta[tt];
				//	condLikeBeta[tt] = controllingProcess.beta[tt];
					newTheta[tt] = parent.sampleSingleTheta(tt);
				}
				StatUtil.dirichlet_rnd(newBeta,condLikeBeta,numTopics+1);
			//	condLike[numClusters] = dataLikeDP(numTopics,dp.clusterNumData, newBeta, dp.alpha, dp.upDownTopicChoice, newTheta);
				condLike[numClusters] = dataLikeDP(numTopics,newBeta,newTheta,dp);
				
				minCondLike = condLike[0];
				for (cc=1;cc<numClusters+1;cc++) {
					if (condLike[cc] > minCondLike) {
						minCondLike = condLike[cc];
					}
				}
				
				for (cc=0;cc<numClusters;cc++) {
					condLike[cc] = condLike[cc] - minCondLike;
					parent = (ParentDP) controllingProcess.children.get(cc);
					if (dp.parent == parent) {
						tmp = ((ParentDP) parent).numChildren() - 1.0;
						if (tmp == 0.0) { tmp = alpha; }
					} else { tmp = ((ParentDP) parent).numChildren(); }
					condLike[cc] = tmp*Math.exp(condLike[cc]);	
				}
				condLike[numClusters] = condLike[numClusters] - minCondLike;
				condLike[numClusters] = alpha*Math.exp(condLike[numClusters]);
				
				for (cc=0;cc<numClusters+1;cc++) {
					if (Double.isNaN(condLike[cc]) | Double.isInfinite(-condLike[cc])) {
						condLike[cc] = 0.0;
						System.out.println("Overflow in group estimation = " + dp.label);
					} else {
						if (Double.isInfinite(condLike[cc])) {
							condLike[cc] = 1.0;
							System.out.println("Overflow in group estimation = " + dp.label);
						}
					}
				}
				
				oldParent = dp.parent;
				cc = StatUtil.multinomial_rnd(condLike,numClusters+1);
				
				if (cc == numClusters) {
					numClusters++;
					condLike = new double[numClusters+2];
					parent = newDirichletProcess(newBeta,newTheta,controllingProcess.children.get(0).concParam,((ParentDP) controllingProcess.children.get(0)).myUpDownConcParam);
				} else {
					parent = (ParentDP) controllingProcess.children.get(cc);
				}
				
				if (oldParent != parent) {			
					dp.parent = parent;
					oldParent.removeChild(dp);
					parent.addChild(dp);
					
					for (tt=0;tt<numTopics;tt++) {
						oldParent.clusterNumData[tt] = oldParent.clusterNumData[tt] - dp.clusterNumTables[tt];
						parent.clusterNumData[tt] = parent.clusterNumData[tt] + dp.clusterNumTables[tt];
						if (dp.upDownTopicChoice[tt]) {
							oldParent.numUpData[tt]--;
							parent.numUpData[tt]++;
						} else {
							oldParent.numDownData[tt]--;
							parent.numDownData[tt]++;
						}
					}
					parent.resampleBeta(numTopics,newBeta);
					parent.resampleTheta(numTopics);
					if (oldParent.numChildren() > 0) {
						oldParent.resampleBeta(numTopics,newBeta);
						oldParent.resampleTheta(numTopics);
					}
				}
		}
		removeEmptyDPs(numTopics);
	}
	
	// computes the data likelihood given that the data is generated by the prior
	// dp_alpha*parent_beta (child_beta is integrated out)
	public double dataLikeDP(int numTopics,int[] clusterNumData, double[] parent_beta, double dp_alpha, boolean[] upDownTopicChoice, double[] theta) {
		double lik = 0.0;
		int tt = 0;
		int j = 0;
		int nd = 0;
		
		for(tt=0;tt<numTopics;tt++) {
			nd = clusterNumData[tt];
			if (nd > 0) {
				for (j=1;j<=nd;j++) {
					lik = lik + Math.log(parent_beta[tt]*dp_alpha + ((double) j));
				}
			}
		
			if (nd > 0) {
				if (upDownTopicChoice[tt]) {
					lik = lik + Math.log(theta[tt]);
				} else {
					lik = lik + Math.log(1.0-theta[tt]);
				}
			}
		}
		
		return lik;
	}
	
//	 computes the data likelihood given that the data is generated by the prior
	// dp_alpha*parent_beta (child_beta is integrated out)
//	dataLikeDP(numTopics,dp.clusterNumData, newBeta, dp.alpha, dp.upDownTopicChoice, newTheta);
	
	public double dataLikeDP(int numTopics, double[] parent_beta, double[] theta, TissueDP dp) {
		double lik = 0.0;
		int tt = 0;
		int j = 0;
		int nd = 0;
		double upChoicePart = 0.0;
		double downChoicePart = 0.0;
		int[] clusterNumData = dp.clusterNumData;
		int[] upTopicCount = dp.upTopicCount;
		double dp_alpha = dp.alpha;
		double matchProb = myModel.upProb;
		
		for(tt=0;tt<numTopics;tt++) {
			nd = dp.clusterNumData[tt];
			if (nd > 0) {
				for (j=1;j<=nd;j++) {
					lik = lik + Math.log(parent_beta[tt]*dp_alpha + ((double) j));
				}
				
				upChoicePart = ((double) upTopicCount[tt])*Math.log(matchProb) + ((double) clusterNumData[tt] - upTopicCount[tt])*Math.log(1.0-matchProb);
				downChoicePart = ((double) upTopicCount[tt])*Math.log(1.0-matchProb) + ((double) clusterNumData[tt] - upTopicCount[tt])*Math.log(matchProb);
				lik = lik + Math.log(theta[tt]*Math.exp(upChoicePart) + (1-theta[tt]*Math.exp(downChoicePart)));
			}
			
		/*	if (nd > 0) {
				if (dp.upDownTopicChoice[tt]) {
					lik = lik + Math.log(theta[tt]);
				} else {
					lik = lik + Math.log(1.0-theta[tt]);
				}
			} */
		}
		
		return lik;
	}
	
	public void removeEmptyDPs(int numTopics) {
		int nc = controllingProcess.numChildren();
		int cc = 0;
		int tt = 0;
		ParentDP dp = null;
		for (cc=nc-1;cc>=0;cc--) {
			dp = (ParentDP) controllingProcess.children.get(cc);
			if (dp.numChildren() == 0) {
				dp.concParam.removeDP(dp);
				dp.myUpDownConcParam.removeDP(dp);
				controllingProcess.removeChild(dp);
				// if this process was initially there, it will need to be removed
				for (tt=0;tt<numTopics;tt++) {
					controllingProcess.clusterNumData[tt] = controllingProcess.clusterNumData[tt] - dp.clusterNumTables[tt];
					((ParentDP) controllingProcess).numUpData[tt] = ((ParentDP) controllingProcess).numUpData[tt] - ((ParentDP) dp).numUpTables[tt];
					((ParentDP) controllingProcess).numDownData[tt] = ((ParentDP) controllingProcess).numDownData[tt] - ((ParentDP) dp).numDownTables[tt];
				}
			}
		}
	}
	
	public ParentDP newDirichletProcess(double[] newBeta,double[] newTheta,HDPConcParam concParam,UpDownConcParam UDconcParam) {
		ParentDP dp = new ParentDP(myModel,controllingProcess,"");
		double[] beta = new double[myModel.allocTopics];
		System.arraycopy(newBeta,0,beta,0,newBeta.length);
		double[] theta = new double[myModel.allocTopics];
		System.arraycopy(newTheta,0,theta,0,newTheta.length);
		
		int[] clusterNumData = new int[myModel.allocTopics];
	  	int[] clusterNumTables = new int[myModel.allocTopics];
	  	
	  	int[] numUpData = new int[myModel.allocTopics];
	  	int[] numDownData = new int[myModel.allocTopics];
	  	int[] numUpTables = new int[myModel.allocTopics];
	  	int[] numDownTables = new int[myModel.allocTopics];
	  	
	  	dp.beta = beta;
	  	dp.theta = theta;
	  	dp.clusterNumData = clusterNumData;
	  	dp.clusterNumTables = clusterNumTables;
	  	dp.numUpData = numUpData;
	  	dp.numDownData = numDownData;
	  	dp.numUpTables = numUpTables;
	  	dp.numDownTables = numDownTables;
  		
  		dp.alpha = concParam.alpha;
  		dp.concParam = concParam;
  		concParam.appendDP(dp);
  		
  		dp.alpha_theta = UDconcParam.alpha_theta;
  		dp.myUpDownConcParam = UDconcParam;
  		UDconcParam.appendDP(dp);
	  	
		return dp;
	}
	
	// sample concentration parameter using aux. variable sampling scheme
	// from Escobar and West (1995)
	public void resampleAlpha(int numIters) {
		int i = 0;
		double eta = 0.0;
		int zz = 0;
		double[] w = new double[2];
		numClusters = controllingProcess.numChildren();
		if (numClusters <= 1) {
			return;
		}
		
		int numPresent = 0;
		
		for (i=0;i<groupableTissues.size();i++) {
			if (groupableTissues.get(i).state != DirichletProcess.HELDOUT) {
				numPresent++;
			}
		}
		
		
		for (i=0;i<numIters;i++) {
			eta = Beta.staticNextDouble(alpha+1.0,(double) numPresent); 
			if (eta > 0) {
				w[0] = (alpha_a + (double) (numClusters - 1))/(alpha_b - Math.log(eta));
				w[1] = (double) numPresent;
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
	
	//	used to increment the pair-wise adjacency matrix
	public void incPairProbs() {
		int cc = 0;
		int ii = 0;
		int jj = 0;
		int nd = 0;
		int ID1 = 0;
		int ID2 = 0;
		ArrayList<DirichletProcess> dps = null;
		
		numClusters = controllingProcess.numChildren();
		
		for (cc=0;cc<numClusters;cc++) {
			dps = ((ParentDP) ((ParentDP) controllingProcess).children.get(cc)).children;
			nd = dps.size();
			
			for (ii=0;ii<nd-1;ii++) {
				ID1 = dps.get(ii).dpID;
				for (jj=(ii+1);jj<nd;jj++) {
					ID2 = dps.get(jj).dpID;
					if (ID2 > ID1) {
						pairProb[ID1][ID2-ID1-1]++;
					} else {
						pairProb[ID2][ID1-ID2-1]++;
					}
				}
			}
		}
		numGoodIter++;
	}
	
	public void normPairProbs() {
		double norm = (double) numGoodIter;
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
	
	// output a file with the actual probabilities _probs.txt
	// and a file with the descriptions _probs_desc.txt
	public void writePairwiseProbs(String fname) throws IOException {
		String f1 = fname + "_probs.txt";
		String f2 = fname + "_probs_desc.txt";
		FileWriter outFile = new FileWriter(f1);
  		String s = "";
  		int ii = 0;
  		int jj = 0;
  		int ng1 = 0;
  		int ng2 = 0;
  		
  		ng1 = pairProb.length;
		
		for (ii=0;ii<ng1;ii++) {
			ng2 = pairProb[ii].length;
			for (jj=0;jj<ng2;jj++) {
				// normalize and convert to a distance
				outFile.write((new Double(pairProb[ii][jj])).toString());
				outFile.write("\n");
			}
		}
		outFile.close();
		
		outFile = new FileWriter(f2);
		
		for (ii=0;ii<numGroupable;ii++) {
			outFile.write(groupableTissues.get(ii).getLabel());
			outFile.write("\n");
		}
		
		outFile.close();
	}
  	
  	// this will build a hierarchical DP based on the concensus model
	public void buildConcensusModel() {
  		HierarchicalCluster hier = new HierarchicalCluster();
  		hier.cluster(pairProb,0.10);
  	//	hier.cluster(pairProb,0.50);
  		
  		numClusters = hier.clusters.size();
  		
  		int cc = 0;
  		int gg = 0;
  		int itemNum = 0;
  		DirichletProcess dp = null;
  		ParentDP parent = null;
  		double[] newBeta = null;
  		double[] newTheta = null;
  		
  		// remove the groups - we'll add them back
  		HDPConcParam cp = controllingProcess.children.get(0).concParam;
  		UpDownConcParam udcp = ((ParentDP) controllingProcess.children.get(0)).myUpDownConcParam;
  		controllingProcess.children.clear();
  	  		
  		// reconstruct the groups
  		for (cc=0;cc<numClusters;cc++) {
  			newBeta = new double[controllingProcess.beta.length];
  			newTheta = new double[controllingProcess.theta.length];
  			parent = newDirichletProcess(newBeta,newTheta,cp,udcp);
  			for (gg=0;gg<hier.clusters.get(cc).size();gg++) {
  				itemNum = hier.clusters.get(cc).get(gg);
  				dp = groupableTissues.get(itemNum);
  				dp.parent = parent;
				parent.addChild(dp);
  			}
  		}
  	}
}
