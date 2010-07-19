package edu.mit.csail.psrg.georg.GeneProgram2;

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
		int j = 0;
		int k = 0;
		int cc = 0;
		int tt = 0;
		numClusters = controllingProcess.numChildren();
		TissueDP dp = null;
		ParentDP parent = null;
		ParentDP oldParent = null;
		double[] condLike = new double[numClusters+2];
		double[] condLikeBeta = new double[myModel.allocExpressionPrograms];
		double[] newBeta = new double[myModel.allocExpressionPrograms];
		double[] newTheta = null;
		double[][][] newThetaMod = null;
		double tmp = 0.0;
		double minCondLike = 0;
		
		if (myModel.useUpDown) {
			newTheta = new double[myModel.allocExpressionPrograms];
		}
		
		if (myModel.modifierLevels != null) {
			newThetaMod = new double[myModel.modifierLevels.length][][];
			for (i=0;i<myModel.modifierLevels.length;i++) {
				newThetaMod[i] = new double[myModel.modifierLevels[i]][myModel.allocExpressionPrograms];
			}
		}
		
		for (i=0;i<numGroupable;i++) {
			dp = groupableTissues.get(i);
			if (dp.state == DirichletProcess.ACTIVE) {
				for (cc=0;cc<numClusters;cc++) {
					parent = (ParentDP) controllingProcess.children.get(cc);	
					condLike[cc] = dataLikeDP(numTopics,parent.beta,parent.theta,parent.thetaMod,dp);
					}
				}
				// sample new beta and theta from the prior
				for (tt=0;tt<numTopics+1;tt++) {
					condLikeBeta[tt] = parent.alpha*controllingProcess.beta[tt];
					
					if (myModel.useUpDown) {
						newTheta[tt] = parent.sampleSingleTheta(tt);
					}
					
					if (myModel.modifierLevels != null) {
						parent.sampleSingleThetaMod(tt);
						for (k=0;k<myModel.modifierLevels.length;k++) {
							for (j=0;j<myModel.modifierLevels[k];j++) {
								newThetaMod[k][j][tt] = parent.mdt[k][j];
							}
						}
					}
				}
				StatUtil.dirichlet_rnd(newBeta,condLikeBeta,numTopics+1);
			//	condLike[numClusters] = dataLikeDP(numTopics,dp.clusterNumData, newBeta, dp.alpha, dp.upDownTopicChoice, newTheta);
				condLike[numClusters] = dataLikeDP(numTopics,newBeta,newTheta,newThetaMod,dp);
				
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
					parent = newDirichletProcess(newBeta,newTheta,newThetaMod,controllingProcess.children.get(0).concParam,((ParentDP) controllingProcess.children.get(0)).myUpDownConcParam,((ParentDP) controllingProcess.children.get(0)).myModifierConcParam);
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
						
						if (myModel.useUpDown) {
							if (dp.upDownExpressionProgramChoice[tt]) {
								oldParent.numUpData[tt]--;
								parent.numUpData[tt]++;
							} else {
								oldParent.numDownData[tt]--;
								parent.numDownData[tt]++;
							}
						}
						
						if (myModel.modifierLevels != null) {
							for (k=0;k<myModel.modifierLevels.length;k++) {
								oldParent.numModifierData[k][dp.modifierExpressionProgramChoice[k][tt]][tt]--;
								parent.numModifierData[k][dp.modifierExpressionProgramChoice[k][tt]][tt]++;
							}
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
	
	public double dataLikeDP(int numTopics, double[] parent_beta, double[] theta, double[][][] thetaMod, TissueDP dp) {
		double lik = 0.0;
		int tt = 0;
		int i = 0;
		int j = 0;
		int nd = 0;
		double upChoicePart = 0.0;
		double downChoicePart = 0.0;
		double modPart = 0.0;
		double accum = 0.0;
		int[] clusterNumData = dp.clusterNumData;
		double dp_alpha = dp.alpha;
		
		double matchProb = myModel.upProb;
		int[] upTopicCount = dp.upExpressionProgramCount;
		
		for(tt=0;tt<numTopics;tt++) {
			nd = dp.clusterNumData[tt];
			if (nd > 0) {
				for (j=1;j<=nd;j++) {
					lik = lik + Math.log(parent_beta[tt]*dp_alpha + ((double) j));
				}
				
				if (myModel.useUpDown) {
					upChoicePart = ((double) upTopicCount[tt])*Math.log(matchProb) + ((double) clusterNumData[tt] - upTopicCount[tt])*Math.log(1.0-matchProb);
					downChoicePart = ((double) upTopicCount[tt])*Math.log(1.0-matchProb) + ((double) clusterNumData[tt] - upTopicCount[tt])*Math.log(matchProb);
					lik = lik + Math.log(theta[tt]*Math.exp(upChoicePart) + (1-theta[tt]*Math.exp(downChoicePart)));
				}
				
				if (thetaMod != null) {
					for (i=0;i<thetaMod.length;i++) {
						accum = 0.0;
						
						for (j=0;j<thetaMod[i].length;j++) {
							// all the expression compatible with this choice
							modPart = ((double) dp.modifierExpressionProgramCount[i][j][tt])*Math.log(myModel.matchProbMod[i]);
							// all the expression incompatible with this choice
							modPart += ((double) (clusterNumData[tt]-dp.modifierExpressionProgramCount[i][j][tt]))*Math.log(myModel.unmatchProbMod[i]);
							accum += thetaMod[i][j][tt]*Math.exp(modPart);
						}
						lik += Math.log(accum);
				//		lik += Math.log(thetaMod[i][dp.modifierExpressionProgramChoice[i][tt]][tt]);
					} 
				}
			}
		}
		
		return lik;
	}
	
	public void removeEmptyDPs(int numTopics) {
		int nc = controllingProcess.numChildren();
		int cc = 0;
		int tt = 0;
		int i = 0;
		int j = 0;
		ParentDP dp = null;
		for (cc=nc-1;cc>=0;cc--) {
			dp = (ParentDP) controllingProcess.children.get(cc);
			if (dp.numChildren() == 0) {
				dp.concParam.removeDP(dp);
				
				if (dp.myUpDownConcParam != null) {
					dp.myUpDownConcParam.removeDP(dp);
				}
				
				if (dp.myModifierConcParam != null) {
					dp.myModifierConcParam.removeDP(dp);
				}
				controllingProcess.removeChild(dp);
				// if this process was initially there, it will need to be removed
				for (tt=0;tt<numTopics;tt++) {
					controllingProcess.clusterNumData[tt] = controllingProcess.clusterNumData[tt] - dp.clusterNumTables[tt];
				
					if (myModel.useUpDown) {
						((ParentDP) controllingProcess).numUpData[tt] = ((ParentDP) controllingProcess).numUpData[tt] - ((ParentDP) dp).numUpTables[tt];
						((ParentDP) controllingProcess).numDownData[tt] = ((ParentDP) controllingProcess).numDownData[tt] - ((ParentDP) dp).numDownTables[tt];
					}
					
					if (myModel.modifierLevels != null) {
						for (i=0;i<myModel.modifierLevels.length;i++) {
							for (j=0;j<myModel.modifierLevels[i];j++) {
								((ParentDP) controllingProcess).numModifierData[i][j][tt] = ((ParentDP) controllingProcess).numModifierData[i][j][tt] - ((ParentDP) dp).numModifierTables[i][j][tt];
							}
						}
					}
				}
			}
		}
	}
	
	public ParentDP newDirichletProcess(double[] newBeta,double[] newTheta,double[][][] newThetaMod,HDPConcParam concParam,UpDownConcParam UDconcParam,ModifierConcParam ModConcParam) {
		ParentDP dp = new ParentDP(myModel,controllingProcess,"");
		double[] beta = new double[myModel.allocExpressionPrograms];
		System.arraycopy(newBeta,0,beta,0,newBeta.length);
		int[] clusterNumData = new int[myModel.allocExpressionPrograms];
	  	int[] clusterNumTables = new int[myModel.allocExpressionPrograms];
		
	  	double[] theta = null;
	  	int[] numUpData = null;
	  	int[] numDownData = null;
	  	int[] numUpTables = null;
	  	int[] numDownTables = null;
	  	double[][] md = null;
		double[][] mdt = null;
		
	  	if (myModel.useUpDown) {
			theta = new double[myModel.allocExpressionPrograms];
			System.arraycopy(newTheta,0,theta,0,newTheta.length);
		  	numUpData = new int[myModel.allocExpressionPrograms];
		  	numDownData = new int[myModel.allocExpressionPrograms];
		  	numUpTables = new int[myModel.allocExpressionPrograms];
		  	numDownTables = new int[myModel.allocExpressionPrograms];
		  }
	  	
	  	double[][][] thetaMod = null;
	  	int[][][] numModifierData = null;
	  	int[][][] numModifierTables = null;
		int i = 0;
		int j = 0;
		if (myModel.modifierLevels != null) {
			thetaMod = new double[myModel.modifierLevels.length][][];
			numModifierData = new int[myModel.modifierLevels.length][][];
			numModifierTables = new int[myModel.modifierLevels.length][][];
			md = new double[myModel.modifierLevels.length][];
			mdt = new double[myModel.modifierLevels.length][];
			for (i=0;i<newThetaMod.length;i++) {
				thetaMod[i] = new double[myModel.modifierLevels[i]][myModel.allocExpressionPrograms];
				numModifierData[i] = new int[myModel.modifierLevels[i]][myModel.allocExpressionPrograms];
				numModifierTables[i] = new int[myModel.modifierLevels[i]][myModel.allocExpressionPrograms];
				md[i] = new double[myModel.modifierLevels[i]];
  				mdt[i] = new double[myModel.modifierLevels[i]];
				for (j=0;j<newThetaMod[i].length;j++) {
					System.arraycopy(newThetaMod[i][j],0,thetaMod[i][j],0,newThetaMod[i][j].length);
				}
			}
		}
	  	
	  	dp.beta = beta;
	  	dp.theta = theta;
	  	dp.thetaMod = thetaMod;
	  	dp.clusterNumData = clusterNumData;
	  	dp.clusterNumTables = clusterNumTables;
	  	dp.numModifierData = numModifierData;
	  	dp.numModifierTables = numModifierTables;
	  	dp.numUpData = numUpData;
	  	dp.numDownData = numDownData;
	  	dp.numUpTables = numUpTables;
	  	dp.numDownTables = numDownTables;
  		
  		dp.alpha = concParam.alpha;
  		dp.concParam = concParam;
  		dp.md = md;
  		dp.mdt = mdt;
  		
  		concParam.appendDP(dp);
  		
  		if (UDconcParam != null) {
  			dp.alpha_theta = UDconcParam.alpha_theta;
  			dp.myUpDownConcParam = UDconcParam;
  			UDconcParam.appendDP(dp);
  		}
  		
  		if (ModConcParam != null) {
  			dp.alpha_thetaMod = new double[ModConcParam.alpha_theta.length];
  			for (i=0;i<ModConcParam.alpha_theta.length;i++) {
  				dp.alpha_thetaMod[i] = ModConcParam.alpha_theta[i];
  			}
  			dp.myModifierConcParam = ModConcParam;
  			ModConcParam.appendDP(dp);
  		}
  		
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
	public void buildConcensusModel(double minMergeProb) {
  		HierarchicalCluster hier = new HierarchicalCluster();
  		hier.cluster(pairProb,1.0-minMergeProb);
  		
  		numClusters = hier.clusters.size();
  		
  		int cc = 0;
  		int gg = 0;
  		int itemNum = 0;
  		DirichletProcess dp = null;
  		ParentDP parent = null;
  		double[] newBeta = null;
  		double[] newTheta = null;
  		double[][][] newThetaMod = null;
  		int i = 0;
  		int j = 0;
  		
  		// remove the groups - we'll add them back
  		HDPConcParam cp = controllingProcess.children.get(0).concParam;
  		UpDownConcParam udcp = ((ParentDP) controllingProcess.children.get(0)).myUpDownConcParam;
  		ModifierConcParam modcp = ((ParentDP) controllingProcess.children.get(0)).myModifierConcParam;
  		controllingProcess.children.clear();
  	  		
  		// reconstruct the groups
  		controllingProcess.myModel.allocExpressionPrograms = controllingProcess.beta.length;
  		for (cc=0;cc<numClusters;cc++) {
  			newBeta = new double[controllingProcess.beta.length];
  			
  			if (myModel.useUpDown) {
  				newTheta = new double[controllingProcess.theta.length];
  			}
  			
  			if (myModel.modifierLevels != null) {
  				newThetaMod = new double[controllingProcess.thetaMod.length][][];
  				for (i=0;i<myModel.modifierLevels.length;i++) {
  					newThetaMod[i] = new double[myModel.modifierLevels[i]][controllingProcess.thetaMod.length];
  				}
  			}
  			parent = newDirichletProcess(newBeta,newTheta,newThetaMod,cp,udcp,modcp);
  			for (gg=0;gg<hier.clusters.get(cc).size();gg++) {
  				itemNum = hier.clusters.get(cc).get(gg);
  				dp = groupableTissues.get(itemNum);
  				dp.parent = parent;
				parent.addChild(dp);
  			}
  		}
  	}
}
