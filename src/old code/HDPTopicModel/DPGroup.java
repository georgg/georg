package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import cern.jet.random.Beta;
import cern.jet.random.Gamma;
import edu.mit.csail.psrg.georg.HierarchicalCluster.HierarchicalCluster;
import edu.mit.csail.psrg.georg.StatUtil.*;

import java.io.Serializable;

// this class implements DP groupings (e.g., the parent of a DP is uncertain) or
// document groupings (e.g., the parent of a document is uncertain)
// in either case, a "controlling" DP is specified, which generates DP's under it
// so, when a new cluster is needed, the "controlling" DP's prior is used to sample
// a new DP
public class DPGroup implements Serializable {
	DirichletProcess controllingProcess = null;
	// concentration parameters for clustering
	double alpha = 1.0;
	double alpha_a = 1.0;
	double alpha_b = 0.25;
	
	// specifies whether to group at the document level or at the DP level
	boolean groupDocuments = true;
	
	// all the document-level DPs in the group (used if the groupDocuments option is not specified)
	ArrayList<DirichletProcess> groupableDPs = null;
	// all the groupable documents (used if the groupDocuments option is specified)
	ArrayList<Document> groupableDocs = null;
	
	int numClusters = 0;
	int numGroupable = 0;
	int numGoodIter = 0;
	
	// used for storing pairwise probabilities
	// for hierarchical clustering
	double[][] pairProb = null;
	
	public DPGroup(DirichletProcess c, boolean gd) {
		controllingProcess = c;
		groupDocuments = gd;
		int i = 0;
		int j = 0;
		DirichletProcess dp = null;
		
		if (groupDocuments) {
			Document doc = null;
			int docID = 0;
			groupableDocs = new ArrayList<Document>();
			for (i=0;i<controllingProcess.numChildren();i++) {
				dp = controllingProcess.children.get(i);
				if (dp.documents != null) {
					for (j=0;j<dp.documents.length;j++) {
						doc = dp.documents[j];
						groupableDocs.add(doc);
						doc.docID = docID;
						docID++;
					}
				}
				numGroupable = docID;
			}
		} else {
			int dpID = 0;
			groupableDPs = new ArrayList<DirichletProcess>();
			for (i=0;i<controllingProcess.numChildren();i++) {
				dp = controllingProcess.children.get(i);
				for (j=0;j<dp.numChildren();j++) {
					groupableDPs.add(dp.children.get(j));
					dp.children.get(j).dpID = dpID;
					dpID++;
				}
			}
			numGroupable = dpID;
		}
		
		initPairProb();
	}
	
	void initPairProb() {
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
		DirichletProcess dp = null;
		DirichletProcess parent = null;
		DirichletProcess oldParent = null;
		double[] condLike = new double[numClusters+1];
		double[] condLikeBeta = new double[numTopics+1];
		double[] newBeta = new double[numTopics+1];
		double tmp = 0.0;
		Document doc = null;
		double minCondLike = 0;
		
		for (i=0;i<numGroupable;i++) {
			if (groupDocuments) {
				doc = groupableDocs.get(i);
			} else{
				dp = groupableDPs.get(i);
			}
			if ((groupDocuments && doc.parent.state == DirichletProcess.ACTIVE) | (!groupDocuments && dp.state == DirichletProcess.ACTIVE)) {
				for (cc=0;cc<numClusters;cc++) {
					parent = controllingProcess.children.get(cc);	
					if (groupDocuments) {
						condLike[cc] = dataLikeDocument(numTopics,doc.clusterNumData, parent.beta);
					} else {
						condLike[cc] = dataLikeDP(numTopics,dp.clusterNumData, parent.beta, dp.alpha);
					//	condLike[cc] = Math.exp(dataLikeDP2(numTopics,dp.clusterNumData, parent.beta, dp.beta, dp.alpha));
					}
				}
				// sample new beta from the prior
				for (tt=0;tt<numTopics+1;tt++) {
					condLikeBeta[tt] = parent.alpha*controllingProcess.beta[tt];
				}
				StatUtil.dirichlet_rnd(newBeta,condLikeBeta,numTopics+1);
				// controllingProcess.resampleBeta(numTopics,condLikeBeta,newBeta);
				if (groupDocuments) {
					condLike[numClusters] = dataLikeDocument(numTopics,doc.clusterNumData, newBeta);
				} else {
					condLike[numClusters] = dataLikeDP(numTopics,dp.clusterNumData, newBeta, dp.alpha);
				//	condLike[numClusters] = alpha*Math.exp(dataLikeDP2(numTopics,dp.clusterNumData, newBeta, dp.beta,dp.alpha));
				}
				
				minCondLike = condLike[0];
				for (cc=1;cc<numClusters+1;cc++) {
					if (condLike[cc] < minCondLike) {
						minCondLike = condLike[cc];
					}
				}
				
				for (cc=0;cc<numClusters;cc++) {
					condLike[cc] = condLike[cc] - minCondLike;
					parent = controllingProcess.children.get(cc);
					if (groupDocuments) {
						if (doc.parent == parent) {
							tmp = parent.numDocuments() - 1.0;
							if (tmp == 0.0) { tmp = alpha; }
						} else { tmp = parent.numDocuments(); }
					} else {
						if (dp.parent == parent) {
							tmp = parent.numChildren() - 1.0;
							if (tmp == 0.0) { tmp = alpha; }
						} else { tmp = parent.numChildren(); }
					}
					condLike[cc] = tmp*Math.exp(condLike[cc]);
				}
				condLike[numClusters] = condLike[numClusters] - minCondLike;
				condLike[numClusters] = alpha*Math.exp(condLike[numClusters]);
				
				if (groupDocuments) {
					oldParent = doc.parent;
				} else {
					oldParent = dp.parent;
				}
			
				cc = StatUtil.multinomial_rnd(condLike,numClusters+1);
				
				if (cc == numClusters) {
					numClusters++;
					condLike = new double[numClusters+1];
					parent = newDirichletProcess(newBeta,controllingProcess.children.get(0).concParam);
				} else {
					parent = controllingProcess.children.get(cc);
				}
				
				if (oldParent != parent) {
				
				//	dp.parent = parent;
				//	oldParent.removeChild(dp);
				//	parent.addChild(dp);
					if (groupDocuments) {
						doc.parent = parent;
						oldParent.removeDocument(doc);
						parent.appendDocument(doc);
					} else {
						dp.parent = parent;
						oldParent.removeChild(dp);
						parent.addChild(dp);
					}
					for (tt=0;tt<numTopics;tt++) {
						if (groupDocuments) {
							oldParent.clusterNumData[tt] = oldParent.clusterNumData[tt] - doc.clusterNumData[tt];
							parent.clusterNumData[tt] = parent.clusterNumData[tt] + doc.clusterNumData[tt];
						} else {
							oldParent.clusterNumData[tt] = oldParent.clusterNumData[tt] - dp.clusterNumTables[tt];
							parent.clusterNumData[tt] = parent.clusterNumData[tt] + dp.clusterNumTables[tt];
						}
					}
					parent.resampleBeta(numTopics,newBeta);
					if (oldParent.numChildren() > 0) {
						oldParent.resampleBeta(numTopics,newBeta);
					}
				}
			}
		}
		removeEmptyDPs(numTopics);
	}
	
//	 computes the data likelihood given that the data is generated by the prior
	// dp_alpha*parent_beta (child_beta is NOT integrated out)
	public double dataLikeDP2(int numTopics,int[] clusterNumData, double[] parent_beta, double[] dp_beta, double dp_alpha) {
		double lik = 0.0;
		int tt = 0;
		int j = 0;
		int nd = 0;
		double tol = 1e-40;
		
		double alpha_0 = 0.0;
		for(tt=0;tt<numTopics;tt++) {
			alpha_0 = alpha_0 + parent_beta[tt];
		}
		alpha_0 = alpha_0*dp_alpha;
		lik = lik + cern.jet.stat.Gamma.logGamma(alpha_0);
		
		for (tt=0;tt<numTopics;tt++) {
			nd = clusterNumData[tt];
			if (dp_beta[tt] > tol) {
			//	lik = lik + (parent_beta[tt]*dp_alpha - 1.0 + ((double) nd))*Math.log(dp_beta[tt]);
				lik = lik + (parent_beta[tt]*dp_alpha - 1.0)*Math.log(dp_beta[tt]);
			}
			if (parent_beta[tt] > tol) {
				lik = lik - cern.jet.stat.Gamma.logGamma(parent_beta[tt]*dp_alpha);
			}
		}
		
		return lik;
	}
	
	// computes the data likelihood given that the data is generated by the prior
	// dp_alpha*parent_beta (child_beta is integrated out)
	public double dataLikeDP(int numTopics,int[] clusterNumData, double[] parent_beta, double dp_alpha) {
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
		}
		
	/*	double alpha_0 = 0.0;
		int td = 0;
		for(tt=0;tt<numTopics;tt++) {
			alpha_0 = alpha_0 + parent_beta[tt];
			td = td + clusterNumData[tt];
		}
		alpha_0 += parent_beta[numTopics];
		alpha_0 = alpha_0*dp_alpha;
		lik = lik + cern.jet.stat.Gamma.logGamma(alpha_0) - cern.jet.stat.Gamma.logGamma(alpha_0 + (double) td); */
		
		return lik;
	}
	
	// computes the data likelihood given that the document is generated by a particular DP w/
	// topic distribution beta
	public double dataLikeDocument(int numTopics,int[] clusterNumData, double[] beta) {
		double lik = 0.0;
		int tt = 0;
		
		for (tt=0;tt<numTopics;tt++) {
			if (clusterNumData[tt] > 0) {
				lik = lik + clusterNumData[tt]*Math.log(beta[tt]);
			}
		}
		
		return lik;
	}
	
	public void removeEmptyDPs(int numTopics) {
		int nc = controllingProcess.numChildren();
		int cc = 0;
		int tt = 0;
		DirichletProcess dp = null;
		for (cc=nc-1;cc>=0;cc--) {
			dp = controllingProcess.children.get(cc);
			if ((dp.numDocuments() == 0 && groupDocuments) || (dp.numChildren() == 0 && !groupDocuments)) {
				dp.concParam.removeDP(dp);
				controllingProcess.removeChild(dp);
				// if this process was initially there, it will need to be removed
				for (tt=0;tt<numTopics;tt++) {
					controllingProcess.clusterNumData[tt] = controllingProcess.clusterNumData[tt] - dp.clusterNumTables[tt];
				}
			}
		}
	}
	
	public DirichletProcess newDirichletProcess(double[] newBeta,HDPConcParam concParam) {
		DirichletProcess dp = new DirichletProcess(controllingProcess,"");
		double[] beta = new double[controllingProcess.beta.length];
		System.arraycopy(newBeta,0,beta,0,newBeta.length);
		int[] clusterNumData = new int[controllingProcess.clusterNumData.length];
	  	int[] clusterNumTables = new int[controllingProcess.clusterNumTables.length];
	  	dp.beta = beta;
	  	dp.clusterNumData = clusterNumData;
	  	dp.clusterNumTables = clusterNumTables;
  		
  		dp.alpha = concParam.alpha;
  		dp.concParam = concParam;
  		concParam.appendDP(dp);
	  	
		return dp;
	}
	
	// sample concentration parameter using aux. variable sampling scheme
	// from Escobar and West (1995)
	void resampleAlpha(int numIters) {
		int i = 0;
		double eta = 0.0;
		int zz = 0;
		double[] w = new double[2];
		numClusters = controllingProcess.numChildren();
		if (numClusters <= 1) {
			return;
		}
		
		int numPresent = 0;
		if (!groupDocuments) {
			for (i=0;i<groupableDPs.size();i++) {
				if (groupableDPs.get(i).state != DirichletProcess.HELDOUT) {
					numPresent++;
				}
			}
		}
		if (groupDocuments) {
			for (i=0;i<groupableDocs.size();i++) {
				if (groupableDocs.get(i).parent.state != DirichletProcess.HELDOUT) {
					numPresent++;
				}
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
	void incPairProbs() {
		int cc = 0;
		int ii = 0;
		int jj = 0;
		int nd = 0;
		int ID1 = 0;
		int ID2 = 0;
		Document[] docs = null;
		ArrayList<DirichletProcess> dps = null;
		
		numClusters = controllingProcess.numChildren();
		
		for (cc=0;cc<numClusters;cc++) {
			if (groupDocuments) {
				docs = controllingProcess.children.get(cc).documents;
				nd = docs.length;
			} else {
				dps = controllingProcess.children.get(cc).children;
				nd = dps.size();
			}
			for (ii=0;ii<nd-1;ii++) {
				if (groupDocuments) { 
					ID1 = docs[ii].docID;
				} else { ID1 = dps.get(ii).dpID; }
				for (jj=(ii+1);jj<nd;jj++) {
					if (groupDocuments) { 
						ID2 = docs[jj].docID;
					} else { ID2 = dps.get(jj).dpID; }
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
	
	void normPairProbs() {
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
	void writePairwiseProbs(String fname) throws IOException {
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
		if (groupDocuments) {
			for (ii=0;ii<numGroupable;ii++) {
				outFile.write(groupableDocs.get(ii).annotation);
				outFile.write("\n");
			}
		} else
		{
			for (ii=0;ii<numGroupable;ii++) {
				outFile.write(groupableDPs.get(ii).getLabel());
				outFile.write("\n");
			}
		}
		
		outFile.close();
	}
	
	// this will find a concensus model from the pair-wise probabilities
	// and write the discovered clusters to a file
  	void writeConcensusModel(String fname) throws IOException {
  		// perform hierarhical clustering, using 1.0-pair-wise probabilities as
  		// distances, cutting clusters with distances greater than 0.50
  		HierarchicalCluster hier = new HierarchicalCluster();
  		hier.cluster(pairProb,0.50);
  		
  		numClusters = hier.clusters.size();
  		
  		int cc = 0;
  		int gg = 0;
  		int itemNum = 0;
  		Document doc = null;
  		DirichletProcess dp = null;
  		
  		FileWriter outFile = new FileWriter(fname);
  		String s = "";
  		
  		// reconstruct the topics
  		for (cc=0;cc<numClusters;cc++) {
  			s = "Cluster #" + (cc+1) + "\n";
  			outFile.write(s);
  			for (gg=0;gg<hier.clusters.get(cc).size();gg++) {
  				itemNum = hier.clusters.get(cc).get(gg);
  				if (groupDocuments) {
  					doc = groupableDocs.get(itemNum);
  					s = doc.annotation + "\n";
  				} else{
  					dp = groupableDPs.get(itemNum);
  					s = dp.getLabel() + "\n";
  				}
  				outFile.write(s);
  			}
  			outFile.write("\n\n");
  		}
  		
  		outFile.close();
  	}
  	
  	// this will build a hierarchical DP based on the concensus model
  	// currently doesn't work if groupDocuments = true
  	void buildConcensusModel() {
  		HierarchicalCluster hier = new HierarchicalCluster();
  	//	hier.cluster(pairProb,0.10);
  		hier.cluster(pairProb,0.50);
  		
  		numClusters = hier.clusters.size();
  		
  		int cc = 0;
  		int gg = 0;
  		int itemNum = 0;
  		Document doc = null;
  		DirichletProcess dp = null;
  		DirichletProcess parent = null;
  		double[] newBeta = null;
  		
  		// remove the groups - we'll add them back
  		HDPConcParam cp = controllingProcess.children.get(0).concParam;
  		controllingProcess.children.clear();
  /*		for (cc=0;cc<controllingProcess.children.size();cc++) {
  			dp = controllingProcess.children.get(cc);
  			controllingProcess.removeChild(dp);
  			dp.concParam.removeDP(dp);
  		} */
  		
  		// reconstruct the groups
  		for (cc=0;cc<numClusters;cc++) {
  			newBeta = new double[controllingProcess.beta.length];
  			parent = newDirichletProcess(newBeta,cp);
  			for (gg=0;gg<hier.clusters.get(cc).size();gg++) {
  				itemNum = hier.clusters.get(cc).get(gg);
  				dp = groupableDPs.get(itemNum);
  				dp.parent = parent;
				parent.addChild(dp);
  			}
  		}
  	}
}
