package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import cern.jet.random.Beta;
import cern.jet.random.Gamma;
import cern.jet.random.Uniform;
import edu.mit.csail.psrg.georg.HierarchicalCluster.*;
import edu.mit.csail.psrg.georg.StatUtil.StatUtil;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.io.ObjectOutputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.FileInputStream;
import java.math.BigInteger;


public class HierDirichletProcess implements Serializable {
	//	maps original gene numbers to re-ordered genes
	int[] geneMap;
	// maps re-ordered genes to original gene numbers
	int[] reverseGeneMap;
	// total unique genes
	int totalGenes = 0;
	
	// topologically sorted list of DP's
	DirichletProcess[] DP;
	
	// a DPGroup is a DP node and children,
	// such that the children of the children (document-level DPs)
	// can reassort among the children nodes (and new children nodes can be created)
	ArrayList<DPGroup> DPGroups = new ArrayList<DPGroup>();
	
	// these items refer to the base distribution
	double priorPos;
	double priorPos_a = 0.5;
	double priorPos_b = 0.5;
	int numClusters = 0;
	
	// the number of clusters across all samples
	ArrayList<Integer> numClustersSamples = new ArrayList<Integer>();
	
	// allocated size of the cluster array
	int allocClusters = 5;
	// amount to grow cluster allocation by
	int growClusterAlloc = 10;
	
	// how many data items (genes) are associated with each cluster
	int[] clusterCounts;
	// how many positive counts (for each gene) are associated with each cluster
	int[][] posCounts;
	
	// used to store intermediate conditional likelihood computations
	double[] condLike;
	
	// used to store concentration parameters
	HDPConcParam[] concParams;
	
	String[] geneNames = null;
	// temporary cluster map used for evaluating document loading on topics
	HashSet[] tempClusterMap;
	// temporary array for evaluating document loading on topics
	ArrayList[] tempClusterArray;
	
	// variables to control snap-shots (saving of the whole HDP object
	// every N # of samples)
	int snapShotInterval = 1000;
	String snapShotFileName = null;
	// variable to control # of iterations after which to start sampling groups
	int groupStartIter = 0;
	
	HashMap<Class,HDPStatCollector> statCollectors = new HashMap<Class,HDPStatCollector>();
	
	HierDirichletProcess(int nc, ArrayList<DirichletProcess> dd, HDPConcParam[] cp,String[] gn) {
		priorPos = 0.001;
		numClusters = nc;
		concParams = cp;
		geneNames = gn;
		
		DP = new DirichletProcess[dd.size()];
		int jj = 0;
		for(jj=0;jj<dd.size();jj++) {
			DP[jj] = dd.get(jj);
		}
	}
	
	public void addStatCollector(HDPStatCollector c) {
		statCollectors.put(c.getClass(),c);
	}
	
	public HDPStatCollector getStatCollector(Class c) {
		return statCollectors.get(c);
	}
	
	public void enableSnapShots(int i,String f) {
		snapShotInterval = i;
		snapShotFileName = f;
	}
	
	public void disableSnapShots() {
		snapShotFileName = null;
	}
	
	void initAllocClusters() {
		int totalAllocate = allocClusters;
		if (allocClusters < numClusters+2) {
			totalAllocate = numClusters+growClusterAlloc;
		}
		
		clusterCounts = new int[totalAllocate];
		posCounts = new int[totalAllocate][totalGenes];
		condLike = new double[totalAllocate];
		tempClusterMap = new HashSet[totalAllocate];
		int i = 0;
		for (i=0;i<totalAllocate;i++) {
			tempClusterMap[i] = new HashSet();
		}
		
		tempClusterArray = new ArrayList[totalAllocate];
		for (i=0;i<totalAllocate;i++) {
			tempClusterArray[i] = new ArrayList();
		}
		
		allocClusters = totalAllocate;
	}
	
	void addDPGroup(DPGroup dg) {
		DPGroups.add(dg);
	}
	
	void addCluster() {
		int jj = 0;
		
		if (clusterCounts.length < numClusters+2) {
			int[] clusterCounts2 = new int[allocClusters+growClusterAlloc];
			System.arraycopy(clusterCounts, 0, clusterCounts2, 0, allocClusters);
			clusterCounts = clusterCounts2;
			
			int[][] posCounts2 = new int[allocClusters+growClusterAlloc][totalGenes];
			System.arraycopy(posCounts, 0, posCounts2, 0, allocClusters);
			posCounts = posCounts2;
			
			double[] condLike2 = new double[allocClusters+growClusterAlloc];
			System.arraycopy(condLike,0,condLike2,0,allocClusters);
			condLike = condLike2;
			
			HashSet[] tempClusterMap2 = new HashSet[allocClusters+growClusterAlloc];
			System.arraycopy(tempClusterMap,0,tempClusterMap2,0,allocClusters);
			for (jj=allocClusters;jj<allocClusters+growClusterAlloc;jj++) {
				tempClusterMap2[jj] = new HashSet();
			}
			tempClusterMap = tempClusterMap2;
			
			ArrayList[] tempClusterArray2 = new ArrayList[allocClusters+growClusterAlloc];
			System.arraycopy(tempClusterArray,0,tempClusterArray2,0,allocClusters);
			for (jj=allocClusters;jj<allocClusters+growClusterAlloc;jj++) {
				tempClusterArray2[jj] = new ArrayList();
			}
			tempClusterArray = tempClusterArray2;
		
			for (jj=0;jj<DP.length;jj++) {
				if (DP[jj].state != DirichletProcess.HELDOUT) {
					DP[jj].growClusters(growClusterAlloc);
				}
			}
			allocClusters = clusterCounts.length;
		}
		
		clusterCounts[numClusters+1] = 0;
		
		int k = 0;
		for (k=0;k<totalGenes;k++) {
			posCounts[numClusters+1][k] = 0;
		}
		
		// init appropriate vectors in each DP
		int dd = 0;
		for (jj=0;jj<DP.length;jj++) {
			if (DP[jj].state != DirichletProcess.HELDOUT) {
				DP[jj].beta[numClusters+1] = 0;
				DP[jj].clusterNumData[numClusters+1] = 0;
				DP[jj].clusterNumTables[numClusters+1] = 0;
				if (DP[jj].documents != null) {
					for (dd=0;dd<DP[jj].documents.length;dd++) {
						DP[jj].documents[dd].clusterNumData[numClusters+1] = 0;
					}
				}
			}
		}
		
		// stick-breaking construction for generating
		// next beta weights
		for (jj=0;jj<DP.length;jj++) {
			if (DP[jj].state != DirichletProcess.HELDOUT) {
				DP[jj].sampleNewBeta(numClusters);
			}
		}
		
		numClusters++;
	}
	
	void deleteCluster(int c) {
		// the cluster to be deleted is zero-indexed, as is standard in Java
		// this method will copy numClusters+1 entries
		
		System.arraycopy(clusterCounts,c+1,clusterCounts,c,numClusters-c);
		
		int k = 0;
		
		int[] pc = posCounts[c];
		for(k=c;k<numClusters;k++) {
			posCounts[k] = posCounts[k+1];
		}
		posCounts[numClusters] = pc;
		
		HashSet cm = tempClusterMap[c];
		for(k=c;k<numClusters;k++) {
			tempClusterMap[k] = tempClusterMap[k+1];
		}
		tempClusterMap[numClusters] = cm;
		//System.arraycopy(posCounts,c+1,posCounts,c,numClusters-c);
		//System.arraycopy(clusterMap,c+1,clusterMap,c,numClusters-c);
		
		ArrayList ca = tempClusterArray[c];
		for(k=c;k<numClusters;k++) {
			tempClusterArray[k] = tempClusterArray[k+1];
		}
		tempClusterArray[numClusters] = ca;
		
		int jj = 0;
		for(jj=0;jj<DP.length;jj++) {
			if (DP[jj].state != DirichletProcess.HELDOUT) {
				DP[jj].deleteCluster(c,numClusters);
			}
		}
		
		numClusters--;
	}
	
  	void sampleDataClusterAssigns(DirichletProcess myDP) {
  		int dd = 0;
  		int gg = 0;
  		Document doc = null;
  		int cc = 0;
  		int oldCC;
  		int gene = 0;
  		int[] genes;
  		int[] geneClusterAssigns;
  		double[] beta = null;
  		double mlt = 0.0;
  		double alpha = myDP.alpha;
  		DirichletProcess parent = myDP.parent;
  		
  		if (parent != null) {
  			beta = parent.beta;
  		}
  		
  		int[] clusterNumData = myDP.clusterNumData;
  		Document[] documents = myDP.documents;
  		
  		for (dd=0;dd<documents.length;dd++) {
  			doc = documents[dd];
  			genes = doc.genes;
  			geneClusterAssigns = doc.geneClusterAssigns;
  			for (gg=0;gg<genes.length;gg++) {
  				oldCC = geneClusterAssigns[gg];
  				gene = genes[gg];
  				
  				// remove the gene from consideration
  				posCounts[oldCC][gene]--;
  				clusterCounts[oldCC]--;
  				clusterNumData[oldCC]--;
  				doc.clusterNumData[oldCC]--;
  				
  				for (cc=0;cc<=numClusters;cc++) {
  					// this is the marginal conditional likelihood that the gene
  					// comes from the cluster (we marginalize out the parameters under
  					// a beta prior)
  					condLike[cc] = (((double) posCounts[cc][gene]) + priorPos)/(((double) clusterCounts[cc]) + priorPos*totalGenes);
  					mlt = clusterNumData[cc];
  					if (parent == null) {
  						if (cc == numClusters) {
  							mlt = mlt + alpha;
  						}
  					}
  					else {
  						mlt = mlt + alpha*beta[cc];
  					}
  					condLike[cc] = condLike[cc]*mlt;
  				}
  				
  				if (!myDP.generateNewClusters) {
  					condLike[numClusters] = 0;
  				}
  				
  				cc = StatUtil.multinomial_rnd(condLike,numClusters+1);
  				
  				// add the gene to the cluster it just joined
  				clusterCounts[cc]++;
  				posCounts[cc][gene]++;
  				clusterNumData[cc]++;
  				doc.clusterNumData[cc]++;
  				geneClusterAssigns[gg] = cc;
  				
  				// if it's a new cluster, we'll have to expand the # of clusters
  				if (cc==numClusters) {
  					addCluster();
  					if (parent != null) {
  			  			beta = parent.beta;
  			  		}
  			  		
  			  		clusterNumData = myDP.clusterNumData;
  				}
  				
  				if (cc != oldCC) {
  				//	clusterMap[oldCC].remove(new Integer(geneID));
  				//	clusterMap[cc].add(new Integer(geneID));
  				}
  			}
  		}
  	}
  	
  	void resampleConcParams(int numiter) {
  		int cp = 0;
  		for (cp=0;cp<concParams.length;cp++) {
  			concParams[cp].resample(numiter);
  		}
  		
  	/*	System.out.print("Conc params=[");
  		for (cp=0;cp<concParams.length;cp++) {
  			System.out.print(" " + concParams[cp].alpha);
  		}
  		System.out.println(" ]"); */
  	}
  	
  	void resamplePriorPos(int numiter) {
  		int iter, cc, gg, nt;
	  	double aa, bb, xx, nd, zz;
	  	priorPos = priorPos*((double) totalGenes);
	  	for (iter = 0 ;iter < numiter ;iter++ ) {
	    	aa = priorPos_a;
	    	bb = priorPos_b;
	    	for ( cc = 0 ; cc < numClusters; cc++ ) {
	    		nd = (double) clusterCounts[cc];
	    		if (nd <= 0) {
	      			xx = 1.0;
	      		} else {
	      			xx = Beta.staticNextDouble(priorPos+1.0, nd);
	      		}
	      		if ( Uniform.staticNextDouble() * (priorPos + nd) < nd ) {
	      			zz = 1;
	      		}
	      		else {
	      			zz = 0;
	      		}
	    		nt = 0;
	    		for (gg=0;gg<totalGenes;gg++) {
	    			if (posCounts[cc][gg] > 0) {
	    				nt++;
	    			}
	    		}
	      		aa += ((double) nt)-zz;
	      		bb -= Math.log(xx);
	    	}
	    	if (aa >= 1.0) {
	    		priorPos = Gamma.staticNextDouble(aa,1.0) / bb;
	    	}
	  	}
	  	priorPos = priorPos/((double) totalGenes);
	//  	System.out.println("priorPos="+priorPos);
  	}
  	
  	void packGenes() {
  		int jj = 0;
  		int kk = 0;
  		int ll = 0;
  		int gene = 0;
  		Document doc = null;
  		ArrayList<Integer> genes = new ArrayList<Integer>();
  		int maxGene = 0;
  		
  		for (jj=0;jj<DP.length;jj++) {
  			if (DP[jj].documents != null) {
	  			for (kk=0;kk<(DP[jj].documents.length);kk++) {
	  				doc = DP[jj].documents[kk];
	  				for (ll=0;ll<(doc.genes.length);ll++) {
	  					gene = doc.genes[ll];
	  					
	  					if (gene > maxGene) {
	  						maxGene = gene;
	  					}
	  					
	  					// the gene was not found
	  					if (StatUtil.find(genes,gene) == -1) {
	  						genes.add(gene);
	  					}
	  				}
	  			}
  			}
  		}
  		
  		totalGenes = genes.size();
  		priorPos = 1.0/((double) totalGenes);
  		priorPos = priorPos*10.0;
  		priorPos_a = 1.0;
  		priorPos_b = 1.0;
  		
  		reverseGeneMap = new int[totalGenes];
  		geneMap = new int[maxGene+1];
  		
  		for (jj=0;jj<genes.size();jj++) {
  			geneMap[genes.get(jj)] = jj;
  			reverseGeneMap[jj] = genes.get(jj);
  		}
  		
  		for (jj=0;jj<(DP.length);jj++) {
  			if (DP[jj].documents != null) {
	  			for (kk=0;kk<(DP[jj].documents.length);kk++) {
	  				doc = DP[jj].documents[kk];
	  				for (ll=0;ll<(doc.genes.length);ll++) {
	  					gene = doc.genes[ll];
	  					doc.genes[ll] = geneMap[gene];
	  				}
	  			}
  			}
  		}
  	}
  	
  	void iterate(int numIters,int burnin) {
  		iterate(numIters,burnin,15);
  	}

  	void iterate(int numIters,int burnin,int concParamIter) {
  		int iters = 0;
  		int jj = 0;
  		int cc = 0;
  		int cp = 0;
  		double likelihood = 0.0;
  		double maxLikelihood = -1.0/0.0;
  		int nc = 0;
  		
  		String snapShotOutName = "";
  		
  		for (iters=0;iters<numIters;iters++) {
  			for (jj=DP.length-1;jj>=0;jj--) {
  				if (DP[jj].documents != null & DP[jj].state == DirichletProcess.ACTIVE) {
  					sampleDataClusterAssigns(DP[jj]);
  				}
  				if (DP[jj].state == DirichletProcess.ACTIVE) {
  					DP[jj].sampleClusterNumTables(numClusters);
  				}
  			}
  			
  			for (cp=0;cp<concParams.length;cp++) {
  				concParams[cp].calcTotals();
  			}
  			
  			for (jj=0;jj<DP.length;jj++) {
  				if (DP[jj].state == DirichletProcess.ACTIVE) {
  					DP[jj].resampleBeta(numClusters,condLike);
  				}
  			} 
  			
  			// now delete any empty clusters
  			nc = numClusters;
  			for (cc=nc-1;cc>=0;cc--) {
  				if (clusterCounts[cc] == 0) {
  					deleteCluster(cc);
  				}
  			}
  			
  			int bicluster = 1;
  			if (!DPGroups.isEmpty() & iters > groupStartIter) {
  				for (cc=0;cc<bicluster;cc++) {
	  				for (jj=0;jj<DPGroups.size();jj++) {
	  					DPGroups.get(jj).sampleDPAssignments(numClusters);
	  					DPGroups.get(jj).resampleAlpha(concParamIter);
	  					System.out.println("grps=" + DPGroups.get(jj).numClusters + " alpha=" + DPGroups.get(jj).alpha);
	  				}
	  				buildDPList();
	  				
	  				for (jj=DP.length-1;jj>=0;jj--) {
	  					if (DP[jj].state == DirichletProcess.ACTIVE) {
	  						DP[jj].sampleClusterNumTables(numClusters);
	  					}
	  	  			}
	  				for (cp=0;cp<concParams.length;cp++) {
	  	  				concParams[cp].calcTotals();
	  	  			}
	  				for (jj=0;jj<DP.length;jj++) {
	  					if (DP[jj].state == DirichletProcess.ACTIVE) {
	  						DP[jj].resampleBeta(numClusters,condLike);
	  					}
		  			}
	  			//	resampleConcParams(concParamIter); 
	  			/*	for (cp=0;cp<concParams.length;cp++) {
	  	  				concParams[cp].calcTotals();
	  	  			} */
  				}
  			}
  			
  			resampleConcParams(concParamIter);
  			resamplePriorPos(concParamIter);
  			System.out.println(iters + " " + numClusters);
  			
  			numClustersSamples.add(numClusters);
  			
  			if (iters >= burnin) {
  				Iterator<HDPStatCollector> statIter = statCollectors.values().iterator();
  				while(statIter.hasNext()) {
  					statIter.next().updateStats();
  				}
  				if (!DPGroups.isEmpty()) {
  					for (jj=0;jj<DPGroups.size();jj++) {
	  					DPGroups.get(jj).incPairProbs();
	  				}
  				}
  			}
  			
  			if (snapShotFileName != null) {
  				double t = ((double) iters+1)/((double) snapShotInterval);
  				if (((double) Math.round(t)) == t) {
  					snapShotOutName = snapShotFileName + "_" + (new Integer(iters)).toString() + ".persist";
  					System.out.println("Snap shot at iteration" + iters);
  					persistSelf(snapShotOutName);
  				}
  			}
  		}
  	}
  	
  	double calcLogLike() {
  		int jj = 0;
  		int dd = 0;
  		int gg = 0;
  		Document doc = null;
  		int cc = 0;
  		int gene = 0;
  		int[] genes;
  		int[] geneClusterAssigns;
  		Document[] documents = null;
  		double lik = 0.0;
  		
  		// remove genes
  		for (jj=0;jj<DP.length;jj++) {
  			if (DP[jj].numDocuments() > 0 & DP[jj].state == DirichletProcess.ACTIVE) {
  				documents = DP[jj].documents;
		  		for (dd=0;dd<documents.length;dd++) {
		  			doc = documents[dd];
		  			genes = doc.genes;
		  			geneClusterAssigns = doc.geneClusterAssigns;
		  			for (gg=0;gg<genes.length;gg++) {
		  				cc = geneClusterAssigns[gg];
		  				gene = genes[gg];	  				
		  				posCounts[cc][gene]--;
		  				clusterCounts[cc]--;
		  			}
		  		}
	  		}
  		}
  		
  		// add genes back and calculate the likelihood
  		for (jj=0;jj<DP.length;jj++) {
  			if (DP[jj].numDocuments() > 0 & DP[jj].state == DirichletProcess.ACTIVE) {
  				documents = DP[jj].documents;
		  		for (dd=0;dd<documents.length;dd++) {
		  			doc = documents[dd];
		  			genes = doc.genes;
		  			geneClusterAssigns = doc.geneClusterAssigns;
		  			for (gg=0;gg<genes.length;gg++) {
		  				cc = geneClusterAssigns[gg];
		  				gene = genes[gg];	  				
		  				posCounts[cc][gene]++;
		  				clusterCounts[cc]++;
		  				lik = lik + Math.log((((double) posCounts[cc][gene]) + priorPos)/(((double) clusterCounts[cc]) + priorPos*totalGenes));
		  			}
		  		}
	  		}
  		}
  		return lik;
  	}

  	void activateAll() {
  		initAllocClusters();
  		int state = 0;
  		
  		int jj = 0;
  		for (jj=0;jj<DP.length;jj++) {
  			if (DP[jj].state == DirichletProcess.ACTIVE) {
  				DP[jj].activate(this);
  			}
  		}
  		
  		for (jj=0;jj<DP.length;jj++) {
  			if (DP[jj].state == DirichletProcess.ACTIVE) {
  				DP[jj].resampleBeta(numClusters,condLike);
  			}
  		}
  		
  		int cp = 0;
  		for (cp=0;cp<concParams.length;cp++) {
  			concParams[cp].calcTotals();
  		}
  	}
  	
  	void output() {
  		int cc = 0;
  		int gg = 0;
  		double pc = 0.0;
  		String s = "";
  		for (cc=0;cc<numClusters;cc++) {
  			s = "Topic #" + cc + "[ ";
  			for (gg=0;gg<totalGenes;gg++) {
  				pc = ((double) posCounts[cc][gg])/((double) clusterCounts[cc]);
  				s = s + pc + " ";
  			}
  			s = s + "]";
  			System.out.println(s);
  		}
  		
  		int jj = 0;
  		for (jj=0;jj<DP.length;jj++) {
  			System.out.println(DP[jj]);
  		}
  	}
  	
  	void outputClustersToFile(String baseFName) throws IOException {
  		String dpFName = baseFName + "_dps.txt";
  		
  		FileWriter outFile = new FileWriter(dpFName);
  		
  		int cc = 0;
  		int gg = 0;
  		int gdx = 0;
  		double pc = 0.0;
  		String s = "";
  		
  		outFile.write((new Integer(DP.length)).toString());
  		outFile.write("\t");
  		outFile.write((new Integer(numClusters)).toString());
  		outFile.write("\t1\t1\t0\t0\n");
  		
  		outFile.write("DP");
  		for (cc=0;cc<numClusters;cc++) {
  			s = "\tTopic #" + (cc+1);
  			outFile.write(s);
  		}
  		outFile.write("\n");
  		
  		int[][] tCt = new int[DP.length][numClusters];
  		HashMap<DirichletProcess,Integer> DPMap = new HashMap<DirichletProcess,Integer>();
  		
  		int jj = 0;
  		for (jj=0;jj<DP.length;jj++) {
  			DPMap.put(DP[jj],jj);
  		}
  		
  		int pIDX = 0;
  		for (jj=DP.length-1;jj>=0;jj--) {
  			for(cc=0;cc<numClusters;cc++) {
  				if (DP[jj].numDocuments() > 0) {
  					tCt[jj][cc] = tCt[jj][cc] + DP[jj].clusterNumData[cc];
  				}
  				if (DP[jj].parent != null) {
  					pIDX = DPMap.get(DP[jj].parent);
  					tCt[pIDX][cc] = tCt[pIDX][cc] + tCt[jj][cc];
  				}
  			}
  		}
  		
  		for (jj=0;jj<DP.length;jj++) {
  			outFile.write(DP[jj].getLabel());
  			for(cc=0;cc<numClusters;cc++) {
  				s = (new Integer(tCt[jj][cc])).toString();
  				outFile.write("\t");
  				outFile.write(s);
  			}
  			outFile.write("\n");
  		}
  		
  	/*	for (jj=0;jj<DP.length;jj++) {
  			outFile.write(DP[jj].getLabel());
  			for(cc=0;cc<numClusters;cc++) {
  				s = (new Double(DP[jj].beta[cc])).toString();
  				outFile.write("\t");
  				outFile.write(s);
  			}
  			outFile.write("\n");
  		} */
  		
  		outFile.close();
  		
  		String topicFName = baseFName + "_topics.txt";
  		outFile = new FileWriter(topicFName);
  		
  		outFile.write((new Integer(totalGenes)).toString());
  		outFile.write("\t");
  		outFile.write((new Integer(numClusters)).toString());
  		outFile.write("\t1\t1\t0\t0\n");
  		
  		outFile.write("Gene");
  		for (cc=0;cc<numClusters;cc++) {
  			s = "\tTopic #" + (cc+1);
  			outFile.write(s);
  		}
  		outFile.write("\n");
  		
  		for (gg=0;gg<totalGenes;gg++) {
  			gdx = reverseGeneMap[gg];
  			outFile.write(geneNames[gdx]);
  			for (cc=0;cc<numClusters;cc++) {
  			//	pc = ((double) posCounts[cc][gg])/((double) clusterCounts[cc]);
  			//	s = (new Double(pc)).toString();
  				s = (new Integer(posCounts[cc][gg])).toString();
  				outFile.write("\t");
  				outFile.write(s);
  			}
  			outFile.write("\n");
  		}
  		
  		outFile.close();
  	}
  	
  	void buildDPList(ArrayList<DirichletProcess> activeList,ArrayList<DirichletProcess> totalList) {
  		int i = 0;
  		int j = 0;
  		DirichletProcess dp = null;
  		ArrayList<DirichletProcess> activeList2 = new ArrayList<DirichletProcess>();
  		
  		for (i=0;i<activeList.size();i++) {
  			dp = activeList.get(i);
  			totalList.add(dp);
  			for(j=0;j<dp.numChildren();j++) {
  				activeList2.add(dp.children.get(j));
  			}
  		}
  		
  		if (!activeList2.isEmpty()) {
			buildDPList(activeList2,totalList);
		}
  	}
  	
  	void buildDPList() {
  		ArrayList<DirichletProcess> activeList = new ArrayList<DirichletProcess>();
  		ArrayList<DirichletProcess> totalList = new ArrayList<DirichletProcess>();
  		activeList.add(DP[0]);
  		buildDPList(activeList,totalList);
  		DirichletProcess[] DP2 = new DirichletProcess[totalList.size()];
  		int i = 0;
  		for (i=0;i<DP2.length;i++) {
  			DP2[i] = totalList.get(i);
  		}
  		DP = DP2;
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
  	
  	public static HierDirichletProcess restoreFromFile(String fName) {
  		HierDirichletProcess myHDP = null;
  		FileInputStream fis = null;
  		ObjectInputStream in = null;
  		try {
  			fis = new FileInputStream(fName);
  		    in = new ObjectInputStream(fis);
  		    myHDP = (HierDirichletProcess) in.readObject();
  		    in.close();
  		} catch(IOException e) {
  			System.out.println(e);
  		} catch(ClassNotFoundException e) {
  			System.out.println(e);
  		}
  		return myHDP;
  	}
  	
  	public void outputNumClustersSamples(String fName) throws IOException {
  		FileWriter outFile = new FileWriter(fName);
  		int i = 0;
  		for (i=0;i<numClustersSamples.size();i++) {
  			outFile.write(numClustersSamples.get(i).toString());
  			outFile.write("\n");
  		}
  		outFile.close();
  	}
}
