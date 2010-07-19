package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.IOException;
import java.util.ArrayList;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

// cluster Macklis data from Arlotta et al paper
public class MacklisModel {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		String fname = "C:\\CPPDataFiles\\macklis_data.txt";
		
		MicroArrayData myData = new MicroArrayData();
		try {
			myData.readFile(fname);
		} catch(IOException e) {
			System.out.println(e);
			System.exit(0);
		}
		
		int i = 0;
		int j = 0;
		int numConditions = 0;
		int cond = 0;
		
		for (j=0;j<myData.numCols;j++) {
			i = myData.experimentCode[j];
			if (i>numConditions) {
				numConditions = i;
			}
		}
		
	numConditions++;
	
	HDPConcParam[] concParams = new HDPConcParam[2+numConditions];
	HDPConcParam gamma_p = new HDPConcParam(1.0,0.5);
	HDPConcParam[] alpha_1 = new HDPConcParam[numConditions];
	HDPConcParam alpha_0 = new HDPConcParam(1.0,1.0);
	
	for (i=0;i<numConditions;i++) {
		alpha_1[i] = new HDPConcParam(1.0,0.5);
	}
	
	concParams[0] = gamma_p;
	for (i=0;i<numConditions;i++) {
		concParams[i+1] = alpha_1[i];
	}
	
	concParams[2+numConditions-1] = alpha_0;
	
	Document doc = null;
	ArrayList<DirichletProcess> DPs = new ArrayList<DirichletProcess>();
	ArrayList<DirichletProcess> alpha_0_DP = new ArrayList<DirichletProcess>();
		
	DirichletProcess root = new DirichletProcess(null,"root");
	gamma_p.addDP(root);
	DPs.add(root);
	DirichletProcess[] condDP = new DirichletProcess[numConditions];
	
	String s = "";
	for(i=0;i<numConditions;i++) {
		if (i == numConditions-1) {
			s = "Common";
		} else {
			s = "Cond#" + (i+1);
		}
		condDP[i] = new DirichletProcess(root,s);
		DPs.add(condDP[i]);
		alpha_1[i].addDP(condDP[i]);
	}
	
	DirichletProcess d = null;
	
	extractCommon(myData,condDP,alpha_0_DP,DPs);
	
	alpha_0.addDPs(alpha_0_DP);
	
	HierDirichletProcess myHierDP = new HierDirichletProcess(100,DPs,concParams,myData.geneNames);
	
	myHierDP.packGenes();
	myHierDP.activateAll();
	
	int iters = 50000;
	int burnin = 5000;
	int targetNumClusters = 119;
	int numConcParamIter = 15;
	
	String mainOutName = "C:\\CPPDataFiles\\macklis_clusters.txt";
//	myHierDP.iterate(iters,burnin,numConcParamIter,targetNumClusters,mainOutName);
	
	try {
	//	myHierDP.outputClustersToFile(mainOutName);
		myHierDP.outputNumClustersSamples("C:\\CPPDataFiles\\macklis_numclusters.txt");
	} catch(IOException e) { System.out.println(e);	}

}

	public static void extractCommon(MicroArrayData myData,DirichletProcess[] condDP,ArrayList<DirichletProcess> alpha_0_DP,ArrayList<DirichletProcess> DPs) {
		int i = 0;
		int j = 0;
		int k = 0;
		int gg = 0;
		int idx = 0;
		int maxTFno = 0;
		int[] commonBind = new int[myData.numRows];
		int[][] uniqueBind = null;
		Document doc = null;
		DirichletProcess d = null;
		int cond = 0;
		
		for (i=0;i<myData.numCols;i++) {
			if (myData.experimentCode2[i] > maxTFno) {
				maxTFno = myData.experimentCode2[i];
			}
		}
		
		ArrayList<Integer> matchedExperiments = new ArrayList<Integer>();
		ArrayList<Integer> glist = new ArrayList<Integer>();
		
		// minimum number of genes that must be bound in common
		// to create a common profile
		int minNumBound = 5;
		int maxCommonBind = 1;
		int numBoundCommon = 0;
		int numBound = 0;
		String s = "";
		for(i=1;i<maxTFno+1;i++) {
			matchedExperiments.clear();
			for(j=0;j<myData.numCols;j++) {
				if (myData.experimentCode2[j] == i) {
					matchedExperiments.add(new Integer(j));
				}
			}
			uniqueBind = new int[myData.numRows][matchedExperiments.size()];
			numBoundCommon = 0;
			if (matchedExperiments.size() == 1) {
				idx = matchedExperiments.get(0);
				for (j=0;j<myData.numRows;j++) {
					uniqueBind[j][0] += (int) myData.values[j][idx];
				}
			} else {
				for (j=0;j<myData.numRows;j++) {
					commonBind[j] = 0;
				}
				for (j=0;j<myData.numRows;j++) {
					maxCommonBind = 10000000;
					for (k=0;k<matchedExperiments.size();k++) {
						idx = matchedExperiments.get(k);
						uniqueBind[j][k] += (int) myData.values[j][idx];
						if (uniqueBind[j][k] < maxCommonBind) {
							maxCommonBind = uniqueBind[j][k];
						}
						if (myData.values[j][idx] > 0) {
							commonBind[j]++;
						}
					}
					if (commonBind[j] == matchedExperiments.size()) {
						commonBind[j] = maxCommonBind;
					} else { commonBind[j] = 0; }
				}
				numBoundCommon = 0;
				for (j=0;j<myData.numRows;j++) {
					if (commonBind[j] > 0) {
						numBoundCommon++;
					}
				}
				// generate common binding
				if (numBoundCommon >= minNumBound) {
					glist.clear();
					for(j=0;j<myData.numRows;j++) {
						if (commonBind[j] > 0) {
							for(gg=0;gg<commonBind[j];gg++) { glist.add(j); }
						}
						for (k=0;k<matchedExperiments.size();k++) {
							uniqueBind[j][k] = uniqueBind[j][k] - commonBind[j];
						}
					}
					s = "common";
					for (k=0;k<matchedExperiments.size();k++) {
						idx = matchedExperiments.get(k);
						s = s + "_" + myData.experimentNames[idx];
					}
					d = new DirichletProcess(condDP[condDP.length-1],s);
					for (k=0;k<matchedExperiments.size();k++) {
						idx = matchedExperiments.get(k);
						doc = new Document("common_" + myData.experimentNames[idx],glist);
						d.addDocument(doc);
					}
					DPs.add(d);
					alpha_0_DP.add(d);
				//	System.out.println(s + " " + numBoundCommon + " genes in common");
				}
			}
			// generate unique binding
			for (k=0;k<matchedExperiments.size();k++) {
				idx = matchedExperiments.get(k);
				numBound = 0;
				glist.clear();
				for (j=0;j<myData.numRows;j++) {
					if (uniqueBind[j][k] > 0) {
						numBound++;
						for (gg=0;gg<uniqueBind[j][k];gg++) { glist.add(j); }
					}
				}
				if (numBound >= minNumBound) {
					if (matchedExperiments.size() > 1 && numBoundCommon >= minNumBound) {
						s = "unique_" + myData.experimentNames[idx];
					} else {
						s = myData.experimentNames[idx];
					}
					cond = myData.experimentCode[idx]-1;
					doc = new Document(myData.experimentNames[i],glist);
					d = new DirichletProcess(condDP[cond],s);
					d.addDocument(doc);
					DPs.add(d);
					alpha_0_DP.add(d);
				//	System.out.println(s + " " + numBound + " unique genes");
				}
			}
		}
	}

}
