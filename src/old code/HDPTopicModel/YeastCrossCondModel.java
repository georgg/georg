package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.IOException;
import java.util.ArrayList;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

public class YeastCrossCondModel {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	//	CondsAsDocuments();
	//	GenesAsDocuments();
		combinedModel();;
	}
	
	public static void combinedModel() {
		String dirPath = "c:\\research_data\\combine_yeast_bind_expression\\";
		String discretizedData_fName = "yeast_combined_discretized_data.txt";
		MicroArrayData myData = new MicroArrayData();
		myData.setDiscrete();
		try {
			myData.readFile(dirPath + discretizedData_fName);
		} catch(IOException e) {
			System.out.println(e);
			System.exit(0);
		}
		
		int i = 0;
		int j = 0;
		int k = 0;
		int intensity = 0;
		int numConditions = 0;
		int cond = 0;
		
		boolean doCommon = true;
		
		for (j=0;j<myData.numCols;j++) {
			i = myData.experimentCode[j];
			if (i>numConditions) {
				numConditions = i;
			}
		}
		
	//	HDPConcParam[] concParams = new HDPConcParam[2+numConditions];
		HDPConcParam[] concParams = new HDPConcParam[3];
		HDPConcParam gamma_p = new HDPConcParam(1.0,0.2);
		HDPConcParam alpha_1 = new HDPConcParam(1.0,0.2);
	//	HDPConcParam[] alpha_1 = new HDPConcParam[numConditions];
		HDPConcParam alpha_0 = new HDPConcParam(1.0,1.0);
		ArrayList<DirichletProcess> alpha_1_DP = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> alpha_0_DP = new ArrayList<DirichletProcess>();
		
	/*	for (i=0;i<numConditions;i++) {
			alpha_1[i] = new HDPConcParam(1.0,0.2);
		} */
		
		concParams[0] = gamma_p;
	/*	for (i=0;i<numConditions;i++) {
			concParams[i+1] = alpha_1[i];
		}
		
		concParams[2+numConditions-1] = alpha_0; */
		concParams[1] = alpha_1;
		concParams[2] = alpha_0;
		
		Document doc = null;
		ArrayList<Integer> glist = new ArrayList<Integer>();
		ArrayList<DirichletProcess> DPs = new ArrayList<DirichletProcess>();
 		
		DirichletProcess root = new DirichletProcess(null,"root");
		gamma_p.addDP(root);
		DPs.add(root);
		DirichletProcess[] condDP = new DirichletProcess[numConditions];
		DirichletProcess[] condDPInitGroup = new DirichletProcess[numConditions];
		DPGroup[] groups = new DPGroup[numConditions];
		boolean groupDocuments = false;
		
		String s = "";
		for(i=0;i<numConditions;i++) {
			s = "Cond#" + (i+1);
			condDP[i] = new DirichletProcess(root,s);
			DPs.add(condDP[i]);
		//	alpha_1[i].addDP(condDP[i]);
			alpha_1_DP.add(condDP[i]);
		}
		
		if (groupDocuments) {
			for(i=0;i<numConditions;i++) {
				condDPInitGroup[i] = new DirichletProcess(condDP[i],"");
				DPs.add(condDPInitGroup[i]);
				alpha_0_DP.add(condDPInitGroup[i]);
			//	alpha_1[i].addDP(condDPInitGroup[i]);
			}
		}
		
		DirichletProcess d = null;
		
		for (j=0;j<myData.numCols;j++) {
			glist.clear();
			for (i=0;i<myData.numRows;i++) {
				if (myData.dvalues[i][j] > 0) {
					intensity = myData.dvalues[i][j];
					for (k=0;k<intensity;k++) {
						glist.add(i);
					}
				}
			}
			doc = new Document(myData.experimentNames[j],glist);
			cond = myData.experimentCode[j]-1;
			if (groupDocuments) {
			//	d = new DirichletProcess(condDP[cond],"");
				d = new DirichletProcess(condDPInitGroup[cond],myData.experimentNames[j]);
			} else {
				d = new DirichletProcess(condDP[cond],myData.experimentNames[j]);
			}
			d.addDocument(doc);
			DPs.add(d);
			alpha_0_DP.add(d);
		}
		
		alpha_1.addDPs(alpha_1_DP);
		alpha_0.addDPs(alpha_0_DP);
		
		HierDirichletProcess myHierDP = new HierDirichletProcess(1,DPs,concParams,myData.geneNames);
		
		if (groupDocuments) {
			for (i=0;i<numConditions;i++) {
			//	groups[i] = new DPGroup(condDP[i],true);
				groups[i] = new DPGroup(condDP[i],false);
				myHierDP.addDPGroup(groups[i]);
			}
		}
		
		myHierDP.packGenes();
		myHierDP.activateAll();
		
		int iters = 150000;
		int burnin = 50000;
		int targetNumClusters = 0;
		int numConcParamIter = 15;
		
		String mainOutName = dirPath + "yeast_combine";
		String snapShotName = mainOutName + "_snap";
		myHierDP.enableSnapShots(1000,snapShotName);
		
	//	myHierDP.iterate(iters,burnin,numConcParamIter,targetNumClusters,mainOutName);
		
		if (groupDocuments) {
			for (i=0;i<numConditions;i++) {
				groups[i].normPairProbs();
			}
		}
		
		try {
			myHierDP.outputClustersToFile(mainOutName);
			myHierDP.outputNumClustersSamples(mainOutName + "_numclusters.txt");
			if (groupDocuments) {
				for (i=0;i<numConditions;i++) {
					s = mainOutName + "_cond_" + (i+1);
					groups[i].writePairwiseProbs(s);
					s = s + ".txt";
					groups[i].writeConcensusModel(s);
				}
			}
		} catch(IOException e) { System.out.println(e);	}

	}
	
	public static void CondsAsDocuments() {
		double bindThresh = 0.005;
		String fname = "C:\\CPPDataFiles\\pvals_filtered.txt";
		MicroArrayData myData = new MicroArrayData();
		try {
			myData.readFile(fname);
		} catch(IOException e) {
			System.out.println(e);
			System.exit(0);
		}
		
		int[] cl = new int[myData.numCols];
		int i = 0;
		int j = 0;
		int numConditions = 0;
		int cond = 0;
		
		boolean doCommon = true;
		
		for (j=0;j<myData.numCols;j++) {
			i = myData.experimentCode[j];
			if (i>numConditions) {
				numConditions = i;
			}
			for (i=0;i<myData.numRows;i++) {
				if (myData.values[i][j] <= bindThresh) {
					cl[j]++;
				}
			}
		}
		
		if (doCommon) {
			numConditions++;
		}
		
		HDPConcParam[] concParams = new HDPConcParam[2+numConditions];
		HDPConcParam gamma_p = new HDPConcParam(1.0,0.1);
		HDPConcParam[] alpha_1 = new HDPConcParam[numConditions];
		HDPConcParam alpha_0 = new HDPConcParam(1.0,1.0);
		
		for (i=0;i<numConditions;i++) {
			alpha_1[i] = new HDPConcParam(1.0,0.1);
		}
		
		concParams[0] = gamma_p;
		for (i=0;i<numConditions;i++) {
			concParams[i+1] = alpha_1[i];
		}
		
		concParams[2+numConditions-1] = alpha_0;
		
		Document doc = null;
		ArrayList<Integer> glist = new ArrayList<Integer>();
		ArrayList<DirichletProcess> DPs = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> alpha_0_DP = new ArrayList<DirichletProcess>();
 		
		DirichletProcess root = new DirichletProcess(null,"root");
		gamma_p.addDP(root);
		DPs.add(root);
		DirichletProcess[] condDP = new DirichletProcess[numConditions];
		DirichletProcess[] condDPInitGroup = new DirichletProcess[numConditions];
		DPGroup[] groups = new DPGroup[numConditions];
		boolean groupDocuments = true;
		
		String s = "";
		for(i=0;i<numConditions;i++) {
			if (doCommon && i == numConditions-1) {
				s = "Common";
			} else {
				s = "Cond#" + (i+1);
			}
			condDP[i] = new DirichletProcess(root,s);
			DPs.add(condDP[i]);
			alpha_1[i].addDP(condDP[i]);
		}
		
		if (groupDocuments) {
			for(i=0;i<numConditions;i++) {
				condDPInitGroup[i] = new DirichletProcess(condDP[i],"");
				DPs.add(condDPInitGroup[i]);
				alpha_0_DP.add(condDPInitGroup[i]);
			//	alpha_1[i].addDP(condDPInitGroup[i]);
			}
		}
		
		DirichletProcess d = null;
		
		if (doCommon) {
			extractCommon(myData,bindThresh,groupDocuments,condDP,condDPInitGroup,alpha_0_DP,DPs);
		} else {
			for (i=0;i<cl.length;i++) {
				if (cl[i] > 0) {
					glist.clear();
					for (j=0;j<myData.numRows;j++) {
						if (myData.values[j][i] <= bindThresh) {
							glist.add(j);
						}
					}
					doc = new Document(myData.experimentNames[i],glist);
					cond = myData.experimentCode[i]-1;
					if (groupDocuments) {
					//	d = new DirichletProcess(condDP[cond],"");
						d = new DirichletProcess(condDPInitGroup[cond],myData.experimentNames[i]);
					} else {
						d = new DirichletProcess(condDP[cond],myData.experimentNames[i]);
					}
					d.addDocument(doc);
					DPs.add(d);
					alpha_0_DP.add(d);
				}
			}
		}
		
		alpha_0.addDPs(alpha_0_DP);
		
		HierDirichletProcess myHierDP = new HierDirichletProcess(100,DPs,concParams,myData.geneNames);
		
		if (groupDocuments) {
			for (i=0;i<numConditions;i++) {
			//	groups[i] = new DPGroup(condDP[i],true);
				groups[i] = new DPGroup(condDP[i],false);
				myHierDP.addDPGroup(groups[i]);
			}
		}
		
		myHierDP.packGenes();
		myHierDP.activateAll();
		
		int iters = 50000;
		int burnin = 5000;
		int targetNumClusters = 157;
		int numConcParamIter = 15;
	//	int iters = 100000;
	//	int burnin = 100000;
		
		String mainOutName = "C:\\CPPDataFiles\\yeast_cross_topic_clusters.txt";
//		myHierDP.iterate(iters,burnin,numConcParamIter,targetNumClusters,mainOutName);
		
		if (groupDocuments) {
			for (i=0;i<numConditions;i++) {
				groups[i].normPairProbs();
			}
		}
		
	//	myHierDP.normPairProbs(iters-burnin);
	//	myHierDP.concensusModel();
		
		String fn = "C:\\CPPDataFiles\\yeast_cross_topic_documents_";
		try {
		//	myHierDP.outputClustersToFile(mainOutName);
			myHierDP.outputNumClustersSamples("C:\\CPPDataFiles\\yeast_cross_topic_numclusters.txt");
			if (groupDocuments) {
				for (i=0;i<numConditions;i++) {
					s = fn + "cond_" + (i+1);
					groups[i].writePairwiseProbs(s);
					s = s + ".txt";
					groups[i].writeConcensusModel(s);
				}
			}
		} catch(IOException e) { System.out.println(e);	}

	}
	
	public static void extractCommon(MicroArrayData myData,double bindThresh,boolean groupDocuments,DirichletProcess[] condDP, DirichletProcess[] condDPInitGroup,ArrayList<DirichletProcess> alpha_0_DP,ArrayList<DirichletProcess> DPs) {
		int i = 0;
		int j = 0;
		int k = 0;
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
					if (myData.values[j][idx] <= bindThresh) {
						uniqueBind[j][0]++;
					}
				}
			} else {
				for (j=0;j<myData.numRows;j++) {
					commonBind[j] = 0;
				}
				for (j=0;j<myData.numRows;j++) {
					for (k=0;k<matchedExperiments.size();k++) {
						idx = matchedExperiments.get(k);
						if (myData.values[j][idx] <= bindThresh) {
							uniqueBind[j][k]++;
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
							glist.add(j);
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
					if (groupDocuments) {
						d = new DirichletProcess(condDPInitGroup[condDPInitGroup.length-1],s);
					} else {
						d = new DirichletProcess(condDP[condDP.length-1],s);
					}
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
						glist.add(j);
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
					if (groupDocuments) {
						d = new DirichletProcess(condDPInitGroup[cond],s);
					} else {
						d = new DirichletProcess(condDP[cond],s);
					}
					d.addDocument(doc);
					DPs.add(d);
					alpha_0_DP.add(d);
				//	System.out.println(s + " " + numBound + " unique genes");
				}
			}
		}
	}
	public static void GenesAsDocuments() {
		String fname = "C:\\CPPDataFiles\\combine_filtered.txt";
		MicroArrayData myData = new MicroArrayData();
		myData.setDiscrete();
		try {
			myData.readFile(fname);
		} catch(IOException e) {
			System.out.println(e);
			System.exit(0);
		}
		int[] gl = new int[myData.numRows];
		int i = 0;
		int j = 0;
		
		for (i=0;i<myData.numRows;i++) {
			for (j=0;j<myData.numCols;j++) {
				if (myData.dvalues[i][j] > 0) {
					gl[i]++;
				}
			}
		}
		
		HDPConcParam[] concParams = new HDPConcParam[2];
		HDPConcParam gamma_p = new HDPConcParam(1.0,0.1);
		HDPConcParam alpha_0 = new HDPConcParam(1.0,1.0);
		
		concParams[0] = gamma_p;
		concParams[1] = alpha_0;
		
		Document doc = null;
		ArrayList<Integer> elist = new ArrayList<Integer>();
		ArrayList<DirichletProcess> DPs = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> alpha_0_DP = new ArrayList<DirichletProcess>();
 		
		DirichletProcess root = new DirichletProcess(null,"root");
		gamma_p.addDP(root);
		DPs.add(root);
		DPGroup[] groups = new DPGroup[1];
		boolean groupDocuments = true;
		
		String s = "";
		
		DirichletProcess d = null;
		int k = 0;
		int v = 0;
		
		for (i=0;i<gl.length;i++) {
			if (gl[i] > 0) {
				elist.clear();
				for (j=0;j<myData.numCols;j++) {
					v = myData.dvalues[i][j];
					if (v > 0) {
						for (k=0;k<v;k++) {
							elist.add(j);
						}
					}
				}
				doc = new Document(myData.geneNames[i],elist);
				d = new DirichletProcess(root,myData.geneNames[i]);
				d.addDocument(doc);
				DPs.add(d);
				alpha_0_DP.add(d);
			}
		}
		
		myData.dvalues = null;
		alpha_0.addDPs(alpha_0_DP);
		HierDirichletProcess myHierDP = new HierDirichletProcess(1,DPs,concParams,myData.experimentNames);
		
		myHierDP.packGenes();
		myHierDP.activateAll();
		
		int iters = 20000;
		int burnin = 10000;
		int targetNumClusters = 0;
		int numConcParamIter = 15;
	//	int iters = 100000;
	//	int burnin = 100000;
		
		String mainOutName = "C:\\CPPDataFiles\\yeast_genedocs";
	//	myHierDP.iterate(iters,burnin,numConcParamIter,targetNumClusters,mainOutName);

		
//		String fn = "C:\\CPPDataFiles\\yeast_cross_topic_documents_";
		try {
			myHierDP.outputClustersToFile(mainOutName);
			myHierDP.outputNumClustersSamples("C:\\CPPDataFiles\\yeast_cross_topic_numclusters.txt");
		} catch(IOException e) { System.out.println(e);	}

	}
}
 