package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.IOException;
import java.util.ArrayList;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.GO.AffyGoAssocReader;
import edu.mit.csail.psrg.georg.GO.GoOBOReader;
import edu.mit.csail.psrg.georg.GO.MouseHumanGoAssocReader;

public class MacklisModel2 {
	public static void main(String[] args) {
	//	cluster();
		outputTopics();
	}
	
	public static void cluster() {
		boolean useGroups = false;
		boolean splitCellTypes = true;
		
		HierDirichletProcess myHierDP = buildHDPWithGroups(useGroups,splitCellTypes);
	//	HierDirichletProcess myHierDP = buildHDPGenesAsDocs();
		myHierDP.activateAll();
		
		int iters = 100000;
		int burnin = 50000;
		int numConcParamIter = 15;
		
		String mainOutName = "C:\\research_data\\macklis\\macklis";
		String snapShotName = mainOutName + "_snap";
		myHierDP.enableSnapShots(2500,snapShotName);
		myHierDP.setCollectStrongModules(false);
		myHierDP.iterate(iters,burnin,numConcParamIter,false);	
	}
	
	public static HierDirichletProcess buildHDPWithGroups(boolean doGroups,boolean splitCellTypes) {
		String fname = "C:\\research_data\\macklis\\macklis_discretized_data.txt";
			
			int i = 0;
			int j = 0;
			int k = 0;
			boolean groupDocuments = false;
			
			MicroArrayData myData = new MicroArrayData();
			myData.setDiscrete();
			try {
				myData.readFile(fname);
			} catch(IOException e) {
				System.out.println(e);
				System.exit(0);
			}
			
			HDPConcParam[] concParams = new HDPConcParam[2];
			
			HDPConcParam root_conc = new HDPConcParam(1.0,1.0);
			HDPConcParam alpha_0 = new HDPConcParam(0.2,2.0);
			concParams[0] = root_conc;
			concParams[1] = alpha_0;
			
			Document doc = null;
			ArrayList<Integer> glist = new ArrayList<Integer>();
			ArrayList<DirichletProcess> DPs = new ArrayList<DirichletProcess>();
			ArrayList<DirichletProcess> alpha_0_DP = new ArrayList<DirichletProcess>();
			ArrayList<DirichletProcess> root_conc_DP = new ArrayList<DirichletProcess>();
			DirichletProcess root = null;
			DirichletProcess[] cellTypeDPs = new DirichletProcess[3];
			
			root = new DirichletProcess(null,"root");
			root_conc_DP.add(root);
			DPs.add(root);
			
			if (splitCellTypes) {
				cellTypeDPs[0] = new DirichletProcess(root,"CPN");
				root_conc_DP.add(cellTypeDPs[0]);
				DPs.add(cellTypeDPs[0]);
				
				cellTypeDPs[1] = new DirichletProcess(root,"CSPN");
				root_conc_DP.add(cellTypeDPs[1]);
				DPs.add(cellTypeDPs[1]);
				
				cellTypeDPs[2] = new DirichletProcess(root,"CTPN");
				root_conc_DP.add(cellTypeDPs[2]);
				DPs.add(cellTypeDPs[2]);
			}
			
			root_conc.addDPs(root_conc_DP);
			
			DPGroup[] groups = null;
			DirichletProcess[] DPInitGroup = null;
			if (doGroups) {
				if (splitCellTypes) {
					groups = new DPGroup[3];
					DPInitGroup = new DirichletProcess[3];
				}
				else {
					groups = new DPGroup[1];
					DPInitGroup = new DirichletProcess[1];
				}
			}
			
			if (doGroups) {
				if (splitCellTypes) {
					for (i=0;i<3;i++) {
						DPInitGroup[i] = new DirichletProcess(cellTypeDPs[i],"");
						alpha_0_DP.add(DPInitGroup[i]);
						DPs.add(DPInitGroup[i]);
					}
				} else {
					DPInitGroup[0] = new DirichletProcess(root,"");
					alpha_0_DP.add(DPInitGroup[0]);
					DPs.add(DPInitGroup[0]);
				}
			}
			
			int numTimes = 0;
			
			DirichletProcess d = null;
			DirichletProcess parent = null;
			for (j=0;j<myData.numCols;j++) {
				glist.clear();
				for (i=0;i<myData.numRows;i++) {
					numTimes = myData.dvalues[i][j];
					if (numTimes > 0) {
						for (k=0;k<numTimes;k++) {
							glist.add(i);
						}
					}
				}
				doc = new Document(myData.experimentNames[j],glist);
				parent = null;
				if (doGroups) {
					if (splitCellTypes) {
						parent = DPInitGroup[myData.experimentCode[j]-1];
					} else {
						parent = DPInitGroup[0];
					}
				} else {
					if (splitCellTypes) {
						parent = cellTypeDPs[myData.experimentCode[j]-1];
					} else {
						parent = root;
					}
						
				}
				
				d = new DirichletProcess(parent,myData.experimentNames[j]);
				d.addDocument(doc);
				DPs.add(d);
				alpha_0_DP.add(d);
			}
			
			alpha_0.addDPs(alpha_0_DP);
			
			int init_numClusters = 1;
			HierDirichletProcess myHierDP = new HierDirichletProcess(init_numClusters,DPs,concParams,myData.geneNames);
			
			if (doGroups) {
				if (splitCellTypes) {
					for (i=0;i<3;i++) {
						groups[i] = new DPGroup(cellTypeDPs[i],groupDocuments);
						myHierDP.addDPGroup(groups[i]);
					}
				} else {
					groups[0] = new DPGroup(root,groupDocuments);
					myHierDP.addDPGroup(groups[0]);
				}
			}
			
			myHierDP.packGenes();
			
			return myHierDP;
	}
	
	public static void transposeData(MicroArrayData myData) {
		int i = 0;
		int j = 0;
		// transpose the matrix
		int[][] values = new int[myData.numCols][myData.numRows];
		for (i=0;i<myData.numRows;i++) {
			for (j=0;j<myData.numCols;j++) {
				values[j][i] = myData.dvalues[i][j];
			}
		}
		myData.dvalues = values;
		myData.numRows = values.length;
		myData.numCols = values[0].length;
		String[] geneNames = new String[myData.experimentNames.length];
		for (i=0;i<geneNames.length;i++) {
			geneNames[i] = myData.experimentNames[i];
		}
		String[] experimentNames = new String[myData.geneNames.length];
		for (i=0;i<experimentNames.length;i++) {
			experimentNames[i] = myData.geneNames[i];
		}
		myData.geneNames = geneNames;
		myData.experimentNames = experimentNames;
	}
	
	public static HierDirichletProcess buildHDPGenesAsDocs() {
		String fname = "C:\\research_data\\macklis\\macklis_discretized_data.txt";
			
			int i = 0;
			int j = 0;
			int k = 0;
			
			MicroArrayData myData = new MicroArrayData();
			myData.setDiscrete();
			try {
				myData.readFile(fname);
			} catch(IOException e) {
				System.out.println(e);
				System.exit(0);
			}
			
			transposeData(myData);
			
			HDPConcParam[] concParams = new HDPConcParam[2];
			
			HDPConcParam root_conc = new HDPConcParam(1.0,0.1);
			HDPConcParam alpha_0 = new HDPConcParam(1.0,5.0);
			concParams[0] = root_conc;
			concParams[1] = alpha_0;
			
			Document doc = null;
			ArrayList<Integer> glist = new ArrayList<Integer>();
			ArrayList<DirichletProcess> DPs = new ArrayList<DirichletProcess>();
			ArrayList<DirichletProcess> alpha_0_DP = new ArrayList<DirichletProcess>();
			ArrayList<DirichletProcess> root_conc_DP = new ArrayList<DirichletProcess>();
			DirichletProcess root = null;
			
			root = new DirichletProcess(null,"root");
			root_conc_DP.add(root);
			DPs.add(root);
			root_conc.addDPs(root_conc_DP);
			
			int numTimes = 0;
			
			DirichletProcess d = null;
			DirichletProcess parent = null;
			for (j=0;j<myData.numCols;j++) {
				glist.clear();
				for (i=0;i<myData.numRows;i++) {
					numTimes = myData.dvalues[i][j];
					if (numTimes > 0) {
						for (k=0;k<numTimes;k++) {
							glist.add(i);
						}
					}
				}
				doc = new Document(myData.experimentNames[j],glist);
				parent = root;
				d = new DirichletProcess(parent,myData.experimentNames[j]);
				d.addDocument(doc);
				DPs.add(d);
				alpha_0_DP.add(d);
			}
			
			alpha_0.addDPs(alpha_0_DP);
			
			int init_numClusters = 1;
			HierDirichletProcess myHierDP = new HierDirichletProcess(init_numClusters,DPs,concParams,myData.geneNames);
			myHierDP.packGenes();
			
			return myHierDP;
	}
	
	public static void outputTopics() {
		String dirPath = "C:\\research_data\\macklis\\";
		String topicOutFileL = dirPath + "macklis_go_long.txt";
		String topicOutFileS = dirPath + "macklis_go_short.txt";
		String normDPsInFile = dirPath + "macklis_dps.txt";
		String normDPsOutFile = dirPath + "macklis_dps_norm_rows.txt";
		String normDPsOutFile2 = dirPath + "macklis_dps_norm_cols.txt";
		String GMLOutFile = dirPath + "macklis.svg";
		String mouseAssocFile = "C:\\research_data\\mouse_human\\Mm.genelist.go.plus";
		
		String snapFileName = dirPath + "macklis_snap_99999.persist";
		HierDirichletProcess myDP = HierDirichletProcess.restoreFromFile(snapFileName);
		try {
			myDP.outputClustersToFile(dirPath + "macklis");
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		int j = 0;
		
		double norm = 0.0;
		MicroArrayData myData = new MicroArrayData();
		try {
			myData.readFile(normDPsInFile);
			
			for (i=0;i<myData.numRows;i++) {
				norm = 0.0;
				for (j=0;j<myData.numCols;j++) {
					norm = norm + myData.values[i][j];
				}
				for (j=0;j<myData.numCols;j++) {
					myData.values[i][j] = myData.values[i][j]/norm;
				}
			}
			myData.writeFile(normDPsOutFile);
			
			myData = new MicroArrayData();
			myData.readFile(normDPsInFile);
			for (j=0;j<myData.numCols;j++) {
				norm = 0.0;
				for (i=3;i<myData.numRows;i++) {
					norm = norm + myData.values[i][j];
				}
				for (i=3;i<myData.numRows;i++) {
					myData.values[i][j] = myData.values[i][j]/norm;
				}
			}
			myData.writeFile(normDPsOutFile2);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		GoOBOReader goReader = null;
		AffyGoAssocReader assocReader = null;
	//	MouseHumanGoAssocReader assocReader = null;
		StrongModuleSummarizer topicSumm = null;
		
		// output annotated topics
		goReader = new GoOBOReader();
	//	assocReader = new MouseHumanGoAssocReader();
	//	assocReader.setFile(mouseAssocFile);
		assocReader = new AffyGoAssocReader();
		assocReader.setUseREFSEQ(true);
		reconstructGroups(myDP);
		topicSumm = new StrongModuleSummarizer(myDP,goReader,assocReader);
	//	topicSumm.outputTopics(topicOutFileL,topicOutFileS,GMLOutFile,false);
	}
	
	public static void reconstructGroups(HierDirichletProcess myDP) {
		if (myDP.DPGroups.size() == 0)
			return;
		
		int i = 0;
		for (i=0;i<myDP.DPGroups.size();i++) {
			myDP.DPGroups.get(i).normPairProbs();
			myDP.DPGroups.get(i).buildConcensusModel();
		}
		myDP.buildDPList();
		System.out.println("reconstructed groups");
	}
}
