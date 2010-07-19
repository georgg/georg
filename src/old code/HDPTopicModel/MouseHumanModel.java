package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;

import edu.mit.csail.psrg.georg.DataAccess.CommandLineParser;
import edu.mit.csail.psrg.georg.DataAccess.CommandLineParserException;
import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.StatUtil.StatUtil;

public class MouseHumanModel {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	//	clusterWithHeldoutMotifs();
	//	crossValidateWithGroups();
		clusterWithGroups(args);
	//	tissuesAsDocs();
	//	reanimate();
	}
	
	public static void reanimate() {
		String dirPath = "C:\\research_data\\mouse_human\\homo_b47_data\\";
		String snapFileName = dirPath + "tissue_as_docs_snap_6999.persist";
		HierDirichletProcess myHDP = HierDirichletProcess.restoreFromFile(snapFileName);
		
		System.out.println("Loaded");
		
		int iters = 3000;
		int burnin = 0;
		int numConcParamIter = 15;
		
		String mainOutName = "C:\\research_data\\mouse_human\\homo_b47_data\\tissue_as_docs";
		String snapShotName = mainOutName + "_snap";
		myHDP.enableSnapShots(1000,snapShotName);
		myHDP.iterate(iters,burnin,numConcParamIter);
	}
	
	public static void tissuesAsDocs() {
		String fname = "C:\\research_data\\mouse_human\\GNF1_discretized_data.txt";
	//	String fname = "C:\\research_data\\mouse_human\\GNF1_merged_PMA.txt";
		
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
		
	//	HDPConcParam[] concParams = new HDPConcParam[4];
		HDPConcParam[] concParams = new HDPConcParam[2];
		HDPConcParam root_conc = new HDPConcParam(1.0,0.2);
	//	HDPConcParam mouse_conc = new HDPConcParam(1.0,0.1);
	//	HDPConcParam human_conc = new HDPConcParam(1.0,0.1);
		HDPConcParam alpha_0 = new HDPConcParam(1.0,1.0);
		concParams[0] = root_conc;
//		concParams[1] = human_conc;
//		concParams[2] = mouse_conc;
//		concParams[3] = alpha_0;
		concParams[1] = alpha_0;
		
		Document doc = null;
		ArrayList<Integer> glist = new ArrayList<Integer>();
		ArrayList<DirichletProcess> DPs = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> alpha_0_DP = new ArrayList<DirichletProcess>();
		
		ArrayList<DirichletProcess> root_conc_DP = new ArrayList<DirichletProcess>();
 		
		DirichletProcess root = new DirichletProcess(null,"root");
	//	root_conc.addDP(root);
		root_conc_DP.add(root);
		DPs.add(root);
		
		DirichletProcess human = new DirichletProcess(root,"human");
	//	human_conc.addDP(human);
		root_conc_DP.add(human);
		DPs.add(human);
		
		DirichletProcess mouse = new DirichletProcess(root,"mouse");
	//	mouse_conc.addDP(mouse);
		root_conc_DP.add(mouse);
		DPs.add(mouse);
		
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
			if (myData.experimentCode[j] == 1) {
				parent = human;
			} else {
				parent = mouse;
			}
			d = new DirichletProcess(parent,myData.experimentNames[j]);
			d.addDocument(doc);
			DPs.add(d);
			alpha_0_DP.add(d);
		}
		
		alpha_0.addDPs(alpha_0_DP);
		
		int init_numClusters = 1;
		HierDirichletProcess myHierDP = new HierDirichletProcess(init_numClusters,DPs,concParams,myData.geneNames);
		
		myHierDP.packGenes();
		myHierDP.activateAll();
		
		int iters = 100000;
		int burnin = 50000;
		int numConcParamIter = 15;
		
		String mainOutName = "C:\\research_data\\mouse_human\\tissue_as_docs";
		String snapShotName = mainOutName + "_snap";
		myHierDP.enableSnapShots(1000,snapShotName);
		myHierDP.iterate(iters,burnin,numConcParamIter);
		
		try {
			myHierDP.outputClustersToFile(mainOutName);
		} catch(IOException e) { System.out.println(e);	}

	}
	
	public static CommandLineParser constructParser() {
		CommandLineParser parser = new CommandLineParser();
		
		try {
			parser.addSwitch("--controlFile","String","file with parameter settings");
			parser.addSwitch("--mainOutFileName","String","base file name for outputing files");
			parser.addSwitch("--dataFileName","String","file name for data file to use");
			parser.addSwitch("--useGroups","NoArgument","use automatic grouping");
			parser.addSwitch("--splitSpecies","NoArgument","split mouse and human tissues into two groups under the root in the hierarchy");
			parser.addSwitch("--humanOnly","NoArgument","cluster only the human tissues");
			parser.addSwitch("--iters","Integer","number of iterations to sample");
			parser.addSwitch("--burnin","Integer","number of iterations to burn in the sampler");
			parser.addSwitch("--groupStartIter","Integer","number of iterations before sampling automatic groups");
			parser.addSwitch("--enableSnapShots","NoArgument","enables snapshots (saving of the whole model in binary format)");
			parser.addSwitch("--snapShotInterval","Integer","interval to save snapshots");
			parser.addSwitch("--collectStrongModules","NoArgument","enable collection of strong modules");
			parser.addSwitch("--collectDistances","NoArgument","enable collection of Hellinger distances between DPs");
			parser.addSwitch("--grpsAlpha_a","Double","alpha_a parameter for group concentration prior");
			parser.addSwitch("--grpsAlpha_b","Double","alpha_b parameter for group concentration prior");
			parser.addSwitch("--init_numClusters","Integer","number of topics to start with");
			parser.addSwitch("--reanimateFileName","String","name of persisted file to load and start iterating on");
		} catch(CommandLineParserException e) {
			System.out.println(e);
		}
		
		return parser;
	}
	
	public static void clusterWithGroups(String[] args) {
	//	boolean useGroups = false;
		boolean useGroups = true;
		boolean splitSpecies = false;
		boolean humanOnly = false;
		boolean useMotifs = false;
		boolean splitMotifs = true;
		int groupStartIter = 29999;
		
		double grpsAlpha_a = 10e-8;
		double grpsAlpha_b = 10e-8;
		int init_numClusters = 1;
		
		int iters = 150000;
		int burnin = 100000;
		int numConcParamIter = 15;
		
		int snapShotInterval = 1000;
		boolean enableSnapShots = true;
		boolean collectStrongModules = true;
		boolean collectDistances = true;
		
	//	String mainOutName = "C:\\research_data\\mouse_human\\homo_b47_data\\tissue_as_docs";
		String mainOutName = "C:\\research_data\\ramaswamy\\GeneProgram\\results\\test\\iter";
		
	//	String mainOutName = "C:\\research_data\\mouse_human\\b47mm7\\geo\\UD\\iter";
		String dataFileName = "C:\\research_data\\mouse_human\\homo_b47_data\\GNF1_discretized_data.txt";
	//	String dataFileName = "C:\\research_data\\mouse_human\\b47mm7\\geo\\GNF1_discretized_data_geo.txt";
		String controlFile = "C:\\research_data\\mouse_human\\homo_b47_data\\tissue_as_docs_control.txt";
	//	String reanimateFileName = null;
		String reanimateFileName = "C:\\research_data\\ramaswamy\\GeneProgram\\results\\test\\iter_snap_3999.persist";
		
	/*	try {
			CommandLineParser parser = constructParser();
			if (args[0].equals("--help")) {
				parser.displayDescriptions();
				return;
			}
			HashMap<String,Object> parseMap = parser.parseCommandLine(args);
			
			if (parseMap.containsKey("--controlFile")) {
				controlFile = (String) parseMap.get("--controlFile");
			}
			
			if (controlFile != null) {
				parseMap = parser.parseCommandFile(controlFile);
				HashMap<String,Object> parseMap2 = parser.parseCommandLine(args);
				parseMap.putAll(parseMap2);
			}
			
			if (parseMap.containsKey("--mainOutFileName")) {
				mainOutName = (String) parseMap.get("--mainOutFileName");
			}
			
			if (parseMap.containsKey("--dataFileName")) {
				dataFileName = (String) parseMap.get("--dataFileName");
			}
			
			if (parseMap.containsKey("--groupStartIter")) {
				groupStartIter = (Integer) parseMap.get("--groupStartIter");
			}
			
			if (parseMap.containsKey("--init_numClusters")) {
				init_numClusters = (Integer) parseMap.get("--init_numClusters");
			}
			
			if (parseMap.containsKey("--snapShotInterval")) {
				snapShotInterval = (Integer) parseMap.get("--snapShotInterval");
			}
			
			if (parseMap.containsKey("--iters")) {
				iters = (Integer) parseMap.get("--iters");
			}
			
			if (parseMap.containsKey("--burnin")) {
				burnin = (Integer) parseMap.get("--burnin");
			}
			
			if (parseMap.containsKey("--useGroups")) {
				useGroups = true;
			} else {
				useGroups = false;
			}
			
			if (parseMap.containsKey("--splitSpecies")) {
				splitSpecies = true;
			} else {
				splitSpecies = false;
			}
			
			if (parseMap.containsKey("--humanOnly")) {
				humanOnly = true;
			} else {
				humanOnly = false;
			}
			
			if (parseMap.containsKey("--enableSnapShots")) {
				enableSnapShots = true;
			} else {
				enableSnapShots = false;
			}
			
			if (parseMap.containsKey("--collectStrongModules")) {
				collectStrongModules = true;
			} else {
				collectStrongModules = false;
			}
			
			if (parseMap.containsKey("--collectDistances")) {
				collectDistances = true;
			} else {
				collectDistances = false;
			}
			
			if (parseMap.containsKey("--grpsAlpha_a")) {
		
				grpsAlpha_a = (Double) parseMap.get("--grpsAlpha_a");
			}
			
			if (parseMap.containsKey("--grpsAlpha_b")) {
				grpsAlpha_b = (Double) parseMap.get("--grpsAlpha_b");
			}
			
			if (parseMap.containsKey("--reanimateFileName")) {
				reanimateFileName = (String) parseMap.get("--reanimateFileName");
			}
			
			parser.outputMap(parseMap);
			
		} catch(CommandLineParserException e) {
			System.out.println(e);
			return;
		} catch(IOException e) {
			System.out.println(e);
			return;
		} */
		
		
		HierDirichletProcess myHierDP = null;
		
		if (reanimateFileName == null) {
			if (!useMotifs) {
				myHierDP = buildHDPWithGroups(humanOnly,useGroups,splitSpecies,init_numClusters,dataFileName);
			} else {
				myHierDP = buildHDPWithMotifs(useGroups,splitMotifs);
			}
			
			if (useGroups) {
				myHierDP.DPGroups.get(0).alpha_a = grpsAlpha_a;
				myHierDP.DPGroups.get(0).alpha_b = grpsAlpha_b;
			}
			
			if (collectStrongModules) {
				myHierDP.addStatCollector(new StrongModuleCollector(myHierDP));
			}
			
			if (collectDistances) {
				int i = 0;
				ArrayList<DirichletProcess> DPs = new ArrayList<DirichletProcess>();
				for (i=0;i<myHierDP.DP.length;i++) {
					if (myHierDP.DP[i].documents != null) {
						DPs.add(myHierDP.DP[i]);
					}
				}
				myHierDP.addStatCollector(new DPDistanceCollector(myHierDP,DPs));
			}
			myHierDP.activateAll();
		} else {
			myHierDP = HierDirichletProcess.restoreFromFile(reanimateFileName);
		}
		
	//	String mainOutName = "C:\\research_data\\mouse_human\\tissue_as_docs";
		String snapShotName = mainOutName + "_snap";
		myHierDP.groupStartIter = groupStartIter;
		if (enableSnapShots) {
			myHierDP.enableSnapShots(1000,snapShotName);
		} else {
			myHierDP.snapShotFileName = null;
		}
		myHierDP.iterate(iters,burnin,numConcParamIter);	
		System.out.println("Done");
	}
	
	public static HierDirichletProcess buildHDPWithGroups(boolean humanOnly, boolean doGroups,boolean splitSpecies) {
//		String fname = "C:\\research_data\\mouse_human\\GNF1_discretized_data.txt";
		String fname = "C:\\research_data\\mouse_human\\homo_b47_data\\GNF1_discretized_data.txt";
		
		return buildHDPWithGroups(humanOnly,doGroups,splitSpecies,1,fname);
	}
		
	public static HierDirichletProcess buildHDPWithGroups(boolean humanOnly, boolean doGroups,boolean splitSpecies,int init_numClusters,String fname) {
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
			
		//	HDPConcParam root_conc = new HDPConcParam(1.0,0.2);
		//	HDPConcParam alpha_0 = new HDPConcParam(1.0,1.0);
			
			HDPConcParam root_conc = new HDPConcParam(10e-8,10e-8);
			HDPConcParam alpha_0 = new HDPConcParam(10e-8,10e-8);
			
			concParams[0] = root_conc;
			concParams[1] = alpha_0;
			
			Document doc = null;
			ArrayList<Integer> glist = new ArrayList<Integer>();
			ArrayList<DirichletProcess> DPs = new ArrayList<DirichletProcess>();
			ArrayList<DirichletProcess> alpha_0_DP = new ArrayList<DirichletProcess>();
			ArrayList<DirichletProcess> root_conc_DP = new ArrayList<DirichletProcess>();
			DirichletProcess root = null;
			DirichletProcess human = null;
			DirichletProcess mouse = null;
			
			root = new DirichletProcess(null,"root");
			root_conc_DP.add(root);
			DPs.add(root);
			
			if (!humanOnly & splitSpecies) {
				human = new DirichletProcess(root,"human");
				root_conc_DP.add(human);
				DPs.add(human);
				
				mouse = new DirichletProcess(root,"mouse");
				root_conc_DP.add(mouse);
				DPs.add(mouse);
			}
			
			root_conc.addDPs(root_conc_DP);
			
			DPGroup[] groups = null;
			if (doGroups) {
				if (humanOnly | !splitSpecies) {
					groups = new DPGroup[1];
				} else {
					groups = new DPGroup[2];
				}
			}
			
			DirichletProcess humanDPInitGroup = null;
			DirichletProcess mouseDPInitGroup = null;
			
			if (doGroups & !groupDocuments) {
				if (!humanOnly & splitSpecies) {
					mouseDPInitGroup = new DirichletProcess(mouse,"");
					alpha_0_DP.add(mouseDPInitGroup);
					DPs.add(mouseDPInitGroup);
					humanDPInitGroup = new DirichletProcess(human,"");
				} else {
					humanDPInitGroup = new DirichletProcess(root,"");
				}
				alpha_0_DP.add(humanDPInitGroup);
				DPs.add(humanDPInitGroup);
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
				doc.code = myData.experimentCode[j];
				doc.docID = j;
				parent = null;
				if (myData.experimentCode[j] == 1) {
					if (doGroups) {
						if (!groupDocuments) {
							parent = humanDPInitGroup;
						} else {
							parent = human;
						}
					} else {
						if (splitSpecies) {
							parent = human;
						} else {
							parent = root;
						}
					}
				} else {
					if (doGroups) {
						if (splitSpecies) {
							if (!groupDocuments) {
								parent = mouseDPInitGroup;
							} else {
								parent = mouse;
							}
						} else {
							if (!groupDocuments) {
								parent = humanDPInitGroup;
							} else {
								parent = human;
							}
						}
					} else {
						if (splitSpecies) {
							parent = mouse;
						} else {
							parent = root;
						}
					}
				}
				if (!humanOnly | (myData.experimentCode[j] == 1)) {
					d = new DirichletProcess(parent,myData.experimentNames[j]);
					d.code = myData.experimentCode[j];
					d.addDocument(doc);
					DPs.add(d);
					alpha_0_DP.add(d);
				}
			}
			
			alpha_0.addDPs(alpha_0_DP);
			
			HierDirichletProcess myHierDP = new HierDirichletProcess(init_numClusters,DPs,concParams,myData.geneNames);
			
			if (doGroups) {
				if (!humanOnly & splitSpecies) {
					groups[0] = new DPGroup(human,groupDocuments);
				} else {
					groups[0] = new DPGroup(root,groupDocuments);
				}
				myHierDP.addDPGroup(groups[0]);
				if (!humanOnly & splitSpecies) {
					groups[1] = new DPGroup(mouse,groupDocuments);
					myHierDP.addDPGroup(groups[1]);
				}
			}
			
			myHierDP.packGenes();
			
			return myHierDP;
			
		/*	try {
				myHierDP.outputClustersToFile(mainOutName);
				myHierDP.outputNumClustersSamples("C:\\research_data\\mouse_human\\numclusters.txt");
			} catch(IOException e) { System.out.println(e);	} */
	}
	
	static void getMotifs(HierDirichletProcess myHierDP,LinkedHashMap<String,DirichletProcess> motifs,LinkedHashMap<String,DirichletProcess> tissues) {
		int i = 0;
		int j = 0;
		
		for (i=0;i<myHierDP.DP.length;i++) {
			if (myHierDP.DP[i].documents != null) {
				j = myHierDP.DP[i].getLabel().indexOf("_mot_M");
				if (j >= 0) {
					motifs.put(myHierDP.DP[i].getLabel(),myHierDP.DP[i]);
				} else {
					tissues.put(myHierDP.DP[i].getLabel(),myHierDP.DP[i]);
				}
			}
		}
	}
	
	public static HierDirichletProcess buildHDPWithMotifs(boolean doGroups,boolean splitMotifs) {
		String fname = "C:\\research_data\\mouse_human\\homo_b47_data\\GNF1_discretized_data.txt";	
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
				
		HDPConcParam[] concParams = null;
		if (splitMotifs) {
			concParams = new HDPConcParam[5];
		} else {
			concParams = new HDPConcParam[2];
		}
				
				
		HDPConcParam root_conc = new HDPConcParam(1.0,0.2);
		HDPConcParam tissue_root_conc = new HDPConcParam(1.0,0.2);
		HDPConcParam motif_root_conc = new HDPConcParam(1.0,0.2);
		HDPConcParam tissue_conc = new HDPConcParam(1.0,0.5);
		HDPConcParam motif_conc = new HDPConcParam(1.0,0.5);
		concParams[0] = root_conc;
		concParams[1] = tissue_root_conc;
		if (splitMotifs) {
			concParams[2] = motif_root_conc;
			concParams[3] = tissue_conc;
			concParams[4] = motif_conc;
		}
				
		Document doc = null;
		ArrayList<Integer> glist = new ArrayList<Integer>();
		ArrayList<DirichletProcess> DPs = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> tissue_DP = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> motif_DP = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> tissue_root_DP = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> motif_root_DP = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> root_DP = new ArrayList<DirichletProcess>();
		DirichletProcess root = null;
		DirichletProcess tissues = null;
		DirichletProcess motifs = null;
				
		root = new DirichletProcess(null,"root");
		root_DP.add(root);
		DPs.add(root);
		
		if (splitMotifs) {
			tissues = new DirichletProcess(root,"tissues");
			tissue_root_DP.add(tissues);
			DPs.add(tissues);
			
			motifs = new DirichletProcess(root,"motifs");
			motif_root_DP.add(motifs);
			DPs.add(motifs);
		}
		
		root_conc.addDPs(root_DP);
		
		if (splitMotifs) {
			tissue_root_conc.addDPs(tissue_root_DP);
			motif_root_conc.addDPs(motif_root_DP);
		}
				
		DPGroup[] groups = null;
		if (doGroups) {
			if (splitMotifs) {
				groups = new DPGroup[2];
			} else {
				groups = new DPGroup[1];
			}
		}
			
		DirichletProcess tissueDPInitGroup = null;
		DirichletProcess motifDPInitGroup = null;
		
		if (doGroups) {
			if (splitMotifs) {
				tissueDPInitGroup = new DirichletProcess(tissues,"");
				tissue_DP.add(tissueDPInitGroup);
				DPs.add(tissueDPInitGroup);
				motifDPInitGroup = new DirichletProcess(motifs,"");
				motif_DP.add(motifDPInitGroup);
				DPs.add(motifDPInitGroup);
			} else {
				tissueDPInitGroup = new DirichletProcess(root,"");
				tissue_DP.add(tissueDPInitGroup);
				DPs.add(tissueDPInitGroup);
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
				if (splitMotifs) {
					if (myData.experimentCode[j] > 2) {
						parent = motifDPInitGroup;
					} else {
						parent = tissueDPInitGroup;
					}
				} else {
					parent = tissueDPInitGroup;
				}
			} else {
				if (splitMotifs) {
					if (myData.experimentCode[j] > 2) {
						parent = motifs;
					} else {
						parent = tissues;
					}
				} else {
					parent = root;
				}
			}
			
			d = new DirichletProcess(parent,myData.experimentNames[j]);
			d.addDocument(doc);
			DPs.add(d);
			if (splitMotifs) {
				if (myData.experimentCode[j] > 2) {
					motif_DP.add(d);
				} else {
					tissue_DP.add(d);
				}
			} else {
				tissue_DP.add(d);
			}
		}
		
		tissue_conc.addDPs(tissue_DP);
		if (splitMotifs) {
			motif_conc.addDPs(motif_DP);
		}
				
		int init_numClusters = 1;
		HierDirichletProcess myHierDP = new HierDirichletProcess(init_numClusters,DPs,concParams,myData.geneNames);
		
		if (doGroups) {
			if (splitMotifs) {
				groups[0] = new DPGroup(tissues,false);
				myHierDP.addDPGroup(groups[0]);
				groups[1] = new DPGroup(motifs,false);
				myHierDP.addDPGroup(groups[1]);
			} else {
				groups[0] = new DPGroup(root,false);
				myHierDP.addDPGroup(groups[0]);
			}
		}
		
		myHierDP.packGenes();
		return myHierDP;
		}
	
	public static void clusterWithHeldoutMotifs() {
		boolean useGroups = false;
		boolean splitMotifs = true;
		
		
		HierDirichletProcess myHierDP = null;
		
		myHierDP = buildHDPWithMotifs(useGroups,splitMotifs);
		myHierDP.activateAll();
		LinkedHashMap<String,DirichletProcess> motifs = new LinkedHashMap<String,DirichletProcess>();
		LinkedHashMap<String,DirichletProcess> tissues = new LinkedHashMap<String,DirichletProcess>();
		getMotifs(myHierDP,motifs,tissues);
		Iterator<DirichletProcess> iter = motifs.values().iterator();
		DirichletProcess dp = null;
		while (iter.hasNext()) {
			dp = iter.next();
			dp.holdout(myHierDP);
		}
		
		int initialIters = 60000;
		int motifBurnin = 10000;
		int motifSamples = 15000;
		int reIter = 1000;
		int motifRounds = 10;
		int numConcParamIter = 15;
		int groupStartIter = 49999;
		
	/*	int initialIters = 5000;
		int motifBurnin = 1000;
		int motifSamples = 2000;
		int reIter = 1000;
		int motifRounds = 3;
		int numConcParamIter = 15;
		int groupStartIter = 49999; */
		
		String mainOutName = "C:\\research_data\\mouse_human\\homo_b47_data\\motif_cluster";
		String snapShotName = mainOutName + "_snap";
	//	String restoreName = snapShotName + "_" + (new Integer(initialIters-1)).toString() + ".persist";
		myHierDP.groupStartIter = groupStartIter;
		myHierDP.enableSnapShots(1000,snapShotName);
	//	myHierDP.setCollectStrongModules(false);
		myHierDP.iterate(initialIters,initialIters+1000,numConcParamIter);
		
		int i = 0;
		double priorPos = myHierDP.priorPos;
		for (i=0;i<motifRounds;i++) {
		/*	myHierDP = HierDirichletProcess.restoreFromFile(restoreName);
			motifs.clear();
			tissues.clear();
			getMotifs(myHierDP,motifs,tissues); */
			
			myHierDP.priorPos = priorPos;
			myHierDP.snapShotFileName = null;
	//		myHierDP.setCollectStrongModules(false);
			myHierDP.iterate(reIter,reIter+1000,numConcParamIter);
			
			iter = tissues.values().iterator();
			while(iter.hasNext()) {
				dp = iter.next();
				dp.freeze();
			}
			iter = motifs.values().iterator();
			while(iter.hasNext()) {
				dp = iter.next();
				dp.activate(myHierDP);
				dp.generateNewClusters = false;
			}
			myHierDP.snapShotFileName = snapShotName + "_round" + (new Integer(i+1)).toString();
	//		myHierDP.setCollectStrongModules(true);
			myHierDP.iterate(motifSamples,motifBurnin,numConcParamIter);
			
			iter = tissues.values().iterator();
			while(iter.hasNext()) {
				dp = iter.next();
				dp.state = DirichletProcess.ACTIVE;
			}
			iter = motifs.values().iterator();
			while(iter.hasNext()) {
				dp = iter.next();
				dp.holdout(myHierDP);
			}
		}
		
	}
}
