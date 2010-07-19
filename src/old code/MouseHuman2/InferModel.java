package edu.mit.csail.psrg.georg.MouseHuman2;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import cern.jet.random.Uniform;

import edu.mit.csail.psrg.georg.DataAccess.CommandLineParser;
import edu.mit.csail.psrg.georg.DataAccess.CommandLineParserException;
import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.HDP2.DirichletProcess;
import edu.mit.csail.psrg.georg.HDP2.HDP;
import edu.mit.csail.psrg.georg.HDP2.HDPConcParam;
import edu.mit.csail.psrg.georg.HDP2.ParentDP;
import edu.mit.csail.psrg.georg.HDP2.REPCollector;
import edu.mit.csail.psrg.georg.HDP2.TissueDP;
import edu.mit.csail.psrg.georg.HDP2.TissueGroup;
import edu.mit.csail.psrg.georg.HDP2.UpDownConcParam;

public class InferModel {
	
	public static void main(String[] args) {
		doInference(args);
	}
	
	public static CommandLineParser constructParser() {
		CommandLineParser parser = new CommandLineParser();
		try {
			parser.addSwitch("--controlFile","String","file with parameter settings");
			parser.addSwitch("--mainOutFileName","String","base file name for outputing snapshot files");
			parser.addSwitch("--dataFileName","String","file name for data file to use");
			parser.addSwitch("--useGroups","NoArgument","use automatic grouping");
			parser.addSwitch("--iters","Integer","number of iterations to sample");
			parser.addSwitch("--burnin","Integer","number of iterations to burn in the sampler");
			parser.addSwitch("--groupStartIter","Integer","number of iterations before sampling automatic groups");
			parser.addSwitch("--enableSnapShots","NoArgument","enables snapshots (saving of the whole model in binary format)");
			parser.addSwitch("--snapShotInterval","Integer","interval to save snapshots");
			parser.addSwitch("--collectREPs","NoArgument","enable collection of REPs");
			parser.addSwitch("--grpsAlpha_a","Double","alpha_a parameter for group concentration prior");
			parser.addSwitch("--grpsAlpha_b","Double","alpha_b parameter for group concentration prior");
			parser.addSwitch("--init_numTopics","Integer","number of topics to start with");
			parser.addSwitch("--reanimateFileName","String","name of persisted file to load and start iterating on");
		} catch(CommandLineParserException e) {
			System.out.println(e);
		}
		return parser;
	}
	
	public static void doInference(String[] args) {
		boolean useGroups = true;
	//	boolean useGroups = false;
		int groupStartIter = 29999;
	//	int groupStartIter = 0;
		
		double grpsAlpha_a = 10e-8;
		double grpsAlpha_b = 10e-8;
		int init_numTopics = 1;
		
		int iters = 150000;
		int burnin = 100000;
		int numConcParamIter = 15;
		
		int snapShotInterval = 1000;
		boolean enableSnapShots = true;
		boolean collectREPs = true;
		
	//	String mainOutName = "C:\\research_data\\mouse_human\\b47mm7\\geo\\UD\\iter";
	//	String dataFileName = "C:\\research_data\\mouse_human\\b47mm7\\geo\\GNF1_geo_merged_UD_discrete_floor5.txt";
		String mainOutName = "C:\\research_data\\ramaswamy\\GeneProgram\\results\\test\\iter";
		String dataFileName = "C:\\research_data\\ramaswamy\\GeneProgram\\RamaGCRMA_hg18_geo_discrete_floor3.txt";
	//	String dataFileName = "C:\\research_data\\mouse_human\\b47mm7\\geo\\GNF1_discretized_data_geo.txt";
	//	String controlFile = "C:\\research_data\\mouse_human\\homo_b47_data\\tissue_as_docs_control.txt";
		String controlFile = "C:\\research_data\\mouse_human\\homo_b47_data\\tissue_as_docs_control.txt";
		String reanimateFileName = null;
	//	String reanimateFileName = "C:\\research_data\\ramaswamy\\GeneProgram\\results\\test\\iter_snap_2999.persist";
		
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
			
			if (parseMap.containsKey("--init_numTopics")) {
				init_numTopics = (Integer) parseMap.get("--init_numTopics");
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
			
			if (parseMap.containsKey("--enableSnapShots")) {
				enableSnapShots = true;
			} else {
				enableSnapShots = false;
			}
			
			if (parseMap.containsKey("--collectREPs")) {
				collectREPs = true;
			} else {
				collectREPs = false;
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
		}
		*/
		
		HDP myHDP = null;
		
		if (reanimateFileName == null) {
			myHDP = buildModel(useGroups,init_numTopics,dataFileName);
			
			if (useGroups) {
				myHDP.DPGroups.get(0).alpha_a = grpsAlpha_a;
				myHDP.DPGroups.get(0).alpha_b = grpsAlpha_b;
			}
			
			if (collectREPs) {
				myHDP.addStatCollector(new REPCollector(myHDP));
			}
			
			myHDP.activateAll();
		} else {
			myHDP = HDP.restoreFromFile(reanimateFileName);
		}
		
		String snapShotName = mainOutName + "_snap";
		myHDP.groupStartIter = groupStartIter;
		if (enableSnapShots) {
			myHDP.enableSnapShots(1000,snapShotName);
		} else {
			myHDP.snapShotFileName = null;
		}
		myHDP.iterate(iters,burnin,numConcParamIter);	
		System.out.println("Done");
	}
	
	public static HDP buildModel(boolean doGroups,int init_numTopics,String fname) {
		int i = 0;
		int j = 0;
		int k = 0;
		
		HDP myHDP = new HDP();
		
		MicroArrayData myData = new MicroArrayData();
		myData.setDiscrete();
		try {
			myData.readFile(fname);
		} catch(IOException e) {
			System.out.println(e);
			System.exit(0);
		}
		
		int num_initGroups = 3;
		
		HDPConcParam[] concParams = new HDPConcParam[2];
	//	HDPConcParam[] concParams = new HDPConcParam[1];
		HDPConcParam root_conc = new HDPConcParam(10e-8,10e-8);
		HDPConcParam alpha_0 = new HDPConcParam(10e-8,10e-8);
		concParams[0] = root_conc;
		concParams[1] = alpha_0;
	//	concParams[0] = alpha_0;
		
		UpDownConcParam[] UDConcParams = null;
		if (doGroups) {
			UDConcParams = new UpDownConcParam[1];
			UDConcParams[0] = new UpDownConcParam(10e-8,10e-8);
		}
		
		DirichletProcess dp = null;
		TissueDP tissue = null;
		ArrayList<Integer> glist = new ArrayList<Integer>();
		ArrayList<Boolean> udlist = new ArrayList<Boolean>();
		ArrayList<DirichletProcess> DPs = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> alpha_0_DP = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> root_conc_DP = new ArrayList<DirichletProcess>();
		ParentDP root = null;
		
		root = new ParentDP(myHDP,null,"root");
		root_conc_DP.add(root);
	//	alpha_0_DP.add(root);
		DPs.add(root);
		
		root_conc.addDPs(root_conc_DP);
		
		TissueGroup group = null;
		ParentDP DPInitGroup = null;
		
		ParentDP[] DPInitGroups = new ParentDP[num_initGroups];
		
	/*	if (doGroups) {
			DPInitGroup = new ParentDP(myHDP,root,"");
			alpha_0_DP.add(DPInitGroup);
			DPs.add(DPInitGroup);
			ArrayList<ParentDP> pdp = new ArrayList<ParentDP>();
			pdp.add(DPInitGroup);
			UDConcParams[0].addDPs(pdp);
		} */
		
		if (doGroups) {
			ArrayList<ParentDP> pdp = new ArrayList<ParentDP>();
			for (i=0;i<num_initGroups;i++) {
				DPInitGroups[i] = new ParentDP(myHDP,root,"");
				alpha_0_DP.add(DPInitGroups[i]);
				DPs.add(DPInitGroups[i]);
				pdp.add(DPInitGroups[i]);
			}
			UDConcParams[0].addDPs(pdp);
		}
		
		int numTimes = 0;
		
		int chooseParent = 0;
		ParentDP parent = null;
		boolean ud = true;
		for (j=0;j<myData.numCols;j++) {
			glist.clear();
			udlist.clear();
			for (i=0;i<myData.numRows;i++) {
				numTimes = myData.dvalues[i][j];
				if (numTimes < 0) {
					ud = false;
					numTimes = -numTimes;
				} else {
					ud = true;
				}
				if (numTimes > 0) {
					for (k=0;k<numTimes;k++) {
						glist.add(i);
						if (ud) {
							udlist.add(true);
						} else {
							udlist.add(false);
						}
					}
				}
			}
			parent = null;
			if (!doGroups) {
				parent = root;
			} else {
				chooseParent = Uniform.staticNextIntFromTo(0,num_initGroups-1);
				parent = DPInitGroups[chooseParent];
			//	parent = DPInitGroup;
			}
			
			tissue = new TissueDP(myHDP,parent,myData.experimentNames[j],glist,udlist);
			if (myData.experimentCode != null) {
				tissue.code = myData.experimentCode[j];
			}
			tissue.dpID = j;
			
			DPs.add(tissue);
			alpha_0_DP.add(tissue);
		}
		
		alpha_0.addDPs(alpha_0_DP);
		
		myHDP.init(init_numTopics,DPs,concParams,UDConcParams,myData.geneNames);
		
		if (doGroups) {
			group = new TissueGroup(myHDP,root);
			myHDP.addDPGroup(group);
		}
		
		myHDP.packGenes();
		
		return myHDP;
	}
}
