package edu.mit.csail.psrg.georg.HDP2;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;

import cern.jet.random.Uniform;

import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;
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

public class InferBasicModel {
	
	public static void main(String[] args) {
		doInference(args);
	}
	
	public static CommandLineParser constructParser() {
		CommandLineParser parser = new CommandLineParser();
		try {
			parser.addSwitch("--controlFile","String","file with parameter settings");
			parser.addSwitch("--mainOutFileName","String","base file name for outputing snapshot files");
			parser.addSwitch("--hierarchyFileName","String","file specifying model hierarchy");
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
			parser.addSwitch("--conc_a","Double","alpha_a parameter for hierarchy concentration prior");
			parser.addSwitch("--conc_b","Double","alpha_b parameter for hierarchy concentration prior");
			parser.addSwitch("--init_numTopics","Integer","number of topics to start with");
			parser.addSwitch("--init_numGroups","Integer","number of groups to start with");
			parser.addSwitch("--reanimateFileName","String","name of persisted file to load and start iterating on");
		} catch(CommandLineParserException e) {
			System.out.println(e);
		}
		return parser;
	}
	
	public static void doInference(String[] args) {
		boolean useGroups = true;
		int groupStartIter = 29999;
		
		double grpsAlpha_a = 10e-8;
		double grpsAlpha_b = 10e-8;
		double conc_a = 10e-8;
		double conc_b = 10e-8;
		int init_numTopics = 1;
		int init_numGroups = 1;
		
		int iters = 150000;
		int burnin = 100000;
		int numConcParamIter = 50;
		
		int snapShotInterval = 1000;
		boolean enableSnapShots = true;
		boolean collectREPs = true;
		
		String mainOutName = "";
		String dataFileName = "";
	//	String controlFile = "";
		String controlFile = "c:\\research_data\\mouse_human\\simulated\\cluster_control_sample.txt";
	//	String controlFile = "c:\\research_data\\golub\\GeneProgram\\amplitude\\combine_fold1_75\\cluster_golub_combine_fold1_75.txt";
		String hierarchyFileName = null;
		String reanimateFileName = null;
	//	String reanimateFileName = "C:\\research_data\\Bild\\Results\\bild_breast_log2_hier\\iter_snap_99999.persist";
		
		try {
			CommandLineParser parser = constructParser();
			if (args.length > 0) {
				if (args[0].equals("--help")) {
					parser.displayDescriptions();
					return;
				}
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
			
			if (parseMap.containsKey("--hierarchyFileName")) {
				hierarchyFileName = (String) parseMap.get("--hierarchyFileName");
			}
			
			if (parseMap.containsKey("--groupStartIter")) {
				groupStartIter = (Integer) parseMap.get("--groupStartIter");
			}
			
			if (parseMap.containsKey("--init_numTopics")) {
				init_numTopics = (Integer) parseMap.get("--init_numTopics");
			}
			
			if (parseMap.containsKey("--init_numGroups")) {
				init_numGroups = (Integer) parseMap.get("--init_numGroups");
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
			
			if (parseMap.containsKey("--conc_a")) {
				conc_a = (Double) parseMap.get("--conc_a");
			}
			
			if (parseMap.containsKey("--conc_b")) {
				conc_b = (Double) parseMap.get("--conc_b");
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
		
		HDP myHDP = null;
		
		if (reanimateFileName == null) {
			if (hierarchyFileName == null) {
				myHDP = buildModel(useGroups,init_numTopics,init_numGroups,dataFileName,conc_a,conc_b);
			} else {
				myHDP = buildModel(init_numTopics,dataFileName,hierarchyFileName,conc_a,conc_b);
			}
			
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
			myHDP.enableSnapShots(snapShotInterval,snapShotName);
		} else {
			myHDP.snapShotFileName = null;
		}
		myHDP.iterate(iters,burnin,numConcParamIter);	
		System.out.println("Done");
	}
	
	public static HDP buildModel(boolean doGroups,int init_numTopics,int init_numGroups,String fname,double conc_a,double conc_b) {
		int i = 0;
		int j = 0;
		int k = 0;
		
	//	double alphaa = 10e-8;
	//	double alphab = 10e-8;
		
	//	double alphaa = 1;
	//	double alphab = 10;
	//	double alphab = 1;
		double alphaa = conc_a;
		double alphab = conc_b;
		
		HDP myHDP = new HDP();
		
		MicroArrayData myData = new MicroArrayData();
		myData.setDiscrete();
		try {
			myData.readFile(fname);
		} catch(IOException e) {
			System.out.println(e);
			System.exit(0);
		}
		
		HDPConcParam[] concParams = new HDPConcParam[2];
		HDPConcParam root_conc = new HDPConcParam(alphaa,alphab);
		HDPConcParam alpha_0 = new HDPConcParam(alphaa,alphab);
		concParams[0] = root_conc;
		concParams[1] = alpha_0;
		
		UpDownConcParam[] UDConcParams = null;
		if (doGroups) {
			UDConcParams = new UpDownConcParam[1];
		//	UDConcParams[0] = new UpDownConcParam(10e-8,10e-8);
			UDConcParams[0] = new UpDownConcParam(conc_a,conc_b);
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
		DPs.add(root);
		
		TissueGroup group = null;
		
		ParentDP DPInitGroup = null;
	/*	if (doGroups) {
			DPInitGroup = new ParentDP(myHDP,root,"");
			alpha_0_DP.add(DPInitGroup);
			DPs.add(DPInitGroup);
			ArrayList<ParentDP> pdp = new ArrayList<ParentDP>();
			pdp.add(DPInitGroup);
			UDConcParams[0].addDPs(pdp);
		} */
		
		ParentDP[] DPInitGroups = new ParentDP[init_numGroups];
			
		if (doGroups) {
			ArrayList<ParentDP> pdp = new ArrayList<ParentDP>();
			for (i=0;i<init_numGroups;i++) {
				DPInitGroups[i] = new ParentDP(myHDP,root,"");
			//	alpha_0_DP.add(DPInitGroups[i]);
				root_conc_DP.add(DPInitGroups[i]);
				DPs.add(DPInitGroups[i]);
				pdp.add(DPInitGroups[i]);
			}
			UDConcParams[0].addDPs(pdp);
		} 
		
		root_conc.addDPs(root_conc_DP);
		
		int numTimes = 0;
		
		ParentDP parent = null;
		boolean ud = true;
		int chooseParent = 0;
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
				chooseParent = Uniform.staticNextIntFromTo(0,init_numGroups-1);
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
	
	public static HDP buildModel(int init_numTopics,String fname,String hierarchyFileName,double conc_a,double conc_b) {
		int i = 0;
		int j = 0;
		int k = 0;
		
	//	double alphaa = 10e-8;
	//	double alphab = 10e-8;
		
	//	double alphaa = 1;
	//	double alphab = 10;
		
		double alphaa = conc_a;
		double alphab = conc_b;
		
		HDP myHDP = new HDP();
		
		ClusterReader reader = new ClusterReader();
		try {
			reader.readFile(hierarchyFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		LinkedHashSet<String> labels = new LinkedHashSet<String>();
		String tempS = null;
		for (i=0;i<reader.clusters.size();i++) {
			tempS = reader.clusters.get(i).get(2);
			if (!tempS.equals("parent")) {
				labels.add(tempS);
			}
		}
		
		HDPConcParam[] concParams = new HDPConcParam[labels.size()];
		LinkedHashMap<String,HDPConcParam> concParamMap = new LinkedHashMap<String,HDPConcParam>();
		LinkedHashMap<HDPConcParam,ArrayList<DirichletProcess> > concParamLists = new LinkedHashMap<HDPConcParam,ArrayList<DirichletProcess> >();
		Iterator<String> iter = labels.iterator();
		
		i = 0;
		while(iter.hasNext()) {
			tempS = iter.next();
			concParams[i] = new HDPConcParam(alphaa,alphab);
			concParamMap.put(tempS,concParams[i]);
			concParamLists.put(concParams[i],new ArrayList<DirichletProcess>());
			i++;
		}
		
		labels.clear();
		for (i=0;i<reader.clusters.size();i++) {
			tempS = reader.clusters.get(i).get(3);
			if (!tempS.equals("NA")) {
				labels.add(tempS);
			}
		}
		
		UpDownConcParam[] UDconcParams = null;
		if (labels.size() > 0) {
			UDconcParams = new UpDownConcParam[labels.size()];
		}
		LinkedHashMap<String,UpDownConcParam> UDconcParamMap = new LinkedHashMap<String,UpDownConcParam>();
		LinkedHashMap<UpDownConcParam,ArrayList<ParentDP> > UDconcParamLists = new LinkedHashMap<UpDownConcParam,ArrayList<ParentDP> >();
		
		if (labels.size() > 0) {
			iter = labels.iterator();
			i = 0;
			while(iter.hasNext()) {
				tempS = iter.next();
				UDconcParams[i] = new UpDownConcParam(alphaa,alphab);
				UDconcParamMap.put(tempS,UDconcParams[i]);
				UDconcParamLists.put(UDconcParams[i],new ArrayList<ParentDP>());
				i++;
			}
		}
		
		DirichletProcess dp = null;
		TissueDP tissue = null;
		ArrayList<Integer> glist = new ArrayList<Integer>();
		ArrayList<Boolean> udlist = new ArrayList<Boolean>();
		ArrayList<DirichletProcess> DPs = new ArrayList<DirichletProcess>();
		
		MicroArrayData myData = new MicroArrayData();
		myData.setDiscrete();
		try {
			myData.readFile(fname);
		} catch(IOException e) {
			System.out.println(e);
			System.exit(0);
		}
		
		LinkedHashMap<String,Integer> dataMap = new LinkedHashMap<String,Integer>();
		LinkedHashMap<String,DirichletProcess> dpMap = new LinkedHashMap<String,DirichletProcess>();
		for (j=0;j<myData.numCols;j++) {
			dataMap.put(myData.experimentNames[j],j);
		}
		
		HDPConcParam hcp = null;
		UpDownConcParam udcp = null;
		ParentDP parent = null;
		ArrayList<DirichletProcess> myCPList = null;
		ArrayList<ParentDP> myUDCPList = null;
		String dpLabel = null;
		int dd = 0;
		boolean ud = true;
		int numTimes = 0;
		for (dd=0;dd<reader.clusters.size();dd++) {
			dp = null;
			dpLabel = reader.clusters.get(dd).get(0);
			
			// get parent
			parent = null;
			tempS = reader.clusters.get(dd).get(1);
			if (!tempS.equals("NA")) {
				parent = (ParentDP) dpMap.get(tempS);
			}
			
			// get concentration parameter
			tempS = reader.clusters.get(dd).get(2);
			hcp = concParamMap.get(tempS);
			
			// get up/down concentration parameter
			tempS = reader.clusters.get(dd).get(3);
			udcp = null;
			if (!tempS.equals("NA")) {
				udcp = UDconcParamMap.get(tempS);
			}
			
			if (dataMap.containsKey(dpLabel)) {
				j = dataMap.get(dpLabel);
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
				tissue = new TissueDP(myHDP,parent,myData.experimentNames[j],glist,udlist);
				if (myData.experimentCode != null) {
					tissue.code = myData.experimentCode[j];
				}
				tissue.dpID = j;
				dp = tissue;
			} else {
				dp = new ParentDP(myHDP,parent,dpLabel);
			}
			
			if (dp != null) {
				DPs.add(dp);
				dpMap.put(dpLabel,dp);
				
				myCPList = concParamLists.get(hcp);
				myCPList.add(dp);
				
				if (udcp != null) {
					myUDCPList = UDconcParamLists.get(udcp);
					myUDCPList.add((ParentDP) dp);
				}
			}
		}
		
		for (i=0;i<concParams.length;i++) {
			myCPList = concParamLists.get(concParams[i]);
			concParams[i].addDPs(myCPList);
		}
		
		if (UDconcParams != null) {
			for (i=0;i<UDconcParams.length;i++) {
				myUDCPList = UDconcParamLists.get(UDconcParams[i]);
				UDconcParams[i].addDPs(myUDCPList);
			}
		}
		
		myHDP.init(init_numTopics,DPs,concParams,UDconcParams,myData.geneNames);
		
		myHDP.packGenes();
		
		return myHDP;
	}
}
