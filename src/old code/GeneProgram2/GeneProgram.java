package edu.mit.csail.psrg.georg.GeneProgram2;

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

public class GeneProgram {
	
	public static void main(String[] args) {
		doInference(args);
	}
	
	public static CommandLineParser constructParser() {
		CommandLineParser parser = new CommandLineParser();
		try {
			parser.addSwitch("--controlFile","String","file with parameter settings");
			parser.addSwitch("--mainOutFileName","String","base file name for outputing snapshot files");
			parser.addSwitch("--hierarchyFileName","String","file specifying model hierarchy");
			parser.addSwitch("--modifierFileName","String","file specifying base file name for modifier data");
			parser.addSwitch("--dataFileName","String","file name for data file to use");
			parser.addSwitch("--useGroups","NoArgument","use automatic grouping");
			parser.addSwitch("--noUpDown","NoArgument","do not model up and down regulation (just induction)");
			parser.addSwitch("--iters","Integer","number of iterations to sample");
			parser.addSwitch("--burnin","Integer","number of iterations to burn in the sampler");
			parser.addSwitch("--groupStartIter","Integer","number of iterations before sampling automatic groups");
			parser.addSwitch("--enableSnapShots","NoArgument","enables snapshots (saving of the whole model in binary format)");
			parser.addSwitch("--snapShotInterval","Integer","interval to save snapshots");
			parser.addSwitch("--grpsAlpha_a","Double","alpha_a parameter for group concentration prior");
			parser.addSwitch("--grpsAlpha_b","Double","alpha_b parameter for group concentration prior");
			parser.addSwitch("--conc_a","Double","alpha_a parameter for hierarchy concentration prior");
			parser.addSwitch("--conc_b","Double","alpha_b parameter for hierarchy concentration prior");
			parser.addSwitch("--init_numPrograms","Integer","number of expression program to start with");
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
		int init_numPrograms = 1;
		int init_numGroups = 1;
		boolean useUpDown = true;
		
		int iters = 150000;
		int burnin = 100000;
		int numConcParamIter = 50;
		
		int snapShotInterval = 1000;
		boolean enableSnapShots = true;
		
		String mainOutName = "";
		String dataFileName = "";
		String controlFile = "";
		String hierarchyFileName = null;
		String reanimateFileName = null;
		String modifierFileName = null;
		
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
			
			if (parseMap.containsKey("--modifierFileName")) {
				modifierFileName = (String) parseMap.get("--modifierFileName");
			}
			
			if (parseMap.containsKey("--groupStartIter")) {
				groupStartIter = (Integer) parseMap.get("--groupStartIter");
			}
			
			if (parseMap.containsKey("--init_numPrograms")) {
				init_numPrograms = (Integer) parseMap.get("--init_numPrograms");
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
			
			if (parseMap.containsKey("--noUpDown")) {
				useUpDown = false;
			}
			
			if (parseMap.containsKey("--enableSnapShots")) {
				enableSnapShots = true;
			} else {
				enableSnapShots = false;
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
				myHDP = buildModel(useGroups,useUpDown,init_numPrograms,init_numGroups,dataFileName,modifierFileName,conc_a,conc_b);
			} else {
				myHDP = buildModel(useUpDown,init_numPrograms,dataFileName,hierarchyFileName,modifierFileName,conc_a,conc_b);
			}
			
			if (useGroups) {
				myHDP.DPGroups.get(0).alpha_a = grpsAlpha_a;
				myHDP.DPGroups.get(0).alpha_b = grpsAlpha_b;
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
	
	public static HDP buildModel(boolean doGroups,boolean useUpDown,int init_numTopics,int init_numGroups,String fname,String modifierFileName,double conc_a,double conc_b) {
		int i = 0;
		int j = 0;
		int k = 0;
		
		double alphaa = conc_a;
		double alphab = conc_b;
		
		HDP myHDP = new HDP();
		
		myHDP.useUpDown = useUpDown;
		
		MicroArrayData myData = new MicroArrayData();
		myData.setDiscrete();
		try {
			myData.readFile(fname);
		} catch(IOException e) {
			System.out.println(e);
			System.exit(0);
		}
		
		int[][][] modifiers = null;
		int[] modifierLevels = null;
		
		int[][][][] modifiersP = new int[1][][][];
		int[][] modifierLevelsP = new int[1][];
		
		if (modifierFileName != null) {
			loadModifierData(modifierFileName,modifierLevelsP,modifiersP);
			modifiers = modifiersP[0];
			modifierLevels = modifierLevelsP[0];
		}
	
		HDPConcParam[] concParams = new HDPConcParam[2];
		HDPConcParam root_conc = new HDPConcParam(alphaa,alphab);
		HDPConcParam alpha_0 = new HDPConcParam(alphaa,alphab);
		concParams[0] = root_conc;
		concParams[1] = alpha_0;
		
		UpDownConcParam[] UDConcParams = null;
		ModifierConcParam[] ModConcParams = null;
		if (doGroups) {
			if (useUpDown) {
				UDConcParams = new UpDownConcParam[1];
				UDConcParams[0] = new UpDownConcParam(conc_a,conc_b);
			}
			
			if (modifierLevels != null) {
				ModConcParams = new ModifierConcParam[1];
				ModConcParams[0] = new ModifierConcParam(conc_a,conc_b,modifierLevels.length);
			}
		}
		
		DirichletProcess dp = null;
		TissueDP tissue = null;
		ArrayList<Integer> glist = new ArrayList<Integer>();
		ArrayList<DirichletProcess> DPs = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> alpha_0_DP = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> root_conc_DP = new ArrayList<DirichletProcess>();
		ParentDP root = null;
		
		ArrayList<Boolean> udlist = null;
		if (useUpDown) {
			udlist = new ArrayList<Boolean>();
		}
		
		ArrayList<ArrayList<Integer> > modlist = null;
		if (modifierLevels != null) {
			modlist = new ArrayList<ArrayList<Integer> >();
			for (i=0;i<modifierLevels.length;i++) {
				modlist.add(new ArrayList<Integer>());
			}
		}
		
		root = new ParentDP(myHDP,null,"root");
		root_conc_DP.add(root);
		DPs.add(root);
		
		TissueGroup group = null;
		
		ParentDP DPInitGroup = null;
		ParentDP[] DPInitGroups = new ParentDP[init_numGroups];
			
		if (doGroups) {
			ArrayList<ParentDP> pdp = new ArrayList<ParentDP>();
			for (i=0;i<init_numGroups;i++) {
				DPInitGroups[i] = new ParentDP(myHDP,root,"");
				root_conc_DP.add(DPInitGroups[i]);
				DPs.add(DPInitGroups[i]);
				pdp.add(DPInitGroups[i]);
			}
			
			if (UDConcParams != null) {
				UDConcParams[0].addDPs(pdp);
			}
			
			if (ModConcParams != null) {
				ModConcParams[0].addDPs(pdp);
			}
		} 
		
		root_conc.addDPs(root_conc_DP);
		
		int numTimes = 0;
		int ml = 0;
		
		ParentDP parent = null;
		boolean ud = true;
		int chooseParent = 0;
		for (j=0;j<myData.numCols;j++) {
			glist.clear();
			
			if (udlist != null)
				udlist.clear();
			
			if (modifierLevels != null) {
				for (ml=0;ml<modifierLevels.length;ml++) {
					modlist.get(ml).clear();
				}
			}
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
						
						if (udlist != null) {
							if (ud) {
								udlist.add(true);
							} else {
								udlist.add(false);
							}
						}
						
						if (modifierLevels != null) {
							for (ml=0;ml<modifierLevels.length;ml++) {
								modlist.get(ml).add(modifiers[ml][i][j]);
							}
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
			}
			
			tissue = new TissueDP(myHDP,parent,myData.experimentNames[j],glist,udlist,modlist);
			if (myData.experimentCode != null) {
				tissue.code = myData.experimentCode[j];
			}
			tissue.dpID = j;
			
			DPs.add(tissue);
			alpha_0_DP.add(tissue);
		}
		
		alpha_0.addDPs(alpha_0_DP);
		
		myHDP.init(init_numTopics,DPs,concParams,UDConcParams,ModConcParams,myData.geneNames,modifierLevels);
		
		if (doGroups) {
			group = new TissueGroup(myHDP,root);
			myHDP.addDPGroup(group);
		}
		
		myHDP.packGenes();
		
		return myHDP;
	}
	
	public static HDP buildModel(boolean useUpDown,int init_numTopics,String fname,String hierarchyFileName,String modifierFileName,double conc_a,double conc_b) {
		int i = 0;
		int j = 0;
		int k = 0;
		
		double alphaa = conc_a;
		double alphab = conc_b;
		
		HDP myHDP = new HDP();
		
		myHDP.useUpDown = useUpDown;
		
		ClusterReader reader = new ClusterReader();
		try {
			reader.readFile(hierarchyFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int[][][] modifiers = null;
		int[] modifierLevels = null;
		
		int[][][][] modifiersP = new int[1][][][];
		int[][] modifierLevelsP = new int[1][];
		
		if (modifierFileName != null) {
			loadModifierData(modifierFileName,modifierLevelsP,modifiersP);
			modifiers = modifiersP[0];
			modifierLevels = modifierLevelsP[0];
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
		if (labels.size() > 0 & useUpDown) {
			UDconcParams = new UpDownConcParam[labels.size()];
		}
		LinkedHashMap<String,UpDownConcParam> UDconcParamMap = new LinkedHashMap<String,UpDownConcParam>();
		LinkedHashMap<UpDownConcParam,ArrayList<ParentDP> > UDconcParamLists = new LinkedHashMap<UpDownConcParam,ArrayList<ParentDP> >();
		
		ModifierConcParam[] ModConcParams = null;
		if (labels.size() > 0 & modifierLevels != null) {
			ModConcParams = new ModifierConcParam[labels.size()];
		}
		LinkedHashMap<String,ModifierConcParam> ModifierParamMap = new LinkedHashMap<String,ModifierConcParam>();
		LinkedHashMap<ModifierConcParam,ArrayList<ParentDP> > ModConcParamLists = new LinkedHashMap<ModifierConcParam,ArrayList<ParentDP> >();
		
		if (labels.size() > 0) {
			iter = labels.iterator();
			i = 0;
			while(iter.hasNext()) {
				tempS = iter.next();
				
				if (useUpDown) {
					UDconcParams[i] = new UpDownConcParam(alphaa,alphab);
					UDconcParamMap.put(tempS,UDconcParams[i]);
					UDconcParamLists.put(UDconcParams[i],new ArrayList<ParentDP>());
				}
				
				if (modifierLevels != null) {
					ModConcParams[i] = new ModifierConcParam(alphaa,alphab,modifierLevels.length);
					ModifierParamMap.put(tempS,ModConcParams[i]);
					ModConcParamLists.put(ModConcParams[i],new ArrayList<ParentDP>());
				}
				
				i++;
			}
		}
		
		DirichletProcess dp = null;
		TissueDP tissue = null;
		ArrayList<Integer> glist = new ArrayList<Integer>();
		ArrayList<DirichletProcess> DPs = new ArrayList<DirichletProcess>();
		
		ArrayList<Boolean> udlist = null;
		if (useUpDown) {
			udlist =  new ArrayList<Boolean>();
		}
		
		ArrayList<ArrayList<Integer> > modlist = null;
		if (modifierLevels != null) {
			modlist = new ArrayList<ArrayList<Integer> >();
			for (i=0;i<modifierLevels.length;i++) {
				modlist.add(new ArrayList<Integer>());
			}
		}
		
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
		ModifierConcParam mcp = null;
		ParentDP parent = null;
		ArrayList<DirichletProcess> myCPList = null;
		ArrayList<ParentDP> myUDCPList = null;
		ArrayList<ParentDP> myModCPList = null;
		String dpLabel = null;
		int dd = 0;
		boolean ud = true;
		int numTimes = 0;
		int ml = 0;
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
			if (!tempS.equals("NA") & !useUpDown) {
				udcp = UDconcParamMap.get(tempS);
			}
			mcp = null;
			if (modifierLevels != null & !tempS.equals("NA")) {
				mcp = ModifierParamMap.get(tempS);
			}
			
			if (dataMap.containsKey(dpLabel)) {
				j = dataMap.get(dpLabel);
				glist.clear();
				
				if (useUpDown) {
					udlist.clear();
				}
				
				if (modifierLevels != null) {
					for (ml=0;ml<modifierLevels.length;ml++) {
						modlist.get(ml).clear();
					}
				}
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
							
							if (useUpDown) {
								if (ud) {
									udlist.add(true);
								} else {
									udlist.add(false);
								}
							}
							
							if (modifierLevels != null) {
								for (ml=0;ml<modifierLevels.length;ml++) {
									modlist.get(ml).add(modifiers[ml][i][j]);
								}
							}
						}
					}
				}
				tissue = new TissueDP(myHDP,parent,myData.experimentNames[j],glist,udlist,modlist);
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
				if (mcp != null) {
					myModCPList = ModConcParamLists.get(mcp);
					myModCPList.add((ParentDP) dp);
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
		
		if (ModConcParams != null) {
			for (i=0;i<ModConcParams.length;i++) {
				myModCPList = ModConcParamLists.get(ModConcParams[i]);
				ModConcParams[i].addDPs(myModCPList);
			}
		}
		
		myHDP.init(init_numTopics,DPs,concParams,UDconcParams,ModConcParams,myData.geneNames,modifierLevels);
		
		myHDP.packGenes();
		
		return myHDP;
	}
	
	public static void loadModifierData(String modifierBaseFileName,int[][] modifierLevelsP,int[][][][] modifiersP) {
		String structureFileName = modifierBaseFileName + "_mod_names.txt";
		String dataBaseName = modifierBaseFileName + "_mod";
		
		ClusterReader struct = new ClusterReader();
		
		try {
			struct.readFile(structureFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		
		int numModifiers = struct.clusters.size();
		int[] modifierLevels = new int[numModifiers];
		for (i=0;i<numModifiers;i++) {
			modifierLevels[i] = struct.clusters.get(i).size();
		}
		
		int[][][] modifiers = new int[numModifiers][][];
		
		String fname;
		MicroArrayData data = null;
		try {
			for (i=0;i<numModifiers;i++) {
				fname = dataBaseName + (new Integer(i+1)).toString() + ".txt";
				data = new MicroArrayData();
				data.setDiscrete();
				data.readFile(fname);
				modifiers[i] = data.dvalues;
			}
		} catch(IOException e) {
			System.out.println(e);
		}
		
		modifierLevelsP[0] = modifierLevels;
		modifiersP[0] = modifiers;
	}
}
