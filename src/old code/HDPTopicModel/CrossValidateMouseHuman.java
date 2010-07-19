package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.StringTokenizer;

import edu.mit.csail.psrg.georg.DataAccess.CommandLineParser;
import edu.mit.csail.psrg.georg.DataAccess.CommandLineParserException;
import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.StatUtil.StatUtil;

public class CrossValidateMouseHuman {
	public static void main(String[] args) {
		createCrossValidationSet();
	//	crossValidateWithGroups(args);
	}
	
	public static CommandLineParser constructParser() {
		CommandLineParser parser = new CommandLineParser();
		
		try {
			parser.addSwitch("--controlFile","String","File with parameter settings");
			parser.addSwitch("--dirPath","String","Directory for saving output");
			parser.addSwitch("--dataFile","String","File containing expression data");
			parser.addSwitch("--crossValidateFile","String","Name of the file that holds cross-validate groups");
			parser.addSwitch("--holdOutBlock","Integer","The number of the block to hold out (0 indexed; -1 to go through all blocks)");
			parser.addSwitch("--useGroups","NoArgument","Use automatic grouping");
			parser.addSwitch("--splitSpecies","NoArgument","Split mouse and human tissues into two groups under the root in the hierarchy");
			parser.addSwitch("--humanOnly","NoArgument","Cluster only the human tissues");
			parser.addSwitch("--initBurnin","Integer","Initial number of iterations to burnin each cross-validation");
			parser.addSwitch("--groupStartIter","Integer","Number of iterations into the burnin to start grouping (only relevant if groups are selected)");
			parser.addSwitch("--sampleSnapshotInterval","Integer","Interval to generate snapshots (e.g., how many iterations between posterior samples)");
			parser.addSwitch("--numPosteriorSamples","Integer","Number of posterior samples to use for each cross-validation");
			parser.addSwitch("--predictBurnin","Integer","Number of iterations to burn-in before capturing prediction likelihoods");
			parser.addSwitch("--numLikelihoodSamples","Integer","Number of predictive likelihood samples to store");
			parser.addSwitch("--grpsAlpha_a","Double","alpha_a parameter for groups gamma prior");
			parser.addSwitch("--grpsAlpha_b","Double","alpha_b parameter for groups gamma prior");
		} catch(CommandLineParserException e) {
			System.out.println(e);
		}
		
		return parser;
	}
	
	public static void createCrossValidationSet() {
		String dirPath = "C:\\research_data\\mouse_human\\b47mm7\\geo\\";
		String fOutName = "cross_validate_all.txt";
		
		boolean uniformSample = true;
		boolean humanOnly = false;
		int numSamples = 140;
		int foldValidate = 10;
		
	//	String fname = "C:\\research_data\\mouse_human\\homo_b47_data\\GNF1_discretized_data.txt";
	//	String fname = "C:\\research_data\\mouse_human\\b47mm7\\geo\\GNF1_discretized_data_geo_common_tissues.txt";
		String fname = "C:\\research_data\\mouse_human\\b47mm7\\geo\\GNF1_discretized_data_geo.txt";
		
		MicroArrayData myData = new MicroArrayData();
		myData.setDiscrete();
		try {
			myData.readFile(fname);
		} catch(IOException e) {
			System.out.println(e);
			System.exit(0);
		}
		
		ArrayList<String> mouseTissues = new ArrayList<String>();
		ArrayList<String> humanTissues = new ArrayList<String>();
		int i = 0;
		int j = 0;
		int k = 0;
		int f = 0;
		
		for (i=0;i<myData.experimentNames.length;i++) {
			if (myData.experimentNames[i].indexOf("h_") == 0) {
				humanTissues.add(myData.experimentNames[i]);
			}
			if (myData.experimentNames[i].indexOf("m_") == 0) {
				mouseTissues.add(myData.experimentNames[i]);
			}
		}
		
		if (uniformSample) {
			humanTissues.addAll(mouseTissues);
		}
		
		ArrayList< ArrayList<String> > holdOuts = new ArrayList< ArrayList <String> >();
		
		ArrayList<String> ho = new ArrayList<String>();
		int[] r_H = StatUtil.randPerm(humanTissues.size());
		int[] r_M = StatUtil.randPerm(mouseTissues.size());
		
		holdOuts.add(ho);
		j = 0;
		k = Math.round(((float) numSamples)/((float) foldValidate));
		for (i=0;i<numSamples;i++) {
			if (j == k) {
				j = 1;
				ho = new ArrayList<String>();
				holdOuts.add(ho);
			} else {
				j++;
			}
			ho.add(humanTissues.get(r_H[i]));
			if (!humanOnly & !uniformSample) {
				ho.add(mouseTissues.get(r_M[i]));
			}
		}
		
		try {
			FileWriter file = new FileWriter(dirPath+fOutName);
			for (i=0;i<holdOuts.size();i++) {
				for (j=0;j<holdOuts.get(i).size();j++) {
					file.write(holdOuts.get(i).get(j));
					if (j<holdOuts.get(i).size()-1) {
						file.write("\t");
					} else {
						file.write("\n");
					}
				}
			}
			file.close();
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void crossValidateWithGroups(String[] args) {
			String dirPath = "C:\\research_data\\mouse_human\\cross_validate_mh_grps_nosplit\\";
			String crossValidateFile = "C:\\research_data\\mouse_human\\cross_validate_tissues.txt";
			String dataFile = "C:\\research_data\\mouse_human\\homo_b47_data\\GNF1_discretized_data.txt";
			String controlFile = null;
			
			int initBurnin = 15000;
			int numConcParamIter = 15;
			int sampleSnapshotInterval = 100;
			int numPosteriorSamples = 10;
			int predictBurnin = 1500;
			int numLikelihoodSamples = 500;
			int groupStartIter = 7500;
			int holdOutBlock = -1;
			
			double grpsAlpha_a = 1.0;
			double grpsAlpha_b = 0.25;
			
			boolean useGroups = true;
			boolean splitSpecies = false;
			boolean humanOnly = false;
			
			try {
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
				
				if (parseMap.containsKey("--dirPath")) {
					dirPath = (String) parseMap.get("--dirPath");
				}
				
				if (parseMap.containsKey("--dataFile")) {
					dataFile = (String) parseMap.get("--dataFile");
				}
				
				if (parseMap.containsKey("--crossValidateFile")) {
					crossValidateFile = (String) parseMap.get("--crossValidateFile");
				}
				
				if (parseMap.containsKey("--initBurnin")) {
					initBurnin = (Integer) parseMap.get("--initBurnin");
				}
				
				if (parseMap.containsKey("--groupStartIter")) {
					groupStartIter = (Integer) parseMap.get("--groupStartIter");
				}
				
				if (parseMap.containsKey("--sampleSnapshotInterval")) {
					sampleSnapshotInterval = (Integer) parseMap.get("--sampleSnapshotInterval");
				}
				
				if (parseMap.containsKey("--numPosteriorSamples")) {
					numPosteriorSamples = (Integer) parseMap.get("--numPosteriorSamples");
				}
				
				if (parseMap.containsKey("--predictBurnin")) {
					predictBurnin = (Integer) parseMap.get("--predictBurnin");
				}
				
				if (parseMap.containsKey("--numLikelihoodSamples")) {
					numLikelihoodSamples = (Integer) parseMap.get("--numLikelihoodSamples");
				}
				
				if (parseMap.containsKey("--holdOutBlock")) {
					holdOutBlock = (Integer) parseMap.get("--holdOutBlock");
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
				
				if (parseMap.containsKey("--grpsAlpha_a")) {
					grpsAlpha_a = (Double) parseMap.get("--grpsAlpha_a");
				}
				
				if (parseMap.containsKey("--grpsAlpha_b")) {
					grpsAlpha_b = (Double) parseMap.get("--grpsAlpha_b");
				}
				
				parser.outputMap(parseMap);
				
			} catch(CommandLineParserException e) {
				System.out.println(e);
				return;
			} catch(IOException e) {
				System.out.println(e);
				return;
			}
			
			ArrayList< ArrayList<String> > holdOuts = readCrossValidateFile(crossValidateFile);
			
			HierDirichletProcess myHierDP = MouseHumanModel.buildHDPWithGroups(humanOnly,useGroups,splitSpecies,1,dataFile);
			
			if (useGroups) {
				myHierDP.DPGroups.get(0).alpha_a = grpsAlpha_a;
				myHierDP.DPGroups.get(0).alpha_b = grpsAlpha_b;
			}
			
			myHierDP.activateAll();
			myHierDP.addStatCollector(new LogLikeCollector(myHierDP));
			
			LogLikeCollector myLogLikeCollector = null; 
			
			if (useGroups){
				myHierDP.groupStartIter = groupStartIter;
			}
	
			myHierDP.persistSelf(dirPath + "cleanHDP.persist");
			
			LinkedHashMap<String,DirichletProcess> tissueMap = null;
			tissueMap = getTissues(myHierDP);
			
			int j = 0;
			DirichletProcess dp = null;
			String mainOutName;
			int k = 0;
			int l = 0;
			int q = 0;
			int n = tissueMap.size();
			int m = 0;
			int mm = 0;
			String s = "";
			double ll = 0.0;
			ArrayList<Double> mySamples = new ArrayList<Double>();
			String tissueName = null;
			ArrayList<Integer> holdOutBlocks = new ArrayList<Integer>();
			if (holdOutBlock == -1) {
				for (j=0;j<holdOuts.size();j++) {
					holdOutBlocks.add(j);
				}
			} else {
				holdOutBlocks.add(holdOutBlock);
			}
			
			for (mm=0;mm<holdOutBlocks.size();mm++) {
				m = holdOutBlocks.get(mm);
				myHierDP = HierDirichletProcess.restoreFromFile(dirPath + "cleanHDP.persist");
				tissueMap = getTissues(myHierDP);
				for (l=0;l<holdOuts.get(m).size();l++) {
					dp = tissueMap.get(holdOuts.get(m).get(l));
					dp.holdout(myHierDP);
				}
				myHierDP.snapShotFileName = null;
				myHierDP.iterate(initBurnin,initBurnin);
				myHierDP.persistSelf(dirPath + "fold_" + (new Integer(m+1)).toString() + ".persist");
				for (l=0;l<holdOuts.get(m).size();l++) {
					myHierDP = HierDirichletProcess.restoreFromFile(dirPath + "fold_" + (new Integer(m+1)).toString() + ".persist");
					tissueMap = getTissues(myHierDP);
					tissueName = holdOuts.get(m).get(l);
					dp = tissueMap.get(tissueName);
					// place the gene count in as the first item, so that perplexity can be calculated
					mySamples.clear();
					mySamples.add((double) dp.documents[0].genes.length);
					mainOutName = dirPath + "cv_sample_" + (new Integer(dp.documents[0].docID+1)).toString();
					myHierDP.enableSnapShots(sampleSnapshotInterval,mainOutName);
					myHierDP.groupStartIter = -1;
					myHierDP.iterate(numPosteriorSamples*sampleSnapshotInterval,numPosteriorSamples*sampleSnapshotInterval);
					k = 0;
					for (j=0;j<numPosteriorSamples;j++) {
						k += sampleSnapshotInterval;
						s = mainOutName + "_" + (new Integer(k-1)).toString() + ".persist";
						myHierDP = HierDirichletProcess.restoreFromFile(s);
						myHierDP.snapShotFileName = null;
						tissueMap = getTissues(myHierDP);
						dp = tissueMap.get(tissueName);
						for (q=0;q<myHierDP.DP.length;q++) {
							if (myHierDP.DP[q].documents != null & myHierDP.DP[q] != dp & myHierDP.DP[q].state != DirichletProcess.HELDOUT) {
								myHierDP.DP[q].freeze();
							}
						}
						dp.activate(myHierDP);
						dp.resampleBeta(myHierDP.numClusters,myHierDP.condLike);
						myHierDP.iterate(numLikelihoodSamples+predictBurnin,predictBurnin,numConcParamIter);
						myLogLikeCollector = (LogLikeCollector) myHierDP.getStatCollector(LogLikeCollector.class);
						for (q=0;q<myLogLikeCollector.logLikeSamples.size();q++) {
							ll = myLogLikeCollector.logLikeSamples.get(q);
							mySamples.add(ll);
						}
					}
					outputDoubleArray(dirPath + "loglike_" + ((new Integer(dp.documents[0].docID+1))).toString() + ".txt",mySamples);
				}
			}
		}
	
	static LinkedHashMap<String,DirichletProcess> getTissues(HierDirichletProcess myHierDP) {
		int i = 0;
		int j = 0;
		LinkedHashMap<String,DirichletProcess> tissueMap = new LinkedHashMap<String,DirichletProcess>();
		
		for (i=0;i<myHierDP.DP.length;i++) {
			if (myHierDP.DP[i].documents != null) {
				tissueMap.put(myHierDP.DP[i].getLabel(),myHierDP.DP[i]);
			}
		} 
		
		return tissueMap;
	}
	
	public static void outputDoubleArray(String fName,ArrayList<Double> x) {
		int i = 0;
		try {
			FileWriter file = new FileWriter(fName);
			for (i=0;i<x.size();i++) {
				file.write(x.get(i).toString());
				file.write("\n");
			}
			file.close();
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static ArrayList< ArrayList<String> > readCrossValidateFile(String fName) {
		ArrayList< ArrayList<String> > holdOut = new ArrayList< ArrayList<String> >();
		ArrayList<String> ho = null;
		
		String[] tissues = null;
		int i = 0;
		try {
			BufferedReader is = new BufferedReader(new FileReader(fName));
			String line = "";
			line = is.readLine();
			while(line != null) {
				tissues = line.split("\t");
				ho = new ArrayList<String>();
				holdOut.add(ho);
				for (i=0;i<tissues.length;i++) {
					ho.add(tissues[i]);
				}
				line = is.readLine();
			}
			is.close();
		} catch(IOException e) {
			System.out.println(e);
		}
		
		return holdOut;
	}
}
