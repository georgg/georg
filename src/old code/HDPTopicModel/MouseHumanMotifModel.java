package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import edu.mit.csail.psrg.georg.DataAccess.CommandLineParser;
import edu.mit.csail.psrg.georg.DataAccess.CommandLineParserException;
import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

public class MouseHumanMotifModel {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		clusterMotifs(args);
	}
	
	public static void clusterMotifs(String[] args) {
		int iters = 150000;
		int burnin = 100000;
		int numConcParamIter = 15;
		int groupStartIter = 29999;
		String controlFile = null;
		int snapShotInterval = 1000;
		boolean enableSnapShots = true;
		boolean collectStrongModules = true;
		int initNumClusters = 25;
		double root_conc_a = 1.0;
		double root_conc_b = 0.2;
		double alpha_1_conc_a = 1.0;
		double alpha_1_conc_b = 0.5;
		double alpha_0_conc_a = 1.0;
		double alpha_0_conc_b = 1.0;
		
		String mainOutName = "C:\\research_data\\mouse_human\\motifs\\genes_as_docs";
		String combinedMotifsFileName = "C:\\research_data\\mouse_human\\motifs\\combined_cons_merged.txt";
		String topGenesFileName = "c:\\research_data\\mouse_human\\homo_b47_clusters_new\\top_genes.txt";
		String reanimateFileName = null;
		
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
			
			if (parseMap.containsKey("--mainOutFileName")) {
				mainOutName = (String) parseMap.get("--mainOutFileName");
			}
			
			if (parseMap.containsKey("--dataFileName")) {
				combinedMotifsFileName = (String) parseMap.get("--dataFileName");
			}
			
			if (parseMap.containsKey("--topGenesFileName")) {
				topGenesFileName = (String) parseMap.get("--topGenesFileName");
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
			
			if (parseMap.containsKey("--initNumClusters")) {
				initNumClusters = (Integer) parseMap.get("--initNumClusters");
			}
			
			if (parseMap.containsKey("--root_conc_a")) {
				root_conc_a = (Double) parseMap.get("--root_conc_a");
			}
			
			if (parseMap.containsKey("--root_conc_b")) {
				root_conc_b = (Double) parseMap.get("--root_conc_b");
			}
			
			if (parseMap.containsKey("--alpha_1_conc_a")) {
				alpha_1_conc_a = (Double) parseMap.get("--alpha_1_conc_a");
			}
			
			if (parseMap.containsKey("--alpha_1_conc_b")) {
				alpha_1_conc_b = (Double) parseMap.get("--alpha_1_conc_b");
			}
			
			if (parseMap.containsKey("--alpha_0_conc_a")) {
				alpha_0_conc_a = (Double) parseMap.get("--alpha_0_conc_a");
			}
			
			if (parseMap.containsKey("--alpha_0_conc_b")) {
				alpha_0_conc_b = (Double) parseMap.get("--alpha_0_conc_b");
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
		
		String snapShotName = mainOutName + "_snap";
		
		double[] concHyperParams = new double[6];
		concHyperParams[0] = root_conc_a;
		concHyperParams[1] = root_conc_b;
		concHyperParams[2] = alpha_1_conc_a;
		concHyperParams[3] = alpha_1_conc_b;
		concHyperParams[4] = alpha_0_conc_a;
		concHyperParams[5] = alpha_0_conc_b;
		HierDirichletProcess myHierDP = null;
		
		if (reanimateFileName == null) {
			myHierDP = constructModel(combinedMotifsFileName,topGenesFileName,initNumClusters,concHyperParams);
			myHierDP.activateAll();
			
			myHierDP.groupStartIter = groupStartIter;
			if (enableSnapShots) {
				myHierDP.enableSnapShots(snapShotInterval,snapShotName);
			} else {
				myHierDP.snapShotFileName = null;
			}
			
			myHierDP.addStatCollector(new NumClustersCollector(myHierDP));
			
			if (collectStrongModules) {
				int i = 0;
				ArrayList<DirichletProcess> clusterDPs = new ArrayList<DirichletProcess>();
				for (i=0;i<myHierDP.DP.length;i++) {
					if (myHierDP.DP[i].label.indexOf("cluster") == 0) {
						clusterDPs.add(myHierDP.DP[i]);
					}
				}
				myHierDP.addStatCollector(new StrongModuleTransposeCollector(myHierDP,clusterDPs));
			}
		} else {
			myHierDP = HierDirichletProcess.restoreFromFile(reanimateFileName);
			
			myHierDP.groupStartIter = groupStartIter;
			if (enableSnapShots) {
				myHierDP.enableSnapShots(snapShotInterval,snapShotName);
			} else {
				myHierDP.snapShotFileName = null;
			}
			
			if (collectStrongModules) {
				int i = 0;
				ArrayList<DirichletProcess> clusterDPs = new ArrayList<DirichletProcess>();
				for (i=0;i<myHierDP.DP.length;i++) {
					if (myHierDP.DP[i].label.indexOf("cluster") == 0) {
						clusterDPs.add(myHierDP.DP[i]);
					}
				}
				myHierDP.addStatCollector(new StrongModuleTransposeCollector(myHierDP,clusterDPs));
			}
		}
		
		myHierDP.iterate(iters,burnin,numConcParamIter);	
		System.out.println("Done");
	}
	
	public static CommandLineParser constructParser() {
		CommandLineParser parser = new CommandLineParser();
		
		try {
			parser.addSwitch("--controlFile","String","file with parameter settings");
			parser.addSwitch("--mainOutFileName","String","base file name for outputing files");
			parser.addSwitch("--dataFileName","String","file name for motif data file to use");
			parser.addSwitch("--topGenesFileName","String","file name for top genes list");
			parser.addSwitch("--iters","Integer","number of iterations to sample");
			parser.addSwitch("--burnin","Integer","number of iterations to burn in the sampler");
			parser.addSwitch("--collectStrongModules","NoArgument","enables collecting strong modules");
			parser.addSwitch("--enableSnapShots","NoArgument","enables snapshots (saving of the whole model in binary format)");
			parser.addSwitch("--snapShotInterval","Integer","interval to save snapshots");
			parser.addSwitch("--initNumClusters","Integer","number of clusters to start with");
			parser.addSwitch("--root_conc_a","Double","concentration parameter a for root");
			parser.addSwitch("--root_conc_b","Double","concentration parameter b for root");
			parser.addSwitch("--alpha_1_conc_a","Double","concentration parameter a for cluster level");
			parser.addSwitch("--alpha_1_conc_b","Double","concentration parameter b for cluster level");
			parser.addSwitch("--alpha_0_conc_a","Double","concentration parameter a for gene level");
			parser.addSwitch("--alpha_0_conc_b","Double","concentration parameter b for gene level");
			parser.addSwitch("--reanimateFileName","String","name of persisted file to load and start iterating on");
		} catch(CommandLineParserException e) {
			System.out.println(e);
		}
		
		return parser;
	}
	
	public static HierDirichletProcess constructModel(String combinedMotifsFileName,String topGenesFileName,int init_numClusters,double[] concHyperParams) {
		ArrayList<ArrayList<String> > modules = null;
		HashMap<String,Integer> geneMap = new HashMap<String,Integer>();
		MicroArrayData motifData = new MicroArrayData();
		motifData.setDiscrete();
		ArrayList<String> genes = null;
		
		try {
			motifData.readFile(combinedMotifsFileName);
			modules = readTopGenes(topGenesFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		int j = 0;
		int k = 0;
		int mm = 0;
		int gg = 0;
		
		for (i=0;i<motifData.numRows;i++) {
			geneMap.put(motifData.geneNames[i],i);
		}
		
		int numMotifs = 0;
		ArrayList<String> keepGenes = new ArrayList<String>();
		for (mm=0;mm<modules.size();mm++) {
			keepGenes.clear();
			genes = modules.get(mm);
			if (genes.size() > 0) {
				for (i=0;i<genes.size();i++) {
					if (geneMap.containsKey(genes.get(i))) {
						numMotifs = 0;
						gg = geneMap.get(genes.get(i));
						for (j=0;j<motifData.numCols;j++) {
							if (motifData.dvalues[gg][j] > 0)
						//	if (motifData.dvalues[gg][j] > 1)
								numMotifs++;
						}
						if (numMotifs > 0)
							keepGenes.add(genes.get(i));
					}
				}
			}
			genes.clear();
			if (keepGenes.size() > 0)
				genes.addAll(keepGenes);
		}
		
		HDPConcParam[] concParams = new HDPConcParam[3];
		
/*		HDPConcParam root_conc = new HDPConcParam(1.0,0.2);
		HDPConcParam alpha_1 = new HDPConcParam(1.0,0.5);
		HDPConcParam alpha_0 = new HDPConcParam(1.0,1.0); */
		
		HDPConcParam root_conc = new HDPConcParam(concHyperParams[0],concHyperParams[1]);
		HDPConcParam alpha_1 = new HDPConcParam(concHyperParams[2],concHyperParams[3]);
		HDPConcParam alpha_0 = new HDPConcParam(concHyperParams[4],concHyperParams[5]);
		
/*		HDPConcParam root_conc = new HDPConcParam(10e-8,10e-8);
		HDPConcParam alpha_1 = new HDPConcParam(10e-8,10e-8);
		HDPConcParam alpha_0 = new HDPConcParam(1.0,1.0); */
		
		concParams[0] = root_conc;
		concParams[1] = alpha_1;
		concParams[2] = alpha_0;
		
		Document doc = null;
		ArrayList<Integer> glist = new ArrayList<Integer>();
		ArrayList<DirichletProcess> DPs = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> alpha_0_DP = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> alpha_1_DP = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> root_conc_DP = new ArrayList<DirichletProcess>();
		DirichletProcess root = null;
		
		root = new DirichletProcess(null,"root");
		root_conc_DP.add(root);
		DPs.add(root);
		root_conc.addDPs(root_conc_DP);
		
		ArrayList<DirichletProcess> moduleDPs = new ArrayList<DirichletProcess>();
		DirichletProcess d = null;
		DirichletProcess parent = null;
		for (mm=0;mm<modules.size();mm++) {
			genes = modules.get(mm);
			if (genes.size() > 0) {
				parent = new DirichletProcess(root,"cluster#"+(new Integer(mm+1)).toString());
				moduleDPs.add(parent);
				DPs.add(parent);
				alpha_1_DP.add(parent);
			} else {
				moduleDPs.add(null);
			}
		}
		alpha_1.addDPs(alpha_1_DP);
		
		int numTimes = 0;
		int docCode = 0;
		for (mm=0;mm<modules.size();mm++) {
			genes = modules.get(mm);
			if (genes.size() > 0) {
				parent = moduleDPs.get(mm);
				for (i=0;i<genes.size();i++) {
					if (geneMap.containsKey(genes.get(i))) {
						gg = geneMap.get(genes.get(i));
						glist.clear();
						for (j=0;j<motifData.numCols;j++) {
							numTimes = motifData.dvalues[gg][j];
							if (numTimes > 0) {
						//	if (numTimes > 1) {
						//		numTimes=2;
								for (k=0;k<numTimes;k++) {
									glist.add(j);
								}
							}
						}
						doc = new Document(genes.get(i) + "_" + (new Integer(mm+1)).toString(),glist);
						doc.code = docCode;
						docCode++;
						d = new DirichletProcess(parent,doc.annotation);
						d.addDocument(doc);
						DPs.add(d);
						alpha_0_DP.add(d);
					}
				}
			}
		}
		alpha_0.addDPs(alpha_0_DP);
		
		HierDirichletProcess myHierDP = new HierDirichletProcess(init_numClusters,DPs,concParams,motifData.experimentNames);
		myHierDP.packGenes();
		
		return myHierDP;
	}
	
	public static ArrayList<ArrayList<String> > readTopGenes(String topGenesFileName) throws IOException {
		ArrayList<ArrayList<String> > modules = new ArrayList<ArrayList<String> >();
		ArrayList<String> genes = null;
		
		int i = 0;
		String[] lineSplit = null;
		BufferedReader is = new BufferedReader(new FileReader(topGenesFileName));
		String line = "";		
		line = is.readLine();
		while(line != null) {
			lineSplit = line.split("\t");
			genes = new ArrayList<String>();
			modules.add(genes);
			if (lineSplit.length > 1) {
				for (i=1;i<lineSplit.length;i++) {
					genes.add(lineSplit[i]);
				}
			}
			line = is.readLine();
		}
		is.close();
		
		return modules;
	}

/*	public static void processMotifs() {
		String humanRatioFileName = "C:\\research_data\\mouse_human\\motifs\\human_2000bp_6ord_ratios.txt";
		String humanRatioFileName2 = "C:\\research_data\\mouse_human\\motifs\\human_2000bp_6ord_ratios2.txt";
		
		MicroArrayData humanRatioData = new MicroArrayData();
		
		try {
			humanRatioData.readFileBlind(humanRatioFileName);
			humanRatioData.writeFile(humanRatioFileName2);
		} catch(IOException e) {
			System.out.println(e);
		}
	} */
}
