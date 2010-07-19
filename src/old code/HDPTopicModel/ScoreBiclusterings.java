package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Vector;

import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;
import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.GO.DAGException;
import edu.mit.csail.psrg.georg.GO.GoOBOReader;
import edu.mit.csail.psrg.georg.GO.MouseHumanGoAssocReader;
import edu.mit.csail.psrg.georg.GO.OntClusters;
import edu.mit.csail.psrg.georg.GO.OntDAG;
import edu.mit.csail.psrg.georg.GO.OntGeneAssociations;
import edu.mit.csail.psrg.georg.GO.OntScore;
import edu.mit.csail.psrg.georg.GO.OntTerm;
import edu.mit.csail.psrg.georg.GO.OntTermScores;
import edu.mit.csail.psrg.georg.StatUtil.StatUtil;
import edu.mit.csail.psrg.georg.StatUtil.StepDownFDR;
import edu.mit.csail.psrg.georg.StatUtil.VectorUtil;

public class ScoreBiclusterings {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	//	outputGOCategoriesSAMBA();
	//	convertBicatFile();
		scoreData();
	}
	
	public static void scoreData() {
	/*	String dirPath = "C:\\research_data\\BicAT_v2.22\\";
		String tissuesFileName = dirPath + "ISA_results_filtered_tissues";
		String genesFileName = dirPath + "ISA_results_filtered_genes"; */
		
	/*	String dirPath = "C:\\research_data\\mouse_human\\b47mm7\\geo\\alt_clusterings\\cluster_rownorm3\\";
		String tissuesFileName = dirPath + "HDP_tissues";
		String genesFileName = dirPath + "HDP_genes"; */
		
	//	String dirPath = "C:\\research_data\\expander\\";
	//	String dirPath = "C:\\research_data\\mouse_human\\b47mm7\\geo\\alt_clusterings\\cluster_MI1_nogrps\\";
	//	String dirPath = "C:\\research_data\\mouse_human\\b47mm7\\geo\\";
		String dirPath = "C:\\research_data\\mouse_human\\bioNMF\\factors_109_nmf\\";
	//	String scoreFilePrefix = "HDP";
		String scoreFilePrefix = "nmf_109";
	//	String scoreFilePrefix = "biclusters_mouse_human_coverall";
		
		scoreTissues(dirPath+scoreFilePrefix+"_tissues.txt",dirPath+scoreFilePrefix+"_tissues_score.txt");
		scoreGenes(dirPath + scoreFilePrefix + "_genes");
		biclusterStats(dirPath,scoreFilePrefix,true);
	}
	
	public static void biclusterStats(String dirPath,String scoreFilePrefix,boolean topicEntropy) {
		String genesFile = dirPath + scoreFilePrefix + "_genes.txt";
		String tissuesFile = dirPath + scoreFilePrefix + "_tissues.txt";
		String genesScoreFile = dirPath + scoreFilePrefix + "_genes_score.txt";
		String tissueLoadFile = dirPath + "tissue_load.txt";
		String outFile = dirPath + "topic_stats.txt";
		String groupOutFile = dirPath + "groups_stats.txt";
		
		ClusterReader tissueClusters = new ClusterReader();
		ClusterReader geneClusters = new ClusterReader();
		ClusterReader genesScoreClusters = new ClusterReader();
		MicroArrayData tissueLoad = new MicroArrayData();
		
		String tissueClassFile = "C:\\research_data\\mouse_human\\b47mm7\\tissue_classifications.txt";
		ClusterReader tissueClasses = new ClusterReader();
		
		try {
			tissueClasses.readFile(tissueClassFile);
			tissueClusters.readFile(tissuesFile);
			geneClusters.readFile(genesFile);
			genesScoreClusters.readFile(genesScoreFile);
			if (topicEntropy) {
				tissueLoad.readFile(tissueLoadFile);
			}
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int numGroups = 0;
		double[] e = null;
		double en = 0.0;
		int i = 0;
		int j = 0;
		int k = 0;
		int gn = 0;
		double[] entropies = null;
		HashMap<String,Integer> groupMap = new HashMap<String,Integer>();
		int[] signifGO = new int[genesScoreClusters.clusters.size()];
		
		ArrayList<String> cluster = null;
		
		for(i=0;i<genesScoreClusters.clusters.size();i++) {
			signifGO[i] = (new Integer(genesScoreClusters.clusters.get(i).get(2))).intValue();
		}
		
		if (topicEntropy) {
			for (i=0;i<tissueLoad.numRows;i++) {
				j = (int) tissueLoad.values[i][0];
				groupMap.put(tissueLoad.geneNames[i],j);
				if (j > numGroups) {
					numGroups = j;
				}
			}
			numGroups++;
			
			e = new double[numGroups];
			entropies = new double[tissueClusters.clusters.size()];
			
			for (j=1;j<tissueLoad.numCols;j++) {
				for (k=0;k<numGroups;k++) {
					e[k] = 0.0;
				}
				for (i=0;i<tissueLoad.numRows;i++) {
					if (tissueLoad.values[i][j] > 0.0) {
						gn = groupMap.get(tissueLoad.geneNames[i]);
						e[gn] += tissueLoad.values[i][j];
					}
				}
				en = 0.0;
				for (k=0;k<numGroups;k++) {
					en += e[k];
				}
				for (k=0;k<numGroups;k++) {
					if (e[k] > 0.0) {
						e[k] = e[k]/en;
						entropies[j-1] += -e[k]*Math.log(e[k])/Math.log(2);
					}
				}
			}
			
			double[][] tissueClassScores = new double[numGroups][tissueClasses.clusters.size()];
			int[][] tissueClassCount = new int[numGroups][tissueClasses.clusters.size()];
			int[] totalTissueClassCount = new int[tissueClasses.clusters.size()];
			int totalTissues = 0;
			HashMap<String,Integer> tissueClassMap = new HashMap<String,Integer>();
			for(i=0;i<tissueClasses.clusters.size();i++) {
				cluster = tissueClasses.clusters.get(i);
				totalTissueClassCount[i] = cluster.size();
				totalTissues += totalTissueClassCount[i];
				for (j=0;j<cluster.size();j++) {
					tissueClassMap.put(cluster.get(j),i);
				}
			}
			
			for (i=0;i<tissueLoad.numRows;i++) {
				tissueClassCount[((int) tissueLoad.values[i][0])][tissueClassMap.get(tissueLoad.geneNames[i])]++;
			}
			
			int totalInGroup = 0;
			try {
				FileWriter gof = new FileWriter(groupOutFile);
				for (i=0;i<tissueClassCount.length;i++) {
					totalInGroup = 0;
					for (j=0;j<tissueClassCount[0].length;j++) {
						totalInGroup += tissueClassCount[i][j];
					}
					for (j=0;j<tissueClassCount[0].length;j++) {
						if (tissueClassCount[i][j] > 0) {
							tissueClassScores[i][j] = 1.0 - StatUtil.hyperGeometricCDF(tissueClassCount[i][j]-1,totalTissues,totalTissueClassCount[j],totalInGroup);
						} else {
							tissueClassScores[i][j] = 1.0;
						}
						gof.write((new Double(tissueClassScores[i][j])).toString());
						if (j < tissueClassCount[0].length - 1) {
							gof.write("\t");
						}
					}
					gof.write("\n");
				}
				gof.close();
			} catch(IOException ep) {
				System.out.println(ep);
			}
		}
		
		boolean[] hasMouse = new boolean[tissueClusters.clusters.size()];
		boolean[] hasHuman = new boolean[tissueClusters.clusters.size()];
		
		for (i=0;i<tissueClusters.clusters.size();i++) {
			hasMouse[i] = false;
			hasHuman[i] = false;
			cluster = tissueClusters.clusters.get(i);
			for (j=0;j<cluster.size();j++) {
				if (cluster.get(j).indexOf("m_") > -1) {
					hasMouse[i] = true;
				} else {
					hasHuman[i] = true;
				}
			}
		}
		
		int totalGO = 0;
		int totalCross = 0;
		int totalCrossGO = 0;
		
		try {
			FileWriter file = new FileWriter(outFile);
			
			file.write("topic\tnum genes\tnum tissues\tspecies\tsignif GO?");
			if (topicEntropy) {
				file.write("\tentropy");
			}
			file.write("\n");
			
			for (i=0;i<tissueClusters.clusters.size();i++) {
				file.write((new Integer(i+1)).toString());
				file.write("\t");
				file.write((new Integer(geneClusters.clusters.get(i).size())).toString());
				file.write("\t");
				file.write((new Integer(tissueClusters.clusters.get(i).size())).toString());
				file.write("\t");
				if (hasMouse[i] & hasHuman[i]) {
					file.write("MH");
					totalCross++;
					if (signifGO[i] == 1) {
						totalCrossGO++;
					}
				} else {
					if (hasMouse[i]) {
						file.write("M");
					} else {
						file.write("H");
					}
				}
				file.write("\t");
				if (signifGO[i] > 0) {
					totalGO++;
				}
				file.write((new Integer(signifGO[i])).toString());
				
				if (topicEntropy) {
					file.write("\t");
					file.write((new Double(entropies[i])).toString());
				}
				file.write("\n");
			}
			
			file.write("\n\n");
			file.write("Percent cross-species\t");
			file.write((new Double(((double) totalCross)/((double) tissueClusters.clusters.size()))).toString());
			file.write("\n");
			file.write("Percent signif. GO\t");
			file.write((new Double(((double) totalGO)/((double) tissueClusters.clusters.size()))).toString());
			file.write("\n");
			file.write("Percent signif. cross-species GO\t");
			file.write((new Double(((double) totalCrossGO)/((double) totalCross))).toString());
			file.write("\n");
			
			file.close();
		} catch(IOException ep) {
			System.out.println(ep);
		}
	}
	
	public static void scoreTissues(String inFile,String outFile) {
		String tissueClassFile = "C:\\research_data\\mouse_human\\b47mm7\\tissue_classifications.txt";
		
		ClusterReader tissueClasses = new ClusterReader();
		ClusterReader biclusters = new ClusterReader();
		
		try {
			tissueClasses.readFile(tissueClassFile);
			biclusters.readFile(inFile);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int[] classCount = new int[tissueClasses.clusters.size()];
		int[] totalClassCount = new int[tissueClasses.clusters.size()];
		double v = 0.0;
		double[] pvals = new double[tissueClasses.clusters.size()];
		int totalTissues = 0;
		double min_pval = 10e-45;
		HashMap<String,Integer> classMap = new HashMap<String,Integer>();
		ArrayList<String> cluster = null;
		HashSet<String> uniqueTissues = new HashSet<String>();
		
		int i = 0;
		int j = 0;
		int classID = 0;
		int mostAbundant = 0;
		double bestPval = 1.0;
		
		for(i=0;i<tissueClasses.clusters.size();i++) {
			cluster = tissueClasses.clusters.get(i);
			totalClassCount[i] = cluster.size();
			totalTissues += totalClassCount[i];
			for (j=0;j<cluster.size();j++) {
				cluster.set(j,removeSpaces(cluster.get(j)));
				classMap.put(cluster.get(j),i);
			}
		}
		
		for (i=0;i<biclusters.clusters.size();i++) {
			cluster = biclusters.clusters.get(i);
			for (j=0;j<cluster.size();j++) {
				uniqueTissues.add(removeSpaces(cluster.get(j)));
			}
		}
		
		Vector<Double> ap = new Vector<Double>();
		Vector<Integer> ab = new Vector<Integer>();
		
		try {
			FileWriter file = new FileWriter(outFile);
			file.write((new Integer(uniqueTissues.size())).toString());
			file.write("\n");
			for (i=0;i<biclusters.clusters.size();i++) {
				cluster = biclusters.clusters.get(i);
				for (j=0;j<classCount.length;j++) {
					pvals[j] = 0.0;
					classCount[j] = 0;
				}
				for (j=0;j<cluster.size();j++) {
					uniqueTissues.add(removeSpaces(cluster.get(j)));
					classID = classMap.get(removeSpaces(cluster.get(j)));
					classCount[classID]++;
				}
				for (j=0;j<classCount.length;j++) {
					v = 1.0 - StatUtil.hyperGeometricCDF(classCount[j]-1,totalTissues,totalClassCount[j],cluster.size());
					if (classCount[j] > 0) {
						ap.add(v);
						ab.add(i);
					}
					if (v < min_pval) {
						v = min_pval;
					}
					pvals[j] = Math.log10(v);
				}
				mostAbundant = classCount[0];
				bestPval = pvals[0];
				for (j=1;j<classCount.length;j++) {
					if (classCount[j] >= mostAbundant) {
						if (classCount[j] == mostAbundant) {
							if (pvals[j] < bestPval) {
								bestPval = pvals[j];
							}
						} else {
							bestPval = pvals[j];
						}
					}
				}
				file.write(new Double(bestPval).toString());
				file.write("\n");
			}
			
			file.close();
		} catch(IOException e) {
			System.out.println(e);
		}
		
		double[] allPvals = new double[ap.size()];
		int[] so = VectorUtil.sortOrder(ap);
		for (i=0;i<ap.size();i++) {
			allPvals[i] = ap.get(so[i]);
		}
		
		boolean[] acceptPvals = StepDownFDR.correctPvals(allPvals,0.05);
		boolean[] signifBicluster = new boolean[biclusters.clusters.size()];
		
		for (i=0;i<signifBicluster.length;i++) {
			signifBicluster[i] = false;
		}
		
		for (i=0;i<ap.size();i++) {
			if (acceptPvals[i]) {
				signifBicluster[ab.get(so[i])] = true;
			}
		}
		
		int totalSignif = 0;
		
		try {
			FileWriter tissueSigFile = new FileWriter(outFile + "_sig.txt");
			for (i=0;i<signifBicluster.length;i++) {
				tissueSigFile.write((new Integer(biclusters.clusters.get(i).size())).toString());
				tissueSigFile.write("\t");
				if (signifBicluster[i]) {
					tissueSigFile.write("1\n");
				} else {
					tissueSigFile.write("0\n");
				}
				if (signifBicluster[i])
					totalSignif++;
			}
			tissueSigFile.close();
		} catch(IOException e) {
			System.out.println(e);
		}
		
		System.out.println("Significant biclusters = " + totalSignif);
		System.out.println("Total biclusters = " + biclusters.clusters.size());
	}
	
	public static void outputGOCategoriesSAMBA() {
		String dirPath = "C:\\research_data\\expander\\";
		String biclusterFileName = dirPath + "biclusters_mouse_human_coverall.txt";
		String GOLongFileName = dirPath + "biclusters_go_long_coverall.txt";
		String GOCountFileName = dirPath + "biclusters_GO_count_coverall.txt";
		String tissuesFileName = dirPath + "biclusters_mouse_human_coverall_tissues.txt";
		String genesFileName = dirPath + "biclusters_mouse_human_coverall_genes.txt";
		
		int minGenes = 5;
		int maxGenes = 200;
		// minimum count of gene in topic
		int minOccurInTopic = 2;
		double FDR = 0.05;
		int namespace = OntTerm.BiologicalProcess;
		double minOccurPercent = 0.00;
		
		String line = "";
		String[] split = null;
		int i = 0;
		int j = 0;
		int f = 0;
		String[] names = null;
		int numBiclusters = -1;
		
		String refDirPath = "C:\\research_data\\mouse_human\\";
		String humanAssocFile = refDirPath + "Hs.genelist.go.plus";
		String mouseAssocFile = refDirPath + "Mm.genelist.go.plus";
		
		GoOBOReader goReader = new GoOBOReader();
		MouseHumanGoAssocReader humanAssocReader = new MouseHumanGoAssocReader();
		humanAssocReader.setFile(humanAssocFile);
		OntDAG humanDAG = null;
		OntGeneAssociations humanAssoc = null;
		
		MouseHumanGoAssocReader mouseAssocReader = new MouseHumanGoAssocReader();
		mouseAssocReader.setFile(mouseAssocFile);
		OntDAG mouseDAG = null;
		OntGeneAssociations mouseAssoc = null;
		
		HashSet<String> humanUseGenes = new HashSet<String>();
		HashSet<String> mouseUseGenes = new HashSet<String>();
		ArrayList<GeneList> geneLists = new ArrayList<GeneList>();
		ArrayList<String> biclusterNames = new ArrayList<String>();
		ArrayList<LinkedHashSet<String>> biclusterTissues = new ArrayList<LinkedHashSet<String>>();
		ArrayList<ArrayList<String>> biclusterGenes = new ArrayList<ArrayList<String>>();
		ArrayList<String> myGenes = new ArrayList<String>();
		GeneList myGeneList = null;
		LinkedHashSet myTissues = null;
		
		int biclusterNum = 0;
		int dataCode = 0;
		String itemName = "";
		double[] geneValues = new double[4];
		for (i=0;i<geneValues.length;i++) {
			geneValues[i] = 1.0;
		}
		
		try {
			BufferedReader readFile = new BufferedReader(new FileReader(biclusterFileName));
			line = readFile.readLine();
			while(readFile.ready()) {
				split = line.split("\t");
				biclusterNum = new Integer(split[0]).intValue();
				dataCode = new Integer(split[1]).intValue();
				itemName = split[2];
				if (biclusterNum > numBiclusters) {
					myGenes = new ArrayList<String>();
					biclusterGenes.add(myGenes);
					myGeneList = new GeneList();
					geneLists.add(myGeneList);
					myTissues = new LinkedHashSet<String>();
					biclusterTissues.add(myTissues);
					numBiclusters++;
				} else {
					myGeneList = geneLists.get(biclusterNum);
					myTissues = biclusterTissues.get(biclusterNum);
				}
				
				if (dataCode == 0) {
					myTissues.add(itemName);
				}
				
				if (dataCode == 1) {
					myGenes.add(itemName);
					names = new String[2];
					f = itemName.indexOf("|");
					names[0] = itemName.substring(0,f);
					names[1] = itemName.substring(f+1);
					humanUseGenes.add(names[0]);
					mouseUseGenes.add(names[1]);
					myGeneList.addGene(names,1.0,geneValues);
				}
				line = readFile.readLine();
			}
			readFile.close();
			ClusterReader gc = new ClusterReader(biclusterGenes);
			gc.writeFile(genesFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		Iterator<String> stringIter = null;
		for (i=0;i<biclusterTissues.size();i++) {
			myTissues = biclusterTissues.get(i);
			itemName = "";
			stringIter = myTissues.iterator();
			while(stringIter.hasNext()) {
				itemName = itemName + stringIter.next() + " | ";
			}
			biclusterNames.add(itemName);
		}
		
		try {
			humanDAG = goReader.readFile();
			System.out.println("Loaded human DAG");
			humanAssoc = humanAssocReader.readFile(humanDAG,humanUseGenes);
			System.out.println("Loaded human associations");
			
			mouseDAG = goReader.readFile();
			System.out.println("Loaded mouse DAG");
			mouseAssoc = mouseAssocReader.readFile(mouseDAG,mouseUseGenes);
			System.out.println("Loaded mouse associations");
		} catch(IOException e) {
			System.out.println(e);
		} catch(DAGException e) {
			System.out.println(e);
		}
		
		OntClusters humanClusters = new OntClusters(humanAssoc,minGenes,maxGenes,FDR,namespace);
		OntClusters mouseClusters = new OntClusters(mouseAssoc,minGenes,maxGenes,FDR,namespace);
		ArrayList<String> humanGeneSetNames = null;
		ArrayList<String> mouseGeneSetNames = null;
		
		for (i=0;i<geneLists.size();i++) {
			humanGeneSetNames = new ArrayList<String>();
			mouseGeneSetNames = new ArrayList<String>();
			myGeneList = geneLists.get(i);
			myGeneList.sort();
			for (j=0;j<myGeneList.size();j++) {
				humanGeneSetNames.add(myGeneList.genes.get(j).names[0]);
				mouseGeneSetNames.add(myGeneList.genes.get(j).names[1]);
			}
			humanClusters.addCluster(humanGeneSetNames,myGeneList.getRanks(),"none",0.0);
			mouseClusters.addCluster(mouseGeneSetNames,myGeneList.getRanks(),"none",0.0);
		}
		
		humanClusters.clusterSignif();
		mouseClusters.clusterSignif();
		
		int[] humanSignif = humanClusters.numSignificantTerms();
		int[] mouseSignif = mouseClusters.numSignificantTerms();
		try {
			FileWriter gf = new FileWriter(GOCountFileName);
			for (i=0;i<humanSignif.length;i++) {
				gf.write((new Integer(humanSignif[i])).toString());
				gf.write("\t");
				gf.write((new Integer(mouseSignif[i])).toString());
				if (humanSignif[i] > 0 | mouseSignif[i] > 0) {
					gf.write("\t1\n");
				} else {
					gf.write("\t0\n");
				}
			}
			gf.close();
		} catch(IOException e) {
			System.out.println(e);
		}
		
		try {
			FileWriter outFile = new FileWriter(GOLongFileName);
			for (i=0;i<geneLists.size();i++) {
				outFile.write(biclusterNames.get(i));
				outFile.write("\n");
				humanClusters.getClusters().get(i).outputScoresLong(outFile,false);
				mouseClusters.getClusters().get(i).outputScoresLong(outFile,false);
				outFile.write("\n");
			}
		} catch(IOException e) {
			System.out.println(e);
		}
		
		try {
			FileWriter tissueWriter = new FileWriter(tissuesFileName);
			for (i=0;i<biclusterTissues.size();i++) {
				myTissues = biclusterTissues.get(i);
				stringIter = myTissues.iterator();
				while(stringIter.hasNext()) {
					tissueWriter.write(stringIter.next());
					if (stringIter.hasNext()) {
						tissueWriter.write("\t");
					}
				}
				tissueWriter.write("\n");
			}
			tissueWriter.close();
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void convertBicatFile() {
		String dirPath = "C:\\research_data\\BicAT_v2.22\\";
		String biclusterFileName = dirPath + "ISA_results_filtered";
		
		String line = null;
		String[] split = null;
		int i = 0;
		
		try {
			BufferedReader readFile = new BufferedReader(new FileReader(biclusterFileName+".txt"));
			FileWriter geneOutFile = new FileWriter(biclusterFileName+"_genes.txt");
			FileWriter tissuesOutFile = new FileWriter(biclusterFileName+"_tissues.txt");
			while(readFile.ready()) {
				line = readFile.readLine();
				line = readFile.readLine();
				split = line.split(" ");
				for (i=0;i<split.length;i++) {
					geneOutFile.write(split[i]);
					if (i < split.length - 1) {
						geneOutFile.write("\t");
					}
				}
				geneOutFile.write("\n");
				
				line = readFile.readLine();
				split = line.split(" ");
				for (i=0;i<split.length;i++) {
					tissuesOutFile.write(split[i]);
					if (i < split.length - 1) {
						tissuesOutFile.write("\t");
					}
				}
				tissuesOutFile.write("\n");
			}
			geneOutFile.close();
			tissuesOutFile.close();
			readFile.close();
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static String removeSpaces(String s) {
		if (s.indexOf(" ") == -1) {
			return s;
		}
		
		String s2 = "";
		String[] split = s.split(" ");
		int i = 0;
		for (i=0;i<split.length;i++) {
			s2 = s2 + split[i];
		}
		
		return s2;
	}
	
	public static void scoreGenes(String genesFileName) {
		
		int minGenes = 5;
		int maxGenes = 200;
		// minimum count of gene in topic
		int minOccurInTopic = 2;
		double FDR = 0.05;
		int namespace = OntTerm.BiologicalProcess;
		double minOccurPercent = 0.00;
		
		int i = 0;
		int j = 0;
		int f = 0;
		
		ClusterReader biclusters = new ClusterReader();
		try {
			biclusters.readFile(genesFileName + ".txt");
		} catch(IOException e) {
			System.out.println(e);
		}
		
		String[] names = null;
		
		String refDirPath = "C:\\research_data\\mouse_human\\";
		String humanAssocFile = refDirPath + "Hs.genelist.go.plus";
		String mouseAssocFile = refDirPath + "Mm.genelist.go.plus";
		
		GoOBOReader goReader = new GoOBOReader();
		MouseHumanGoAssocReader humanAssocReader = new MouseHumanGoAssocReader();
		humanAssocReader.setFile(humanAssocFile);
		OntDAG humanDAG = null;
		OntGeneAssociations humanAssoc = null;
		
		MouseHumanGoAssocReader mouseAssocReader = new MouseHumanGoAssocReader();
		mouseAssocReader.setFile(mouseAssocFile);
		OntDAG mouseDAG = null;
		OntGeneAssociations mouseAssoc = null;
		
		HashSet<String> humanUseGenes = new HashSet<String>();
		HashSet<String> mouseUseGenes = new HashSet<String>();
		ArrayList<GeneList> geneLists = new ArrayList<GeneList>();
		GeneList myGeneList = null;
		
		HashSet<String> uniqueGenes = new HashSet<String>();
		
		double[] geneValues = new double[4];
		for (i=0;i<geneValues.length;i++) {
			geneValues[i] = 1.0;
		}
		
		ArrayList<String> genes = null;
		String geneNames = null;
		
		for (i=0;i<biclusters.clusters.size();i++) {
			genes = biclusters.clusters.get(i);
			myGeneList = new GeneList();
			geneLists.add(myGeneList);
			for (j=0;j<genes.size();j++) {
				names = new String[2];
				geneNames = genes.get(j);
				uniqueGenes.add(geneNames);
				f = geneNames.indexOf("|");
				names[0] = geneNames.substring(0,f);
				names[1] = geneNames.substring(f+1);
				humanUseGenes.add(names[0]);
				mouseUseGenes.add(names[1]);
				myGeneList.addGene(names,1.0,geneValues);
			}
		}
		
		try {
			humanDAG = goReader.readFile();
			System.out.println("Loaded human DAG");
			humanAssoc = humanAssocReader.readFile(humanDAG,humanUseGenes);
			System.out.println("Loaded human associations");
			
			mouseDAG = goReader.readFile();
			System.out.println("Loaded mouse DAG");
			mouseAssoc = mouseAssocReader.readFile(mouseDAG,mouseUseGenes);
			System.out.println("Loaded mouse associations");
		} catch(IOException e) {
			System.out.println(e);
		} catch(DAGException e) {
			System.out.println(e);
		}
		
		OntClusters humanClusters = new OntClusters(humanAssoc,minGenes,maxGenes,FDR,namespace);
		OntClusters mouseClusters = new OntClusters(mouseAssoc,minGenes,maxGenes,FDR,namespace);
		ArrayList<String> humanGeneSetNames = null;
		ArrayList<String> mouseGeneSetNames = null;
		
		for (i=0;i<geneLists.size();i++) {
			humanGeneSetNames = new ArrayList<String>();
			mouseGeneSetNames = new ArrayList<String>();
			myGeneList = geneLists.get(i);
			myGeneList.sort();
			for (j=0;j<myGeneList.size();j++) {
				humanGeneSetNames.add(myGeneList.genes.get(j).names[0]);
				mouseGeneSetNames.add(myGeneList.genes.get(j).names[1]);
			}
			humanClusters.addCluster(humanGeneSetNames,myGeneList.getRanks(),"none",0.0);
			mouseClusters.addCluster(mouseGeneSetNames,myGeneList.getRanks(),"none",0.0);
		}
		
		humanClusters.clusterSignif();
		mouseClusters.clusterSignif();
		
		int[] humanSignif = humanClusters.numSignificantTerms();
		int[] mouseSignif = mouseClusters.numSignificantTerms();
		try {
			FileWriter gf = new FileWriter(genesFileName + "_score.txt");
			for (i=0;i<humanSignif.length;i++) {
				gf.write((new Integer(humanSignif[i])).toString());
				gf.write("\t");
				gf.write((new Integer(mouseSignif[i])).toString());
				if (humanSignif[i] > 0 | mouseSignif[i] > 0) {
					gf.write("\t1\n");
				} else {
					gf.write("\t0\n");
				}
			}
			gf.close();
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int mostAbundant = 0;
		double bestPVal = 1.0;
		
		OntTermScores otsH = null;
		OntTermScores otsM = null;
		Iterator<OntScore> scoreIter = null;
		Iterator<OntTerm> termIter = null;
		OntTerm term = null;
		OntScore score = null;
		double[] bestScores = new double[humanClusters.clusters.size()];
		double minPVal = 10e-45;
		
		try {
			FileWriter gfb = new FileWriter(genesFileName + "_bestcats.txt");
			for (i=0;i<humanClusters.clusters.size();i++) {
				otsH = humanClusters.clusters.get(i);
				otsM = humanClusters.clusters.get(i);
				mostAbundant = 0;
				bestPVal = 1.0;
				scoreIter = otsH.termMap.values().iterator();
				termIter = otsH.termMap.keySet().iterator();
				String bestGOTerm = "";
				String bestGODesc = "";
				boolean isSignif = false;
				
				while(scoreIter.hasNext()) {
					score = scoreIter.next();
					term = termIter.next();
					if (score.genes.size() >= mostAbundant & score.acceptPVal) {
						if (score.genes.size() == mostAbundant) {
							if (score.pval < bestPVal) {
								mostAbundant = score.genes.size();
								bestPVal = score.pval;
								bestGOTerm = term.ID;
								bestGODesc = term.name;
								isSignif = score.acceptPVal;
							}
						} else {
							mostAbundant = score.genes.size();
							bestPVal = score.pval;
							bestGOTerm = term.ID;
							bestGODesc = term.name;
							isSignif = score.acceptPVal;
						}
						
					}
				}
				
				scoreIter = otsM.termMap.values().iterator();
				termIter = otsM.termMap.keySet().iterator();
				while(scoreIter.hasNext()) {
					score = scoreIter.next();
					term = termIter.next();
					if (score.genes.size() >= mostAbundant & score.acceptPVal) {
						if (score.genes.size() == mostAbundant) {
							if (score.pval < bestPVal) {
								mostAbundant = score.genes.size();
								bestPVal = score.pval;
								bestGOTerm = term.ID;
								bestGODesc = term.name;
								isSignif = score.acceptPVal;
							}
						} else {
							mostAbundant = score.genes.size();
							bestPVal = score.pval;
							bestGOTerm = term.ID;
							bestGODesc = term.name;
							isSignif = score.acceptPVal;
						}
						
					}
				}
				
				if (bestPVal < minPVal) {
					bestPVal = minPVal;
				}
				
				if (isSignif) {
					gfb.write(bestGOTerm + "\t" + bestGODesc + "\t" + ((new Double(bestPVal)).toString()) + "\n");
				} else {
					gfb.write("none" + "\t" + "none" + "\t" + "N/A" + "\n");
				}
 				
				
				bestScores[i] = Math.log10(bestPVal);
			}
			gfb.close();
		} catch(IOException e) {
			System.out.println(e);
		}
		
		try {
			FileWriter gf = new FileWriter(genesFileName + "_cumscore.txt");
			gf.write((new Integer(uniqueGenes.size())).toString());
			gf.write("\n");
			for (i=0;i<bestScores.length;i++) {
				gf.write((new Double(bestScores[i])).toString());
				gf.write("\n");
			}
			gf.close();
		} catch(IOException e) {
			System.out.println(e);
		}
	}
}
