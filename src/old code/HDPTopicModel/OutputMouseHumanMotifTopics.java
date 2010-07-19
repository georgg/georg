package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.Vector;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.StatUtil.StatUtil;
import edu.mit.csail.psrg.georg.StatUtil.VectorUtil;

public class OutputMouseHumanMotifTopics {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	//	outputFiles();
		outputStrongModules();
	}
	
	// construct strong modules (transposed) for existing topics in HDP
	public static ArrayList<StrongModuleTranspose> constructStrongModules(HierDirichletProcess myHDP) {
		ArrayList<StrongModuleTranspose> modules = new ArrayList<StrongModuleTranspose>();
		
		ArrayList<DirichletProcess> clusterDPs = new ArrayList<DirichletProcess>();
		
		int i = 0;
		int j = 0;
		for (i=0;i<myHDP.DP.length;i++) {
			if (myHDP.DP[i].label.indexOf("cluster") == 0) {
				clusterDPs.add(myHDP.DP[i]);
			}
		}
		
		int docNum = 0;
		DirichletProcess dp = null;
		DirichletProcess clusterDP = null;
		for (i=0;i<clusterDPs.size();i++) {
			clusterDP = clusterDPs.get(i);
			for (j=0;j<clusterDP.children.size();j++) {
				dp = clusterDP.children.get(j);
				dp.documents[0].code = docNum;
				docNum++;
			}
		}
		
		StrongModuleTransposeCollector collector = new StrongModuleTransposeCollector(myHDP,clusterDPs);
		myHDP.addStatCollector(collector);
		collector.updateStats();
		
		modules = collector.collapseModules();
		
	/*	Iterator<StrongModuleTranspose> iter = collector.strongModuleMap.values().iterator();
		while (iter.hasNext()) {
			modules.add(iter.next());
		} */
		
		return modules;
	}
	
	public static void outputGeneMotifMatrix(ArrayList<StrongModuleTranspose> modules,StrongModuleTransposeCollector myCollector,String[] clusterNames,HierDirichletProcess myHDP,String outFileName) {
		ArrayList<LinkedHashMap<Integer,ArrayList<Double>>> clm = new ArrayList<LinkedHashMap<Integer,ArrayList<Double>>>();
		LinkedHashMap<Integer,ArrayList<Double>> clmMap = null;
		ArrayList<Double> geneCounts = null;
		
		int i = 0;
		int j = 0;
		int mm = 0;
		int clusterID = 0;
		int geneID = 0;
		int count = 0;
		Iterator<Integer> clusterIter = null;
		Iterator<Integer> geneIter = null;
		StrongModuleTranspose module = null;
		HashMap<Integer,Integer> geneMap = null;
		
		for (i=0;i<myCollector.clusterDPList.size();i++) {
			clm.add(new LinkedHashMap<Integer,ArrayList<Double>>());
		}
		
		for (mm=0;mm<modules.size();mm++) {
			module = modules.get(mm);
			clusterIter = module.clusterTotal.keySet().iterator();
			while (clusterIter.hasNext()) {
				clusterID = clusterIter.next();
				clmMap = clm.get(clusterID);
				geneMap = module.clusterTotal.get(clusterID);
				geneIter = geneMap.keySet().iterator();
				while(geneIter.hasNext()) {
					geneID = geneIter.next();
					if (!clmMap.containsKey(geneID)) {
						geneCounts = new ArrayList<Double>();
						for (i=0;i<modules.size();i++) {
							geneCounts.add(0.0);
						}
						clmMap.put(geneID,geneCounts);
					} else {
						geneCounts = clmMap.get(geneID);
					}
					count = geneMap.get(geneID);
					geneCounts.set(mm,count+geneCounts.get(mm));
				}
			}
		}
		
		Vector<Double> v = new Vector<Double>();
		double norm = 0.0;
		int[] order = null;
		String geneName = null;
		HashMap<Integer,String> geneNames = new HashMap<Integer,String>();
		for (i=0;i<myHDP.DP.length;i++) {
			if (myHDP.DP[i].documents != null) {
				geneNames.put(myHDP.DP[i].documents[0].code,myHDP.DP[i].documents[0].annotation);
			}
		}
		
		
		try {
			FileWriter file = new FileWriter(outFileName);
			for (clusterID=0;clusterID<clm.size();clusterID++) {
				file.write(clusterNames[clusterID]);
				file.write("\n");
				clmMap = clm.get(clusterID);
				geneIter = clmMap.keySet().iterator();
				while(geneIter.hasNext()) {
					geneID = geneIter.next();
					geneCounts = clmMap.get(geneID);
					v.clear();
					norm = 0.0;
					for (i=0;i<geneCounts.size();i++) {
						v.add(geneCounts.get(i));
						norm += geneCounts.get(i);
					}
					if (norm > 0.0) {
						for (i=0;i<v.size();i++) {
							v.set(i,-v.get(i)/norm);
						}
					}
					order = VectorUtil.sortOrder(v);
					file.write(geneNames.get(geneID));
					for (i=0;i<order.length;i++) {
						file.write("\t");
						file.write((new Integer(order[i]+1)).toString());
					}
					file.write("\n");
					file.write("--");
					for (i=0;i<order.length;i++) {
						file.write("\t");
						file.write((new Double(-v.get(order[i]))).toString());
					}
					file.write("\n");
				}
				file.write("\n");
			}
			file.close();
		} catch(IOException e) {
			
		}
	}
	
	public static void outputStrongModules() {
		boolean constructStrongModules = false;
	//	String dirPath = "C:\\research_data\\mouse_human\\b47mm7\\geo\\motifs_cons_expression\\6\\";
		String dirPath = "C:\\research_data\\mouse_human\\b47mm7\\geo\\motifs\\50_6\\";
		String snapFileName = dirPath + "genes_as_docs_strong_snap_49999.persist";
		String motifTopicsFileName = dirPath + "motif_topics.txt";
		String motifLoadMatrixFileName = dirPath + "motif_load.txt";
		String clusterLoadMatrixFileName = dirPath + "cluster_load.txt";
		String clusterRankMatrixFileName = dirPath + "cluster_rank.txt";
		String geneMotifMatrixFileName = dirPath + "gene_load.txt";
		
		HierDirichletProcess myHDP = HierDirichletProcess.restoreFromFile(snapFileName);
		int i = 0;
		int j = 0;
		Vector<Double> v = new Vector<Double>();
		
		MicroArrayData dataFile = new MicroArrayData();
		
		StrongModuleTransposeCollector myCollector = null;
		ArrayList<StrongModuleTranspose> modules = null;
		
		if (constructStrongModules) {
			modules = constructStrongModules(myHDP);
		}
		
		myCollector = (StrongModuleTransposeCollector) myHDP.getStatCollector(StrongModuleTransposeCollector.class);
		myCollector.finalFilterModuleOccurPercent = 0.05;
		myCollector.minClusterOverlap = 0.50;
	/*	myCollector.minModuleOccurPercent = 0.0005;
		myCollector.finalFilterModuleOccurPercent = 0.005;
		myCollector.minClusterOverlap = 0.80; */
		
		if (!constructStrongModules) {
			modules = myCollector.collapseModules();
		}
		
		String[] motifNames = new String[myHDP.totalGenes];
		for (i=0;i<myHDP.totalGenes;i++) {
			motifNames[i] = myHDP.geneNames[myHDP.reverseGeneMap[i]];
		}
		
		dataFile.values = myCollector.generateMotifLoadMatrix(modules);
		dataFile.numCols = dataFile.values[0].length;
		dataFile.numRows = dataFile.values.length;
		dataFile.geneNames = motifNames;
		
		try {
			dataFile.writeFile(motifLoadMatrixFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		double[] modulePercents = new double[modules.size()];
		for (i=0;i<modules.size();i++) {
			modulePercents[i] = ((double) modules.get(i).numOccur)/((double) myCollector.numSamples);
		}
		
		String[] clusterNames = new String[myCollector.clusterDPList.size()];
		for (i=0;i<clusterNames.length;i++) {
			clusterNames[i] = myCollector.clusterDPList.get(i).label;
		}
		
		outputGeneMotifMatrix(modules,myCollector,clusterNames,myHDP,geneMotifMatrixFileName);
		
		double[][] clusterLoadValues = myCollector.generateClusterLoadMatrix(modules);
		double norm = 0.0;
		
		int[] sortOrder = null;
		try {
			FileWriter writer = new FileWriter(motifTopicsFileName);
			for (j=0;j<dataFile.values[0].length;j++) {
				writer.write("Topic#" + (new Integer(j+1)).toString());
				writer.write("\t");
				writer.write((new Double(modulePercents[j])).toString());
				writer.write("\n");
				v.clear();
				for (i=0;i<dataFile.values.length;i++) {
					v.add(-dataFile.values[i][j]);
				}
				sortOrder = VectorUtil.sortOrder(v);
				for (i=0;i<dataFile.values.length;i++) {
					if (dataFile.values[sortOrder[i]][j] > 0) {
						writer.write(dataFile.geneNames[sortOrder[i]]);
						writer.write("\t");
						writer.write((new Double(dataFile.values[sortOrder[i]][j])).toString());
						writer.write("\n");
					}
				}
			/*	v.clear();
				norm = 0.0;
				for (i=0;i<clusterLoadValues.length;i++) {
					v.add(-clusterLoadValues[i][j]);
					norm += clusterLoadValues[i][j];
				}
				writer.write("------------\n");
				sortOrder = VectorUtil.sortOrder(v);
				for (i=0;i<sortOrder.length;i++) {
					if (v.get(sortOrder[i]) < 0.0) {
						writer.write(clusterNames[sortOrder[i]]);
						writer.write("\t");
						writer.write((new Double(-v.get(sortOrder[i]))).toString());
						writer.write("\n");
					}
				} */
				
				writer.write("\n\n");
			}
			writer.close();
		} catch (IOException e) {
			System.out.println(e);
		}
		
	/*	for (i=0;i<clusterLoadValues.length;i++) {
			norm = 0;
			for (j=0;j<clusterLoadValues[0].length;j++) {
				if (!Double.isNaN(clusterLoadValues[i][j])) {
					norm += clusterLoadValues[i][j];
				}
			}
			if (norm > 0.0) {
				for (j=0;j<clusterLoadValues[0].length;j++) {
					clusterLoadValues[i][j] = clusterLoadValues[i][j]/norm;
				}
			}
		} */
		
		dataFile.values = clusterLoadValues;
		dataFile.numCols = dataFile.values[0].length;
		dataFile.numRows = dataFile.values.length;
		dataFile.geneNames = clusterNames;
		
		try {
			dataFile.writeFile(clusterLoadMatrixFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int[][] clusterRank = new int[dataFile.values.length][dataFile.values[0].length];
		for (i=0;i<clusterRank.length;i++) {
			v.clear();
			for (j=0;j<clusterRank[0].length;j++) {
				v.add(-dataFile.values[i][j]);
			}
			clusterRank[i] = VectorUtil.sortOrder(v);
			// correct for the fact that it's zero indexed
			for (j=0;j<clusterRank[0].length;j++) {
				clusterRank[i][j]++;
			}
		}
		dataFile.dvalues = clusterRank;
		dataFile.discreteData = true;
		dataFile.values = null;
		
		try {
			dataFile.writeFile(clusterRankMatrixFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	
	public static void outputFiles() {
		String dirPath = "C:\\research_data\\mouse_human\\motifs\\";
		String snapFileName = dirPath + "genes_as_docs_snap_98999.persist";
		String dpFileName = dirPath + "genes_as_docs_dps.txt";
		String topicFileName = dirPath + "genes_as_docs_topics.txt";
		String dpOutFileName = dirPath + "genes_as_docs_dps_rank.txt";
		String topicOutFileName = dirPath + "genes_as_docs_topics_rank.txt";
		HierDirichletProcess myDP = HierDirichletProcess.restoreFromFile(snapFileName);
		
		MicroArrayData dpData = new MicroArrayData();
		MicroArrayData topicData = new MicroArrayData();
		
		try {
			myDP.outputClustersToFile(dirPath + "genes_as_docs");
			dpData.readFile(dpFileName);
			topicData.readFile(topicFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		ArrayList<Integer> clusters = new ArrayList<Integer>();
		int k = 0;
		int i = 0;
		int j = 0;
		
		for (i=0;i<dpData.numRows;i++) {
			if (dpData.geneNames[i].indexOf("cluster") > -1) {
		//	if (dpData.geneNames[i].indexOf("NM_") > -1) {
				clusters.add(i);
			}
		}
		
		double norm = 0.0;
		Vector<Double> v = new Vector<Double>();
		int[] order = null;
		FileWriter outFile = null;
		double[][] normDP = new double[clusters.size()][dpData.numCols];
		double[] clusterTotal = new double[clusters.size()];
		
		Vector<Integer> v2 = new Vector<Integer>();
		try {
			outFile = new FileWriter(dpOutFileName);
			for (k=0;k<clusters.size();k++) {
				i = clusters.get(k);
				norm = 0.0;
				v.clear();
				v2.clear();
				for (j=0;j<dpData.numCols;j++) {
					norm += dpData.values[i][j];
				}
				clusterTotal[k] = norm;
				for (j=0;j<dpData.numCols;j++) {
					normDP[k][j] = dpData.values[i][j]/norm;
					if (dpData.values[i][j] >= 5.0) {
						v.add(dpData.values[i][j]/norm);
						v2.add(j);
					}
				}
				order = VectorUtil.sortOrder(v);
				outFile.write(dpData.geneNames[i]);
				for (j=order.length-1;j>=0;j--) {
					outFile.write("\t");
					outFile.write((new Integer(v2.get(order[j])+1)).toString());
				}
				outFile.write("\n");
				outFile.write("---");
				for (j=order.length-1;j>=0;j--) {
					outFile.write("\t");
					outFile.write((new Double(v.get(order[j])).toString()));
				}
				outFile.write("\n");
			}
			outFile.close();
		} catch(IOException e) {
			System.out.println(e);
		}
		
		double[] motifNorm = new double[topicData.numRows];
		double totalCount = 0.0;
		double pval;
		for(i=0;i<topicData.numRows;i++) {
			for (j=0;j<topicData.numCols;j++) {
				motifNorm[i] += topicData.values[i][j];
				totalCount += topicData.values[i][j];
			}
		}
		
		try {
			outFile = new FileWriter(topicOutFileName);
			for (j=0;j<topicData.numCols;j++) {
				norm = 0.0;
				for (i=0;i<topicData.numRows;i++) {
					norm += topicData.values[i][j];
				}
				v.clear();
				v2.clear();
				for (i=0;i<topicData.numRows;i++) {
					pval = 1.0-StatUtil.hyperGeometricCDF((int)topicData.values[i][j],(int)totalCount,(int)motifNorm[i],(int)norm);
				//	if (topicData.values[i][j] >= 5.0) {
					if (pval <= 0.001) {
						v2.add(i);
						v.add(topicData.values[i][j]/norm);
					//	v.add(pval);
					}
				}
				order = VectorUtil.sortOrder(v);
				outFile.write("topic#" + (j+1));
				for (i=order.length-1;i>=0;i--) {
			//	for (i=0;i<order.length;i++) {
					outFile.write("\t");
					outFile.write(topicData.geneNames[v2.get(order[i])]);
				}
				outFile.write("\n");
				outFile.write("---");
				for (i=order.length-1;i>=0;i--) {
			//	for (i=0;i<order.length;i++) {
					outFile.write("\t");
					outFile.write((new Double(v.get(order[i])).toString()));
				}
				outFile.write("\n");
				
				norm = 0.0;
				v.clear();
				v2.clear();
				
				for (k=0;k<clusters.size();k++) {
					i = clusters.get(k);
					norm += dpData.values[i][j];
				}
				for (k=0;k<clusters.size();k++) {
					i = clusters.get(k);
					pval = 1.0-StatUtil.hyperGeometricCDF((int)dpData.values[i][j],(int)totalCount,(int)clusterTotal[k],(int)norm);
				//	if (dpData.values[i][j] >= 5.0) {
					if (pval <= 0.001) {
						v2.add(i);
					//	v.add(normDP[k][j]*dpData.values[i][j]/norm);
						v.add(normDP[k][j]);
					//	v.add(dpData.values[i][j]/norm);
					//	v.add(pval);
					}
				}
				order = VectorUtil.sortOrder(v);
				for (k=order.length-1;k>=0;k--) {
			//	for (k=0;k<order.length;k++) {
					outFile.write((new Double(v.get(order[k]))).toString());
					outFile.write("\t");
					outFile.write(dpData.geneNames[v2.get(order[k])]);
					outFile.write("\n");
				}
				outFile.write("\n");
			}
			outFile.close();
		} catch(IOException e) {
			System.out.println(e);
		}
	}

}
