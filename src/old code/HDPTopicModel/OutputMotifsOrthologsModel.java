package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Vector;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.GO.DAGException;
import edu.mit.csail.psrg.georg.GO.GoOBOReader;
import edu.mit.csail.psrg.georg.GO.MouseHumanGoAssocReader;
import edu.mit.csail.psrg.georg.GO.OntClusters;
import edu.mit.csail.psrg.georg.GO.OntDAG;
import edu.mit.csail.psrg.georg.GO.OntGene;
import edu.mit.csail.psrg.georg.GO.OntGeneAssociations;
import edu.mit.csail.psrg.georg.GO.OntTerm;
import edu.mit.csail.psrg.georg.StatUtil.StatUtil;
import edu.mit.csail.psrg.georg.StatUtil.VectorUtil;

public class OutputMotifsOrthologsModel {
	
	public static void main(String[] args) {
		outputRankedFiles();
	}
	
	public static void outputRankedFiles() {
		String dirPath = "C:\\research_data\\mouse_human\\motifs_orthologs\\";
		String snapFileName = dirPath + "genes_as_docs_snap_49999.persist";
		String dpFileName = dirPath + "genes_as_docs_dps.txt";
		String topicFileName = dirPath + "genes_as_docs_topics.txt";
		String dpOutFileName = dirPath + "genes_as_docs_dps_rank.txt";
		String topicOutFileName = dirPath + "genes_as_docs_topics_rank.txt";
		String humanGOFileName = dirPath + "dps_rank_human_GO.txt";
		String mouseGOFileName = dirPath + "dps_rank_mouse_GO.txt";
		
		String refDirPath = "C:\\research_data\\mouse_human\\";
		String humanAssocFile = refDirPath + "Hs.genelist.go.plus";
		String mouseAssocFile = refDirPath + "Mm.genelist.go.plus";
		
		OrthoHierDirichletProcess myDP = (OrthoHierDirichletProcess) HierDirichletProcess.restoreFromFile(snapFileName);
		
		MicroArrayData dpData = new MicroArrayData();
		MicroArrayData topicData = new MicroArrayData();
		
		try {
			myDP.outputClustersToFile(dirPath + "genes_as_docs");
			dpData.readFile(dpFileName);
			topicData.readFile(topicFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		GoOBOReader goReader = new GoOBOReader();
		MouseHumanGoAssocReader humanAssocReader = new MouseHumanGoAssocReader();
		humanAssocReader.setFile(humanAssocFile);
		OntDAG humanDAG = null;
		OntGeneAssociations humanAssoc = null;
		
		MouseHumanGoAssocReader mouseAssocReader = new MouseHumanGoAssocReader();
		mouseAssocReader.setFile(mouseAssocFile);
		OntDAG mouseDAG = null;
		OntGeneAssociations mouseAssoc = null;
		
		ArrayList<Integer> clusters = new ArrayList<Integer>();
		int k = 0;
		int i = 0;
		int j = 0;
		
		HashSet<String> humanGenes = new HashSet<String>();
		HashSet<String> mouseGenes = new HashSet<String>();
		
		for (i=0;i<myDP.myOrthologStats.humanDPs.length;i++) {
			humanGenes.add(myDP.myOrthologStats.humanDPs[i].label);
			mouseGenes.add(myDP.myOrthologStats.mouseDPs[i].label);
		}
		
		for (i=0;i<dpData.numRows;i++) {
			if (dpData.geneNames[i].indexOf("NM_") > -1) {
				clusters.add(i);
			}
		}
		
		try {
			humanDAG = goReader.readFile();
			System.out.println("Loaded human DAG");
			humanAssoc = humanAssocReader.readFile(humanDAG,humanGenes);
			System.out.println("Loaded human associations");
			
			mouseDAG = goReader.readFile();
			System.out.println("Loaded mouse DAG");
			mouseAssoc = mouseAssocReader.readFile(mouseDAG,mouseGenes);
			System.out.println("Loaded mouse associations");
		} catch(IOException e) {
			System.out.println(e);
		} catch(DAGException e) {
			System.out.println(e);
		}
		
		int minGenes = 5;
		int maxGenes = 200;
		double FDR = 0.05;
		int namespace = OntTerm.BiologicalProcess;
		OntClusters humanClusters = new OntClusters(humanAssoc,minGenes,maxGenes,FDR,namespace);
		OntClusters mouseClusters = new OntClusters(mouseAssoc,minGenes,maxGenes,FDR,namespace);
		
		double norm = 0.0;
		Vector<Double> v = new Vector<Double>();
		int[] order = null;
		FileWriter outFile = null;
		double[][] normDP = new double[clusters.size()][dpData.numCols];
		double[] clusterTotal = new double[clusters.size()];
		String geneName;
		String commonGeneName;
		String geneDesc;
		OntGene oGene = null;
		double pValThresh = 10e-6;
		
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
		
		ArrayList<String> humanList = new ArrayList<String>();
		ArrayList<Double> humanRank = new ArrayList<Double>();
		ArrayList<String> mouseList = new ArrayList<String>();
		ArrayList<Double> mouseRank = new ArrayList<Double>();
		String clusterDescription;
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
					if (pval <= pValThresh & norm > 0.0) {
						v2.add(i);
						v.add(topicData.values[i][j]/norm);
					//	v.add(pval);
					}
				}
				order = VectorUtil.sortOrder(v);
				outFile.write("topic#" + (j+1));
				clusterDescription = "";
				for (i=order.length-1;i>=0;i--) {
			//	for (i=0;i<order.length;i++) {
					outFile.write("\t");
					outFile.write(topicData.geneNames[v2.get(order[i])]);
					clusterDescription = clusterDescription + topicData.geneNames[v2.get(order[i])] + "(" + (new Double(v.get(order[i])).toString()) + ")\t";
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
					if (pval <= pValThresh & normDP[k][j] > 0.0) {
						v2.add(i);
					//	v.add(normDP[k][j]*dpData.values[i][j]/norm);
						v.add(normDP[k][j]);
					//	v.add(dpData.values[i][j]/norm);
					//	v.add(pval);
					}
				}
				order = VectorUtil.sortOrder(v);
				humanList.clear();
				mouseList.clear();
				humanRank.clear();
				mouseRank.clear();
				for (k=order.length-1;k>=0;k--) {
			//	for (k=0;k<order.length;k++) {
					outFile.write((new Double(v.get(order[k]))).toString());
					outFile.write("\t");
					geneName = dpData.geneNames[v2.get(order[k])];
					oGene = null;
					if (humanGenes.contains(geneName)) {
						humanList.add(geneName);
						humanRank.add(v.get(order[k]));
						oGene = humanAssoc.getGene(geneName);
						geneName = "h_" + geneName;
					}
					if (mouseGenes.contains(geneName)) {
						mouseList.add(geneName);
						mouseRank.add(v.get(order[k]));
						oGene = mouseAssoc.getGene(geneName);
						geneName = "m_" + geneName;
					}
					outFile.write(geneName);
					outFile.write("\t");
					outputGeneDescription(outFile,oGene);
					outFile.write("\n");
				}
				if (humanList.size() > 0) {
					humanClusters.addCluster(humanList,humanRank,clusterDescription,0.0);
				}
				if (mouseList.size() > 0)  {
					mouseClusters.addCluster(mouseList,mouseRank,clusterDescription,0.0);
				}
				outFile.write("\n");
			}
			outFile.close();
			
			humanClusters.clusterSignif();
			mouseClusters.clusterSignif();
			humanClusters.outputClustersLong(humanGOFileName);
			mouseClusters.outputClustersLong(mouseGOFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void outputGeneDescription(FileWriter outFile,OntGene gene) throws IOException {
		if (gene == null) {
			return;
		}
		if (gene.Name != null) {
			outFile.write(gene.Name);
		} else {
			outFile.write("none");
		}
		outFile.write("\t");
		if (gene.Description != null) {
			outFile.write(gene.Description);
		} else {
			outFile.write("none");
		}
	}
}
