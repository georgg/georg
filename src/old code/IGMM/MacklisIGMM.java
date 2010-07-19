package edu.mit.csail.psrg.georg.IGMM;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.GO.AffyGoAssocReader;
import edu.mit.csail.psrg.georg.GO.DAGException;
import edu.mit.csail.psrg.georg.GO.GoOBOReader;
import edu.mit.csail.psrg.georg.GO.OntDAG;
import edu.mit.csail.psrg.georg.GO.OntGene;
import edu.mit.csail.psrg.georg.GO.OntGeneAssociations;

public class MacklisIGMM {
	public static void main(String[] args) {
	//	cluster();
	//	outputClusters();
		displayClusters();
	}
	
	public static void cluster() {
		MicroArrayData myData = processData();
		System.out.println("Processed data");
		
		GaussianDP DP = new GaussianDP(myData,1);
		
		int iter = 20000;
		int burnin = 5000;
		
		String mainOutName = "C:\\research_data\\macklis\\macklis_IGMM";
		String snapShotName = mainOutName + "_snap";
		DP.enableSnapShots(500,snapShotName);
		
		DP.iterate(iter,burnin,15);
	}
	
	public static void outputClusters() {
		String mainOutName = "C:\\research_data\\macklis\\macklis_IGMM";
		String snapShotName = mainOutName + "_snap_19999.persist";
		String clustersName = mainOutName + "_clusters.txt";
		String clusterExpName = mainOutName + "_clusters_expr.txt";
		GaussianDP DP = null;
		
		DP = GaussianDP.restoreFromFile(snapShotName);
		DP.normPairProbs(DP.numGoodIters);
		DP.concensusModel();
		
		try {
			DP.writeClusterStatsToFile(clusterExpName);
			DP.writeClustersToFile(clustersName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static MicroArrayData processData() {
		String dirName = "C:\\research_data\\macklis\\";
		String intensity_fName = "macklis_intensity.txt";
		String pvals_fName = "macklis_pvals.txt";
		String normalized_fName = "macklis_normalized_expression.txt";
		double pvalThresh = 0.0001;
		
		MicroArrayData intensityData = new MicroArrayData();
		MicroArrayData pvalData = new MicroArrayData();
		
		try {
			intensityData.readFile(dirName + intensity_fName);
			pvalData.readFile(dirName + pvals_fName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		int j = 0;
		
		// floor all intensities to 10
		for (i=0;i<intensityData.numRows;i++) {
			for (j=0;j<intensityData.numCols;j++) {
				if (intensityData.values[i][j] <= 0.0) {
					intensityData.values[i][j] = 10.0;
				}
			}
		}
		
		transformRelativeAbundance(intensityData.values);
		combineReplicateProbes(intensityData,pvalData);
		
		int[] filterList = new int[intensityData.numRows];
		int numFiltered = 0;
		int numExpressed = 0;
		int numChange = 0;
		
		// filter genes that are expressed everywhere
		for (i=0;i<intensityData.numRows;i++) {
			numExpressed = 0;
			numChange = 0;
			for (j=0;j<intensityData.numCols;j++) {
				if (pvalData.values[i][j] < pvalThresh) {
					numExpressed++;
				}
				if ((intensityData.values[i][j] >= 2.0 || intensityData.values[i][j] <= 0.5)) {
					numChange++; 
				}
			}
			
			if (numExpressed < 2 || numChange < 2) {
				filterList[i] = 1;
				numFiltered++;
			}
		}
		
		int k = 0;
		if (numFiltered > 0) {
			k = 0;
			double[][] intensity2 = new double[intensityData.numRows-numFiltered][intensityData.numCols];
			double[][] pvals2 = new double[intensityData.numRows-numFiltered][intensityData.numCols];
			String[] geneNames2 = new String[intensityData.numRows-numFiltered];
			for (i=0;i<filterList.length;i++) {
				if (filterList[i] == 0) {
					intensity2[k] = intensityData.values[i];
					pvals2[k] = pvalData.values[i];
					geneNames2[k] = intensityData.geneNames[i];
					k++;
				}
			}
			intensityData.values = intensity2;	intensityData.numRows = intensity2.length;	intensityData.geneNames = geneNames2;
			pvalData.values = pvals2;	pvalData.numRows = pvals2.length;	pvalData.geneNames = geneNames2;
		}
		
		for (i=0;i<intensityData.numRows;i++) {
			for (j=0;j<intensityData.numCols;j++) {
				intensityData.values[i][j] = Math.log(intensityData.values[i][j]);
			}
		}
		
		intensityData.normRows();
		
		try {
			intensityData.writeFile(dirName + normalized_fName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		return intensityData;
	}
	
	public static int[] uniqueGenes(String[] geneNames) {
		GoOBOReader goReader = new GoOBOReader();
		AffyGoAssocReader assocReader = new AffyGoAssocReader();
		OntDAG DAG = null;
		HashSet<String> useGenes = new HashSet<String>();
		OntGeneAssociations assoc = null;
		
		int i = 0;
		for (i=0;i<geneNames.length;i++) {
			useGenes.add(geneNames[i]);
		}
		
		try {
			DAG = goReader.readFile();
			assoc = assocReader.readFile(DAG,useGenes);
			
		} catch(IOException e) {
			System.out.println(e);
		} catch(DAGException e) {
			System.out.println(e);
		}
		
		LinkedHashMap<String,ArrayList<Integer>> gmap = assocReader.uniqueGeneMap(assoc,geneNames);
		int[] ugenes = new int[geneNames.length];
		OntGene gene = null;
		ArrayList<Integer> geneList = null;
		for (i=0;i<geneNames.length;i++) {
			gene = assoc.getGene(geneNames[i]);
			geneList = gmap.get(gene.REFSEQ);
			if (geneList.size() <= 1) {
				ugenes[i] = 1;
			}
		}
		return ugenes;
	}
	
	public static void combineReplicateProbes(MicroArrayData intensityValues,MicroArrayData pvalsData) {
		GoOBOReader goReader = new GoOBOReader();
		AffyGoAssocReader assocReader = new AffyGoAssocReader();
		OntDAG DAG = null;
		HashSet<String> useGenes = new HashSet<String>();
		OntGeneAssociations assoc = null;
		
		int i = 0;
		for (i=0;i<intensityValues.geneNames.length;i++) {
			useGenes.add(intensityValues.geneNames[i]);
		}
		
		try {
			DAG = goReader.readFile();
			assoc = assocReader.readFile(DAG,useGenes);
			
		} catch(IOException e) {
			System.out.println(e);
		} catch(DAGException e) {
			System.out.println(e);
		}
		
		LinkedHashMap<String,ArrayList<Integer>> gmap = assocReader.uniqueGeneMap(assoc,intensityValues.geneNames);
		double[][] values = new double[gmap.size()][intensityValues.values[0].length];
		double[][] pvals = new double[gmap.size()][intensityValues.values[0].length];
		String[] geneNames = new String[gmap.size()];
		String geneName;
		Iterator<String> iter = gmap.keySet().iterator();
		
		double[] mu = new double[intensityValues.values[0].length];
		double[] norm = new double[intensityValues.values[0].length];
		double[] pv = new double[intensityValues.values[0].length];
		ArrayList<Integer> genes = null;
		int j = 0;
		int k = 0;
		int idx = 0;
		while(iter.hasNext()) {
			geneName = iter.next();
			geneNames[k] = geneName;
			genes = gmap.get(geneName);
			for (j=0;j<intensityValues.numCols;j++) {
				mu[j] = 0.0;
				norm[j] = 0.0;
				pv[j] = 100.0;
			}
			for (idx=0;idx<genes.size();idx++) {
				i = genes.get(idx);
				for (j=0;j<intensityValues.numCols;j++) {
					if (pvalsData.values[i][j] < pv[j]) {
						pv[j] = pvalsData.values[i][j];
					}
					if (!Double.isNaN(intensityValues.values[i][j])) {
						mu[j] += intensityValues.values[i][j];
						norm[j]++;
					}
				}
			}
			for (j=0;j<intensityValues.numCols;j++) {
				values[k][j] = mu[j]/norm[j];
			}
			
			k++;
		}
		
		intensityValues.values = values;
		intensityValues.geneNames = geneNames;
		intensityValues.numRows = values.length;
		
		pvalsData.values = pvals;
		pvalsData.geneNames = geneNames;
		pvalsData.numRows = values.length;
	}
	
	public static void transformRelativeAbundance(double[][] expression) {
		double norm = 0.0;
		int i = 0;
		int j = 0;
		
		for (i=0;i<expression.length;i++) {
			norm = 0.0;
			for (j=0;j<expression[i].length;j++) {
				norm = norm + expression[i][j];
			}
			for (j=0;j<expression[0].length;j++) {
				expression[i][j] = ((double) expression[0].length)*expression[i][j]/norm;
			}
		}
	}
	
	public static void displayClusters() {
		String mainOutName = "C:\\research_data\\macklis\\macklis_IGMM";
		String clustersName = mainOutName + "_clusters.txt";
		String clusterExpName = mainOutName + "_clusters_expr.txt";
		String GOName = mainOutName + "_clusters_GO.txt";
		String heatMapName = mainOutName + "_cluster_means.svg";
		
		MicroArrayData clusterAssignData = new MicroArrayData();
		clusterAssignData.setDiscrete();
		MicroArrayData clusterExprData = new MicroArrayData();
		
		try {
			clusterAssignData.readFile(clustersName);
			clusterExprData.readFile(clusterExpName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		GoOBOReader goReader = new GoOBOReader();
		AffyGoAssocReader assocReader = new AffyGoAssocReader();
		assocReader.setUseREFSEQ(true);
		
		int[] clusterAssigns = new int[clusterAssignData.numRows];
		int i = 0;
		for (i=0;i<clusterAssigns.length;i++) {
			clusterAssigns[i] = clusterAssignData.dvalues[i][0];
		}
		
		IGMMClusterSummarizer summarizer = new IGMMClusterSummarizer(goReader,assocReader,clusterAssignData.geneNames);
		summarizer.outputClusters(GOName,heatMapName,clusterAssigns,clusterAssignData.geneNames,clusterExprData.experimentNames,clusterExprData.values);
	}
}
