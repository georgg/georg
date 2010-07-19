package edu.mit.csail.psrg.georg.BiclusterAnalysis;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;

import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;
import edu.mit.csail.psrg.georg.StatUtil.StatUtil;
import edu.mit.csail.psrg.georg.StatUtil.StepDownFDR;
import edu.mit.csail.psrg.georg.StatUtil.VectorUtil;

// compare two different clusterings by using one as the "gold standard"
// and testing the other for significance of enrichments
public class CompareClusterings {
	
	public static void main(String[] args) {
				
		String clusterFile1 = "c:\\research_data\\Shyamsundar\\nmf\\MHS_nmf_bicluster_genes";
		String clusterFile2 = "c:\\research_data\\Shyamsundar\\nmf\\SHY_nmf_bicluster_genes";
		
	//	String clusterFile1 = "c:\\research_data\\Shyamsundar\\MHS_2lvls_UP\\MHS_2lvls_UP_bicluster_genes";
	//	String clusterFile2 = "c:\\research_data\\Shyamsundar\\SHY_2lvls\\SHY_2lvls_bicluster_genes";
		
		scoreClusterings(clusterFile1,clusterFile2);
		scoreClusterings(clusterFile2,clusterFile1);
	}
	
	public static void scoreClusterings(String clusterFile1,String clusterFile2) {
		
		ClusterReader clusters2 = new ClusterReader();
		ClusterReader clusters1 = new ClusterReader();
		
		int minGenes = 25;
		int maxGenes = 2000;
		
		try {
			clusters2.readFile(clusterFile2+".txt");
			clusters1.readFile(clusterFile1+".txt");
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		int j = 0;
		ArrayList<ArrayList<String>> clusters1New = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> clusters2New = new ArrayList<ArrayList<String>>();
		
		for (i=0;i<clusters1.clusters.size();i++) {
			if (clusters1.clusters.get(i).size() >= minGenes & clusters1.clusters.get(i).size() <= maxGenes) {
				clusters1New.add(clusters1.clusters.get(i));
			}
		}
		clusters1.clusters = clusters1New;
		
		for (i=0;i<clusters2.clusters.size();i++) {
			if (clusters2.clusters.get(i).size() >= minGenes & clusters2.clusters.get(i).size() <= maxGenes) {
				clusters2New.add(clusters2.clusters.get(i));
			}
		}
		clusters2.clusters = clusters2New;
		
		int[] classCount = new int[clusters2.clusters.size()];
		int[] totalClassCount = new int[clusters2.clusters.size()];
		double v = 0.0;
		double[] pvals = new double[clusters2.clusters.size()];
		int totalGenes = 0;
		double min_pval = 10e-45;
		HashMap<String,Integer> classMap = new HashMap<String,Integer>();
		ArrayList<String> cluster = null;
		HashSet<String> uniqueGenes = new HashSet<String>();
		
		int classID = 0;
		int mostAbundant = 0;
		double bestPval = 1.0;
		
		for(i=0;i<clusters2.clusters.size();i++) {
			cluster = clusters2.clusters.get(i);
			totalClassCount[i] = cluster.size();
			totalGenes += totalClassCount[i];
			for (j=0;j<cluster.size();j++) {
				classMap.put(cluster.get(j),i);
			}
		}
		
		for (i=0;i<clusters1.clusters.size();i++) {
			cluster = clusters1.clusters.get(i);
			for (j=0;j<cluster.size();j++) {
				uniqueGenes.add(cluster.get(j));
			}
		}
		
		Vector<Double> ap = new Vector<Double>();
		Vector<Integer> ab = new Vector<Integer>();
		
		try {
			FileWriter file = new FileWriter(clusterFile1+"_compare_cumscore.txt");
			file.write((new Integer(uniqueGenes.size())).toString());
			file.write("\n");
			for (i=0;i<clusters1.clusters.size();i++) {
				cluster = clusters1.clusters.get(i);
				for (j=0;j<classCount.length;j++) {
					pvals[j] = 0.0;
					classCount[j] = 0;
				}
				for (j=0;j<cluster.size();j++) {
					uniqueGenes.add(cluster.get(j));
					if (classMap.containsKey(cluster.get(j))) {
						classID = classMap.get(cluster.get(j));
						classCount[classID]++;
					}
				}
				for (j=0;j<classCount.length;j++) {
					v = 1.0 - StatUtil.hyperGeometricCDF(classCount[j]-1,totalGenes,totalClassCount[j],cluster.size());
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
		boolean[] signifBicluster = new boolean[clusters1.clusters.size()];
		
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
			FileWriter clusterSigFile = new FileWriter(clusterFile1 + "_compare_signif.txt");
			for (i=0;i<signifBicluster.length;i++) {
				clusterSigFile.write((new Integer(clusters1.clusters.get(i).size())).toString());
				clusterSigFile.write("\t");
				if (signifBicluster[i]) {
					clusterSigFile.write("1\n");
				} else {
					clusterSigFile.write("0\n");
				}
				if (signifBicluster[i])
					totalSignif++;
			}
			clusterSigFile.close();
		} catch(IOException e) {
			System.out.println(e);
		}
		
		System.out.println("Significant clusters = " + totalSignif);
		System.out.println("Total clusters = " + clusters1.clusters.size());
	}
}
