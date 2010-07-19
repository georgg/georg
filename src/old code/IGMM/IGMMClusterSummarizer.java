package edu.mit.csail.psrg.georg.IGMM;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Vector;

import edu.mit.csail.psrg.georg.GO.DAGException;
import edu.mit.csail.psrg.georg.GO.DAGReader;
import edu.mit.csail.psrg.georg.GO.OntAssociationReader;
import edu.mit.csail.psrg.georg.GO.OntClusters;
import edu.mit.csail.psrg.georg.GO.OntDAG;
import edu.mit.csail.psrg.georg.GO.OntGeneAssociations;
import edu.mit.csail.psrg.georg.GO.OntTerm;
import edu.mit.csail.psrg.georg.Graph.ExpressionHeatMap;
import edu.mit.csail.psrg.georg.HierarchicalCluster.HierarchicalCluster;
import edu.mit.csail.psrg.georg.StatUtil.VectorUtil;

public class IGMMClusterSummarizer {
	OntGeneAssociations assoc = null;
	OntDAG DAG = null;
	int minGenes = 5;
	int maxGenes = 200;
	// minimum count of gene in topic
	double FDR = 0.05;
//	int namespace = OntTerm.AnyCategory;
	int namespace = OntTerm.BiologicalProcess;
	HashSet<String> useGenes = new HashSet<String>();
	
	public IGMMClusterSummarizer(DAGReader init_DAGReader, OntAssociationReader init_ontAssociationReader,String[] geneNames) {
		processGenes(geneNames);
		readDAG(init_DAGReader);
		readGeneAssocations(init_ontAssociationReader);
	}
	
	void processGenes(String[] geneNames) {
		int i = 0;
		for (i=0;i<geneNames.length;i++) {
			useGenes.add(geneNames[i]);
		}
	}
	
	void readDAG(DAGReader reader) {
		try {
				DAG = reader.readFile();
				System.out.println("Loaded DAG");
			} catch(IOException e) {
				System.out.println(e);
			} catch(DAGException e) {
				System.out.println(e);
			}
	}
	
	void readGeneAssocations(OntAssociationReader reader) {
		try {
			assoc = reader.readFile(DAG,useGenes);
			System.out.println("Loaded associations");
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public void outputClusters(String ontOutName,String heatMapOutName,int[] clusterAssigns,String[] geneNames,String[] experimentNames,double[][] mu) {
		ArrayList<String> genes = new ArrayList<String>();
		OntClusters clusters = new OntClusters(assoc,minGenes,maxGenes,FDR,namespace);
		double[] ranks = null;
		
		ArrayList<ArrayList<Integer>> clusterMap = new ArrayList<ArrayList<Integer>>();
		int i = 0;
		int j = 0;
		int k = 0;
		
		int numClusters = -1;
		for (i=0;i<clusterAssigns.length;i++) {
			if (clusterAssigns[i] > numClusters) {
				numClusters = clusterAssigns[i];
			}
		}
		numClusters++;
		
		ArrayList<Integer> geneList = null;
		for (i=0;i<numClusters;i++) {
			geneList = new ArrayList<Integer>();
			clusterMap.add(geneList);
		}
		
		
		for (j=0;j<clusterAssigns.length;j++) {
			i = clusterAssigns[j];
			geneList = clusterMap.get(i);
			geneList.add(j);
		}
		
		int[] useCluster = new int[numClusters];
		int numUseClusters = 0;
		int geneIDX = 0;
		ArrayList<ArrayList<Integer>> clusterMap2 = new ArrayList<ArrayList<Integer>>();
		for (i=0;i<numClusters;i++) {
			geneList = clusterMap.get(i);
			if (geneList.size() >= minGenes) {
				clusterMap2.add(geneList);
				numUseClusters++;
				useCluster[i] = 1;
			}
		}
		clusterMap = clusterMap2;
			
		double[][] mu2 = new double[numUseClusters][mu[0].length-1];
		double[] std = new double[numUseClusters];
		String[] cn = new String[numUseClusters];
		String[] en = new String[experimentNames.length-1];
		k = 0;
		
		for (i=0;i<numClusters;i++) {
			if (useCluster[i] > 0) {
				for (j=1;j<mu[0].length;j++) {
					mu2[k][j-1] = mu[i][j];
					std[k] = mu[i][0];
				}
				cn[k] = (new Integer(k+1)).toString();
				k++;
			}
		}
		for (j=1;j<experimentNames.length;j++) {
			en[j-1] = experimentNames[j];
		}
		
		ArrayList<Integer> order = orderDataRows(mu2);
		double[][] mu3 = new double[mu2.length][mu2[0].length];
		double[] std2 = new double[numUseClusters];
				
		for (i=0;i<mu3.length;i++) {
			geneIDX = order.get(i);
			System.arraycopy(mu2[i],0,mu3[geneIDX],0,mu2[i].length);
			std2[geneIDX] = std[i];
		}
		mu2 = mu3;
		std = std2;
		
		ExpressionHeatMap hm = new ExpressionHeatMap(mu2,cn,en);
		hm.outputHeatMap(heatMapOutName);
		
		ArrayList<Integer> order2 = new ArrayList<Integer>();
		for (i=0;i<order.size();i++) {
			order2.add(0);
		}
		for (i=0;i<order.size();i++) {
			geneIDX = order.get(i);
			order2.set(geneIDX,i);
		}
		
		String muString = "";
		int qq = 0;
		for (qq=0;qq<clusterMap.size();qq++) {
			i = order2.get(qq);
			geneList = clusterMap.get(i);
			genes.clear();
			for (j=0;j<geneList.size();j++) {
				geneIDX = geneList.get(j);
				genes.add(geneNames[geneIDX]);
			}
			muString = "";
			for (k=0;k<mu2[i].length;k++) {
				muString = muString + (new Double(mu2[qq][k])).toString();
				if (k < mu[i].length - 1) {
					muString += "\t";
				}
			}
			clusters.addCluster(genes,ranks,muString,std[qq]);
		}
		
		clusters.clusterSignif();
		
		try {
			clusters.outputClustersLong(ontOutName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public ArrayList<Integer> orderDataRows(double[][] data) {
		int nr = data.length;
		double[][] distances = new double[nr-1][];
		int i = 0;
		int j = 0;
		int k = 0;
		Vector<Double> x = new Vector<Double>();
		Vector<Double> y = new Vector<Double>();
		double v = 0.0;
		double maxDistance = 0.0;
		
		for (i=0;i<nr-1;i++) {
			x.clear();
			for (k=0;k<data[i].length;k++) {
				x.add(data[i][k]);
			}
			distances[i] = new double[nr-i-1];
			for (j=(i+1);j<nr;j++) {
				y.clear();
				for (k=0;k<data[j].length;k++) {
					y.add(data[j][k]);
				}
				v = (1.0-VectorUtil.correlation(x,y))/2.0;
			//	v = VectorUtil.euclidean(x,y);
				distances[i][j-i-1] = v;
				if (v > maxDistance) {
					maxDistance = v;
				}
			}
		}
		HierarchicalCluster hierCluster = new HierarchicalCluster();
		hierCluster.cluster(distances,maxDistance + 1.0);
		
		ArrayList<Integer> rowPosition = new ArrayList<Integer>();
		for (i=0;i<nr;i++) {
			rowPosition.add(0);
		}
		for (i=0;i<hierCluster.clusters.get(0).size();i++) {
			rowPosition.set(hierCluster.clusters.get(0).get(i),i);
		}
		
		return rowPosition;
	}
}
