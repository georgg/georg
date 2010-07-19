package edu.mit.csail.psrg.georg.BiclusterAnalysis;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Vector;

import edu.mit.csail.psrg.georg.Annotation.Annotations;
import edu.mit.csail.psrg.georg.Annotation.DAGException;
import edu.mit.csail.psrg.georg.Annotation.GeneSetCollection;
import edu.mit.csail.psrg.georg.Annotation.OntGene;
import edu.mit.csail.psrg.georg.Annotation.OntTermScore;
import edu.mit.csail.psrg.georg.Annotation.OntTerm;
import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;
import edu.mit.csail.psrg.georg.StatUtil.StatUtil;
import edu.mit.csail.psrg.georg.StatUtil.StepDownFDR;
import edu.mit.csail.psrg.georg.StatUtil.VectorUtil;

public class ScoreBiclusters {
	//	 parameters for GO categories
	public int minGOGenes = 5;
	public int maxGOGenes = 200;
	public double FDR = 0.05;
	public String namespace = "biological_process";
	public HashSet<String> useGenes = new HashSet<String>();
	public String GeneAssocFileName = "C:\\research_data\\mouse_human\\b47mm7\\Hs.genelist.go.plus";
	public String GOHierarchyFileName = "C:\\research_data\\go_data\\gene_ontology_plus.obo";
	
	public static void main(String[] args) {
		ScoreBiclusters me = new ScoreBiclusters();
		me.GeneAssocFileName = "c:\\research_data\\Shyamsundar\\gene_annotations.txt";
		
	//	me.GeneAssocFileName = "c:\\research_data\\mouse_human2\\gene_annotations.txt";
		
		String baseFile = "c:\\research_data\\Shyamsundar\\MHS_1lvls_UP\\MHS_1lvls_UP";
	//	String baseFile = "c:\\research_data\\mouse_human2\\UD2lvls\\MHUD2lvls";
		
		me.scoreGOCategories(baseFile);
		me.scoreTissues(baseFile,"c:\\research_data\\Shyamsundar\\MHS_tissue_classifications.txt");
	//	me.scoreTissues(baseFile,"c:\\research_data\\Shyamsundar\\SHY_tissue_classifications.txt");
	}
	
	public void scoreGOCategories(String baseFileName) {
		Annotations myAnnotations = new Annotations();
		
		ClusterReader genesFile = new ClusterReader();
		try {
			genesFile.readFile(baseFileName + "_bicluster_genes.txt");
		} catch(IOException e) {
			System.out.println(e);
		}
		
		HashSet<String> useGenes = new HashSet<String>();
		int i = 0;
		int j = 0;
		ArrayList<String> genes = null;
		for (i=0;i<genesFile.clusters.size();i++) {
			genes = genesFile.clusters.get(i);
			for (j=0;j<genes.size();j++) {
				useGenes.add(genes.get(j));
			}
		}
		
		try {
			myAnnotations.readFiles(GeneAssocFileName,GOHierarchyFileName,useGenes);
			System.out.println("Loaded annotations");
		} catch(IOException e) {
			System.out.println(e);
		} catch(DAGException e) {
			System.out.println(e);
		}
		
		GeneSetCollection geneSets = new GeneSetCollection(myAnnotations,minGOGenes,maxGOGenes,FDR,namespace);
		OntGene oGene = null;
		
		for (i=0;i<genesFile.clusters.size();i++) {
			genes = genesFile.clusters.get(i);
			geneSets.addGeneSet(genes);
		}
		
		geneSets.geneSetSignif();
		int[] numSignif = geneSets.numSignificantTerms();
		
		try {
			FileWriter signifCountsFile = new FileWriter(baseFileName + "_bicluster_genes_signif_counts.txt");
			for (i=0;i<numSignif.length;i++) {
				signifCountsFile.write((new Integer(numSignif[i])).toString());
				signifCountsFile.write("\n");
			}
			signifCountsFile.close();
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int mostAbundant = 0;
		double bestPVal = 1.0;
		double[] bestScores = new double[geneSets.geneSets.size()];
		double minPVal = 10e-45;
		Iterator<OntTermScore> scoreIter = null;
		Iterator<OntTerm> termIter = null;
		OntTerm term = null;
		OntTermScore score = null;
		boolean isSignif = false;
		
		for (i=0;i<geneSets.geneSets.size();i++) {
			mostAbundant = 0;
			bestPVal = 1.0;
			isSignif = false;
			scoreIter = geneSets.geneSets.get(i).termMap.values().iterator();
			termIter = geneSets.geneSets.get(i).termMap.keySet().iterator();
			while(scoreIter.hasNext()) {
				score = scoreIter.next();
				term = termIter.next();
				
				if (score.genes.size() >= mostAbundant & score.acceptPVal) {
					if (score.genes.size() == mostAbundant) {
						if (score.pval < bestPVal) {
							mostAbundant = score.genes.size();
							bestPVal = score.pval;
							isSignif = score.acceptPVal;
						}
					} else {
						mostAbundant = score.genes.size();
						bestPVal = score.pval;
						isSignif = score.acceptPVal;
					}	
				}
			}
			if (bestPVal < minPVal) {
				bestPVal = minPVal;
			}	
			bestScores[i] = Math.log10(bestPVal);
		}
		
		try {
			FileWriter gf = new FileWriter(baseFileName + "_bicluster_genes_cumsum.txt");
			gf.write((new Integer(useGenes.size())).toString());
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
	
	public void scoreTissues(String baseFile,String tissueClassFile) {
		
		ClusterReader tissueClasses = new ClusterReader();
		ClusterReader biclusters = new ClusterReader();
		
		try {
			tissueClasses.readFile(tissueClassFile);
			biclusters.readFile(baseFile+"_bicluster_tissues.txt");
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
			FileWriter file = new FileWriter(baseFile+"_bicluster_tissues_cumscore.txt");
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
					try {
						v = 1.0 - StatUtil.hyperGeometricCDF(classCount[j]-1,totalTissues,totalClassCount[j],cluster.size());
					} catch(Exception e) {
						v = 1.0;
					}
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
			FileWriter tissueSigFile = new FileWriter(baseFile + "_bicluster_tissues_signif.txt");
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
	
	public String removeSpaces(String s) {
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
}
