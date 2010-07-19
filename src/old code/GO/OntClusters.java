package edu.mit.csail.psrg.georg.GO;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import edu.mit.csail.psrg.georg.StatUtil.StepDownFDR;

public class OntClusters {
	public ArrayList<OntTermScores> clusters = new ArrayList<OntTermScores>();
	int minGenes = 5;
	int maxGenes = 200;
	double FDR = 0.05;
	int namespace = OntTerm.AnyCategory;
	OntGeneAssociations assoc = null;
	
	public OntClusters(OntGeneAssociations i_assoc) { assoc = i_assoc; }
	
	public OntClusters(OntGeneAssociations i_assoc,int i_minGenes, int i_maxGenes, double i_FDR, int i_namespace) {
		minGenes = i_minGenes;
		maxGenes = i_maxGenes;
		FDR = i_FDR;
		namespace = i_namespace;
		assoc = i_assoc;
	}
	
	public ArrayList<OntTermScores> getClusters() {
		return clusters;
	}
	
	public void addCluster(ArrayList<String> geneIDs, double[] geneRanks,String description,double p) {
		OntTermScores scores = new OntTermScores();
		scores.description = description;
		scores.probability = p;
		clusters.add(scores);
		scores.scoreGeneSet(geneIDs,geneRanks,assoc,minGenes,maxGenes,FDR,namespace);
	}
	
	public void addCluster(ArrayList<String> geneIDs, ArrayList<Double> gr,String description,double p) {
		OntTermScores scores = new OntTermScores();
		scores.description = description;
		scores.probability = p;
		clusters.add(scores);
		
		double[] geneRanks = new double[gr.size()];
		int i = 0;
		for (i=0;i<geneRanks.length;i++) {
			geneRanks[i] = gr.get(i);
		}
		
		scores.scoreGeneSet(geneIDs,geneRanks,assoc,minGenes,maxGenes,FDR,namespace);
	}
	
	public void addCluster(ArrayList<String> geneIDs, double[] geneRanks) {
		addCluster(geneIDs,geneRanks,"",0.0);
	}
	
	// determine whether enrichments are significant, correcting for multiple
	// hypotheses across all clusters
	public void clusterSignif() {
		int i = 0;
		int totalScores = 0;
		
		for (i=0;i<clusters.size();i++) {
			totalScores += clusters.get(i).termMap.size();
		}
		
		OntScore[] scores = new OntScore[totalScores];
		Iterator<OntScore> iter = null;
		int j = 0;
		for (i=0;i<clusters.size();i++) {
			iter = clusters.get(i).termMap.values().iterator();
			while(iter.hasNext()) {
				scores[j] = iter.next();
				j++;
			}
		}
		
		i = 0;
		for (i=0;i<scores.length;i++) {
			scores[i].rank = i;
		}
		
		List<OntScore> scoreList = Arrays.asList(scores);
		Collections.sort(scoreList);
		
		boolean[] acceptPvals = null;
		double[] pvals = new double[scores.length];
		OntScore score = null;
		for (i=0;i<scores.length;i++) {
			score = scoreList.get(i);
			pvals[i] = score.pval;
		}
		acceptPvals = StepDownFDR.correctPvals(pvals,FDR);
		for (i=0;i<scores.length;i++) {
			score = scoreList.get(i);
			score.acceptPVal = acceptPvals[i];
		}
	}
	
	public int[] numSignificantTerms() {
		int[] sigTerms = new int[clusters.size()];
		int i = 0;
		for (i=0;i<clusters.size();i++) {
			sigTerms[i] = clusters.get(i).numSignificantTerms();
		}
		return sigTerms;
	}
	
	public void outputClustersLong(String fName) throws IOException {
		FileWriter outFile = new FileWriter(fName);
		int i = 0;
		for (i=0;i<clusters.size();i++) {
			outFile.write("Cluster#");
			outFile.write((new Integer(i+1)).toString());
			outFile.write("\t");
			outFile.write((new Integer(clusters.get(i).numGenesInSet)).toString());
			outFile.write("\t");
			outFile.write((new Double(clusters.get(i).probability)).toString());
			outFile.write("\n");
			if (clusters.get(i).description != null) {
				outFile.write(clusters.get(i).description);
				outFile.write("\n");
			}
			clusters.get(i).outputScoresLong(outFile,true);
		}
		outFile.close();
	}
	
	public void outputClustersShort(String fName) throws IOException {
		FileWriter outFile = new FileWriter(fName);
		int i = 0;
		for (i=0;i<clusters.size();i++) {
			if (clusters.get(i).description != null) {
				outFile.write(clusters.get(i).description);
				outFile.write("\t");
			}
			outFile.write((new Double(clusters.get(i).probability)).toString());
			outFile.write("\t");
			clusters.get(i).outputBestScore(outFile);
		}
		outFile.close();
	}
}
