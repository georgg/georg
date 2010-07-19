package edu.mit.csail.psrg.georg.Annotation;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import edu.mit.csail.psrg.georg.StatUtil.StepDownFDR;

public class GeneSetCollection {
	public ArrayList<GeneSet> geneSets = new ArrayList<GeneSet>();
	int minGenes = 5;
	int maxGenes = 200;
	double FDR = 0.05;
	String namespace = "AnyCategory";
	Annotations myAnnotations = null;
	
	public GeneSetCollection(Annotations a) { myAnnotations = a; }
	
	public GeneSetCollection(Annotations a,int i_minGenes, int i_maxGenes, double i_FDR, String i_namespace) {
		minGenes = i_minGenes;
		maxGenes = i_maxGenes;
		FDR = i_FDR;
		namespace = i_namespace;
		myAnnotations = a;
	}
	
	public ArrayList<GeneSet> getGeneSets() {
		return geneSets;
	}
	
	public void addGeneSet(ArrayList<String> geneIDs) {
		double[] ranks = new double[geneIDs.size()];
		addGeneSet(geneIDs,ranks);
	}
	
	public void addGeneSet(ArrayList<String> geneIDs, double[] geneRanks) {
		GeneSet set = new GeneSet();
		geneSets.add(set);
		set.scoreGeneSet(geneIDs,geneRanks,myAnnotations,minGenes,maxGenes,FDR,namespace);
	}
	
	public void addGeneSet(ArrayList<String> geneIDs, ArrayList<Double> gr) {
		GeneSet set = new GeneSet();
		geneSets.add(set);
		
		double[] geneRanks = new double[gr.size()];
		int i = 0;
		for (i=0;i<geneRanks.length;i++) {
			geneRanks[i] = gr.get(i);
		}
		
		set.scoreGeneSet(geneIDs,geneRanks,myAnnotations,minGenes,maxGenes,FDR,namespace);
	}
	
	// determine whether enrichments are significant, correcting for multiple
	// hypotheses across all clusters
	public void geneSetSignif() {
		int i = 0;
		int totalTerms = 0;
		
		for (i=0;i<geneSets.size();i++) {
			totalTerms += geneSets.get(i).termMap.size();
		}
		
		OntTermScore[] scores = new OntTermScore[totalTerms];
		Iterator<OntTermScore> iter = null;
		int j = 0;
		for (i=0;i<geneSets.size();i++) {
			iter = geneSets.get(i).termMap.values().iterator();
			while(iter.hasNext()) {
				scores[j] = iter.next();
				j++;
			}
		}
		
		i = 0;
		for (i=0;i<scores.length;i++) {
			scores[i].rank = i;
		}
		
		List<OntTermScore> scoreList = Arrays.asList(scores);
		Collections.sort(scoreList);
		
		boolean[] acceptPvals = null;
		double[] pvals = new double[scores.length];
		OntTermScore score = null;
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
		int[] sigTerms = new int[geneSets.size()];
		int i = 0;
		for (i=0;i<geneSets.size();i++) {
			sigTerms[i] = geneSets.get(i).numSignificantTerms();
		}
		return sigTerms;
	}
}
