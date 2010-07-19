package edu.mit.csail.psrg.georg.Annotation;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;

import edu.mit.csail.psrg.georg.StatUtil.*;

public class GeneSet {
	public LinkedHashMap<OntTerm,OntTermScore> termMap = new LinkedHashMap<OntTerm,OntTermScore>();
	public int numGenesInSet = 0;
	int numGenesInPopulation = 0;
	ArrayList<OntGene> genes = null;
	double[] geneRanks = null;
	static String outSep1 = "--------------------------------\n";
	static String outSep2 = "================================\n";
	
	void addGene(OntTerm term,OntGene gene) {
		OntTermScore s = null;
		
		if (termMap.containsKey(term)) {
			s = termMap.get(term);
		} else {
			s = new OntTermScore();
			termMap.put(term,s);
		}
		s.addGene(gene);
	}
	
	public ArrayList<OntGene> getGenes() {
		return genes;
	}
	
	public double[] getGeneRanks() {
		return geneRanks;
	}
	
	public void scoreGeneSet(ArrayList<String> geneIDs,double[] gr, Annotations myAnnotations,int minGenes,int maxGenes,double FDR,String namespace) {
		int i = 0;
		int j = 0;
		HashSet<OntTerm> parents = new HashSet<OntTerm>();
		HashSet<OntTerm> terms = null;
		OntTerm term = null;
		OntGene gene = null;
		Iterator<OntTerm> iter = null;
		int totalGenes = 0;
		genes = new ArrayList<OntGene>();
		geneRanks = gr;
		
		for (i=0;i<geneIDs.size();i++) {
			gene = myAnnotations.getGene(geneIDs.get(i));
			if (gene == null) {
				gene = new OntGene();
				gene.Name = geneIDs.get(i);
				gene.ID = geneIDs.get(i);
			//	totalGenes++;
			}
			if (gene != null) {
				genes.add(gene);
				terms = gene.getTerms();
				parents.clear();
				if (terms != null) {
					iter = terms.iterator();
					while(iter.hasNext()) {
						term = iter.next();
						term.getAllParents(parents,minGenes,maxGenes,namespace);
					}
					iter = parents.iterator();
					while(iter.hasNext()) {
						term = iter.next();
					//	System.out.println(term.ID);
						addGene(term,gene);
					}
					totalGenes++;
				}
			}
		}
		
		OntTerm[] terms2 = new OntTerm[termMap.size()];
		termMap.keySet().toArray(terms2);
		OntTermScore score = null;
		for(i=0;i<terms2.length;i++) {
			score = termMap.get(terms2[i]);
			if (score != null) {
				if (score.genes.size() < minGenes) {
					termMap.remove(terms2[i]);
				}
			}
		}
		
		numGenesInSet = totalGenes;
		numGenesInPopulation = myAnnotations.genes.size();
		computePvals(FDR);
	}
	
	void computePvals(double FDR) {
		Iterator<OntTerm> iter = termMap.keySet().iterator();
		OntTermScore score = null;
		OntTerm term = null;
		while(iter.hasNext()) {
			term = iter.next();
			score = termMap.get(term);
			score.pval = 1.0-StatUtil.hyperGeometricCDF(score.genes.size()-1,numGenesInPopulation,term.numGenes,numGenesInSet);
		//	score.pval = 1.0-StatUtil.hyperGeometricCDF(score.genes.size(),numGenesInPopulation,term.numGenes,numGenesInSet);
			score.fraction = ((double) score.genes.size())/((double) numGenesInSet);
		}
		sortScores();
		// now correct the p-values
		Iterator<OntTermScore> iter2 = termMap.values().iterator();
		boolean[] acceptPvals = null;
		double[] pvals = new double[termMap.size()];
		int i = 0;
		while(iter2.hasNext()) {
			score = iter2.next();
			pvals[i] = score.pval;
			i++;
		}
		acceptPvals = StepDownFDR.correctPvals(pvals,FDR);
		iter2 = termMap.values().iterator();
		i = 0;
		while(iter2.hasNext()) {
			score = iter2.next();
			score.acceptPVal = acceptPvals[i];
			i++;
		}
	}
	
	void sortScores() {
		if (termMap.isEmpty()) {
			return;
		}
		OntTermScore[] scores = new OntTermScore[termMap.size()];
		OntTerm[] terms = new OntTerm[termMap.size()];
		termMap.values().toArray(scores);
		termMap.keySet().toArray(terms);
		int i = 0;
		for (i=0;i<scores.length;i++) {
			scores[i].rank = i;
		}
		
		List<OntTermScore> scoreList = Arrays.asList(scores);
		Collections.sort(scoreList);
		
		termMap.clear();
		OntTerm term = null;
		OntTermScore score = null;
		for(i=0;i<scores.length;i++) {
			score = scoreList.get(i);
			term = terms[score.rank];
			termMap.put(term,score);
		}
	}
	
	void outputBestScore(FileWriter outFile) throws IOException {
		Iterator<OntTerm> iter = termMap.keySet().iterator();
		OntTerm term = null;
		OntTermScore score = null;
		OntTerm bestTerm = null;
		OntTermScore bestScore = null;
		double bestFraction = 0.0;
		double bestPVal = 1.0;
		
		while(iter.hasNext()) {
			term = iter.next();
			score = termMap.get(term);
			if (score.acceptPVal) {
				if (score.fraction > bestFraction) {
			//	if (score.pval < bestPVal) {
					bestScore = score;
					bestFraction = score.fraction;
				//	bestPVal = score.pval;
					bestTerm = term;
				}
			}
		}
		
		if (bestScore != null) {
			outFile.write(bestTerm.ID);
			outFile.write("\t");
			outFile.write((new Double(bestScore.pval)).toString());
			outFile.write("\t");
			outFile.write((new Integer(bestScore.genes.size())).toString());
			outFile.write("\\");
			outFile.write((new Integer(numGenesInSet)).toString());
			outFile.write("\t");
			outFile.write(bestTerm.name);
			outFile.write("\t");
			outFile.write((new Integer(bestTerm.numGenes)).toString());
			outFile.write("\n");
		} else {
			outFile.write("none");
			outFile.write("\t");
			outFile.write("1.0");
			outFile.write("\t");
			outFile.write("0");
			outFile.write("\\");
			outFile.write((new Integer(numGenesInSet)).toString());
			outFile.write("\t");
			outFile.write("none");
			outFile.write("\t");
			outFile.write("0");
			outFile.write("\n");
		}
	}
	
	public int numSignificantTerms() {
		int numSigTerms = 0;
		Iterator<OntTerm> iter2 = termMap.keySet().iterator();
		OntTerm term = null;
		OntTermScore score = null;
		
		while(iter2.hasNext()) {
			term = iter2.next();
			score = termMap.get(term);
			if (score.acceptPVal) {
				numSigTerms++;
			}
		}
		
		return numSigTerms;
	}
	
	public void outputScoresOneLine(FileWriter outFile) throws IOException {
		Iterator<OntTerm> iter2 = termMap.keySet().iterator();
		OntTerm term = null;
		OntTermScore score = null;
		int numSig = 0;
		
		while(iter2.hasNext()) {
			term = iter2.next();
			score = termMap.get(term);
			if (score.acceptPVal) {
				numSig++;
			}
		}
		
		if (numSig == 0) {
			outFile.write("none\n");
			return;
		}
		
		iter2 = termMap.keySet().iterator();
		
		while(iter2.hasNext()) {
			term = iter2.next();
			score = termMap.get(term);
			if (score.acceptPVal) {
				outFile.write(term.name);
				outFile.write(" ");
				outFile.write(term.ID);
				outFile.write(", p=");
				outFile.write((new Double(score.pval)).toString());
				outFile.write(" (");
				outFile.write((new Integer(score.genes.size())).toString());
				outFile.write("\\");
				outFile.write((new Integer(numGenesInSet)).toString());
				outFile.write(")");
			//	outFile.write((new Integer(term.numGenes)).toString());
				outFile.write("\t");
			}
		}
		outFile.write("\n");
	}
	
	public void outputScoresLong(FileWriter outFile,boolean outputGeneInfo) throws IOException {
		Iterator<OntGene> iter = genes.iterator();
		OntGene gene = null;
		int i = 0;
		
		if (outputGeneInfo) {
			outFile.write(GeneSet.outSep1);
			while(iter.hasNext()) {
				gene = iter.next();
				if (geneRanks != null) {
					outFile.write((new Double(geneRanks[i]).toString()));
					outFile.write("\t");
				}
				outFile.write(gene.ID);
				outFile.write("\t");
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
				outFile.write("\n");
				i++;
			}
			outFile.write("\n");
		}
		
		outFile.write(GeneSet.outSep2);
		
		Iterator<OntTerm> iter2 = termMap.keySet().iterator();
		OntTerm term = null;
		OntTermScore score = null;
		
		while(iter2.hasNext()) {
			term = iter2.next();
			score = termMap.get(term);
			if (score.acceptPVal) {
				outFile.write(term.ID);
				outFile.write("\t");
				outFile.write((new Double(score.pval)).toString());
				outFile.write("\t");
				outFile.write((new Integer(score.genes.size())).toString());
				outFile.write("\\");
				outFile.write((new Integer(numGenesInSet)).toString());
				outFile.write("\t");
			//	outFile.write((new Double(score.fraction)).toString());
			//	outFile.write("\t");
				outFile.write(term.name);
				outFile.write("\t");
				outFile.write((new Integer(term.numGenes)).toString());
				outFile.write("\n");
			}
		}
		outFile.write(GeneSet.outSep2);
		outFile.write("\n");
	}
}
