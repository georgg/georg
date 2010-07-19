package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Vector;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.StatUtil.StatUtil;
import edu.mit.csail.psrg.georg.StatUtil.StepDownFDR;
import edu.mit.csail.psrg.georg.StatUtil.VectorUtil;

public class MotifSignificance {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		outputFiles();
	}
	
	public static void outputFiles() {
		String dirPath = "C:\\research_data\\mouse_human\\b47mm7\\geo\\";
		String topGenesFileName = dirPath + "top_genes_50.txt";
		String motifFileName = dirPath + "\\motifs\\50_6\\combined_cons_merged_50_6.txt";
		String outputFileName = dirPath + "\\motifs\\50_6\\motif_signif.txt";
		
		double FDR = 0.05;
		
		ArrayList<String> topicNums = new ArrayList<String>();
		ArrayList<ArrayList<String> > geneSets = null;
		MicroArrayData motifData = new MicroArrayData();
		motifData.setDiscrete();
		
		try {
			geneSets = readTopGenes(topGenesFileName,topicNums);
			motifData.readFile(motifFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		HashMap<String,Integer> geneMap = new HashMap<String,Integer>();
		
		int i = 0;
		int j = 0;
		
		for (i=0;i<motifData.numRows;i++) {
			geneMap.put(motifData.geneNames[i],i);
		}
		
		int[] motifCounts = new int[motifData.numCols];
		
		int N = 0;
		int s = 0;
		int x = 0;
		int n = 0;
		
		for (j=0;j<motifData.numCols;j++) {
			for (i=0;i<motifData.numRows;i++) {
				if (motifData.dvalues[i][j] > 0) {
					motifCounts[j]++;
				}
			}
		}
		
		ArrayList<Vector<Double>> pvalSets = new ArrayList<Vector<Double>>();
		Vector<Double> pvals = null;
		
		int cc = 0;
		int gg = 0;
		ArrayList<String> genes = null;
		String gene = null;
		
		HashSet<String> tempHash = new HashSet<String>();
		for (cc=0;cc<geneSets.size();cc++) {
			pvalSets.add(new Vector<Double>());
			genes = geneSets.get(cc);
			for (gg=0;gg<genes.size();gg++) {
				tempHash.add(genes.get(gg));
			}
		}
		N = tempHash.size();
		double pval = 0.0;
		
		for (cc=0;cc<geneSets.size();cc++) {
			genes = geneSets.get(cc);
			n = genes.size();
			pvals = pvalSets.get(cc);
			for (j=0;j<motifData.numCols;j++) {
				s = motifCounts[j];
				x = 0;
				for (gg=0;gg<genes.size();gg++) {
					gene = genes.get(gg);
					if (geneMap.containsKey(gene)) {
						i = geneMap.get(gene);
						if (motifData.dvalues[i][j] > 0) {
							x++;
						}
					}
				}
				if (x >= 3) {
					try {
						pval = 1.0 - StatUtil.hyperGeometricCDF(x-1,N,s,n);
					} catch(ArithmeticException e) {
						pval = 1.0;
					}
				} else {
					pval = 1.0;
				}
				pvals.add(pval);
			}
		}
		
		int[] order = null;
		double[] sortPvals = null;
		boolean[] accept = null;
		
		String currentTopic = "";
		
		try {
			FileWriter file = new FileWriter(outputFileName);
			for (cc=0;cc<pvalSets.size();cc++) {
				pvals = pvalSets.get(cc);
				currentTopic = topicNums.get(cc);
				file.write(currentTopic);
				file.write("\n");
				order = VectorUtil.sortOrder(pvals);
				Collections.sort(pvals);
				sortPvals = new double[pvals.size()];
				for (i=0;i<pvals.size();i++) {
					sortPvals[i] = pvals.get(i);
				}
				accept = StepDownFDR.correctPvals(sortPvals,FDR);
				for (j=0;j<pvals.size();j++) {
					if (accept[j]) {
						file.write(motifData.experimentNames[order[j]]);
						file.write("\t");
						file.write((new Double(sortPvals[j])).toString());
						file.write("\n");
					}
				}
				file.write("\n");
			}
			
			file.close();
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static ArrayList<ArrayList<String> > readTopGenes(String topGenesFileName,ArrayList<String> topicNums) throws IOException {
		ArrayList<ArrayList<String> > modules = new ArrayList<ArrayList<String> >();
		ArrayList<String> genes = null;
		
		int i = 0;
		String[] lineSplit = null;
		BufferedReader is = new BufferedReader(new FileReader(topGenesFileName));
		String line = "";		
		line = is.readLine();
		while(line != null) {
			lineSplit = line.split("\t");
			genes = new ArrayList<String>();
			modules.add(genes);
			topicNums.add(lineSplit[0]);
			if (lineSplit.length > 1) {
				for (i=1;i<lineSplit.length;i++) {
					genes.add(lineSplit[i]);
				}
			}
			line = is.readLine();
		}
		is.close();
		
		return modules;
	}

}
