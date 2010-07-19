package edu.mit.csail.psrg.georg.DiscretizeData;

import java.io.IOException;
import java.util.ArrayList;
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
import edu.mit.csail.psrg.georg.IGMM.MacklisIGMM;

public class DiscretizeMacklis {
	public static void main(String[] args) {
		processData();
	}
	
	public static void processData() {
		String dirName = "C:\\research_data\\macklis\\";
		String intensity_fName = "macklis_intensity.txt";
		String pvals_fName = "macklis_pvals.txt";
		String discretizedMI_fName = "macklis_discretized_MI.txt";
		String discretizedLevels_fName = "macklis_discretized_levels.txt";
		String discretizedData_fName = "macklis_discretized_data.txt";
		double pvalThresh = 0.0001;
		
		MicroArrayData intensityData = new MicroArrayData();
		MicroArrayData pvalData = new MicroArrayData();
		
		try {
			intensityData.readFile(dirName + intensity_fName);
			pvalData.readFile(dirName + pvals_fName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		transformRelativeAbundance(intensityData.values,pvalData.values,pvalThresh);
		combineReplicates(intensityData);
		
		int i = 0;
		int j = 0;
		
		int[] filterList = new int[intensityData.numRows];
		int numFiltered = 0;
		double maxExpressedPercent = 0.75;
		int maxExpressed = (int) Math.floor(maxExpressedPercent*((double) intensityData.numCols));
		int numExpressed = 0;
		int numUp = 0;
		
		// filter genes that are expressed everywhere
		for (i=0;i<intensityData.numRows;i++) {
			numExpressed = 0;
			numUp = 0;
			for (j=0;j<intensityData.numCols;j++) {
			//	if (pvalData.values[i][j] < pvalThresh) {
				if (!Double.isNaN(intensityData.values[i][j])) {
					numExpressed++;
				}
				if (!Double.isNaN(intensityData.values[i][j])) {
					numUp++;
				}
			}
			
			if (numExpressed >= maxExpressed || numUp == 0) {
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
		
		DiscretizeExpression discretizer = new DiscretizeExpression(intensityData.values,50);
		intensityData.values = discretizer.expression;
		
		try {
			discretizer.discretizeDownTree(dirName+discretizedMI_fName);
			discretizer.mergeDownLevels(5,dirName+discretizedLevels_fName);
			discretizer.transposeDExpression();
			intensityData.setDiscrete();
			intensityData.values = null;
			intensityData.dvalues = discretizer.dExpression;
			intensityData.writeFile(dirName+discretizedData_fName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
	}
	
	public static void combineReplicates(MicroArrayData intensityValues) {
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
		String[] geneNames = new String[gmap.size()];
		String geneName;
		Iterator<String> iter = gmap.keySet().iterator();
		
		double[] mu = new double[intensityValues.values[0].length];
		double[] norm = new double[intensityValues.values[0].length];
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
			}
			for (idx=0;idx<genes.size();idx++) {
				i = genes.get(idx);
				for (j=0;j<intensityValues.numCols;j++) {
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
	}
	
	public static void transformRelativeAbundance(double[][] expression,double[][] pvals,double pvalThresh) {
		double norm = 0.0;
		int i = 0;
		int j = 0;
		int numPresent = 0;
		
		for (i=0;i<expression.length;i++) {
			norm = 0.0;
			numPresent = 0;
			for (j=0;j<expression[i].length;j++) {
				if (expression[i][j] > 0.0) {
					norm = norm + expression[i][j];
				}
				numPresent++;
			}
			for (j=0;j<expression[0].length;j++) {
				expression[i][j] = expression[i][j]/norm;
				if (pvals[i][j] > pvalThresh || expression[i][j] <= 1.0/((double)numPresent)) {
					expression[i][j] = 0.0/0.0;
				}
			}
		}
	}
}
