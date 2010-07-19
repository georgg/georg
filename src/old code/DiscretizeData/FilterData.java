package edu.mit.csail.psrg.georg.DiscretizeData;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;

import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;
import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

// Do simple filtering on data: remove genes w/ insufficient
// fold-change, etc.
// also, log2 transform the data
public class FilterData {
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String inFName = "C:\\research_data\\ramaswamy\\GeneProgram\\RamaGCRMA.geo.bygene";
		String outFName = "C:\\research_data\\ramaswamy\\GeneProgram\\RamaGCRMA_log2_filtered.txt";
	//	String selectFName = "C:\\research_data\\ramaswamy\\GeneProgram\\select_genes.txt";
		String selectFName = null;
		
	/*	String inFName = "C:\\research_data\\Bild\\Lung.geo.bygene";
		String outFName = "C:\\research_data\\Bild\\Lung_geo_filter_log2.txt";
		String selectFName = null; */
		
		// change in log space
		double foldChange = 1.0;
		double logBase = 2.0;
		filterFile(inFName,outFName,selectFName,foldChange,logBase);
	}
	
	public static void filterFile(String inFName,String outFName,String selectFName,double foldChange,double logBase) {
		int minChangeNum = 5;
		
		MicroArrayData data = new MicroArrayData();
		
		try {
			data.readFileBlindContinuous(inFName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		int j = 0;
		
		for (i=0;i<data.numRows;i++) {
			for (j=0;j<data.numCols;j++) {
				data.values[i][j] = Math.log(data.values[i][j])/Math.log(logBase);
			}
		}
		
		int numChange = 0;
		HashSet<Integer> drows = new HashSet<Integer>();
		
		for (i=0;i<data.numRows;i++) {
			numChange = 0;
			for (j=0;j<data.numCols;j++) {
				if (Math.abs(data.values[i][j]) >= foldChange) {
					numChange++;
				}
			}
			if (numChange < minChangeNum) {
				drows.add(i);
			}
		}
		
		ClusterReader selectReader = new ClusterReader();
		
		if (selectFName != null) {
			try {
				selectReader.readFile(selectFName);
			} catch(IOException e) {
				System.out.println(e);
			}
			
		
			HashSet<String> selectGenes = new HashSet<String>();
			ArrayList<String> genes = selectReader.clusters.get(0);
			for (i=0;i<genes.size();i++) {
				selectGenes.add(genes.get(i));
			}
		
			for (i=0;i<data.numRows;i++) {
				if (!selectGenes.contains(data.geneNames[i])) {
					drows.add(i);
				}
			}
		}
	
		data.deleteRows(drows);
		
		try {
			data.writeFile(outFName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}

}
