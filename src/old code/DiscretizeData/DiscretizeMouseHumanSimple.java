package edu.mit.csail.psrg.georg.DiscretizeData;

import java.io.IOException;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

public class DiscretizeMouseHumanSimple {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		discretizeData();
	//	doPolish();
	}

	public static void discretizeData() {
		String dirName = "C:\\research_data\\mouse_human\\b47mm7\\";
		String humanEName = dirName + "GNF1H.gcRMA.geo.txt";
		String mouseEName = dirName + "GNF1M.gcRMA.geo.txt";
		
		String outEName = dirName + "\\geo\\GNF1_geo_merged_log2.txt";
		boolean doRows = true;
		
		double threshold = 1.0;
		
		MicroArrayData humanData = new MicroArrayData();
		MicroArrayData mouseData = new MicroArrayData();
		
		try {
			humanData.readFileBlindContinuous(humanEName);
			mouseData.readFileBlindContinuous(mouseEName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int[][] humanDiscrete = new int[humanData.numRows][humanData.numCols];
		int[][] mouseDiscrete = new int[mouseData.numRows][mouseData.numCols];
		int[] humanCode = new int[humanData.numCols];
		int[] mouseCode = new int[mouseData.numCols];
		
		int i = 0;
		int j = 0;
		
		for (j=0;j<humanCode.length;j++) {
			humanCode[j] = 1;
			humanData.experimentNames[j] = "h_" + humanData.experimentNames[j];
		}
		for (j=0;j<mouseCode.length;j++) {
			mouseCode[j] = 2;
			mouseData.experimentNames[j] = "m_" + mouseData.experimentNames[j];
		}
		
		for (i=0;i<humanData.numRows;i++) {
			humanData.geneNames[i] = humanData.geneNames[i] + "|" + mouseData.geneNames[i];
		}
		
		humanData.experimentCode = humanCode;
		mouseData.experimentCode = mouseCode;
		
		humanData.log2Data();
		mouseData.log2Data();
		
	/*	if (doRows) {
			humanData.normRows();
			mouseData.normRows();
		} else {
			humanData.normCols();
			mouseData.normCols();
		} */
		
	/*	for (i=0;i<humanData.numRows;i++) {
			for (j=0;j<humanData.numCols;j++) {
				if (humanData.values[i][j] > threshold) {
					humanDiscrete[i][j] = (int) Math.round(humanData.values[i][j]);
				}
			}
		}
		
		for (i=0;i<mouseData.numRows;i++) {
			for (j=0;j<mouseData.numCols;j++) {
				if (mouseData.values[i][j] > threshold) {
					mouseDiscrete[i][j] = (int) Math.floor(mouseData.values[i][j]);
				}
			}
		}
		
		humanData.values = null;
		humanData.dvalues = humanDiscrete;
		humanData.setDiscrete();
		
		mouseData.values = null;
		mouseData.dvalues = mouseDiscrete;
		mouseData.setDiscrete(); */
		
		MicroArrayData merged = MicroArrayData.mergeCols(humanData,mouseData);
		
		try {
			merged.writeFile(outEName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void doPolish() {
		String dirName = "C:\\research_data\\mouse_human\\b47mm7\\";
		String humanEName = dirName + "GNF1H.gcRMA.geo.txt";
		String mouseEName = dirName + "GNF1M.gcRMA.geo.txt";
		
		String outEName = dirName + "\\geo\\GNF1_polished_merged.txt";
		
		double threshold = 1.0;
		
		MicroArrayData humanData = new MicroArrayData();
		MicroArrayData mouseData = new MicroArrayData();
		
		try {
			humanData.readFileBlindContinuous(humanEName);
			mouseData.readFileBlindContinuous(mouseEName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int[][] humanDiscrete = new int[humanData.numRows][humanData.numCols];
		int[][] mouseDiscrete = new int[mouseData.numRows][mouseData.numCols];
		int[] humanCode = new int[humanData.numCols];
		int[] mouseCode = new int[mouseData.numCols];
		
		int i = 0;
		int j = 0;
		
		for (j=0;j<humanCode.length;j++) {
			humanCode[j] = 1;
			humanData.experimentNames[j] = "h_" + humanData.experimentNames[j];
		}
		for (j=0;j<mouseCode.length;j++) {
			mouseCode[j] = 2;
			mouseData.experimentNames[j] = "m_" + mouseData.experimentNames[j];
		}
		
		for (i=0;i<humanData.numRows;i++) {
			humanData.geneNames[i] = humanData.geneNames[i] + "|" + mouseData.geneNames[i];
		}
		
		humanData.experimentCode = humanCode;
		mouseData.experimentCode = mouseCode;
		
		humanData.log2Data();
		mouseData.log2Data();
		
	//	humanData.normCols();
	//	mouseData.normCols();
		
		MicroArrayData mergeData = MicroArrayData.mergeCols(humanData,mouseData);
		
		polishValues(mergeData);
		
		try {
			mergeData.writeFile(outEName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
	}
	
	public static void polishValues(MicroArrayData data) {
		int i = 0;
		int j = 0;
		double res_old = 100.0;
		double res_new = 50.0;
		int max_Iter = 1000;
		double tol = 10e-8;
		double v = 0.0;
		int iter = 0;
		
		double[][] values2 = data.copyContinousData();
		
		boolean flipFlop = false;
		
		while(Math.abs(res_old-res_new) > tol & iter < max_Iter) {
			if (flipFlop) {
				data.normRows();
				flipFlop = false;
			} else {
				data.normCols();
				flipFlop = true;
			}
			res_old = res_new;
			res_new = 0.0;
			for (i=0;i<data.numRows;i++) {
				for (j=0;j<data.numCols;j++) {
					v = Math.abs(data.values[i][j] - values2[i][j]);
					if (v > res_new) {
						res_new = v;
					}
				}
			}
			iter++;
			System.out.println(iter);
		}
	}
	
}
