package edu.mit.edu.csail.psrg.georg.AgglomerativeCluster;

import java.io.IOException;
import java.util.ArrayList;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

public class MouseHumanAgglom {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		clusterTissues();
	}
	
	public static void clusterTissues() {
		String dirPath = "C:\\research_data\\mouse_human\\b47mm7\\geo\\";
		String dataFile = dirPath + "GNF1_discretized_data_geo.txt";
		String outFile = dirPath + "agglom_tissues_geo.txt";
		
		MicroArrayData data = new MicroArrayData();
		data.setDiscrete();
		
		try {
			data.readFile(dataFile);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		int j = 0;
	/*	int k = 0;
		int kt = 0;
		ArrayList<Integer> dt = new ArrayList<Integer>();
		
		int[][] dv = new int[data.numCols][];
		for (j=0;j<data.numCols;j++) {
			dt.clear();
			for (i=0;i<data.numRows;i++) {
				kt = data.dvalues[i][j];
				if (kt > 0) {
					for (k=0;k<kt;k++) {
						dt.add(i);
					}
				}
			}
			dv[j] = new int[dt.size()];
			for (i=0;i<dt.size();i++) {
				dv[j][i] = dt.get(i);
			}
		} */ 
		int[][] dv = new int[data.numCols][data.numRows];
		for (i=0;i<data.numRows;i++) {
			for (j=0;j<data.numCols;j++) {
				dv[j][i] = data.dvalues[i][j];
			}
		}
		
		AgglomerativeCluster myClusterer = new AgglomerativeCluster(dv,data.numRows,data.experimentNames);
		myClusterer.iterate();
		
		try {
			myClusterer.outputClusters(outFile);
		} catch(IOException e) {
			System.out.println(e);
		}
	}

}
