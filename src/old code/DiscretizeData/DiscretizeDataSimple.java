package edu.mit.csail.psrg.georg.DiscretizeData;

import java.io.IOException;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

// takes as input log2 transformed data
// the values are column normalized and then rounded off to produce discrete values (including negative values)
// values with absolute fold changes < threshold are floored to zero
public class DiscretizeDataSimple {
	
	public static void main(String args[]) {
	//	String inFName = "C:\\research_data\\mouse_human\\b47mm7\\geo\\GNF1_geo_merged_log2.txt";
	//	String outFName = "C:\\research_data\\mouse_human\\b47mm7\\geo\\GNF1_geo_merged_UD_discrete_nonorm_floor5.txt";
		
		String inFName = "C:\\research_data\\ramaswamy\\GeneProgram\\RamaGCRMA_log2_filtered.txt";
		String outFName = "C:\\research_data\\ramaswamy\\GeneProgram\\RamaGCRMA_log2_discrete.txt";
		
	//	String inFName = "C:\\research_data\\Bild\\Lung_geo_filter_log2.txt";
	//	String outFName = "C:\\research_data\\Bild\\Lung_geo_discrete_log2.txt";
		
		double threshold = 1.0;
		
		doDiscretization(inFName,outFName,threshold);
	}
	
	public static void doDiscretization(String inFName,String outFName,double threshold) {
		MicroArrayData data = new MicroArrayData();
		
		try {
			data.readFile(inFName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		int j = 0;
		// restrict values to +/ floorValue
		double floorValue = 3.0;
		
		int[][] discrete = new int[data.numRows][data.numCols];
		double[][] orgValues = data.copyContinousData();
		
	//	data.normCols();
		
		for (i=0;i<data.numRows;i++) {
			for (j=0;j<data.numCols;j++) {
				
				if (data.values[i][j] > floorValue) {
					data.values[i][j] = floorValue;
				} else {
					if (data.values[i][j] < -floorValue) {
						data.values[i][j] = -floorValue;
					}
				}
				if (Math.abs(orgValues[i][j]) >= threshold) {
			//	if (orgValues[i][j] >= threshold) {
					discrete[i][j] = (int) Math.round(data.values[i][j]);
				}
			}
		}
		
		data.values = null;
		data.dvalues = discrete;
		data.setDiscrete();
		
		try {
			data.writeFile(outFName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
}
