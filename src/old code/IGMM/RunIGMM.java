package edu.mit.csail.psrg.georg.IGMM;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import java.io.IOException;


public class RunIGMM {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		MicroArrayData myData = new MicroArrayData();
		String fnameIn = "C:\\CPPDataFiles\\cdc28_expression_java.txt";
		String fnameOut = "C:\\CPPDataFiles\\cdc28_expression_clusters.txt";
	//	String fnameIn = "C:\\CPPDataFiles\\cell_cycle_expression_java.txt";
	//	String fnameOut = "C:\\CPPDataFiles\\cell_cycle_expression_clusters.txt";
		
		try {
			myData.readFile(fnameIn);
		} catch(Exception e) {
			System.out.println(e);
		}
		
		GaussianDP DP = new GaussianDP(myData,1);
		
		int iter = 20000;
		int burnin = 10000;
		
		DP.iterate(iter,burnin,15);
		DP.concensusModel();
		
		try {
			DP.writeClustersToFile(fnameOut);
		} catch(Exception e) {
			System.out.println(e);
		}
	}

}
