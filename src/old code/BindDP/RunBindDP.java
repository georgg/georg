package edu.mit.csail.psrg.georg.BindDP;

import java.io.IOException;
import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

public class RunBindDP {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		MicroArrayData myData = new MicroArrayData();
		try {
		//	myData.readFile("C:\\CPPDataFiles\\pvals_test.txt");
			myData.readFile("C:\\CPPDataFiles\\pvals_filtered.txt");
			BindDP myDP = new BindDP(myData.values,0.005,50);
			myDP.iterate(10000,5000,1);
			myDP.expandMeanInferredBinding();
			myData.values = myDP.meanInferredBinding;
		//	myData.writeFile("C:\\CPPDataFiles\\test_output.txt");
			myData.writeFile("C:\\CPPDataFiles\\test_output_full.txt");
			myDP.outputClusters();
		} catch(IOException e) {
			System.out.println(e);
		}

	}

}
