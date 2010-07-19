package com.newmacondo.georg.DynamicDUO;

import java.io.IOException;

public class DynamicDUO {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	//	SynthesizeData ms = new SynthesizeData();
	//	DP myDP = ms.init();
		DP myDP = new DP();
		String fname = "c:\\DynamicDUO\\test.dat";
		String clusterName = "c:\\DynamicDUO\\clusters.dat";
		DataReader dr = new DataReader(myDP);
		try {
		//	dr.writeFile(fname);
			dr.readFile(fname);
		} catch(IOException e) {
			System.out.print(e);
		}
		myDP.init();
		myDP.iterate(10000);
		try {
			myDP.writeClusters(clusterName);
		} catch(IOException e) {
			System.out.print(e);
		}
	}

}
