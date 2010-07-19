package edu.mit.csail.psrg.georg.SampleUtil;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Vector;
import edu.mit.csail.psrg.georg.StatUtil.*;

public class TestSample {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	//	testAdaptiveRejectionSampler();
		testAdaptiveRejectionSampleCR();
	}
	
	public static void testAdaptiveRejectionSampleCR() {
		int numSamples = 500000;
		GaussAdaptSampleFn2 g_fn = new GaussAdaptSampleFn2(0.0,1.0);
		double[] samples = new double[numSamples];
		int i = 0;
		AdaptiveRejectionSampleCR sampler = new AdaptiveRejectionSampleCR();
		
		for(i=0;i<numSamples;i++) {
			samples[i] = sampler.ars(0.0,0.1,g_fn);
			System.out.println(i);
		}
		
		String fname = "C:\\CPPDataFiles\\gauss_samples.txt";
		
		try {
			FileWriter outFile = new FileWriter(fname);
			for(i=0;i<numSamples;i++) {
				outFile.write((new Double(samples[i])).toString());
				outFile.write("\n");
			}
			outFile.close();
		} catch(Exception e) {
			System.out.println(e);
		}
	}
	
	public static void testAdaptiveRejectionSampler() {
		GaussianAdaptSampleFn g_fn = new GaussianAdaptSampleFn(0.0,1.0);
		
		Vector<Double> x_init = new Vector<Double>();
		VectorUtil.zeros(x_init,4);
		x_init.set(0,-2.0);
		x_init.set(1,-1.0);
		x_init.set(2,1.0);
		x_init.set(3,2.0);
		
		AdaptiveRejectionSampler sampler = new AdaptiveRejectionSampler(g_fn,x_init,-1.0/0.0,1.0/0.0,100000);
		try {
			sampler.writeSamplesToFile("C:\\CPPDataFiles\\gauss_samples.txt");
		} catch(Exception e) {
			System.out.println(e);
		}
	}
}
