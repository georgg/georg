package com.newmacondo.georg.DynamicDUO;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.StringTokenizer;

public class DataReader {
	DP myDP = null;
	
	public DataReader(DP d) {
		myDP = d;
	}
	
	public void readFile(String fileName) throws IOException {
		BufferedReader is = new BufferedReader(new FileReader(fileName));
		String line = "";
		
		line = is.readLine();
		StringTokenizer st = new StringTokenizer(line, "\t");
		
		myDP.numOTUs = new Integer(st.nextToken());
		myDP.numLocations = new Integer(st.nextToken());
		myDP.numTimepoints = new Integer(st.nextToken());
		myDP.numReplicates = new Integer(st.nextToken());
		
		myDP.OTUCounts = new int[myDP.numOTUs][myDP.numLocations][myDP.numTimepoints][myDP.numReplicates];
		myDP.OTUNames = new String[myDP.numOTUs];
		
		String OTUName;
		
		int i = 0; int j = 0; int k = 0; int l = 0;
		for (i=0;i<myDP.numOTUs;i++) {
			for (j=0;j<myDP.numLocations;j++) {
				line = is.readLine();
				st = new StringTokenizer(line,"\t");
				OTUName = st.nextToken();
				if (j == 0) {
					myDP.OTUNames[i] = OTUName;
				}
				st.nextToken();
				for (k=0;k<myDP.numTimepoints;k++) {
					for (l=0;l<myDP.numReplicates;l++) {
						myDP.OTUCounts[i][j][k][l] = new Integer(st.nextToken());
					}
				}
			}
		}
		is.close();
		System.out.print("OTUs=" + Integer.toString(myDP.numOTUs) + "\n");
		System.out.print("Locations=" + Integer.toString(myDP.numLocations) + "\n");
		System.out.print("Timepoints=" + Integer.toString(myDP.numTimepoints) + "\n");
		System.out.print("Replicates=" + Integer.toString(myDP.numReplicates) + "\n");
	}
	
	public void writeFile(String fname) throws IOException {
		int i = 0;
		int j = 0;
		int k = 0;
		int l = 0;
		FileWriter outFile = new FileWriter(fname);
		String s;
		int value = 0;
		
		s = Integer.toString(myDP.numOTUs) + "\t" + Integer.toString(myDP.numLocations) + "\t" + Integer.toString(myDP.numTimepoints) + "\t" + Integer.toString(myDP.numReplicates) + "\n";
		outFile.write(s);
		
		for(i=0;i<myDP.numOTUs;i++) {
			for(j=0;j<myDP.numLocations;j++) {
				s = myDP.OTUNames[i] + "\t" + Integer.toString(j);
				for (k=0;k<myDP.numTimepoints;k++) {
					for (l=0;l<myDP.numReplicates;l++) {
						value = myDP.OTUCounts[i][j][k][l];
						s = s + "\t" + Integer.toString(value);
					}
				}
				outFile.write(s+"\n");
			}
		}
		
		outFile.close();
	}
}
