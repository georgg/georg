package edu.mit.csail.psrg.georg.DataAccess;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

// reads a data file consisting of multiple lines
// each line is tab delimited and contains a variable number of items
public class ClusterReader {
	public ArrayList<ArrayList<String>> clusters = new ArrayList<ArrayList<String>>();
	
	public ClusterReader() {
		
	}
	
	public ClusterReader(ArrayList<ArrayList<String>> a) {
		clusters = a;
	}
	
	public void writeFile(String fname) throws IOException {
		int i = 0;
		int j = 0;
		ArrayList<String> cluster = null;
		
		FileWriter file = new FileWriter(fname);
		for (i=0;i<clusters.size();i++) {
			cluster = clusters.get(i);
			for (j=0;j<cluster.size();j++) {
				file.write(cluster.get(j));
				if (j < cluster.size() - 1) {
					file.write("\t");
				}
			}
			file.write("\n");
		}
		file.close();
	}
	
	public void readFile(String fname) throws IOException {
		String[] items = null;
		int i = 0;
		ArrayList<String> members = null;
		
		BufferedReader is = new BufferedReader(new FileReader(fname));
		String line = "";
		line = is.readLine();
		while(line != null) {
			items = line.split("\t");
			members = new ArrayList<String>();
			clusters.add(members);
			
			for (i=0;i<items.length;i++) {
				members.add(items[i]);
			}
			line = is.readLine();
		}
		is.close();
	}
}
