package edu.mit.csail.psrg.georg.BiclusterAnalysis;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedHashSet;

import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;
import edu.mit.csail.psrg.georg.HDPTopicModel.GeneList;

public class ConvertFiles {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
	
		convertSamba("c:\\research_data\\Shyamsundar\\MHS_samba");
	}
	
	public static void convertSamba(String baseFile) {
		ArrayList<ArrayList<String>> biclusterTissues = new ArrayList<ArrayList<String>>();
		ArrayList<ArrayList<String>> biclusterGenes = new ArrayList<ArrayList<String>>();
		String line = "";
		String[] split = null;
		int dataCode = 0;
		int biclusterNum = 0;
		String itemName = null;
		int numBiclusters = -1;
		ArrayList<String> myGenes = null;
		ArrayList<String> myTissues = null;
		
		try {
			BufferedReader readFile = new BufferedReader(new FileReader(baseFile + ".txt"));
			line = readFile.readLine();
			while(readFile.ready()) {
				split = line.split("\t");
				biclusterNum = new Integer(split[0]).intValue();
				dataCode = new Integer(split[1]).intValue();
				itemName = split[2];
				if (biclusterNum > numBiclusters) {
					myGenes = new ArrayList<String>();
					biclusterGenes.add(myGenes);
					myTissues = new ArrayList<String>();
					biclusterTissues.add(myTissues);
					numBiclusters++;
				} else {
					myGenes = biclusterGenes.get(biclusterNum);
					myTissues = biclusterTissues.get(biclusterNum);
				}
				
				if (dataCode == 0) {
					myTissues.add(itemName);
				}
				
				if (dataCode == 1) {
					myGenes.add(itemName);
				}
				line = readFile.readLine();
			}
			readFile.close();
		} catch(IOException e) {
			System.out.println(e);
		}
		
		String geneOutFile = baseFile + "_bicluster_genes.txt";
		String tissueOutFile = baseFile + "_bicluster_tissues.txt";
		int i = 0;
		int j = 0;
		
		try {
			FileWriter geneFile = new FileWriter(geneOutFile);
			FileWriter tissueFile = new FileWriter(tissueOutFile);
			
			for (i=0;i<biclusterTissues.size();i++) {
				myTissues = biclusterTissues.get(i);
				for (j=0;j<myTissues.size();j++) {
					tissueFile.write(myTissues.get(j));
					if (j < myTissues.size() - 1) {
						tissueFile.write("\t");
					}
				}
				tissueFile.write("\n");
			}
			tissueFile.close();
			
			for (i=0;i<biclusterGenes.size();i++) {
				myGenes = biclusterGenes.get(i);
				for (j=0;j<myGenes.size();j++) {
					geneFile.write(myGenes.get(j));
					if (j < myGenes.size() - 1) {
						geneFile.write("\t");
					}
				}
				geneFile.write("\n");
			}
			geneFile.close();
			
		} catch(IOException e) {
			System.out.println(e);
		}
	}

}
