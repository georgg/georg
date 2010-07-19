package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.StatUtil.StatUtil;

public class MotifsOrthologsModel {
	public static void main(String[] args) {
	//	divideGenes();
	//	clusterMotifs();
		outputSimilarities();
	}
	
	public static void clusterMotifs() {
		OrthoHierDirichletProcess myHierDP = constructModel(1);
		myHierDP.activateAll();
		
		int iters = 50000;
		int burnin = 25000;
		int numConcParamIter = 15;
		int groupStartIter = 29999;
		
		String mainOutName = "C:\\research_data\\mouse_human\\motifs_orthologs\\genes_as_docs";
		String snapShotName = mainOutName + "_snap";
		myHierDP.groupStartIter = groupStartIter;
		myHierDP.enableSnapShots(1000,snapShotName);
		myHierDP.iterate(iters,burnin,numConcParamIter,false);	
		System.out.println("Done");
	}
	
	public static void outputSimilarities() {
		String dirPath = "C:\\research_data\\mouse_human\\motifs_orthologs\\";
		String snapFileName = dirPath + "genes_as_docs_snap_49999.persist";
		String dpFileName = dirPath + "genes_as_docs_dps.txt";
		String similaritiesFileName = dirPath + "similarities.txt";
		OrthoHierDirichletProcess myDP = (OrthoHierDirichletProcess) OrthoHierDirichletProcess.restoreFromFile(snapFileName);
		
		MicroArrayData dpData = new MicroArrayData();
		MicroArrayData topicData = new MicroArrayData();
		
		try {
		//	myDP.outputClustersToFile(dirPath + "genes_as_docs");
		//	dpData.readFile(dpFileName);
			myDP.myOrthologStats.outputSimilarities(similaritiesFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static OrthoHierDirichletProcess constructModel(int divisionNum) {
		String combinedMotifFileName = "C:\\research_data\\mouse_human\\motifs_orthologs\\combined_2000bp_6ord_discretized.txt";
		String divisionFileName = "C:\\research_data\\mouse_human\\motifs_orthologs\\divide_genes.txt";
	
		ArrayList<ArrayList<Integer> > divisions = null;
		MicroArrayData motifData = new MicroArrayData();
		motifData.setDiscrete();
		
		try {
			motifData.readFile(combinedMotifFileName);
			divisions = readDivisions();
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		int j = 0;
		int k = 0;
		int mm = 0;
		int gg = 0;
		
		HDPConcParam[] concParams = new HDPConcParam[2];
		HDPConcParam root_conc = new HDPConcParam(1.0,0.1);
		HDPConcParam alpha_0 = new HDPConcParam(1.0,1.0);
		concParams[0] = root_conc;
		concParams[1] = alpha_0;
		
		Document doc = null;
		ArrayList<Integer> glist = new ArrayList<Integer>();
		ArrayList<DirichletProcess> DPs = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> alpha_0_DP = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> root_conc_DP = new ArrayList<DirichletProcess>();
		DirichletProcess root = null;
		
		root = new DirichletProcess(null,"root");
		root_conc_DP.add(root);
		DPs.add(root);
		root_conc.addDPs(root_conc_DP);
		
		int numTimes = 0;
		String[] humanOrthologs = new String[divisions.get(divisionNum).size()];
		String[] mouseOrthologs = new String[divisions.get(divisionNum).size()];
		
		ArrayList<Integer> geneList = new ArrayList<Integer>();
		DirichletProcess d = null;
		for (i=0;i<divisions.get(divisionNum).size();i++) {
			gg = divisions.get(divisionNum).get(i);
			geneList.add(gg);
			geneList.add(gg+motifData.numRows/2);
			humanOrthologs[i] = motifData.geneNames[gg];
			mouseOrthologs[i] = motifData.geneNames[gg+motifData.numRows/2];
		}
		
		for (i=0;i<geneList.size();i++) {
			gg = geneList.get(i);
			glist.clear();
			for (j=0;j<motifData.numCols;j++) {
				numTimes = motifData.dvalues[gg][j];
				if (numTimes > 0) {
					for (k=0;k<numTimes;k++) {
						glist.add(j);
					}
				}
			}
			doc = new Document(motifData.geneNames[gg],glist);
			d = new DirichletProcess(root,doc.annotation);
			d.addDocument(doc);
			DPs.add(d);
			alpha_0_DP.add(d);
		}
		alpha_0.addDPs(alpha_0_DP);

		int init_numClusters = 1;
		OrthoHierDirichletProcess myHierDP = new OrthoHierDirichletProcess(init_numClusters,DPs,concParams,motifData.experimentNames);
		myHierDP.packGenes();
		myHierDP.myOrthologStats.init(myHierDP,humanOrthologs,mouseOrthologs);
		myHierDP.setCollectStrongModules(false);
		
		return myHierDP;
	}
	
	public static ArrayList<ArrayList<Integer> > readDivisions() {
		String divisionFileName = "C:\\research_data\\mouse_human\\motifs_orthologs\\divide_genes.txt";	
		ArrayList<ArrayList<Integer> > divisions = new ArrayList<ArrayList<Integer> >();
		ArrayList<Integer> geneSet = null;
		int totalGenes = 0;
		
		int i = 0;
		String[] lineSplit = null;
		try {
			BufferedReader is = new BufferedReader(new FileReader(divisionFileName));
			String line = "";
			line = is.readLine();
			while(line != null) {
				lineSplit = line.split("\t");
				geneSet = new ArrayList<Integer>();
				divisions.add(geneSet);
				for (i=0;i<lineSplit.length;i++) {
					geneSet.add(new Integer(lineSplit[i]));
					totalGenes++;
				}
				line = is.readLine();
			}
			is.close();
		} catch(IOException e) {
			System.out.println(e);
		}
		
		System.out.println(totalGenes);
		return divisions;
	}
	
	public static void divideGenes() {
		int numDivisions = 6;
		String combinedMotifFileName = "C:\\research_data\\mouse_human\\motifs_orthologs\\combined_2000bp_6ord_discretized.txt";
		String divisionFileName = "C:\\research_data\\mouse_human\\motifs_orthologs\\divide_genes.txt";	
		
		MicroArrayData motifData = new MicroArrayData();
		motifData.setDiscrete();
		
		try {
			motifData.readFile(combinedMotifFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int numGenes = motifData.numRows/2;
		
		int i = 0;
		int j = 0;
		int k = 0;
		ArrayList<ArrayList<Integer> > divisions = new ArrayList<ArrayList<Integer> >();
		ArrayList<Integer> geneSet = new ArrayList<Integer>();
		int[] r = StatUtil.randPerm(numGenes);
		
		divisions.add(geneSet);
		j = 0;
		k = Math.round(((float) numGenes)/((float) numDivisions));
		for (i=0;i<numGenes;i++) {
			if (j == k) {
				j = 1;
				geneSet = new ArrayList<Integer>();
				divisions.add(geneSet);
			} else {
				j++;
			}
			geneSet.add(r[i]);
		}
		
		try {
			FileWriter file = new FileWriter(divisionFileName);
			for (i=0;i<divisions.size();i++) {
				for (j=0;j<divisions.get(i).size();j++) {
					file.write((new Integer(divisions.get(i).get(j))).toString());
					if (j<divisions.get(i).size()-1) {
						file.write("\t");
					} else {
						file.write("\n");
					}
				}
			}
			file.close();
		} catch(IOException e) {
			System.out.println(e);
		}
	}
}
