package edu.mit.csail.psrg.georg.GeneProgram2;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Vector;

import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;

// this class will take a series of .persist files,
// extract EPs, and then merge EPs from subsequent samples
// to create a coherent set of REPs
public class REPMatrixBuilder {
	public ArrayList<REPOverlap> overlaps = new ArrayList<REPOverlap>();
	
	public int minGenesInREP = 5;
	// the percentage of the REP's genes that must be used by the tissue
	// in a sample to be deemed significant
	public double tissueUseGenesThreshold = 0.0;
	
	public double minSimilarity = 0.50;
	// the minimum similarity to merge tissues into the same consensus group
	public double minGroupMerge = 0.90;
	
	public void processSamples(String baseFile,String outfName,String modifierFileName, int sampleStart,int sampleInterval,int numSamples) {
		int i = 0;
		String fName = null;
		int sampleID = 0;
		HDP myHDP = null;
		
		fName = baseFile + ((Integer) sampleStart).toString();
		fName = fName + ".persist";
		myHDP = HDP.restoreFromFile(fName);
		LinkedHashMap<String,Integer> tissueMap = new LinkedHashMap<String,Integer>();
		int tissueNum = 0;
		for (i=0;i<myHDP.DP.length;i++) {
			if (myHDP.DP[i].getClass() == TissueDP.class) {
				tissueMap.put(myHDP.DP[i].label,tissueNum);
				tissueNum++;
			}
		} 
		
		for (i=0;i<numSamples;i++) {
			sampleID = sampleStart + i*sampleInterval;
			fName = baseFile + ((Integer) sampleID).toString();
			fName = fName + ".persist";
			
			myHDP = HDP.restoreFromFile(fName);
			System.out.println("Processing sample " + ((Integer) sampleID).toString());
			addSample(myHDP,tissueMap);
		}
		
		try {
			outputREPMatrices(outfName,myHDP,numSamples,tissueMap,modifierFileName);
			// use last HDP to output number of expression programs for samples
			myHDP.outputNumExpressionProgramSamples(outfName+"_num_EPs.txt");
			// output consensus tissue groups if groups have been collected
			if (myHDP.DPGroups.size() > 0) {
			//	if (myHDP.DPGroups.get(0).numGoodIter > 0) {
					reconstructGroups(myHDP,outfName+"_groups.txt");
			//	}
			}
		} catch(IOException e) {
			System.out.println(e);
		}
		
		System.out.println("Done");
	}
	
	public void addSample(HDP myHDP,LinkedHashMap<String,Integer> tissueMap) {
		int i = 0;
  		int j = 0;
  		int k = 0;
  		double p = 0;
  		DirichletProcess dp = null;
  		REPOverlap[] tempREPs = new REPOverlap[myHDP.numExpressionPrograms];
  		int[] modifiers = null;
  		boolean UDChoice = true;
  		
  		if (myHDP.modifierLevels != null) {
  			modifiers = new int[myHDP.modifierLevels.length];
  		}
  		
  		for (i=0;i<myHDP.numExpressionPrograms;i++) {
  			tempREPs[i] = new REPOverlap(myHDP.totalGenes);
  			tempREPs[i].addGeneCounts(myHDP.expressionPrograms.get(i).posCounts);
  		}
  		
  		for (j=0;j<myHDP.DP.length;j++) {
  			dp = myHDP.DP[j];
  			if (dp.state != DirichletProcess.HELDOUT & dp.getClass() == TissueDP.class) {
				((TissueDP) dp).allocateGenesToExpressionProgramMap();
				for (i=0;i<myHDP.numExpressionPrograms;i++) {
					if (modifiers != null) {
						for (k=0;k<modifiers.length;k++) {
							modifiers[k] = ((TissueDP) dp).modifierExpressionProgramChoice[k][i];
						}
					}
					if (myHDP.useUpDown) {
						UDChoice = ((TissueDP) dp).upDownExpressionProgramChoice[i];
					}
					tempREPs[i].addTissue(myHDP,UDChoice,modifiers,myHDP.expressionPrograms.get(i).expressionProgramMap,myHDP.expressionPrograms.get(i).expressionProgramArray,(TissueDP) dp,minGenesInREP,tissueUseGenesThreshold);
				}
			}
		}
  		
  		for (i=0;i<myHDP.numExpressionPrograms;i++) {
  			tempREPs[i].normCounts();
  		}
  		
  		double[] bestSimilarity = new double[tempREPs.length];
  		int[] bestSimilarityIDX = new int[tempREPs.length];
  		double similarity = 0.0;
  		
  		for (i=0;i<tempREPs.length;i++) {
  			bestSimilarity[i] = 0;
  			bestSimilarityIDX[i] = -1;
  			for (j=0;j<overlaps.size();j++) {
  				similarity = overlaps.get(j).similarity(tempREPs[i]);
  				if (similarity > bestSimilarity[i] & similarity >= minSimilarity) {
  					bestSimilarity[i] = similarity;
  					bestSimilarityIDX[i] = j;
  				}
  			}
  		}
  		
  		for (i=0;i<tempREPs.length;i++) {
  			if (bestSimilarityIDX[i] >= 0) {
  				overlaps.get(bestSimilarityIDX[i]).add(tempREPs[i]);
  			} else {
  				overlaps.add(tempREPs[i]);
  			}
  		}
	}
	
	public void outputREPMatrices(String fOutName,HDP myHDP,int numSamples,LinkedHashMap<String,Integer> tissueMap,String modifierFileName) throws IOException {
		int i = 0;
		int j = 0;
		int k = 0;
		
		Collections.sort(overlaps);
		
		String[] humanGeneNames = myHDP.geneNames;
		
		double norm = 0.0;
		REPOverlap rep = null;
		double v1 = 0.0;
		double v2 = 0.0;
	//	Iterator<OverlapTissue> tissueIter = null;
		OverlapTissue tissue = null;
		
		FileWriter file = new FileWriter(fOutName+"_REPs_occur.txt");
		
		for (i=0;i<overlaps.size();i++) {
			rep = overlaps.get(i);
			v1 = ((double) rep.numOccur)/((double) numSamples);
			file.write((new Double(v1)).toString());
			if (i < overlaps.size() - 1) {
				file.write("\t");
			}
		}
		file.close();
		
		file = new FileWriter(fOutName+"_REPs_generality_all.txt");
		
		for (i=0;i<overlaps.size();i++) {
			rep = overlaps.get(i);
			v1 = ((double) rep.generality)/((double) rep.numOccur);
			file.write((new Double(v1)).toString());
			if (i < overlaps.size() - 1) {
				file.write("\t");
			}
		}
		file.close();
		
		file = new FileWriter(fOutName+"_REPs_geneMatrix.txt");
		int gene = 0;
		
		file.write("genes");
		for (i=0;i<overlaps.size();i++) {
			file.write("\t");
			file.write("r" + (new Integer(i+1)).toString());
		}
		file.write("\n");
		
		for(j=0;j<myHDP.totalGenes;j++) {
			gene = myHDP.reverseGeneMap[j];
			file.write(humanGeneNames[gene]);
			for (i=0;i<overlaps.size();i++) {
				rep = overlaps.get(i);
				v1 = rep.geneOccurNormalized[j]/((double) rep.numOccur);
				v2 = rep.geneIntensityNormalized[j]/((double) rep.numOccur);
				file.write("\t");
				file.write((new Double(v2)).toString());
			}
			file.write("\n");
		}
		file.close();
		
		file = new FileWriter(fOutName+"_REPs_tissueMatrix.txt");
		Iterator<String> tNames = tissueMap.keySet().iterator();
		file.write("REPs");
		while(tNames.hasNext()) {
			file.write("\t");
			file.write(tNames.next());
		}
		file.write("\n");
		
		ArrayList<ArrayList<String> > modifierNames = null;
		int[][] modifiers = null;
		int mi = 0;
		int mj = 0;
		if (modifierFileName != null) {
			ClusterReader modifierReader = new ClusterReader();
			modifierReader.readFile(modifierFileName);
			modifierNames = modifierReader.clusters;
		}
		
		for (i=0;i<overlaps.size();i++) {
			rep = overlaps.get(i);
			if (myHDP.useUpDown) {
				file.write("UD");
				tNames = tissueMap.keySet().iterator();
				while(tNames.hasNext()) {
					tissue = rep.tissues.get(tNames.next());
					file.write("\t");
					if (tissue != null) {
						v1 = (double) tissue.numUp;
						v1 = v1 / ((double) rep.numOccur);
						if (v1 < 0.50) {
							v1 = v1 - 1.0;
						}
					} else {
						v1 = 0.0;
					}
					file.write(String.format("%.3f",v1));
				}
				file.write("\n");
			}
			
			file.write("topicUse");
			tNames = tissueMap.keySet().iterator();
			while(tNames.hasNext()) {
				tissue = rep.tissues.get(tNames.next());
				file.write("\t");
				if (tissue != null) {
					v1 = tissue.topicUse;
					v1 = v1 / ((double) rep.numOccur);
				} else {
					v1 = 0.0;
				}
				file.write(String.format("%.3f",v1));
			}
			file.write("\n");
			file.write("tissueUse");
			tNames = tissueMap.keySet().iterator();
			while(tNames.hasNext()) {
				tissue = rep.tissues.get(tNames.next());
				file.write("\t");
				if (tissue != null) {
					v1 = tissue.tissueUse;
					v1 = v1 / ((double) rep.numOccur);
				} else {
					v1 = 0.0;
				}
				file.write(String.format("%.3f",v1));
			}
			file.write("\n");
			file.write("signif");
			tNames = tissueMap.keySet().iterator();
			while(tNames.hasNext()) {
				tissue = rep.tissues.get(tNames.next());
				file.write("\t");
				if (tissue != null) {
					v1 = (double) tissue.numSignificant;
					v1 = v1 / ((double) rep.numOccur);
				} else {
					v1 = 0.0;
				}
				file.write(String.format("%.3f",v1));
			}
			file.write("\n");
			
			if (modifierFileName != null) {
				for (mi=0;mi<modifierNames.size();mi++) {
					for (mj=0;mj<modifierNames.get(mi).size();mj++) {
						file.write(modifierNames.get(mi).get(mj));
						tNames = tissueMap.keySet().iterator();
						while(tNames.hasNext()) {
							tissue = rep.tissues.get(tNames.next());
							if (tissue != null) {
								v1 = (double) tissue.numModifier[mi][mj];
								v1 = v1 / ((double) rep.numOccur);
							} else {
								v1 = 0.0;
							}
							file.write("\t");
							file.write(String.format("%.3f",v1));
						}
						file.write("\n");
					}
				}
			}
		}
		
		file.close();
	}
	
	public void reconstructGroups(HDP myDP,String groupFName) {
		if (myDP.DPGroups.size() == 0)
			return;
		
		int i = 0;
		int j = 0;
		for (i=0;i<myDP.DPGroups.size();i++) {
			if (myDP.DPGroups.get(i).numGoodIter > 0) {
				myDP.DPGroups.get(i).normPairProbs();
				myDP.DPGroups.get(i).buildConcensusModel(minGroupMerge);
			}
		}
		myDP.buildDPList();
		
		ParentDP pdp = null;
		TissueDP dp = null;
		TissueGroup group = myDP.DPGroups.get(0);
		
		try {
			FileWriter file = new FileWriter(groupFName);
			
			int dpID = 0;
			for (i=0;i<group.controllingProcess.numChildren();i++) {
				pdp = (ParentDP) group.controllingProcess.children.get(i);
				for (j=0;j<pdp.numChildren();j++) {
					dp = (TissueDP) pdp.children.get(j);
					file.write(dp.label);
					if (j < pdp.numChildren() - 1) {
						file.write("\t");
					}
				}
				file.write("\n");
			}
			file.close();
		} catch(IOException e) {
			System.out.println(e);
		}
		
		System.out.println("reconstructed groups");
	}
}
