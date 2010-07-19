package edu.mit.csail.psrg.georg.HDP2;

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

import edu.mit.csail.psrg.georg.GO.DAGException;
import edu.mit.csail.psrg.georg.GO.GoOBOReader;
import edu.mit.csail.psrg.georg.GO.MouseHumanGoAssocReader;
import edu.mit.csail.psrg.georg.GO.OntClusters;
import edu.mit.csail.psrg.georg.GO.OntDAG;
import edu.mit.csail.psrg.georg.GO.OntGeneAssociations;
import edu.mit.csail.psrg.georg.GO.OntTerm;

// this class will take a series of .persist files,
// extract EPs, and then merge EPs from subsequent samples
// to create a coherent set of REPs
public class REPOverlapBuilder {
	public ArrayList<REPOverlap> overlaps = new ArrayList<REPOverlap>();
	public ArrayList<REPOverlapTrack> overlapTrack = new ArrayList<REPOverlapTrack>();
	public int minGenesInREP = 5;
	public double pvalThreshold = 0.05;
	public double minSimilarity = 0.50;
	public int minGenes = 5;
	public int maxGenes = 200;
	// minimum count of genes in topic
	public int minOccurInTopic = 5;
	public double FDR = 0.05;
	// minimum percentage occurance for a gene in a topic
	public double minGeneOccurPercent = 0.05;
//	public double minTopicOccurPercent = 0.95;
	public double minTopicOccurPercent = 0.51;
	public double minPercentSignificance = 0.05;

	public int namespace = OntTerm.BiologicalProcess;
	public HashSet<String> useGenes = new HashSet<String>();
	
	public String refDirPath = "C:\\research_data\\mouse_human\\b47mm7\\";
	public String humanAssocFile = refDirPath + "Hs.genelist.go.plus";
	
	private class geneHolder implements Comparable {
		String name = "";
		double rank = 0.0;
		double occur = 0.0;
		
		public geneHolder(String n,double r,double o) {
			name = n;
			rank = r;
			occur = o;
		}
		
		public int compareTo(Object o) {
			geneHolder holder = (geneHolder) o;
			if (rank == holder.rank)
				return 0;
			if (rank < holder.rank)
				return 1;
			return -1;
		}
	}
	
	private class geneList {
		ArrayList<geneHolder> genes = new ArrayList<geneHolder>();
		
		void addGene(String name,double rank,double occur) {
			genes.add(new geneHolder(name,rank,occur));
		}
		
		int size() {
			return genes.size();
		}
		
		double[] sort(ArrayList<String> geneNames) {
			geneHolder[] garray = new geneHolder[genes.size()];
			garray = genes.toArray(garray);
			List<geneHolder> glist = Arrays.asList(garray);
			Collections.sort(glist);
			
			double[] ranks = new double[genes.size()];
			int i = 0;
			Iterator<geneHolder> iter = glist.iterator();
			geneHolder gene = null;
			while(iter.hasNext()) {
				gene = iter.next();
				ranks[i] = gene.rank;
				geneNames.add(gene.name);
				i++;
			}
			
			return ranks;
		}
	}
	
	public void processSamples(String baseFile,String outfName,int sampleStart,int sampleInterval,int numSamples) {
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
	/*	for (i=0;i<40;i++) {
			tissueMap.put("sample_" + ((new Integer(i+1))).toString(),i);
		} */
		
		for (i=0;i<numSamples;i++) {
			sampleID = sampleStart + i*sampleInterval;
			fName = baseFile + ((Integer) sampleID).toString();
			fName = fName + ".persist";
			myHDP = HDP.restoreFromFile(fName);
			System.out.println("Processing sample " + ((Integer) sampleID).toString());
			addSample(myHDP,tissueMap);
		}
		
	/*	try {
			outputOverlapREPs(outfName,myHDP,numSamples,tissueMap);
		} catch(IOException e) {
			System.out.println(e);
		} */
		
		System.out.println("Done");
	}
	
	public void addSample(HDP myHDP,LinkedHashMap<String,Integer> tissueMap) {
		int i = 0;
  		int j = 0;
  		double p = 0;
  		DirichletProcess dp = null;
  		REPOverlap[] tempREPs = new REPOverlap[myHDP.numTopics];
  	//	REPOverlapTrack[] tempREPTrack = new REPOverlapTrack[myHDP.numTopics];
  		
  		for (i=0;i<myHDP.numTopics;i++) {
  			tempREPs[i] = new REPOverlap(myHDP.totalGenes);
  			tempREPs[i].addGeneCounts(myHDP.topics.get(i).posCounts);
  	//		tempREPTrack[i] = new REPOverlapTrack(myHDP.geneNames.length,tissueMap.size());
  		}
  		
  		for (j=0;j<myHDP.DP.length;j++) {
  			dp = myHDP.DP[j];
  			if (dp.state != DirichletProcess.HELDOUT & dp.getClass() == TissueDP.class) {
				((TissueDP) dp).allocateGenesToTopicMap();
				for (i=0;i<myHDP.numTopics;i++) {
					tempREPs[i].addTissue(((TissueDP) dp).upDownTopicChoice[i],myHDP.topics.get(i).topicMap,myHDP.topics.get(i).topicArray,(TissueDP) dp,minGenesInREP,pvalThreshold);
		//			tempREPTrack[i].addTissue(myHDP,((TissueDP) dp),tissueMap,myHDP.topics.get(i).topicArray);
				}
			}
		}
  		
  		for (i=0;i<myHDP.numTopics;i++) {
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
  			//	overlapTrack.get(bestSimilarityIDX[i]).add(tempREPTrack[i]);
  			} else {
  				overlaps.add(tempREPs[i]);
  			//	overlapTrack.add(tempREPTrack[i]);
  			}
  		}
	}
	
	public void outputOverlapREPs(String fOutName,HDP myHDP,int numSamples,LinkedHashMap<String,Integer> tissueMap) throws IOException {
		ArrayList<String> genes = new ArrayList<String>();
		int i = 0;
		int j = 0;
		int k = 0;
		GoOBOReader goReader = new GoOBOReader();
		MouseHumanGoAssocReader humanAssocReader = new MouseHumanGoAssocReader();
		humanAssocReader.setFile(humanAssocFile);
		OntDAG humanDAG = null;
		OntGeneAssociations humanAssoc = null;
		
		String[] humanGeneNames = myHDP.geneNames;
		
		HashSet<String> humanUseGenes = new HashSet<String>();
		for (i=0;i<humanGeneNames.length;i++) {
			humanUseGenes.add(humanGeneNames[i]);
		}
		
		try {
			humanDAG = goReader.readFile();
			System.out.println("Loaded human GO DAG");
			humanAssoc = humanAssocReader.readFile(humanDAG,humanUseGenes);
			System.out.println("Loaded human GO associations");
		} catch(IOException e) {
			System.out.println(e);
		} catch(DAGException e) {
			System.out.println(e);
		}
		
		OntClusters clusters = new OntClusters(humanAssoc,minGenes,maxGenes,FDR,namespace);
		
		geneList glist = null;
		double[] ranks = null;
		double norm = 0.0;
		REPOverlap rep = null;
		double v1 = 0.0;
		double v2 = 0.0;
		Vector<OverlapTissue> tissues = new Vector<OverlapTissue>();
		Iterator<OverlapTissue> tissueIter = null;
		OverlapTissue tissue = null;
		
		FileWriter file = new FileWriter(fOutName+".txt");
		ArrayList<REPOverlap> useREPs = new ArrayList<REPOverlap>();
		int gene = 0;
		
		int numGoodREPs = 0;
		for (i=0;i<overlaps.size();i++) {
			rep = overlaps.get(i);
			v1 = ((double) rep.numOccur)/((double) numSamples);
			if (v1 >= minTopicOccurPercent) {
				numGoodREPs++;
			}
		}
		
		double[][] REPMatrix = new double[humanGeneNames.length][numGoodREPs];
		int REPIDX = 0;
		
		FileWriter REPMatrixFiles = null;
		
		for (i=0;i<overlaps.size();i++) {
			rep = overlaps.get(i);
			v1 = ((double) rep.numOccur)/((double) numSamples);
			if (v1 >= minTopicOccurPercent) {
				tissues.clear();
				tissueIter = rep.tissues.values().iterator();
				while (tissueIter.hasNext()) {
					tissue = tissueIter.next();
					v2 = ((double) tissue.numSignificant)/((double) rep.numOccur);
					if (v2 >= minPercentSignificance) {
						tissues.add(tissue);
					}
				}
				if (tissues.size() >= 1) {
					genes.clear();
					glist = new geneList();
					for (j=0;j<rep.geneOccurNormalized.length;j++) {
						v1 = rep.geneOccurNormalized[j]/((double) rep.numOccur);
						v2 = rep.geneIntensityNormalized[j]/((double) rep.numOccur);
					//	if (v1 >= 0.50) {
						if (v1 >= minGeneOccurPercent) {
							gene = myHDP.reverseGeneMap[j];
							REPMatrix[gene][REPIDX] = v2;
							glist.addGene(humanGeneNames[gene],v2,v1);
						//	glist.addGene(humanGeneNames[gene],v1,v2);
						}
					}
					
					if (glist.size() >= minGenesInREP) {
						useREPs.add(rep);
						ranks = glist.sort(genes);
						clusters.addCluster(genes,ranks,null,1.0);
						
						REPMatrixFiles = new FileWriter(fOutName + "_REPMatrix_" + (new Integer(REPIDX+1)).toString()+".txt");
						REPMatrixFiles.write("genes");
						Iterator<String> sIter = tissueMap.keySet().iterator();
						while(sIter.hasNext()) {
							REPMatrixFiles.write("\t");
							REPMatrixFiles.write(sIter.next());
						}
						REPMatrixFiles.write("\n");
						for (k=0;k<myHDP.geneNames.length;k++) {
							REPMatrixFiles.write(myHDP.geneNames[k]);
							for (j=0;j<tissueMap.size();j++) {
								REPMatrixFiles.write("\t");
								v1 = (double) overlapTrack.get(i).numOccur;
								v2 = (double) overlapTrack.get(i).rMatrix[k][j];
								if (v1 > 0) {
									v1 = v2/v1;
								}
								REPMatrixFiles.write((new Double(v1)).toString());
							}
							REPMatrixFiles.write("\n");
						}
					}
					REPMatrixFiles.close();
				}
				REPIDX++;
			}
		}
		
		clusters.clusterSignif();
		
		for (i=0;i<useREPs.size();i++) {
			rep = useREPs.get(i);
			v1 = ((double) rep.numOccur)/((double) numSamples);
			
			tissues.clear();
			tissueIter = rep.tissues.values().iterator();
			while (tissueIter.hasNext()) {
				tissue = tissueIter.next();
				v2 = ((double) tissue.numSignificant)/((double) rep.numOccur);
				if (v2 >= minPercentSignificance) {
					tissues.add(tissue);
				}
			}
			Collections.sort(tissues);
			file.write("REP#" + ((new Integer(i+1))).toString());
			file.write("\t" + String.format("%.3f",v1));
			file.write("\n");
			file.write("tissues");
			for (j=0;j<tissues.size();j++) {
				file.write("\t");
				file.write(tissues.get(j).tissueName);
			}
			file.write("\n");
			file.write("UD");
			for (j=0;j<tissues.size();j++) {
				file.write("\t");
				v1 = (double) tissues.get(j).numUp;
				v1 = v1 / ((double) rep.numOccur);
				if (v1 < 0.50) {
					v1 = v1 - 1.0;
				}
				file.write(String.format("%.3f",v1));
			}
			file.write("\n");
			file.write("topicUse");
			for (j=0;j<tissues.size();j++) {
				file.write("\t");
				v1 = tissues.get(j).topicUse;
				v1 = v1 / ((double) rep.numOccur);
				file.write(String.format("%.3f",v1));
			}
			file.write("\n");
			file.write("tissueUse");
			for (j=0;j<tissues.size();j++) {
				file.write("\t");
				v1 = tissues.get(j).tissueUse;
				v1 = v1 / ((double) rep.numOccur);
				file.write(String.format("%.3f",v1));
			}
			file.write("\n");
			file.write("signif");
			for (j=0;j<tissues.size();j++) {
				file.write("\t");
				v1 = (double) tissues.get(j).numSignificant;
				v1 = v1 / ((double) rep.numOccur);
				file.write(String.format("%.3f",v1));
			}
			file.write("\n");
			clusters.clusters.get(i).outputScoresLong(file,true);
		}	
		file.close();
		
		file = new FileWriter(fOutName+"_matrix.txt");
		for (i=0;i<REPMatrix.length;i++) {
			for (j=0;j<REPMatrix[0].length;j++) {
				file.write((new Double(REPMatrix[i][j])).toString());
				if (j<REPMatrix[0].length-1) {
					file.write("\t");
				}
			}
			file.write("\n");
		}
		file.close();
	}
}
