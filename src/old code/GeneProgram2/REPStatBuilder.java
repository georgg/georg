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

import edu.mit.csail.psrg.georg.Annotation.Annotations;
import edu.mit.csail.psrg.georg.Annotation.GeneSetCollection;
import edu.mit.csail.psrg.georg.Annotation.OntGene;
import edu.mit.csail.psrg.georg.Annotation.DAGException;
import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;
import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.StatUtil.VectorUtil;

// this class will take a series of .persist files,
// extract EPs, and then merge EPs from subsequent samples
// to create a coherent set of REPs
public class REPStatBuilder {
	public int minGenesInREP = 5;
	// the percentage of the REP's genes that must be used by the tissue
	// to be deemed significant
	public double tissueUseGenesThreshold = 0.05;
	// minimum percentage occurance for a gene in a REP
	public double minGeneOccurPercent = 0.05;
	// the minimum percentage of samples the REP must occur in
	public double minREPOccurPercent = 0.50;
	
	// parameters for GO categories
	public int minGOGenes = 5;
	public int maxGOGenes = 200;
	public double FDR = 0.05;
	public String namespace = "biological_process";
	public HashSet<String> useGenes = new HashSet<String>();
	public String GeneAssocFileName = "C:\\research_data\\mouse_human\\b47mm7\\Hs.genelist.go.plus";
	public String GOHierarchyFileName = "C:\\research_data\\go_data\\gene_ontology_plus.obo";
	public boolean outputAnnotations = false;
	
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
	
	public void outputREPStats(String baseFileName) throws IOException {
		ArrayList<String> genes = new ArrayList<String>();
		int i = 0;
		int j = 0;
		int k = 0;
		
		Annotations myAnnotations = new Annotations();
		
		String geneMatrixFileName = baseFileName + "_REPs_geneMatrix.txt";
		String tissueMatrixFileName = baseFileName + "_REPs_tissueMatrix.txt";
		String occurFileName = baseFileName + "_REPs_occur.txt";
		String generalityFileName = baseFileName + "_REPs_generality_all.txt";
		
		MicroArrayData gmd = new MicroArrayData();
		gmd.readFileBlindContinuous(geneMatrixFileName);
		String[] humanGeneNames = gmd.geneNames;
		double[][] geneMatrix = gmd.values;
		
		MicroArrayData tmd = new MicroArrayData();
		tmd.readFileBlindContinuous(tissueMatrixFileName);
		double[][] tissueMatrix = tmd.values;
		int entryInc = 0;
		i = 0;
		while(entryInc == 0 & i < tmd.numRows) {
			if (tmd.geneNames[i].equals("topicUse")) {
				entryInc = i;
			}
			i++;
		}
		
		ClusterReader omd = new ClusterReader();
		omd.readFile(occurFileName);
		double[] REPOccur = new double[omd.clusters.get(0).size()];
		for (i=0;i<REPOccur.length;i++) {
			REPOccur[i] = new Double(omd.clusters.get(0).get(i));
		}
		
		ClusterReader genmd = new ClusterReader();
		genmd.readFile(generalityFileName);
		ArrayList<Double> generalities = new ArrayList<Double>();
		
		HashSet<String> humanUseGenes = new HashSet<String>();
		for (i=0;i<humanGeneNames.length;i++) {
			humanGeneNames[i] = humanGeneNames[i].toUpperCase();
			humanUseGenes.add(humanGeneNames[i]);
		}
		
		try {
			myAnnotations.readFiles(GeneAssocFileName,GOHierarchyFileName,humanUseGenes);
			System.out.println("Loaded annotations");
		} catch(IOException e) {
			System.out.println(e);
		} catch(DAGException e) {
			System.out.println(e);
		}
		
		GeneSetCollection geneSets = new GeneSetCollection(myAnnotations,minGOGenes,maxGOGenes,FDR,namespace);
		OntGene oGene = null;
		
		geneList glist = null;
		double[] ranks = null;
		double norm = 0.0;
		double v1 = 0.0;
		double v2 = 0.0;
		Vector<Double> tissueUse = new Vector<Double>();
		Vector<Integer> tissueSelect = new Vector<Integer>();
		int[] sortIdx = null;
		
		FileWriter file = new FileWriter(baseFileName + "REPs_use.txt");
		int baseIdx = 0;
		int useREPs = 0;
		int tissueNum = 0;
		
		FileWriter file2 = null;
		double[] modAvg = null;
		double[][] tissueModAvg = null;
		
		ArrayList<Double[]> reviseTissueMatrix = new ArrayList<Double[]>();
		Double[] entriesDouble = null;
		ArrayList<Integer[]> reviseModMatrix = new ArrayList<Integer[]>();
		Integer[] entriesInteger = null;
		double maxMod = -1.0;
		int maxIdx = -1;
		String sterms;
		
		int numMod = entryInc - 3;
		if (numMod > 0) {
			file2 = new FileWriter(baseFileName + "REPs_mod_avg.txt");
			file2.write("rep");
			for (i=0;i<numMod;i++) {
				file2.write("\t");
				file2.write(tmd.geneNames[i+3]);
			}
			file2.write("\n");
			modAvg = new double[numMod];
			tissueModAvg = new double[tmd.numCols][numMod];
		}
		
		FileWriter file3 = new FileWriter(baseFileName + "REPs_num_tissues_genes.txt");
		
		FileWriter fileBTissues = new FileWriter(baseFileName + "_bicluster_tissues.txt");
		FileWriter fileBGenes = new FileWriter(baseFileName + "_bicluster_genes.txt");
		
		for (i=0;i<gmd.numCols;i++) {
			baseIdx = i*entryInc;
			if (REPOccur[i] >= minREPOccurPercent) {
				tissueUse.clear();
				tissueSelect.clear();
				for (j=0;j<tissueMatrix[0].length;j++) {
					v1 = tissueMatrix[baseIdx][j];
					if (v1 >= tissueUseGenesThreshold) {
						tissueSelect.add(j);
						tissueUse.add(-v1);
					}
				}
				if (tissueSelect.size() > 0) {
					genes.clear();
					glist = new geneList();
					for (j=0;j<geneMatrix.length;j++) {
						v1 = geneMatrix[j][i];
						if (v1 >= minGeneOccurPercent) {
							glist.addGene(humanGeneNames[j],v1,1.0);
						}
					}
					
					if (glist.size() >= minGenesInREP) {
						ranks = glist.sort(genes);
						geneSets.addGeneSet(genes,ranks);
						
						generalities.add(new Double(genmd.clusters.get(0).get(i)));
						entriesDouble = new Double[tissueMatrix[0].length];
						for (j=0;j<tissueSelect.size();j++) {
							entriesDouble[tissueSelect.get(j)] = -tissueUse.get(j);
						}
						reviseTissueMatrix.add(entriesDouble);
						
						entriesInteger = new Integer[tissueMatrix[0].length];
						for (j=0;j<tissueSelect.size();j++) {
							maxMod = -1.0;
							maxIdx = -1;
							tissueNum = tissueSelect.get(j);
							if (numMod > 0) {
								for (k=0;k<numMod;k++) {
									v2 = tissueMatrix[baseIdx+k+3][tissueNum];
									if (v2 > maxMod & v2 > 0.0) {
										maxMod = v2;
										maxIdx = k;
									}
								}
							}
							entriesInteger[tissueNum] = maxIdx;
						}
						reviseModMatrix.add(entriesInteger);
					
						sortIdx = VectorUtil.sortOrder(tissueUse);
						file.write("REP#" + ((new Integer(useREPs+1))).toString());
						file.write("\t" + String.format("%.3f",REPOccur[i]));
						file.write("\t" + String.format("%.3f",new Double(genmd.clusters.get(0).get(i))));
						file.write("\t" + ((new Integer(tissueUse.size()))).toString());
						file.write("\t" + ((new Integer(glist.size()))).toString());
						
						file3.write(String.format("%.3f",new Double(genmd.clusters.get(0).get(i))));
						file3.write("\t" + ((new Integer(tissueUse.size()))).toString());
						file3.write("\t" + ((new Integer(glist.size()))).toString());
						file3.write("\n");
						
						file.write("\n");
						file.write("tissues");
						for (k=0;k<sortIdx.length;k++) {
							file.write("\t");
							tissueNum = tissueSelect.get(sortIdx[k]);
							file.write(tmd.experimentNames[tissueNum]);
							fileBTissues.write(tmd.experimentNames[tissueNum]);
							if (k<sortIdx.length - 1) {
								fileBTissues.write("\t");
							}
						}
						file.write("\n");
						fileBTissues.write("\n");
						for (k=0;k<entryInc;k++) {
							file.write(tmd.geneNames[baseIdx + k]);
							for (j=0;j<sortIdx.length;j++) {
								tissueNum = tissueSelect.get(sortIdx[j]);
								file.write("\t" + String.format("%.3f",tissueMatrix[baseIdx+k][tissueNum]));
							}
							file.write("\n");
						}
						file.write("\n");
						
						for (j=0;j<glist.size();j++) {
							file.write(String.format("%.3f",ranks[j]));
							file.write("\t");
							file.write(genes.get(j));
							oGene = myAnnotations.getGene(genes.get(j));
							if (oGene != null) {
								if (outputAnnotations) {
									file.write("\t");
									sterms = oGene.getTermString(namespace);
									if (sterms != null) {
										file.write(oGene.getTermString(namespace));
									}
								}
								if (oGene.Description != null) {
									file.write("\t");
									file.write(oGene.Description);
								}
							}
							file.write("\n");
							
							fileBGenes.write(genes.get(j));
							if (j < glist.size()-1) {
								fileBGenes.write("\t");
							}
						}
						file.write("\n");
						fileBGenes.write("\n");
						
						if (file2 != null) {
							norm = 0.0;
							for (j=0;j<numMod;j++) {
								modAvg[j] = 0.0;
							}
							for (k=0;k<sortIdx.length;k++) {
								tissueNum = tissueSelect.get(sortIdx[k]);
								for (j=0;j<numMod;j++) {
									v1 = tissueMatrix[baseIdx+3+j][tissueNum]*tissueMatrix[baseIdx][tissueNum];
									norm += v1;
									modAvg[j] += v1;
									tissueModAvg[tissueNum][j] += v1;
								}
							}
							file2.write("REP#" + ((new Integer(useREPs+1))).toString());
							for (j=0;j<numMod;j++) {
								modAvg[j] = modAvg[j]/norm;
								file2.write("\t");
								file2.write(String.format("%.3f",modAvg[j]));
							}
							file2.write("\n");
						}
						useREPs++;
					}
				}
			}
		}
		file.close();
		file3.close();
		fileBGenes.close();
		fileBTissues.close();
		
		if (file2 != null) {
			file2.close();
			
			file2 = new FileWriter(baseFileName + "REPs_tissue_mod_avg.txt");
			file2.write("tissue");
			for (i=0;i<numMod;i++) {
				file2.write("\t");
				file2.write(tmd.geneNames[i+3]);
			}
			file2.write("\n");
			for (i=0;i<tmd.numCols;i++) {
				file2.write(tmd.experimentNames[i]);
				norm = 0.0;
				for (j=0;j<numMod;j++) {
					norm += tissueModAvg[i][j];
				}
				for (j=0;j<numMod;j++) {
					file2.write("\t");
					file2.write(String.format("%.3f",tissueModAvg[i][j]/norm));
				}
				file2.write("\n");
			}
			file2.close();
		}
		
		file2 = new FileWriter(baseFileName + "REPs_generality.txt");
		for (i=0;i<generalities.size();i++) {
			file2.write(generalities.get(i).toString());
			if (i < generalities.size() - 1) {
				file2.write("\t");
			}
		}
		file2.close();
		
		file = new FileWriter(baseFileName + "REPs_tissueMatrix_revise.txt");
		file2 = new FileWriter(baseFileName + "REPs_modMatrix_revise.txt");
		
		file.write("tissue");
		file2.write("tissue");
		for (i=0;i<reviseTissueMatrix.size();i++) {
			file.write("\t");
			file.write((new Integer(i+1)).toString());
			file2.write("\t");
			file2.write((new Integer(i+1)).toString());
		}
		file.write("\n");
		file2.write("\n");
		
		for (j=0;j<tmd.numCols;j++) {
			file.write(tmd.experimentNames[j]);
			file2.write(tmd.experimentNames[j]);
			for (i=0;i<reviseTissueMatrix.size();i++) {
				entriesDouble = reviseTissueMatrix.get(i);
				entriesInteger = reviseModMatrix.get(i);
				file.write("\t");
				if (entriesDouble[j] != null) {
					file.write(entriesDouble[j].toString());
				} else {
					file.write("0.0");
				}
				file2.write("\t");
				if (entriesInteger[j] != null) {
					file2.write(entriesInteger[j].toString());
				} else {
					file2.write("-1");
				}
			}
			file.write("\n");
			file2.write("\n");
		}
		file.close();
		file2.close();
		
		geneSets.geneSetSignif();
		file = new FileWriter(baseFileName + "REPs_"+namespace+".txt");
		
		for (i=0;i<geneSets.geneSets.size();i++) {
			geneSets.geneSets.get(i).outputScoresOneLine(file);
		}
		file.close();
		
	}
}
