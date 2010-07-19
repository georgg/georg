package edu.mit.csail.psrg.georg.HDPTopicModel;

import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;
import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.GO.DAGException;
import edu.mit.csail.psrg.georg.GO.GoOBOReader;
import edu.mit.csail.psrg.georg.GO.MouseHumanGoAssocReader;
import edu.mit.csail.psrg.georg.GO.OntClusters;
import edu.mit.csail.psrg.georg.GO.OntDAG;
import edu.mit.csail.psrg.georg.GO.OntGeneAssociations;
import edu.mit.csail.psrg.georg.GO.OntTerm;
import edu.mit.csail.psrg.georg.StatUtil.VectorUtil;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Vector;

public class OutputMouseHumanTopics {
	
	public static void main(String[] args) {
	//	outputFiles();
	//	outputTopGenes();
		correctTopicFile("C:\\research_data\\mouse_human\\b47mm7\\geo\\");
	//	correctGOCategories("C:\\research_data\\mouse_human\\b47mm7\\geo\\");
	}
	
	public static void outputFiles() {
	//	String refDirPath = "C:\\research_data\\mouse_human\\b47mm7\\";
		String dirPath = "C:\\research_data\\mouse_human\\b47mm7\\geo\\alt_clusterings\\cluster_MI1_nogrps\\";
	//	String humanAssocFile = refDirPath + "Hs.genelist.go.plus";
	//	String mouseAssocFile = refDirPath + "Mm.genelist.go.plus";
		String topicInFile = dirPath + "tissue_as_docs_topics.txt";
		String humanTopicOutFileL = dirPath + "tissue_as_docs_topics_go_human_long.txt";
		String humanTopicOutFileS = dirPath + "tissue_as_docs_topics_go_human_short.txt";
		String humanMotifsOutFile = dirPath + "tissue_as_docs_human_motifs.txt";
		String mouseTopicOutFileL = dirPath + "tissue_as_docs_topics_go_mouse_long.txt";
		String mouseTopicOutFileS = dirPath + "tissue_as_docs_topics_go_mouse_short.txt";
		String mouseMotifsOutFile = dirPath + "tissue_as_docs_mouse_motifs.txt";
		String numTopicsOutFile = dirPath + "tissue_as_docs_num_topics.txt";
		String normDPsInFile = dirPath + "tissue_as_docs_dps.txt";
		String normDPsOutFile = dirPath + "tissue_as_docs_dps_norm_rows.txt";
		String normDPsOutFile2 = dirPath + "tissue_as_docs_dps_norm_cols.txt";
//		String GMLOutFile = dirPath + "tissue_as_docs.gml";
		String GMLOutFile = dirPath + "tissue_as_docs.svg";
		String geneCrossOutFile = dirPath + "gene_cross.txt";
		String topGenesOutFile = dirPath + "top_genes.txt";
		String tissueMatrixLoadFile = dirPath + "tissue_load.txt";
		String topicGOCountFile = dirPath + "topic_GO_count.txt";
		String topicCountFile = dirPath + "numtopics.txt";
		
		String snapFileName = dirPath + "tissue_as_docs_snap_149999.persist";
	//	String snapFileName = dirPath + "motif_cluster_snap_round1_5999.persist";
		HierDirichletProcess myDP = HierDirichletProcess.restoreFromFile(snapFileName);
		
	/*	MicroArrayData tempReader = new MicroArrayData();
		try {
			tempReader.readFile(dirPath + "ESNovartis_combined_discretized.txt");
			myDP.geneNames = tempReader.geneNames;
			myDP.persistSelf(snapFileName);
			tempReader = null;
		} catch(IOException e) {
			System.out.println(e);
		} */
		
		try {
			myDP.outputClustersToFile(dirPath + "tissue_as_docs");
		} catch(IOException e) {
			System.out.println(e);
		}
		
		int i = 0;
		int j = 0;
		
		double norm = 0.0;
		MicroArrayData myData = new MicroArrayData();
		try {
			myData.readFile(normDPsInFile);
			
			for (i=0;i<myData.numRows;i++) {
				norm = 0.0;
				for (j=0;j<myData.numCols;j++) {
					norm = norm + myData.values[i][j];
				}
				for (j=0;j<myData.numCols;j++) {
					myData.values[i][j] = myData.values[i][j]/norm;
				}
			}
			myData.writeFile(normDPsOutFile);
			
			myData = new MicroArrayData();
			myData.readFile(normDPsInFile);
			for (j=0;j<myData.numCols;j++) {
				norm = 0.0;
				for (i=3;i<myData.numRows;i++) {
					norm = norm + myData.values[i][j];
				}
				for (i=3;i<myData.numRows;i++) {
					myData.values[i][j] = myData.values[i][j]/norm;
				}
			}
			myData.writeFile(normDPsOutFile2);
			
			
		} catch(IOException e) {
			System.out.println(e);
		}
		
		
	/*	MicroArrayData myData = new MicroArrayData();
		try {
			myData.setDiscrete();
			myData.readFile(topicInFile);
		} catch(IOException e) {
			System.out.println(e);
		} */
		
		String[] humanGeneNames = new String[myDP.geneNames.length];
		String[] mouseGeneNames = new String[myDP.geneNames.length];
		String s1 = "";
		String s2 = "";
		int f = 0;
		
		for (i=0;i<myDP.geneNames.length;i++) {
			f = myDP.geneNames[i].indexOf("|");
			s1 = myDP.geneNames[i].substring(0,f);
			s2 = myDP.geneNames[i].substring(f+1);
			humanGeneNames[i] = s1;
			mouseGeneNames[i] = s2;
		}
		
		if (myDP.DPGroups.size() > 0) {
			if (myDP.DPGroups.get(0).numGoodIter > 0) {
				reconstructGroups(myDP);
			}
		}
		CrossSpeciesStrongModuleSummarizer topicSumm = new CrossSpeciesStrongModuleSummarizer(myDP);
		
		boolean collapseModules = false;
		
		if (myDP.statCollectors.containsKey(StrongModuleCollector.class)) {
			StrongModuleCollector ms = (StrongModuleCollector) myDP.getStatCollector(StrongModuleCollector.class);
			if (ms.numSamples > 0) {
				collapseModules = true;
			}
		}
		
	//	topicSumm.outputTopics(humanTopicOutFileL,humanTopicOutFileS,GMLOutFile,geneCrossOutFile,topGenesOutFile,tissueMatrixLoadFile,topicGOCountFile,true);
		topicSumm.outputTopics(humanTopicOutFileL,humanTopicOutFileS,GMLOutFile,geneCrossOutFile,topGenesOutFile,tissueMatrixLoadFile,topicGOCountFile,collapseModules);
		
	//	topicSumm.outputTopics(humanTopicOutFileL,humanTopicOutFileS,GMLOutFile,geneCrossOutFile,topGenesOutFile,false);
		
/*		GoOBOReader goReader = null;
		MouseHumanGoAssocReader assocReader = null;
		StrongModuleSummarizer topicSumm = null;
		
		// output annotated human topics
		goReader = new GoOBOReader();
		assocReader = new MouseHumanGoAssocReader();
		assocReader.setFile(humanAssocFile);
		myDP.geneNames = humanGeneNames;
		
		reconstructGroups(myDP);
		topicSumm = new StrongModuleSummarizer(myDP,goReader,assocReader);
		topicSumm.outputTopics(humanTopicOutFileL,humanTopicOutFileS,GMLOutFile,humanMotifsOutFile,true,true); */ 
		
		// output annotated mouse topics
/*		goReader = new GoOBOReader();
		assocReader = new MouseHumanGoAssocReader();
		assocReader.setFile(mouseAssocFile);
		myDP = HierDirichletProcess.restoreFromFile(snapFileName);
		myDP.geneNames = mouseGeneNames;
	
		reconstructGroups(myDP);
		topicSumm = new StrongModuleSummarizer(myDP,goReader,assocReader);
		topicSumm.outputTopics(mouseTopicOutFileL,mouseTopicOutFileS,GMLOutFile,mouseMotifsOutFile,true,false);
		*/
		
		try {
			myDP.outputNumClustersSamples(topicCountFile);
		} catch(IOException e) {
			System.out.println(e);
		}
		
		correctTopicFile(dirPath);
	}
	
	public static void reconstructGroups(HierDirichletProcess myDP) {
		if (myDP.DPGroups.size() == 0)
			return;
		
		int i = 0;
		for (i=0;i<myDP.DPGroups.size();i++) {
			myDP.DPGroups.get(i).normPairProbs();
			myDP.DPGroups.get(i).buildConcensusModel();
		}
		myDP.buildDPList();
		System.out.println("reconstructed groups");
	}
	
	public static void outputTopGenes() {
		String dirPath = "C:\\research_data\\mouse_human\\b47mm7\\geo\\";
		String topicFileName = dirPath + "corrected_topics_go.txt";
		String topGenesFileName = dirPath + "top_genes_50.txt";
		
		LinkedHashMap<String,ArrayList<String>> topicMap = null;
		String[] fields = null;
		ArrayList<String> geneList = null;
		String gene = null;
		String topicName = null;
		
		try {
			topicMap = parseTopicFile(topicFileName);
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		Iterator<String> topicIter = topicMap.keySet().iterator();
		int i = 0;
		double humanPercent = 0.0;
		double mousePercent = 0.0;
		
		try {
			FileWriter file = new FileWriter(topGenesFileName);
			while(topicIter.hasNext()) {
				topicName = topicIter.next();
				geneList = topicMap.get(topicName);
				fields = topicName.split("\t");
				topicName = fields[0];
				file.write(topicName);
				for (i=0;i<geneList.size();i++) {
					fields = geneList.get(i).split("\t");
					humanPercent = (new Double(fields[3])).doubleValue();
					mousePercent = (new Double(fields[5])).doubleValue();
				/*	if (humanPercent > 0.0 & mousePercent > 0.0) {
						// write human gene
						file.write("\t");
						file.write(fields[0]);
						// write mouse gene
						file.write("\t");
						file.write(fields[1]);
					} */
					if (humanPercent >= 0.50) {
						file.write("\t");
						file.write(fields[0]);
					}
					if (mousePercent >= 0.50) {
						file.write("\t");
						file.write(fields[1]);
					}
				}
				file.write("\n");
			}
			file.close();
		} catch(IOException e) {
			e.printStackTrace();
		}
	}
	
	public static LinkedHashMap<String,ArrayList<String>> parseTopicFile(String fileName) throws IOException {
		BufferedReader is = new BufferedReader(new FileReader(fileName));
		int f = 0;
		String[] fields = null;
		String line = null;
		String line2 = null;
		int readingGenes = 1;
		int readingGO = 2;
		int state = readingGO;
		LinkedHashMap<String,ArrayList<String>> genesMap = new LinkedHashMap<String,ArrayList<String>>();
		ArrayList<String> geneList = null;
		
		while(is.ready()) {
			line = is.readLine();
			if (state == readingGO) {
				f = line.indexOf("Cluster#");
				if (f == 0) {
					state = readingGenes;
					geneList = new ArrayList<String>();
					line2 = is.readLine();
					genesMap.put(line + "\t" + line2,geneList);
				}
			} else {
				if (state == readingGenes) {
					f = line.indexOf("===");
					if (f==0) {
						state = readingGO;
					} else {
						geneList.add(line);
					}
				}
			}
		}
		
		is.close();
		return genesMap;
	}
	
	public static void correctGOCategories(String dirPath) {
		String topicFileName = dirPath + "corrected_topics_go.txt";
		String GOLongFileName = dirPath + "topics_go_long.txt";
		String GOCountFileName = dirPath + "topic_GO_count.txt";
		
		int minGenes = 5;
		int maxGenes = 200;
		// minimum count of gene in topic
		int minOccurInTopic = 2;
		double FDR = 0.05;
		int namespace = OntTerm.BiologicalProcess;
		double minOccurPercent = 0.0;
		
		String line = "";
		String[] split = null;
		int i = 0;
		int j = 0;
		int f = 0;
		boolean nextTopic = true;
		double[] geneValues = null;
		double v = 0;
		String[] names = null;
		
		String refDirPath = "C:\\research_data\\mouse_human\\";
		String humanAssocFile = refDirPath + "Hs.genelist.go.plus";
		String mouseAssocFile = refDirPath + "Mm.genelist.go.plus";
		
		GoOBOReader goReader = new GoOBOReader();
		MouseHumanGoAssocReader humanAssocReader = new MouseHumanGoAssocReader();
		humanAssocReader.setFile(humanAssocFile);
		OntDAG humanDAG = null;
		OntGeneAssociations humanAssoc = null;
		
		MouseHumanGoAssocReader mouseAssocReader = new MouseHumanGoAssocReader();
		mouseAssocReader.setFile(mouseAssocFile);
		OntDAG mouseDAG = null;
		OntGeneAssociations mouseAssoc = null;
		
		HashSet<String> humanUseGenes = new HashSet<String>();
		HashSet<String> mouseUseGenes = new HashSet<String>();
		ArrayList<GeneList> geneLists = new ArrayList<GeneList>();
		ArrayList<String> topicNames = new ArrayList<String>();
		GeneList myGeneList = null;
		
		try {
			BufferedReader readFile = new BufferedReader(new FileReader(topicFileName));
			line = readFile.readLine();
			while(readFile.ready()) {
				split = line.split("\t");
				topicNames.add(split[0]);
				myGeneList = new GeneList();
				geneLists.add(myGeneList);
			
				nextTopic = false;
				while(!nextTopic & readFile.ready()) {
					line = readFile.readLine();
					f = line.indexOf("Cluster#");
					if (f == 0) {
						nextTopic = true;
					} else {
						f = line.indexOf("NM_");
						// parse genes
						if (f == 0) {
							split = line.split("\t");
							humanUseGenes.add(split[0]);
							mouseUseGenes.add(split[1]);
							v = (new Double(split[2])).doubleValue();
							geneValues = new double[4];
							for (i=3;i<7;i++) {
								geneValues[i-3] = (new Double(split[i])).doubleValue();
							}
							if (geneValues[0] >= minOccurPercent | geneValues[2] >= minOccurPercent) {
								names = new String[2];
								names[0] = split[0];
								names[1] = split[1];
								myGeneList.addGene(names,v,geneValues);
							}
						}
					}
				}
			}
			readFile.close();
		} catch(IOException e) {
			System.out.println(e);
		}
		
		try {
			humanDAG = goReader.readFile();
			System.out.println("Loaded human DAG");
			humanAssoc = humanAssocReader.readFile(humanDAG,humanUseGenes);
			System.out.println("Loaded human associations");
			
			mouseDAG = goReader.readFile();
			System.out.println("Loaded mouse DAG");
			mouseAssoc = mouseAssocReader.readFile(mouseDAG,mouseUseGenes);
			System.out.println("Loaded mouse associations");
		} catch(IOException e) {
			System.out.println(e);
		} catch(DAGException e) {
			System.out.println(e);
		}
		
		OntClusters humanClusters = new OntClusters(humanAssoc,minGenes,maxGenes,FDR,namespace);
		OntClusters mouseClusters = new OntClusters(mouseAssoc,minGenes,maxGenes,FDR,namespace);
		ArrayList<String> humanGeneSetNames = null;
		ArrayList<String> mouseGeneSetNames = null;
		
		for (i=0;i<geneLists.size();i++) {
			humanGeneSetNames = new ArrayList<String>();
			mouseGeneSetNames = new ArrayList<String>();
			myGeneList = geneLists.get(i);
			myGeneList.sort();
			for (j=0;j<myGeneList.size();j++) {
				humanGeneSetNames.add(myGeneList.genes.get(j).names[0]);
				mouseGeneSetNames.add(myGeneList.genes.get(j).names[1]);
			}
			humanClusters.addCluster(humanGeneSetNames,myGeneList.getRanks(),"none",0.0);
			mouseClusters.addCluster(mouseGeneSetNames,myGeneList.getRanks(),"none",0.0);
		}
		
		humanClusters.clusterSignif();
		mouseClusters.clusterSignif();
		
		int[] humanSignif = humanClusters.numSignificantTerms();
		int[] mouseSignif = mouseClusters.numSignificantTerms();
		try {
			FileWriter gf = new FileWriter(GOCountFileName);
			for (i=0;i<humanSignif.length;i++) {
				gf.write((new Integer(humanSignif[i])).toString());
				gf.write("\t");
				gf.write((new Integer(mouseSignif[i])).toString());
				if (humanSignif[i] > 0 | mouseSignif[i] > 0) {
					gf.write("\t1\n");
				} else {
					gf.write("\t0\n");
				}
			}
			gf.close();
		} catch(IOException e) {
			System.out.println(e);
		}
		
		try {
			FileWriter outFile = new FileWriter(GOLongFileName);
			for (i=0;i<geneLists.size();i++) {
				outFile.write(topicNames.get(i));
				outFile.write("\n");
				humanClusters.getClusters().get(i).outputScoresLong(outFile,false);
				mouseClusters.getClusters().get(i).outputScoresLong(outFile,false);
				outFile.write("\n");
			}
			outFile.close();
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	public static void correctTopicFile(String dirPath) {
		String oldTopicFileName = dirPath + "tissue_as_docs_topics_go_human_long.txt";
		String newTopicFileName = dirPath + "corrected_topics_go.txt";
		String tissueOutputName = dirPath + "HDP_tissues.txt";
		String geneOutputName = dirPath + "HDP_genes.txt";
		String line = "";
		double topicPercent = 0.0;
		String[] split = null;
		ArrayList<String> tissueNames = new ArrayList<String>();
		String tempString1 = "";
		String tempString2 = "";
		int i = 0;
		int f = 0;
		int f2 = 0;
		Vector<Double> tissueLoads = new Vector<Double>();
		int[] tissueOrder = null;
		boolean nextTopic = true;
		double geneValue = 0.0;
		double v = 0;
		
		ArrayList<ArrayList<String>> geneClusters = new ArrayList<ArrayList<String>>();
		ArrayList<String> myGeneCluster = new ArrayList<String>();
		geneClusters.add(myGeneCluster);
		
		try {
			BufferedReader readFile = new BufferedReader(new FileReader(oldTopicFileName));
			FileWriter writeFile = new FileWriter(newTopicFileName);
			FileWriter tissuesFile = new FileWriter(tissueOutputName);
			line = readFile.readLine();
			while(readFile.ready()) {
				split = line.split("\t");
				topicPercent = (new Double(split[2])).doubleValue();
				if (topicPercent > 1.0)
				{
					split[2] = "1.0";
				} else {
					split[2] = (new Double(topicPercent)).toString();
				}
				line = split[0] + "\t" + split[1] + "\t" + split[2] + "\n";
				writeFile.write(line);
				line = readFile.readLine();
				
				tissueLoads.clear();
				tissueNames.clear();
				
				f = line.indexOf(" | ");
				while(f > 1) {
					tempString1 = line.substring(0,f-1);
					line = line.substring(f+3);
					f2 = tempString1.lastIndexOf(" (");
					tissueNames.add(tempString1.substring(0,f2));
					tempString2 = tempString1.substring(f2+2);
					tempString2 = tempString2.substring(0,tempString2.length()-1);
					tissueLoads.add(new Double(tempString2));
					f = line.indexOf(" | ");
				}
				for (i=0;i<tissueLoads.size();i++) {
					tissueLoads.set(i,-tissueLoads.get(i));
					tissuesFile.write(tissueNames.get(i));
					if (i < tissueLoads.size() - 1) {
						tissuesFile.write("\t");
					}
				}
				tissuesFile.write("\n");
				tissueOrder = VectorUtil.sortOrder(tissueLoads);
				for (i=0;i<tissueOrder.length;i++) {
					writeFile.write(tissueNames.get(tissueOrder[i]));
					v = -tissueLoads.get(tissueOrder[i]);
					writeFile.write(" (" + (new Double(v)).toString() + ")");
					if (i < tissueOrder.length - 1) {
						writeFile.write(" | ");
					}
				}
				writeFile.write("\n");
			
				nextTopic = false;
				while(!nextTopic & readFile.ready()) {
					line = readFile.readLine();
					f = line.indexOf("Cluster#");
					if (f == 0) {
						nextTopic = true;
						myGeneCluster = new ArrayList<String>();
						geneClusters.add(myGeneCluster);
					} else {
						f = line.indexOf("NM_");
						// parse genes
						if (f == 0) {
							split = line.split("\t");
							myGeneCluster.add(split[0] + "|" + split[1]);
							if (topicPercent > 1.0) {
								for (i=2;i<7;i++) {
									geneValue = (new Double(split[i])).doubleValue();
									geneValue = geneValue * topicPercent;
									split[i] = (new Double(geneValue)).toString();
								}
							}
							for (i=0;i<split.length;i++) {
							//	if (i < 4 | i > 6) {
									writeFile.write(split[i]);
									if (i < split.length - 1) {
										writeFile.write("\t");
									}
							//	}
							}
							writeFile.write("\n");
						} else {
							writeFile.write(line);
							writeFile.write("\n");
						}
					}
				}
			}
			readFile.close();
			writeFile.close();
			tissuesFile.close();
			
			ClusterReader gc = new ClusterReader(geneClusters);
			gc.writeFile(geneOutputName);
		} catch(IOException e) {
			System.out.println(e);
		}
		correctGOCategories(dirPath);
	}
}
