package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.GO.DAGException;
import edu.mit.csail.psrg.georg.GO.DAGReader;
import edu.mit.csail.psrg.georg.GO.GoOBOReader;
import edu.mit.csail.psrg.georg.GO.MIPSDAG;
import edu.mit.csail.psrg.georg.GO.OntAssociationReader;
import edu.mit.csail.psrg.georg.GO.OntClusters;
import edu.mit.csail.psrg.georg.GO.OntDAG;
import edu.mit.csail.psrg.georg.GO.OntGeneAssociations;
import edu.mit.csail.psrg.georg.GO.OntScore;
import edu.mit.csail.psrg.georg.GO.OntTerm;

public class TopicSummarizer {
		OntGeneAssociations assoc = null;
		OntDAG DAG = null;
		int minGenes = 5;
		int maxGenes = 200;
		// minimum count of gene in topic
		int minOccurInTopic = 2;
		double FDR = 0.05;
	//	int namespace = OntTerm.AnyCategory;
		int namespace = OntTerm.BiologicalProcess;
		HashSet<String> useGenes = new HashSet<String>();
		MicroArrayData myData = null;
		
		private class geneHolder implements Comparable {
			String name = "";
			double rank = 0.0;
			
			public geneHolder(String n,double r) {
				name = n;
				rank = r;
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
			
			void addGene(String name,double rank) {
				genes.add(new geneHolder(name,rank));
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
		
		public TopicSummarizer(MicroArrayData init_MicroArrayData,DAGReader init_DAGReader, OntAssociationReader init_ontAssociationReader) {
			processGenes(init_MicroArrayData);
			readDAG(init_DAGReader);
			readGeneAssocations(init_ontAssociationReader);
		}
		
		public void setMinOccurInTopic(int m) {
			minOccurInTopic = m;
		}
		
		void processGenes(MicroArrayData init_data) {
			myData = init_data;
			int i = 0;
			for (i=0;i<myData.numRows;i++) {
				useGenes.add(myData.geneNames[i]);
			}
		}
		
		void readDAG(DAGReader reader) {
			try {
					DAG = reader.readFile();
					System.out.println("Loaded DAG");
				} catch(IOException e) {
					System.out.println(e);
				} catch(DAGException e) {
					System.out.println(e);
				}
		}
		
		void readGeneAssocations(OntAssociationReader reader) {
			try {
				assoc = reader.readFile(DAG,useGenes);
				System.out.println("Loaded associations");
			} catch(IOException e) {
				System.out.println(e);
			}
		}
		
		public void outputTopics(String fOutNameLong,String fOutNameShort) {
			int j = 0;
			int i = 0;
			ArrayList<String> genes = new ArrayList<String>();
			OntClusters clusters = new OntClusters(assoc,minGenes,maxGenes,FDR,namespace);
			geneList glist = null;
			double[] ranks = null;
			double norm = 0.0;
			
			for(j=0;j<myData.numCols;j++) {
				genes.clear();
				glist = new geneList();
				norm = 0.0;
				for (i=0;i<myData.numRows;i++) {
					if (myData.dvalues[i][j] > minOccurInTopic) {
						norm = norm + (double) myData.dvalues[i][j];
					}
				}
				for (i=0;i<myData.numRows;i++) {
					if (myData.dvalues[i][j] > minOccurInTopic) {
					//	genes.add(myData.geneNames[i]);
					//	glist.addGene(myData.geneNames[i],((double) myData.dvalues[i][j])/norm);
						glist.addGene(myData.geneNames[i],myData.dvalues[i][j]);
					}
				}
				ranks = glist.sort(genes);
				clusters.addCluster(genes,ranks);
			}
			clusters.clusterSignif();
			
			try {
				clusters.outputClustersLong(fOutNameLong);
				clusters.outputClustersShort(fOutNameShort);
			} catch(IOException e) {
				System.out.println(e);
			}
		}
		
}
