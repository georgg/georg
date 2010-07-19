package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.IOException;
import java.util.ArrayList;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;

public class TestModelGroups {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		String fname = "C:\\CPPDataFiles\\test_docs.txt";
		
		MicroArrayData myData = new MicroArrayData();
		try {
			myData.readFile(fname);
		} catch(IOException e) {
			System.out.println(e);
			System.exit(0);
		}
		
		double[] cl = myData.meanCols();
		
		int i = 0;
		int j = 0;
		
		HDPConcParam[] concParams = new HDPConcParam[2];
		HDPConcParam gamma_p = new HDPConcParam(1.0,1.0);
	//	HDPConcParam alpha_1 = new HDPConcParam(1.0,1.0);
		HDPConcParam alpha_0 = new HDPConcParam(1.0,1.0);
		concParams[0] = gamma_p;
	//	concParams[1] = alpha_1;
		concParams[1] = alpha_0;
		
		Document doc = null;
		ArrayList<Integer> glist = new ArrayList<Integer>();
		ArrayList<DirichletProcess> DPs = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> alpha_0_DP = new ArrayList<DirichletProcess>();
	//	ArrayList<DirichletProcess> alpha_1_DP = new ArrayList<DirichletProcess>();
 		
		DirichletProcess root = new DirichletProcess(null,"root");
		gamma_p.addDP(root);
		DPs.add(root);
		
		DirichletProcess d0 = null;
		DirichletProcess d1 = null;
		
		boolean groupDocuments = true;
		
		if (!groupDocuments) {
			d1 = new DirichletProcess(root,"");
			alpha_0_DP.add(d1);
		}
		
		ArrayList<Document> docs = new ArrayList<Document>();
		
		for (i=0;i<cl.length;i++) {
			if (cl[i] > 0.0) {
				glist.clear();
				for (j=0;j<myData.numRows;j++) {
					if (myData.values[j][i] == 1.0) {
						glist.add(j);
					}
				}
				doc = new Document(myData.experimentNames[i],glist);
				if (groupDocuments) {
					docs.add(doc);
				}
				if (!groupDocuments) {
					d0 = new DirichletProcess(d1,"DP_" + myData.experimentNames[i]);
					d0.addDocument(doc);
					alpha_0_DP.add(d0);
				}
			}
		}
		
		if (groupDocuments) {
			d0 = new DirichletProcess(root,"");
			d0.addDocuments(docs);
			alpha_0_DP.add(d0);
		}
		
	/*	for (i=0;i<alpha_1_DP.size();i++) {
			DPs.add(alpha_1_DP.get(i));
		} */
		
		for (i=0;i<alpha_0_DP.size();i++) {
			DPs.add(alpha_0_DP.get(i));
		}
		
		alpha_0.addDPs(alpha_0_DP);
	//	alpha_1.addDPs(alpha_1_DP);
		
		HierDirichletProcess myHierDP = new HierDirichletProcess(2,DPs,concParams,myData.geneNames);
		DPGroup group1 = new DPGroup(root,groupDocuments);
		myHierDP.addDPGroup(group1);
		
		myHierDP.packGenes();
		myHierDP.activateAll();
		
		int iters = 20000;
		int burnin = 5000;
		
		myHierDP.iterate(iters,burnin);
		group1.normPairProbs();
		
		try {
			myHierDP.outputClustersToFile("C:\\CPPDataFiles\\test_cluster_output.txt");
			group1.writeConcensusModel("C:\\CPPDataFiles\\test_doc_groups.txt");
		} catch(IOException e) { System.out.println(e);	}
	}
}
