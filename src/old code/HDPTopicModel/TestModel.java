package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.IOException;
import java.util.ArrayList;

import edu.mit.csail.psrg.georg.DataAccess.*;

public class TestModel {
	public static void main(String[] args) {
	//	String fname = args[0];
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
		HDPConcParam alpha_0 = new HDPConcParam(1.0,1.0);
		concParams[0] = gamma_p;
		concParams[1] = alpha_0;
		
		Document doc = null;
		ArrayList<Integer> glist = new ArrayList<Integer>();
		ArrayList<DirichletProcess> DPs = new ArrayList<DirichletProcess>();
		ArrayList<DirichletProcess> alpha_0_DP = new ArrayList<DirichletProcess>();
 		
		DirichletProcess root = new DirichletProcess(null,"root");
		gamma_p.addDP(root);
		DPs.add(root);
		
		DirichletProcess d = null;
		
		for (i=0;i<cl.length;i++) {
			if (cl[i] > 0.0) {
				glist.clear();
				for (j=0;j<myData.numRows;j++) {
					if (myData.values[j][i] == 1.0) {
						glist.add(j);
					}
				}
				doc = new Document(myData.experimentNames[i],glist);
				d = new DirichletProcess(root,myData.experimentNames[i]);
				d.addDocument(doc);
				DPs.add(d);
				alpha_0_DP.add(d);
			}
		}
		
		alpha_0.addDPs(alpha_0_DP);
		
		HierDirichletProcess myHierDP = new HierDirichletProcess(1,DPs,concParams,myData.geneNames);
		
		myHierDP.packGenes();
		myHierDP.activateAll();
		
		int iters = 10000;
		int burnin = 5000;
		
		myHierDP.iterate(iters,burnin);
		
	//	myHierDP.normPairProbs(iters-burnin);
	//	myHierDP.concensusModel();
	//	myHierDP.output();
		try {
			myHierDP.outputClustersToFile("C:\\CPPDataFiles\\test_cluster_output.txt");
		} catch(IOException e) { System.out.println(e);	}
	}
}
