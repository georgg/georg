package com.newmacondo.georg.DynamicDUO;

import umontreal.iro.lecuyer.randvar.NegativeBinomialGen;
import umontreal.iro.lecuyer.rng.MRG32k3a;

import com.newmacondo.georg.StatUtil.StatUtil;

public class SynthesizeData {
	
	public DP init() {
		int numClusters = 5;
		DP myDP = new DP();
		myDP.numOTUs = 50;
		myDP.numLocations = 3;
		myDP.numTimepoints = 5;
		myDP.numReplicates = 5;
		DynamicsCluster cluster = null;
		
		myDP.OTUCounts = new int[myDP.numOTUs][myDP.numLocations][myDP.numTimepoints][myDP.numReplicates];
		myDP.OTUNames = new String[myDP.numOTUs];
		myDP.states = new int[myDP.numOTUs][myDP.numLocations][myDP.numTimepoints];
		myDP.OTUAssignment = new DynamicsCluster[myDP.numOTUs][myDP.numLocations];
		myDP.negBinomR = new double[myDP.numStates];
		myDP.negBinomP = new double[myDP.numStates];
		
		int i = 0;
		int j = 0;
		int k = 0;
		int l = 0;
		
		int state = 0;
		
		double[] clusterProb = new double[numClusters];
		double[] alpha = new double[numClusters];
		for (i=0;i<numClusters;i++) {
			alpha[i] = 1.0;
			cluster = new DynamicsCluster(myDP);
			myDP.clusters.add(cluster);
			cluster.sampleTransitionProbsPosterior();
		}
		StatUtil.dirichlet_rnd(clusterProb, alpha, numClusters);
		int assignment = 0;
		
		for (i=0;i<myDP.numOTUs;i++) {
			assignment = StatUtil.multinomial_rnd(clusterProb, numClusters);
			for (j=0;j<myDP.numLocations;j++) {
				myDP.OTUAssignment[i][j] = myDP.clusters.get(assignment);
			}
		}
		
		myDP.negBinomR[0] = 5.0;	myDP.negBinomP[0] = 5.0/(10.0+5.0);
		myDP.negBinomR[1] = 20.0;	myDP.negBinomP[1] = 20.0/(20.0+50.0);
		myDP.negBinomR[2] = 40.0; 	myDP.negBinomP[2] = 40.0/(40.0+100.0);
		myDP.negBinomR[3] = 10.0; 	myDP.negBinomP[3] = 10.0/(200.0+10.0);
		
		MRG32k3a rStream = new MRG32k3a();
		
		for (i=0;i<myDP.numOTUs;i++) {
			myDP.OTUNames[i] = "OTU " + Integer.toString(i+1);
			for (j=0;j<myDP.numLocations;j++) {
				cluster = myDP.OTUAssignment[i][j];
				state = StatUtil.multinomial_rnd(cluster.initProbs, myDP.numStates);
				myDP.states[i][j][0] = state;
				for (k=1;k<myDP.numTimepoints;k++) {
					state = StatUtil.multinomial_rnd(cluster.transitionProbs[state], myDP.numStates);
					myDP.states[i][j][k] = state;
				}
				for (k=0;k<myDP.numTimepoints;k++) {
					state = myDP.states[i][j][k];
					for (l=0;l<myDP.numReplicates;l++) {
						myDP.OTUCounts[i][j][k][l] = NegativeBinomialGen.nextInt(rStream, myDP.negBinomR[state], myDP.negBinomP[state]);
					}
				}
			}
		}
		
		return myDP;
	}
	
}
