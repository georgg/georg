package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.util.ArrayList;

public class OrthoHierDirichletProcess extends HierDirichletProcess {
	public OrthologStats myOrthologStats = new OrthologStats();
	
	OrthoHierDirichletProcess(int nc, ArrayList<DirichletProcess> dd, HDPConcParam[] cp, String[] gn) {
		super(nc, dd, cp, gn);
	}
	
	void iterate(int numIters,int burnin,int concParamIter,boolean calcLikelihood) {
  		int iters = 0;
  		int jj = 0;
  		int cc = 0;
  		int cp = 0;
  		double likelihood = 0.0;
  		double maxLikelihood = -1.0/0.0;
  		int nc = 0;
  		
  		String snapShotOutName = "";
  		
  		for (iters=0;iters<numIters;iters++) {
  			for (jj=DP.length-1;jj>=0;jj--) {
  				if (DP[jj].documents != null & DP[jj].state == DirichletProcess.ACTIVE) {
  					sampleDataClusterAssigns(DP[jj]);
  				}
  				if (DP[jj].state == DirichletProcess.ACTIVE) {
  					DP[jj].sampleClusterNumTables(numClusters);
  				}
  			}
  			
  			for (cp=0;cp<concParams.length;cp++) {
  				concParams[cp].calcTotals();
  			}
  			
  			for (jj=0;jj<DP.length;jj++) {
  				if (DP[jj].state == DirichletProcess.ACTIVE) {
  					DP[jj].resampleBeta(numClusters,condLike);
  				}
  			} 
  			
  			// now delete any empty clusters
  			nc = numClusters;
  			for (cc=nc-1;cc>=0;cc--) {
  				if (clusterCounts[cc] == 0) {
  					deleteCluster(cc);
  				}
  			}
  			
  			int bicluster = 1;
  			if (!DPGroups.isEmpty() & iters > groupStartIter) {
  				for (cc=0;cc<bicluster;cc++) {
	  				for (jj=0;jj<DPGroups.size();jj++) {
	  					DPGroups.get(jj).sampleDPAssignments(numClusters);
	  					DPGroups.get(jj).resampleAlpha(concParamIter);
	  					System.out.println("grps=" + DPGroups.get(jj).numClusters + " alpha=" + DPGroups.get(jj).alpha);
	  				}
	  				buildDPList();
	  				
	  				for (jj=DP.length-1;jj>=0;jj--) {
	  					if (DP[jj].state == DirichletProcess.ACTIVE) {
	  						DP[jj].sampleClusterNumTables(numClusters);
	  					}
	  	  			}
	  				for (cp=0;cp<concParams.length;cp++) {
	  	  				concParams[cp].calcTotals();
	  	  			}
	  				for (jj=0;jj<DP.length;jj++) {
	  					if (DP[jj].state == DirichletProcess.ACTIVE) {
	  						DP[jj].resampleBeta(numClusters,condLike);
	  					}
		  			}
	  			//	resampleConcParams(concParamIter); 
	  			/*	for (cp=0;cp<concParams.length;cp++) {
	  	  				concParams[cp].calcTotals();
	  	  			} */
  				}
  			}
  			
  			resampleConcParams(concParamIter);
  			resamplePriorPos(concParamIter);
  			System.out.println(iters + " " + numClusters);
  			
  			if (iters >= burnin) {
  				if (calcLikelihood) {
  					logLikeSamples.add(calcLogLike());
  				}
  				numClustersSamples.add(numClusters);
  				if (collectStrongModules) {
  					updateClusterStatistics();
  				}
  				if (!DPGroups.isEmpty()) {
  					for (jj=0;jj<DPGroups.size();jj++) {
	  					DPGroups.get(jj).incPairProbs();
	  				}
  				}
  				if (myOrthologStats != null) {
  					myOrthologStats.increment();
  				}
  			}
  			
  			if (snapShotFileName != null) {
  				double t = ((double) iters+1)/((double) snapShotInterval);
  				if (((double) Math.round(t)) == t) {
  					snapShotOutName = snapShotFileName + "_" + (new Integer(iters)).toString() + ".persist";
  					System.out.println("Snap shot at iteration" + iters);
  					persistSelf(snapShotOutName);
  				}
  			}
  		}
  	}

}
