package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;

import cern.jet.random.Uniform;

public class OrthologStats implements Serializable {
	public DirichletProcess[] humanDPs = null;
	public DirichletProcess[] mouseDPs = null;
	public HierDirichletProcess myHDP = null;
	public double[][] similarity = null;
	int numIters = 0;
	
	public void init(HierDirichletProcess HDP,String[] humanOrthologs,String[] mouseOrthologs) {
		humanDPs = new DirichletProcess[humanOrthologs.length];
		mouseDPs = new DirichletProcess[mouseOrthologs.length];
		myHDP = HDP;
		similarity = new double[humanOrthologs.length][4];
		
		HashMap<String,DirichletProcess> DPMap = new HashMap<String,DirichletProcess>();
		int i = 0;
		
		for (i=0;i<myHDP.DP.length;i++) {
			DPMap.put(myHDP.DP[i].label,myHDP.DP[i]);
		}
		
		DirichletProcess dp = null;
		for (i=0;i<humanOrthologs.length;i++) {
			dp = DPMap.get(humanOrthologs[i]);
			humanDPs[i] = dp;
			dp = DPMap.get(mouseOrthologs[i]);
			mouseDPs[i] = dp;
		}
	}
	
	public void increment() {
		numIters++;
		int i = 0;
		DirichletProcess mouseDP = null;
		DirichletProcess humanDP = null;
		int match = 0;
		double dist = 0.0;
		double symDist = 0.0;
		
		for (i=0;i<humanDPs.length;i++) {
			humanDP = humanDPs[i];
			mouseDP = mouseDPs[i];
			dist = hellingerDistance(humanDP,mouseDP);
			similarity[i][0] += hellingerDistance(humanDP,mouseDP);
			similarity[i][1] += Math.pow(dist,2.0);
			
			if (!Double.isNaN(similarity[i][0])) {
				dist = 0.0/0.0;
				while (Double.isNaN(dist)) {
					match = Uniform.staticNextIntFromTo(0,similarity.length-1);
					dist = hellingerDistance(humanDP,mouseDPs[match]);
				}
				symDist = dist;
				
				dist = 0.0/0.0;
				while (Double.isNaN(dist)) {
					match = Uniform.staticNextIntFromTo(0,similarity.length-1);
					dist = hellingerDistance(humanDPs[match],mouseDP);
				}
				symDist += symDist;
				dist = symDist/2.0;
				similarity[i][2] += dist;
				similarity[i][3] += Math.pow(dist,2.0);
			}
		}
	}
	
	public double hellingerDistance(DirichletProcess humanDP,DirichletProcess mouseDP) {
		int cc = 0;
		int numTopics = myHDP.numClusters;
		int humanNorm = 0;
		int mouseNorm = 0;
		double dist = 0.0;
		double norm = 0.0;
		
		dist = 0.0;
		for (cc=0;cc<numTopics;cc++) {
			humanNorm += humanDP.clusterNumData[cc];
			mouseNorm += mouseDP.clusterNumData[cc];
			// compute the Hellinger distance between the distributions
			dist += Math.sqrt((double) (humanDP.clusterNumData[cc]*mouseDP.clusterNumData[cc]));
		}
		norm = Math.sqrt((double) (humanNorm*mouseNorm));
		dist = dist/norm;
		
		return dist;
	}
	
	public void normalize() {
		int i = 0;
		int j = 0;
		double norm = (double) numIters;
		
		for (i=0;i<similarity.length;i++) {
			for (j=0;j<4;j++) {
				similarity[i][j] = similarity[i][j]/norm;
			}
		}
	}
	
	public void outputSimilarities(String fName) throws IOException {
		normalize();
		FileWriter file = new FileWriter(fName);
		int i = 0;
		int j = 0;
		for (i=0;i<similarity.length;i++) {
			file.write(humanDPs[i].label);
			file.write("\t");
			file.write(mouseDPs[i].label);
			for (j=0;j<4;j++) {
				file.write("\t");
				file.write((new Double(similarity[i][j])).toString());
			}
			file.write("\n");
		}
		file.close();
	}
}
