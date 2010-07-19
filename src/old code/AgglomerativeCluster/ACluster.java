package edu.mit.edu.csail.psrg.georg.AgglomerativeCluster;

import java.util.ArrayList;

public class ACluster {
	public ArrayList<Integer> dataItems = new ArrayList<Integer>();
	public double lik = 0.0;
	public double priorCount = 1.0;
	public AgglomerativeCluster myClusterer = null;
	
	public ACluster(AgglomerativeCluster a) {
		myClusterer = a;
	}
	
	public ACluster(AgglomerativeCluster a,int i) {
		myClusterer = a;
		addItem(i);
	}
	
	public static ACluster merge(ACluster a1,ACluster a2) {
		ACluster a3 = new ACluster(a1.myClusterer);
		
		a3.dataItems.addAll(a1.dataItems);
		a3.dataItems.addAll(a2.dataItems);
		a3.updateLikelihood();
		
		return a3;
	}
	
	public void addItem(int i) {
		dataItems.add(i);
	}
/*	
	public void updateLikelihood() {
		double[] v = new double[myClusterer.numCols];
		double vt = 0.0;
		int i = 0;
		int idx = 0;
		int j = 0;
		int[][] data = myClusterer.data;
		
		lik = 0.0;
		
		for (idx=0;idx<dataItems.size();idx++) {
			i = dataItems.get(idx);
			for (j=0;j<data[i].length;j++) {
				v[data[i][j]] += 1.0;
				vt += 1.0;
			}
		}
		
		for (j=0;j<v.length;j++) {
			v[j] = v[j]/vt;
		}
		
		for (idx=0;idx<dataItems.size();idx++) {
			i = dataItems.get(idx);
			for (j=0;j<data[i].length;j++) {
				lik += Math.log(v[data[i][j]]);
			}
		}
	} 
 */
	public void updateLikelihood() {
		double[] v = new double[myClusterer.maxValue];
		double vt = 0.0;
		int i = 0;
		int idx = 0;
		int j = 0;
		int k = 0;
		int[][] data = myClusterer.data;
		
		lik = 0.0;
		
		for (j=0;j<data[0].length;j++) {
			vt = 0.0;
			for (k=0;k<myClusterer.maxValue;k++) {
				v[k] = 0.0;
			}
			for (idx=0;idx<dataItems.size();idx++) {
				i = dataItems.get(idx);
				v[myClusterer.data[i][j]]++;
				vt++;
			}
			for (k=0;k<myClusterer.maxValue;k++) {
				v[k] = v[k]/vt;
			}
			
			for (idx=0;idx<dataItems.size();idx++) {
				i = dataItems.get(idx);
				lik += Math.log(v[myClusterer.data[i][j]]);
			}
		}
	} 
}
