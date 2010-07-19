package edu.mit.csail.psrg.georg.HDPTopicModel;

import java.awt.Color;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Vector;

import edu.mit.csail.psrg.georg.DataAccess.MicroArrayData;
import edu.mit.csail.psrg.georg.Graph.TabText;
import edu.mit.csail.psrg.georg.Graph.TabularComponent;
import edu.mit.csail.psrg.georg.HierarchicalCluster.HierarchicalCluster;
import edu.mit.csail.psrg.georg.StatUtil.VectorUtil;

public class OutputStrongModulesTabular {
	HierDirichletProcess myHDP = null;
	ArrayList<StrongModule> modules = null;
	String fileName = "";
	TabularComponent sheet = null;
	double[][] data = null;
	String matrixFileName = "";
	
	public OutputStrongModulesTabular(HierDirichletProcess dp,ArrayList<StrongModule> m,String fn,String mfn) {
		myHDP = dp;
		modules = m ;
		fileName = fn;
		matrixFileName = mfn;
		
		initDataMatrix();
		
		try {
			sheet.outputSVG(fileName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	// determine if the module is exclusively mouse, human, or mixed
	int moduleType(StrongModule module) {
		boolean hasMouse = false;
		boolean hasHuman = false;
		Document doc = null;
		Iterator<Document> docIter = module.controllingDocuments.iterator();
		while (docIter.hasNext()) {
			doc = docIter.next();
			if (doc.annotation.indexOf("h_") == 0) {
				hasHuman = true;
			}
			if (doc.annotation.indexOf("m_") == 0) {
				hasMouse = true;
			}
		}
		
		if (hasHuman & !hasMouse)
			return 1;
		
		if (!hasHuman & hasMouse)
			return 2;
		
		if (!hasHuman & !hasMouse)
			return -1;
		
		// has both mouse and human
		return 3;
	}
	
	void initDataMatrix() {
		int i = 0;
		int j = 0;
		int k = 0;
		DirichletProcess DP = null;
		LinkedHashMap<Document,Integer> documents = new LinkedHashMap<Document,Integer>();
		HashMap<DirichletProcess,Integer> processes = new HashMap<DirichletProcess,Integer>();
		int[] moduleCode = new int[modules.size()];
		int numDocs = 0;
		Document doc = null;
		StrongModule module = null;
		Iterator<Document> docIter = null;
		for (i=0;i<myHDP.DP.length;i++) {
			DP = myHDP.DP[i];
			processes.put(DP,i);
			if (DP.documents != null) {
				for (j=0;j<DP.documents.length;j++) {
					doc = DP.documents[j];
					documents.put(doc,numDocs);
					numDocs++;
				}
			}
		}
		
		data = new double[documents.size()][modules.size()];
		int relative = 0;
		boolean crossesGroup = false;
		for (j=0;j<modules.size();j++) {
			module = modules.get(j);
			docIter = module.controllingDocuments.iterator();
			relative = closestRelative(module,processes);
			crossesGroup = false;
			while(docIter.hasNext()) {
				doc = docIter.next();
				i = documents.get(doc);
				data[i][j] = ((double) module.documentLoad.get(doc))/((double) module.numOccur);
				if (relative != processes.get(doc.parent.parent)) {
					crossesGroup = true;
				}
			}
			if (crossesGroup) {
				// the closest relative is the root
				if (relative == 0) {
					moduleCode[j] = 2;
				} else {
					moduleCode[j] = 1;
				}
			}
		}
		
		// this species splitting will override everything else - get rid of
		// this to use more of the tree information
		int mt = 0;
		for (j=0;j<modules.size();j++) {
			mt = moduleType(modules.get(j));
			if (mt == 3)
				moduleCode[j] = 2;
			if (mt == 1)
				moduleCode[j] = 1;
			if (mt == 2)
				moduleCode[j] = 0;
		}
		
		// normalize rows of the data matrix (loading on each tissue)
		double norm = 0.0;
		for (i=0;i<data.length;i++) {
			norm = 0.0;
			for (j=0;j<data[0].length;j++) {
				norm += data[i][j];
			}
			if (norm > 0.0) {
				for (j=0;j<data[0].length;j++) {
					data[i][j] = data[i][j]/norm;
				}
			}
		}
		
		int headerRows = 2;
		sheet = new TabularComponent(documents.size()+headerRows*2,modules.size()+2,30);
		int[] colWidths = new int[modules.size()+2];
		for (i=0;i<colWidths.length;i++) {
			colWidths[i] = 20;
		}
		colWidths[1] = 150;
		sheet.setColWidths(colWidths);
		sheet.setGrid(true);
		
		double maxLoad = 0.0;
		double minLoad = 1000.0;
		for (i=0;i<data.length;i++) {
			for (j=0;j<data[0].length;j++) {
				if (data[i][j] > maxLoad) {
					maxLoad = data[i][j];
				}
				if (data[i][j] < minLoad & data[i][j] > 0.0) {
					minLoad = data[i][j];
				}
			}
		}
		
		ArrayList<Integer> colOrder = orderDataColumns();
		ArrayList<Integer> rowPosition = orderDataRows();
		int rp = 0;
		
		docIter = documents.keySet().iterator();
		TabText myTab = null;
		float color = 0.0f;
		double minIntensity = 0.25;
		double maxIntensity = 1.0;
		double v = 0.0;
		double a = (maxIntensity-minIntensity)/(maxLoad-minLoad);
		double b = minIntensity - minLoad*a;
		
		int[] groupIDs = new int[documents.keySet().size()];
		String[] tissueNames = new String[documents.keySet().size()];
		
		while(docIter.hasNext()) {
			doc = docIter.next();
			i = documents.get(doc);
		//	rp = rowPosition.get(i);
			rp = i;
			myTab = new TabText(processHDPLabel(doc.annotation));
			tissueNames[i] = processHDPLabel(doc.annotation);
			sheet.setCell(rp+headerRows,1,myTab);
			j = processes.get(doc.parent.parent);
			groupIDs[i] = j-1;
			tissueNames[i] = tissueNames[i] + "\t" + (new Integer(groupIDs[i])).toString();
			// assumes there are 3 processes in front (root, human, mouse)
			myTab = new TabText((new Integer(j-2)).toString());
			myTab.setFill(true);
			myTab.setFillColor(colorCycle(j-2));
			sheet.setCell(rp+headerRows,0,myTab);
			for (j=0;j<data[0].length;j++) {
				k = colOrder.get(j);
				if (data[i][k] > 0.0) {
					v = data[i][k]*a + b;
					color = (float) v;
					if (color > 1.0f) {
						color = 1.0f;
					}
					if (color < 0.0f) {
						color = 0.0f;
					}
					myTab = new TabText("");
					myTab.setFillColor(1.0f-color,1.0f-color,1.0f);
					myTab.setFill(true);
					sheet.setCell(rp+headerRows,j+2,myTab);
				}
			}
		}
		
		for (j=2;j<modules.size()+2;j++) {
			k = colOrder.get(j-2);
			myTab = new TabText((new Integer(k+1)).toString());
			myTab.setTextAlign(TabText.CENTER);
			sheet.setCell(headerRows-1,j,myTab);
			sheet.setCell(sheet.getNumRows()-2,j,myTab);
			if (moduleCode[k] > 0) {
				if (moduleCode[k] == 2) {
					myTab = new TabText("");
					myTab.setFill(true);
					myTab.setFillColor(1.0f,0.0f,0.0f);
				} else {
					myTab = new TabText("");
					myTab.setFill(true);
					myTab.setFillColor(0.0f,0.5f,0.0f);
				}
				sheet.setCell(headerRows-2,j,myTab);
				sheet.setCell(sheet.getNumRows()-1,j,myTab);
			}
		}
		
		MicroArrayData outputData = new MicroArrayData();
		outputData.values = data;
		outputData.numCols = data[0].length;
		outputData.numRows = data.length;
		outputData.geneNames = tissueNames;
		
		try {
			outputData.writeFile(matrixFileName);
		} catch(IOException e) {
			System.out.println(e);
		}
	}
	
	String processHDPLabel(String label1) {
		String label2 = null;
		int f = label1.indexOf("human_");
		if (f == 0) {
			label2 = label1.substring(f+6);
			label2 = "h_" + label2;
		} else {
			f = label1.indexOf("mouse_");
			if (f == 0) {
				label2 = label1.substring(f+6);
				label2 = "m_" + label2;
			} else {
				label2 = label1;
			}
		}
		return label2;
	}
	
	int closestRelative(StrongModule module,HashMap<DirichletProcess,Integer> DPMap) {
		int relative = -1;
		
		Iterator<Document> iter = module.controllingDocuments.iterator();
		Document doc = null;
		HashSet<Integer> relatives = new HashSet<Integer>();
		HashSet<Integer> myAncestors = null;
		
		while(iter.hasNext()) {
			doc = iter.next();
			myAncestors = ancestors(doc,DPMap);
			if (relatives.isEmpty()) {
				relatives.addAll(myAncestors);
			} else {
				relatives.retainAll(myAncestors);
			}
		}
		
		Iterator<Integer> relativeIter = relatives.iterator();
		int curRelative = 0;
		while (relativeIter.hasNext()) {
			curRelative = relativeIter.next();
			if (curRelative > relative) {
				relative = curRelative;
			}
		}
		
		return relative;
	}
	
	HashSet<Integer> ancestors(Document doc,HashMap<DirichletProcess,Integer> DPMap) {
		HashSet<Integer> myAncestors = new HashSet<Integer>();
		
		DirichletProcess parent = null;
		parent = doc.parent.parent;
		
		while(parent != null) {
			myAncestors.add(DPMap.get(parent));
			parent = parent.parent;
		}
		
		return myAncestors;
	}
	
	public ArrayList<Integer> orderDataRows() {
		int nr = data.length;
		double[][] distances = new double[nr-1][];
		int i = 0;
		int j = 0;
		int k = 0;
		Vector<Double> x = new Vector<Double>();
		Vector<Double> y = new Vector<Double>();
		double v = 0.0;
		double maxDistance = 0.0;
		
		for (i=0;i<nr-1;i++) {
			x.clear();
			for (k=0;k<data[i].length;k++) {
				x.add(data[i][k]);
			}
			distances[i] = new double[nr-i-1];
			for (j=(i+1);j<nr;j++) {
				y.clear();
				for (k=0;k<data[j].length;k++) {
					y.add(data[j][k]);
				}
				v = (1.0-VectorUtil.correlation(x,y))/2.0;
			//	v = VectorUtil.euclidean(x,y);
				distances[i][j-i-1] = v;
				if (v > maxDistance) {
					maxDistance = v;
				}
			}
		}
		HierarchicalCluster hierCluster = new HierarchicalCluster();
		hierCluster.cluster(distances,maxDistance + 1.0);
		
		ArrayList<Integer> rowPosition = new ArrayList<Integer>();
		for (i=0;i<nr;i++) {
			rowPosition.add(0);
		}
		for (i=0;i<hierCluster.clusters.get(0).size();i++) {
			rowPosition.set(hierCluster.clusters.get(0).get(i),i);
		}
		
		return rowPosition;
	}
	
	public ArrayList<Integer> orderDataColumns() {
		int nc = data[0].length;
		double[][] distances = new double[nc-1][];
		int i = 0;
		int j = 0;
		int k = 0;
		Vector<Double> x = new Vector<Double>();
		Vector<Double> y = new Vector<Double>();
		double v = 0.0;
		double maxDistance = 0.0;
		
		for (i=0;i<nc-1;i++) {
			x.clear();
			for (k=0;k<data.length;k++) {
				x.add(data[k][i]);
			}
			distances[i] = new double[nc-i-1];
			for (j=(i+1);j<nc;j++) {
				y.clear();
				for (k=0;k<data.length;k++) {
					y.add(data[k][j]);
				}
				v = (1.0-VectorUtil.correlation(x,y))/2.0;
				distances[i][j-i-1] = v;
				if (v > maxDistance) {
					maxDistance = v;
				}
			}
		}
		HierarchicalCluster hierCluster = new HierarchicalCluster();
		hierCluster.cluster(distances,1.5);
		
		return hierCluster.clusters.get(0);
	}
	
	public Color colorCycle(int c) {
		int r = (int) (12.0*Math.floor( ((double) c)/12.0 ));
		r = c - r;
		
		switch(r) {
		case 0:
			return Color.WHITE;
		case 1:
			return Color.BLUE;
		case 2:
			return Color.CYAN;
		case 3:
			return Color.DARK_GRAY;
		case 4:
			return Color.GREEN;
		case 5:
			return Color.GRAY;
		case 6:
			return Color.MAGENTA;
		case 7:
			return Color.ORANGE;
		case 8:
			return Color.LIGHT_GRAY;
		case 9:
			return Color.PINK;
		case 10:
			return Color.YELLOW;
		case 11:
			return Color.RED;
		}
		return Color.BLACK;
	}
}
