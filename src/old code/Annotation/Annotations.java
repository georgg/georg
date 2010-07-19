package edu.mit.csail.psrg.georg.Annotation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;

import edu.mit.csail.psrg.georg.DataAccess.ClusterReader;

public class Annotations {

	LinkedHashMap<String,OntTerm> DAG = new LinkedHashMap<String,OntTerm>();
	HashMap<String,OntGene> genes = new HashMap<String,OntGene>();
	
	public void addGene(OntGene gene) {
		genes.put(gene.ID,gene);
	}
	
	public boolean checkGene(OntGene gene) {
		return genes.containsKey(gene.ID);
	}
	
	public OntGene getGene(String name) {
		return genes.get(name);
	}
	
	public void addTerm(OntTerm term) {
		DAG.put(term.ID,term);
	}
	
	public OntTerm getTerm(String ID) {
		return DAG.get(ID);
	}
	
	public void associateParents() {
		OntTerm[] mv = new OntTerm[DAG.size()];
		DAG.values().toArray(mv);
		associateParentsUtil(mv);
	}
	
	void associateParentsUtil(OntTerm[] mv) {
		int i = 0;
		for(i=0;i<mv.length;i++) {
			mv[i].associateParents(this);
		}
	}
	
	void topologicalSort() throws DAGException {
		ArrayList<OntTerm> list = new ArrayList<OntTerm>();
		Iterator<OntTerm> iter = DAG.values().iterator();
		OntTerm term = null;
		while(iter.hasNext()) {
			term = iter.next();
			if (term.color == OntTerm.White) {
				DFSVisit(term,list);
			}
		}
		DAG.clear();
		iter = list.iterator();
		while(iter.hasNext()) {
			term = iter.next();
			addTerm(term);
		//	System.out.println(term.ID);
		}
	}
	
	void DFSVisit(OntTerm term,ArrayList<OntTerm> list) throws DAGException {
		if (term.color != OntTerm.White) {
			System.out.println("DAG error detected in " + term.ID);
			throw new DAGException();
		}
		OntTerm child = null;
		term.color = OntTerm.Gray;
		if (term.children != null) {
			Iterator<OntTerm> iter = term.children.iterator();
			while(iter.hasNext()) {
				child = iter.next();
				if (child.color == OntTerm.White) {
					DFSVisit(child,list);
				}
			}
		}
		term.color = OntTerm.Black;
		list.add(0,term);
	}
	
	void propagateGeneCounts() {
		OntTerm[] list = new OntTerm[DAG.size()];
		DAG.values().toArray(list);
		int i = 0;
		int j = 0;
		
		for (i=list.length-1;i>=0;i--) {
			if (list[i].genes != null) {
				list[i].numGenes = list[i].numGenes + list[i].genes.size();
			}
			if (list[i].parents != null) {
				for (j=0;j<list[i].parents.length;j++) {
					list[i].parents[j].numGenes = list[i].parents[j].numGenes + list[i].numGenes;
				}
			}
		}
 	}
	
	public void readOBOFile(String OBOFileName) throws IOException, DAGException {
		BufferedReader is = new BufferedReader(new FileReader(OBOFileName));
		String line = "";

		line = is.readLine();
		OntTerm term = null;
		String tag = "";
		String content = "";
		int fi = 0;
		int terms = 0;
		
		while(line != null) {
			if (line.equals("[Term]")) {
				if (term != null) {
					addTerm(term);
				}
				term = new OntTerm();
				terms++;
			} else {
				if (terms > 0) {
					fi = line.indexOf(":");
					if (fi >= 0) {
						tag = line.substring(0,fi);
						content = line.substring(fi+2);
					}
					term.setField(tag,content);
				}
			}
			line = is.readLine();
		}
		
		is.close();
	}
	
	public void readFiles(String geneAssocFileName,String OBOFileName,HashSet<String> useGenes) throws IOException, DAGException {
		ClusterReader reader = new ClusterReader();
		
		readOBOFile(OBOFileName);
		reader.readFile(geneAssocFileName);
		
		// the first row contains column names
		// the first four columns must at least be gene ID, gene name, gene description, and GO
		// other annotation types can be specified optionally
		
		int i = 0;
		int j = 0;
		OntGene gene = null;
		boolean useGene = true;
		int numOtherTypes = 0;
		
		numOtherTypes = reader.clusters.get(0).size() - 4;
		String temp = null;
		
		for (i=1;i<reader.clusters.size();i++) {
			gene = new OntGene();
			
			gene.ID = reader.clusters.get(i).get(0).toUpperCase();
			gene.Name = reader.clusters.get(i).get(1);
			gene.Description = reader.clusters.get(i).get(2);
			
			// see if the gene is in the filter list
			useGene = true;
			if (useGenes != null) {
				useGene = useGenes.contains(gene.ID);
			}
			
			if (!checkGene(gene) && useGene) {
				addGene(gene);
			}
			
			if (useGene) {
				// process GO terms
				if (reader.clusters.get(i).get(3).length() > 1) {
					addToGO(gene,reader.clusters.get(i).get(3));
				}
				
				if (numOtherTypes > 0) {
					for (j=0;j<numOtherTypes;j++) {
						temp = reader.clusters.get(i).get(j+4);
						if (temp.length() > 1) {
							addToOtherType(gene,temp,reader.clusters.get(0).get(j+4));
						}
					}
				}
			}
		}
		associateParents();
		topologicalSort();
		propagateGeneCounts();
	}
	
	void addToOtherType(OntGene gene,String cats,String type) {
	//	String[] terms = cats.split(",");
		String[] terms = cats.split(" ");
		int i = 0;
		OntTerm parent = null;
		for (i=0;i<terms.length;i++) {
			parent = getTerm(terms[i]);
			
			if (parent == null) {
				parent = new OntTerm();
				parent.ID = terms[i];
				parent.name = terms[i];
				parent.def = terms[i];
				parent.namespace = type;
				addTerm(parent);
			}
			
			parent.addOntTerm(gene);
		}
	}
	
	void addToGO(OntGene gene,String cats) {
	//	String[] terms = cats.split(",");
		String[] terms = cats.split(" ");
		int i = 0;
		OntTerm parent = null;
		for (i=0;i<terms.length;i++) {
			parent = getTerm(terms[i]);
			if (parent != null) {
				parent.addOntTerm(gene);
			}
		}
	}
}
