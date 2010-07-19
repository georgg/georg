package edu.mit.csail.psrg.georg.Annotation;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

public class OntTerm {
	ArrayList<String> parentNames = null;
	OntTerm[] parents = null;
	HashSet<OntGene> genes = null;
	HashSet<OntTerm> children = null;
	public String ID;
	public String name;
	public String def;
	String namespace = "AnyCategory";
	int numGenes = 0;
	int color = 0;
	static int White = 0;
	static int Gray = 1;
	static int Black = 2;
	
	public void addParentName(String pname) {
		if (parentNames == null) {
			parentNames = new ArrayList<String>();
		}
		parentNames.add(pname);
	}
	
	public void associateParents(Annotations DAG) {
		if (parentNames == null) { return; }
		
		int i = 0;
		OntTerm parent = null;
		parents = new OntTerm[parentNames.size()];
		for (i=0;i<parentNames.size();i++) {
			parent = DAG.getTerm(parentNames.get(i));
			if (parent != null) {
				parents[i] = parent;
				parent.addChild(this);
			}
		}
		parentNames = null;
	}
	
	void addChild(OntTerm g) {
		if (children == null) {
			children = new HashSet<OntTerm>();
		}
		children.add(g);
	}
	
	public void setField(String tag,String content) {
		if (tag.contentEquals("id")) {
			ID = content;
		}
		
		if (tag.contentEquals("name")) {
			name = content;
		}
		
		if (tag.contentEquals("is_a") || tag.contentEquals("part_of")) {
			int fi = content.indexOf(" ");
			String pID;
			if (fi >= 0) {
				pID = content.substring(0,fi);
			} else{
				pID = content;
			}
			addParentName(pID);
		}
		
		if (tag.contentEquals("def")) {
			def = content;
		}
		
		if (tag.contentEquals("namespace")) {
			namespace = content;
		/*	if (content.contentEquals("biological_process")) {
				namespace = "BiologicalProcess";
			}
			if (content.contentEquals("molecular_function")) {
				namespace = "MolecularFunction";
			}
			if (content.contentEquals("cellular_component")) {
				namespace = "CellularComponent";
			} */
		}
	}
	
	public void addGene(OntGene gene) {
		if (genes == null) {
			genes = new HashSet<OntGene>();
		}
		genes.add(gene);
	}
	
	public void addGenes(HashSet<OntGene> g2) {
		if (genes == null) {
			genes = new HashSet<OntGene>();
		}
		genes.addAll(g2);
	}
	
	public void getAllParents(HashSet<OntTerm> list,int minGenes,int maxGenes,String useNamespace) {
		if (numGenes >= minGenes && numGenes <= maxGenes && (namespace.equals(useNamespace) || useNamespace.equals("AnyCategory"))) {
			list.add(this);
		}
		if (parents == null || numGenes > maxGenes || (!namespace.equals(useNamespace) && !useNamespace.equals("AnyCategory"))) {
			return;
		}
		
		int i = 0;
		for (i=0;i<parents.length;i++) {
			parents[i].getAllParents(list,minGenes,maxGenes,useNamespace);
		}
	}

	public void addOntTerm(OntGene gene) {
		if (gene.terms == null) {
			gene.terms = new HashSet<OntTerm>();
		}
		gene.terms.add(this);
		addGene(gene);
	}
}
