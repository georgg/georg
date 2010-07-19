package edu.mit.csail.psrg.georg.DataAccess;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.StringTokenizer;

public class CommandLineParser {
	LinkedHashMap<String,String> typeMap = new LinkedHashMap<String,String>();
	LinkedHashMap<String,String> descriptionMap = new LinkedHashMap<String,String>();
	
	private class NoArgument {	
	}
	
	public void displayDescriptions() {
		String s1;
		String s2;
		Iterator<String> iter = typeMap.keySet().iterator();
		
		while(iter.hasNext()) {
			s1 = iter.next();
			s2 = s1 + " " + descriptionMap.get(s1);
			System.out.println(s2);
		}
	}
	
	public HashMap<String,Object> parseCommandFile(String fName)throws CommandLineParserException,IOException {
		HashMap<String,Object> parseMap = new HashMap<String,Object>();
		String type;
		Object o = null;
		
		BufferedReader is = new BufferedReader(new FileReader(fName));
		String[] args = null;
		String line = "";
		line = is.readLine();
		while(line != null) {
			args = line.split("\t");
			
			if (typeMap.containsKey(args[0])) {
				type = typeMap.get(args[0]);
				if (!type.equals("NoArgument")) {
					o = getType(type,args[1]);
				} else {
					o = getType(type,null);
				}
				parseMap.put(args[0],o);
			} else {
				throw new CommandLineParserException(args[0] + " not a valid switch");
			}
			line = is.readLine();
		}
		return parseMap;
	}
	
	public HashMap<String,Object> parseCommandLine(String[] args) throws CommandLineParserException {
		HashMap<String,Object> parseMap = new HashMap<String,Object>();
		String arg;
		String type;
		Object o = null;
		
		if (args.length == 0) {
			return parseMap;
		}
		
		int k = 0;
		while(k < args.length) {
			arg = args[k];
			if (typeMap.containsKey(arg)) {
				type = typeMap.get(arg);
				if (!type.equals("NoArgument")) {
					k++;
					o = getType(type,args[k]);
				} else {
					o = getType(type,null);
				}
				parseMap.put(arg,o);
			} else {
				throw new CommandLineParserException(arg + " not a valid switch");
			}
			k++;
		}
		return parseMap;
	}
	
	public void addSwitch(String switchName,String type,String description) throws CommandLineParserException {
		Object o = getType(type,null);
		if (o == null) {
			throw new CommandLineParserException(type);
		}
		typeMap.put(switchName,type);
		descriptionMap.put(switchName,description);
	}
	
	public Object getType(String type,String contents) {
		Object o = null;
		
		if (type.equals("String")) {
			if (contents != null) {
				o = new String(contents);
			} else {
				o = new String();
			}
		}
		if (type.equals("Integer")) {
			if (contents != null) {
				o = new Integer(contents);
			} else {
				o = new Integer(0);
			}
		}
		if (type.equals("Boolean")) {
			if (contents != null) {
				o = new Boolean(contents);
			} else {
				o = new Boolean(false);
			}
		}
		if (type.equals("Double")) {
			if (contents != null) {
				o = new Double(contents);
			} else {
				o = new Double(0.0);
			}
		}
		if (type.equals("NoArgument")) {
			o = new NoArgument();
		}
		
		return o;
	}
	
	public void outputMap(HashMap<String,Object> map) {
		Iterator<String> iter = map.keySet().iterator();
		String s = null;
		Object o = null;
		while(iter.hasNext()) {
			s = iter.next();
			o = map.get(s);
			System.out.print(s);
			System.out.print(" ");
			System.out.println(o);
		}
	}
}
