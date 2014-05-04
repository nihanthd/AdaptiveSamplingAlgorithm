package com.utils;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import com.algorithm.SamplingAlgortihm;

public class InputReader {

	private static BufferedReader br;
//	static String fileName = "E:\\workspace\\PGM\\SamplingAlgorithm2\\input\\17_4_s.binary.uai";
	

	/**
	 * @author Nihanth
	 * @param sa 
	 * @param filepath 
	 * @param networkType 
	 * @param functionNodesArray 
	 * @param networkNodesArray 
	 * @throws NumberFormatException
	 * @throws IOException
	 */
	public static void readInput(SamplingAlgortihm sa, String fileName) throws NumberFormatException, IOException {
		br = new BufferedReader(new FileReader(fileName));
		BufferedReader br2 = new BufferedReader(new FileReader(fileName));

		String line = null;
		String line2 = null;
		int noOfVariables = 0;
		String domainSizes[] = null;
		NetworkNode[] nodesArray = null;
		int noOfFunctions = 0;
		FunctionNode[] functionsArray = null;
		int noOfEntriesInFunction;
		String network_Type = null;
		//reading network type
		if((line = br.readLine()) != null)
		{
			network_Type = line;
		}

		//reading no of variables
		if((line = br.readLine()) != null)
		{
			noOfVariables = Integer.parseInt(line);
		}
		
		//setting nodes array
		nodesArray = new NetworkNode[noOfVariables];
		
		//reading domain sixes of each variable.
		if((line = br.readLine()) != null)
		{
			domainSizes = line.split("[\\s]+");
		}
		
		//creating each node with the domain size and variable.
		//and setting it to the array of node for tracking
		for(int i = 0; i < noOfVariables; i++)
		{
			NetworkNode node = new NetworkNode(i, domainSizes[i]);
			nodesArray[i] = node;
		}
		
		//reading no of functions
		if((line = br.readLine()) != null)
		{
			noOfFunctions = Integer.parseInt(line);
		}
		
		functionsArray = new FunctionNode[noOfFunctions];
		
		//traversing the buffered reader until the functions.
		while((line2 = br2.readLine())!= null && !line2.isEmpty())
		{
		}
		int i = 0;
		
		//reading functions
		while(((line = br.readLine())!= null && !line.isEmpty()) && (line2 = br2.readLine())!= null)
		{
			String[] lineArray = line.split("[\\s]+");
			while(line2.isEmpty())
			{
				line2 = br2.readLine();
			}
			noOfEntriesInFunction = Integer.parseInt(line2);
			double[] valuesInFunction = getValuesInFunction(br2, noOfEntriesInFunction);
			FunctionNode functionNode = 
					new FunctionNode(nodesArray, lineArray, valuesInFunction);
			functionsArray[i] = functionNode;
			i++;
		}
		sa.setNetworkType(network_Type);
		sa.setNetworkNodesArray(nodesArray);
		sa.setFunctionNodesArray(functionsArray);
	}
	
	
	/**
	 * Reads all the values of the function and returns an array of the values.
	 * @param br2
	 * @param noOfEntriesInFunction 
	 * @return
	 * @throws IOException
	 */
	private static double[] getValuesInFunction(BufferedReader br2, int noOfEntriesInFunction) throws IOException {
		double[] values = new double[noOfEntriesInFunction];
		String line = null;
		String[] lineArray = null;
		int i = 0;
		while(i < noOfEntriesInFunction && (line = br2.readLine()) != null && !line.isEmpty())
		{
			line = line.trim();
			lineArray = line.split("[\\s]+");
			for(int j = 0; i < noOfEntriesInFunction && j < lineArray.length; j++,i++)
			{
				values[i] = Double.parseDouble(lineArray[j]);
			}
		}
		return values;
	}


	public static void readEvidence(SamplingAlgortihm sa, String fileName) throws IOException {
		String evidenceFileName = fileName + ".evid";
		br = new BufferedReader(new FileReader(evidenceFileName));
		String line = null;
		int noOfEvidenceVariables = 0;
		if((line = br.readLine()) != null)
		{
			line = line.trim();
			noOfEvidenceVariables = Integer.parseInt(line);
		}
		int[] evidenceVariables = new int[noOfEvidenceVariables];
		int[] evidenceVariablesValues = new int[noOfEvidenceVariables];
		int i = 0;
		while(noOfEvidenceVariables > 0 && (line = br.readLine()) != null)
		{
			line = line.trim();
			String[] lineArray = line.split("[\\s]+");
			evidenceVariables[i] = Integer.parseInt(lineArray[0]);
			evidenceVariablesValues[i] = Integer.parseInt(lineArray[1]);
			i++;
		}
		sa.setEvidenceVariables(evidenceVariables);
		sa.setEvidenceVariablesValues(evidenceVariablesValues);
	}
}
