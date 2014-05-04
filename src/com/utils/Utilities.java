package com.utils;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

import com.algorithm.SamplingAlgortihm;

public class Utilities {
	
	public static HashSet<Integer> getHashSetForVariablesInFunctions(
			FunctionNode node1, FunctionNode node2) {
		HashSet<Integer> hashSetOfCommonVariables = new HashSet<Integer>();
		for(int i = 0; i <= node1.nodesValuesArray.length-1; i++)
		{
			hashSetOfCommonVariables.add(node1.nodesValuesArray[i]);
		}
		for(int i = 0; i <= node2.nodesValuesArray.length-1; i++)
		{
			hashSetOfCommonVariables.add(node2.nodesValuesArray[i]);
		}
		return hashSetOfCommonVariables;
	}
	
	/**
	 * This method returns true if the nodeElement is present in the function array
	 * @param nodeValue
	 * @param 
	 * @return
	 */
	public static boolean isVariablePresentInFunction(int nodeValue,
			int[] nodesValuesArray) {
		int length = nodesValuesArray.length;
                for(int i = 0; i < length; i++)
                {
                        if(nodeValue == nodesValuesArray[i])
                        {
                                return true;
                        }
                }
		return false;
	}
	
	/**
	 * This method adds the elements in nodesValuesArray to hashset.
	 * @param adjacencySet 
	 * @param nodesValuesArray
	 */
	public static void addVariablesIntoAdjacencySet(HashSet<Integer> adjacencySet, int[] nodesValuesArray) {
		int length = nodesValuesArray.length;
		for(int i = 0; i < length; i++)
		{
			adjacencySet.add(nodesValuesArray[i]);
		}
	}
	
	/**
	 * This method updates the adjacency nodes.
	 * @param node
	 * @param networkNodesArray
	 */
	public static void updateAdjacentNodes(NetworkNode node,
			NetworkNode[] networkNodesArray) {
		for(Integer k: node.adjacentNodes)
		{
			networkNodesArray[k].adjacentNodes.remove(node.nodeValue);
			for(Integer j : node.adjacentNodes)
			{
				if(j.intValue() != k.intValue())
				{
					networkNodesArray[k].adjacentNodes.add(networkNodesArray[j].nodeValue);
				}
			}
		}
	}
	
	/**
	 * This method returns the new array of variables in the function
	 * @param variableToBeSummed
	 * @param nodesValuesArray
	 * @return
	 */
	public static int[] getNewVariablesArray(int variableToBeSummed,
			int[] nodesValuesArray) {
		int[] newNodeValuesArray = new int[nodesValuesArray.length -1];
		for(int i = 0, j = 0; i < nodesValuesArray.length; i++)
		{
			if(nodesValuesArray[i] != variableToBeSummed)
			{
				newNodeValuesArray[j] = nodesValuesArray[i];
				j++;
			}
		}
		return newNodeValuesArray;
	}
	
	
	/**
	 * 
	 * @param nodeVariableStrideMap
	 * @param node
	 * @return 
	 */
	public static HashMap<Integer, Integer> constructVariableStrideMap(FunctionNode node) {
		HashMap<Integer, Integer> variableStrideMap = new HashMap<Integer, Integer>();
		for(int i = 0; i <= node.nodesValuesArray.length - 1; i++)
		{
			variableStrideMap.put(node.nodesValuesArray[i], node.stepsForVariables[i]);
		}
		return variableStrideMap;
	}
	
	/**
	 * This method returns the stride of the variables in a function/factor.
	 * @param node
	 * @param variable
	 * @return
	 */
	public static int getStrideForVariableInFunction(FunctionNode node,
			int variable) {
		for(int i = 0; i < node.nodesValuesArray.length; i++)
		{
			if(node.nodesValuesArray[i] == variable)
			{
				return node.stepsForVariables[i];
			}
		}
		return 0;
	}

	public static LinkedList<FunctionNode> removeElementsFromList(
			LinkedList<FunctionNode> removalElements,
			LinkedList<FunctionNode> temporaryFunctions) {
//		System.out.println("Before - " + temporaryFunctions.size());
		for(FunctionNode removalNode: removalElements)
		{
			temporaryFunctions.remove(removalNode);
		}
//		System.out.println("After - " + temporaryFunctions.size());
		return temporaryFunctions;
	}

	public static int getStrideForPreviousVariableInFunction(
			FunctionNode node, int variable) {
		for(int i = 0; i < node.nodesValuesArray.length; i++)
		{
			if(node.nodesValuesArray[i] == variable)
			{
				if(i == 0)
				{
					return 0; 
				}
				else
				{
					return node.stepsForVariables[i-1];
				}
			}
			
		}
		return 0;
	}

	public static int getDomainSizeOfPreviousVariableInFunction(
			FunctionNode node, int variable, NetworkNode[] networkNodesArray) {
		
		for(int i = 0; i < node.nodesValuesArray.length; i++)
		{
			if(node.nodesValuesArray[i] == variable)
			{
				if(i == 0)
				{
					return 0; 
				}
				else
				{
					return networkNodesArray[node.nodesValuesArray[i-1]].domainSize;
				}
			}
			
		}
		return 0;
	}
	
	/*
	 * This method returns true if variable is present in the evidence variables or else false.
	 * @param networkNode
	 * @param evidenceVariables
	 * @return
	 */
	public static boolean isVariablePresentInEvidenceVariables(
			int variable, int[] evidenceVariables) {
		for(int i = 0; i < evidenceVariables.length; i++)
		{
			if(variable == evidenceVariables[i])
			{
				return true;
			}
		}
		return false;
	}

	/**
	 * This method builds adjacency list for variables in a function.
	 * @param functionNode
	 * @param networkNodesArray
	 * @param evidenceVariables
	 */
	public static void buildAdjacencyListForEachVariableInFunction(
			FunctionNode functionNode, NetworkNode[] networkNodesArray,
			int[] evidenceVariables) {
		int[] nodesValuesArray = functionNode.nodesValuesArray;
		for(int i = 0; i < nodesValuesArray.length; i++)
		{
			if(!isVariablePresentInEvidenceVariables(nodesValuesArray[i], evidenceVariables))
			{
				for(int j = 0; j < nodesValuesArray.length; j++)
				{
					if(i != j && !isVariablePresentInEvidenceVariables(nodesValuesArray[j], evidenceVariables))
					{
						networkNodesArray[nodesValuesArray[i]].adjacentNodes.add(nodesValuesArray[j]);
					}
				}
			}
		}
	}

	public static void printElementsInIntegerArray(int[] integerArray) {
		for(int i = 0; i < integerArray.length; i++)
		{
			System.out.print(integerArray[i] + ", ");
		}
		System.out.println("");
	}
	
	/**
	 * This method sets the variables of a bucket from nodesValuesArray into VariablesSet 
	 * @param variablesSet
	 * @param nodesValuesArray
	 */
	public static void setVariablesOfBucket(HashSet<Integer> variablesSet,
			int[] nodesValuesArray) {
		for(int variable : nodesValuesArray)
		{
			variablesSet.add(variable);
		}
	}

	public static void setVariablesOfBucket(HashSet<Integer> variablesSet,
			HashSet<Integer> tempVariablesSet) {
		for(int variable : tempVariablesSet)
		{
			variablesSet.add(variable);
		}
		
	}

	public static HashSet<Integer> getBucketWithClusterSizeGreaterThan(int w,
			LinkedList<HashSet<Integer>> treeDecomposition) {
		HashSet<Integer> bucket = null;
		for(HashSet<Integer> temp : treeDecomposition)
		{
			if(w+1 < temp.size())
			{
				bucket = temp;
			}
		}
		return bucket;
	}
	
	/**
	 * This method returns the variable with max count in the bucket. 
	 */
	public static int getVariableFromBucketWithMaximumCountInClusters(
			HashSet<Integer> bucket, int[] variablesCountInClusters) {
		int variable = 0, count = 0;
		for(Integer temp : bucket)
		{
			if(count < variablesCountInClusters[temp.intValue()])
			{
				count = variablesCountInClusters[temp.intValue()];
				variable = temp.intValue();
			}
		}
		return variable;
	}

	public static void generateCopyOfNetworkNodes(SamplingAlgortihm sa) {
		sa.networkNodesArrayCopy = new NetworkNode[sa.networkNodesArray.length];
		for(int i = 0; i < sa.networkNodesArray.length; i++)
		{
			sa.networkNodesArrayCopy[i] = new NetworkNode(sa.networkNodesArray[i]);
		}
	}

	public static void generateCopyOfFunctionNodes(SamplingAlgortihm sa) {
		sa.functionNodesArrayCopy = new FunctionNode[sa.functionNodesArray.length];
		for(int i = 0; i < sa.functionNodesArray.length; i++)
		{
			sa.functionNodesArrayCopy[i] = new FunctionNode(sa.functionNodesArray[i], sa.networkNodesArray);
		}
	}

	public static void generateCopyOfEvidence(SamplingAlgortihm sa) {
		sa.evidenceVariablesCopy = new int[sa.evidenceVariables.length];
		sa.evidenceVariablesValuesCopy = new int[sa.evidenceVariablesValues.length];
		System.arraycopy(sa.evidenceVariables, 0, sa.evidenceVariablesCopy, 0, sa.evidenceVariables.length);
		System.arraycopy(sa.evidenceVariablesValues, 0, sa.evidenceVariablesValuesCopy, 0, sa.evidenceVariablesValues.length);
	}

	public static double sumUpAdaptiveValues(double[] adaptiveValueNew) {
		double value = 0;
		for(int j = 0; j < adaptiveValueNew.length;j++)
		{
        	value += adaptiveValueNew[j]; 
		}
		if(value > 0)
		{
			return value;
		}
		else
		{
			return 1;
		}
	}
	
	/**
	 * This method resets the Q arrays in network nodes and network nodes copy
	 * @param sa
	 */
	public static void resetQInNetworkNodes(SamplingAlgortihm sa) {
		for(int i = 0; i < sa.networkNodesArray.length; i++)
		{
			NetworkNode node = sa.networkNodesArray[i];
			NetworkNode nodeCopy = sa.networkNodesArrayCopy[i];
			node.adaptiveValueNew = new double[node.domainSize];
			node.adaptiveValueOld = new double[node.domainSize];
			nodeCopy.adaptiveValueNew = new double[nodeCopy.domainSize];
			nodeCopy.adaptiveValueOld = new double[nodeCopy.domainSize];
		}
	}

	/**
	 * This method returns the assignment from domain size for a variable when we know its value.
	 * @param networkNode
	 * @param assignmentValue
	 * @return
	 */
	public static int getValueForAVariableWithAssignmentValue(
			NetworkNode networkNode, double assignmentValue) {
		int assignment = 0;
		for(int i = 0; i < networkNode.assignmentRange.length; i++)
		{
			if(assignmentValue < networkNode.assignmentRange[i])
			{
				assignment = i;
				break;
			}
		}
		return assignment;
	}
	
	/**
	 * This method updates the adaptive old values with the new value after computing them
	 * and also calculates the domain ranges.
	 * @param adaptiveValueOld
	 * @param adaptiveValueNew
	 * @param temp
	 * @param valueSum 
	 */
	public static void updateAdaptiveOldQ(double[] adaptiveValueOld,
			double[] adaptiveValueNew, NetworkNode temp, double valueSum) {
		double range = 0.0;
		temp.assignmentRange = new double[adaptiveValueNew.length];
		for(int i = 0; i < adaptiveValueNew.length; i++)
		{
			adaptiveValueOld[i] = (double)adaptiveValueNew[i]/valueSum;
			range += adaptiveValueOld[i];
			if(adaptiveValueOld[i] <= 0.0)
			{
				temp.assignmentRange[i] = 0.0;
			}
			else
			{
				temp.assignmentRange[i] = range;
			}
		}
	}

	public static void resetAdaptiveQInNetworkNodes(SamplingAlgortihm sa) {
		for(int i = 0; i < sa.networkNodesArray.length; i++)
		{
			sa.networkNodesArray[i].adaptiveValueNew = new double[sa.networkNodesArray[i].domainSize];
			sa.networkNodesArray[i].adaptiveValueOld = new double[sa.networkNodesArray[i].domainSize];
			sa.networkNodesArray[i].assignmentRange = new double[sa.networkNodesArray[i].domainSize];
			sa.networkNodesArrayCopy[i].adaptiveValueNew = new double[sa.networkNodesArray[i].domainSize];
			sa.networkNodesArrayCopy[i].adaptiveValueOld = new double[sa.networkNodesArray[i].domainSize];
			sa.networkNodesArrayCopy[i].assignmentRange = new double[sa.networkNodesArray[i].domainSize];
		}
	}
}
