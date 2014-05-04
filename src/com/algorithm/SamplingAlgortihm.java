package com.algorithm;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Comparator;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.PriorityQueue;

import com.utils.AdjacencyListComparator;
import com.utils.FunctionNode;
import com.utils.InputReader;
import com.utils.NetworkNode;
import com.utils.Utilities;

public class SamplingAlgortihm {
	
	String networkType = null;
	public NetworkNode[] networkNodesArray = null;
	public FunctionNode[] functionNodesArray = null;
	public int[] evidenceVariables = null;
	public int[] evidenceVariablesValues = null;
	public NetworkNode[] networkNodesArrayCopy = null;
	public FunctionNode[] functionNodesArrayCopy = null;
	public int[] evidenceVariablesCopy = null;
	public int[] evidenceVariablesValuesCopy = null;
	LinkedList<HashSet<Integer>> treeDecomposition = null;
	int[] variablesCountInClusters;
	int w = 1;
	int N = 1;
//	int[] Narray = {1000};
	int[] Narray = {100,1000,10000,20000};
	HashSet<Integer> wCutset = null;
	int[] wCutsetArray = null;
	int[] wCutsetVariablesValues = null;
	double z = 0;
	double q = 1;
	double evidence;
	double[] zValues;
	float[] timeValues;
	double zMean;
	float timeMean;
	double zDeviation;
	double timeDeviation;
	
	public static void main(String[] args) throws NumberFormatException, IOException {
		
		SamplingAlgortihm sa = new SamplingAlgortihm();
		InputReader.readInput(sa, args[0]);
		InputReader.readEvidence(sa, args[0]);

		Utilities.generateCopyOfNetworkNodes(sa);
		Utilities.generateCopyOfFunctionNodes(sa);
		Utilities.generateCopyOfEvidence(sa);
		
		instantiateEvidences(sa, sa.evidenceVariables, sa.evidenceVariablesValues, sa.functionNodesArray, sa.networkNodesArray);

		generateInteractionGraph(sa.functionNodesArray, sa.networkNodesArray, sa.evidenceVariables);
		
		int[] bucketsMinDegreeOrder = generateMinDegree(sa.networkNodesArray, sa.evidenceVariables);
		
//		Utilities.printElementsInIntegerArray(bucketsMinDegreeOrder);
		
		BufferedWriter output = new BufferedWriter(new FileWriter(new File(args[1])));
		
		for(int i = 1; i <= 5; i++)
		{

			sa.zValues = new double[10];
			sa.timeValues = new float[10];
			
			sa.w = i;
			sa.getTreeDecomposition(sa, bucketsMinDegreeOrder);
			sa.getWCutSet(sa);
			
			sa.wCutsetArray = new int[sa.wCutset.size()];
			int p = 0;
			for(Integer temp : sa.wCutset)
			{
				sa.wCutsetArray[p] =  temp.intValue();
				p++;
			}
			
			for(int j = 0;j < sa.Narray.length; j++)
			{
				sa.N = sa.Narray[j];
				for(int d = 1; d <=  10; d++)
				{
					System.out.println("For W = " + sa.w + " N = " + sa.N + " Run = " + d +"\n");
					output.write("For W = " + sa.w + " N = " + sa.N + "\n");
					long startTime = System.currentTimeMillis();
					calculateUniformQ(sa);
					
					Utilities.resetAdaptiveQInNetworkNodes(sa);
					
					for(int n = 0; n < sa.N; n++)
					{
						if(n != 0 && n%100 == 0)
						{
							assignNewAdaptiveToOldAdaptive(sa);
//							System.out.print("Updated Q  ---  ");
						}
						sa.wCutsetVariablesValues = new int[sa.wCutset.size()];
						if(n <= 100)
						{
							for(int k = 0; k < sa.wCutsetArray.length; k++)
							{
								int domainSize = sa.networkNodesArray[sa.wCutsetArray[k]].domainSize;
								int assignment = 0 + (int)(Math.random() * domainSize);
								sa.wCutsetVariablesValues[k] =  assignment;
//								System.out.print(assignment);
							}
//							System.out.println("");
						}
						else
						{
							for(int k = 0; k < sa.wCutsetArray.length; k++)
							{
								double assignmentValue = Math.random();
								int assignment = Utilities.getValueForAVariableWithAssignmentValue(sa.networkNodesArray[sa.wCutsetArray[k]], assignmentValue);
								sa.wCutsetVariablesValues[k] =  assignment;
							}
						}
						
						if(n > 100)
						{
							calculateAdaptiveQ(sa);
						}
						
//						System.out.println(sa.q);
//						Utilities.printElementsInIntegerArray(sa.wCutsetVariablesValues);
						instantiateEvidences(sa, sa.wCutsetArray, sa.wCutsetVariablesValues, sa.functionNodesArray, sa.networkNodesArray);
						bucketElimination(bucketsMinDegreeOrder, sa, sa.networkType);

						//generating a copy.
						sa.networkNodesArray = sa.networkNodesArrayCopy;
						sa.functionNodesArray = sa.functionNodesArrayCopy;
						sa.evidenceVariables = sa.evidenceVariablesCopy;
						sa.evidenceVariablesValues = sa.evidenceVariablesValuesCopy;
						
						//Update the adaptive Q is newQ in every network node.
						updateAdaptiveNewQ(sa);
						
						Utilities.generateCopyOfNetworkNodes(sa);
						Utilities.generateCopyOfFunctionNodes(sa);
						Utilities.generateCopyOfEvidence(sa);
					}
					
					//Reset the old and new Q values in NetworkNodes
					Utilities.resetQInNetworkNodes(sa);
					
					sa.z = sa.z / (sa.N);
					System.out.println("sampled probability is : " + sa.z + "\n");
					output.write("sampled probability is : " + sa.z + "\n");
					long endTime = System.currentTimeMillis();
					output.write("Execution Time : " + (endTime - startTime)+ "\n");
					sa.zValues[d-1] = Math.abs((Math.log10(Double.parseDouble(args[2])) - Math.log10(sa.z))/Math.log10(Double.parseDouble(args[2])));
					sa.timeValues[d-1] = (float)(endTime - startTime)/60000;
					sa.z = 0;
				}
				
				sa.calculateMeanAndStandardDeviation(sa.zValues, sa.timeValues, Math.log10(Double.parseDouble(args[2])), Float.parseFloat(args[3]));
				System.out.println("Mean Z : " + sa.zMean);
				System.out.println("Deviation Z : " + sa.zDeviation);
				System.out.println("Mean time : " + sa.timeMean);
				System.out.println("Deviation Time : " + sa.timeDeviation);
				System.out.println("");
				output.write("\n\nMean Z : " + sa.zMean+ "\n");
				output.write("Deviation Z : " + sa.zDeviation+ "\n");
				output.write("Mean time : " + sa.timeMean+ "\n");
				output.write("Deviation Time : " + sa.timeDeviation+ "\n"+ "\n");
			}
		}
		output.close();
	}
	
	/**
	 * This method calculates the mean and standard deviation
	 * @param timeValues 
	 * @param zValues
	 * @param originalZ 
	 * @param originalTime 
	 */
	private void calculateMeanAndStandardDeviation(double[] zValues, float[] timeValues, double originalZ, float originalTime) {
		double zvalue = 0.0d;
		float time = 0.0f;
		float timeDeviation = 0.0f;
		double deviation = 0.0d;
		for(int i = 0; i < zValues.length; i ++)
		{
			zvalue += zValues[i];
			time += timeValues[i];
		}
		this.zMean = (double)zvalue/10;
		this.timeMean = (float)time/10;
		for(int i = 0; i < zValues.length; i ++)
		{
			deviation += Math.pow((this.zMean- zValues[i]), 2);
			timeDeviation += Math.pow((this.timeMean - timeValues[i]), 2);
		}
		this.zDeviation = Math.sqrt(deviation/10);
		this.timeDeviation = Math.sqrt(timeDeviation/10);
	}

	
	
	/**
	 * This method updates the adaptive Q
	 * @param sa
	 */
	private static void updateAdaptiveNewQ(SamplingAlgortihm sa) {
		NetworkNode temp = null;
		NetworkNode tempcopy = null;
		for(int i = 0; i < sa.wCutsetArray.length; i++)
		{
			temp = sa.networkNodesArray[sa.wCutsetArray[i]];
			tempcopy = sa.networkNodesArrayCopy[sa.wCutsetArray[i]];
			temp.adaptiveValueNew[sa.wCutsetVariablesValues[i]] = temp.adaptiveValueNew[sa.wCutsetVariablesValues[i]] + sa.evidence;
			tempcopy.adaptiveValueNew[sa.wCutsetVariablesValues[i]] = tempcopy.adaptiveValueNew[sa.wCutsetVariablesValues[i]] + sa.evidence;
		}
	}
	
	/**
	 * This method assigns values in new into the old after calculating the probabilities and also the range assignment.
	 * @param sa
	 */
	private static void assignNewAdaptiveToOldAdaptive(SamplingAlgortihm sa) {
		double valueSum = 1.0;
		if(sa.wCutsetArray.length > 0)
		{
			valueSum = Utilities.sumUpAdaptiveValues(sa.networkNodesArray[sa.wCutsetArray[0]].adaptiveValueNew);
		}
		for(int i = 0; i < sa.wCutsetArray.length; i++)
		{
			NetworkNode temp = sa.networkNodesArray[sa.wCutsetArray[i]];
			NetworkNode tempCopy = sa.networkNodesArrayCopy[sa.wCutsetArray[i]];
			Utilities.updateAdaptiveOldQ(temp.adaptiveValueOld, temp.adaptiveValueNew, temp, valueSum);
			Utilities.updateAdaptiveOldQ(tempCopy.adaptiveValueOld, tempCopy.adaptiveValueNew, tempCopy, valueSum); 
//			temp.adaptiveValueOld = new double[temp.adaptiveValueNew.length]; 
//			tempCopy.adaptiveValueOld = new double[temp.adaptiveValueNew.length]; 
//			System.arraycopy(temp.adaptiveValueNew, 0, temp.adaptiveValueOld, 0, temp.adaptiveValueNew.length);
//			System.arraycopy(tempCopy.adaptiveValueNew, 0, tempCopy.adaptiveValueOld, 0, tempCopy.adaptiveValueNew.length);
			temp.adaptiveValueNew = new double[temp.domainSize];
			tempCopy.adaptiveValueNew = new double[tempCopy.domainSize];
		}
	}
	
	/**
	 * This method calculates the adaptive Q
	 * @param sa
	 */
	private static void calculateAdaptiveQ(SamplingAlgortihm sa) {
		sa.q = 1.0;
		for(int i = 0; i < sa.wCutsetArray.length; i++)
		{
			NetworkNode temp = sa.networkNodesArray[sa.wCutsetArray[i]];
			sa.q = (Double)(sa.q * temp.adaptiveValueOld[sa.wCutsetVariablesValues[i]]);
		}
	}

	/**
	 * This method calculates uniform Q
	 * @param sa
	 */
	private static void calculateUniformQ(SamplingAlgortihm sa) {
		sa.q = 1.0;
		for(Integer i : sa.wCutset)
		{
			sa.q = (Double)sa.q/sa.networkNodesArray[i.intValue()].domainSize;
		}
	}

	/**
	 * This method returns the tree decomposition of the graph that would be passed to it.
	 * @param sa
	 * @param bucketsMinDegreeOrder 
	 */
	private void getTreeDecomposition(SamplingAlgortihm sa, int[] bucketsMinDegreeOrder) {
		sa.treeDecomposition = new LinkedList<HashSet<Integer>>();
		LinkedList<HashSet<Integer>> temporaryVariablesSetList = new LinkedList<HashSet<Integer>>();
		for(int b: bucketsMinDegreeOrder)
		{
//			System.out.println("Processing : " + b);
			setVariablesForBucket(sa.functionNodesArray, b, sa.networkNodesArray, temporaryVariablesSetList, sa);
		}
	}


	private void setVariablesForBucket(FunctionNode[] functionNodesArray2,
			int b, NetworkNode[] networkNodesArray2,
			LinkedList<HashSet<Integer>> temporaryVariablesSetList, SamplingAlgortihm sa) {
		HashSet<Integer> variablesSet = new HashSet<Integer>();
		for(int i = 0; i < functionNodesArray.length; i++)
		{
			if(!functionNodesArray[i].isVisitedForTreeDecomposition && Utilities.isVariablePresentInFunction(b, functionNodesArray[i].nodesValuesArray))
			{
				Utilities.setVariablesOfBucket(variablesSet, functionNodesArray[i].nodesValuesArray);
				functionNodesArray[i].isVisitedForTreeDecomposition = true;
			}
		}
		HashSet<Integer> tempVariablesSet = null;
		for(int i = 0; i < temporaryVariablesSetList.size();)
		{
			tempVariablesSet = temporaryVariablesSetList.get(i);
			if(tempVariablesSet.contains(b))
			{
				Utilities.setVariablesOfBucket(variablesSet, tempVariablesSet);
				temporaryVariablesSetList.remove(i);
			}
			else
			{
				 i++;
			}
		}
		sa.treeDecomposition.add(variablesSet);
		if(variablesSet != null)
		{
			variablesSet.remove(b);
			temporaryVariablesSetList.add(variablesSet);
		}
	}

	/**
	 * This method returns the Wcutset for the network
	 * @param sa
	 */
	private void getWCutSet(SamplingAlgortihm sa) {
		sa.wCutset = new HashSet<Integer>();
//		System.out.println(sa.treeDecomposition.size());
		sa.variablesCountInClusters = new int[sa.networkNodesArray.length];
		int max = 0;
		for(int i = 0; i < sa.treeDecomposition.size(); i++)
		{
			HashSet<Integer> tempCluster = sa.treeDecomposition.get(i);
			if(max < tempCluster.size())
			{
				max = tempCluster.size();
			}
			for(Integer temp : tempCluster)
			{
				sa.variablesCountInClusters[temp] = sa.variablesCountInClusters[temp]+1; 
			}
		}
		HashSet<Integer> bucket = null;  
		while((bucket = Utilities.getBucketWithClusterSizeGreaterThan(sa.w, sa.treeDecomposition)) != null)
		{
			int variable = Utilities.getVariableFromBucketWithMaximumCountInClusters(bucket, sa.variablesCountInClusters);
			sa.wCutset.add(variable);
			for(HashSet<Integer> temp : sa.treeDecomposition)
			{
				temp.remove(variable);
			}
		}
//		System.out.println("Max cluster size is : " + max);
	}
	
	/**
	 * This method does the bucket elimination.
	 * @param bucketsMinDegreeOrder
	 * @param sa
	 * @param networkType 
	 */
	private static void bucketElimination(int[] bucketsMinDegreeOrder,
			SamplingAlgortihm sa, String networkType) {
		LinkedList<FunctionNode> temporaryFunctions = new LinkedList<FunctionNode>();
		for(int b: bucketsMinDegreeOrder)
		{
			updateFunctionsForBucket(sa.functionNodesArray, b, sa.networkNodesArray, temporaryFunctions);
		}
		sa.evidence = getProductOfEmptyVariableFunction(sa.functionNodesArray, temporaryFunctions);
		if(sa.q != 0)
		{
			sa.z += (sa.evidence/sa.q);
		}
		else
		{
			calculateUniformQ(sa);
			sa.z += (sa.evidence/sa.q);
		}
//		System.out.println(sa.z + " - " +sa.q);
	}
	
	public static void updateFunctionsForBucket(FunctionNode[] functionNodesArray, int b, NetworkNode[] networkNodesArray, LinkedList<FunctionNode> temporaryFunctions) {
		FunctionNode node = null;
		for(int i = 0; i < functionNodesArray.length; i++)
		{
			if(!functionNodesArray[i].isVisited && functionNodesArray[i].nodesValuesArray.length > 0 && Utilities.isVariablePresentInFunction(b, functionNodesArray[i].nodesValuesArray))
			{
				if(node != null)
				{
					node = factorProductOfTwoFunctions(node, functionNodesArray[i], networkNodesArray);
				}
				else
				{
					node = functionNodesArray[i];
				}
				functionNodesArray[i].isVisited = true;
			}
		}
		for(FunctionNode tempNode : temporaryFunctions)
		{
			if(!tempNode.isVisited && tempNode.nodesValuesArray.length > 0 && Utilities.isVariablePresentInFunction(b, tempNode.nodesValuesArray))
			{
				if(node != null)
				{
					node = factorProductOfTwoFunctions(node, tempNode, networkNodesArray);
				}
				else
				{
					node = tempNode;
				}
				tempNode.isVisited = true;
			}
		}
		if(node != null)
		{
			node = factorSum(node, b, networkNodesArray);
			temporaryFunctions.add(node);
		}
	}
	
	/**
	 * This function returns factor product of two functions.
	 * @param node2 
	 * @param node1 
	 * @param networkNodesArray2 
	 * @return valuesinFunction
	 */
	private static FunctionNode factorProductOfTwoFunctions(FunctionNode node1, FunctionNode node2, NetworkNode[] networkNodesArray2) {
		int j = 0;
		int k = 0;
		HashSet<Integer> variablesSet = Utilities.getHashSetForVariablesInFunctions(node1, node2);
		int noOfValuesInFunction = getNoOfValuesinFunction(variablesSet, networkNodesArray2);
		double[] valuesInFunction = new double[noOfValuesInFunction];
		int[] variablesArray = new int[variablesSet.size()];
		int y = 0;
		int maxElement = 0;
		for(int l : variablesSet)
		{
			if(l > maxElement)
			{
				maxElement = l;
			}
			variablesArray[y] = l;
			y++;
		}
		int[] assignment = new int[maxElement+1];
		int temp = -100;
		for(int i = 0; i <= noOfValuesInFunction - 1; i++)
		{
			valuesInFunction[i] = node1.functionValues[j] * node2.functionValues[k];
//  			System.out.println( node1.functionValues[j] + " * " + node2.functionValues[k] + " ---> " + valuesInFunction[i]);
//  			System.out.println(i + " - " + valuesInFunction[i]);
  			for(int x = variablesArray.length-1; x >= 0; x--)
			{
  				int l = variablesArray[x];
//				System.out.println(i);
				assignment[l] = assignment[l]+1;
				if(assignment[l] == networkNodesArray2[l].domainSize)
				{
					assignment[l] = 0;
					temp = -100;
					if((temp = Utilities.getStrideForVariableInFunction(node1, l)) > 0 && j > 0)
					{
						j = j - (networkNodesArray2[l].domainSize - 1 ) * temp;
					}
					temp = -100;
					if((temp = Utilities.getStrideForVariableInFunction(node2, l)) > 0  && k > 0)
					{
						k = k - (networkNodesArray2[l].domainSize - 1 ) * temp;
					}
				}
				else
				{
					temp = -100;
					if((temp = Utilities.getStrideForVariableInFunction(node1, l)) > 0 )
					{
						j = j + temp;
					}
					temp = -100;
					if((temp = Utilities.getStrideForVariableInFunction(node2, l)) > 0)
					{
						k = k + temp;
					}
					break;
				}
			}
//			System.out.println(i + " - " + j + " - " + k);
		}
		return new FunctionNode(networkNodesArray2, variablesArray, valuesInFunction);
	}

	
	/**
	 * This function returns the size of the function table calculated for product of two factors
	 * @param variablesSet
	 * @param networkNodesArray2
	 * @return
	 */
	private static int getNoOfValuesinFunction(HashSet<Integer> variablesSet,
			NetworkNode[] networkNodesArray2) {
		int noOfValuesInFunction = 1;
		for(int i : variablesSet)
		{
			noOfValuesInFunction = noOfValuesInFunction * networkNodesArray2[i].domainSize;
		}
		return noOfValuesInFunction;
	}
		
	/**
	 * This function gives product of the functions that have no variables in it and are not visited..
	 * @param functionNodesArray
	 * @param temporaryFunctions
	 * @return
	 */
	private static double getProductOfEmptyVariableFunction(
			FunctionNode[] functionNodesArray,
			LinkedList<FunctionNode> temporaryFunctions) {
		double product = 1.0;
		for(int i = 0; i < functionNodesArray.length; i++)
		{
			if(!functionNodesArray[i].isVisited && functionNodesArray[i].nodesValuesArray.length == 0 && functionNodesArray[i].functionValues.length == 1)
			{
				product = product * functionNodesArray[i].functionValues[0];
			}
		}
		for(FunctionNode tempFunction : temporaryFunctions)
		{
			if(!tempFunction.isVisited && tempFunction.nodesValuesArray.length == 0 && tempFunction.functionValues.length == 1)
			{
				product = product * tempFunction.functionValues[0];
			}
		}
		return product;
	}
	
	/**
	 * This method generates the min order computation of the given graph
	 * @param networkNodesArray 
	 * @param evidenceVariables 
	 * @return 
	 */
	private static int[] generateMinDegree(NetworkNode[] networkNodesArray, int[] evidenceVariables) {
		int length = networkNodesArray.length;
		Comparator<NetworkNode> comparator = new AdjacencyListComparator();
        PriorityQueue<NetworkNode> queue = 
            new PriorityQueue<NetworkNode>(length, comparator);
        for(int i = 0; i < length; i++)
        {
        	if(!Utilities.isVariablePresentInEvidenceVariables(networkNodesArray[i].nodeValue, evidenceVariables))
        	{
        		queue.add(networkNodesArray[i]);
        	}
        }
        int[] bucketsMinDegreeOrder = new int[queue.size()];
        int i = 0;
        while(!queue.isEmpty())
        {
        	NetworkNode node = queue.poll();
        	Utilities.updateAdjacentNodes(node, networkNodesArray);
        	for(int k: node.adjacentNodes)
        	{
        		if(queue.remove(networkNodesArray[k]))
        		{
        			queue.add(networkNodesArray[k]);
        		}
        	}
        	bucketsMinDegreeOrder[i++] = node.nodeValue;
        }
        return bucketsMinDegreeOrder;
	}

	
	/**
	 * This method generates interaction graph for given nodes and factors
	 * @param evidenceVariables 
	 * @param functionNodesArray2
	 * @param networkNodesArray2
	 */
	private static void generateInteractionGraph(
			FunctionNode[] functionNodesArray, NetworkNode[] networkNodesArray, int[] evidenceVariables) {

		for(int j = 0; j < functionNodesArray.length; j++)
		{
			Utilities.buildAdjacencyListForEachVariableInFunction(functionNodesArray[j], networkNodesArray, evidenceVariables);
		}
	}
	
	/**
	 * This function instantiates all the evidences for a given function.
	 * @param sa
	 * @param evidenceVariables
	 * @param evidenceVariablesValues
	 * @param functionNodesArray
	 * @param networkNodesArray
	 */
	private static void instantiateEvidences(SamplingAlgortihm sa, int[] evidenceVariables,
			int[] evidenceVariablesValues, FunctionNode[] functionNodesArray, NetworkNode[] networkNodesArray) {
		for(int i = 0; i < evidenceVariables.length; i++)
		{
			for(int j = 0; j < functionNodesArray.length; j++)
			{
				int domainSizeOfVariable = networkNodesArray[evidenceVariables[i]].domainSize;
				if(Utilities.isVariablePresentInFunction(evidenceVariables[i], functionNodesArray[j].nodesValuesArray))
				{
					sa.functionNodesArray[j] = instantiateEvidenceForVariableAndFunction(evidenceVariables[i], evidenceVariablesValues[i], functionNodesArray[j], domainSizeOfVariable, networkNodesArray);
				}
			}
		}
	}
	
	/**
	 * This method instantiates the variable on a function for a variable.
	 * @param be 
	 * @param evidenceVariable
	 * @param evidenceVariableValue
	 * @param functionNode
	 * @param domainSizeOfVariable 
	 * @param networkNodesArray 
	 * @return 
	 */
	private static FunctionNode instantiateEvidenceForVariableAndFunction(int evidenceVariable, int evidenceVariableValue,
			FunctionNode functionNode, int domainSizeOfVariable, NetworkNode[] networkNodesArray) {
		int strideofVariable = Utilities.getStrideForVariableInFunction(functionNode, evidenceVariable);
		int stepIncrement = 0;
		int domainIncrementer = 0;
		for(int i = 0; i < functionNode.functionValues.length; i++)
		{
			if(strideofVariable == stepIncrement)
			{
				stepIncrement = 1;
				domainIncrementer++;
			}
			else
			{
				stepIncrement++;
			}
			if(domainIncrementer == domainSizeOfVariable)
			{
				domainIncrementer = 0;
			}
			if(domainIncrementer != evidenceVariableValue)
			{
				functionNode.functionValues[i] = 0.0;
			}
			else
			{
				i = i+ strideofVariable-1;
				stepIncrement = stepIncrement + strideofVariable - 1;
			}
		}
		functionNode = factorSum(functionNode, evidenceVariable, networkNodesArray);
		return functionNode;
	}
	
	/**
	 * This method returns the summation of a factor over a given range of variables
	 * @param node
	 * @param variableToBeSummed
	 * @param networkNodesArray
	 * @return
	 */
	private static FunctionNode factorSum(FunctionNode node, int variableToBeSummed, NetworkNode[] networkNodesArray) {
		int j = 0;
		int noOfValuesInFunction = node.functionValues.length/ networkNodesArray[variableToBeSummed].domainSize;
		double[] valuesInFunction = new double[noOfValuesInFunction];
		int[] newVariablesArray = Utilities.getNewVariablesArray(variableToBeSummed, node.nodesValuesArray);
		int variableStride = Utilities.getStrideForVariableInFunction(node, variableToBeSummed);
		boolean[] visitedIndexes = new boolean[node.functionValues.length];
		for(int i = 0; i <= node.functionValues.length - 1; i++)
		{
			if(visitedIndexes[i] == false)
			{
				for(int k = 0; k < networkNodesArray[variableToBeSummed].domainSize; k++)
				{
					int stepIncrementer = i + (k * variableStride);
					valuesInFunction[j] = valuesInFunction[j] + node.functionValues[stepIncrementer];
					visitedIndexes[stepIncrementer] = true;
				}
//				System.out.println(j + "   ----->   " +valuesInFunction[j]);
				j++;
			}
		}
		return new FunctionNode(networkNodesArray, newVariablesArray, valuesInFunction);
	}
	
	public void setNetworkType(String networkType) {
		this.networkType = networkType;
	}

	public void setNetworkNodesArray(NetworkNode[] networkNodesArray) {
		this.networkNodesArray = networkNodesArray;
	}

	public void setFunctionNodesArray(FunctionNode[] functionNodesArray) {
		this.functionNodesArray = functionNodesArray;
	}


	public void setEvidenceVariables(int[] evidenceVariables) {
		this.evidenceVariables = evidenceVariables;
	}

	public void setEvidenceVariablesValues(int[] evidenceVariablesValues) {
		this.evidenceVariablesValues = evidenceVariablesValues;
	}

}
