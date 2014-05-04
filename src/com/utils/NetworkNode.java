package com.utils;

import java.util.HashSet;


public class NetworkNode{
	public int nodeValue;
	public int domainSize;
	public int domainValues[];
	public HashSet<Integer> adjacentNodes;
	public int[] functionsArray;
	public double[] adaptiveValueOld;
	public double[] adaptiveValueNew;
	public double[] assignmentRange;
	
	
	public NetworkNode(int i, String domainSize) {
		this.nodeValue = i;
		this.domainSize = Integer.parseInt(domainSize);
		this.domainValues = new int[this.domainSize];
		this.adjacentNodes = new HashSet<Integer>();
		for(int j = 0; j < this.domainSize; j++)
		{
			domainValues[j] = j;
		}
		this.adaptiveValueOld = new double[this.domainSize];
		this.adaptiveValueNew = new double[this.domainSize];
		this.assignmentRange = new double[this.domainSize];
	}

	/**
	 * this constructor generates a copy of the network node 
	 * @param networkNode
	 */
	public NetworkNode(NetworkNode networkNode) {
		this.nodeValue = networkNode.nodeValue;
		this.domainSize = networkNode.domainSize;
		this.domainValues = networkNode.domainValues;
		this.adjacentNodes = new HashSet<Integer>();
		for(Integer temp: networkNode.adjacentNodes)
		{
			Integer temp1 = new Integer(temp.intValue());
			this.adjacentNodes.add(temp1);
		}
		/*for(int j = 0; j < this.domainSize; j++)
		{
			this.domainValues[j] = j;
		}*/
		this.adaptiveValueOld = networkNode.adaptiveValueOld;
		this.adaptiveValueNew = networkNode.adaptiveValueNew;
		this.assignmentRange = networkNode.assignmentRange;
//		this.adaptiveValueOld = new double[this.domainSize];
//		this.adaptiveValueNew = new double[this.domainSize];
//		System.arraycopy(networkNode.adaptiveValueOld, 0, this.adaptiveValueOld, 0, this.domainSize);
//		System.arraycopy(networkNode.adaptiveValueNew, 0, this.adaptiveValueNew, 0, this.domainSize);
	}
}
