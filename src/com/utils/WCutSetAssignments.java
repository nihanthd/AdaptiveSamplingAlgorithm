package com.utils;

public class WCutSetAssignments {
	public int[] wCutSetAssignment;
	public double value;
	
	public WCutSetAssignments(int[] wCutsetVariablesValues, double evidence) {
		this.wCutSetAssignment = wCutsetVariablesValues;
		this.value = evidence;
	}
}
