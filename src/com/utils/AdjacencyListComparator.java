package com.utils;

import java.util.Comparator;

public class AdjacencyListComparator implements Comparator<NetworkNode> {

	@Override
	public int compare(NetworkNode node1, NetworkNode node2) {
		
		if (node1.adjacentNodes.size() < node2.adjacentNodes.size())
        {
            return -1;
        }
		if (node1.adjacentNodes.size() > node2.adjacentNodes.size())
        {
            return 1;
        }
		return 0;
	}

}
