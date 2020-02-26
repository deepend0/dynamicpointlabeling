/*
 * AlleleGroupers.cpp
 *
 *  Created on: Aug 19, 2018
 *      Author: oakile
 */

#include <random>

#include "AlleleGroupers.h"
#include "BoundDiameterMinCutClustering.h"

BoundedDiameterMinCutClusterGrouper::BoundedDiameterMinCutClusterGrouper(int depth) {
	this->depth = depth;
}

GroupsT BoundedDiameterMinCutClusterGrouper::operator()(graph::Graph* graph) {
	GroupsT groups = __null;
	graph::clustering::BoundDiameterMinCutClustering clustering;
	std::vector<graph::Graph*>* clusters = clustering.cluster(
			graph, depth);
	int numberOfGroups = clusters->size();
	groups = new std::vector<std::vector<int>*>();
	for (int i = 0; i < numberOfGroups; i++) {
		graph::Graph* cluster = clusters->at(i);
		groups->push_back(new std::vector<int>(*cluster->getVertices()));
		delete cluster;
	}
	delete clusters;
	return groups;
}


RandomPivotDiameterClusterGrouper::RandomPivotDiameterClusterGrouper(int depth, int maxBranchingFactor) {
	this->depth = depth;
	this->maxBranchingFactor = maxBranchingFactor;
}

GroupsT RandomPivotDiameterClusterGrouper::clusterGraph(graph::Graph* graph) {
	const int vertexNumber = graph->getVertexNumber();
	std::vector<std::vector<int>*>* clusters = new std::vector<std::vector<int>*>();

	if (depth > -1) {
		std::vector<int> verticesToSelectLimited;
		std::vector<int> verticesToSelectOverLimit;
		for (int i = 0; i < vertexNumber; i++) {
			if (maxBranchingFactor == 0) {
				verticesToSelectLimited.push_back(i);
			} else {
				int vertexDegree =
						graph->getAdjacencyList()[i].size();
				if (vertexDegree > maxBranchingFactor) {
					verticesToSelectOverLimit.push_back(i);
				} else {
					verticesToSelectLimited.push_back(i);
				}
			}
		}
		std::vector<int> verticesSelected;
		bool verticesSelectedIndicator[vertexNumber];
		for (int i = 0; i < vertexNumber; i++) {
			verticesSelectedIndicator[i] = false;
		}
		int lastDepthCurStart = -1;
		int lastDepthCurEnd = -1;
		int groupNumber = 0;
		for (int currentClusterSize = 0; currentClusterSize < vertexNumber;
				) {
			std::vector<int>* cluster = new std::vector<int>();
			clusters->push_back(cluster);
			groupNumber++;
			std::vector<int>* verticesToSelect;
			if (!verticesToSelectLimited.empty()) {
				verticesToSelect =
						&verticesToSelectLimited;
			} else {
				verticesToSelect =
						&verticesToSelectOverLimit;
			}

			std::random_device rd;
			std::mt19937 gen(rd());
			std::uniform_int_distribution<> vertexIxDist(0, verticesToSelect->size() - 1);
			int vertexToSelectIx = vertexIxDist(gen);
			int vertexToSelect = verticesToSelect->at(vertexToSelectIx);

			verticesToSelect->erase(
					verticesToSelect->begin() + vertexToSelectIx);
			verticesSelected.push_back(vertexToSelect);
			verticesSelectedIndicator[vertexToSelect] = true;
			cluster->push_back(vertexToSelect);
			currentClusterSize++;

			lastDepthCurStart = lastDepthCurEnd = verticesSelected.size() - 1;
			for (int curDepth = 0;
					curDepth < depth && currentClusterSize < vertexNumber;
					curDepth++) {
				for (;
						lastDepthCurStart <= lastDepthCurEnd
								&& currentClusterSize < vertexNumber;
						lastDepthCurStart++) {
					int vertexToFindRivals =
							verticesSelected[lastDepthCurStart];
					std::vector<int> adjacentVertices =
							graph->getAdjacencyList()[vertexToFindRivals];
					if (adjacentVertices.size() > maxBranchingFactor
							&& !verticesToSelectLimited.empty()) {
						continue;
					}
					for (std::vector<int>::iterator it = adjacentVertices.begin();
							it != adjacentVertices.end()
									&& currentClusterSize < vertexNumber;
							it++) {
						int adjacentVertex = *it;
						if (!verticesSelectedIndicator[adjacentVertex]) {
							for (std::vector<int>::iterator it2 =
									verticesToSelect->begin();
									it2 != verticesToSelect->end();
									it2++) {
								if (*it2 == adjacentVertex) {
									verticesToSelect->erase(it2);
									break;
								}
							}
							verticesSelected.push_back(adjacentVertex);
							verticesSelectedIndicator[adjacentVertex] = true;
							cluster->push_back(adjacentVertex);
							currentClusterSize++;
						}
					}
				}
				lastDepthCurEnd = verticesSelected.size() - 1;
			}
		}
	}
	return clusters;
}

GroupsT RandomPivotDiameterClusterGrouper::operator()(graph::Graph* graph) {
	GroupsT groups = __null;
	std::vector<std::vector<int>*>* clusters = clusterGraph(graph);
	int numberOfGroups = clusters->size();
	groups = new std::vector<std::vector<int>*>();
	for (int i = 0; i < numberOfGroups; i++) {
		groups->push_back(new std::vector<int>(*clusters->at(i)));
		delete clusters->at(i);
	}
	delete clusters;

	return groups;
}
