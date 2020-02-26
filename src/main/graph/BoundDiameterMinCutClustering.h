/*
 * BoundDiamMinCutClustering.h
 *
 *  Created on: Jan 31, 2018
 *      Author: oakile
 */

#ifndef _BOUNDDIAMMINCUTCLUSTERING_H_
#define _BOUNDDIAMMINCUTCLUSTERING_H_

#include <vector>

#include "Graph.h"
namespace graph {
namespace clustering {

class BoundDiameterMinCutClustering {
public:
	BoundDiameterMinCutClustering();
	virtual ~BoundDiameterMinCutClustering();
	Graph* getGraph() const;
	void setGraph(Graph* graph);
	std::vector<Graph*>* findComponentsBFS(Graph* graph, int vertexSizeLimit);
	void floydWarshallAllPairsShortestPaths(Graph* component,
			int **shortestPaths);
	int getDiameterOfComponent(Graph* component);
	std::vector<std::pair<int, int>*>* randomCut(Graph* component,
			Graph* subComponent1, Graph* subComponent2);
	Graph** boundedDiameterMinCut(Graph* component, int diameter);
	std::vector<Graph*>* cluster(graph::Graph* graph, int diameter);
private:
	int diameter;
	Graph* graph;
	int nextVertex(int* vertices, int size);
	void buildGraph(std::vector<int>* vertexNos, std::vector<int>* vertices,
			std::vector<std::pair<int, int>*>* edges, Graph *graph);
};

}
}

#endif /* _BOUNDDIAMMINCUTCLUSTERING_H_ */
