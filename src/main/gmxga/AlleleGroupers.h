/*
 * AlleleGroupers.h
 *
 *  Created on: Aug 20, 2018
 *      Author: oakile
 */

#ifndef _GMXGA_ALLELEGROUPERS_H_
#define _GMXGA_ALLELEGROUPERS_H_

#include <vector>

#include "Graph.h"

typedef std::vector<std::vector<int>*>* GroupsT;
typedef GroupsT AlleleGrouper(graph::Graph* graph);

struct BoundedDiameterMinCutClusterGrouper {

	int depth;
	BoundedDiameterMinCutClusterGrouper(int depth);
	GroupsT operator()(graph::Graph* graph);
};

struct RandomPivotDiameterClusterGrouper {

	int depth;
	int maxBranchingFactor;
	GroupsT operator()(graph::Graph* graph);
	GroupsT clusterGraph(graph::Graph* graph);
};

#endif /* _GMXGA_ALLELEGROUPERS_H_ */
