/*
 * ConflictGraph.h
 *
 *  Created on: Dec 17, 2017
 *      Author: oakile
 */

#ifndef _CONFLICTGRAPH_H_
#define _CONFLICTGRAPH_H_

#include <vector>

#include "Graph.h"

namespace labelplacement {

using namespace graph;

class ConflictGraph {
public:
	ConflictGraph(int positionNumberPerPoint, int numberOfPositions, std::vector<int>* adjacencyListOfPositions);
	ConflictGraph();
	virtual ~ConflictGraph();
	Graph* getConflictGraphOfPoints() const;
	void setConflictGraphOfPoints(Graph* conflictGraphOfPoints);
	Graph* getConflictGraphOfPositions() const;
	void setConflictGraphOfPositions(Graph* conflictGraphOfPositions);
	int getPositionNumberPerPoint() const;
	void setPositionNumberPerPoint(int positionNumberPerPoint);
	void initializeOf(int positionNumberPerPoint, int numberOfPositions, std::vector<int>* adjacencyListOfPositions);

private:
	Graph* conflictGraphOfPoints;
	Graph* conflictGraphOfPositions;
	int positionNumberPerPoint;
};

} /* namespace labelplacement */

#endif /* SRC_CONFLICTGRAPH_H_ */
