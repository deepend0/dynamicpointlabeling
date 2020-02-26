/*
 * ConflictGraph.h
 *
 *  Created on: Dec 17, 2017
 *      Author: oakile
 */

#ifndef _CONFLICTGRAPH_H_
#define _CONFLICTGRAPH_H_

#include <string>
#include <vector>

#include "Graph.h"

namespace labelplacement {

using namespace graph;

class ConflictGraph {
public:
	static void writeConflictGraphIntoFile(std::vector<int>* adjacencyListOfPositions, int pointSize, int positionsPerPoint, std::string filename);

	ConflictGraph(int positionNumberPerPoint, int numberOfPositions, std::vector<int>* adjacencyListOfPositions);
	ConflictGraph();
	static int compareConflictGraphs(int size, std::vector<int>* cg1, std::vector<int>* cg2);
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

#endif /* _CONFLICTGRAPH_H_ */
