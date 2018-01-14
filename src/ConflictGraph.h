/*
 * ConflictGraph.h
 *
 *  Created on: Dec 17, 2017
 *      Author: oakile
 */

#ifndef _CONFLICTGRAPH_H_
#define _CONFLICTGRAPH_H_

#include <vector>

namespace labelplacement {

class ConflictGraph {
public:
	ConflictGraph();
	virtual ~ConflictGraph();

	int getPointNumber() const;

	void setPointNumber(int pointNumber);

	int getPositionNumber() const;

	void setPositionNumber(int positionNumber);

	bool** getAdjacencyMatrix() const;

	void setAdjacencyMatrix(bool** adjacencyMatrix);

	std::vector<int>* getConflictingPositions() const;

	void setConflictingPositions(std::vector<int>* conflictingPositions);

	std::vector<int>* getConflictingPoints() const;

	void setConflictingPoints(std::vector<int>* conflictingPoints);

private:
	int pointNumber;
	int positionNumber;
	bool** adjacencyMatrix;
	std::vector<int>* conflictingPositions;
	std::vector<int>* conflictingPoints;

};

} /* namespace labelplacement */

#endif /* SRC_CONFLICTGRAPH_H_ */
