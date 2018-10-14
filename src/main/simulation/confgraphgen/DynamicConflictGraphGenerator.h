/*
 * DynamicConflictGraphGenerator.h
 *
 *  Created on: May 26, 2018
 *      Author: oakile
 */

#ifndef _DYNAMICCONFLICTGRAPHGENERATOR_H_
#define _DYNAMICCONFLICTGRAPHGENERATOR_H_

#include <string>
#include <vector>

#include "MovingPointsGenerator.h"
#include "ConflictGraph.h"

class Rectangle {
public:
	Rectangle();
	Rectangle(Position upRightPos, int width, int height);
	Position upRight = Position(0,0);
	Position upLeft = Position(0,0);
	Position downRight = Position(0,0);
	Position downLeft = Position(0,0);

	void print();
};
struct DynamicConflictGraphGeneratorParameters {
	DynamicConflictGraphGeneratorParameters(int areaWidth, int areaHeight, int numberOfPoints, int positionsPerPoint,
			int labelWidth, int labelHeight);
	int areaWidth;
	int areaHeight;
	int numberOfPoints;
	int labelWidth;
	int labelHeight;
	int positionsPerPoint;
};
class DynamicConflictGraphGenerator {
public:
	DynamicConflictGraphGenerator();
	DynamicConflictGraphGenerator(DynamicConflictGraphGeneratorParameters parameters);
	virtual ~DynamicConflictGraphGenerator();
	virtual labelplacement::ConflictGraph* generate(int forwardStates);
private:
	labelplacement::ConflictGraph* generateConflictGraph(std::vector<MovingPoint*>* movingPoints);
	bool intersectionOccurs(Rectangle& r1, Rectangle& r2);
	int compareConflictGraphs(int size, std::vector<int>* cg1, std::vector<int>* cg2);
	void writeConflictGraphIntoFile(std::vector<int>* adjacencyListOfPositions, int pointSize, std::string filename);
	int labelWidth;
	int labelHeight;
	int positionsPerPoint;
	MovingPointsGenerator mpGen;
	std::vector<MovingPoint*>* movingPoints;
};

#endif /* _DYNAMICCONFLICTGRAPHGENERATOR_H_ */
