/*
 * RandomDynamicConflictGraphGenerator.h
 *
 *  Created on: May 26, 2018
 *      Author: oakile
 */

#ifndef _RANDOMDYNAMICCONFLICTGRAPHGENERATOR_H_
#define _RANDOMDYNAMICCONFLICTGRAPHGENERATOR_H_

#include <string>
#include <vector>

#include "ConflictGraph.h"
#include "ConflictGraphGenerator.h"
#include "MovingPointsGenerator.h"

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
struct RandomDynamicConflictGraphGeneratorParameters {
	RandomDynamicConflictGraphGeneratorParameters(int areaWidth, int areaHeight, int numberOfPoints, int positionsPerPoint,
			int labelWidth, int labelHeight);
	~RandomDynamicConflictGraphGeneratorParameters();
	int areaWidth;
	int areaHeight;
	int numberOfPoints;
	int labelWidth;
	int labelHeight;
	int positionsPerPoint;
	bool mpParamsIndicator = false;
	double maxSpd;
	double meanSpd;
	double stdSpd;
	double maxDir;
	double meanDir;
	double stdDir;
};
class RandomDynamicConflictGraphGenerator : public ConflictGraphGenerator {
public:
	RandomDynamicConflictGraphGenerator(RandomDynamicConflictGraphGeneratorParameters parameters);
	virtual ~RandomDynamicConflictGraphGenerator();
	virtual labelplacement::ConflictGraph* generate(int forwardStates);
private:
	labelplacement::ConflictGraph* generateConflictGraph(std::vector<MovingPoint*>* movingPoints);
	bool intersectionOccurs(Rectangle& r1, Rectangle& r2);
	int labelWidth;
	int labelHeight;
	int positionsPerPoint;
	MovingPointsGenerator* mpGen;
	std::vector<MovingPoint*>* movingPoints;
};

#endif /* _RANDOMDYNAMICCONFLICTGRAPHGENERATOR_H_ */
