/*
 * DynamicConflictGraphGenerator.cpp
 *
 *  Created on: May 26, 2018
 *      Author: oakile
 */

#include <iostream>
#include <cmath>
#include "RandomDynamicConflictGraphGenerator.h"

using namespace std;


Rectangle::Rectangle() {
}

Rectangle::Rectangle(Position upRight, int width, int height) {
	this->upRight = upRight;
	this->upLeft = Position(upRight.x-width, upRight.y);
	this->downRight = Position(upRight.x, upRight.y-height);
	this->downLeft = Position(upRight.x-width, upRight.y-height);
}

void Rectangle::print() {
	cout<<"UR("<<this->upRight.x<<","<<this->upRight.y<<"), ";
	cout<<"UL("<<this->upLeft.x<<","<<this->upLeft.y<<"), ";
	cout<<"DR("<<this->downRight.x<<","<<this->downRight.y<<"), ";
	cout<<"DL("<<this->downLeft.x<<","<<this->downLeft.y<<")"<<endl;
}

RandomDynamicConflictGraphGeneratorParameters::RandomDynamicConflictGraphGeneratorParameters(int areaWidth, int areaHeight,
		int numberOfPoints, int positionsPerPoint, int labelWidth, int labelHeight) {
	this->areaWidth = areaWidth;
	this->areaHeight = areaHeight;
	this->numberOfPoints = numberOfPoints;
	this->labelWidth = labelWidth;
	this->labelHeight = labelHeight;
	this->positionsPerPoint = positionsPerPoint;
	this->meanSpd = 0;
	this->stdSpd = 0;
	this->maxSpd = 0;
	this->meanDir = 0;
	this->stdDir = 0;
	this->maxDir = 0;
}

RandomDynamicConflictGraphGeneratorParameters::~RandomDynamicConflictGraphGeneratorParameters() {
}

RandomDynamicConflictGraphGenerator::RandomDynamicConflictGraphGenerator(RandomDynamicConflictGraphGeneratorParameters parameters){
	this->labelWidth = parameters.labelWidth;
	this->labelHeight = parameters.labelHeight;
	this->positionsPerPoint = parameters.positionsPerPoint;
	if(parameters.mpParamsIndicator) {
		MovingPointsGenerator::Parameters mpParams;
		mpParams.areaHeight = parameters.areaHeight;
		mpParams.areaWidth = parameters.areaWidth;
		mpParams.numPoints = parameters.numberOfPoints;
		mpParams.meanSpd = parameters.meanSpd;
		mpParams.stdSpd = parameters.stdSpd;
		mpParams.maxSpd = parameters.maxSpd;
		mpParams.meanDir = parameters.meanDir;
		mpParams.stdDir = parameters.stdDir;
		mpParams.maxDir = parameters.maxDir;
		mpGen = new MovingPointsGenerator(mpParams);
	} else {
		mpGen = new MovingPointsGenerator(parameters.areaWidth, parameters.areaHeight, parameters.numberOfPoints);
	}
	movingPoints = NULL;
}

RandomDynamicConflictGraphGenerator::~RandomDynamicConflictGraphGenerator() {
	if(mpGen!= NULL) {
		delete mpGen;
	}
	if(movingPoints != NULL) {
		for(int j=0; j<movingPoints->size(); j++) {
			delete movingPoints->at(j);
		}
	}
}

labelplacement::ConflictGraph* RandomDynamicConflictGraphGenerator::generate(int forwardStates) {
	for(int i=0; i<forwardStates; i++) {
		if(movingPoints == NULL) {
			movingPoints = mpGen->generate();
		} else {
			for(int j=0; j<movingPoints->size(); j++)
			{
				MovingPoint* mp = movingPoints->at(j);
				mp->updatePosition();
			}
		}
	}
	return generateConflictGraph(movingPoints);
}

labelplacement::ConflictGraph* RandomDynamicConflictGraphGenerator::generateConflictGraph(vector<MovingPoint*>* movingPoints) {
	int numberOfPoints = movingPoints->size();
	int numberOfPositions = numberOfPoints*positionsPerPoint;

	vector<int>* adjacencyListOfPositions = new vector<int>[numberOfPositions];

	int pointIx1=0;
	for(vector<MovingPoint*>::iterator it = movingPoints->begin(); it!=movingPoints->end(); it++) {
		MovingPoint* mp1 = *it;
		Rectangle positionsMp1[4];
		for(int i=0; i<4; i++) {
			int xCoef = i&0x02;
			int yCoef = i&0x01;
			Position upRightPos(mp1->position.x + xCoef*labelWidth,mp1->position.y + yCoef*labelHeight);
			positionsMp1[i] = Rectangle(upRightPos, labelWidth, labelHeight);
		}
		int pointIx2=0;
		for(vector<MovingPoint*>::iterator it2 = movingPoints->begin(); it2!=movingPoints->end(); it2++) {
			if(pointIx1 != pointIx2) {
				MovingPoint* mp2 = *it2;
				Rectangle positionsMp2[4];
				for(int i=0; i<4; i++) {
					int xCoef = i&0x02;
					int yCoef = i&0x01;
					positionsMp2[i] = Rectangle(Position(mp2->position.x + xCoef*labelWidth,mp2->position.y + yCoef*labelHeight), labelWidth, labelHeight);
				}
				for(int i=0; i<4; i++) {
					int positionIx1 = pointIx1*positionsPerPoint + i;
					for(int j=0; j<4; j++) {
						int positionIx2 = pointIx2*positionsPerPoint + j;

						//std::cout<<positionIx1<<" "<<positionIx2<<endl;
						if(intersectionOccurs(positionsMp1[i], positionsMp2[j])) {
							//std::cout<<"intersection"<<endl;
							adjacencyListOfPositions[positionIx1].push_back(positionIx2);
						}
					}
				}
			}
			pointIx2++;
		}
		pointIx1++;
	}
	labelplacement::ConflictGraph* conflictGraph = new labelplacement::ConflictGraph();
	conflictGraph->initializeOf(positionsPerPoint, numberOfPositions, adjacencyListOfPositions);
	return conflictGraph;
}

bool RandomDynamicConflictGraphGenerator::intersectionOccurs(Rectangle& r1, Rectangle& r2) {
	int r1XUpper = r1.upRight.x;
	int r1XLower = r1.upLeft.x;
	int r1YUpper = r1.upRight.y;
	int r1YLower = r1.downRight.y;
	int r2XUpper = r2.upRight.x;
	int r2XLower = r2.upLeft.x;
	int r2YUpper = r2.upRight.y;
	int r2YLower = r2.downRight.y;
	//cout<<r1XUpper<<","<<r2XLower<<","<<r1XLower<<","<<r2XUpper<<","<<r1YUpper<<","<<r2YLower<<","<<r1YLower<<","<<r2YUpper<<endl;
	if(r1XUpper >= r2XLower && r1XLower <= r2XUpper
			&& r1YUpper >= r2YLower && r1YLower <= r2YUpper) {
		return true;
	}
	return false;
}

