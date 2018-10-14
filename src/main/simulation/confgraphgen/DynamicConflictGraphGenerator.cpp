/*
 * DynamicConflictGraphGenerator.cpp
 *
 *  Created on: May 26, 2018
 *      Author: oakile
 */

#include "DynamicConflictGraphGenerator.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

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

DynamicConflictGraphGeneratorParameters::DynamicConflictGraphGeneratorParameters(int areaWidth, int areaHeight,
		int numberOfPoints, int positionsPerPoint, int labelWidth, int labelHeight) {
	this->areaWidth = areaWidth;
	this->areaHeight = areaHeight;
	this->numberOfPoints = numberOfPoints;
	this->labelWidth = labelWidth;
	this->labelHeight = labelHeight;
	this->positionsPerPoint = positionsPerPoint;
}
DynamicConflictGraphGenerator::DynamicConflictGraphGenerator() {
	mpGen = MovingPointsGenerator();
	labelWidth = 0;
	labelHeight = 0;
	positionsPerPoint = 0;
	movingPoints = NULL;
}
DynamicConflictGraphGenerator::DynamicConflictGraphGenerator(DynamicConflictGraphGeneratorParameters parameters) {
	mpGen = MovingPointsGenerator(parameters.areaWidth, parameters.areaHeight, parameters.numberOfPoints);
	this->labelWidth = parameters.labelWidth;
	this->labelHeight = parameters.labelHeight;
	this->positionsPerPoint = parameters.positionsPerPoint;
	movingPoints = NULL;
}

DynamicConflictGraphGenerator::~DynamicConflictGraphGenerator() {
	if(movingPoints != NULL) {
		for(int j=0; j<movingPoints->size(); j++) {
			delete movingPoints->at(j);
		}
	}
}

labelplacement::ConflictGraph* DynamicConflictGraphGenerator::generate(int forwardStates) {
	for(int i=0; i<forwardStates; i++) {
		if(movingPoints == NULL) {
			movingPoints = mpGen.generate();
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

labelplacement::ConflictGraph* DynamicConflictGraphGenerator::generateConflictGraph(vector<MovingPoint*>* movingPoints) {
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

void DynamicConflictGraphGenerator::writeConflictGraphIntoFile(std::vector<int>* adjacencyListOfPositions, int pointSize, string filename) {
	int positionSize = pointSize * positionsPerPoint;
	ofstream ofs(filename, ios_base::out);
	ofs<<endl;
	ofs<<pointSize<<endl;
	ofs<<positionsPerPoint<<endl;
	for(int i=0; i<positionSize; i++) {
		vector<int>* adjacencyListOfPosition = adjacencyListOfPositions + i;
		ofs<<adjacencyListOfPosition->size()<<endl;
		for(vector<int>::iterator it=adjacencyListOfPosition->begin(); it!=adjacencyListOfPosition->end(); it++) {
			ofs<<*it<<" ";
		}
		ofs<<endl;
	}
	ofs.close();
}
bool DynamicConflictGraphGenerator::intersectionOccurs(Rectangle& r1, Rectangle& r2) {
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


int DynamicConflictGraphGenerator::compareConflictGraphs(int size, vector<int>* cg1, vector<int>* cg2) {
	int diff=0;
	int conflictSize1=0;
	int conflictSize2=0;
	for(int i=0; i<size; i++) {
		vector<int> adjList1 = vector<int>(*(cg1 + i));
		vector<int> adjList2 = vector<int>(*(cg2 + i));
		conflictSize1+=adjList1.size();
		conflictSize2+=adjList2.size();
		for(int j=0; j<adjList1.size();) {
			bool isFound = false;
			int k=0;
			for(; k < adjList2.size(); k++) {
				if(adjList1[j] == adjList2[k]) {
					isFound = true;
					break;
				}
			}
			if(isFound) {
				adjList1.erase(adjList1.begin()+j);
				adjList2.erase(adjList2.begin()+k);
			}
			else
			{
				j++;
			}
		}
		diff+= adjList1.size() + adjList2.size();
	}
	diff/=2;
	conflictSize1/=2;
	conflictSize2/=2;
	int diffConflictSize = (int)abs(conflictSize1-conflictSize2);
	cout<<conflictSize1<<" "<<conflictSize2<<" "<<diffConflictSize<<endl;
	return diff;
}

/*int main(int argc, char** argv)
{
	int numMovements = 25;
	int numPoints[5] = {100, 250, 500, 750, 1000};

	DynamicConflictGraphGenerator cggenerator;
	cggenerator.positionsPerPoint = 4;
	cggenerator.labelWidth = 20;
	cggenerator.labelHeight = 10;
	for(int h=0; h<5; h++) {
		MovingPointsGenerator generator(800, 600, numPoints[h]);
		std::vector<MovingPoint*>* movingPoints = generator.generate();
		string filePath = "output/d" + to_string(numPoints[h]) + "/";
		vector<int>* adjacencyListOfPositions = NULL;
		for(int i=0; i<numMovements; i++)
		{
			//cout<<i<<endl;
			//for(vector<MovingPoint*>::iterator it=movingPoints->begin();it!=movingPoints->end(); it++){
			//	MovingPoint* mp = *it;
			//	cout<<"("<<mp->position.x<<" "<<mp->position.y<<") ";
			//}
			//cout<<endl;
			vector<int>* adjacencyListOfPositionsNew = cggenerator.generateConflictGraph(movingPoints);
			int confGraphDiff = 0;
			if(adjacencyListOfPositions != NULL) {
				cout<<"Num Points: "<<numPoints[h]<<" Mvt: "<<i<<" Diff: "<<cggenerator.compareConflictGraphs(cggenerator.positionsPerPoint*numPoints[h] , adjacencyListOfPositions, adjacencyListOfPositionsNew)<<endl;
			}
			adjacencyListOfPositions = adjacencyListOfPositionsNew;
	//		for(int j=0; j<numPoints*4; j++) {
	//			cout<<j+1<<" : ";
	//			for(int k = 0; k<adjacencyListOfPositions[j].size(); k++) {
	//				cout<<adjacencyListOfPositions[j][k]<<" ";
	//			}
	//			cout<<endl;
	//		}
			//cout<<endl;

			string filename = "d" + to_string(numPoints[h]) +"_" + (i+1<10?"0":"") + to_string(i+1) + string(".dat");
			cggenerator.writeConflictGraphIntoFile(adjacencyListOfPositions, numPoints[h], filePath + filename);
			for(int j=0; j<numPoints[h]; j++)
			{
				MovingPoint* mp = movingPoints->at(j);
				mp->updatePosition();
			}
		}
	}
	return 0;
}*/
