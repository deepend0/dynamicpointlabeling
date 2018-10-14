#ifndef _MOVINGPOINTSGENERATOR_HPP_
#define _MOVINGPOINTSGENERATOR_HPP_

#include <cmath>
#include <vector>

class Position
{
	public:
		int x;
		int y;
		Position(int x, int y);
};

class MovingPoint;

class MovingPointsGenerator
{
	public:
		static const double pi;
		MovingPointsGenerator();
		MovingPointsGenerator(int areaWidth, int areaHeight, int numPoints);
		
		int areaWidth;
		int areaHeight;
		int numPoints;
		
		double maxSpd;
		double meanSpd;
		double stdSpd;
		double maxDir;
		double meanDir;
		double stdDir;
		
		std::vector<MovingPoint*>* generate();
};

class MovingPoint
{
	public:
	MovingPoint();
	
	MovingPointsGenerator* parent;
	Position position;
	double speed;
	double direction;
	
	double maxAcc;
	double meanAcc;
	double stdAcc;
	double accFreq;
	double meanDirChg;
	double stdDirChg;
	double dirChgFreq;
	
	void updatePosition();
};
#endif




