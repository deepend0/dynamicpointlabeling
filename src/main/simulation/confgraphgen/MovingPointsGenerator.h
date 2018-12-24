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
		struct ParametersCoef {
			int areaWidth;
			int areaHeight;
			int numPoints;
			double meanSpdCoef;
			double stdSpdCoef;
			double maxSpdCoef;
			double meanDirDivergence;
			double stdDirCoef;
			double maxDirDivergence;
		};

		struct Parameters {
			int areaWidth;
			int areaHeight;
			int numPoints;
			double meanSpd;
			double stdSpd;
			double maxSpd;
			double meanDir;
			double stdDir;
			double maxDir;
		};

		static const double pi;
		MovingPointsGenerator(int areaWidth, int areaHeight, int numPoints);
		MovingPointsGenerator(ParametersCoef& parameters);
		MovingPointsGenerator(Parameters& parameters);
		
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
	private:
		void calculateMovementParameters(ParametersCoef parameters);
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




