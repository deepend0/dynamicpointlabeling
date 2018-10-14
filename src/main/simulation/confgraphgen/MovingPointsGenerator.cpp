#include "MovingPointsGenerator.h"

#include <cmath>
#include <random>
#include <vector>

#include <iostream>

const double MovingPointsGenerator::pi = acos(-1);

Position::Position(int x, int y)
{
	this->x = x;
	this->y = y;
}

MovingPointsGenerator::MovingPointsGenerator()
{
}

MovingPointsGenerator::MovingPointsGenerator(int areaWidth, int areaHeight, int numPoints)
{
	this->areaWidth = areaWidth;
	this->areaHeight = areaHeight;
	this->numPoints = numPoints;
}

std::vector<MovingPoint*>* MovingPointsGenerator::generate()
{

	double diagonLength = sqrt(areaWidth*areaWidth + areaHeight*areaHeight);
	meanSpd = diagonLength/100;
	stdSpd = meanSpd/100;
	maxSpd = meanSpd*2;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> freqDist(0.0, 1.0);
	std::uniform_real_distribution<double> cosDist(-1.0, 1.0);

	meanDir = acos((double)cosDist(gen));
	stdDir = meanDir/200;
	maxDir = 2*pi;

	std::uniform_int_distribution<> posWidthDist(0, areaWidth);
	std::uniform_int_distribution<> posHeightDist(0, areaHeight);
	std::normal_distribution<double> spdDist(meanSpd, stdSpd);
	std::normal_distribution<double> dirDist(meanDir, stdDir);

	std::vector<MovingPoint*>* movingPoints = new std::vector<MovingPoint*>();

	for(int i=0; i<numPoints; i++)
	{
		MovingPoint* mp = new MovingPoint();
		mp->parent = this;
		mp->position = Position(posWidthDist(gen),posHeightDist(gen));
		do
		{
			mp->speed = spdDist(gen);
		}
		while(mp->speed > maxSpd && mp->speed < 0);
		mp->accFreq = freqDist(gen);

		do
		{
			mp->direction = dirDist(gen);
		}
		while(mp->direction > maxDir && mp->direction < 0);
		mp->dirChgFreq = freqDist(gen);

		movingPoints->push_back(mp);
	}

	return movingPoints;
}

MovingPoint::MovingPoint() : position(0,0)
{
}

void MovingPoint::updatePosition()
{
	meanAcc = 0;
	stdAcc = speed/100;
	maxAcc = speed/20;
	meanDirChg = 0;
	stdDirChg = MovingPointsGenerator::pi/360;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::normal_distribution<double> accDist(meanAcc, stdAcc);
	std::normal_distribution<double> dirChgDist(meanDirChg, stdDirChg);
	std::uniform_real_distribution<double> freqDist(0.0, 1.0);

	double accFreq = freqDist(gen);
	double acc = 0;

	if(this->accFreq >= accFreq)
	{
		do
		{
			acc = accDist(gen);
		}
		while(acc>maxAcc);
	}
	else
	{
		//std::cout<<"Acc 0"<<std::endl;
	}

	double dirChgFreq = freqDist(gen);
	double dirChg = 0;
	if(this->dirChgFreq >= dirChgFreq) {
		do
		{
			dirChg=dirChgDist(gen);
		}
		while(dirChg>MovingPointsGenerator::pi*2 && dirChg<0);
	}
	else
	{
		//std::cout<<"Dir Chg 0"<<std::endl;
	}
	//acc = 0;
	//dirChg = 0;
	speed += acc;
	direction = direction + dirChg;
	
	if(speed > parent->maxSpd)
		speed = parent->maxSpd;
	else if(speed < 0)
		speed = 0;
	/*if(direction > parent->maxDirection)
		direction = parent->maxDirection;
	else if(direction < 0)
		direction = 0;*/
	
	position.x += speed*sin(direction);
	position.y += speed*cos(direction);
	
	if(position.x < 0) {
		position.x = (-1*position.x/parent->areaWidth + 1) * parent->areaWidth;
	}
	else if(position.x > parent->areaWidth) {
		position.x %= parent->areaWidth;
	}

	if(position.y < 0) {
			position.y = (-1*position.y/parent->areaHeight + 1) * parent->areaHeight;
	}
	else if(position.y > parent->areaHeight) {
		position.y %= parent->areaHeight;
	}
	//std::cout<<"SPD: "<<speed<<" DIR: "<<direction<<std::endl;
	//std::cout<<"HPOS: "<<position.x<<" VPOS: "<<position.y;
	
	//std::cout<<std::endl;
}
