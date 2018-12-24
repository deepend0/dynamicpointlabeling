/*
 * PresentFileSetDCGG.cpp
 *
 *  Created on: Nov 4, 2018
 *      Author: oakile
 */

#include "PresentFileSetDCGG.h"

#include <iostream>

#include "ConflictGraphLoader.h"

PresentFileSetDCGG::Parameters::Parameters(std::string path) {
	this->path = path;
}

PresentFileSetDCGG::PresentFileSetDCGG(Parameters parameters) {
	this->path = parameters.path;
	currentState = 0;

}

PresentFileSetDCGG::~PresentFileSetDCGG() {
}

labelplacement::ConflictGraph* PresentFileSetDCGG::generate(int forwardStates) {
	//std::cout<<"Current State: "<<currentState<<std::endl;
	labelplacement::ConflictGraphLoader cgl;
	currentState+=forwardStates;
	std::string filename = (currentState < 10 ? "0" : "") + std::to_string(currentState) + ".dat";
	std::string instanceFullPath = path + filename;
	labelplacement::ConflictGraph* cg = cgl.load(instanceFullPath);
	return cg;
}
