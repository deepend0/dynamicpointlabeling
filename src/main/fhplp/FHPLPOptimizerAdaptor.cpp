/*
 * FHPLPOptimizerAdaptor.cpp
 *
 *  Created on: Jun 7, 2019
 *      Author: oakile
 */

#include "FHPLPOptimizerAdaptor.h"

#include <vector>

namespace fhplp {

FHPLPOptimizerAdaptor::FHPLPOptimizerAdaptor(FastHeuristicPLP* fhplp) {
	this->fhplp = fhplp;
}

FHPLPOptimizerAdaptor::~FHPLPOptimizerAdaptor() {
}


labelplacement::Solution& FHPLPOptimizerAdaptor::optimize(labelplacement::ConflictGraph& conflictGraph) {
	fhplp::ConflictGraph* fhplpConflictGraph = convertConflictGraph(conflictGraph);
	fhplp::Solution* fhplpSolution = fhplp->optimize(fhplpConflictGraph);
	int numPoints = conflictGraph.getConflictGraphOfPoints()->getVertexNumber();
	int numPositions = conflictGraph.getConflictGraphOfPositions()->getVertexNumber();
	int numPositionsPerPoint = numPositions / numPoints;
	labelplacement::Solution& solution = convertSolution(fhplpSolution, numPositionsPerPoint);
	return solution;
}

fhplp::ConflictGraph* FHPLPOptimizerAdaptor::convertConflictGraph(labelplacement::ConflictGraph& conflictGraph) {
	int numPoints = conflictGraph.getConflictGraphOfPoints()->getVertexNumber();
	int numPositions = conflictGraph.getConflictGraphOfPositions()->getVertexNumber();
	int numPositionsPerPoint = numPositions / numPoints;

	fhplp::PointFeature::numPositions = numPositionsPerPoint;

	fhplp::ConflictGraph* fhplpConflictGraph = new fhplp::ConflictGraph();
	std::vector<fhplp::PointFeature*>* pointFeatures = new std::vector<PointFeature*>();
	fhplpConflictGraph->setPointFeatures(pointFeatures);
	fhplpConflictGraph->setNumPositions(numPositionsPerPoint);

	fhplp::PositionNode** allPositions = new fhplp::PositionNode*[numPositions];
	for(int i=0; i<numPositions; i++) {
		allPositions[i] = new fhplp::PositionNode();
	}

	std::vector<int>* adjacencyListOfPositions = conflictGraph.getConflictGraphOfPositions()->getAdjacencyList();
	for(int i=0; i<numPoints; i++) {
		fhplp::PointFeature* pointFeature = new fhplp::PointFeature(i+1);
		pointFeatures->push_back(pointFeature);
		std::vector<fhplp::PositionNode*>* positions = new std::vector<fhplp::PositionNode*>();
		pointFeature->setPositions(positions);
		for(int j=0; j<numPositionsPerPoint; j++) {
			int positionIndex = i*numPositionsPerPoint + j;
			fhplp::PositionNode* position = allPositions[positionIndex];
			positions->push_back(position);
			position->setNo(positionIndex+1);
			position->setPointFeature(pointFeature);
			std::vector<fhplp::PositionNode*>* conflictingPositions = new std::vector<fhplp::PositionNode*>();
			position->setConflictingPositions(conflictingPositions);

			std::vector<int>* adjacencyListOfPosition = adjacencyListOfPositions + positionIndex;
			for(int k=0; k<adjacencyListOfPosition->size(); k++) {
				int conflictingPositionIx = adjacencyListOfPosition->at(k);
				conflictingPositions->push_back(allPositions[conflictingPositionIx]);
			}
		}
	}
	delete [] allPositions;
	return fhplpConflictGraph;
}

labelplacement::Solution& FHPLPOptimizerAdaptor::convertSolution(fhplp::Solution* fhplpSolution, int numPositionsPerPoint) {
	labelplacement::Solution* solution = new labelplacement::Solution(NULL);
	labelplacement::Solution& solutionRef = *solution;
	int numPositions = fhplpSolution->getPositions()->size();
	int* labelPlacements = new int[numPositions];
	solutionRef.setLabelPlacements(labelPlacements);
	for(int i=0; i<numPositions; i++) {
		labelPlacements[i] = (fhplpSolution->getPositions()->at(i)->getNo()-1) % numPositionsPerPoint;
	}
	delete fhplpSolution;
	return solutionRef;
}

} /* namespace fhplp */
