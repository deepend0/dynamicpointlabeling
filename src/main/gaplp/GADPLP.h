/*
 * GADPLP.h
 *
 *  Created on: Aug 20, 2018
 *      Author: oakile
 */

#ifndef _GADPLP_GADPLP_H_
#define _GADPLP_GADPLP_H_

#include <ga/GA1DArrayGenome.h>

#include "ConflictGraph.h"
#include "Solution.h"
#include "RandomImmigrantSGMXGA.h"

namespace GADPLP {

struct GADPLPParameters {
	int numPoints;
	int numPositionsPerPoint;
	int populationSize;
	int selector;
	bool elitist;
	int period;
	int numInitGen;
	int numImmigrateGen;
	double immigrationRate;
	int crossoverType;
	double crossoverRate;
	double mutationRate;
	double groupProportion;
	double groupProportionMargin;
	double individualProportion;
	double individualProportionMargin;
	double individualScoreLimit;
	GADPLPParameters();
	GADPLPParameters(int numPoints, int numPositionsPerPoint, int populationSize, double crossoverRate, double mutationRate, int selector,
			bool elitist, int period, int numInitGen, int numImmigrateGen, double immigrationRate,
			int crossoverType, double groupProportion, double groupProportionMargin,
			double individualProportion, double individualProportionMargin, double individualScoreLimit);
};

class GADPLP {
public:
	struct Objective {
		Objective();
		labelplacement::ConflictGraph** conflictGraph;
		float operator()(GAGenome & g);
	};

	struct AlleleEvaluator {
		AlleleEvaluator();
		labelplacement::ConflictGraph* conflictGraph;
		float operator()(const GA1DArrayGenome<int> &, int);
	};

	GADPLP();
	virtual ~GADPLP();
	void init(GADPLPParameters& parameters);
	labelplacement::Solution& optimize(labelplacement::ConflictGraph& conflictGraph);
private:
	Objective objective;
	gmxga::RandomImmigrantSGMXGA* ga;
	GADPLPParameters parameters;
	void updateProblem(labelplacement::ConflictGraph& conflictGraph);
	labelplacement::Solution& buildSolution(labelplacement::ConflictGraph& conflictGraph, const GA1DArrayAlleleGenome<unsigned int> & bestIndividual);
};

} /* namespace GADPLP */

#endif /* _GADPLP_GADPLP_H_ */
