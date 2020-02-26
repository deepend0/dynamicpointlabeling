/*
 * FHPLPOptimizerAdaptor.h
 *
 *  Created on: Jun 7, 2019
 *      Author: oakile
 */

#ifndef _FHPLP_FHPLPOPTIMIZERADAPTOR_H_
#define _FHPLP_FHPLPOPTIMIZERADAPTOR_H_

#include "ConflictGraph.h"
#include "Solution.h"
#include "fhplp/FastHeuristicPLP.hpp"
#include "fhplp/ConflictGraphElements.hpp"
#include "fhplp/Solution.hpp"

namespace fhplp {

class FHPLPOptimizerAdaptor {
public:
	FastHeuristicPLP* fhplp;
	FHPLPOptimizerAdaptor(FastHeuristicPLP* fhplp);
	virtual ~FHPLPOptimizerAdaptor();
	labelplacement::Solution& optimize(labelplacement::ConflictGraph& conflictGraph);
	fhplp::ConflictGraph* convertConflictGraph(labelplacement::ConflictGraph& conflictGraph);
	labelplacement::Solution& convertSolution(fhplp::Solution* solution, int numPositionsPerPoint);
};

} /* namespace fhplp */

#endif /* _FHPLP_FHPLPOPTIMIZERADAPTOR_H_ */
