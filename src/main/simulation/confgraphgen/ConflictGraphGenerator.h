/*
 * ConflictGraphGenerator.h
 *
 *  Created on: Nov 4, 2018
 *      Author: oakile
 */

#ifndef SIMULATION_CONFGRAPHGEN_CONFLICTGRAPHGENERATOR_H_
#define SIMULATION_CONFGRAPHGEN_CONFLICTGRAPHGENERATOR_H_

class ConflictGraphGenerator {
public:
	virtual ~ConflictGraphGenerator(){};
	virtual labelplacement::ConflictGraph* generate(int forwardStates)=0;
};


#endif /* SIMULATION_CONFGRAPHGEN_CONFLICTGRAPHGENERATOR_H_ */
