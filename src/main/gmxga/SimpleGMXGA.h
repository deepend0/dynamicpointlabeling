/*
 * SimpleGMXGA.h
 *
 *  Created on: Aug 13, 2018
 *      Author: oakile
 */

#ifndef _GMXGA_SimpleGMXGA_H_
#define _GMXGA_SimpleGMXGA_H_

#include "GroupMaskingCrossover.h"

#include <ga/GA1DArrayGenome.h>
#include <ga/GAGenome.h>
#include <ga/GASimpleGA.h>

namespace gmxga {
class SimpleGMXGA : public GASimpleGA{
public:
	SimpleGMXGA(const GAGenome& gaGenome);
	virtual ~SimpleGMXGA();
	GAGenome::SexualCrossoverFn crossover(GAGenome::SexualCrossoverFn f);
	GroupMaskingCrossover* getCrossover();
	void setCrossover(GroupMaskingCrossover* gmx);

private:
	GroupMaskingCrossover* gmx;
};
}
#endif /* _GMXGA_SimpleGMXGA_H_ */
