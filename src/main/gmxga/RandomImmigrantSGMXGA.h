/*
 * RandomImmigrantSGMXGA.h
 *
 *  Created on: Aug 19, 2018
 *      Author: oakile
 */

#ifndef _RANDOMIMMIGRANTSGMXGA_H_
#define _RANDOMIMMIGRANTSGMXGA_H_

#include "gmxga/SimpleGMXGA.h"

namespace gmxga {

class RandomImmigrantSGMXGA: public SimpleGMXGA {
public:
	RandomImmigrantSGMXGA(const GAGenome& gaGenome);
	virtual ~RandomImmigrantSGMXGA();

	int getPeriod() const;

	void setPeriod(int period);

	int getNumImmigrateGen() const;

	void setNumImmigrateGen(int numImmigrateGen);

	int getNumInitGen() const;

	void setNumInitGen(int numInitGen);

	double getImmigrationRate() const;

	void setImmigrationRate(double immigrationRate);

	int nGenerations(unsigned int n);

	virtual void initialize(unsigned int seed=0);

	virtual void evolve(unsigned int seed=0);

private:
	int period;
	int numInitGen;
	int numImmigrateGen;
	double immigrationRate;

	int curStage;
};

} /* namespace gmxga */

#endif /* _RANDOMIMMIGRANTSGMXGA_H_ */
