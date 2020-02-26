/*
 * GroupMaskingCrossover.h
 *
 *  Created on: Aug 6, 2018
 *      Author: oakile
 */

#ifndef _GROUPMASKINGCROSSOVER_H_
#define _GROUPMASKINGCROSSOVER_H_

#include "Graph.h"
#include "AlleleGroupers.h"

#include <ga/GA1DArrayGenome.h>
#include <ga/GAGenome.h>
#include <ga/GASimpleGA.h>

#include <functional>
#include <vector>

class GroupMaskingCrossover {
public:

	GroupMaskingCrossover();
	virtual ~GroupMaskingCrossover();
	virtual int operator()(const GAGenome&, const GAGenome&, GAGenome*, GAGenome*)=0;
	virtual void initialize(graph::Graph* graph);
	virtual int* selectGroupsOfIndividual(const GA1DArrayGenome<int>& individual)=0;
	void setGrouper(std::function<AlleleGrouper> grouper);

protected:
	std::function<AlleleGrouper> grouper;
	GroupsT groups;
};

struct GMXDelegate {
	GroupMaskingCrossover* gmx;
	GMXDelegate(GroupMaskingCrossover* gmx);
	int operator()(const GAGenome& p1, const GAGenome& p2, GAGenome* c1, GAGenome* c2);
};

class UniformRankedGroupMaskingCrossover : public GroupMaskingCrossover {
public:
	typedef float AlleleEvaluator(const GA1DArrayGenome<int> &, int);

	UniformRankedGroupMaskingCrossover();
	~UniformRankedGroupMaskingCrossover();
	void setParams(double groupProportion, double groupProportionMargin,
			double individualProportion, double individualProportionMargin, double individualScoreLimit);
	std::function<AlleleEvaluator>* getEvaluator();
	void setEvaluator(std::function<AlleleEvaluator>* evaluator);
	void setGA(GASimpleGA* ga);
	int operator()(const GAGenome&, const GAGenome&, GAGenome*, GAGenome*);
	void initialize(graph::Graph* graph);
	int* selectGroupsOfIndividual(const GA1DArrayGenome<int>& individual);

	struct ComparePairIntDoubleMin {
		bool operator()(std::pair<int, double> p1, std::pair<int, double> p2) {
			return p1.second < p2.second;
		}
	};
	struct ComparePairIntDoubleMax {
		bool operator()(std::pair<int, double> p1, std::pair<int, double> p2) {
			return p1.second > p2.second;
		}
	};
	
private:
	//EXTERNAL PARAMETERS
	double groupProportion;
	double groupProportionMargin;
	double individualProportion;
	double individualProportionMargin;
	double individualScoreLimit;
	//INTERNAL VARIABLES
	int alleleNumber;
	int numberOfGroups;
	int curNumberOfGroups;

	GASimpleGA* ga;
	std::function<AlleleEvaluator>* evaluator;
};

#endif /* _GROUPMASKINGCROSSOVER_H_ */
