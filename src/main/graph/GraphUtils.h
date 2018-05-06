/*
 * GraphUtils.h
 *
 *  Created on: Feb 17, 2018
 *      Author: oakile
 */

#ifndef SRC_GRAPHUTILS_H_
#define SRC_GRAPHUTILS_H_

#include "Graph.h"

namespace graph {

class GraphUtils {
public:
	GraphUtils();
	virtual ~GraphUtils();
	static void printGraph(Graph* graph);
};

} /* namespace labelplacement */

#endif /* _GRAPHUTILS_H_ */
