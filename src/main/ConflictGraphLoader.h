/*
 * ConflictGraphLoader.h
 *
 *  Created on: Dec 16, 2017
 *      Author: oakile
 */

#ifndef _CONFLICTGRAPHLOADER_H_
#define _CONFLICTGRAPHLOADER_H_

#include "ConflictGraph.h"
#include <string>

namespace labelplacement {

class ConflictGraphLoader {
public:
	ConflictGraphLoader();
	virtual ~ConflictGraphLoader();
	ConflictGraph* load(std::string filePath);
};
}
#endif /* _CONFLICTGRAPHLOADER_H_ */
