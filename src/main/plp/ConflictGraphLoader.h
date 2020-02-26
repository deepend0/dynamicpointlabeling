/*
 * ConflictGraphLoader.h
 *
 *  Created on: Dec 16, 2017
 *      Author: oakile
 */

#ifndef _CONFLICTGRAPHLOADER_H_
#define _CONFLICTGRAPHLOADER_H_

#include <string>
#include "ConflictGraph.h"

namespace labelplacement {

class ConflictGraphLoader {
public:
	ConflictGraphLoader();
	virtual ~ConflictGraphLoader();
	ConflictGraph* load(std::string filePath);
};
}
#endif /* _CONFLICTGRAPHLOADER_H_ */
