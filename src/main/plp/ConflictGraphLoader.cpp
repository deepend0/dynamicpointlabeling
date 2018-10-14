/*
 * ConflictGraphLoader.cpp
 *
 *  Created on: Dec 16, 2017
 *      Author: oakile
 */


#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/regex_find_format.hpp>

#include "ConflictGraph.h"
#include "ConflictGraphLoader.h"
#include "Graph.h"

using namespace boost;
using namespace std;

namespace labelplacement {

ConflictGraphLoader::ConflictGraphLoader() {
	// TODO Auto-generated constructor stub
}

ConflictGraphLoader::~ConflictGraphLoader() {
	// TODO Auto-generated destructor stub
}

ConflictGraph* ConflictGraphLoader::load(string filePath) {
	ifstream ifs;
	ifs.open(filePath.c_str(), ifstream::in);
	if (ifs.fail()) {
		cout << "Opening instance data file is failed." << endl;
		return NULL;
	}

	string line;
	istringstream iss;

	getline(ifs, line); //Skip first line as it is empty.

	getline(ifs, line); //Line of point feature count
	iss.str(line);
	getline(iss, line, ' ');
	int pointNumber = atoi(line.c_str());

	getline(ifs, line); //Line of position count
	iss.str(line);
	getline(iss, line, ' ');
	int positionNumber = atoi(line.c_str());
	const int totalPositionCount = pointNumber * positionNumber;

	vector<int>* adjacencyListOfPositions = new vector<int> [totalPositionCount];

	for (int i = 0; i < totalPositionCount; i++) {
		getline(ifs, line); //Line of potential conflict count
		iss.str(line);
		getline(iss, line, ' ');
		int positionConflictCount = atoi(line.c_str());
		//int conflictSize = positionConflictCount - positionNumber + 1;

		getline(ifs, line); //Line of potential conflict label numbers

		typedef find_iterator<string::iterator> string_find_iterator;
		regex numRegex("[0-9]+");
		string_find_iterator sfit = make_find_iterator(line,
				regex_finder(numRegex));
		int currentPointIx = i / positionNumber;
		for (int j = 0; j < positionConflictCount; j++) {
			string positionNoStr;
			positionNoStr = copy_range<string>(*sfit++);
			int positionNo = atoi(positionNoStr.c_str());
			int positionIx = positionNo - 1;

			int pointIx = positionIx / positionNumber;
			if (currentPointIx != pointIx) {
				adjacencyListOfPositions[i].push_back(positionIx);
			}
		}
	}
	ifs.close();
	ConflictGraph* conflictGraph = new ConflictGraph();
	conflictGraph->initializeOf(positionNumber, totalPositionCount, adjacencyListOfPositions);
	return conflictGraph;
}
}
