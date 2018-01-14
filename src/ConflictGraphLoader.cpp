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

#include "ConflictGraphLoader.h"
#include "ConflictGraph.h"

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
	ConflictGraph* conflictGraph = new ConflictGraph();
	ifstream ifs;
	ifs.open(filePath.c_str(), ifstream::in);
	if(ifs.fail())
	{
		cout<<"Opening instance data file is failed."<<endl;
		return NULL;
	}

	string line;
	istringstream iss;

	getline(ifs, line); //Skip first line as it is empty.

	getline(ifs, line); //Line of point feature count
	iss.str(line);
	getline(iss, line, ' ');
	int pointNumber = atoi(line.c_str());
	conflictGraph->setPointNumber(pointNumber);

	getline(ifs, line); //Line of position count
	iss.str(line);
	getline(iss, line, ' ');
	int positionNumber = atoi(line.c_str());
	conflictGraph->setPositionNumber(positionNumber);

	const int totalPositionCount = pointNumber * positionNumber;

	bool** adjacencyMatrix = new bool*[totalPositionCount];
	for(int i=0; i<totalPositionCount; i++)
	{
		adjacencyMatrix[i] = new bool[totalPositionCount];
		for(int j=0; j<totalPositionCount; j++)
		{
			adjacencyMatrix[i][j] = false;
		}
	}

	vector<int>* conflictingPositions = new vector<int>[totalPositionCount];
	vector<int>* conflictingPoints = new vector<int>[pointNumber];

	conflictGraph->setAdjacencyMatrix(adjacencyMatrix);
	conflictGraph->setConflictingPositions(conflictingPositions);
	conflictGraph->setConflictingPoints(conflictingPoints);

	for(int i=0; i<totalPositionCount; i++)
	{
		getline(ifs, line); //Line of potential conflict count
		iss.str(line);
		getline(iss, line, ' ');
		int positionConflictCount = atoi(line.c_str());
		int conflictSize = positionConflictCount-positionNumber+1;

		getline(ifs, line); //Line of potential conflict label numbers

		typedef find_iterator<string::iterator> string_find_iterator;
		regex numRegex("[0-9]+");
		string_find_iterator sfit=make_find_iterator(line, regex_finder(numRegex));
		int currentPointIx = i/positionNumber;
		for(int j=0; j<positionConflictCount; j++)
		{
			string positionNoStr;
			positionNoStr = copy_range<string>(*sfit++);
			int positionNo = atoi(positionNoStr.c_str());
			int positionIx = positionNo - 1;

			int pointIx = positionIx/positionNumber;
			if(currentPointIx != pointIx)
			{
				adjacencyMatrix[i][positionIx]=true;
				conflictingPositions[i].push_back(positionIx);

				vector<int>* currentConflictingPoints = conflictingPoints+currentPointIx;
				bool isConflictingPointFound = false;
				for(vector<int>::iterator it=currentConflictingPoints->begin(); it!=currentConflictingPoints->end(); it++) {
					int conflictingPoint = *it;
					if(conflictingPoint==pointIx) {
						isConflictingPointFound = true;
						break;
					}
				}
				if(!isConflictingPointFound) {
					currentConflictingPoints->push_back(pointIx);
				}
			}
		}
	}
	return conflictGraph;
}
}