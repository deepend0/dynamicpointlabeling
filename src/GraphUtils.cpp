/*
 * GraphUtils.cpp
 *
 *  Created on: Feb 17, 2018
 *      Author: oakile
 */

#include "GraphUtils.h"
#include <iostream>
#include <vector>

using namespace std;
namespace graph {

GraphUtils::GraphUtils() {
	// TODO Auto-generated constructor stub

}

GraphUtils::~GraphUtils() {
	// TODO Auto-generated destructor stub
}

void GraphUtils::printGraph(Graph* graph) {
	int vertexNumber = graph->getVertexNumber();
	vector<int>* vertices = graph->getVertices();

	cout<<"   ";
	for(vector<int>::iterator it2=vertices->begin(); it2!=vertices->end(); it2++) {
		cout<<*it2<<" ";
	}
	cout<<endl;
	cout<<"   ";
	for(int i=0; i<vertexNumber; i++) {
		cout<<"__";
	}
	cout<<endl;
	bool** adjacencyMatrix = graph->getAdjacencyMatrix();
	for(int i=0; i<vertexNumber; i++) {
		cout<<vertices->at(i)<<"| ";
		for(int j=0; j<vertexNumber; j++) {
			cout << (adjacencyMatrix[i][j]?"1":"0") << " ";
		}
		cout<<endl;
	}
	cout<<endl;
}
} /* namespace labelplacement */
