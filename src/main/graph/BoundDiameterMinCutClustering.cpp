/*
 * BoundDiamMinCutClustering.cpp
 *
 *  Created on: Jan 31, 2018
 *      Author: oakile
 */

#include "BoundDiameterMinCutClustering.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <limits.h>
#include <map>
#include <iostream>
#include <ga/garandom.h>

#include "GraphUtils.h"

using namespace std;
namespace graph {
namespace clustering {

BoundDiameterMinCutClustering::BoundDiameterMinCutClustering() {
	this->graph = NULL;
	this->diameter = 0;
}

BoundDiameterMinCutClustering::~BoundDiameterMinCutClustering() {

}
int sumWithInfinity(int a, int b) {
	if (a == -1 || b == -1) {
		return -1;
	} else {
		return a + b;
	}
}

int compareWithInfinity(int a, int b) {
	if (a == b) {
		return 0;
	} else if (a == -1 || (b != -1 && a > b)) {
		return -1;
	} else {
		return 1;
	}
}
Graph* BoundDiameterMinCutClustering::getGraph() const {
	return this->graph;
}
void BoundDiameterMinCutClustering::setGraph(Graph* graph) {
	this->graph = graph;
}
vector<Graph*>* BoundDiameterMinCutClustering::findComponentsBFS(Graph* graph,
		int vertexSizeLimit) {
	int vertexNumber = graph->getVertexNumber();

	int verticesToBeVisited[vertexNumber];
	for (int i = 0; i < vertexNumber; i++) {
		verticesToBeVisited[i] = i;
	}
	int vertex;
	vector<Graph*>* components = new vector<Graph*>();
	while ((vertex = nextVertex(verticesToBeVisited, vertexNumber)) >= 0) {
		Graph* component = new Graph();
		vector<int>* componentVertices = new vector<int>();
		vector<int>* componentVertexIndices = new vector<int>();
		int currentVertexIx = 0;
		map<int, int> vertexIndicesMap;
		component->setVertices(componentVertices);
		components->push_back(component);

		componentVertexIndices->push_back(vertex);
		componentVertices->push_back(graph->getVertices()->at(vertex));
		vertexIndicesMap[vertex] = currentVertexIx++;
		verticesToBeVisited[vertex] = -1;

		for (int i = 0;
				i < componentVertexIndices->size()
						&& compareWithInfinity(currentVertexIx, vertexSizeLimit)
								> 0; i++) {
			int vertexToIterate = (*componentVertexIndices)[i];
			vector<int>* adjacencyList = graph->getAdjacencyList()
					+ vertexToIterate;
			for (vector<int>::iterator it2 = adjacencyList->begin();
					it2 != adjacencyList->end()
							&& compareWithInfinity(currentVertexIx,
									vertexSizeLimit) > 0; it2++) {
				int adjacentVertex = *it2;
				bool adjacentVertexIsFound = false;
				for (int j = 0; j < vertexNumber; j++) {
					if (verticesToBeVisited[j] == adjacentVertex) {
						adjacentVertexIsFound = true;
						break;
					}
				}
				if (adjacentVertexIsFound) {
					componentVertexIndices->push_back(adjacentVertex);
					componentVertices->push_back(
							graph->getVertices()->at(adjacentVertex));
					vertexIndicesMap[adjacentVertex] = currentVertexIx++;
					verticesToBeVisited[adjacentVertex] = -1;
				}
			}
		}
		int componentVertexSize = componentVertexIndices->size();
		component->setVertexNumber(componentVertexSize);
		vector<int>* componentAdjacencyList =
				new vector<int> [componentVertexSize];

		for (int i = 0; i < componentVertexSize; i++) {
			int vertex = componentVertexIndices->at(i);
			vector<int>* adjacencyList = &graph->getAdjacencyList()[vertex];
			for (vector<int>::iterator it2 = adjacencyList->begin();
					it2 != adjacencyList->end(); it2++) {
				int adjacentVertex = *it2;
				bool adjacentVertexIsFound = false;
				for (int j = 0; j < componentVertexSize; j++) {
					if (componentVertexIndices->at(j) == adjacentVertex) {
						adjacentVertexIsFound = true;
						break;
					}
				}
				if (adjacentVertexIsFound) {
					componentAdjacencyList[i].push_back(
							vertexIndicesMap[adjacentVertex]);
				}
			}
		}
		component->setAdjacencyList(componentAdjacencyList);
		component->generateAdjacencyMatrixFromAdjacencyList();
		delete componentVertexIndices;
	}
	return components;
}

void BoundDiameterMinCutClustering::floydWarshallAllPairsShortestPaths(
		Graph* graph, int** shortestPaths) {
	int vertexNumber = graph->getVertexNumber();

	for (int i = 0; i < vertexNumber; i++) {
		for (int j = 0; j < vertexNumber; j++) {
			if (i == j) {
				shortestPaths[i][i] = 0;
			} else {
				if (graph->getAdjacencyMatrix()[i][j]) {
					shortestPaths[i][j] = 1;
				} else {
					shortestPaths[i][j] = -1;
				}
			}
		}
	}
	for (int i = 0; i < vertexNumber; i++) {
		for (int j = 0; j < vertexNumber; j++) {
			for (int k = 0; k < vertexNumber; k++) {
				int sumTwoOther = sumWithInfinity(shortestPaths[j][i],
						shortestPaths[i][k]);
				if (compareWithInfinity(shortestPaths[j][k], sumTwoOther)
						== -1) {
					shortestPaths[j][k] = sumTwoOther;
				}
			}
		}
	}
}
int BoundDiameterMinCutClustering::getDiameterOfComponent(Graph* component) {
	int vertexNumber = component->getVertexNumber();
	int** shortestPaths = new int*[vertexNumber];
	for (int i = 0; i < vertexNumber; i++) {
		shortestPaths[i] = new int[vertexNumber];
	}
	floydWarshallAllPairsShortestPaths(component, shortestPaths);
	int maxLength = 0;
	for (int i = 1; i < vertexNumber; i++) {
		for (int j = 0; j < i; j++) {
			if (maxLength < shortestPaths[i][j]) {
				maxLength = shortestPaths[i][j];
			}
		}
	}
	for (int i = 0; i < vertexNumber; i++) {
		delete[] shortestPaths[i];
	}
	delete[] shortestPaths;

	return maxLength;
}

vector<pair<int, int>*>* BoundDiameterMinCutClustering::randomCut(
		Graph* component, Graph* subComponent1, Graph* subComponent2) {
	int componentSize = component->getVertexNumber();
	if (componentSize < 2) {
		return NULL;
	}
	//cout<<"randomCut>>component size: "<<componentSize<<endl;
	vector<int>* componentVertices = component->getVertices();

	// Init Multi Adjacency Matrix, Adjacency Indicator List And Adjacency Indicator Map
	vector<pair<int, int>*>** multiGraphAdjMatrixWithSourceVertex = new vector<
			pair<int, int>*>*[componentSize];
	for (int i = 0; i < componentSize; i++) {
		multiGraphAdjMatrixWithSourceVertex[i] =
				new vector<pair<int, int>*> [componentSize];
	}

	//cout<<"randomCut>>init"<<endl;
	vector<pair<int, vector<int>*>*> adjacencyIndicatorList;
	vector<int>* adjacencyIndicatorMap = new vector<int> [componentSize];
	for (int vertex = 0; vertex < componentSize; vertex++) {
		pair<int, vector<int>*>* adjacencyIndicatorListOfVertexPair = new pair<
				int, vector<int>*>();
		adjacencyIndicatorList.push_back(adjacencyIndicatorListOfVertexPair);
		vector<int>* adjacencyIndicatorListOfVertex =
				&adjacencyIndicatorMap[vertex];
		adjacencyIndicatorListOfVertexPair->first = vertex;
		adjacencyIndicatorListOfVertexPair->second =
				adjacencyIndicatorListOfVertex;
		for (vector<int>::iterator it2 =
				component->getAdjacencyList()[vertex].begin();
				it2 != component->getAdjacencyList()[vertex].end(); it2++) {
			int adjacentVertex = *it2;
			adjacencyIndicatorListOfVertex->push_back(adjacentVertex);
			multiGraphAdjMatrixWithSourceVertex[vertex][adjacentVertex].push_back(
					new pair<int, int>(vertex, adjacentVertex));
		}
	}

	//cout<<"randomCut>>contract"<<endl;
	//Contraction
	vector<int>* mergedVertices[componentSize];
	vector<int>* mergedVertexNos[componentSize];
	vector<pair<int, int>*>* mergedEdges[componentSize];
	for (int i = 0; i < componentSize; i++) {
		mergedVertices[i] = new vector<int>();
		mergedVertexNos[i] = new vector<int>();
		mergedEdges[i] = new vector<pair<int, int>*>();
		mergedVertices[i]->push_back(i);
		mergedVertexNos[i]->push_back(componentVertices->at(i));
	}
	int currentVertexSize;
	while ((currentVertexSize = adjacencyIndicatorList.size()) > 2) {
		srand(time(0));
		int contractVertexIx1 = GARandomInt(0, currentVertexSize - 1);
		pair<int, vector<int>*>* adjacencyIndicatorListOfVertexPair =
				adjacencyIndicatorList[contractVertexIx1];
		int contractVertex = adjacencyIndicatorListOfVertexPair->first;
		vector<int>* adjacencyIndicatorListOfVertex =
				adjacencyIndicatorListOfVertexPair->second;
		int targetVertexIx1 = GARandomInt(0,
				adjacencyIndicatorListOfVertex->size() - 1);
		int targetVertex = adjacencyIndicatorListOfVertex->at(targetVertexIx1);
		//cout<<"randomCut>>contract vertex: "<<contractVertex<<" target vertex: "<<targetVertex<<endl;
		//cout<<"randomCut>> mergedVertices[contractVertex] size: "<<mergedVertices[contractVertex]->size()<<" mergedEdges[contractVertex] size: "<<mergedEdges[contractVertex]->size()<<endl;
		for (int i = 0; i < mergedVertices[contractVertex]->size(); i++) {
			mergedVertices[targetVertex]->push_back(
					mergedVertices[contractVertex]->at(i));
			mergedVertexNos[targetVertex]->push_back(
					mergedVertexNos[contractVertex]->at(i));
		}
		for (int i = 0; i < mergedEdges[contractVertex]->size(); i++) {
			mergedEdges[targetVertex]->push_back(
					mergedEdges[contractVertex]->at(i));
		}
		delete mergedVertices[contractVertex];
		delete mergedVertexNos[contractVertex];
		delete mergedEdges[contractVertex];

		for (vector<int>::iterator it = adjacencyIndicatorListOfVertex->begin();
				it != adjacencyIndicatorListOfVertex->end(); it++) {
			int adjacentVertex = *it;
			//cout<<"randomCut>>adjacent vertex: "<<adjacentVertex<<endl;

			vector<int>* adjacencyOfAdjacentVertex =
					&adjacencyIndicatorMap[adjacentVertex];
			vector<pair<int, int>*>* edgesContractAdjacent =
					&multiGraphAdjMatrixWithSourceVertex[contractVertex][adjacentVertex];
			vector<pair<int, int>*>* edgesAdjacentContract =
					&multiGraphAdjMatrixWithSourceVertex[adjacentVertex][contractVertex];
			if (adjacentVertex != targetVertex) {
				vector<pair<int, int>*>* edgesTargetAdjacent =
						&multiGraphAdjMatrixWithSourceVertex[targetVertex][adjacentVertex];
				vector<pair<int, int>*>* edgesAdjacentTarget =
						&multiGraphAdjMatrixWithSourceVertex[adjacentVertex][targetVertex];

				//Update Adjacency Indicator List
				//1-If adjacent is not target and no previous edge exists between them
				//	add target to adjacency list of adjacent, vice versa.
				if (edgesTargetAdjacent->empty()) {
					//cout<<"randomCut>>indicate adjacency"<<endl;
					adjacencyOfAdjacentVertex->push_back(targetVertex);
					vector<int>* adjacencyOfTargetVertex =
							&adjacencyIndicatorMap[targetVertex];
					adjacencyOfTargetVertex->push_back(adjacentVertex);
				}
				for (int i = 0; i < edgesContractAdjacent->size(); i++) {
					edgesTargetAdjacent->push_back(
							edgesContractAdjacent->at(i));
					edgesAdjacentTarget->push_back(
							edgesAdjacentContract->at(i));
				}

			} else {
				//cout<<"randomCut>>merge edges"<<endl;
				for (int i = 0; i < edgesContractAdjacent->size(); i++) {
					mergedEdges[targetVertex]->push_back(
							edgesContractAdjacent->at(i));
					delete edgesAdjacentContract->at(i);
				}
			}
			//Update Adjacency Indicator List
			//2-Remove contract vertex from its adjacent's adjacency list.

			for (vector<int>::iterator it3 = adjacencyOfAdjacentVertex->begin();
					it3 != adjacencyOfAdjacentVertex->end(); it3++) {
				if (*it3 == contractVertex) {
					adjacencyOfAdjacentVertex->erase(it3);
					break;
				}
			}
		}

		//cout<<"randomCut>>reduce adj ind list"<<endl;
		adjacencyIndicatorList.erase(
				adjacencyIndicatorList.begin() + contractVertexIx1);
		delete adjacencyIndicatorListOfVertexPair;
	}

	int lastVertex1 = adjacencyIndicatorList[0]->first;
	int lastVertex2 = adjacencyIndicatorList[1]->first;
	//cout<<"randomCut>>adj ind list size: "<<adjacencyIndicatorList.size()<<" last vertex 1: "<<lastVertex1<<" last vertex 2: "<<lastVertex2<<endl;
	vector<int>* verticesOfSubComponent = mergedVertices[lastVertex1];
	vector<int>* vertexNosOfSubComponent = mergedVertexNos[lastVertex1];
	vector<pair<int, int>*>* edgesOfSubComponent = mergedEdges[lastVertex1];
	//cout<<"randomCut>>verticesOfSubComponent1 size: "<<verticesOfSubComponent->size()<<" edgesOfSubComponent size: "<<edgesOfSubComponent->size()<<endl;
	for (vector<int>::iterator it = verticesOfSubComponent->begin();
			it != verticesOfSubComponent->end(); it++) {
		//cout<<*it<<",";
	}
	//cout<<endl;
	buildGraph(vertexNosOfSubComponent, verticesOfSubComponent,
			edgesOfSubComponent, subComponent1);
	verticesOfSubComponent = mergedVertices[lastVertex2];
	vertexNosOfSubComponent = mergedVertexNos[lastVertex2];
	edgesOfSubComponent = mergedEdges[lastVertex2];
	//cout<<"randomCut>>verticesOfSubComponent2 size: "<<verticesOfSubComponent->size()<<" edgesOfSubComponent size: "<<edgesOfSubComponent->size()<<endl;
	for (vector<int>::iterator it = verticesOfSubComponent->begin();
			it != verticesOfSubComponent->end(); it++) {
		//cout<<*it<<",";
	}
	//cout<<endl;
	buildGraph(vertexNosOfSubComponent, verticesOfSubComponent,
			edgesOfSubComponent, subComponent2);

	vector<pair<int, int>*>* cut = new vector<pair<int, int>*>(
			multiGraphAdjMatrixWithSourceVertex[lastVertex1][lastVertex2]);
	vector<pair<int, int>*>* copyOfCut =
			&multiGraphAdjMatrixWithSourceVertex[lastVertex2][lastVertex1];

	delete[] adjacencyIndicatorMap;
	for (int i = 0; i < copyOfCut->size(); i++) {
		delete copyOfCut->at(i);
	}

	for (int i = 0; i < componentSize; i++) {
		delete[] multiGraphAdjMatrixWithSourceVertex[i];
	}
	delete[] multiGraphAdjMatrixWithSourceVertex;

	//cout<<"randomCut>>end"<<endl;

	return cut;
}
Graph** BoundDiameterMinCutClustering::boundedDiameterMinCut(Graph* component,
		int diameter) {
	int componentVertexNumber = component->getVertexNumber();
	int componentEdgeNumber = component->getEdgeNumber();
	double avgIncidenceFactor = 2 * (double) componentEdgeNumber
			/ (double) componentVertexNumber;
	//cout<<"boundedDiameterMinCut>>diameter limit: "<<diameter<<" vertex number: "<<componentVertexNumber<<" avg inc: "<<avgIncidenceFactor<<endl;
	//int numIterations = (int)ceil((componentVertexNumber-1.0)*componentVertexNumber/2.0*log(componentVertexNumber));
	int numIterations = componentVertexNumber;
	//cout<<"boundedDiameterMinCut>>numIterations: "<<numIterations<<endl;
	double minEval = INT_MAX;
	Graph** minSubComponents = new Graph*[2];
	minSubComponents[0] = NULL;
	minSubComponents[1] = NULL;
	for (int i = 0; i < numIterations; i++) {
		//cout<<"boundedDiameterMinCut>>i: "<<i<<endl;
		Graph* subComponent1 = new Graph();
		Graph* subComponent2 = new Graph();
		vector<pair<int, int>*>* cut = randomCut(component, subComponent1,
				subComponent2);
		int cutSize = 0;
		if (cut != NULL) {
			//cout<<"boundedDiameterMinCut>>cut: ";
			cutSize = cut->size();
			for (vector<pair<int, int>*>::iterator it = cut->begin();
					it != cut->end(); it++) {
				//cout<<"("<<(*it)->first<<","<<(*it)->second<<"),";
				delete *it;
			}
			//cout<<endl;
			delete cut;
		}

		//GraphUtils::printGraph(subComponent1);
		//GraphUtils::printGraph(subComponent2);
		int subComponent1Diameter = getDiameterOfComponent(subComponent1);
		int subComponent2Diameter = getDiameterOfComponent(subComponent2);
		//cout<<"boundedDiameterMinCut>>subComponent1Diameter: "<<subComponent1Diameter<<" subComponent2Diameter: "<<subComponent2Diameter<<endl;
		if (subComponent1Diameter < diameter) {
			subComponent1Diameter = diameter;
		}
		if (subComponent2Diameter < diameter) {
			subComponent2Diameter = diameter;
		}
		double subComponent1IncidenceFactor = 2.0
				* (double) subComponent1->getEdgeNumber()
				/ (double) subComponent1->getVertexNumber();
		double subComponent2IncidenceFactor = 2.0
				* (double) subComponent2->getEdgeNumber()
				/ (double) subComponent2->getVertexNumber();
		//cout<<"boundedDiameterMinCut>>subComponent1IncidenceFactor: "<<subComponent1IncidenceFactor<<" subComponent2IncidenceFactor: "<<subComponent2IncidenceFactor<<endl;
		double expectedCutSize1 = ceil(
				log2((double) subComponent1Diameter / (double) diameter))
				* subComponent1IncidenceFactor;
		double expectedCutSize2 = ceil(
				log2((double) subComponent2Diameter / (double) diameter))
				* subComponent2IncidenceFactor;
		double eval = cutSize + expectedCutSize1 + expectedCutSize2;
		//cout<<"boundedDiameterMinCut>>cut size: "<<cutSize<<" expected cut size 1: "<<expectedCutSize1<<" expected cut size 2: "<<expectedCutSize2<<" eval: "<<eval<<endl;

		if (eval < minEval) {
			//cout<<"boundedDiameterMinCut>>update min eval"<<endl;
			minEval = eval;
			if (minSubComponents[0] != NULL && minSubComponents[1] != NULL) {
				delete minSubComponents[0];
				delete minSubComponents[1];
			}
			minSubComponents[0] = subComponent1;
			minSubComponents[1] = subComponent2;
		} else {
			//cout<<"boundedDiameterMinCut>>delete current cuts"<<endl;
			delete subComponent1;
			delete subComponent2;
		}
	}
	//cout<<"boundedDiameterMinCut>>end: "<<endl;
	return minSubComponents;
}
vector<Graph*>* BoundDiameterMinCutClustering::cluster(graph::Graph* graph,
		int diameterLimit) {

	vector<Graph*>* components = findComponentsBFS(graph, 50);
	//cout<<"cluster>>components size: "<<components->size()<<endl;
	vector<Graph*>* clusters = new vector<Graph*>();
	for (int i = 0; i < components->size(); i++) {
		//cout<<"cluster>>component no: "<<i<<endl;
		Graph* component = components->at(i);
		//cout<<"cluster>>component size: "<<component->getVertexNumber();
		int diameter = getDiameterOfComponent(component);
		//cout<<" diameter: "<<diameter<<endl;
		//GraphUtils::printGraph(component);
		if (diameter <= diameterLimit) {
			//cout<<"cluster>>diameter acceptable"<<endl;
			clusters->push_back(component);
		} else {
			//cout<<"cluster>>diameter not acceptable"<<endl;
			Graph** subComponents = boundedDiameterMinCut(component,
					diameterLimit);
			if (subComponents[0] == NULL || subComponents[1] == NULL) {
				//cout<<"cluster>>one of new components null"<<endl;
				return NULL;
			}
			components->push_back(subComponents[0]);
			components->push_back(subComponents[1]);
			//cout<<"cluster>>subComponents[0] size: "<<subComponents[0]->getVertexNumber()<<endl;
			//cout<<"cluster>>subComponents[1] size: "<<subComponents[1]->getVertexNumber()<<endl;
			//GraphUtils::printGraph(component);
			delete component;
			delete[] subComponents;
			//cout<<"cluster>>component deleted"<<endl;
		}
	}
	delete components;
	//cout<<"cluster>>end"<<endl;
	return clusters;
}

void BoundDiameterMinCutClustering::buildGraph(vector<int>* vertexNos,
		vector<int>* vertices, vector<pair<int, int>*>* edges, Graph *graph) {
	int vertexNumber = vertices->size();
	map<int, int> vertexIndexMap;
	for (int i = 0; i < vertices->size(); i++) {
		vertexIndexMap[vertices->at(i)] = i;
	}
	vector<int>* adjacencyList = new vector<int> [vertexNumber];
	for (vector<pair<int, int>*>::iterator it = edges->begin();
			it != edges->end(); it++) {
		pair<int, int>* edge = *it;
		int vertex1OfEdge = edge->first;
		int vertex2OfEdge = edge->second;
		int indexOfVertex1 = vertexIndexMap[vertex1OfEdge];
		int indexOfVertex2 = vertexIndexMap[vertex2OfEdge];
		adjacencyList[indexOfVertex1].push_back(indexOfVertex2);
		adjacencyList[indexOfVertex2].push_back(indexOfVertex1);
		delete edge;
	}
	delete vertices;
	delete edges;
	graph->setVertexNumber(vertexNumber);
	graph->setVertices(vertexNos);
	graph->setAdjacencyList(adjacencyList);
	graph->generateAdjacencyMatrixFromAdjacencyList();
}
int BoundDiameterMinCutClustering::nextVertex(int* vertices, int size) {
	int vertexFound = -1;
	for (int i = 0; i < size; i++) {
		if (vertices[i] >= 0) {
			vertexFound = vertices[i];
			break;
		}
	}
	return vertexFound;
}
}
} /* namespace labelplacement */
