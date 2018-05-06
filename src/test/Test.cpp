#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>

#include "BoundDiameterMinCutClustering.h"
#include "ConflictGraphLoader.h"
#include "GAPointLabelPlacement.h"
#include "Graph.h"
#include "GraphUtils.h"
#include "Solution.h"

using namespace std;
using namespace labelplacement;
using namespace graph;
using namespace clustering;

void calculateTotalNumberOfConflictingLabels(bool** adjacencyMatrix,
		int pointNumber, int positionNumber);
Solution& randomPlacement(ConflictGraph& conflictGraph);
void testFindComponents();
void testAllPairsShortestPaths();
void testCheckDiameter();
void testRandomCut();
void testBoundDiameterMinCutCluster();
void testMemoryLeakBoundDiameterMinCutCluster();
void testGAPointLabelPlacementOptimize();
void testClusterConflictGraphs();
int main(int argc, char **argv) {
	testGAPointLabelPlacementOptimize();
	return 0;
}

void calculateTotalNumberOfConflictingLabels(bool** adjacencyMatrix,
		int pointNumber, int positionNumber) {
	int totalConflictSize = 0;
	for (int i = 0; i < pointNumber * positionNumber; i++) {
		for (int j = 0; j < pointNumber * positionNumber; j++) {
			if (i / positionNumber != j / positionNumber)
				if (adjacencyMatrix[i][j]) {
					totalConflictSize++;
					break;
				}
		}
	}
	cout << totalConflictSize << endl;
}

Solution& randomPlacement(ConflictGraph& conflictGraph) {
	Solution* solution = new Solution(&conflictGraph);
	Solution& solutionRef = *solution;
	int* labelPlacements =
			new int[conflictGraph.getConflictGraphOfPoints()->getVertexNumber()];
	solution->setLabelPlacements(labelPlacements);
	std::srand(std::time(0)); //use current time as seed for random generator

	for (int i = 0;
			i < conflictGraph.getConflictGraphOfPoints()->getVertexNumber();
			i++) {
		int random_variable = std::rand();
		labelPlacements[i] = random_variable
				% conflictGraph.getPositionNumber();
	}
	return solutionRef;
}

void testGAPointLabelPlacementOptimize() {
	ConflictGraphLoader cgl;
	GAPointLabelPlacement ga;
	string rootpath =
			string(
					"/home/oakile/Workspace/FastStaticPointLabelPlacement/problem_instance/");
	int sets[5] = { 100, 250, 500, 750, 1000 };
	for (int c = 0; c < 5; c++) {
		string pathofset = "d" + to_string(sets[c]) + "/";
		for (int i = 1; i <= 25; i++) {
			string filename = "d" + to_string(sets[c]) + "_"
					+ (i < 10 ? "0" : "") + to_string(i) + ".dat";
			cout << "I: " << filename;
			string path = rootpath + pathofset + filename;
			ConflictGraph* cg = cgl.load(path);
			//calculateTotalNumberOfConflictingLabels(cg->getConflictGraphOfPositions()->getAdjacencyMatrix(), cg->getConflictGraphOfPoints()->getVertexNumber(), cg->getConflictGraphOfPositions()->getVertexNumber());
			clock_t startG = clock();
			Solution& solution = ga.optimize(*cg);
			clock_t endG = clock();
			cout << "\tC: " << solution.getConflictSize();
			//Solution& solution = randomPlacement(*cg);
			int gaTimeCost = 1000 * (endG - startG) / CLOCKS_PER_SEC;
			cout << "\tT: " << gaTimeCost << endl;
			/*for(int i=0; i<cg->getPointNumber(); i++)
			 {
			 cout<<solution.getLabelPlacements()[i]<<" ";
			 }
			 cout<<endl;*/

			delete &solution;

			/*GAPointLabelPlacement::setConflictGraph(cg);
			 GAPointLabelPlacement::generateGroupedMask(2,0.70);
			 cout<<"Number of groupes: "<<ga.getNumberOfMaskGroupes()<<endl;
			 cout<<"Mask:: "<<endl;
			 int* mask = GAPointLabelPlacement::getGroupedMask();
			 for(int i=0; i<cg->getPointNumber();i++) {
			 cout<<mask[i]<<" ";
			 }*/
			delete cg;
		}
	}
}

void testRandomPlacement() {

}
Graph* getGraphInstance1() {
	Graph* graph = new Graph();
	graph->setVertexNumber(7);
	vector<int>* vertices = new vector<int>();
	vertices->push_back(1);
	vertices->push_back(2);
	vertices->push_back(3);
	vertices->push_back(4);
	vertices->push_back(5);
	vertices->push_back(6);
	vertices->push_back(7);
	graph->setVertices(vertices);

	vector<int>* adjacencyList = new vector<int> [7];
	adjacencyList[0].push_back(4);
	adjacencyList[1].push_back(2);
	adjacencyList[1].push_back(3);
	adjacencyList[2].push_back(1);
	adjacencyList[2].push_back(3);
	adjacencyList[2].push_back(6);
	adjacencyList[3].push_back(1);
	adjacencyList[3].push_back(2);
	adjacencyList[3].push_back(6);
	adjacencyList[4].push_back(0);
	adjacencyList[6].push_back(2);
	adjacencyList[6].push_back(3);
	graph->setAdjacencyList(adjacencyList);
	graph->generateAdjacencyMatrixFromAdjacencyList();
	return graph;
}

Graph* getGraphInstance2() {
	Graph* graph = new Graph();
	graph->setVertexNumber(8);
	vector<int>* vertices = new vector<int>();
	vertices->push_back(1);
	vertices->push_back(2);
	vertices->push_back(3);
	vertices->push_back(4);
	vertices->push_back(5);
	vertices->push_back(6);
	vertices->push_back(7);
	vertices->push_back(8);
	graph->setVertices(vertices);

	vector<int>* adjacencyList = new vector<int> [8];
	adjacencyList[0].push_back(2);
	adjacencyList[1].push_back(2);
	adjacencyList[2].push_back(0);
	adjacencyList[2].push_back(1);
	adjacencyList[2].push_back(3);
	adjacencyList[2].push_back(4);
	adjacencyList[3].push_back(2);
	adjacencyList[3].push_back(5);
	adjacencyList[4].push_back(2);
	adjacencyList[4].push_back(5);
	adjacencyList[5].push_back(3);
	adjacencyList[5].push_back(4);
	adjacencyList[5].push_back(6);
	adjacencyList[5].push_back(7);
	adjacencyList[6].push_back(5);
	adjacencyList[7].push_back(5);
	graph->setAdjacencyList(adjacencyList);
	graph->generateAdjacencyMatrixFromAdjacencyList();

	return graph;
}
void testFindComponents() {
	Graph* graph = getGraphInstance1();

	BoundDiameterMinCutClustering clustering;
	vector<Graph*>* components = clustering.findComponentsBFS(graph, -1);
	for (vector<Graph*>::iterator it = components->begin();
			it != components->end(); it++) {

		Graph* component = *it;
		GraphUtils::printGraph(component);
	}
}

void testAllPairsShortestPaths() {
	Graph* graph = getGraphInstance1();

	BoundDiameterMinCutClustering clustering;
	vector<Graph*>* clusters = clustering.findComponentsBFS(graph, -1);
	for (vector<Graph*>::iterator it = clusters->begin(); it != clusters->end();
			it++) {
		Graph* component = *it;
		GraphUtils::printGraph(component);
		int componentSize = component->getVertexNumber();
		int **shortestPaths = new int*[componentSize];
		for (int i = 0; i < componentSize; i++) {
			shortestPaths[i] = new int[componentSize];
		}
		clustering.floydWarshallAllPairsShortestPaths(component, shortestPaths);

		for (int i = 0; i < componentSize; i++) {
			for (int j = 0; j < componentSize; j++) {
				cout << shortestPaths[i][j] << " ";
			}
			cout << endl;
		}
	}
}
void testCheckDiameter() {
	Graph* graph = getGraphInstance1();
	BoundDiameterMinCutClustering clustering;
	vector<Graph*>* components = clustering.findComponentsBFS(graph, -1);
	for (vector<Graph*>::iterator it = components->begin();
			it != components->end(); it++) {
		Graph* component = *it;
		GraphUtils::printGraph(component);
		cout << "Diameter: " << clustering.getDiameterOfComponent(component)
				<< endl;
	}
}

void testRandomCut() {

	Graph* graph = getGraphInstance1();
	BoundDiameterMinCutClustering clustering;
	vector<Graph*>* components = clustering.findComponentsBFS(graph, -1);
	for (vector<Graph*>::iterator it = components->begin();
			it != components->end(); it++) {
		Graph* component = *it;
		GraphUtils::printGraph(component);
		Graph subComponent1, subComponent2;

		vector<pair<int, int>*>* cut = clustering.randomCut(component,
				&subComponent1, &subComponent2);
		cout << "SubComponent 1 " << subComponent1.getAdjacencyList()->size()
				<< endl;
		GraphUtils::printGraph(&subComponent1);
		cout << "SubComponent 2 " << subComponent2.getAdjacencyList()->size()
				<< endl;
		GraphUtils::printGraph(&subComponent2);
		if (cut != NULL) {
			cout << "Cut Edges Size: " << cut->size() << endl;
			for (vector<pair<int, int>*>::iterator it2 = cut->begin();
					it2 != cut->end(); it2++) {
				pair<int, int>* cutEdge = *it2;
				cout << "(" << component->getVertices()->at(cutEdge->first)
						<< "," << component->getVertices()->at(cutEdge->second)
						<< "), ";
				delete cutEdge;
			}
			delete cut;
			cout << endl << "next component" << endl;
			;
		} else {
			cout << "NO CUT" << endl;
		}
	}
}
void testBoundDiameterMinCutCluster() {

	Graph* graph = getGraphInstance1();
	BoundDiameterMinCutClustering clustering;
	vector<Graph*>* components = clustering.cluster(graph, 1);
	for (vector<Graph*>::iterator it = components->begin();
			it != components->end(); it++) {
		Graph* component = *it;
		GraphUtils::printGraph(component);
		delete component;
	}
	delete components;
	delete graph;

	graph = getGraphInstance2();
	components = clustering.cluster(graph, 2);
	for (vector<Graph*>::iterator it = components->begin();
			it != components->end(); it++) {
		Graph* component = *it;
		GraphUtils::printGraph(component);
		delete component;
	}
	delete components;
	delete graph;
}

void testMemoryLeakBoundDiameterMinCutCluster() {
	for (int i = 0; i < 1000000; i++) {
		testBoundDiameterMinCutCluster();
	}
}

void testClusterConflictGraphs() {
	ConflictGraphLoader cgl;

	string rootpath =
			string(
					"/home/oakile/Workspace/FastStaticPointLabelPlacement/problem_instance/");
	int sets[1] = { 1000 };
	for (int c = 0; c < 1; c++) {
		string pathofset = "d" + to_string(sets[c]) + "/";
		for (int i = 1; i <= 1; i++) {
			string filename = "d" + to_string(sets[c]) + "_"
					+ (i < 10 ? "0" : "") + to_string(i) + ".dat";
			cout << "I: " << filename << endl;
			cout
					<< "================================================================================"
					<< endl;
			string path = rootpath + pathofset + filename;
			ConflictGraph* cg = cgl.load(path);
			BoundDiameterMinCutClustering clustering;
			vector<Graph*>* components = clustering.cluster(
					cg->getConflictGraphOfPoints(), 2);
			for (vector<Graph*>::iterator it = components->begin();
					it != components->end(); it++) {
				Graph* component = *it;
				GraphUtils::printGraph(component);
				delete component;
			}
			delete cg;
		}
	}
}
