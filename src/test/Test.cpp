#include <iostream>
#include <cstdlib>
#include <ctime>
#include <string>

#include "ConflictGraphLoader.h"
#include "GAPointLabelPlacement.h"
#include "Solution.h"

using namespace std;
using namespace labelplacement;

void calculateTotalNumberOfConflictingLabels(bool** adjacencyMatrix, int pointNumber, int positionNumber);
Solution& randomPlacement(ConflictGraph& conflictGraph);

int main(int argc, char **argv)
{
	ConflictGraphLoader cgl;
	GAPointLabelPlacement ga;

	string rootpath = string("/home/oakile/Workspace/FastStaticPointLabelPlacement/problem_instance/");
	int sets[5] = {100,250,500,750,1000};
	for(int c=0; c<5; c++)
	{
		string pathofset = "d" + to_string(sets[c]) + "/";
		for (int i=1; i<=25; i++)
		{
			string filename = "d" + to_string(sets[c]) + "_" + (i<10?"0":"") + to_string(i) + ".dat";
			cout<<"I: "<<filename;
			string path = rootpath + pathofset + filename;
			ConflictGraph* cg = cgl.load(path);
			//calculateTotalNumberOfConflictingLabels(cg->getAdjacencyMatrix(), cg->getPointNumber(), cg->getPositionNumber());
			clock_t startG = clock();
			Solution& solution = ga.optimize(*cg);
			clock_t endG = clock();
			cout<<"\tC: "<<solution.getConflictSize();
			//Solution& solution = randomPlacement(*cg);
			int gaTimeCost = 1000 * (endG - startG) / CLOCKS_PER_SEC;
			cout<<"\tT: "<<gaTimeCost<<endl;
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
	return 0;
}

void calculateTotalNumberOfConflictingLabels(bool** adjacencyMatrix, int pointNumber, int positionNumber)
{
	int totalConflictSize = 0;
	for(int i=0; i<pointNumber*positionNumber; i++)
	{
		for(int j=0; j<pointNumber*positionNumber; j++)
		{
			if(i/positionNumber!=j/positionNumber)
			if(adjacencyMatrix[i][j])
			{
				totalConflictSize ++;
				break;
			}
		}
	}
	cout<<totalConflictSize<<endl;
}

Solution& randomPlacement(ConflictGraph& conflictGraph)
{
	Solution* solution = new Solution(&conflictGraph);
	Solution& solutionRef = *solution;
	int* labelPlacements = new int[conflictGraph.getPointNumber()];
	solution->setLabelPlacements(labelPlacements);
	std::srand(std::time(0)); //use current time as seed for random generator

	for(int i=0; i<conflictGraph.getPointNumber(); i++)
	{
		int random_variable = std::rand();
		labelPlacements[i] = random_variable%conflictGraph.getPositionNumber();
	}
	return solutionRef;
}
