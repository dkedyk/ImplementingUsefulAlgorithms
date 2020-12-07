#include "Graph.h"
#include "GraphsTestAuto.h"
#include "../Utils/Debug.h"
using namespace igmdk;

void DDDGraph()
{
    typedef GraphAA<double> G;
    G Graph05;
	for(int i = 0; i < 6; ++i)
	{
		Graph05.addVertex();
	}
	Graph05.addEdge(0,1,6);
	Graph05.addEdge(0,2,8);
	Graph05.addEdge(0,3,18);
	Graph05.addEdge(1,4,11);
	Graph05.addEdge(2,3,9);
	Graph05.addEdge(4,5,3);
	Graph05.addEdge(5,2,7);
	Graph05.addEdge(5,3,4);

	cout << "breakpoint" << endl;
}

int main()
{
	testAllAutoGraphs();
    DDDGraph();
	return 0;
}
