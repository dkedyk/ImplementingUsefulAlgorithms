#ifndef IGMDK_OPT_TEST_AUTO_H
#define IGMDK_OPT_TEST_AUTO_H
#include <string>
using namespace std;
#include "SearchAlgorithms.h"
#include "../Graphs/Graph.h"

namespace igmdk{

struct GraphProblem
{
    typedef int STATE_ID;
    int nullState()const{return -1;}
    typedef EHash<BUHash> HASHER;
    GraphAA<double> graph;
    int from, to;
    int start()const{return from;}
    template<typename DUMMY>
    bool isGoal(int i, DUMMY const& d)const{return i == to;}
    template<typename DUMMY>
    Vector<int> nextStates(int j, DUMMY const& d)const
    {
        Vector<int> result;
        for(GraphAA<double>::AdjacencyIterator i = graph.begin(j);
            i != graph.end(j); ++i) result.append(i.to());
        return result;
    }
    template<typename DUMMY>
    double remainderLowerBound(int dummy, int i, DUMMY const& d)const{return 0;}
    template<typename DUMMY>
    double distance(int k, int j, DUMMY const& d)const
    {
        for(GraphAA<double>::AdjacencyIterator i = graph.begin(j);
            i != graph.end(j); ++i) if(i.to() == k) return i.data();
        return 0;
    }
};

void testAStarAuto()
{
    typedef GraphAA<double> G;
    G sp;
    GraphProblem Gp;
	for(int i = 0; i < 6; ++i)
	{
		Gp.graph.addVertex();
	}
	Gp.graph.addEdge(0,1,6);
	Gp.graph.addEdge(0,2,8);
	Gp.graph.addEdge(0,3,18);
	Gp.graph.addEdge(1,4,11);
	Gp.graph.addEdge(2,3,9);
	Gp.graph.addEdge(4,5,3);
	Gp.graph.addEdge(5,2,7);
	Gp.graph.addEdge(5,3,4);
	Gp.from = 0;
	Gp.to = 5;

	pair<Vector<int>, bool> result = AStar<GraphProblem>::solve(Gp);
	assert(result.second);
	Vector<int> expected;
	expected.append(0);
	expected.append(1);
	expected.append(4);
	expected.append(5);
	assert(result.first == expected);
	DEBUG("testAStartAuto passed");
}

void testRBFSAuto()
{
    typedef GraphAA<double> G;
    G sp;
    GraphProblem Gp;
	for(int i = 0; i < 6; ++i)
	{
		Gp.graph.addVertex();
	}
	Gp.graph.addEdge(0,1,6);
	Gp.graph.addEdge(0,2,8);
	Gp.graph.addEdge(0,3,18);
	Gp.graph.addEdge(1,4,11);
	Gp.graph.addEdge(2,3,9);
	Gp.graph.addEdge(4,5,3);
	Gp.graph.addEdge(5,2,7);
	Gp.graph.addEdge(5,3,4);
	Gp.from = 0;
	Gp.to = 5;

	RecursiveBestFirstSearch<GraphProblem> dk(Gp);
	assert(dk.foundGoal);
	assert(dk.pred.pop() == 5);
	assert(dk.pred.pop() == 4);
	assert(dk.pred.pop() == 1);
	assert(dk.pred.pop() == 0);
	assert(dk.pred.isEmpty());
	DEBUG("testRBFSAuto passed");
}

void testAllAutoOpt()
{
    DEBUG("testAllAutoOpt");
    testAStarAuto();
    testRBFSAuto();
}

}//end namespace
#endif
