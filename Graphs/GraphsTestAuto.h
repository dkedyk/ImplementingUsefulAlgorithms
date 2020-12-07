#ifndef IGMDK_GRAPHS_TEST_AUTO_H
#define IGMDK_GRAPHS_TEST_AUTO_H
using namespace std;
#include "Graph.h"
#include "NetworkFlowTestAuto.h"

namespace igmdk{

void testDFSAuto()
{
    DEBUG("testDFSAuto");
    typedef GraphAA<bool> G;
    G sp;
    for(int i = 0; i < 6; ++i)
	{
		sp.addVertex();
	}
	sp.addEdge(0,1,6);
	sp.addEdge(0,2,8);
	sp.addEdge(0,3,18);
	sp.addEdge(1,4,11);
	sp.addEdge(2,3,9);
	sp.addEdge(4,5,3);
	sp.addEdge(5,2,7);
	sp.addEdge(5,3,4);

    Vector<Vector<int> > components = connectedComponents(sp);
    assert(components.getSize() == 1);
    for(int i = 0; i < components.getSize(); ++i)
    {
        Vector<int>& c = components[i];
        assert(c.getSize() == 6);
        assert(c[0] == 0);
        assert(c[1] == 1);
        assert(c[2] == 4);
        assert(c[3] == 5);
        assert(c[4] == 2);
        assert(c[5] == 3);
    }

    Vector<int> ranks = topologicalSort(sp);
    assert(ranks.getSize() == 6);
    assert(ranks[0] == 0);
    assert(ranks[1] == 1);
    assert(ranks[2] == 4);
    assert(ranks[3] == 5);
    assert(ranks[4] == 2);
    assert(ranks[5] == 3);
    DEBUG("testDFSAuto passed");
}

void testMSTAuto()
{
    DEBUG("testMSTAuto");
    typedef GraphAA<double> G;
    G sp;
    for(int i = 0; i < 5; ++i)
	{
		sp.addVertex();
	}
	sp.addEdge(0,1,2);
	sp.addEdge(0,3,8);
	sp.addEdge(0,4,4);
	sp.addEdge(1,0,2);
	sp.addEdge(1,2,3);
	sp.addEdge(2,1,3);
	sp.addEdge(2,3,5);
	sp.addEdge(2,4,1);
	sp.addEdge(3,0,8);
	sp.addEdge(3,2,5);
	sp.addEdge(3,4,7);
	sp.addEdge(4,0,4);
	sp.addEdge(4,2,1);
	sp.addEdge(4,3,7);
    Vector<int> parents = MST(sp);
    assert(parents[0] == -1);
    assert(parents[1] == 0);
    assert(parents[2] == 1);
    assert(parents[3] == 2);
    assert(parents[4] == 2);
    DEBUG("testMSTAuto passed");
}

void testShortestPathAuto()
{
    DEBUG("testShortestPathAuto");
    typedef GraphAA<double> G;
    G sp;
	for(int i = 0; i < 6; ++i)
	{
		sp.addVertex();
	}
	sp.addEdge(0,1,6);
	sp.addEdge(0,2,8);
	sp.addEdge(0,3,18);
	sp.addEdge(1,4,11);
	sp.addEdge(2,3,9);
	sp.addEdge(4,5,3);
	sp.addEdge(5,2,7);
	sp.addEdge(5,3,4);

	Vector<Vector<int> > preds;
	preds.append(ShortestPath(sp, 0));
	BellmanFord<G> dk(sp, 0);
	preds.append(dk.pred);
    for(int i = 0; i < preds.getSize(); ++i)
    {
        Vector<int>& pred = preds[i];
        assert(pred.getSize() == 6);
        assert(pred[0] == -1);
        assert(pred[1] == 0);
        assert(pred[2] == 0);
        assert(pred[3] == 2);
        assert(pred[4] == 1);
        assert(pred[5] == 4);
    }
    DEBUG("testShortestPathAuto passed");
}

void testStableMatchingAuto()
{
    DEBUG("testStableMatchingAuto");
    Vector<Vector<int> > womenOrder, menRanks;
    Vector<int> order1, order2, order3, rank1, rank2, rank3;
    order1.append(0);
    order1.append(1);
    order1.append(2);

    order2.append(1);
    order2.append(2);
    order2.append(0);

    order3.append(2);
    order3.append(1);
    order3.append(0);
    womenOrder.append(order1);
    womenOrder.append(order2);
    womenOrder.append(order3);
    rank1.append(0);
    rank1.append(1);
    rank1.append(2);

    rank2.append(1);
    rank2.append(0);
    rank2.append(2);

    rank3.append(2);
    rank3.append(1);
    rank3.append(0);
    menRanks.append(rank1);
    menRanks.append(rank2);
    menRanks.append(rank3);
    Vector<int> womenResult = stableMatching(womenOrder, menRanks);
    assert(womenResult[0] == 2);
    assert(womenResult[1] == 0);
    assert(womenResult[2] == 1);
    DEBUG("testStableMatchingAuto passed");
}

void testAllAutoGraphs()
{
    DEBUG("testAllAutoGraphs");
    testDFSAuto();
    testMSTAuto();
    testShortestPathAuto();
    testAllAutoNetworkFlow();
    testStableMatchingAuto();
}

}//end namespace
#endif
