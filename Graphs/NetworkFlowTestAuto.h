#ifndef IGMDK_NETWORK_FLOW_TEST_AUTO_H
#define IGMDK_NETWORK_FLOW_TEST_AUTO_H
using namespace std;
#include "NetworkFlow.h"

namespace igmdk{

void testMaxFlowAuto()
{
    DEBUG("testMaxFlowAuto");
    GraphAA<int> sp(6);
    Vector<FlowData> data;
    data.append(FlowData(0, 2));//to 1
	sp.addUndirectedEdge(0,1,0);
	data.append(FlowData(0, 3));//to 2
	sp.addUndirectedEdge(0,2,1);
	data.append(FlowData(1, 3));//to 3
	sp.addUndirectedEdge(1,3,2);
	data.append(FlowData(1, 1));//to 4
	sp.addUndirectedEdge(1,4,3);
	data.append(FlowData(2, 1));//to 3
	sp.addUndirectedEdge(2,3,4);
	data.append(FlowData(2, 1));//to 4
	sp.addUndirectedEdge(2,4,5);
	data.append(FlowData(3, 2));//to 5
	sp.addUndirectedEdge(3,5,6);
	data.append(FlowData(4, 3));//to 5
	sp.addUndirectedEdge(4,5,7);
	ShortestAugmentingPath<GraphAA<int> > dk(sp, data, 0, 5);
    assert(dk.getTotalFlow() == 4);
	assert(data[0].flow == 2);
	assert(data[1].flow == 2);
	assert(data[2].flow == 1);
	assert(data[3].flow == 1);
	assert(data[4].flow == 1);
	assert(data[5].flow == 1);
	assert(data[6].flow == 2);
	assert(data[7].flow == 2);
    DEBUG("testMaxFlowAuto passed");
}

void testMinCostAuto()
{
    DEBUG("testMinCostAuto");
    GraphAA<int> sp(6);
    Vector<FlowData> data;
    data.append(FlowData(0, 300, 0));//to 1
	sp.addUndirectedEdge(0,1,0);
	data.append(FlowData(0, 300, 0));//to 2
	sp.addUndirectedEdge(0,2,1);
	data.append(FlowData(1, 200, 7));//to 3
	sp.addUndirectedEdge(1,3,2);
	data.append(FlowData(1, 200, 6));//to 4
	sp.addUndirectedEdge(1,4,3);
	data.append(FlowData(2, 280, 4));//to 3
	sp.addUndirectedEdge(2,3,4);
	data.append(FlowData(2, 350, 6));//to 4
	sp.addUndirectedEdge(2,4,5);
	data.append(FlowData(3, 300, 0));//to 5
	sp.addUndirectedEdge(3,5,6);
	data.append(FlowData(4, 300, 0));//to 5
	sp.addUndirectedEdge(4,5,7);

	ShortestAugmentingPath<GraphAA<int> > dk(sp, data, 0, 5, 600);
	assert(dk.getTotalFlow() == 600);
	assert(data[0].flow == 300);
	assert(data[1].flow == 300);
	assert(data[2].flow == 100);
	assert(data[3].flow == 200);
	assert(data[4].flow == 200);
	assert(data[5].flow == 100);
	assert(data[6].flow == 300);
	assert(data[7].flow == 300);
    DEBUG("testMinCostAuto passed");
}

void testBipartiteAuto()
{
    DEBUG("testBipartiteAuto");
    Vector<pair<int, int> > allowed;
    allowed.append(make_pair(0, 5));
    allowed.append(make_pair(1, 4));
    allowed.append(make_pair(1, 3));
    allowed.append(make_pair(2, 4));
    allowed = bipartiteMatching(3, 3, allowed);

    assert(allowed.getSize() == 3);
    assert(allowed[0].first == 0);
    assert(allowed[0].second == 5);
    assert(allowed[1].first == 1);
    assert(allowed[1].second == 3);
    assert(allowed[2].first == 2);
    assert(allowed[2].second == 4);
    DEBUG("testBipartiteAuto passed");
}

void testAssignmentAuto()
{
    DEBUG("testAssignmentAuto");
    Vector<pair<pair<int, int>, double> > allowed;
    allowed.append(make_pair(make_pair(0, 5), 0));
    allowed.append(make_pair(make_pair(1, 4), 0));
    allowed.append(make_pair(make_pair(1, 3), 0));
    allowed.append(make_pair(make_pair(2, 4), 0));
    Vector<pair<int, int> > allowed2 = assignmentProblem(3, 3, allowed);
    assert(allowed2.getSize() == 3);
    assert(allowed2[0].first == 0);
    assert(allowed2[0].second == 5);
    assert(allowed2[1].first == 1);
    assert(allowed2[1].second == 3);
    assert(allowed2[2].first == 2);
    assert(allowed2[2].second == 4);
    DEBUG("testAssignmentAuto passed");
}

void testAllAutoNetworkFlow()
{
    DEBUG("testAllAutoNetworkFlow");
    testMaxFlowAuto();
    testMinCostAuto();
    testBipartiteAuto();
    testAssignmentAuto();
}

}//end namespace
#endif
