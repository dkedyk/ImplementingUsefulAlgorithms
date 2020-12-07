#ifndef IGMDK_NETWORK_FLOW_H
#define IGMDK_NETWORK_FLOW_H
#include "Graph.h"
#include "../Utils/Vector.h"
#include "../Utils/Debug.h"
#include "../RandomNumberGeneration/Random.h"
#include "../Utils/Queue.h"
#include "../Utils/UnionFind.h"
#include <cassert>
#include "../Heaps/Heap.h"
using namespace std;
namespace igmdk{

struct FlowData
{
    int from;
    double flow, capacity, cost;//cost used only for min flow
    FlowData(int theFrom, double theCapacity, double theCost = 0):
        from(theFrom), capacity(theCapacity), flow(0), cost(theCost) {}
    double capacityTo(int v)const{return v == from ? flow : capacity - flow;}
    //flow can step out of (0, capacity) numerically but OK
    void addFlowTo(int v, double change)
        {flow += change * (v == from ? -1 : 1);}
};
template<typename GRAPH> class ShortestAugmentingPath
{
    int v;
    Vector<int> path, pred;
    double totalFlow;
    bool hasAugmentingPath(GRAPH const& g, Vector<FlowData>& fedges, int from,
        int to)
    {
        for(int i = 0; i < v; ++i) pred[i] = -1;
        Queue<int> queue;
        queue.push(pred[from] = from);
        while(!queue.isEmpty())
        {
            int i = queue.pop();
            for(typename GRAPH::AdjacencyIterator j = g.begin(i);
                j != g.end(i); ++j)
                if(pred[j.to()] == -1 &&//unvisited with capacity
                    fedges[j.data()].capacityTo(j.to()) > 0)
                {
                    pred[j.to()] = i;
                    path[j.to()] = j.data();
                    queue.push(j.to());
                }
        }
        return pred[to] != -1;
    }
    bool hasMinCostAugmentingPath(GRAPH const& g, Vector<FlowData>& fedges,
        int from, int to, double neededFlow)
    {//if need more flow, make a graph from edges with available capacity
        if(totalFlow >= abs(neededFlow)) return false;
        GraphAA<double> costGraph(v);
        for(int i = 0; i < v; ++i)
            for(typename GRAPH::AdjacencyIterator j = g.begin(i);
                j != g.end(i); ++j)
                if(fedges[j.data()].capacityTo(j.to()) > 0)
                    costGraph.addEdge(i, j.to(), fedges[j.data()].cost);
        if(neededFlow > 0) pred = ShortestPath(costGraph, from, to);
        else
        {//negative costs
            BellmanFord<GraphAA<double> > bf(costGraph, from);
            if(bf.hasNegativeCycle)
            {
                totalFlow = numeric_limits<double>::infinity();
                return false;
            }
            pred = bf.pred;
        };
        //extract edges for the path
        for(int i = to; pred[i] != -1; i = pred[i])
            for(typename GRAPH::AdjacencyIterator j = g.begin(pred[i]);
                j != g.end(pred[i]); ++j)
                if(j.to() == i)
                {
                    path[i] = j.data();
                    break;
                }
        return pred[to] != -1;
    }
public:
    double getTotalFlow()const{return totalFlow;}
    ShortestAugmentingPath(GRAPH const& g, Vector<FlowData>& fedges, int from,
        int to, double neededFlow = 0): v(g.nVertices()), totalFlow(0),
        path(v, -1), pred(v, -1)
    {//iteratively, first find a path
        assert(from >= 0 && from < v && to >= 0 && to < v);
        while(neededFlow == 0 ? hasAugmentingPath(g, fedges, from, to) :
            hasMinCostAugmentingPath(g, fedges, from, to, neededFlow))
        {//then from it the amount of flow to add
            double increment = numeric_limits<double>::max();
            for(int j = to; j != from; j = pred[j])
                increment = min(increment, fedges[path[j]].capacityTo(j));
            if(neededFlow != 0)//only relevant to min cost flow
                increment = min(increment, abs(neededFlow) - totalFlow);
            //then add to all edges
            for(int j = to; j != from; j = pred[j])
                fedges[path[j]].addFlowTo(j, increment);
            totalFlow += increment;
        }
    }
};

Vector<pair<int, int> > bipartiteMatching(int n, int m,
    Vector<pair<int, int> > const& allowedMatches)
{//v = n + m + 2, e = n + m + allowedMatches.getSize(), time is O(ve)
    GraphAA<int> sp(n + m + 2);//setup graph and flow edges
    Vector<FlowData> data;
    for(int i = 0; i < allowedMatches.getSize(); ++i)
    {
        data.append(FlowData(allowedMatches[i].first, 1));
        sp.addUndirectedEdge(allowedMatches[i].first,
            allowedMatches[i].second, i);
    }
    int source = n + m, sink = source + 1;//setup source and sink groups
    for(int i = 0; i < source; ++i)
    {
        int from = i, to = sink;
        if(i < n)
        {
            from = source;
            to = i;
        }
        data.append(FlowData(from, 1));
        sp.addUndirectedEdge(from, to, i + allowedMatches.getSize());
    }//calculate the matching
    ShortestAugmentingPath<GraphAA<int> > dk(sp, data, source, sink);
    //return edges with positive flow
    Vector<pair<int, int> > result;
    for(int i = 0; i < allowedMatches.getSize(); ++i)
        if(data[i].flow > 0) result.append(allowedMatches[i]);
    return result;
}

Vector<pair<int, int> > assignmentProblem(int n, int m,
    Vector<pair<pair<int, int>, double> > const& allowedMatches)
{//v = n + m + 2, e = n + m + allowedMatches.getSize()
    GraphAA<int> sp(n + m + 2);//setup graph and flow edges
    Vector<FlowData> data;
    for(int i = 0; i < allowedMatches.getSize(); ++i)
    {
        data.append(FlowData(allowedMatches[i].first.first, 1,
            allowedMatches[i].second));
        sp.addUndirectedEdge(allowedMatches[i].first.first,
            allowedMatches[i].first.second, i);
    }
    int source = n + m, sink = source + 1;//setup source and sink groups
    for(int i = 0; i < source; ++i)
    {
        int from = i, to = sink;
        if(i < n)
        {
            from = source;
            to = i;
        }
        data.append(FlowData(from, 1, 0));
        sp.addUndirectedEdge(from, to, i + allowedMatches.getSize());
    }//calculate the matching
    ShortestAugmentingPath<GraphAA<int> > dummy(sp, data, source, sink,
        min(n, m));
    //return edges with positive flow
    Vector<pair<int, int> > result;
    for(int i = 0; i < allowedMatches.getSize(); ++i)
        if(data[i].flow > 0) result.append(allowedMatches[i].first);
    return result;
}

}//end namespace
#endif
