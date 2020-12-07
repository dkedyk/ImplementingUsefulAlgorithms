#ifndef IGMDK_TSP_H
#define IGMDK_TSP_H
#include "../RandomNumberGeneration/Random.h"
#include "../ComputationalGeometry/Point.h"
#include "../Graphs/Graph.h"
#include "SearchAlgorithms.h"
#include "MetaHeuristics.h"
namespace igmdk{

struct TSPRandomInstance
{
    Vector<Vector<double> > points;
    TSPRandomInstance(int n)
    {
        for(int i = 0; i < n; ++i)
        {
            Vector<double> point;
            for(int j = 0; j < 2; ++j) point.append(GlobalRNG().uniform01());
            points.append(point);
        }
    }
    double evalStep(int from, int to, bool isReturn = false)const
    {//evaluate single step
        if(isReturn) return 0;//no cost to return
        assert(from >= 0 && from < points.getSize() && to >= 0 &&
            to < points.getSize());
        EuclideanDistance<Vector<double> >::Distance d;
        return d(points[from], points[to]);
    }
    double operator()(Vector<int> const& permutation)const
    {//evaluate complete solution
        assert(permutation.getSize() == points.getSize());
        double result = 0;
        for(int i = 1; i < points.getSize(); ++i)
            result += evalStep(permutation[i - 1], permutation[i]);
        return result;
    }
};

template<typename PROBLEM> class TSPProxy
{
    PROBLEM const& p;
    double evalStep(Vector<int> const& permutation, int from, int to)const
        {return p.evalStep(permutation[from], permutation[to], to == 0);}
public:
    typedef double INCREMENTAL_STATE;//current score
    TSPProxy(PROBLEM const& theP): p(theP) {}
    INCREMENTAL_STATE updateIncrementalStateReverse(int i, int j,
        Vector<int> const& permutation, INCREMENTAL_STATE const& is)const
    {//take out edges i-1 to i and j to j + 1; add i-1 to j and i to j + 1
        int n = permutation.getSize();
        double im1Factor = i > 0 ? evalStep(permutation, i - 1, j) -
            evalStep(permutation, i - 1, i) : 0,
            jp1Factor = j + 1 < n ? evalStep(permutation, i, j + 1) -
            evalStep(permutation, j, j + 1) : 0;
        return is + im1Factor + jp1Factor;
    }
    double evalReverse(int i, int j, Vector<int> const& permutation,
        INCREMENTAL_STATE const& is)const
        {return updateIncrementalStateReverse(i, j, permutation, is) - is;}
    INCREMENTAL_STATE updateIncrementalStateSwap(int i, int j,
        Vector<int> const& permutation, INCREMENTAL_STATE const& is)const
    {//take out edges i-1 to i and i to i + 1, j-1 to j and j to j + 1
     //add edges i-1 to j and j to i + 1, j-1 to i and i to j + 1
        int n = permutation.getSize();
        double im1Factor = i > 0 ? evalStep(permutation, i - 1, j) -
                evalStep(permutation, i - 1, i) : 0,
            ip1Factor = i + 1 < n ? evalStep(permutation, j, i + 1) -
                evalStep(permutation, i, i + 1) : 0,
            jm1Factor = j > 0 ? evalStep(permutation, j - 1, i) -
                evalStep(permutation, j - 1, j) : 0,
            jp1Factor = j + 1 < n ? evalStep(permutation, i, j + 1) -
                evalStep(permutation, j, j + 1) : 0;
        return is + im1Factor + ip1Factor + jm1Factor + jp1Factor;
    }
    double evalSwap(int i, int j, Vector<int> const& permutation,
        INCREMENTAL_STATE const& is)const
        {return updateIncrementalStateSwap(i, j, permutation, is) - is;}
    INCREMENTAL_STATE getIncrementalState(Vector<int> const& permutation)const
        {return p(permutation);}
    double operator()(Vector<int> const& permutation)const
        {return getIncrementalState(permutation);}
};

template<typename PROBLEM> double findMSTCost(PROBLEM const& instance,
    Vector<int> const& remPoints)
{
    int n = remPoints.getSize();
    if(n <= 1) return 0;
    GraphAA<double> g(n);
    for(int i = 0; i < n; ++i)
        for(int j = i + 1; j < n; ++j) g.addUndirectedEdge(i, j,
            instance.evalStep(remPoints[i], remPoints[j]));
    assert(validateGraph(g));
    Vector<int> parents = MST(g);
    double sum = 0;
    for(int i = 0; i < parents.getSize(); ++i)
    {
        int parent = parents[i];
        if(parent != -1)
            sum += instance.evalStep(remPoints[i], remPoints[parent]);
    }
    return sum;
}
template<typename PROBLEM> class BranchAndBoundPermutation
{
    PROBLEM const& p;
    int const n;
public:
    typedef Vector<int> INCREMENTAL_STATE;
    INCREMENTAL_STATE getInitialState()const{return INCREMENTAL_STATE();}
    typedef Vector<int> X;
    BranchAndBoundPermutation(int theN, PROBLEM const& theProblem): n(theN),
        p(theProblem){}
    typedef int Move;
    bool isSolutionComplete(INCREMENTAL_STATE const& permutation)const
        {return permutation.getSize() == n;}
    double evaluate(INCREMENTAL_STATE const& permutation)const
    {
        assert(isSolutionComplete(permutation));
        return p(permutation);
    }
    X extractSolution(INCREMENTAL_STATE const& permutation)const
    {
        assert(isSolutionComplete(permutation));
        return permutation;
    }
    Vector<pair<double, Move> > generateMoves(
        INCREMENTAL_STATE const& permutation)const
    {
        double sumNow = 0;
        Vector<bool> isIncluded(n, false);
        for(int i = 0; i < permutation.getSize(); ++i)
        {
            isIncluded[permutation[i]] = true;
            if(i > 0) sumNow += p.evalStep(permutation[i - 1], permutation[i]);
        }
        Vector<pair<double, Move> > result;
        for(int i = 0; i < n; ++i)
            if(!isIncluded[i])
            {
                Vector<int> remainder;
                for(int j = 0; j < n; ++j) if(j != i && !isIncluded[j])
                    remainder.append(j);
                double lb = sumNow + (permutation.getSize() > 0 ?
                    p.evalStep(permutation.lastItem(), i) : 0) +
                    findMSTCost(p, remainder);
                result.append(pair<double, Move>(lb, i));
            }
        return result;
    }
    void move(INCREMENTAL_STATE& permutation, Move m)const
        {permutation.append(m);}
    void undoMove(INCREMENTAL_STATE& permutation, Move m)const
        {permutation.removeLast();}
};
template<typename INSTANCE> pair<Vector<int>, bool> solveTSPBranchAndBound(
    INSTANCE const& instance, int maxLowerBounds)
{
    return branchAndBound(BranchAndBoundPermutation<INSTANCE>(
        instance.points.getSize(), instance), maxLowerBounds);
}
template<typename INSTANCE> Vector<int> solveTSPRTAS(INSTANCE const&
    instance)
{
    return realtimeAStar(BranchAndBoundPermutation<INSTANCE>(
        instance.points.getSize(), instance));
}

template<typename PROBLEM> struct AStartTSPProblem
{
    PROBLEM const& problem;
    typedef pair<Bitset<>, int> STATE_ID;
    STATE_ID nullState()const
        {return make_pair(Bitset<>(problem.points.getSize()), -1);}
    struct HASHER
    {
        DataHash<> h;
        HASHER(unsigned long long m): h(m){}
        unsigned long long operator()(STATE_ID const& s)const
        {
            Vector<unsigned long long> storage = s.first.getStorage();
            storage.append(s.second);
            return h(storage);
        }
    };//don't know a legitimate best first node
    STATE_ID start()const{return nullState();}
    Vector<int> convertStatePath(Vector<STATE_ID> const& path)const
    {//first state start
        Vector<int> result;
        for(int i = 1; i < path.getSize(); ++i)
            result.append(path[i].second);
        return result;
    }
    Vector<int> findRemainder(STATE_ID const& id)const
    {
        Vector<int> result;
        for(int i = 0; i < id.first.getSize(); ++i)
            if(!id.first[i]) result.append(i);
        return result;
    }
    template<typename VISITOR>
    bool isGoal(STATE_ID id, VISITOR const& dummy)const
    {
        id.first.flip();
        return id.first.isZero();
    }
    template<typename VISITOR>
    Vector<STATE_ID> nextStates(STATE_ID const& id, VISITOR const& dummy)const
    {
        Vector<STATE_ID> result;
        Vector<int> remainder = findRemainder(id);
        for(int i = 0; i < remainder.getSize(); ++i)
        {
            STATE_ID to = id;
            to.first.set(remainder[i], true);
            to.second = remainder[i];
            result.append(to);
        }
        return result;
    }
    template<typename VISITOR> double remainderLowerBound(STATE_ID const&dummy,
        STATE_ID const& to, VISITOR const& dummy2)const
        {return findMSTCost(problem, findRemainder(to));}
    template<typename VISITOR> double distance(STATE_ID const& j,
        STATE_ID const& to, VISITOR const& dummy)const
        {return j == start() ? 0 : problem.evalStep(j.second, to.second);}
};

template<typename INSTANCE> pair<Vector<int>, bool> solveTSPAStar(INSTANCE
    const& instance, int maxSetSize)
{
    AStartTSPProblem<INSTANCE> p = {instance};
    pair<Vector<typename AStartTSPProblem<INSTANCE>::STATE_ID>, bool> result =
        AStar<AStartTSPProblem<INSTANCE> >::solve(p, maxSetSize);
    return make_pair(p.convertStatePath(result.first), result.second);
}

template<typename INSTANCE> pair<Vector<int>, bool> solveTSPRBFS(INSTANCE
    const& instance, int maxLowerBounds)
{
    AStartTSPProblem<INSTANCE> p = {instance};
    RecursiveBestFirstSearch<AStartTSPProblem<INSTANCE> > rbfs(p,
        maxLowerBounds);
    return make_pair(p.convertStatePath(rbfs.bestPath), rbfs.foundGoal);
}

}//end namespace
#endif
