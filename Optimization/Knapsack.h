#ifndef IGMDK_KNAPSACK_H
#define IGMDK_KNAPSACK_H
#include "../RandomNumberGeneration/Random.h"
#include "SearchAlgorithms.h"
#include "MetaHeuristics.h"
namespace igmdk{

class KnapsackRandomInstance
{
    Vector<double> profits, weights;
    double capacity;
public:
    double getProfit(int i)const
    {
        assert(i >= 0 && i < profits.getSize());
        return profits[i];
    }
    double getWeight(int i)const
    {
        assert(i >= 0 && i < weights.getSize());
        return weights[i];
    }
    int getSize()const{return profits.getSize();}
    double getCapacity()const{return capacity;}
    KnapsackRandomInstance(int n): capacity(0)
    {//uniformly random weights and profits and capacity half of weights
        for(int i = 0; i < n; ++i)
        {
            profits.append(GlobalRNG().uniform01());
            weights.append(GlobalRNG().uniform01());
            capacity += weights[i];
        }
        capacity /= 2;
    }
    double operator()(Vector<bool> const& subset)const
    {
        assert(subset.getSize() == getSize());
        double totalProfit = 0, totalWeight = 0;
        for(int i = 0; i < subset.getSize(); ++i) if(subset[i])
        {
             totalProfit += getProfit(i);
             totalWeight += getWeight(i);
        }//infeasible marked by infinite score
        return totalWeight <= capacity ? -totalProfit :
            numeric_limits<double>::infinity();
    }
};

template<typename PROBLEM> class KnapsackPenaltyProxy
{
    PROBLEM const& p;
    double profitLimit;
    double getScore(pair<double, double> const& is)const
    {//evaluate complete solution as minimization problem
        double score = -is.first, capacity = p.getCapacity();
        if(is.second > capacity)//simple penalty for exceeding capacity
            score += profitLimit * exp((is.second - capacity)/capacity);
        return score;
    }
public:
    typedef pair<double, double> INCREMENTAL_STATE;//first is profit
    KnapsackPenaltyProxy(PROBLEM const& theP): p(theP), profitLimit(0)
        {for(int i = 0; i < p.getSize(); ++i) profitLimit += p.getProfit(i);}
    INCREMENTAL_STATE updateIncrementalState(int iToFlip,
        Vector<bool> const& subset, INCREMENTAL_STATE const& is)const
    {
        bool selection = !subset[iToFlip];
        return INCREMENTAL_STATE(is.first + p.getProfit(iToFlip) *
            (selection ? 1 : -1), is.second + p.getWeight(iToFlip) *
            (selection ? 1 : -1));
    }
    double evalStep(int iToFlip, Vector<bool> const& subset,
        INCREMENTAL_STATE const& is)const
    {//evaluate single step difference
        return getScore(updateIncrementalState(iToFlip, subset, is)) -
            getScore(is);
    }
    INCREMENTAL_STATE getIncrementalState(Vector<bool> const& subset)const
    {
        assert(subset.getSize() == p.getSize());
        double totalProfit = 0, totalWeight = 0;
        for(int i = 0; i < subset.getSize(); ++i) if(subset[i])
        {
             totalProfit += p.getProfit(i);
             totalWeight += p.getWeight(i);
        }
        return INCREMENTAL_STATE(totalProfit, totalWeight);
    }
    double operator()(Vector<bool> const& subset)const
        {return getScore(getIncrementalState(subset));}
};

template<typename PROBLEM>
Vector<pair<double, int> > getKnapsackRatios(PROBLEM const& p)
{
    int n = p.getSize();
    double capacity = p.getCapacity();
    Vector<pair<double, int> > ratios(n);
    for(int i = 0; i < n; ++i)
        ratios[i] = make_pair(p.getProfit(i)/p.getWeight(i), i);
    quickSort(ratios.getArray(), 0, n - 1,
        PairFirstComparator<double, int, ReverseComparator<double> >());
    return ratios;
}
template<typename PROBLEM>
Vector<bool> approximateKnapsackGreedy(PROBLEM const& p)
{
    double remainingCapacity = p.getCapacity();
    Vector<pair<double, int> > ratios = getKnapsackRatios(p);
    Vector<bool> solution(ratios.getSize(), false);
    for(int i = 0; i < ratios.getSize(); ++i)
    {
        int j = ratios[i].second;
        if(remainingCapacity > p.getWeight(j))
        {
            remainingCapacity -= p.getWeight(j);
            solution[j] = true;
        }
    }
    return solution;
}

template<typename PROBLEM> class BranchAndBoundKnapsack
{
    PROBLEM const& p;//order 1
    int const n;//order 2
    Vector<pair<double, int> > const ratios;
public:
    struct INCREMENTAL_STATE
    {
        Vector<bool> subset;
        int current;
        double currentTotalProfit, currentTotalWeight;
    };
    INCREMENTAL_STATE getInitialState()const
    {
        INCREMENTAL_STATE is = {Vector<bool>(n, false), 0, 0, 0};
        return is;
    }
    typedef Vector<bool> X;
    BranchAndBoundKnapsack(PROBLEM const& theProblem): p(theProblem),
        n(p.getSize()), ratios(getKnapsackRatios(p)) {}
    typedef pair<int, bool> Move;
    bool isSolutionComplete(INCREMENTAL_STATE const& is)const
        {return is.current == n;}
    double evaluate(INCREMENTAL_STATE const& is)const
    {
        assert(isSolutionComplete(is));
        return p(is.subset);
    }
    X extractSolution(INCREMENTAL_STATE const& is)const
    {
        assert(isSolutionComplete(is));
        return is.subset;
    }
    Vector<pair<double, Move> > generateMoves(INCREMENTAL_STATE const& is)const
    {
        Vector<pair<double, Move> > result;
        for(int i = 0; i < 2; ++i)
        {
            double lb = is.currentTotalProfit,
                remainingCapacity = p.getCapacity() - is.currentTotalWeight;
            int j = ratios[is.current].second;
            if(i)//current move
            {
                remainingCapacity -= p.getWeight(j);
                lb += p.getProfit(j);
            }
            if(remainingCapacity >= 0)//remaining moves if any
            {
                for(int k = is.current + 1; k < n; ++k)
                {
                    int j = ratios[k].second;
                    if(remainingCapacity >= p.getWeight(j))
                    {
                        remainingCapacity -= p.getWeight(j);
                        lb += p.getProfit(j);
                    }
                    else
                    {//reached non-fitting item, use fractional
                        lb += remainingCapacity * ratios[k].first;
                        break;
                    }
                }
                result.append(pair<double, Move>(-lb, make_pair(j, i)));
            }
        }
        return result;
    }
    void move(INCREMENTAL_STATE& is, Move const& m)const
    {
        if(m.second)
        {
            is.subset[m.first] = m.second;
            is.currentTotalProfit += p.getProfit(m.first);
            is.currentTotalWeight += p.getWeight(m.first);
        }
        ++is.current;
    }
    void undoMove(INCREMENTAL_STATE& is, Move const& m)const
    {
        if(m.second)
        {
            is.subset[m.first] = false;
            is.currentTotalProfit -= p.getProfit(m.first);
            is.currentTotalWeight -= p.getWeight(m.first);
        }
        --is.current;
    }
};
template<typename INSTANCE> Vector<bool> solveKnapsackBranchAndBound(
    INSTANCE const& instance, int maxLowerBounds)
{
    return branchAndBound(BranchAndBoundKnapsack<INSTANCE>(instance),
        maxLowerBounds).first;
}

}//end namespace
#endif
