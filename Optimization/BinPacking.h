#ifndef IGMDK_BINPACKING_H
#define IGMDK_BINPACKING_H
#include "../RandomNumberGeneration/Random.h"
#include "SearchAlgorithms.h"
#include "MetaHeuristics.h"
namespace igmdk{

class BinPackingRandomInstance
{
    Vector<double> weights;
    double binSize;
public:
    double getWeight(int i)const
    {
        assert(i >= 0 && i < weights.getSize());
        return weights[i];
    }
    double getBinSize()const{return binSize;}
    BinPackingRandomInstance(int n): binSize(5)//uniformly random weights
        {for(int i = 0; i < n; ++i) weights.append(GlobalRNG().uniform01());}
};

template<typename PROBLEM> class BinPackingProxy
{
    PROBLEM const& p;
    double getBinScore(Vector<int> const& bin)const
    {//evaluate complete solution as minimization problem
        if(bin.getSize() == 0) return 0;
        double sum = 0;
        for(int i = 0; i < bin.getSize(); ++i)
            sum += p.getWeight(bin[i]);
        double score = 1, binSize = p.getBinSize();
        if(sum > binSize)//simple penalty for exceeding capacity
            score += exp((sum - binSize)/binSize);
        return score;
    }
public:
    typedef double INCREMENTAL_STATE;//current score
    BinPackingProxy(PROBLEM const& theP): p(theP) {}
    INCREMENTAL_STATE updateIncrementalState(int index, int from, int to,
        Vector<Vector<int> >const& bins, INCREMENTAL_STATE const& is)const
    {
        assert(index < bins[from].getSize());
        Vector<int> &bin1 = (Vector<int>&)bins[from],
            &bin2 = (Vector<int>&)bins[to];
        double partialScore = is - getBinScore(bin1) - getBinScore(bin2);
        //do move
        bin2.append(bin1[index]);
        bin1[index] = bin1.lastItem();
        bin1.removeLast();
        double score = partialScore+ getBinScore(bin1) + getBinScore(bin2);
        //undo move
        bin1.append(bin2.lastItem());
        bin2.removeLast();
        swap(bin1[index], bin1.lastItem());
        return score;
    }
    double evalStep(int index, int from, int to, Vector<Vector<int> >const&
        bins, INCREMENTAL_STATE const& is)const
        {return updateIncrementalState(index, from, to, bins, is) - is;}
    INCREMENTAL_STATE getIncrementalState(Vector<Vector<int> >const& bins)
        const
    {
        double sum = 0;
        for(int i = 0; i < bins.getSize(); ++i) sum += getBinScore(bins[i]);
        return sum;
    }
    double operator()(Vector<Vector<int> > const& bins)const
        {return getIncrementalState(bins);}
};

template<typename PROBLEM>
Vector<Vector<int> > firstFitDecreasing(PROBLEM const& p, int n)
{
    assert(n > 0);
    Vector<pair<double, int> > items(n);
    double binSize = p.getBinSize();
    for(int i = 0; i < n; ++i)
    {
        items[i].first = p.getWeight(i);
        assert(items[i].first <= binSize);//else bad problem
        items[i].second = i;
    }
    quickSort(items.getArray(), 0, n - 1,
        PairFirstComparator<double, int, ReverseComparator<double> >());
    Vector<Vector<int> > result(1);//start with one empty bin
    Vector<double> sums(1, 0);
    for(int i = 0; i < n; ++i)
    {
        bool fit = false;
        for(int bin = 0; bin < result.getSize(); ++bin)
            if(sums[bin] + items[i].first <= binSize)
            {
                fit = true;
                result[bin].append(items[i].second);
                sums[bin] += items[i].first;
                break;
            }
            else if(bin == result.getSize() - 1)
            {//doesn't fit in last bin, add a new one
                result.append(Vector<int>());
                sums.append(0);
            }
    }
    return result;
}

}//end namespace
#endif
