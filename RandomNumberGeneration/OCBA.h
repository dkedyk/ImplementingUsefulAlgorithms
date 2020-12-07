#ifndef IGMDK_OCBA_H
#define IGMDK_OCBA_H

#include "../Utils/Vector.h"
#include "Statistics.h"
#include "MultipleComparison.h"
#include <cmath>
namespace igmdk{

bool isNormal0BestBonf(Vector<NormalSummary> const& data, double aLevel)
{//smallest is best with precision meanPrecision
    int k = data.getSize();
    assert(k > 1);
    double z = find2SidedConfZBonf(k, 1 - aLevel),
        upper = data[0].mean + z * data[0].stddev();
    for(int i = 1; i < k; ++i)
    {
        double lower = data[i].mean - z * data[i].stddev();
        if(lower <= upper) return false;
    }
    return true;
}

template<typename MULTI_FUNCTION> struct OCBA
{
    MULTI_FUNCTION const& f;
    Vector<IncrementalStatistics> data;
    int nDone;
    OCBA(MULTI_FUNCTION const& theF = MULTI_FUNCTION(), int initialSims = 30):
        f(theF), data(theF.getSize())
    {
        int k = f.getSize();
        for(int i = 0; i < k; ++i)
            for(int j = 0; j < initialSims; ++j) data[i].addValue(f(i));
        nDone = k * initialSims;
    }
    pair<Vector<NormalSummary>, int> findBest()
    {
        int k = f.getSize();
        Vector<NormalSummary> s;
        for(int i = 0; i < k; ++i) s.append(data[i].getStandardErrorSummary());
        int bestI = 0, bestRatioI = -1;
        double bestMean = s[0].mean, ratioSum = 0, bestRatio;
        for(int i = 1; i < k; ++i)
            if(s[i].mean < bestMean) bestMean = s[bestI = i].mean;
        swap(s[0], s[bestI]);
        return make_pair(s, bestI);
    }
    void simulateNext()
    {
        pair<Vector<NormalSummary>, int> best = findBest();
        int k = f.getSize(), bestI = best.second, bestRatioI = -1;;
        Vector<NormalSummary> s = best.first;
        //compute the largest OCBA ratio
        double bestMean = s[0].mean, ratioSum = 0, bestRatio;
        for(int i = 1; i < k; ++i)
        {
            double meanDiff = s[i].mean - bestMean, ratio =
                s[i].variance/(meanDiff * meanDiff);
            ratioSum += ratio * ratio/s[i].variance;
            if(bestRatioI == -1 || ratio > bestRatio)
            {
                bestRatio = ratio;
                bestRatioI = i;
            }
        }
        double ratioBest = sqrt(ratioSum * s[0].variance);
        if(ratioBest > bestRatio) bestRatioI = bestI;
        else if(bestRatioI == bestI) bestRatioI = 0;
        //simulate the largest ratio alternative
        data[bestRatioI].addValue(f(bestRatioI));
        ++nDone;
    }
    int simulateTillBest(int simBudget = 100000, double aLevel = 0.05)
    {
        assert(nDone < simBudget);
        int k = f.getSize(), nTests = lgCeiling(simBudget) - lgFloor(nDone);
        while(nDone < simBudget)
        {
            simulateNext();
            if(isPowerOfTwo(nDone) || nDone == simBudget - 1)
            {
                Vector<NormalSummary> s = findBest().first;
                if(isNormal0BestBonf(s, aLevel/nTests)) break;
            }
        }
        return nTests;
    }
};

}//end namespace
#endif
