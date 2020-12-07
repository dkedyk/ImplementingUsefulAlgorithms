#ifndef IGMDK_REGRESSION_COMMON_H
#define IGMDK_REGRESSION_COMMON_H
#include "LearningCommon.h"
#include "../Utils/Vector.h"
#include "../RandomNumberGeneration/Random.h"
#include <cmath>

namespace igmdk{

struct RegressionStats
{
    double expStd, rmse, l1Err, lInfErr;
    void debug()const
    {
        DEBUG(expStd);
        DEBUG(rmse);
        DEBUG(l1Err);
        DEBUG(lInfErr);
    }
};
RegressionStats evaluateRegressor(
    Vector<pair<double, double> > const& testResult)
{
    IncrementalStatistics yStats, l2Stats, l1Stats;
    for(int i = 0; i < testResult.getSize(); ++i)
    {
        yStats.addValue(testResult[i].first);
        double diff = testResult[i].second - testResult[i].first;
        l1Stats.addValue(abs(diff));
        l2Stats.addValue(diff * diff);
    }
    RegressionStats result;
    result.lInfErr = l1Stats.maximum;
    result.l1Err = l1Stats.getMean();
    result.rmse = sqrt(l2Stats.getMean());
    result.expStd = 1 - result.rmse/yStats.stdev();
    return result;
}
template<typename LEARNER, typename DATA, typename PARAMS> double
    crossValidateReg(PARAMS const& p, DATA const& data, int nFolds = 5)
{
    return evaluateRegressor(crossValidateGeneral<LEARNER,
        typename DATA::Y_TYPE>(p, data, nFolds)).rmse;
}
template<typename LEARNER, typename PARAM, typename DATA>
struct RRiskFunctor
{
    DATA const& data;
    RRiskFunctor(DATA const& theData): data(theData) {}
    double operator()(PARAM const& p)const
        {return crossValidateReg<LEARNER>(p, data);}
};

template<typename LEARNER, typename DATA, typename PARAMS> double
    repeatedCVReg(PARAMS const& p, DATA const& data, int nFolds = 5,
    int nRepeats = 5)
{
    return evaluateRegressor(repeatedCVGeneral<double>(
        LEARNER(data, p), data, nFolds, nRepeats)).rmse;
}
template<typename LEARNER, typename PARAM, typename DATA>
struct RRCVRiskFunctor
{
    DATA const& data;
    RRCVRiskFunctor(DATA const& theData): data(theData) {}
    double operator()(PARAM const& p)const
        {return repeatedCVReg<LEARNER>(p, data);}
};

}//end namespace
#endif

