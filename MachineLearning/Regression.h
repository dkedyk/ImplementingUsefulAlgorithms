#ifndef IGMDK_REGRESSION_H
#define IGMDK_REGRESSION_H
#include "LearningCommon.h"
#include "RandomForestRegression.h"
#include "Lasso.h"
#include <cmath>

namespace igmdk{

template<typename SUBSET_LEARNER = RandomForestReg> struct SmartFSLearnerReg
{
    typedef FeatureSubsetLearner<SUBSET_LEARNER> MODEL;
    MODEL model;
public:
    template<typename DATA> SmartFSLearnerReg(DATA const& data,
        int subsampleLimit = 20): model(data, selectFeaturesSmart(
        RRiskFunctor<MODEL, Bitset<>, DATA>(data), getD(data),
        subsampleLimit)) {}
    double predict(NUMERIC_X const& x)const{return model.predict(x);}
};

class SimpleBestCombinerReg
{
    BestCombiner<double> c;
public:
    template<typename DATA> SimpleBestCombinerReg(DATA const& data)
    {
        c.addNoParamsClassifier<RandomForestReg>(data, RRiskFunctor<
            NoParamsLearner<RandomForestReg, double>, EMPTY, DATA>(data));
        c.addNoParamsClassifier<SLasso>(data, RRiskFunctor<
            NoParamsLearner<SLasso, double>, EMPTY, DATA>(data));
    }
    double predict(NUMERIC_X const& x)const{return c.predict(x);}
};

}//end namespace
#endif

