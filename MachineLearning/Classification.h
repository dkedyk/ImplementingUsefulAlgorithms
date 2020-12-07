#ifndef IGMDK_MACHINELEARNING_H
#define IGMDK_MACHINELEARNING_H
#include "ClassificationCommon.h"
#include "RandomForest.h"
#include "KernelSVM.h"
namespace igmdk{

template<typename SUBSET_LEARNER = RandomForest> struct SmartFSLearner
{
    typedef FeatureSubsetLearner<SUBSET_LEARNER> MODEL;
    MODEL model;
public:
    template<typename DATA> SmartFSLearner(DATA const& data, int limit = 20):
        model(data, selectFeaturesSmart(SCVRiskFunctor<MODEL, Bitset<>,DATA>(
        data), getD(data), limit)) {}
    int predict(NUMERIC_X const& x)const{return model.predict(x);}
};

class SimpleBestCombiner
{
    BestCombiner<int> c;
public:
    template<typename DATA> SimpleBestCombiner(DATA const& data)
    {
        c.addNoParamsClassifier<RandomForest>(data, SCVRiskFunctor<
            NoParamsLearner<RandomForest, int>, EMPTY, DATA>(data));
        c.addNoParamsClassifier<SSVM>(data, SCVRiskFunctor<
            NoParamsLearner<SSVM, int>, EMPTY, DATA>(data));
    }
    int predict(NUMERIC_X const& x)const{return c.predict(x);}
};

}//end namespace
#endif

