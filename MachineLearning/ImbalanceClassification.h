#ifndef IGMDK_IMBALANCE_CLASSIFICATION_H
#define IGMDK_IMBALANCE_CLASSIFICATION_H
#include "ClassificationCommon.h"
#include "RandomForest.h"
#include "KernelSVM.h"
#include "../Utils/Vector.h"
#include "../RandomNumberGeneration/Statistics.h"
#include <cmath>
namespace igmdk{

class WeightedRF
{
    Vector<DecisionTree> forest;
    int nClasses;
public:
    template<typename DATA> WeightedRF(DATA const& data, Vector<double> const
        & weights, int nTrees = 300): nClasses(findNClasses(data))
    {
        assert(data.getSize() > 1);
        AliasMethod sampler(weights);
        for(int i = 0; i < nTrees; ++i)
        {
            PermutedData<DATA> resample(data);
            for(int j = 0; j < data.getSize(); ++j)
                resample.addIndex(sampler.next());
            forest.append(DecisionTree(resample, 0, true));
        }
    }
    int predict(NUMERIC_X const& x)const
        {return RandomForest::classifyWork(x, forest, nClasses);}
};

template<typename DATA>
Vector<double> findImbalanceWeights(DATA const& data)
{
    int n = data.getSize(), properK = 0, nClasses = findNClasses(data);
    Vector<double> counts(nClasses);
    for(int i = 0; i < n; ++i) ++counts[data.getY(i)];
    for(int i = 0; i < nClasses; ++i) if(counts[i] > 0) ++properK;
    Vector<double> dataWeights(n, 0);
    for(int i = 0; i < data.getSize(); ++i)
        dataWeights[i] = 1.0/properK/counts[data.getY(i)];
    return dataWeights;
}
class ImbalanceRF
{
    WeightedRF model;
public:
    template<typename DATA> ImbalanceRF(DATA const& data, int nTrees = 300):
        model(data, findImbalanceWeights(data), nTrees) {}
    int predict(NUMERIC_X const& x)const{return model.predict(x);}
};

template<typename LEARNER, typename PARAMS = EMPTY>
    class WeightedBaggedLearner
{
    Vector<LEARNER> models;
    int nClasses;
public:
    template<typename DATA> WeightedBaggedLearner(DATA const& data,
        Vector<double> weights, PARAMS const& p = PARAMS(), int nBags = 15):
        nClasses(findNClasses(data))
    {
        assert(data.getSize() > 1);
        AliasMethod sampler(weights);
        for(int i = 0; i < nBags; ++i)
        {
            PermutedData<DATA> resample(data);
            for(int j = 0; j < data.getSize(); ++j)
                resample.addIndex(sampler.next());
            models.append(LEARNER(resample, p));
        }
    }
    int predict(NUMERIC_X const& x)const
        {return RandomForest::classifyWork(x, models, nClasses);}
};

class ImbalanceSVM
{
    WeightedBaggedLearner<MulticlassSVM<>,
        pair<GaussianKernel, double> > model;
public:
    template<typename DATA> ImbalanceSVM(DATA const& data): model(data,
        findImbalanceWeights(data), NoParamsSVM::gaussianMultiClassSVM(data))
        {}
    int predict(NUMERIC_X const& x)const{return model.predict(x);}
};
typedef ScaledLearner<NoParamsLearner<ImbalanceSVM, int>, int> SImbSVM;

}//end namespace
#endif

