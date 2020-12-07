#ifndef IGMDK_SAMME_H
#define IGMDK_SAMME_H

#include "../MachineLearning/Classification.h"

namespace igmdk{

template<typename LEARNER = NoParamsLearner<DecisionTree, int>,
    typename PARAMS = EMPTY, typename X = NUMERIC_X> class AdaBoostSamme
{
    Vector<LEARNER> classifiers;
    Vector<double> weights;
    int nClasses;
public:
    template<typename DATA> AdaBoostSamme(DATA const& data, PARAMS const&
        p = PARAMS(), int nClassifiers = 300): nClasses(findNClasses(data))
    {
        int n = data.getSize();
        assert(n > 0 && nClassifiers > 0 && nClasses > 0);
        Vector<double> dataWeights(n, 1.0/n);
        for(int i = 0; i < nClassifiers; ++i)
        {
            AliasMethod sampler(dataWeights);
            PermutedData<DATA> resample(data);
            for(int j = 0; j < n; ++j) resample.addIndex(sampler.next());
            classifiers.append(LEARNER(resample, p));
            double error = 0;
            Bitset<> isWrong(n);
            for(int j = 0; j < n; ++j) if(classifiers.lastItem().predict(
                data.getX(j)) != data.getY(j))
                {
                    isWrong.set(j);
                    error += dataWeights[j];
                }
            if(error >= 1 - 1.0/nClasses) classifiers.removeLast();
            else if(error == 0)
            {//replace ensemble by classifier
                Vector<LEARNER> temp;
                temp.append(classifiers.lastItem());
                classifiers = temp;
                weights = Vector<double>(1, 1);
                break;
            }
            else
            {
                double expWeight = (nClasses - 1) * (1 - error)/error;
                weights.append(log(expWeight));
                for(int j = 0; j < n; ++j)
                    if(isWrong[j]) dataWeights[j] *= expWeight;
                normalizeProbs(dataWeights);
            }
        }
    }
    int predict(X const& x)const
    {
        Vector<double> counts(nClasses, 0);
        for(int i = 0; i < classifiers.getSize(); ++i)
            counts[classifiers[i].predict(x)] += weights[i];
        return argMax(counts.getArray(), counts.getSize());
    }
};

}//end namespace
#endif

