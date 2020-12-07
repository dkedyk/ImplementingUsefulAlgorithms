#ifndef IGMDK_RANDOM_FOREST_H
#define IGMDK_RANDOM_FOREST_H
#include "ClassificationCommon.h"
#include "DecisionTree.h"
#include "../Utils/Vector.h"
#include "../RandomNumberGeneration/Random.h"
#include <cmath>
namespace igmdk{

class RandomForest
{
    Vector<DecisionTree> forest;
    int nClasses;
public:
    template<typename DATA> RandomForest(DATA const& data, int nTrees = 300):
        nClasses(findNClasses(data))
    {
        assert(data.getSize() > 1);
        for(int i = 0; i < nTrees; ++i)
        {
            PermutedData<DATA> resample(data);
            for(int j = 0; j < data.getSize(); ++j)
                resample.addIndex(GlobalRNG().mod(data.getSize()));
            forest.append(DecisionTree(resample, 0, true));
        }
    }
    template <typename ENSEMBLE> static int classifyWork(NUMERIC_X const& x,
        ENSEMBLE const& e, int nClasses)
    {
        Vector<int> counts(nClasses, 0);
        for(int i = 0; i < e.getSize(); ++i) ++counts[e[i].predict(x)];
        return argMax(counts.getArray(), counts.getSize());
    }
    int predict(NUMERIC_X const& x)const
        {return classifyWork(x, forest, nClasses);}
    Vector<double> classifyProbs(NUMERIC_X const& x)const
    {
        Vector<double> counts(nClasses, 0);
        for(int i = 0; i < forest.getSize(); ++i)
            ++counts[forest[i].predict(x)];
        normalizeProbs(counts);
        return counts;
    }
};

}//end namespace
#endif

