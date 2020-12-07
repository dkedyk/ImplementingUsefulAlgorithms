#ifndef IGMDK_RANDOM_FOREST_REGRESSION_H
#define IGMDK_RANDOM_FOREST_REGRESSION_H
#include "LearningCommon.h"
#include "../Utils/Vector.h"
#include "../RandomNumberGeneration/Statistics.h"
#include "RegressionTree.h"

namespace igmdk{

class RandomForestReg
{
    Vector<RegressionTree> forest;
public:
    template<typename DATA> RandomForestReg(DATA const& data,
        int nTrees = 300){addTrees(data, nTrees);}
    template<typename DATA> void addTrees(DATA const& data, int nTrees)
    {
        assert(data.getSize() > 1);
        for(int i = 0, D = getD(data); i < nTrees; ++i)
        {
            PermutedData<DATA> resample(data);
            for(int j = 0; j < data.getSize(); ++j)
                resample.addIndex(GlobalRNG().mod(data.getSize()));
            forest.append(RegressionTree(resample, 0, 50, true));
        }
    }
    double predict(NUMERIC_X const& x)const
    {
        IncrementalStatistics s;
        for(int i = 0; i < forest.getSize(); ++i)
            s.addValue(forest[i].predict(x));
        return s.getMean();
    }
};

}//end namespace
#endif

