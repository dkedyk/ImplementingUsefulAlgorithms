#ifndef IGMDK_KNN_REGRESSION_H
#define IGMDK_KNN_REGRESSION_H
#include "LearningCommon.h"
#include "../ComputationalGeometry/KDTree.h"
#include "../RandomNumberGeneration/Statistics.h"
#include <cmath>

namespace igmdk{

template<typename X = NUMERIC_X, typename INDEX = VpTree<X, double, typename
    EuclideanDistance<X>::Distance> > class KNNReg
{
    mutable INDEX instances;
    int k;
public:
    template<typename DATA> KNNReg(DATA const& data, int theK = -1): k(theK)
    {
        assert(data.getSize() > 0);
        if(k == -1) k = 2 * int(log(data.getSize())/2) + 1;
        for(int i = 0; i < data.getSize(); ++i)
            learn(data.getY(i), data.getX(i));
    }
    void learn(double label, X const& x){instances.insert(x, label);}
    double predict(X const& x)const
    {
        Vector<typename INDEX::NodeType*> neighbors = instances.kNN(x, k);
        IncrementalStatistics s;
        for(int i = 0; i < neighbors.getSize(); ++i)
            s.addValue(neighbors[i]->value);
        return s.getMean();
    }
};
typedef ScaledLearner<NoParamsLearner<KNNReg<>, double>, double> SKNNReg;

}//end namespace
#endif

