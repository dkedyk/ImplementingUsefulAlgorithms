#ifndef IGMDK_KNN_H
#define IGMDK_KNN_H
#include "ClassificationCommon.h"
#include "../ComputationalGeometry/KDTree.h"
#include <cmath>
namespace igmdk{

template<typename X = NUMERIC_X, typename INDEX = VpTree<X, int, typename
    EuclideanDistance<X>::Distance> > class KNNClassifier
{
    mutable INDEX instances;
    int n, nClasses;
public:
    KNNClassifier(int theNClasses): nClasses(theNClasses), n(0) {}
    template<typename DATA> KNNClassifier(DATA const& data): n(0),
        nClasses(findNClasses(data))
    {
        for(int i = 0; i < data.getSize(); ++i)
            learn(data.getY(i), data.getX(i));
    }
    void learn(int label, X const& x){instances.insert(x, label); ++n;}
    int predict(X const& x)const
    {
        Vector<typename INDEX::NodeType*> neighbors =
            instances.kNN(x, 2 * int(log(n))/2 + 1);
        Vector<int> votes(nClasses);
        for(int i = 0; i < neighbors.getSize(); ++i)
            ++votes[neighbors[i]->value];
        return argMax(votes.getArray(), votes.getSize());
    }
};

}//end namespace
#endif

