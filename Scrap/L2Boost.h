#ifndef IGMDK_L2_BOOST_H
#define IGMDK_L2_BOOST_H
#include "../MachineLearning/Regression.h"
#include <cmath>

namespace igmdk{

template<typename LEARNER = NoParamsLearner<RegressionTree, double>,
    typename PARAMS = EMPTY, typename X = NUMERIC_X> class L2Boost
{
    Vector<LEARNER> classifiers;
    Vector<double> weights;
    double getWeight(int i)const{return 1/pow(i + 1, 0.501);}
    struct L2Loss
    {
        Vector<double> F;
        L2Loss(int n): F(n, 0) {}
        double getNegGrad(int i, double y){return 2 * (y - F[i]);}
        double loss(int i, double y){return (F[i] - y) * (F[i] - y);}
    };
public:
    template<typename DATA> L2Boost(DATA const& data,
        PARAMS const& p = PARAMS(), int nClassifiers = 300)
    {
        int n = data.getSize();
        assert(n > 0 && nClassifiers > 0);
        L2Loss l(n);
        RelabeledData<DATA> regData(data);
        for(int j = 0; j < n; ++j) regData.addLabel(data.getY(j));
        for(int i = 0; i < nClassifiers; ++i)
        {//find gradients of relabel data and fit learner
            for(int j = 0; j < n; ++j)
                regData.labels[j] = l.getNegGrad(j, data.getY(j));
            classifiers.append(LEARNER(regData, p));
            Vector<double> h;
            for(int j = 0; j < n; ++j)
                h.append(classifiers.lastItem().predict(data.getX(j)));
            //calculate weights and update F
            double sumH2 = 0, weight = 0;
            for(int j = 0; j < n; ++j)
            {
                sumH2 += h[j] * h[j];
                weight += (data.getY(j) - l.F[j]) * h[j];
            }
            if(weight > 0 && isfinite(weight/sumH2))
            {
                weights.append(weight/sumH2);
                for(int j = 0; j < n; ++j)
                    l.F[j] += weights.lastItem() * h[j];
            }
            else
            {
                classifiers.removeLast();
                break;
            }
        }
    }
    double predict(X const& x)const
    {
        double sum = 0;
        for(int i = 0; i < classifiers.getSize(); ++i)
            sum += classifiers[i].predict(x) * weights[i];
        return sum;
    }
};

}//end namespace
#endif

