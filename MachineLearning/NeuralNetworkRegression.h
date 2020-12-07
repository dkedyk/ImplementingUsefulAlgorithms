#ifndef IGMDK_NEURAL_NETWORK_REGRESSION_H
#define IGMDK_NEURAL_NETWORK_REGRESSION_H
#include "RegressionCommon.h"
#include "../Utils/Utils.h"
#include "../NumericalOptimization/NumericalOptimization.h"
#include "../NumericalOptimization/DiscreteGlobalOptimization.h"
#include "../RandomNumberGeneration/Statistics.h"
#include <cmath>

namespace igmdk{

class HiddenLayerNNReg
{
    Vector<NeuralNetwork> nns;
public:
    template<typename DATA> HiddenLayerNNReg(DATA const& data,
        Vector<double>const& p, int nGoal = 100000, int nNns = 5):
        nns(nNns, NeuralNetwork(getD(data), true, p[0]))
    {//structure
        int nHidden = p[1], D = getD(data),
            nRepeats = ceiling(nGoal, data.getSize());
        double a = sqrt(3.0/D);
        for(int l = 0; l < nns.getSize(); ++l)
        {
            NeuralNetwork& nn = nns[l];
            nn.addLayer(nHidden);
            for(int j = 0; j < nHidden; ++j)
                for(int k = -1; k < D; ++k)
                    nn.addConnection(0, j, k, k == -1 ? 0 :
                        GlobalRNG().uniform(-a, a));
            nn.addLayer(1);
            for(int k = -1; k < nHidden; ++k)
                nn.addConnection(1, 0, k, 0);
        }
        //training
        for(int j = 0; j < nRepeats; ++j)
            for(int i = 0; i < data.getSize(); ++i)
                learn(data.getX(i), data.getY(i));
    }
    void learn(NUMERIC_X const& x, double label)
    {
        for(int l = 0; l < nns.getSize(); ++l)
            nns[l].learn(x, Vector<double>(1, label));
    }
    double evaluate(NUMERIC_X const& x)const
    {
        double result = 0;
        for(int l = 0; l < nns.getSize(); ++l)
            result += nns[l].evaluate(x)[0];
        return result/nns.getSize();
    }
    int predict(NUMERIC_X const& x)const{return evaluate(x);}
};
struct NoParamsNNReg
{
    HiddenLayerNNReg model;
    template<typename DATA> static Vector<double> findParams(DATA const&
        data, int rLow = -15, int rHigh = 5, int hLow = 0, int hHigh = 6)
    {
        Vector<Vector<double> > sets(2);
        for(int i = rLow; i <= rHigh; i += 2) sets[0].append(pow(2, i));
        for(int i = hLow; i <= hHigh; i += 2) sets[1].append(pow(2, i));
        return gridMinimize(sets,
            RRiskFunctor<HiddenLayerNNReg, Vector<double>, DATA>(data));
    }
    template<typename DATA> NoParamsNNReg(DATA const& data):
        model(data, findParams(data)) {}
    double predict(NUMERIC_X const& x)const{return model.predict(x);}
};
typedef ScaledLearner<NoParamsLearner<NoParamsNNReg, double>, double, EMPTY,
    ScalerMQ> SNNReg;

}//end namespace
#endif

