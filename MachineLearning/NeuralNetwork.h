#ifndef IGMDK_NEURAL_NETWORK_H
#define IGMDK_NEURAL_NETWORK_H
#include "ClassificationCommon.h"
#include "../Utils/Vector.h"
#include "../NumericalOptimization/NumericalOptimization.h"
#include "../RandomNumberGeneration/Random.h"
#include <cmath>
namespace igmdk{

class BinaryNN
{
    Vector<NeuralNetwork> nns;
    void setupStructure(int D, int nHidden)
    {
        double a = sqrt(3.0/D);
        for(int l = 0; l < nns.getSize(); ++l)
        {
            NeuralNetwork& nn = nns[l];
            nn.addLayer(nHidden);
            for(int j = 0; j < nHidden; ++j)
                for(int k = -1; k < D; ++k) nn.addConnection(0, j, k,
                    k == -1 ? 0 : GlobalRNG().uniform(-a, a));
            nn.addLayer(1);
            for(int k = -1; k < nHidden; ++k) nn.addConnection(1, 0, k, 0);
        }
    }
public:
    BinaryNN(int D, int nHidden = 5, int nNns = 5):
        nns(nNns, NeuralNetwork(D)){setupStructure(D, nHidden);}
    template<typename DATA> BinaryNN(DATA const& data, int nHidden = 5, int
        nGoal = 100000, int nNns = 5): nns(nNns, NeuralNetwork(getD(data)))
    {
        int D = getD(data), nRepeats = ceiling(nGoal, data.getSize());
        setupStructure(D, nHidden);
        for(int j = 0; j < nRepeats; ++j)
            for(int i = 0; i < data.getSize(); ++i)
                learn(data.getX(i), data.getY(i));
    }
    void learn(NUMERIC_X const& x, int label)
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
    int predict(NUMERIC_X const& x)const{return evaluate(x) > 0.5;}
};

class MulticlassNN
{
    MulticlassLearner<NoParamsLearner<BinaryNN, int>, EMPTY> model;
public:
    template<typename DATA> MulticlassNN(DATA const& data): model(data) {}
    int predict(NUMERIC_X const& x)const {return model.classifyByProbs(x);}
};
typedef ScaledLearner<NoParamsLearner<MulticlassNN, int>, int, EMPTY,
    ScalerMQ> SNN;

class SOnlineNN
{
    ScalerMQ s;
    OnlineMulticlassLearner<BinaryNN, int> model;
public:
    template<typename DATA> SOnlineNN(DATA const& data): model(getD(data)),
        s(getD(data))
    {
        for(int j = 0; j < 1000000; ++j)
        {
            int i = GlobalRNG().mod(data.getSize());
            learn(data.getX(i), data.getY(i));
        }
    }
    SOnlineNN(int D): model(D), s(D) {}
    void learn(NUMERIC_X const& x, int label)
    {
        s.addSample(x);
        model.learn(s.scale(x), label);
    }
    int predict(NUMERIC_X const& x)const{return model.predict(s.scale(x));}
};

/*
class DeepNN
{
    int D, nClasses;
    DeepNeuralNetwork nn;
public:
    void learnUnsupI(NUMERIC_X const& x, int i, int n){nn.learnUnsupI(x, i, n);}
    template<typename MODEL> static int findNHidden(
        typename DATA<>::LABELED_DATA const& data)
    {
        int nLow = 1, nHigh = 6;
        Vector<int> sizes;
        for(int j = nLow; j <= nHigh; ++j)
        {
            int size = pow(2, j) + 1;
            sizes.append(size);
        }
        return valMinFunc(sizes.getArray(), sizes.getSize(),
            SCVRiskFunctor<MODEL, int>(data));
    }
    DeepNN(DATA<>::LABELED_DATA const& data, int nHidden,
        int nGoal = 1000000, double l = 0): D(data[0].first.getSize()),
        nClasses(findNClasses<NUMERIC_X>(data)), nn(D, nClasses, nHidden, l)
    {
        int nTrains = max(nGoal, data.getSize());
        for(int i = 0; i < nn.nHiddenLayers(); ++i)
        {
            for(int k = 0; k < nTrains; ++k)
            {
                int k = GlobalRNG().mod( data.getSize());
                learnUnsupI(data[k].first, i, data.getSize());
            }
        }
        while(nTrains--)
        {
            int i = GlobalRNG().mod( data.getSize());
            learn(data[i].first, data[i].second, data.getSize());
        }
    }
    void learn(NUMERIC_X const& x, int label, int n)
    {
        assert(label >= 0 && label < nClasses);
        Vector<double> result(nClasses);
        result[label] = 1;
        nn.learn(x, result, n);
    }
    int predict(NUMERIC_X const& x)const
    {
        Vector<double> result = nn.evaluate(x);
        return argMax(result.getArray(), result.getSize());
    }
};
class NoParamDeepNN
{
    DeepNN model;
public:
    NoParamDeepNN(DATA<>::LABELED_DATA const& data): model(data, DeepNN::findNHidden<DeepNN>(data)){}
    int predict(NUMERIC_X const& x)const {return model.predict(x);}
};*/

}//end namespace
#endif

