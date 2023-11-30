#ifndef IGMDK_LEARNINGCOMMON_H
#define IGMDK_LEARNINGCOMMON_H
#include "../Utils/Utils.h"
#include "../HashTable/LinearProbingHashTable.h"
#include "../ComputationalGeometry/KDTree.h"
#include "../ComputationalGeometry/Point.h"
#include "../NumericalMethods/Matrix.h"
#include "../NumericalMethods/NumericalMethods.h"
#include "../NumericalOptimization/NumericalOptimization.h"
#include "../RandomTreap/LCPTreap.h"
#include "../RandomNumberGeneration/Statistics.h"
#include "../ExternalMemoryAlgorithms/EMVector.h"
#include <cmath>
#include <tgmath.h>
namespace igmdk{

typedef Vector<double> NUMERIC_X;
template<typename DATA> int getD(DATA const& data)
{
    assert(data.getSize() > 0);
    return data.getX(0).getSize();
}
/*
template<typename X> struct InMemoryUnlabelledData
{
    typedef X X_TYPE;
    typedef X Y_TYPE;
};

template<typename X> struct DiskUnlabelledData
{
    typedef X X_TYPE;
    typedef X Y_TYPE;
};*/

template<typename X, typename Y> struct InMemoryData
{
    Vector<pair<X, Y> > data;
    typedef X X_TYPE;
    typedef Y Y_TYPE;
    typedef X const& X_RET;
    InMemoryData(){}
    template<typename DATA> InMemoryData(DATA const& theData):
        data(theData.getSize())
    {
        for(int i = 0; i < data.getSize(); ++i)
        {
            data[i].first = theData.getX(i);
            data[i].second = theData.getY(i);
        }
    }
    void addZ(X const& x, Y const& y){data.append(make_pair(x, y));}
    X_RET getX(int i)const
    {
        assert(i >= 0 && i < data.getSize());
        return data[i].first;
    }
    double getX(int i, int feature)const
    {
        assert(i >= 0 && i < data.getSize() && feature >= 0 &&
            feature < data[i].first.getSize());
        return data[i].first[feature];
    }
    Y_TYPE const& getY(int i)const
    {
        assert(i >= 0 && i < data.getSize());
        return data[i].second;
    }
    int getSize()const{return data.getSize();}
};

template<typename X, typename Y> struct DiskData
{
    EMVector<pair<X, Y> > data;
    DiskData(string const& filename): data(filename) {}
    template<typename DATA> DiskData(DATA const& theData,
        string const& filename): data(filename)
    {
        for(int i = 0; i < theData.getSize(); ++i)
            addZ(theData.getX(i), theData.getY(i));
    }
    void addZ(X const& x, Y const& y){data.append(make_pair(x, y));}
    typedef X X_TYPE;
    typedef Y Y_TYPE;
    typedef X const& X_RET;
    X_RET getX(int i)const
    {
        assert(i >= 0 && i < data.getSize());
        return data[i].first;
    }
    double getX(int i, int feature)const
    {
        assert(i >= 0 && i < data.getSize() &&
            feature >= 0 && feature < data[i].first.getSize());
        return data[i].first[feature];
    }
    Y_TYPE const& getY(int i)const
    {
        assert(i >= 0 && i < data.getSize());
        return data[i].second;
    }
    int getSize()const{return data.getSize();}
};

template<typename LEARNER, typename BUFFER, typename PARAMS =
    EMPTY> class BufferLearner
{
    LEARNER model;
public:
    template<typename DATA> BufferLearner(DATA const& data,
        PARAMS const& p = PARAMS()): model(BUFFER(data), p) {}
    typename BUFFER::Y_TYPE predict(typename BUFFER::X_TYPE const& x)const
        {return model.predict(x);}
};

template<typename DATA> struct RelabeledData
{
    DATA const& data;
    typedef typename DATA::X_TYPE X_TYPE;
    typedef typename DATA::Y_TYPE Y_TYPE;
    typedef typename DATA::X_RET X_RET;
    Vector<Y_TYPE> labels;
    RelabeledData(DATA const& theData): data(theData) {}
    int getSize()const{return data.getSize();}
    void addLabel(Y_TYPE y){labels.append(y);}
    void checkI(int i)const
    {
        assert(i >= 0 && i < data.getSize() &&
            labels.getSize() == data.getSize());
    }
    X_RET getX(int i)const
    {
        checkI(i);
        return data.getX(i);
    }
    double getX(int i, int feature)const
    {
        checkI(i);
        return data.getX(i, feature);
    }
    Y_TYPE const& getY(int i)const
    {
        checkI(i);
        return labels[i];
    }
};

template<typename DATA> struct JoinedData//same type only
{
    DATA const &data1, &data2;
    typedef typename DATA::X_TYPE X_TYPE;
    typedef typename DATA::Y_TYPE Y_TYPE;
    typedef typename DATA::X_RET X_RET;
    JoinedData(DATA const& theData1, DATA const& theData2): data1(theData1),
        data2(theData2){}
    int getSize()const{return data1.getSize() + data2.getSize();}
    X_RET getX(int i)const
    {
        assert(i >= 0 && i < getSize());
        return i < data1.getSize() ? data1.getX(i) :
            data2.getX(i - data1.getSize());
    }
    double getX(int i, int feature)const
    {
        assert(i >= 0 && i < getSize());
        return i < data1.getSize() ? data1.getX(i, feature) :
            data2.getX(i - data1.getSize(), feature);
    }
    Y_TYPE const& getY(int i)const
    {
        assert(i >= 0 && i < getSize());
        return i < data1.getSize() ? data1.getY(i) :
            data2.getY(i - data1.getSize());
    }
};

template<typename DATA> struct PermutedData
{
    DATA const& data;
    Vector<int> permutation;
    typedef typename DATA::X_TYPE X_TYPE;
    typedef typename DATA::Y_TYPE Y_TYPE;
    typedef typename DATA::X_RET X_RET;
    PermutedData(DATA const& theData): data(theData) {}
    int getSize()const{return permutation.getSize();}
    void addIndex(int i){permutation.append(i);}
    void checkI(int i)const
    {
        assert(i >= 0 && i < permutation.getSize() &&
            permutation[i] >= 0 && permutation[i] < data.getSize());
    }
    X_RET getX(int i)const
    {
        checkI(i);
        return data.getX(permutation[i]);
    }
    double getX(int i, int feature)const
    {
        checkI(i);
        return data.getX(permutation[i], feature);
    }
    Y_TYPE const& getY(int i)const
    {
        checkI(i);
        return data.getY(permutation[i]);
    }
};

template<typename DATA> pair<PermutedData<DATA>, PermutedData<DATA> >
    createTrainingTestSetsDetPerm(DATA const& data,
    double relativeTestSize = 0.8)
{
    int n = data.getSize(), m = n * relativeTestSize;
    assert(m > 0 && m < n);
    pair<PermutedData<DATA>, PermutedData<DATA> > result(data, data);
    Vector<int> perm(n);
    for(int i = 0; i < n; ++i) perm[i] = i;
    permuteDeterministically(perm.getArray(), n);
    for(int i = 0; i < n; ++i)
    {
        if(i < m) result.first.addIndex(perm[i]);
        else result.second.addIndex(perm[i]);
    }
    return result;
}

class FeatureSelector
{
public:
    Vector<int> fMap;
public:
    FeatureSelector(Bitset<> const& selection)
    {
        for(int i = 0; i < selection.getSize(); ++i)
            if(selection[i]) fMap.append(i);
    }
    NUMERIC_X select(NUMERIC_X const& x)const
    {
        NUMERIC_X result;
        for(int i = 0; i < fMap.getSize(); ++i) result.append(x[fMap[i]]);
        return result;
    }
    double select(NUMERIC_X const& x, int feature)const
    {
        assert(feature >= 0 && feature < fMap.getSize());
        return x[fMap[feature]];
    }
};

template<typename DATA> struct FSData
{
    DATA const& data;
    FeatureSelector const & f;
    typedef typename DATA::X_TYPE X_TYPE;
    typedef typename DATA::Y_TYPE Y_TYPE;
    typedef X_TYPE X_RET;
    FSData(DATA const& theData, FeatureSelector const & theF): data(theData),
        f(theF) {}
    int getSize()const{return data.getSize();}
    X_RET getX(int i)const{return f.select(data.getX(i));}
    double getX(int i, int feature)const
        {return f.select(data.getX(i), feature);}
    Y_TYPE const& getY(int i)const{return data.getY(i);}
};

template<typename DATA, typename SCALER> struct ScaledData
{
    DATA const& data;
    SCALER const& s;
    typedef typename DATA::X_TYPE X_TYPE;
    typedef typename DATA::Y_TYPE Y_TYPE;
    typedef X_TYPE X_RET;
    ScaledData(DATA const& theData, SCALER const& theS):
        data(theData), s(theS){}
    int getSize()const{return data.getSize();}
    X_RET getX(int i)const
    {
        assert(i >= 0 && i < data.getSize());
        return s.scale(data.getX(i));
    }
    double getX(int i, int feature)const
    {
        assert(i >= 0 && i < data.getSize());
        return s.scaleI(data.getX(i, feature), feature);
    }
    Y_TYPE const& getY(int i)const{return data.getY(i);}
};

template<typename LEARNER, typename Y, typename X = NUMERIC_X>
struct NoParamsLearner
{
    LEARNER model;
    template<typename DATA> NoParamsLearner(DATA const& data,
        EMPTY const& p): model(data) {}
    void learn(X const& x, Y const& y){model.learn(x, y);}
    Y predict(X const& x)const{return model.predict(x);}
    double evaluate(X const& x)const{return model.evaluate(x);}
};

template<typename Y, typename DATA, typename LEARNER> Vector<pair<Y, Y> >
    evaluateLearner(LEARNER const& l, DATA const& test)
{
    Vector<pair<Y, Y> > result;
    for(int i = 0; i < test.getSize(); ++i)
        result.append(pair<Y, Y>(test.getY(i), l.predict(test.getX(i))));
    return result;
}

template<typename LEARNER, typename Y, typename DATA, typename PARAM>
    Vector<pair<Y, Y> > crossValidateGeneral(PARAM const& p,
    DATA const& data, int nFolds = 5)
{
    assert(nFolds > 1 && nFolds <= data.getSize());
    Vector<pair<Y, Y> > result;
    int testSize = data.getSize()/nFolds;//roundoff goes to training
    for(int i = 0; i < nFolds; ++i)
    {
        PermutedData<DATA> trainData(data), testData(data);
        int testStart = i * testSize, testStop = (i + 1) * testSize;
        for(int j = 0; j < data.getSize(); ++j)
        {
            if(testStart <= j && j < testStop) testData.addIndex(j);
            else trainData.addIndex(j);
        }
        LEARNER l(trainData, p);
        result.appendVector(evaluateLearner<Y>(l, testData));
    }
    return result;
}

template<typename PARAM, typename Y, typename DATA, typename LEARNER>
    Vector<pair<Y, Y> > repeatedCVGeneral(PARAM const& p, DATA const& data,
    int nFolds = 5, int nRepeats = 5)
{
    Vector<pair<Y, Y> > result;
    PermutedData<DATA> pData(data);
    for(int i = 0; i < data.getSize(); ++i) pData.addIndex(i);
    for(int i = 0; i < nRepeats; ++i)
    {
        GlobalRNG().randomPermutation(pData.permutation.getArray(),
            data.getSize());
        result.appendVector(crossValidateGeneral<PARAM, Y>(p, data, nFolds));
    }
    return result;
}

struct ScalerMinMax
{
    NUMERIC_X minX, maxX;
public:
    ScalerMinMax(int D): minX(D, numeric_limits<double>::infinity()),
        maxX(D, -numeric_limits<double>::infinity()){}
    template<typename DATA> ScalerMinMax(DATA const& data)
    {
        assert(data.getSize() > 0);
        minX = maxX = data.getX(0);
        for(int i = 1; i < data.getSize(); ++i) addSample(data.getX(i));
    }
    void addSample(NUMERIC_X const& x)
    {
        assert(minX.getSize() == x.getSize());
        for(int j = 0; j < x.getSize(); ++j)
        {
            minX[j] = min(minX[j], x[j]);
            maxX[j] = max(maxX[j], x[j]);
        }
    }
    double scaleI(double xi, int i)const
    {
        double delta = maxX[i] - minX[i];
        return delta > 0 ? (xi - minX[i])/delta : 0;
    }
    NUMERIC_X scale(NUMERIC_X x)const
    {
        for(int i = 0; i < x.getSize(); ++i) x[i] = scaleI(x[i], i);
        return x;
    }
};

template<typename LEARNER, typename Y, typename PARAMS = EMPTY,
    typename SCALER = ScalerMinMax> class ScaledLearner
{
    SCALER s;
    LEARNER l;
public:
    template<typename DATA> ScaledLearner(DATA const& data, PARAMS const& p =
        PARAMS()): s(data), l(ScaledData<DATA, SCALER>(data, s), p){}
    Y predict(NUMERIC_X const& x)const{return l.predict(s.scale(x));}
};

class ScalerMQ
{
    Vector<IncrementalStatistics> ic;
public:
    ScalerMQ(int D): ic(D) {}
    template<typename DATA> ScalerMQ(DATA const& data): ic(getD(data))
        {for(int i = 0; i < data.getSize(); ++i) addSample(data.getX(i));}
    void addSample(NUMERIC_X const& x)
        {for(int j = 0; j < x.getSize(); ++j) ic[j].addValue(x[j]);}
    double scaleI(double xi, int i)const
    {
        double q = ic[i].stdev();
        return q > 0 ? (xi - ic[i].getMean())/q : 0;
    }
    NUMERIC_X scale(NUMERIC_X x)const
    {
        for(int i = 0; i < x.getSize(); ++i) x[i] = scaleI(x[i], i);
        return x;
    }
};

typedef Vector<int> CATEGORICAL_X;
class DiscretizerEqualWidth
{
    ScalerMinMax s;
    int nBins;
    int discretize(double x, int i)const
    {
        x = s.scaleI(x, i);
        if(x < 0) return 0;
        if(x >= 1) return nBins - 1;
        return nBins * x;
    }
public:
    template<typename DATA> DiscretizerEqualWidth(DATA const& data,
        int theNBins = -1): s(data), nBins(theNBins)
    {
        assert(data.getSize() > 1);
        if(nBins == -1) nBins = lgCeiling(data.getSize());
    }
    CATEGORICAL_X operator()(NUMERIC_X const& x)const
    {
        CATEGORICAL_X result;
        for(int i = 0; i < x.getSize(); ++i)
            result.append(discretize(x[i], i));
        return result;
    }
};

class NeuralNetwork
{
    long long learnedCount;
    double initialLearningRate;
    double learningRate()
        {return initialLearningRate * RMRate(learnedCount++);}
    struct Neuron
    {
        Vector<int> sources;
        Vector<double> weights;
        double output, error;
    };
    mutable NUMERIC_X inputs;
    bool isContinuousOutput;
    mutable Vector<Vector<Neuron> > layers;
    double actOutput(double x)const
        {return isContinuousOutput ? x : 1/(1 + exp(-x));}
    double actOutputDeriv(double fx)const
        {return isContinuousOutput ? 1 : fx * (1 - fx);}
    double actInner(double x)const{return tanh(x);}
    double actInnerDeriv(double fx)const{return 1 - fx * fx;}
    void propagateInputs(NUMERIC_X const& x)const
    {
        inputs = x;
        for(int i = 0; i < layers.getSize(); ++i)
            for(int j = 0; j < layers[i].getSize(); ++j)
            {
                Neuron& n = layers[i][j];
                double sum = n.error = 0;
                for(int k = 0; k < n.sources.getSize(); ++k)
                    sum += n.weights[k] * getInput(i, n.sources[k]);
                n.output = i == layers.getSize() - 1 ?
                    actOutput(sum) : actInner(sum);
            }
    }
    double getInput(int layer, int source)const
    {
        return source == -1 ? 1 : layer == 0 ?
            inputs[source] : layers[layer - 1][source].output;
    }
public:
    NeuralNetwork(int D, bool theIsContinuousOutput = false,
        double theInitialLearningRate = 1): inputs(D), learnedCount(0),
        initialLearningRate(theInitialLearningRate),
        isContinuousOutput(theIsContinuousOutput) {}
    void addLayer(int nNeurons){layers.append(Vector<Neuron>(nNeurons));}
    void addConnection(int layer, int from, int to, double weight)
    {
        Vector<Neuron>& last = layers[layer];
        last[from].sources.append(to);
        last[from].weights.append(weight);
    }
    void learn(NUMERIC_X const& x, Vector<double> const& results)
    {
        assert(results.getSize() == layers.lastItem().getSize());
        propagateInputs(x);
        Vector<Neuron>& last = layers.lastItem();
        for(int j = 0; j < last.getSize(); ++j) last[j].error =
            last[j].output - results[j];
        double r = learningRate();
        for(int i = layers.getSize() - 1; i >= 0; --i)
            for(int j = 0; j < layers[i].getSize(); ++j)
            {
                Neuron& n = layers[i][j];
                double temp = n.error * (i == layers.getSize() - 1 ?
                    actOutputDeriv(n.output) : actInnerDeriv(n.output));
                for(int k = 0; k < n.sources.getSize(); ++k)
                {//update weights and prev layer output errors
                    int source = n.sources[k];
                    if(i > 0 && source != -1)
                        layers[i - 1][source].error += n.weights[k] * temp;
                    n.weights[k] -= r * temp * getInput(i, source);
                }
            }
    }
    Vector<double> evaluate(NUMERIC_X const& x)const
    {
        propagateInputs(x);
        Vector<double> result;
        for(int i = 0; i < layers.lastItem().getSize(); ++i)
            result.append(layers.lastItem()[i].output);
        return result;
    }
};

/*
class DeepNeuralNetwork
{
    Vector<int> learnedCounts;
    double initialLearningRate, l;
    double learningRate(int layer)
        {return initialLearningRate * RMRate(learnedCounts[layer]++);}
    struct Neuron
    {
        Vector<int> sources;
        Vector<double> weights;
        double output, error;
    };
    mutable NUMERIC_X inputs;
    bool isContinuousOutput;
    mutable Vector<Vector<Neuron> > layers, layersUnsup;
    double actOutput(double x)const
        {return isContinuousOutput ? x : 1/(1 + exp(-x));}
    double actOutputDeriv(double fx)const
        {return isContinuousOutput ? 1 : fx * (1 - fx);}
    double actInner(double x)const{return tanh(x);}
    double actInnerDeriv(double fx)const{return 1 - fx * fx;}
    double getInput(int layer, int source)const
    {
        return source == -1 ? 1 : layer == 0 ?
            inputs[source] : layers[layer - 1][source].output;
    }
    double collectSum(Neuron& n, int layer)const
    {
        double sum = n.error = 0;
        for(int k = 0; k < n.sources.getSize(); ++k)
            sum += n.weights[k] * getInput(layer, n.sources[k]);
        return sum;
    }
    void propagateUnsup(int i)
    {
        for(int j = 0; j < layers[i].getSize(); ++j)
        {
            Neuron& n = layers[i][j];
            n.output = actInner(collectSum(n, i));
        }
        for(int j = 0; j < layersUnsup[i].getSize(); ++j)
        {
            Neuron& n = layersUnsup[i][j];
            n.output = actOutput(collectSum(n, i + 1));
        }
    }
    void propagateInputs(NUMERIC_X const& x)const
    {
        inputs = x;
        for(int i = 0; i < layers.getSize(); ++i)
            for(int j = 0; j < layers[i].getSize(); ++j)
            {
                Neuron& n = layers[i][j];
                double sum = collectSum(n, i);
                n.output = i == layers.getSize() - 1 ?
                    actOutput(sum) : actInner(sum);
            }
    }
    void addLayer(int nNeurons)
        {layers.append(Vector<Neuron>(nNeurons));}
    void addConnection(Vector<Neuron>& last, int from, int to, double weight)
    {
        last[from].sources.append(to);
        last[from].weights.append(weight);
    }
    void updateNeuron(Neuron& n, int i, double temp, double r, double decay)
    {
        //DEBUG("update");
        for(int k = 0; k < n.sources.getSize(); ++k)
        {//update weights and prev layer output errors
            int source = n.sources[k];
            if(source != -1)
            {
                if(i > 0) layers[i - 1][source].error += n.weights[k] * temp;
                n.weights[k] *= 1 - r * decay;
            }
            n.weights[k] -= r * temp * getInput(i, source);
            //DEBUG(r * temp);
            //DEBUG(getInput(i, source));
            //DEBUG(n.weights[k]);
        }
    }
public:
    DeepNeuralNetwork(int D, int nOutputs, int nHiddenNeurons, double theL = 0, int nHiddenLayers = 3,
        bool theIsContinuousOutput = false,
        double theInitialLearningRate = 1): inputs(D), learnedCounts(nHiddenLayers + 1),
        initialLearningRate(theInitialLearningRate), l(theL),
        isContinuousOutput(theIsContinuousOutput)
    {
        //nHiddenNeurons = pow(D, 0.75);
        assert(D > 0 && nHiddenLayers > 0 && nOutputs > 0);
        Vector<Neuron> decoder(D);
        for(int i = 0; i < nHiddenLayers; ++i)
        {
            addLayer(nHiddenNeurons);
            for(int j = 0; j < nHiddenNeurons; ++j)
                for(int k = -1; k < (i ? nHiddenNeurons : D); ++k)
                    addConnection(layers[i], j, k, GlobalRNG().normal01());
            layersUnsup.append(decoder);
            for(int j = 0; j < D; ++j)
                for(int k = -1; k < nHiddenNeurons; ++k)
                    addConnection(layersUnsup[i], j, k, 0);
        }
        addLayer(nOutputs);
        for(int j = 0; j < nOutputs; ++j)
            for(int k = -1; k < nHiddenNeurons; ++k)
                addConnection(layers[nHiddenLayers], j, k, 0);
    }
    int nHiddenLayers(){return layersUnsup.getSize();}
    void learnUnsupI(NUMERIC_X const& x, int i, int m)
    {
        //DEBUG(i);
        inputs = x;
        for(int j = 0; j <= i; ++j) propagateUnsup(j);
        double r = learningRate(i);
        Vector<Neuron>& last = layersUnsup[i];
        for(int j = 0; j < last.getSize(); ++j) last[j].error =
            last[j].output - x[j];
        //DEBUG("encoder");
        for(int j = 0; j < last.getSize(); ++j)
        {
            Neuron& n = last[j];
            double temp = n.error * actOutputDeriv(n.output);
            updateNeuron(n, i + 1, temp, r, l/m);
        }
        for(int j = 0; j < layers[i].getSize(); ++j)
        {
            Neuron& n = layers[i][j];
            double temp = n.error * actInnerDeriv(n.output);
            updateNeuron(n, i, temp, r, l/m);
        }
    }
    void learnUnsup(NUMERIC_X const& x, int m)
    {
        inputs = x;
        for(int i = 0; i < layersUnsup.getSize(); ++i)
        {
            double r = learningRate(i);
            propagateUnsup(i);
            //DEBUG(i);
            Vector<Neuron>& last = layersUnsup[i];
            for(int j = 0; j < last.getSize(); ++j){ last[j].error =
                last[j].output - x[j]; //DEBUG(last[j].error);
                }
            //system("PAUSE");
            for(int j = 0; j < last.getSize(); ++j)
            {
                Neuron& n = last[j];
                updateNeuron(n, i + 1, n.error, r, l/m);
            }
            for(int j = 0; j < layers[i].getSize(); ++j)
            {
                Neuron& n = layers[i][j];
                double temp = n.error * actInnerDeriv(n.output);
                updateNeuron(n, i, temp, r, l/m);
            }
        }
    }
    void learn(NUMERIC_X const& x, Vector<double> const& results, int m)
    {
        assert(results.getSize() == layers.lastItem().getSize());
        propagateInputs(x);
        Vector<Neuron>& last = layers.lastItem();
        for(int j = 0; j < last.getSize(); ++j) last[j].error =
            last[j].output - results[j];
        for(int lastI = layers.getSize() - 1, i = lastI; i >= 0; --i)
        {
            double r = learningRate(i);
            for(int j = 0; j < layers[i].getSize(); ++j)
            {
                Neuron& n = layers[i][j];
                double temp = n.error * (i == lastI ?
                    actOutputDeriv(n.output) : actInnerDeriv(n.output));
                updateNeuron(n, i, temp, r, i < lastI - 1 ? 0 : l/m);
            }
        }
    }
    Vector<double> evaluate(NUMERIC_X const& x)const
    {
        propagateInputs(x);
        Vector<double> result;
        for(int i = 0; i < layers.lastItem().getSize(); ++i)
            result.append(layers.lastItem()[i].output);
        return result;
    }
};*/

struct BinaryLoss
{
    double operator()(int predicted, int actual)const
        {return predicted != actual;}
};
template<typename LEARNER, typename PARAMS = EMPTY, typename Y = int,
    typename LOSS = BinaryLoss, typename X = NUMERIC_X> class RaceLearner
{
    Vector<LEARNER> learners;
    Vector<double> losses;
    LOSS l;
    int n;
public:
    RaceLearner(Vector<PARAMS> const& p): losses(p.getSize(), 0), n(0)
        {for(int i = 0; i < p.getSize(); ++i)learners.append(LEARNER(p[i]));}
    void learn(X const& x, Y y)
    {
        for(int i = 0; i < learners.getSize(); ++i)
        {
            if(n > 30) losses[i] += l(learners[i].predict(x), y);
            learners[i].learn(x, y);
        }
        ++n;
    }
    Y predict(NUMERIC_X const& x)const
    {
        return learners[argMin(losses.getArray(), losses.getSize())].
            predict(x);
    }
};

template<typename LEARNER> struct FeatureSubsetLearner
{
    FeatureSelector f;
    LEARNER l;
public:
    template<typename DATA> FeatureSubsetLearner(DATA const& data,
        Bitset<>const& selection): f(selection), l(FSData<DATA>(data, f)) {}
    int predict(NUMERIC_X const& x)const{return l.predict(f.select(x));}
};

/*
template<typename VECTOR> struct MaxDistance
{
    static double iDistanceIncremental(
        VECTOR const& lhs, VECTOR const& rhs, int i)
        {return abs(lhs[i] - rhs[i]);}
    static double distanceIncremental(VECTOR const& lhs,
        VECTOR const& rhs, double bound = numeric_limits<double>::max())
    {
        assert(lhs.getSize() == rhs.getSize());
        double sum = 0;
        for(int i = 0; i < lhs.getSize() && sum < bound; ++i)
            sum = max(sum, iDistanceIncremental(lhs, rhs, i));
        return sum;
    }
    struct Distance
    {
        double operator()(VECTOR const& lhs, VECTOR const& rhs)const
            {return distanceIncremental(lhs, rhs);}
    };
    struct DistanceIncremental
    {
        double operator()(VECTOR const& lhs, VECTOR const& rhs)const
            {return distanceIncremental(lhs, rhs);}
        double operator()(VECTOR const& lhs, VECTOR const& rhs, int i)const
            {return iDistanceIncremental(lhs, rhs, i);}
        double operator()(double bound, VECTOR const& lhs, VECTOR const& rhs)
            const{return distanceIncremental(lhs, rhs, bound);}
    };
};
double digamma(double x)
    {return log(x) - 1/(2 * x) - 1/(12 * x * x);}

double estimateMI(Vector<pair<double, double> > const& data, int k = 6)
{
    assert(data.getSize() > k);
    double minX = data[0].first, maxX = data[0].first, minY = data[0].second, maxY = data[0].second;
    for(int i = 1; i < data.getSize(); ++i)
    {
        minX = min(minX, data[i].first);
        maxX = max(maxX, data[i].first);
        minY = min(minY, data[i].second);
        maxY = max(maxY, data[i].second);
    }
    double deltaX = maxX - minX, deltaY = maxY - minY;
    if(deltaX == 0) deltaX = 1;
    if(deltaY == 0) deltaY = 1;
    KDTree<Point2, bool, 2> index;

    for(int i = 0; i < data.getSize(); ++i)
    {
        index.insert(Point2((data[i].first - minX)/deltaX, (data[i].second - minY)/deltaY), true);
    }
    MaxDistance<Point2>::DistanceIncremental d;
    double sum = 0;
    for(int i = 0; i < data.getSize(); ++i)
    {
        Point2 p((data[i].first - minX)/deltaX, (data[i].second - minY)/deltaY);

        Vector<KDTree<Point2, bool, 2>::NodeType*> neighbors = index.kNN(p, k + 1, d);
        double dist = d(p, neighbors.lastItem()->key);
        int nx = 1, ny = 1;
        for(int j = 0; j < data.getSize(); ++j) if(i != j)
        {
            if(abs(data[j].first - data[i].first) < dist) ++nx;
            if(abs(data[j].second - data[i].second) < dist) ++ny;
        }
        sum += digamma(nx) + digamma(ny);
    }
    return digamma(k) + digamma(data.getSize()) - sum/data.getSize();
}*/

template<typename RISK_FUNCTOR> Vector<Bitset<> > selectFeatures1F(
    RISK_FUNCTOR const &r, int D)
{
    Vector<Bitset<> > selections;
    Vector<double> risks(D);
    for(int i = 0; i < D; ++i)
    {
        Bitset<> temp(D);
        temp.set(i);
        risks[i] = r(temp);
    }
    Vector<int> indices(D);
    for(int i = 0; i < D; ++i) indices[i] = i;
    IndexComparator<double> c(risks.getArray());
    quickSort(indices.getArray(), 0, D - 1, c);
    Bitset<> resultI(D);
    for(int i = 0; i < D; ++i)
    {
        resultI.set(indices[i]);
        selections.append(resultI);
    }
    return selections;
}

template<typename RISK_FUNCTOR> Bitset<> pickBestSubset(
    RISK_FUNCTOR const &r, Vector<Bitset<> >const& subsets)
{
    Bitset<> best = valMinFunc(subsets.getArray(), subsets.getSize(), r);
    DEBUG(best.popCount());
    best.debug();
    return best;
}

template<typename RISK_FUNCTOR> Bitset<> pickBestSubsetGreedy(
    RISK_FUNCTOR const &r, Vector<Bitset<> >const& subsets)
{
    int best = subsets.getSize() - 1;
    double fullRisk = r(subsets[best]);
    for(int i = 0; i < best; ++i) if(r(subsets[i]) <= fullRisk) best = i;
    DEBUG(subsets[best].popCount());
    subsets[best].debug();
    return subsets[best];
}

Vector<Bitset<> > subSampleSubsets(Vector<Bitset<> >const& subsets, int limit)
{
    assert(subsets.getSize() > 0 && limit > 0);
    Vector<Bitset<> > result;
    int skip = ceiling(subsets.getSize(), limit);
    for(int i = subsets.getSize() - 1; i >= 0; i -= skip)
        result.append(subsets[i]);
    result.reverse();
    return result;
}

template<typename RISK_FUNCTOR>
Bitset<> selectFeaturesForwardGreedy(RISK_FUNCTOR const &r, int D)
{
    Bitset<> resultI(D);
    resultI.setAll();
    double fullRisk = r(resultI);
    resultI.setAll(0);
    for(int i = 0; i < D; ++i)
    {
        double bestRisk;
        int bestJ = -1;
        for(int j = 0; j < D; ++j) if(!resultI[j])
        {
            if(i == D - 1)
            {
                bestJ = j;
                break;
            }
            resultI.set(j, true);
            double risk = r(resultI);
            resultI.set(j, false);
            if(bestJ == -1 || risk < bestRisk)
            {
                bestRisk = risk;
                bestJ = j;
            }
        }
        resultI.set(bestJ, true);
        if(r(resultI) <= fullRisk)
        {
            DEBUG(resultI.popCount());
            resultI.debug();
            return resultI;
        }
    }
    resultI.setAll();
    DEBUG(resultI.popCount());
    resultI.debug();
    return resultI;
}

struct SubsetLengthComparator
{
    bool operator()(Bitset<> const& lhs, Bitset<> const& rhs)const
        {return lhs.popCount() < rhs.popCount();}
    bool isEqual(Bitset<> const& lhs, Bitset<> const& rhs)const
        {return lhs.popCount() == rhs.popCount();}
};
Vector<Bitset<> > selectFeaturesAllSubsets(int D)
{
    assert(D <= 20);//computational safety
    int n = pow(2, D) - 1;
    Vector<Bitset<> > result(n, Bitset<>(D));
    for(int i = 0; i < n; ++i)
    {
        int rank = i + 1;
        for(int j = 0; rank > 0; ++j, rank /= 2)
            if(rank % 2) result[i].set(j);
    }
    quickSort(result.getArray(), 0, result.getSize() - 1,
        SubsetLengthComparator());
    return result;
}

template<typename RISK_FUNCTOR> Bitset<> selectFeaturesSmart(
    RISK_FUNCTOR const& r, int D, int subsampleLimit = 20)
{
    if(D <= 12) return pickBestSubset(r, selectFeaturesAllSubsets(D));
    else if(D <= 40) return selectFeaturesForwardGreedy(r, D);
    else return pickBestSubsetGreedy(r, subSampleSubsets(
        selectFeatures1F(r, D), subsampleLimit));
}

template<typename X, typename Y> struct LearnerInterface
{
    virtual Y predict(X const& x)const = 0;
    virtual LearnerInterface* clone()const = 0;
};
template<typename LEARNER, typename X, typename Y>
struct TypeFreeLearner: public LearnerInterface<X, Y>
{
    LEARNER model;
    template<typename DATA> TypeFreeLearner(DATA const& data): model(data) {}
    Y predict(X const& x)const{return model.predict(x);}
    LearnerInterface<X, Y>* clone()const{return new TypeFreeLearner(*this);}
};
template<typename Y, typename X = NUMERIC_X>
class BestCombiner
{
    LearnerInterface<X, Y>* model;
    double risk;
public:
    BestCombiner(): model(0) {}
    BestCombiner(BestCombiner const& rhs): model(rhs.model->clone()) {}
    BestCombiner& operator=(BestCombiner const& rhs)
        {return genericAssign(*this, rhs);}
    template<typename LEARNER, typename DATA, typename RISK_FUNCTOR>
        void addNoParamsClassifier(DATA const& data, RISK_FUNCTOR const& r)
    {
        double riskNew = r(EMPTY());
        if(!model || riskNew < risk)
        {
            delete model;
            model = new TypeFreeLearner<LEARNER, X, Y>(data);
            risk = riskNew;
        }
    }
    Y predict(X const& x)const{assert(model); return model->predict(x);}
    ~BestCombiner(){delete model;}
};

}//end namespace
#endif

