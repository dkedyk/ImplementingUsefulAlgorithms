#ifndef IGMDK_COST_CLASSIFICATION_H
#define IGMDK_COST_CLASSIFICATION_H
#include "ClassificationCommon.h"
#include "RandomForest.h"
#include "KernelSVM.h"
#include "ImbalanceClassification.h"
#include "../Utils/Utils.h"
#include "../HashTable/ChainingHashTable.h"
#include "../HashTable/LinearProbingHashTable.h"
#include "../NumericalMethods/Matrix.h"
#include "../NumericalMethods/NumericalMethods.h"
#include "../NumericalOptimization/NumericalOptimization.h"
#include "../RandomNumberGeneration/Statistics.h"
#include <cmath>
namespace igmdk{

int costClassify(Vector<double> const& probs, Matrix<double> const& cost)
{
    int k = probs.getSize();
    assert(k == cost.getRows());
    Vector<double> losses(k);
    for(int i = 0; i < k; ++i)
        for(int j = 0; j < k; ++j) losses[i] += probs[j] * cost(j, i);
    return argMin(losses.getArray(), k);
}

template<typename LEARNER = RandomForest> class CostLearner
{
    Matrix<double> cost;
    LEARNER model;
public:
    template<typename DATA> CostLearner(DATA const& data,
        Matrix<double>const& costMatrix): model(data), cost(costMatrix) {}
    int predict(NUMERIC_X const& x)const
        {return costClassify(model.classifyProbs(x), cost);}
};

void scaleCostMatrix(Matrix<double>& cost)
{
    double maxCost = 0;
    for(int r = 0; r < cost.getRows(); ++r)
        for(int c = 0; c < cost.getRows(); ++c)
            maxCost = max(maxCost, cost(r, c));
    cost *= 1/maxCost;
}

template<typename LEARNER = NoParamsLearner<DecisionTree, int>,
    typename PARAMS = EMPTY, typename X = NUMERIC_X> class RMBoost
{
    Vector<LEARNER> classifiers;
    int nClasses;
    struct BinomialLoss
    {
        Vector<Vector<double> > F;
        BinomialLoss(int n, int nClasses): F(n, Vector<double>(nClasses, 0))
            {}
        int findBestFalse(int i, int label)
        {
            double temp = F[i][label];
            F[i][label] = -numeric_limits<double>::infinity();
            double result = argMax(F[i].getArray(), F[i].getSize());
            F[i][label] = temp;
            return result;
        }
        double getNegGrad(int i, int label, Matrix<double>const& costMatrix)
        {
            int bestFalseLabel = findBestFalse(i, label);
            double margin = F[i][label] - F[i][bestFalseLabel];
            return costMatrix(label, bestFalseLabel)/(exp(margin) + 1);
        }
    };
public:
    template<typename DATA> RMBoost(DATA const& data, Matrix<double>
        costMatrix, PARAMS const& p = PARAMS(),
        int nClassifiers = 100): nClasses(findNClasses(data))
    {//initial weights are based on ave cost
        int n = data.getSize();
        assert(n > 0 && nClassifiers > 0);
        BinomialLoss l(n, nClasses);
        Vector<double> dataWeights(n), classWeights(nClasses);
        for(int i = 0; i < nClasses; ++i)
            for(int j = 0; j < nClasses; ++j)
                classWeights[i] += costMatrix(i, j);
        for(int i = 0; i < n; ++i)
            dataWeights[i] = classWeights[data.getY(i)];
        for(int i = 0; i < nClassifiers; ++i)
        {
            normalizeProbs(dataWeights);
            AliasMethod sampler(dataWeights);
            PermutedData<DATA> resample(data);
            for(int j = 0; j < n; ++j) resample.addIndex(sampler.next());
            classifiers.append(LEARNER(resample, p));
            for(int j = 0; j < n; ++j)
            {
                l.F[j][classifiers.lastItem().predict(data.getX(j))] +=
                    RMRate(i);
                dataWeights[j] = l.getNegGrad(j, data.getY(j), costMatrix);
            }
        }
    }
    int predict(X const& x)const
    {
        Vector<double> counts(nClasses, 0);
        for(int i = 0; i < classifiers.getSize(); ++i)
            counts[classifiers[i].predict(x)] += RMRate(i);
        return argMax(counts.getArray(), counts.getSize());
    }
};

class BoostedCostSVM
{
    RMBoost<MulticlassSVM<>, pair<GaussianKernel, double> > model;
public:
    template<typename DATA> BoostedCostSVM(DATA const& data,
        Matrix<double> const& cost = Matrix<double>(1, 1)):
        model(data, cost, NoParamsSVM::gaussianMultiClassSVM(data), 15) {}
    int predict(NUMERIC_X const& x)const{return model.predict(x);}
};
typedef ScaledLearner<BoostedCostSVM, int, Matrix<double> > SBoostedCostSVM;

template<typename LEARNER, typename PARAMS = EMPTY,
    typename X = NUMERIC_X> class AveCostLearner
{
    LEARNER model;
    template<typename DATA> static Vector<double> findWeights(
        DATA const& data, Matrix<double> const& costMatrix)
    {//init with average weights
        int k = costMatrix.getRows(), n = data.getSize();
        assert(k > 1 && k == findNClasses(data));
        Vector<double> classWeights(k), result(n);
        for(int i = 0; i < k; ++i)
            for(int j = 0; j < k; ++j)
                classWeights[i] += costMatrix(i, j);
        for(int i = 0; i < n; ++i) result[i] = classWeights[data.getY(i)];
        normalizeProbs(result);
        return result;
    }
public:
    template<typename DATA> AveCostLearner(DATA const& data,
        Matrix<double> const& costMatrix, PARAMS const& p = PARAMS()):
        model(data, findWeights(data, costMatrix), p) {}
    int predict(X const& x)const{return model.predict(x);}
};

class AveCostSVM
{
    typedef pair<GaussianKernel, double> P;
    AveCostLearner<WeightedBaggedLearner<MulticlassSVM<>, P>, P> model;
public:
    template<typename DATA> AveCostSVM(DATA const& data,
        Matrix<double>const& cost = Matrix<double>(1, 1)):
        model(data, cost, NoParamsSVM::gaussianMultiClassSVM(data)) {}
    int predict(NUMERIC_X const& x)const{return model.predict(x);}
};
typedef ScaledLearner<AveCostSVM, int, Matrix<double> > SAveCostSVM;

}//end namespace
#endif

