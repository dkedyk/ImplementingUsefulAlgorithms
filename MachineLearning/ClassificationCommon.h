#ifndef IGMDK_CLASSIFICATION_COMMON_H
#define IGMDK_CLASSIFICATION_COMMON_H
#include "LearningCommon.h"
#include "../Utils/Debug.h"
#include "../NumericalMethods/Matrix.h"
#include "../HashTable/ChainingHashTable.h"
#include "../RandomNumberGeneration/Statistics.h"
#include <cmath>
namespace igmdk{

template<typename DATA> int findNClasses(DATA const& data)
{
    int maxClass = -1;
    for(int i = 0; i < data.getSize(); ++i)
        maxClass = max(maxClass, data.getY(i));
    return maxClass + 1;
}

template<typename DATA> pair<PermutedData<DATA>, PermutedData<DATA> >
    createTrainingTestSetsStatified(DATA const& data,
    double relativeTestSize = 0.8)
{
    int n = data.getSize(), m = n * relativeTestSize;
    assert(m > 0 && m < n);
    pair<PermutedData<DATA>, PermutedData<DATA> > result(data, data);
    Vector<int> counts(findNClasses(data)), p(n);//need p for legacy only
    for(int i = 0; i < n; ++i){++counts[data.getY(i)]; p[i] = i;}
    for(int i = 0; i < counts.getSize(); ++i) counts[i] *= relativeTestSize;
    for(int i = 0; i < p.getSize(); ++i)
    {
        int label = data.getY(p[i]);
        if(counts[label]){--counts[label]; result.first.addIndex(p[i]);}
        else
        {
            result.second.addIndex(p[i]);
            p[i--] = p.lastItem();
            p.removeLast();
        }
    }
    return result;
}

Matrix<int> evaluateConfusion(Vector<pair<int, int> > const& testResult,
    int nClasses = -1)
{
    if(nClasses == -1)
    {//calculate nClasses if unknown
        int maxClass = 0;
        for(int i = 0; i < testResult.getSize(); ++i) maxClass =
            max(maxClass, max(testResult[i].first, testResult[i].second));
        nClasses = maxClass + 1;
    }
    Matrix<int> result(nClasses, nClasses);
    for(int i = 0; i < testResult.getSize(); ++i)
        ++result(testResult[i].first, testResult[i].second);
    return result;
}

double evalConfusionCost(Matrix<int>const& confusion,
    Matrix<double>const& cost)
{
    int k = confusion.rows, total = 0;
    assert(k == confusion.columns && k == cost.rows && k == cost.columns);
    double sum = 0;
    for(int r = 0; r < k; ++r)
        for(int c = 0; c < k; ++c)
        {
            total += confusion(r, c);
            sum += confusion(r, c) * cost(r, c);
        }
    return sum/total;
}

struct ClassifierStats
{
    double acc, bac;
    pair<double, double> accConf, bacConf;
    Vector<double> accByClass, confByClass;
    int total;
    ClassifierStats(Matrix<int> const& confusion): total(0)
    {//same row = same label, same column = same prediction
        Vector<int> confTotal, accTotal;
        int k = confusion.getRows(), nBac = 0, actualK = 0;
        IncrementalStatistics accS, basSW;
        Vector<IncrementalStatistics> precS(k);
        Vector<double> weights(k);
        for(int r = 0; r < k; ++r)
        {
            int totalR = 0;
            for(int c = 0; c < k; ++c)
            {
                totalR += confusion(r, c);
                weights[r] += confusion(r, c);
                total += confusion(r, c);
            }
            accTotal.append(totalR);
            actualK += (totalR > 0);
        }
        double M = 0;
        for(int r = 0; r < k; ++r)
        {
            weights[r] = total/weights[r]/actualK;
            IncrementalStatistics bacS;
            for(int c = 0; c < k; ++c)
            {
                int count = confusion(r, c);
                bool correct = r == c;
                while(count--)
                {
                    accS.addValue(correct);
                    basSW.addValue(correct * weights[r]);
                    M += weights[r] * weights[r];
                    bacS.addValue(correct);
                    precS[c].addValue(correct);
                }
            }
            accByClass.append(bacS.getMean());
        }
        M = sqrt(M/total);
        for(int c = 0; c < k; ++c)
        {
            int totalC = 0;
            for(int r = 0; r < k; ++r) totalC += confusion(r, c);
            confTotal.append(totalC);
            confByClass.append(precS[c].getMean());
        }
        acc = accS.getMean();
        accConf = wilsonScoreInterval(acc, accS.n);
        bac = basSW.getMean();
        bacConf = HoefFunctor::conf(bac, basSW.n);
        bac *= M;
        bacConf.first *= M;
        bacConf.second *= M;
    }
    void debug()const
    {
        DEBUG(acc * total);
        DEBUG(total);
        cout << "Accuracy: ";
        cout << acc << " ";
        cout << "95% interval: ";
        cout << accConf.first << " ";
        cout << accConf.second;
        cout << endl;
        cout << "Balanced Accuracy: ";
        cout << bac << " ";
        cout << "95% interval: ";
        cout << bacConf.first << " ";
        cout << bacConf.second;
        cout << endl;
        cout << "Accuracy by class: " << endl;
        for(int i = 0; i < accByClass.getSize(); ++i)
        {
            cout << accByClass[i] << " ";
            cout << endl;
        }
        cout << "Confidence by class: " << endl;
        for(int i = 0; i < confByClass.getSize(); ++i)
        {
            cout << confByClass[i] << " ";
            cout << endl;
        }
    }
};

template<typename LEARNER, typename DATA, typename PARAMS>
    Vector<pair<int, int> > crossValidationStratified(PARAMS const& p,
    DATA const& data, int nFolds = 5)
{
    assert(nFolds > 1 && nFolds <= data.getSize());
    int nClasses = findNClasses(data), testSize = 0;
    Vector<int> counts(nClasses, 0), starts(nClasses, 0);
    PermutedData<DATA> pData(data);
    for(int i = 0; i < data.getSize(); ++i)
    {
        pData.addIndex(i);
        ++counts[data.getY(i)];
    }
    for(int i = 0; i < counts.getSize(); ++i)
        counts[i] /= nFolds;//roundoff goes to training
    for(int i = 0; i < counts.getSize(); ++i) testSize += counts[i];
    Vector<pair<int, int> > result;
    for(int i = 0;; ++i)
    {//create list of included test examples in increasing order
        Vector<int> includedCounts(nClasses, 0), includedIndices;
        for(int j = valMin(starts.getArray(), starts.getSize());
            includedIndices.getSize() < testSize; ++j)
        {
            int label = data.getY(j);
            if(starts[label] <= j && includedCounts[label] < counts[label])
            {
                ++includedCounts[label];
                includedIndices.append(j);
                starts[label] = j + 1;
            }
        }
        PermutedData<DATA> testData(data);
        for(int j = testSize - 1; j >= 0; --j)
        {
            testData.addIndex(includedIndices[j]);
            pData.permutation[includedIndices[j]] =
                pData.permutation.lastItem();
            pData.permutation.removeLast();
        }
        result.appendVector(evaluateLearner<int>(LEARNER(pData, p),
            testData));
        //put test data back into data in correct places
        if(i == nFolds - 1) break;
        for(int j = 0; j < testSize; ++j)
        {
            pData.addIndex(includedIndices[j]);
            pData.permutation[includedIndices[j]] =
                testData.permutation[testSize - 1 - j];
        }
    }
    return result;
}
template<typename LEARNER, typename DATA, typename PARAMS> double
    crossValidation(PARAMS const& p, DATA const& data, int nFolds = 5)
{
    return ClassifierStats(evaluateConfusion(
        crossValidationStratified<LEARNER>(p, data, nFolds))).acc;
}

template<typename LEARNER, typename PARAM, typename DATA>
struct SCVRiskFunctor
{
    DATA const& data;
    SCVRiskFunctor(DATA const& theData): data(theData) {}
    double operator()(PARAM const& p)const
        {return 1 - crossValidation<LEARNER>(p, data);}
};

template<typename LEARNER, typename PARAMS = EMPTY, typename X = NUMERIC_X>
class OnlineMulticlassLearner
{
    mutable Treap<int, LEARNER> binaryLearners;
    int nClasses;
    PARAMS p;
    int makeKey(short label1, short label2) const
        {return label1 * numeric_limits<short>::max() + label2;}
public:
    OnlineMulticlassLearner(PARAMS const& theP = PARAMS(),
        int initialNClasses = 0): nClasses(initialNClasses), p(theP) {}
    void learn(X const& x, int label)
    {
        nClasses = max(nClasses, label + 1);
        for(int j = 0; j < nClasses; ++j) if(j != label)
        {
            int key = j < label ? makeKey(j, label) : makeKey(label, j);
            LEARNER* s = binaryLearners.find(key);
            if(!s)
            {
                binaryLearners.insert(key, LEARNER(p));
                s = binaryLearners.find(key);
            }
            s->learn(x, int(j < label));
        }
    }
    int predict(X const& x)const
    {
        assert(nClasses > 0);
        Vector<int> votes(nClasses, 0);
        for(int j = 0; j < nClasses; ++j)
            for(int k = j + 1; k < nClasses; ++k)
            {
                LEARNER* s = binaryLearners.find(makeKey(j, k));
                if(s) ++votes[s->predict(x) ? k : j];
            }
        return argMax(votes.getArray(), votes.getSize());
    }
};

template<typename LEARNER, typename PARAMS = EMPTY, typename X = NUMERIC_X>
class MulticlassLearner
{//if params not passed, uses default value!
    mutable ChainingHashTable<int, LEARNER> binaryLearners;
    int nClasses;
public:
    Vector<LEARNER const*> getLearners()const
    {
        Vector<LEARNER const*> result;
        for(typename ChainingHashTable<int, LEARNER>::Iterator i =
            binaryLearners.begin(); i != binaryLearners.end(); ++i)
            result.append(&i->value);
        return result;
    };
    template<typename DATA> MulticlassLearner(DATA const& data,
        PARAMS const&p = PARAMS()): nClasses(findNClasses(data))
    {
        Vector<Vector<int> > labelIndex(nClasses);
        for(int i = 0; i < data.getSize(); ++i)
            labelIndex[data.getY(i)].append(i);
        for(int j = 0; j < nClasses; ++j) if(labelIndex[j].getSize() > 0)
            for(int k = j + 1; k < nClasses; ++k)
                if(labelIndex[k].getSize() > 0)
                {
                    PermutedData<DATA> twoClassData(data);
                    RelabeledData<PermutedData<DATA> >
                        binaryData(twoClassData);
                    for(int l = 0, m = 0; l < labelIndex[j].getSize() ||
                        m < labelIndex[k].getSize(); ++l, ++m)
                    {
                        if(l < labelIndex[j].getSize())
                        {
                            twoClassData.addIndex(labelIndex[j][l]);
                            binaryData.addLabel(0);
                        }
                        if(m < labelIndex[k].getSize())
                        {
                            twoClassData.addIndex(labelIndex[k][m]);
                            binaryData.addLabel(1);
                        }
                    }
                    binaryLearners.insert(j * nClasses + k,
                        LEARNER(binaryData, p));
                }
    }
    int predict(X const& x)const
    {
        Vector<int> votes(nClasses, 0);
        for(int j = 0; j < nClasses; ++j)
            for(int k = j + 1; k < nClasses; ++k)
            {
                LEARNER* s = binaryLearners.find(j * nClasses + k);
                if(s) ++votes[s->predict(x) ? k : j];
            }
        return argMax(votes.getArray(), votes.getSize());
    }
    int classifyByProbs(X const& x)const
    {//for probability-output learners like neural network
        Vector<double> votes(nClasses, 0);
        for(int j = 0; j < nClasses; ++j)
            for(int k = j + 1; k < nClasses; ++k)
            {
                LEARNER* s = binaryLearners.find(j * nClasses + k);
                if(s)
                {
                    double p = s->evaluate(x);
                    votes[k] += p;
                    votes[j] += 1 - p;
                }
            }
        return argMax(votes.getArray(), votes.getSize());
    }
};

}//end namespace
#endif

