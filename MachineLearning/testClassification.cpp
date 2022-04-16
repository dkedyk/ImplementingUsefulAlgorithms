#include "CostClassification.h"
#include "ImbalanceClassification.h"
#include "Classification.h"
#include "NaiveBayes.h"
#include "KNN.h"
#include "DecisionTree.h"
#include "RandomForest.h"
#include "LSVM.h"
#include "KernelSVM.h"
#include "NeuralNetwork.h"
#include "ReadClassificationData.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "../RandomNumberGeneration/Statistics.h"
#include "../ExternalMemoryAlgorithms/File.h"
using namespace igmdk;

/*
template<typename DISTANCE = EuclideanDistance<NUMERIC_X>::Distance>
struct MeanNN
{
    Vector<NUMERIC_X> means;
    DISTANCE d;
public:
    template <typename DATA> MeanNN(DATA const& data,
        DISTANCE const& theD = DISTANCE()): d(theD)
    {
        int D = getD(data), nClasses = findNClasses(data);
        NUMERIC_X zero(D);
        Vector<int> counts(nClasses);
        means = Vector<NUMERIC_X>(nClasses, zero);
        for(int i = 0; i < data.getSize(); ++i)
        {
            int y = data.getY(i);
            means[y] += data.getX(i);
            ++counts[y];
        }
        for(int i = 0; i < nClasses; ++i)
            if(counts[i] > 0) means[i] *= 1.0/counts[i];
    }
public:
    int predict(NUMERIC_X const& x)const
    {
        int best = -1;
        double bestD;
        for(int i = 0; i < means.getSize(); ++i)
        {
            double dist = d(means[i], x);
            if(best == -1 || dist < bestD)
            {
                best = i;
                bestD = dist;
            }
        }
        return best;
    }
};*/

Matrix<double> sampleCostDeterministic(int k)
{
    assert(k > 0);
    int count = 0;
    Matrix<double> result(k, k);
    for(int r = 0; r < k; ++r)
        for(int c = 0; c < k; ++c) if(r != c) result(r, c) = count++ % 2 ? 0.01 : 1;
    scaleCostMatrix(result);
    return result;
}

template<typename LEARNER> int testNumericalClassifier()
{
    DEBUG("Started Reading");
    typedef InMemoryData<NUMERIC_X, int> T;
    Vector<T> dataM(50);//make many enough to avoid ref realloc
    Vector<pair<PermutedData<T>, PermutedData<T> > > data;

    dataM.append(T());
    readIrisData(dataM.lastItem());
    data.append(makeData<T>(dataM));

    /*dataM.append(T());
    readDigitData(dataM.lastItem(), true);
    dataM.append(T());
    readDigitData(dataM.lastItem(), false);
    data.append(makeDataDivided<T>(dataM));*/

    /*dataM.append(T());
    readCNEAData(dataM.lastItem());
    data.append(makeData<T>(dataM));
    dataM.append(T());
    readBanknoteData(dataM.lastItem());
    data.append(makeData<T>(dataM));
    dataM.append(T());
    readWDBCData(dataM.lastItem());
    data.append(makeData<T>(dataM));
    dataM.append(T());
    readGlassData(dataM.lastItem());
    data.append(makeData<T>(dataM));
    dataM.append(T());
    readIonosphereData(dataM.lastItem());
    data.append(makeData<T>(dataM));

    dataM.append(T());
    readLetterData(dataM.lastItem());
    data.append(makeData<T>(dataM));
    dataM.append(T());
    readPimaData(dataM.lastItem());
    data.append(makeData<T>(dataM));
    dataM.append(T());
    readSpamData(dataM.lastItem());
    data.append(makeData<T>(dataM));
    dataM.append(T());
    readSpectData(dataM.lastItem(), true);
    dataM.append(T());
    readSpectData(dataM.lastItem(), false);
    data.append(makeDataDivided<T>(dataM));

    dataM.append(T());
    readStatlogData(dataM.lastItem(), true);
    dataM.append(T());
    readStatlogData(dataM.lastItem(), false);
    data.append(makeDataDivided<T>(dataM));

    dataM.append(T());
    readWineData(dataM.lastItem());
    data.append(makeData<T>(dataM));

    dataM.append(T());
    readArceneData(dataM.lastItem(), true);
    dataM.append(T());
    readArceneData(dataM.lastItem(), false);
    data.append(makeDataDivided<T>(dataM));

    dataM.append(T());
    readMadelonData(dataM.lastItem(), true);
    dataM.append(T());
    readMadelonData(dataM.lastItem(), false);
    data.append(makeDataDivided<T>(dataM));*/

    DEBUG("Done Reading");
    int reportNumber = time(0);
    string fAcc = "reportAcc" + to_string(reportNumber) + ".csv";
    string fAveAcc = "reportAveAcc" + to_string(reportNumber) + ".csv";
    string fTimer = "reportTimer" + to_string(reportNumber) + ".csv";
    string fCount = "reportFeature" + to_string(reportNumber) + ".csv";
    string fCost = "reportCost" + to_string(reportNumber) + ".csv";
    ++reportNumber;
    for(int i = 0; i < data.getSize(); ++i)
    {
        if(false)//cost
        {
            /*int start = clock();
            int k = findNClasses(data[i].first);
            Matrix<double> c = sampleCostDeterministic(k);
            for(int z = 0; z < k; ++z)
            {
                for(int j = 0; j < k; ++j)
                    cout << c(z, j) << " ";
                cout << endl;
            }
            //LEARNER s(T(data[i].first).data);
            LEARNER s(T(data[i].first).data, c);
            Matrix<int> confusion = evaluateConfusion(evaluateLearner<int>(s, T(data[i].second).data));
            for(int z = 0; z < k; ++z)
            {
                for(int j = 0; j < k; ++j)
                    cout << confusion(z, j) << " ";
                cout << endl;
            }
            ClassifierStats cs(confusion);
            double timediff = 1.0 * (clock() - start)/CLOCKS_PER_SEC;
            double cost = evalConfusionCost(confusion, c);
            DEBUG(cost);
            cs.debug();
            addToCSV(Vector<string>(1, to_string(cost)), fCost.c_str());
            addToCSV(Vector<string>(1, to_string(timediff)), fTimer.c_str());*/
        }
        else
        {
            int start = clock();
            LEARNER l(data[i].first);
            ClassifierStats cs(evaluateConfusion(evaluateLearner<int>(l, data[i].second)));
            double timediff = 1.0 * (clock() - start)/CLOCKS_PER_SEC;
            cs.debug();

            /*addToCSV(Vector<string>(1, to_string(cs.acc.mean)), fAcc.c_str());
            addToCSV(Vector<string>(1, to_string(cs.bac.mean)), fAveAcc.c_str());
            addToCSV(Vector<string>(1, to_string(timediff)), fTimer.c_str());*/
            //addToCSV(Vector<string>(1, to_string(s.model.f.fMap.getSize())), fCount.c_str());
        }
        //system("PAUSE");
    }
    return 0;
}

Matrix<double> getEqualCostMatrix(int nClasses)
{//for testing RMBoost
    Matrix<double> result(nClasses, nClasses);
    for(int i = 0; i < nClasses; ++i)
        for(int j = 0; j < nClasses; ++j) if(i != j) result(i, j) = 1;
    return result;
}
void testNumericalClassifiers()
{
//    testNumericalClassifier<SSVM>();
//    testNumericalClassifier<DecisionTree>();
//
//
//    testNumericalClassifier<SNN>();
    testNumericalClassifier<SLSVM>();
//    testNumericalClassifier<RMBoost<> >();
//    testNumericalClassifier<RandomForest>();
//
//
//    testNumericalClassifier<SImbSVM>();
//
//    testNumericalClassifier<SRaceLSVM>();
//    testNumericalClassifier<SOnlineNN>();
//
//    testNumericalClassifier<ScaledLearner<NoParamsLearner<KNNClassifier<>, int>, int> >();
//
//
//
    //testNumericalClassifier<SimpleBestCombiner>();

    //inactive learners
    //testNumericalClassifier<ScaledLearner<MeanNN<>, int> >();
    //testNumericalClassifier<NumericalBayes>();

    //feature selection
    //testNumericalClassifier<SmartFSLearner<> >();

    //cost learning
    //testNumericalClassifier<SBoostedCostSVM>();
    //testNumericalClassifier<SAveCostSVM>();
    //testNumericalClassifier<CostLearner<> >();

}

int main(int argc, char *argv[])
{
    for(int i = 0; i < 1; ++i) testNumericalClassifiers();
	return 0;
}


