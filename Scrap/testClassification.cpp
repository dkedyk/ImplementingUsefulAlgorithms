#include "SAMME.h"
#include "MRMR.h"
#include "LSHLearner.h"
#include "../MachineLearning/ReadClassificationData.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "../RandomNumberGeneration/Statistics.h"
#include "../ExternalMemoryAlgorithms/File.h"
using namespace igmdk;

template<typename LEARNER> int testNumericalClassifier()
{
    DEBUG("Started Reading");
    typedef InMemoryData<NUMERIC_X, int> T;
    Vector<T> dataM(50);//make many enough to avoid ref realloc
    Vector<pair<PermutedData<T>, PermutedData<T> > > data;

    dataM.append(T());
    readIrisData(dataM.lastItem());
    data.append(makeData<T>(dataM));
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
            addToCSV(Vector<string>(1, toStringDouble(cost)), fCost.c_str());
            addToCSV(Vector<string>(1, toStringDouble(timediff)), fTimer.c_str());*/
        }
        else
        {
            int start = clock();
            LEARNER NeuralNetworkIris(data[i].first);
            ClassifierStats cs(evaluateConfusion(evaluateLearner<int>(NeuralNetworkIris, data[i].second)));
            double timediff = 1.0 * (clock() - start)/CLOCKS_PER_SEC;
            cs.debug();
        }
        //system("PAUSE");
    }
    return 0;
}

int main(int argc, char *argv[])
{
    testNumericalClassifier<AdaBoostSamme<> >();
    testNumericalClassifier<MRMRLearner<DecisionTree> >();
    testNumericalClassifier<LSHRange01>();
	return 0;
}


