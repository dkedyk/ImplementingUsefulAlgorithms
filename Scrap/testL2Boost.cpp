#include "L2Boost.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "../RandomNumberGeneration/Statistics.h"
#include "../ExternalMemoryAlgorithms/File.h"
using namespace igmdk;

template<typename DATA> void readEnergyHeatData(DATA& result)
{
    //the data format is 5 space-separated values per line
    //the last value is y
    ifstream fin("../MachineLearning/RegDatasets/ENB2012_dataConverted.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        for(int i = 0; i < 8; i++)
        {
            double currentDigit;
            data >> currentDigit;
            x.append( currentDigit );
        }
        double label;
        data >> label;
        result.addZ(x, label);
    }
}

template<typename LEARNER> int testRegressor()
{
    DEBUG("Started Reading");
    typedef InMemoryData<NUMERIC_X, double> T;
    Vector<T> data;
    data.append(T());
    readEnergyHeatData(data.lastItem());
    DEBUG("Done Reading");
    for(int i = 0; i < data.getSize(); ++i)
    {
        pair<PermutedData<T>, PermutedData<T> > tt(createTrainingTestSetsDetPerm(data[i]));
        int start = clock();
        LEARNER NeuralNet(tt.first);
        RegressionStats cs = evaluateRegressor(evaluateLearner<T::Y_TYPE>(NeuralNet, tt.second));
        double timediff = 1.0 * (clock() - start)/CLOCKS_PER_SEC;
        cs.debug();
    }
    return 0;
}

void testRegressors()
{
    DEBUG("L2Boost");
    testRegressor<L2Boost<> >();
}

int main(int argc, char *argv[])
{
    for(int i = 0; i < 1; ++i) testRegressors();
	return 0;
}


