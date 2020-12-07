#include "Regression.h"
#include "KNNRegression.h"
#include "NeuralNetworkRegression.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "../RandomNumberGeneration/Statistics.h"
#include "../ExternalMemoryAlgorithms/File.h"
using namespace igmdk;

template<typename DATA> void readAirfoilData(DATA& result)
{
    //the data format is 6 space-separated values per line
    //the last value is y
    ifstream fin("RegDatasets/airfoil_self_noise.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        for(int i = 0; i < 5; i++)
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

template<typename DATA> void readHardwareData(DATA& result)
{
    //the data format is 10 comma-separated values per line
    //the last value is y, first 2 irrelevant
    ifstream fin("RegDatasets/machine.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        getline(data, line, ',');
        getline(data, line, ',');
        for(int i = 0; i < 9; i++)
        {
            double currentDigit;
            getline(data, line, ',');
            stringstream data2;
            data2 << line;
            data2 >> currentDigit;
            x.append( currentDigit );
        }
        double label;
        getline(data, line, ',');
        stringstream data2;
        data2 << line;
        data2 >> label;
        result.addZ(x, label);
    }
}

template<typename DATA> void readConcreteData(DATA& result)
{
    //the data format is 9 space-separated values per line
    //the last value is y
    ifstream fin("RegDatasets/Concrete_Data.data");
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

template<typename DATA> void readCASPData(DATA& result)
{
    //the data format is 10 space-separated values per line
    //the last value is y
    ifstream fin("RegDatasets/CASP.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        double label;
        data >> label;
        for(int i = 0; i < 9; i++)
        {
            double currentDigit;
            data >> currentDigit;
            x.append( currentDigit );
        }
        result.addZ(x, label);
    }
}

template<typename DATA> void  readYachtData(DATA& result)
{
    //the data format is 10 space-separated values per line
    //the last value is y
    ifstream fin("RegDatasets/yacht_hydrodynamics.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        for(int i = 0; i < 6; i++)
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

template<typename DATA> void readSkillData(DATA& result)
{
    //the data format is 10 space-separated values per line
    //the last value is y
    ifstream fin("RegDatasets/SkillCraft1Converted.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        double label;
        data >> label;
        for(int i = 0; i < 18; i++)
        {
            double currentDigit;
            data >> currentDigit;
            x.append( currentDigit );
        }
        result.addZ(x, label);
    }
}

template<typename DATA> void readCCPPData(DATA& result)
{
    //the data format is 5 space-separated values per line
    //the last value is y
    ifstream fin("RegDatasets/Folds5x2_ppConverted.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        for(int i = 0; i < 4; i++)
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

template<typename DATA> void  readWineData(DATA& result)
{
    //the data format is 5 space-separated values per line
    //the last value is y
    ifstream fin("RegDatasets/winequality-whiteConverted.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        for(int i = 0; i < 11; i++)
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

template<typename DATA> void readEnergyHeatData(DATA& result)
{
    //the data format is 5 space-separated values per line
    //the last value is y
    ifstream fin("RegDatasets/ENB2012_dataConverted.data");
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
    /*data.append(T());
    readAirfoilData(data.lastItem());
    data.append(T());
    readHardwareData(data.lastItem());
    data.append(T());
    readConcreteData(data.lastItem());
    data.append(T());
    readCASPData(data.lastItem());
    data.append(T());
    readYachtData(data.lastItem());
    data.append(T());
    readSkillData(data.lastItem());
    data.append(T());
    readCCPPData(data.lastItem());
    data.append(T());
    readWineData(data.lastItem());*/
    data.append(T());
    readEnergyHeatData(data.lastItem());
    DEBUG("Done Reading");
    int reportNumber = time(0);
    string suffix = to_string(reportNumber) + ".csv";
    string fExpStd = "RegReportExpStd" + suffix;
    string fL1 = "RegReportL1" + suffix;
    string fTimer = "RegReportTimer" + suffix;
    string fCount = "RegeportFeature" + suffix;
    ++reportNumber;
    for(int i = 0; i < data.getSize(); ++i)
    {
        pair<PermutedData<T>, PermutedData<T> > tt(createTrainingTestSetsDetPerm(data[i]));
        int start = clock();
        LEARNER NeuralNet(tt.first);
        RegressionStats cs = evaluateRegressor(evaluateLearner<T::Y_TYPE>(NeuralNet, tt.second));
        double timediff = 1.0 * (clock() - start)/CLOCKS_PER_SEC;
        cs.debug();
        /*addToCSV(Vector<string>(1, toStringDouble(cs.expStd)), fExpStd.c_str());
        addToCSV(Vector<string>(1, toStringDouble(cs.l1Err)), fL1.c_str());
        addToCSV(Vector<string>(1, toStringDouble(timediff)), fTimer.c_str());*/
        //addToCSV(Vector<string>(1, to_string(s.model.f.fMap.getSize())), fCount.c_str());
    }
    return 0;
}

void testRegressors()
{
    DEBUG("SLasso>");
    testRegressor<SLasso>();
    DEBUG("RegressionTree");
    testRegressor<RegressionTree>();
    DEBUG("RandomForestReg");
    for(int i = 0; i < 5; ++i)
        testRegressor<RandomForestReg>();
    DEBUG("SKNNReg>");
    testRegressor<SKNNReg>();
    DEBUG("SNNReg>");
    testRegressor<SNNReg>();
    DEBUG("SmartFSLearnerReg>");
    testRegressor<SmartFSLearnerReg<> >();
    DEBUG("SimpleBestCombinerReg");
    testRegressor<SimpleBestCombinerReg>();
    DEBUG("SRaceLasso");
    testRegressor<SRaceLasso>();
}

int main(int argc, char *argv[])
{
    for(int i = 0; i < 1; ++i) testRegressors();
	return 0;
}


