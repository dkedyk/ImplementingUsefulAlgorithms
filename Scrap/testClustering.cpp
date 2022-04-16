#include "../MachineLearning/Classification.h"
#include "ScrapClustering.h"
#include "../MachineLearning/ReadClassificationData.h"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "../RandomNumberGeneration/Statistics.h"
#include "../ExternalMemoryAlgorithms/File.h"
using namespace igmdk;



template<typename CLUSTERER> int testNumericalClusterer()
{
    DEBUG("Started Reading");
    typedef InMemoryData<NUMERIC_X, int> T;
    Vector<T> dataM(50);//make many enough to avoid ref realloc
    Vector<pair<PermutedData<T>, PermutedData<T> > > data;

    dataM.append(T());
    readIrisData(dataM.lastItem());//iris data duplicated here from ML directory to use same code
    data.append(makeData<T>(dataM));

    DEBUG("Done Reading");
    for(int i = 0; i < data.getSize(); ++i)
    {
        int start = clock();
        CLUSTERER c;
        //Vector<int> result = c(data[i].first, findNClasses(data[i].first));
        //ScalerMinMax s(data[i].first);
        //ScalerMQ s(data[i].first);
        //ScaledData<PermutedData<T>, ScalerMQ> sd(data[i].first, s);
        ScalerMinMax s(data[i].first);
        ScaledData<PermutedData<T>, ScalerMinMax> sd(data[i].first, s);
        //Vector<int> result = c(sd).assignments;
        Vector<int> result = c(sd, findNClasses(sd)).assignments;

        Matrix<int> counts = clusterContingencyMatrix(result, sd);
        double purity = clusterPurity(counts);
        double aRand = AdjustedRandIndex(counts);
        double relkDiff = (counts.rows - counts.columns)*1.0/counts.columns;
        double cAcc = clusterClassificationAccuracy(counts);

        DEBUG(purity);
        DEBUG(aRand);
        DEBUG(relkDiff);
        DEBUG(cAcc);

        double timediff = 1.0 * (clock() - start)/CLOCKS_PER_SEC;
        DEBUG(timediff);
    }
    return 0;
}

int main(int argc, char *argv[])
{
    testNumericalClusterer<EMSmart>();
    testNumericalClusterer<DBSCAN<> >();
    //testNumericalClusterer<HierarchicalClustering<>::Functor>();
	return 0;
}


