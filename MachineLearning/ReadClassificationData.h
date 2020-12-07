#ifndef IGMDK_READCLASSIFICATIONDATA_H
#define IGMDK_READCLASSIFICATIONDATA_H
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "../ExternalMemoryAlgorithms/File.h"
namespace igmdk{

template<typename DATA> void readWineData(DATA& result)
{
    ifstream fin("Datasets/wine.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        int label;
        data >> label;
        --label;
        if(label < 0 || label > 2) continue;
        data >> comma;
        for(int i = 0; i < 12; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        result.addZ(x, label);
    }
}

template<typename DATA> void readStatlogData(DATA& result, bool isTrain)
{
    ifstream fin(isTrain ? "Datasets/sat.trn" : "Datasets/sat.tst");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        for(int i = 0; i < 36; i++)
        {
            double currentDigit;
            data >> currentDigit;
            x.append( currentDigit );
        }
        int label;
        data >> label;
        if(label == 7) label = 6;
        if(label < 0 || label > 6) continue;
        result.addZ(x, label);
    }
}

template<typename DATA> void readSpectData(DATA& result, bool isTrain)
{
    ifstream fin(isTrain ? "Datasets/SPECT.train" : "Datasets/SPECT.test");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        int label;
        data >> label;
        data >> comma;
        for(int i = 0; i < 22; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        result.addZ(x, label);
    }
}

template<typename DATA> void readSpamData(DATA& result)
{
    ifstream fin("Datasets/spambase.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        for(int i = 0; i < 57; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        int label;
        data >> label;
        if(label != 0 && label != 1) continue;
        result.addZ(x, label);
    }
}

template<typename DATA> void readPimaData(DATA& result)
{
    ifstream fin("Datasets/pima-indians-diabetes.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        for(int i = 0; i < 8; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        int label;
        data >> label;
        if(label != 0 && label != 1) continue;
        result.addZ(x, label);
    }
}

template<typename DATA> void readLetterData(DATA& result)
{
    ifstream fin("Datasets/letter-recognition.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        char sLabel;
        data >> sLabel;
        int label = sLabel - 'A';
        if(label < 0 || label > 26) continue;
        data >> comma;
        for(int i = 0; i < 16; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        result.addZ(x, label);
    }
}

template<typename DATA> void readIonosphereData(DATA& result)
{
    ifstream fin("Datasets/ionosphere.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        for(int i = 0; i < 34; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        string sLabel;
        data >> sLabel;
        int label;
        if(sLabel == "g") label = 0;
        else if(sLabel == "b") label = 1;
        else continue;
        result.addZ(x, label);
    }
}

template<typename DATA> void readGlassData(DATA& result)
{
    ifstream fin("Datasets/glass.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        for(int i = 0; i < 10; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        int label;
        data >> label;
        if(label > 7 || label < 1) continue;
        result.addZ(x, label - 1);
    }
}

template<typename DATA> void readMadelonData(DATA& result, bool isTrain)
{;
    ifstream fin(isTrain ? "Datasets/madelon_train.data" : "Datasets/madelon_valid.data");
    ifstream fin2(isTrain ? "Datasets/madelon_train.labels" : "Datasets/madelon_valid.labels");
    while(!fin.eof() && !fin2.eof())
    {
        string line, line2;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        for(int i = 0; i < 500; i++)
        {
            double currentDigit;
            data >> currentDigit;

            x.append( currentDigit );
        }
        assert(x.getSize() == 500);
        getline(fin2, line2, '\n');
        stringstream data2;
        data2 << line2;
        int dummy, label;
        data2 >> dummy;
        if(dummy == -1) label = 0;
        else if(dummy == 1) label = 1;
        else continue;
        result.addZ(x, label);
    }
}

template<typename DATA> void readArceneData(DATA& result, bool isTrain)
{
    ifstream fin(isTrain ? "Datasets/arcene_train.data" : "Datasets/arcene_valid.data");
    ifstream fin2(isTrain ? "Datasets/arcene_train.labels" : "Datasets/arcene_valid.labels");
    while(!fin.eof() && !fin2.eof())
    {
        string line, line2;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        for(int i = 0; i < 10000; i++)
        {
            double currentDigit;
            data >> currentDigit;

            x.append( currentDigit );
        }
        assert(x.getSize() == 10000);
        getline(fin2, line2, '\n');
        stringstream data2;
        data2 << line2;
        int dummy, label;
        data2 >> dummy;
        if(dummy == -1) label = 0;
        else if(dummy == 1) label = 1;
        else continue;
        result.addZ(x, label);
    }
}

template<typename DATA> void readWDBCData(DATA& result)
{
    ifstream fin("Datasets/wdbc.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        int dummy;
        data >> dummy;
        data >> comma;
        char sLabel;
        data >> sLabel;
        int label;
        if(sLabel == 'B') label = 0;
        else if(sLabel == 'M') label = 1;
        else continue;
        data >> comma;
        for(int i = 0; i < 30; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        assert(x.getSize() == 30);
        result.addZ(x, label);
    }
}

template<typename DATA> void readBanknoteData(DATA& result)
{
    ifstream fin("Datasets/data_banknote_authentication.txt");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        for(int i = 0; i < 4; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        int label;
        data >> label;
        assert(x.getSize() == 4);
        result.addZ(x, label);
    }
}

template<typename DATA> void readCNEAData(DATA& result)
{
    ifstream fin("Datasets/CNAE-9.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        int label;
        data >> label;
        data >> comma;
        for(int i = 0; i < 856; i++)
        {
            double currentDigit;
            data >> currentDigit;
            data >> comma;
            x.append( currentDigit);
        }
        assert(x.getSize() == 856);
        result.addZ(x, label - 1);
    }
}

template<typename DATA> void readIrisData(DATA& result)
{
    ifstream fin("Datasets/iris.data");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        for(int i = 0; i < 4; i++)
        {
            double currentDigit;
            data >> currentDigit;

            data >> comma;
            x.append( currentDigit );
        }
        string sLabel;
        data >> sLabel;
        int label;
        if(sLabel == "Iris-setosa") label = 0;
        else if(sLabel == "Iris-versicolor") label = 1;
        else if(sLabel == "Iris-virginica") label = 2;
        else continue;
        result.addZ(x, label);
    }
}


template<typename DATA> void readDigitData(DATA& result, bool isTrain)
{
    ifstream fin(isTrain? "optdigits.tra" : "optdigits.tes");
    while(!fin.eof())
    {
        string line;
        getline(fin, line, '\n');
        stringstream data;
        data << line;
        NUMERIC_X x;
        char comma;
        for(int i = 0; i < 64; i++)
        {
            double currentDigit;
            data >> currentDigit;
            data >> comma;
            x.append( currentDigit );
        }
        int label;
        data >> label;
        result.addZ(x, label);
    }
}

template<typename T> pair<PermutedData<T>, PermutedData<T> >
makeData(Vector<T>& dataM)
{
    return createTrainingTestSetsStatified(dataM.lastItem());
}

template<typename T> pair<PermutedData<T>, PermutedData<T> >
makeDataDivided(Vector<T>& dataM)
{
    PermutedData<T> p1(dataM[dataM.getSize() - 2]);
    for(int i = 0; i < dataM[dataM.getSize() - 2].getSize(); ++i) p1.addIndex(i);
    PermutedData<T> p2(dataM.lastItem());
    for(int i = 0; i < dataM.lastItem().getSize(); ++i) p2.addIndex(i);
    return make_pair(p1, p2);
}

}//end namespace
#endif
