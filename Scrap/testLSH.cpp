#include <iostream>
#include <cmath>
#include "LSH.h"
#include "../ComputationalGeometry/Point.h"
#include "../RandomNumberGeneration/Random.h"
#include "../NumericalMethods/NumericalMethods.h"
using namespace igmdk;

void testLSH()
{
    //DEBUG(1/exp(1));
    //DEBUG(E2LSHHasher::p(1, 1));
    int D = 2;
    LSH<E2LSHHasher> tree = buildE2LSH(D, 1, 1, D * 5);
    int N = 1000000;
    for(int i = 0; i < N; ++i)
    {
        tree.insert(Vector<double>(2, i));
    }
    for(int i = 0; i < 1; ++i)
    {
        Vector<Vector<double> > neighbors = tree.cNeighbors(Vector<double>(2, i));
        DEBUG(neighbors.getSize());
        for(int j = 0; j < neighbors.getSize(); ++j)
        {
            DEBUG(j);
            for(int k = 0; k < neighbors[j].getSize(); ++k) DEBUG(neighbors[j][k]);
        }
    }
}

void testLSH2()
{
    //DEBUG(1/exp(1));
    //DEBUG(E2LSHHasher::p(1, 1));
    int D = 2;
    NearestNeighborLSH<E2LSHHasher> tree = buildE2NNLSH(D, 1, 10, D * 5);
    int N = 100000;
    for(int i = 0; i < N; ++i)
    {
        tree.insert(Vector<double>(2, i));
    }
    int noneCount = 0;
    for(int i = 0; i < N; ++i)
    {
        pair<Vector<double>, bool> neighbor = tree.cNeighbor(Vector<double>(2, i));
        //DEBUG(neighbor.second);
        if(neighbor.second)
        {
            //for(int k = 0; k < neighbor.first.getSize(); ++k) DEBUG(neighbor.first[k]);
        }
        else ++noneCount;
    }
    DEBUG(noneCount);
}

void testLSH3()
{
    //DEBUG(1/exp(1));
    //DEBUG(E2LSHHasher::p(1, 1));
    int D = 100;
    NearestNeighborLSH<E2LSHHasher> tree = buildE2NNLSH(D, 1, 10, D * 0.5);
    int N = 100000;
    for(int i = 0; i < N; ++i)
    {
        Vector<double> x;
        for(int j = 0; j < D; ++j) x.append(i);
        tree.insert(x);
    }
    int noneCount = 0;
    for(int i = 0; i < N; ++i)
    {
        Vector<double> x;
        for(int j = 0; j < D; ++j) x.append(i);
        pair<Vector<double>, bool> neighbor = tree.cNeighbor(x);
        //DEBUG(neighbor.second);
        if(neighbor.second)
        {
            //for(int k = 0; k < neighbor.first.getSize(); ++k) DEBUG(neighbor.first[k]);
        }
        else ++noneCount;
    }//20 secs
    DEBUG(noneCount);
}

int main()
{
    testLSH();
    testLSH2();
    testLSH3();
	return 0;
}
