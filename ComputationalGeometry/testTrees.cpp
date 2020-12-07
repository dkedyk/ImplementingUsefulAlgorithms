#include <iostream>
#include <cmath>
#include "KDTree.h"
#include "Point.h"
#include "ComputationalGeometryTestAuto.h"
#include "../RandomNumberGeneration/Random.h"
#include "../NumericalMethods/NumericalMethods.h"
#include "../HashTable/ChainingHashTable.h"
#include "../RandomNumberGeneration/Random.h"
#include "../RandomNumberGeneration/Statistics.h"
using namespace igmdk;

template<typename KEY, typename VALUE, typename DISTANCE> class KNNBruteForce
{
    DISTANCE distance;
    typedef KVPair<KEY, VALUE> Node;
    Vector<Node> nodes;
    struct QNode
    {
        double distance;
        int result;
        bool operator<(QNode const& rhs)const
            {return distance > rhs.distance;}
    };
public:
    KNNBruteForce(DISTANCE const& theDistance = DISTANCE()):
        distance(theDistance){}
    typedef Node NodeType;
    void insert(KEY const& key, VALUE const& value)
        {nodes.append(Node(key, value));}
    Vector<NodeType*> kNN(KEY const& key, int k)
    {
        Heap<QNode> q;
        for(int i = 0; i < nodes.getSize(); ++i)
        {
            QNode node = {distance(key, nodes[i].key), i};
            if(q.getSize() < k) q.insert(node);
            else if(node.distance < q.getMin().distance)
                q.changeKey(0, node);
        }
        Vector<NodeType*> result;
        while(!q.isEmpty()) result.append(&nodes[q.deleteMin().result]);
        result.reverse();
        return result;
    }
    NodeType* nearestNeighbor(KEY const& key){return kNN(key, 1)[0];}
};

void testKD3()
{
    KDTree<Point<double, 100>, bool> kdtree(100);
    int D = 100;
    int N = 100000;
    for(int i = 0; i < N; ++i)
    {
        Point<double, 100> x;
        for(int j = 0; j < D; ++j) x[j] = j;
        kdtree.insert(x, true);
    }
    for(int i = 0; i < N; ++i)
    {
        Point<double, 100> x;
        for(int j = 0; j < D; ++j) x[j] = j + GlobalRNG().uniform01();
        assert(kdtree.nearestNeighbor(x, EuclideanDistance<Point<double, 100> >::DistanceIncremental()));
    }
}

void testKNNBF()
{
    KNNBruteForce<Point<double, 100>, bool, EuclideanDistance<Point<double, 100> >::Distance> kdtree;
    int D = 100;
    int N = 100000;
    for(int i = 0; i < N; ++i)
    {
        Point<double, 100> x;
        for(int j = 0; j < D; ++j) x[j] = j;
        kdtree.insert(x, true);
    }
    for(int i = 0; i < N; ++i)
    {
        Point<double, 100> x;
        for(int j = 0; j < D; ++j) x[j] = j + GlobalRNG().uniform01();
        assert(kdtree.nearestNeighbor(x));
    }
}

void DDDVPTree()
{
    Random<> rng(0);
    VpTree<Point<double>, int, EuclideanDistance<Point<double> >::DistanceIncremental> VPTree0to9;
    int D = 2;
    for(int i = 0; i < 10; ++i)
    {
        Point<double, 2> x;
        for(int j = 0; j < D; ++j) x[j] = rng.uniform01();
        VPTree0to9.insert(x, i);
    }

    cout << "breakpoint" << endl;
}

void DDDKDTree()
{
    Random<> rng(0);
    KDTree<Point<double, 2>, int> KDTree0to9(2);
    int D = 2;
    for(int i = 0; i < 10; ++i)
    {
        Point<double, 2> x;
        for(int j = 0; j < D; ++j) x[j] = rng.uniform01();
        KDTree0to9.insert(x, i);
    }

    cout << "breakpoint" << endl;
}

int main()
{
    DDDVPTree();
    DDDKDTree();
    testAllAutoComputationalGeometry();
    /*testKD3();//very fast
    //testKNNBF();//very slow*/
	return 0;
}
