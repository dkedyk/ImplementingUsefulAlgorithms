#ifndef IGMDK_COMPUTATIONAL_GEOMETRY_TEST_AUTO_H
#define IGMDK_COMPUTATIONAL_GEOMETRY_TEST_AUTO_H
#include "KDTree.h"
#include "Point.h"
#include <cassert>
using namespace std;

namespace igmdk{


template<typename MAP2D> void testPointMapAutoHelper(MAP2D& trie)
{
    int n = 100000;
    Vector<int> permutation(n);
    for(int i = 0; i < n; ++i) permutation[i] = i;
    GlobalRNG().randomPermutation(permutation.getArray(), n);
    for(int j = 0; j < n; ++j)
    {
        int i = permutation[j];
        Vector<int> key(2, i);
        trie.insert(key, -i);
        assert(trie.find(key));
        assert(*trie.find(key) == -i);
    }
}
template<typename N> void testAutoCheckRange(Vector<N*> const& v, int from,
    int to)
{
    for(int i = from; i < to; ++i)
    {
        bool found = false;
        for(int j = 0; j < v.getSize(); ++j)
            if(v[j]->key[0] == i)
            {
                found = true;
                break;
            }
        assert(found);
    }
}
void testKDTreeAuto()
{
    DEBUG("testKDTreeAuto");
    typedef KDTree<Vector<int>, int> TREE;
    TREE tree(2);
    testPointMapAutoHelper(tree);
    EuclideanDistance<Vector<int> >::DistanceIncremental di;
    Vector<int> p1(2, 100);
    //radius test; 7 * 7 + 7 * 7 = 98, so just makes it in 10^2
    testAutoCheckRange(tree.distanceQuery(p1, 100, di), 93, 107);
    //nn tests
    typename TREE::NodeType* nn = tree.nearestNeighbor(p1, di);
    assert(nn && nn->key == p1);
    testAutoCheckRange(tree.kNN(p1, 5, di), 98, 102);
    //range test
    testAutoCheckRange(tree.rangeQuery(Vector<int>(2, 97), Vector<int>(2, 103),
        Vector<bool>(2, true)), 97, 103);
    DEBUG("testKDTreeAuto passed");
}

void testVPTreeAuto()
{
    DEBUG("testVPTreeAuto");
    typedef VpTree<Vector<int>, int, EuclideanDistance<Vector<int> >::Distance>
        TREE;
    TREE tree;
    testPointMapAutoHelper(tree);
    Vector<int> p1(2, 100);
    //radius test; 7 * 7 + 7 * 7 = 98, so just makes it in 10^2
    testAutoCheckRange(tree.distanceQuery(p1, 10), 93, 107);
    //nn tests
    typename TREE::NodeType* nn = tree.nearestNeighbor(p1);
    assert(nn && nn->key == p1);
    testAutoCheckRange(tree.kNN(p1, 5), 98, 102);
    DEBUG("testVPTreeAuto passed");
}

void testAllAutoComputationalGeometry()
{
    DEBUG("testAllAutoComputationalGeometry");
    testKDTreeAuto();
    testVPTreeAuto();
}

}//end namespace
#endif
