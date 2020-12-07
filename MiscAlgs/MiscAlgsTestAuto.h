#ifndef IGMDK_MISC_ALGORITHMS_TEST_AUTO_H
#define IGMDK_MISC_ALGORITHMS_TEST_AUTO_H
#include "LRUCache.h"
#include "PrimeTable.h"
#include "CombinatorialGeneration.h"
#include "KBitWordVector.h"

using namespace std;

namespace igmdk{

void testLRUAuto()
{
    DEBUG("testLRUAuto");
    int n = 1000; int k = 5;
    LRUCache<int, int> c(k);
    for(int i = 0; i < n; ++i) c.write(i, i);
    for(LRUCache<int, int>::Iterator i = c.begin(); i != c.end(); ++i)
        assert(i->value >= n - k && i->value < n);
    for(int i = n - k; i < n; ++i)
    {
        assert(c.read(i));
        assert(*c.read(i) == i);
    }
    for(int i = 0; i < n - k; ++i)
    {
        assert(!c.read(i));
    }
    DEBUG("testLRUAuto passed");
}

void testPrimeTableAuto()
{
    DEBUG("testPrimeTableAuto");
    int smallPrimes[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47};
    Vector<bool> primality(50, false);
    for(int i = 0; i < sizeof(smallPrimes)/sizeof(smallPrimes[0]); ++i)
        primality[smallPrimes[i]] = true;
    PrimeTable pt(primality.getSize());
    for(int i = 1; i < primality.getSize(); ++i)
        assert(pt.isPrime(i) == primality[i]);
    DEBUG("testPrimeTableAuto passed");
}

void testKBitVectorAuto()
{
    DEBUG("testKBitWordVectorAuto");
    KBitWordVector<5> x;
    int N = 100000;
    for(int i = 0; i < N; ++i) x.append(i);
    for(int i = 0; i < N; ++i) x.set(i, i);
    for(int i = 0; i < N; ++i) assert(x[i] == i % 32);
    DEBUG("testKBitVectorAuto passed");
}

void testAllAutoMiscAlgorithms()
{
    DEBUG("testAllAutoMicAlgorithms");
    testLRUAuto();
    testPrimeTableAuto();
    testKBitVectorAuto();
}

}//end namespace
#endif
