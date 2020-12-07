#ifndef IGMDK_BLOOM_FILTER_H
#define IGMDK_BLOOM_FILTER_H

#include <cassert>
#include "../Utils/Bitset.h"
#include "HashFunction.h"
using namespace std;

namespace igmdk{

template<typename KEY, typename HASHER = EHash<BUHash> >
class BloomFilter
{
    Bitset<unsigned char> items;//must be before h1, h2
    HASHER h1, h2;
    int nHashes;
    int hash(int hash1, int hash2, int i)
    {
        if(i == 0) return hash1;
        if(i == 1) return hash2;
        return (hash1 + i * hash2) % items.getSize();
    }
public:
    BloomFilter(int m, int theNHashes = 7): nHashes(theNHashes), items(
        nextPowerOfTwo(m)), h1(items.getSize()), h2(items.getSize())
        {assert(m > 0 && theNHashes > 0);}
    void insert(KEY const& key)
    {
        int hash1 = h1(key), hash2 = h2(key);
        for(int i = 0; i < nHashes; ++i) items.set(hash(hash1, hash2, i));
    }
    bool isInserted(KEY const& key)
    {
        int hash1 = h1.hash(key), hash2 = h2.hash(key);
        for(int i = 0; i < nHashes; ++i)
            if(!items[hash(hash1, hash2, i)]) return false;
        return true;
    }
};

}//end namespace
#endif
