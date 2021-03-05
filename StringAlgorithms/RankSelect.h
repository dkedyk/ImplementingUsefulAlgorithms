#ifndef IGMDK_RANKSELECT_H
#define IGMDK_RANKSELECT_H

#include <cassert>
#include "../Utils/Bitset.h"
#include "../Sorting/Sort.h"
namespace igmdk{

int rank64(unsigned long long x, int i)
{//ith bit not included
    assert(0 <= i && i < numeric_limits<unsigned long long>::digits);
    //set bits >= i to 0 and popCount
    return popCountWord(x & ((1ull << i) - 1));
    //to include ith bit use i+1 instead of i
    //overflow is fine because result is -1
}
int select64(unsigned long long x, int i, bool is0 = false)
{
    assert(0 <= i && i < numeric_limits<unsigned long long>::digits);
    static PopCount8 p8;
    int result = 0, byteBits = numeric_limits<unsigned char>::digits;
    //use popCount to get to the byte with the desired bit position
    while(x)
    {
        int temp = p8(x & 0xff);
        if(is0) temp = byteBits - temp;
        if(i - temp < 0) break;
        result += byteBits;
        x >>= byteBits;
        i -= temp;
    }//scan the bits in the found byte for the desired bit position
    for(int j = 0; j <= byteBits; ++j)
    {
        bool temp = x & 1 << j;
        if(is0) temp = !temp;
        if(temp && i-- == 0) return result;
        ++result;
    }
    return -1;
}
class RankSelect
{
    enum{B = numeric_limits<unsigned long long>::digits};
    Bitset<unsigned long long> bitset;
    //cumulative counts of every last bit in a word
    Vector<unsigned long long> counts;
    long long counts0(long long i){return (i + 1) * B - counts[i];}
    long long getCount(long long i, bool is0)
        {return is0 ? (i + 1) * B - counts[i] : counts[i];}
public:
    RankSelect(unsigned long long initialSize = 0): bitset(initialSize){}
    Bitset<unsigned long long>& getBitset(){return bitset;}
    void finalize()
    {
        counts = bitset.getStorage();
        for(long long i = 0; i < bitset.wordSize(); ++i) counts[i] =
            popCountWord(bitset.getStorage()[i]) + (i == 0 ? 0 : counts[i-1]);
    }
    long long rank(long long i)
    {
        assert(0 <= i && i < bitset.getSize());
        long long index = i / B;
        return (index > 0 ? counts[index - 1] : 0) +
            rank64(bitset.getStorage()[index], i % B);
    }
    long long rank0(long long i){return i - rank(i);}
    long long select(long long i, bool is0 = false)
    {
        assert(0 <= i && i < bitset.getSize());
        long long left = 0, right = bitset.wordSize() - 1;
        while(left < right)
        {
            long long middle = (left + right)/2;
            if(getCount(middle, is0) <= i) left = middle + 1;
            else right = middle - 1;
        }
        long long result = select64(bitset.getStorage()[left],
            i - (left == 0 ? 0 : getCount(left - 1, is0)), is0);
        if(result != -1) result += left * B;
        return result;
    }
    long long select0(long long i){return select(i, true);}
    long long prev(long long i){return select(rank(i) - 1);}
    long long prev0(long long i){return select0(rank0(i) - 1);}
    long long next(long long i){return select(rank(i));}
    long long next0(long long i){return select0(rank0(i));}
};

class BinaryTree
{//this does not extend to d-ary trees unlike the heaps
    //which need at least log(d)n bits to be distinguished
    RankSelect rs;
    int convert(int i){return rs.getBitset()[i] ? rs.rank(i) : -1;}
public:
    void addNodeInLevelOrder(bool isNotExternalDummyLeaf)
        {rs.getBitset().append(isNotExternalDummyLeaf);}
    void finalize(){rs.finalize();}
    int parent(int i){return (rs.select(i)+1)/2-1;}
    int leftChild(int i){return convert(2 * rs.rank(i) + 1);}
    int rightChild(int i){return convert(2 * rs.rank(i) + 2);}
};
class OrdinalTree
{
    RankSelect rs;
    int convert(int i){return rs.getBitset()[i] ? rs.rank(i) : -1;}
public:
    OrdinalTree(){rs.getBitset().append(1);rs.getBitset().append(0);}
    void addNodeInLevelOrder(int nChildren)
    {
        for(int i = 0; i < nChildren; ++i) rs.getBitset().append(1);
        rs.getBitset().append(0);
    }
    void finalize(){rs.finalize();}
    int parent(int i){return rs.select(i) - i - 1;}
    int firstChild(int i){return convert(rs.select0(i) + 1);}
    int nextChild(int i) {return convert(rs.select(i) + 1);}
};

}
#endif
