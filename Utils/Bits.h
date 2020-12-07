#ifndef IGMDK_BITS_H
#define IGMDK_BITS_H
#include <cassert>
#include <limits>
#include "../Utils/Utils.h"
using namespace std;
namespace igmdk{

unsigned long long twoPower(int x){return 1ull << x;}
//e.g., for 8 have !(0100 & 0011); careful: returns true for 0
bool isPowerOfTwo(unsigned long long x){return !(x & (x - 1));}
int lgFloor(unsigned long long x)
{//shift until 0 and count; e.g., lgFloor(2) = lgFloor(3) = 1
    assert(x);//log of 0 is undefined
    int result = 0;
    while(x >>= 1) ++result;
    return result;
}
//E.g., lgCeiling(3) = 2
int lgCeiling(unsigned long long x){return lgFloor(x) + !isPowerOfTwo(x);}
//E.g., nextPowerOfTwo(7) = nextPowerOfTwo(8) = 8
unsigned long long nextPowerOfTwo(unsigned long long x)
    {return isPowerOfTwo(x) ? x : twoPower(lgFloor(x) + 1);}

namespace Bits
{//below useful due to type because literals are int
unsigned long long const ZERO = 0, FULL = ~ZERO;
bool get(unsigned long long x, int i){return x & twoPower(i);}//and
bool flip(unsigned long long x, int i){return x ^ twoPower(i);}//xor
template<typename WORD> void set(WORD& x, int i, bool value)
{
    assert(!numeric_limits<WORD>::is_signed);
    if(value) x |= twoPower(i);
    else x &= ~twoPower(i);
}
unsigned long long upperMask(int n){return FULL << n;}//11110000
unsigned long long lowerMask(int n){return ~upperMask(n);}//00001111
unsigned long long middleMask(int i, int n){return lowerMask(n)<<i;}//00111000
unsigned long long getValue(unsigned long long x, int i, int n)
    {return (x >> i) & lowerMask(n);}//shift down and mask
template<typename WORD>
void setValue(WORD& x, unsigned long long value, int i, int n)
{
    assert(!numeric_limits<WORD>::is_signed);
    WORD mask = middleMask(i, n);//note the cast
    x &= ~mask;//erase
    x |= mask & (value << i);//set
}
}//end namespace

class PopCount8
{
    char table[256];
public:
    static int popCountBruteForce(unsigned long long x)
    {//for computing the table
        int count = 0;
        while(x)
        {//count rightmost bits one-by-one
            count += x & 1;
            x >>= 1;
        }
        return count;
    }
    PopCount8(){for(int i = 0; i < 256; ++i) table[i] = popCountBruteForce(i);}
    int operator()(unsigned char x)const{return table[x];}
};
int popCountWord(unsigned long long x)//loses no efficiency over template
{//initialization on first call
    static PopCount8 p8;
    int result = 0;
    for(; x; x >>= 8) result += p8(x);//equal performance for every word type
    return result;
}

int rightmost0Count(unsigned long long x){return popCountWord(~x & (x - 1));}

class ReverseBits8
{
    unsigned char table[256];
public:
    template<typename WORD> static WORD reverseBitsBruteForce(WORD x)
    {
        assert(!numeric_limits<WORD>::is_signed);
        WORD result = 0;
        for(int i = 0; i < numeric_limits<WORD>::digits; ++i)
        {
            result = (result << 1) + (x & 1);
            x >>= 1;
        }
        return result;
    }
    ReverseBits8()
    {
        for(int i = 0; i < 256; ++i)
            table[i] = reverseBitsBruteForce<unsigned char>(i);
    }
    unsigned char operator()(unsigned char x)const{return table[x];}
};
template<typename WORD> WORD reverseBits(WORD x)
{
    assert(!numeric_limits<WORD>::is_signed);
    static ReverseBits8 r8;
    WORD result = 0;
    for(int i = 0; i < sizeof(x); ++i, x >>= 8) result = (result << 8) + r8(x);
    return result;
}
template<typename WORD> WORD reverseBits(WORD x, int n)
{
    int shift = sizeof(x) * 8 - n;
    assert(!numeric_limits<WORD>::is_signed && n > 0 && shift >= 0);
    return reverseBits<WORD>(x & Bits::lowerMask(n)) >> shift;
}

}//end namespace
#endif
