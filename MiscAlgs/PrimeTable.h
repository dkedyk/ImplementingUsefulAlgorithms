#ifndef IGMDK_PRIME_TABLE_H
#define IGMDK_PRIME_TABLE_H
#include "../Utils/Bitset.h"

namespace igmdk{

class PrimeTable
{
    long long maxN;
    Bitset<> table;//marks odd numbers starting from 3
    long long nToI(long long n)const{return (n - 3)/2;}
public:
    PrimeTable(long long primesUpto): maxN(primesUpto - 1),
        table(nToI(maxN) + 1)
    {
        assert(primesUpto > 1);
        table.setAll(true);
        for(long long i = 3; i <= sqrt(maxN); i += 2)
            if(isPrime(i))//set every odd multiple i <= k <= maxN/i to false
                for(long long k = i; i * k <= maxN; k += 2)
                    table.set(nToI(i * k), false);
    }
    bool isPrime(long long n)const
    {
        assert(n > 0 && n <= maxN);
        return n == 2 || (n > 2 && n % 2 && table[nToI(n)]);
    }
};

}//end namespace
#endif
