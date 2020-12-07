#ifndef IGMDK_K_BIT_VECTOR_H
#define IGMDK_K_BIT_VECTOR_H
#include "../Utils/Bitset.h"
using namespace std;
namespace igmdk{

template<int N, typename WORD = unsigned long long> class KBitWordVector
{
    Bitset<WORD> b;
public:
    unsigned long long getSize()const{return b.getSize()/N;}
    KBitWordVector(){};
    KBitWordVector(int n, WORD item = 0): b(n * N)
    {
        if(Bits::getValue(item, 0, N) == 0) b.setAll(0);
        else for(unsigned long long i = 0; i < getSize(); ++i) set(item, i);
    }
    WORD operator[](unsigned long long i)const
        {assert(i < getSize()); return b.getValue(i * N, N);}
    void set(WORD value, unsigned long long i)
        {assert(i < getSize()); b.setValue(value, i * N, N);}
    void append(WORD value){b.appendValue(value, N);}
};

}//end namespace
#endif
