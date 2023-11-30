#ifndef IGMDK_HASH_FUNCTION_H
#define IGMDK_HASH_FUNCTION_H
#include <string>
#include "../RandomNumberGeneration/Random.h"
#include "../Utils/Bits.h"
#include "../Utils/Bitset.h"
namespace igmdk{

class PrimeHash
{
    static uint32_t const PRIME = (1ull << 32) - 5;
    uint32_t seed;
public:
    PrimeHash(): seed((GlobalRNG().next() | 1) % PRIME) {}//ensure non-0
    typedef uint32_t WORD_TYPE;
    unsigned long long max()const{return PRIME - 1;}
    unsigned long long operator()(WORD_TYPE x)const
        {return (unsigned long long)seed * x % PRIME;}
    class Builder
    {
        unsigned long long sum;
        WORD_TYPE a;
        friend PrimeHash;
        Builder(WORD_TYPE theSeed): sum(0), a(theSeed) {}
    public:
        void add(WORD_TYPE xi)
        {//unlikely possible overflow from adding but that's ok
            sum += (unsigned long long)a * xi;
            a = xorshiftTransform(a);
        }
    };
    Builder makeBuilder()const{return Builder(seed);}
    unsigned long long operator()(Builder b)const{return b.sum % PRIME;}
};

class PrimeHash2
{
    static uint32_t const PRIME = (1ull << 32) - 5;
    uint32_t seed;
public:
    PrimeHash2(): seed((GlobalRNG().next() | 1) % PRIME) {}//ensure non-0
    typedef uint32_t WORD_TYPE;
    unsigned long long max()const{return PRIME - 1;}
    unsigned long long operator()(WORD_TYPE x)const
        {return (unsigned long long)seed * x % PRIME;}
    class Builder
    {
        unsigned long long sum;
        WORD_TYPE seed;
        friend PrimeHash2;
        Builder(WORD_TYPE theSeed): sum(0), seed(theSeed) {}
    public://unlikely possible overflow from add & mult but that's ok
        void add(WORD_TYPE xi){sum = seed * (sum + xi) % PRIME;}
    };
    Builder makeBuilder()const{return Builder(seed);}
    unsigned long long operator()(Builder b)const{return b.sum;}
};

template<typename HASHER> class MHash
{
    unsigned long long m;
    HASHER h;
public:
    MHash(unsigned long long theM): m(theM)
        {assert(theM > 0 && theM <= h.max());}
    typedef typename HASHER::WORD_TYPE WORD_TYPE;
    unsigned long long max()const{return m - 1;}
    unsigned long long operator()(WORD_TYPE const& x)const{return h(x) % m;}
    typedef typename HASHER::Builder Builder;
    Builder makeBuilder()const{return h.makeBuilder();}
    unsigned long long operator()(Builder b)const{return h(b) % m;}
};
template<typename HASHER> class BHash
{
    unsigned long long mask;
    HASHER h;
public:
    BHash(unsigned long long m): mask(m - 1){assert(m > 0 && isPowerOfTwo(m));}
    typedef typename HASHER::WORD_TYPE WORD_TYPE;
    unsigned long long max()const{return mask;}
    unsigned long long operator()(WORD_TYPE const& x)const{return h(x) & mask;}
    typedef typename HASHER::Builder Builder;
    Builder makeBuilder()const{return h.makeBuilder();}
    unsigned long long operator()(Builder b)const{return h(b) & mask;}
};

class BUHash
{
    uint32_t a, wLB;
    BHash<PrimeHash> h;
public:
    BUHash(unsigned long long m): a(GlobalRNG().next() | 1),//ensure non-0
        wLB(32 - lgCeiling(m)), h(m) {assert(m > 0 && isPowerOfTwo(m));}
    typedef uint32_t WORD_TYPE;
    unsigned long long max()const{return h.max();}
    uint32_t operator()(WORD_TYPE const& x)const{return (a * x) >> wLB;}
    typedef BHash<PrimeHash>::Builder Builder;
    Builder makeBuilder()const{return h.makeBuilder();}
    unsigned long long operator()(Builder b)const{return h(b);}
};

template<typename HASHER> class EHash
{//takes special care to avoid template substitution compile errors
    HASHER h;
    template<typename WORD> unsigned long long hashWord(WORD x, true_type)const
    {//integral type - hash as word if possible
        if(sizeof(WORD) <= sizeof(WORD_TYPE)) return h(x);
        return operator()(&x, 1);//word to big, will break in chunks
    }
public:
    EHash(){}//for h that use no m
    EHash(unsigned long long m): h(m) {}
    typedef typename HASHER::WORD_TYPE WORD_TYPE;
    unsigned long long max()const{return h.max();}
    template<typename WORD> unsigned long long operator()(WORD x)const
        {return hashWord(x, is_integral<WORD>());}//integral words only
    class Builder
    {
        enum{K = sizeof(WORD_TYPE)};
        union
        {
            WORD_TYPE xi;
            unsigned char bytes[K];
        };
        int byteIndex;
        typename HASHER::Builder b;
        friend EHash;
        Builder(EHash const& eh): xi(0), byteIndex(0), b(eh.h.makeBuilder()) {}
        template<typename WORD> void add(WORD const& xi, true_type)
        {//only word type supported for safety
            if(sizeof(WORD) % sizeof(WORD_TYPE) == 0)
                //exact multiple, add as word_type
                for(int i = 0; i < sizeof(WORD)/sizeof(WORD_TYPE); ++i)
                    b.add(((WORD_TYPE*)&xi)[i]);
            else//as char sequence
                for(int i = 0; i < sizeof(xi); ++i)
                    add(((unsigned char*)&xi)[i]);
        }
    public:
        void add(unsigned char bi)
        {
            bytes[byteIndex++] = bi;
            if(byteIndex >= K)
            {
                byteIndex = 0;
                b.add(xi);
                xi = 0;
            }
        }
        template<typename WORD> void add(WORD const& xi)
            {add(xi, is_integral<WORD>());}//integral words only
        typename HASHER::Builder operator()()
        {//finalize remaining xi if any
            if(byteIndex > 0) b.add(xi);
            return b;
        }
    };
    Builder makeBuilder()const{return Builder(*this);}
    unsigned long long operator()(Builder b)const{return h(b());}
    template<typename WORD>
    unsigned long long operator()(WORD* array, int size)const
    {
        Builder b(makeBuilder());
        for(int i = 0; i < size; ++i) b.add(array[i]);
        return operator()(b);
    }
};

class TableHash
{
    enum{N = 1 << numeric_limits<unsigned char>::digits};
    unsigned int table[N];
public:
    TableHash(){for(int i = 0; i < N; ++i) table[i] = GlobalRNG().next();}
    typedef unsigned char WORD_TYPE;
    unsigned long long max()const{return numeric_limits<unsigned int>::max();}
    unsigned int operator()(WORD_TYPE const& x)const
    {
        Builder b(makeBuilder());
        b.add(x);
        return b.sum;
    }
    unsigned int update(unsigned int currentHash, unsigned char byte)
        const{return currentHash ^ table[byte];}//for both add and remove
    class Builder
    {
        unsigned long long sum;
        TableHash const& h;
        friend TableHash;
        Builder(TableHash const& theH): sum(0), h(theH) {}
    public://unlikely possible overflow from add & mult but that's ok
        void add(unsigned char xi){sum ^= h.table[xi];}
    };
    Builder makeBuilder()const{return Builder(*this);}
    unsigned long long operator()(Builder b)const{return b.sum;}
};

struct FNVHash
{
    typedef unsigned char WORD_TYPE;
    unsigned long long max()const{return numeric_limits<uint32_t>::max();}
    uint32_t operator()(WORD_TYPE const& x)const
    {
        Builder b(makeBuilder());
        b.add(x);
        return b.sum;
    }
    class Builder
    {
        uint32_t sum;
        friend FNVHash;
        Builder(): sum(2166136261u) {}
    public://unlikely possible overflow from add & mult but that's ok
        void add(WORD_TYPE xi){sum = (sum * 16777619) ^ xi;}
    };
    Builder makeBuilder()const{return Builder();}
    uint32_t operator()(Builder b)const{return b.sum;}
};
struct FNVHash64
{
    typedef unsigned char WORD_TYPE;
    unsigned long long max()const{return numeric_limits<uint64_t>::max();}
    uint64_t operator()(WORD_TYPE const& x)const
    {
        Builder b(makeBuilder());
        b.add(x);
        return b.sum;
    }
    class Builder
    {
        uint64_t sum;
        friend FNVHash64;
        Builder(): sum(14695981039346656037ull) {}
    public:
        void add(WORD_TYPE xi){sum = (sum * 1099511628211ull) ^ xi;}
    };
    Builder makeBuilder()const{return Builder();}
    uint64_t operator()(Builder b)const{return b.sum;}
};

struct Xorshift64Hash
{
    typedef uint64_t WORD_TYPE;
    unsigned long long max()const{return numeric_limits<WORD_TYPE>::max();}
    uint64_t operator()(WORD_TYPE x)const
        {return QualityXorshift64::transform(x);}
    class Builder
    {
        uint64_t sum;
        friend Xorshift64Hash;
        Builder(): sum(0) {}
    public:
        void add(WORD_TYPE xi){sum = QualityXorshift64::transform(sum + xi);}
    };
    Builder makeBuilder()const{return Builder();}
    uint64_t operator()(Builder b)const{return b.sum;}
};

template<typename HASHER = EHash<BUHash> > class DataHash
{
    HASHER h;
public:
    DataHash(unsigned long long m): h(m){}
    typedef typename HASHER::WORD_TYPE WORD_TYPE;
    unsigned long long max()const{return h.max();}
    unsigned long long operator()(string const& item)const
        {return h(item.c_str(), item.size());}
    unsigned long long operator()(Bitset<> const& item)
        const{return h(item.getStorage().getArray(), item.wordSize());}
    template<typename VECTOR> unsigned long long operator()(VECTOR const& item)
        const{return h(item.getArray(), item.getSize());}
    typedef EMPTY Builder;
};

}//end namespace
#endif
