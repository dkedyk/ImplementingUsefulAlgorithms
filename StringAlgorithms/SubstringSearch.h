#ifndef IGMDK_SUBSTRING_SEARCH_H
#define IGMDK_SUBSTRING_SEARCH_H
#include "../Utils/Utils.h"
#include "../Utils/Vector.h"
#include "../Utils/Bits.h"
namespace igmdk{

template<typename VECTOR, typename VECTOR2> bool matchesAt(int position,
    VECTOR2 text, VECTOR pattern, int patternSize)
{//allows different text and pattern types
    int i = 0;
    while(i < patternSize && pattern[i] == text[i + position]) ++i;
    return i == patternSize;
}

struct LecroqHash//identity for q = 1
{//ignore the size parameter - the matcher will have enough table size
    LecroqHash(int dummy){}
    struct Builder
    {
        unsigned char result;
        Builder(): result(0){}
        void add(unsigned char c){result = result << 1 + c;}
    };
    Builder makeBuilder(){return Builder();}
    unsigned char operator()(Builder b){return b.result;}
};
template<typename VECTOR, typename HASHER = LecroqHash>
class HashQ
{
    enum{CHAR_ALPHABET_SIZE = 1 << numeric_limits<unsigned char>::digits};
    int patternSize, q;
    Vector<int> shifts;//size is a power of 2 for fast hashing
    VECTOR const &pattern;
    HASHER h;
    typedef typename HASHER::Builder B;
public:
    HashQ(VECTOR const& thePattern, int thePatternSize, int theQ = 1): q(theQ),
        pattern(thePattern), patternSize(thePatternSize), shifts(max<int>(
        CHAR_ALPHABET_SIZE, nextPowerOfTwo(ceiling(patternSize, q)))),
        h(shifts.getSize())
    {//precompute shifts
        assert(patternSize >= q);
        int temp = patternSize - q;
        for(int i = 0; i < shifts.getSize(); ++i) shifts[i] = temp + 1;
        for(int i = 0; i < temp; ++i)
        {
            B b(h.makeBuilder());
            for(int j = 0; j < q; ++j) b.add(pattern[i + j]);
            shifts[h(b)] = temp - i;
        }
    }//return match position and next start position of (-1, -1)
    template<typename VECTOR2> pair<int, int> findNext(VECTOR2 const& text,
        int textSize, int start = 0)//allow different text and pattern types
    {
        while(start + patternSize <= textSize)
        {
            int result = start, hStart = start + patternSize - q;
            B b(h.makeBuilder());
            for(int j = 0; j < q; ++j) b.add(text[hStart + j]);
            start += shifts[h(b)];
            if(matchesAt(result, text, pattern, patternSize))
                return make_pair(result, start);
        }
        return make_pair(-1, -1);
    }
};

template<typename VECTOR, typename HASHER = LecroqHash> class WuManber
{
    enum{CHAR_ALPHABET_SIZE = 1 << numeric_limits<unsigned char>::digits};
    int q, minPatternSize;
    Vector<pair<VECTOR, int> > const& patterns;
    Vector<int> shifts;//size is a power of 2 for fast hashing
    Vector<Vector<int> > candidates;
    HASHER h;
    typedef typename HASHER::Builder B;
public:
    WuManber(Vector<pair<VECTOR, int> > const& thePatterns, int theQ = 1,
        double avePatternSize = 1): q(theQ), patterns(thePatterns), shifts(
        max<int>(CHAR_ALPHABET_SIZE, nextPowerOfTwo(avePatternSize *
        patterns.getSize()/q))), candidates(shifts.getSize()),
        h(shifts.getSize()), minPatternSize(numeric_limits<int>::max())
    {//precompute shifts
        for(int i = 0; i < patterns.getSize(); ++i)
            minPatternSize = min(patterns[i].second, minPatternSize);
        assert(patterns.getSize() > 0 && minPatternSize >= q);
        int temp = minPatternSize - q;
        for(int i = 0; i < shifts.getSize(); ++i) shifts[i] = temp + 1;
        for(int j = 0; j < patterns.getSize(); ++j)
            for(int i = 0; i < temp + 1; ++i)
            {
                B b(h.makeBuilder());
                for(int k = 0; k < q; ++k) b.add(patterns[j].first[i + k]);
                int hi = h(b);
                if(i == temp) candidates[hi].append(j);
                else shifts[hi] = min(temp - i, shifts[hi]);
            }
    }//return match position and next start position of (-1, -1) and indices
    //of patterns that match; allow different text and pattern types
    template<typename VECTOR2> pair<Vector<int>, int> findNext(
        VECTOR2 const& text, int textSize, int start = 0)
    {
        Vector<int> matches(patterns.getSize(), -1);
        while(start + minPatternSize <= textSize)
        {
            B b(h.makeBuilder());
            for(int j = 0; j < q; ++j)
                b.add(text[start + minPatternSize - q + j]);
            int hValue = h(b);
            bool foundAMatch = false;
            for(int i = 0; i < candidates[hValue].getSize(); ++i)
            {
                int j = candidates[hValue][i];
                if(start + patterns[j].second <= textSize && matchesAt(
                    start, text, patterns[j].first, patterns[j].second))
                {
                    foundAMatch = true;
                    matches[j] = start;
                }
            }
            start += shifts[hValue];
            if(foundAMatch) return make_pair(matches, start);
        }
        return make_pair(matches, -1);
    }
};

}//end namespace
#endif

