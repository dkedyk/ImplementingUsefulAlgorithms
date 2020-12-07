#ifndef IGMDK_SUFFIX_ARRAY_H
#define IGMDK_SUFFIX_ARRAY_H
#include "../Sorting/Sort.h"
#include "../Utils/Vector.h"
namespace igmdk{

struct SARank
{
    int* ranks;
    int n, k;
    int operator()(int i)const{i += k; return i < n ? ranks[i] + 1 : 0;}
};
struct BWTRank
{
    int* ranks;
    int n, k;
    int operator()(int i)const{return ranks[(i + k) % n];}
};
template<typename RANKER, typename ITEM>
Vector<int> suffixArray(ITEM* const vector, int n)
{
    Vector<int> ranks(n, 0), sa(n, 0);
    for(int i = 0; i < n; ++i) sa[i] = i;
    quickSort(sa.getArray(), 0, n - 1, IndexComparator<ITEM>(vector));
    ranks[sa[0]] = 0;//set ranks based on first char
    for(int i = 1, r = 0; i < n; ++i)
    {
        if(vector[sa[i]] != vector[sa[i - 1]]) ++r;
        ranks[sa[i]] = r;
    }
    for(int k = 1; k < n; k *= 2)
    {//double
        RANKER r1 = {ranks.getArray(), n, k}, r2 = {ranks.getArray(), n, 0};
        KSort(sa.getArray(), n, n + 1, r1);
        KSort(sa.getArray(), n, n + 1, r2);
        if(k * 2 < n)//else ranks already unique after sort
        {//set doubled ranks based on the tuples
            Vector<int> ranks2(n, 0);
            ranks2[sa[0]] = 0;
            for(int i = 1, r = 0; i < n; ++i)
            {//advance rank if either tuple different
                if(r1(sa[i]) != r1(sa[i - 1]) || r2(sa[i]) != r2(sa[i - 1]))
                    ++r;
                ranks2[sa[i]] = r;
                if(r == n - 1) return sa;//ranks already unique
            }
            ranks.swapWith(ranks2);//more efficient than assignment
        }
    }
    return sa;
}

template<typename ITEM> Vector<int> LCPArray(ITEM* text, int size, int* sa)
{
    Vector<int> pred(size, 0), PLCP(size, 0);
    for(int i = 0; i < size; ++i) pred[sa[i]] = sa[(i ? i : size) - 1];
    for(int i = 0, p = 0; i < size; ++i)
    {
        while(text[i + p] == text[pred[i] + p]) ++p;
        PLCP[i] = p;
        p = max(p - 1, 0);
    }//pred becomes the LCP array now, 0 has lcp(sa[0], sa[n - 1])
    for(int i = 0; i < size; ++i) pred[i] = PLCP[sa[i]];
    return pred;
}

template<typename ITEM> struct SuffixIndex
{
    Vector<ITEM> const& text;
    Vector<int> sa, lcpa;
    SuffixIndex(Vector<ITEM> const& theText): text(theText),
        sa(suffixArray<SARank>(text.getArray(), text.getSize())),
        lcpa(LCPArray(text.getArray(), text.getSize(), sa.getArray())) {}
    bool isKLess(ITEM const* a, int aSize, ITEM const* b, int bSize, int k)
        const
    {
        for(int i = 0; i < min(aSize, bSize); ++i)
            if(a[i] > b[i]) return false;
        return max(aSize, bSize) < k;
    }
    pair<int, int> interval(ITEM* pattern, int size)
    {
        int left = 0, right = sa.getSize();
        while(left < right)
        {
            int i = (left + right)/2;
            if(isKLess(&text[sa[i]], sa.getSize() - sa[i], pattern, size,
                size)) left = i + 1;
            else right = i - 1;
        }
        int left2 = left - 1, right2 = sa.getSize() - 1;
        while(left2 < right2)
        {
            int i = (left2 + right2)/2;
            if(isKLess(pattern, size, &text[sa[i]], sa.getSize() - sa[i],
                size)) right2 = i - 1;
            else left2 = i + 1;
        }
        return make_pair(left, right2);
    }
};

}//end namespace
#endif
