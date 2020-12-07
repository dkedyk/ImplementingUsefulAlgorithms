#ifndef IGMDK_COMBINATORIAL_GENERATION_H
#define IGMDK_COMBINATORIAL_GENERATION_H
#include "../Sorting/Sort.h"

namespace igmdk{

struct Permutator
{
    Vector<int> p;
    Permutator(int size){for(int i = 0; i < size; ++i) p.append(i);}
    bool next()
    {//find largest i such that p[i] < p[i + 1]
        int j = p.getSize() - 1, i = j - 1;//start with one-before-last
        while(i >= 0 && p[i] >= p[i + 1]) --i;
        bool backToIdentity = i == -1;
        if(!backToIdentity)
        {//find j such that p[j] is next largest element after p[i]
            while(i < j && p[i] >= p[j]) --j;
            swap(p[i], p[j]);
        }
        p.reverse(i + 1, p.getSize() - 1);
        return backToIdentity;//true if returned to smallest permutation
    }
    bool advance(int i)
    {
        assert(i >= 0 && i < p.getSize());
        quickSort(p.getArray(), i + 1, p.getSize() - 1,
            ReverseComparator<int>());
        return next();
    }
};

struct Combinator
{
    int n;
    Vector<int> c;
    Combinator(int m, int theN): n(theN), c(m, -1)
    {
        assert(m <= n && m > 0);
        skipAfter(0);
    }
    void skipAfter(int i)
    {//increment c[i], and reset all c[j] for j > i
        assert(i >= 0 && i < c.getSize());
        ++c[i];
        for(int j = i + 1; j < c.getSize(); ++j) c[j] = c[j - 1] + 1;
    }
    bool next()
    {//find rightmost c[i] which can be increased
        int i = c.getSize() - 1;
        while(i >= 0 && c[i] == n - c.getSize() + i) --i;
        bool finished = i == -1;
        if(!finished) skipAfter(i);
        return finished;
    }
};

struct Partitioner
{
    Vector<int> p;
    Partitioner(int n): p(n, 0) {assert(n > 0);}
    bool skipAfter(int k)
    {//set trailing elements to maximum values and call next
        assert(k >= 0 && k < p.getSize());
        for(int i = k; i < p.getSize(); ++i) p[i] = i;
        return next();
    }
    bool next()
    {//find rightmost p[j] which can be increased
        int m = 0, j = -1;
        for(int i = 0; i < p.getSize(); ++i)
        {
            if(p[i] < m) j = i;
            m = max(m, p[i] + 1);
        }
        bool finished = j == -1;
        if(!finished)
        {//increase it and reset the tail
            ++p[j];
            for(int i = j + 1; i < p.getSize(); ++i) p[i] = 0;
        }
        return finished;
    }
};

}//end namespace
#endif
