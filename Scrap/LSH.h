#ifndef IGMDK_LSH_H
#define IGMDK_LSH_H
#include "../Utils/Utils.h"
#include "../Utils/Debug.h"
#include "../Utils/Vector.h"
#include "../Heaps/Heap.h"
#include "../HashTable/ChainingHashTable.h"
#include "../RandomNumberGeneration/Random.h"
#include "../RandomNumberGeneration/Statistics.h"
#include "../Utils/GCFreelist.h"
#include <cmath>
namespace igmdk{

class E2LSHHasher
{
    Vector<EHash<Xorshift64Hash> > mappers;
    struct Hasher
    {
        Vector<double> a;
        double w, b;
        Hasher(int D, int r): w(r), b(GlobalRNG().uniform01() * w)
            {for(int i = 0; i < D; ++i) a.append(GlobalRNG().normal01());}
        int operator()(Vector<double> const& x)const
            {return 1 + int((dotProduct(a, x) + b)/w);}
    };
    Vector<Hasher> h;
public:
    typedef unsigned long long RESULT_TYPE;
    typedef Vector<double> ITEM_TYPE;
    E2LSHHasher(int k, int l, int D, double w): mappers(l)
        {for(int i = 0; i < k * l; ++i) h.append(Hasher(D, w));}
    RESULT_TYPE operator()(ITEM_TYPE const& x, int bucket)const
    {
        Vector<int> result;
        int k = h.getSize()/mappers.getSize();
        for(int i = 0; i < k; ++i) result.append(h[bucket * k + i](x));
        return mappers[bucket](result.getArray(), result.getSize());
    }
    static double p(double w, double r)
    {
        double z = r/w/1;
        return 2 * approxNormalCDF(z) - 1 -
            2/sqrt(2 * PI()) * (1 - exp(-z * z/2));
    }
    static double p1(double r){return p(r, r);}
    static double p2(double r, double c){return p(r, r * c);}
    static double distance(ITEM_TYPE const& x1, ITEM_TYPE const& x2)
    {
        EuclideanDistance<ITEM_TYPE>::Distance ed;
        return ed(x1, x2);
    }
};

namespace LSHKLFinder
{
    int LSHGetL(int k, double p1, double e)
    {
        double l = log(e)/log(1 - pow(p1, k));
        return (!isfinite(l) && l > numeric_limits<int>::max()) ? -1 : 1 + int(l);
    }
    double LSHCost(int k, double e, double p1, double p2, int n)
    {
        int l = LSHGetL(k, p1, e);
        return (10 * k + pow(p2, k) * n) * l;
    }
    int minL(double p1, double e){return LSHGetL(1, p1, e);};
    int LSHFindK(double e, double p1, double p2, int n, int maxL)
    {
        int bestK = -1;
        double bestV;
        for(int k = 1;; ++k)
        {
            DEBUG(k);
            int l = LSHGetL(k, p1, e);
            DEBUG(l);
            double v = LSHCost(k, e, p1, p2, n);
            DEBUG(v);
            if(v < 0) break;

            if(bestK == -1 || (l > 0 && l < maxL && v < bestV)) {bestK = k; bestV = v;}
            if(l < 0 || l >= maxL) break;
        }
        DEBUG(bestV);
        //DEBUG(LSHGetL(bestK, p1, e));
        return bestK;
    }
}

template<typename HASHER> class LSH
{
    typedef typename HASHER::ITEM_TYPE ITEM;
    typedef typename HASHER::RESULT_TYPE RESULT_TYPE;
    Vector<ChainingHashTable<RESULT_TYPE, Vector<int> > > buckets;
    Vector<ITEM> items;
    HASHER g;
    double r2;
public:
    LSH(HASHER const& theG, int l, double theR2): buckets(l, l), g(theG), r2(theR2){}
    void insert(ITEM const& x)
    {
        for(int i = 0; i < buckets.getSize(); ++i)
        {
            typename HASHER::RESULT_TYPE hash = g(x, i);
            //DEBUG(i);
            //DEBUG(hash);
            Vector<int>* xBucket = buckets[i].find(hash);
            if(!xBucket)
            {
                buckets[i].insert(hash, Vector<int>());
                xBucket = buckets[i].find(hash);
            }
            xBucket->append(items.getSize());//have linear probing return chain instead?
        }
        items.append(x);
    }
    Vector<ITEM> cNeighbors(ITEM const& x)
    {
        Vector<ITEM> result;
        ChainingHashTable<int, bool> retrievedItems;
        int hitItems = 0;
        for(int i = 0; i < buckets.getSize(); ++i)
        {
            typename HASHER::RESULT_TYPE hash = g(x, i);
            //DEBUG(i);
            //DEBUG(hash);
            Vector<int>* xBucket = buckets[i].find(hash);
            if(xBucket)
                for(int i = 0; i < xBucket->getSize(); ++i)
                {
                    int itemIndex = (*xBucket)[i];
                    ++hitItems;
                    if(!retrievedItems.find(itemIndex))
                    {
                        retrievedItems.insert(itemIndex, true);
                        if(HASHER::distance(x, items[itemIndex]) < r2)
                            result.append(items[itemIndex]);
                    }
                }
        }
        //DEBUG(hitItems);
        return result;
    }
};

LSH<E2LSHHasher> buildE2LSH(int D, double r, double c, int maxL, double e = 10e-6, int maxN = 1000000)
{
    double p1 = E2LSHHasher::p(1, 1), r2 = r * (1 + c);
    int k = LSHKLFinder::LSHFindK(e, p1, E2LSHHasher::p(r, r2), maxN, maxL);
    //DEBUG(k);
    int l = LSHKLFinder::LSHGetL(k, p1, e);
    //DEBUG(l);
    return LSH<E2LSHHasher>(E2LSHHasher(k, l, D, r), l, r2);
}

template<typename HASHER> class NearestNeighborLSH
{
    typedef typename HASHER::ITEM_TYPE ITEM;
    Vector<LSH<HASHER> > lshs;//items are duplicated dont store them!
public:
    void addLSH(LSH<HASHER> const& lsh){lshs.append(lsh);}
    void insert(ITEM const& x)
        {for(int i = 0; i < lshs.getSize(); ++i) lshs[i].insert(x);}
    pair<ITEM, bool> cNeighbor(ITEM const& x)
    {
        for(int i = 0; i < lshs.getSize(); ++i)
        {
            Vector<ITEM> items = lshs[i].cNeighbors(x);
            if(items.getSize() > 0)
            {
                int best = -1, bestD;
                for(int j = 0; j < items.getSize(); ++j)
                {
                    double d = HASHER::distance(x, items[j]);
                    if(best == -1 || d < bestD)
                    {
                        best = j;
                        bestD = d;
                    }
                }
                return pair<ITEM, bool>(items[best], true);
            }
        }
        return pair<ITEM, bool>(ITEM(), false);
    }
};

NearestNeighborLSH<E2LSHHasher> buildE2NNLSH(int D, double rMin, double rMax, int maxL, double c = 1, double e = 10e-6, int maxN = 1000000)
{
    NearestNeighborLSH<E2LSHHasher> result;
    for(double r = rMin; r < rMax; r *= (1 + c))
    {
        result.addLSH(buildE2LSH(D, r, c, maxL, e, maxN));
    }
    return result;
}

}//end namespace
#endif
