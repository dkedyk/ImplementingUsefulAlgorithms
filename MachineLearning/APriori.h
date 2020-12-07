#ifndef IGMDK_APRIORI_H
#define IGMDK_APRIORI_H

#include "../Sorting/Sort.h"
#include "../RandomTreap/LCPTreap.h"
#include "../MiscAlgs/CombinatorialGeneration.h"
namespace igmdk{

struct APriori
{
    LCPTreap<Vector<int>, int> counts;
    int processBasket(Vector<int> const& basket, int round,
        int rPrevMinCount = 0, int r1MinCount = 0)
    {
        int addedCount = 0;
        if(basket.getSize() > round)
        {
            Combinator c(round, basket.getSize());
            do//prepare the current combination of ids, needn't sort if each
            {//basket is already sorted
                Vector<int> key, single;
                for(int i = 0; i < round; ++i) key.append(basket[c.c[i]]);
                quickSort(key.getArray(), key.getSize());
                int* count = counts.find(key);
                if(count) ++*count;//combination is frequent if already
                else if(round == 1)//frequent or round is 1
                {
                    counts.insert(key, 1);
                    ++addedCount;
                }
                else//combination is frequent if the last item and
                {//combination without the last item are both frequent
                    single.append(key.lastItem());
                    if(*counts.find(single) >= r1MinCount)
                    {
                        key.removeLast();
                        if(*counts.find(key) >= rPrevMinCount)
                        {
                            key.append(single[0]);
                            counts.insert(key, 1);
                            ++addedCount;
                        }
                    }
                }
            }while(!c.next());
        }
        return addedCount;
    }
    void noCutProcess(Vector<Vector<int> >const& baskets, int nRounds)
    {
        for(int k = 1; k <= nRounds; ++k)
            for(int i = 0; i < baskets.getSize(); ++i)
                processBasket(baskets[i], k);
    }
};

}//end namespace
#endif

