#ifndef IGMDK_RMQ_H
#define IGMDK_RMQ_H

#include "../Utils/Vector.h"
namespace igmdk{

template<typename ITEM> int LCA(ITEM* array, int i, int j)
{//O(h) performance
    int hDiff = 0;
    for(int k = i; array[k] != -1; k = array[k]) --hDiff;
    for(int k = j; array[k] != -1; k = array[k]) ++hDiff;
    if (hDiff < 0) {swap(i, j); hDiff = -hDiff;}
    while(hDiff--) j = array[j];
    while(i != j) {i = array[i]; j = array[j];}
    return i;
}

template<typename ITEM> struct LRMTree
{//previous smaller value tree
    Vector<int> parents;
    LRMTree(ITEM* array, int size): parents(size, -1)
    {
        for(int i = 1; i < size; ++i)
        {//expand the rightmost branch unless a smaller value causes another
            //rightmost branch to be created
            parents[i] = i - 1;
            while(parents[i] != -1 && array[parents[i]] >= array[i])
                parents[i] = parents[parents[i]];
        }
    }
    int RMQ(int left, int right)
    {
        int l = LCA(parents.getArray(), left, right);
        if(l == left) return l;
        while(parents[right] != l) right = parents[right];
        return right;
    }
};

}
#endif
