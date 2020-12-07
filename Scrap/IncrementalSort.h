#ifndef IGMDK_INCREMENTAL_SORT_H
#define IGMDK_INCREMENTAL_SORT_H
#include "../Utils/Stack.h"
#include "../Sorting/Sort.h"
namespace igmdk{

template<typename ITEM, typename COMPARATOR> ITEM incrementalQuickSelect(
    ITEM* vector, int left, Stack<int>& s, COMPARATOR const& c)
{
    for(int right, i, j; left < (right = s.getTop()); s.push(j))
        partition3(vector, left, right, i, j, c);
    s.pop();
    return vector[left];
}
template<typename ITEM, typename COMPARATOR> void incrementalSort(
    ITEM* vector, int n, COMPARATOR const& c)
{
    Stack<int> s;
    s.push(n - 1);
    for(int i = 0; i < n; ++i) incrementalQuickSelect(vector, i, s, c);
}

}//end namespace
#endif
