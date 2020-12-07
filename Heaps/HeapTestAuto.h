#ifndef IGMDK_HEAP_TEST_AUTO_H
#define IGMDK_HEAP_TEST_AUTO_H
#include <string>
using namespace std;
#include "Heap.h"
#include "IndexedHeap.h"
#include "../Sorting/Sort.h"

namespace igmdk{

void testHeapAuto()
{
    DEBUG("testHeapAuto");
    int N = 1000000;
    Heap<int> heap;
    for(int i = 0; i < N; ++i)
	{
		heap.insert(i);
	}
    Vector<int> nums;
	for(int i = 0; i < N; ++i)
	{
		nums.append(heap.deleteMin());
	}
    assert(isSorted(nums.getArray(), 0, N - 1, DefaultComparator<int>()));
    DEBUG("testHeapAuto passed");
}

void testAllAutoHeaps()
{
    DEBUG("testAllAutoHeaps");
    testHeapAuto();
}

}//end namespace
#endif
