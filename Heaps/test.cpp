#include "HeapTestAuto.h"
#include <cassert>
#include <iostream>
#include <cstdlib>
#include <ctime>
using namespace igmdk;

void timeSRT()
{
	IndexedHeap<int> heap;
	int N = 1500000;
	//IndexedArrayHeap<int> heap;
	for(int i = 0; i < N; ++i)
	{
		heap.insert(rand()%10, i);
	}
	for(int i = 0; i < N; ++i)
	{
		heap.deleteMin();
	}
}

void DDDIndexedHeap()
{
    IndexedHeap<int> IndexedHeap0to3;
    for(int i = 0; i < 4; ++i)
	{
		IndexedHeap0to3.insert(rand(), i);
	}
	cout << "breakpoint" << endl;
}

void DDDIndexedArrayHeap()
{
    IndexedArrayHeap<int> IndexedArrayHeap0to3;
    for(int i = 0; i < 4; ++i)
	{
		IndexedArrayHeap0to3.insert(rand(), i);
	}
	cout << "breakpoint" << endl;
}

int main()
{
    testAllAutoHeaps();
    DDDIndexedHeap();
    DDDIndexedArrayHeap();
	timeSRT();
	return 0;
}
