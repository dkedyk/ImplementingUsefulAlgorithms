#include "IncrementalSort.h"
#include <iostream>
#include <cstdlib>
using namespace igmdk;

void testIncremental()
{
    int simple[] = {5, 1, 3 ,0,4,2};
    incrementalSort(simple, 6, DefaultComparator<int>());
    for(int i = 0; i < 6; ++i)
    {
        DEBUG(simple[i]);
    }
}

int main()
{
	testIncremental();
	return 0;
}
