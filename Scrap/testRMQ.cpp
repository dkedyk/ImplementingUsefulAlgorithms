#include <cassert>
#include <iostream>
#include "../Utils/Debug.h"
#include "RMQ.h"
using namespace std;
using namespace igmdk;

void testLRMTree()
{
    int array[] = {15,8,13,7,11,16,1,10,9,14,2,12,3,6,5,4},
        size = sizeof(array)/sizeof(*array);
    LRMTree<int> lrm(array, size);
    for(int i = 0; i < size; ++i) DEBUG(lrm.parents[i]);
    DEBUG(lrm.RMQ(4,8));
    DEBUG(lrm.RMQ(0,5));
    DEBUG(lrm.RMQ(0,size-1));
    DEBUG(lrm.RMQ(11,size-1));
}

int main()
{
    testLRMTree();
    return 0;
}
