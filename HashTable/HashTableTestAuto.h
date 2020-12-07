#ifndef IGMDK_HASH_TABLE_TEST_AUTO_H
#define IGMDK_HASH_TABLE_TEST_AUTO_H

#include "MapTestAutoHelper.h"
#include "ChainingHashTable.h"
#include "LinearProbingHashTable.h"

namespace igmdk{

void testChainingHashTableAuto()
{
    DEBUG("testChainingHashTableAuto");
    testMapAutoHelper<ChainingHashTable<int, int> >();
    DEBUG("testChainingHashTableAuto passed");
}

void testLinearProbingHashTableAuto()
{
    DEBUG("testLinearProbingHashTableAuto");
    testMapAutoHelper<LinearProbingHashTable<int, int> >();
    DEBUG("testLinearProbingHashTableAuto passed");
}

void testAllAutoHashTable()
{
    DEBUG("testAllAutoHashTable");
    testChainingHashTableAuto();
    testLinearProbingHashTableAuto();
}

}//end namespace
#endif
