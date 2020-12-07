#ifndef IGMDK_EXTERNAL_MEMORY_ALGORITHMS_TEST_AUTO_H
#define IGMDK_EXTERNAL_MEMORY_ALGORITHMS_TEST_AUTO_H

#include "File.h"
#include "EMVector.h"
#include "EMBTree.h"
#include "../Utils/Utils.h"

namespace igmdk{

void testFileBasicsAuto()
{
    DEBUG("testFileBasicsAuto");
    string filename = "m1_13.igmdk";
    File::remove(filename.c_str());
    {
        File file(filename.c_str(), true);
        int m1 = -1;
        file.append((unsigned char*)&m1, sizeof(m1));
        double oneover3 = 1.0/3;//not exactly representable
        file.append((unsigned char*)&oneover3, sizeof(oneover3));
        assert(file.getSize() == sizeof(m1) + sizeof(oneover3));
        file.setPosition(0);
        file.read((unsigned char*)&m1, sizeof(m1));
        assert(m1 == -1);
        file.read((unsigned char*)&oneover3, sizeof(oneover3));
        assert(abs(oneover3 - 1.0/3)/(1.0/3) <//check rel error
            numeric_limits<double>::epsilon());
    }
    File::remove(filename.c_str());
    assert(!File::exists(filename.c_str()));
    DEBUG("testFileBasicsAuto passed");
}

void testVectorAuto()
{
    DEBUG("testVectorAuto");
    File::remove("EMVector.igmdk");
    int n = 100000;
    {//this pass write
        EMVector<int> v("EMVector.igmdk");
        for(int i = 0; i < n; ++i)
        {
            v.append(-1);
        }
    }//force destructor
    {//this pass read random access
        EMVector<int> v("EMVector.igmdk");
        int n = v.getSize();
        int sum = 0;
        for(int i = 0; i < n; ++i) sum += v[GlobalRNG().mod(n)];
        assert(sum == -n);
    }
    File::remove("EMVector.igmdk");
    DEBUG("testVectorAuto passed");
}

void testIOSortAuto()
{
    DEBUG("testIOSortAuto");
    File::remove("111.igmdk");
    {
        EMVector<int> vec("111.igmdk");
        int K = 100000;
        for(int i = 0; i < K; ++i)
        {
            vec.append(-i);
        }
        IOSort(vec);
        assert(vec.getSize() == K);
        for(int i = 1; i < vec.getSize(); ++i)
        {
            assert(vec[i] == i + 1 - K);
            assert(vec[i - 1] <= vec[i]);
        }
    }
    File::remove("111.igmdk");
    DEBUG("testIOSortAuto passed");
}

void testBTreeAuto()
{
    DEBUG("testBTreeAuto");
    File::remove("HeaderBPlusTree.igmdk");
    File::remove("KeysBPlusTree.igmdk");
    File::remove("NodesRecordsBPlusTree.igmdk");
    File::remove("ReturnedRecordsBPlusTree.igmdk");
    int N = 10000;
    {
        EMBPlusTree<int, int> trie("BPlusTree.igmdk");
        for(int i = 0; i < N; ++i)
        {
            trie.insert(-i, -i);
        }
    }//force destructor
    {
        EMBPlusTree<int, int> trie("BPlusTree.igmdk");
        for(int i = 0; i < N; ++i)
        {
            bool status;
            int item = trie.find(-i, status);
            assert(status);
            assert(item == -i);
            trie.remove(-i);
            trie.find(-i, status);
            assert(!status);
        }
    }
    File::remove("HeaderBPlusTree.igmdk");
    File::remove("KeysBPlusTree.igmdk");
    File::remove("NodesRecordsBPlusTree.igmdk");
    File::remove("ReturnedRecordsBPlusTree.igmdk");
    DEBUG("testBTreeAuto passed");
}

void testAllAutoExternalMemoryAlgorithms()
{
    DEBUG("testAllAutoExternalMemoryAlgorithms");
    testVectorAuto();
    testIOSortAuto();
    testBTreeAuto();
}

}//end namespace
#endif
