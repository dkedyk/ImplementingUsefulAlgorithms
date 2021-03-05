#ifndef IGMDK_SORT_TEST_AUTO_H
#define IGMDK_SORT_TEST_AUTO_H

using namespace std;
#include "../Utils/Vector.h"
#include "Sort.h"

namespace igmdk{

Vector<int> makeRandVector(int n = 1000000)
{
    Vector<int> v(n);
    for(int i = 0; i < n; ++i) v[i] = GlobalRNG().mod(n);
    return v;
}

Vector<int> makeNegatedVector(int n = 1000000)
{
    Vector<int> v(n);
    for(int i = 0; i < n; ++i) v[i] = -i;
    return v;
}

void verifySortedNegatedVector(Vector<int> const& v)
{
    int n = v.getSize();
    for(int i = 0; i < n; ++i) assert(v[i] == i + 1 - n);
}

void testQuicksortAuto()
{
    DEBUG("testQuicksortAuto");
    Vector<int> nums = makeRandVector();
    quickSort(nums.getArray(), nums.getSize());
    assert(isSorted(nums.getArray(), 0, nums.getSize() - 1,
        DefaultComparator<int>()));
    nums = makeNegatedVector();
    quickSort(nums.getArray(), nums.getSize());
    verifySortedNegatedVector(nums);
    DEBUG("testQuicksortAuto passed");
}

void testMergesortAuto()
{
    DEBUG("testMergesortAuto");
    Vector<int> nums = makeRandVector();
    mergeSort(nums.getArray(), nums.getSize());
    assert(isSorted(nums.getArray(), 0, nums.getSize() - 1,
        DefaultComparator<int>()));
    nums = makeNegatedVector();
    mergeSort(nums.getArray(), nums.getSize());
    verifySortedNegatedVector(nums);
    DEBUG("testMergesortAuto passed");
}

void testPermutationSortAuto()
{
    DEBUG("testPermutationSortAuto");
    char v[4] = {'a', 'b', 'c', 'd'};
    int p[4] = {3, 2, 1, 0};
    permutationSort(v, p, 4);
    assert(v[0] == 'd');
    assert(v[1] == 'c');
    assert(v[2] == 'b');
    assert(v[3] == 'a');
    DEBUG("testPermutationSortAuto passed");
}

void testCountingsortAuto()
{
    DEBUG("testCountingsortAuto");
    Vector<int> nums = makeRandVector();
    countingSort(nums.getArray(), nums.getSize(), nums.getSize());
    assert(isSorted(nums.getArray(), 0, nums.getSize() - 1,
        DefaultComparator<int>()));
    DEBUG("testCountingsortAuto passed");
}

void testKsortAuto()
{
    DEBUG("testKsortAuto");
    Vector<int> nums = makeRandVector();
    KSort(nums.getArray(), nums.getSize(), nums.getSize(),
        IdentityHash<int>());
    assert(isSorted(nums.getArray(), 0, nums.getSize() - 1,
        DefaultComparator<int>()));
    DEBUG("testKsortAuto passed");
}

Vector<int> makeRandPermVector(int n = 1000000)
{
    Vector<int> v(n);
    for(int i = 0; i < n; ++i) v[i] = i;
    GlobalRNG().randomPermutation(v.getArray(), n);
    return v;
}
void testQuickselect()
{
    DEBUG("testQuickselectAuto");
    Vector<int> nums = makeRandPermVector(),
        goal = GlobalRNG().randomCombination(100, nums.getSize());
    for(int i = 0; i < goal.getSize(); ++i)
    {
        quickSelect(nums.getArray(), nums.getSize(), goal[i]);
        assert(nums[goal[i]] == goal[i]);
    }
    countingSort(nums.getArray(), nums.getSize(), nums.getSize());
    DEBUG("testQuickselectAuto passed");
}
void testMultipleQuickselect()
{
    DEBUG("testMultipleQuickselectAuto");
    Vector<int> nums = makeRandPermVector(),
        goal = GlobalRNG().randomCombination(100, nums.getSize());
    Vector<bool> selected(nums.getSize(), false);
    for(int i = 0; i < goal.getSize(); ++i) selected[goal[i]] = true;
    multipleQuickSelect(nums.getArray(), selected.getArray(), 0,
        nums.getSize() - 1, DefaultComparator<int>());
    for(int i = 0; i < goal.getSize(); ++i) assert(nums[goal[i]] == goal[i]);
    DEBUG("testMultipleQuickselectAuto passed");
}

void testMultikeyAuto()
{
    DEBUG("testMultikeyAuto");
    Vector<string> s;
    s.append("fsdlfjl");
    s.append("wejk");
    s.append("iosufrwhrew");
    s.append("wqjklhdsaiohd");
    s.append("wioeurksd");
    s.append("");
    s.append("w");
    Vector<int> perm;
    for(int i = 0; i < s.getSize(); ++i) perm.append(i);
    typedef IndexedAccessor<string, char, StringAccessor> T;
    MultikeyQuicksortVectorComparator<int, char, T> c(T(s.getArray()));
    int k = 3;
    multikeyQuickselect(perm.getArray(), 0, s.getSize() - 1, k, c);
    int kth = perm[k];
    multikeyQuicksortNR(perm.getArray(), 0, s.getSize() - 1, c);
    assert(isSorted(perm.getArray(), 0, s.getSize() - 1,
        IndexComparator<string>(s.getArray())));
    assert(kth == perm[k]);
    DEBUG("testMultikeyAuto passed");
}

void testAllAutoSort()
{
    DEBUG("testAllAutoSort");
    testQuicksortAuto();
    testMergesortAuto();
    testPermutationSortAuto();
    testCountingsortAuto();
    testKsortAuto();
    testQuickselect();
    testMultipleQuickselect();
    testMultikeyAuto();
}

}//end namespace
#endif
