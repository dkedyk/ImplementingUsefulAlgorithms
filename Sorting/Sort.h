#ifndef IGMDK_SORT_H
#define IGMDK_SORT_H
#include "../Utils/Utils.h"
#include "../Utils/Stack.h"
#include "../RandomNumberGeneration/Random.h"
namespace igmdk{

template<typename ITEM, typename COMPARATOR>
void insertionSort(ITEM* vector, int left, int right, COMPARATOR const& c)
{//allow more general left != 0
    for(int i = left + 1; i <= right; ++i)
    {
        ITEM e = vector[i];
        int j = i;
        for(;j > left && c(e, vector[j - 1]); --j) vector[j] = vector[j - 1];
        vector[j] = e;
    }
}

template<typename ITEM, typename COMPARATOR>
int pickPivot(ITEM* vector, int left, int right, COMPARATOR const& c)
{//OK if pivots are same by chance occasionally
    int i = GlobalRNG().inRange(left, right), j =
        GlobalRNG().inRange(left, right), k = GlobalRNG().inRange(left, right);
    if(c(vector[j], vector[i])) swap(i, j);
    //i <= j, decide where k goes
    return c(vector[k], vector[i]) ? i : c(vector[k], vector[j]) ? k : j;
}

template<typename ITEM, typename COMPARATOR> void partition3(ITEM* vector,
    int left, int right, int& i, int& j, COMPARATOR const& c)
{//i/j are the current left/right pointers
    ITEM p = vector[pickPivot(vector, left, right, c)];
    int lastLeftEqual = i = left - 1, firstRightEqual = j = right + 1;
    for(;;)//the pivot is the sentinel for the first pass
    {//after one swap swapped items act as sentinels
        while(c(vector[++i], p));
        while(c(p, vector[--j]));
        if(i >= j) break;//pointers crossed
        swap(vector[i], vector[j]);//both pointers found swappable items
        //swap equal items to the sides
        if(c.isEqual(vector[i], p))//i to the left
            swap(vector[++lastLeftEqual], vector[i]);
        if(c.isEqual(vector[j], p))//j to the right
            swap(vector[--firstRightEqual], vector[j]);
    }
    //invariant: i == j if they stop at an item = pivot
    //and this can happen at both left and right item
    //or they cross over and i = j + 1
    if(i == j){++i; --j;}//don't touch pivot in the middle
    //swap side items to the middle; left with "<" section and right with ">"
    for(int k = left; k <= lastLeftEqual; ++k) swap(vector[k], vector[j--]);
    for(int k = right; k >= firstRightEqual; --k)
        swap(vector[k], vector[i++]);
}

template<typename ITEM, typename COMPARATOR>
void quickSort(ITEM* vector, int left, int right, COMPARATOR const& c)
{
    while(right - left > 16)
    {
        int i, j;
        partition3(vector, left, right, i, j, c);
        if(j - left < right - i)//pick smaller
        {
            quickSort(vector, left, j, c);
            left = i;
        }
        else
        {
            quickSort(vector, i, right, c);
            right = j;
        }
    }
    insertionSort(vector, left, right, c);
}
template<typename ITEM> void quickSort(ITEM* vector, int n)
    {quickSort(vector, 0, n - 1, DefaultComparator<ITEM>());}

template<typename ITEM, typename COMPARATOR> ITEM quickSelect(ITEM* vector,
    int left, int right, int k, COMPARATOR const& c)
{
    assert(k >= left && k <= right);
    for(int i, j; left < right;)
    {
        partition3(vector, left, right, i, j, c);
        if(k >= i) left = i;
        else if(k <= j) right = j;
        else break;
    }
    return vector[k];
}
template<typename ITEM> ITEM quickSelect(ITEM* vector, int size, int k)
    {return quickSelect(vector, 0, size - 1, k, DefaultComparator<ITEM>());}

template<typename ITEM, typename COMPARATOR> void multipleQuickSelect(ITEM*
    vector, bool* selected, int left, int right, COMPARATOR const& c)
{
    while(right - left > 16)
    {
        int i, j;
        for(i = left; i <= right && !selected[i]; ++i);
        if(i == right+1) return;//none are selected
        partition3(vector, left, right, i, j, c);
        if(j - left < right - i)//smaller first
        {
            multipleQuickSelect(vector, selected, left, j, c);
            left = i;
        }
        else
        {
            multipleQuickSelect(vector, selected, i, right, c);
            right = j;
        }
    }
    insertionSort(vector, left, right, c);
}

template<typename ITEM, typename COMPARATOR> void merge(ITEM* vector,
    int left, int middle, int right, COMPARATOR const& c, ITEM* storage)
{//i for the left half, j for the right, merge until fill up vector
    for(int i = left, j = middle + 1; left <= right; ++left)
    {//either i or j can get out of bounds
        bool useRight = i > middle || (j <= right &&
            c(storage[j], storage[i]));
        vector[left] = storage[(useRight ? j : i)++];
    }
}
template<typename ITEM, typename COMPARATOR> void mergeSortHelper(
    ITEM* vector, int left, int right, COMPARATOR const& c, ITEM* storage)
{
    if(right - left > 16)
    {//sort storage using vector as storage, then merge into vector
        int middle = (right + left)/2;
        mergeSortHelper(storage, left, middle, c, vector);
        mergeSortHelper(storage, middle + 1, right, c, vector);
        merge(vector, left, middle, right, c, storage);
    }
    else insertionSort(vector, left, right, c);
}
template<typename ITEM, typename COMPARATOR>
void mergeSort(ITEM* vector, int n, COMPARATOR const& c)
{//copy vector to storage first
    if(n <= 1) return;
    Vector<ITEM> storage(n, vector[0]);//reserve space for n with 1st item
    for(int i = 1; i < n; ++i) storage[i] = vector[i];
    mergeSortHelper(vector, 0, n - 1, c, storage.getArray());
}
template<typename ITEM> void mergeSort(ITEM* vector, int n)
    {mergeSort(vector, n, DefaultComparator<ITEM>());}

template<typename VECTOR, typename ITEM> struct DefaultVectorAccessor
{//for sorting vectors
    ITEM const& getItem(VECTOR const& v, int i)const{return v[i];}
    int getSize(VECTOR const& v)const{return v.getSize();}
};
struct StringAccessor
{//for sorting strings
    char getItem(string const& s, int i)const{return s[i];}
    int getSize(string const& s)const{return s.length();}
};
template<typename VECTOR, typename ITEM, typename ACCESSOR =
    DefaultVectorAccessor<VECTOR, ITEM> > struct IndexedAccessor
{
    VECTOR const*const v;
    ACCESSOR a;
    IndexedAccessor(VECTOR const*const theV, ACCESSOR const& theA = ACCESSOR())
        : v(theV), a(theA){}
    ITEM getItem(int i, int j)const{return a.getItem(v[i], j);}
    int getSize(int i)const{return a.getSize(v[i]);}
};
template<typename VECTOR, typename ITEM, typename ACCESSOR =
    DefaultVectorAccessor<VECTOR, ITEM>, typename COMPARATOR =
    DefaultComparator<ITEM> > struct MultikeyQuicksortVectorComparator
{
    ACCESSOR s;
    COMPARATOR c;
    mutable int depth;
    MultikeyQuicksortVectorComparator(ACCESSOR const& theS = ACCESSOR(),
        COMPARATOR const& theC = COMPARATOR()): s(theS), c(theC), depth(0){}
    bool operator()(VECTOR const& lhs, VECTOR const& rhs)const
    {
        return depth < s.getSize(lhs) ?
            depth < s.getSize(rhs) && c(s.getItem(lhs, depth),
            s.getItem(rhs, depth)) : depth < s.getSize(rhs);
    }
    bool isEqual(VECTOR const& lhs, VECTOR const& rhs)const
    {
        return depth < s.getSize(lhs) ?
            depth < s.getSize(rhs) && c.isEqual(s.getItem(lhs, depth),
            s.getItem(rhs, depth)) : depth >= s.getSize(rhs);
    }
};

template<typename VECTOR, typename COMPARATOR> void multikeyQuicksortNR(
    VECTOR* vector, int left, int right, COMPARATOR const& c,
    int maxDepth = numeric_limits<int>::max())
{
    Stack<int> stack;
    stack.push(left);
    stack.push(right);
    stack.push(0);
    while(!stack.isEmpty())
    {
        c.depth = stack.pop();
        right = stack.pop();
        left = stack.pop();
        if(right - left > 0 && c.depth < maxDepth)
        {
            int i, j;
            partition3(vector, left, right, i, j, c);
            //left
            stack.push(left);
            stack.push(j);
            stack.push(c.depth);
            //right
            stack.push(i);
            stack.push(right);
            stack.push(c.depth);
            //middle
            stack.push(j + 1);
            stack.push(i - 1);
            stack.push(c.depth + 1);
        }
    }
}

template<typename VECTOR, typename COMPARATOR> void multikeyQuickselect(
    VECTOR* vector, int left, int right, int k, COMPARATOR const& c)
{
    assert(k >= left && k <= right);
    for(int d = 0, i, j; right - left >= 1;)
    {
        partition3(vector, left, right, i, j, c);
        if(k <= j) right = j;
        else if (k < i)//equal case j < k < i
        {
            left = j + 1;
            right = i - 1;
            ++c.depth;
        }
        else left = i;
    }
}

template<typename ITEM, typename COMPARATOR> int binarySearch(ITEM const*
    vector, int left, int right, ITEM const& key, COMPARATOR const& c)
{
    while(left <= right)
    {//careful to avoid overflow in middle calculation
        int middle = left + (right - left)/2;
        if(c.isEqual(key, vector[middle])) return middle;
        c(key, vector[middle]) ? right = middle - 1 : left = middle + 1;
    }
    return -1;//not found
}

template<typename ITEM> void permutationSort(ITEM* a, int* permutation, int n)
{//need permutation to be valid, else get an infinite loop
    for(int i = 0; i < n; ++i) if(permutation[i] != i)
        {//start cycle
            ITEM temp = a[i];
            int to = i;
            do
            {
                a[to] = a[permutation[to]];//put element in right place
                swap(permutation[to], to);//mark to done, and advance cycle
            }while(permutation[to] != i);//until find what goes to i
            a[to] = temp;//complete cycle
            permutation[to] = to;//becomes identity
        }
}

void countingSort(int* vector, int n, int N)
{
    Vector<int> counter(N, 0);
    for(int i = 0; i < n; ++i) ++counter[vector[i]];//count
    for(int i = 0, index = 0; i < N; ++i)//create in order
        while(counter[i]-- > 0) vector[index++] = i;
}

template<typename ITEM> struct IdentityHash
    {int operator()(ITEM const& item)const{return item;}};
template<typename ITEM, typename ORDERED_HASH> void KSort(ITEM* a, int n,
    int N, ORDERED_HASH const& h)
{
    ITEM* temp = rawMemory<ITEM>(n);
    Vector<int> count(N + 1, 0);
    for(int i = 0; i < n; ++i) ++count[h(a[i]) + 1];
    for(int i = 0; i < N; ++i) count[i + 1] += count[i];//accumulate counts
    //rearrange items
    for(int i = 0; i < n; ++i) new(&temp[count[h(a[i])]++])ITEM(a[i]);
    for(int i = 0; i < n; ++i) a[i] = temp[i];
    rawDelete(temp);
}

template<typename ITEM, typename COMPARATOR> bool isSorted(ITEM const*
    vector, int left, int right, COMPARATOR const& c)
{
    for(int i = left + 1; i <= right; ++i)
        if(c(vector[i], vector[i - 1])) return false;
    return true;
}

}//end namespace
#endif
