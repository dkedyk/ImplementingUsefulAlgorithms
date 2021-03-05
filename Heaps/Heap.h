#ifndef IGMDK_HEAP_H
#define IGMDK_HEAP_H
#include "../Utils/Vector.h"
namespace igmdk{

template<typename ITEM>
struct ReportDefault{void operator()(ITEM& item, int i){}};
template<typename ITEM, typename COMPARATOR = DefaultComparator<ITEM>,
    typename REPORTER = ReportDefault<ITEM> > class Heap
{
    REPORTER r;
    int getParent(int i)const{return (i - 1)/2;}
    int getLeftChild(int i)const{return 2 * i + 1;}
    Vector<ITEM> items;
    void moveUp(int i)
    {
        ITEM temp = items[i];
        for(int parent; i > 0 && c(temp, items[parent = getParent(i)]);
            i = parent) r(items[i] = items[parent], i);
        r(items[i] = temp, i);
    }
    void moveDown(int i)
    {
        ITEM temp = items[i];
        for(int child; (child = getLeftChild(i)) < items.getSize(); i = child)
        {//find smaller child
            int rightChild = child + 1;
            if(rightChild < items.getSize() && c(items
                [rightChild], items[child])) child = rightChild;
            //replace with the smaller child if any
            if(!c(items[child], temp)) break;
            r(items[i] = items[child], i);
        }
        r(items[i] = temp, i);
    }
public:
    COMPARATOR c;
    Heap(COMPARATOR const& theC = COMPARATOR(), REPORTER const&
        theReporter = REPORTER()): r(theReporter), c(theC) {}
    bool isEmpty()const{return items.getSize() == 0;}
    int getSize()const{return items.getSize();}
    ITEM const& getMin()const
    {
        assert(!isEmpty());
        return items[0];
    }
    void insert(ITEM const& item)
    {
        items.append(item);
        moveUp(items.getSize() - 1);
    }
    ITEM const& operator[](int i)const
    {//random access is useful with item handles
        assert(i >= 0 && i < items.getSize());
        return items[i];
    }
    void changeKey(int i, ITEM const& item)
    {
        assert(i >= 0 && i < items.getSize());
        bool decrease = c(item, items[i]);
        items[i] = item;
        decrease ? moveUp(i) : moveDown(i);
    }
    ITEM deleteMin(){return remove(0);}
    ITEM remove(int i)
    {
        assert(i >= 0 && i < items.getSize());
        ITEM result = items[i];
        r(result, -1);
        if(items.getSize() > i)
        {//not last item
            items[i] = items.lastItem();
            r(items[i], i);//report move
            moveDown(i);//won't touch the last item
        }
        items.removeLast();
        return result;
    }
};

}//end namespace
#endif
