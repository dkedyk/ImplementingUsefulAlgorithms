#ifndef IGMDK_SYMMETRIC_MIN_MAX_HEAP_H
#define IGMDK_SYMMETRIC_MIN_MAX_HEAP_H

#include "../Utils/Vector.h"

namespace igmdk{

template<typename ITEM, typename COMPARATOR = DefaultComparator<ITEM> >
class SymmetricMinMaxHeap
{
	Vector<ITEM> items;
	COMPARATOR c;
	int getParent(int i){return i/2 - 1;}
	int getLeftChild(int i){return 2*i+2;}
	int getRightChild(int i){return 2*i+3;}
	int getLeftSibling(int i){return getLeftChild(getParent(i));}
	void moveUp(int i)
	{
	    for(;;)
	    {
	        int next = getParent(getParent(i)), ppl, ppr,
                pl = getLeftSibling(i);
	        if(c(items[i], items[pl])) next = pl;//p1
            else if(c(items[i], items[ppl = getLeftChild(next)]))
                next = ppl;//p2
            else if(c(items[ppr = getRightChild(next)], items[i]))
                next = ppr;//p3
            else return;
            swap(items[next], items[i]);
	        i = next;
	    }
	}
	void moveDownMin(int i)
	{
        for(;;)
	    {
            int next = getParent(i), l = getLeftChild(i),
                pr = getRightChild(next);
            if(l < items.getSize())
            {
                int prl = getLeftChild(pr);
                next = prl < items.getSize() && c(items[prl], items[l]) ?
                    prl : l;
                if(c(items[i], items[next])) return;
            }
            else if(pr < items.getSize() && c(items[pr], items[i])) next = pr;
            else return;
            swap(items[next], items[i]);
	        i = next;
	    }
	}
	void moveDownMax(int i)
	{
        for(;;)
	    {
            int next, pl = getLeftChild(getParent(i)), plr = getRightChild(pl);
            if(plr < items.getSize())
            {
                int r = getRightChild(i),
                    child = r < items.getSize() ? r : getLeftChild(i);
                next = child < items.getSize() && c(items[plr], items[child]) ?
                    child : plr;
                if(c(items[next], items[i])) return;
            }
            else if(c(items[i], items[pl])) next = pl;
            else return;
            swap(items[next], items[i]);
	        i = next;
	    }
	}
public:
	bool isEmpty(){return items.getSize() <= 0;}
	ITEM const& getMin(){assert(!isEmpty());return items[0];}
	ITEM const& getMax(){assert(!isEmpty());return items[items.getSize()!= 1];}
	void insert(ITEM const& item)
	{
		items.append(item);
		if(items.getSize() > 1) moveUp(items.getSize()-1);
	}
	void deleteMin()
	{
		assert(!isEmpty());
		items[0] = items.lastItem();
		moveDownMin(0);
		items.removeLast();
	}
	void deleteMax()
	{
		assert(!isEmpty());
		int index = items.getSize() != 1;
		items[index] = items[items.getSize()-1];
		moveDownMax(index);
		items.remove(items.getSize()-1);
	}
};

}
#endif
