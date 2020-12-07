#ifndef IGMDK_SCRAP_HEAP_H
#define IGMDK_SCRAP_HEAP_H

#include "../Utils/Vector.h"
#include "../Utils/GCFreeList.h"

namespace igmdk{

template<typename ITEM, typename COMPARATOR = DefaultComparator<ITEM> >
class IndexedPointerHeap
{
    COMPARATOR c;
    struct Item
	{
		ITEM item;
		int index;
		Item(ITEM const& theItem, int theIndex):item(theItem),index(theIndex){}
	};
	Freelist<Item> freelist;
    Vector<Item*> items;
    int getParent(int i){return (i-1)/2;}
	int getLeftChild(int i){return 2*i+1;}
	void report(Item* item, int i){item->index = i;}
	void moveUp(int i)
	{
		Handle temp = items[i];
		for(int parent; i > 0 && c(temp->item,
            items[parent = getParent(i)]->item); i = parent)
            report(items[i] = items[parent], i);
		report(items[i] = temp, i);
	}
	void moveDown(int i)
	{
		Handle temp = items[i];
		for(int child; (child = getLeftChild(i)) <
            items.getSize(); i = child)
		{
			int rightChild = child + 1;
			if(rightChild < items.getSize() && c(items
                [rightChild]->item, items[child]->item)) child = rightChild;
			if(!c(items[child]->item, temp->item)) break;
			report(items[i] = items[child], i);
		}
		report(items[i] = temp, i);
	}
	void heapify()
        {for(int i = getParent(items.getSize()-1); i >= 0; --i) moveDown(i);}
    void remove(int i)
    {
        freelist.remove(items[i]);
        if(items.getSize() > 1)
        {
            items[i] = items.lastItem();
            moveDown(i);
        }
        items.removeLast();
    }
public:
    typedef Item* Handle;
    bool isEmpty(){return items.getSize() <= 0;}
	ITEM const& getMin(){assert(!isEmpty());return items[0]->item;}
    Handle insert(ITEM const& item)
    {
        Handle result = new(freelist.allocate())Item(item, items.getSize());
        items.append(result);
		moveUp(items.getSize()-1);
		return result;
    }
	ITEM deleteMin(){ITEM result = getMin();remove(0);return result;}
	void changeKey(Handle pointer, ITEM const& item)
	{//assert(pointer && pointer is not garbage);
		bool decrease = c(item, pointer->item);
		pointer->item = item;
		decrease ? moveUp(pointer->index) : moveDown(pointer->index);
	}
	void decreaseKey(Handle index, ITEM const& item){changeKey(index, item);}
    void deleteKey(Handle pointer){remove(pointer->index);}
};

template<typename ITEM> class BucketQueue
{
	int capacity, minIndex;
	struct Node
	{
		ITEM item;
		Node* next;
		Node(ITEM const& theItem, Node* theNext):item(theItem), next(theNext){}
	}** buckets;
	Freelist<Node> freelist;
public:
	BucketQueue(int maxN): capacity(maxN + 1), buckets(new Node*[capacity]),
        minIndex(capacity){for(int i = 0; i < capacity; ++i)buckets[i] = 0;}
	void insert(int priority, ITEM const& item)
	{
		assert(priority < capacity);
		buckets[priority] =
            new(freelist.allocate())Node(item, buckets[priority]);
		if(priority < minIndex) minIndex = priority;
	}
	bool isEmpty(){return minIndex == capacity;}
	ITEM findMin(){assert(!isEmpty()); return buckets[minIndex]->item;}
	void deleteMin()
	{
		assert(!isEmpty());
		Node* temp = buckets[minIndex];
		buckets[minIndex] = temp->next;
		freelist.remove(temp);
		while(minIndex < capacity && !buckets[minIndex]) ++minIndex;
	}
	~BucketQueue(){delete[] buckets;}
};
}
#endif
