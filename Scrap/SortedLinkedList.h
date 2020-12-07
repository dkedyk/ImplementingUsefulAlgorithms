#ifndef IGMDK_SORTED_LINKED_LIST_H
#define IGMDK_SORTED_LINKED_LIST_H

#include "../Utils/GCFreeList.h"
#include "../Utils/Utils.h"
namespace igmdk{

template<typename ITEM, typename COMPARATOR = DefaultComparator<ITEM> >
struct LinkedList
{
	struct Node
	{
		ITEM item;
		Node* next;
		Node(ITEM const& theItem, Node* theNext)
		:	item(theItem), next(theNext)	{}
	}*root;
	int size;
	Freelist<Node> freelist;
	COMPARATOR c;
	LinkedList(COMPARATOR theC = COMPARATOR()):	root(0), size(0), c(theC) {}
	void prepend(ITEM const& item)
	{
		root = new(freelist.allocate())Node(item, root);
		++size;
	}
	Node* advanceSmaller(Node*& a, Node*& b)
	{
		Node*& smaller = c(b->item, a->item) ? b : a;
		Node* nodeToAppend = smaller;
		smaller = smaller->next;
		return nodeToAppend;
	}
	Node* merge(Node* a, Node* b)
	{//pick head
		Node* head = advanceSmaller(a, b), *tail = head;
		//append from smaller until one runs out
		while(a && b) tail = tail->next = advanceSmaller(a, b);
		//append the rest of the remaining list
		tail->next = a ? a : b;
		return head;
	}
	Node* mergesort(Node* list, int n, Node*& nextAfterLast)
	{
		if(n==1)
		{
			nextAfterLast = list->next;
			list->next = 0;
			return list;
		}
		int middle = n/2;
		Node *secondHalf, *m1 = mergesort(list, middle, secondHalf);
		return merge(m1, mergesort(secondHalf, n - middle, nextAfterLast));
	}
	void sort()
	{
		if(size > 1)
		{
			Node* dummy;
			root = mergesort(root, size, dummy);
		}
	}
};

}//end namespace
#endif
