#ifndef IGMDK_PAIRING_HEAP_H
#define IGMDK_PAIRING_HEAP_H

#include "../Utils/GCFreeList.h"
#include "../Utils/Utils.h"
#include "../HashTable/LinearProbingHashTable.h"

namespace igmdk{

template<typename ITEM, typename COMPARATOR = DefaultComparator<ITEM> >
class PairingHeap
{
    /*Delete in pairing heap is unlink tree with the element, delete min on it,
    then if the deleted element was not the root merge with the root. To
    increase key need to cut children and put in siblings list. These
    additional operations are not as useful and have not been implemented.*/
    COMPARATOR c;
    int size;
	struct Node
	{
		ITEM element;
		Node* elder, *oldestChild, *youngerSibling;
		Node(ITEM const& theItem)
		:element(theItem), elder(0), oldestChild(0), youngerSibling(0){}
	} *root;
	Freelist<Node> freeList;
	PairingHeap(PairingHeap<ITEM> const&);
	PairingHeap& operator=(PairingHeap<ITEM> const&);
	void insertNode(Node* node)
	{
		if(c(root->element, node->element)) linkRoots(root, node);
		else{linkRoots(node, root);root = node;}
	}
	void linkRoots(Node* parent, Node* child)
	{
		//current youngerSiblings are ignored
		Node* oldestChild = parent->oldestChild;
		child->youngerSibling = oldestChild;
		if(oldestChild) oldestChild->elder = child;
		parent->oldestChild = child;
		child->elder = parent;
	}
	void pairUp()
	{
		for(Node* currentRoot = root->youngerSibling; currentRoot; )
		{
			Node* nextRoot = currentRoot->youngerSibling;
			if(!nextRoot)
            {
                insertNode(currentRoot);
                break;
            }
			Node* youngerSibling = nextRoot->youngerSibling;
			if(c(currentRoot->element, nextRoot->element))
                swap(nextRoot, currentRoot);
            linkRoots(nextRoot, currentRoot);
            insertNode(nextRoot);
			currentRoot = youngerSibling;
		}
		correctRoot();
	}
	void correctRoot(){root->elder = root->youngerSibling = 0;}
public:
	typedef Node* Handle;
	PairingHeap():	root(0), size(0){}
	Handle insert(ITEM const& theItem)
	{
	    ++size;
		Node* node = new(freeList.allocate())Node(theItem);
		if(isEmpty()) root = node;
		else insertNode(node);
		return node;
		//assert(client will not corrupt the heap through newRoot)
	}
	Handle replaceMin(ITEM const& theItem)
	{
	    deleteMin();
	    return insert(theItem);
	}
	void decreaseKey(Handle node, ITEM const& newItem)
	{
		assert(node && !c(node->element, newItem) && !isEmpty());
		node->element = newItem;
		if(node != root)
		{
			Node *elder = node->elder, *youngerSibling = node->youngerSibling;
			if(youngerSibling) youngerSibling->elder = elder;
			(elder->oldestChild == node ?
				elder->oldestChild : elder->youngerSibling) = youngerSibling;
			insertNode(node);
			correctRoot();
		}
	}
	bool isEmpty(){return !root;}
	ITEM const& getMin(){assert(!isEmpty());return root->element;}
	ITEM deleteMin()
	{
	    --size;
		ITEM result = getMin();
		Node* oldRoot = root;
		root = root->oldestChild;
		freeList.remove(oldRoot);
		if(root) pairUp();
		return result;
	}
	int getSize(){return size;}
	//make sure external freelist is used when merging is needed!
	void merge(Node* otherRoot){insertNode(otherRoot);}
};

template<typename ITEM> struct IndexedPaHeap
{
    typedef typename PairingHeap<ITEM>::Handle HANDLE;
    LinearProbingHashTable<int, HANDLE> map;
    PairingHeap<ITEM> heap;
public:
    bool isEmpty(){return heap.isEmpty();}
    void insert(ITEM const& item, int handle)
        {map.insert(handle, heap.insert(item));}
    ITEM const& getMin(){return heap.getMin();}
    ITEM deleteMin(){return heap.deleteMin();}
    void changeKey(ITEM const& item, int handle)
	{
	    HANDLE h = map.find(handle);
        if(h) heap.changeKey(h->item, Item(item, h));
        else insert(item, handle);
	}
	void deleteKey(int handle)
    {
        HANDLE h = map.find(handle);
        assert(h);
        heap.remove(*h);
    }
};

}//end namespace
#endif
