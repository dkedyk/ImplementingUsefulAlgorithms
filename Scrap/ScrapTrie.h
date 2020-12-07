#ifndef IGMDK_SCRAP_TRIE_H
#define IGMDK_SCRAP_TRIE_H

#include "../Utils/GCFreeList.h"
#include <cassert>
#include "../Utils/Debug.h"
using namespace std;
namespace igmdk{

template<typename KEY> struct DefaultRank
{
    unsigned char* array;
    int size;
	DefaultRank(KEY const& key)
	{//works with pod types only
		array = (unsigned char*)&key;
		size = sizeof(key);
	}
};

template<typename KEY, typename ITEM, typename RANK = DefaultRank<KEY> >
class PatriciaTrie
{
	enum{B = 8, TABLE_SIZE = 2, S = B-1};
	int extract(unsigned char* key, int i)
        {return (unsigned char)(key[i/B] << (i%B)) >> S;}
	struct Node
	{
		KEY key;
		ITEM item;
		int index;
		Node *next[2];
		Node(KEY const& theKey, ITEM const& theItem, int theTo, bool bit,
            Node* other): key(theKey), item(theItem), index(theTo)
		{
		    next[!bit] = other;
		    next[bit] = this;
		}
	}*root;
	Freelist<Node> freelist;
    template<typename ACTION>
	void forEachHelper(Node* t, ACTION& action)
	{
        for(int i = 0; t && i < 2; ++i)
        {
            Node* child = t->next[i];
            if(child)
            {
                if(t->index >= child->index) action(child);
                else forEachHelper(child, action);
            }
        }
	}
	struct AppendAction
	{
	    Vector<Node*>& result;
	    AppendAction(Vector<Node*>& theResult):result(theResult){}
	    void operator()(Node* node){result.append(node);}
	};
	Node** findForwardPointer(Node* query, Node** tree)
	{
	    RANK rank(query->key);
		Node* node;
		while((node = *tree) != query)
		{
		    tree = &node->next[extract(rank.array, node->index)];
		}
		return tree;
	}
	Node** findBackwardPointer(Node* query, Node** tree)
	{
	    RANK rank(query->key);
		int prevIndex = -1;
		Node** pointer = tree;
		Node* node;
		while((node = *pointer) && prevIndex < node->index)
		{
			prevIndex = node->index;
            pointer = &node->next[extract(rank.array, prevIndex)];
		}
		return pointer;
	}
	void removeSingleChildNode(Node** forwardPointer)
	{
	    Node* node = *forwardPointer;
        *forwardPointer = node->next[!node->next[0]];
        freelist.remove(node);
	}
public:
	PatriciaTrie():	root(0){}
	//acts on nodes in lexicographic bit order
	template<typename ACTION> void forEach(ACTION& action)
        {forEachHelper(root, action);}
	Node* findLongestMatch(KEY const& key)
	{//will crash if key is prefix of another
	    if(!root) return 0;
	    RANK rank(key);
		int prevIndex = -1;
		Node* node = root;
		while(prevIndex < node->index)
		{
			prevIndex = node->index;
			Node* next = node->next[extract(rank.array, prevIndex)];
			if(!next) break;
            node = next;
		}
		return node;
		/*
		Longest match because all nodes of untested bits on the path have the
		same bits and this is the one that can match furtherst. Other nodes
		can't be eliminated by untested bits because they are the same their
		because the bits are not tested or by tested bits because if so search
		would take to different path
		*/
	}
	ITEM* find(KEY const& key)
	{
	    Node* node = findLongestMatch(key);
	    return node && key == node->key ? &node->item : 0;
	}
	void insert(KEY const& key, ITEM const& item)
	{//1. Find index of key to be inserted
	    Node* lcpNode = findLongestMatch(key);
	    RANK rank(key);
	    int theIndex = 0;
	    if(lcpNode)
	    {
	        if(key == lcpNode->key){lcpNode->item = item; return;}
            RANK rank2(lcpNode->key);
            while(extract(rank.array, theIndex) ==
                extract(rank2.array, theIndex)) ++theIndex;
            if(theIndex == lcpNode->index && theIndex < rank.size * B)
                ++theIndex;
	    }//2. Create and insert the node
	    Node **pointer = &root, *node;
		int prevIndex = -1;
		//new node goes before an existing one if node->index >= theIndex or
		//after a self pointing node if prevIndex >= node->index
		while((node = *pointer) && node->index < theIndex &&
            prevIndex < node->index)
		{
		    prevIndex = node->index;
		    pointer = &node->next[extract(rank.array, prevIndex)];
		}
		*pointer = new(freelist.allocate())
            Node(key, item, theIndex, extract(rank.array, theIndex), node);
	}
	Vector<Node*> prefixFind(KEY const& key, int minLCP)
	{
	    RANK rank(key);
		int prevIndex = -1;
		Node* node = root;
		while(node && node->index < minLCP && prevIndex < node->index)
		{
			prevIndex = node->index;
			node = node->next[extract(rank.array, prevIndex)];
		}
		Vector<Node*> result;
		if(node && node->index >= minLCP - 1)//&& findlcp is good
		{
		    if(prevIndex >= node->index) result.append(node);
		    else
		    {
		        AppendAction action(result);
		        forEachHelper(node, action);
		    }
		}
        return result;
	}
	void remove(KEY const& key)
    {
        RANK rank(key);
		int prevIndex = -1;
		Node** pointer = &root, **parentForwardPointer = &root;
		Node* node, *parent = 0;
		while((node = *pointer) && prevIndex < node->index)
		{
			prevIndex = node->index;
			parentForwardPointer = pointer;
            pointer = &node->next[extract(rank.array, prevIndex)];
		}
		if(node && key == node->key)
		{
		    Node* parent = *parentForwardPointer;
		    *pointer = 0;
		    //if self pointer then remove and link forward pointer to the other
		    //child
            if(node == parent)
                removeSingleChildNode(findForwardPointer(node, &root));
            else
            {//if not replace item by that of the parent, redirect parent's
                //backward pointer to it and remove parent as single child node
                node->key = parent->key;
                node->item = parent->item;
                *findBackwardPointer(parent, parentForwardPointer) = node;
                removeSingleChildNode(parentForwardPointer);
            }
		}
    }
};

}
#endif
