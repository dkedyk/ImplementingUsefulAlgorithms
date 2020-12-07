#ifndef IGMDK_LEFT_LEANING_RED_BLACK_TREE_H
#define IGMDK_LEFT_LEANING_RED_BLACK_TREE_H

#include "../RandomNumberGeneration/Random.h"
#include "../Utils/GCFreeList.h"
using namespace std;

namespace igmdk{

template<typename KEY, typename ITEM>
class LeftLeaningRedBlackTree
{
private:
	enum{RED = 1, BLACK = 0};
	struct Node
	{
		KEY key;
		ITEM item;
		Node *left, *right;
		bool color;
		Node(KEY const& theKey, ITEM const& theItem, bool theColor):
		    key(theKey), item(theItem), left(0), right(0), color(theColor){}

	}*root;
	Freelist<Node> freelist;
	bool isRed(Node* x){return x && x->color == RED;}
	Node* rotateLeft(Node* h)
	{
		Node* x = h->right;
		h->right = x->left;
		x->left = h;
		x->color = x->left->color;
		x->left->color = RED;
		return x;
	}
	Node* rotateRight(Node* h)
	{
		Node* x = h->left;
		h->left = x->right;
		x->right = h;
		x->color = x->right->color;
		x->right->color = RED;
		return x;
	}
	Node* colorFlip(Node* x)
	{
		x->color = !x->color;
		x->left->color = !x->left->color;
		x->right->color = !x->right->color;
		return x;
	}
	Node* fixUp(Node* h)
	{
		if(isRed(h->right)) h = rotateLeft(h);
		if(isRed(h->left) && isRed(h->left->left)) h = rotateRight(h);
		if(isRed(h->left) && isRed(h->right)) colorFlip(h);
		return h;
	}
	Node* moveRedLeft(Node* h)
	{
		colorFlip(h);
		if (isRed(h->right->left))
		{
			h->right = rotateRight(h->right);
			h = rotateLeft(h);
			colorFlip(h);
		}
		return h;
	}
	Node* moveRedRight(Node* h)
	{
		colorFlip(h);
		if(isRed(h->left->left))
		{
			h = rotateRight(h);
			colorFlip(h);
		}
		return h;
	}
	Node* insert(Node* h, KEY const& key, ITEM const& item)
	{
		if(!h) return new(freelist.allocate())Node(key, item, RED);
		if(isRed(h->left) && isRed(h->right)) colorFlip(h);
		Node*& child = key < h->key ? h->left : h->right;
		child = insert(child, key, item);
		if(isRed(h->right)) h = rotateLeft(h);
		if(isRed(h->left) && isRed(h->left->left)) h = rotateRight(h);
		return h;
	}
	Node* deleteMin(Node* h)
	{
		if(h->left == 0) return 0;
		if(!isRed(h->left) && !isRed(h->left->left)) h = moveRedLeft(h);
		h->left = deleteMin(h->left);
		return fixUp(h);
	}
	KEY* min(Node* x)
	{
		Node* result = x;
		if(result)
		{
			while(result = result->left);
			return &result->key;
		}
		return 0;
	}
	Node* remove(Node* h, KEY const& key)
	{
		if(key < h->key)
		{
			if (!isRed(h->left) && !isRed(h->left->left)) h = moveRedLeft(h);
			h->left = remove(h->left, key);
		}
		else
		{
			if (isRed(h->left)) h = rotateRight(h);
			if (key == h->key && (h->right == 0)) return 0;
			if (!isRed(h->right) && !isRed(h->right->left))h = moveRedRight(h);
			if (key == h->key)
			{
				h->key = *min(h->right);
				h->item = *find(h->right, h->key);
				h->right = deleteMin(h->right);
			}
			else h->right = remove(h->right, key);
		}
		return fixUp(h);
	}
	ITEM* find(Node* node, KEY const& key)
	{
		Node* result = root;
		while(result && key != result->key)
            {result = key < result->key ? result->left : result->right;}
		return result ? &result->item : 0;
	}
public:
	LeftLeaningRedBlackTree(): root(0){}
	void insert(KEY const& key, ITEM const& item)
        {root = insert(root, key, item);}
	void remove(KEY const& key){root = remove(root, key);}
	ITEM* find(KEY const& key){return find(root, key);}
};

}//end namespace
#endif
