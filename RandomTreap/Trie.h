#ifndef IGMDK_TRIE_H
#define IGMDK_TRIE_H
#include "../Utils/GCFreeList.h"
#include "../Utils/Vector.h"
#include "../RandomNumberGeneration/Random.h"
namespace igmdk{

template<typename ITEM, typename KEY_OBJECT = unsigned char, typename
    COMPARATOR = DefaultComparator<ITEM> > class TernaryTreapTrie
{
    COMPARATOR c;
    struct Node
    {
        KEY_OBJECT pivot;
        unsigned int priority;
        Node *next, *left, *right;
        ITEM* item;
        Node(KEY_OBJECT const& thePivot): next(0), left(0), right(0),
            item(0), pivot(thePivot), priority(GlobalRNG().next()){}
    }* root;
    Freelist<ITEM> itemF;
    Freelist<Node> nodeF;
    Node* rotateRight(Node* tree)
    {
        Node* goingUp = tree->left;
        tree->left = goingUp->right;
        goingUp->right = tree;
        return goingUp;
    }
    Node* rotateLeft(Node* tree)
    {
        Node* goingUp = tree->right;
        tree->right = goingUp->left;
        goingUp->left = tree;
        return goingUp;
    }

    Node* insertNode(Node* node, KEY_OBJECT* key, int keySize,
        ITEM const& item, int i)
    {
        if(!node) node = new(nodeF.allocate())Node(key[i]);
        if(c.isEqual(key[i], node->pivot))//go down
            if(i == keySize - 1)//found the wanted node
            {
                if(node->item) *node->item = item;
                else node->item = new(itemF.allocate())ITEM(item);
            }
            else node->next = insertNode(node->next, key, keySize, item,
                i + 1);//keep going
        else
        {//go sideways
            bool goLeft = c(key[i], node->pivot);
            Node*& chosenChild = goLeft ? node->left : node->right;
            chosenChild = insertNode(chosenChild, key, keySize, item, i);
            if(chosenChild->priority < node->priority)
                node = goLeft ? rotateRight(node) : rotateLeft(node);
        }
        return node;
    }

    Node* join(Node* left, Node* right)
    {
        if(!left) return right;
        if(!right) return left;
        if(left->priority < right->priority)//lower priority goes up
        {
            left->right = join(left->right, right);
            return left;
        }
        else
        {
            right->left = join(left, right->left);
            return right;
        }
    }
    Node* removeR(Node* node, KEY_OBJECT* key, int keySize, int i)
    {
        if(node)
        {
            bool isEqual = c.isEqual(key[i], node->pivot);
            if(isEqual && i == keySize - 1)
            {//remove found item if it exists
                if(!node->item) return node;
                itemF.remove(node->item);
                node->item = 0;
            }
            else
            {//go to next node
                Node** child;
                if(isEqual){child = &node->next; ++i;}
                else if(c(key[i], node->pivot)) child = &node->left;
                else child = &node->right;
                *child = removeR(*child, key, keySize, i);
            }
            if(!node->item && !node->next)
            {//remove empty node
                Node* left = node->left, *right = node->right;
                nodeF.remove(node);
                node = (left || right) ? join(left, right) : 0;
            }
        }
        return node;
    }
    template<typename ACTION> void forEachNode(Node* node, ACTION& action)
    {
        if(node)
        {
            action(node);
            forEachNode(node->left);
            forEachNode(node->next);
            forEachNode(node->right);
        }
    }
    struct CopyAction
    {
        Vector<ITEM*>& result;
        CopyAction(Vector<ITEM*>& theResult): result(theResult){}
        void operator()(Node* node){if(node->item)result.append(node->item);}
    };
    Node* constructFrom(Node* node)
    {
        Node* result = 0;
        if(node)
        {
            result = new(nodeF.allocate())Node(node->pivot);
            result->priority = node->priority;
            if(node->item)result->item = new(itemF.allocate())ITEM(node->item);
            result->left = constructFrom(node->left);
            result->next = constructFrom(node->next);
            result->right = constructFrom(node->right);
        }
        return result;
    }
public:
    TernaryTreapTrie(COMPARATOR const& theC = COMPARATOR()): root(0), c(theC){}
    TernaryTreapTrie(TernaryTreapTrie const& other):
        c(other.c){root = constructFrom(other.root);}
    TernaryTreapTrie& operator=(TernaryTreapTrie const& rhs)
        {return genericAssign(*this, rhs);}

    void insert(unsigned char* key, int keySize, ITEM const& item)
    {
        assert(keySize > 0);
        root = insertNode(root, key, keySize, item, 0);
    }
    Vector<ITEM*> prefixFind(KEY_OBJECT* key, int lcp)
    {
        Vector<ITEM*> result;
        CopyAction action(result);
        forEachNode(findNode(key, lcp), action);
        return result;
    }
    ITEM* longestMatch(KEY_OBJECT* key, int keySize)
    {
        Node* node = root, *result = 0;
        for(int i = 0; node && i < keySize;)
            if(c.isEqual(key[i], node->pivot))
            {
                result = node;
                if(i == keySize - 1) break;//reached last key object
                else {node = node->next; ++i;}
            }
            else if(c(key[i], node->pivot)) node = node->left;
            else node = node->right;
        return result ? result->item : 0;
    }
    struct Handle
    {
        Node* node;
        int i;
        Handle(Node* theNode = 0, int theI = 0): node(theNode), i(theI){}
    };
    Node* findNodeIncremental(KEY_OBJECT* key, int keySize, Handle& h)
    {
        while(h.node && h.i < keySize)
        {
            if(c.isEqual(key[h.i], h.node->pivot))
            {
                if(h.i == keySize - 1) return h.node;
                else{h.node = h.node->next; ++h.i;}
            }
            else if(c(key[h.i], h.node->pivot)) h.node = h.node->left;
            else h.node = h.node->right;
        }
        h = Handle();//not found
        return 0;
    }
    ITEM* findIncremental(KEY_OBJECT* key, int keySize, Handle& h)
    {
        if(!h.node) h.node = root;
        Node* result = findNodeIncremental(key, keySize, h);
        return result ? result->item : 0;
    }
    Node* findNode(KEY_OBJECT* key, int keySize)
    {
        Handle h(root, 0);
        return findNodeIncremental(key, keySize, h);
    }
    ITEM* find(KEY_OBJECT* key, int keySize)
    {
        Node* result = findNode(key, keySize);
        return result ? result->item : 0;
    }

    void remove(KEY_OBJECT* key, int keySize)
    {
        assert(keySize > 0);
        root = removeR(root, key, keySize, 0);
    }
};

}//end namespace
#endif
