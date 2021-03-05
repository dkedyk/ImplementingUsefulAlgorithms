#ifndef IGMDK_LCP_TREAP_H
#define IGMDK_LCP_TREAP_H
#include "../RandomNumberGeneration/Random.h"
#include "../Utils/GCFreeList.h"
#include "Treap.h"

namespace igmdk{

template<typename KEY, typename VALUE, typename INDEXED_COMPARATOR =
    LexicographicComparator<KEY> > class LCPTreap
{
    INDEXED_COMPARATOR c;
    struct Node
    {
        KEY key;
        VALUE value;
        Node *left, *right, *parent;
        unsigned int priority, nodeCount;
        unsigned short predLcp, succLcp;
        Node(KEY const& theKey, VALUE const& theValue): key(theKey),
            value(theValue), left(0), right(0), priority(GlobalRNG().next()),
            predLcp(0), succLcp(0), parent(0), nodeCount(1) {}
    }* root;
    Freelist<Node> f;
    Node* rotateRight(Node* node)
    {
        Node *goingUp = node->left, *&movedChild = goingUp->right;
        node->predLcp = goingUp->succLcp;
        goingUp->succLcp = min(node->succLcp, goingUp->succLcp);
        node->left = movedChild;
        TreapGeneric<Node>::rotateHelper(node, goingUp, movedChild);
        return goingUp;
    }
    Node* rotateLeft(Node* node)
    {
        Node *goingUp = node->right, *&movedChild = goingUp->left;
        node->succLcp = goingUp->predLcp;
        goingUp->predLcp = min(node->predLcp, goingUp->predLcp);
        node->right = movedChild;
        TreapGeneric<Node>::rotateHelper(node, goingUp, movedChild);
        return goingUp;
    }
    int findLCP(KEY const& key, Node* node, int predM, int& m)
    {
        int lcp = predM == m ? node->predLcp : node->succLcp;//get l or r
        if(lcp >= m)
        {
            while(m < c.getSize(key) && c.isEqual(key, node->key, m)) ++m;
            lcp = m;
        }
        return lcp;
    }
    Node* constructFrom(Node* node)
    {
        Node* tree = 0;
        if(node)
        {
            tree = new(f.allocate())Node(node->key, node->value);
            tree->priority = node->priority;
            tree->nodeCount = node->nodeCount;
            tree->predLcp = node->predLcp;
            tree->succLcp = node->succLcp;
            tree->left = constructFrom(node->left);
            if(tree->left) tree->left->parent = tree;
            tree->right = constructFrom(node->right);
            if(tree->right) tree->left->parent = tree;
        }
        return tree;
    }
public:
    typedef Node NodeType;
    LCPTreap(INDEXED_COMPARATOR theC = INDEXED_COMPARATOR()):root(0), c(theC){}
    LCPTreap(LCPTreap const& other): c(other.c)
        {root = constructFrom(other.root);}
    LCPTreap& operator=(LCPTreap const&rhs){return genericAssign(*this,rhs);}
    unsigned int getSize(){return root ? root->nodeCount : 0;}
    Node* findNode(KEY const& key)
    {
        Node* node = root;
        int m = 0, predM = 0;
        while(node)
        {
            int lcp = findLCP(key, node, predM, m);
            if(c.getSize(key) == lcp && c.getSize(node->key) == lcp) break;
            if(c(key, node->key, lcp)) node = node->left;
            else
            {
                node = node->right;
                predM = lcp;
            }
        }
        return node;
    }
    VALUE* find(KEY const& key)
    {
        Node* node = findNode(key);
        return node ? &node->value : 0;
    }
    Node* insertNode(Node* newNode, Node* node, int m = 0)
    {
        if(!node) return newNode;
        int lcp = findLCP(newNode->key, node, newNode->predLcp, m);
        bool goLeft = c(newNode->key, node->key, lcp);
        (goLeft ? newNode->succLcp : newNode->predLcp) = lcp;
        Node*& chosenChild = goLeft ? node->left : node->right;
        chosenChild = insertNode(newNode, chosenChild, m);
        chosenChild->parent = node;
        ++node->nodeCount;
        if(chosenChild->priority < node->priority)
            node = goLeft ? rotateRight(node) : rotateLeft(node);
        return node;
    }
    Node* insert(KEY const& key, VALUE const& value)
    {
        Node* node = findNode(key);
        if(node) node->value = value;
        else
        {
            node = new(f.allocate())Node(key, value);
            root = insertNode(node, root);
        }
        return node;
    }
    Node* removeFound(Node* node)
    {
        Node *left = node->left, *right = node->right;
        if(left && right)
        {
            bool goRight = left->priority < right->priority;
            node = goRight ? rotateRight(node) : rotateLeft(node);
            Node*& child = goRight ? node->right : node->left;
            child = removeFound(child);
            if(child) child->parent = node;
            --node->nodeCount;
        }
        else
        {
            f.remove(node);
            node = left ? left : right;
        }
        return node;
    }
    void remove(KEY const& key)
    {
        Node* node = findNode(key);
        if(node)
        {
            Node* parent = node->parent;
            bool wasLeft = parent && node == parent->left;
            node = removeFound(node);
            if(node) node->parent = parent;
            (parent ? (wasLeft ? parent->left : parent->right) : root) = node;
            for(; parent; parent = parent->parent) --parent->nodeCount;
        }
    }
    NodeType* predecessor(KEY const& key)
    {
        int m = 0, predM = 0;
        Node* pred = 0;
        for(Node* node = root; node;)
        {
            int lcp = findLCP(key, node, predM, m);
            if(c(node->key, key, lcp))//found pred set member
            {
                pred = node;
                node = node->right;
                predM = lcp;
            }
            else node = node->left;
        }
        return pred;
    }
    NodeType* inclusivePredecessor(KEY const& key)
    {
        Node* node = findNode(key);
        return node ? node : predecessor(key);
    }
    NodeType* successor(KEY const& key)
    {
        int m = 0, predM = 0;
        Node* succ = 0;
        for(Node* node = root; node;)
        {
            int lcp = findLCP(key, node, predM, m);
            if(c(key, node->key, lcp))//found succ set member
            {
                succ = node;
                node = node->left;
                predM = lcp;
            }
            else node = node->right;
        }
        return succ;
    }
    NodeType* inclusiveSuccessor(KEY const& key)
    {
        Node* node = findNode(key);
        return node ? node : successor(key);
    }
    NodeType* prefixSuccessor(KEY const& prefix)
    {
        int m = 0, predM = 0;
        for(Node* node = root; node;)
        {
            int lcp = findLCP(prefix, node, predM, m);
            if(c.getSize(prefix) == lcp || !c(prefix, node->key, lcp))
            {
                node = node->right;
                predM = lcp;
            }
            else if(!node->left) return node;
            else node = node->left;
        }
        return 0;
    }
    NodeType* findMin(){return TreapGeneric<Node>::findMin(root);}
    NodeType* findMax(){return TreapGeneric<Node>::findMax(root);}
    NodeType* nthElement(int n)
        {return TreapGeneric<Node>::nthElement(n, root);}
    typedef TreeIterator<Node> Iterator;
    Iterator begin(){return Iterator(findMin());}
    Iterator end(){return Iterator(0);}
    Iterator rBegin(){return Iterator(findMax());}
    Iterator rEnd(){return Iterator(0);}
};

}//end namespace
#endif
