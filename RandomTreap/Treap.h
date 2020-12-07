#ifndef IGMDK_TREAP_H
#define IGMDK_TREAP_H
#include "../RandomNumberGeneration/Random.h"
#include "../Utils/GCFreeList.h"
#include "../Utils/Stack.h"
namespace igmdk{

template<typename NODE> class TreeIterator
{
    NODE* current;
public:
    TreeIterator(NODE* node){current = node;}
    TreeIterator& operator++()
    {
        assert(current);
        if(current->right)
        {//if have right child go there, then maximally left
            current = current->right;
            while(current->left) current = current->left;
        }
        else
        {//parent if came from left child, else keep going up
            while(current->parent && current != current->parent->left)
                current = current->parent;
            current = current->parent;
        }
        return *this;
    }
    TreeIterator& operator--()
    {
        assert(current);
        if(current->left)
        {//if have left child go there, then maximally right
            current = current->left;
            while(current->right) current = current->right;
        }
        else
        {//parent if came from right child, else keep going up
            while(current->parent && current != current->parent->right)
                current = current->parent;
            current = current->parent;
        }
        return *this;
    }
    NODE& operator*(){assert(current); return *current;}
    NODE* operator->(){assert(current); return current;}
    bool operator==(TreeIterator const& rhs)const
        {return current == rhs.current;}
};

template<typename NODE> struct TreapGeneric
{
    static void rotateHelper(NODE* node, NODE* goingUp, NODE*& movedChild)
    {
        if(movedChild) movedChild->parent = node;
        movedChild = node;
        goingUp->nodeCount = node->nodeCount;
        goingUp->parent = node->parent;
        node->parent = goingUp;
        node->nodeCount = 1 + (node->left ? node->left->nodeCount: 0) +
            (node->right ? node->right->nodeCount: 0);
    }
    static NODE* findMin(NODE* root)
    {
        NODE* node = root;
        if(node) while(node->left) node = node->left;
        return node;
    }
    static NODE* findMax(NODE* root)
    {
        NODE* node = root;
        if(node) while(node->right) node = node->right;
        return node;
    }
    static NODE* nthElement(int n, NODE* root)
    {
        assert(n >= 0 && root && n < root->nodeCount);
        NODE* node = root;
        for(;;)
        {
            unsigned int lc = node->left ? node->left->nodeCount : 0;
            if(n == lc) break;
            if(n < lc) node = node->left;
            else
            {
                n -= lc + 1;
                node = node->right;
            }
        }
        return node;
    }
};

template<typename KEY, typename VALUE, typename COMPARATOR =
    DefaultComparator<KEY> > class Treap
{
    COMPARATOR c;
    struct Node
    {
        KEY key;
        VALUE value;
        Node *left, *right, *parent;
        unsigned int priority, nodeCount;
        Node(KEY const& theKey, VALUE const& theValue): key(theKey),
            value(theValue), left(0), right(0), parent(0),
            priority(GlobalRNG().next()), nodeCount(1){}
    }* root;
    Freelist<Node> f;
    Node* rotateRight(Node* node)
    {
        Node *goingUp = node->left, *&movedChild = goingUp->right;
        node->left = movedChild;
        TreapGeneric<Node>::rotateHelper(node, goingUp, movedChild);
        return goingUp;
    }
    Node* rotateLeft(Node* node)
    {
        Node *goingUp = node->right, *&movedChild = goingUp->left;
        node->right = movedChild;
        TreapGeneric<Node>::rotateHelper(node, goingUp, movedChild);
        return goingUp;
    }
    Node* insertNode(Node* newNode, Node* node)
    {
        if(!node) return newNode;
        bool goLeft = c(newNode->key, node->key);
        Node*& chosenChild = goLeft ? node->left : node->right;
        chosenChild = insertNode(newNode, chosenChild);
        chosenChild->parent = node;
        ++node->nodeCount;
        if(chosenChild->priority < node->priority)
            node = goLeft ? rotateRight(node) : rotateLeft(node);
        return node;
    }
    Node* constructFrom(Node* node)
    {
        Node* tree = 0;
        if(node)
        {
            tree = new(f.allocate())Node(node->key, node->value);
            tree->priority = node->priority;
            tree->nodeCount = node->nodeCount;
            tree->left = constructFrom(node->left);
            if(tree->left) tree->left->parent = tree;
            tree->right = constructFrom(node->right);
            if(tree->right) tree->right->parent = tree;
        }
        return tree;
    }
public:
    typedef Node NodeType;
    unsigned int getSize(){return root ? root->nodeCount : 0;}
    Treap(COMPARATOR const& theC = COMPARATOR()): root(0), c(theC){}
    Treap(Treap const& other): c(other.c)
    {
        root = constructFrom(other.root);
        if(root) root->parent = 0;
    }
    Treap& operator=(Treap const& rhs){return genericAssign(*this, rhs);}

    Node* findNode(KEY const& key)
    {
        Node* node = root;
        while(node && !c.isEqual(key, node->key)) node =
            c(key, node->key) ? node->left : node->right;
        return node;
    }
    VALUE* find(KEY const& key)
    {
        Node* node = findNode(key);
        return node ? &node->value : 0;
    }
//CALC USING ITERATORS? + FIND BOTH AT SAME TIME?
//BRACKET IS MORE COMPLEX HERE?
//JUST COMMENT THAT CAN DO ITERATORS INSTEAD????
    NodeType* predecessor(KEY const& key)
    {
        Node* pred = 0;
        for(Node* node = root; node;)
            if(c(node->key, key))//found pred set member
            {
                pred = node;
                node = node->right;
            }
            else node = node->left;
        return pred;
    }
    NodeType* inclusivePredecessor(KEY const& key)
    {
        Node* node = findNode(key);
        return node ? node : predecessor(key);
    }
    NodeType* successor(KEY const& key)
    {
        Node* succ = 0;
        for(Node* node = root; node;)
            if(c(key, node->key))//found succ set member
            {
                succ = node;
                node = node->left;
            }
            else node = node->right;
        return succ;
    }
    NodeType* inclusiveSuccessor(KEY const& key)
    {
        Node* node = findNode(key);
        return node ? node : successor(key);
    }

    NodeType* insert(KEY const& key, VALUE const& value)
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
    NodeType* findMin(){return TreapGeneric<Node>::findMin(root);}
    NodeType* findMax(){return TreapGeneric<Node>::findMax(root);}
    NodeType* nthElement(int n)
        {return TreapGeneric<Node>::nthElement(n, root);}
    void debugHelper(Node* node)
    {
        if(node)
        {
            DEBUG(node);
            DEBUG(node->key);
            DEBUG(node->left);
            DEBUG(node->right);
            DEBUG(node->parent);
            DEBUG(node->nodeCount);
            debugHelper(node->left);
            debugHelper(node->right);
        }
    }
    void debug()
    {
        debugHelper(root);
    }
    typedef TreeIterator<Node> Iterator;
    Iterator begin(){return Iterator(findMin());}
    Iterator end(){return Iterator(0);}
    Iterator rBegin(){return Iterator(findMax());}
    Iterator rEnd(){return Iterator(0);}
};

template<typename WORD = unsigned long long> class Key2DBuilder
{
    unsigned int n;
    bool firstNotSmaller;
public:
    typedef WORD WORD_TYPE;
    Key2DBuilder(unsigned int theN = numeric_limits<unsigned int>::max(),
        bool theFirstNotSmaller = true): n(theN),
        firstNotSmaller(theFirstNotSmaller){}
    WORD to1D(unsigned int n1, unsigned int n2)const
    {
        assert(n1 < n && n2 < n);
        return firstNotSmaller ? n1 * n + n2 : n2 * n + n1;
    }
    pair<unsigned int, unsigned int> to2D(WORD key)const
    {
        pair<unsigned int, unsigned int> n1n2(key/n, key % n);
        assert(n1n2.first < n && n1n2.second < n);
        if(!firstNotSmaller) swap(n1n2.first, n1n2.second);
        return n1n2;
    }
};

}//end namespace
#endif
