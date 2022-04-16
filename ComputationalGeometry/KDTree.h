#ifndef IGMDK_KDTREE_H
#define IGMDK_KDTREE_H
#include "../Utils/Utils.h"
#include "../Utils/Debug.h"
#include "../Utils/Vector.h"
#include "../Heaps/Heap.h"
#include "../Utils/GCFreeList.h"
#include "../NumericalMethods/Matrix.h"//for eLess
#include <cmath>
namespace igmdk{

template<typename NODE> struct QNode
{
    NODE* n;
    double d;
    bool operator<(QNode const& rhs)const{return d > rhs.d;}
    static double dHeap(Heap<QNode>& heap, int k)
    {
        return heap.getSize() < k ?
            numeric_limits<double>::max() : heap.getMin().d;
    }
};

template<typename KEY, typename VALUE, typename DISTANCE> class VpTree
{
    DISTANCE lowerBound;
    static double bound(double keyDistance, double rLow, double rHigh)
        {return max(0., max(keyDistance - rHigh, rLow - keyDistance));}
    struct Node
    {
        KEY key;
        VALUE value;
        double leftChildDistance, radius;
        Node *left, *right;
        Node(KEY const& theKey, VALUE const& theValue): key(theKey), left(0),
            right(0), value(theValue), leftChildDistance(0), radius(0) {}
        double leftChildBound(double keyDistance)
            {return bound(keyDistance, 0, leftChildDistance);}
        double rightChildBound(double keyDistance)
            {return bound(keyDistance, leftChildDistance, radius);}
    }* root;
    Freelist<Node> f;
    Node* constructFrom(Node* n)
    {
        Node* tree = 0;
        if(n)
        {
            tree = new(f.allocate())Node(n->key, n->value);
            tree->leftChildDistance = n->leftChildDistance;
            tree->radius = n->radius;
            tree->left = constructFrom(n->left);
            tree->right = constructFrom(n->right);
        }
        return tree;
    }
    void distanceQuery(KEY const& key, double radius, Vector<Node*>& result,
        Node* n)const
    {
        if(!n) return;
        double d = lowerBound(n->key, key);
        if(d <= radius) result.append(n);
        if(n->leftChildBound(d) <= radius)//first go left if not pruned
            distanceQuery(key, radius, result, n->left);
        if(n->rightChildBound(d) <= radius)//then right if not pruned
            distanceQuery(key, radius, result, n->right);
    }
    typedef QNode<Node> HEAP_ITEM;
    void kNN(Node* n, KEY const& key, Heap<HEAP_ITEM>& heap, int k)const
    {
        if(!n) return;
        //replace furthest n in heap with the current n if it's closer
        HEAP_ITEM x = {n, lowerBound(key, n->key)};
        if(heap.getSize() < k) heap.insert(x);
        else if(x.d < HEAP_ITEM::dHeap(heap, k)) heap.changeKey(0, x);
        //expand closer child first
        double lb = n->leftChildBound(x.d), rb = n->rightChildBound(x.d);
        Node* l = n->left, *r = n->right;
        if(lb > rb)//go to smaller-lower-bound node first to reduce the chance
        {//of going to the other one by placing closer nodes on the heap
            swap(lb, rb);
            swap(l, r);
        }
        if(lb <= HEAP_ITEM::dHeap(heap, k)) kNN(l, key, heap, k);
        if(rb <= HEAP_ITEM::dHeap(heap, k)) kNN(r, key, heap, k);
    }
public:
    typedef DISTANCE DISTANCE_TYPE;//update doc!
    DISTANCE const& getDistance(){return lowerBound;}
    typedef Node NodeType;
    bool isEmpty()const{return !root;}
    VpTree(DISTANCE const& theDistance = DISTANCE()): root(0),
        lowerBound(theDistance) {}
    VpTree(VpTree const& rhs): lowerBound(rhs.lowerBound)
        {root = constructFrom(rhs.root);}
    VpTree& operator=(VpTree const& rhs){return genericAssign(*this, rhs);}

    Vector<NodeType*> distanceQuery(KEY const& key, double radius)const
    {
        Vector<NodeType*> result;
        distanceQuery(key, radius, result, root);
        return result;
    }

    VALUE* find(KEY const& key)const
    {
        Node* n = root;
        while(n && key != n->key) n = !isELess(n->leftChildDistance,
            lowerBound(key, n->key)) ? n->left : n->right;
        return n ? &n->value : 0;
    }
    void insert(KEY const& key, VALUE const& value)
    {
        Node **pointer = &root, *n;
        while((n = *pointer) && key != n->key)
        {
            double d = lowerBound(key, n->key);
            n->radius = max(n->radius, d);
            if(!n->left) n->leftChildDistance = d;//will make left child
            pointer = &(!isELess(n->leftChildDistance, d)? n->left : n->right);
        }
        if(n) n->value = value;//equality--assign new value
        else *pointer = new(f.allocate())Node(key, value);
    }

    Vector<NodeType*> kNN(KEY const& key, int k)const
    {
        Heap<HEAP_ITEM> heap;
        kNN(root, key, heap, k);
        Vector<Node*> result;//heap-sort found nodes in by distance
        while(!heap.isEmpty()) result.append(heap.deleteMin().n);
        result.reverse();
        return result;
    }
    NodeType* nearestNeighbor(KEY const& key)const
    {
        assert(!isEmpty());
        return kNN(key, 1)[0];
    }
};

template<typename KEY, typename VALUE,
    typename INDEXED_COMPARATOR = LexicographicComparator<KEY> > class KDTree
{
    INDEXED_COMPARATOR c;
    struct Node
    {
        KEY key;
        VALUE value;
        Node *left, *right;
        Node(KEY const& theKey, VALUE const& theValue): key(theKey),
            value(theValue), left(0), right(0) {}
    }* root;
    Freelist<Node> f;
    int D;
    Node* constructFrom(Node* n)
    {
        Node* tree = 0;
        if(n)
        {
            tree = new(f.allocate())Node(n->key, n->value);
            tree->left = constructFrom(n->left);
            tree->right = constructFrom(n->right);
        }
        return tree;
    }
    Node** findPointer(KEY const& key, Node*& parent)const
    {
        Node* n, **pointer = (Node**)&root;//cast for const
        parent = 0;
        for(int i = 0; (n = *pointer) && !c.isEqual(key, n->key);
            i = (i + 1) % D)
        {
            parent = n;
            pointer = &(c(key, n->key, i) ?
                n->left : n->right);
        }
        return pointer;
    }
    void rangeQuery(KEY const& l, KEY const& u, Vector<bool> const& dimensions,
        Vector<Node*>& result, Node* n, int i)const
    {
        if(!n) return;
        bool inRange = true;//check if current node in range
        for(int j = 0; j < D; ++j)
            if(dimensions[j] && (c(n->key, l, j) ||
                c(u, n->key, i))) inRange = false;
        if(inRange) result.append(n);
        int j = (i + 1) % D;//only check range for the wanted dimensions
        if(!(dimensions[i] && c(n->key, l, i)))
            rangeQuery(l, u, dimensions, result, n->left, j);
        if(!(dimensions[i] && c(u, n->key, i)))
            rangeQuery(l, u, dimensions, result, n->right, j);
    }
    template<typename DISTANCE> void distanceQuery(KEY const& x,
        double partialRadius, Vector<Node*>& result, Node* n, int i,
        DISTANCE const& distanceIncremental, KEY& partial,
        double partialDistance)const
    {//first try to prune subtree
        if(!n || partialDistance > partialRadius) return;
        if(distanceIncremental(n->key, x) <= partialRadius)
            result.append(n);
        i = (i + 1) % D;
        Node* nodes[] = {n->left, n->right};
        for(int j = 0; j < 2; ++j)
        {//apply partial to right subtree if x[i] on the left side of n and
            //to left if on the right; equality not a problem
            bool applyPartial = c(x, n->key, i) == (j == 1);
            double dDelta = 0;
            if(applyPartial)
            {
                dDelta = distanceIncremental(x, n->key, i) -
                    distanceIncremental(x, partial, i);
                swap(partial[i], n->key[i]);//use n as temp storage
            }
            distanceQuery(x, partialRadius, result, nodes[j], i,
                distanceIncremental, partial, partialDistance + dDelta);
            if(applyPartial) swap(partial[i], n->key[i]);
        }
    }
    typedef QNode<Node> HEAP_ITEM;
    template<typename DISTANCE> void kNN(Node* n, KEY const& key,
        Heap<HEAP_ITEM>& heap, int k, int i, KEY& partial,
        double partialDistance, DISTANCE const& distance)const
    {
        double best = HEAP_ITEM::dHeap(heap, k);
        if(n && partialDistance < best)
        {//update partial distance
            double newPartialDistance = distance(key, n->key, i) -
                distance(key, partial, i);
            if(heap.getSize() < k)
            {
                HEAP_ITEM x = {n, distance(key, n->key)};
                heap.insert(x);
            }
            //use new partial distance to check for a cut again
            else if(newPartialDistance < best)
            {//incremental calculate-compare
                double d = distance(best, key, n->key);
                if(d < best)
                {
                    HEAP_ITEM x = {n, d};
                    heap.changeKey(0, x);
                }
            }
            int j = (i + 1) % D;
            //swap children for best order
            Node *l = n->left, *r = n->right;
            if(!c(key, n->key, i)) swap(l, r);
            kNN(l, key, heap, k, j, partial, partialDistance, distance);
            //set partial component to the n component, use the n
            //as temporary storage
            swap(partial[i], n->key[i]);
            kNN(r, key, heap, k, j, partial, newPartialDistance, distance);
            swap(partial[i], n->key[i]);
        }
    }
public:
    typedef Node NodeType;
    bool isEmpty()const{return !root;}
    KDTree(int theD, INDEXED_COMPARATOR const& theC = INDEXED_COMPARATOR()):
        root(0), c(theC), D(theD) {}
    KDTree(KDTree const& rhs): c(rhs.c){root = constructFrom(rhs.root);}
    KDTree& operator=(KDTree const& rhs){return genericAssign(*this, rhs);}
    VALUE* find(KEY const& key)const
    {
        Node *n = *findPointer(key, n);
        return n ? &n->value : 0;
    }
    void insert(KEY const& key, VALUE const& value)
    {
        Node *dummy, **pointer = findPointer(key, dummy);
        if(*pointer) (*pointer)->value = value;
        else *pointer = new(f.allocate())Node(key, value);
    }
    Vector<NodeType*> rangeQuery(KEY const& l, KEY const& u,
        Vector<bool> const& dimensions)const
    {
        Vector<Node*> result;
        rangeQuery(l, u, dimensions, result, root, 0);
        return result;
    }
    template<typename DISTANCE> Vector<NodeType*> distanceQuery(KEY const& x,
        double partialRadius, DISTANCE const& distanceIncremental)const
    {
        Vector<Node*> result;
        KEY partial = x;
        distanceQuery(x, partialRadius, result, root, 0, distanceIncremental,
            partial, 0);
        return result;
    }
    template<typename DISTANCE> Vector<NodeType*> kNN(KEY const& key, int k,
        DISTANCE const& distance)const
    {
        Heap<HEAP_ITEM> heap;
        KEY partial = key;
        kNN(root, key, heap, k, 0, partial, 0, distance);
        Vector<Node*> result;//heap-sort found nodes in by distance
        while(!heap.isEmpty()) result.append(heap.deleteMin().n);
        result.reverse();
        return result;
    }
    template<typename DISTANCE> NodeType* nearestNeighbor(KEY const& key,
        DISTANCE const& distance)const
    {
        assert(!isEmpty());
        Node* parent, *result = *findPointer(key, parent);
        if(result) return result;//found equal-value node, d = 0
        Heap<HEAP_ITEM> heap;//put parent on heap
        HEAP_ITEM x = {parent, distance(key, parent->key)};
        heap.insert(x);
        KEY partial = key;
        kNN(root, key, heap, 1, 0, partial, 0, distance);
        return heap.getMin().n;
    }
};

}//end namespace
#endif
