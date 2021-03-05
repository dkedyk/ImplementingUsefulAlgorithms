#ifndef IGMDK_SKIPLIST_H
#define IGMDK_SKIPLIST_H
#include "../RandomNumberGeneration/Random.h"
#include "../Utils/GCFreeList.h"
namespace igmdk{

template<typename KEY, typename VALUE, typename COMPARATOR =
    DefaultComparator<VALUE> > class SkipList
{
    COMPARATOR c;
    enum{MAX_HEIGHT = 32};
    struct Node
    {
        KEY key;
        VALUE value;
        Node** tower;
        Node(KEY const& theKey, VALUE const& theValue, int height):
            key(theKey), value(theValue), tower(new Node*[height]) {}
        ~Node(){delete[] tower;}
    }* head[MAX_HEIGHT];
    Freelist<Node> f;
    int currentLevel;
public:
    SkipList(COMPARATOR const& theC = COMPARATOR()): currentLevel(0), c(theC)
        {for(int i = 0; i < MAX_HEIGHT; ++i) head[i] = 0;}
    SkipList(SkipList const& rhs): currentLevel(0), c(rhs.c)
    {//order of items with nonunique keys in copy is reversed
        for(int i = 0; i < MAX_HEIGHT; ++i) head[i] = 0;
        for(Node* node = rhs.head[0]; node; node = node->tower[0])
            insert(node->key, node->value, false);
    }
    SkipList& operator=(SkipList const&rhs){return genericAssign(*this,rhs);}

    class Iterator
    {
        Node* current;
        Iterator(Node* node): current(node){}
        friend SkipList;
    public:
        Iterator& operator++()
        {
            assert(current);
            current = current->tower[0];
            return *this;
        }
        Node& operator*()
        {
            assert(current);
            return *current;
        }
        Node* operator->()
        {
            assert(current);
            return current;
        }
        bool operator==(Iterator const& rhs)const
            {return current == rhs.current;}
    };
    Iterator begin(){return Iterator(findMin());}
    Iterator end(){return Iterator(0);}

    Iterator predecessor(KEY const& key)
    {
        Node **tower = head, *pred = 0;
        for(int level = currentLevel; level >= 0; --level)
            for(Node* node; (node = tower[level]) && c(node->key, key);
                tower = node->tower) pred = node;
        return Iterator(pred);
    }
    Iterator inclusiveSuccessor(KEY const& key)
    {//next(pred) = inc succ
        Iterator pred = predecessor(key);
        assert(pred == end() || c(pred->key, key));
        return pred == end() ? begin() : Iterator(pred->tower[0]);
    }
    Iterator findNode(KEY const& key)
    {//general pattern – return a pointer for reuse in other operations
        Iterator node = inclusiveSuccessor(key);//match if not larger than key
        assert(node == end() || !c(node->key, key));
        return node == end() || c(key, node->key) ? end() : node;
    }
    VALUE* find(KEY const& key)
    {
        Iterator result = findNode(key);
        return result == end() ? 0 : &result->value;
    }
    Iterator successor(KEY const& key)
    {//is inclusiveSuccessor with nonunique keys, else loop over equal keys
        Node* node = inclusiveSuccessor(key)->current;
        while(node && !c(node->key, key)) node = node->tower[0];
        return Iterator(node);
    }

    Iterator insert(KEY const& key, VALUE const& value, bool unique = true)
    {
        if(unique)
        {//for unique keys check if one already exists
            Iterator result = findNode(key);
            if(result != end())
            {
                result->value = value;
                return result;
            }
        }//level = height - 1
        int newLevel = min<int>(MAX_HEIGHT, GlobalRNG().geometric(0.632)) - 1;
        Node* newNode = new(f.allocate())Node(key, value, newLevel + 1);
        if(currentLevel < newLevel) currentLevel = newLevel;
        Node** tower = head;
        for(int level = currentLevel; level >= 0; --level)
        {
            for(Node* node; (node = tower[level]) &&
                c(node->key, key); tower = node->tower);
            if(level <= newLevel)
            {//relink pointers
                newNode->tower[level] = tower[level];
                tower[level] = newNode;
            }
        }
        return Iterator(newNode);
    }
    void remove(KEY const& key)
    {//with nonunique items will remove first found
        Node **prevTower = head, *result = 0;
        for(int level = currentLevel; level >= 0; --level)
        {//go down if node->key < key (when result ==0), else keep moving right
            for(Node* node; (node = prevTower[level]) && !c(key, node->key);
                prevTower = node->tower)
                //found if hit remembered node or nothing remembered and ==
                if(node == result || (!result && !c(node->key, key)))
                {//unlink the node from current level
                    prevTower[level] = node->tower[level];
                    node->tower[level] = 0;
                    if(!head[currentLevel])--currentLevel;//if removed top node
                    result = node;//remember node
                    break;//go down
                }
        }
        if(result) f.remove(result);
    }
    Iterator findMin(){return Iterator(head[0]);}
    Iterator findMax()
    {
        Node *result = 0, **tower = head;
        for(int level = currentLevel; level >= 0; --level)
            for(Node* node; node = tower[level]; tower = node->tower)
                result = node;
        return Iterator(result);
    }
};

}//end namespace
#endif
