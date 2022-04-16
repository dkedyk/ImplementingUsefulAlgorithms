#ifndef IGMDK_CHAINING_HASH_TABLE_H
#define IGMDK_CHAINING_HASH_TABLE_H
#include "HashFunction.h"
#include "../Utils/GCFreeList.h"
namespace igmdk{

template<typename KEY, typename VALUE, typename HASHER = EHash<BUHash>,
    typename COMPARATOR = DefaultComparator<KEY> > class ChainingHashTable
{
    int capacity, size;//capacity must be before h
    struct Node
    {
        KEY key;
        VALUE value;
        Node* next;
        Node(KEY const& theKey, VALUE const& theValue): key(theKey),
            value(theValue), next(0) {}
    }** table;
    Freelist<Node> f;
    HASHER h;
    COMPARATOR c;
    enum{MIN_CAPACITY = 8};//for efficiency require at least size 8
    void allocateTable()
    {
        h = HASHER(capacity);
        table = new Node*[capacity];
        for(int i = 0; i < capacity; ++i) table[i] = 0;
    }
    void resize()
    {
        int oldCapacity = capacity;
        Node** oldTable = table;
        capacity = nextPowerOfTwo(size * 2);
        allocateTable();
        for(int i = 0; i < oldCapacity; ++i)
            for(Node* j = oldTable[i], *tail; j; j = tail)
            {
                tail = j->next;
                j->next = 0;
                *findPointer(j->key) = j;//insert node
            }
        delete[] oldTable;
    }
    Node** findPointer(KEY const& key)
    {//for code reuse get pointer to node pointer
        Node** pointer = &table[h(key)];
        for(;*pointer && !c.isEqual((*pointer)->key, key);
            pointer = &(*pointer)->next);
        return pointer;//if not found return pointer to next of last node
    }
public:
    typedef Node NodeType;
    int getSize(){return size;}
    ChainingHashTable(int initialCapacity = 8, COMPARATOR const&
        theC = COMPARATOR()): capacity(nextPowerOfTwo(max<int>(initialCapacity,
        MIN_CAPACITY))), c(theC), h(capacity), size(0) {allocateTable();}
    ChainingHashTable(ChainingHashTable const& rhs): capacity(rhs.capacity),
        size(rhs.size), h(rhs.h), table(new Node*[capacity]), c(rhs.c)
    {//copy just mirrors the source, without trying to compact
        for(int i = 0; i < capacity; ++i)
        {
            table[i] = 0;
            Node** target = &table[i];
            for(Node* j = rhs.table[i]; j; j = j->next)
            {
                *target = new(f.allocate())Node(*j);
                target = &(*target)->next;
            }
        }
    }
    ChainingHashTable& operator=(ChainingHashTable const& rhs)
        {return genericAssign(*this, rhs);}
    ~ChainingHashTable(){delete[] table;}
    Node* insert(KEY const& key, VALUE const& value)
    {
        Node** pointer = findPointer(key);
        if(*pointer)
        {//already exists, just update value
            (*pointer)->value = value;
            return *pointer;
        }
        else
        {
            Node* node = *pointer = new(f.allocate())Node(key, value);
            if(++size >= capacity) resize();//not > where will have x4 size
            return node;
        }
    }

    //chaining has node persistence, so allow pointer return
    Node* findNode(KEY const& key){return *findPointer(key);}
    VALUE* find(KEY const& key)
    {
        Node* next = findNode(key);
        return next ? &next->value : 0;
    }
    void remove(KEY const& key)
    {
        Node** pointer = findPointer(key);
        Node* i = *pointer;
        if(i)
        {//found
            *pointer = i->next;
            f.remove(i);
            if(--size < capacity * 0.1 && size * 2 >= MIN_CAPACITY) resize();
        }
    }
    class Iterator
    {
        int i;//current cell index
        Node* node;//node in cell
        ChainingHashTable& t;
        friend ChainingHashTable;
        void advanceCell()//if at null node and not at end, try next cell
            {if(!node) while(i + 1 < t.capacity && !(node = t.table[++i]));}
        Iterator(ChainingHashTable& theHashTable, int theI = -1): i(theI),
            node(0), t(theHashTable) {advanceCell();}
    public:
        Iterator& operator++()
        {
            assert(node);
            node = node->next;
            advanceCell();
            return *this;
        }
        NodeType& operator*(){assert(node); return *node;}
        NodeType* operator->(){assert(node); return node;}
        bool operator==(Iterator const& rhs)const{return node == rhs.node;}
    };
    Iterator begin(){return Iterator(*this);}
    Iterator end()
    {
        Iterator result(*this, capacity);
        return result;
    }
};

}//end namespace
#endif
