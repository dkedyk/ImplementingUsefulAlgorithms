#ifndef IGMDK_LINEAR_PROBING_HASH_TABLE_H
#define IGMDK_LINEAR_PROBING_HASH_TABLE_H
#include <new>
#include "HashFunction.h"
namespace igmdk{

template<typename KEY, typename VALUE, typename HASHER = EHash<BUHash>,
typename COMPARATOR = DefaultComparator<KEY> > class LinearProbingHashTable
{
    int capacity, size;//capacity must be before h
    typedef KVPair<KEY, VALUE> Node;
    Node* table;
    bool* isOccupied;
    HASHER h;
    COMPARATOR c;
    enum{MIN_CAPACITY = 8};//for efficiency require at least size 8
    void allocateTable()
    {//create an unoccupied table of size capacity
        h = HASHER(capacity);
        size = 0;
        table = rawMemory<Node>(capacity);
        isOccupied = new bool[capacity];
        for(int i = 0; i < capacity; ++i) isOccupied[i] = false;
    }//helper to remove an unused table
    static void cleanUp(Node* theTable, int theCapacity, bool* isOccupied)
    {//destruct the occupied nodes, and deallocate the arrays
        for(int i = 0; i < theCapacity; ++i)
            if(isOccupied[i]) theTable[i].~Node();
        rawDelete(theTable);
        delete[] isOccupied;
    }
    void destroy(int cell)
    {
        table[cell].~Node();
        isOccupied[cell] = false;
        --size;
    }
    void resize()
    {
        int oldCapacity = capacity;
        Node* oldTable = table;
        bool* oldIsOccupied = isOccupied;
        capacity = nextPowerOfTwo(size * 2);
        allocateTable();
        for(int i = 0; i < oldCapacity; ++i)//reinsert
            if(oldIsOccupied[i]) insert(oldTable[i].key, oldTable[i].value);
        cleanUp(oldTable, oldCapacity, oldIsOccupied);//remove old table
    }
    int findNode(KEY const& key)
    {//find the cell where the key would belong if inserted
        int cell = h(key);
        for(; isOccupied[cell] && !c.isEqual(key, table[cell].key);
            cell = (cell + 1) % capacity);
        return cell;
    }
public:
    typedef Node NodeType;
    int getSize(){return size;}
    LinearProbingHashTable(int initialCapacity = 8, COMPARATOR const& theC =
        COMPARATOR()): capacity(nextPowerOfTwo(max<int>(initialCapacity,
        MIN_CAPACITY))), c(theC), h(capacity) {allocateTable();}
    LinearProbingHashTable(LinearProbingHashTable const& rhs):
        capacity(rhs.capacity), h(rhs.h), size(rhs.size), c(rhs.c),
        isOccupied(new bool[capacity]), table(rawMemory<Node>(capacity))
    {//copy just mirrors the source, without trying to compact
        for(int i = 0; i < capacity; ++i)
            if(isOccupied[i] = rhs.isOccupied[i]) table[i] = rhs.table[i];
    }
    LinearProbingHashTable& operator=(LinearProbingHashTable const& rhs)
        {return genericAssign(*this, rhs);}
    ~LinearProbingHashTable(){cleanUp(table, capacity, isOccupied);}

    VALUE* find(KEY const& key)
    {
        int cell = findNode(key);
        return isOccupied[cell] ? &table[cell].value : 0;
    }
    void insert(KEY const& key, VALUE const& value)
    {
        int cell = findNode(key);
        if(isOccupied[cell]) table[cell].value = value;//update
        else
        {//insert
            new(&table[cell])Node(key, value);
            isOccupied[cell] = true;
            if(++size > capacity * 0.8) resize();//resize if reach a
        }
    }
    void remove(KEY const& key)
    {
        int cell = findNode(key);
        if(isOccupied[cell])
        {//reinsert subsequent nodes in the found value's chain
            destroy(cell);//remove item
            if(size < capacity * 0.1 && size * 2 >= MIN_CAPACITY) resize();
            else//reinsert chain
                while(isOccupied[cell = (cell + 1) % capacity])
                {
                    Node temp = table[cell];
                    destroy(cell);//destroy item
                    insert(temp.key, temp.value);//reinsert it
                }
        }
    }//below optimized algorithm has bug somewhere - will debug later
    /*void remove(KEY const& key)
    {
        int cell = findNode(key);
        if(isOccupied[cell])
        {//reinsert subsequent nodes in the found value's chain
            destroy(cell);//no need to compact if will resize
            if(size < capacity * 0.1 && size * 2 >= MIN_CAPACITY) resize();
            else
            {//compact chain
                int last = cell;//first find last cell
                while(isOccupied[(last + 1) % capacity])
                    last = (last + 1) % capacity;
                while(cell < last)//search from last
                {
                    int next = last;
                    for(; next > cell; --next)
                    {
                        if(findNode(table[next].key) == cell)
                        {//found cell to use as fill
                            new(&table[cell])Node(table[next]);
                            isOccupied[cell] = true;
                            ++size;
                            destroy(next);
                            cell = next;//now compact from next cell
                        }
                    }
                    if(next >= cell) break;//nothing to compact
                }
            }
        }
    }*/
    class Iterator
    {
        int i;//current cell index
        LinearProbingHashTable& t;
        friend LinearProbingHashTable;
        void advance(){while(i < t.capacity && !t.isOccupied[i]) ++i;}
        Iterator(LinearProbingHashTable& theHashTable, int theI = 0): i(theI),
            t(theHashTable) {advance();}
    public:
        Iterator& operator++()
        {
            ++i;
            advance();
            return *this;
        }
        NodeType& operator*()const{assert(i < t.capacity); return t.table[i];}
        NodeType* operator->()const{assert(i < t.capacity);return &t.table[i];}
        bool operator==(Iterator const& rhs)const{return i == rhs.i;}
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
