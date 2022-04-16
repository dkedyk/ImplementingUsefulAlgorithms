#ifndef IGMDK_CUCKOO_HASH_TABLE_H
#define IGMDK_CUCKOO_HASH_TABLE_H

#include "../HashTable/HashFunction.h"
#include "../Utils/Utils.h"

namespace igmdk{

template<typename KEY, typename VALUE, typename HASHER = EHash<BUHash> >
class CuckooHashTable
{
    enum{NOT_FOUND = -1};
	int capacity, size;
	typedef KVPair<KEY, VALUE> Node;
	Node* table;
	bool* isOccupied;
	HASHER h1, h2;
	int hash2(KEY const& key, int hash1)
	{
	    //to remember which hash function was used to hash the key
	    //the trick is to use h1 = hash1(key) and
	    //h2 = (hash2(key) - h1) & capacity because this way
	    //given cell, the other cell is (hash2(key) - cell) & capacity
	    int result = h2(key) - hash1;
	    if(result < 0) result += capacity;
	    return result;
	}
	int findNode(KEY const& key)
	{
		for(int i = 0, cell = h1(key); i < 2; ++i)
		{
			if(isOccupied[cell] && table[cell].key == key) return cell;
			if(i == 0) cell = hash2(key, cell);
        }
		return NOT_FOUND;
	}
	void allocateTable()
	{
	    table = rawMemory<Node>(capacity);
	    isOccupied = new bool[capacity];
		setGoodState();
    }
    void setGoodState()
	{
        h1 = HASHER(capacity);
        h2 = HASHER(capacity);
        size = 0;
        for(int i = 0; i < capacity; ++i) isOccupied[i] = false;
	}
    static void cleanUp(Node* theTable, int theCapacity, bool* isOccupied)
	{
		for(int i = 0; i < theCapacity; ++i)
            if(isOccupied[i]) theTable[i].~Node();
		rawDelete(theTable);
		delete[] isOccupied;
	}
    bool insertHelper(KEY const& key, VALUE const& value, bool duringResize)
	{
        //load factor must be < 0.5, which is phase transition for a universal
	    //hash unction, for kickout phase to take
	    if(size > capacity * 0.4) resize(true);
	    //can try to first check if both locations are empty rather then just
	    //the first one. This gives about 17% insertion speedup at expense of
	    //more code, so not worth is since Cuckoo Hashing is not the
	    //method of choice
        for(Node node(key, value);;resize(false))
	    {
	        int cell = h1(node.key);
	        //best max limit choice is not clear, but 50 works fine
            for(int limit = 0; limit < 50; ++limit)
            {
                if(!isOccupied[cell])
                {
                    ++size;
                    new(&table[cell])Node(node);
                    isOccupied[cell] = true;
                    return true;
                }
                if(limit < 2 && table[cell].key == node.key)
                {
                    table[cell].value = node.value;
                    return true;
                }
                swap(node, table[cell]);
                cell = hash2(node.key, cell);
            }
            if(duringResize) return false;
	    }
	}
	void resize(bool increaseSize)
	{
	    Node* oldTable = table;
		int oldCapacity = capacity;
		bool* oldIsOccupied = isOccupied;
		if(increaseSize) capacity = nextPowerOfTwo(max(4 * size, 8));
		allocateTable();
        for(int i = 0; i < oldCapacity; ++i)
        {
            if(oldIsOccupied[i] && !insertHelper(oldTable[i].key,
                oldTable[i].value, true))
            {
                setGoodState();
                i = -1;
            }
        }
		cleanUp(oldTable, oldCapacity, oldIsOccupied);
	}
public:
    CuckooHashTable(int initialSize = 8): capacity(nextPowerOfTwo(max(
        initialSize, 8))), h1(capacity), h2(capacity) {allocateTable();}
    VALUE* find(KEY const& key)
	{
		int result = findNode(key);
		return result == NOT_FOUND ? 0 : &table[result].value;
	}
	~CuckooHashTable(){cleanUp(table, capacity, isOccupied);}
	void insert(KEY const& key, VALUE const& value)
        {insertHelper(key, value, false);}
	void remove(KEY const& key)
	{
		int result = findNode(key);
		if(result != NOT_FOUND)
		{
			table[result].~Node();
			isOccupied[result] = false;
			if(--size < capacity * 0.1) resize(true);
		}
	}
	typedef Node NodeType;
	int getSize(){return size;}
    template<typename ACTION> void forEach(ACTION& action)
	{
        for(int i = 0; i < capacity; ++i)
            if(table[i].isOccupied) action(&table[i]);
	}
	class Iterator
	{
	    int i;
	    CuckooHashTable& hashTable;
	    void advance()
            {while(i < hashTable.capacity && !hashTable.isOccupied[i])++i;}
    public:
        Iterator(CuckooHashTable& theHashTable): i(0), hashTable(theHashTable)
            {advance();}
        bool hasNext()
            {return i < hashTable.capacity && hashTable.isOccupied[i];}
        NodeType* next()
        {
            Node* result = hasNext() ? &hashTable.table[i++] : 0;
            advance();
            return result;
        }
	};
};

}
#endif
