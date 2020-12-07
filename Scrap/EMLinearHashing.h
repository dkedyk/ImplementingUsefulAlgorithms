#ifndef IGMDK_EM_LINEAR_HASH_TABLE_H
#define IGMDK_EM_LINEAR_HASH_TABLE_H

#include "../ExternalMemoryAlgorithms/EMFreelist.h"
#include "../Utils/Utils.h"
#include "../Utils/Debug.h"
#include "../HashTable/HashFunction.h"
using namespace std;
namespace igmdk{

template<typename KEY, typename ITEM, typename HASH_FUNCTION = EHash<BUHash> >
struct EMLinearHashTable
{
    typedef KVPair<KEY, ITEM> Item;
    enum{BLOCK_SIZE = 2048, NULL_IO_POINTER = -1, N = (BLOCK_SIZE - 16)/sizeof(Item)};
    struct Cell
    {
		long long next;
		Item item[N];
		int size;
		Cell():	next(NULL_IO_POINTER), size(0){}
    };
    int counts[30];

	EMFreelist<Cell> freelist;
	EMVector<Cell> table;
	long long split, bitSize, n, ioFindCount, listCount, maxLength;
	HASH_FUNCTION hashFunction;
	long long sizeBit(){return 1ull << bitSize;}
	long long hashValue(KEY const& key)
	{
	    long long hash = hashFunction(key), size = sizeBit();
	    //DEBUG(key);
	    //DEBUG(hash);
	    if(hash % size < split) size *= 2;
	    return hash % size;
	}
	EMLinearHashTable():freelist("EMLHFreelist.igmdk"), table("EMLHTable.igmdk"), bitSize(0), split(0), n(0),ioFindCount(0),listCount(0),maxLength(0),
        hashFunction(1ull << 30)
        {
            for(int i = 0; i < 30; ++i) counts[i] = 0;
            table.append(Cell());}
	ITEM find(KEY const& key, bool& status)
	{
	    long long length = 0;
		status = true;
		//DEBUG(hashValue(key));
		for(Cell cell = table[hashValue(key)];;cell = freelist[cell.next])
		{
		    ++length;
		    ++ioFindCount;
		    for(int i = 0; i < cell.size; ++i)
            {
                if(cell.item[i].key == key)
                {
                    maxLength = max(maxLength, length);
                    ++counts[length];
                    return cell.item[i].value;
                }
            }
            if(cell.next == NULL_IO_POINTER) break;

		}
		status = false;
	}

	void saveItem(long long index, long long hash, Cell const& cell)
	{
	    if(index == NULL_IO_POINTER) table[hash] = cell;
        else freelist.set(cell, index);
	}
	void splitCell()
	{
	    long long oldSplit = split, pageOld = NULL_IO_POINTER,
            pageNew = NULL_IO_POINTER, pageSpace = NULL_IO_POINTER;
        Cell cell = table[oldSplit], result, spaceCell;
        table.append(result);
        if(++split >= sizeBit())
        {
            split = 0;
            ++bitSize;
        }
        for(;;cell = freelist[pageOld = cell.next])
		{
		    for(int i = 0; i < cell.size; ++i)
            {
                if(hashValue(cell.item[i].key) != oldSplit)
                {
                    if(result.size >= N)
                    {
                        ++listCount;
                        result.next = freelist.allocate();
                        saveItem(pageNew, table.getSize() - 1, result);
                        result = freelist[pageNew = result.next];
                    }
                    result.item[result.size++] = cell.item[i];
                    cell.item[i--] = cell.item[--cell.size];
                }
            }
            saveItem(pageOld, oldSplit, cell);
            if(cell.next == NULL_IO_POINTER) break;
		}
		saveItem(pageNew, table.getSize() - 1, result);

        /*
        //optional optimization - compact the lists when done, this gives
        //maybe 3% speedup
        cell = spaceCell = table[oldSplit];
        pageSpace = pageOld = NULL_IO_POINTER;
		for(;;)
		{
            if(pageSpace == pageOld || cell.size == 0)
            {
                saveItem(pageOld, oldSplit, cell);
                if(cell.next == NULL_IO_POINTER) break;
                cell = freelist[pageOld = cell.next];
            }
            else if(spaceCell.size >= N)
            {
                saveItem(pageSpace, oldSplit, spaceCell);
                assert(spaceCell.next != NULL_IO_POINTER);
                spaceCell = freelist[pageSpace = spaceCell.next];
            }
            else spaceCell.item[spaceCell.size++] = cell.item[--cell.size];
		}
        saveItem(pageSpace, oldSplit, spaceCell);*/
	}
	void insert(KEY const& key, ITEM const& item)
	{
        long long hash = hashValue(key), index = NULL_IO_POINTER;
		for(Cell cell = table[hash];;cell = freelist[index = cell.next])
		{
		    if(cell.size < N)
		    {
		        cell.item[cell.size++] = Item(key, item);
		        saveItem(index, hash, cell);
		        break;
		    }
            if(cell.next == NULL_IO_POINTER)
            {
                ++listCount;
                cell.next = freelist.allocate();
                saveItem(index, hash, cell);
            }
		}
		if(++n > sizeBit() * N) splitCell();//faster but resizes every every
		//insertion and then after none, alternating, the latter is slower
		//but resizes uniformly, this may be better because hash function has
		//a bias
		//if(++n > table.getSize() * N) splitCell();
	}
	long long getSize(){return n;}
	void remove(KEY const& key)
	{
        /*Find the item in a block and
        swap it with the last item in the block, decrease size
        if this is the last item, delete the block, replace the
        pointer to it with its next pointer. Can also join last
        unjoined blocks and merge the buckets to shrink the file*/

        //Will take the lazy approach here and not shrink the cells
        //or merge blocks
		for(Cell cell = table[hashValue(key)];;cell = freelist[cell.next])
		{
		    for(int i = 0; i < cell.size; ++i)
            {
                if(cell.item[i].key == key)
                {
                    cell.item[i] = cell.item[--cell.size];
                    --n;
                    return;
                }
            }
            if(cell.next == NULL_IO_POINTER) break;
		}
	}
};


template<typename KEY, typename ITEM, typename HASHER = EHash<BUHash> >
class EMLinearProbingHashTable
{
	int capacity, size;
	struct Node
    {
        KEY key;
        ITEM item;
        bool isOccupied;
        Node():isOccupied(false){}
        Node(KEY const& theKey, ITEM const& theITEM)
        :	key(theKey), item(theITEM), isOccupied(true)	{}
    };
	EMVector<Node> *table;
	HASHER hashFunction;
	void allocateTable(int requestedSize)
	{
        capacity = nextPowerOfTwo(max(requestedSize, 1));
		hashFunction = HASHER(capacity);
		size = 0;
		//inefficient but good-for-testing allocation strategy
		string newFileName = string("EMLPTHTable") + to_string(capacity) + ".igmdk";
		File::remove(newFileName.c_str());
		table = new EMVector<Node>(newFileName);
		for(int i = 0; i < capacity; ++i) table->append(Node());
	}
	void resize()
	{
		int oldCapacity = capacity;
		EMVector<Node>* oldTable = table;
		allocateTable(4*size);
		for(int i = 0; i < oldCapacity; ++i)
		{
		    Node oldItem = (*oldTable)[i];
		    if(oldItem.isOccupied) insert(oldItem.key, oldItem.item);
		}
		delete oldTable;
	}
	int findNode(KEY const& key)
	{
	    int cell = hashFunction(key);
		for(;; cell = (cell + 1) % capacity)
		{
		    Node result = (*table)[cell];
		    if(!result.isOccupied || key == result.key) break;
		}
		return cell;
	}
public:
	EMLinearProbingHashTable(int initialCapacity = 8):hashFunction(8){allocateTable(initialCapacity);}
    ITEM find(KEY const& key, bool& status)
	{
	    Node result = (*table)[findNode(key)];
	    status = result.isOccupied;
		return result.item;
	}
	void insert(KEY const& key, ITEM const& item)
	{
	    (*table)[findNode(key)] = Node(key, item);
        if(++size * 5 > 4 * capacity) resize();
	}
	void remove(KEY const& key)
	{
        for(int cell = findNode(key);;cell = (cell + 1) % capacity)
        {
            Node result = (*table)[cell];
            if(!result.isOccupied) break;
            --size;
            result.isOccupied = false;
            (*table)[cell] = result;
            //reinsert subsequent nodes in the found item's chain
            if(key != result.key) insert(result.key, result.item);
        }
        //no shrinking
	}
	~EMLinearProbingHashTable(){delete table;}
};

}
#endif
