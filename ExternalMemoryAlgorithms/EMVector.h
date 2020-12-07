#ifndef IGMDK_EMVECTOR_H
#define IGMDK_EMVECTOR_H
#include "File.h"
#include "../Utils/Vector.h"
#include "../Utils/Stack.h"
#include "../Sorting/Sort.h"
#include "../Heaps/Heap.h"
#include "../Utils/Queue.h"
#include <cmath>
using namespace std;

namespace igmdk{

template<typename POD> struct CastSerializer//unportable - only for convenience
{//uncomment when this is supported in a few years
    //CastSerializer(){assert(is_trivially_copyable<POD>::value);}
    constexpr static int byteSize(){return sizeof(POD);}
    POD operator()(Vector<unsigned char> const& bytes)
    {
        assert(bytes.getSize() == byteSize());
        POD item;
        for(int i = 0; i < byteSize(); ++i)
            ((unsigned char*)&item)[i] = bytes[i];
        return item;
    }
    Vector<unsigned char> operator()(POD const& item)
    {
        Vector<unsigned char> bytes(byteSize());
        for(int i = 0; i < byteSize(); ++i)
            bytes[i] = ((unsigned char*)&item)[i];
        return bytes;
    }
};
template<typename POD> struct IntegralSerializer
{//common use case
    IntegralSerializer(){assert(is_integral<POD>::value);}
    constexpr static int byteSize(){return sizeof(POD);}
    POD operator()(Vector<unsigned char> const& bytes)
    {
        assert(bytes.getSize() == byteSize());
        return ReinterpretDecode(bytes);
    }
    Vector<unsigned char> operator()(POD const& item)
        {return ReinterpretEncode(item, byteSize());}
};

template<typename POD, typename SERIALIZER = CastSerializer<POD> >
class EMVector
{
    BlockFile blockFile;//must be first
    long long size;
    int itemsPerBlock()const{return blockFile.getBlockSize()/sizeof(POD);}
    long long block(long long i){return i/itemsPerBlock();}
    long long index(long long i){return i % itemsPerBlock();}
    SERIALIZER s;
    enum{HEADER_SIZE = 4};
    int extraItems()const{return blockFile.getSize() * itemsPerBlock() - size;}
    static int calculateBlockSize()
    {//if not exact try to go under, if not go over
        int result = BlockFile::targetBlockSize()/SERIALIZER::byteSize();
        if(result == 0) ++result;//can't go under, must go over
        return result * SERIALIZER::byteSize();
    }
    EMVector(EMVector const&);//no copying allowed
    EMVector& operator=(EMVector const&);
public:
    long long getSize(){return size;}
    EMVector(string const& filename, int cacheSize = 2): size(0),
        blockFile(filename, calculateBlockSize(), cacheSize, HEADER_SIZE)
    {
        assert(blockFile.getBlockSize() % SERIALIZER::byteSize() == 0);
        //check if file already exists - header is number of extra items
        if(blockFile.getSize() > 0) size = blockFile.getSize() *
            itemsPerBlock() - ReinterpretDecode(blockFile.readHeader());
    }
    ~EMVector()
    {//write number of extra items to header
        Vector<unsigned char> header =
            ReinterpretEncode(extraItems(), HEADER_SIZE);
        blockFile.writeHeader(header);
    }
    void append(POD const& item)
    {
        ++size;
        if(extraItems() < 0) blockFile.appendEmptyBlock();
        set(item, size - 1);
    }
    void set(POD const& item, long long i)
    {
        assert(i >= 0 && i < size);
        blockFile.set(s(item), block(i), index(i) * SERIALIZER::byteSize());
    }
    POD operator[](long long i)
    {
        assert(i >= 0 && i < size);
        return s(blockFile.get(block(i), index(i) * SERIALIZER::byteSize(),
            SERIALIZER::byteSize()));
    }
    void removeLast()
    {
        assert(size > 0);
        --size;
    }

    friend void IOSort(EMVector& vector)
    {
        {//scope to remove temp vector before file deletion
            long long n = vector.getSize(), C = sqrt(n *
                vector.itemsPerBlock()), Q = n/C, lastQSize = n % C;
            File::remove("IOSortTempFile.igmdk");//in case exists already
            EMVector temp("IOSortTempFile.igmdk");//potentially different
            //block size if old file and BUFSIZ changed, but ok
            typedef pair<POD, long long> HeapItem;
            Heap<HeapItem, PairFirstComparator<POD, long long> > merger;
            for(long long q = 0, i = 0; q < Q + 1; ++q)
            {
                long long m = q == Q ? lastQSize : C;
                if(m > 0)
                {//sort each block
                    Vector<POD> buffer;
                    for(long long j = 0; j < m; ++j)buffer.append(vector[i++]);
                    quickSort(buffer.getArray(), buffer.getSize());
                    //put smallest item of each block on the heap
                    merger.insert(HeapItem(buffer[0], q));
                    //and the rest to the temp vector
                    for(long long j = 1; j < m; ++j) temp.append(buffer[j]);
                }
            }
            Vector<Queue<POD> > buffers(Q + 1);
            Vector<long long> pointers(Q + 1, 0);
            for(long long i = 0; i < n; ++i)
            {//merge, remember that temp blocks are 1 less
                long long q = merger.getMin().second;
                vector.set(merger.deleteMin().first, i);
                if(buffers[q].isEmpty())//refill if needed
                    while(pointers[q] < (q == Q ? lastQSize : C) - 1)
                        buffers[q].push(temp[q * (C - 1) + pointers[q]++]);
                if(!buffers[q].isEmpty())//check if done with block
                    merger.insert(HeapItem(buffers[q].pop(), q));
            }
        }
        File::remove("IOSortTempFile.igmdk");
    }
};

}//end namespace
#endif
