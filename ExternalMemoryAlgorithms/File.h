#ifndef IGMDK_FILE_H
#define IGMDK_FILE_H

#include <cassert>
#include <fstream>
#include <sstream>
#include <cstdio>//for remove and rename
#include "../Utils/Vector.h"
#include "../Compression/Stream.h"//for reinterpret code
#include "../MiscAlgs/LRUCache.h"//for LRU cache
using namespace std;
namespace igmdk{

string toStringDouble(double x)
{
    stringstream s;
    s << setprecision(17) << x;
    string result;
    s >> result;
    return result;
}

class File
{
    fstream f;//C++ IO object
    long long size;//cached for efficiency
    void create(const char* filename){ofstream dummy(filename, ios::trunc);}
    File(File const&);//no copying
    File& operator=(File const&);
    void goToEnd()
    {
        f.seekg(0, ios::end);
        assert(f);//watch out for external issues
    }
public:
    static bool exists(const char* filename){return bool(ifstream(filename));}
    static void remove(const char* filename)
    {
        if(exists(filename))
        {
            int returnCode = std::remove(filename);
            assert(returnCode == 0);//watch out for external issues
        }
    }
    File(const char* filename, bool truncate)
    {
        if(truncate || !exists(filename)) create(filename);
        f.open(filename, ios::binary | ios::in | ios::out);
        assert(f);//make sure f not locked, etc.
        //calculate size
        goToEnd();
        size = getPosition();
        setPosition(0);//come back to the beginning

    }
    long long getPosition(){return f.tellg();}
    long long getSize()const{return size;}
    long long bytesToEnd(){return getSize() - getPosition();}
    void setPosition(long long position)
    {
        assert(0 <= position && position <= getSize());
        f.seekg(position);
        assert(f);//watch out for external issues
    }
    void read(unsigned char* buffer, long long n)
    {
        assert(n > 0 && n <= bytesToEnd());
        f.read((char*)buffer, n);
        assert(f);//watch out for external issues
    }
    void write(unsigned char* buffer, long long n)
    {
        assert(n > 0 && n <= bytesToEnd());//to prevent errors, not needed
        f.write((char*)buffer, n);
        f.flush();
        assert(f);//watch out for external issues
    }
    void append(unsigned char* buffer, long long n)
    {
        goToEnd();
        size += n;//do this first to prevent write assertion fail
        write(buffer, n);//write from one-past-last
    }
};

//The below version also caches position - doesn't seem worth it but need
//further study
/*
class File
{
    fstream f;//C++ IO object
    long long size, position;//cached here for efficiency
    void create(const char* filename){ofstream dummy(filename, ios::trunc);}
    File(File const&);//no copying
    File& operator=(File const&);
    void goToEnd()
    {
        f.seekg(0, ios::end);
        assert(f);//watch out for external issues
        position = f.tellg();
    }
public:
    static bool exists(const char* filename){return ifstream(filename);}
    static void remove(const char* filename)
    {
        if(exists(filename))
        {
            int returnCode = std::remove(filename);
            assert(returnCode == 0);//watch out for external issues
        }
    }
    File(const char* filename, bool truncate)
    {//open an existing f or start with a blank one
        if(truncate || !exists(filename)) create(filename);
        f.open(filename, ios::binary | ios::in | ios::out);
        assert(f);//make sure f not locked, etc.
        //calculate size
        goToEnd();
        size = position;
        setPosition(0);//come back to the beginning
    }
    long long getSize()const{return size;}
    long long getPosition()const{return position;}
    //how many can be consumed one-by-one
    long long bytesToEnd()const{return getSize() - getPosition();}
    void setPosition(long long thePosition)
    {//likely to flush buffer
        assert(0 <= thePosition && thePosition <= getSize());
        if(thePosition != position)
        {
            position = thePosition;
            f.seekg(position);
            assert(f);//watch out for external issues
        }
    }
    void read(unsigned char* block, long long n)
    {
        assert(n > 0 && n <= bytesToEnd());
        f.read((char*)block, n);
        assert(f);//watch out for external issues
        position += n;
    }
    void write(unsigned char* block, long long n)
    {//after write choose to commit immediately
        assert(position == f.tellg());
        assert(n > 0 && n <= bytesToEnd());//to prevent errors, not needed
        f.write((char*)block, n);
        f.flush();
        assert(f);//watch out for external issues
        position += n;
    }
    void append(unsigned char* block, long long n)
    {
        goToEnd();
        size += n;//do this first to prevent write assertion fail
        write(block, n);//write from one-past-last
    }
};*/

class BlockFile
{
    File f;
    long long size;//the number of blocks, excluding header ones
    enum{SELF_HEADER_SIZE = 4};
    int headerSize, blockSize;
    int getNHeaderBlocks()const
        {return ceiling(SELF_HEADER_SIZE + headerSize, blockSize);}
    void setBlock(long long blockId){f.setPosition(blockId * blockSize);}
    void write(long long blockId, Vector<unsigned char> const& block)
    {
        assert(block.getSize() == blockSize);
        setBlock(blockId);
        f.write(block.getArray(), blockSize);
    }
    Vector<unsigned char> read(long long blockId)
    {
        Vector<unsigned char> block(blockSize);
        setBlock(blockId);
        f.read(block.getArray(), blockSize);
        return block;
    }
    typedef DelayedCommitLRUCache<long long, Vector<unsigned char>, BlockFile>
        CACHE;
    friend CACHE;//to allow access to read and write
    CACHE cache;//declared last to be destructed first
    Vector<unsigned char> getHelper(long long blockId, int start, int n)
    {
        Vector<unsigned char> data(n);
        Vector<unsigned char> const& block = cache.read(blockId);
        for(int i = 0; i < n; ++i) data[i] = block[start + i];
        return data;
    }
    void setHelper(Vector<unsigned char> const& data, long long blockId,
        int start)
    {
        Vector<unsigned char> block = cache.read(blockId);
        for(int i = 0; i < data.getSize(); ++i) block[start + i] = data[i];
        cache.write(blockId, block);
    }
public:
    constexpr static int targetBlockSize(){return max<int>(BUFSIZ, 4096);}
    int getBlockSize()const{return blockSize;}
    //header blocks not included in size
    long long getSize()const{return size - getNHeaderBlocks();}
    BlockFile(string const& filename, int theBlockSize, int cacheSize,
        int theHeaderSize = 0): f(filename.c_str(), false), size(0),
        headerSize(theHeaderSize), blockSize(theBlockSize),
        cache(*this, cacheSize)
    {
        assert(blockSize > 0);
        long long fileSize = f.getSize();
        if(fileSize > 0)
        {//already exists, importing settings, can't use getHelper because
            //don't know blockSize yet
            Vector<unsigned char> selfHeader(SELF_HEADER_SIZE);
            f.setPosition(0);
            f.read(selfHeader.getArray(), SELF_HEADER_SIZE);
            blockSize = ReinterpretDecode(selfHeader);
            assert(blockSize > 0 && fileSize % blockSize == 0);//basic check
            size = fileSize/blockSize;
        }
        else
        {//append header blocks and write blockSize to own header
            for(int i = 0; i < getNHeaderBlocks(); ++i) appendEmptyBlock();
            setHelper(ReinterpretEncode(blockSize, SELF_HEADER_SIZE), 0, 0);
        }
    }
    void appendEmptyBlock()
    {
        ++size;
        Vector<unsigned char> block(blockSize, 0);
        f.append(block.getArray(), blockSize);
    }
    Vector<unsigned char> get(long long blockId, int start, int n)
    {//for non-header blocks
        assert(0 <= blockId && blockId < getSize() && n > 0 && start >= 0 &&
            start + n <= blockSize);
        return getHelper(blockId + getNHeaderBlocks(), start, n);
    }
    void set(Vector<unsigned char> const& data, long long blockId, int start)
    {//for non-header blocks
        assert(0 <= blockId && blockId < getSize() && start >= 0 &&
            start + data.getSize() <= blockSize);
        setHelper(data, blockId + getNHeaderBlocks(), start);
    }
    void writeHeader(Vector<unsigned char> const& header)
    {//caller header may span several blocks
        assert(header.getSize() == headerSize);
        for(int i = 0, toWrite = headerSize; i < getNHeaderBlocks(); ++i)
        {
            int start = i == 0 ? SELF_HEADER_SIZE : 0,
                n = min(toWrite, blockSize - start);
            Vector<unsigned char> headerBlockData(n);
            for(int j = 0; j < n; ++j) headerBlockData[j] =
                header[(headerSize - toWrite) + j];
            setHelper(headerBlockData, i, start);
            toWrite -= n;
        }
    }
    Vector<unsigned char> readHeader()
    {//caller header may span several blocks
        assert(headerSize > 0);
        Vector<unsigned char> header;
        for(int i = 0, toRead = headerSize; i < getNHeaderBlocks(); ++i)
        {
            int start = i == 0 ? SELF_HEADER_SIZE : 0,
                n = min(toRead, blockSize - start);
            header.appendVector(getHelper(i, start, n));
            toRead -= n;
        }
        return header;
    }
};

}//end namespace
#endif
