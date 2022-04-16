
#include "File.h"
#include "EMVector.h"
#include "EMBTree.h"
#include "CSV.h"
#include "../Utils/Debug.h"
#include "ExternalMemoryAlgorithmsTestAuto.h"
using namespace std;

using namespace igmdk;

void testCSV()
{
    Vector<Vector<string> > matrix;
    Vector<string> row1;
    row1.append("Animal Recognition");
    row1.append("Child");
    row1.append("1");
    row1.append("Computer");
    row1.append("5");
    matrix.append(row1);
    Vector<string> row2;
    row2.append("Simple Multiplication");
    row2.append("Child");
    row2.append("2");
    row2.append("Computer");
    row2.append("1");
    matrix.append(row2);

    Vector<string> names;
    names.append("Quality");
    createAugmentedCSVFiles(matrix, names, "CVSTest.csv");
}

void DDDEMVector()
{
    {
        EMVector<int> EMVector0to9B16("EMVector.igmdk", 16);
        int K = 10;
        for(int i = 0; i < K; ++i)
        {
            EMVector0to9B16.append(i);
        }
        cout << "breakpoint" << endl;
    }
    File::remove("EMVector.igmdk");
}

//for testing only, setting buffer doesn't seem worth it
class FileWithBuffSet
{
    Vector<char> buffer;//must be first to be destroyed last
    fstream file;//C++ IO object
    long long size, position;//cached here for efficiency
    void create(const char* filename){ofstream dummy(filename, ios::trunc);}
    FileWithBuffSet(FileWithBuffSet const&);//no copying
    FileWithBuffSet& operator=(FileWithBuffSet const&);
    void goToEnd()
    {
        file.seekg(0, ios::end);
        position = file.tellg();
    }
public:
    static bool exists(const char* filename){return bool(ifstream(filename));}
    static int remove(const char* filename){return std::remove(filename);}
    FileWithBuffSet(const char* filename, bool truncate, int B = -1)
    {//open an existing file or start with a blank one
        if(truncate || !exists(filename)) create(filename);
        if(B == 0) file.rdbuf()->pubsetbuf(0, 0);
        else if(B > 0)
        {//must be done before open file to have effect
            buffer = Vector<char>(B);
            file.rdbuf()->pubsetbuf(buffer.getArray(), B);
        }
        file.open(filename, ios::binary | ios::in | ios::out);
        assert(file);//make sure file not locked, etc.
        //calculate size
        goToEnd();
        size = position;
        position = 0;
        file.seekg(0);
    }
    long long getSize()const{return size;}
    long long getPosition()const{return position;}
    long long bytesLeft()const{return getSize() - getPosition();}
    void setPosition(long long thePosition)
    {//likely to flush buffer
        assert(0 <= thePosition && thePosition < getSize());
        if(thePosition != position)
        {
            position = thePosition;
            file.seekg(position);
        }
    }
    void read(unsigned char* buffer, long long n)
    {
        assert(n <= bytesLeft());
        file.read((char*)buffer, n);
    }
    void write(unsigned char* buffer, long long n)
    {//after write choose to commit immediately
        assert(n <= bytesLeft());//to prevent errors, not actually needed
        position += n;
        file.write((char*)buffer, n);
        file.flush();
    }
    void append(unsigned char* buffer, long long n)
    {
        goToEnd();
        size += n;
        write(buffer, n);
    }
};

void testBufferSizeSlowRead(int B)
{
    {
        int sum = 0;
        FileWithBuffSet f("File.igmdk", true, B);
        int n = 100000;
        Vector<unsigned char> buffer(n, -1);
        f.append(buffer.getArray(), n);
        //write is fast due to internal buffer, now check read
        for(int j = 0; j < 100; ++j)
        {
            f.setPosition(0);
            for(int i = 0; i < n; ++i)
            {
                //f.setPosition(i);
                Vector<unsigned char> buffer2(1);
                f.read(buffer2.getArray(), 1);
                sum += buffer2[0];
            }
        }
        DEBUG(sum);
    }
    File::remove("File.igmdk");
}

void testBufferSize(int B, int BSet = -1)
{
    if(BSet == -1)
    {
        DEBUG("Default");
        File f("File.igmdk", true);
        int n = 10000000;
        Vector<unsigned char> buffer(n, -1);
        f.append(buffer.getArray(), n);
        //write is fast due to internal buffer, now check read
        for(int i = 0; i < 100000; ++i)
        {
            int start = GlobalRNG().mod(n);
            if(start + B < n)
            {
                Vector<unsigned char> buffer2(B);
                f.setPosition(start);
                f.read(buffer2.getArray(), B);
            }
        }
    }
    else
    {
        DEBUG(BSet);
        FileWithBuffSet f("File.igmdk", true, BSet);
        int n = 10000000;
        Vector<unsigned char> buffer(n, -1);
        f.append(buffer.getArray(), n);
        //write is fast due to internal buffer, now check read
        for(int i = 0; i < 100000; ++i)
        {
            int start = GlobalRNG().mod(n);
            if(start + B < n)
            {
                Vector<unsigned char> buffer2(B);
                f.setPosition(start);
                f.read(buffer2.getArray(), B);
            }
        }
    }
    File::remove("File.igmdk");
}

void testBufferSizeDriver()
{
    Vector<Vector<string> > sensitivityMatrix, bufferMatrix;
    DEBUG(BUFSIZ);
    Vector<string> titles;
    titles.append("B");
    titles.append("Seconds");
    sensitivityMatrix.append(titles);
    titles.append("Seconds Bset = 0");//guaranteed small
    titles.append("Seconds Bset = 4096");//just right
    titles.append("Seconds Bset = 2^16");//too large
    titles.append("Seconds Bset = B Slow Read");//too large
    bufferMatrix.append(titles);
    //for(int b = 11; b <= 12; ++b)
    for(int b = 2; b <= 16; ++b)
    {
        int B = twoPower(b);
        Vector<string> row;
        DEBUG(B);
        row.append(to_string(B));

        int now = clock();
        testBufferSize(B);
        double time = (clock() - now) * 1.0/CLOCKS_PER_SEC;
        DEBUG(time);
        row.append(to_string(time));

        sensitivityMatrix.append(row);

        now = clock();
        testBufferSize(B, 0);
        time = (clock() - now) * 1.0/CLOCKS_PER_SEC;
        DEBUG(time);
        row.append(to_string(time));

        now = clock();
        testBufferSize(B, 4096);
        time = (clock() - now) * 1.0/CLOCKS_PER_SEC;
        DEBUG(time);
        row.append(to_string(time));

        now = clock();
        testBufferSize(B, twoPower(16));
        time = (clock() - now) * 1.0/CLOCKS_PER_SEC;
        DEBUG(time);
        row.append(to_string(time));

        now = clock();
        testBufferSizeSlowRead(B);
        time = (clock() - now) * 1.0/CLOCKS_PER_SEC;
        DEBUG(time);
        row.append(to_string(time));

        bufferMatrix.append(row);
    }
    createCSV(sensitivityMatrix, "BlockSensitivity.csv");
    createCSV(bufferMatrix, "BufferSize.csv");
}

void testBlockSizeSeq(int size)
{
    {
        EMVector<int> v("EMVector.igmdk", size);
        int n = 100000, sum = 0;
        DEBUG("start append");
        for(int i = 0; i < n; ++i)
        {
            //DEBUG(i);
            v.append(i);
        }
        DEBUG("append done");
        for(int i = 0; i < n; ++i) sum ^= v[i];
    }
    File::remove("EMVector.igmdk");
}
void testBlockSizeRand(int size)
{
    {
        EMVector<int> v("EMVector.igmdk", size);
        int n = 100000, sum = 0;
        for(int i = 0; i < n; ++i)
        {
            v.append(i);
        }
        for(int i = 0; i < n; ++i) sum ^= v[GlobalRNG().mod(n)];
    }
    File::remove("EMVector.igmdk");
}

void testVectorBlockSizeDriver()
{
    Vector<Vector<string> > matrix;
    DEBUG(BUFSIZ);
    Vector<string> titles;
    titles.append("B");
    titles.append("Seconds vec int sequential 10^5");
    titles.append("Seconds vec int rand 10^5");
    matrix.append(titles);
    for(int b = 2; b <= 2; ++b)//16
    {
        Vector<string> row;
        int B = twoPower(b);
        DEBUG(B);
        row.append(to_string(B));
        int now = clock();
        testBlockSizeSeq(B);
        double time = (clock() - now) * 1.0/CLOCKS_PER_SEC;
        DEBUG(time);
        row.append(to_string(time));
        now = clock();
        testBlockSizeRand(B);
        time = (clock() - now) * 1.0/CLOCKS_PER_SEC;
        DEBUG(time);
        row.append(to_string(time));
        matrix.append(row);
    }
    //createCSV(matrix, "VectorBChoice.csv");
}




int main()
{
    testCSV();
    //return 0;
    testVectorBlockSizeDriver();
    //return 0;
    //testBufferSizeDriver();
    //return 0;

    DDDEMVector();
    testAllAutoExternalMemoryAlgorithms();

	return 0;
}
