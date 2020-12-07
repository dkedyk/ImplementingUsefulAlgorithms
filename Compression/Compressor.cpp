#include "../ExternalMemoryAlgorithms/File.h"
#include "../ExternalMemoryAlgorithms/CSV.h"
#include "../Utils/Debug.h"
#include <string>
#include "Compression.h"
using namespace std;

using namespace igmdk;

int compressor(File& in, File&out, bool compress, string const& smethod)
{
    char method;
    enum{HUF, BWT, LZW};
    if(smethod == "Huffman") method = HUF;
    else if(smethod == "BWT") method = BWT;
    else if(smethod == "LZW") method = LZW;
    else{DEBUG("Method Unknown"); return 0;}

    enum{N = 8096};
    unsigned char buffer[N];
    Vector<unsigned char> original, v;
    for(;;)
    {
        int size = min<long long>(N, in.bytesToEnd());
        in.read(buffer, size);
        for(int i = 0; i < size; ++i)
        {
            original.append(buffer[i]);
        }
        if(size < N) break;
    }
    if(compress)
    {
        if(method == LZW)
        {
            BitStream result;
            BitStream in(original);
            LZWCompress(in, result);
            v = ExtraBitsCompress(result.bitset);
        }
        else if(method == BWT)
        {
            v = BWTCompress(original);
        }
        else if(method == HUF)
        {
            v = HuffmanCompress(original);
        }
    }
    else
    {
        if(method == LZW)
        {
            BitStream in(ExtraBitsUncompress(original));
            BitStream result;
            LZWUncompress(in, result);
            v = result.bitset.getStorage();
        }
        else if(method == BWT)
        {
            v = BWTUncompress(original);
        }
        else if(method == HUF)
        {
            v = HuffmanUncompress(original);
        }
    }
    out.append(v.getArray(), v.getSize());
    return out.getSize();
}

void testAllMethods()
{
    //AAR decomp has bug for all
    string methods[] = {"Huf", "BWT", "LZW"}, files[] = {"a.txt", "bible.txt",
        "dickens.txt", "ecoli.txt", "mobydick.txt", "pi10mm.txt",//
        "world192.txt"};
    Vector<Vector<string> > matrix;
    Vector<string> titles;
    titles.append("File");
    titles.append("Size");
    for(int j = 0; j < sizeof(methods)/sizeof(methods[0]); ++j)
        titles.append(methods[j]);
    matrix.append(titles);
    for(int i = 0; i < sizeof(files)/sizeof(files[0]); ++i)
    {
        File in(files[i].c_str(), false);
        Vector<string> row;
        DEBUG(files[i]);
        row.append(files[i]);
        int oriSize = in.getSize();
        row.append(to_string(oriSize));
        for(int j = 0; j < sizeof(methods)/sizeof(methods[0]); ++j)
        {
            in.setPosition(0);
            DEBUG(methods[j]);
            int size;
            string outName = files[i] + "." + methods[j],
                backName = outName + ".ori";
            {
                File out(outName.c_str(), true);
                int start = clock();
                size = compressor(in, out, true, methods[j]);
                row.append(to_string(size));
                int elapsed = clock()-start;
            }
            {
                File out(outName.c_str(), false), back(backName.c_str(), true);
                int start = clock();
                int size2 = compressor(out, back, false, methods[j]);
                assert(oriSize == size2);
                int elapsed = clock()-start;
            }
            File::remove(outName.c_str());
            File::remove(backName.c_str());
        }
        matrix.append(row);
    }
    createCSV(matrix, "CompressionResult.csv");
}

int main(int argc, char *argv[])
{
    testAllMethods();
    return 0;
}
