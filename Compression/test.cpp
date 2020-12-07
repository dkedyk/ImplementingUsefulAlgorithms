#include "Compression.h"
#include "CompressionTestAuto.h"
#include <iostream>
using namespace igmdk;

void compress(Vector<unsigned char> const& byteArray)
{
    BitStream result;
    FibonacciEncode(1597, result);//first large value that loses to byte code
    result.bitset.debug();
    DEBUG(FibonacciDecode(result));

    GammaEncode(32, result);
    result.bitset.debug();
    DEBUG(GammaDecode(result));

    byteEncode(128 * 128, result);
    result.bitset.debug();
    DEBUG(byteDecode(result));

    HuffmanTree HuffTree(byteArray);
    Vector<unsigned char> MTF_Mississippi = MoveToFrontTransform(true, byteArray);

    cout << "breakpoint" << endl;//if have to recompute HUFFMAN do MISSISSIPPI not "large ascii text"
}

void timeRT()
{
    char text[] = //"abbabab";
    "mississippi";
    Vector<unsigned char> uncompressed;
    for(int i = 0; i < sizeof(text)-1; ++i) uncompressed.append(text[i]);

    cout << "uncompressed.size" << uncompressed.getSize() << endl;
    compress(uncompressed);
}

int main()
{
    timeRT();//fail Huffman has bug???
    testAllAutoCompression();
	return 0;
}
