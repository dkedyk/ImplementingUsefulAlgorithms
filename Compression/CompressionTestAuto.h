#ifndef IGMDK_COMPRESSION_TEST_AUTO_H
#define IGMDK_COMPRESSION_TEST_AUTO_H
#include <string>
using namespace std;
#include "Compression.h"

namespace igmdk{

void testGammaCodeAuto()
{
    DEBUG("testGammaCodeAuto");
    BitStream result;
    for(int i = 1; i < 1000; ++i) GammaEncode(i, result);
    for(int i = 1; i < 1000; ++i) assert(GammaDecode(result) == i);
    DEBUG("testGammaCodeAuto passed");
}

void testFibonacciCodeAuto()
{
    DEBUG("testFibonacciCodeAuto");
    BitStream result;
    for(int i = 1; i < 1000; ++i) FibonacciEncode(i, result);
    for(int i = 1; i < 1000; ++i) assert(FibonacciDecode(result) == i);
    DEBUG("testFibonacciCodeAuto passed");
}

void testByteCodeAuto()
{
    DEBUG("testGammaCodeAuto");
    BitStream result;
    for(int i = 0; i < 1000; ++i) byteEncode(i, result);
    for(int i = 0; i < 1000; ++i) assert(byteDecode(result) == i);
    DEBUG("testGammaCodeAuto passed");
}

Vector<unsigned char> getRandomBytes(int n = 10000)
{
    Vector<unsigned char> w(n, 0);
    for(int i = 0; i < n; ++i) w[i] = GlobalRNG().next();
    return w;
}
void testBWTCompressAuto()
{
    DEBUG("testBWTCompressAuto");
    Vector<unsigned char> byteArray = getRandomBytes();
    assert(byteArray == BWTUncompress(BWTCompress(byteArray)));
    DEBUG("testBWTCompressAuto passed");
}

void testLZWAuto()
{
    DEBUG("testLZWAuto");
    Vector<unsigned char> byteArray = getRandomBytes(), code;
    {
        BitStream in(byteArray);
        BitStream out;
        LZWCompress(in, out);
        code = ExtraBitsCompress(out.bitset);
    }
    {
        BitStream in(ExtraBitsUncompress(code));
        BitStream out;
        LZWUncompress(in, out);
        assert(byteArray == out.bitset.getStorage());
    }
    DEBUG("testLZWAuto passed");
}

void testAllAutoCompression()
{
    DEBUG("testAllAutoCompression");
    testGammaCodeAuto();
    testFibonacciCodeAuto();
    testByteCodeAuto();
    testBWTCompressAuto();
    testLZWAuto();
}

}//end namespace
#endif
