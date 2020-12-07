#ifndef IGMDK_CRYPTOGRAPHY_H
#define IGMDK_CRYPTOGRAPHY_H
#include "../RandomNumberGeneration/Random.h"
#include "../Utils/Utils.h"
#include "../ErrorCorrectingCodes/CRC.h"
#include "../Compression/Compression.h"
using namespace std;

namespace igmdk{

void applyARC4(uint32_t seed, Vector<unsigned char> temp,
    Vector<unsigned char>& data)
{
    for(int i = 0; i < temp.getSize(); ++i)
        temp[i] ^= (seed = xorshiftTransform(seed));
    ARC4 arc4(temp.getArray(), temp.getSize());
    for(int i = 0; i < data.getSize(); ++i) data[i] ^= arc4.nextByte();
}
Vector<unsigned char> simpleEncrypt(Vector<unsigned char> data,
    Vector<unsigned char> const& key)
{
    uint32_t seed = time(0), s = sizeof(int);
    CRC32 crc32;
    Vector<unsigned char> theSeed = ReinterpretEncode(seed, s), crc =
        ReinterpretEncode(crc32.hash(data.getArray(), data.getSize()), s);
    for(int i = 0; i < s; ++i) data.append(crc[i]);
    applyARC4(seed, key, data);
    for(int i = 0; i < s; ++i) data.append(theSeed[i]);
    return data;
}
pair<Vector<unsigned char>, bool> simpleDecrypt(Vector<unsigned char> code,
    Vector<unsigned char> const& key)
{
    assert(code.getSize() >= 8);
    enum{s = sizeof(uint32_t)};
    Vector<unsigned char> seed, crc;
    for(int i = 0; i < s; ++i) seed.append(code[code.getSize() + i - 4]);
    for(int i = 0; i < s; ++i) code.removeLast();
    applyARC4(ReinterpretDecode(seed), key, code);
    for(int i = 0; i < s; ++i) crc.append(code[code.getSize() + i - 4]);
    for(int i = 0; i < s; ++i) code.removeLast();
    CRC32 crc32;
    return make_pair(code, crc32.hash(code.getArray(), code.getSize()) ==
        ReinterpretDecode(crc));
}

}//end namespace
#endif
