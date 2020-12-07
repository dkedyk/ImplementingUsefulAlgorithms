#ifndef IGMDK_STREAM_H
#define IGMDK_STREAM_H
#include "../Utils/Bitset.h"
#include "../Utils/Vector.h"
namespace igmdk{

Vector<unsigned char> ReinterpretEncode(unsigned long long n, int size)
{
    assert(size > 0);
    enum{M = 1 << numeric_limits<unsigned char>::digits};
    Vector<unsigned char> result;
    while(size-- > 0)
    {
        result.append(n % M);
        n /= M;
    }
    return result;
}
unsigned long long ReinterpretDecode(Vector<unsigned char> const& code)
{
    assert(code.getSize() > 0);
    unsigned long long n = 0, base = 1;
    enum{M = 1 << numeric_limits<unsigned char>::digits};
    for(int i = 0; i < code.getSize(); ++i)
    {
        n += base * code[i];
        base *= M;
    }
    return n;
}

struct Stream
{
    unsigned long long position;
    Stream(): position(0) {}
};
struct BitStream : public Stream
{
    Bitset<unsigned char> bitset;//unsigned char for portability
    enum{B = numeric_limits<unsigned char>::digits};
    BitStream() {}
    BitStream(Bitset<unsigned char> const& aBitset): bitset(aBitset) {}
    BitStream(Vector<unsigned char> const& vector): bitset(vector) {}
    void writeBit(bool value){bitset.append(value);}
    bool readBit()
    {
        assert(bitsLeft());
        return bitset[position++];
    }
    void writeByte(unsigned char byte){writeValue(byte, B);}
    void writeBytes(Vector<unsigned char> const& bytes)
        {for(int i = 0; i < bytes.getSize(); ++i) writeByte(bytes[i]);}
    unsigned char readByte(){return readValue(B);}
    Vector<unsigned char> readBytes(int n)
    {
        assert(n <= bytesLeft());
        Vector<unsigned char> result(n);
        for(int i = 0; i < n; ++i) result[i] = readByte();
        return result;
    }
    void debug()const{bitset.debug();}
    void writeValue(unsigned long long value, int bits)
        {bitset.appendValue(value, bits);}
    unsigned long long readValue(int bits)
    {
        assert(bits <= bitsLeft());
        position += bits;
        return bitset.getValue(position - bits, bits);
    }
    unsigned long long bitsLeft()const{return bitset.getSize() - position;}
    unsigned long long bytesLeft()const{return bitsLeft()/B;}
};

}//end namespace
#endif
