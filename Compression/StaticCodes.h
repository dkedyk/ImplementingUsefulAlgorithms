#ifndef IGMDK_STATIC_CODES_H
#define IGMDK_STATIC_CODES_H
#include "Stream.h"
#include <cstdlib>
namespace igmdk{

Vector<unsigned char> ExtraBitsCompress(Bitset<unsigned char> const& bitset)
{
    assert(bitset.getSize() > 0);//makes no sense otherwise
    Vector<unsigned char> result = bitset.getStorage();
    result.append(bitset.garbageBits());
    return result;
}
Bitset<unsigned char> ExtraBitsUncompress(Vector<unsigned char> byteArray)
{
    assert(byteArray.getSize() > 1 && byteArray.lastItem() < BitStream::B);
    int garbageBits = byteArray.lastItem();
    byteArray.removeLast();
    Bitset<unsigned char> result(byteArray);
    while(garbageBits--) result.removeLast();
    return result;
}

void byteEncode(unsigned long long n, BitStream& result)
{
    enum{M05 = 1 << (numeric_limits<unsigned char>::digits - 1)};
    do
    {
        unsigned char r = n % M05;
        n /= M05;
        if(n) r += M05;
        result.writeByte(r);
    }while(n);
}
unsigned long long byteDecode(BitStream& stream)
{
    unsigned long long n = 0, base = 1;
    enum{M05 = 1 << (numeric_limits<unsigned char>::digits - 1)};
    for(;; base *= M05)
    {
        unsigned char code = stream.readByte(), value = code % M05;
        n += base * value;
        if(value == code) break;
    }
    return n;
}

void UnaryEncode(int n, BitStream& result)
{
    while(n--) result.writeBit(true);
    result.writeBit(false);
}
int UnaryDecode(BitStream& code)
{
    int n = 0;
    while(code.readBit()) ++n;
    return n;
}

void GammaEncode(unsigned long long n, BitStream& result)
{
    assert(n > 0);
    int N = lgFloor(n);
    UnaryEncode(N, result);
    if(N > 0) result.writeValue(n - twoPower(N), N);
}
unsigned long long GammaDecode(BitStream& code)
{
    int N = UnaryDecode(code);
    return twoPower(N) + (N > 0 ? code.readValue(N) : 0);
}

void advanceFib(unsigned long long& f1, unsigned long long& f2)
{
    unsigned long long temp = f2;
    f2 += f1;
    f1 = temp;
}
void FibonacciEncode(unsigned long long n, BitStream& result)
{
    assert(n > 0);
    //find largest fib number f1 <= n
    unsigned long long f1 = 1, f2 = 2;
    while(f2 <= n) advanceFib(f1, f2);
    //mark the numbers from highest to lowest
    Bitset<unsigned char> reverse;
    while(f2 > 1)
    {
        reverse.append(n >= f1);
        if(n >= f1) n -= f1;
        unsigned long long temp = f1;
        f1 = f2 - f1;
        f2 = temp;
    }//change order to lowest to highest and add terminator
    reverse.reverse();
    result.bitset.appendBitset(reverse);
    result.writeBit(true);
}
unsigned long long FibonacciDecode(BitStream& code)
{
    unsigned long long n = 0, f1 = 1, f2 = 2;
    for(bool prevBit = false;; advanceFib(f1, f2))
    {//add on the next Fibonacci number until see 11
        bool bit = code.readBit();
        if(bit)
        {
            if(prevBit) break;
            n += f1;
        }
        prevBit = bit;
    }
    return n;
}

}//end namespace
#endif
