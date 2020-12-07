#ifndef IGMDK_ARITHMETIC_CODING_H
#define IGMDK_ARITHMETIC_CODING_H
#include "../Utils/Debug.h"
#include "../Compression/Compression.h"
#include "../Compression/Stream.h"
#include <cstdlib>
namespace igmdk{

struct Arithmetic
{
    typedef unsigned long long WORD;
    enum{END = 256, EXTENDED_SIZE = END + 1};
    static WORD const N05 = 1 + (~0ull)/2, N025 = N05/2;
    WORD L, I, D, r, undecidedBits;
    void encodeSymbol(WORD l, WORD h, WORD t, BitStream& out)
    {
        WORD r = I / t, temp = r * l;
        if(h < t) I = r * h;
        for(L += temp, I -= temp; I <= N025; L *= 2, I *= 2)
        {
            if(L - 1 + I < N05) sendBit(0, out);
            else if(N05 <= L)
            {
                sendBit(1, out);
                L -= N05;
            }
            else
            {
                ++undecidedBits;
                L -= N025;
            }
        }
    }
    void sendBit(bool bit, BitStream& out)
    {
        out.writeBit(bit);
        while(undecidedBits)
        {
            out.writeBit(!bit);
            --undecidedBits;
        }
    }
    void encodeEOF(BitStream& out)
    {
        int b = numeric_limits<WORD>::digits - 1;
        sendBit(L >> b, out);
        out.writeValue(L, b);
    }
    void startEncoding()
    {
        L = 0;
        I = N05;
        undecidedBits = 0;
    }
    WORD decodeTarget(WORD t)
    {
        r = I / t;
        return min(t-1, D / r);
    }
    void adjustInterval(WORD l, WORD h, WORD t, BitStream& in)
    {
        WORD temp = r * l;
        D -= temp;
        if(h < t) I = r * h;
        I -= temp;
        while(I <= N025)
        {
            I *= 2;
            D = 2 * D + in.readBit();
        }
    }
    WORD fp[EXTENDED_SIZE];
    WORD size(WORD s){return s & -s;}
    WORD forw(WORD s){return s + size(s);}
    WORD back(WORD s){return s - size(s);}
    WORD getL(WORD s)
    {
        WORD l = 0;
        for(int i = s; i != 0; i = back(i)) l += fp[i-1];
        return l;
    }
    WORD getAndUpdateCount(WORD s)
    {
        WORD count = fp[s];
        for(int i = s; i != back(s+1); i = back(i)) count -= fp[i-1];
        for(int i = s+1; i <= EXTENDED_SIZE; i = forw(i)) ++fp[i-1];
        return count;
    }
    void setupFP(){for(int i = 0; i < EXTENDED_SIZE; ++i) fp[i] = i+1 - back(i+1);}
    Vector<unsigned char> adaptiveEncode(Vector<unsigned char>const& in)
    {
        BitStream out;
        setupFP();
        startEncoding();
        for(int i = 0; i <= in.getSize(); ++i)
        {
            int c = i < in.getSize() ? in[i] : END;
            WORD l = getL(c);
            encodeSymbol(l, l + getAndUpdateCount(c), i + EXTENDED_SIZE, out);
        }
        encodeEOF(out);
        return out.bitset.getStorage();
    }
    static Vector<unsigned char> AdaptiveCompress(Vector<unsigned char>const& inV)
    {
        Arithmetic en;
        return en.adaptiveEncode(inV);
    }
    int findCFP(WORD target)
    {
        int s = 0;
        for(int mid = 1 << lgCeiling(EXTENDED_SIZE); mid > 0; mid /= 2)
        {
            int i = s + mid - 1;
            if(i < EXTENDED_SIZE && fp[i] <= target)
            {
                target -= fp[i];
                s += mid;
            }
        }
        return s;
    }
    Vector<unsigned char> adaptiveDecode(Vector<unsigned char>const& inV)
    {
        BitStream in(inV);
        setupFP();
        Vector<unsigned char> result;
        I = N05;
        D = in.readValue(numeric_limits<WORD>::digits);
        for(int i = 0;; ++i)
        {
            WORD target = decodeTarget(i + EXTENDED_SIZE);
            int c = findCFP(target);
            WORD l = getL(c);
            adjustInterval(l, l + getAndUpdateCount(c), i + EXTENDED_SIZE, in);
            if(c == END) break;
            result.append(c);
        }
        return result;
    }
    static Vector<unsigned char> AdaptiveUncompress(Vector<unsigned char>const& inV)
    {
        Arithmetic ad;
        return ad.adaptiveDecode(inV);
    }
};

}
#endif
