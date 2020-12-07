#ifndef IGMDK_BITSET_H
#define IGMDK_BITSET_H
#include "../Utils/Vector.h"
#include "../Utils/Utils.h"
#include "Bits.h"
using namespace std;
namespace igmdk{

//update
template<typename WORD = unsigned long long> class Bitset
{
    enum{B = numeric_limits<WORD>::digits};
    unsigned long long bitSize;//must come before storage
    Vector<WORD> storage;
    void zeroOutRemainder()
    {
        if(bitSize > 0)
            storage.lastItem() &= Bits::lowerMask(lastWordBits());
    }
    bool get(int i)const
    {
        assert(i >= 0 && i < bitSize);
        return Bits::get(storage[i/B], i % B);
    }
    unsigned long long wordsNeeded()const{return ceiling(bitSize, B);}
public:
    Bitset(unsigned long long initialSize = 0): bitSize(initialSize),
        storage(wordsNeeded(), 0){}//set all bits to 0
    Bitset(Vector<WORD> const& vector): bitSize(B * vector.getSize()),
        storage(vector) {}//direct construction from storage vector is useful
    int lastWordBits()const
    {//1 to B if > 1 bit
        assert(bitSize > 0);
        int result = bitSize % B;
        if(result == 0) result = B;
        return result;
    }
    int garbageBits()const{return bitSize > 0 ? B - lastWordBits() : 0;}
    Vector<WORD> const& getStorage()const{return storage;}
    unsigned long long getSize()const{return bitSize;}
    unsigned long long wordSize()const{return storage.getSize();}
    bool operator[](int i)const{return get(i);}
    void set(int i, bool value = true)
    {
        assert(i >= 0 && i < bitSize);
        Bits::set(storage[i/B], i % B, value);
    }
    void append(bool value)
    {//increase size if needed, and get last bit
        ++bitSize;
        if(wordSize() < wordsNeeded()) storage.append(0);
        set(bitSize - 1, value);
    }
    void removeLast()
    {
        assert(bitSize > 0);
        if(lastWordBits() == 1) storage.removeLast();//shrink if can
        --bitSize;
        zeroOutRemainder();//removed bit might have been 1
    }
    void setAll(bool value = true)
    {
        for(int i = 0; i < wordSize(); ++i)
            storage[i] = value ? Bits::FULL : Bits::ZERO;
        zeroOutRemainder();
    }
    bool isZero()const
    {
        for(int i = 0; i < wordSize(); ++i) if(storage[i])return false;
        return true;
    }
    unsigned long long getValue(int i, int n)const
    {
        assert(n <= numeric_limits<unsigned long long>::digits && i >= 0 &&
            i + n <= bitSize && n > 0);
        unsigned long long result = 0;
        for(int word = i/B, bit = i % B, shift = 0; n > 0; bit = 0)
        {//get lower bits first
            int m = min(n, B - bit);//all bits or as much as the word has
            result |= Bits::getValue(storage[word++], bit, m) << shift;
            shift += m;
            n -= m;
        }
        return result;
    }
    void setValue(unsigned long long value, int i, int n)
    {
        assert(n <= numeric_limits<unsigned long long>::digits && i >= 0 &&
            i + n <= bitSize && n > 0);
        for(int word = i/B, bit = i % B, shift = 0; n > 0; bit = 0)
        {//set lower bits first
            int m = min(n, B - bit);//all bits or as much as the word has
            Bits::setValue(storage[word++], value >> shift, bit, m);
            shift += m;
            n -= m;
        }
    }
    void appendValue(unsigned long long value, int n)
    {
        int start = bitSize;
        bitSize += n;
        int k = wordsNeeded() - wordSize();
        for(int i = 0; i < k; ++i) storage.append(0);
        setValue(value, start, n);
    }
    void appendBitset(Bitset const& rhs)
    {//append storage words and remove extra words if any
        if(rhs.getSize() > 0)
        {
            for(int i = 0; i < rhs.wordSize(); ++i)
                appendValue(rhs.storage[i], B);
            bitSize -= B - rhs.lastWordBits();
            if(wordSize() > wordsNeeded()) storage.removeLast();
        }
    }

    bool operator==(Bitset const& rhs)const{return storage == rhs.storage;}
    Bitset& operator&=(Bitset const& rhs)
    {//only makes sense for equal sizes
        assert(bitSize == rhs.bitSize);
        for(int i = 0; i < wordSize(); ++i) storage[i] &= rhs.storage[i];
        return *this;
    }
    void flip()
    {
        for(int i = 0; i < wordSize(); ++i) storage[i] = ~storage[i];
        zeroOutRemainder();
    }

    Bitset& operator|=(Bitset const& rhs)
    {//only makes sense for equal sizes
        assert(bitSize == rhs.bitSize);
        for(int i = 0; i < wordSize(); ++i) storage[i] |= rhs.storage[i];
        return *this;
    }
    Bitset& operator^=(Bitset const& rhs)
    {//only makes sense for equal sizes
        assert(bitSize == rhs.bitSize);
        for(int i = 0; i < wordSize(); ++i) storage[i] ^= rhs.storage[i];
        return *this;
    }
    Bitset& operator>>=(int shift)
    {//shift by 0 no-op
        if(shift < 0) return operator<<=(-shift);
        int normalShift = shift % bitSize, wordShift = normalShift/B,
            bitShift = normalShift % B;
        if(wordShift > 0)//shift words
            for(int i = 0; i + wordShift < wordSize(); ++i)
            {
                storage[i] = storage[i + wordShift];
                storage[i + wordShift] = 0;
            }
        if(bitShift > 0)//shift bits
        {//for word layout 00000101|00111000 >>= 4 -> 10000000|00000011
            WORD carry = 0;
            for(int i = wordSize() - 1 - wordShift; i >= 0; --i)
            {
                WORD tempCarry = storage[i] << (B - bitShift);
                storage[i] >>= bitShift;
                storage[i] |= carry;
                carry = tempCarry;
            }
        }
        return *this;
    }
    Bitset& operator<<=(int shift)
    {
        if(shift < 0) return operator>>=(-shift);
        int normalShift = shift % bitSize, wordShift = normalShift/B,
            bitShift = normalShift % B;
        if(wordShift > 0)//shift words
            for(int i = wordSize() - 1; i - wordShift >= 0; --i)
            {
                storage[i] = storage[i - wordShift];
                storage[i - wordShift] = 0;
            }
        if(bitShift > 0)//shift bits
        {//for word layout 10000000|00000011 <<= 4 -> 00000000|00111000
            WORD carry = 0;
            for(int i = wordShift; i < wordSize(); ++i)
            {
                WORD tempCarry = storage[i] >> (B - bitShift);
                storage[i] <<= bitShift;
                storage[i] |= carry;
                carry = tempCarry;
            }
        }//some 1 bits could have shifted into the remainder
        zeroOutRemainder();
        return *this;
    }
    void debug()const
    {
        DEBUG(bitSize);
        for(int i = 0; i < bitSize; ++i) cout << get(i);
        cout << endl;
    }
    void reverse()
    {//fill up garbage bits
        int nFill = garbageBits();
        bitSize += nFill;
        (*this)<<=(nFill);
        //reverse storage words
        storage.reverse();
        for(int i = 0; i < wordSize(); ++i)
            storage[i] = reverseBits(storage[i]);
        //delete the garbage
        bitSize -= nFill;
        zeroOutRemainder();
    }
    int popCount()const
    {
        int sum = 0;
        for(int i = 0; i < wordSize(); ++i) sum += popCountWord(storage[i]);
        return sum;
    }
};

}//end namespace
#endif
