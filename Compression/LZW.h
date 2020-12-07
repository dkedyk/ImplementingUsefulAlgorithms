#ifndef IGMDK_LZW_H
#define IGMDK_LZW_H
#include "../RandomTreap/Trie.h"
#include "Stream.h"
#include <cstdlib>
namespace igmdk{

void LZWCompress(BitStream& in, BitStream& out, int maxBits = 16)
{
    assert(in.bytesLeft());
    byteEncode(maxBits, out);//store as config
    TernaryTreapTrie<int> dictionary;
    TernaryTreapTrie<int>::Handle h;
    int n = 0;
    while(n < (1 << numeric_limits<unsigned char>::digits))
    {//initialize with all bytes
        unsigned char letter = n;
        dictionary.insert(&letter, 1, n++);
    }
    Vector<unsigned char> word;
    while(in.bytesLeft())
    {
        unsigned char c = in.readByte();
        word.append(c);
        //if found keep appending
        if(!dictionary.findIncremental(word.getArray(), word.getSize(), h))
        {//word without the last byte guaranteed to be in the dictionary
            out.writeValue(*dictionary.find(word.getArray(),
                word.getSize() - 1), lgCeiling(n));
            if(n < twoPower(maxBits))//add new word if have space
                dictionary.insert(word.getArray(), word.getSize(), n++);
            word = Vector<unsigned char>(1, c);//set to read byte
        }
    }
    out.writeValue(*dictionary.find(word.getArray(), word.getSize()),
        lgCeiling(n));
}

void LZWUncompress(BitStream& in, BitStream& out)
{
    int maxBits = byteDecode(in), size = twoPower(maxBits), n = 0,
        lastIndex = -1;
    assert(maxBits >= numeric_limits<unsigned char>::digits);
    Vector<Vector<unsigned char> > dictionary(size);
    for(; n < (1 << numeric_limits<unsigned char>::digits); ++n)
        dictionary[n].append(n);
    while(in.bitsLeft())
    {
        int index = in.readValue(lastIndex == -1 ? 8 :
            min(maxBits, lgCeiling(n + 1)));
        if(lastIndex != -1 && n < size)
        {
            Vector<unsigned char> word = dictionary[lastIndex];
            word.append((index == n ? word : dictionary[index])[0]);
            dictionary[n++] = word;
        }
        for(int i = 0; i < dictionary[index].getSize(); ++i)
            out.writeByte(dictionary[index][i]);
        lastIndex = index;
    }
}

}//end namespace
#endif
