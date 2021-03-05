#ifndef IGMDK_COMPRESSION_H
#define IGMDK_COMPRESSION_H
#include "../StringAlgorithms/SuffixArray.h"
#include "Stream.h"
#include "StaticCodes.h"
#include "HuffmanTree.h"
#include "LZW.h"
#include <cstdlib>
namespace igmdk{

enum {RLE_E1 = (1 << numeric_limits<unsigned char>::digits) - 1,
    RLE_E2 = RLE_E1 - 1};
Vector<unsigned char> RLECompress(Vector<unsigned char> const& byteArray)
{
    Vector<unsigned char> result;
    for(int i = 0; i < byteArray.getSize();)
    {
        unsigned char byte = byteArray[i++];
        result.append(byte);
        int count = 0;
        while(count < RLE_E2 - 1 && i + count < byteArray.getSize() &&
            byteArray[i + count] == byte) ++count;
        if(count > 1 || (byte == RLE_E1 && count == 1))
        {
            result.append(RLE_E1);
            result.append(count);
            i += count;
        }
        else if(byte == RLE_E1) result.append(RLE_E2);
    }
    return result;
}
Vector<unsigned char> RLEUncompress(Vector<unsigned char> const& byteArray)
{
    Vector<unsigned char> result;
    for(int i = 0; i < byteArray.getSize();)
    {
        unsigned char byte = byteArray[i++];
        if(byte == RLE_E1 && byteArray[i] != RLE_E1)
        {
            unsigned char count = byteArray[i++];
            if(count == RLE_E2) count = 1;
            else byte = result.lastItem();//need temp if vector reallocates
            while(count--) result.append(byte);
        }
        else result.append(byte);
    }
    return result;
}

Vector<unsigned char> MoveToFrontTransform(bool compress,
    Vector<unsigned char> const& byteArray)
{
    unsigned char list[1 << numeric_limits<unsigned char>::digits], j, letter;
    for(int i = 0; i < sizeof(list); ++i) list[i] = i;
    Vector<unsigned char> resultArray;
    for(int i = 0; i < byteArray.getSize(); ++i)
    {
        if(compress)
        {//find and output rank
            j = 0;
            letter = byteArray[i];
            while(list[j] != letter) ++j;
            resultArray.append(j);
        }
        else
        {//rank to byte
            j = byteArray[i];
            letter = list[j];
            resultArray.append(letter);
        }//move list back to make space for front item
        for(; j > 0; --j) list[j] = list[j - 1];
        list[0] = letter;
    }
    return resultArray;
}

Vector<unsigned char> BurrowsWheelerTransform(
    Vector<unsigned char> const& byteArray)
{
    int original = 0, size = byteArray.getSize();
    Vector<int> BTWArray = suffixArray<BWTRank>(byteArray.getArray(), size);
    Vector<unsigned char> result;
    for(int i = 0; i < size; ++i)
    {
        int suffixIndex = BTWArray[i];
        if(suffixIndex == 0)
        {//found the original string
            original = i;
            suffixIndex = size;//avoid the % size in next step
        }
        result.append(byteArray[suffixIndex - 1]);
    }//assume that 4 bytes is enough
    Vector<unsigned char> code = ReinterpretEncode(original, 4);
    for(int i = 0; i < code.getSize(); ++i) result.append(code[i]);
    return result;
}

Vector<unsigned char> BurrowsWheelerReverseTransform(
     Vector<unsigned char> const& byteArray)
{
    enum{M = 1 << numeric_limits<unsigned char>::digits};
    int counts[M], firstPositions[M], textSize = byteArray.getSize() - 4;
    for(int i = 0; i < M; ++i) counts[i] = 0;
    Vector<int> ranks(textSize);//compute ranks
    for(int i = 0; i < textSize; ++i) ranks[i] = counts[byteArray[i]]++;
    firstPositions[0] = 0;//compute first positions
    for(int i = 0; i < M - 1; ++i)
        firstPositions[i + 1] = firstPositions[i] + counts[i];
    Vector<unsigned char> index, result(textSize);//extract original rotation
    for(int i = 0; i < 4; ++i) index.append(byteArray[i + textSize]);
    //construct in reverse order
    for(int i = textSize - 1, ix = ReinterpretDecode(index); i >= 0; --i)
        ix = ranks[ix] + firstPositions[result[i] = byteArray[ix]];
    return result;
}

Vector<unsigned char> BWTCompress(Vector<unsigned char> const& byteArray)
{
    return HuffmanCompress(RLECompress(MoveToFrontTransform(true,
        BurrowsWheelerTransform(byteArray))));
}
Vector<unsigned char> BWTUncompress(Vector<unsigned char> const& byteArray)
{
    return BurrowsWheelerReverseTransform(MoveToFrontTransform(false,
       RLEUncompress(HuffmanUncompress(byteArray))));
}

}//end namespace
#endif
