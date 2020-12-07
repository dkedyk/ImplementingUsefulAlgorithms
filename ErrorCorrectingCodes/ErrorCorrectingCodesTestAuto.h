#ifndef IGMDK_ERROR_CORRECTING_CODES_TEST_AUTO_H
#define IGMDK_ERROR_CORRECTING_CODES_TEST_AUTO_H
#include "CRC.h"
#include "ReedSolomon.h"
#include "LDPC.h"
#include <cassert>
using namespace std;

namespace igmdk{

void testReedSolomonAuto()
{
    DEBUG("testReedSolomonAuto");
    Vector<unsigned char> message;
    string text = "Four score and seven years ago our fathers brought forth"
        " on this continent a new nation ...";
    for(int i = 0; i < text.size(); ++i) message.append(text[i]);
    int n = 255, p = (n - 223)/2;
    for(int j = 0; j < 1000; ++j)
    {
        ReedSolomon rs;
        Vector<unsigned char> code = rs.encodeBlock(rs.lengthPadBlock(message));
        for(int i = 0; i < p; ++i)//introduce upto p random errors
        {
            int location = GlobalRNG().mod(n);
            code[location] = GlobalRNG().mod(256);
        }
        pair<Vector<unsigned char>, bool> result = rs.decodeBlock(code);
        assert(result.second);
        result = rs.lengthUnpadBlock(result.first);
        assert(result.second);
        Vector<unsigned char> messageDecoded = result.first;
        assert(message == messageDecoded);
    }
    DEBUG("testReedSolomonAuto passed");
}

void testAllAutoErrorCorrectingCodes()
{
    DEBUG("testAllAutoErrorCorrectingCodes");
    testReedSolomonAuto();
}

}//end namespace
#endif
