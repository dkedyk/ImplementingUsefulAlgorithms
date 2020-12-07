#include <cassert>
#include <iostream>
#include "ErrorCorrectingCodesTestAuto.h"
using namespace std;
using namespace igmdk;

void testLDPCAuto()
{
    int n = 20, k = 5;
    LDPC l(n, k);
    Bitset<> message(l.getNewK());
    for(int i = 0; i < message.getSize(); ++i)
        message.set(i, GlobalRNG().mod(2));//random message
    int nFailed = 0;
    for(int j = 0; j < 10; ++j)
    {
        Bitset<> code = l.encode(message);

        for(int i = 0; i < (n - l.getNewK()/2); ++i)
            code.set(GlobalRNG().mod(code.getSize()), GlobalRNG().mod(2));

        pair<Bitset<>, bool> result = l.decode(code);
        if(!result.second) ++nFailed;
        else assert(message == result.first);
    }
    DEBUG(nFailed);
}

void testLDPC()
{
    LDPC l(20, 5);
    Bitset<> message(l.getNewK());
    message.setAll();
    message.set(1, 0);
    message.set(0, 0);
    DEBUG("message");
    message.debug();
    Bitset<> code = l.encode(message);
    DEBUG("code");
    code.debug();
    for(int i = 0; i < 3; ++i)
        code.set(GlobalRNG().mod(20), GlobalRNG().mod(2));
    DEBUG("code");
    code.debug();
    pair<Bitset<>, bool> result = l.decode(code);
    message = result.first;
    DEBUG(result.second);
    DEBUG("message");
    message.debug();
}

int main()
{
    testAllAutoErrorCorrectingCodes();
    testLDPC();
    //testLDPCAuto();
    return 0;
}
