#include <cassert>
#include <iostream>
#include "ErrorCorrectingCodesTestAuto.h"
using namespace std;
using namespace igmdk;

void testLDPCAuto()//takes ~100 seconds on my pc
{//not quite auto need better statistical tests here
    int n = 20, k = 5;
    int nFailed = 0, nFalseSuccess = 0, nCodes = 100, nTests = 100;
    for(int m = 0; m < nCodes; ++m)
    {
        LDPC l(n, k);
        Bitset<> message(l.getNewK());
        for(int i = 0; i < message.getSize(); ++i)
            message.set(i, GlobalRNG().mod(2));//random message
        for(int j = 0; j < nTests; ++j)
        {
            Bitset<> code = l.encode(message);
            //below use the worst-case bound but need to try other values
            for(int i = 0; i < (n - l.getNewK())/2; ++i)
                code.set(GlobalRNG().mod(code.getSize()), GlobalRNG().mod(2));

            pair<Bitset<>, bool> result = l.decode(code);
            if(!result.second) ++nFailed;
            else if(message != result.first) ++nFalseSuccess;
        }
    }
    DEBUG(nFailed * 1.0/nTests/nCodes);//0.34 particular run
    DEBUG(nFalseSuccess * 1.0/nTests/nCodes);//0.09 particular run
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
        code.set(GlobalRNG().mod(code.getSize()), GlobalRNG().mod(2));
    DEBUG("code");
    code.debug();
    pair<Bitset<>, bool> result = l.decode(code);
    message = result.first;
    DEBUG(result.second);
    DEBUG("message");
    message.debug();//fails very occasionally
}

int main()
{
    testAllAutoErrorCorrectingCodes();
    testLDPC();
    testLDPCAuto();
    return 0;
}
