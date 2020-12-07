#ifndef IGMDK_CRYPTOGRAPHY_TEST_AUTO_H
#define IGMDK_CRYPTOGRAPHY_TEST_AUTO_H
#include "Cryptography.h"
#include <string>
#include <cassert>
using namespace std;

namespace igmdk{

void testSimpleEncryptAuto()
{
    DEBUG("testSimpleEncryptAuto");
    string s = "top secret info", key = "123456", sDecrypted;
    enum{B = 128};
    Vector<unsigned char> data, password;
    for(int i = 0; i < B; ++i)
    {
        data.append(i < s.length() ? s[i] : 0);
        password.append(i < key.length() ? key[i] : 0);
    }
    assert(data == simpleDecrypt(simpleEncrypt(data, password),
        password).first);
    DEBUG("testSimpleEncryptAuto passed");
}

void testAllAutoCryptography()
{
    DEBUG("testAllAutoCryptography");
    testSimpleEncryptAuto();
}

}//end namespace
#endif
