#ifndef IGMDK_LARGE_NUMBER_TEST_AUTO_H
#define IGMDK_LARGE_NUMBER_TEST_AUTO_H
#include "LargeNumber.h"
#include "LargeRational.h"

namespace igmdk{

void testNumberAuto()
{
    assert((-Number(0) - Number(2)) == (Number(0) - Number(2)));
    assert((Number(2) - Number(0)) == (Number(2) - -Number(0)));

    Number m = power(Number(2), Number(128));
    m >>= 125;
    assert(m == Number(8));
    m <<= 125;
    assert(power(Number(2), Number(2)) == Number(4));
    assert(Number(4) % Number(3) == Number(1));
    assert(modInverse(Number(4), Number(7)) == Number(2));

    assert((Number(-11) % Number(103)) == Number(-11));
    assert(gcd(m, Number(3)) == Number(1));
    assert(sqrtInt(Number(99)) == Number(9));
    assert(modPower(Number(2), Number(2), Number(3)) == Number(1));
    assert(isPrime(Number((3))));
    assert(isPrime(Number((53))));
    assert(!isPrime(Number((616460792))));
    assert(Number(-23).toDecimalString() == "-23");
    assert(Number("-23") == Number(-23));
}

void testAllAutoLargeNumber()
{
    DEBUG("testAllAutoLargeNumber");
    testNumberAuto();
}

}//end namespace
#endif
