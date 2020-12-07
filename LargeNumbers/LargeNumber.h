#ifndef IGMDK_LARGE_NUMBER_H
#define IGMDK_LARGE_NUMBER_H
#include "../Utils/Bits.h"//for lgFloor
#include "../Utils/Utils.h"
#include "../Utils/Debug.h"
#include "../RandomNumberGeneration/Random.h"//for prime generation
#include "../Utils/Vector.h"
#include <string>//for decimal conversion
#include <algorithm>//for string reverse
using namespace std;
namespace igmdk{

class Number : public ArithmeticType<Number>
{
    bool isMinus;
    typedef uint32_t DIGIT;
    typedef unsigned long long LARGE_DIGIT;
    enum{BASE_RADIX = numeric_limits<DIGIT>::digits};
    Vector<DIGIT> digits;
    DIGIT getDigit(int i)const{return i < nDigits() ? digits[i] : 0;}
    Number(int size, DIGIT fill): digits(size, fill), isMinus(false) {}
    void trim(){while(nDigits() > 1 && isZero()) digits.removeLast();}

    static DIGIT fullAdder(DIGIT a, DIGIT b, bool& carry)
    {
        LARGE_DIGIT sum = LARGE_DIGIT(a) + b + carry;
        carry = sum >> BASE_RADIX;
        return sum;
    }
    static Number add(Number const& a, Number const& b)
    {//O(|a|+|b|)
        int n = max(a.nDigits(), b.nDigits());
        Number result(n + 1, 0);
        bool carry = 0;
        for(int i = 0; i < n; ++i)
            result[i] = fullAdder(a.getDigit(i), b.getDigit(i), carry);
        result[n] = carry;
        result.trim();
        return result;
    }

    static void sub(Number& a, Number const& b)
    {//O(|a| + |b|)
        bool carry = 0;
        for(int i = 0; i < a.nDigits(); ++i)
        {
            LARGE_DIGIT digit = LARGE_DIGIT(b.getDigit(i)) + carry;
            carry = a[i] < digit;
            a[i] -= digit;//implicitly mod B
        }
        a.trim();
    }

    static DIGIT digitMult(DIGIT a, DIGIT b, DIGIT& carry)
    {
        LARGE_DIGIT prod = LARGE_DIGIT(a) * b;
        carry = prod >> BASE_RADIX;
        return prod;
    }
    static Number mult(Number const& a, DIGIT const& b)
    {
        Number result(a.nDigits() + 1, 0);
        bool addCarry = 0;
        DIGIT multCarry = 0, newMultCarry;
        for(int i = 0; i < a.nDigits(); ++i)
        {
            result[i] = fullAdder(digitMult(a[i], b, newMultCarry), multCarry,
                addCarry);
            multCarry = newMultCarry;
        }
        result[a.nDigits()] = fullAdder(0, multCarry, addCarry);
        result.trim();
        return result;
    }

    static DIGIT findK(Number const& s, Number  const& b)
    {//O(|s|), find k such that 0 <= k < BASE and kb <= s < (k + 1)b;
        DIGIT guess = s.digits.lastItem()/b.digits.lastItem();
        if(s.nDigits() > b.nDigits()) guess = (s.digits.lastItem() *
            (1ull << BASE_RADIX) + s[s.nDigits() - 2])/b.digits.lastItem();
        while(s < mult(b, guess)) --guess;//executes <= 2 times
        return guess;
    }
public:
    typedef DIGIT DIGIT_TYPE;
    int nDigits()const{return digits.getSize();}
    DIGIT& operator[](unsigned int i){return digits[i];}
    DIGIT const& operator[](unsigned int i)const{return digits[i];}
    void appendDigit(DIGIT const& digit){digits.append(digit);}
    bool isZero()const{return digits.lastItem() == 0;}
    bool isPositive()const{return !isMinus && !isZero();}
    bool isNegative()const{return isMinus && !isZero();}
    void negate(){isMinus = !isMinus;}
    Number abs()const{return isMinus ? -*this : *this;}
    bool isOdd()const{return digits[0] % 2;}
    bool isEven()const{return !isOdd();}
    Number(): isMinus(false), digits(1, 0) {}//default is 0
    //convenience constructor for single-digit numbers
    explicit Number(long long x): isMinus(x < 0), digits(1, std::abs(x))
        {assert(std::abs(x) <= numeric_limits<DIGIT_TYPE>::max());}

    bool absLess(Number const& rhs)const
    {//more digits larger, else MSD digit decides
        if(nDigits() != rhs.nDigits()) return nDigits() < rhs.nDigits();
        for(int i = nDigits() - 1; i >= 0; --i)
            if(digits[i] != rhs[i]) return digits[i] < rhs[i];
        return false;
    }
    bool absEqual(Number const& rhs)const
    {//more digits larger, else MSD digit decides
        if(nDigits() != rhs.nDigits()) return false;
        for(int i = 0; i < nDigits(); ++i)
            if(digits[i] != rhs[i]) return false;
        return true;
    }
    bool operator<(Number const& rhs)const//handle 0 = -0
        {return (isMinus && !rhs.isMinus && !isZero()) || absLess(rhs);}
    bool operator==(Number const& rhs)const//handle 0 = -0
        {return (isMinus == rhs.isMinus || isZero()) && absEqual(rhs);}

    Number operator-()const
    {
        Number result = *this;
        result.negate();
        return result;
    }
    Number& operator+=(Number const& rhs)
    {
        if(isMinus == rhs.isMinus)
        {//if same sign add
            digits = add(*this, rhs).digits;
            return *this;
        }
        else return *this -= -rhs;//else subtract
    }

    Number& operator-=(Number const& rhs)
    {
        if(isMinus == rhs.isMinus)
        {//if same sign subtract
            if(absLess(rhs))
            {
                Number temp = rhs;
                sub(temp, *this);
                temp.negate();
                *this = temp;
            }
            else sub(*this, rhs);
            trim();
            return *this;
        }
        else return *this += -rhs;//else add
    }

    Number& operator>>=(unsigned int k)
    {
        int skipDigits = k/BASE_RADIX, last = nDigits() - skipDigits,
            skipBits = k % BASE_RADIX;
        if(skipDigits > 0)//first apply whole-digit shifts
            for(int i = 0; i < nDigits(); ++i)
                digits[i] = i < last ? digits[i + skipDigits] : 0;
        if(skipBits > 0)//then bit shifts
        {
            DIGIT carry = 0, tempCarry;
            for(int i = last - 1; i >= 0; --i)
            {
                tempCarry = digits[i] << (BASE_RADIX - skipBits);
                digits[i] = (digits[i] >> skipBits) | carry;
                carry = tempCarry;
            }
        }
        trim();//in case introduced 0's
        return *this;
    }
    Number operator<<=(unsigned int k)
    {//first make space for extra digits
        int skipDigits = k/BASE_RADIX, skipBits = k % BASE_RADIX;
        for(int i = 0; i < skipDigits + 1; ++i) digits.append(0);
        if(skipDigits > 0)//apply whole-digit shifts
            for(int i = nDigits() - 1; i >= 0; --i)
                digits[i] = i < skipDigits ? 0 : digits[i - skipDigits];
        if(skipBits > 0)//then bit shifts
        {
            DIGIT carry = 0;
            for(int i = skipDigits; i < nDigits(); ++i)
            {
                DIGIT tempCarry = digits[i] >> (BASE_RADIX - skipBits);
                digits[i] = (digits[i] << skipBits) | carry;
                carry = tempCarry;
            }
        }
        trim();//in case bit shift not large enough to fill in all the space
        return *this;
    }

    Number& operator*=(Number const& rhs)
    {//O(|a| * |b|); multiply by one digit at a time
        Number const& a = *this, b = rhs;
        Number product(a.nDigits() + b.nDigits(), 0);
        for(int j = 0; j < b.nDigits(); ++j)
            product += mult(a, b[j]) << BASE_RADIX * j;
        product.isMinus = a.isMinus != b.isMinus;//negative if different signs
        product.trim();
        return *this = product;
    }

    static Number divide(Number const& a, Number const& b1, Number& q)
    {//O(|a| * |b|)
        assert(!b1.isZero());
        q = Number(0);
        Number b = b1.abs(), r = a.abs();//first normalize
        int norm = BASE_RADIX - lgFloor(b.digits.lastItem()) - 1;
        r <<= norm;
        b <<= norm;
        for(int i = r.nDigits() - b.nDigits(); i >= 0; --i)
        {
            int shift = i * BASE_RADIX;
            Number s = b << shift;
            DIGIT k = findK(r, s);
            q += mult(Number(1) << shift, k);//q += pk
            r -= mult(s, k);
        }//both q and r negative if different signs
        q.isMinus = r.isMinus = a.isMinus != b1.isMinus;
        return r >>= norm;//renormalize
    }
    Number& operator%=(Number const& rhs)
    {
        Number quotient(0);
        return *this = divide(*this, rhs, quotient);
    }
    Number& operator/=(Number const& rhs)
    {
        Number quotient(0);
        divide(*this, rhs, quotient);
        return *this = quotient;
    }

    string toDecimalString()const
    {
        string result;
        Number r = *this;
        while(!r.isZero())
        {
            Number q(0);
            result.push_back('0' + divide(r, Number(10), q)[0]);
            r = q;
        }
        if(isMinus) result.push_back('-');
        reverse(result.begin(), result.end());
        return result;
    }
    Number(string const& decimalS)
    {
        assert(decimalS.length() > 0);
        int firstDigit = 0;
        if(decimalS[0] == '-')//only accept n-dash for minus
        {
            assert(decimalS.length() > 1);
            isMinus = true;
            firstDigit = 1;
        }
        Number result(0);
        for(int i = firstDigit; i < decimalS.length(); ++i)
        {
            if(i == firstDigit) assert(decimalS[i] != '0');//disallow MSD = 0
            else result *= Number(10);
            assert(decimalS[i] >= '0' && decimalS[i] <= '9');
            result += Number(decimalS[i] - '0');

        }
        digits = result.digits;
    }

    int lg()const
        {return BASE_RADIX * (nDigits() - 1) + lgFloor(digits.lastItem());}

    void debug()const
    {
        DEBUG("begin");
        for(int i = nDigits() - 1; i >= 0; --i) DEBUG(digits[i]);
        DEBUG(isMinus);
        DEBUG("end");
    }
};

Number power(Number const& t, Number n)
{
    Number x = t, result(1);
    for(;;)
    {
        if(n.isOdd()) result *= x;
        n >>= 1;//cheap division by 2
        if(n.isZero()) break;
        x *= x;
    }
    return result;
}
Number modPower(Number const& t, Number n, Number const& modulus)
{
    assert(!modulus.isZero());
    Number x = t, result(1);
    for(;;)
    {
        if(n.isOdd())
        {
            result *= x;
            result %= modulus;
        }
        n >>= 1;//cheap division by 2
        if(n.isZero()) break;
        x *= x;
        x %= modulus;
    }
    return result;
}

Number sqrtInt(Number const& t)
{//start with a good guess
    Number x(Number(1) << (1 + t.lg()/2));
    for(;;)
    {
        Number y = (x + t/x)/Number(2);
        if(y < x) x = y;
        else return x;
    }
}

Number extendedGcdR(Number const& a, Number const& b, Number& x, Number& y)
{
    if(!b.isPositive())
    {
        x = Number(1);
        y = Number(0);
        return a;
    }
    Number q, r = Number::divide(a, b, q), gcd = extendedGcdR(b, r, y, x);
    y -= q * x;
    return gcd;
}
Number extendedGcd(Number const& a, Number const& b, Number& x, Number& y)
{
    assert(a.isPositive() && b.isPositive());
    return a < b ? extendedGcdR(b, a, y, x) : extendedGcdR(a, b, x, y);
}
Number gcd(Number const& a, Number const& b)
{
    Number x, y;
    return extendedGcd(a, b, x, y);
}

Number modInverse(Number const& a, Number const& n)
{
    assert(a.isPositive() && a < n);
    Number x, y;
    extendedGcd(a, n, x, y);
    if(x.isNegative()) x += n;//adjust range if needed
    return x;
}

bool provenComposite(Number const& a, Number const& n)
{
    Number ONE = Number(1), oddPart = n - ONE;
    int nSquares = 0;
    while(oddPart.isEven())
    {
        oddPart >>= 1;
        ++nSquares;
    }
    Number x = modPower(a, oddPart, n);
    for(int i = 0; i < nSquares; ++i)
    {//if x2 is 1 x must have been 1 or -1 if n is prime
        Number x2 = modPower(x, Number(2), n);
        if(x2 == ONE && x != ONE && x != n - ONE) return true;
        x = x2;
    }
    return x != ONE;
}

bool isPrime(Number const& n)
{
    n.debug();
    if(n.isEven() || n < Number(2)) return false;
    int smallPrimes[] = {3,5,7,11,13,17,19,23,29,31,37,41,43,47};
    for(int i = 0; i < sizeof(smallPrimes)/sizeof(int); ++i)
    {
        Number p = Number(smallPrimes[i]);
        if(n == p) return true;
        if((n % p).isZero()) return false;
    }//Miller-Rabin if trial division was inconclusive
    int nTrials = 1;
    int sizes[] = {73,105,132,198,223,242,253,265,335,480,543,627,747,927,
        1233,1854,4096}, nTests[] = {47,42,35,29,23,20,18,17,16,12,8,7,6,5,4,
        3,2};
    for(int i = 0; i < sizeof(sizes)/sizeof(*sizes); ++i)
        if(n.lg() < sizes[i])
        {
            nTrials = nTests[i];
            break;
        }
    while(nTrials--)
    {//use single-digit exponents for efficiency
        Number::DIGIT_TYPE max = numeric_limits<Number::DIGIT_TYPE>::max();
        if(provenComposite(Number(GlobalRNG().inRange(2, (Number(max) < n ?
            max : int(n[0])) - 1)), n)) return false;
    }
    return true;
}

}//end namespace
#endif
