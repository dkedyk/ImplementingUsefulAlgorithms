#ifndef IGMDK_LARGE_RATIONAL_H
#define IGMDK_LARGE_RATIONAL_H
#include "LargeNumber.h"
#include "../Utils/Utils.h"
#include "../Utils/Debug.h"
#include "../Utils/Vector.h"
using namespace std;
namespace igmdk{

pair<long long, int> rationalize(double x)
{//only support the usual binary (not decimal)
    assert(numeric_limits<double>::radix == 2 &&
        numeric_limits<double>::digits <= numeric_limits<long long>::digits);
    int w = numeric_limits<double>::digits, e;
    x = frexp(x, &e);//normalize x into [0.5, 1)
    long long mantissa = ldexp(x, w);//find x^53
    return make_pair(mantissa, e - w);
}

struct Rational: public ArithmeticType<Rational>
{
    Number numerator, denominator;
    Rational(Number const& theNumerator = Number(0),
        Number const& theDenominator = Number(1)): numerator(theNumerator),
        denominator(theDenominator)
    {
        assert(!denominator.isZero());
        reduce();
    }
    Rational(double x): denominator(1), numerator(1)
    {
        pair<long long, int> mantissaExponent = rationalize(x);
        numerator = Number(mantissaExponent.first);
        int e = mantissaExponent.second;
        if(e < 0) denominator <<= -e;
        else if(e > 0) numerator <<= e;
    }
    void reduce()
    {
        Number g = gcd(numerator, denominator);
        numerator /= g;
        denominator /= g;
    }
    bool isZero()const{return numerator.isZero();}
    bool isMinus()const
        {return numerator.isNegative() != denominator.isNegative();}

    Rational operator-()const
    {
        Rational result = *this;
        result.numerator.negate();
        return result;
    }
    Rational& operator+=(Rational const& rhs)
    {
        numerator = numerator * rhs.denominator + rhs.numerator *
            denominator;
        denominator *= rhs.denominator;
        reduce();
        return *this;
    }
    Rational& operator-=(Rational const& rhs){return *this += -rhs;}

    Rational& operator*=(Rational const& rhs)
    {
        numerator *= rhs.numerator;
        denominator *= rhs.denominator;
        reduce();
        return *this;
    }
    Rational& operator/=(Rational const& rhs)
    {
        assert(!rhs.isZero());
        numerator *= rhs.denominator;
        denominator * rhs.numerator;
        reduce();
        return *this;
    }

    int lg()const{return numerator.lg() - denominator.lg();}
    Number evaluate(Number const& scale = Number(1))
        {return numerator * scale / denominator;}
};

}//end namespace
#endif
