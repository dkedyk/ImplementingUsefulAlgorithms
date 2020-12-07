#ifndef IGMDK_INTERVAL_NUMBER_H
#define IGMDK_INTERVAL_NUMBER_H
#include <cfenv>
#include "../Utils/Utils.h"
namespace igmdk{
//beware below, though part of C99 standard, not supported by some compilers
//such as clang
#pragma STDC FENV_ACCESS ON
template<typename ITEM = double>
class IntervalNumber: public ArithmeticType<IntervalNumber<ITEM> >
{//all operations must restore nearest rounding for other calculations
    ITEM left, right;
public:
    IntervalNumber(ITEM rounded)
    {//epsilon addition exact; change this to ulp?
        std::fesetround(FE_DOWNWARD);
        left = rounded * (1 + numeric_limits<ITEM>::epsilon());
        std::fesetround(FE_UPWARD);
        right = rounded * (1 - numeric_limits<ITEM>::epsilon());
        std::fesetround(FE_TONEAREST);
    }
    IntervalNumber(long long numerator, long long denominator)
    {
        std::fesetround(FE_DOWNWARD);
        left = numerator/denominator;
        std::fesetround(FE_UPWARD);
        right = numerator/denominator;
        std::fesetround(FE_TONEAREST);
    }
    bool isfinite()const{return std::isfinite(left) && std::isfinite(right);}
    bool isnan()const{return std::isnan(left) && std::isnan(right);}
    IntervalNumber& operator+=(IntervalNumber const& rhs)
    {
        std::fesetround(FE_DOWNWARD);
        left += rhs.left;
        std::fesetround(FE_UPWARD);
        right += rhs.right;
        std::fesetround(FE_TONEAREST);
        return *this;
    }
    IntervalNumber operator-()const//minus exact in floating represenation
        {return IntervalNumber(-right, -left);}
    IntervalNumber& operator-=(IntervalNumber const& rhs)
        {return *this += -rhs;}
    IntervalNumber& operator*=(IntervalNumber const& rhs)
    {
        std::fesetround(FE_DOWNWARD);
        ITEM temp1 = min(left * rhs.left, left * rhs.right);
        ITEM temp2 = min(right * rhs.left, right * rhs.right);
        ITEM newLeft = min(temp1, temp2);
        std::fesetround(FE_UPWARD);
        temp1 = max(left * rhs.left, left * rhs.right);
        temp2 = max(right * rhs.left, right * rhs.right);
        ITEM newRight = max(temp1, temp2);
        std::fesetround(FE_TONEAREST);
        left = newLeft;
        right = newRight;
        return *this;
    }
    IntervalNumber& operator*=(long long a)
    {//exact constant multiplication
        std::fesetround(FE_DOWNWARD);
        left *= a;
        std::fesetround(FE_UPWARD);
        right *= a;
        std::fesetround(FE_TONEAREST);
        return *this;
    }
    bool contains(ITEM a)const{return left <= a && a <= right;}
    IntervalNumber& operator/=(IntervalNumber rhs)
    {
        if(rhs.contains(ITEM(0)))
        {
            rhs.left = -numeric_limits<ITEM>::infinity();
            rhs.right = numeric_limits<ITEM>::infinity();
        }
        else
        {
            std::fesetround(FE_DOWNWARD);
            rhs.left = ITEM(1)/rhs.right;
            std::fesetround(FE_UPWARD);
            rhs.right = ITEM(1)/rhs.left;
        }
        return (*this) *= rhs;
    }
    void debug()const
    {
        cout << left << " ";
        cout << right << " ";
        cout << endl;
    }
};

}// end namespace
#endif
