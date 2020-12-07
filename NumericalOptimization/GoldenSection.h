#ifndef IGMDK_GOLDEN_SECTION_H
#define IGMDK_GOLDEN_SECTION_H
#include <cmath>
#include "../NumericalMethods/NumericalCommon.h"
namespace igmdk{

template<typename FUNCTION> pair<double, double> minimizeGS(
    FUNCTION const& f, double xLeft, double xRight,
    double relAbsXPrecision = defaultPrecEps)
{//don't want precision too low
    assert(isfinite(xLeft) && isfinite(xRight) && xLeft <= xRight &&
        relAbsXPrecision >= numeric_limits<double>::epsilon());
    double GR = 0.618, xMiddle = xLeft * GR + xRight * (1 - GR),
        yMiddle = f(xMiddle);
    while(isELess(xLeft, xRight, relAbsXPrecision))
    {
        bool chooseR = xRight - xMiddle > xMiddle - xLeft;
        double xNew = GR * xMiddle + (1 - GR) *
            (chooseR ? xRight : xLeft), yNew = f(xNew);
        if(yNew < yMiddle)
        {
            (chooseR ? xLeft : xRight) = xMiddle;
            xMiddle = xNew;
            yMiddle = yNew;
        }
        else (chooseR ? xRight : xLeft) = xNew;
    }
    return make_pair(xMiddle, yMiddle);
}
int roundToNearestInt(double x)
{
    return int(x > 0 ? x + 0.5 : x - 0.5);
}
template<typename FUNCTION> pair<int, double> minimizeGSDiscrete(
    FUNCTION const& f, int xLeft, int xRight)
{
    assert(isfinite(xLeft) && isfinite(xRight) && xLeft <= xRight);
    double GR = 0.618;
    int xMiddle = roundToNearestInt(xLeft * GR + xRight * (1 - GR));
    double yMiddle = f(xMiddle);
    while(xLeft < xRight)
    {
        bool chooseR = xRight - xMiddle > xMiddle - xLeft;
        int xNew = roundToNearestInt(GR * xMiddle + (1 - GR) *
            (chooseR ? xRight : xLeft));
        double yNew = xNew == xMiddle ? yMiddle : f(xNew);
        if(yNew < yMiddle)
        {
            (chooseR ? xLeft : xRight) = xMiddle;
            xMiddle = xNew;
            yMiddle = yNew;
        }
        else (chooseR ? xRight : xLeft) = xNew;
    }
    return make_pair(xMiddle, yMiddle);
}

template<typename FUNCTION> pair<double, double> unimodalMinBracket(FUNCTION
    const& f, double x0, double fx, bool twoSided, double d, int maxEvals)
{
    assert(abs(d) > 0 && isfinite(x0) && isfinite(d));// && maxEvals > 2?
    pair<double, double> best(x0, x0 + d);
    double fMin = f(x0 + d);
    maxEvals -= 2;
    if(fx < fMin && twoSided)//check decrease direction if 2-sided
    {
        d *= -1;//maximal pattern must be in the other direction
        fMin = fx;
        swap(best.first, best.second);
    }
    if(!(fx < fMin))//if 1-sided can't increase current bracket
        while(d * 2 != d && maxEvals-- > 0)//handle d = 0, inf, and NaN
        {
            d *= 2;
            double xNext = x0 + d, fNext = f(xNext);
            if(fNext < fMin)
            {//shift
                best.first = best.second;
                best.second = xNext;
                fMin = fNext;
            }
            else
            {//found 3-point pattern, form interval
                best.second = xNext;
                break;
            }
        }//ensure sorted interval
    if(best.first > best.second) swap(best.first, best.second);
    return best;
}

template<typename FUNCTION> pair<double, double> minimizeGSBracket(
    FUNCTION const& f, double x, double fx = numeric_limits<double>::
    quiet_NaN(), bool twoSided = true, double step = 0.1,
    double relAbsXPrecision = defaultPrecEps, int bracketMaxEvals = 100)
{
    if(!isfinite(x)) return make_pair(x, fx);
    if(!isfinite(fx))
    {
        fx = f(x);
        --bracketMaxEvals;
    }
    pair<double, double> bracket = unimodalMinBracket(f, x, fx, twoSided,
        step * max(1.0, abs(x)), bracketMaxEvals), result = minimizeGS(f,
        bracket.first, bracket.second, relAbsXPrecision);//ensure nondecrease
    return result.second < fx ? result : make_pair(x, fx);
}

}//end namespace igmdk
#endif
