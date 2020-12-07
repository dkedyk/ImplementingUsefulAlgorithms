#ifndef IGMDK_SPSA_H
#define IGMDK_SPSA_H
#include "../RandomNumberGeneration/Random.h"
#include "NumericalOptimization.h"
namespace igmdk{

template<typename POINT, typename FUNCTION> POINT SPSA(POINT x,
    FUNCTION const& f, int maxEvals = 10000, double initialStep = 1)
{
    POINT direction = x;
    for(int i = 0, D = x.getSize(); i < maxEvals/2; ++i)
    {
        for(int j = 0; j < D; ++j) direction[j] =
            GlobalRNG().next() % 2 ? 1 : -1;
        double step = initialStep/pow(i + 1, 0.101), temp = RMRate(i) *
            (f(x + direction * step) - f(x - direction * step))/2;
        if(!isfinite(temp)) break;
        for(int j = 0; j < D; ++j) x[j] -= temp/direction[j];
    }
    return x;
}
template<typename POINT, typename FUNCTION> pair<POINT, double> metaSPSA(
    POINT x, FUNCTION const& f, int spsaEvals = 100000, int estimateEvals =
    100, double step = pow(2, 10), double minStep = pow(2, -20))
{
    pair<POINT, double> xy(x, numeric_limits<double>::infinity());
    for(; step > minStep; step /= 2)
    {
        if(isfinite(xy.second)) x = SPSA(xy.first, f, spsaEvals, step);
        double sum = 0;
        for(int i = 0; i < estimateEvals; ++i) sum += f(x);
        if(sum/estimateEvals < xy.second)
        {
            xy.first = x;
            xy.second = sum/estimateEvals;
        }
    }
    return xy;
}

}//end namespace igmdk
#endif
